"""

.. module:: rcsb
   :synopsis: Defines the PDBParser class
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import urllib.request
import os
import logging
import yaml
import numpy as np
from mmcif.io.IoAdapterCore import IoAdapterCore
from . import resources
from .baseparsers import ListParsers, ListParser
from .baserecord import BaseRecordParser
from .pdbrecord import PDBRecord
from .mmcif_parse import MMCIF_Parser
from pidibble import resources
logger=logging.getLogger(__name__)
import importlib.metadata
import json
from .hex import str2atomSerial, hex_reset

def safe_float(x):
    if x=='nan':
        return 0.0
    return float(x)

__version__ = importlib.metadata.version("pidibble")

def str2int_sig(arg:str):
    if not arg.strip().isnumeric():
        if arg.strip()[0]=='-':
            return int(arg)
        else:
            return -1
    return int(arg)

class PDBParser:
    # mappers={'Integer':int,'String':str,'Float':float}
    mappers={'HxInteger':str2atomSerial,'Integer':str2int_sig,'String':str,'Float':safe_float}
    mappers.update(ListParsers)
    comment_lines=[]
    comment_chars=['#']
    default_input_format='PDB'
    def __init__(self,**options):
        logger.debug(f'Pidibble v. {__version__}')
        # loglevel=options.get('loglevel','info')
        # logfile=options.get('logfile','pidibble.log')
        # loglevel_numeric=getattr(logging,loglevel.upper())
        # logging.basicConfig(filename=logfile,filemode='w',format='%(asctime)s %(name)s.%(funcName)s %(levelname)s> %(message)s',level=loglevel_numeric)
        self.parsed={}
        self.input_format=options.get('input_format','PDB')
        self.pdb_code=options.get('PDBcode','')
        self.overwrite=options.get('overwrite',False)
        self.alphafold=options.get('alphafold','')
        self.pdb_format_file=options.get('pdb_format_file',os.path.join(
            os.path.dirname(resources.__file__),
            'pdb_format.yaml'))
        self.mmcif_format_file=options.get('mmcif_format_file',os.path.join(
            os.path.dirname(resources.__file__),
            'mmcif_format.yaml'))
        if os.path.exists(self.pdb_format_file):
            with open(self.pdb_format_file,'r') as f:
                self.pdb_format_dict=yaml.safe_load(f)
                logger.debug(f'Pidibble uses the installed config file:')
                logger.debug(self.pdb_format_file)
        else:
            logger.error(f'{self.pdb_format_file}: not found. You have a bad installation of pidibble.')
        if os.path.exists(self.mmcif_format_file):
            with open(self.mmcif_format_file,'r') as f:
                self.mmcif_format_dict=yaml.safe_load(f)
                logger.debug(f'Pidibble uses the installed config file:')
                logger.debug(self.mmcif_format_file)
        else:
            logger.error(f'{self.mmcif_format_file}: not found. You have a bad installation of pidibble.')
        delimiter_dict=self.pdb_format_dict.get('delimiters',{})
        for map,d in delimiter_dict.items():
            if not map in self.mappers:
                self.mappers[map]=ListParser(d).parse
        # define mappers for custom formats of substrings
        cformat_dict=self.pdb_format_dict.get('custom_formats',{})
        for cname,cformat in cformat_dict.items():
            if not cname in self.mappers:
                self.mappers[cname]=BaseRecordParser(cformat,PDBParser.mappers).parse
            
    def fetch(self):
        # TODO: allow for back-defaulting to mmCIF format
        # if PDB is not available
        assert self.pdb_code!='' or self.alphafold!=''
        if self.pdb_code!='':
            if self.input_format=='PDB':
                self.filepath=f'{self.pdb_code}.pdb'
            elif self.input_format=='mmCIF':
                self.filepath=f'{self.pdb_code}.cif'
            else:
                logger.warning(f'Input format {self.input_format} not recognized; using PDB')
                self.filepath=f'{self.pdb_code}.pdb'
            BASE_URL=self.pdb_format_dict['BASE_URL']
            target_url=os.path.join(BASE_URL,self.filepath)
            if not os.path.exists(self.filepath) or self.overwrite:
                try:
                    urllib.request.urlretrieve(target_url,self.filepath)
                except:
                    logger.warning(f'Could not fetch {self.filepath}')
                    return False
            return True
        elif self.alphafold!='':
            self.filepath=f'{self.alphafold}.pdb'
            BASE_URL=self.pdb_format_dict['ALPHAFOLD_API_URL']
            target_url=os.path.join(BASE_URL,self.alphafold)
            try:
                urllib.request.urlretrieve(target_url,f'{self.alphafold}.json')
            except:
                logger.warning(f'Could not fetch metadata for entry with accession code {self.alphafold} from AlphaFold')
                return False
            with open(f'{self.alphafold}.json') as f:
                result=json.load(f)
            try:
                urllib.request.urlretrieve(result[0]['pdbUrl'],self.filepath)
            except:
                logger.warning(f'Could not retrieve {result[0]["pdbUrl"]}')
                return False
            return True
        else:
            pass # assert statement at top of method body should suppress this branch

    def read_PDB(self):
        self.pdb_lines=[]
        with open(self.filepath,'r') as f:
            self.pdb_lines=f.read().split('\n')
            if self.pdb_lines[-1]=='':
                self.pdb_lines=self.pdb_lines[:-1]
    def read_mmCIF(self):
        io=IoAdapterCore()
        l_dc=io.readFile(self.filepath)
        self.cif_data=l_dc[0]
    def read(self):
        self.pdb_lines=[]
        self.cif_data=None
        if self.input_format=='mmCIF':
            self.read_mmCIF()
        else:
            self.read_PDB()

    def parse_base(self):
        if self.input_format=='mmCIF':
            self.parse_mmCIF()
        else:
            self.parse_PDB()

    def parse_mmCIF(self):
        mmcif_parser=MMCIF_Parser(self.mmcif_format_dict,self.pdb_format_dict['record_formats'],self.cif_data)
        self.parsed=mmcif_parser.parse()

    def parse_PDB(self):
        hex_reset()
        record_formats=self.pdb_format_dict['record_formats']
        key=''
        record_format={}
        group_open_record=None
        for i,pdbrecord_line in enumerate(self.pdb_lines):
            tc=pdbrecord_line[0]
            if tc in PDBParser.comment_chars:
                continue
            pdbrecord_line+=' '*(80-len(pdbrecord_line))
            base_key=pdbrecord_line[:6].strip()
            assert base_key in record_formats,f'{base_key} is not found in among the available record formats'
            base_record_format=record_formats[base_key]
            record_type=base_record_format['type']
            new_record=PDBRecord.newrecord(base_key,pdbrecord_line,base_record_format,self.mappers)
            key=new_record.key
            record_format=new_record.format
            if record_type in [1,2,6]:
                if not key in self.parsed:
                    self.parsed[key]=new_record
                else:
                    # this must be a continuation record
                    assert record_type!=1,f'{key} may not have continuation records'
                    root_record=self.parsed[key]
                    root_record.continue_record(new_record,record_format,all_fields=('REMARK' in key))
            elif record_type in [3,4,5]:
                if not key in self.parsed:
                    # this is necessarily the first occurance of a record with this key, but since there can be multiple instances this must be a list of records
                    if 'groupuntil' in record_format:
                        group_open_record=new_record
                        logger.debug(f'opening group {group_open_record.serial} until {group_open_record.format["groupuntil"]}')
                    if group_open_record!=None and key==group_open_record.format['groupuntil']:
                        logger.debug(f'closing group {group_open_record.serial}')
                        group_open_record=None
                    if 'groupby' in record_format:
                        tok=new_record.format['groupby'].split('.')
                        if group_open_record!=None:
                            if tok[0]==group_open_record.key:
                                groupid=getattr(group_open_record,tok[1])
                                setattr(new_record,group_open_record.key.lower(),groupid)
                    self.parsed[key]=[new_record]
                else:
                    # this is either
                    # (a) a continuation record of a given key.(determinants)
                    # or
                    # (b) a new set of (determinants) on this key
                    # note (b) is only option if there are no determinants
                    # first, look for key.(determinants)
                    root_record=None
                    if 'determinants' in record_format:
                        nrd=[new_record.__dict__[k] for k in record_format['determinants']]
                        for r in self.parsed[key]:
                            td=[r.__dict__[k] for k in record_format['determinants']]
                            if nrd==td:
                                root_record=r
                                break
                    if root_record:
                        # case (a)
                        assert root_record.continuation<new_record.continuation,f'continuation parsing error {record_type}'
                        root_record.continue_record(new_record,record_format)
                    else:
                        # case (b)
                        if 'groupuntil' in record_format:
                            group_open_record=new_record
                            logger.debug(f'opening group {group_open_record.serial} until {group_open_record.format["groupuntil"]}')
                        if group_open_record!=None and key==group_open_record.format['groupuntil']:
                            logger.debug(f'closing group {group_open_record.serial}')
                            group_open_record=None
                        if 'groupby' in record_format:
                            tok=new_record.format['groupby'].split('.')
                            if group_open_record!=None:
                                if tok[0]==group_open_record.key:
                                    groupid=getattr(group_open_record,tok[1])
                                    setattr(new_record,group_open_record.key.lower(),groupid)
                        self.parsed[key].append(new_record)

    def post_process(self):
        if self.input_format!='mmCIF':
            self.parse_embedded_records()
            self.parse_tokens()
            self.parse_tables()

    # def parse_models(self):
    #     n_models=self.parsed.get('NUMMDL',1)
    #     for i in range(n_models):
    #         self.parsed['MODEL'][i+1]={}
    #         # in progress
        
    def parse_embedded_records(self):
        new_parsed_records={}
        for key,p in self.parsed.items():
            if type(p)==PDBRecord:
                rf=p.format
                if 'embedded_records' in rf:
                    new_parsed_records.update(p.parse_embedded(self.pdb_format_dict['record_formats'],self.mappers))
            elif type(p)==list:
                for q in p:
                    rf=q.format
                    if 'embedded_records' in rf:
                        new_parsed_records.update(q.parse_embedded(self.pdb_format_dict['record_formats'],self.mappers))
        self.parsed.update(new_parsed_records)

    def parse_tokens(self):
        for key,p in self.parsed.items():
            if type(p)==PDBRecord:
                rf=p.format
                if 'token_formats' in rf:
                    p.parse_tokens(self.mappers)
            elif type(p)==list:
                for q in p:
                    rf=q.format
                    if 'token_formats' in rf:
                        q.parse_tokens(self.mappers)

    def parse_tables(self):
        for key,p in self.parsed.items():
            if type(p)==list:
                continue # don't expect to read a table from a multiple-record entry
            rf=p.format
            if 'tables' in rf:
                p.parse_tables(self.mappers)                        

    def parse(self):
        if self.fetch():
            self.read()
            self.parse_base()
            self.post_process()
        else:
            logger.warning(f'No data.')
        return self
            
def get_symm_ops(rec:PDBRecord):
    M=np.identity(3)
    T=np.array([0.,0.,0.])
    assert len(rec.row)==3,f'a transformation matrix record should not have more than 3 rows'
    for c,r in zip(rec.coordinate,rec.row):
        row=c-1
        M[row][0]=r.m1
        M[row][1]=r.m2
        M[row][2]=r.m3
        T[row]=r.t
    return M,T
