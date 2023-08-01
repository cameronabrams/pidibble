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
from . import resources
from .baseparsers import StringParser, ListParsers, ListParser
from .pdbrecord import PDBRecord

logger=logging.getLogger(__name__)

class PDBParser:
    mappers={'Integer':int,'String':str,'Float':float}
    mappers.update(ListParsers)
    comment_lines=[]
    comment_chars=['#']
    def __init__(self,**options):
        self.parsed={}
        self.pdb_code=options.get('PDBcode','')
        # print(self.pdb_code)
        self.overwrite=options.get('overwrite',False)
        self.pdb_format_file=options.get('pdb_format_file',os.path.join(
            os.path.dirname(resources.__file__),
            'pdb_format.yaml'))
        if os.path.exists(self.pdb_format_file):
            with open(self.pdb_format_file,'r') as f:
                self.pdb_format_dict=yaml.safe_load(f)
        else:
            logger.fatal(f'{self.pdb_format_file}: not found.')
        delimiter_dict=self.pdb_format_dict.get('delimiters',{})
        for map,d in delimiter_dict.items():
            if not map in self.mappers:
                self.mappers[map]=ListParser(d).parse
        # define mappers for custom formats of substrings
        cformat_dict=self.pdb_format_dict.get('custom_formats',{})
        for cname,cformat in cformat_dict.items():
            if not cname in self.mappers:
                self.mappers[cname]=StringParser(cformat,PDBParser.mappers).parse
            
    def fetch(self):
        self.filename=f'{self.pdb_code}.pdb'
        BASE_URL=self.pdb_format_dict['BASE_URL']
        target_url=os.path.join(BASE_URL,self.filename)
        if not os.path.exists(self.filename) or self.overwrite:
            try:
                urllib.request.urlretrieve(target_url,self.filename)
            except:
                logger.warning(f'Could not fetch {self.filename}')

    def read(self):
        self.pdb_lines=[]
        with open(self.filename,'r') as f:
            self.pdb_lines=f.read().split('\n')
            if self.pdb_lines[-1]=='':
                self.pdb_lines=self.pdb_lines[:-1]
        
    def parse_base(self):
        record_formats=self.pdb_format_dict['record_formats']
        key=''
        record_format={}
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
                        assert root_record.continuation<new_record.continuation,f'continuation parsing error'
                        root_record.continue_record(new_record,record_format)
                    else:
                        # case (b)
                        self.parsed[key].append(new_record)


    def post_process(self):
        self.parse_embedded_records()
        self.parse_tokens()
        self.parse_tables()

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
            # print(key)
            if type(p)==PDBRecord:
                rf=p.format
                # if key=='REMARK.300':
                #     print('nonlist remark300',rf)
                if 'token_formats' in rf:
                    # print('non-list',p.key,rf)
                    p.parse_tokens(self.mappers)
            elif type(p)==list:
                for q in p:
                    rf=q.format
                    if 'token_formats' in rf:
                        # if key=='REMARK.300':
                        #     print('list remark300',rf)
                        # print('list',q.key,rf)
                        q.parse_tokens(self.mappers)

    def parse_tables(self):
        for key,rec in self.parsed.items():
            if type(rec)==list:
                continue # don't expect to read a table from a multiple-record entry
            fmt=rec.format
            if 'tables' in fmt:
                rec.tables={}
                scanbegin=0
                for tname,table in fmt['tables'].items():
                    # print(f'{key} will acquire a table {tname} from line {scanbegin}')
                    sigparser=StringParser({'signal':table['signal']},self.mappers).parse
                    sigval=table['value']
                    skiplines=table.get('skiplines',0)
                    rowparser=StringParser(table['fields'],self.mappers).parse
                    rec.tables[tname]=[]
                    scanfield=table['from']
                    triggered=False
                    capturing=False
                    lskip=0
                    for i in range(scanbegin,len(rec.__dict__[scanfield])):
                        # check for signal
                        l=rec.__dict__[scanfield][i]
                        if not triggered and sigparser(l).signal==sigval:
                            # this is a signal-line
                            triggered=True
                            if not skiplines:
                                capturing=True
                        elif triggered and not capturing:
                            if skiplines:
                                lskip+=1
                                if lskip==skiplines:
                                    capturing=True
                        elif capturing:
                            if sigparser(l).signal=='':
                                # print(f'Terminate table {tname}')
                                scanbegin=i+1
                                break
                            parsedrow=rowparser(l)
                            if not all([x=='' for x in parsedrow.__dict__.values()]):
                                rec.tables[tname].append(parsedrow)
                        

    def parse(self):
        self.fetch()
        self.read()
        self.parse_base()
        self.post_process()
        return self
            
def get_symm_ops(rec:PDBRecord):
    # assert rec.key=='REMARK.290.CRYSTSYMMTRANS'
    # TODO: infer attribute names from rec.format
    M=np.identity(3)
    T=np.array([0.,0.,0.])
    Mlist=[]
    Tlist=[]
    for r,i,c1,c2,c3,t in zip(rec.rowName,rec.replNum,rec.m1,rec.m2,rec.m3,rec.t):
        row=int(r[-1])-1
        M[row][0]=c1
        M[row][1]=c2
        M[row][2]=c3
        T[row]=t
        if row==2:
            Mlist.append(M.copy())
            Tlist.append(T.copy())
            M=np.identity(3)
    return Mlist,Tlist
