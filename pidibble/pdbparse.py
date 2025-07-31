# Author: Cameron F. Abrams <cfa22@drexel.edu>
"""

.. module:: pdbparse
   :synopsis: Defines the PDBParser class
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import urllib.request
import os
import logging
import yaml
import numpy as np
from pathlib import Path
from mmcif.io.IoAdapterCore import IoAdapterCore
from typing import List, Dict
from . import resources
from .baseparsers import ListParsers, ListParser, str2int_sig, safe_float
from .baserecord import BaseRecordParser
from .pdbrecord import PDBRecord, PDBRecordDict, PDBRecordList
from .mmcif_parse import MMCIF_Parser
from pidibble import resources
logger=logging.getLogger(__name__)
import importlib.metadata
import json
from .hex import str2atomSerial, hex_reset

__version__ = importlib.metadata.version("pidibble")

class PDBParser:
    """
    A class for parsing PDB files and extracting structured data.
    This class handles fetching PDB files, reading them, and parsing their contents into structured records
    based on predefined formats.

    Attributes
    ----------

    parsed : PDBRecordDict
        A dictionary containing parsed records, where keys are record types and values are :class:`.pdbrecord.PDBRecord` instances or lists of instances.
        This dictionary is populated after parsing the PDB or mmCIF file.

    mappers : dict
        A dictionary of mappers for parsing different data types, including custom formats and delimiters.

    pdb_lines : list
        A list of lines read from the PDB file. Empty if no input file is provided.

    cif_data : dict
        A dictionary containing the parsed mmCIF data. Empty if no input file is provided.
    """
    def __init__(self,
                 input_format='PDB',
                 overwrite=False,
                 alphafold='',
                 filepath: Path = None,
                 mappers={'HxInteger':str2atomSerial,'Integer':str2int_sig,'String':str,'Float':safe_float},
                 pdbcode_synonyms=['PDBcode','pdb_code','pdbCode','pdbcode','pdbID','pdb_id','pdbid','PDBID','PDBId','PDBid'],
                 comment_chars=['#'],
                 pdb_format_file='pdb_format.yaml',
                 mmcif_format_file='mmcif_format.yaml',
                 **kwargs):
        logger.debug(f'Pidibble v. {__version__}')
        self.input_format=input_format
        self.overwrite=overwrite
        self.alphafold=alphafold
        self.filepath=filepath
        self.mappers=mappers
        self.mappers.update(ListParsers)
        self.pdb_code=''
        for k in pdbcode_synonyms:
            if k in kwargs:
                self.pdb_code=kwargs[k]
                break
        self.comment_chars=comment_chars
        self.pdb_lines=[]
        self.cif_data={}

        self.parsed=PDBRecordDict()
        self.pdb_format_file=pdb_format_file
        if not os.path.isfile(self.pdb_format_file):
            # if pdb_format_file is not a file in the CWD, assume it is a relative path to the resources directory
            # this is useful for testing
            self.pdb_format_file=os.path.join(os.path.dirname(resources.__file__),pdb_format_file)
        self.mmcif_format_file=mmcif_format_file
        if not os.path.isfile(self.mmcif_format_file):
            # if mmcif_format_file is not a file in the CWD, assume it is a relative path to the resources directory
            # this is useful for testing
            self.mmcif_format_file=os.path.join(os.path.dirname(resources.__file__),mmcif_format_file)
        if os.path.exists(self.pdb_format_file):
            with open(self.pdb_format_file,'r') as f:
                self.pdb_format_dict=yaml.safe_load(f)
                logger.debug(f'Pidibble uses the installed config file {self.pdb_format_file}')
        else:
            raise FileNotFoundError(f'{self.pdb_format_file} not found, either locally ({os.getcwd()}) or in resources ({os.path.dirname(resources.__file__)})')
        if os.path.exists(self.mmcif_format_file):
            with open(self.mmcif_format_file,'r') as f:
                self.mmcif_format_dict=yaml.safe_load(f)
                logger.debug(f'Pidibble uses the installed config file {self.mmcif_format_file}')
        else:
            raise FileNotFoundError(f'{self.mmcif_format_file} not found, either locally ({os.getcwd()}) or in resources ({os.path.dirname(resources.__file__)})')

        # update mappers with delimiters and custom formats
        delimiter_dict=self.pdb_format_dict.get('delimiters',{})
        for map,d in delimiter_dict.items():
            if not map in self.mappers:
                self.mappers[map]=ListParser(d).parse
        cformat_dict=self.pdb_format_dict.get('custom_formats',{})
        for cname,cformat in cformat_dict.items():
            if not cname in self.mappers:
                self.mappers[cname]=BaseRecordParser(cformat,self.mappers).parse
            
    def fetch(self):
        """
        Fetch the PDB file based on the provided PDB code or AlphaFold ID.
        This method checks if the PDB code or AlphaFold ID is provided, constructs the appropriate file path,
        and attempts to download the file from the PDB or AlphaFold API.
        
        Returns
        -------
        bool
            True if the file was successfully fetched, False otherwise.
        """
        assert self.pdb_code!='' or self.alphafold!='' or self.filepath is not None, "PDB code or AlphaFold ID must be provided, or a path to a PDB file must be specified."
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
        elif self.filepath is not None:
            if not os.path.exists(self.filepath):
                logger.warning(f'File {self.filepath} does not exist.')
                return False
            return True
        else:
            pass # assert statement at top of method body should suppress this branch

    def read_PDB(self):
        """
        Read the PDB file and store its lines in :attr:`PDBParser.pdb_lines`.
        This method opens the PDB file, reads its contents, and splits it into lines.
        If the last line is empty, it removes it from the list of lines.
        """
        with open(self.filepath,'r') as f:
            self.pdb_lines=f.read().split('\n')
            if self.pdb_lines[-1]=='':
                self.pdb_lines=self.pdb_lines[:-1]

    def read_mmCIF(self):
        """
        Read the mmCIF file and store its data in :attr:`PDBParser.cif_data`.
        This method uses the :class:`mmcif.io.IoAdapterCore.IoAdapterCore` to read the
        mmCIF file and store the data in :attr:`PDBParser.cif_data`.
        """
        io=IoAdapterCore()
        l_dc=io.readFile(self.filepath)
        self.cif_data=l_dc[0]

    def read(self):
        """
        Read the PDB or mmCIF file based on the input format.
        This method checks the input format and calls the appropriate read method.
        """
        if self.input_format=='mmCIF':
            self.read_mmCIF()
        else:
            self.read_PDB()

    def parse_base(self):
        """
        Parse the base records from the PDB or mmCIF file.
        This method initializes the parsing process based on the input format.
        If the input format is mmCIF, it uses the :class:`.mmcif_parse.MMCIF_Parser` to parse the mmCIF data.
        If the input format is PDB, it uses the :class:`.pdbrecord.PDBRecord` class to parse the PDB lines.
        """
        if self.input_format=='mmCIF':
            self.parse_mmCIF()
        else:
            self.parse_PDB()

    def parse_mmCIF(self):
        """
        Parse the mmCIF data and generate a dictionary of :class:`.pdbrecord.PDBRecord` instances.
        This method uses the :class:`.mmcif_parse.MMCIF_Parser` to parse the mmCIF data and store the parsed records
        in :attr:`PDBParser.parsed`.
        """
        mmcif_parser=MMCIF_Parser(self.mmcif_format_dict,self.pdb_format_dict['record_formats'],self.cif_data)
        self.parsed=mmcif_parser.parse()

    def parse_PDB(self):
        """
        Parse the PDB lines and generate a dictionary of :class:`.pdbrecord.PDBRecord` instances.
        This method iterates through the PDB lines, identifies the record type based on the first character,
        and creates a new :class:`.pdbrecord.PDBRecord` instance for each record.
        It handles different record types, including continuation records and grouped records.
        """
        hex_reset()
        record_formats=self.pdb_format_dict['record_formats']
        key=''
        record_format={}
        group_open_record=None
        for i,pdbrecord_line in enumerate(self.pdb_lines):
            tc=pdbrecord_line[0]
            if tc in self.comment_chars:
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
                    self.parsed[key]=PDBRecordList([new_record])
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
        """
        Post-process the parsed records to handle embedded records, tokens, and tables.
        This method checks if the input format is mmCIF and processes the records accordingly.
        If the input format is PDB, it processes the records to handle embedded records, tokens, and tables.
        """
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
        """
        Parse embedded records within the parsed records.
        This method iterates through the parsed records and checks if any record has embedded records.
        If an embedded record is found, it calls the :meth:`.pdbrecord.PDBRecord.parse_embedded` method to parse the embedded records.
        It updates the :attr:`PDBParser.parsed` dictionary with the new parsed records.
        """
        new_parsed_records={}
        for key,p in self.parsed.items():
            if type(p)==PDBRecord:
                rf=p.format
                if 'embedded_records' in rf:
                    new_parsed_records.update(p.parse_embedded(self.pdb_format_dict['record_formats'],self.mappers))
            elif type(p)==PDBRecordList:
                for q in p:
                    rf=q.format
                    if 'embedded_records' in rf:
                        new_parsed_records.update(q.parse_embedded(self.pdb_format_dict['record_formats'],self.mappers))
        self.parsed.update(new_parsed_records)

    def parse_tokens(self):
        """
        Parse tokens within the parsed records.
        This method iterates through the parsed records and checks if any record has token formats.
        If a token format is found, it calls the :meth:`.pdbrecord.PDBRecord.parse_tokens` method to parse the tokens.
        It updates the :attr:`PDBParser.parsed` dictionary with the new parsed records.
        """
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
        """
        Parse tables within the parsed records.
        This method iterates through the parsed records and checks if any record has table formats.
        If a table format is found, it calls the :meth:`.pdbrecord.PDBRecord.parse_tables` method to parse the tables.
        It updates the :attr:`PDBParser.parsed` dictionary with the new parsed records.
        """
        for key,p in self.parsed.items():
            if type(p)==PDBRecordList:
                continue # don't expect to read a table from a multiple-record entry
            rf=p.format
            if 'tables' in rf:
                p.parse_tables(self.mappers)                        

    def parse(self):
        """
        Parse the PDB or mmCIF file and generate a dictionary of :class:`.pdbrecord.PDBRecord` instances.
        This method first fetches the PDB or mmCIF file based on the provided PDB code or AlphaFold ID.
        It then reads the file and parses its contents into structured records.
        If the input format is mmCIF, it uses the :class:`.mmcif_parse.MMCIF_Parser` to parse the mmCIF data.
        If the input format is PDB, it uses the :class:`.pdbrecord.PDBRecord` class to parse the PDB lines.

        Returns
        -------
        self : PDBParser
            The instance of :class:`.pdbrecord.PDBRecord` containing the parsed records.
        """
        if self.fetch():
            self.read()
            self.parse_base()
            self.post_process()
        else:
            logger.warning(f'No data.')
        return self
            
def get_symm_ops(rec:PDBRecord):
    """
    Extract the symmetry operations from a PDB record.
    This function processes the symmetry operations from a PDB record and returns the transformation matrix and translation vector.

    Parameters
    ----------
    rec : :class:`.pdbrecord.PDBRecord`
        The PDBRecord instance containing the symmetry operations.

    Returns
    -------
    M : :class:`numpy.ndarray`
        The 3x3 transformation matrix.
    T : :class:`numpy.ndarray`
        The 3x1 translation vector.
    """
    M = np.identity(3)
    T = np.array([0., 0., 0.])
    if not (hasattr(rec, 'row') and hasattr(rec, 'coordinate')):
        raise ValueError('Invalid PDBRecord: missing row or coordinate attributes')
    assert len(rec.row) == 3, f'a transformation matrix record should not have more than 3 rows'
    for c, r in zip(rec.coordinate, rec.row):
        row = c - 1
        M[row][0] = r.m1
        M[row][1] = r.m2
        M[row][2] = r.m3
        T[row] = r.t
    return M, T
