"""

.. module:: rcsb
   :synopsis: Manages all downloading and parsing of data from RCSB
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import urllib.request
import os
import logging
import yaml

from . import resources

logger=logging.getLogger(__name__)

BASE_URL='https://files.rcsb.org/download'

class Listparser:
    def __init__(self,d=','):
        self.d=d
    def parse(self,string):
        if self.d==None:
            return [x for x in string.split() if x.strip()!='']
        else:
            return [x.strip() for x in string.split(self.d) if x.strip()!='']
    
def list_parse(obj,d):
    return obj(d).parse

ListParsers={
    'CList':list_parse(Listparser,','),
    'SList':list_parse(Listparser,';'),
    'WList':list_parse(Listparser,None),
    'DList':list_parse(Listparser,':')
}

class Stringparser:
    def __init__(self,fmtdict,typemap):
        self.typemap=typemap
        self.fields={k:v for k,v in fmtdict.items()}
    def parse(self,record):
        input_dict={}
        total_bytes=0
        for k,v in self.fields.items():
            typestring,byte_range=v
            total_bytes+=byte_range[1]-byte_range[0]+1
        while len(record)<=total_bytes:
                record+=' '
        # print(f'({record})',len(record))
        for k,v in self.fields.items():
            typestring,byte_range=v
            typ=self.typemap[typestring]
            assert byte_range[1]<=len(record)
            # using columns beginning with "1" not "0"
            fieldstring=record[byte_range[0]-1:byte_range[1]]
            # print(k,f'({fieldstring})')
            fieldstring=fieldstring.strip()
            # print(typestring,typ)
            input_dict[k]='' if fieldstring=='' else typ(fieldstring)
            if typ==str:
                input_dict[k]=input_dict[k].strip()
        return BaseRecord(input_dict)

# def isempty(adict):
#     isempty=True
#     for v in adict.values():
#         isempty&=(v=='')
#     return isempty

def stringparse(fmtdict,typemap):
    return Stringparser(fmtdict,typemap).parse

class BaseRecord:
    def __init__(self,input_dict):
        self.__dict__.update(input_dict)
    def empty(self):
        isempty=True
        for v in self.__dict__.values():
            isempty&=(v=='')
        return isempty
    def __str__(self):
        return ';'.join([f'{k}:[{v}]' for k,v in self.__dict__.items()])

class PDBRecord(BaseRecord):
    def __init__(self,input_dict):
        for k,v in input_dict.items():
            self.__dict__[k]=v
    
    @classmethod
    def base_parse(cls,pdbrecordline:str,record_format:dict,typemap:dict):
        while len(pdbrecordline)<80:
            pdbrecordline+=' '
        input_dict={'key':pdbrecordline[0:6].strip()}
        fields=record_format.get('fields',{})
        subrecords=record_format.get('subrecords',{})
        allowed_values=record_format.get('allowed',{})
        concats=record_format.get('concatenate',{})
        for k,v in fields.items():
            typestring,byte_range=v
            typ=typemap[typestring]
            assert byte_range[1]<=len(pdbrecordline)
            # using columns beginning with "1" not "0"
            fieldstring=pdbrecordline[byte_range[0]-1:byte_range[1]]
            # print(typestring,typ)
            input_dict[k]='' if fieldstring.strip()=='' else typ(fieldstring)
            if type(input_dict[k])==str:
                input_dict[k]=input_dict[k].strip()
            # if typ==str:
            #     input_dict[k]=input_dict[k].strip()
            if typestring in allowed_values:
                assert input_dict[k] in allowed_values[typestring]
        for cfield,subf in concats.items():
            if not cfield in input_dict:
                input_dict[cfield]=[]
                for f in subf:
                    assert f in input_dict,f'{input_dict["key"]} specifies a field for concatenation ({f}) that is not found'
                    # print(f,input_dict[f])
                    if input_dict[f]:
                        input_dict[cfield].append(input_dict[f])
                    # del input_dict[f]
        if subrecords:
            assert 'formats' in subrecords,f'{input_dict["key"]} is missing formats from its subrecords specification'
            assert 'branchon' in subrecords,f'{input_dict["key"]} is missing specification of base key from its subrecords specification'
            assert subrecords['branchon'] in input_dict,f'{input_dict["key"]} specifies a base record that is not found'
            assert input_dict[subrecords['branchon']] in subrecords['formats'],f'{input_dict["key"]} is missing specification of a subrecord format for {subrecords["branchon"]} from its subrecords specification'
            subrecord_format=subrecords['formats'][input_dict[subrecords['branchon']]]
            subr_dict=PDBRecord.base_parse(cls,pdbrecordline,subrecord_format,typemap)
            input_dict.update(subr_dict)
        return input_dict
    
    @classmethod
    def newrecord(cls,pdbrecordline,record_format,typemap):
        input_dict=PDBRecord.base_parse(cls,pdbrecordline,record_format,typemap)
        return cls(input_dict)

    def update_2(self,pdbrecord,record_format,typemap):
        if not 'continuation' in record_format['fields']:
            logging.fatal(f'A type-2 pdb record must have a continuation field')
        input_dict=PDBRecord.base_parse(pdbrecord,record_format,typemap)
        for k in input_dict.keys():
            # just add to the end of the string
            if type(self.__dict__[k])==str:
                self.__dict__[k]+=' '+input_dict[k]
            else:
                if not type(self.__dict__[k])==list:
                    self.__dict__[k]=[self.__dict__[k]]
                if type(input_dict[k])==list:
                    self.__dict__[k].extend(input_dict[k])
                else:
                    self.__dict__[k].append(input_dict[k])

    def parse_tokens(self,record_format,typemap):
        if not 'tokens' in record_format:
            return
        attr_w_tokens=record_format['tokens']
        # print(attr_w_tokens)
        self.tokens={}
        current_parent_toknames={}
        for a in attr_w_tokens.keys():
            obj=self.__dict__[a] # expect to be a list
            assert type(obj)==list
            self.tokens[a]={}
            tdict=attr_w_tokens[a]
            for pt in self.__dict__[a]:
                toks=[x.strip() for x in pt.split(':')]
                assert len(toks)==2,f'Malformed token string {pt}'
                tokname,tokvalue=[x.strip() for x in pt.split(':')]
                assert tokname in tdict.keys(),f'Unrecognized token {tokname}'
                typ=typemap[tdict[tokname]['type']]
                if 'associated_to' in tdict[tokname]:
                    parent_tokname=tdict[tokname]['associated_to']
                    asstokname=current_parent_toknames[parent_tokname]
                    self.tokens[a][asstokname][tokname]=typ(tokvalue)

                else:
                    tokvalue=typ(tokvalue)
                    ntokname=f'{tokname}.{tokvalue}'
                    self.tokens[a][ntokname]={}
                    current_parent_toknames[tokname]=ntokname
    def merge(self,other,record_format):
        continues=self.__dict__ if (not 'continues' in record_format or not record_format['continues']) else record_format['continues']
        # print(continues)
        # print(self.__dict__)
        for cfield in continues:
            if type(self.__dict__[cfield])==list:
                if type(other.__dict__[cfield])==list:
                    self.__dict__[cfield].extend(other.__dict__[cfield])
                else:
                    self.__dict__[cfield].append(other.__dict__[cfield])
            elif type(self.__dict__[cfield])==str: # assume a string
                # print(self.__dict__[cfield],other.__dict__[cfield])
                self.__dict__[cfield]+=' '+other.__dict__[cfield]
            elif type(self.__dict__[cfield])==int:
                self.__dict__[cfield]=[self.__dict__[cfield],other.__dict__[cfield]]
            else:
                self.__dict__[cfield]=[self.__dict__[cfield],other.__dict__[cfield]]
        # concats=record_format.get('concatenate',{})
        # for concatfield,subfields in concats.items():
        #     if not concatfield in self.__dict__.keys():
        #         self.__dict__[concatfield]=[]
            
class PDBParser:
    previous_key=None
    previous_record_format={}
    last_parsed_3=None
    last_parsed_5=None
    # keys_encountered=[]
    parsed={}
    mappers={'Integer':int,'String':str,'Float':float}
    mappers.update(ListParsers)
    comment_lines=[]
    comment_chars=['#']
    def __init__(self,**options):
        self.pdb_code=options.get('PDBcode','')
        # print(self.pdb_code)
        self.overwrite=options.get('overwrite',False)
        pdb_format_file=options.get('pdb_format_file',os.path.join(
            os.path.dirname(resources.__file__),
            'pdb_format.yaml'))
        if os.path.exists(pdb_format_file):
            with open(pdb_format_file,'r') as f:
                self.pdb_format_dict=yaml.safe_load(f)
        delimiter_dict=self.pdb_format_dict.get('delimiters',{})
        for map,d in delimiter_dict.items():
            if not map in self.mappers:
                self.mappers[map]=Listparser(d).parse
        cformat_dict=self.pdb_format_dict.get('custom_formats',{})
        for cname,cformat in cformat_dict.items():
            if not cname in self.mappers:
                self.mappers[cname]=Stringparser(cformat,PDBParser.mappers).parse
        # print(self.mappers)
    # def get_record_format(self,key):
    #     for rf in self.pdb_format_dict['record_formats']:
    #         if rf['key']==key:
    #             return rf
    #     return {}
            
    def fetch(self):
        self.filename=f'{self.pdb_code}.pdb'
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

    def parse(self,**options):
        for i,l in enumerate(self.pdb_lines):
            key=l[:6].strip()
            if key[0] in PDBParser.comment_chars:
                self.comment_lines.append([i,l])
                continue
            assert key in self.pdb_format_dict['record_formats'],f'{key} is not found in among the available record formats'
            record_format=self.pdb_format_dict['record_formats'][key]
            # print(key,record_format)
            record_type=record_format['type']
            new_record=PDBRecord.new(l,record_format,self.mappers)

            if record_type==1: # one-time-single-line
                # print(key,record_format)
                if key in self.parsed:
                    logger.fatal(f'PDB record type {key} can appear only once in a PDB file.')
                self.parsed[key]=new_record
            elif record_type==2: # one-time-multiple-lines
                if not key in self.parsed:
                    self.parsed[key]=new_record
                else:
                    continuation=new_record.continuation
                    assert continuation.isdigit()
                    base_record=self.parsed[key]
                    base_record.extend(new_record,record_format,self.mappers)
            elif record_type in [3,4]: # multiple-times-single-line, except for some reason REVDAT, so let's make it all the same
                if not key in self.parsed:
                    self.parsed[key]=[new_record]
                else:
                    if 'continuation' in new_record.__dict__:
                        assert key=='REVDAT'
                        if new_record.continuation.isdigit():
                            base_record=None
                            nrd=[new_record.__dict__[k] for k in record_format['determinants']]
                            for r in self.parsed[key]:
                                td=[r.__dict__[k] for k in record_format['determinants']]
                                if nrd==td:
                                    base_record=r
                                    break
                            if base_record:
                                base_record.extend(new_record,record_format,self.mappers)
                            else:
                                self.parsed[key].append(new_record)
                    else:
                        self.parsed[key].append(new_record)
            elif record_type==5:
                if not key in self.parsed:
                    self.parsed[key]=[]
                    self.last_parsed_5=None
                parsed_record=PDBRecord.new(l,record_format,self.mappers,continue_on=self.last_parsed_5)
                if parsed_record!=self.last_parsed_5:
                    self.parsed[key].append(parsed_record)
                self.last_parsed_5=self.parsed[key][-1]
            # elif record_type==6:
            #     if not key in self.parsed:
            #         self.parsed[key]=[]
            #     if 'subrecords' in record_format:
            #         parsed_subrecord=PDBRecord.newsubrecord(l,record_format,self.mappers)
            #     else:
            #         parsed_record=PDBRecord.new(l,record_format,self.mappers)
            #         self.parsed[key].append(parsed_record)
            

