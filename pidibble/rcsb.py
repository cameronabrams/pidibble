"""

.. module:: rcsb
   :synopsis: Manages all downloading and parsing of data from RCSB
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import urllib.request
import os
import logging
import yaml
import numpy as np
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
    'DList':list_parse(Listparser,':'),
    'LList':list_parse(Listparser,'\n')
}

class Stringparser:
    def __init__(self,fmtdict,typemap):
        self.typemap=typemap
        self.fields={k:v for k,v in fmtdict.items()}
    def parse(self,record):
        input_dict={}
        record+=' '*(80-len(record)) # pad
        for k,v in self.fields.items():
            typestring,byte_range=v
            typ=self.typemap[typestring]
            assert byte_range[1]<=len(record),f'{record} {byte_range}'
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

# def stringparse(fmtdict,typemap):
#     return Stringparser(fmtdict,typemap).parse

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

class tokengroup:
    def __init__(self,tokname,tokval,determinant=True):
        if determinant:
            self.label=f'{tokname}.{tokval}'
        else:
            self.label=f'{tokname}'
            self.add_token(tokname,tokval)
    def add_token(self,tokname,tokval):
        self.__dict__[tokname]=tokval

class PDBRecord(BaseRecord):
    continuation='0'
    @classmethod
    def base_parse(cls,current_key,pdbrecordline:str,current_format:dict,typemap:dict):
        local_record_format=current_format.copy()
        fields=local_record_format.get('fields',{})
        subrecords=local_record_format.get('subrecords',{})
        allowed_values=local_record_format.get('allowed',{})
        concats=local_record_format.get('concatenate',{})
        input_dict={}
        for k,v in fields.items():
            typestring,byte_range=v
            typ=typemap[typestring]
            assert byte_range[1]<=len(pdbrecordline)
            # columns begin at "1" not "0"
            fieldstring=pdbrecordline[byte_range[0]-1:byte_range[1]]
            input_dict[k]='' if fieldstring.strip()=='' else typ(fieldstring)
            if type(input_dict[k])==str:
                input_dict[k]=input_dict[k].strip()
            if typestring in allowed_values:
                assert input_dict[k] in allowed_values[typestring]
        for cfield,subf in concats.items():
            if not cfield in input_dict:
                input_dict[cfield]=[]
                for f in subf:
                    assert f in input_dict,f'{current_key} specifies a field for concatenation ({f}) that is not found'
                    if input_dict[f]:
                        input_dict[cfield].append(input_dict[f])
        if subrecords:
            assert 'formats' in subrecords,f'{current_key} is missing formats from its subrecords specification'
            assert 'branchon' in subrecords,f'{current_key} is missing specification of base key from its subrecords specification'
            assert subrecords['branchon'] in input_dict,f'{current_key} specifies a base record that is not found'
            required=subrecords.get('required',True)
            if required or input_dict[subrecords['branchon']] in subrecords['formats']:
                assert input_dict[subrecords['branchon']] in subrecords['formats'],f'Key "{current_key}" is missing specification of a required subrecord format for field "{subrecords["branchon"]}" value "{input_dict[subrecords["branchon"]]}" from its subrecords specification'
                subrecord_format=subrecords['formats'][input_dict[subrecords['branchon']]]
                new_key=f'{current_key}.{input_dict[subrecords["branchon"]]}'
                input_dict,current_key,current_format=PDBRecord.base_parse(new_key,pdbrecordline,subrecord_format,typemap)
        return input_dict,current_key,current_format
    
    @classmethod
    def newrecord(cls,base_key:str,pdbrecordline:str,record_format:dict,typemap:dict):
        while len(pdbrecordline)<80:
            pdbrecordline+=' '
        input_dict,current_key,current_format=PDBRecord.base_parse(base_key,pdbrecordline,record_format,typemap)
        continuation_custom_fieldname=current_format.get('continuation',None)
        if continuation_custom_fieldname:
            input_dict['continuation']=str(input_dict[continuation_custom_fieldname])
        if input_dict.get('continuation','')=='':
            input_dict['continuation']='0'
        inst=cls(input_dict)
        inst.key=current_key
        inst.format=current_format
        return inst

    def continue_record(self,other,record_format,**kwargs):
        all_fields=kwargs.get('all_fields',False)
        continuing_fields=record_format.get('continues',record_format['fields'].keys() if all_fields else {})
        # print(f'{self.key} {continuing_fields}')
        for cfield in continuing_fields:
            if type(self.__dict__[cfield])==str:
                if type(other.__dict__[cfield])==str:
                    self.__dict__[cfield]+=' '+other.__dict__[cfield]
                elif type(other.__dict__[cfield])==list:
                    self.__dict__[cfield]=[self.__dict__[cfield]]
                    self.__dict__[cfield].extend(other.__dict__[cfield])
            elif type(self.__dict__[cfield])==list:
                if type(other.__dict__[cfield])!=list:
                    assert type(self.__dict__[cfield][0])==type(other.__dict__[cfield])
                    self.__dict__[cfield].append(other.__dict__[cfield])
                else:
                    self.__dict__[cfield].extend(other.__dict__[cfield])
            else:
                self.__dict__[cfield]=[self.__dict__[cfield],other.__dict__[cfield]]

    def parse_tokens(self,typemap):
        record_format=self.format
        if not 'token_formats' in record_format:
            return
        attr_w_tokens=record_format['token_formats']
        # print(self.key,list(attr_w_tokens.keys()))
        self.tokengroups={}
        for a in attr_w_tokens.keys():
            obj=self.__dict__[a] # expect to be a list
            assert type(obj)==list
            tdict=attr_w_tokens[a]['tokens']
            determinants=attr_w_tokens[a].get('determinants',[])
            assert len(determinants) in [0,1],f'Token group for field {a} of {self.key} may not have more than one determinant'
            # print('token names',list(tdict.keys()),'determinants',determinants)
            self.tokengroups[a]={}
            current_tokengroup=None
            for pt in self.__dict__[a]:
                toks=[x.strip() for x in pt.split(':')]
                if len(toks)!=2: # this is not a token-bearing string
                    # print('ignoring tokenstring:',toks)
                    continue
                tokkey=None
                tokname,tokvalue=[x.strip() for x in pt.split(':')]
                if not tokname in tdict.keys():
                    for k,v in tdict.items():
                        if 'key' in v:
                            if tokname==v['key']:
                                tokkey=k
                else:
                    tokkey=tokname
                if not tokkey:
                    # logger.info(f'Ignoring token {tokname} in record {self.key}')
                    continue
                typ=typemap[tdict[tokkey]['type']]
                tokvalue=typ(tokvalue)
                if tokkey in determinants:
                    detrank=determinants.index(tokkey)
                    if detrank==0:
                        # print('new det tokgroup',tokkey,tokvalue)
                        new_tokengroup=tokengroup(tokkey,tokvalue)
                        self.tokengroups[a][new_tokengroup.label]=new_tokengroup
                        current_tokengroup=self.tokengroups[a][new_tokengroup.label]
                    else:
                        assert False,'should never have a detrank>0'
                        pass # should never happen
                else: # assume we are adding tokens to the last group
                    if not current_tokengroup:
                        # we have not encoutered the determinant token 
                        # so we assume there is not one
                        # print('new nondet tokgroup',tokkey,tokvalue)
                        new_tokengroup=tokengroup(tokkey,tokvalue,determinant=False)
                        self.tokengroups[a][new_tokengroup.label]=new_tokengroup
                    else:
                        current_tokengroup.add_token(tokkey,tokvalue)
    def parse_embedded(self,format_dict,typemap):
        new_records={}
        record_format=self.format 
        if not 'embedded_records' in record_format:
            return
        base_key=self.key
        embedspec=record_format.get('embedded_records',{})
        for ename,espec in embedspec.items():
            embedfrom=espec['from']
            assert embedfrom in self.__dict__,f'Record {self.key} references an invalid base field [{embedfrom}] from which to extract embeds'
            assert 'signal' in espec,f'Record {self.key} has an embed spec {ename} for which no signal is specified'
            sigparser=Stringparser({'signal':espec['signal']},typemap).parse
            assert 'value' in espec,f'Record {self.key} has an embed spec {ename} for which no value for signal {espec["signal"]} is specified'
            idxparser=None
            if 'record_index' in espec:
                idxparser=Stringparser({'record_index':espec['record_index']},typemap).parse
            if type(espec['record_format'])==str:
                embedfmt=format_dict.get(espec['record_format'],{})
                assert embedfmt!={},f'Record {self.key} contains an embedded_records specification with an invalid record format [{espec["record_format"]}]'
            else:
                assert type(espec['record_format'])==dict,'Record {self.key} has an embed spec {ename} for which no format is specified'
                embedfmt=espec['record_format']
            skiplines=espec.get('skiplines',0)
            tokenize=espec.get('tokenize',{})
            if tokenize:
                self.tokens={}
                tokenparser=Stringparser({'token':tokenize['from']},typemap).parse
            key=base_key
            embedkey=base_key
            idx=-1
            lskip=0
            triggered=False
            capturing=False
            for record in self.__dict__[embedfrom]:
                # check for signal
                if not triggered and sigparser(record).signal==espec['value']:
                    # this is a signal-line
                    triggered=True
                    if not skiplines and not tokenize:
                        capturing=True
                    embedkey=f'{base_key}.{ename}'
                    if idxparser:
                        idx=idxparser(record).record_index
                        embedkey=f'{base_key}.{ename}.{idx}'
                    print(f'initiating capture for key {key} using {embedkey}')
                elif triggered and not capturing:
                    if skiplines:
                        print(f'Skipping {record}')
                        lskip+=1
                        if lskip==skiplines:
                            capturing=True
                    elif tokenize:
                        tokenrec=tokenparser(record)
                        tokenstr=tokenrec.token
                        print(f'Checking non-data line for tokenization "{tokenstr}"')
                        if tokenize['d'] in tokenstr:
                            k,v=tokenstr.split(tokenize['d'])
                            self.tokens[k]=v
                        else:
                            capturing=True
                elif capturing:
                    if sigparser(record).signal=='':
                        print(f'Terminate embed capture for {embedkey}')
                        break # blank line ends the subrecord
                    print(f'Parsing "{record}"')
                    new_record=PDBRecord.newrecord(embedkey,record,embedfmt,typemap)
                    key=new_record.key
                    record_format=new_record.format
                    print(f'new record has key {key}')
                    if not key in new_records:
                # print(f'new record for {key}')
                        new_records[key]=new_record
                    else:
                # this must be a continuation record
                # print(f'continuing record for {key}')
                        root_record=new_records[key]
                        root_record.continue_record(new_record,record_format)
                else:
                    print(f'Ingoring {record}')
                        
        # print(f'embed rec new keys',new_records)
        return new_records

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
                self.mappers[map]=Listparser(d).parse
        # define mappers for custom formats of substrings
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
                    print(f'{key} will acquire a table {tname} from line {scanbegin}')
                    sigparser=Stringparser({'signal':table['signal']},self.mappers).parse
                    sigval=table['value']
                    skiplines=table.get('skiplines',0)
                    rowparser=Stringparser(table['fields'],self.mappers).parse
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
                                print(f'Terminate table {tname}')
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
