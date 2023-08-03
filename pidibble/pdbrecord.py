"""

.. module:: pdbrecord
   :synopsis: defines the PDBRecrrd class
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""

from .baserecord import BaseRecord, BaseRecordParser
from .baseparsers import StringParser

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
        input_dict=StringParser(fields,typemap,allowed=allowed_values).parse(pdbrecordline)
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

    def get_token(self,key):
        if not hasattr(self,'tokengroups'):
            return None
        values={}
        for k,tg in self.tokengroups.items():
            for kk,tl in tg.items():
                if key in tl.__dict__:
                    values[kk]=tl.__dict__[key]
        if len(values)==1:
            return list(values.values())[0]
        else:
            # print(f'Warning: token-key {key} occurs more than once in the {len(self.tokengroups)} tokengroup(s) of {self.key}')
            return values

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
            for i in range(len(self.__dict__[a])):
                pt=self.__dict__[a][i]
                toks=[x.strip() for x in pt.split(':')]
                if len(toks)!=2: # this is not a token-bearing string
                    # print('ignoring tokenstring:',toks)
                    continue
                tokkey=None
                tokname,tokvalue=[x.strip() for x in pt.split(':')]
                # print(f'Found {tokname} : {tokvalue}')
                if not tokname in tdict.keys():
                    for k,v in tdict.items():
                        if 'key' in v:
                            # print(f'comparing {tokname} to {v["key"]}')
                            if tokname==v['key']:
                                tokkey=k
                else:
                    tokkey=tokname
                if not tokkey:
                    # print(f'Ignoring token {tokname} in record {self.key}')
                    continue
                typ=typemap[tdict[tokkey]['type']]
                multiline=tdict[tokkey].get('multiline',False)
                tokvalue=typ(tokvalue)
                if multiline:
                    i+=1
                    while i<len(self.__dict__[a]) and self.__dict__[a][i]!='':
                        tokvalue+=' '+self.__dict__[a][i].strip()
                        i+=1
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
            sigparser=BaseRecordParser({'signal':espec['signal']},typemap)
            assert 'value' in espec,f'Record {self.key} has an embed spec {ename} for which no value for signal {espec["signal"]} is specified'
            if 'record_index' in espec:
                sigparser.add_fields({'record_index':espec['record_index']})
            sigparse=sigparser.parse
            terparser=BaseRecordParser({'blank':['String',[12,80]]},typemap).parse
            if type(espec['record_format'])==str:
                embedfmt=format_dict.get(espec['record_format'],{})
                assert embedfmt!={},f'Record {self.key} contains an embedded_records specification with an invalid record format [{espec["record_format"]}]'
            else:
                assert type(espec['record_format'])==dict,'Record {self.key} has an embed spec {ename} for which no format is specified'
                embedfmt=espec['record_format']
            skiplines=espec.get('skiplines',0)
            tokenize=espec.get('tokenize',{})
            header=embedfmt.get('header',{})
            if tokenize:
                self.tokens={}
                tokenparser=BaseRecordParser({'token':tokenize['from']},typemap).parse
            # if header:
            #     headerparser=
            key=base_key
            embedkey=base_key
            idx=-1
            lskip=0
            triggered=False
            capturing=False
            for record in self.__dict__[embedfrom]:
                # check for signal

                sigrec=sigparse(record)
                if not triggered and sigrec.signal==espec['value']:
                    # this is a signal-line
                    triggered=True
                    if not skiplines and not tokenize:
                        capturing=True
                    embedkey=f'{base_key}.{ename}'
                    if hasattr(sigrec,'record_index'):
                        idx=sigrec.record_index
                        embedkey=f'{base_key}.{ename}.{idx}' # this separates biomolecules
                    # print(f'initiating capture for key {key} using {embedkey}')
                elif triggered and not capturing:
                    if skiplines:
                        # print(f'Skipping {record}')
                        lskip+=1
                        if lskip==skiplines:
                            capturing=True
                    elif tokenize:
                        tokenrec=tokenparser(record)
                        tokenstr=tokenrec.token
                        # print(f'Checking non-data line for tokenization "{tokenstr}"')
                        if tokenize['d'] in tokenstr:
                            k,v=tokenstr.split(tokenize['d'])
                            self.tokens[k]=v
                        else:
                            capturing=True
                elif capturing:
                    if(terparser(record).blank==''):
                        # print(f'Terminate embed capture for {embedkey} from record  {record}')
                        break # blank line ends the subrecord
                    # print(f'Parsing "{record}"')
                    # we may hit a token line in the middle of the table
                    if tokenize:
                        tokenrec=tokenparser(record)
                        tokenstr=tokenrec.token
                        # print(f'Checking non-data line for tokenization "{tokenstr}"')
                        if tokenize['d'] in tokenstr:
                            k,v=tokenstr.split(tokenize['d'])
                            # print(f'while capturing encountered token line: {k} {v}')
                            if k in self.tokens:
                                if type(self.tokens[k])!=list:
                                    self.tokens[k]=[self.tokens[k]]
                                self.tokens[k].append(v)
                            else:
                                self.tokens[k]=v
                            # we are done with this record
                            # print(f'We are done with {record}')
                            continue
                    new_record=PDBRecord.newrecord(embedkey,record,embedfmt,typemap)
                    key=new_record.key
                    record_format=new_record.format
                    # print(f'new record has key {key}')
                    if not key in new_records:
                # print(f'new record for {key}')
                        new_records[key]=new_record
                    else:
                # this must be a continuation record
                # print(f'continuing record for {key}')
                        root_record=new_records[key]
                        root_record.continue_record(new_record,record_format)
                # else:
                    # print(f'Ingoring {record}')
                        
        # print(f'embed rec new keys',new_records)
        return new_records

    def parse_tables(self,typemap):
        fmt=self.format
        self.tables={}
        scanbegin=0
        for tname,table in fmt['tables'].items():
            # print(f'{key} will acquire a table {tname} from line {scanbegin}')
            sigparser=BaseRecordParser({'signal':table['signal']},typemap).parse
            sigval=table['value']
            skiplines=table.get('skiplines',0)
            rowparser=BaseRecordParser(table['fields'],typemap).parse
            self.tables[tname]=[]
            scanfield=table['from']
            triggered=False
            capturing=False
            lskip=0
            for i in range(scanbegin,len(self.__dict__[scanfield])):
                # check for signal
                l=self.__dict__[scanfield][i]
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
                        self.tables[tname].append(parsedrow)
