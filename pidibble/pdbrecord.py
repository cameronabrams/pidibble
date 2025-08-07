# Author: Cameron F. Abrams <cfa22@drexel.edu>
"""

.. module:: pdbrecord
   :synopsis: defines the PDBRecord class
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
from __future__ import annotations
from collections import UserList, UserDict
from .baserecord import BaseRecord, BaseRecordParser
from .baseparsers import StringParser
import logging

logger=logging.getLogger(__name__)

class tokengroup:
    """
    A class to represent a group of tokens with a label and a method to add tokens.
    """
    def __init__(self,tokname,tokval,determinant=True):
        if determinant:
            self.label=f'{tokname}.{tokval}'
        else:
            self.label=f'{tokname}'
            self.add_token(tokname,tokval)

    def add_token(self,tokname,tokval):
        """
        Add a token to the token group.
        
        Parameters
        ----------
        tokname : str
            The name of the token.
        tokval : str
            The value of the token.
        """
        self.__dict__[tokname]=tokval

class PDBRecord(BaseRecord):
    """
    A class representing a PDB record, inheriting from :class:`.baserecord.BaseRecord`.
    It provides methods for parsing and handling PDB records, including embedded records and tokens.
    """
    continuation='0'
    @classmethod
    def base_parse(cls,current_key,pdbrecordline:str,current_format:dict,typemap:dict):
        """
        Parse a PDB record line based on the provided format and type mapping.
            
        This method handles the parsing of the PDB record line according to the specified format.
        It extracts fields, subrecords, allowed values, and concatenated fields based on the format
        and type mapping provided. It also checks for subrecords and handles them accordingly.

        Parameters
        ----------
        current_key : str
            The key for the current record being parsed.
        pdbrecordline : str
            The line from the PDB file to be parsed.
        current_format : dict
            The format dictionary defining the structure of the PDB record.
        typemap : dict
            A dictionary mapping field names to their types.

        Returns
        -------
        tuple
            A tuple containing the parsed input dictionary, the current key, and the current format.
        """
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
        """
        Create a new PDBRecord instance from a PDB record line and its format.
        This method parses the PDB record line according to the specified format and type mapping,
        and returns a new instance of the PDBRecord class with the parsed data.
        
        Parameters
        ----------
        base_key : str
            The base key for the PDB record.
        pdbrecordline : str
            The line from the PDB file to be parsed.
        record_format : dict
            The format dictionary defining the structure of the PDB record.
        typemap : dict
            A dictionary mapping field names to their types.
        
        Returns
        -------
        :class:`PDBRecord`
            A new instance of the PDBRecord class containing the parsed data.
        """
        # logger.debug(f'newrecord pdbrecordline "{pdbrecordline}"')
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
        """
        Retrieve a token value from the PDBRecord instance based on the provided key.
        
        Parameters
        ----------
        key : str
            The key for the token to retrieve.
        
        Returns
        -------
        str or None
            The value of the token if found, or None if the token does not exist.
        """
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
            return values

    def continue_record(self,other,record_format,**kwargs):
        """
        Continue a PDBRecord instance with another PDBRecord instance.
        This method merges the attributes of the other PDBRecord instance into the current instance,
        handling continuation fields and concatenating values as necessary.
        
        Parameters
        ----------
        other : PDBRecord
            The other PDBRecord instance to merge with.
        record_format : dict
            The format dictionary defining the structure of the PDB record.
        kwargs : dict, optional
            Additional keyword arguments, such as 'all_fields' to specify whether to include all fields.
        all_fields : bool, optional
            If True, all fields from the record format will be considered for continuation.
            If False, only the fields specified in the record format will be considered.
            
        Returns
        -------
        None
            This method modifies the current instance in place.
        """
        all_fields=kwargs.get('all_fields',False)
        continuing_fields=record_format.get('continues',record_format['fields'].keys() if all_fields else {})
        logger.debug(f'{self.key} {continuing_fields}')
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
        """
        Parse tokens from the PDBRecord instance based on the record format.
        This method checks if the record format contains token formats and parses them accordingly.
        
        Parameters
        ----------
        typemap : dict
            A dictionary mapping field names to their types.
            
        Returns
        -------
        None
            This method modifies the PDBRecord instance in place, adding a `tokengroups` attribute
            that contains the parsed tokens grouped by their labels.
        """
        record_format=self.format
        if not 'token_formats' in record_format:
            return
        attr_w_tokens=record_format['token_formats']
        logger.debug(f'{self.key} {list(attr_w_tokens.keys())}')
        self.tokengroups={} # one tokengroup per attribute in attr_w_tokens
        for a in attr_w_tokens.keys():
            obj=self.__dict__[a] # expect to be a list
            assert type(obj)==list,f'Invalid type {type(obj)} for {obj} for token parsing; expecting a list of token-strings'
            tdict=attr_w_tokens[a]['tokens']
            determinants=attr_w_tokens[a].get('determinants',[])
            assert len(determinants) in [0,1],f'Token group for field {a} of {self.key} may not have more than one determinant'
            logger.debug(f'token names {list(tdict.keys())} determinants {determinants}')
            self.tokengroups[a]={}
            current_tokengroup=None
            for i in range(len(self.__dict__[a])):
                pt=self.__dict__[a][i]
                toks=[x.strip() for x in pt.split(':')]
                if len(toks)!=2: # this is not a token-bearing string
                    logger.debug(f'ignoring tokenstring: {toks}')
                    continue
                tokkey=None
                try:
                    tokname,tokvalue=[x.strip() for x in pt.split(':')]
                except:
                    logger.warning(f'Invalid format for token-string {pt}')
                    continue
                logger.debug(f'Found {tokname} : {tokvalue}')
                if not tokname in tdict.keys():
                    for k,v in tdict.items():
                        if 'key' in v:
                            logger.debug(f'comparing {tokname} to {v["key"]}')
                            if tokname==v['key']:
                                tokkey=k
                else:
                    tokkey=tokname
                if not tokkey:
                    logger.debug(f'Ignoring token {tokname} in record {self.key}')
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
                        logger.debug(f'new det tokgroup {tokkey} {tokvalue}')
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
                        logger.debug(f'new nondet tokgroup {tokkey} {tokvalue}')
                        new_tokengroup=tokengroup(tokkey,tokvalue,determinant=False)
                        self.tokengroups[a][new_tokengroup.label]=new_tokengroup
                    else:
                        current_tokengroup.add_token(tokkey,tokvalue)

    def parse_embedded(self,format_dict,typemap):
        """
        Parse embedded records within the PDBRecord instance based on the record format.
        This method checks if the record format contains embedded records and parses them accordingly.

        Parameters
        ----------
        format_dict : dict
            A dictionary mapping field names to their formats.
        typemap : dict
            A dictionary mapping field names to their types.
        """
        logger.debug(f'Parsing embedded')
        new_records={}
        record_format=self.format 
        if not 'embedded_records' in record_format:
            return
        base_key=self.key
        embedspec=record_format.get('embedded_records',{})
        for ename,espec in embedspec.items():
            logger.debug(f'Embedded {ename}')
            embedfrom=espec['from']
            assert embedfrom in self.__dict__,f'Record {self.key} references an invalid base field [{embedfrom}] from which to extract embeds'
            assert 'signal' in espec,f'Record {self.key} has an embed spec {ename} for which no signal is specified'
            sigparse=BaseRecordParser({'signal':espec['signal']},typemap).parse
            assert 'value' in espec,f'Record {self.key} has an embed spec {ename} for which no value for signal {espec["signal"]} is specified'
            idxparse=None
            if 'record_index' in espec:
                idxparse=BaseRecordParser({'record_index':espec['record_index']},typemap).parse
            terparse=BaseRecordParser({'blank':['String',[12,80]]},typemap).parse
            if type(espec['record_format'])==str:
                embedfmt=format_dict.get(espec['record_format'],{})
                assert embedfmt!={},f'Record {self.key} contains an embedded_records specification with an invalid record format [{espec["record_format"]}]'
            else:
                assert type(espec['record_format'])==dict,'Record {self.key} has an embed spec {ename} for which no format is specified'
                embedfmt=espec['record_format']
            skiplines=espec.get('skiplines',0)
            tokenize=espec.get('tokenize',{})
            headers=espec.get('headers',{})
            token_hold={}
            header_hold=[]
            if tokenize:
                tokenparser=BaseRecordParser({'token':tokenize['from']},typemap).parse
                if headers:
                    headertokenparser=BaseRecordParser({k:v['format'] for k,v in headers['formats'].items()},typemap).parse
            embedkey=base_key
            idx=-1
            lskip=0
            triggered=False # True when an embedded_record signal is encountered
            capturing=False # True once we have skipped the desired number of lines
                            # after the signal OR if we are tokenizing and encounter
                            # the FIRST non-tokenizable line
            current_division=0
            for record in self.__dict__[embedfrom]:
                # check for signal
                sigrec=sigparse(record)
                if not triggered and sigrec.signal==espec['value']:
                    idx=None if not idxparse else idxparse(record).record_index
                    # this is a signal-line
                    triggered=True
                    if not skiplines and not tokenize and not headers:
                        capturing=True
                    embedkey=f'{base_key}.{ename}'
                    if idx:
                        embedkey=f'{base_key}.{ename}{idx}'
                    # savkey=embedkey
                    continue # go to next record
                if triggered and not capturing:
                    # we can skip lines or we can gather tokens
                    if skiplines:
                        logger.debug(f'Skipping {record}')
                        lskip+=1
                        if lskip==skiplines:
                            capturing=True
                        continue
                    logger.debug(f'Parsing "{record}"')
                    if tokenize:
                        is_ht=header_or_token(record,tokenize['d'],headers,tokenparser,headertokenparser,token_hold,header_hold)
                        if is_ht:
                            continue # go to next record
                        else:
                            # if it not a token, it must be a record described by record_format
                            capturing=True
                            new_div=capture_record(record,embedfmt,typemap,embedkey,headers,header_hold,token_hold,current_division,new_records)
                            if new_div:
                                current_division+=1
                                logger.debug(f'First capture into division {current_division}')
                            continue
                elif capturing:
                    # if we are capturing, the first occurrence of a blank line
                    # terminates the search for embedded records
                    if(terparse(record).blank==''):
                        logger.debug(f'Terminate embed capture for {embedkey} from record {record}')
                        break # finished!
                    logger.debug(f'Parsing "{record}"')
                    # capturing can capture embedded records or tokens
                    if tokenize:
                        is_ht=header_or_token(record,tokenize['d'],headers,tokenparser,headertokenparser,token_hold,header_hold)
                        if is_ht:
                            continue
                    new_div=capture_record(record,embedfmt,typemap,embedkey,headers,header_hold,token_hold,current_division,new_records)
                    if new_div:
                        current_division+=1
                    continue
                else:
                    logger.debug(f'Ingoring {record}')
                        
        logger.debug(f'embed rec new keys {new_records}')
        return new_records

    def parse_tables(self,typemap):
        """
        Parse tables from the PDBRecord instance based on the record format.
        This method checks if the record format contains table formats and parses them accordingly.
        
        Parameters
        ----------
        typemap : dict
            A dictionary mapping field names to their types.
            
        Returns
        -------
        None. This method modifies the PDBRecord instance in place, adding a `tables` attribute
        that contains the parsed tables, where each table is a list of :class:`.pdbrecord.PDBRecord` instances.
        """
        fmt=self.format
        self.tables={}
        scanbegin=0
        for tname,table in fmt['tables'].items():
            logger.debug(f'{self.key} will acquire a table {tname} from line {scanbegin}')
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
                        logger.debug(f'Terminate table {tname}')
                        scanbegin=i+1
                        break
                    parsedrow=rowparser(l)
                    if not all([x=='' for x in parsedrow.__dict__.values()]):
                        self.tables[tname].append(parsedrow)

class PDBRecordList(UserList):
    """
    A class representing a list of PDBRecord instances, inheriting from UserList.
    It provides methods for parsing and handling multiple PDB records.
    """
    def __init__(self, initlist=None):
        if initlist is not None:
            self._validate_all(initlist)
        super().__init__(initlist or [])

    def _validate(self, item):
        if not isinstance(item, PDBRecord):
            raise TypeError(f"All items must be instances of PDBRecord, got {type(item)}")

    def _validate_all(self, iterable):
        for item in iterable:
            self._validate(item)

    def __setitem__(self, index, item):
        # Support slice assignment
        if isinstance(index, slice):
            self._validate_all(item)
        else:
            self._validate(item)
        super().__setitem__(index, item)

    def append(self, item):
        self._validate(item)
        super().append(item)

    def insert(self, index, item):
        self._validate(item)
        super().insert(index, item)

    def extend(self, other):
        self._validate_all(other)
        super().extend(other)

    def __add__(self, other):
        self._validate_all(other)
        return PDBRecordList(super().__add__(other))

    def __iadd__(self, other):
        self._validate_all(other)
        return super().__iadd__(other)
    
class PDBRecordDict(UserDict):
    """
    A class representing a dictionary of PDBRecord or PDBRecordList instances, inheriting from UserDict.
    It provides methods for parsing and handling multiple PDB records stored in a dictionary.
    """
    def __init__(self, *args, **kwargs):
        super().__init__()
        self.update(*args, **kwargs)

    def _validate(self, value):
        if not isinstance(value, (PDBRecord, PDBRecordList)):
            raise TypeError(f"Values must be PDBRecord or PDBRecordList, got {type(value)}")

    def __setitem__(self, key, value):
        self._validate(value)
        super().__setitem__(key, value)

    def update(self, *args, **kwargs):
        other = dict(*args, **kwargs)
        for key, value in other.items():
            self[key] = value  # Triggers __setitem__

    def setdefault(self, key, default=None):
        if key not in self:
            self[key] = default  # Triggers __setitem__
        return self[key]

def header_check(record,headers,parse,hold=[]):
    """
    Check if a record is a header line and parse it accordingly.
    
    Parameters
    ----------
    record : str
        The record line to check.
    headers : dict
        A dictionary containing header formats and their specifications.
    parse : function
        A function to parse the header line.
    hold : list, optional
        A list to hold the parsed header values. Defaults to an empty list.
    """
    r=parse(record)
    if r.mainline==headers['formats']['mainline']['signalvalue']:
        assert len(hold)==0
        hold.extend([x.strip() for x in r.value.strip().split(',')])
        if '' in hold:
            hold.remove('')
    elif r.andline==headers['formats']['andline']['signalvalue']:
        hold.extend([x.strip() for x in r.value.strip().split(',')])

def gather_token(k,v,hold={}):
    """
    Gather a token into a holder dictionary.
    If the key already exists in the holder, it appends the value to the list.
    
    Parameters
    ----------
    k : str
        The key for the token.
    v : str
        The value of the token.
    hold : dict, optional
        A dictionary to hold the tokens. Defaults to an empty dictionary.
    """
    if k in hold:
        if not type(hold[k])==list:
            hold[k]=[hold[k],v]
        else:
            hold[k].append(v)
    else:
        hold[k]=v

def header_or_token(rec,d,hdrs,tp,htp,th,hh):
    """
    Check if a record is a token or a header line and parse it accordingly.
    This function checks if the record can be tokenized or if it matches a header format.
    If it matches a header format, it updates the header holder. If it is a token,
    it gathers the token into the token holder.
    
    Parameters
    ----------
    rec : str
        The record line to check.
    d : str
        The delimiter used to separate tokens in the record.
    hdrs : dict
        A dictionary containing header formats and their specifications.
    tp : function
        A function to parse the token line.
    htp : function
        A function to parse the header line.
    th : dict
        A dictionary to hold the tokens.
    hh : list
        A list to hold the header values.
    
    Returns
    -------
    bool
        True if the record was parsed as a token or header, False otherwise.
    """
    tokenstr=tp(rec).token
    if d in tokenstr:
        k,v=tokenstr.split(d)
        # check to see if this a special "header" line
        header_check(rec,hdrs,htp,hh)
        logger.debug(f'header_hold {hh}')
        if not hh:
            gather_token(k,v,th)
        return True
    return False

def capture_record(rec,fmt,typemap,key,hdrs,hh,th,divno,rh):
    """
    Capture a record from the PDB file and create a new PDBRecord instance.
    This function checks if the record is a continuation of an existing record or a new record.
    If it is a continuation, it updates the existing record. If it is a new record,
    it creates a new PDBRecord instance and adds it to the record holder.
    
    Parameters
    ----------
    rec : str
        The record line to capture.
    fmt : dict
        The format dictionary defining the structure of the PDB record.
    typemap : dict
        A dictionary mapping field names to their types.
    key : str
        The key for the current record being captured.
    hdrs : dict
        A dictionary containing header formats and their specifications.
    hh : list
        A list to hold the header values.
    th : dict
        A dictionary to hold the tokens.
    divno : int
        The current division number.
    rh : dict
        A dictionary to hold the records, where keys are record keys and values are PDBRecord instances.    
    
    Returns
    -------
    bool
        True if a new division was detected, False otherwise.
    """
    new_division=False
    embedkey=key
    if hh:
        divno+=1
        logger.debug(f'Capture inherits a header hold; currdivno {divno}')
        new_division=True
    # if we are not holding headers, we still could encounter a new division
    # of the data if the record's divnumber is not the current divnumber
    if hdrs:
        embedkey=f'{key}.{hdrs["divlabel"]}{divno}'
    logger.debug(f'record to {embedkey}')
    new_record=PDBRecord.newrecord(embedkey,rec,fmt,typemap)
    if hasattr(new_record,'divnumber'):
        mydivno=new_record.divnumber
        if mydivno==divno+1:
            logger.debug(f'New division detected {divno} -> {mydivno}')
            new_division=True
        divlabel=hdrs.get('divlabel','')
        embedkey=f'{key}.{divlabel}{mydivno}'
        new_record.key=embedkey
    thiskey=new_record.key
    record_format=new_record.format
    if hh:
        new_record.header=hh.copy()
        while hh:
            hh.pop(0)
    if th:
        new_record.tokens=th.copy()
        keys=list(th.keys())
        for k in keys:
            del th[k]
    logger.debug(f'new record has key {thiskey}')
    if not thiskey in rh:
        logger.debug(f'new record for {thiskey}')
        rh[thiskey]=new_record
    else:
        logger.debug(f'continuing record for {thiskey}')
        root_record=rh[thiskey]
        root_record.continue_record(new_record,record_format)
        if hasattr(new_record,'tokens'):
            if hasattr(root_record,'tokens'):
                root_record.tokens.update(new_record.tokens)
            else:
                root_record.tokens=new_record.tokens
    return new_division
