"""

.. module:: baseparsers
   :synopsis: defines some basic string and list parsing functions
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import logging
logger=logging.getLogger(__name__)

class ListParser:
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
    'CList':list_parse(ListParser,','),
    'SList':list_parse(ListParser,';'),
    'WList':list_parse(ListParser,None),
    'DList':list_parse(ListParser,':'),
    'LList':list_parse(ListParser,'\n')
}

cols="""
         1         2         3         4         5         6         7         8
12345678901234567890123456789012345678901234567890123456789012345678901234567890
"""
class StringParser:
    def __init__(self,fmtdict,typemap,allowed={}):
        self.typemap=typemap
        self.fields={k:v for k,v in fmtdict.items()}
        self.allowed=allowed
    def parse(self,record):
        if len(record)<=80:
            self.report_record_error(record)
        assert len(record)<=80,f'Record is too long; something wrong with your PDB file?'
        input_dict={}
        record+=' '*(80-len(record)) # pad
        for k,v in self.fields.items():
            typestring,byte_range=v
            typ=self.typemap[typestring]
            assert byte_range[1]<=len(record),f'{record} {byte_range}'
            # using columns beginning with "1" not "0"
            fieldstring=record[byte_range[0]-1:byte_range[1]]
            fieldstring=fieldstring.strip()
            try:
                input_dict[k]='' if fieldstring=='' else typ(fieldstring)
            except:
                self.report_field_error(record,k)
                input_dict[k]=''
            if typ==str:
                input_dict[k]=input_dict[k].strip()
            if fieldstring in self.allowed:
                assert input_dict[k] in self.allowed[fieldstring],f'Value {input_dict[k]} is not allowed for field {k}; allowed values are {self.allowed[fieldstring]}'
        return input_dict
    def report_record_error(self,record):
        logger.info(cols)
        logger.info(record)
        
    def report_field_error(self,record,k):
        logger.info(f'ERROR: Could not parse field {k}:')
        self.report_record_error(record)
        byte_range=self.fields[k][1]
        logger.info(f'{" "*byte_range[0]}{"-"*(byte_range[1]+1-byte_range[0])}')