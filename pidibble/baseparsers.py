"""

.. module:: baseparsers
   :synopsis: defines some basic string and list parsing functions
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
from .baserecord import BaseRecord

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

class StringParser:
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