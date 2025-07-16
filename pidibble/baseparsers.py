"""

.. module:: baseparsers
   :synopsis: defines some basic string and list parsing functions
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import logging
logger=logging.getLogger(__name__)

class ListParser:
    """
    A simple parser for lists of strings, with a customizable delimiter.
    """

    def __init__(self,d=','):
        self.d=d
    
    def parse(self,string):
        """
        Parse a string into a list of strings, using the specified delimiter.
        If no delimiter is specified, it splits on whitespace.
        
        Parameters
        ----------
        string : str
            The string to parse.

        Returns
        -------
        list
            A list of strings parsed from the input string.
        """
        if self.d==None:
            return [x for x in string.split() if x.strip()!='']
        else:
            return [x.strip() for x in string.split(self.d) if x.strip()!='']
    
def list_parse(obj,d):
    """
    A factory function to create a ListParser with a specific delimiter.

    Parameters
    ----------
    obj : type
        The class to instantiate (should be ListParser).
    d : str or None
        The delimiter to use for parsing. If None, it will split on whitespace.

    Returns
    -------
    function
        A function that takes a string and returns a list of parsed strings.
    """
    return obj(d).parse

"""
Define a dictionary of parsers for different list formats
"""
ListParsers={
    'CList':list_parse(ListParser,','),
    'SList':list_parse(ListParser,';'),
    'WList':list_parse(ListParser,None),
    'DList':list_parse(ListParser,':'),
    'LList':list_parse(ListParser,'\n')
}

_cols="""
         1         2         3         4         5         6         7         8
12345678901234567890123456789012345678901234567890123456789012345678901234567890"""
class StringParser:
    """
    A parser for fixed-width strings, with a customizable field map.
    
    Parameters
    ----------
    fmtdict : dict
        A dictionary mapping field names to tuples of (type, byte_range).
    typemap : dict
        A dictionary mapping type names to Python types.
    allowed : dict, optional
        A dictionary mapping field values to allowed values, for validation.
    """
    def __init__(self,fmtdict,typemap,allowed={}):
        self.typemap=typemap
        self.fields={k:v for k,v in fmtdict.items()}
        self.allowed=allowed

    def parse(self,record):
        """
        Parse a fixed-width string record into a dictionary of fields.
        
        Parameters
        ----------
        record : str
            The fixed-width string record to parse.

        Returns
        -------
        dict
            A dictionary of fields parsed from the input record.
        """
        if len(record)>80:
            logger.warning('The following record exceeds 80 bytes in length:')
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
            fieldstring=fieldstring.rstrip()
            try:
                # if len(fieldstring)>0 and not typ==str:
                #     fieldstring=''
                input_dict[k]='' if fieldstring=='' else typ(fieldstring)
            except:
                self.report_field_error(record,k)
                input_dict[k]=''
            if typ==str:
                input_dict[k]=input_dict[k].strip()
            if fieldstring in self.allowed:
                assert input_dict[k] in self.allowed[fieldstring],f'Value {input_dict[k]} is not allowed for field {k}; allowed values are {self.allowed[fieldstring]}'
        return input_dict
    
    def report_record_error(self,record,byte_range=[]):
        """
        Report an error in parsing a fixed-width string record.
        
        Parameters
        ----------
        record : str
            The fixed-width string record that caused the error.
        byte_range : list, optional
            A list of byte ranges to highlight in the error message.
            If empty, the entire record is reported. 
        """
        if byte_range:
            record=record[:byte_range[0]-1]+'\033[91m'+record[byte_range[0]:byte_range[1]+1]+'\033[0m'+record[byte_range[1]+1:]
        repstr=_cols+'\n'+record
        logger.warning(repstr)
        
    def report_field_error(self,record,k):
        """
        Report an error in parsing a specific field from a fixed-width string record.
        
        Parameters
        ----------
        record : str
            The fixed-width string record that caused the error.
        k : str
            The field name that caused the error.
        """
        byte_range=self.fields[k][1]
        logger.warning(f'Could not parse field {k} from bytes {byte_range}:')
        self.report_record_error(record,byte_range=byte_range)

def safe_float(x):
    """
    Convert a string to a float, returning 0.0 if the string is 'nan'.
    """
    if x=='nan':
        return 0.0
    return float(x)

def str2int_sig(arg:str):
    """
    Convert a string to an integer, returning -1 if the string is not numeric.
    If the string starts with a '-', it is returned as an integer.
    """
    if not arg.strip().isnumeric():
        if arg.strip()[0]=='-':
            return int(arg)
        else:
            return -1
    return int(arg)