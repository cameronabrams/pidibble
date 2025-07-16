"""

.. module:: baserecord
   :synopsis: defines the BaseRecord class
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
from .baseparsers import StringParser
import logging
logger=logging.getLogger(__name__)

def rstr(d,excludes,pad):
    """
    Generate a formatted string representation of a dictionary, excluding specified keys.
    
    Parameters
    ----------
    d : dict
        The dictionary to format.
    excludes : list
        A list of keys to exclude from the output.
    pad : int
        The padding for the keys in the output string.
    
    Returns
    -------
    str
        A formatted string representation of the dictionary.
    """
    retstr=''
    kfstr=r'{:>'+str(pad)+r's}:'
    for k,v in d.items():
        if not k in excludes:
            retstr+=kfstr.format(k)
            if type(v)==dict:
                retstr+='\n'
                retstr+=rstr(v,excludes,pad+5)
            elif hasattr(v,'__len__') and not type(v)==str:
                ch=['','']
                if hasattr(v[0],'__dict__'):
                    ch=['[',']']
                retstr+=' '+', '.join([f'{ch[0]}{str(x)}{ch[1]}' for x in v])+'\n'
            else: # type(v)==str:
                retstr+=f' {str(v)}'+'\n'
    return retstr    

class BaseRecord:
    """
    A class representing a base record with fields and methods for parsing and displaying.
    """
    def __init__(self,input_dict):
        self.__dict__.update(input_dict)

    def empty(self):
        """
        Check if the record is empty, meaning all fields are empty strings.

        Returns
        -------
        bool
            True if all fields are empty strings, False otherwise.
        """
        isempty=True
        for v in self.__dict__.values():
            isempty&=(v=='')
        return isempty
    
    def __str__(self):
        """
        Generate a string representation of the BaseRecord instance.
        
        Returns
        -------
        str
            A string representation of the BaseRecord instance, showing its attributes and values.
        """
        return '; '.join([f'{k}: {v}' for k,v in self.__dict__.items()])
    
    def pstr(self,excludes=['key','format','continuation'],pad=20):
        """
        Generate a formatted string representation of the BaseRecord instance, excluding specified keys.

        Parameters
        ----------
        excludes : list, optional
            A list of keys to exclude from the output (default is ['key', 'format', 'continuation']).
        pad : int, optional
            The padding for the keys in the output string (default is 20).

        Returns
        -------
        str
            A formatted string representation of the BaseRecord instance, excluding specified keys.
        """
        retstr=f'{self.key}'+'\n'
        retstr+=rstr(self.__dict__,excludes,pad)
        return retstr

class BaseRecordParser(StringParser):
    """
    A parser for fixed-width string records that generates BaseRecord instances.  Inherits from :class:`StringParser`.
    """

    def add_fields(self,fields):
        """
        Add fields to the parser's field map.

        Parameters
        ----------
        fields : dict
            A dictionary of fields to add to the parser's field map.
        """
        self.fields.update(fields)

    def parse(self,record):
        """
        Parse a fixed-width string record into a BaseRecord instance.
        
        Parameters
        ----------
        record : str
            The fixed-width string record to parse.
            
        Returns
        -------
        BaseRecord
            A BaseRecord instance containing the parsed fields.
        """
        input_dict=super().parse(record)
        return BaseRecord(input_dict)