"""

.. module:: baserecord
   :synopsis: defines the BaseRecord class
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
from .baseparsers import StringParser
import logging
logger=logging.getLogger(__name__)

def rstr(d,excludes,pad):
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
    def __init__(self,input_dict):
        self.__dict__.update(input_dict)
    def empty(self):
        isempty=True
        for v in self.__dict__.values():
            isempty&=(v=='')
        return isempty
    def __str__(self):
        return '; '.join([f'{k}: {v}' for k,v in self.__dict__.items()])
    def pstr(self,excludes=['key','format','continuation'],pad=20):
        """pstr generates a pretty string for this record
        """
        retstr=f'{self.key}'+'\n'
        retstr+=rstr(self.__dict__,excludes,pad)
        return retstr

class BaseRecordParser(StringParser):
    def add_fields(self,fields):
        self.fields.update(fields)
    def parse(self,record):
        input_dict=super().parse(record)
        return BaseRecord(input_dict)