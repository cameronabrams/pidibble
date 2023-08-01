"""

.. module:: baserecord
   :synopsis: defines the BaseRecord class
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""

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
    def pstr(self):
        """pstr generates a pretty string for this record
        """
        retstr=f'{self.key}:\n'
        for k,v in self.__dict__.items():
            if k!='format':
                retstr+=f'{k:>20s}: {v}'
        return retstr