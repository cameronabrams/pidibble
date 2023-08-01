"""

.. module:: baserecord
   :synopsis: defines the BaseRecord class
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""

def rstr(d,excludes,pad):
    retstr=''
    kfstr=r'{:>'+str(pad)+r's}:'
    for k,v in d.items():
        if not k in excludes:
            retstr+=kfstr.format(k)
            if type(v)==dict:
                retstr+='\n'
                retstr+=rstr(v,excludes,pad+5)
            else:
                retstr+=f' {v}'+'\n'
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