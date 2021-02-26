# -*- coding: utf-8 -*-
"""
===============================================================================
    Copyright © 2021  Kouadio K.Laurent
    
    This file is part of pyCSAMT.
    
    pyCSAMT is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    
    pyCSAMT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.
    
    You should have received a copy of the GNU Lesser General Public License
    along with pyCSAMT.  If not, see <https://www.gnu.org/licenses/>.

===============================================================================  

    .Metaclasses of geostrata . Takes conventional strata geocodes 
    look in _geocodes folder . 
    
Created on Wed Nov 25 16:45:33 2020

@author: KouaoLaurent K. alias @Daniel03
"""

from csamtpy.ff.processing.callffunc import set_stratum_on_dict as setdico
try : 
    from csamtpy.utils._csamtpylog import csamtpylog
    _logger=csamtpylog.get_csamtpy_logger(__name__)
except :
    pass


traces_classes={}

class MetaStrata(type):
    "control the list of classes "
    def __init__(cls, name, bases, dict):
        type.__init__(cls, name, bases, dict)
        traces_classes[name]=cls
        
    def __repr__(self):
        return '<{0}.{1}: {2}, {3}>'.format(type(self).__name__, self.name,
                                            self.value, self.traces_classes)
        
class Strata(metaclass=MetaStrata):
    def __new__(cls, name, bases, dicstrata):
        dicstrata=setdico()[0]
        for keys , value in dicstrata.items():
            setattr(cls, keys, value)
    
class SedimentDeposits(metaclass=Strata) :
    """
    Parents of all metasediments classes
    """
    pass 

############### PROPRIETY META ################
    
class argillite(metaclass=MetaStrata):
    def __new__(cls, name, bases, dicstratum):
        dicstratum=Strata.__dict__[cls.__name__]
        for keys, values  in dicstratum.items():
            if values =='':values=None
            setattr(cls,keys,values )
        return cls.__dict__[keys]

class AGLT(metaclass=argillite):
    pass

class alluvium(metaclass=MetaStrata):
    def __new__(cls, name, bases, dicstratum):
        dicstratum=Strata.__dict__[cls.__name__]
        for keys, values  in dicstratum.items():
            if values =='':values=None
            setattr(cls,keys,values )
        return cls.__dict__[keys] 
    
class ALUV(metaclass=alluvium):
    pass

class amphibolite(metaclass=MetaStrata):
    def __new__(cls, name, bases, dicstratum):
        dicstratum=Strata.__dict__[cls.__name__]
        for keys, values  in dicstratum.items():
            if values =='':values=None
            setattr(cls,keys,values )
        return cls.__dict__[keys] 
    
class AMP(metaclass=amphibolite):
    pass 

class anorthosite(metaclass=MetaStrata):
    def __new__(cls, name, bases, dicstratum):
        dicstratum=Strata.__dict__[cls.__name__]
        for keys, values  in dicstratum.items():
            if values =='':values=None
            setattr(cls,keys,values )
        return cls.__dict__[keys] 

class ANS(metaclass=anorthosite):
    pass

class andesite(metaclass=MetaStrata):
    def __new__(cls, name, bases, dicstratum):
        dicstratum=Strata.__dict__[cls.__name__]
        for keys, values  in dicstratum.items():
            if values =='':values=None
            setattr(cls,keys,values )
        return cls.__dict__[keys]
    
class ANT(metaclass=andesite):
    pass

class aplite(metaclass=MetaStrata):
    def __new__(cls, name, bases, dicstratum):
        dicstratum=Strata.__dict__[cls.__name__]
        for keys, values  in dicstratum.items():
            if values =='':values=None
            setattr(cls,keys,values )
        return cls.__dict__[keys] 
    
class APL(metaclass=aplite):
    pass

class arkose(metaclass=MetaStrata):
    def __new__(cls, name, bases, dicstratum):
        dicstratum=Strata.__dict__[cls.__name__]
        for keys, values  in dicstratum.items():
            if values =='':values=None
            setattr(cls,keys,values )
        return cls.__dict__[keys] 
    
class ARKS(metaclass=arkose):
    pass

class arenite(metaclass=MetaStrata):
    def __new__(cls, name, bases, dicstratum):
        dicstratum=Strata.__dict__[cls.__name__]
        for keys, values  in dicstratum.items():
            if values =='':values=None
            setattr(cls,keys,values )
        return cls.__dict__[keys] 
      
class ARNT(metaclass=arenite):
    pass

class altered_rock(metaclass=MetaStrata):
    def __new__(cls, name, bases, dicstratum):
        dicstratum=Strata.__dict__[cls.__name__]
        for keys, values  in dicstratum.items():
            if values =='':values=None
            setattr(cls,keys,values )
        return cls.__dict__[keys] 
    
class ATRK(metaclass=altered_rock):
    pass

class augen_gneiss(metaclass=MetaStrata):
    def __new__(cls, name, bases, dicstratum):
        dicstratum=Strata.__dict__[cls.__name__]
        for keys, values  in dicstratum.items():
            if values =='':values=None
            setattr(cls,keys,values )
        return cls.__dict__[keys] 
    
class AUGN (metaclass=augen_gneiss):
    pass

class banded_iron_formation(metaclass=MetaStrata):
    def __new__(cls, name, bases, dicstratum):
        dicstratum=Strata.__dict__[cls.__name__]
        for keys, values  in dicstratum.items():
            if values =='':values=None
            setattr(cls,keys,values )
        return cls.__dict__[keys] 

class BIF(metaclass=banded_iron_formation):
    pass

class black_shale(metaclass=MetaStrata):
    def __new__(cls, name, bases, dicstratum):
        dicstratum=Strata.__dict__[cls.__name__]
        for keys, values  in dicstratum.items():
            if values =='':values=None
            setattr(cls,keys,values )
        return cls.__dict__[keys] 
    
class BLSH(metaclass=black_shale):
    pass

class basalt(metaclass=MetaStrata):
    def __new__(cls, name, bases, dicstratum):
        dicstratum=Strata.__dict__[cls.__name__]
        for keys, values  in dicstratum.items():
            if values =='':values=None
            setattr(cls,keys,values )
        return cls.__dict__[keys] 
    
class BLT(metaclass=basalt):
    pass

class breccia(metaclass=MetaStrata):
    def __new__(cls, name, bases, dicstratum):
        dicstratum=Strata.__dict__[cls.__name__]
        for keys, values  in dicstratum.items():
            if values =='':values=None
            setattr(cls,keys,values )
        return cls.__dict__[keys] 
    
class BX(metaclass=breccia):
    pass
class calcarenite(metaclass=MetaStrata):
    def __new__(cls, name, bases, dicstratum):
        dicstratum=Strata.__dict__[cls.__name__]
        for keys, values  in dicstratum.items():
            if values =='':values=None
            setattr(cls,keys,values )
        return cls.__dict__[keys] 
    
class CALR(metaclass=calcarenite):
    pass

class cavity(metaclass=MetaStrata):
    def __new__(cls, name, bases, dicstratum):
        dicstratum=Strata.__dict__[cls.__name__]
        for keys, values  in dicstratum.items():
            if values =='':values=None
            setattr(cls,keys,values )
        return cls.__dict__[keys] 
    
class CAV(metaclass=cavity):
    pass

class carbonate_iron_formation(metaclass=MetaStrata):
    def __new__(cls, name, bases, dicstratum):
        dicstratum=Strata.__dict__[cls.__name__]
        for keys, values  in dicstratum.items():
            if values =='':values=None
            setattr(cls,keys,values )
        return cls.__dict__[keys] 
    
class CBIF(metaclass=carbonate_iron_formation):
    pass

class carbonate_rock(metaclass=MetaStrata):
    def __new__(cls, name, bases, dicstratum):
        dicstratum=Strata.__dict__[cls.__name__]
        for keys, values  in dicstratum.items():
            if values =='':values=None
            setattr(cls,keys,values )
        return cls.__dict__[keys] 
    
class CBRK(metaclass=carbonate_rock):
    pass

class carbonatite(metaclass=MetaStrata):
    def __new__(cls, name, bases, dicstratum):
        dicstratum=Strata.__dict__[cls.__name__]
        for keys, values  in dicstratum.items():
            if values =='':values=None
            setattr(cls,keys,values )
        return cls.__dict__[keys] 
    
class CBT(metaclass=carbonatite):
    pass

class chalk(metaclass=MetaStrata):
    def __new__(cls, name, bases, dicstratum):
        dicstratum=Strata.__dict__[cls.__name__]
        for keys, values  in dicstratum.items():
            if values =='':values=None
            setattr(cls,keys,values )
        return cls.__dict__[keys] 
    
class CHLK(metaclass=chalk):
    pass

class chert(metaclass=MetaStrata):
    def __new__(cls, name, bases, dicstratum):
        dicstratum=Strata.__dict__[cls.__name__]
        for keys, values  in dicstratum.items():
            if values =='':values=None
            setattr(cls,keys,values )
        return cls.__dict__[keys] 

class CHRT(metaclass=chert):
    pass

class clast(metaclass=MetaStrata):
    def __new__(cls, name, bases, dicstratum):
        dicstratum=Strata.__dict__[cls.__name__]
        for keys, values  in dicstratum.items():
            if values =='':values=None
            setattr(cls,keys,values )
        return cls.__dict__[keys] 
    
class CLAS(metaclass=clast):
    pass

class clast_supported_breccia(metaclass=MetaStrata):
    def __new__(cls, name, bases, dicstratum):
        dicstratum=Strata.__dict__[cls.__name__]
        for keys, values  in dicstratum.items():
            if values =='':values=None
            setattr(cls,keys,values )
        return cls.__dict__[keys] 
    
class CLBX(metaclass=clast_supported_breccia):
    pass

class calcrete(metaclass=MetaStrata):
    def __new__(cls, name, bases, dicstratum):
        dicstratum=Strata.__dict__[cls.__name__]
        for keys, values  in dicstratum.items():
            if values =='':values=None
            setattr(cls,keys,values )
        return cls.__dict__[keys] 
    
class CLCR(metaclass=calcrete):
    pass

class clay(metaclass=MetaStrata):
    def __new__(cls, name, bases, dicstratum):
        dicstratum=Strata.__dict__[cls.__name__]
        for keys, values  in dicstratum.items():
            if values =='':values=None
            setattr(cls,keys,values )
        return cls.__dict__[keys] 
    
class CLY(metaclass=clay):
    pass

class conglomerate(metaclass=MetaStrata):
    def __new__(cls, name, bases, dicstratum):
        dicstratum=Strata.__dict__[cls.__name__]
        for keys, values  in dicstratum.items():
            if values =='':values=None
            setattr(cls,keys,values )
        return cls.__dict__[keys] 

class CNGL(metaclass=conglomerate):
    pass

class coal(metaclass=MetaStrata):
    def __new__(cls, name, bases, dicstratum):
        dicstratum=Strata.__dict__[cls.__name__]
        for keys, values  in dicstratum.items():
            if values =='':values=None
            setattr(cls,keys,values )
        return cls.__dict__[keys] 
    
class COAL(metaclass=coal):
    pass

class colluvium(metaclass=MetaStrata):
    def __new__(cls, name, bases, dicstratum):
        dicstratum=Strata.__dict__[cls.__name__]
        for keys, values  in dicstratum.items():
            if values =='':values=None
            setattr(cls,keys,values )
        return cls.__dict__[keys] 
      
class COLV(metaclass=colluvium):
    pass

class calc_silicate_rock(metaclass=MetaStrata):
    def __new__(cls, name, bases, dicstratum):
        dicstratum=Strata.__dict__[cls.__name__]
        for keys, values  in dicstratum.items():
            if values =='':values=None
            setattr(cls,keys,values )
        return cls.__dict__[keys] 
      
class CSRK(metaclass=calc_silicate_rock):
    pass

class carbonate_vein(metaclass=MetaStrata):
    def __new__(cls, name, bases, dicstratum):
        dicstratum=Strata.__dict__[cls.__name__]
        for keys, values  in dicstratum.items():
            if values =='':values=None
            setattr(cls,keys,values )
        return cls.__dict__[keys] 

class CVN(metaclass=carbonate_vein):
    pass



        

if __name__=='__main__':
    
    # print(amphibolite.name)
    # print(type(argillite))
    # print(Strata.amphibolite)
    print(amphibolite.color)
    # # # print(basalt.size)
    # print(dir(wood))
    # for keys, values in set_stratum_on_dict()[0].items():
    #     if keys =='argillite':
    #         print(values)
