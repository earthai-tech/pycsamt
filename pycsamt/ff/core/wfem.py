# -*- coding: utf-8 -*-
#       Created on Mon Jan 11 11:37:51 2021
#       Author: Kouadio K.Laurent<etanoyau@gmail.com>
#       Licence: LGPL

"""
Created on Wed Aug 10 21:32:15 2022

@author: Daniel
"""

from pycsamt.ff.core.edi import Edi 




class WFEM (object): 
    """
    Class of wide field electromagnetic method(WFEM).
    
    The WFEM based on the pseudo-random artificial field source was developed 
    to obtain reliable resistivity information in areas with strong interferences
    The WFEM has many advantages. First, it can transmit multiple excitation 
    signals of different frequencies simultaneously, and obtain multiple frequency 
    geoelectric information at the same time, which avoids the electromagnetic 
    interference caused by random signals. Second, WFEM is simple and efficient,
    since it only  needs to measure the electric field component. Furthermore,
    with great transmission current,  WFEM has strong anti-interference ability
    and has been widely used in metal ore and shale gas exploration. The WFEM 
    can capture the low resistivity of carbonaceous shale and identify the 
    electrical distribution patterns. The WFEM  could provide information about
    strata distribution, geological structures, and fracture zones to reveal 
    the prospect of hydrocarbons, which may contribute to the shale gas exploration. 
    In addition, WFEM can be used to find some hidden gold deposits and delineate
    the favorable areas, which  has been confirmed by excavation. Thus, the 
    method is effective and feasible in exploring deep metal mineral resources.
    
    """
    def __new__(cls, data_fn, *args, **kwargs): 
        
        return super(WFEM, cls).__new__(cls, *args, **kwargs)
    
    def __init__(self, data_fn, *args, **kwargs): 
        
        self.data_fn=data_fn 
        
        