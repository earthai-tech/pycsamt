# -*- coding: utf-8 -*-
"""
    Script to convert j to edi .
Created on Mon Oct 25 13:43:02 2021
"""


from pycsamt.ff.core  import CSAMT
data_fn = 'data/j' 
savepath =None 

CSAMT().j2edi(data_fn= data_fn, savepath = savepath)
