# -*- coding: utf-8 -*-
"""
    Script to convert j to edi .
Created on Mon Oct 25 13:43:02 2021

@author: @Daniel03
"""


from pycsamt.ff.core  import CSAMT
path2j = 'data/j' 
savepath =None 

CSAMT().j2edi(jfn= path2j , savepath = savepath)
# JObjs().j2edi(jfn=path2j, 
#             savepath =savepath )
