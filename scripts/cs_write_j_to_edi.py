# -*- coding: utf-8 -*-
"""
    Script to convert j to edi .
Created on Mon Oct 25 13:43:02 2021

@author: @Daniel03
"""


from pycsamt.ff.core.j import J_collection as JObjs
path2j = 'data/j' 
savepath =None 

JObjs().j2edi(jfn=path2j, 
            savepath =savepath )