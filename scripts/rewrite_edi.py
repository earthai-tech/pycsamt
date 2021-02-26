# -*- coding: utf-8 -*-
"""
    . Script to rewrite edifiles .
Created on Tue Feb 23 13:15:59 2021

@author: @Daniel03
"""

import os 
from csamtpy.ff.core.edi import Edi


# Directory where edifiles are located 

edipath =os.path.join(os.environ['pyCSAMT'], 'data', 'edi') 

# save new edipath 
savepath =r'C:\Users\Administrator\Desktop\scripts\edinew'
       #) # change the current directory     
    
# list of Edi files 
stationlst=[os.path.join(edipath,file) for 
            file in os.listdir(edipath) if file.endswith(".edi")]

for file in stationlst: 
    Edi().write_edifile(edi_fn= file ,
                        savepath =savepath )
