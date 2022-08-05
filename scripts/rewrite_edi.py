# -*- coding: utf-8 -*-
"""
    . Script to rewrite edifiles .
Created on Tue Feb 23 13:15:59 2021

@author:K.L ~ @Daniel03

"""

import os 
from pycsamt.ff.core.edi import Edi


# Directory where edifiles are located 
edipath = r'C:\Users\Daniel\Desktop\Data\AMT\E1\edi_ss'#'data/edi' 

# save new edipath 
savepath = r'C:\Users\Daniel\Desktop\Data\AMT\E1\edi_ss2' # None 
       #) # change the current directory     
# list of Edi files 
stationlst=[os.path.join(edipath,file) for 
            file in os.listdir(edipath) if file.endswith(".edi")]

for file in stationlst: 
    Edi().write_edifile(edi_fn= file ,
                        savepath =savepath )
