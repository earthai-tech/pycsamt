# -*- coding: utf-8 -*-
"""
    script to write Avg file ASTATIC to plain File 
    
Created on Fri Dec 11 20:26:52 2020

@author: @Daniel03

"""

import os 
from pycsamt.ff.core.avg import Avg 

#set your directory to your astatic file and your savepath  
path = os.path.join(os.path.abspath('.'), 'data', 'avg', 'K2.AVG')# to file

savepath =None

# for mutltiread 
    #path = os.path.join(os.path.abspath('.'), 'data', 'avg')
    #avg_files =[file for file in  os.listdir(path) if file.endswith('.AVG')]
    
avg_files = [path]

for file in avg_files :
    Avg().avg_write_2_to_1(data_fn= file , #os.path.join(path, file),
                           savepath =savepath)
