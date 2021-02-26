# -*- coding: utf-8 -*-
"""
    script to write Avg file ASTATIC to plain File 
    
Created on Fri Dec 11 20:26:52 2020

@author: @Daniel03

"""

import os 
from csamtpy.ff.core.avg import Avg 

#set your directory to your astatic file and your savepath  

path = r'F:\__main__csamt__\可控源资料\汝城资料_Infos sur la ville\CSAMT原始数据_Donnees brutes de CSAMT\avg_F2'
savepath =r'C:\Users\Administrator\Desktop'

avg_files =[file for file in  os.listdir(path) if file.endswith('.AVG')]
# print(avg_files)

for file in avg_files :
    Avg().avg_write_2_to_1(data_fn=os.path.join(path, file),
                           savepath =savepath)
