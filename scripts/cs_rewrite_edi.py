# -*- coding: utf-8 -*-
"""
    .Script to rewrite EDI file. Can Rewrite MT edifile or EMAP edifile . 
    User can force  the programm to write file into different section
     by setting argument"datatype" either to 'mtsect' or 'emapsect'.
    If set to 'emap' programm automatically will compute according SEG 
    "emap section" resistivity and phase by considering that 
    it's electromagnetic array profile. If CSAMT edifile will provide , 
    it will compute different resistivities  according to different components . 
   
    
Created on Thu Jan 14 22:00:21 2021
@author: @Daniel03

"""
import os 

from pycsamt.ff.core.edi import Edi 

#---> set edipath 
path = 'data/edi' 
# path to save edifile # if None , save edi in your current work directory 
save_path =  None 

#--> set the type of datasection 
data_section  ='emap'  #  may be "mt" or "emap"

# add new output edi name
new_edi_name = None # 'k6'
# get edilist 
edilist = [os.path.join(path,edifile) 
           for edifile in os.listdir(path) 
           if edifile.endswith('.edi')]

# print(edilist)
for edi in edilist :
    edi_obj =Edi(edi_filename=edi )
    edi_obj.write_edifile(savepath =save_path,
                          datatype=data_section,
                          new_edifilename =new_edi_name )
    

