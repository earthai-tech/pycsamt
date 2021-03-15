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

from csamtpy.ff.core.edi import Edi 

#---> set edipath 
# exmple : r'F:\__main__csamt__\paper2_data_old\data_edifiles - numStations\K1_edi\new_EDI'
path = os.path.join(os.environ['pyCSAMT'], 'data','edi') #'_outputAVG2EDI_')#

# path = r'C:\Users\Administrator\OneDrive\Python\project\pyCSAMT\csamtpy\data'

# path to save edifile # if None , save edi in your current work directory 
save_path =  None #r'C:\Users\Administrator\Desktop\test\edirewrite'   

#--> set the type of datasection 
data_section  =None  #  may be "mt" or "emap"

# get edilist 
edilist = [os.path.join(path,edifile) for edifile in os.listdir(path) if edifile.endswith('.edi')]

# print(edilist)
for edi in edilist :
    edi_obj =Edi(edi_filename=edi )
    edi_obj.write_edifile(savepath =save_path,
                          datatype=data_section )
    

