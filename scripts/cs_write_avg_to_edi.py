# -*- coding: utf-8 -*-
"""
    .Script to write edifiles from Zonge Engineering AVG format . Script accept 3 filters "
    "TMA" : Triming Moving Average , "AMA": Adaptative Moving Average of Torres abd verding 
    "FLMA": Fixed Length dipole moving average . Actuallyfilters "TMA & FLMA" work, 
    AMA filter will be add soon. Once applied it , it will generate Rho corrected into 
    your output edifiles. If filter is None , just write edifile , no filter will be applied. 
    If  zonge astatic file is given , program will take account to write both resistivites . 
    the non corrected apparent resistivities RHO and the static shift corrected FRHO. 
   
Created on Sat Jan 16 19:55:58 2021

@author: @Daniel03
"""

import os 
from csamtpy.ff.core import avg 


#--> path to avg file : example where avg file  is located
# avg file and station profile file must be located in the same path  
path_to_avgfile =os.path.join(os.environ ['pyCSAMT'], 'data', 'avg') 

#--savepath :
save_edipath = None         #r'C:\Users\Administrator\Desktop\test\edi_from_avg'

#path to Zonge AVG file
avgfile = 'K1.AVG'
#--- > add profile file : 
station_profile_file = 'K1.stn'

#---Appply filter . 
add_filter = None              # can be "tma" or "lma" filters , if filter is not None
                                # will compute FRHO
                               
# add reference frequency 
reference_frequency = None      #8192  # if reference frequency is None AND add filter is not None , 
                                # reference frequency  will be compute automatically  
#-- Call avgobject 

avg_obj= avg.Avg()

avg_obj.avg_to_edifile(data_fn= os.path.join(path_to_avgfile, avgfile) , 
                       profile_fn = os.path.join(path_to_avgfile, station_profile_file), 
                       savepath =save_edipath, 
                       reference_frequency= reference_frequency, 
                       apply_filter=add_filter ) 

