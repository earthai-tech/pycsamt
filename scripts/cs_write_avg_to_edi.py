# -*- coding: utf-8 -*-
"""
    .Script to write edifiles from Zonge Engineering AVG format . Script accept 3 filters "
    "TMA" : Triminf Moving Average , "AMA": Adaptative Moving Average of Torres abd verding 
    "FLMA": Fixed Lend dipole Moving average . Actually only filter "TMA" works , for the next upgrade , 
    will input the AMA and FLMA. Once applied it , it will generate Rho corrected into 
    your output edifiles. If filter is None , just write edifile , no filter will be applied. 
    iF astatic filed is given , programm will take account to write both resistivites . 
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
add_filter =None            # can be "TMA", "AMA" or "FLMA" filters :
                            #actually only TMA is  available, work for AMA and FLMA are 
                            # on progress


#-- Call avgobject 

avg_obj= avg.Avg()

avg_obj.avg_to_edifile(data_fn= os.path.join(path_to_avgfile, avgfile) , 
                       profile_fn = os.path.join(path_to_avgfile, station_profile_file), 
                       savepath =save_edipath, 
                       apply_filter=add_filter ) 

