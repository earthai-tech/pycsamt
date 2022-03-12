# -*- coding: utf-8 -*-
"""
    .Script to write edifiles from Zonge Engineering AVG format .
    Script accept 3 filters , TMA" : Triming Moving Average ,
     "AMA": Adaptative Moving Average of Torres-Verdin
    "FLMA": Fixed Length dipole moving average . If filter is applied 
    Rho corrected shoud be in output edifiles as FRHO. If filter
    is None , just write edifile. If  zonge astatic file is given , it assumes 
    a filter is already applied, then will read `avg` file and set filter 
    resistivities values as FRHO, means , `edi`outputfile 
    will contain static shift corrected FRHO and Uncorrected resistivities. 

Created on Sat Jan 16 19:55:58 2021
@author: @Daniel03
"""

import os 
from pycsamt.ff.core import CSAMT

#--> path to avg file : example where avg file  is located
# avg file and station profile file must be located in the same path  

path_to_avgfile= 'data/avg'
#--savepath :
#save_edipath =  None         
save_edipath ='data/K1_edi'
#path to Zonge AVG file
avgfile = 'K1.AVG'
#--- > add profile file : 
station_profile_file = 'K1.stn'

#---Appply filter . 
add_filter = None              # can be "tma", 'flma' or "ama" filters , if filter is not None
                                # will compute FRHO
#------------------------------------------------------------------------------
# when applied filter , can set the number of points for FLMA or TMA filter and number of skin depth 
# for AMA filters. if add filter is None , numberof points and skin depth will turn off. 
number_of_points =7.            # default is 7. can be change to any value you want 
number_of_skin_depth =7.        # default is 3. 
#------------------------------------------------------------------------------                       
# add reference frequency 
reference_frequency = None      #8192  # if reference frequency is None AND add filter is not None , 
                                # reference frequency  will be compute automatically  
#-- Call avgobject 

avg_obj= CSAMT().avg2edi(data_fn= os.path.join(path_to_avgfile, avgfile) , 
                       profile_fn = os.path.join(path_to_avgfile, 
                                                 station_profile_file), 
                       savepath =save_edipath, 
                       reference_frequency= reference_frequency, 
                       apply_filter=add_filter, 
                       number_of_points=number_of_points, 
                       number_of_skin_depth=number_of_skin_depth) 

