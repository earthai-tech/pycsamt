# -*- coding: utf-8 -*-
"""
    . script to rewrite station coordinates , correct coordinates values for UTM projection 

Created on Tue Jan  5 17:38:54 2021

@author: @Daniel03
"""
import os 
from csamtpy.ff.core.cs import Profile as prof_obj

# --- > path to  your stn file
path =os.path.join(os.environ['pyCSAMT'], 'data', 'stn_profiles') 
#path =r'C:\Users\Administrator\Desktop\ThesisImp\raw_stn_file'

#--> original stn file . It could be the Station file. 
file_stn='K6.stn'

# For UTM projection .Sometimes need to adjust X -easting and Y-northing , If None : value will be 0. 

utm_X_cor=0.                                #-300238.702 
utm_Y_cor=0.                                #-2369.252

# to see documentation , set to True
see_doc =False
# set to False if you just want to see documentation . 
REWRITE=True 

# savepath 
savepath = None                         #r'C:\Users\Administrator\Desktop\ThesisImp\scaled_stn'
# link obj with path 
path_fn = os.path.join(path, file_stn)

# read multifile stations files 

lst=[os.path.join(path, file) for file in os.listdir(path) if file.endswith('.stn')]

# call  profile_obj
if REWRITE :
    for stn_file in lst:
        profile =prof_obj.reajust_coordinates_values(x=utm_X_cor ,
                                                    y=utm_Y_cor , 
                                                    stn_fn=stn_file ,
                                                    rewrite = REWRITE, 
                                                    savepath=savepath )
        
if see_doc: help(prof_obj.reajust_coordinates_values)