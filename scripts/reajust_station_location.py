# -*- coding: utf-8 -*-
"""
    . script to rewrite station coordinates , correct coordinates values for UTM projection 

Created on Tue Jan  5 17:38:54 2021

@author: @Daniel03
"""
import os 
from csamtpy.ff.core.cs import Profile as prof_obj

# --- > path to  your stn file 
path =r'C:\Users\Administrator\OneDrive\Python\project\pyCSAMT\csamtpy\data'
#--> original stn file . It could be the Station file. 
file_stn='K1.stn'

# For UTM projection .Sometimes need to adjust X -easting and Y-northing , If None : value will be 0. 

utm_X_cor=None
utm_Y_cor=None

# to see documentation , set to True
see_doc =False
# set to False if you just want to see documentation . 
REWRITE=True 

# link obj with path 
path_fn = os.path.join(path, file_stn)



# call  profile_obj
if REWRITE :
    profile =prof_obj.reajust_coordinates_values(x=utm_X_cor ,
                                                y=utm_Y_cor , 
                                                stn_fn=path_fn ,
                                                rewrite = REWRITE, 
                                                savepath=None )
if see_doc: help(prof_obj.reajust_coordinates_values)