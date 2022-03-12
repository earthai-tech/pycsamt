# -*- coding: utf-8 -*-
"""
    . script to rewrite station coordinates, 
    correct coordinates values for UTM projection 

Created on Tue Jan  5 17:38:54 2021

@author:K.L ~ @Daniel03

"""
import os 
from pycsamt.ff.site import Profile

# --- > path to  your stn file
path ='data/stn_profiles' 

#--> original stn file . It could be the Station file. 
file_stn='K6.stn'

# For UTM projection .Sometimes need to adjust 
#X -easting and Y-northing , If None : value will be 0. 

utm_X_cor=0.                                #-300238.702 
utm_Y_cor=0.                                #-2369.252

# output new station location file. 
REWRITE=True 

# savepath 
savepath = None                       
# link obj with path 
path_fn = os.path.join(path, file_stn)

# read multifile stations files 

lst=[os.path.join(path, file) 
     for file in os.listdir(path) 
     if file.endswith('.stn')]

# call  profile_obj
for stn_file in lst:
    profile =Profile.reajust_coordinates_values(x=utm_X_cor ,
                                                y=utm_Y_cor , 
                                                stn_fn=stn_file ,
                                                rewrite = REWRITE, 
                                                savepath=savepath )
