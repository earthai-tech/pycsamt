# -*- coding: utf-8 -*-
"""
Created on Fri Apr 16 15:11:29 2021

@author: @Daniel03
"""

# import required modules
import os 
from pycsamt.ff.core.avg import Avg 
from pycsamt.viewer.plot import Plot1d
from pycsamt.ff.core.cs import Profile 
# create profile_obj 
profile_obj = Profile(station_profile)
# to get 
profile_obj.east    # get easting coordinate  
profile_obj.north   # get northing coordinates 
profile_obj.lon 
profile_obj.lat 

profile_obj.dipole_length 
profile_obj.stn_interval 
profile_obj.stn_interval 
profile_obj.azim 
profile_obj.Site.stn_name

profile_obj.Site.east['S00']
profile_obj.Site.north['S00']


# Plot1d().plot_curves(fn = avg_file, selected_stations=['S00', 'S22', 'S46'])
# convert avg file to SEG EDI file .
# edipath ='data/edi' # path to edi-files 
# selected_frequencies =[1024, 4096., 8000.] # selected frequencies to visualize
# from pycsamt.viewer.plot import Plot2d
# Plot1d().plot2d_obj.penetration2D(fn = edipath, doi ='2km') # can be doi=2000. 
