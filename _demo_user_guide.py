# -*- coding: utf-8 -*-
"""
Created on Fri Apr 16 15:11:29 2021

@author: @Daniel03
"""

# import required modules
import os 
from pycsamt.ff.core.avg import Avg 
from pycsamt.viewer.plot import Plot1d
from pycsamt.viewer.plot import Plot2d
from pycsamt.ff.core.cs import Profile
from pycsamt.ff.core.cs import Site
# create profile_obj 
profile_obj = Profile(station_profile)
# to get 
profile_obj.east    # get easting coordinate  
profile_obj.north   # get northing coordinates 
profile_obj.lon     # longitude value 
profile_obj.lat     # latitue value 
# to get dipole length in meters 
profile_obj.dipole_length 

profile_obj.stn_interval  # interval between stations 
profile_obj.stn_position  # scaled position of each stations 
profile_obj.azimuth          #  azimuth of profile_line 


profile_obj.Site.stn_name # station id 

# to get single value of latitude and longitude or easing northing at each station 
# call Site obj 
# 
straighten_out_mode ='classic' # can be 'natural/distord' or equidistant
# contrinute value for x,y coordinates hidden. 
adjust__x_utm_coordinates = -300238.702 
adjust__y_utm_coordinates = -2369.252
get_new_station_profile =True 
Plot1d().plot_station_profile(fn = station_profile_file, 
                                 reajust_coordinates=(adjust__x_utm_coordinates,
                                                      adjust__y_utm_coordinates),
                                 straighten_type =straighten_out_mode , 
                                 outputfile =get_new_station_profile, 
                                 savefig=savepath)

