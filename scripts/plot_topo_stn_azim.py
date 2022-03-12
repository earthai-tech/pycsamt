# -*- coding: utf-8 -*-
"""
    . Script to plot Topography _station_separation and azimuth . 
    Note : Most of matplotlib property can be handled on plot function to
    get the the kind of profile you want. 
    
Created on Wed Dec 30 19:27:04 2020

@author:K.L ~ @Daniel03

"""

from pycsamt.viewer.plot import Plot1d

# set the stn profile file or EDI-file
# set path to your profile file 

file_stn='K6.stn'           # name of zonge station profile file 
# uncomment `path_to_stn_profile_file is zonge station file is used 
path_to_stn_profile_file = None #  os.path.join('data', 'stn_profiles',file_stn)
            # OR 
            
# provided edipath or jpath when used edifiles or jfiles and 
edipath_or_jpath = 'data/edi' # None 
#edipath_or_jpath =None                     # uncomment section if station stn file is provided 

plot_type ='*'          # could be |"topography or "topo" |"speration" or "stn"| "azimuth or "az|
                        # could use only the two first letter or more for ploting or 1,2,3 or 123|*
set_stnNames =True      # set it to True if you want to see station names appear on xaxis 


# create profile_obj 
plot_1d_obj= Plot1d()
plot_1d_obj.plot_topo_sep_azim(fn = edipath_or_jpath,
                               profile_fn= path_to_stn_profile_file ,
                               plot=plot_type,
                               set_station_names=set_stnNames,
                               )

