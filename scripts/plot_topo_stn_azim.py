# -*- coding: utf-8 -*-
"""
    . Script to plot Topography _station_separation and azimuth . 
    Note : Most of matplotlib property can be handled on plot function to
    get the the kind of profile you want. To see documentation of Plot topography-
    station_separation and Azimuth . set to True "see_documentation"  
    
Created on Wed Dec 30 19:27:04 2020

@author: @Daniel03

"""
import os 
from viewer.plot import Plot1d

# set the pyCSAMT environment : 
os.environ['pyCSAMT']=r'C:/Users\Administrator\OneDrive\Python\project\pyCSAMT'

# set the stn profile file or EDI-file  
file_stn='K6.stn'
edipath_or_jpath = os.path.join(os.environ["pyCSAMT"],'csamtpy','data')

# set path to your profile file 
path =  os.path.join(os.environ["pyCSAMT"],'csamtpy','data', file_stn)

# to see documentation of that function , set "see_documentation " to 'true'. 
see_documentation=False
# to plot profile : set PLOT to "True"
PLOT=True 
# set "*" for three profile. 

plot_type ='3'          # could be |"topography or "topo" |"speration" or "stn"| "azimuth or "az|
                        # could use only the two first letter or more for ploting or 1,2,3 or 123|*
set_stnNames =True      # set it to True if you want to see station names appear on xaxis 


# create profile_obj 
plot_1d_obj= Plot1d()

if PLOT: 
    plot_1d_obj.plot_topo_sep_azim(fn = path , profile_fn= path , 
                                   plot=plot_type,
                                   set_station_names=set_stnNames,
                                   )
    
if see_documentation : help(plot_1d_obj.plot_topo_sep_azim)
    #print (plot_1d_obj.plot_topo_sep_azim.__doc__)