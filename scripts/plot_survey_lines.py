# -*- coding: utf-8 -*-
"""
    .Script to plot survey lines . Deal with `stn` files directly from Zonge 
    Engineering file. If such file is not at available , use 
    module  `rewrite_station_profile` from Profile class  to write *.stn file
    so to call it directly. 
    
        >>> from pycsamt.ff.site import Profile 
        >>> Profile().rewrite_station_profile(easting , northing, **kws)
        
    User can provide also provide  `X`, `Y` coordinates
    for each survey  lines  on a list of eastings and northings.
    Can customize plot using differents matplotlib properties 

Created on Mon Feb 22 22:14:59 2021
"""
 
from pycsamt.viewer.plot import Plot1d

# path to station profiles files 
path_to_profiles = 'data/stn_profiles'

# profile_lines : specify the different lines , you want to plot 
profile_lines = ['K9.stn', 'K8.stn']           # if Will plot all survey lines 
                                                # located on path_to_profiles 
                                                #directory 
# path to save figure 
savefig =None

# write multiple lines                                          
# profile_lines =[
#                 # 'K1_reaj.stn', 'K2_reaj.stn', 'K3_reaj.stn', 
#                 # 'K4_reaj.stn', 'K5_reaj.stn', 
#                 'K6_reaj.stn', 
#                 'K7_reaj.stn', 'K8_reaj.stn', 'K9_reaj.stn'
#                 ]
#scaled the lines : scale :
scale ='km'                     # can be `m` or `km` . Default is `m`

#set to False if you dont want to see stations labels 
show_station_labels = True                  # default is TRue  
 
figsize  =[5,3]                 # figure size 
    
# call object 
plot1d_obj = Plot1d( fig_size =figsize)
plot1d_obj.plot_multiStations(path = path_to_profiles, 
                               # profile_lines =profile_lines, 
                              scale =scale, 
                              savefig =savefig, 
                              show_station_labels = show_station_labels )

