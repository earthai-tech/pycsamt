# -*- coding: utf-8 -*-
"""
    . Script  to straighten profile .
    The purpose of this script is to straighten out the coordinates taken during the data collection.
    Due to the aleas of the terrain being in mountainous or in forest regions, the land harvesting
    is not an exact distance. This script offers several possibilities for recovery: 
        
    1- Classical reajustment {"classic"} in linear regime. Here all the coordinates are straightened in the form
        of a function y = f (dl) : dl = dipole length . 
    2. Natural or distorted recovery {"natural|distored"}. In this type of reajustment the staion, seperation respects 
        the distance  of the dipole, but the poinst do not form a line regime.
        f (dl) = a dl^m+ b dl^s+...+z^0 
       
    3. Equidistant recovery {"equidistant"}. It doesnt based on the computation of dipole length, just use derived 
        from distanciation of coordinates, compute the mean between each station and and reajusted cordinates. 
        **Default** is 'classic'
        User can also take advantage at the same time to correct easting and northing coordinates using reajust
        arguments (x, y) ; index  
        
Created on Fri Jan 29 18:57:32 2021

@author: @Daniel03

"""
import os 
from viewer.plot import Plot1d
from csamtpy.ff.core.cs import Profile

# full path to your profile file : stn 
path =  os.path.join(os.environ["pyCSAMT"],# 'STN-28800_reaj.stn')
                      'data', 'avg', 'K1.stn')            # change your path 

#saveyour figure path 
savepath =None 

# if stn file is not available , set CREATE to True and REAJUST to False
# Then you will generate your own stn file r=for reajusting coordinates
CREATE =False
REAJUST=True 
# choose the type to straithen your profile 
straigthen_out_mode ='classic'          # can be 'distored|natural, or "equisistant".Default is "classic".
# correct coordinates : 
X = 0.                                  # correct easting value : eg : X =2800               
Y = 0.                                  # correct northing coordinates : eg Y = 3330.

# set to True if ou want to get new coordinates values 
ouputnew_stnfile = True 
#===============================================================================================
# if Zonge *.stn file is available , dont need to fill this part , set only REAJUST to True.
easting = None                         # (ndarray, 1) or list # station easting coordinates(required) 
northing =None                         # (ndarray,1) or list # station northing coordinates (required)
elevation = None                       # (ndarray,1) or list  # station elevation value (optional)
compute_azimuth = False                # if set to True , will compute azimuth (optional)
username =None                         # user name will use for customizing output (optional)
location_name =None                    # name of location of area .(optional )
output_profile_name =None              # name of output file  (optional)

#NB : easting, northing and elevation if given MUST have the same size|length. 
#=====================================================================================================


if CREATE :
    #create profile object
    profile_obj =Profile()
    profile_obj.rewrite_station_profile (easting =easting ,
                                         northing=northing ,
                                         elevation =elevation , 
                                         area_name =location_name, 
                                         username =username, 
                                         add_azimuth =compute_azimuth, 
                                         savepath =savepath , 
                                         output_name =output_profile_name )
if REAJUST : 
    #create plot object 
    plot_1d_obj= Plot1d()    
    plot_1d_obj.plot_station_profile(fn = path, 
                                     reajust_coordinates=(X,Y),
                                     straighten_type =straigthen_out_mode, 
                                     savefig=savepath  , 
                                     outputfile =ouputnew_stnfile  
                                     )
