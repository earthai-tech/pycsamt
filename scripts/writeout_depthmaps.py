# -*- coding: utf-8 -*-
"""
    .Script to generate the depth maps. 
    Deal with Oasis montaj software and Goldenn softwares like `surfer`. 
    Please refer to the module `pycsamt.geodrill.geocore.geodrill` 
    to build  the output files before using this script.
    Once the files are built, it is possible to visualize the resistivity 
    distribution of different lines at each gievn depth.
    
Created on Wed Nov 25 11:39:11 2020
@author:K.L ~ @Daniel03
"""

from pycsamt.geodrill.geocore import Geosurface 

#Oasis output files are input files  to use geosurface module 
path_to_oasisfiles = 'data/InputOas'
# path to save geosurfaces outputfiles 
savepath = None                   
# depth values in meters 
values_for_imaging = [38, 100]      # mean surface at 38 m deep and 100 m deep 
values_for_imaging =494.

# file format is default output format : 
output_format = 'xlsx'              # could be '`xlsx` or `csv
                                    # default is `csv`
                                
# call geosurface object 
geo_surface_obj = Geosurface( path =path_to_oasisfiles, 
                             depth_values = values_for_imaging, 
                             )
geo_surface_obj.write_file(fileformat = output_format, 
                            savepath =savepath )

