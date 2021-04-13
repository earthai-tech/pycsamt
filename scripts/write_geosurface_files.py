# -*- coding: utf-8 -*-
"""
    .Script to generate the geosurface data for geosurface map. 
    Deal with Oasis montaj software and Goldenn softwares like `surfer`. 
    From occam2D model files to oasis outfiles model files . please see the module 
    `geodrill` to build output files before using this script. Once files are built
    can visualize the resistivity distribution of different lines at that depth.
    
Created on Wed Nov 25 11:39:11 2020


@author: @Daniel03

"""
import os 
from pycsamt.geodrill.geoCore.geodrill import Geosurface 

#Oasis output files are input files  to use geosurface module 
path_to_oasisfiles = os.path.join(os.environ['pyCSAMT'], 'data',
                                    'InputOas')

# path to save geosurfaces outputfiles 
savepath =None                      # if None , will create a directory 
# depth values in meters 
values_for_imaging = [38, 100]      # mean surface at 38 m deep and 100 m deep 

# file format is default output format : 
output_format = 'csv'               # could be '`xlsx` or `csv
                                    # default is `csv`
                                    

# call geosurface object 
geo_surface_obj = Geosurface( path =path_to_oasisfiles, 
                             depth_values = values_for_imaging, 
                             )

geo_surface_obj.write_file(fileformat = output_format, 
                            savepath =savepath )

