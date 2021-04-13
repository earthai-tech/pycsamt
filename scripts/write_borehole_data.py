# -*- coding: utf-8 -*-
"""
    . Script to Build Drill data . Deal only with OasisMonj Drillhole modules 
    User can build manually it own drill files by enter layer names as well has 
    the depth . program will recohnize the top and the bottom of the stratigraphy 
    sequences. User can also provided a typical file  parser : User must folow 
    how file is arranged .Once  the script is called , it will recognize the file 
    to build all the Drill data.' An examples of `parser file` are located in 
    ::r'//pyCSAMT/geodrill/data/drill_example_files'::
    Later , an other simple parser file will propose . 

Created on Sun Feb 21 19:26:21 2021

@author: @Daniel03
"""

import os 
from pycsamt.geodrill.geoCore.geodrill import Drill 


path_to_parser_files =os.path.join(os.environ['pyCSAMT'], 
                                   'geodrill', 'data', 'drill_example_files')

#name of parserfile : eg: data collected from `location `nble`. 

parser_file ='nbleDH.csv'

#savepath : path to save outfile borehole files 
savepath = None 
# data2zrite : whick kind of data do you want to output ?
# borehole geology? borehole geochemistry sample ? or borehole survey elevation ? 
kind_of_data2output = '*'       # can be '5',"all",
                                # `Collar`,`Geology`,`Sample`,`Elevation`,`Azimuth`
                                #            or   `*`.
                                #`*` is a joker , mean output all data .

# if set to True , user will add step by step all data with the layer thicknesses 
buid_borehole_manually =False 
# if elevation of all borehole is there , upload , will take account 
add_elevation =None 
# if azimuth of all borehole is available , call it 
add_azimuth = None  

borehole_obj = Drill (well_filename= os.path.join(path_to_parser_files, parser_file)  , 
                   build_manually_welldata= buid_borehole_manually)
borehole_obj._collar()
borehole_obj.dhGeology()
borehole_obj.dhSample()
borehole_obj.dhSurveyElevAz(add_azimuth=add_azimuth, 
                            add_elevation=add_elevation)
borehole_obj.writeDHData(data2write=kind_of_data2output, 
                         savepath  = savepath )

