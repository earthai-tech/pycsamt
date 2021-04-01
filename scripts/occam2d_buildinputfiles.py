# -*- coding: utf-8 -*-
"""
    .Script to generate occam2d building inputfiles from 'MTpy' module. 
    MTpy is a magnetotelluric toolbox from authors Alison Louise Kirkby1,Fei Zhang1, Jared Peacock2, 
    Rakib Hassan1, and Jingming Duan is already available. Get the documentation here :.
        - DOI : 10.21105/joss.01358 
        - https://mtpy2.readthedocs.io/en/develop/core.html
        
Created on Wed Apr 15 09:51:02 2015
@author: Alison Kirkby
sets up input files for running 2d occam inversions using the occam2d_rewrite module    
    
    if you are already MTpy installed on your computer you may use this script to 
    build occam 2D input files.  if 'MTpy' is not intalled 'pyCSAMT' will try tp 
    install mtpy with its dependancies. To avoid some packages conflits , better 
    approch is to create a virtual environnement to use 'pyCSAMT' . If automatic
    downloading failed ,  you can download MTpy  floowing steps in MTpy wiki page : 
        https://github.com/MTgeophysics/mtpy/wiki : see intallation guide and dependancies.
    
Edited  on Mon Feb 15 16:27:49 2021
by: @Daniel03

"""

from  csamtpy.modeling.occam2d import occam2d_write

import os

#path where edi files are located
edipath = os.path.join(os.environ['pyCSAMT'], 'data', 'edi') # specify the path where edi is located 

 # specify the path to save the Occam 2D inputfiles
savepath = os.path.join(os.path.abspath('.'),
                        'data', 'occam2dBuildInputfiles')

# occam_output_dataname 
OccamDataFile= 'OccamDataFile.dat'
StartupFile_name= 'Startup'

# if interpolate frequency is set to true , bring limit of interpolatation il logspace frequency 

interpolate_frequency =True 
number_of_frequency = 17            # number of frequency for interpolated 

frequency_interpolate_range = (-1,4,number_of_frequency ) # in logspace log10. last item is number of frequency 

# ---> create model 
# number of layers
number_of_model_layers = 31.
# investigation depth 

expected_investigation_depth_or_z_target = 1100. 
# exaggeration depth , must be enough as possible. Around 5* time  the investigation depth 

exaggerate_vertical_z_bottom =5000.      # exagerate bottom 

model_first_layer_thickness_z1 = 5.     # first layer value 

# starting resistivity model value in log 10 res 

starting_res_model_value = 1. 

#number_of_iteration to run 
number_of_iteration_to_run = 100.

# bring geoelectrik strike if provided 
geoelectric_strike =34.                     # if not given set to 0.
# error floors 
    #----TM mode 
occam_mode ='6'         # must be string 
TM_phase_error= 20.      # in percentage 
TM_res_error =10.        # in percentage 

    #TE mode 
#occam_mode = '5'
TE_phase_error= 20.      # in percentage 
TE_res_error =10.        # in percentage 
 
# Configuration mesh features 

model_cell_width = 5 # cell width to aim for, note 
                      # this is the mesh size (2 mesh # blocks per model block)
horizontal_node_x_pad_multiplier = 1.7 # controls size of padding
brick_trigger = 1.12            # controls aspect ratio of blocks


occam2d_write.buildingInputfiles(edi_fn =edipath,
                           geoelectric_strike= geoelectric_strike,
                           interpolate_freq= interpolate_frequency, 
                           intp_freq_logspace =frequency_interpolate_range, 
                            iteration_to_run= number_of_iteration_to_run, 
                            resistivity_start = starting_res_model_value, 
                            startup_basename=StartupFile_name,
                            res_tm_err= TM_res_error, 
                            phase_tm_err= TM_phase_error , 
                            occam_mode= occam_mode, 
                            n_layers = number_of_model_layers , 
                            cell_width = model_cell_width , 
                            x_pad_multiplier = horizontal_node_x_pad_multiplier, 
                            trigger =brick_trigger, 
                            z_bottom =exaggerate_vertical_z_bottom , 
                            z1_layer =model_first_layer_thickness_z1, 
                            z_target = expected_investigation_depth_or_z_target, 
                            occamDataFile = OccamDataFile, 
                        )
    

