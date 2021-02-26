# -*- coding: utf-8 -*-
"""
    .Script is MTpy property. a MT tool  of Alison Louise Kirkby1,Fei Zhang1, Jared Peacock2, 
    Rakib Hassan1, and Jingming Duan is already available. Get the documentation here :.
        - DOI : 10.21105/joss.01358 
        - https://mtpy2.readthedocs.io/en/develop/core.html
        
Created on Wed Apr 15 09:51:02 2015
@author: Alison Kirkby
sets up input files for running 2d occam inversions using the occam2d_rewrite module    
    
    if you are already MTpy installed on your computer with different environment variables , you may use this script to 
    build occam 2D input files. In other case , you can download MTpy through : 
        https://github.com/MTgeophysics/mtpy/wiki : see intallation guide and dependancies.
    
Edited  on Mon Feb 15 16:27:49 2021
by: @Daniel03

"""


import mtpy.modeling.occam2d_rewrite as occam2d
import os
import os.path as op
import numpy as np

# change environnment variables to MTpy 
os.environ['MTpy'] = '/mtpy_path'

# path where edi files are located
edipath = r"//data/edi"                # specify the path where edi is located 

# path to save to
savepath = r"//data/edi_output"        # specify the path to save the Occam 2D inputfiles

if not op.exists(savepath):
    os.mkdir(savepath)


# list of stations
slst=[edi[0:-4] for edi in os.listdir(edipath) if edi.find('.edi')>0]



# create an occam data object
ocd = occam2d.Data(edi_path=edipath,
                   station_list=slst,
                   interpolate_freq=True,
                   freq=np.logspace(-1,4,17)
                   )
ocd.save_path = savepath
ocd.freq_num = 17 # number of frequencies to invert for

#### make data file
# geoelectric strike for rotation
# if not specified will calculate from the data
ocd.geoelectric_strike = None

# error floors
#ocd.res_te_err = 10
ocd.res_tm_err = 10
#ocd.phase_te_err = 5
ocd.phase_tm_err = 20
ocd.model_mode= 6
ocd.write_data_file()


# make model and mesh files
ocr = occam2d.Regularization(ocd.station_locations)
# number of layers
ocr.n_layers = 31

ocr.cell_width = 5 # cell width to aim for, note 
                      # this is the mesh size (2 mesh 
                      # blocks per model block)
ocr.x_pad_multiplier = 1.7 # controls size of padding
ocr.trigger= 1.12 # controls aspect ratio of blocks
ocr.z_bottom =5000

# z1 layer and target depth in metres
ocr.z1_layer = 5
ocr.z_target_depth = 1100
ocr.save_path=ocd.save_path
ocr.build_mesh()
ocr.build_regularization()
ocr.write_mesh_file()
ocr.write_regularization_file()
ocr.plot_mesh()

# make startup file
ocs=occam2d.Startup()
ocs.iterations_to_run=100
ocs.data_fn=op.join(ocd.save_path,'OccamDataFile.dat')
ocs.resistivity_start=1.0
ocr.get_num_free_params()
ocs.param_count=ocr.num_free_param
ocs.save_path=ocd.save_path
ocs.model_fn=ocr.reg_fn
ocs.write_startup_file()

