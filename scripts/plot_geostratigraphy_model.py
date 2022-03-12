# -*- coding: utf-8 -*-
"""
    Plot different models
    - Occam Model by setting `kind`=``crm`` 
    - strata model by setting `kind`=``nm`` 
    - Misfit by setting kind ='strata' and `plot_misfit` to ``True``
Created on Wed Sep 29 09:50:19 2021

"""

import os
from pycsamt.geodrill.geocore import GeoStratigraphy 

# path to OCCAM2D model files
occamPath = 'data/Occam2d'
# type of model plots
kindOfPlot ='nm'      # for New strata model, `crm` forOccamResistivity model` plot
                # `strata` and `plot_misfit` to True plot misfitG 
# set to True to plot error between CRM and NM
plotMisfitG=False  
inversion_files = {'model_fn':'Occam2DModel', 
                    'mesh_fn': 'Occam2DMesh',
                    "iter_fn":'ITER17.iter',
                    'data_fn':'OccamDataFile.dat'
                    }
inversion_files = {key:os.path.join(occamPath , vv) for key,
                    vv in inversion_files.items()}
# input_True_resistivities (TRES)
TRES= [10,70, 180, 1000,   3000, 7000]     
#[10, 60, 70, 180, 1000,  3000, 7000]                                 
# Input layers names (LN) 
LN =['MWG2', 'MWG1', 
     'FG', 'LWG', 'Igneous rocks','Basement rocks' ]
# ['river water','sedimentary rocks', 'fracture zone', 
#     'gravel', 'granite','igneous rocks','basement rocks' ]
geosObj = GeoStratigraphy(**inversion_files,
                      input_resistivities=TRES, 
                      input_layers=LN)

_= geosObj.strataModel(kind=kindOfPlot ,
                       misfit_G =plotMisfitG)
