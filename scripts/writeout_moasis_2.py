# -*- coding: utf-8 -*-
"""
    .Script to write geo_outputfiles for oasis montaj of Geosoft corporation
    To write oasis montaj file, station profile (*stn) or another similar file 
    MUST be supplied otherwise add easting and northing arrays and/or 
    elevation array(Optional). The script can generate 3 additionals  outputs 
    files if INPUT true resistivity values  are supplied.
        1. model for rho averaged (*_aver), transitory data betwern the
            calculated rho and true rho.
        2. second for the rho valuereplaced. Replaced calcualted resistivities
            model structures
        3. the most important files: the  cut out resistivities value 
            with step descent (._sd) 
            show most dominant stratigraphy sequences .
    If input resistivities are not provided, will generate OCCAM2D model files
        with station coordinates files
 
Created on Sat Feb 20 16:07:54 2021
"""
import os
from pycsamt.geodrill.geocore import Geodrill 

# path to OCCAM 2D folder 
path ='data/occam2D'
# if you want to write 
path2 = 'data/iter2dat'
# path to occam 2d model files 
oc2dkws= {k: os.path.join(path, v) for k, v in 
                dict(
                data_fn='OccamDataFile.dat',
                mesh_fn = 'Occam2DMesh',
                model_fn = 'Occam2DModel',
                iter_fn='ITER17.iter',
                ).items(
                    )
}
#-------------
#---
#if you want to write model fom 
#x, y, z model file and set the oc2dkws as empt dict
i2dkws= dict (
    iter2dat_fn = 'data/iter2dat/iter17.2412.dat',
    bln_file = 'data/iter2dat/iter17.2412.bln',
)
                                                          
# station coordinates profile files 
profile_fn = 'data/avg/K1.stn'
# name of outputfile 
filename = 'ybkro'
   
savepath =None
#  Maximum depth investigation  for CSAMT,
# if not provided , will set to 1km 
#  can be float like 1000 = 1km 
DOI = '1km'                 

# output resistivity to log10 
log10rho =True
# write file to negative depth:
# set to True if you want to keep your depth as positive value 
to_negative =False  
# ouput scalled file 
#  if set to False , scalled profile file will not generate .
out_put_scalled_file =True             

#================== optional params ============
# easting 
easting =None                           
# Input easting in the case where station profile is not provided 
# northing coordinate 
northing =None                          
# Input northinge in the case where station profile is not provided 
# correct coordinates : 
X = 0.                   # correct easting value : eg : X =2800               
Y = 0.                   # correct northing coordinates : eg Y = 3330.
# if elevation is provided set it 
elevation =None 

# Main parameters : value must be less or 
step_descent = 200.         
input_resistivity_values =[66, 70, 180, 1000, 3000 ]
# optional parameters 
input_layer_names =['river water', 'fracture zone', 'granite']
#=============================================================

#------Create geodrill object ----------------
geo_obj = Geodrill( 
                input_resistivities=input_resistivity_values, 
                input_layers =input_layer_names ,
                step_descent = step_descent,
                doi =DOI, 
                **oc2dkws, 
                # **i2dkws
                )    
                               
geo_obj.to_oasis_montaj (profile_fn = profile_fn , 
                         to_negative_depth = to_negative, 
                         elevation =elevation , 
                         savepath =savepath , 
                         filename =filename, 
                         easting =easting  , 
                         northing =northing , 
                         scalled_east_north=(X,Y), 
                         output_s_XY= out_put_scalled_file, 
                         to_log10 =log10rho)


