# -*- coding: utf-8 -*-
"""
    .Script to write geo_outputfiles to Golden software plots 
    will generate 3 outputs : 
        1. One model for rho averaged (*_aver), transitory data betwen the 
            calculated rho and true rho.
        2. second for rho value replaced . Replaced calcualted model 
            structures resistivities 
        to their closest resistivities as resistivities reference from input 
            resistivities.
        3. the most important files: the  cut out resistivities value with
            step descent (._sd) 
            show most dominant stratigraphy sequences .
        4. the station location file (.bln)
        
        Plot the 3 files in Golden software to see transition from model
        calculation to truth sequence detail models
        which is most closest to reality.

Created on Mon Feb 15 14:56:18 2021
"""
import os
from pycsamt.geodrill.geocore import Geodrill 

# path to OCCAM 2D folder 
path = 'data/occam2d'
# prefix outputfile
filename = 'k8' #'ybkro'
savepath ='data/saveGS'
#  Maximum depth investigation  for CSAMT,
# if not provided , will set to 1km 
 #  can be float like 1000 = 1km 
DOI = '1km'                

# out put file either "m" or "km". Default is "meter"
scale =None 

# write file to negative depth: 
#set to True if you want to keep 
#your depth as positive value 
to_negative =True                       

# if elevation is provided set it 
elevation =None 

# Main parameters : value must be less or
STEP_DESCENT = 200.         
# list of resistivites values : order is insensitives  
INPUT_RESISTIVITIES = [66,70, 180, 1000, 3000, 10000, 20000] 
# Truth layer names if given must match the input resistivities 
# if nname of layers not provided, program will seek in dataBase to find
# the name of layer whom its resistivities is much closer to abose 
#INPUT_LAYERS = ['alluvium', 'amphibolite','altered rock',
#'augen gneiss', 'granite'] #
INPUT_LAYERS = ['river water', 'fracture zone' , 'granite ',
                'Most Weathered', 'Less Weathered'] # 

#-----Read Occam 2D output files  ---------
# path to occam Data file
path_to_occam_data='OccamDataFile.dat'

# path to occam Mesh file 
path_to_occam_mesh = 'Occam2DMesh'

# path to occam Model file 
path_to_occam_model = 'Occam2DModel'
# path to Occam Iteration file 
path_to_occam_iter='ITER17.iter'

# ---------------Read with Iter2DAT FILE ------------------

# see Occam_module Iter2Dat to see what file is it 

#x, y, z model file 
iter2dat_fn = 'iter17.2412.dat'
# station location file 

bln_file = 'iter17.2412.bln'
#------Create geodrill object ----------------

geo_obj = Geodrill(  input_resistivities=INPUT_RESISTIVITIES, 
                            input_layers =INPUT_LAYERS ,
                            step_descent =STEP_DESCENT,
                            doi =DOI, 
                            mesh_fn = os.path.join(path , path_to_occam_mesh),
                            iter_fn = os.path.join(path , path_to_occam_iter), 
                            model_fn =os.path.join(path, path_to_occam_model) , 
                            data_fn =os.path.join(path, path_to_occam_data ),
                            # iter2dat_fn = os.path.join(path2 , iter2dat_fn),
                            # bln_fn = os.path.join(path2 , bln_file),
                            )    
geo_obj.to_golden_software(filename =filename , 
                               elevation =elevation, 
                               savepath = savepath, 
                               scale =scale, 
                               to_negative_depth=to_negative)     
                                  