# -*- coding: utf-8 -*-
"""
    .Script to write geo_outputfiles for oasis montaj of geosoft 
    To write oasis montaj file , station profile (*stn) or another similar file MUST 
    be profiles. If not file is provided , Add easting -northing arrays and if necessary an elevation array 
    Script can generate 3 additionals  outputs files if INPUT truth resistivity values  is provided. 
        of the some layers of the area.  : 
        1. One model for rho averaged (*_aver), transitory data betwen the calculated rho and true rho.
        2. second for rho value replaced . Replaced calcualted model structures resistivities 
        3. the most important files: the  cut out resistivities value with step descent (._sd) 
            show most dominant stratigraphy sequences .
    If input resistivity is not provided , will generate OCCAM2D model files with station coordinates files
 
Created on Sat Feb 20 16:07:54 2021

@author: @Daniel03

"""
import os
from geodrill.geoCore.geodrill import Geodrill 

# path to OCCAM 2D folder 
path =os.path.join(os.environ ['pyCSAMT'],  'data', 'occam2D')
path2 =os.path.join(os.environ ['pyCSAMT'],  'data', '_iter2dat_2')
profile_path = profile_fn =os.path.join(os.environ['pyCSAMT'], 'data', 'avg')
                                                                
# station coordinates profile files 
profile_fn = 'K1.stn'

# name of outputfile 
filename = 'ybkro'

#save path# path to save output files 

savepath  = None                #os.path.join(os.environ ['pyCSAMT'], 'data','_output2Oasis_2')     

#  Maximum depth investigation  for CSAMT , if not provided , will set to 1km 
DOI = '1km'                 #  can be float like 1000 = 1km 

# output resistivity to log10 
log10rho =True
# write file to negative depth : set to True if you want to keep your depth as positive value 
to_negative =True   
# ouput scalled file 
out_put_scalled_file =True             #  if set to False , scalled profile file will not generate .

#-----Read Occam 2D output files  ---------
# path to occam Data file
path_to_occam_data='OccamDataFile.dat'

# path to occam Mesh file 
path_to_occam_mesh = 'Occam2DMesh'

# path to occam Model file 
path_to_occam_model = 'Occam2DModel'
# path to Occam Iteration file 
path_to_occam_iter='ITER17.iter'

#========================================= optional params =========================================
# easting 
easting =None                           # Input easting in the case where station profile is not provided 
# northing coordinate 
northing =None                          # Input northinge in the case where station profile is not provided 
# correct coordinates : 
X = 0.                                  # correct easting value : eg : X =2800               
Y = 0.                                  # correct northing coordinates : eg Y = 3330.
# if elevation is provided set it 
elevation =None 

# Main parameters : value must be less or eagl of DOI. if step =DOi WILL return your Occam model files.
step_descent = 200.         # float value. Step to cut out data and to  force resistivites calcualted 
                            # to match the reference data as input resistivities 
                            # if not provided the step will be 20% of D0I
# Truth resistivities values otained on the sites or from other companies
#optional parameters 
input_resistivity_values = [312, 525, 1235., 2202., 4000, 7000.]   # list of resistivites values : order is insensitives 

# Truth layer names if given must match the input resistivities 
# if nname of layers not provided, program will seek in dataBase to find
# the closet resistivity to the layers.
# optional parameters 
input_layer_names = ['alluvium', 'amphibolite','altered rock','augen gneiss', 'granite'] 

# ---------------Read with Iter2DAT FILE ------------------

# see Occam_module Iter2Dat to see what file is it 

#x, y, z model file 
iter2dat_fn = 'iter17.2412.dat'
# station location file 

bln_file = 'iter17.2412.bln'

#===========================================================================================


#------Create geodrill object ----------------

geo_obj = Geodrill(   mesh_fn = os.path.join(path , path_to_occam_mesh),
                            iter_fn = os.path.join(path , path_to_occam_iter), 
                            model_fn =os.path.join(path, path_to_occam_model) , 
                            data_fn =os.path.join(path, path_to_occam_data ),

                            # iter2dat_fn = os.path.join(path2 , iter2dat_fn),
                            # bln_fn = os.path.join(path2 , bln_file),
                            input_resistivities=input_resistivity_values, 
                            input_layers =input_layer_names ,
                            step_descent = step_descent,
                            doi =DOI, 
                            )    
                               
geo_obj.to_oasis_montaj (profile_fn =os.path.join(profile_path, profile_fn) , 
                         to_negative_depth = to_negative, 
                         elevation =elevation , 
                         savepath =savepath , 
                         filename =filename, 
                         easting =easting  , 
                         northing =northing , 
                         scalled_east_north=(X,Y), 
                         output_s_XY= out_put_scalled_file, 
                         to_log10 =log10rho)