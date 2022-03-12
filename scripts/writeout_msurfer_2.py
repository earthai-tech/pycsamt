# -*- coding: utf-8 -*-
"""
    .Script to write models outputfiles to Golden software plots 
        1. One model for rho averaged (*_aver), transitory data betwen
            the calculated rho and true rho.
        2. second for rho value replaced . Replaced calcualted model 
            structures resistivities to their closest resistivities as
                resistivities reference from input resistivities.
        3. the most important files: the  cut out resistivities 
            value with step descent (._sd) 
            show most dominant stratigraphy sequences .
        4. the station location file (.bln)
        
        Plot the 3 files in Golden software to see the difference 
        the between the calculated model with the sequence detail models
        alow to conclude the efficienty of the proposed approach.
        
Created on Mon Feb 15 14:56:18 2021
@author:K.L ~ @Daniel03

"""
import os
from pycsamt.geodrill.geocore import Geodrill 

# path to OCCAM 2D folder 
path = 'data/occam2D'
#-----Read Occam 2D output files  ---------
# path to occam model files
oc2d_data= {k: os.path.join(path, v) for k, v in 
                dict(
                data_fn='OccamDataFile.dat',
                mesh_fn = 'Occam2DMesh',
                model_fn = 'Occam2DModel',
                iter_fn='ITER17.iter',
                ).items(
                    )
}
# prefix output file
filename = 'k8' #'ybkro'
savepath = 'data/saveGS'
#  Maximum depth investigation  for CSAMT,
# if not provided , will set to 1km 
DOI = '1km'                 #  can be float like 1000 = 1km 

# out put file either
# "m" or "km". Default is "meter"
scale = None

# write file to negative depth:
    # set to True if you want to keep your depth as positive value 
to_negative =True                        

# if elevation is provided set it 
elevation =None 

# Main parameters : value must be less or eagl of DOI.
# if step =DOi WILL return your Occam model files.
# float value. Step to cut out data and to force the calcualted rho
# to match the reference data as input resistivities  
# if not provided the step will be 20% of D0I
STEP_DESCENT = 200.         
# Truth resistivities values otained
 #on the sites or from other companies
#COMPULSORY parameter  

INPUT_RESISTIVITIES = [66,70, 180, 1000, 3000, 10000, 20000] 
# Truth layer names if given must match the input resistivities 

INPUT_LAYERS = ['river water', 'fracture zone' ,
                'granite ', 'Most Weathered', 'Less Weathered'] # 

# ---------------Read with Iter2DAT FILE ------------------
# see Occam_module Iter2Dat to see what the file is . 
# to run this file set the oc2_data to empty dict. 
#x, y, z model file 
i2dkws={'iter2dat_fn ': 'data/iter2dat/iter17.2412.dat',
        'bln_file' :'data/iter2dat/iter17.2412.bln'}
#------Create geodrill object ----------------

geo_obj = Geodrill(  input_resistivities=INPUT_RESISTIVITIES, 
                            input_layers =INPUT_LAYERS ,
                            step_descent =STEP_DESCENT,
                            doi =DOI, 
                            # **i2dkws,
                            **oc2d_data,
                            )    
geo_obj.to_golden_software(filename =filename , 
                               elevation =elevation, 
                               savepath = savepath, 
                               scale =scale, 
                               to_negative_depth=to_negative
                               )     
                                  