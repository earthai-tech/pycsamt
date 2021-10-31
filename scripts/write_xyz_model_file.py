# -*- coding: utf-8 -*-
"""
    . Script to write YangBo Model data file
        called 'x,y,z' *.data file for post-processing.
        Script use output Occam 2D files such as Model, Iteration, Mesh and data file .

Created on Fri Feb  5 12:45:05 2021

@author: @Daniel03

"""
import os 
from pycsamt.modeling.occam2d import Iter2Dat as i2d


# path to OCCAM 2D folder 
#path =os.path.join(os.environ ['pyCSAMT'], 'data', 'occam2D')
path ='data/occam2D'
#savepath folder 
#savepath =None              # if None , will create a folder to hold differents output files

savepath =None

outputfilename =None        # if None , will create automatically 
# scale the output data 
scale_output =None          # if None : default is "km" .can be [m|km]

# imaging depth : Maximum depth investigation 
doi = '1km'                 #  can be float like 1000 = 1km 

#set elevation if you need 
elevation =None             # provided elevation if you need on list or array_like.

#-----Bring OCCAM data files ---------
# path to occam Data file
path_to_occam_data='OccamDataFile.dat'

# path to occam Mesh file 
path_to_occam_mesh = 'Occam2DMesh'

# path to occam Model file 
path_to_occam_model = 'Occam2DModel'
# path to Occam Iteration file 
path_to_occam_iter='ITER17.iter'


# call Iter2Dat object 

occam_iter2dat_obj =i2d(mesh_fn=os.path.join(path, path_to_occam_mesh), 
                    iter_fn = os.path.join(path, path_to_occam_iter), 
                    model_fn =os.path.join(path, path_to_occam_model ), 
                    data_fn =os.path.join(path, path_to_occam_data))
occam_iter2dat_obj.write_iter2dat_file(filename =outputfilename,
                                       scale=scale_output, 
                                       doi=doi, 
                                       elevation =elevation, 
                                       savepath=savepath)
                
