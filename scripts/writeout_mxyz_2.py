# -*- coding: utf-8 -*-
"""
    . Script to write YangBo Model data file
        called 'x,y,z' *.data file for post-processing.
        Script use output Occam 2D files such as Model,
        Iteration, Mesh and data file .

Created on Fri Feb  5 12:45:05 2021

"""
from pycsamt.modeling.occam2d import Iter2Dat as i2d


# path to OCCAM 2D folder 
i2dkws = dict(
    data_fn='data/occam2D/OccamDataFile.dat',
    mesh_fn= 'data/occam2D/Occam2DMesh',
    model_fn= 'data/occam2D/Occam2DModel',
    iter_fn='data/occam2D/ITER17.iter',
)

#savepath folder 
savepath =None          

# give an output file name 
# if None, will create automatically 
outputfilename =None        
# scale the output data 
# if None : default is "km" .can be [m|km]
scale_output =None          

# imaging depth : Maximum depth investigation 
#  can be float like 1000 = 1km 
doi = '1km'                 

#set elevation if you need 
 # provided elevation if you need on list or array_like.
elevation =None            

# call Iter2Dat object 
occam_iter2dat_obj =i2d(**i2dkws, savepath = savepath )
occam_iter2dat_obj.write_iter2dat_file(filename =outputfilename,
                                       scale=scale_output, 
                                       doi=doi, 
                                       elevation =elevation, 
                                       )
                
