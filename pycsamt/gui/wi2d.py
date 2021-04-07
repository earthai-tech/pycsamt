# -*- coding: utf-8 -*-
"""
.Console scripts. 
    $ python pycsamt write_occam2golden
    
    write model occam2d file to write YangBo Model data file
        called 'x,y,z' *.data file for post-processing.
        Script use output Occam 2D files such as Model, Iteration, Mesh and data file .
        
Created on Tue Apr  6 19:56:09 2021

@author: @Daniel03
"""
import os 
from pycsamt.modeling.occam2d import Iter2Dat as i2d
from pycsamt.gui.wrap_console_scripts import wrap_cscripts as wrs


default_kwargs = {'Path to occam2d files - str - ' : None , 
                  'Occam2d mesh file name - str - ' : None , 
                  'Occam2d model file - str - ': None , 
                  'Occam2d data file - str - ' : None, 
                  'Occam2d iteration file - str - ' : None, 
                  'Depth of investigation in meter** - float - ': 1000., 
                  'Scale  in [m|km] - str - ': 'm',
                  'Path to iter2dat files - str -':None,
                  'Output filename - str -': None , 
                  'Output negative depth - bool -' :True,
                  }

def main(): 
    """
    Can give input_kwargs to test consoles scripts.
    
    :Example:
        
        >>> input_kwargs ={'Path to occam2d files': r'C:/Users\Administrator\Desktop\ThesisImp\'\
        ...                       'occam2D\invers+files\inver_res\K1', 
        ...            'data':'OccamDataFile.dat',
        ...            'mesh' : 'Occam2DMesh',
        ...            'model' : 'Occam2DModel',
        ...            'iter':'ITER17.iter' }      
        >>> wrs_obj = wrs(default_kwargs =default_kwargs, input_kwargs =input_kwargs)
    
    """

    wrs_obj = wrs(default_kwargs =default_kwargs)
    
    for kwey, kwvalues in wrs_obj.sanitize_kwargs.items():
        if kwey.find('Path to occam2d') >=0 : path = kwvalues
        if kwey.find('mesh')>=0 : path_to_occam_mesh = kwvalues
        if kwey.find('data') >= 0 : path_to_occam_data = kwvalues 
        if kwey.find('iteration') >= 0 : path_to_occam_iter = kwvalues 
        if kwey.find('model') >=0: path_to_occam_model = kwvalues
        if kwey.find('investigation') >=0 :DOI = kwvalues
        if kwey.find('Scale') >= 0 : scale_output= kwvalues 
        if kwey.find('Output filename')>=0 : outputfilename = kwvalues

    occam_iter2dat_obj =i2d(mesh_fn=os.path.join(path, path_to_occam_mesh), 
                        iter_fn = os.path.join(path, path_to_occam_iter), 
                        model_fn =os.path.join(path, path_to_occam_model ), 
                        data_fn =os.path.join(path, path_to_occam_data))
    occam_iter2dat_obj.write_iter2dat_file(filename =outputfilename,
                                            scale=scale_output, 
                                            doi=DOI)
    
if __name__=='__main__':
    main()