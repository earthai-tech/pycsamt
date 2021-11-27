# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 10:19:19 2021
    .Console scripts. 
    $ python pycsamt write_avg2edi 
    
    .. note:: Try to use dict_kwargs to collect input arguments values, 
            as well as defaults values 

"""

from pycsamt.ff.core import avg 
from pycsamt.gui.wrap_console_scripts import wrap_cscripts as wrs

default_kwargs = {'Path to avg file - type Path like' : None , 
                  'Save ouput edipath - type Path_like': None , 
                  'Add zonge station profile (.stn) - os.path' : None, 
                  'Name of filter - str' : 'tma', 
                  'Number of filter points or dipoles * - int': 7,
                  'number of skin depth * - int' :3, 
                  'Reference frequency** - float' : None  
                  }

def main(): 
    wrs_obj = wrs(default_kwargs =default_kwargs)
    for kwey, kwvalues in wrs_obj.sanitize_kwargs.items():
        if kwey.find('avg') >=0 : avgpath = kwvalues
        if kwey.find('frequency')>=0 : reference_frequency = kwvalues
        if kwey.find('edipath') >= 0 : save_edipath = kwvalues 
        if kwey.find('profile')>=0 : profile_fn = kwvalues 
        if kwey.find('points or dipoles') >=0: number_of_points = kwvalues
        if kwey.find('skin depth')>=0 : number_of_skin_depth = kwvalues
        if kwey.find('Name of filter') >=0 :add_filter = kwvalues
    

    avg.Avg().avg_to_edifile(data_fn= avgpath , 
                        profile_fn = profile_fn, 
                        savepath =save_edipath, 
                        reference_frequency= reference_frequency, 
                        apply_filter=add_filter, 
                        number_of_points=number_of_points, 
                        number_of_skin_depth=number_of_skin_depth
                        ) 
    
if __name__=='__main__':

    main()

    
    
    
    
    
    
    