# -*- coding: utf-8 -*-
"""
 .Console scripts. 
 
    $ python pycsamt plot_pseudolog
    
 `**` before the default keys means to convrt input values into float
    `*` before the default keys means to convrt into interger 
    Should compolsory to specify , if not will let inbput values as str 

    Jokers `**` and `*` can be anywhere else in the string default_kwargs key.
    
Created on Tue Apr  6 18:03:21 2021

@author: @Daniel03

"""

import os 
from pycsamt.viewer.plot import Plot2d
from pycsamt.gui.wrap_console_scripts import wrap_cscripts as wrs


default_kwargs = {'Path to occam2d files - str - ' : None , 
                  'occam2d mesh file name - str - ' : None , 
                  'occam2d model file - str - ': None , 
                  'occam2d data file - str - ' : None, 
                  'occam2d iteration file - str - ' : None, 
                  'occam2d save oudir - str - ' : None,
                  'Station id to plot - str|int -': 1,
                  'Input resistivities  - list - '  : None, 
                  'Input layers  - list - ': None,
                  'Step descent in meter** - float - ' :200., 
                  'Depth of investigation in meter** - float - '  : 1000., 
                  'Scale  in [m|km] - str - ': 'm',
                  'Plot style  - str - ': 'pcolormesh',
                  'Path to iter2dat files - str -':None,
                  'Iter2dat `x,y,z` model filename - str -': None , 
                  'Iter2dat `.bln` filename - str - ': None
                  }

def main(): 
    
    wrs_obj = wrs(default_kwargs =default_kwargs)
    for kwey, kwvalues in wrs_obj.sanitize_kwargs.items():
        if kwey.find('Path to') >=0 : path = kwvalues
        if kwey.find('mesh')>=0 : path_to_occam_mesh = kwvalues
        if kwey.find('data') >= 0 : path_to_occam_data = kwvalues 
        if kwey.find('iteration') >= 0 : path_to_occam_iter = kwvalues 
        if kwey.find('model') >=0: path_to_occam_model = kwvalues
        if kwey.find('save outdir')>=0 : savefigure= kwvalues
        if kwey.find('investigation') >=0 :DOI = kwvalues
        if kwey.find('style') >= 0 : plotStyle= kwvalues 
        if kwey.find('Scale') >= 0 : scale= kwvalues 
        if kwey.find('Station id') >=0 :station_to_plot =kwvalues 
        if kwey.find('Input resistivities') >=0: INPUT_RESISTIVITIES= kwvalues
        if kwey.find('Input layers')>=0 : INPUT_LAYERS= kwvalues
        if kwey.find('descent') >=0: STEP_DESCENT= kwvalues
        if kwey.find('to iter2dat files')>=0 : path2 = kwvalues
        if kwey.find('Iter2dat `x,y,z`') >=0: iter2dat_fn= kwvalues
        if kwey.find('Iter2dat `.bln`')>=0 : bln_file = kwvalues

    # Matplotlib properties to customize plots 
    
    average_curve_color = (0.5, 0.8, 0.)
    sequence_curve_color = 'blue'
    font_dict_site = {'size': 8, 
                      'color':'saddlebrown', 
                      'weight': 'bold', 
                      'style': 'italic'}
    # can use other matplotlib to customize your plot 

    Plot2d(station_label_rotation=None , 
            show_grid= True , 
            font_size =8., 
            lc='r' , 
            fig_size=[5,5], 
            markerfacecolor='k', 
            markeredgecolor='k' ).plot_Pseudolog( station_id= station_to_plot, 
                                    input_resistivities=INPUT_RESISTIVITIES, 
                                    input_layers =INPUT_LAYERS ,
                                    step_descent =STEP_DESCENT,
                                    doi =DOI, 
                                    mesh_fn = os.path.join(path , path_to_occam_mesh),
                                    iter_fn = os.path.join(path , path_to_occam_iter), 
                                    model_fn =os.path.join(path, path_to_occam_model) , 
                                    data_fn =os.path.join(path, path_to_occam_data ),
                                    plot_style= plotStyle ,
                                    scale =scale, 
                                    iter2dat_fn = os.path.join(path2 , iter2dat_fn),
                                    bln_fn = os.path.join(path2 , bln_file),
                                    savefig =savefigure,
                                     lc_AD_curves= (average_curve_color, sequence_curve_color),
    
                                    font_dict_site=font_dict_site, 
                                )

    
if __name__=='__main__':

    main()
