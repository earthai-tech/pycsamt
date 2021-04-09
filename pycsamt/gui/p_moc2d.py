# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 17:43:47 2021
    .Console scripts. 
    
    $ python pycsamt plot_model_oc2d
    
 `**` before the default keys means to convrt input values into `float`
    `*` before the default keys means to convrt into `integer`. 
    Should compolsory to specify , if not will let inbput values as str 

    Jokers `**` and `*` can be anywhere else in the string default_kwargs key.
    
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
                  'Depth of investigation in meter** - float - '  : 1000., 
                  'Plot style  - str - ': 'pcolormesh',
                  'See groundwater report - bool - ' :False, 
                  'Change figure size - list -':[9,9],
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
        if kwey.find('investigation') >=0 :doi = kwvalues
        if kwey.find('style') >= 0 : plotStyle= kwvalues 
        if kwey.find('report') >=0: see_report= kwvalues
        if kwey.find('figure')>=0 : figsize = kwvalues

        
    plot2d_obj = Plot2d(fig_size =figsize)#,fig_dpi = 600. )
    plot2d_obj.plot_occam2dModel(mesh_fn=os.path.join(path, path_to_occam_mesh), 
                        iter_fn = os.path.join(path, path_to_occam_iter), 
                        model_fn =os.path.join(path, path_to_occam_model ), 
                        data_fn =os.path.join(path, path_to_occam_data), 
                        doi= doi, 
                        savefig =savefigure, 
                        plot_style =plotStyle, 
                        show_report = see_report )
    
if __name__=='__main__':

    main()
