# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 17:19:54 2021
    .Console scripts. 
    $ python pycsamt corrected_edi
    
@author: @Daniel03

    `**` before the default keys means to convrt input values into float
    `*` before the default keys means to convrt into interger 
    Should compolsory to specify , if not will let inbput values as str 

    Jokers `**` and `*` can be anywhere else in the string default_kwargs key.
    
"""

from pycsamt.ff.processing.corr import shifting
from pycsamt.gui.wrap_console_scripts import wrap_cscripts as wrs

default_kwargs = {'Path to edi file - str - ' : None , 
                  'Save  corrected edi outdir - str - ': None , 
                  'New edi filename - str - ' : None, 
                  'FILTER - str - '  : 'ama', 
                  'Number of filter points* - int - ': 3,
                  'Number of skin depth * - int - ' :3, 
                  'Dipole length in meter** - float -' : 50., 
                  'Reference frequency in Hz** - float - ': None, 
                  'Data type - str - '  : None, 
                  'Reduce rho factor x** - float -': 1.,
                  'Reduce rho factor y** - float -':1.,
                  'Distortion tensor - ndarray(2,2, dtype=real)- ':None, 
                  'Distortion  error tensor - ndarray(2,2, dtype=real) -':None
                  }

def main(): 
    wrs_obj = wrs(default_kwargs =default_kwargs)
    for kwey, kwvalues in wrs_obj.sanitize_kwargs.items():
        if kwey.find('to edi file') >=0 : edipath = kwvalues
        if kwey.find('outdir')>=0 : savepath = kwvalues
        if kwey.find('filename') >= 0 : new_edifilename = kwvalues 
        if kwey.find('FILTER') >= 0 : FILTER = kwvalues 
        if kwey.find('of filter point') >=0: number_of_filter_points = kwvalues
        if kwey.find('of skin depth')>=0 : number_of_skin_depth= kwvalues
        if kwey.find('length') >=0 :dipoleLength = kwvalues
        if kwey.find('frequency in') >= 0 : reference_frequency= kwvalues 
        if kwey.find('type') >=0: datatype= kwvalues
        if kwey.find('rho factor x')>=0 : reduce_res_factor_x = kwvalues
        if kwey.find('rho factor y') >=0 :reduce_res_factor_y = kwvalues
        if kwey.find('Distortion tensor')>=0 : distortion_tensor = kwvalues
        if kwey.find('Distortion  error tensor') >=0 :distortion_err_tensor = kwvalues
        
    
    shifting().write_corrected_edi(data_fn = edipath, 
                             number_of_points =number_of_filter_points,
                             reference_frequency=reference_frequency,
                             number_of_skin_depth=number_of_skin_depth, 
                             dipole_length =dipoleLength, 
                             FILTER=FILTER, 
                             filename = new_edifilename, 
                             datatype =datatype, 
                             reduce_res_factor_x=reduce_res_factor_x, 
                             reduce_res_factor_y = reduce_res_factor_y, 
                             distortion_tensor= distortion_tensor, 
                             distortion_err_tensor = distortion_err_tensor, 
                             savepath =savepath )
    
if __name__=='__main__':

    main()
