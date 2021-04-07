# -*- coding: utf-8 -*-
"""
    .Console scripts. 
    $ python pycsamt occam2d_build_in 
    
    write  occam2d input files with edifiles.
    expected files are ('OccamDataFile.dat', 'Occam2DMesh', 'Occam2DModel', 
                        'ITERxx.iter', 'Startup')
    `xx` is= number of iteration.
    
Created on Tue Apr  6 20:53:41 2021

@author: @Daniel03
"""
from  pycsamt.modeling.occam2d import occam2d_write
from pycsamt.gui.wrap_console_scripts import wrap_cscripts as wrs

default_kwargs = {'Path to edifiles - str - ' : None , 
                  'Interpolate frequencies  - bool- '  : False, 
                  'Interpolate frequencies range (freqmin, freqmax, number of frequency) - tuple - ':(-1, 4, 17),
                  'Number of model layers*  - int - ': 31,
                  'Investigation depth (Z target) in meter** - float - ' :1100., 
                  'Z bottom in meter**  - float - ': 5000.,
                  'First layer thickness  in meter**  - float - ': 5.,
                  'Start model resistivity in log10**  - float -':1.,
                  'Number of iteration to run * - int -': 100 , 
                  'Geoelectric strike in degree E-N** - float - ': 0.,
                  'Occam2d model mode - str -' :'6',
                  'TM phase error in % - float** -': 20., 
                  'TM resistivity error in % **- float -': 10., 
                  'TE phase error in % **- float -': 20., 
                  'TE resistivity error in %** - float -': 10., 
                  'Model cell width in meter **- float - ':5, 
                  'Horizontal x pad multiplier** -float -':1.7, 
                  'Brick trigger **- float -': 1.12, 
                  
                  }
      
def main(): 
    
    wrs_obj = wrs(default_kwargs =default_kwargs)
    for kwey, kwvalues in wrs_obj.sanitize_kwargs.items():
        if kwey.find('Path to edifiles') >=0 : edipath = kwvalues
        if kwey.find('frequencies  - bool')>=0 : interpolate_frequency = kwvalues
        if kwey.find('frequencies range') >= 0 : frequency_interpolate_range = kwvalues 
        if kwey.find('of model layers') >= 0 : number_of_model_layers = kwvalues 
        if kwey.find('depth (z target)') >=0: expected_investigation_depth_or_z_target = kwvalues
        if kwey.find('Z bottom in')>=0 : exaggerate_vertical_z_bottom  = kwvalues
        if kwey.find('First layer thickness') >=0 :model_first_layer_thickness_z1 = kwvalues
        if kwey.find('Start model resistivity') >= 0 : starting_res_model_value= kwvalues 
        if kwey.find('iteration to run') >= 0 : number_of_iteration_to_run= kwvalues 
        if kwey.find('strike in degree') >= 0 : geoelectric_strike= kwvalues 
        if kwey.find('Occam2d model mode') >=0: occam_mode= kwvalues
        if kwey.find('TM phase')>=0 : TM_phase_error = kwvalues
        if kwey.find('TM resistivity') >=0: TM_res_error= kwvalues
        # if kwey.find('TE phase')>=0 : path2 = kwvalues
        # if kwey.find('TE resistivity') >=0: iter2dat_fn= kwvalues
        if kwey.find('cell width in')>=0 : model_cell_width = kwvalues
        if kwey.find('x pad multiplier')>=0 : horizontal_node_x_pad_multiplier = kwvalues
        if kwey.find('Brick trigger  ')>=0 : brick_trigger = kwvalues

  

    occam2d_write.buildingInputfiles(edi_fn =edipath,
                               geoelectric_strike= geoelectric_strike,
                               interpolate_freq= interpolate_frequency, 
                               intp_freq_logspace =frequency_interpolate_range, 
                                iteration_to_run= number_of_iteration_to_run, 
                                resistivity_start = starting_res_model_value, 
                                res_tm_err= TM_res_error, 
                                phase_tm_err= TM_phase_error , 
                                occam_mode= occam_mode, 
                                n_layers = number_of_model_layers , 
                                cell_width = model_cell_width , 
                                x_pad_multiplier = horizontal_node_x_pad_multiplier, 
                                trigger =brick_trigger, 
                                z_bottom =exaggerate_vertical_z_bottom , 
                                z1_layer =model_first_layer_thickness_z1, 
                                z_target = expected_investigation_depth_or_z_target, 
                            )
if __name__ == '__main__': 
    main()