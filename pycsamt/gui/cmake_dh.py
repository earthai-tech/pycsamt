# -*- coding: utf-8 -*-
"""
    .Console scripts. 
    $ python pycsamt write_drillhole 
    
    write drill hole from location drill parser file . To see how parser file is 
    organized , consult ~/ 'drill_example_files' directory.
    exemple of parser file is 'nbleDH.csv'
    to write drillhole individually either `Collar`,`Geology`,
    `Sample`,`Elevation`, or `Azimuth` and to also add elevation or azimuth
    please use the scripts in ~/scripts directory. 
    
    default`*` means write all data. 
    
Created on Tue Apr  6 21:33:23 2021

@author: @Daniel03

"""
import os 
from pycsamt.geodrill.geocore import Drill 
from pycsamt.gui.wrap_console_scripts import wrap_cscripts as wrs

default_kwargs = {'Path to parser drill files - str - ' : None , 
                  'Name of parser file  - str - '  : False, 
                  'Type of data to write - str - ':'*',
                  'Save borehole outdir  - str - ':None,
                  }
      
def main(): 
    
    wrs_obj = wrs(default_kwargs =default_kwargs)
    for kwey, kwvalues in wrs_obj.sanitize_kwargs.items():
        if kwey.find('Path to parser drill') >=0 : path_to_parser_files = kwvalues
        if kwey.find('Name of parser file ')>=0 : parser_file = kwvalues
        if kwey.find('Type of data to write') >= 0 : kind_of_data2output = kwvalues 
        if kwey.find('Save borehole outdir') >= 0 : savepath = kwvalues 
 

    borehole_obj = Drill (well_filename= os.path.join(path_to_parser_files, parser_file)  , 
                       auto= True)
    borehole_obj._collar()
    borehole_obj.dhGeology()
    borehole_obj.dhSample()
    borehole_obj.dhSurveyElevAz(add_azimuth=None, 
                                add_elevation=None)
    borehole_obj.writeDHData(data2write=kind_of_data2output, 
                             savepath  = savepath )
if __name__ == '__main__': 
    main()
