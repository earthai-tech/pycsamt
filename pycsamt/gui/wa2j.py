# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 17:01:00 2021
    .Console scripts. 
    $ python pycsamt write_avg2j
    
    .. note:: If  defaults key change . Specify a shortcut name(key) to let 
        sanitize dict to fing exactly the input arguments values. For instance 
        default_kwargs['Path to avg file'] --> shortcut key =='to avg'
        shortcut name should the unique `str` name that could be found in 
        sanitize__dict_keys(). Indeed sanitize dict is  the combinaison of 
        default_kwargs.keys + default_keyward.values. 
        like  `Path to avg file - str - [None]`
        ....


@author: @Daniel03

"""
from pycsamt.ff.core import avg 
from pycsamt.gui.wrap_console_scripts import wrap_cscripts as wrs

default_kwargs = {'Path to avg file - str - ' : None , 
                  'Save  j outdir - str -': None , 
                  'Add Zonge Station profile (.stn) - str - ' : None, 
                  'j extension - str - '  : '.dat', 
                  ' Survey name - str -': None,
                  'Write avg file infos * - bool - ' :False, 
                  'UTM zone - str  -' : '49N',  
                  }

def main(): 
    wrs_obj = wrs(default_kwargs =default_kwargs)
    for kwey, kwvalues in wrs_obj.sanitize_kwargs.items():
        if kwey.find('to avg') >=0 : path_to_avg = kwvalues
        if kwey.find('outdir')>=0 : savepath = kwvalues
        if kwey.find('profile') >= 0 : path_to_station_file = kwvalues 
        if kwey.find('extension')>=0 : j_extension = kwvalues 
        if kwey.find('Survey name') >=0: surveyName = kwvalues
        if kwey.find('file infos')>=0 : write_avg_file_infos = kwvalues
        if kwey.find('UTM') >=0 :utmZone = kwvalues
    
    avg.Avg().avg_to_jfile(avg_data_fn=path_to_avg, 
                                profile_fn=path_to_station_file,
                                j_extension=j_extension,
                                savepath=savepath  ,
                                writeInfos=write_avg_file_infos, 
                                utm_zone =utmZone , 
                                survey_name =surveyName)
    
if __name__=='__main__':

    main()
