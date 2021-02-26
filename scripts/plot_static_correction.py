# -*- coding: utf-8 -*-
"""
     .Script to plot static correction . Deal with [EDI|J|AVG]files. If user provided 
     Zonge Engineering file "*avg", it may provided also profile station file '*.stn'
     Add filter for corrected apparent resistivities AMA , TMA or FLMA. 
     If reference value set is not in the freq array , it will be interpolated 
     to a frequency value with highest clean data . 
     
     
Created on Tue Jan 19 16:57:08 2021

@author: @Daniel03

"""
import os 
from viewer.plot import Plot1d 

#--- > path to your file 

path_to_file =os.path.join(os.environ['pyCSAMT'],
                           'data', 'edi' )

# stn station profile file 
profile_stn =None 

#---> Filter points 
TMApoints = 7
#reference frequency at that station  . for multipleplot , add reference frequency into a list eg : [8, 1024, 2012]

reference_frequency =8000.
# can fill between to delineate the correction effect 
fillBetween =True 
#---call object --- 

rhoplot_obj = Plot1d().plot_static_correction(data_fn =path_to_file , 
                                              profile_fn=profile_stn, 
                                              frequency_id= reference_frequency,
                                              number_of_TMA_points =TMApoints, 
                                              fill_between =fillBetween )



