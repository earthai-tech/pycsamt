# -*- coding: utf-8 -*-
"""
     .Script to plot static correction . Deal with [EDI|J|AVG]files. If 
     Zonge Engineering file "*avg" is given, user should  provided 
     also profile station file '*.stn'. Availbale Filters to correct apparent 
     resistivities are  `AMA` , `TMA` and `FLMA`. 
     AMA: Adapatative moving average based on the idea of Torres-Verdin, 
     FLMA: fixed length dipole moving average ,
     TMA: Trimming moving average most used by Zonge Engineering company.
     If reference value set is not in frequency range, it shoul be interpolated.
     reference frequency value is highest frequency with clean data. 
     
     
Created on Tue Jan 19 16:57:08 2021

@author: @Daniel03

"""
import os 
# from pycsamt.viewer import plot 
from pycsamt.viewer.plot import Plot1d 

#--- > path to your file 

path_to_file = os.path.join(os.environ['pyCSAMT'],'data','avg','K1.AVG')
path_to_file= r'F:\ThesisImp\edis\K1_edi'

#path_to_file = r'C:\Users\Administrator\Desktop\ThesisImp\edis\K9_edi' 
#path_to_file =os.path.join(os.environ['pyCSAMT'],'data', 'edi' )
                           
#save figure 
savefigure = r'C:\Users\Administrator\Desktop\ThesisImp\plots\correction\K9.png'
savefigure=None 
# stn station profile file 
#profile_stn = os.path.join(os.environ['pyCSAMT'],'data','avg','K1.stn')
profile_stn =None  #  

FILTER='ama'                   # Can be `tma` or `flma`

#fipole length in meter if `flma filter is provided 
dipole_length =50.

# When plot AMA filter add number of filter : default is 1 , can be 1 to 10
number_of_skinDepth=3.

#---> Filter points 
FILERpoints = 7.
#reference frequency at that station  . for multipleplot , 
#add reference frequency into a list eg : [8, 1024, 2012]
reference_frequency =8192.

# can fill between to delineate the correction effect 
fillBetween =False 

#customize plots 
fill_between_color= 'thistle'
tma_color= 'blue'
flma_color='aqua'
   
marker=  'x'
ms =  3.
ls =  '-'
markeredgecolor=  'k'
fs = 2.
lw= 1.5

#---call object --- 
rhoplot_obj = Plot1d(mstyle= marker, 
                     ms=ms,
                     fs=fs, 
                     lw=lw, 
                     markeredgecolor= markeredgecolor,
                     )

rhoplot_obj.plot_static_correction(data_fn =path_to_file , 
                                    profile_fn=profile_stn, 
                                    frequency_id= reference_frequency,
                                    number_of_points=FILERpoints, 
                                    fill_between =fillBetween,
                                    ADD_FILTER =FILTER, 
                                    dipole_length=dipole_length,
                                    fill_between_color= fill_between_color, 
                                    tma_color=tma_color, 
                                    flma_color= flma_color,
                                    number_of_skin_depth = number_of_skinDepth, 
                                    savefig =savefigure)



