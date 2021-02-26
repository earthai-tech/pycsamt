# -*- coding: utf-8 -*-
"""
     . Script to plot differents curves from avg file .Deal only with Zonge Engineering file.  
         plot will add some informations of provenance file . 
     
Created on Wed Jan 20 20:23:39 2021

@author: @Daniel03

"""

import os 
from viewer.plot import Plot1d 


#--- avg filename 
avgfile = 'K1.AVG' 

# path to avgfile 
avgPath  = os.path.join(os.environ["pyCSAMT"], 
                        'data','avg',  avgfile) # can change the path 

#see errorbar plot 
see_errobar = True              # set to False if you dont want to see errorbar on the plot 

#station to plot : can reas station name liKE "S00".station number start by one . 
station_id= [1,10, 20]          # can be [S00, S01 , 18 ] etc.. for multiples stations plot , put stations on list
                          
# create plot_obj 
plot_1d_obj =Plot1d()
plotcurves = plot_1d_obj .plot_curves(fn = avgPath  , 
                                      selected_stations=station_id, 
                                      error_bar=see_errobar)