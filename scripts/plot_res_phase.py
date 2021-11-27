# -*- coding: utf-8 -*-
"""
    . Script to plot the resistivity and phase at each station.
    Deal with [EDI|J|AVG].*AVG is Zonge Engineering file.If provided, 
    please add your profile (*.stn) file.
    
Created on Thu Jan 21 08:57:20 2021
"""

from pycsamt.viewer.plot import Plot1d

#--- path to files [EDI|J|AVG] 
path ='data/avg/K1.AVG' #  can change the path

#---> profile *stn file                   
profile_stn ='data/avg/K1.stn' 
#---selected station to plot 
station_reference =[ 24]             # can be only string  like ["S04,"S25"]
#---> rename your station 
new_stations = None 
#--> see error bar 
showError_bar=True


#create plotObj 
plot1D_obj = Plot1d()
plot1D_obj.plot_freqVSRhoPhase(fn =path , 
                               profile_fn=profile_stn, 
                               station_id =station_reference,
                               rename_stations= new_stations,
                               show_error=showError_bar )
