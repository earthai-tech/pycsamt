# -*- coding: utf-8 -*-
"""
    . Script to plotskin depth . 
        Deal with [EDI|AVG|J]. AVG is Zong Engineering file . If you bring Zonge engenering file, 
        you can also add its station file know as profile file  (*.stn)
        
Created on Wed Jan 20 15:00:33 2021

@author: @Daniel03
"""
import os 
from viewer.plot import Plot1d 


#-----> path to your files , may be ED|J|OR AVG file 
filespath = os.path.join(os.environ ['pyCSAMT'], 'data', 'edi')

#--> profile station file . Necessary if you bring AVG file 
profile_stn =None
# savefigure path 
savefig =None 
# -- > selected frequencies 
selected_frequencies =[ 1024, 3010, 8000 ]            # for single frequency , you dont need to put on list 

#rotate station name in angle 
rotstn =90 
#---> bring a new station-name on list with the same lenght with the default one. 
new_station_names = None                        # if the defalut station name is not convenient user can change the name of station 
#--- > figure orientation 
orient = 'landscape'

#-- define plot_obj 

plot_1d_obj =Plot1d()
plot1d_depth = plot_1d_obj.penetrated1D(fn =filespath ,
                                        profile_fn= profile_stn, 
                                        selected_frequency =selected_frequencies, 
                                        rename_station =new_station_names, 
                                        rotate_station_names = rotstn , 
                                        orientation=orient, 
                                        savefigure =savefig)

