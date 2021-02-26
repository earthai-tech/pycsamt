# -*- coding: utf-8 -*-
"""

    .script to plot Penetration 2D . Deal with[AVG|EDI|J]
    User can  customize plot by using matplotlib properties in Plot2d class. 
    
Created on Thu Jan 21 21:20:39 2021

@author: @Daniel03
"""

import os 
from viewer.plot import Plot2d


#----path to your EDI|J|AVG 
#--- avg filename 
avgfile = 'K1.AVG' # can plot directly from abg path using : 
    # path_to_avgfile =os.path.join(os.environ["pyCSAMT"], 'data', 'avg', avgfile)

# path to avgfile 
pathfile  = os.path.join(os.environ["pyCSAMT"], 'data', 'j') # can change your path 
#---If pathfile point to AVG FILE . bring  also profile file . 
profile_stn = None 

#---depth to image : ImageDepth units is in Meter (m)
imageDepth = 2000               # 2000 means 2km . user can use the str ["2km"|2000m"]

#plot style : can be [pycolormesh|imshow]
plotStyle ='imshow'              # default is pcolormesh
#Customize your stations name by bring a new list of station names . 
rename_station =None 

#---ceate 2D plot obj 
plot2d_obj= Plot2d()

plot2d_obj.penetration2D(fn = pathfile, 
                         profile_fn= profile_stn ,
                             plot_style=plotStyle,  
                             doi=imageDepth , 
                             )