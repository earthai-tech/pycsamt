# -*- coding: utf-8 -*-
"""

    .script to plot Penetration 2D . Deal with[AVG|EDI|J]
    User can  customize plot by using matplotlib properties 
    in Plot2d class. 
    
Created on Thu Jan 21 21:20:39 2021
@author: @Daniel03
"""

from pycsamt.viewer.plot import Plot2d

#----path to your EDI|J|AVG 
#--- avg filename 
# pathfile = 'data/avg/K1.AVG' 

# path to EDifile or jfile 
pathfile  = 'data/j'

#save figure 
savefigure = None 

#---If path points to avg file . bring  also the profile file *.stn. 
# profile_stn = 'data/avg/K1.stn' 
# otherwise set to None 
profile_stn = None          

#---depth to image : ImageDepth units is in Meter (m)
# 2000 means 2km . user can use the str ["2km"|2000m"]
imageDepth = 2000               

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
                        savefig=savefigure )
