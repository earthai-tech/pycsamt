# -*- coding: utf-8 -*-
"""

    .script to plot Penetration 2D . Deal with[AVG|EDI|J]
    User can  customize plot by using matplotlib properties in Plot2d class. 
    
Created on Thu Jan 21 21:20:39 2021

@author: @Daniel03
"""

import os 
from pycsamt.viewer.plot import Plot2d


#----path to your EDI|J|AVG 
#--- avg filename 
avgfile = 'K1.AVG' # can plot directly from avg using path path below : 
#pathfile =os.path.join(os.environ["pyCSAMT"], 'data', 'avg', avgfile)

# path to EDifile or jfile 
pathfile  = os.path.join(os.environ["pyCSAMT"], 'data', 'j') 

#pathfile  = os.path.join(os.environ["pyCSAMT"],'data', 'correctedEDI') # test with corrected edi

#save figure 
savefigure = None # r'C:\Users\Administrator\Desktop\ThesisImp\plots\penetration2D\K6_2.png'

#---If path points to avg file . bring  also profile file . 
#profile_file = 'K1.stn' 

#profile_stn =os.path.join(os.environ["pyCSAMT"], 'data', 'avg', profile_file)
# otherwise set to None 
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
                            savefig=savefigure )