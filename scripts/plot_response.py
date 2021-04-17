# -*- coding: utf-8 -*-
"""
    . Script to plot Occam 2D Forwad response and Residual 
    
Created on Tue Feb  9 07:53:16 2021

@author: @Daniel03
"""

import os 
from pycsamt.viewer.plot import Plot2d


# par=th to occam folder 

path = os.path.join(os.environ ['pyCSAMT'],  'data', 'occam2D' )

# path to Occam response file 
resp= 'RESP17.resp'
# path to occam data file 
data='OccamDataFile.dat'

# delineate _resistivity cure 
delineate= 700.                 # value or resistivity to delineate is on ohm m not in log10 resistivity
                                # for multiple value of contour, use the list like , [500, 700] means 
                                # contour500ohm and 7000 ohm meter 
# path to save figure 
savefigure =None                

# show a report 
see_report =False               # if see report is True , then program will generate a samll report about 
                                # resistivites value in the survey area and propose in the case of groundwater exploration 
                                # an overview of drilling verification point . If your pupose is not a groundater exploration 
                                # you can ignore it . 
# plot style 
plotStyle ='imshow'                 # if None , default is 'imshow' can be 'pcolormesh'.


# create plot_object 
 
plot2d_obj = Plot2d(fig_size = [6,6], 
                    font_size =6, 
                    show_grid =True)       # can use lot of matplotlib properties to customize your plot 

    
plot2d_obj.plot_Response(data_fn =os.path.join(path, data) , 
                         response_fn =  os.path.join(path, resp) , 
                         show_report =see_report  ,
                         delineate_resistivity=delineate,
                         plot_style =plotStyle, 
                         savefig =savefigure)