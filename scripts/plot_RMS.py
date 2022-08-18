# -*- coding: utf-8 -*-
"""
    .Script to plot RMS .from occam 2D or other software. 
    Browse to see others params to customize your plot.
    
Created on Sat Jan 23 18:31:02 2021

@author:K.L ~ @Daniel03

"""
import os 
from pycsamt.viewer.plot import Plot1d

#---path to logile if exists 
occampath = r'C:\Users\Daniel\Desktop\Data\AMT\E1\oci_20m' #'data/occam2D'
occam2dLogfile = 'LogFile.logfile'          # Occam2D logfile 
logPath =os.path.join(occampath, occam2dLogfile )

#R.M.S Target value 
RMS_target = 1.

savefigure = None
#--- set to True to show_target line 
showTargetLine =True 
# set to True to see Grid 
showGrid =False

#---create obj -----------

plot_1d_obj =Plot1d()
plot_1d_obj.plotRMS(fn =logPath, 
                    target=RMS_target, 
                    show_grid =showGrid ,
                    show_target_line = showTargetLine , 
                    savefig = savefigure )

# display figure in non-interactive mode
# see  http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.show
import matplotlib as mpl
mpl.pyplot.show() 
