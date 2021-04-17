# -*- coding: utf-8 -*-
"""
    .Script to plot RMS .from occam 2D or other software. 
    Browse to see others params to customize your plot.
    
Created on Sat Jan 23 18:31:02 2021

@author: @Daniel03

"""
import os 
from pycsamt.viewer.plot import Plot1d

#---path to logile if exists 
occam2dLogfile = 'LogFile.logfile'          # Occam2D logfile 
logPath = os.path.join(r'C:\Users\Administrator\Desktop\ThesisImp\occam2D\invers+files\inver_res\K9', 
                       occam2dLogfile)
# logPath = os.path.join(os.environ['pyCSAMT'], 
#                        'data', 'occam2D', occam2dLogfile)

#R.M.S Target value 
RMS_target = 1.

savefigure = r'C:\Users\Administrator\Desktop\ThesisImp\plots\plotRMS\K9.png'
#--- set to True to show_target line 
showTargetLine =True 
# set to True to see Grid 
showGrid =False

#---create obj -----------

plot_1d_obj =Plot1d()
plot_1d_obj.plotRMS(fn =logPath, 
                    target=RMS_target, 
                    show_grid =showGrid ,
                    show_target_line = showTargetLine , savefig = savefigure )
