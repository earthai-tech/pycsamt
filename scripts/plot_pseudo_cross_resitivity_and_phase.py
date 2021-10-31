# -*- coding: utf-8 -*-
"""
    .Script to plot pseudo-cross-section of apparent resistivity and phase.
    Deals with[AVG|J|EDI] file. *avg file is Zonge Engineering file. If brought
    if , add the station profile file *stn.
    Explore to customize  your plot.
    
Created on Fri Jan 22 18:44:21 2021

@author: @Daniel03

"""
import os 
from pycsamt.viewer.plot import Plot2d

#path file 

pathfile =  os.path.join(os.environ["pyCSAMT"],'data','j' )
#pathfile = r'C:\Users\Administrator\Desktop\ThesisImp\edis\_special_K6_edi\k6_AMA'
#pathfile= os.path.join(os.environ["pyCSAMT"],'data', 'correctedEDI_AMA') # test with corrected edi
# if user user avg file , add the profile fn 
profile_fn =None                # not necessary when use [EDI|J] file. 

#savefigure 
savefigure =None

#delineate resitivity : delineate resistivity contour. Resitivities are not on logarithm scale (ohm m ) . 
contourRes =[1000]                # for multiple contour delineation , put value on list eg: [500, 7000]

#delineate phase : delineate phase contour  
contourPhase =[45]              # can be 45 degree or else :None 
    
# plot style : [pcolormesh |imshow] . Default is pcolormesh 
plotStyle =None
# create objet 
#define contout line style 
contour_lines_style='-'

plot2d_obj = Plot2d()
plot2d_obj.pseudocrossResPhase(fn=pathfile, 
                                profile_fn=profile_fn, 
                                delineate_resistivity=contourRes,
                                delineate_phase=contourPhase,
                                plot_style =plotStyle, 
                                savefig = savefigure, 
                                contour_lines_style=contour_lines_style)


