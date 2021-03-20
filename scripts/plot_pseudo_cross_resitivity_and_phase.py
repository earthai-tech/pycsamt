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
from viewer.plot import Plot2d

#path file 

pathfile =  os.path.join(os.environ["pyCSAMT"],'data','j' )
# pathfile= os.path.join(os.environ["pyCSAMT"],'data', 'correctedEDI') # test with corrected edi
# if user user avg file , add the profile fn 
profile_fn =None                # not necessary when use [EDI|J] file. 

#delineate resitivity : delineate resistivity contour. Resitivities are not on logarithm scale (ohm m ) . 
contourRes =[1000]                # for multiple contour delineation , put value on list eg: [500, 7000]

#delineate phase : delineate phase contour  
contourPhase =[45]              # can be 45 degree or else :None 
    
# plot style : [pcolormesh |imshow] . Default is pcolormesh 
plotStyle =None
# create objet 

plot2d_obj = Plot2d()
plot2d_obj.pseudocrossResPhase(fn=pathfile, 
                                profile_fn=profile_fn, 
                                delineate_resistivity=contourRes,
                                delineate_phase=contourPhase,
                                plot_style =plotStyle)


