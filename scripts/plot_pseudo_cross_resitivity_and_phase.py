# -*- coding: utf-8 -*-
"""
    .Script to plot pseudo-cross-section of apparent
    resistivity and phase.Deals with[AVG|J|EDI] file. *avg file 
    is Zonge Engineering file. If brought
    if , add the station profile file *stn.
    Explore to customize  your plot.
    
Created on Fri Jan 22 18:44:21 2021
@author:K.L ~ @Daniel03

"""

from pycsamt.viewer.plot import Plot2d

#path file 
pathfile =  r"C:\Users\Daniel\Desktop\Data\AMT\E1\edi_test" 
pathfile = r'C:\Users\Daniel\Desktop\Data\AMT\E1\edi_rw'#ediout_batch1"#'data/K1_edi'
# if user user avg file , add the profile fn 
profile_fn =None                # not necessary when use [EDI|J] file. 

#savefigure 
savefigure =None

#delineate resitivity : delineate resistivity contour.
#Resitivities are not on logarithm scale (ohm m ) . 
# for multiple contour delineation , put value on list eg: [500, 7000]
contourRes =None #[1000]               

#delineate phase : delineate phase contour  
contourPhase = None              # can be 45 degree or else :None 
    
# plot style : [pcolormesh |imshow] . Default is pcolormesh 
plotStyle = 'imshow'#'pcolormesh'
# create objet 
#define contout line style 
contour_lines_style='-'

plot2d_obj = Plot2d(plot_comp = 'yx', fig_size =(20, 9), #(12,6),
                    station_names ='edi',  # showing edifinesnames as stations, 
                    station_label_rotation =90, 
                    cmap =  'jet',
                    
                    )
plot2d_obj.pseudocrossResPhase(fn=pathfile, 
                                profile_fn=profile_fn, 
                                delineate_resistivity=contourRes,
                                delineate_phase=contourPhase,
                                plot_style =plotStyle, 
                                savefig = savefigure, 
                                contour_lines_style=contour_lines_style,
                               )

