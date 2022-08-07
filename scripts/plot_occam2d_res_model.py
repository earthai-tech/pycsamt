# -*- coding: utf-8 -*-
"""
    .Script to plot Occam Resistivity Model .
    More arguments can be set to customize the model plot 
       
Created on Fri Feb  5 21:53:08 2021
@author: @Daniel03
"""
import os  

from pycsamt.viewer.plot import Plot2d

# path to OCCAM 2D folder 

path = r'C:\Users\Daniel\Desktop\Data\AMT\E1\oci_20m' #'data/occam2d'#
path =r'C:\Users\Daniel\Desktop\Data\AMT\E1\oci_res100e20'

#savefigure 
savefigure =r'C:\Users\Daniel\Desktop\Data\AMT\E1\model2im.png' #None                   

# scale the output data 
scale= None                # if None : default is "m" .can be [m|km]

# imaging depth : Maximum depth investigation 
doi = '500m'                 #  can be float like 1000 = 1km 

#plot style 
plotStyle =     "imshow" # "pcolormesh"#  # if None Default is 'imshow', can be 
                            #["pcolormesh"]

#-----OCCAM 2D output data files ---------
# path to occam Data file
path_to_occam_data='OccamDataFile.dat'

# path to occam Mesh file 
path_to_occam_mesh = 'Occam2DMesh'

# path to occam Model file 
path_to_occam_model = 'Occam2DModel'

# path to Occam Iteration file 
path_to_occam_iter='ITER26.iter'#17 #

figsize =[10,4]
# call plot obj 

curve =[ 1  ]
plot2d_obj = Plot2d(fig_size =figsize, station_label_rotation=90, ms=.75, fig_dpi =600)
plot2d_obj.plot_occam2dModel(mesh_fn=os.path.join(path, path_to_occam_mesh), 
                    iter_fn = os.path.join(path, path_to_occam_iter), 
                    model_fn =os.path.join(path, path_to_occam_model ), 
                    data_fn =os.path.join(path, path_to_occam_data), 
                    doi= doi, 
                    savefig =savefigure, 
                    plot_style =plotStyle, 
                    #delineate_rho= curve,
                    #show_contour = True, 
                    )


