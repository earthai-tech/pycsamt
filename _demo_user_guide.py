# -*- coding: utf-8 -*-
"""
Created on Fri Apr 16 15:11:29 2021

@author: @Daniel03
"""

# import required modules
import os 
from pycsamt.ff.core.avg import Avg 
from pycsamt.viewer.plot import Plot1d
from pycsamt.viewer.plot import Plot2d
from pycsamt.ff.core.cs import Profile
from pycsamt.ff.core.cs import Site
# create profile_obj 
profile_obj = Profile(station_profile)
# to get 
profile_obj.east    # get easting coordinate  
profile_obj.north   # get northing coordinates 
profile_obj.lon     # longitude value 
profile_obj.lat     # latitue value 
# to get dipole length in meters 
profile_obj.dipole_length 

profile_obj.stn_interval  # interval between stations 
profile_obj.stn_position  # scaled position of each stations 
profile_obj.azimuth          #  azimuth of profile_line 


profile_obj.Site.stn_name # station id 

# to get single value of latitude and longitude or easing northing at each station 
# call Site obj 
# 
straighten_out_mode ='classic' # can be 'natural/distord' or equidistant
# contrinute value for x,y coordinates hidden. 
adjust__x_utm_coordinates = -300238.702 
adjust__y_utm_coordinates = -2369.252
get_new_station_profile =True 
Plot1d().plot_station_profile(fn = station_profile_file, 
                                 reajust_coordinates=(adjust__x_utm_coordinates,
                                                      adjust__y_utm_coordinates),
                                 straighten_type =straighten_out_mode , 
                                 outputfile =get_new_station_profile, 
                                 savefig=savepath)

from pycsamt.viewer.plot import Plot2d
contouRes = 1000.       # resistivity in ohm.meters 
for path2edi_obj in [
                    'data/edi', 
                    'data/correctedEDI_TMA', # edi corrected with tma filter
                    'data/correctedEDI_FLMA',  # edi corrected with flma filter 
                    'data/correctedEDI_AMA',  # edi corrected with ama filters 
                    ]:
    Plot2d().pseudocrossResPhase(fn=path2edi_obj,
                                 delineate_resistivity=[1000])
    
from pycsamt.modeling.occam2d.import occam2d_write 
occam2d_write.buildingInputfiles(edi_fn =edipath,
                           geoelectric_strike= 34.,
                           interpolate_freq= True, # interpolate if you want
                           intp_freq_logspace =(-1, 4, 17), 
                            iteration_to_run= 100., 
                            resistivity_start = 1.,  # in log10 resistivity 
                            res_tm_err= 10.,        # in %
                            phase_tm_err= 20. ,     # in %
                            occam_mode= '6',        # see documentation
                            n_layers =31. , 
                            cell_width = 5 , 
                            x_pad_multiplier = 1.7, 
                            trigger =1.12, 
                            z_bottom =5000., 
                            z1_layer =5., 
                            z_target = 1100., 
                        )

#-> first option : use occam2d modelfiles   ---------

oc2d_inversion_kwargs={'data_fn':'OccamDataFile.dat', 
                      'mesh_fn': 'Occam2DMesh', 
                      'model_fn':'Occam2DModel' , 
                      'iter_fn': 'ITER17.iter', 
                      }
# join the occam2d directory with occam inversion files
# to create a realpath for each file. 
for oc2dkey, inversion_file  in oc2d_inversion_kwargs.items(): 
    oc2d_inversion_kwargs[oc2dkey]= os.path.join('data/occam2D', inversion_file)
    
from pycsamt.viewer.plot import Plot1d
RMS_target =1. #  Root Mean-square target , default is 1.0
show_grid =False  #set to let grid to be visible
showTargetLine= True # show rms target line . 
savefigure = os.path.abspath('./test_rmsplot.png') # savepath 
Plot1d().plotRMS(fn ='data/occam2D/logFile.logfile', 
                    target=RMS_target , 
                    show_grid =show_grid,
                    show_target_line = showTargetLine,
                    savefig =savefigure )

from pycsamt.modeling.occam2d import Iter2Dat as i2d
# give an output file name , if None , will create automatically 
outputfilename ='testi2d_area'       
# scale the output data 
scale_output =None          # if None,  default is "km" .can be "m" 
# imaging depth : Maximum depth investigation 
doi = '1km'                 #  can be float like 1000 = 1km 
elevation =None             # provided elevation if you need on list or array_like.
# create i2d or modelxyz object and entering aruments 
occam_iter2dat_obj =i2d(**oc2d_inversion_kwargs, 
                        savepath =savepath)
occam_iter2dat_obj.write_iter2dat_file(filename =outputfilename,
                                       scale=scale_output, 
                                       doi=doi, 
                                       elevation=elevation)


contourRes =None  #  value or resistivity in on ohm m, can be 700 or 100.
showReport = True 
savefigure = os.path.abspath('./test_fwdresp.png') # savepath            
Plot2d().plot_Response(data_fn ='data/occam2D/OccamDataFile.dat', 
                         response_fn = 'data/occam2D/RESP17.resp' , 
                         delineate_resistivity =contourRes   , 
                         show_report =showReport , 
                         savefig =savefigure)


