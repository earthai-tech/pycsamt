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

from pycsamt.viewer.plot import Plot2d


# additional geological informations collected 
INPUT_RESISTIVITIES = [66.,70., 180.,
                       1000., 3000., 10000., 20000.] 
INPUT_LAYERS = ['river zone', 'fracture zone' , 
                'granite ', 'Most Weathered', 
                'Less Weathered'] 
STEP_DESCENT =200.      # step descent in meters 
DOI  ='1km'             # investigation depath 
additional_geological_infos={
        'doi':DOI,
        'step_descent': STEP_DESCENT, 
        'input_resistivities' : INPUT_RESISTIVITIES, 
        'input_layers' : INPUT_LAYERS
                    }

from pycsamt.viewer.plot import Plot2d
Plot2d().plot_Pseudolog(station_id = 'S43', # station to visualize
                        **additional_geological_infos, 
                        **oc2d_inversion_kwargs)

i2d_files_kwargs={
                 'iter2dat_fn' : 'data/iter2dat/K1.iter.dat',
                 'bln_fn':'data/iter2dat/K1.bln'
                 }
Plot2d().plot_Pseudolog(station_id = 'S43', # station to visualize
                        **additional_geological_infos, 
                        **i2d_files_kwargs)       
       


from pycsamt.geodrill.geoCore.geodrill import Geodrill 
# create a geological object from geodrill module 
geo_obj = Geodrill(**oc2d_inversion_kwargs , **additional_geological_infos)
# ca export to golden software 
geo_obj.to_golden_software(filename= 'test_area', 
                           to_negative_depth =True , # default output
                           savepath=savepath) 

from pycsamt.geodrill.geoCore.geodrill import Geodrill 
geo_obj.to_oasis_montaj(profile_fn ='data/avg/K1.stn', # can be create if coordinates doesnt exists) 
                           to_negative_depth =True , # default output
                           to_log10=True,   #output resistivity to log10 values 
                           filename ='test_area',
                           savepath=savepath) 


from pycsamt.geodrill.geoCore.geodrill import Geosurface
path_to_oasisfiles ='data/InputOas' # loaction of oasis output files
# section depth map assumed to be  40  m and 100m .
output_format ='.csv'
values_for_imaging = [40.,100.]  
# we create self container of geosurface object        
geo_surface_obj = Geosurface( path =path_to_oasisfiles, 
                             depth_values = values_for_imaging, 
                             )
geo_surface_obj.write_file(fileformat = output_format, 
                            savepath =savepath )

parser_file ='nbleDH.csv'

#savepath : path to save outfile borehole files 
savepath = None 

# if set to False , user will add step by step all data with the layer thicknesses 
build_borehole_auto=True
# create a borehole object 
borehole_obj = Drill (well_filename= 'data/drill_example_files/nbleDH.csv', 
                   auto= build_borehole_auto)
# data2write : which kind of data do you want to output ?
# borehole geology? borehole geochemistry sample ? or borehole survey elevation ? 
kind_of_data2output = '*'   # can be 'collar' `geology` , `sample`
borehole_obj.writeDHData(data2write=kind_of_data2output, 
                         savepath  = savepath )


# set the stn profile file or EDI-file
# set path to your profile file 

file_stn='K6.stn'           # name of zonge station profile file 
# uncomment `path_to_stn_profile_file is zonge station file is used 
path_to_stn_profile_file = None #  os.path.join(os.environ["pyCSAMT"],'csamtpy','data', file_stn) 
            # OR 
            
# provided edipath or jpath when used edifiles or jfiles and 
edipath_or_jpath = os.path.join(os.environ["pyCSAMT"],'data','edi') # None 
#edipath_or_jpath =None                     # uncomment section if station stn file is provided 

# to see documentation of that function , set "see_documentation " to 'true'. 
see_documentation=False
# to plot profile : set PLOT to "True"
PLOT=True 
# set "*" for three profile. 

plot_type ='*'          # could be |"topography or "topo" |"speration" or "stn"| "azimuth or "az|
                        # could use only the two first letter or more for ploting or 1,2,3 or 123|*
set_stnNames =True      # set it to True if you want to see station names appear on xaxis 


# create profile_obj 
plot_1d_obj= Plot1d()

from pycsamt.viewer.plot import Plot1d 
set_stnNames =True # show staion names labels 
plot_type ='*' #or 123  or can be [Top|az|sep] or [1|2|3] for individually
path_to_stn_profile_file = 'data/avg/K1.stn'
plot1d().plot_topo_sep_azim(fn = None, # set edipath or jpath 
                               profile_fn= path_to_stn_profile_file ,
                               plot=plot_type,
                               set_station_names=set_stnNames,
                               )
plot1d_obj = Plot1d( fig_size =figsize)

from pycsamt.viewer.plot import Plot1d 
# path to profiles stn files location 
path_to_profiles = 'data/stn_profiles'
# if you want to plot some specific profiles betwen many profiles 
# specify the profile lines  LIKE ['K9.stn', 'K8.stn']
profile_lines = ['K9.stn', 'K8.stn'] #profile_lines = ['K9.stn', 'K8.stn'] 
#scaled the line scale # can be `m` or `km` . Default is `m`
scale ='km'                     
#save figure 
savefig ='test_multisites.png'
#set to False if you dont want to see stations labels 
show_station_labels = True  
Plot1d().plot_multiStations(path = path_to_profiles, 
                              profile_lines =profile_lines,
                              scale =scale, 
                              savefig =savefig, 
                              show_station_labels = show_station_labels )