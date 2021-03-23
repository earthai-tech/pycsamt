# -*- coding: utf-8 -*-
"""
   .Script to write corrected edi by applying filter. Data can be Electromagnetic 
   Array Profiling(EMAP) or Magnetotelluric data (MT). Available filters are
   `tma`, `flma`, `ss` and `dist`.'Trimming moving-average' `tma` or Adaptative
   moving average `ama` or Fixed-length dipole moving average `flma` are mainly 
   used to corrected EMAP data , however can use to force processing MT data.
   static shift removal `ss` and distorsion removal `dist` are usually use
   to process MT data.

Created on Fri Mar 19 19:15:46 2021

@author: @Daniel03
"""
import os 
from csamtpy.ff.processing.corr import shifting


# from csamtpy.ff.processing.corr import shifting 

# profile edipath : full path to edifiles or single edifile 
edipath =os.path.join(os.environ['pyCSAMT'], 'data','_outputEDIFiltered_AMA') #'edi')#, 'new_csa000.edi' )

# Applied filter 
FILTER = 'ama'                        # availables filters [`tma`, `flma`,`ama`, `ss`, `dist`]
                                        # default is `tma`
                                        
# number of points : to computed the window width  
number_of_filter_points = 7.            # default is 7. set to 01. when use si ngle edifiles 
# number of skin depth : specially provided to compute rho ith AMA filter 
number_of_skin_depth =7.                # default is 3. can be 1 to 10 skin depths 

# reference frequency 
reference_frequency = 8192.             # frequency at clean data , usefull when data is EMAP data 
                                        # not use for MT data 
                                        
#dipole length  in meter : provided to integrate on one segment of dipole 
dipoleLength = 50.                      # default is 50m for CSAMT survey 

#datatype correspond either EMAP data or MT data , if None , will detect automatically       
datatype = None                     # Type of edifile , can be `mt` or`emap` 
    
# path to hold edi outputs files 
savepath =  None

# new edi output filenames 
new_edifilename = None 
#------------------------------------------------------------------------------
# Optional params but usefull when used MT data 

reduce_res_factor_x= 1.             #  static shift factor to be applied to x
                                    # components (ie z[:, 0, :]).  This is
                                    # assumed to be in resistivity scale
reduce_res_factor_y = 1.            #static shift factor to be applied to y
                                    #components (ie z[:, 1, :]).  This is
                                    # assumed to be in resistivity scale
distortion_tensor = None            # real distortion tensor error as a 2x2
                                    #np.ndarray(2, 2, dtype=real) 
distortion_err_tensor = None
#------------------------------------------------------------------------------
# call correction object 

corr_obj= shifting().write_corrected_edi(data_fn = edipath, 
                             number_of_points =number_of_filter_points,
                             reference_frequency=reference_frequency,
                             number_of_skin_depth=number_of_skin_depth, 
                             dipole_length =dipoleLength, 
                             FILTER=FILTER, 
                             edi_newname = new_edifilename, 
                             datatype =datatype, 
                             reduce_res_factor_x=reduce_res_factor_x, 
                             reduce_res_factor_y = reduce_res_factor_y, 
                             distortion_tensor= distortion_tensor, 
                             distortion_err_tensor = distortion_err_tensor)


