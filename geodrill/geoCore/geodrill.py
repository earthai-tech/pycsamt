# -*- coding: utf-8 -*-
"""
===============================================================================
    Copyright © 2021  Kouadio K.Laurent
    
    This file is part of pyCSAMT.
    
    pyCSAMT is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    
    pyCSAMT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.
    
    You should have received a copy of the GNU Lesser General Public License
    along with pyCSAMT.  If not, see <https://www.gnu.org/licenses/>.

===============================================================================
  
.. _module-geodrill::`geodrill.geoCore.geodrill`
    
        :synopsis: Deal with geological strata , geological structural info and 
            geological dataBase. Module to deal with inversion file ,and
            geologica data and Borehole data  Can read , create an output file
            directly use by Oasis Montaj of Geosoft corporation and Goldern Software
            Can also build borehode data andgenerate a report of survey location 
            Deal with geostrata   module to incorporate conventional
            geological codes , pattern and colors. 
    
Created on Sat Sep 19 12:37:42 2020
@author: Daniel03 
"""

import os 
import copy 
import warnings
import datetime
import numpy as np 
import pandas as pd 
import scipy as sp

import  csamtpy.utils.exceptions as CSex
from csamtpy.modeling import occam2d
# from csamtpy.utils.avg_utils import ReadFile as rfi
from csamtpy.ff.core.cs import Profile

from csamtpy.utils import func_utils as func
from csamtpy.utils import plot_utils as punc
from csamtpy.utils import Agso
from geodrill.geoCore import structural as STRL
from geodrill.geoDB.sql_recorder import GeoDataBase 

try : 
    from csamtpy.utils._csamtpylog import csamtpylog
    _logger=csamtpylog.get_csamtpy_logger(__name__)
except :
    pass
try : 
      import scipy.stats as spSTAT
      scipy_version = [int(vers) for vers in sp.__version__.split('.')] 
      if scipy_version [0] == 1 : 
          if scipy_version [1] < 4 :
              warnings.warn('Note: need scipy version 1.4.0 or more . '
                            'It may probably  get a trouble when import "stats" attribute'
                            'under such version.It may probably move from ', ImportWarning)
              _logger.warning('Note: need scipy version 0.14.0 or higher'
                              ' or for stats.linearegress. Under such version'
                              'it might not work.')
      # from sp import stats 
      stats_import =True 
          
except :
    warnings.warn('Could not find scipy.stats, cannot use method linearregression'
                  'check installation you can get scipy from scipy.org.')
    _logger.warning('Could not find scipy.stats, cannot use method linearegress'
                    'check installation you can get scipy from scipy.org.')
    
    stats_import =False 


class Geodrill (object): 
    """
    Class to manage data from Occam2D and Model so to create section of each  
    sites and it depth.
    
    Each station constitues an attribute framed by two closest  
    point of station offsets from model resistivities.
    Deal with  True resistivities get on the survey or with others firms. 
    In fact , input truth resistivities values into our model , produce an accuracy 
    underground map.The challenge to build pseudolog allow to know how layers are 
    disposal in underground so to emphasize the large conductive zone  especially 
    in the case of groundwater exploration.  Program works in combinaison with
    geophysic data especially Occam 2D inversion data,  and  geological data.
    Actually the program deal only with  Occam 2D inversion files or  Bo Yang
    (x,y,z) files. We intend to extend later with other external softares but
    can generate output directly see use with Golder sofware('surfer'). If user
     has a golder software installed on its computer , It can use  output files 
     generated here  to produce 2D map so to compare both maps to see how far 
    is the difference between  Model map and detail-sequences map ) 
    "pseudosequences model" could match better the reality of underground).
    Details sequences map is most closest to the reality  When {step descent}
    parameter is not too  small at all. Indeed True geological data allow
    to harmonize the value of resistivity produced by Occam2D model so to force 
    the pogramm to make a correlation  between data fromtruth layers and the model
    values.


    Arguments
    ----------
        **model_fn** : str,  
                    full path to Occam model  file .                             
        **iter_fn** :  str,                
                    full path to occam iteration  file
                                                        
        **data_fn** :  str,  
                    full path to occam_data file 
        **input_resistivities** :   array_like, 
                            Truth values of resistivities 
        **step_descent**: float,   
                        step to enforce the model resistivities to
                         keep the  truth layers values as reference
                         Step to cut out data and to  force resistivites 
                         calcualted to match the reference data as input 
                         resistivities if not provided the step will be
                         20% of D0I
        **input_layers** :  array_like,  
                        True input_layers names : geological 
                        informations of encountered layers 

    =============  ===========  ===============================================
    Attributes     Type         Explanation 
    =============  ===========  ===============================================
    doi            str          depth of investigation might
                                be float or str like "1km" =1000
    depth_scale    str          scale of imaging depth can be 
                                `km` or `m`.  Default is`m` 
    lc_AD_curves   dict         customize line color of average curve 
                                and details sequences 
    elevation      array_like   elevation of survey area 
    iter2dat_fn    str          full path to Bo Yang (x, y, z) model_file     
    bln_file       str          full path to station location file  additional 
                                file issue from Bo Yang (x, y, z file)
    =============  ===========  ===============================================
    
    .. note:: In this module, all resistivites are in ohm.meter not in log10.
    
    :Example: 
        
        >>> from geodrill.geoCore.geodrill import Geodrill 
        >>> path =os.path.join(os.environ ['pyCSAMT'], 'csamtpy', 
        ...                       'data', 'occam2D')
        >>> geo_obj = Geodrill( input_resistivities=[300, 500, 
        ...                                             1000, 2000, 4000, 6000],
        ...                       input_layers =['alluvium', 
        ...                                      'amphibolite','altered rock',
        ...                                                    'augen gneiss', 'granite'],
        ...                       mesh_fn=os.path.join(path, 'Occam2DMesh')
        ...                       iter_fn = os.path.join(path, 'ITER17.iter'), 
        ...                       model_fn =os.path.join(path, 'Occam2DModel') , 
        ...                       data_fn =os.path.join(path, 'OccamDataFile.dat'),
        ...                       doi='1km', 
        ...                       step_descent=200., 
        ...                          )
        >>> geo_obj.geo_build_strata_logs()
        >>> geo_obj.geo_name_S01
        >>> geo_obj.geo_name_S00
        >>> geo_obj.geo_name_S46
        >>> geo_obj.geo_d
        >>> geo_obj.geo_drr
        >>> geo_obj.geo_daver
        >>> geo_obj.geo_depth 
        >>> geo_obj.geo_dstep_descent
    """
    
    geo_rocks_properties={
                    "basement rocks" :          [1e99,1e6 ],
                    "igneous rocks":            [1e6, 1e3], 
                    "duricrust"   :             [5.1e3 , 5.1e2],
                    "gravel/sand" :             [1e4  , 7.943e0],
                    "conglomerate"    :         [1e4  , 8.913e1],
                    "dolomite/limestone" :      [1e5 ,  1e3],
                   "permafrost"  :              [1e5  , 4.169e2],
                    "metamorphic rocks" :       [5.1e2 , 1e1],
                    "tills"  :                   [8.1e2 , 8.512e1],
                    "standstone conglomerate" :  [1e4 , 8.318e1],
                    "lignite/coal":               [7.762e2 , 1e1],
                    "shale"   :                   [5.012e1 , 3.20e1],
                    "clay"   :                    [1e2 ,  5.012e1],
                    "saprolite" :                 [6.310e2 , 3.020e1],
                    "sedimentary rocks":          [1e4 , 1e0],
                    "fresh water"  :              [3.1e2 ,1e0],
                    "salt water"   :              [1e0 , 1.41e0],
                    "massive sulphide" :           [1e0   ,  1e-2],
                    "sea water"     :             [1.231e-1 ,1e-1],
                    "ore minerals"  :             [1e0   , 1e-4],
                    "graphite"    :               [3.1623e-2, 3.162e-3]
        
        }
    
    geo_params =['model_x_nodes', 'model_z_nodes', 'model_res', 'station_names', 
                    'station_location', 'elevation']
    geoparams =['geo_offsets', 'geo_depth', 'geo_res', 'geo_name', 
                    'geo_location', 'geo_elevation']
    
    def __init__(self, iter2dat_fn=None , model_fn =None , data_fn=None ,
                 iter_fn=None , mesh_fn=None , bln_fn =None , **kwargs):
                
        self._logging = csamtpylog.get_csamtpy_logger(self.__class__.__name__)
        
        self.iter2dat_fn = iter2dat_fn
        self.data_fn =data_fn 
        self.iter_fn =iter_fn 
        self.mesh_fn =mesh_fn 
        self.model_fn =model_fn 
        self.bln_fn =bln_fn  
        self.input_layers =kwargs.pop('input_layers', None)
        self.input_resistivities =kwargs.pop('input_resistivities', None)
        self.doi =kwargs.pop('doi', '1km')
        
        self.elevation =kwargs.pop('elevation', None)
        self.step_descent =kwargs.pop('step_descent', None)
    
        
        
        for key in self.geo_params :  # Initialise attributes to None value to recover it later  
            setattr(self, key, None)
        
        for key in list(kwargs.keys()):
            setattr(self, key, kwargs[key])
            
        if self.model_fn is not None and self.data_fn is not None or self.iter2dat_fn is not None : 
            self.set_geodata()
    
    def set_geodata(self,iter2dat_fn =None ,  model_fn=None , data_fn=None , 
                    iter_fn=None , mesh_fn=None, **kwargs):
        """
        Readmodel data and collected from each site its value from surface to depth 
        
        1. Read with Occam 2D outputs files
        
        :Example: 
            
            >>> form csamtpy.geodrill import Geodrill
            >>> path_occam2d = os.path.join(os.environ ['pyCSAMT'],
            ...                                'csamtpy', 'data', 'occam2D')
            >>> path_i2d =os.path.join(os.environ ['pyCSAMT'], 
            ...                           'csamtpy', 'data', '_iter2dat_2')
            >>> geo_obj =Geodrill(model_fn  = os.path.join(path_occam2d,'Occam2DModel'), 
            ...                  mesh_fn = os.path.join(path_occam2d, 'Occam2DMesh'), 
            ...                  iter_fn = os.path.join(path_occam2d, 'ITER17.iter'), 
            ...                  data_fn = os.path.join(path_occam2d, 'OccamDataFile.dat'))
            
        2. Read with only iter2dat file and bln file 
        
        :Example:
            
            >>> form csamtpy.geodrill import Geodrill
            >>> geo_obj =Geodrill(iter2dat_fn  = os.path.join(path_i2d, 'K1.iter.dat'), 
            ...                  bln_fn = os.path.join(path_i2d, 'K1.bln'),)
            >>> geo_model_off = geo_obj.model_x_nodes
            >>> geo_model_res= geo_obj.model_res
            >>> geoS01 = geo_obj.geo_name_S01
            >>> geoS00= geo_obj.geo_name_S00
            >>> geoS46=geo_obj.geo_name_S46
        """
        
        print('{0:-^77}'.format('GeoDrill *Data* info'))
        
        def frame_each_site_into_3_offsets(site_name, model_res, 
                                           model_offset, site_offset): 
            """
            Fonction to frame each site into 3 offsets and put array 

            :param site_names:  name of site , eg S01
            :type site_names: str  
               
            :param model_res: ndarray(depth, station_offset) resistivity 
            :type model_res: ndarray(depth, station_offset)

            :param model_offset: the offset generated by mesh 
            :type model_offset: array_like  
                 
            :site _offset:  offset of that site  : 15.0 m 
            :site _offset: float 
                                 
            """
            # fist check if site offset in in model offset
            model_offset=np.array([np.around(ss, 2) for ss in model_offset])
            site_offset =round(site_offset,2)
 
            if site_offset in model_offset : 
         
                # get new index 
                new_index = int(np.where(model_offset==site_offset)[0])
                new_off = model_offset[new_index]
        
            elif site_offset not in model_offset :
     
                if site_offset > model_offset.max (): # the case where offset is much larger than the mesh (very rare case )
                    new_index = int( np.where(model_offset==model_offset.max())[0])
                    new_off = model_offset.max ()
                else : 
                    for ii, di in enumerate(model_offset): 
                        if di > site_offset : # if value of model is greater than check the closest site offset 
                            dx0 = round(abs(site_offset -di ),2) # compute the distance betewn di and offset (dxmax )
                        # try co compute the previous and the next point distance 
                        # use try because sometimes , previous point is none if model value is at index 0 
                        # if model value at index 0 , then dx0 = dmin 
                            try : dxmin =abs(model_offset[ii-1]-site_offset)
                            except : dxmin =round(dx0 , 2)
                            if dx0 <=  dx0 : 
                                new_off , new_index = di, ii  
                                break 
                            elif dxmin < dx0 : 
                                 new_off, new_index = model_offset[ii-1], ii-1 # take the preious 
                                 break 
     
            # now get the list of resistivity between this index 
          
            if new_off == model_offset[0] : # the first site index then take ii+1, ii+2 
                dm0 = np.concatenate ((model_res[:,new_index].reshape(model_res[:, new_index].shape[0], 1), 
                                       model_res[:,new_index+1].reshape(model_res[:, new_index+1].shape[0], 1), 
                                       model_res[:,new_index +2 ].reshape(model_res[:, new_index+2].shape[0], 1)), 
                                       
                                       axis =1)
            elif new_off == model_offset[-1] : 
                dm0 = np.concatenate ((model_res[:,new_index-2].reshape(model_res[:, new_index-2].shape[0], 1), 
                                       model_res[:,new_index-1].reshape(model_res[:, new_index-1].shape[0], 1), 
                                       model_res[:,new_index ].reshape(model_res[:, new_index].shape[0], 1)),
                                       axis =1)
            
            else : 
                dm0 = np.concatenate ((model_res[:,new_index-1].reshape(model_res[:, new_index-1].shape[0], 1), 
                                       model_res[:,new_index ].reshape(model_res[:, new_index].shape[0], 1), 
                                       model_res[:,new_index +1].reshape(model_res[:, new_index+1].shape[0], 1)),
                                       axis =1)
                
            return site_name, dm0 

 
        # ----> attributes statements 
        if model_fn is not None : self.model_fn =model_fn 
        if data_fn is not None : self.data_fn =data_fn 
        if iter_fn is not None : self.iter_fn =iter_fn 
        if mesh_fn is not None : self.mesh_fn =mesh_fn 
        
        if iter2dat_fn is not None : self.iter2dat_fn = iter2dat_fn 
        
        # then assert all input files 
        if self.iter2dat_fn is not None : 
            self._logging.info ('Read Geodata from  Bo Yang Iter2Dat files ')
            if self.bln_fn is None : 
                mess =''.join(['No (*bln) station file found. Need sites names and sites locations.',
                               ' Use mode "Iter2Dat"  :<from csamtpy.modeling.occam2d import Iter2Dat> : to write a *bln file ', 
                               ' use a Profile module : < from csamtpy.ff.core.cs import Profile > .', 
                               ' You can also use occam2D output files like model, mesh, data, and iteration files.'])
                warnings.warn('! Error reading *.bln flile !' + mess)
                self._logging.error(mess)
                raise CSex.pyCSAMTError_geodrill_inputarguments(
                    '! Error reading *.bln flile ! No (*bln) station file found.'
                     ' Need sites names and sites locations.')
            i2d_obj = occam2d.Iter2Dat(iter2dat_fn=self.iter2dat_fn,
                                       bln_fn = self.bln_fn)

            self.model_x_nodes= i2d_obj.model_x_nodes
            self.model_z_nodes = i2d_obj.model_z_nodes
            self.model_res = i2d_obj.model_res
            self.station_names = i2d_obj.station_names
            self.station_location = i2d_obj.station_location
   
            
        else : 
            self._logging.info ('Read Geodata from  Occam2D files  ')
            # check if all attributes are not None 
            for attr , ins  in zip(['model_fn', 'iter_fn', 'data_fn', 'mesh_fn'],
                                   [self.model_fn, self.iter_fn, self.data_fn, self.mesh_fn]): 
                if ins ==None : 
                    msg =' !No {0} file found ! Please provide the Occam 2D {0} file'.\
                        format(str(attr).replace('_fn', ''))
                    self._logging.error(msg)
                    warnings.warn(msg)
                    raise CSex.pyCSAMTError_plot_geoinputargument(msg)
                
            # Recreate object and get ncessaries attributes 
            model_obj= occam2d.Model(model_fn =self.model_fn, iter_fn =self.iter_fn, 
                                    mesh_fn =self.mesh_fn )
            data_obj =occam2d.Data(data_fn= self.data_fn)
            
            # then get important attributes and re
            self.model_x_nodes =model_obj.model_station_offsets
            self.model_z_nodes =model_obj.model_depth_offsets
            self.model_res = model_obj.model_resistivity 
            # get from data sites names and station location and elevation 
            self.station_names =data_obj.data_sites 
            self.station_location = data_obj.data_offsets
            
        # build for each site an array of three values framed into the main sites 
        #    get the the offset of the site 
        # build a dictionnary of geoname
        
        #------Mange DOI (investigation depth) --------
        
        # check the depth of investigation and ge the new_nodel matrix 
        
        if self.doi is None : self.doi = 1000.
        elif self.doi is not None :
            self.doi =punc.depth_of_investigation(self.doi) # return meter value

        # for consistency let check how the depth is ranged : Depth is  minimum to max depth , if not flip data 
        if self.model_z_nodes [0] > self.model_z_nodes[-1]: 
            self.model_z_nodes= self.model_z_nodes[::-1]
        
        # chech the doi data and compare to max depth data 
        
        if self.doi > self.model_z_nodes.max (): 
            mess='Maximum Input doi value = {0} m is larger'\
                  ' than model depth = {1} m. Investigation depth "doi" will be resetting at = "{1}"m.'.\
                              format(self.doi, self.model_z_nodes.max())
            
            self.doi = self.model_z_nodes.max()         # resseting doi to maximum value in model depth 
            
            warnings.warn(mess)
            self._logging.debug (mess)
        
        if self.doi in self.model_z_nodes: 
            dep_index = int(np.where(self.model_z_nodes== self.doi)[0])
            
            if self.doi == self.model_z_nodes[-1]:  # if doi is the  max depth , 
                self.__setattr__('geo_depth', self.model_z_nodes)
                
            elif self.doi  <     self.model_z_nodes[-1]:         # resize model _resistivity 
                self.model_res = np.resize(self.model_res, (dep_index +1 , self.model_res.shape[1]))
                self.__setattr__('geo_depth', self.model_z_nodes[:dep_index+1]) # resize the large depth to doi max  
            
            print('---> resetting model doi to = {0} m depth !'.format(self.doi))
        else :
            for index, dep  in enumerate(self.model_z_nodes) : 
                if dep  > self.doi  :       # seek the depth index that match better the  input doi 1014>1000 :index=23,
                    dep_index = index-1     # get the previous index dep_index = 22
                    self.model_res = np.resize(self.model_res, (index, self.model_res.shape[1])) # index because of python count (add+1) to depth index
                    self.doi = self.model_z_nodes[dep_index]
                    print('---> resetting model doi to ={0} m depth !'.format(self.doi))
                    
                    self.__setattr__('geo_depth', self.model_z_nodes[:index])
                    break 
                
        # build geo_station_names and geodict rho 
        self.geo_d ={}
        if isinstance(self.station_names , (tuple, np.ndarray)): 
                self.station_names=self.station_names.tolist()
    
        for name in self.station_names : 
    
            site_offset_value = self.station_location[self.station_names.index(name)]
            
            geo_name , geodat = frame_each_site_into_3_offsets(site_name= name, 
                                                               model_res = self.model_res, 
                                                               model_offset=self.model_x_nodes, 
                                                               site_offset=site_offset_value)
            
            self.__setattr__('geo_name_{0}'.format(geo_name), geodat)
            # convert resistivities value from log10 rho to ohm meter 
            self.geo_d[name]= np.power(10, getattr(self, 'geo_name_{0}'.format(geo_name))) # 
            
            # self.geo_dict[name]=  getattr(self, 'geo_name_{0}'.format(geo_name))
            
  
    def geo_replace_rho (self, input_resistivity_range =None, 
                         input_layers=None,  **kwargs ): 
        
        """
       Allow ro replace the calculated resistivities from model
       to real resistivites obtained on survey area with other companies. 
       The accuracy depend the significant values of input  resistivites.
       More input resitivities are,more accuracy in the design of underground model
       geostratigrapyhy model should be.

       :param input_resistivity_range: an array of resistivity on the site
       :type input_resistivity_range: array_like 
                      
       :param input_layers: layers_names , eg .`granite`, `fault`, `river `
       :param input_layers: list or arrays 
                            
       :param depth_range: array _of depth of specific layer 
       :type depth_range: array _like
                    
   
        Example of table  of geoinformation collected on the field somewhere 
        
        ====================  ====================  ===========================
        Structure	          Rho mean value (Ω.m)	Rho range  (Ω.m)
        ====================  ====================  ===========================
        Granite	                2000	               1000-3000
        Fault fracture zone     120	                   60-180
        River water 	        68	                   66-70
        ====================  ====================  ===========================
        
        """
  
        if self.iter2dat_fn is None  or self.model_fn is None : 
            self.set_geodata() # loead attribute 
        
        def find_and_replace_rho(rowlines, value_range): 
            """
            Small function replacement function, the model's resistivities are 
            replaced by the true value of the survey. T higher the replacement is,
            the more prominent the specificities of the layers become.
 
            Parameters
            ------------
                * rowlines :  array_like ,
                        array of resistivities values 
                * value_range : array_like 
                            array of resistivities from survey area.
            Returns 
            ---------
                array_like 
                      rowlines , new data range                    
            """
            # for consistency convert to float value 
            try : 
                value_range= np.array([float(ss) for ss in value_range])
            except : 
                raise CSex.pyCSAMTError_geodrill_inputarguments(
                    'Value provided as resistivity range must be an array of float number.')
                                                                
            tem=[]
            for mm in value_range :     # loop the value range 
                op = rowlines-mm        # substract rowlines from value range to seek the minimum close
                op=np.array([abs(aj) for aj in op]) # keep absolute value for difference 
                tem.append(op)                      # keep it on temporray list  
            ts = func.concat_array_from_list(tem)   # concat list to axis =0 order 
            
            for ii in range(len(rowlines)):         # loop the rowlines line now 
                us=ts[:,ii].argmin()                # keep minimum index 
                rowlines[ii]=value_range[us]        # change the resistivity values 
                
            return rowlines
        
        
        def ascertain_input_layers(input_layers, input_rho): 
            """
            Function to assert the length of input layers and the input resistivities 
            to avoid miscomputation .
            
            Parameters
            -----------
                * input_layers : array_like |list, 
                            list of input
                * input_rho : list ,
                    layer list of resistivities provided
            Returns 
            -------
                list 
                     new list of input layers 
            """
            ilay =[str(ly) for ly in input_layers]      # for consistency put on string if user provide a digit
            if len(input_rho) ==len(ilay): 
                return ilay
            
            elif len(input_rho) > len(ilay): 
                sec_res = input_rho[len(ilay):]         # get the last value of resistivities  to find the structres names ans structures sresistivities 
                geos =Geodrill.get_structure(resistivities_range=sec_res) # get the name of structure as possible 
                if len(geos)>1 : tm = 's'
                else :tm =''
                print('---> !We added other {0} geological struture{1}. You may ignore it.'.format(len(geos), tm))
                ilay.extend(geos)                       # then , extend the list 
                return ilay 
            elif len(ilay) > len(input_rho): 
                ilay = ilay[:len(input_rho)]            # truncated the list 
                return ilay 
        
        
        # if range of resistivity is provided then use it for average 
        if input_resistivity_range is not None : 
            self.input_resistivities = input_resistivity_range 
        if self.input_resistivities is None and self.input_layers is not None : 
            warnings.warn(' !Without any input resistivities values, we can not set only your layer names.'\
                          ' Bring more details about your layers by adding their corresponding resistivities values.')
            # print('-->!')
            
            self.input_layers=None       #abort the input layers by renitializing to NoneType values
            
         
        auto_mess=''            # keep automatic message 
        if self.input_resistivities is None :
            if self.input_layers is None : 
                n_layers =7.  # number of slices layer juged by the program 7 
            else : n_layers =len(self.input_layers)
            
            mess="".join(["Resistivity range is not provided . Sites Depth will",
                          " be cout out into {0} slices as possible layers ", 
                          "  below the site. If the number of slices doesnt",
                          " suit the purpose , please change the number ",
                          "  of slices using argument <input_layers> ",
                          "to provided the real layer's names."])
            
            warnings.warn(mess.format(int( n_layers))) 
            self._logging.info(mess.format(int( n_layers)))
           
          
            # get the model resistivities minimum and maximum from selected doi, it is much  better 
            # than  to select min max into the  global resistivities models.  
            
            minmax=[(res_values.min(), res_values.max()) for stn, res_values  in self.geo_d.items()] #
            maxres= max([mm[1] for mm in minmax])
            minres= min([mm[0] for mm in minmax])
            
            self.input_resistivities =np.linspace( minres, maxres, int( n_layers))
            self.input_resistivities  = np.around(self.input_resistivities,2) # rounded to 2 , because resistiviris are in ohm-m.
            auto_mess ='Automatic'
 
            
        if isinstance (self.input_resistivities  , (tuple, list)): 
            self.input_resistivities  = np.array(self.input_resistivities )

        elif isinstance(self.input_resistivities , (float, int, str)): 
            try : 
                self.input_resistivities = float(self.input_resistivities )
            except : 
                raise CSex.pyCSAMTError_plot_geoinputargument(
                    'Can not converted value <%s> into float number.'\
                    ' Please provide a float number.'%self.input_resistivities )
                
        # Display infos 
        print('**{0:<37} {1} {2}'.format('{0} Layers sliced'.format(auto_mess ),
                                         '=' , len(self.input_resistivities  )))
        print('**{0:<37} {1} {2} (Ω.m)'.format('{0} Rho range'.format(auto_mess ),
                                               '=' , tuple(self.input_resistivities  ) ))
        print('**{0:<37} {1} {2} {3}'.format(' Minimum rho ',
                                             '=' ,self.input_resistivities.min(), 'Ω.m' ))
        print('**{0:<37} {1} {2} {3}'.format(' Minimum rho ',
                                             '=' , self.input_resistivities.max(), 'Ω.m' ))
       
        # so to get the structures for each input resistivities 
                 # ascertain unput layers first if no misleading value is inputted 
        if self.input_layers is not None : 
            # rebuild new list of layer by adding necessary strutures or substracting unecessaries structures
            formations= ascertain_input_layers(input_layers= self.input_layers,
                                                  input_rho= self.input_resistivities)
        #     self.depth_range =depth_range 
        else : formations = Geodrill.get_structure(self.input_resistivities)
        
        self.input_layers=formations 
        

        # self.__setattr__('geo_dict_rho', nOne)
        
        self.geo_drr=copy.deepcopy(self.geo_d)              # use deep copy to make difference between two dict data set 
        
        
        tem =[] # dictionnary opf  replacement rho
        for key , geovalue in self.geo_drr.items() : 
            for ii, rowlines in enumerate(geovalue): 
               out_aver =find_and_replace_rho(rowlines= rowlines,
                                     value_range= self.input_resistivities)
               tem.append(out_aver)
            self.geo_drr[key]=func.concat_array_from_list(list_of_array=tem)
            tem=[]
            

            
    @staticmethod 
    def get_structure (resistivities_range):
        """
        function to get according the range of resistivities values , 
        the corresponding associated geological rocks 
        The list of electrical properties of rocks is not exhaustive ,
        can be fill by others 

        :param resistivities_range:   array of input_resistivities 
        :type resistivities_range: array_like,                      
        
        :returns: the list of geological structures
                    form go_electrical _rocks properties 
        :rtype: list 
        
        """
        if isinstance(resistivities_range, (float, str, int)): # only single value provided than put on list.
            try : 
                resistivities_range=[float(resistivities_range)]
            except : 
                raise CSex.pyCSAMTError_geodrill_inputarguments(
                    'Can not convert <%s> to float number !Input resistivity '
                     'must be a float number not <%s>!'% (resistivities_range,
                                                          type(resistivities_range)))
        
        if isinstance(resistivities_range, (tuple, list, np.ndarray)):  # for consistency , chek again input values 
            try :  resistivities_range =np.array([float(ss) for ss in resistivities_range])
            except : 
                raise CSex.pyCSAMTError_geodrill_inputarguments(
                    'Input argument provided is wrong.Please check your resistivities '
                    ' range, values must be a float number.')
             
        geo_structures=[]                   # in fact append the idde of resistivities values located 
                                            # so to be sure that the resistivities will match exactly the layer found in geo_electricl_property of rocks 
        f=False 
        for resr in resistivities_range :
            f=False
            for rocks, eprops in Geodrill.geo_rocks_properties.items(): 
                if eprops[-1]<resr < eprops[0] : 
                    f=True
                    geo_structures.append(rocks.capitalize())
                    
                    break 
            if f==False : 
                geo_structures.append('*! Struture not found')
                
        return geo_structures
 
    @staticmethod 
    def get_average_rho(data_array, transpose =False ): 
        """
        Function to average rho to one point to onother . 
        It show the lowest point and the maximum point averaged . Function averaged rho 
        value between local maximum and local minima values . 
        if data  values of station are located on columnlines , set transpose to True 
        then rotate the matrix to find minima and maxima  locals value then 
        calculated averaged rho after will return matrix transpose
         as the same shape as inputted . .Defaut is **False** . 
     
        :param data_array: data of resistivities collected at the site point  
        :type data_array: ndarray  
        
        """
        if data_array.dtype not in ['float', 'int']:  # be sure data type is on float value
            try :
                data_array= data_array.astype('float64')
            except : 
                warnings.warn("It seems somethings wrong"
                              " happened during data conversion to float values.")
                raise CSex.pyCSAMTError_inputarguments('Could not convert value'
                                                       ' to float numbers , Please check your data!')
            
      
        if transpose is True :          # transpose the data to loop the rownlines in first times 
            data_array =data_array.T
        
        exem =[punc.average_rho_with_locals_minmax(rowlines)
               for rowlines in data_array ]         # build list of array transposed  and concat vaue 
        new_data = func.concat_array_from_list(exem)
        if transpose is True :          # transpose the data to keep finally the same shape  
            new_data = new_data.T
            
        return new_data
        
  
    def geo_build_strata_logs (self, input_resistivities=None, input_layers =None, 
                               step_descent = None,  **kwargs):
        
        """"
        Read resistivits data  got on survey area and build geological strata block. 
        If constrained_electrical_properties_of_rocks is False , will build a log 
        by considering the  conventional electrical property of rocks 
        can be find through this link.
            
        Parameters
        ------------
            * input_resistivity :array_like
                            an array of resistivity on the site
            * step_descent : float 
                            depth value to averaged rho . Must be smaller as possible 
                            if None , it take the 2% times the investigation depth  
            * input_layers : list or array 
                        layers_names , eg `granite`, `fault`, `river` 
        
        Returns 
        ---------
            obj , 
                Build stratigraphy log obj 
            
            
        A Sample of electrical_properties_of_rocks is below
            
        =========================  ===================  =======================
        Rocks                           Max Rho             Min Rho (ohm-m)
        =========================  ===================  =======================
        igneous rocks                   10^6                10^3 
        duricrust                       5.10^3              5.10^2
        gravel and sand                 10^4                10^2.90(800)
        conglomerate                    10^4                10^1.95(90)
        dolomite/limestone              10^5                10^3
        permafrost                      10^5                10^2.62(750)
        metamorphic rocks               5.10^2              10.^1
        tills                           8.10^2              10^1.93(85)
        standstone conglomerate         10^4                10^1.92 (80)
        lignite/coal                    10^2.89(790)        10^1
        shales                          10^1.7(50)          10^1.48(30)
        clays                           10^2                10^1.7(50)
        saprolite                       10^2.08(120)        10^1.48(30)
        sedimentary rocks               10^4                10^0
        fresh water                     3.10^2              10^0
        salt water                      10^0                10^-0.15
        massive sulfide                 10^0                10^-2
        sea water                       10^-0.09(0.8)       10^-1
        ore minerals                    10^0                10^-4
        Graphite                        10^-2.5             10^-3.5
        =========================  ===================  =======================
        
        .. seealso:: https://www.eoas.ubc.ca/ubcgif/iag/foundations/properties/resistivity.htm
                list is not exhaustive and depend of the geological formations of 
                survey area. 
                
        .. note:: list is not Exhaustive, use the data base script to populate 
                    most of goeological electrical properties.
        """

        def get_conductive_zone (dep_array, rho_array, step_in_deep):
            
            """
            Get a conductive zone is important for many purposes .In the case of 
            groundwater exploration for instance, sometimes in deeper, because of 
            heterogeneities of structures in underground , it is much difficult 
            to take some minima rho as a conductive area or conductiv productive 
            veification drill point. To be sure that this point is really among
            a good saturated zone, it is better to get   averaged rho in some 
            distance and to see how the layer resistivity value goes on on this range,
            whether it's an effective conductive zone or overlapping zone. In addition , 
            fixing resistivities values at specific distance depth allow us to detect 
            the unsatured zone as weel as the satured zone at the set of 
            investigation depth imaged. This technique allow us to build a pseudo
            specific strata that could match the zone .
          
            Parameters
            ----------
                * dep_array : array_like 
                             the imaged depth 
                * step_in_deep : float, 
                                value to averaged resistivities in deeper
                *  rho_array :  array_like ,
                                the resistivities at the sites depth {Top to bottom}
            Returns 
            --------
                array_like 
                    rho average for each station dep_averaged for each station
                                 
            
            .. note:: Much larger is the step , the accuracy becomes weak , 
                        it must be as possible 1/20 of investigation depth (doi)
                        
            :Example: 
                
                >>> import numpy as np
                >>> pseudo_depth = np.arange(0, 1201.,100)
                >>> print(pseudo_depth)
                ... dep_aver= dep(dep_array=pseudo_depth, value=100)
                ... print(dep_aver)
                ... print(pseudo_depth.shape)
                ... print(dep_aver.shape)
            """
   
            step_in_deep  =punc.depth_of_investigation(step_in_deep)
            
            v,r , dm, rm=[[] for i in range(4)]
            
            if isinstance(step_in_deep, (str, int)): 
                try : 
                    step_in_deep =float(step_in_deep)
                except : 
                    raise CSex.pyCSAMTError_plot_geoinputargument(
                        'Could not convert depth value ={} to float.'
                        ' Please check your value.'.format(step_in_deep))

            if step_in_deep< dep_array.min(): 
                raise CSex.pyCSAMTError_plot_geoinputargument(
                    'Value provided ={0} m is less than the minimum'
                    ' depth ={1} m.'.format(step_in_deep, dep_array.min()))
            
            if step_in_deep > dep_array.max(): 
                raise CSex.pyCSAMTError_plot_geoinputargument(
                    'Value provided is = {0} m is greater than '
                    'maximum depth ={1}m.'.format(step_in_deep, dep_array.max()))
                
            _init_depth =step_in_deep
            
            for index , depth in enumerate(dep_array):
                if depth <= step_in_deep :      # value less than step descent must be averaged 
                    v.append(depth)             # keep resistivities values onto list 
                    r.append(rho_array[index])
        
                if depth > step_in_deep :       # if the next value is greater than the previous one 
                    if v !=[]:                  # ccheck if kist is not enmpty 
                    
                        dm.append(np.repeat(np.array(v).mean(), len(v)))  # rebuild resistivities values with rho averaged 
                        rm.append(np.repeat(np.array(r).mean(), len(r)))
                        step_in_deep += _init_depth              #increment the next descent to step of descent 
                        v=[depth]       # initialise new list by adding the index value greater one 
                        r=[rho_array[index]]
                      
                if depth ==dep_array[-1]:
                    if len(v)==1 :                  # it length last value ==1 , means is the last value of depth
                        dm.append(dep_array[index])
                        rm.append(rho_array[index])
                    elif len(v) !=1 :               # averaged the reamin rho values  
                        dm.append(np.repeat(np.array(v).mean(), len(v)))
                        rm.append(np.repeat(np.array(r).mean(), len(r)))
              
        
            return np.hstack(tuple(rm)),np.hstack(tuple(dm))
            

        self.__setattr__('qc', .0)          # set quality control in the data 
        
        #------  STATEMENT OF INPUT ARGUMENTS------------------ 
        # check data files if provided (Dont need to check all )
        if input_resistivities is not None : self.input_resistivities= input_resistivities 
        
        if self.model_fn is None or self.iter2dat_fn is  None :  # read to populate attributes 
            self.geo_replace_rho()
            
        if input_layers is not None :
            # rebuild input_layers 
            self.input_layers = ascertain_layers_with_its_resistivities(real_layer_names= input_layers, 
                                                                        real_layer_resistivities=self.input_resistivities) 
            
          

        if step_descent is not None : self.step_descent = step_descent
        if self.step_descent is None : self.step_descent = .2 * self.doi 
        

        # ---------------get the geodictionnary from from geodata average ---------------------
        self.geo_daver ={}
        for stn, vrho in self.geo_d.items(): 
            # transpose data so to read rowlines S01,S02 etc 

            self.geo_daver[stn] =  Geodrill.get_average_rho( data_array= vrho,
                                                            transpose =True)
            
        self.qc += .25 
        #----- set attribute of geo_dstep_descent ---------------------------------------------------
        self.geo_dstep_descent={}
        for stn , geo_sd_values in self.geo_d.items(): 
            # for ii in range(3) :
                
            mv= [get_conductive_zone(dep_array= self.geo_depth,
                                     rho_array=geo_sd_values[:,ii], 
                                     step_in_deep= self.step_descent)[0] for ii in range(3)] # concat three array 
                
            self.geo_dstep_descent[stn]=func.concat_array_from_list(mv, concat_axis=1)
                
        #------------------build a dict of pseudosequences layer and resistivities ----------------------
        
        # build at dictionnary of strata from resistivities in the model and themoel depth
        self._logging.info ('Build the pseudosequences of strata.')
        
        # if constrained_electrical_properties_of_rocks is False : 
            # build layer according to geo_drr (geo_replacedrho)
            
        self.geo_dpseudo_sequence, self.geo_dpseudo_sequence_rho =[{} for ii in range(2)]
        sc=[]
        
        self.__setattr__('geo_secure_pseudo_sequence', None)
         
        for stn , georr in self.geo_drr.items(): 
            # for jj in range(3): 
            if stn =='S00' :            # first station id framed into three started at the left  
                svm = punc.build_resistivity_barplot(depth_values=self.geo_depth,
                                           res_values=georr[:,0])
            elif stn == self.station_names[-1] : # lst station id names start is lopcated at the right 
                svm = punc.build_resistivity_barplot(depth_values=self.geo_depth,
                                           res_values=georr[:,2])
            else :          # station id located at the midle framed into three stations 
                svm = punc.build_resistivity_barplot(depth_values=self.geo_depth,
                                           res_values=georr[:,1])
                

            sc.append(svm[-1])

            self.geo_dpseudo_sequence[stn]= svm[0]      # the thickness of layers 
            self.geo_dpseudo_sequence_rho[stn]=svm[1]   # the resistivities of layers 
            
        self.geo_secure_pseudo_sequence = np.array(sc)
        
        #---------------- Quality control --------------------------------.--------------------
        mess =' !Your data passes safety the Quality Control of Geodrill'
        
        if len(self.geo_secure_pseudo_sequence) == len(self.station_location):
            self.qc +=.25 
        if np.all(self.geo_secure_pseudo_sequence== self.doi ): #chek whether total thickness 
                                                                #cover the total depth
            self.qc += .25 
        else :
            warnings.warn('Data provided are inconsistencies ,'\
                          ' Please you may try to undestand carefuuly the code.')
            self._logging.debug('Data provided are inconsistencies ,'\
                          ' Please you may try to undestand carefuuly the code.')
            
        if self.qc ==1. : 
            print('---> {}.'.format(mess))
            
        print('**{0:<37} {1} {2} {3}'.format(' QC flux rate','=' , 100. * self.qc, '%' ))
        
        #------------------------rewrite info with real layers resistivities if provided and layers names ---
        # self.input_layers, _  , _= Geodrill.get_geo_formation_properties(structures_resistivities= self.input_resistivities,)
        
        
        #-------------------Print Info -------------------------------------------------
        
        print('**{0:<37} {1} {2}'.format(' Number of layers','=' , len(self.input_layers) ))
        
        print('-'*77)
        print('{0:<25}{1:<25} {2:<25}'.format('Structure', 'Rho mean value (Ω.m)', 
                                             'Rho range  (Ω.m)'))
        print('-'*77)
  
    
        for ij  , (ires, geos) in enumerate(zip (self.input_resistivities,self.input_layers )): 
            if len(self.input_resistivities ) >1 :
                if ij==0 : 
                    rho_rg = '{0}-{1}'.format(ires, self.input_resistivities [ij+1])
                    rho_mean = np.around((ires+ self.input_resistivities [ij+1])/2,2)
                elif ires == self.input_resistivities[-1]: 
                     rho_rg = '{0}-{1}'.format(self.input_resistivities[ij-1], self.input_resistivities[ij])
                     rho_mean =np.around( (self.input_resistivities[ij-1]+ self.input_resistivities[ij])/2,2)
                else : 
                    rho_rg = '{0}-{1}'.format(self.input_resistivities[ij-1], self.input_resistivities[ij+1])
                    rho_mean = np.around( (self.input_resistivities[ij-1]+ self.input_resistivities[ij+1])/2,2)
            if len(self.input_resistivities)==1:
                rho_rg = '{0}'.format(ires)
                rho_mean = ires
                
            print('{0:<25}{1:<25} {2:>25}'.format(geos, rho_mean,rho_rg ))
        print('-'*77) 
        
        
    @staticmethod 
    def get_geo_formation_properties (structures_resistivities, real_layer_names=None,
                                      constrained_electrical_properties_of_rocks=True, 
                                       **kwargs
                                      ):
        """
        Get the list of stuctures and their names of after replaced , flexible tools.  
        wherever structures provided, the name, color ,  as well as the pattern .
        if constrained electrical properties if True , will keep the resistivities with their 
        corresponding layers as reference. If the layer names if found on the data Base then , 
        will return its pattern and color else defaultcolor is black and patter is "+.-".
        if constrained_electrical_properties_of_rocks is False , will check under data base 
        to find the resistivities that match better the layers 
        
        Parameters
        ------------
            * structures resistivities  : array_like, 
                                       resistivities of structures 
                                        
            * real_layer_names : array_like |list 
                               names of layer of survey area 
                               if not provided , will use resistivities to find the closet 
                               layer that match the best the resistivities
                                
            * constrained_electrical_properties_of_rocks: bool 
                        set to True mean the realy_layer is provided. if not program 
                        will enforce to False , will use default conventional layers
                        Default is False, assume to povided layer names for  accuracy
        
        Returns
        --------
             f_name: array_like 
                  names of formations find with their corresponding rho 
             f_pattern: array_like 
                 patter of differents geological formations 
             f_color: array_like   
                 color of differents geological formations 
        """

        unknow_layer_color = kwargs.pop('default_layer_color', 'white')
        unknow_layer_pattern =  kwargs.pop('default_layer_pattern', ".//")
        
        find_pattern , __f_db= False, -1             # flag to find pattern , flag to raise error for unsuccessufful connection to database 
        pattern, geo_color = [[] for i in range(2)] 
        
        geoformation_obj =STRL.Geo_formation()
        geof_names = geoformation_obj.names 
        
        
        if constrained_electrical_properties_of_rocks is True : 
            if real_layer_names is not None : #check weither the layer of resistivities
                                    #are the same dimension and return new layer with same dimension with rho
                real_layer_names = ascertain_layers_with_its_resistivities(real_layer_names =real_layer_names,
                                                                       real_layer_resistivities= structures_resistivities)
            else : 
                msg =''.join(["{constrained_electrical_properties_of_rocks} argument is set to <True> as default value ",
                              "  unfortunately , you did not provide any layer's names with its resistivity values.", 
                              " We can not set layer' resistivities and layer' names as reference data.",
                              " We will ressetting {constrained_electrical_properties_of_rocks} to False. However", 
                              " be sure that layer's provided automatically don't match exactly the underground informations.", 
                              "To accurate informations, you need ABSOLUTELY to provided layer' resistivities as well as its names get on the field or any other firms."])
                
                print('---> Important Note :!' )
                text = punc.fmt_text(data_text=msg, fmt='+',return_to_line=90)
                print(text)
                warnings.warn(msg)
                
                constrained_electrical_properties_of_rocks = False 
           
        if constrained_electrical_properties_of_rocks is False : 
            if real_layer_names is None : 
                real_layer_names = Geodrill.get_structure(structures_resistivities)
            elif real_layer_names is not None :  # enforce to use the layer's names provided by the user (constrained is True) 
                real_layer_names = ascertain_layers_with_its_resistivities(real_layer_names =real_layer_names,
                                                                       real_layer_resistivities= structures_resistivities) 

        #-------------------------------------CONNECT TO GEODATABASE -----------------------------
        # connect to geodataBase 
        try : 
            geodatabase_obj = GeoDataBase ()
    
        except :
            __f_db= 0 
        #-----------------------SEARCH ON DATABASE -------------------*
        if geodatabase_obj.success == 1 :
            mess = '----> Successfull connexion to geoDataBase !'
            print(mess)
            _logger.info(mess)
            
            for res,  lnames  in zip ( structures_resistivities, real_layer_names ) :
                    geodatabase_obj. _get_geo_structure(structure_name =lnames) 
                    if geodatabase_obj.geo_structure_exists is True :  # ckek geostructure existence
                    
                        if geodatabase_obj.colorMPL != 'none':  
                            # rebuild the tuple of mpl colors blocked as string in database  
                            lnames_color = tuple([ float(ss) for ss in 
                                                    geodatabase_obj.colorMPL.replace('(', '').replace(')', '').split(',')])
                            
                            geo_color.append(lnames_color) # append rgb colors not include alpha 
                        
                        else : 
                            mess = 'Sorry ! {0} matplotlib color is  not filled yet in our DataBase.'.format(lnames)
                            warnings.warn(mess)
                            geo_color.append(unknow_layer_color) # then append default colors 
                            
                        if geodatabase_obj.hatch != 'none': 
                            
                            pattern.append(geodatabase_obj.hatch.replace('(', '').replace(')', ''))
                        
                        else : 
                            mess = 'Sorry ! {0} matplotlib pattern is not filled yet in our DataBase.'.format(lnames)
                            warnings.warn(mess)
                            pattern.append(unknow_layer_pattern)
                        
                    elif geodatabase_obj.geo_structure_exists is False   : 
                        # then aborted the process and go to search into strata and strutral modules 
                        __f_db=0
                        
                        mess =' Actually {0} does not exist in our dataBase,'\
                            ' we alternatively use other ways to suitable respond to your request.'.format(lnames)
                        warnings.warn(mess)
                        _logger.debug(mess)
                        break
   
            geodatabase_obj.manage_geoDataBase.closeDB() # close de DB
            
        #-------------------------SEARCH ON STRATA AND STRUCTURAL MODULES ----------------------*
        if geodatabase_obj.success ==0 or  __f_db == 0 : 
            geo_color , pattern = [[] for i in range(2)] # clean initialise  and restart searching 
            
            for res,  lnames  in zip ( structures_resistivities, real_layer_names ) :
                if lnames in geof_names :  # if exist on geoformation array names then get the index and the coorsponding values 
                   
                   #-------------------------------------------------------------------
                   if ' ' in lnames  : new_lnames = lnames.replace(' ', '_') # be sure to find attribute names
                   else : new_lnames =lnames
                   
                   col = getattr(geoformation_obj, new_lnames)['color'] # for consistency 
                   geo_color.append(col)
    
                elif lnames not in geof_names : #check wether the names belongs to structural pattern 
                   if lnames.lower() in STRL.geo_pattern.pattern.keys():
                       geo_color.append(STRL.geo_pattern.pattern[lnames.lower()][1])
                   else :
                       geo_color.append( unknow_layer_color)
      
                 # actually will use the default pattern : Then check the resistivities inot geoelectrical 
                 #    propertty. if structures found , then took it pattern 
                  
                for ii, ( keyprop , resprops)  in enumerate( Geodrill.geo_rocks_properties.items ()) : 
                    if min(resprops) <= res <= max(resprops) : #resprops)[-1] <= res <= resprops[0] :
                        if '/' in keyprop:  # split and search for their appropriate pattern
                            if res <= np.array(resprops).mean() : keyprop=keyprop.split('/')[1] # take the first na
                            else : keyprop=keyprop.split('/')[0]
                            
                        pat =STRL.geo_pattern.pattern[keyprop][0]       # keep it pattern
                        pattern.append(pat)
                        find_pattern =True                           # switch on flag 
                        break  # dont need to continue to take other pattern 
                        
                    if ii == len(Geodrill.geo_rocks_properties.items ()) and find_pattern is False : # mean no pattern 
                    # is found after looping all items of dictionnary 
                        pattern.append( unknow_layer_pattern)
                        
                find_pattern =False      # resetting flag to False (switch off )

        return real_layer_names, geo_color , pattern      
                       
    def to_golden_software(self, input_resistivities=None, input_layers =None, 
                               step_descent = None,  filename =None,
                               savepath =None , **kwargs ): 
        """
        Output average files, rehoreplaced files, spseudosequennces
        files and station locations files  will generate 3 outputs for 
        Golden software plots 
        
        1. One model for rho averaged (*_aver), transitory data betwen 
            the calculated rho and true rho.
            
        2. second for rho value replaced . Replaced calcualted model 
            structures resistivities 
            
        to their closest resistivities as resistivities reference from
            input resistivities.
            
        3. the most important files: the  cut out resistivities value 
            with step descent (._sd) show most dominant stratigraphy sequences .
            
        4. the station location file (.bln)
    
        Plot the 3 files in Golden software to see transition from model 
        calculation to truthsequence detail models which is most closest to reality.

        Parameters
        ------------
            * input_resistivities :   array_like 
                Truth values of resistivities 
                    
            * filename :str 
                name of output file 
                    
            * step_descent :  float,
                Step to cut out data and to  force resistivites calcualted 
                 to match the reference data as input resistivities 
                 if not provided the step will be 20% of D0I.
                 
            * input_layers  : array_like 
                    True input_layers names ( geological informations of
                     encountered layers )
                    
            * savepath :  str, 
                    full path to the savepath ,
                    if None , will create folder name to savepath 
           
        Returns 
        ----------
            obj , str 
             golden software outputfiles `_rr.dat`, `_aver.dat`, `_sd.dat`. 
         
        Holding a followings informations:
            
        ===================  ================  ================================
        Optional params          Type                   Explanation 
        ===================  ================  ================================
        elevation              array_like       elevation of survey area 
        to_negative_depth      bool             export deth in negative 
                                                value or positive 
                                                default is "negative".
        scale                  str              scale to export data .Must be
                                                *m* or `km`. *default* is **m**       
        ===================  ================  ================================
        
        """
        elevation = kwargs.pop('elevation', None)
        write_negative_depth=kwargs.pop('to_negative_depth', True)
        scale =kwargs.pop('scale', 'm')
        if scale is None : scale= 'm'
        
        # rescale data 
        if scale =='m': 
            df =1. 
        elif scale =='km': 
            df =1000. 
        else : df =1.
        
        if input_layers is not None : self.input_layers = input_layers 
        if input_resistivities is not None : 
            self.input_resistivities  = input_resistivities 
        if step_descent is not None : self.step_descent= step_descent 
        if elevation is not None : self.elevation = elevation 
        
        
        if self.elevation is None : 
            mess='!Elevation is not provided. We gonna set to 0.'
            print('-->'+mess)
            self._logging.debug(mess)
            
            if self.station_location is not None :  # assume station location exist
                self.elevation  = np.repeat(0., len(self.station_location))
        elif self.elevation is not None : # check elevation length with station location 
            
            self.elevation = geo_length_checker(main_param= self.station_location,
                                                optional_param = self.elevation, 
                param_names = ('station location', 'elevation'), fill_value=0.)
        
        # now read files 
        if  self.input_resistivities is not None : 
            self.geo_build_strata_logs()
        else : 
            mess ="".join(["Need ABSOLUTELY an input resistivities get"
                           " on the field or other companies ", 
                           " Input resistivities MUST provided in order"
                           " to take data as reference resistivities.", 
                           " please provided a least a truth resitivities of one layer."])
            
            warnings.warn(mess)
            self._logging.error(mess)
            raise CSex.pyCSAMTError_plot_geoinputargument("Can not write details sequences log files "
                                                          "! Please provided at least one"
                                                          " truth layer resistivity.")
        
        matrix_rhoaver ,matrix_rhorr, matrix_rho_stepdescent =[[]for i in range(3)]
        for site in self.station_names : 
            if site ==self.station_names[0]: # station id to read is the first index =0
                index =0 
            if site ==self.station_names [-1] : index =2 # station id  to read is the 2index 
            else : index =1  # read th middle array 
    
            matrix_rhorr.append(self.geo_drr[site][:, index])
            matrix_rhoaver.append(self.geo_daver[site][:, index])
            matrix_rho_stepdescent.append(self.geo_dstep_descent[site][:, index ])
                
        # build a matrix of all data collected 

        matrix_rhorr= func.concat_array_from_list(matrix_rhorr, concat_axis=1)
        matrix_rhoaver= func.concat_array_from_list(matrix_rhoaver, concat_axis=1)
        matrix_rho_stepdescent= func.concat_array_from_list( matrix_rho_stepdescent, concat_axis=1)
        
        # create empty list to collect infos of read values
        write_averlines , write_rrlines, write_stepdescent_lines , write_blnfiles =[[]for i in range(4)]
        
        # writes multilines in the same times 
        
        if write_negative_depth  is True : self.geo_depth  *=-1
                # scalled station location and geo_depth 
        self.geo_depth /=df 
        self.station_location /=df
        
        for writes_lines , matrix in zip ( [write_averlines , write_rrlines, write_stepdescent_lines], 
                                          [matrix_rhoaver, matrix_rhorr, matrix_rho_stepdescent]):
            
            for ii in range(len(self.geo_depth)) : 
                for jj in range(len(self.station_names)): 
                    writes_lines.append(''.join([
                        '{0:>15.7f}'.format(self.station_location[jj]),
                        '{0:>15.7f}'.format(self.geo_depth[ii]), 
                        '{0:>15.7f}'.format(matrix[ii,jj]),
                        '\n',
                        ]))
                   
                
      # create fileame if not provided 
        if filename is None : 
            filename = '.{0}'.format(datetime.datetime.now().month)
        else : 
            filename ='{0}.{1}'.format(os.path.basename(filename), datetime.datetime.now().month)
         
        for  ikey, stn in enumerate(self.station_names) : 
            write_blnfiles .append(''.join(['{0:<7.6f},'.format(self.station_location[ikey]), 
                                                '{0:<7.6f},'.format(self.elevation[ikey]),
                                                '{0:<4}'.format(stn), '\n']))
            
        if savepath is None : # create a folder in your current work directory
            try :
                savepath = os.path.join(os.getcwd(), '_outputGeoSD_')
                if not os.path.isdir(savepath):
                    os.mkdir(savepath)#  mode =0o666)
            except : 
                warnings.warn("It seems the path already exists !")
                
        #writes files 
        for ii , (tfiles , wfiles) in enumerate(zip(['_aver', '_rr', '_sd', '_yb'], 
                                                    [write_averlines, write_rrlines, 
                                                     write_stepdescent_lines,  write_blnfiles])): 
            if ii == 3 : mm='bln'
            else :mm='dat'
            with open(''.join([filename,'{0}.{1}'.format(tfiles, mm)]), 'w') as fw : 
                fw.writelines(wfiles)
                
        #savefile inot it savepath 
        if savepath is not None : 
            import shutil 
            try :
                for jj, file in enumerate( [filename + '_aver.dat', filename +'_rr.dat',
                                            filename +'_sd.dat', filename + '_yb.bln']):
                    shutil.move(os.path.join(os.getcwd(),file) ,
                        os.path.join(savepath , file))
            except : 
                warnings.warn("It seems the files already exists !")
                
        filenames =[filename + '_aver.dat', filename +'_rr.dat',
                                            filename +'_sd.dat', filename + '_yb.bln']   

        print('---> geo output files {0}, {1}, {2}  & {3}'
              ' have been successfully written to  <{4}>.'.format(*filenames, savepath))
           
    def to_oasis_montaj (self,  model_fn=None , iter_fn=None ,profile_fn =None, 
                         mesh_fn=None  , data_fn=None , filename=None,
                         savepath =None , **kwargs) : 
        
        """
        write output to oasis montaj when station loocation and profile 
        coordinates are provided . We assune before using this method , you are
        already the coordinates files at disposal(*stn), if not  use the method 
        `to_golden_software`.Coordinated files are Easting Northing value  
        not in degree decimals . It uses occam 2D outputfiles or Bo Yang iter2Dat file 
        also add profile XY coordinates (utm_zone ).
        
        Parameters
        ------------
            * model_fn :  str  
                    full path to Occam2D model file 
                    
            * mesh_fn :  str   
                    full path to Occam2D mesh file 
                    
            * data_fn :   str    
                    full path to Occam2D data file 
                    
            * iter_fn :   str    
                    full path to Occam2D iteration file 
                    OR Bo Yang (x, y, z ) files  
                                
            * iter2dat_fn : str      
                    full path to iter2dat file 
                    (see _occam2d module to see which file is it. or call  
                     occam2d.Iter2Dat.__doc__)
                        
            * bln : str  
                    full path to satation location file 
                    profile files (Easting , Northing Coordinates) , 
                    "see :ref:module-cs"
                                        
            * profile_fn: str   
                    full path to station profile file . You can   
                    useProfile module to rewrite _coordinate files 
                     
            * write_negative_depth:  bool 
                            output negative depth. 
                            *Default* is True
                                
            * scalled_east_north:  tuple
                            scalled the easting and northing. 
                            Substract or add value to easting
                            or northing values. first `index1` 
                            equal easting and  `index2` equal 
                            Nothing.*Default* is (0,0).

         Returns 
         --------
            obj , str 
             oasis outputfiles `.xlxs`, `.csv`
         
        Holding a followings informations:
            
        =================  ==============  ====================================
        Optional params         Type            Explanation 
        =================  ==============  ====================================
        filename            str             New name of output file .
                                            *Default* is None 
        normalize_depth     bool            Set the depth ega spacing depth .
                                            In fact for oasis montaj ,depth 
                                            must be equidistant  the same
                                            spacing when image in deeper.
                                            *default* is True. if false , 
                                            will generate depth as Occam2D 
                                            mesh z nodes. 
        easting             array_like      Easting UTM coordinates .Profile
                                            *Stn* file is  provided,
                                            no need to input. *default* is None 
        northing            array_like      Northing coordinates .  If station
                                            coordinates is provided from *stn 
                                            file .It will detect  automatically.
                                            *default* is None.
        elevation           array_like      Elevation area . *Default* is None.
        output_s_XY         bool            Scalled coordinates  output files. 
                                            *Default* is True 
        input_rho           array_like      input resistivities of True
                                            geological formations (optional)
                                            if input resistivities is provided 
                                            will generated output step descent 
                                            files (_sd),  roughness_rho (_rr),
                                            rho_averaged (_aver). *Default* is 
                                            output the main resistivty model. 
        input_layers         array_like     list of true names of layers(opt.) 
        step_descent         float          Step to cut out data and to  
                                            force resistivites calcualted 
                                            to match the reference data as 
                                            input resistivities if not provided 
                                            the step will be 20% of D0I.(opt)          
        writeType            str            Writer format *.csv or *.xlsx .
                                            *Default* is *.xlsx
        add_header           bool           Add head on exported sheet. 
                                            set False to mask heads.
                                            *Default* is True.
        csv_separateType     str            Indicated for csv exported files.
                                            The type of comma delimited.
                                            *Defaut* is ','
        to_log10             bool           if True : will ouput all 
                                            resistivities data to  log10 values         
        =================  ==============  ====================================
        
        """
        def create_model_matrix (station_names, geo_dict_rho, transpose =True): 
            """
            simple function to generate matrix from geo_dictionnary data 
            main site is frames into3 stations 
            
            Parameters
            ----------
                * station _names: list  
                            name of stations 
                * geo_dict_rho: dict 
                            dictionnary of model_resistivities 
                             can  be a step _descent dict , 
                             geo_average_or geo_replace rho 
                * transpose:  bool 
                             if transpose is True , depth will be a xoffset 
                             and and station names as y offsets
                                
            Returns 
            --------
                ndarray 
                    model_rho_matrix , matrix , depth & station locations 
            """

            model_rho_matrix =[]
            for site in station_names : 
                if site ==station_names[0]: # station id to read is the first index =0
                    index =0 
                if site ==station_names [-1] : index =2 # station id  to read is the 2index 
                else : index =1  # read th middle array 
                model_rho_matrix.append(geo_dict_rho[site][:, index])
                
            model_rho_matrix= func.concat_array_from_list(model_rho_matrix , concat_axis=1)
            if transpose is True : model_rho_matrix= model_rho_matrix.T
            
            return model_rho_matrix
            
            
        print('-'*77)
        normalize_depth =kwargs.pop('normalize_depth', True )
        iter2dat_fn =kwargs.pop('iter2dat_fn', None)
        bln_fn =kwargs.pop('bln_fn', None)
        easting =kwargs.pop('easting', None)
        northing =kwargs.pop('northing', None)
        elevation = kwargs.pop('elevation', None)
        scalled_east_north =kwargs.pop('scalled_east_north', (0,0))
        output_sprofile =kwargs.pop('output_s_XY', True)
        write_negative_depth=kwargs.pop('to_negative_depth', True)
        to_log10 = kwargs.pop('to_log10', True)
        
        # write -externals files : Step descent , average rho 

        input_rho =kwargs.pop('input_resistivities', None)
        input_layers =kwargs.pop('input_layers', None)
        step_descent =kwargs.pop('step_descent', None)
 
        if input_rho is not None :  self.input_resistivities  = input_rho 
        
        if self.input_resistivities is not None : 
            write_external_files=True 
        else : write_external_files = False
            

        # write sheet 
        writeType =kwargs.pop('writeType', 'xlsx')
        csvsetHeader =kwargs.pop('add_header',True)
        csvsep =kwargs.pop('csv_separateType',',')
        writeIndex=kwargs.pop('write_index_on_sheet',False)
    
 
        # create filename if not provided 
        if filename is None and profile_fn is not None : 
            filename = os.path.basename(profile_fn)
            
        if filename is None : 
            filename = '.{0}'.format(datetime.datetime.now().month)
        else: 
            filename ='{0}.{1}'.format(os.path.basename(filename),
                                       datetime.datetime.now().month)

        # statements of arguments ----

        for mattar, mfile in zip(['model_fn', 'mesh_fn', 'data_fn',
                                  'iter_fn', 'iter2dat_fn', 'bln_fn' ],   
                                  [model_fn, mesh_fn, data_fn, iter_fn, iter2dat_fn, bln_fn ]):
            if mfile is not None : 
                setattr(self, mattar, mfile)
        
        # -----------------------------------Buil profile object ------------------------------------------
        if elevation is not None : self.elevation = elevation 
        profile_obj = Profile ()
        profile_obj.read_stnprofile(profile_fn , easting=easting, northing =northing , 
                                                 elevation = self.elevation)
        profile_obj.straighten_profileline ( reajust=scalled_east_north , output =output_sprofile)
        
        print('** {0:<37} {1} {2} {3}'.format('Dipole length','=', profile_obj.dipole_length, 'm.'  ))
        
        profile_angle, geo_electric_strike = profile_obj.get_profile_angle()
        
        print('** {0:<37} {1} {2} {3}'.format(' Profile angle ','=',
                                              round(profile_angle, 2), 'deg N.E' ))
        print('** {0:<37} {1} {2} {3}'.format('Geo_electric strike','=',
                                              round(geo_electric_strike , 2), 'deg N.E' ))
        
        # In the case where elevation is contain on the profile stn file , get elevation 
        if profile_obj.elev is not None : self.elevation = profile_obj.elev 
        
        # build section Infos 
        #Stations	Easting_X_m	Northing_Y_m	Elev_H_m	x_m	, Norm_h_m	offsets_m	DOI_max_m
        info_list = ['station', 'easting_X_m', 'northing_Y_m', 'elev_H_m', 'position_x_m', 
                     'nivelize_h_m', 'offset_m', 'doi']
        #normalH 
        nivelize_elevation = np.around(profile_obj.elev - profile_obj.elev.min(), 2)
        # stn position 
        if write_negative_depth is True : doi =-1* self.doi 
        else : doi = self.doi 
        infos_oasis_montaj = func.concat_array_from_list(list_of_array= [self.station_names ,
                                                                  profile_obj.east , 
                                                                  profile_obj.north , 
                                                                  profile_obj.elev , 
                                                                  profile_obj.stn_position , 
                                                                  nivelize_elevation , 
                                                                  self.station_location , 
                                                                  np.full((self.station_location.shape[0],), doi),
                                                                  ], 
                                                                 concat_axis=1)
        #-----------------------Build main data for oasis -----------------------------------------------------
        
        # manage the depth info and create matrix of depth 
        spacing_depth =abs(self.geo_depth.max()-self.geo_depth.min())/ (len(self.geo_depth)-1)
        print('** {0:<37} {1} {2} {3}'.format('Spacing depth','~=', round(spacing_depth,2), 'm.' ))
        
        if normalize_depth is True : 
            # get the step_deth 
            geo_depth = np.around(np.linspace(self.geo_depth [0], self.geo_depth[-1], len(self.geo_depth )))
            spacing_depth =abs(self.geo_depth.max()- self.geo_depth.min())/ (len(self.geo_depth)-1)
    
            print('---> Depth normalized !')
            print('---> {0:<37} {1} {2} {3}.'.format('new spacing depth',
                                                     '=', round(spacing_depth), 'm' ))
        else :
            print('---> UNnormalized  depth!')
           
        depth_oasis_montaj  = np.repeat(geo_depth.reshape(geo_depth.shape[0], 1),
                              len(self.station_location), axis =1 ) # build a mtrix of depth
        depth_oasis_montaj = depth_oasis_montaj .T          # Transpose the depth 
        
        if write_negative_depth is True :depth_oasis_montaj *= -1
            
        # build data 
        
        model_oasis_montaj = create_model_matrix (station_names = self.station_names, 
                                                  geo_dict_rho= self.geo_d)
        
        if to_log10: np.log10( model_oasis_montaj )
        
        data_oasis = np.concatenate((infos_oasis_montaj,
                                     depth_oasis_montaj, 
                                     model_oasis_montaj), axis =1)
        # buidd pandas columns 
        depres_pandas_columns = info_list + ['dep_{0}'.format(int(dd)) for dd in geo_depth] +\
            ['res_{0}'.format(int(dd)) for dd in geo_depth] 
            
            
        oas_pandas = pd.DataFrame(data=data_oasis, columns=depres_pandas_columns)
        
        fex=0  # falg to write external files 
        if write_external_files is True :
            if input_layers is not None : self.input_layers = input_layers
            if step_descent is not None : self.step_descent = step_descent 
            if self.input_resistivities is None : 
                mess = ''.join(['! No input resistivities provided! Could not write external geo-files.', 
                                ' If you expect to write external GEO {step_descent  - roughness  and average_rho} files,', 
                                ' You ABSOLUTELY need  to provided at least a truth resistivity values', 
                                ' of some geological formation of the area.'])
                print(punc.fmt_text(mess,fmt='+'))
                
                warnings.warn(mess)
                self._logging.debug (mess)
            if self.input_resistivities is not None : 
                self.geo_build_strata_logs(input_resistivities= self.input_resistivities, 
                                           input_layers=self.input_layers,
                                           step_descent=step_descent)
             
            # now build matrix of all external files 
            model_step_descent = create_model_matrix (station_names = self.station_names, 
                                                  geo_dict_rho= self.geo_dstep_descent)   
            model_geo_roughness = create_model_matrix (station_names = self.station_names, 
                                                  geo_dict_rho= self.geo_drr)
            model_average_rho = create_model_matrix (station_names = self.station_names, 
                                                  geo_dict_rho= self.geo_daver)
            
            if to_log10 : 
                model_step_descent= np.log10(model_step_descent)
                model_geo_roughness=np.log10(model_geo_roughness)
                model_average_rho =np.log10( model_average_rho)
            
            data_oas_stepDescent = np.concatenate((infos_oasis_montaj,
                                                   depth_oasis_montaj,
                                                   model_step_descent), axis =1)
            data_oas_roughness = np.concatenate((infos_oasis_montaj,
                                                 depth_oasis_montaj,
                                                 model_geo_roughness), axis =1)
            data_oas_averageRho = np.concatenate((infos_oasis_montaj,
                                                  depth_oasis_montaj,
                                                  model_average_rho), axis =1)
            # create pandas dataFrame 
            
            STD_pandas = pd.DataFrame(data=data_oas_stepDescent,
                                      columns=depres_pandas_columns)
            ROUGH_pandas = pd.DataFrame(data=data_oas_roughness,
                                        columns=depres_pandas_columns)
            AVER_pandas = pd.DataFrame(data=data_oas_averageRho,
                                       columns=depres_pandas_columns)
            
            
            
            fex =1 # everything where ok then take one 
        # write files 
        # if f== 1 : 
            # create filename 
        if writeType.lower().find('exc')>= 0 or writeType.lower().find('xls')>=0 :
            writeType = '.xlsx'
            if fex ==1 :  # write all files into one worksheet 
                with pd.ExcelWriter(''.join([filename + '_cor_oas',writeType])) as writer :
                    for  sname , df_ in zip(['.main.', '_sd', '_rr', '_aver'],
                                            [oas_pandas, STD_pandas,ROUGH_pandas, AVER_pandas ]):
                        df_.to_excel(writer,sheet_name=sname, index =writeIndex)
                        
                    fex =2 
            else : 
     
                oas_pandas.to_excel(''.join([filename,'.main._cor_oas', writeType]),
                                    sheet_name=filename +'main_cor_oas', index=writeIndex)  

            
        elif writeType.lower().find('csv')>= 0 : 
            writeType ='.csv'
            oas_pandas.to_csv(''.join([filename,'.main._cor_oas',writeType]),
                              header=csvsetHeader,
                  index =writeIndex,sep=csvsep) 
            if fex==1 : 
                for fnx, _df in zip( ['_sd', '_rr', '_aver'] ,
                                    [STD_pandas,ROUGH_pandas, AVER_pandas ]):
                    
                    _df.to_excel(''.join([filename + fnx,'_cor_oas', writeType]),
                                sheet_name=filename +fnx +'_cor_oas', index=writeIndex) 
                    
        else :
            mess = 'The type you provide is wrong. Support only *.xlsx and *.csv format'
            self._logging.error (mess)
            warnings.warn (mess)
            raise CSex.pyCSAMTError_geodrill_inputarguments('wrong format ={0} !'\
                                                            ' Could not write geo file to oasis.'\
                                                                ' Support only *.xlsx or *.csv format ')

        # export to savepath 
        if savepath is None : # create a folder in your current work directory
            try :
                savepath = os.path.join(os.getcwd(), '_output2Oasis_')
                if not os.path.isdir(savepath):
                    os.mkdir(savepath)#  mode =0o666)
            except : 
                warnings.warn("It seems the path already exists !")
                
        #savefile inot it savepath 
        if savepath is not None : 
            import shutil 
            try :
                if fex ==2 : 
                   shutil.move(os.path.join(os.getcwd(),''.join([filename,'_cor_oas', writeType])) ,
                   os.path.join(savepath , ''.join([filename,'_cor_oas', writeType])))
                else : 
                    for jj, file in enumerate(  ['.main.', '_sd', '_rr', '_aver']):
                        if fex ==0 : 
                            if jj ==1 : break # donc continue because other files does not exists  
                        shutil.move(os.path.join(os.getcwd(),
                                                 ''.join([filename, file, 
                                                          '_cor_oas', writeType])) ,
                            os.path.join(savepath , 
                                         ''.join([filename, file, 
                                                  '_cor_oas', writeType])))
            except : 
                warnings.warn("It seems the files already exists !")
        
        # write Infos 
        print('** {0:<37} {1} {2} '.format('number of stations',
                                           '=', len(self.station_names)))
        print('** {0:<37} {1} {2} {3}'.format('minimum offset',
                                              '=', self.station_location.min(), 'm' ))
        print('** {0:<37} {1} {2} {3}'.format('maximum offset',
                                              '=', self.station_location.max(), 'm' ))
        print('** {0:<37} {1} {2} {3}'.format('maximum depth',
                                              '=', geo_depth.max(), 'm' ))
        print('** {0:<37} {1} {2} {3}'.format('spacing depth ',
                                              '=', round(spacing_depth,2), 'm' ))
        if self.elevation is None or np.all(self.elevation ==0.) : 
            print('--->  Elevation no added !')
        else : 

            print('** {0:<37} {1} {2} {3}'.format('minumum elevation',
                                                  '=', self.elevation.min(), 'm' ))
            print('** {0:<37} {1} {2} {3}'.format('maximum elevation ',
                                                  '=', self.elevation.max(), 'm' ))
        
        print('** {0:<37} {1} {2} {3}'.format('minumum resistivity value','=',
                                              round(model_oasis_montaj.min(), 2), 'Ω.m' ))
        print('** {0:<37} {1} {2} {3}'.format('maximum resistivity value','=',
                                              round(model_oasis_montaj.max(), 2), 'Ω.m' ))
   
        print('** {0:<37} {1} {2} '.format('Lowest sation','=',self.station_names[ 
                                    int(np.where(nivelize_elevation==nivelize_elevation.min())[0])] ))
        
        print('** {0:<37} {1} {2}'.format('Highest site','=', self.station_names[ 
                                    int(np.where(nivelize_elevation==nivelize_elevation.max())[0])]))
        print('** {0:<37} {1} {2} {3}'.format('Altitude gap','=', nivelize_elevation.max(), 'm' ))
        print('** {0:<37} {1} {2}'.format('Number of running ','=', 
                                          len(depres_pandas_columns) * len(self.station_location)))
        
        
        if fex ==0 : 
            print('---> geo output file {0},  has been successfully written to  <{1}>.'.\
                  format(''.join([filename, '.main._cor_oas', writeType]), savepath))
        elif fex ==1 : #read external files 

            filenames =[ filename + file +'_cor_oas'+ writeType for file ,
                        exten in zip(['main.', '_sd', '_rr', '_aver'], 
                                     [writeType for i in range(4)])]
                        

            print('---> geo output files {0}, {1}, {2}  & {3} have been successfully '\
                  'written to  <{4}>.'.format(*filenames, savepath))
                

        print('-'*77)
        
          
class Geosurface :
    """
    Read Multidata from oasis montaj output files  generated by `geodrill`
    module . Class to Build  a depth surface map for imaging . 
 
    Arguments
    -----------
        **path**: string 
                path to oasis ouput files , frequently the files generated by
                `geodrill`modules.
       
    Holding a followings informations:
        
    ==============  =============  ============================================
    KeyWords        Type                        Description    
    ==============  =============  ============================================
    depth _values   array_like      Values of depth for imaging .
                                    can be a list of depth like [38,912]
    fileformat      str             output format . actually geosurface can 
                                    only generate ouput in `csv` and `xls`.   
                                    *default* is   `csv`. 
    savepath        str             full path to sve directory . If  none 
                                    will create an a new directory 
    ==============  =============  ============================================
        
    :Example:
        
        >>> from geodrill.geoCore.geodrill import Geosurface 
        >>> gs_obj = geosurface (path =os.path.join(os.environ['pyCSAMT'], 
        ...                                   'geodrill', 'data', 
        ...                                   InputOas), 
        ...                            depth_values = [28, 100])
        >>> gs_obj.write_file()
    """
    geo_surface_format=["csv",  "xlsx", "json","html","sql"] 
    
    read_dico = {".csv":pd.read_csv,  ".xlsx":pd.read_excel,".json":pd.read_json,
                 ".html":pd.read_json,".sql" : pd.read_sql}  
    
    
    def __init__(self , path=None, **kwargs):
        self._logging =csamtpylog.get_csamtpy_logger(self.__class__.__name__)
        self.path=path
        self.depth_values =kwargs.pop('depth_values', None)
        self.export_format =kwargs.pop ('fileformat', 'xlsx')
        self.savepath =kwargs.pop('savepath', None)

    
        if self.path is not None : 
            self._validate_oasis_path() 
               
                        
    def _validate_oasis_path(self, path =None ) : 
        """
        Validate oasis montaj files and get extension files
        
        :param path: full path to oasis montaj files
        :type path: str 
        
        """  
        if path is not None : self.path =path 
        if self.path is not None : 
            if os.path.isdir(self.path):  # get the list of files in folder 
                # get file and be sure that format exist in 
                self.oasis_data = [os.path.join(self.path, file) for file in os.listdir(self.path)
                                   if ( file.split('.')[-1]  in self.geo_surface_format)]
                
                # print(self.oasis_data)
                if self.oasis_data is None or len(self.oasis_data)==0 : 
                    mess ='No files detected!. Please provided right path to oasis models files.'
                    warnings.warn(mess)
                    self._logging.error(mess)
                    raise CSex.pyCSAMTError_inputarguments(mess)
      
        elif self.path is None :
            raise CSex.pyCSAMTError_inputarguments('No path found ! Please provided a right path.')
        
        # get extension file 
        extension_files =[os.path.splitext(pathfile)[1] for pathfile in self.oasis_data]
        self.extension_file =extension_files
        if self.extension_file.replace('.','') not in self.geo_surface_format : 
            mess = 'Unacceptable format = {0}. Could read format ={1}'.\
                format(self.extension_file, tuple(self.geo_surface_format))
            self._logging.warning(mess)
            raise CSex.pyCSAMTError_file_handling(mess)
        
        
    @property 
    def extension_file (self):
        return self._extension_file
    
    @extension_file.setter 
    def extension_file(self,  get_extension): 
        # put on array and use np.unique to get the most frequent format
        if isinstance(get_extension, (list, tuple)): 
            get_extension = np.array(get_extension)
        ex , ex_count =np.unique(get_extension , return_counts= True)
        self._extension_file  = ex[int(np.where (ex_count== ex_count.max())[0])]
        
    def read_oasis_files (self, path =None ):
        """
        Method to get depth spacing , station info  data infos, 
        Each line becomes it own attributes components
        of info values are  `Stations`,	`Easting_X_m`,`Northing_Y_m`,
        `v_H_m`,`x_m`	`Norm_h_m`,`sets_m`	`DOI_max_m`.

        :param path: full path to the oasis files .  
        :type path: str  
 
        :Example:
            
             >>> from geodrill.geoCore.geodrill  import Geosurface 
             >>> path =  r'F:/__main__csamt__\oasis data\OASISWORKS\all_data'
             >>> geo_surface_obj = Geosurface( path =path )
             >>> geo_surface_obj.read_oasis_files()
             >>> geofilenames = geo_surface_obj.filenames 
         
        ... note::To get the values of line K1_cor_oas "K1_cor_oas.csv
                do ``k1_obj = geo_surface_obj.K1_cor_oas`` 
        """
        print('{0:-^77}'.format('GeoSurface * Data * info'))
        
        def get_depth_values(df ): 
            """
            get specific depth values from geomodel dataframe
            
            Parameters
            ----------
                * df : pandas.Core.DataFrame 
                        Dataframe pandas 
            Returns
            --------
                 info_names: list 
                          header names info 
                 depth_values: array_like 
                            spacing depth offsets 
                 depth_spacing: float 
                          value of spacing depth   
            """
            info_names = [name.lower() for name in df.columns if not
                          (name.find('dep_')>=0  or name.find('res_')>= 0) ]
            # The depth spacing value  is the same like the res  
            depth_values = np.array([float(name.replace('dep_', '')) for 
                                     name in df.columns  if name.find('dep_')>=0  ])
    
            return info_names, depth_values , np.abs(depth_values.max() - 
                                                     depth_values.min())/ (len(depth_values)-1)
        
        def get_infohead_index(fnames, info_matrix ,  name, dtype ='float'): 
            """
            If oasis file is litthel messy about the head , better way to avoid 
            miscomputation 
            if to get the real index of colums names you are seeking so to
            get corresponding values
            
            Parameters
            -----------
                * fnames: list 
                        head of info names except the depth matrix
                        (dep_) and ressitivity matrix (res_)  
                * name: str 
                        name of head user are looking for 
                * info_matrix: ndarray 
                            ndarray array of oasis montag info .
                            start by ]station] to end normally  [doi max] 
            """
            for head in fnames : 
                if head.find(name)>=0 : 
                    name =head
                    index = fnames.index(name)
                    break 
            if dtype is not None : 
                indexmatrix = info_matrix[:, index].astype(dtype)
            return  indexmatrix
        
        
        if path is not None :  self.path = path 
        if self.path is not None : 
            self._validate_oasis_path()
      
        self.global_dico = {} # initiliase gloabal dico to takes all files with key as survey lines 
        for keys, vitems in self.read_dico.items(): 
            if keys == self.extension_file : 
                for files in self.oasis_data : 
                     self.global_dico [os.path.basename(files).\
                                  replace('{0}'.format(self.extension_file), '')] =\
                         vitems(files)
        
        self.info_depthvalues = { filename : get_depth_values(df_) for filename ,
                            df_ in self.global_dico.items()}
        # get the input filenames and set attributes 
    
        self.filenames = list(self.global_dico.keys())
    
            # set main attributes for selects key
           
        for oasnames in self.filenames :
            # convert dataframe to numpy : numpy very fast 
            tem_array = self.global_dico[oasnames].to_numpy() 
            # get the length (of head infos)
            len_info_names = len(self.info_depthvalues[oasnames][0]) 
            # get the length of depth infos 
            len_depth_offset = len(self.info_depthvalues[oasnames][1]) 
                # now set temporary attributes for each limes 
            self.__setattr__(oasnames, (tem_array[::, :len_info_names],  # info head array 
                                        tem_array[::, len_info_names:len_info_names +len_depth_offset], # depth array 
                                        tem_array[::, len_info_names +len_depth_offset:]) ) # get the resistivity array 
           
        # initiliase dictionary to hold  
        self.geos_profile_strike_angles ={}
        for gskeys,gsvalues in self.info_depthvalues.items(): 
            print('** ----- {0} : --|>{1:<37} :'.format('file',gskeys))
            print('** {0:<37} {1} {2} {3}'.format('depth spacing ','=', gsvalues[2], 'm' ))
            print('** {0:<37} {1} {2} {3}'.format('maximum depth','=', gsvalues[1].max(), 'm' ))
            # compute profile ange and geolectrike strike 
            easting = get_infohead_index(fnames=gsvalues[0], 
                                         info_matrix= getattr(self, gskeys)[0],
                                         name='east')
          
            northing= get_infohead_index(fnames=gsvalues[0], 
                                         info_matrix= getattr(self, gskeys)[0],
                                         name='north')
            try : 
                gstrike , profile_angle,_ = geostrike.\
                    compute_geoelectric_strike(easting=easting , 
                                                      northing = northing)
                print('** {0:<37} {1} {2} {3}'.format('profile angle ','=',
                                                      profile_angle, 'degrees E of N.' ))
                print('** {0:<37} {1} {2} {3}'.format('geoelectric strike','=',
                                                      gstrike, 'degrees E of N.' ))
            except : 
                warnings.warn('Trouble occurs when computing profile'\
                              ' angle and geo electric strike.')
                self._logging.debug('Trouble occurs when computing profile'\
                              ' angle and geo electric strike.')
            else : 
                self.geos_profile_strike_angles[gskeys]= (profile_angle, gstrike)
                 # pass when something wrong 
 
    def get_depth_surfaces(self, path =None , depth_values =None): 
        """
        get the depth surfaces for multi-lines and build the numpy 
        corresponding array  at that depth .

        :param depth_values: array of depth 
        :type  depth_values: float  or array_like  
                           
        :param path: full path to oasis outputfiles .
        :type path: str 
        
        """
        if path is not None : 
            self.path = path 
        if self.path is not None : self.read_oasis_files()
        
        if depth_values is not None : self.depth_values= depth_values    

        if self.depth_values is None : 
            raise CSex.pyCSAMTError_plot_geoinputargument(
                'NoneType could not be computed. Please provided'
                '  right values !')
        
        if isinstance(self.depth_values, (float, str)): 
            try : 
                self.depth_values = float(self.depth_values)
            except : 
                raise CSex.pyCSAMTError_float('Could not convert {0} to float'.\
                                              format(self.depth_values))
            else : self.depth_values =np.array([self.depth_values])
  
        elif isinstance(self.depth_values , (list, tuple, np.ndarray)): 
            self.depth_values = np.array(self.depth_values)
            try : 
                self.depth_values = self.depth_values.astype('float')
            except : raise CSex.pyCSAMTError_float('Could not convert values to float!')

        def get_single_surface_from_one_line(site_name , depth_value)  : 
            """
            Fonction to get single surface from depth value fron one line .
            
            :param depth value: depth value to image
            :type depth value: float 
                          
            :param site_name: str name of site
            :type site_name: str
            """
            # get attributes array : depth and resistivity
            #get depth index from info_depthvalues dictionnary (key  :
                    # values = info_names, depth_values, depth spacing )
            new_depth_value, index =  get_closest_value(
                values_range= self.info_depthvalues[site_name][1],
                              input_value= depth_value)
            # now get especially array or depth and resistivy at that value
            # from attribute name  value attr = infomatrix , depth matrix and res matrix 
            depth =  getattr(self, site_name)[1][:, index]
            res = getattr(self, site_name)[2][:, index]
            return new_depth_value , depth , res 
        
        # for loop to concat array 
        self.geosurface_dico ={}         # dict to keep all surface matrix  
        tem=[] 
        for dvalue in self.depth_values : 
            for names in self.filenames : 
                new_depth_value, sdepth , sres =\
                    get_single_surface_from_one_line(site_name=names , 
                                                     depth_value = dvalue)
                #build matrix array info depth res 
                globalmatrix =np.concatenate((getattr(self, names)[0], 
                                              sdepth.reshape(sdepth.shape[0], 1), 
                                              sres.reshape(sres.shape[0], 1)), 
                                             axis =1 )
                
                tem.append(globalmatrix)
                
            self.geosurface_dico[new_depth_value] =tem
            tem=[]
        for gkey, gvalue in self.geosurface_dico.items(): 
            # buid the array of each depth 
            self.geosurface_dico[gkey]=func.concat_array_from_list(gvalue)
            
        # create pandas dataframe 
        for hnames, hvalues  in self.info_depthvalues.items() : 
            # Generally for heads from oasis  are the same for all lines 
            # the let get one: 
            headnames = self.info_depthvalues[hnames][0] # break 
            break 
        
        for gnkey, gnvalue in self.geosurface_dico.items(): 
            columnames = headnames + [ 'dep_{0}'.format(int(gnkey))] + \
                [ 'rho_{0}'.format(int(gnkey))]
            self.geosurface_dico[gnkey] = pd.DataFrame(data=gnvalue,
                                                       columns = columnames)
            
    def write_file(self, path =None, depth_values =None, 
                   fileformat='.csv', **kwargs):
        """
        Write output files. Output files are `.xlsx` or `.csv` . 
     
        :param path:  full path to  `geodrill` ouput files
        :type path: str 
        
        :param depth_value:  depth values for imaging  
        :type depth_value: float or array_like
                   
        :param fileformat:`xlsx` or `csv` are actually the acceptable formats.
        :type fileformat: str
        
        """
        savepath = kwargs.pop('savepath', None)
        writeIndex= kwargs.pop('write_index', False)
        csvsetHeader =kwargs.pop('add_header',True)
        csvsep =kwargs.pop('csv_separateType',',')
        if savepath is not None : self.savepath =savepath 

        if path is not None : self.path = path 
        if depth_values is not None : self.depth_values =depth_values 
        if fileformat is not None : self.export_format = fileformat.lower() 
        if self.export_format is not None : 
            self.export_format= self.export_format.lower() # for consistency 
            
        if self.export_format.find('exc')>=0 or \
            self.export_format.find('xls')>=0 : 
                self.export_format = 'xlsx'
        elif self.export_format.find('csv')>=0 : 
            self.export_format ='csv'
        else : 
            mess=''.join(['-->Sorry! Geosurface actually does not support the {0} ',
                          'format provided. Only support `xlsx` or `csv` format.', 
                          ' Please provided the rigth format !'])
            self._logging.error(mess.format(self.export_format))
            warnings.warn(mess.format(self.export_format))
            raise CSex.pyCSAMTError_file_handling(
                'Format provided = {0} is wrong ! Could output only `csv` or `xlsx`.'.\
                    format(self.export_format))
            
        
        if self.depth_values is None : 
            warnings.warn('Need to specify the value of depth for imaging !')
            self._logging.warn('Need to specify the value of depth for imaging !')
            raise CSex.pyCSAMTError_inputarguments(
                '! Need to specify the depth value for imaging.')
        if self.path is None : 
            warnings.warn('Need to provide the rigth path `geodrill` model output files.')
            self._logging.warning('Need to provide the rigth path `geodrill` model output files.')
            raise CSex.pyCSAMTError_inputarguments(
                'No path detected ! Please provide a right path to `geodrill`\
                    oasis montaj outputfiles ')
  
        self.get_depth_surfaces() # call methods to create 
        
        # get values from dictionnary and write file 
        #make output filename 
        filenames =[]
        filename =''.join([file.replace('_cor_oas', '') 
                                for file in self.filenames])
        if self.export_format =='xlsx':
            with pd.ExcelWriter(''.join([filename+ 
                                         '_gs{0}.'.format(datetime.datetime.now().month ),
                                         self.export_format])) as writer :
                for file , df_ in self.geosurface_dico.items(): 
                        df_.to_excel(writer,
                                     sheet_name='dep_{0}'.format(file), 
                                     index =writeIndex)
            filenames.append(''.join(
                            [filename+ 
                            '_gs{0}.'.format(datetime.datetime.now().month ),
                                      self.export_format]))
        if self.export_format =='csv': 
            
            for file , df_ in self.geosurface_dico.items():
                df_.to_csv(''.join([filename + str(file),
                                    '_gs{0}.'.format(datetime.datetime.now().month),
                                    self.export_format]),
                              header=csvsetHeader,
                              index =writeIndex,sep=csvsep)
                
                filenames.append(
                    ''.join([filename + str(file) ,
                             '_gs{0}.'.format(datetime.datetime.now().month ), 
                             self.export_format]))
                
        # export to savepath 
        if self.savepath is None : # create a folder in your current work directory
            try :
                self.savepath = os.path.join(os.getcwd(), '_outputGS_')
                if not os.path.isdir(self.savepath):
                    os.mkdir(self.savepath)#  mode =0o666)
            except : 
                warnings.warn("It seems the path already exists !")
                
        #savefile inot it savepath 
        if self.savepath is not None : 
            import shutil 
            try :
                for nfile in filenames : 
                   shutil.move(os.path.join(
                       os.getcwd(),  nfile),           
                       os.path.join(self.savepath , nfile))

            except : 
                warnings.warn("It seems files already exist !")
                
        # fmt =''.join(['{0}'.format(ifile) for ifile in filenames])
        print('---> Geosurfaces outputfiles :{0} : have been successfully '\
                  'written to  <{1}>.'.format(','.join(
                      ['{0}'.format(ifile) for ifile in filenames]),
                                              self.savepath))
        print('-'*77)
           
class geostrike : 
    """
    Class to deal with computation with profile angle and geo_electrical strike.
    Compute Profile angle and strike angle 
     Need to import scipy.stats as one module . Sometimes import scipy 
    differently  with stats may not work . 
    either `import scipy.stats` rather than `import scipy as sp` to use : 
        `sp.stats.linregress` .
        
    """

    @staticmethod 
    def compute_profile_angle (easting=None, northing=None): 
        """
        Essentially dedicated to compute geoprofile angle. 
        
        Parameters 
        -----------
            * easting : array_like 
                    easting coordiantes values 
            * northing : array_like 
                    northing coordinates values
                
        Returns 
        ---------
            float
                 profile_angle 
            float 
                geo_electric_strike 
            str 
                message of return 
        """
        _logger.info('Computing  profile angle from Easting and Nothing coordinates.')
        if easting is None or northing is None : 
            raise CSex.pyCSAMTError_inputarguments('NoneType can not be computed !')
            
            # use the one with the lower standard deviation
        try :
            easting = easting.astype('float')
            northing = northing.astype('float')
        except : 
            raise CSex.pyCSAMTError_float('Could not convert input argument to float!')
        
        if stats_import is True : 
            profile1 = spSTAT.linregress(easting, northing)
    
            profile2 =spSTAT.linregress(northing, easting)
        else :
            warnings.warn('Could not find scipy.stats, cannot use method linearRegression '
                  'check installation you can get scipy from scipy.org.')
            _logger.warning('Could not find scipy.stats, cannot use method lineaRegress '
                    'check installation you can get scipy from scipy.org.')
            raise ImportError('Could not find scipy.stats, cannot use method lineaRegress '
                    'check installation you can get scipy from scipy.org.')
            
        profile_line = profile1[:2]
        # if the profile is rather E=E(N), the parameters have to converted  into N=N(E) form:
        
        if profile2[4] < profile1[4]:
            profile_line = (1. / profile2[0], -profile2[1] / profile2[0])

        # if self.profile_angle is None:
        profile_angle = (90 - (np.arctan(profile_line[0]) * 180 / np.pi)) % 180

        # otherwise: # have 90 degree ambiguity in strike determination# choose strike which offers larger angle with profile
        # if profile azimuth is in [0,90].
        
        return np.around(profile_angle,2) , 'Profile angle is {0:+.2f} degrees E of N'.format(profile_angle)
         
        
      
    @staticmethod 
    def compute_geoelectric_strike (profile_angle = None , easting =None, 
                                    northing=None, **kws):
        """
        Compute geoelectric strike
        
        Parameters
        -------------
            *  profile_angle : float 
                      If not provided , will comput with 
                      easting and northing coordinates 
            * easting : array_like 
                       Easting coordiantes values 
            * northing : array_like 
                       Northing coordinates values 
            * geo_strike : float 
                       strike value , if provided, will 
                       recomputed geo_electric strike .
        Returns 
        --------
            float
                 profile_angle in degree E of N 
            float 
                geo_electric_strike in degrees E of N
            str 
                message of return 
        
        """
        gstrike =kws.pop('geo_strike', None)
        
        if profile_angle is None and  easting is None and northing is None : 
            mess =''.join(['None type detected for profile angle , easting and'\
                           ' northing. NoneType can not be computed', 
                           ' will check whether geostrike is provided .'])
            _logger.debug(mess)
            if gstrike is None :
                
                _logger.warning('NoneType found.Could not compute geo-electrike strike!')
                raise CSex.pyCSAMTError_inputarguments(
                    'NoneType found. Could not compute geo-electrike strike!')
        
        if profile_angle is None : 
            if easting is not None and northing is not None : 
                profile_angle ,_ = geostrike.compute_profile_angle(easting, northing)
        
        if gstrike is None : 
            if 0<= profile_angle < 90 :
                geo_electric_strike  = profile_angle + 90  
            elif 90<=profile_angle < 180 :
                geo_electric_strike = profile_angle -90
            elif 180 <= profile_angle <270 :
                geo_electric_strike = - profile_angle +90 
            else :
                geo_electric_strike  = - profile_angle -90 
            
            geo_electric_strike  %= 180   
    
        if gstrike is not None : # recomputed geo_electrike strike 
            if 0 <= profile_angle < 90:
                if np.abs(profile_angle - gstrike) < 45:
                    geo_electric_strike  = gstrike+ 90
     
            elif 90 <= profile_angle < 135:
                if profile_angle - gstrike < 45:
                    geo_electric_strike = gstrike - 90
            else:
                if profile_angle - gstrike >= 135:
                   geo_electric_strike = gstrike+ 90
            geo_electric_strike %=  180         # keep value of geoelectrike strike less than 180 degree
            
        geo_electric_strike =np.floor(geo_electric_strike)
        
        return  geo_electric_strike, profile_angle , 'Profile angle is {0:+.2f} degrees E of N'.format(geo_electric_strike)
        
      
class Drill(object):
    """
    This class is focus on well logs . How to generate well Log for Oasis:
        
    Arguments
    -----------
        **well_filename** : string ,
                     The  well filename. 02 options is set : 
                     1rst option is to build well data manually and the program will  
                     generate a report.  2nd option is to send to
                     the program a typical file type to be parsed . the programm parses
                     only the typical well datafile. If None ,  the program will 
                    redirect to build mannually option . 
                    
        **build_manually_welldata** : bool  
                     option to build manually the well data . set to True 
                      if you want to build manually a well data .
                     *default* is False

    ====================  ==========  =========================================    
    Key Words/Attributes  Type          Description    
    ====================  ==========  =========================================   
    utm_zone                str         utm WGS84 zone. should be N or S.   
                                        *default* is 49N .
    compute_azimuth         bool        if no azimuth is provided. 
                                        set to True to letprogram to compute
                                        azimuth .*Default* is False.
    Drill_dip               float       The dip of drill hole.*default* is 90
    Drill_buttom            float       The average bottom of drill , 
                                        can be filled during the well
                                        buiding . *default* is  None
    mask                    int         the mask of DrillHole(DH) data. 
                                        *Default * is 1.
    ====================  ==========  =========================================
    
    ==================  =======================================================
    Methods                   Description
    ==================  =======================================================
    _collar             build _collar data *return*  collar log dataframe 
                        format
     dhGeology          build DH log geology *return* geology log dataframe.        
    dhSample            build DH Geochemistry-Strutural sample, *return* Sample
                        log dataframe    
    dhSurveyElevAz      build DH Elevation & Azimuth logs.*return * Elevation
                        & Azimuth dataframes
    writeDHDATA          output log :* return *  the right log to output for
                        Oasis Montaj 
    ==================  =======================================================
        
    :Example: 
        
        >>> from geodrill.geoCore.geodrill import Drill 
        >>> parser_file ='nbleDH.csv'
        >>> drill_obj=Drill(well_filename=os.path.join(os.environ['pyCSAMT'],
        ...                                           parser_file),
        ...      build_manually_welldata=False)
        >>>  scollar=drill._collar(DH_Top=None)
        >>> sgeo=drill.dhGeology()
        >>> ssam=drill.dhSample()
        >>> selevaz=drill.dhSurveyElevAz( add_elevation=None, 
        ...                             add_azimuth=None)
        >>> swrite=drill.writeDHData(data2write ="*",
                                 savepath =None)
    """
    
    def __init__(self, well_filename=None , build_manually_welldata=False, **kwargs):
        
        self._logging=csamtpylog.get_csamtpy_logger(self.__class__.__name__)
        
        self.wfilename=well_filename
        self.buildmanuel=build_manually_welldata
        
        self.mask=kwargs.pop("mask",1)
        self.utm_zone=kwargs.pop("utm_zone","49N")
        self.compute_azimuth=kwargs.pop("compute_azimuth",False)
        self.dip =kwargs.pop("Drill_dip",90)
        self.buttom=kwargs.pop("Drill_buttom", None)
        # self.elevfile=kwargs.pop("elevation_file", None)
        
        self.easts=None
        self.norths= None
        self.wellnames= None
        self._f=None 
        
        #populate attribute later 
        self.wdico={"DH_Hole" :None, 
                    "DH_East":None, 
                    "DH_North":None, 
                    "Mask": None,
                    "DH_RH":None, 
                    'DH_From':None , 
                    "DH_To": None , 
                    "Rock": None , 
                    "DH_Azimuth":None , 
                    'DH_Top':None, 
                    'DH_Bottom':None,
                    'DH_PlanDepth':None, 
                    'DH_Decr':None, 
                    'Sample':None,
                    'DH_Dip': None, 
                    'Elevation':None,
                    'DH_RL':None,
                    }        
        
        
        if self.buildmanuel is True and self.wfilename is None :
            
            self.daTA=func.build_wellData (add_azimuth=False, 
                                            utm_zone=self.utm_zone,
                                            report_path=self.compute_azimuth)
            self.wdata=self.daTA[1]
            
            self.wdico["DH_East"]   =   self.wadata[:,1]
            self.wdico["DH_North"]  =   self.wdata[:,2]
            self.wdico["DH_Hole"]   =   self.wadata[:,0]
            self.wdico['DH_Dip']    =   self.wadata[:,4]
            self.wdico['DH_Bottom'] =   self.wdata[:,3]
            self.wdico['DH_Decr'] =   self.wdata[:,7]
            self.wdico['DH_PlanDepth'] =   self.wdata[:,6]
            self.wdico["DH_Azimuth"] =   self.wdata[:,5]
            
            self._f=0


        elif  self.wfilename is not None :
            
            self.daTA=func.parse_wellData(filename=self.wfilename,
                                          include_azimuth=False,
                                          utm_zone=self.utm_zone)
            self.wdata=self.daTA[1]
            self.wdico.__setitem__("DH_East", self.wdata[:,1])
            self.wdico.__setitem__("DH_North", self.wdata[:,2])
            self.wdico.__setitem__("DH_Hole", self.wdata[:,0])
            self.wdico.__setitem__('DH_Dip', self.wdata[:,3])
            self.wdico.__setitem__('DH_PlanDepth', self.wdata[:,8])
            self.wdico.__setitem__("DH_Azimuth", self.wdata[:,5])
            self.wdico.__setitem__('DH_Decr', self.wdata[:,9])
            self.wdico.__setitem__('DH_Bottom', self.wdata[:,7])
            
            self._f=1
            

        #set Mask and set dr_rh
        self.mask=np.full((self.wdata.shape[0]),self.mask,dtype='<U12')
        # print(self.mask.shape)
        self.wdico.__setitem__("Mask", self.mask)
        self.dh_rh=np.zeros((self.wdata.shape[0]))
        self.wdico.__setitem__("DH_RH", self.dh_rh)  

        for keys in kwargs.keys():
            self.__setattr__(keys, kwargs[keys])
            
 
  
    def _collar(self, DH_Top=None,add_elevation =None ):
        """
        Method to build Collar Data 
        
        Parameters 
        ----------
            * DH_Top  : np.ndarray ,
                    it's the Top of data for each Hole Name. ndaray (number of DH , 1) 
                    *Default* is None.
        Returns
        -------
            pd.DataFrme 
                collar Drillhole log
        """

        if DH_Top is None :
            DH_Top=np.zeros((self.wdata.shape[0]))
        elif type(DH_Top) is float or type(DH_Top) is int :
            DH_Top=np.full((self.wdata.shape[0]),DH_Top,dtype='<U12')
            
        elif DH_Top is not None :
            if type(DH_Top)==list:
                DH_Top=np.array(DH_Top)
                
            assert DH_Top.shape[0]==self.wdata.shape[0],'the input DH_Top '\
                'shape doesnt match. The convenience  shape is %d.'%self.wdata.shape[0]
        
        # print(DH_Top)
        self.wdico.__setitem__('DH_Top',DH_Top)
        
        if self._f == 0 :
            if add_elevation is None :
                #No topography is added , set to 0 
                add_elevation=np.full((len(self.wdico['DH_EAST']),1),0,dtype='<U12')
            elif add_elevation is not None :
                if type(add_elevation ) is list :
                    add_elevation =np.array(add_elevation)
                assert add_elevation.shape[0]==self.wdico['DH_EAST'].shape[0],"INDEXERROR:"\
                    " The the current dimention of Elevation data is {0}.It's must be"\
                        " the size {1}.".format(add_elevation.shape[0],self.wdico['DH_EAST'].shape[0])
            
            self.wdico.__setitem__("Elevation", add_elevation)
                    
        elif self._f == 1 :
            
            if add_elevation is not None:
                
                if type(add_elevation ) is list :
                    add_elevation =np.array(add_elevation)
                try :
                    np.concat((add_elevation,self.wdico['DH_East']))
                except Exception : 
                    mess =''.join(['SIZEERROR! Try to set the elevation dimentional as ', 
                                    'same size like the collar data'])
                    self._logging.error(mess)
                    warnings.warn(mess)
                    
            elif add_elevation is None :
                add_elevation=self.daTA [1][:,4]
        
            self.wdico.__setitem__("Elevation", add_elevation)
        
        collarKeys=["DH_Hole",	"DH_East",	"DH_North",	"DH_RH",
                    "DH_Dip", "Elevation", "DH_Azimuth","DH_Top", "DH_Bottom",
                    "DH_PlanDepth",	"DH_Decr",	"Mask"] 
        
        # print(self.wdico)
        collar=self.wdico[collarKeys[0]]
        collar=collar.reshape((collar.shape[0],1))
        for ss, collk in enumerate(collarKeys[1:]):  
            # print(collk)
            for key , value in self.wdico.items():
                if key == collk :
                    value=value.reshape((value.shape[0],1))
                    collar=np.concatenate((collar,value), axis=1)
        
        
        self.coLLAR=pd.DataFrame(data=collar, columns=collarKeys)

        return self.coLLAR
    
    
    def dhGeology (self, dh_geomask=None):
        """
        Method to build geology drillhole log. The name of input rock must feell exaction accordinag to
        a convention AGSO file . If not sure for the name of rock and Description and label .
        you may consult the geocode folder before building the well_filename. If the entirely
        rock name is given , program will search on the AGSO file the corresponding Label and code . 
        If the rock name is  founc then it will take its CODE else it will generate exception. 
 
        Parameters
        ----------
            * dh_geomask : np.ndarray, optional
                        geology mask. send mask value can take exactly
                        the np.ndarray(num_of_geology set ,). The better way 
                        to set geology maskis to fill on the wellfilename.
                        if not , programm will take the general mask value. 
                        The *default* is None.

        Returns
        -------
            pd.DataFrame 
                geology drillhole log.
        """
        
        
        geolKeys=["DH_Hole","DH_From",	"DH_To","Rock",	"Sample",
                  "East",	"DH_North",	"DH_RH",	"Mask"]
        
        wgeo=self.daTA[2]
        # print(wgeo)
        
        self.wdico.__setitem__('DH_From', wgeo[:,1])
        self.wdico.__setitem__('DH_To', wgeo[:,2])
        self.wdico.__setitem__("Rock",wgeo[:,3])
        dhgeopseudosamp=np.zeros((wgeo.shape[0]))
        
        ###### FIND AGSO MODULE #######
        #Try to check the name of rocks and their acronym
        geoelm=Agso._agso_on_dict_(set_agsoDataFrame=False, return_orientation="SERIES")
            # #extract elem with their acronym 
        geolemDico_AGSO={key:value for key , value in \
                         zip (geoelm["CODE"],geoelm['__DESCRIPTION'])}
        # elemgeo_AGSO=sorted(geolemDico.items())
    
        
        for ii, elm in enumerate (self.wdico['Rock']):
            if elm.upper() in geolemDico_AGSO.keys():
                pass 
            elif elm.upper() not in geolemDico_AGSO.keys():
                if elm.lower() in geolemDico_AGSO.values():
                    for key, values in geolemDico_AGSO.items():
                        if elm.lower() == values :
                            self.wdico['Rock'][ii]=key
                else  :
                    mess=''.join(['The Geological Name ({0}) given in is wrong'.format(elm),
                                'Please provide a right name the right Name.', 
                                'Please consult the AGSO file in _geocodes folder', 
                                'without changing anything.'])
                    self._logging.warn(mess)
                    warnings.warn(mess)

        ######END AGS0 ########
        
        self.dh_geoleast=np.zeros((wgeo.shape[0]))
        self.dh_geol_norths=np.zeros((wgeo.shape[0]))
        
        for ss , value in enumerate(self.dh_geoleast):
            for indix, val in enumerate(self.wdico["DH_East"]):
                if wgeo[:,0][ss] in self.wdico["DH_Hole"]:
                    value=val
                    self.dh_geoleast[ss] =value
                    self.dh_geol_norths[ss]=self.wdico["DH_North"][indix]
                    
        dhgeopseudosamp=np.zeros((wgeo.shape[0]))

        if dh_geomask == None :
            dh_geomask =self.mask[0]
        maskgeo= np.full((wgeo.shape[0]),dh_geomask,dtype='<U12')
        dhrhgeo=np.array([ -1* np.float(ii) for ii in self.wdico['DH_From']])
        dhGeol=np.concatenate((wgeo[:,0].reshape(wgeo[:,0].shape[0],1),
                              self.wdico['DH_From'].reshape((self.wdico['DH_From'].shape[0],1)),
                              self.wdico['DH_To'].reshape((self.wdico['DH_To'].shape[0],1)),
                              self.wdico['Rock'].reshape((self.wdico['Rock'].shape[0],1)),
                              dhgeopseudosamp.reshape((dhgeopseudosamp.shape[0],1)),
                              self.dh_geoleast.reshape((self.dh_geoleast.shape[0],1)),
                              self.dh_geol_norths.reshape((self.dh_geol_norths.shape[0],1)),
                              dhrhgeo.reshape((dhrhgeo.shape[0],1)),
                              maskgeo.reshape((maskgeo.shape[0],1))),axis=1)
        self.geoDHDATA=pd.DataFrame(data=dhGeol, columns=geolKeys)
        
        return self.geoDHDATA
    
           
    def dhSample (self,path_to_agso_codefile=None, dh_sampmask=None):
        """
        Method to build Sample log. This method focuses on the sample obtained 
        during the DH trip.it may georeferenced as the well_filename needed. 
        A main thing is to set the AGSO_STCODES file. AGSO_STCODES is the 
        conventional code of structurals sample. If you have an own AGSO_STCODES ,
        you may provide the path * kwargs=path_to_ags_codefile * . 
        the program will read and generate logs according to the  DESCRIPTION 
        and STCODES figured. if None, the program will take it STCODES  and set
        the samplelogs. When you set the Sample code aor sample name , 
        make sur that the name match the same name on STCODES. If not ,
        program will raises an error. 

        Parameters
        ----------
            * path_to_agso_codefile : str, optional
                                path to conventional
                                AGSO_STRUCTURAL CODES.
                                The *default* is None.
                                
            * dh_sampmask : np.ndarray, optional
                            Structural mask. The default is None.

        Returns
        -------
            pd.DataFrame 
                Sample DH log.
        """
        
        sampKeys=["DH_Hole","DH_From",	"DH_To","Rock",	"Sample",
                  "East",	"DH_North",	"DH_RH",	"Mask"]
        
        wsamp=self.daTA[3]
        # print(wgeo)
        
        self.wdico.__setitem__('DH_From', wsamp[:,1])
        self.wdico.__setitem__('DH_To', wsamp[:,2])
        self.wdico.__setitem__("Sample",wsamp[:,3])
        dhsampseudorock=np.zeros((wsamp.shape[0]))
        
        ###### FIND AGSO MODULE (AGSO_STCODES) #######
        #Try to check the name of sample and their acronym
        
        if path_to_agso_codefile is None:
            path_to_agso_codefile= os.path.join(os.environ ['pyCSAMT'],
                                                'geodrill','_geocodes' )
       
            agsofilename=[file for file in os.listdir(path_to_agso_codefile) 
                          if file =='AGSO_STCODES.csv' ]
            if agsofilename is not None :
                sampelm=Agso._agso_on_dict_(set_agsoDataFrame=True, 
                                            return_orientation="series", 
                            agso_codefile=os.path.join(
                                path_to_agso_codefile,
                                agsofilename[0]))
            elif agsofilename is None :

                self._logging.warning(
                'None AGSO_STCODES.csv file is found. Please provide the right path '
                                      )
                warnings.warn('None AGSO_STCODES.csv file is found. Please provide the right path ')
                                      
        elif path_to_agso_codefile is not None : 
            #os.chdir(os.path.dirname(path_to_agso_codefile))
            sampelm=Agso._agso_on_dict_(set_agsoDataFrame=True, 
                                        return_orientation="series", 
                            agso_codefile=path_to_agso_codefile)
        
            # #extrcat elem with their acronym 
        sampelmDico_AGSO={key:value for key , value in \
                         zip (sampelm["CODE"],sampelm['__DESCRIPTION'])}
        # elemgeo_AGSO=sorted(geolemDico.items())

        
        for ii, elm in enumerate (self.wdico['Sample']):
            if elm.lower() in sampelmDico_AGSO.keys():
                pass 
            elif elm.lower() not in sampelmDico_AGSO.keys():
                if elm in sampelmDico_AGSO.values():
                    for key, values in sampelmDico_AGSO.items():
                        if elm  == values :
                            self.wdico['Sample'][ii]=key
                else  :
                    mess=''.join(['The Sample Name({0}) given in is wrong'.format(elm),
                                'Please provide a right name the right Name.', 
                                'Please consult the AGSO_STCODES.csvfile in _geocodes folder', 
                                'without changing anything.'])
                    self._logging.warn(mess)
                    warnings.warn(mess)

        ######END AGS0_STCODES ########
        
        dh_sampeast=np.zeros((wsamp.shape[0]))
        dh_sampnorths=np.zeros((wsamp.shape[0]))
        
        for ss , value in enumerate(dh_sampeast):
            for indix, val in enumerate(self.wdico["DH_East"]):
                if wsamp[:,0][ss] in self.wdico["DH_Hole"]:
                    value=val
                    dh_sampeast[ss] =value
                    dh_sampnorths[ss]=self.wdico["DH_North"][indix]
                    
        dhsampseudorock=np.zeros((wsamp.shape[0]))

        if dh_sampmask == None :
            dh_sampmask =self.mask[0]
        masksamp= np.full((wsamp.shape[0]),dh_sampmask,dtype='<U12')
        dhrhsamp=np.array([ -1* np.float(ii) for ii in self.wdico['DH_From']])
        dhSample=np.concatenate((wsamp[:,0].reshape(wsamp[:,0].shape[0],1),
                              self.wdico['DH_From'].reshape((self.wdico['DH_From'].shape[0],1)),
                              self.wdico['DH_To'].reshape((self.wdico['DH_To'].shape[0],1)),
                              dhsampseudorock.reshape((dhsampseudorock.shape[0],1)),
                              self.wdico['Sample'].reshape((self.wdico['Sample'].shape[0],1)),
                              dh_sampeast.reshape((dh_sampeast.shape[0],1)),
                              dh_sampnorths.reshape((dh_sampnorths.shape[0],1)),
                              dhrhsamp.reshape((dhrhsamp.shape[0],1)),
                              masksamp.reshape((masksamp.shape[0],1))),axis=1)
        self.sampleDHDATA=pd.DataFrame(data=dhSample, columns=sampKeys)
        
        return self.sampleDHDATA
    
    def dhSurveyElevAz(self, add_elevation=None, add_azimuth=None, **kwargs):
        """
        Method to build Elevation & Azimuth DH logs. if add_elevation and . 
        add_azimuth are set . The programm will ignore the computated azimuth,
        and it will replace to the new azimuth   provided . all elevation will 
        be ignore and set by the new elevation . *kwargs arguments 
        {add_elevation , add-azimuth }  must match the same size like the 
        number of Drillholes . Each one must be on ndarray(num_of_holes, 1). 
        
        Parameters
        ----------
            * add_elevation : np.nadarray , optional
                    elevation data (num_of_holes, 1) 
                    The *default* is None.
                    
            * add_azimuth : np.ndarray , optional
                    azimuth data (num_of_holes,1). 
                    The *default* is None.
                    
            * DH_RL :np.float or np.ndarray(num_of_hole,1),
                    if not provided , it's set to 0. means No topography is added'.
                
        Returns
        -------
            pd.Dataframe 
                Elevation DH log .
            pd.DataFrame 
                Azimuth DH log.
        """
        
        
        dh_rl=kwargs.pop("DH_RL",None)
        
        sizep=self.wdico['DH_East'].shape[0]
        if self._f == 0 :
            if add_elevation is None :
                #No topography is added , set to 0 
                add_elevation=np.full((len(self.wdico['DH_EAST']),1),0,dtype='<U12')
            elif add_elevation is not None :
                if type(add_elevation ) is list :
                    add_elevation =np.array(add_elevation)
                assert add_elevation.shape[0]==self.wdico['DH_EAST'].shape[0],"INDEXERROR:"\
                    " The the current dimention of Elevation data is {0}.It's must be"\
                        " the size {1}.".format(add_elevation.shape[0],self.wdico['DH_EAST'].shape[0])
            
            self.wdico.__setitem__("Elevation", add_elevation)
                    
        elif self._f == 1 :
            
            if add_elevation is not None:
                
                if type(add_elevation ) is list :
                    add_elevation =np.array(add_elevation)
                try :
                    np.concat((add_elevation,self.wdico['DH_East']))
                except :
                    mess= ''.join(['SIZEERROR! Try to set the elevation dimentional. ', 
                                   'same like the collar data '])
                    self._logging.error(mess)
                    warnings.warn(mess)
            elif add_elevation is None :
                add_elevation=self.daTA [1][:,4]
        
            self.wdico.__setitem__("Elevation", add_elevation)
            
        #set DH_RL
        if dh_rl is not None : 
            if type (dh_rl) is list : 
                dh_rl=np.array (dh_rl)
            assert dh_rl.shape[0]==self.data.shape[0]," DH_RL data size is out"\
                " of the range.Must be {0}".format(self.data.shape[0])
                
            self.wdico.__setitem__("DH_RL",dh_rl)
            
        elif dh_rl is None :
            #if None set DH_RL to None :
            self.wdico.__setitem__("DH_RL",np.full((self.daTA[1].shape[0]),0,dtype='<U12'))
        
        #set azimuth 
        if add_azimuth  is not None : 
            if type(add_azimuth) ==list : 
                add_azimuth=np.array(add_azimuth)
            assert add_azimuth.shape[0]==self.data.shape[0]," Azimuth data size is out"\
                " of the range.Must be {0}".format(self.data.shape[0])
                
            self.wdico.__setitem__("DH_Azimuth",add_azimuth) 
            
        elif add_azimuth is None : 
            pass 
                
        elevazKeys=['DH_Hole', 'Depth','DH_East',
                    'DH_North','Elevation','DH_RL','DH_Dip']
        
        self.wdico.__setitem__("DH_RL",np.full((self.daTA[1].shape[0]),0,dtype='<U12'))
        # add Hole and Depth 
        
        surveyELEV =np.concatenate((self.wdico['DH_Hole'].reshape((self.wdico['DH_Hole'].shape[0],1)),
                                    self.wdico["DH_Bottom"].reshape((self.wdico["DH_Bottom"].shape[0],1))),
                                       axis=1)
        surveyAZIM=np.concatenate((self.wdico['DH_Hole'].reshape((self.wdico['DH_Hole'].shape[0],1)),
                                    self.wdico["DH_Bottom"].reshape((self.wdico["DH_Bottom"].shape[0],1))),
                                      axis=1)
        
        for ss , elm in enumerate (elevazKeys[2:]):
            for key, values in self.wdico.items():
                if elm==key :
                    values=values.reshape((values.shape[0],1))
                    if elm =='DH_RL'or elm=='DH_Dip':
                        # print(values)
                        surveyAZIM=np.concatenate((surveyAZIM,values),axis=1)
                    elif  elm=='Elevation':
                        surveyELEV =np.concatenate((surveyELEV,values),axis=1)
                    else:
                        surveyAZIM=np.concatenate((surveyAZIM,values),axis=1)
                        if ss < elevazKeys.index('Elevation')-1: 
                            surveyELEV =np.concatenate((surveyELEV,values),axis=1)
                            
        
        self.surveyDHELEV=pd.DataFrame(data=surveyELEV, columns=elevazKeys[:5])
        # pop the elevation elm on the list 
        [elevazKeys.pop(ii) for ii, elm in enumerate(elevazKeys) if elm=='Elevation']
        
        self.surveyDHAZIM=pd.DataFrame(data=surveyAZIM, columns=elevazKeys)
        
        return (self.surveyDHELEV, self.surveyDHAZIM)
        
                    
    def writeDHData (self, data2write =None ,**kwargs):
        """ 
        Method to write allDH logs. It depends to the users to sort which data 
        want to export and which format. the program support only two format 
        (.xlsx and .csv) if one is set , it will ouptput the convenience format.
        Users can give a list of  the name of log he want to export.
        Program is dynamic and flexible. It tolerates quite symbols number to
         extract data logs. 
        
        Parameters
        ----------
            * data2write : str or list , optional
                        the search key. The default is None.
            
            * datafn :str
                    savepath to exported file 
                    *Default* is current work directory.
                    
            * write_index_on_sheet : bool, 
                    choice to write the sheet with pandas.Dataframe index. 
                    
            * writeType : str , 
                    file type . its may *.csv or *.xlsx .
                    *Default* is *.xlsx
                    
            * add_header : bool, 
                    add head on exported sheet. set False to mask heads. 
                    *Default* is True. 
                    
            * csv_separateType : str , 
                    Indicated for csv exported files , 
                    the type of comma delimited . defaut is ','.
        """
    
        writepath =kwargs.pop("savepath",None )
        writeIndex=kwargs.pop('write_index_on_sheet',False)
        writeType =kwargs.pop('writeType', 'xlsx')
        csvencoding =kwargs.pop('encoding','utf-8')
        csvsetHeader =kwargs.pop('add_header',True)
        csvsep =kwargs.pop('csv_separateType',',')
        
        
        wDATA ={"collar": self._collar,
                 "geology": self.dhGeology,
                 'sample':self.dhSample,
                 'elevazim':self.dhSurveyElevAz}
        
        _all=['5',"all","__all__",'CollGeoSampElevAz','CGSAZ','cgsaz',
              ['Collar','Geology','Sample','Elevation','Azimuth'],
              'colgeosamelevaz','alldata','*']
        
        df_collar=wDATA['collar']()
        df_geology=wDATA['geology']()
        df_sample=wDATA['sample']()
        df_elevation,df_azimuth=wDATA['elevazim']()
        
        # for df_ in  [df_collar, df_geology, df_sample, df_elevation,df_azimuth]: 
        #         df_.set_index(setIndex) # this is unnecessary 
        
        _dHDico ={'collar': [['1','c'], df_collar],
                 'geology':[['2','g'],df_geology],
                 'sample': [['3','s'],df_sample],
                 'survey_elevation':[['4','elev', 'topo','topography','e'], df_elevation],
                 'survey_azimuth': [['5','-1','azim','a'],df_azimuth]}
        
            
        if data2write is None or data2write in _all :  # write all 

            with pd.ExcelWriter(''.join([self.daTA[0][:-1],'.xlsx'])) as writer :
                for keys, df_ in _dHDico.items():
                    df_[1].to_excel(writer,sheet_name=keys, index =writeIndex)

                                
        elif data2write is not None :
            
            if type(data2write) is not list:
                data2write=str(data2write)
                
                try :

                    if writeType in ['xlsx','.xlsx', 'excell','Excell','excel','Excel','*.xlsx']:
                        for keys, df in _dHDico.items():
                            if data2write ==keys or data2write.lower() in keys or  data2write in df[0]:
                              df[1].to_excel('.'.join([self.daTA[0][:-1],'xlsx']),sheet_name=keys,index =writeIndex)  

                        
                    elif writeType in ['csv','.csv', 'comma delimited','*.csv',
                                       'comma-separated-value','comma seperated value',
                                       'comsepval']:
                        # print('passed')
                        for keys, df_ in _dHDico.items():
                            if data2write == keys or data2write.lower() in keys or data2write in df_[0]:
                              df_[1].to_csv(''.join([self.daTA[0][:-1],'.csv']), header=csvsetHeader,
                                    index =writeIndex,sep=csvsep, encoding=csvencoding)  

                except Exception as error :
                    self._logging.error ('The type you provide as WriteType argument is wrong.'
                                          ' Support only *.xlsx and *.csv format',error)
                    warnings.warn ('Argument writeType support only [xlsx or csv] format.'
                                    ' Must change your *.{0} format'.format(writeType))

                
            elif type(data2write) is list :
                data2write=[str(elm) for elm in data2write] # check the string format
                with pd.ExcelWriter(''.join([self.daTA[0][:-1],'.xlsx'])) as writer :
                    for ii, df in enumerate (data2write):
                        for keys, df__ in _dHDico.items():
                            if df.lower() in keys or df in df__[0] : 
                                df__[1].to_excel(writer,sheet_name=keys, index =writeIndex)
            else :
                self._logging.error ('The key you provide  as agrument of data2write is wrong. '
                                     'the data2write argument should be either [collar, geology,'
                                         ' sample, elevation, azimuth] or all (*). ')
                warnings.warn ('Wrong format of input data2write ! Argument dataType is str,'
                               ' or list of string element choosen among [collar, geology,'
                                   'sample, elevation, azimuth] or all (*), not {0}'.format(data2write))
        
         # export to savepath 
        if writepath is None : # create a folder in your current work directory
            try :
                writepath = os.path.join(os.getcwd(), '_outputDH_')
                if not os.path.isdir(writepath):
                    os.mkdir(writepath)#  mode =0o666)
            except : 
                warnings.warn("It seems the path already exists !")
        
        
        if writepath is not None  :
            import shutil
            
            if writeType in ['csv','.csv', 'comma delimited',
                             'comma-separated-value','comma sperated value',
                                       'comsepval']:
                shutil.move ( os.path.join(os.getcwd(),
                                           ''.join([self.daTA[0][:-1],'csv'])),writepath)
                print('---> Borehole output <{0}> has been written to {1}.'.\
                      format(os.path.basename(
                    ''.join([self.daTA[0][:-1],'.csv'])), writepath))
                
            elif writeType in ['xlsx','.xlsx', 'excell','Excell','excel','Excel']:
                
                shutil.move ('.'.join([self.daTA[0][:-1],'xlsx']),writepath)
                
                print('---> Borehole output <{0}> has been written to {1}.'.\
                      format(os.path.basename(
                      ''.join([self.daTA[0][:-1],'.xlsx'])), writepath))
                
#------Usefull functions --------------------------
def get_closest_value (values_range, input_value): 
    """
    Fonction  to get closest values when input values is not in the values range
    we assume that values are single on array. if the same value is repeated 
    will take the first index and the value at that index 
    
    :param values_range: values to get 
    :type values_range: array_like   
                 
    :param input_value: specific value
    :type input_value: float,

    :returns: the closest value and its index
    :rtye: float 
    
    """

    values_range = np.array(values_range)
    if np.all(values_range <0) : # if all element less than zero than convert 
                                    #inoput value to negative value 
        if input_value >0 : input_value *=-1  # for depth select purpose 
        
    if input_value < values_range.min(): 
        print('--> ! Input value ={0} is out  the range, min value = {1}. Value'\
              ' should be reset to ={1}.'.format(input_value, values_range.min()))
        
        warnings.warn('Input value ={0} is out  the range ! min value = {1}'.\
                      format(input_value, values_range.min()))
        _logger.debug('Input value ={0} is out  the range ! min value = {1}'.\
                      format(input_value, values_range.min()))  
            
        input_value = values_range.min()
    elif input_value > values_range.max(): 
        
        warnings.warn('Input value ={0} is out  the range ! max value = {1}'.\
                      format(input_value, values_range.max()))
        _logger.debug('Input value ={0} is out  the range ! max value = {1}'.\
                      format(input_value, values_range.max()))
            
        input_value = values_range.max()
        print('!--> Input value ={0} is out  the range , min value = {1}. Value'\
              ' should be reset to ={1}.'.format(input_value, values_range.max()))
        
    if input_value in values_range : 
        indexes,*_= np.where(values_range==input_value)
        if len(indexes)==1 : 
            index =int(indexes )
        elif len(indexes)>1 : # mean element is repeated then take the first index
            index = int(indexes)[0]
            
        value = values_range[index]
        return value , index 
    
    elif values_range.min() < input_value < values_range.max(): 
        # values_range = sorted(values_range)
        for ii, xx in enumerate(values_range): 
            if xx > input_value : 
                #compute distance : 
                d0 = abs(input_value-xx) # make the diffence between distances 
                d1= abs(values_range[ii-1]-input_value) # nad take the short
                if d0 < d1 : 
                    return  xx , ii
                elif d0> d1 : 
                    return   values_range[ii-1], ii-1
                

 
def ascertain_layers_with_its_resistivities(real_layer_names , 
                                            real_layer_resistivities): 
    """
    Assert the length of of real resistivites with their corresponding layers.
    If length of resistivities of larger that the layer's names, then will add
    layer that match the best the remained resistivities. If the length 
    of layer is larger than resistivities , to avoid miscomputation , will
    cut out this more layer and work only the length of resistivities provided.
    
    Parameters
    ----------
        * real_layer_names: array_like , list 
                    list of input layer names as real 
                    layers names encountered in area 
                    
        * real_layer_resistivities :array_like , list 
                    list of resistivities get on survey area
                
    Returns 
    --------
        list 
            real_layer_names,  new list of input layers 
    """
    # for consistency put on string if user provide a digit
    real_layer_names =[str(ly) for ly in real_layer_names]      
    
    if len(real_layer_resistivities) ==len(real_layer_names): 
        return real_layer_names
    
    elif len(real_layer_resistivities) > len(real_layer_names): 
         # get the last value of resistivities  to find the structres names ans structures sresistivities 
        sec_res = real_layer_resistivities[len(real_layer_names):]        
       
        geos =Geodrill.get_structure(resistivities_range=sec_res) # get the name of structure as possible 
        if len(geos)>1 : tm = 's'
        else :tm =''
        print('---> !We added other {0} geological struture{1}.'\
              ' You may ignore it.'.format(len(geos), tm))
        real_layer_names.extend(geos)                       # then , extend the list 
        return real_layer_names 
    elif len(real_layer_names) > len(real_layer_resistivities): 
        real_layer_names = real_layer_names[:len(real_layer_resistivities)]        # truncated the list 
        return real_layer_names

def geo_length_checker(main_param, optional_param, force =False, 
                param_names =('input_resistivities', 'input_layers'), **kws): 
    """
    Geo checker is a function to check differents length of different params.
    
    the length of optional params  depend of the length of main params . if
    the length of optional params is larger than the length of main 
    params, the length of optional params will be reduce to the length of main params .
    if the optional params  length is shorther than the length of
    main params, will filled it either with "None" if dtype param is string
    of 0.if dtype params is float or 0 if integer.if Force  is set True , 
    it will absolutely check if the main params and the optional params have
    the same length. if not the case ,   will generate an error occurs.
  
    Parameters 
    ------------
        * main_param : array_like, list 
                 main parameter that must took 
                 its length as reference length 
                 
        * optional params : array_like, list 
                 optional params, whom length depend 
                 to the length of main params
                 
        * param_names : tuple or str 
                 names of main params and optional params 
                 so to generate error if exits.
                 
        * fill_value: str, float, optional  
                Default value to fill thearray in the case where 
                the length of optional param is 
                less than the length of the  main param .If None ,
                will fill according to array dtype
            
    Returns 
    --------
        array_like 
           optional param truncated according to the man params 
    """
    add_v =kws.pop('fill_value', None)
    

    if isinstance(main_param, (str, float, int)):
        main_param = np.array([main_param])
    if isinstance(optional_param, (str, float, int)):
       optional_param = np.array([optional_param])
    if isinstance(main_param, (list, tuple)):
        main_param =np.array(main_param)
    if isinstance(optional_param, (list, tuple)):
        optional_param =np.array(optional_param)
            

    mes=''
    if len(optional_param) > len(main_param): 
        mes ="---> Note ! {0} will be truncated to length = {1} as the same length of {2} .".\
            format(param_names[1], len(main_param),param_names[0] )
        warnings.warn(mes)
        
        optional_param= optional_param[:len(main_param)]
        if force is True : 
            mess = '--> force argument  is set <True>, Can not truncate {0}'\
                ' = {1} to match the length of {2} = {3}.'\
                .format(param_names[1], len(param_names[1]), param_names[0],
                        len(param_names[0]))
            raise CSex.pyCSAMTError_parameter_number(mess)

    elif len(optional_param) < len(main_param) : 
        if force is True : 
            mess = '--> force argument  is set <True>, Can not fill value of {0} '\
                'to match the length of {1} = {2}.'\
                .format(param_names[1], param_names[0], len(param_names[0]))
            raise CSex.pyCSAMTError_parameter_number(mess)
            
        if add_v is not None : 
            add_v =[add_v for vv in range( len(main_param)-len(optional_param))] # repeat the value to add 
            add_v = np.array(add_v)
        if add_v is None :
             if optional_param.dtype not in [ 'float', 'int'] : 
                 add_v =['None' for i in range(len(main_param)-len(optional_param))]
                 
             else : 
                 for type_param, fill_value in zip([ 'float', 'int'],[ 0., 0] ): 
                     if  type_param  == optional_param.dtype :
                         add_v =[fill_value for i in range(len(main_param)-len(optional_param))]

        mes ="--> !Note Length of {0} is ={1} which  length of {2} is ={3}.'\
            We'll add {4} to fill {5} value.".\
            format(param_names[1], len(optional_param), param_names[0],
                   len(main_param), add_v[0], param_names[1])
            
        warnings.warn(mes)
        optional_param= optional_param.tolist()
        optional_param.extend(add_v)


    
    return np.array(optional_param)
    
                                                                                 


# if __name__=="__main__" :

#     path = r'F:/__main__csamt__\oasis data\OASISWORKS\all_data'
#     geo_surface_obj = Geosurface( path =path )
    
#     geo_surface_obj.read_oasis_files()
    # print(geo_surface_obj.oasis_data)
    # DICO= geo_surface_obj.global_dico
    
    # TEST0 = geo_surface_obj.K1_cor_oas
    

    # from scipy.signal import argrelextrema as extrema 
    
    # maxi = extrema(geoDict['S07'][:, 1], np.greater )
    # mini = extrema(geoDict['S07'][:, 1], np.less )
    # print(geoDict2['S00'][:,0])
    


    
    
    
    
    

                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    


