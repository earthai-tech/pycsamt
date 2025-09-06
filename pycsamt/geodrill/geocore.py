# -*- coding: utf-8 -*-
#       Create:on Sat Sep 19 12:37:42 2020
#       Author: Kouadio K.Laurent<etanoyau@gmail.com>
#       Licence: LGPL
"""  
Module Geocore 
===============
Module deals with Occam2D inversion files, the geological data, the geological 
structural infos and borehole data. Can plot the stratigraphy model and 
pseudostratigraphic logs under each station. Module can also build borehole
data collected in exploration area. It also generates a report of each geological 
structure colllected on survey area and compare it with its corresponding 
into a conventional geological codes < USGS(US Geological Survey ), pattern and colors
for plotting purposes or else. This is usefull for module `GeoStratigraphic_`
for topping 'fake' or 'unknow' structures into a pseudostratigraphic model.
  
            ...
"""
from __future__ import division 
import os 
import copy 
import warnings
import datetime
import numpy as np 
import pandas as pd 

import pycsamt.bases as BS
import pycsamt.utils.geo_utils as GU
import pycsamt.geodrill.structural as STRL
from pycsamt.geodrill.geodatabase import GeoDataBase 
from pycsamt.modeling import occam2d
from pycsamt.site import Profile
from pycsamt.utils import func_utils as func
from pycsamt.utils import plot_utils as punc
from pycsamt.utils.plotdecorator import geoplot2d
from  pycsamt.utils.exceptions import ( 
	GeoError, 
	FileHanglingError, 
	GeoMemoryError
) 
try : 
    from pycsamt._csamtpylog import csamtpylog
    _logger=csamtpylog.get_csamtpy_logger(__name__)
except :
    pass

try : 
      import scipy 
      import scipy.stats as spSTAT
      scipy_version = [int(vers) for vers in scipy.__version__.split('.')] 
      if scipy_version [0] == 1 : 
          if scipy_version [1] < 4 :
              warnings.warn('Note: need scipy version 1.4.0 or more . '
                'It may probably get a trouble when import `stats` attribute'
                f'under version {scipy_version}. It may probably deprecated ',
                ImportWarning)
              
              _logger.warning('Note: need scipy version 0.14.0 or higher'
                       ' or for stats.linearegress. Under such version'
                          'it might not work.')
      stats_import =True 
          
except :
    warnings.warn(
        'Could not find `scipy.stats`, cannot use method linearregression.'
        'Please check installation you can get scipy from https://scipy.org/')
    _logger.warning('Could not find scipy.stats, cannot use method linearegress.'
                    'check installation you may get scipy from https://scipy.org/')
    
    stats_import =False 
    
try : 
    from pycsamt.__init__ import itqdm 
    if itqdm : 
        import tqdm
except: 
    itqdm =False 
    
    pass 

class Geodrill (object): 
    """
    Class to manage occam2D model files comnbined with geoogical informations  
    to create a stratigraphic log of exploration area.
    
    Each station is condidered as an attribute and  framed by two closest  
    points from the station offsets. The class deals with the true resistivity
    values collected on exploration area from drilling or geolgical companies. 
    Indeed, the input true resistivity values into the occam2d inversion data 
    could yield an accuracy underground map. The challenge to build a pseudolog 
    framed between two stations allow to know the layers disposal or supperposition
    from to top to the investigation depth. The aim is to  emphasize a  large 
    conductive zone usefful in of groundwater exploration. 
    The class `Geodrill` deals at this time with Occam 2D inversion files or 
    Bo Yang model (x,y,z) files. We intend to extend later with other external 
    softwares like the Modular System EM (MODEM) and else.
    It's also possible to generate output straighforwardly for others external 
    softwares  like Golder sofwares('surfer')or Oasis Montaj:
        
        - Surfer: :https://www.goldensoftware.com/products/surfer) 
        - Oasis: http://updates.geosoft.com/downloads/files/how-to-guides/Oasis_montaj_Gridding.pdf
         <https://www.seequent.com/products-solutions/geosoft-oasis-montaj/>
       
    Note: 
        If the user has a golder software installed on its computer, it 's 
        possible to use the output files generated here  to yield a 2D map so 
        to compare both maps to see the difference between model map (
            only inversion files and detail-sequences map after including the 
            input true resistivity values and layer names) 
    Futhermore, the "pseudosequences model" could match and describe better 
    the layers disposal (thcikness and contact) in underground than the 
    raw model map which seems to be close to the reality  when `step descent` 
    parameter is not too small at all. 


    Arguments
    ----------
        **model_fn** : str,  
                    full path to Occam model  file .                             
        **iter_fn** :  str,                
                    full path to occam iteration  file
                                                        
        **data_fn** :  str,  
                    full path to occam_data file 
        **input_resistivities** : list, array_like, 
                            True resistivity values
        **step_descent**: float,   
                         Param to force the model resistivity to
                         fit the truth input resistivity values. As a 
                         reference, the param is 20% of D0I.
        **input_layers**:  array_like,  
                        True input layers names- geological 
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
        
        >>> from pycsamt.geodrill.geocore import Geodrill 
        >>> path =os.path.join(os.environ ['pyCSAMT'],  
        ...                       'data', 'occam2D')
        >>> geo_obj = Geodrill( input_resistivities=[300, 500, 
        ...                                             1000, 2000, 4000, 6000],
        ...                       input_layers =['alluvium', 
        ...                                      'amphibolite','altered rock',
        ...                                                    'augen gneiss',
        ...                                                        'granite'],
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
    
    def __init__(self,
                 iter2dat_fn=None,
                 model_fn =None ,
                 data_fn=None ,
                 iter_fn=None ,
                 mesh_fn=None ,
                 bln_fn =None ,
                 verbose =0,
                 **kwargs):
                
        self._logging = csamtpylog.get_csamtpy_logger(self.__class__.__name__)
        
        self.iter2dat_fn = iter2dat_fn
        self.data_fn =data_fn 
        self.iter_fn =iter_fn 
        self.mesh_fn =mesh_fn 
        self.model_fn =model_fn 
        self.bln_fn =bln_fn  
        self.verbose =verbose
        self.input_layers =kwargs.pop('input_layers', None)
        self.input_resistivities =kwargs.pop('input_resistivities', None)
        self.doi =kwargs.pop('doi', '1km')
        
        self.elevation =kwargs.pop('elevation', None)
        self.etacoeff = kwargs.pop('etaCoef', 5)
        self.step_descent =kwargs.pop('step_descent', None)
        self.savepath =kwargs.pop('savepath', None)
        
        self.model_res =None 
        
        for key in self.geo_params :  
            setattr(self, key, None)

        # for key in list(kwargs.keys()):
        #     setattr(self, key, kwargs[key])
            
        if self.model_fn is not None and self.data_fn is not None\
            or self.iter2dat_fn is not None : 
            self.set_geodata()
    
    def set_geodata(self,iter2dat_fn =None ,  model_fn=None , data_fn=None , 
                    iter_fn=None , mesh_fn=None, **kwargs):
        """
        Readmodel data and collected from each site its value from surface to depth 
        
        1. Read with Occam 2D outputs files
        
        :Example: 
            
            >>> from pycsamt.geodrill.geocore import Geodrill
            >>> path_occam2d = os.path.join(os.environ ['pyCSAMT'],
            ...                                 'data', 'occam2D')
            >>> path_i2d =os.path.join(os.environ ['pyCSAMT'], 
            ...                            'data', '_iter2dat_2')
            >>> geo_obj =Geodrill(
            ...      model_fn  = os.path.join(path_occam2d,'Occam2DModel'), 
            ...      mesh_fn = os.path.join(path_occam2d, 'Occam2DMesh'), 
            ...      iter_fn = os.path.join(path_occam2d, 'ITER17.iter'), 
            ...     data_fn = os.path.join(path_occam2d, 'OccamDataFile.dat'))
            
        2. Read with only iter2dat file and bln file 
        
        :Example:
            
            >>> from pycsamt.geodrill.geocore import Geodrill
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
     
                if site_offset > model_offset.max (): # the case where 
                #offset is much larger than the mesh (very rare case )
                    new_index = int( np.where(model_offset==model_offset.max())[0])
                    new_off = model_offset.max ()
                else : 
                    for ii, di in enumerate(model_offset): 
                        if di > site_offset : # if value of model is 
                        #greater than check the closest site offset 
                            dx0 = round(abs(site_offset -di ),2) # compute the 
                            #distance betewn di and offset (dxmax )
                        # try co compute the previous and the next point distance 
                        # use try because sometimes , previous point is none 
                        #if model value is at index 0 
                        # if model value at index 0 , then dx0 = dmin 
                            try : dxmin =abs(model_offset[ii-1]-site_offset)
                            except : dxmin =round(dx0 , 2)
                            if dx0 <=  dx0 : 
                                new_off , new_index = di, ii  
                                break 
                            elif dxmin < dx0 : # take the preious 
                                 new_off, new_index = model_offset[ii-1], ii-1 
                                 break 
     
            # now get the list of resistivity between this index 
          
            if new_off == model_offset[0] : # the first site index then take ii+1, ii+2 
                dm0 = np.concatenate ((model_res[:,new_index].reshape(
                                          model_res[:, new_index].shape[0], 1), 
                                       model_res[:,new_index+1].reshape(
                                           model_res[:, new_index+1].shape[0], 1), 
                                       model_res[:,new_index +2 ].reshape(
                                           model_res[:, new_index+2].shape[0], 1)), 
                                       
                                       axis =1)
            elif new_off == model_offset[-1] : 
                dm0 = np.concatenate ((model_res[:,new_index-2].reshape(
                    model_res[:, new_index-2].shape[0], 1), 
                                       model_res[:,new_index-1].reshape(
                                           model_res[:, new_index-1].shape[0], 1), 
                                       model_res[:,new_index ].reshape(
                                           model_res[:, new_index].shape[0], 1)),
                                       axis =1)
            
            else : 
                dm0 = np.concatenate ((model_res[:,new_index-1].reshape(
                    model_res[:, new_index-1].shape[0], 1), 
                                       model_res[:,new_index ].reshape(
                                           model_res[:, new_index].shape[0], 1), 
                                       model_res[:,new_index +1].reshape(
                                           model_res[:, new_index+1].shape[0], 1)),
                                       axis =1)
                
            return site_name, dm0 

 
        # ----> attributes statements 
        if model_fn is not None :
            self.model_fn =model_fn 
        if data_fn is not None :
            self.data_fn =data_fn 
        if iter_fn is not None : 
            self.iter_fn =iter_fn 
        if mesh_fn is not None :
            self.mesh_fn =mesh_fn 
        
        if iter2dat_fn is not None :
            self.iter2dat_fn = iter2dat_fn 
        
        # then assert all input files 
        if self.iter2dat_fn is not None : 
            self._logging.info ('Read Geodata from  Bo Yang Iter2Dat files ')
            if self.bln_fn is None : 
                mess =''.join(['No (*bln) station file found. Need sites',
                               ' names and sites locations.',
                               ' Use mode "Iter2Dat"  :',
                               '<from pycsamt.modeling.occam2d import Iter2Dat>',
                               ' : to write a *bln file ', 
                               ' use a Profile module :',
                               ' < from pycsamt.core.cs import Profile > .', 
                               ' You can also use occam2D output files ',
                               'like model, mesh, data, and iteration files.'])
                warnings.warn('! Error reading *.bln flile !' + mess)
                self._logging.error(mess)
                raise GeoError(
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
            for attr , ins  in zip([
                    'model_fn', 'iter_fn', 'data_fn', 'mesh_fn'],
                    [self.model_fn, self.iter_fn, self.data_fn, self.mesh_fn]): 
                if ins ==None : 
                    msg =' !No {0} file found ! Please provide the Occam 2D {0} file'.\
                        format(str(attr).replace('_fn', ''))
                    self._logging.error(msg)
                    warnings.warn(msg)
                    raise GeoError(msg)
                
            # Recreate object and get ncessaries attributes 
            model_obj= occam2d.Model(model_fn =self.model_fn,
                                     iter_fn =self.iter_fn, 
                                    mesh_fn =self.mesh_fn )
            data_obj =occam2d.Data(data_fn= self.data_fn)
            
            # then get important attributes and re
            self.model_x_nodes =model_obj.model_station_offsets
            self.model_z_nodes =model_obj.model_depth_offsets
            self.model_res = model_obj.model_resistivity 
            # get from data sites names and station location and elevation 
            self.station_names =data_obj.data_sites 
            self.station_location = data_obj.data_offsets
            
            #get rms and roughness value
            self.model_rms = model_obj.model_rms
            self.model_roughness =model_obj.model_roughness
 
        # build for each site an array of three values framed into the main sites 
        #    get the the offset of the site 
        # build a dictionnary of geoname
        
        #------Mange DOI (investigation depth) --------
        
        # check the depth of investigation and ge the new_nodel matrix 
        
        if self.doi is None : self.doi = 1000.
        elif self.doi is not None :
            self.doi =punc.depth_of_investigation(self.doi) # return meter value

        # for consistency let check how the depth is ranged :
            #Depth is  minimum to max depth , if not flip data 
        if self.model_z_nodes [0] > self.model_z_nodes[-1]: 
            self.model_z_nodes= self.model_z_nodes[::-1]
        
        # chech the doi data and compare to max depth data 
        
        if self.doi > self.model_z_nodes.max (): 
            mess='Maximum Input doi value = {0} m is larger'\
                  ' than model depth = {1} m. Investigation'\
                      ' depth "doi" will be resetting at = "{1}"m.'.\
                              format(self.doi, self.model_z_nodes.max())
             # resseting doi to maximum value in model depth 
            self.doi = self.model_z_nodes.max()        
            
            warnings.warn(mess)
            self._logging.debug (mess)
        
        if self.doi in self.model_z_nodes: 
            dep_index = int(np.where(self.model_z_nodes== self.doi)[0])
            
            if self.doi == self.model_z_nodes[-1]:  # if doi is the  max depth , 
                self.__setattr__('geo_depth', self.model_z_nodes)
                 # resize model _resistivity 
            elif self.doi  <     self.model_z_nodes[-1]:        
                self.model_res = np.resize(self.model_res, (
                    dep_index +1 , self.model_res.shape[1]))
                # resize the large depth to doi max  
                self.__setattr__('geo_depth', self.model_z_nodes[:dep_index+1]) 
            
            print('---> resetting model doi to = {0} m depth !'.format(self.doi))
        else :
            for index, dep  in enumerate(self.model_z_nodes) : 
                # seek the depth index that match better the  input 
                if dep  > self.doi  :       #doi 1014>1000 :index=23,
                    dep_index = index-1     # get the previous index dep_index = 22
                    self.model_res = np.resize(self.model_res, (
                        index, self.model_res.shape[1])) 
                    self.doi = self.model_z_nodes[dep_index]
                    print('---> resetting model doi to ={0} m depth !'.format(
                        self.doi))
                    
                    self.__setattr__('geo_depth', self.model_z_nodes[:index])
                    break 
                
        # build geo_station_names and geodict rho 
        self.geo_d ={}
        if isinstance(self.station_names , (tuple, np.ndarray)): 
                self.station_names=self.station_names.tolist()
    
        for name in self.station_names : 
    
            site_offset_value = self.station_location[
                self.station_names.index(name)]
            
            geo_name , geodat = frame_each_site_into_3_offsets(
                site_name= name, 
                model_res = self.model_res, 
                model_offset=self.model_x_nodes, 
                site_offset=site_offset_value)
            
            self.__setattr__('geo_name_{0}'.format(geo_name), geodat)
            # convert resistivities value from log10 rho to ohm meter 
            self.geo_d[name]= np.power(10, getattr(
                self, 'geo_name_{0}'.format(geo_name))) # 
            
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
            self.set_geodata() # load attributes 
        
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
                raise GeoError(
                    'Value provided as resistivity '
                    'range must be an array of float number.')
                                                                
            tem=[]                   # loop the value range 
            for mm in value_range :   # substract rowlines from 
                op = rowlines-mm        #value range to seek the minimum close
                op=np.array([abs(aj) for aj in op]) # keep absolute value for difference 
                tem.append(op)                      # keep it on temporray list  
            ts = func.concat_array_from_list(tem)   # concat list to axis =0 order 
            
            for ii in range(len(rowlines)):         # loop the rowlines line now 
                us=ts[:,ii].argmin()                # keep minimum index 
                rowlines[ii]=value_range[us]        # change the resistivity values 
                
            return rowlines
        
        
        def ascertain_input_layers(input_layers, input_rho): 
            """
            Function to assert the length of input layers and the 
            input resistivities to avoid miscomputation .
            
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
            # for consistency put on string if user provide a digit
            ilay =[str(ly) for ly in input_layers]      
            if len(input_rho) ==len(ilay): 
                return ilay
            #get the last value of resistivities  to find
            elif len(input_rho) > len(ilay): #the structres names ans structures
                sec_res = input_rho[len(ilay):]         #  sresistivities  
                geos =Geodrill.get_structure(resistivities_range=sec_res) 
                if len(geos)>1 : tm = 's'
                else :tm =''
                print('---> !We added other {0} geological '
                      'structure{1}. You may ignore it.'.format(len(geos), tm))
                ilay.extend(geos)                       # then , extend the list 
                return ilay 
            elif len(ilay) > len(input_rho): 
                ilay = ilay[:len(input_rho)]            # truncated the list 
                return ilay 
        
        
        # if range of resistivity is provided then use it for average 
        if input_resistivity_range is not None : 
            self.input_resistivities = input_resistivity_range 
        if self.input_resistivities is None and self.input_layers is not None : 
            warnings.warn(
                ' !Without any input resistivities values,'
                ' we can not set only your layer names.'
                ' Bring more details about your layers by adding their '
                'corresponding resistivities values.')
            # print('-->!')
             #abort the input layers by renitializing to NoneType values
            self.input_layers=None      
            
         
        auto_mess=''            # keep automatic message 
        if self.input_resistivities is None :
            if self.input_layers is None : 
                n_layers =7.  # number of slices layer juged by the program 7 
            else : n_layers =len(self.input_layers)
            
            mess="".join([
                "Resistivity range is not provided . Sites Depth will",
                " be cout out into {0} slices as possible layers ", 
                "  below the site. If the number of slices doesnt",
                " suit the purpose , please change the number ",
                "  of slices using argument <input_layers> ",
                "to provided the real layer's names."])
            
            warnings.warn(mess.format(int( n_layers))) 
            self._logging.info(mess.format(int( n_layers)))
           
          
            # get the model resistivities minimum and maximum from selected doi,  
            # it is much  betterthan  to select min max into the  
            #global resistivities models.  
            
            minmax=[(res_values.min(), res_values.max()) 
                    for stn, res_values  in self.geo_d.items()] #
            maxres= max([mm[1] for mm in minmax])
            minres= min([mm[0] for mm in minmax])
            
            self.input_resistivities =np.linspace(
                minres, maxres, int( n_layers))
            self.input_resistivities  = np.around(
                self.input_resistivities,2) 
            auto_mess ='Automatic'
 
            
        if isinstance (self.input_resistivities  , (tuple, list)): 
            self.input_resistivities  = np.array(self.input_resistivities )

        elif isinstance(self.input_resistivities , (float, int, str)): 
            try : 
                self.input_resistivities = float(self.input_resistivities )
            except : 
                raise GeoError(
                    'Can not converted value <%s> into float number.'\
                    ' Please provide a float number.'%self.input_resistivities )
                
        # Display infos
        if self.verbose >0 :
            print('**{0:<37} {1} {2}'.format('{0} Layers sliced'.format(
                auto_mess ), '=' , len(self.input_resistivities  )))
        if self.verbose>3:                                
            print('**{0:<37} {1} {2} (Ω.m)'.format('{0} Rho range'.format(
                auto_mess ),'=' ,  tuple(self.input_resistivities  ) ))
            print('**{0:<37} {1} {2} {3}'.format(
                ' Minimum rho ','=' ,self.input_resistivities.min(), 'Ω.m' ))
        if self.verbose >0 :
            print('**{0:<37} {1} {2} {3}'.format(' Maximum rho ',
                                                 '=' ,
                                                 self.input_resistivities.max(), 
                                                 'Ω.m' ))
       
        # so to get the structures for each input resistivities 
                 # ascertain unput layers first if no misleading value is inputted 
        if self.input_layers is not None : 
            # rebuild new list of layer by adding 
            #necessary strutures or substracting unecessaries structures
            formations= ascertain_input_layers(
                input_layers= self.input_layers,
                input_rho= self.input_resistivities)
        #     self.depth_range =depth_range 
        else :
            formations = Geodrill.get_structure(
                self.input_resistivities)
        
        self.input_layers=formations 

        # self.__setattr__('geo_dict_rho', nOne)
        self.geo_drr=copy.deepcopy(self.geo_d)             
        
        
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
        # only single value provided than put on list.
        if isinstance(resistivities_range, (float, str, int)): 
            try : 
                resistivities_range=[float(resistivities_range)]
            except : 
                raise GeoError(
                    'Can not convert <%s> to float number !Input resistivity '
                     'must be a float number not <%s>!'% (resistivities_range,
                                                          type(resistivities_range)))
         # for consistency , check again input values 
        if isinstance(resistivities_range, (tuple, list, np.ndarray)): 
            try :  resistivities_range =np.array(
                [float(ss) for ss in resistivities_range])
            except : 
                raise GeoError(
                    'Input argument provided is wrong. '
                    'Please check your resistivities '
                    ' range, values must be a float number.')
          # in fact append the idde of resistivities values located  
        # so to be sure that the resistivities will match exactly the
        # layer found in geo_electricl_property of rocks   
        geo_structures=[]                  
                                            
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
        It show the lowest point and the maximum point averaged . Function 
        averaged rho value between local maximum and local minima values . 
        if data  values of station are located on columnlines , set transpose 
        to True then rotate the matrix to find minima and maxima  locals value 
        then calculated averaged rho after will return matrix transpose
         as the same shape as inputted . .Defaut is **False** . 
     
        :param data_array: data of resistivities collected at the site point  
        :type data_array: ndarray  
        
        """
        if data_array.dtype not in ['float', 'int']:  
            try :
                data_array= data_array.astype('float64')
            except : 
                warnings.warn("It seems somethings wrong"
                              " happened during data conversion to float values.")
                raise GeoError(
                    'Could not convert value'
                    ' to float numbers , Please check your data!')
            
      
        if transpose is True :  
            data_array =data_array.T
        
        exem =[punc.average_rho_with_locals_minmax(rowlines)
               for rowlines in data_array ]         
        new_data = func.concat_array_from_list(exem)
        if transpose is True :          
            new_data = new_data.T
            
        return new_data
        
  
    def geo_build_strata_logs (self, input_resistivities=None, 
                               input_layers =None, 
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
        
        .. seealso:: https://www.eoas.ubc.ca/ubcgif/iag/foundations/properties/resistivity.html
                list is not exhaustive and depend of the geological formations of 
                survey area. 
                
        .. note:: list is not Exhaustive, use the data base script to populate 
                    most of goeological electrical properties.
        """
        etaCoef = kwargs.pop('etaCoef', None)
        def get_conductive_zone (dep_array, rho_array, step_in_deep):
            
            """
            Get a conductive zone is important for many purposes .In the case 
            of groundwater exploration for instance, sometimes in deeper, 
            because of heterogeneities of structures in underground , it is 
            much difficult to take some minima rho as a conductive area or 
            conductiv productive veification drill point. To be sure that this 
            point is really among a good saturated zone, it is better to get   
            averaged rho in some distance and to see how the layer resistivity 
            value goes on on this range, whether it's an effective conductive
            zone or overlapping zone. In addition, fixing resistivities values
             at specific distance depth allow us to detect  the unsatured
            zone as weel as the satured zone at the set of investigation depth 
            imaged. This technique allow us to build a pseudo specific strata
             that could match the zone .
          
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
                    raise GeoError(
                        'Could not convert depth value ={} to float.'
                        ' Please check your value.'.format(step_in_deep))

            if step_in_deep< dep_array.min(): 
                raise GeoError(
                    'Value provided ={0} m is less than the minimum'
                    ' depth ={1} m.'.format(step_in_deep, dep_array.min()))
            
            if step_in_deep > dep_array.max(): 
                raise GeoError(
                    'Value provided is = {0} m is greater than '
                    'maximum depth ={1}m.'.format(step_in_deep, dep_array.max()))
                
            _init_depth =step_in_deep
            
            for index , depth in enumerate(dep_array):
                if depth <= step_in_deep :      # value less than step descent 
                    v.append(depth)            
                    r.append(rho_array[index])
        
                if depth > step_in_deep :       
                    if v !=[]:                 
                        # rebuild resistivities values with rho averaged 
                        dm.append(np.repeat(np.array(v).mean(), len(v)))  
                        rm.append(np.repeat(np.array(r).mean(), len(r)))
                        step_in_deep += _init_depth              
                        v=[depth]       # initialise new list by adding 
                        r=[rho_array[index]] #the index value greater one 
                      
                if depth ==dep_array[-1]:
                    # it length last value ==1, means is the last value of depth        
                    if len(v)==1 :         
                        dm.append(dep_array[index])
                        rm.append(rho_array[index])
                    elif len(v) !=1 :    # averaged the reamin rho values  
                        dm.append(np.repeat(np.array(v).mean(), len(v)))
                        rm.append(np.repeat(np.array(r).mean(), len(r)))
              
        
            return np.hstack(tuple(rm)),np.hstack(tuple(dm))
            

        self.__setattr__('qc', .0)          # set quality control in the data 
        
        #------  STATEMENT OF INPUT ARGUMENTS------------------ 
        # check data files if provided (Dont need to check all )
        if input_resistivities is not None :
            self.input_resistivities= input_resistivities 
        
        if self.model_fn is None or self.iter2dat_fn is  None :  
            self.geo_replace_rho()
            
        if input_layers is not None :
            # rebuild input_layers 
            self.input_layers = assert_len_layers_with_resistivities(
                real_layer_names= input_layers, 
                real_layer_resistivities=self.input_resistivities) 

        # if self.etacoeff is not None : 
            # self.step_descent = self.doi / int(self.etacoeff)
        if etaCoef is not None : 
            self.etacoeff = etaCoef 

        if step_descent is not None : 
            self.step_descent = step_descent

        if self.step_descent is None : 
            self.step_descent = self.doi / int(self.etacoeff)

        if self.step_descent > self.doi : 
            self._logging.info('Step descent is = {0} much larger than doi ={1}.'
                               'Value is ressetting to 20% of DOI = {2}m.'
                               .format(round(self.step_descent,3),
                                       round(self.doi,3),
                                       .2 * self.doi ))
            
            self.step_descent = .2 * self.doi 
        # recompute a new block coefficient if step descent is given 
        
        if self.step_descent is not None : 
            self.etacoeff = int(self.doi/ self.step_descent)
        

        # ---------------get the geodictionnary from from geodata average ----
        self.geo_daver ={}
        for stn, vrho in self.geo_d.items(): 
            # transpose data so to read rowlines S01,S02 etc 

            self.geo_daver[stn] =  Geodrill.get_average_rho( data_array= vrho,
                                                            transpose =True)
            
        self.qc += .25 
        #----- set attribute of geo_dstep_descent -----------------------------
        self.geo_dstep_descent={}
        for stn , geo_sd_values in self.geo_d.items(): 
            # for ii in range(3) :
                
            mv= [get_conductive_zone(dep_array= self.geo_depth,
                                     rho_array=geo_sd_values[:,ii], 
                                     step_in_deep= self.step_descent)[0] 
                 for ii in range(3)] # concat three array 
                
            self.geo_dstep_descent[stn]=func.concat_array_from_list(
                mv, concat_axis=1)
                
        #----build a dict of pseudosequences layer and resistivities ----------
        
        # build at dictionnary of strata from resistivities
        self._logging.info ('Build the pseudosequences of strata.')
        
        # if constrained_electrical_properties_of_rocks is False : 
            # build layer according to geo_drr (geo_replacedrho)
            
        self.geo_dpseudo_sequence, self.geo_dpseudo_sequence_rho =[
            {} for ii in range(2)]
        sc=[]
        
        self.__setattr__('geo_secure_pseudo_sequence', None)
         
        for stn , georr in self.geo_drr.items(): 
            # for jj in range(3): 
            if stn =='S00' :            
                # first station id framed into three started at the left  
                svm = punc.build_resistivity_barplot(depth_values=self.geo_depth,
                                           res_values=georr[:,0])
            elif stn == self.station_names[-1] :  
                svm = punc.build_resistivity_barplot(depth_values=self.geo_depth,
                                           res_values=georr[:,2])
            else :         
                svm = punc.build_resistivity_barplot(depth_values=self.geo_depth,
                                           res_values=georr[:,1])

            sc.append(svm[-1])

            self.geo_dpseudo_sequence[stn]= svm[0]     
            self.geo_dpseudo_sequence_rho[stn]=svm[1]  
            
        self.geo_secure_pseudo_sequence = np.array(sc)
        
        #---------------- Quality control --------------------------------.---
        mess =' Your data pass safety the Quality Control !'
        
        if len(self.geo_secure_pseudo_sequence) == len(self.station_location):
            self.qc +=.25 
        if np.all(self.geo_secure_pseudo_sequence== self.doi ):
                                                               
            self.qc += .25 
        else :
            warnings.warn('Data provided are inconsistencies ,'
                          ' Please you may try to undestand carefuuly the code.')
            self._logging.debug('Data provided are inconsistencies ,'
                          ' Please you may try to understand carefully the code.')
            
        if self.qc ==1. : 
            print('---> {}.'.format(mess))
            
        print('**{0:<37} {1} {2} {3}'.format(' QC flux rate','=' , 
                                             100. * self.qc, '%' ))
        
        #-rewrite info with real layers resistivities 
        # if provided and layers names ---
        # self.input_layers, _ , _= pycsamt.geodrill.get_geo_formation_properties(
        # structures_resistivities= self.input_resistivities,)

        #-------------------Print Info ----------------------------------------
        print('**{0:<37} {1} {2}'.format(' Number of layers','=' ,
                                         len(self.input_layers) ))
        print('**{0:<37} {1} {2} {3}'.format(' Step descent','=' ,
                                         round(self.step_descent,3), 'm' ))
        print('**{0:<37} {1} {2}'.format(' Block coefficient','=' ,
                                         int(self.etacoeff)))
        
        if self.verbose >7:
            print('-'*77)
            print('{0:<25}{1:<25} {2:<25}'.format('Structure', 'Rho mean value (Ω.m)', 
                                                  'Rho range  (Ω.m)'))
            print('-'*77)
      
            for ij  , (ires, geos) in enumerate(zip (
                    self.input_resistivities,self.input_layers )): 
                if len(self.input_resistivities ) >1 :
                    if ij==0 : 
                        rho_rg = '{0}-{1}'.format(
                            ires, self.input_resistivities [ij+1])
                        rho_mean = np.around(
                            (ires+ self.input_resistivities [ij+1])/2,2)
                    elif ires == self.input_resistivities[-1]: 
                          rho_rg = '{0}-{1}'.format(self.input_resistivities[ij-1],
                                                    self.input_resistivities[ij])
                          rho_mean =np.around( (self.input_resistivities[ij-1]+\
                                                self.input_resistivities[ij])/2,2)
                    else : 
                        rho_rg = '{0}-{1}'.format(self.input_resistivities[ij-1],
                                                  self.input_resistivities[ij+1])
                        rho_mean = np.around( (self.input_resistivities[ij-1]+\
                                                self.input_resistivities[ij+1])/2,2)
                if len(self.input_resistivities)==1:
                    rho_rg = '{0}'.format(ires)
                    rho_mean = ires
                    
                print('{0:<25}{1:<25} {2:>25}'.format(geos, rho_mean,rho_rg ))
            print('-'*77) 
        
        
    @staticmethod 
    def get_geo_formation_properties (structures_resistivities, 
                            real_layer_names=None,
                            constrained_electrical_properties_of_rocks=True, 
                             **kwargs):         
        """
        Get the list of stuctures and their names of after replaced , 
        flexible tools.  wherever structures provided, the name, color , 
         as well as the pattern .if constrained electrical properties if True , 
         will keep the resistivities with their  corresponding layers
        as reference. If the layer names if found on the data Base then , 
        will return its pattern and color else defaultcolor is black and 
        patter is "+.-". If `constrained_electrical_properties_of_rocks`
         is False , will check under data base to find the 
        resistivities that match better the layers. 
        
        Parameters
        ------------
            * structures resistivities  : array_like, 
                                       resistivities of structures 
                                        
            * real_layer_names : array_like |list 
                               names of layer of survey area 
                               if not provided , will use resistivities 
                               to find the closet layer that match
                                the best the resistivities
                                
            * constrained_electrical_properties_of_rocks: bool 
                        set to True mean the realy_layer is provided.
                        if not program  will enforce to False ,
                        will use default conventional layers
                        Default is False, assume to povided layer 
                        names for  accuracy.
        
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
        
        find_pattern , __f_db= False, -1             
        pattern, geo_color = [[] for i in range(2)] 
        
        geoformation_obj =STRL.Geo_formation()
        geof_names = geoformation_obj.names 

        if constrained_electrical_properties_of_rocks is True : 
            #check weither the layer of resistivities
            if real_layer_names is not None : 
               #are the same dimension and return new layer with same dimension
                real_layer_names = assert_len_layers_with_resistivities(
                    real_layer_names =real_layer_names,
                    real_layer_resistivities= structures_resistivities)
            else : 
                msg =''.join(["{constrained_electrical_properties_of_rocks} ",
                    "argument is set to <True> as default value ",
                    "  unfortunately , you did not provide any layer's",
                    " names with its resistivity values.", 
                    " We can not set layer' resistivities and layer'",
                    " names as reference data.",
                    " We will ressetting",
                    " {constrained_electrical_properties_of_rocks} to False.", 
                    "  However be sure that layer's provided automatically",
                    " don't match exactly the underground informations.", 
                    "To accurate informations, you need ABSOLUTELY",
                    " to provided layer' resistivities as well as ", 
                    "its names get on the field or any other firms."])
                
                print('---> Important Note :!' )
                text = punc.fmt_text(data_text=msg, fmt='+',return_to_line=90)
                print(text)
                warnings.warn(msg)
                
                constrained_electrical_properties_of_rocks = False 
           
        if constrained_electrical_properties_of_rocks is False : 
            if real_layer_names is None : 
                real_layer_names = Geodrill.get_structure(
                    structures_resistivities)
            elif real_layer_names is not None :  
                real_layer_names = assert_len_layers_with_resistivities(
                    real_layer_names =real_layer_names,
                     real_layer_resistivities= structures_resistivities) 

        #-------CONNECT TO GEODATABASE -----------------------------
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
            
            for res,  lnames  in zip ( structures_resistivities,
                                      real_layer_names ) :
                    geodatabase_obj. _get_geo_structure(
                        structure_name =lnames) 
                    if geodatabase_obj.geo_structure_exists is True : 
                    
                        if geodatabase_obj.colorMPL != 'none':  
                            # rebuild the tuple of mpl colors blocked 
                            lnames_color = tuple([ float(ss) 
                                   for ss in  geodatabase_obj.colorMPL.replace(
                                             '(', '').replace(')', ''
                                                              ).split(',')])

                            geo_color.append(lnames_color) # append rgb colors 
                        
                        else : 
                            mess = 'Sorry ! {0} matplotlib color is '\
                                ' not filled yet in our DataBase.'.format(
                                    lnames)
                            warnings.warn(mess)
                            geo_color.append(unknow_layer_color) 
                            
                        if geodatabase_obj.hatch != 'none': 
                            
                            pattern.append(geodatabase_obj.hatch.replace(
                                '(', '').replace(')', ''))
                        
                        else : 
                            mess = 'Sorry ! {0} matplotlib pattern'\
                                ' is not filled yet in our DataBase.'.format(
                                    lnames)
                            warnings.warn(mess)
                            pattern.append(unknow_layer_pattern)
                        
                    elif geodatabase_obj.geo_structure_exists is False   : 
                        # then aborted the process and go to search 
                        #into strata and strutral modules 
                        __f_db=0
                        
                        mess =' Actually {0} does not exist in our dataBase,'\
                            ' we alternatively use other ways '\
                                'to suitable respond to your request.'.format(
                                    lnames)
                        warnings.warn(mess)
                        _logger.debug(mess)
                        break
   
            geodatabase_obj.manage_geoDataBase.closeDB() # close de DB
            
        #----SEARCH ON STRATA AND STRUCTURAL MODULES ----------------------*
        if geodatabase_obj.success ==0 or  __f_db == 0 : 
            # clean initialise  and restart searching 
            geo_color , pattern = [[] for i in range(2)] 
            
            for res,  lnames  in zip ( structures_resistivities,
                                      real_layer_names ) :
                # if exist on geoformation array names then
                # get the index and the coorsponding values 
                if lnames in geof_names :  
                   
                   if ' ' in lnames  : new_lnames = lnames.replace(' ', '_') 
                   else : new_lnames =lnames
                   
                   col = getattr(geoformation_obj, new_lnames)['color'] 
                   geo_color.append(col)
    
                elif lnames not in geof_names : 
                   if lnames.lower() in STRL.geo_pattern.pattern.keys():
                       geo_color.append(STRL.geo_pattern.pattern[lnames.lower()][1])
                   else :
                       geo_color.append( unknow_layer_color)
      
                 # actually will use the default pattern : Then check the 
                 #resistivities inot geoelectrical 
                 #    propertty. if structures found , then took it pattern 
                for ii, ( keyprop , resprops)  in enumerate( 
                        Geodrill.geo_rocks_properties.items ()) :
                    #resprops)[-1] <= res <= resprops[0] :
                    if min(resprops) <= res <= max(resprops) : 
                        # split and search for their appropriate pattern
                        if '/' in keyprop:  
                            if res <= np.array(resprops).mean() :
                                keyprop=keyprop.split('/')[1]
                            else : keyprop=keyprop.split('/')[0]
                            
                        pat =STRL.geo_pattern.pattern[keyprop][0]     
                        pattern.append(pat)
                        find_pattern =True                 
                        break  # dont need to continue to take other pattern 
                        
                    if ii == len(Geodrill.geo_rocks_properties.items ()) \
                        and find_pattern is False : # mean no pattern 
                    # is found after looping all items of dictionnary 
                        pattern.append( unknow_layer_pattern)
                        
                find_pattern =False  # resetting flag to False (switch off )

        return real_layer_names, geo_color , pattern      
                       
    def to_golden_software(self, input_resistivities=None,input_layers =None, 
                            step_descent = None, filename =None,
                            savepath =None ,**kwargs ): 
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
        etaCoef                int              number of sperated model blocks
        ===================  ================  ================================
        
        """
        elevation = kwargs.pop('elevation', None)
        write_negative_depth=kwargs.pop('to_negative_depth', True)
        etaCoef = kwargs.pop('etaCoef', None)
        scale =kwargs.pop('scale', 'm')
    
        if scale is None : scale= 'm'
        
        # rescale data 
        if scale =='m': 
            df =1. 
        elif scale =='km': 
            df =1000. 
        else : df =1.
        
        if savepath is not None: 
            self.savepath =savepath 
        if input_layers is not None :
            self.input_layers = input_layers 
        if input_resistivities is not None : 
            self.input_resistivities  = input_resistivities 
            
        if etaCoef is not None : 
            self.etacoeff = etaCoef     
        if step_descent is not None :
            self.step_descent= step_descent 
        
        if elevation is not None : 
            self.elevation = elevation 
        
        
        if self.elevation is None : 
            mess='!Elevation is not provided. We gonna set to 0.'
            print('-->'+mess)
            self._logging.debug(mess)
            
            if self.station_location is not None : 
                self.elevation  = np.repeat(0., len(self.station_location))
        elif self.elevation is not None : 
            
            self.elevation = func.geo_length_checker(
                main_param= self.station_location,
                optional_param = self.elevation, 
                param_names = ('station location', 'elevation'),
                fill_value=0.)
        
        # now read files 
        if  self.input_resistivities is not None : 
            self.geo_build_strata_logs()
        else : 
            mess ="".join(["Need ABSOLUTELY an input resistivities get"
                           " on the field or other companies ", 
                           " Input resistivities MUST provided in order"
                           " to take data as reference resistivities.", 
                           " please provided a least",
                           " a truth resitivities of one layer."])
            
            warnings.warn(mess)
            self._logging.error(mess)
            raise GeoError(
                "Can not write details sequences log files "
                "! Please provided at least one"
                " truth layer resistivity.")
        
        matrix_rhoaver ,matrix_rhorr,\
            matrix_rho_stepdescent =[[]for i in range(3)]
        for site in self.station_names : 
            #station id to read is the first index =0
            if site ==self.station_names[0]: 
                index =0 
            if site ==self.station_names [-1] : 
                index =2 # station id  to read is the 2index 
            else : index =1  # read th middle array 
    
            matrix_rhorr.append(self.geo_drr[site][:, index])
            matrix_rhoaver.append(self.geo_daver[site][:, index])
            matrix_rho_stepdescent.append(
                self.geo_dstep_descent[site][:, index ])
                
        # build a matrix of all data collected 

        matrix_rhorr= func.concat_array_from_list(
            matrix_rhorr, concat_axis=1)
        matrix_rhoaver= func.concat_array_from_list(
            matrix_rhoaver, concat_axis=1)
        matrix_rho_stepdescent= func.concat_array_from_list(
            matrix_rho_stepdescent, concat_axis=1)
        
        # create empty list to collect infos of read values
        write_averlines , write_rrlines, write_stepdescent_lines ,\
            write_blnfiles =[[]for i in range(4)]
        
        # writes multilines in the same times 
        
        if write_negative_depth  is True : self.geo_depth  *=-1
                # scalled station location and geo_depth 
        self.geo_depth /=df 
        self.station_location /=df
        
        for writes_lines , matrix in zip ( [
                write_averlines , write_rrlines, write_stepdescent_lines], 
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
            filename ='{0}.{1}'.format(os.path.basename(
                filename), datetime.datetime.now().month)
         
        for  ikey, stn in enumerate(self.station_names) : 
            write_blnfiles .append(''.join([
                '{0:<7.6f},'.format(
                self.station_location[ikey]), 
                '{0:<7.6f},'.format(
                    self.elevation[ikey]),
                '{0:<4}'.format(stn), '\n']))
   
        
        #writes files 
        for ii , (tfiles , wfiles) in enumerate(zip(
                ['_aver', '_rr', '_sd', '_yb'], 
                [write_averlines, write_rrlines, 
                 write_stepdescent_lines,  write_blnfiles])): 
            if ii == 3 : mm='bln'
            else :mm='dat'
            with open(''.join([filename,'{0}.{1}'.format(tfiles, mm)]),
                      'w') as fw : 
                fw.writelines(wfiles)
                
        #savefile or create one if not provide
        self.savepath = func.cpath(self.savepath , '_outputGeoSD_' )
        if self.savepath is not None : 
            import shutil 
            try :
                for jj, file in enumerate( [
                        filename + '_aver.dat', filename +'_rr.dat',
                        filename +'_sd.dat', filename + '_yb.bln']):
                    shutil.move(os.path.join(os.getcwd(),file) ,
                        os.path.join(self.savepath , file))
            except : 
                warnings.warn("It seems the files already exists !")
                
        filenames =[filename + '_aver.dat', filename +'_rr.dat',
                                            filename +'_sd.dat',
                                            filename + '_yb.bln']   

        print('---> geo output files {0}, {1}, {2}  & {3}'
              ' have been successfully written to  <{4}>.'.format(
                  *filenames, self.savepath))
           
    def to_oasis_montaj (self, 
                         model_fn=None,
                         iter_fn=None,
                         profile_fn =None, 
                         mesh_fn=None,
                         data_fn=None,
                         filename=None,
                         savepath =None,
                         **kwargs) : 
        """
        write output to oasis montaj when station loocation and profile 
        coordinates are provided . We assune before using this method , you are
        already the coordinates files at disposal(*stn), if not  use the method 
        `to_golden_software`.Coordinated files are Easting Northing value  
        not in degree decimals . It uses occam 2D outputfiles or Bo Yang 
        iter2Dat file also add profile XY coordinates (utm_zone ).
        
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
        etaCoef             int             number to separate the model blocks 
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
                if site ==station_names[0]: 
                    index =0 
                if site ==station_names [-1] : index =2 
                else : index =1  # read th middle array 
                model_rho_matrix.append(geo_dict_rho[site][:, index])
                
            model_rho_matrix= func.concat_array_from_list(model_rho_matrix ,
                                                          concat_axis=1)
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
        etaCoef = kwargs.pop('etaCoef', None)
        
        # write -externals files : Step descent , average rho 

        input_rho =kwargs.pop('input_resistivities', None)
        input_layers =kwargs.pop('input_layers', None)
        step_descent =kwargs.pop('step_descent', None)
        
        if savepath is not None: 
            self.savepath =savepath 
        if input_rho is not None : 
            self.input_resistivities  = input_rho 
        
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
                                  [model_fn, mesh_fn, data_fn, iter_fn,
                                   iter2dat_fn, bln_fn ]):
            if mfile is not None : 
                setattr(self, mattar, mfile)
        
        # --Buil profile object ----------------------------
        if elevation is not None : self.elevation = elevation 
        profile_obj = Profile ()
        profile_obj.read_stnprofile(profile_fn , easting=easting,
                                    northing =northing , 
                                    elevation = self.elevation)
        profile_obj.straighten_profileline ( reajust=scalled_east_north , 
                                            output =output_sprofile)
        
        print('** {0:<37} {1} {2} {3}'.format('Dipole length','=',
                                              profile_obj.dipole_length, 
                                              'm.'  ))
        
        profile_angle, geo_electric_strike = profile_obj.get_profile_angle()
        
        print('** {0:<37} {1} {2} {3}'.format(' Profile angle ','=',
                                              round(profile_angle, 2), 
                                              'deg N.E' ))
        print('** {0:<37} {1} {2} {3}'.format('Geo_electric strike','=',
                                              round(geo_electric_strike , 2), 
                                              'deg N.E' ))
        
        # In the case where elevation is 
        #contain on the profile stn file , get elevation 
        if profile_obj.elev is not None : self.elevation = profile_obj.elev 
        
        # build section Infos 
        #Stations	Easting_X_m	Northing_Y_m	Elev_H_m	
        #x_m	, Norm_h_m	offsets_m	DOI_max_m
        info_list = ['station', 'easting_X_m', 'northing_Y_m', 'elev_H_m',
                     'position_x_m', 'nivelize_h_m', 'offset_m', 'doi']
        #normalH 
        nivelize_elevation = np.around(
            profile_obj.elev - profile_obj.elev.min(), 2)
        # stn position 
        if write_negative_depth is True : doi =-1* self.doi 
        else : doi = self.doi 
        infos_oasis_montaj = func.concat_array_from_list(list_of_array= [
                            self.station_names ,
                             profile_obj.east , 
                             profile_obj.north , 
                             profile_obj.elev , 
                             profile_obj.stn_position , 
                             nivelize_elevation , 
                             self.station_location , 
                             np.full((self.station_location.shape[0],), doi),
                             ], 
                            concat_axis=1)
        #---Build main data for oasis ------------------------------------
        
        # manage the depth info and create matrix of depth 
        spacing_depth =abs(self.geo_depth.max()-self.geo_depth.min())/ (
            len(self.geo_depth)-1)
        print('** {0:<37} {1} {2} {3}'.format('Spacing depth','~=', round(
            spacing_depth,2), 'm.' ))
        
        if normalize_depth is True : 
            # get the step_deth 
            geo_depth = np.around(np.linspace(self.geo_depth [0],
                                              self.geo_depth[-1], 
                                              len(self.geo_depth )))
            spacing_depth =abs(self.geo_depth.max()- self.geo_depth.min())/ (
                len(self.geo_depth)-1)
    
            print('---> Depth normalized !')
            print('---> {0:<37} {1} {2} {3}.'.format('new spacing depth',
                                                     '=', round(spacing_depth),
                                                     'm' ))
        else :
            print('---> UNnormalized  depth!')
        # build a mtrix of depth  
        depth_oasis_montaj  = np.repeat(geo_depth.reshape(geo_depth.shape[0], 1),
                              len(self.station_location), axis =1 ) 
        depth_oasis_montaj = depth_oasis_montaj .T      # Transpose the depth 
        
        if write_negative_depth is True :depth_oasis_montaj *= -1
            
        # build data 
        
        model_oasis_montaj = create_model_matrix (
                                        station_names = self.station_names, 
                                        geo_dict_rho= self.geo_d)
        
        if to_log10: np.log10( model_oasis_montaj )
        
        data_oasis = np.concatenate((infos_oasis_montaj,
                                     depth_oasis_montaj, 
                                     model_oasis_montaj), axis =1)
        # buidd pandas columns 
        depres_pandas_columns = info_list + ['dep_{0}'.format(int(dd)) 
                                             for dd in geo_depth] +\
            ['res_{0}'.format(int(dd)) for dd in geo_depth] 
            
            
        oas_pandas = pd.DataFrame(data=data_oasis,
                                  columns=depres_pandas_columns)
        
        fex=0  # flag to write external files 
        if write_external_files is True :
            if input_layers is not None : self.input_layers = input_layers
            if step_descent is not None : self.step_descent = step_descent 
            if etaCoef is not None : self.etacoeff =etaCoef
            
            if self.input_resistivities is None : 
                mess = ''.join(['! No input resistivities provided! ',
                                'Could not write external geo-files.', 
                                ' If you expect to write external GEO ',
                                '{step_descent  - roughness  and average_rho} ', 
                                ' files, you ABSOLUTELY need  to provided at',
                                ' least a truth resistivity values', 
                                ' of some geological formation of the area.'])
                print(punc.fmt_text(mess,fmt='+'))
                
                warnings.warn(mess)
                self._logging.debug (mess)
            if self.input_resistivities is not None : 
                self.geo_build_strata_logs(
                    input_resistivities= self.input_resistivities, 
                    input_layers=self.input_layers,
                    step_descent=self.step_descent,
                    etaCoef =self.etacoeff)
             
            # now build matrix of all external files 
            model_step_descent = create_model_matrix (
                station_names = self.station_names, 
                geo_dict_rho= self.geo_dstep_descent)   
            model_geo_roughness = create_model_matrix (
                station_names = self.station_names, 
                 geo_dict_rho= self.geo_drr)
            model_average_rho = create_model_matrix (
                station_names = self.station_names, 
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
                with pd.ExcelWriter(''.join([
                        filename + '.main._cor_oas',writeType])) as writer :
                    for  sname , df_ in zip([
                            '.main.', '_sd', '_rr', '_aver'],
                            [oas_pandas, STD_pandas,ROUGH_pandas, AVER_pandas ]):
                        df_.to_excel(writer,sheet_name=sname, index =writeIndex)
                        
                    fex =2 
            else : 
     
                oas_pandas.to_excel(''.join([
                    filename,'.main._cor_oas', writeType]),
                    sheet_name=filename +'main_cor_oas', index=writeIndex)  

            
        elif writeType.lower().find('csv')>= 0 : 
            writeType ='.csv'
            oas_pandas.to_csv(''.join([filename,'.main._cor_oas',writeType]),
                              header=csvsetHeader,index =writeIndex,sep=csvsep) 
                  
            if fex==1 : 
                for fnx, _df in zip( ['_sd', '_rr', '_aver'] ,
                                    [STD_pandas,ROUGH_pandas, AVER_pandas ]):
                    
                    _df.to_csv(''.join([filename + fnx,'_cor_oas', writeType]),
                                header=csvsetHeader,index =writeIndex,sep=csvsep) 
                                         

                    
        else :
            mess = 'The type you provide is wrong.'\
                ' Support only *.xlsx and *.csv format'
            self._logging.error (mess)
            warnings.warn (mess)
            raise GeoError(
                'wrong format ={0} !'\
                ' Could not write geo file to oasis.'\
                    ' Support only *.xlsx or *.csv format ')

        # export to savepath  
        self.savepath = func.cpath(self.savepath , '_output2Oasis_')

        if self.savepath is not None : 
            import shutil 
            try :
                if fex ==2 : 
                   shutil.move(os.path.join(os.getcwd(),''.join(
                       [filename,'.main._cor_oas', writeType])) ,
                   os.path.join(self.savepath , ''.join([
                       filename,'.main._cor_oas', writeType])))
                else : 
                    for jj, file in enumerate(  
                            ['.main.', '_sd', '_rr', '_aver']):
                        if fex ==0 : 
                            if jj ==1 : 
                                break # dont continue, other files dont exists  
                        shutil.move(os.path.join(os.getcwd(),
                                                 ''.join([filename, file, 
                                                '_cor_oas', writeType])) ,
                            os.path.join(self.savepath , 
                                ''.join([filename, file, 
                                         '_cor_oas', writeType])))
            except : 
                warnings.warn("It seems the files already exists !")
        
        # write Infos 
        print('** {0:<37} {1} {2} '.format('number of stations',
                                           '=', len(self.station_names)))
        print('** {0:<37} {1} {2} {3}'.format('minimum offset',
                                              '=', self.station_location.min(), 
                                              'm' ))
        print('** {0:<37} {1} {2} {3}'.format('maximum offset',
                                              '=', self.station_location.max(),
                                              'm' ))
        print('** {0:<37} {1} {2} {3}'.format('maximum depth',
                                              '=', geo_depth.max(), 'm' ))
        print('** {0:<37} {1} {2} {3}'.format('spacing depth ',
                                              '=', round(spacing_depth,2), 
                                              'm' ))
        if self.elevation is None or np.all(self.elevation ==0.) : 
            print('--->  Elevation no added !')
        else : 

            print('** {0:<37} {1} {2} {3}'.format('minumum elevation',
                                                  '=', self.elevation.min(), 
                                                  'm' ))
            print('** {0:<37} {1} {2} {3}'.format('maximum elevation ',
                                                  '=', self.elevation.max(),
                                                  'm' ))
        
        print('** {0:<37} {1} {2} {3}'.format('minumum resistivity value','=',
                                              round(model_oasis_montaj.min(), 2),
                                              'Ω.m' ))
        print('** {0:<37} {1} {2} {3}'.format('maximum resistivity value','=',
                                              round(model_oasis_montaj.max(), 2), 
                                              'Ω.m' ))
   
        print('** {0:<37} {1} {2} '.format('Lowest station','=',
                                           self.station_names[ 
                                               int(np.where(
                                nivelize_elevation==nivelize_elevation.min(
                                    ))[0])] ))
        
        print('** {0:<37} {1} {2}'.format('Highest station','=', 
                    self.station_names[ 
                    int(np.where(
                        nivelize_elevation==nivelize_elevation.max())[0])]))
        print('** {0:<37} {1} {2} {3}'.format('Altitude gap','=',
                                              nivelize_elevation.max(), 'm' ))
        print('** {0:<37} {1} {2}'.format('Number of running ','=', 
                                          len(depres_pandas_columns)
                                          * len(self.station_location)))

        if fex ==0 or fex==2: 
            print('---> geo output file {0}, '
                  ' has been successfully written to  <{1}>.'.
                  format(''.join([filename, '.main._cor_oas',
                                  writeType]), self.savepath))
        elif fex ==1 : #read external files 

            filenames =[ filename + file +'_cor_oas'+ writeType for file ,
                        exten in zip(['.main.', '_sd', '_rr', '_aver'], 
                                     [writeType for i in range(4)])]
                        

            print('---> geo output files {0}, {1}, {2} '
                  ' & {3} have been successfully '\
                  'written to  <{4}>.'.format(*filenames, self.savepath))
                

        print('-'*77)

class Geosurface :
    """
    Read Multidata from oasis montaj output files  generated by `geodrill`
    module. Class to Build  a depth surface map for depth imaging . 
 
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
        
        >>> from pycsamt.geodrill.geocore import Geosurface 
        >>> gs_obj = geosurface (path =os.path.join(os.environ['pyCSAMT'], 
        ...                                   'geodrill', 'data', 
        ...                                   InputOas), 
        ...                            depth_values = [28, 100])
        >>> gs_obj.write_file()
    """
    geo_surface_format=["csv",  "xlsx", "json","html","sql"] 
    
    read_dico = {
                ".csv":pd.read_csv, 
                 ".xlsx":pd.read_excel,
                 ".json":pd.read_json,
                 ".html":pd.read_json,
                 ".sql" : pd.read_sql
                 }  
    
    
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
        if path is not None : 
            self.path =path 
        if self.path is not None : 
            if os.path.isdir(self.path):  # get the list of files in folder 
                # get file and be sure that format exist in 
                self.oasis_data = [os.path.join(self.path, file)
                                   for file in os.listdir(self.path)
                                   if ( file.split('.')[-1] 
                                       in self.geo_surface_format)]
                
                # print(self.oasis_data)
                if self.oasis_data is None or len(self.oasis_data)==0 : 
                    mess ='No files detected!. Please provided'\
                        ' right path to oasis models files.'
                    warnings.warn(mess)
                    self._logging.error(mess)
                    raise GeoError(mess)
      
        elif self.path is None :
            raise GeoError(
                'No path found ! Please provided a right path.')
        
        # get extension file 
        extension_files =[os.path.splitext(pathfile)[1] 
                          for pathfile in self.oasis_data]
        self.extension_file =extension_files
        if self.extension_file.replace('.','') not in self.geo_surface_format : 
            mess = 'Unacceptable format = {0}. Could read format ={1}'.\
                format(self.extension_file, tuple(self.geo_surface_format))
            self._logging.warning(mess)
            raise FileHanglingError(mess)
        
        
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
            
             >>> from pycsamt.geodrill.geocore import Geosurface 
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
                                     name in df.columns 
                                     if name.find('dep_')>=0  ])
    
            return info_names, depth_values , np.abs(depth_values.max() - 
                                                     depth_values.min())/ (
                                                         len(depth_values)-1)
        
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
        # initiliase gloabal dico to takes all files with key as survey lines 
        self.global_dico = {} 
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
            self.__setattr__(oasnames, (tem_array[::, :len_info_names],  
            tem_array[::, len_info_names:len_info_names +len_depth_offset], 
            tem_array[::, len_info_names +len_depth_offset:]) )
           
        # initiliase dictionary to hold  
        self.geos_profile_strike_angles ={}
        for gskeys,gsvalues in self.info_depthvalues.items(): 
            print('** ----- {0} : --|>{1:<37} :'.format('file',gskeys))
            print('** {0:<37} {1} {2} {3}'.format('depth spacing ','=',
                                                  gsvalues[2], 'm' ))
            print('** {0:<37} {1} {2} {3}'.format('maximum depth','=',
                                                  gsvalues[1].max(), 'm' ))
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
                                                      profile_angle, 
                                                      'degrees E of N.' ))
                print('** {0:<37} {1} {2} {3}'.format('geoelectric strike','=',
                                                      gstrike, 
                                                      'degrees E of N.' ))
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
            raise GeoError(
                'NoneType could not be computed. Please provided'
                '  right values !')
        
        if isinstance(self.depth_values, (float, str)): 
            try : 
                self.depth_values = float(self.depth_values)
            except : 
                raise TypeError('Could not convert {0} to float'.\
                                              format(self.depth_values))
            else : self.depth_values =np.array([self.depth_values])
  
        elif isinstance(self.depth_values , (list, tuple, np.ndarray)): 
            self.depth_values = np.array(self.depth_values)
            try : 
                self.depth_values = self.depth_values.astype('float')
            except : raise TypeError(
                    'Could not convert values to float!')

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
            new_depth_value, index =  func.get_closest_value(
                values_range= self.info_depthvalues[site_name][1],
                              input_value= depth_value)
            # now get especially array or depth and resistivy at that value
            # from attribute name  value attr = infomatrix ,
            #depth matrix and res matrix 
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
        if savepath is not None : 
            self.savepath =savepath 

        if path is not None :
            self.path = path 
        if depth_values is not None : 
            self.depth_values =depth_values 
        if fileformat is not None : 
            self.export_format = fileformat.lower() 
        if self.export_format is not None : 
            self.export_format= self.export_format.lower() # for consistency 
            
        if self.export_format.find('exc')>=0 or \
            self.export_format.find('xls')>=0 : 
                self.export_format = 'xlsx'
        elif self.export_format.find('csv')>=0 : 
            self.export_format ='csv'
        else : 
            mess=''.join([
                '-->Sorry! Depth map actually does not support the {0} ',
                    'format provided. Only support `xlsx` or `csv` format.', 
                    ' Please provide the rigth format !'])
            self._logging.error(mess.format(self.export_format))
            warnings.warn(mess.format(self.export_format))
            raise FileHanglingError(
                'Format provided = {0} is wrong ! '
                'Could output only `csv` or `xlsx`.'.
                    format(self.export_format))
            
        
        if self.depth_values is None : 
            warnings.warn(
                'Need to specify the value of depth for imaging !')
            self._logging.warn(
                'Need to specify the value of depth for imaging !')
            raise GeoError(
                'Need to specify the depth value for imaging.')
        if self.path is None : 
            warnings.warn(
                'Need to provide the rigth path `geodrill` model output files.')
            self._logging.warning(
                'Need to provide the rigth path `geodrill` model output files.')
            raise GeoError(
                'No path detected! Please provide a right path to `geodrill`\
                    oasis montaj outputfiles ')
  
        self.get_depth_surfaces() # call methods to create 
        
        # get values from dictionnary and write file 
        #make output filename 
        filenames =[]
        filename =''.join([file.replace('_cor_oas', '') 
                                for file in self.filenames])
        if self.export_format =='xlsx':
            with pd.ExcelWriter(''.join([filename+ 
                                         '_gs{0}.'.format(
                                             datetime.datetime.now().month ),
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
                                    '_gs{0}.'.format(
                                        datetime.datetime.now().month),
                                    self.export_format]),
                              header=csvsetHeader,
                              index =writeIndex,sep=csvsep)
                
                filenames.append(
                    ''.join([filename + str(file) ,
                             '_gs{0}.'.format(datetime.datetime.now().month ), 
                             self.export_format]))
                
        # export to savepath 
        self.savepath = func.cpath(self.savepath , '_outputGS_')
        if self.savepath is not None : 
            import shutil 
            try :
                for nfile in filenames : 
                   shutil.move(os.path.join(
                       os.getcwd(),  nfile),           
                       os.path.join(self.savepath , nfile))
            except : 
                warnings.warn("It seems files already exist !")

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
        _logger.info(
            'Computing  profile angle from Easting and Nothing coordinates.')
        if easting is None or northing is None : 
            raise GeoError(
                'NoneType can not be computed !')
            
            # use the one with the lower standard deviation
        try :
            easting = easting.astype('float')
            northing = northing.astype('float')
        except : 
            raise TypeError(
                'Could not convert input argument to float!')
        
        if stats_import is True : 
            profile1 = spSTAT.linregress(easting, northing)
    
            profile2 =spSTAT.linregress(northing, easting)
        else :
            warnings.warn(
                'Could not find scipy.stats, cannot use method linearRegression '
                  'check installation you can get scipy from scipy.org.')
            _logger.warning(
                'Could not find scipy.stats, cannot use method lineaRegress '
                    'check installation you can get scipy from scipy.org.')
            raise ImportError(
                'Could not find scipy.stats, cannot use method lineaRegress '
                    'check installation you can get scipy from scipy.org.')
            
        profile_line = profile1[:2]
        # if the profile is rather E=E(N),
        # the parameters have to converted  into N=N(E) form:
        
        if profile2[4] < profile1[4]:
            profile_line = (1. / profile2[0], -profile2[1] / profile2[0])

        # if self.profile_angle is None:
        profile_angle = (90 - (np.arctan(profile_line[0]) * 180 / np.pi)) % 180

        # otherwise: # have 90 degree ambiguity in 
        #strike determination# choose strike which offers larger
        #  angle with profile if profile azimuth is in [0,90].
        
        return np.around(
            profile_angle,2) , 'Profile angle is {0:+.2f} degrees E of N'.format(
                profile_angle)
         
        
      
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
                
                _logger.warning(
                    'NoneType found.Could not compute geo-electrike strike!')
                raise GeoError(
                    'NoneType found. Could not compute geo-electrike strike!')
        
        if profile_angle is None : 
            if easting is not None and northing is not None : 
                profile_angle ,_ = geostrike.compute_profile_angle(
                                    easting, northing)
        
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
            geo_electric_strike %=  180         # keep value of
            #geoelectrike strike less than 180 degree
            
        geo_electric_strike =np.floor(geo_electric_strike)
        
        return  geo_electric_strike, profile_angle ,\
            'Profile angle is {0:+.2f} degrees E of N'.format(geo_electric_strike)
        
      
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
                    
        **auto** : bool  
                     option to automatically well data . set to True 
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
        
        >>> from pycsamt.geodrill.geocore import Drill 
        >>> parser_file ='nbleDH.csv'
        >>> drill_obj=Drill(well_filename=os.path.join(os.environ['pyCSAMT'],
        ...                  'data', 'drill_example_files',parser_file),
        ...      build_manually_welldata=False)
        >>>  scollar=drill._collar(DH_Top=None)
        >>> sgeo=drill.dhGeology()
        >>> ssam=drill.dhSample()
        >>> selevaz=drill.dhSurveyElevAz( add_elevation=None, 
        ...                             add_azimuth=None)
        >>> swrite=drill.writeDHData(data2write ="*",
                                 savepath =None)
    """
    try:
        import openpyxl
    except ImportError: 
        func.subprocess_module_installation('openpyxl')
        
    def __init__(self, well_filename=None , auto=True, **kwargs):
        
        self._logging=csamtpylog.get_csamtpy_logger(self.__class__.__name__)
        
        self.wfilename=well_filename
        self.auto=auto
        
        self.mask=kwargs.pop("mask",1)
        self.utm_zone=kwargs.pop("utm_zone","49N")
        self.compute_azimuth=kwargs.pop("compute_azimuth",False)
        self.dip =kwargs.pop("Drill_dip",90)
        self.buttom=kwargs.pop("Drill_buttom", None)
        self.savepath =kwargs.pop('savepath', None )

        
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
        
        if self.auto is False and self.wfilename is None :
            
            self.daTA=func.build_wellData (add_azimuth=self.compute_azimuth, 
                                            utm_zone=self.utm_zone,
                                            report_path = self.savepath, 
                                            )
            self.wdata=self.daTA[1]
            
            self.wdico["DH_East"]   =   self.wdata[:,1]
            self.wdico["DH_North"]  =   self.wdata[:,2]
            self.wdico["DH_Hole"]   =   self.wdata[:,0]
            self.wdico['DH_Dip']    =   self.wdata[:,4]
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
                    it's the Top of data for each Hole Name. 
                    ndaray (number of DH , 1) 
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
                'shape doesnt match. The convenience '\
                    ' shape is %d.'%self.wdata.shape[0]
        
        # print(DH_Top)
        self.wdico.__setitem__('DH_Top',DH_Top)
        
        if self._f == 0 :
            if add_elevation is None :
                #No topography is added , set to 0 
                add_elevation=np.full((len(self.wdico['DH_East']),1),0,
                                      dtype='<U12')
            elif add_elevation is not None :
                if type(add_elevation ) is list :
                    add_elevation =np.array(add_elevation)
                assert add_elevation.shape[0]==\
                    self.wdico['DH_East'].shape[0],"INDEXERROR:"\
                    " The the current dimention of Elevation data is {0}.It's must be"\
                        " the size {1}.".format(
                            add_elevation.shape[0],self.wdico['DH_East'].shape[0])
            
            self.wdico.__setitem__("Elevation", add_elevation)
                    
        elif self._f == 1 :
            
            if add_elevation is not None:
                
                if type(add_elevation ) is list :
                    add_elevation =np.array(add_elevation)
                try :
                    np.concat((add_elevation,self.wdico['DH_East']))
                except Exception : 
                    mess =''.join([
                        'SIZEERROR! Try to set the elevation dimentional as ', 
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
        Method to build geology drillhole log. The name of input rock must
        feell exaction accordinag to a convention AGSO file . If not sure
        for the name of rock and Description and label. You may consult
        the geocode folder before building the well_filename. If the entirely
        rock name is given , program will search on the AGSO file the 
        corresponding Label and code . If the rock name is  founc then 
        it will take its CODE else it will generate exception. 
 
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
        geoelm= GU.get_agso_properties()
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
                    mess=''.join(['The Geological Name ({0})'
                                  ' given in is wrong'.format(elm),
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
                              self.wdico['DH_From'].reshape((
                                  self.wdico['DH_From'].shape[0],1)),
                              self.wdico['DH_To'].reshape((
                                  self.wdico['DH_To'].shape[0],1)),
                              self.wdico['Rock'].reshape((
                                  self.wdico['Rock'].shape[0],1)),
                              dhgeopseudosamp.reshape((
                                  dhgeopseudosamp.shape[0],1)),
                              self.dh_geoleast.reshape((
                                  self.dh_geoleast.shape[0],1)),
                              self.dh_geol_norths.reshape((
                                  self.dh_geol_norths.shape[0],1)),
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
        if wsamp is None :
            self.sampleDHDATA = None 
            return  # mean no geochemistry sample is provided 
        
        self.wdico.__setitem__('DH_From', wsamp[:,1])
        self.wdico.__setitem__('DH_To', wsamp[:,2])
        self.wdico.__setitem__("Sample",wsamp[:,3])
        dhsampseudorock=np.zeros((wsamp.shape[0]))
        
        ###### FIND AGSO MODULE (AGSO_STCODES) #######
        #Try to check the name of sample and their acronym
        
        if path_to_agso_codefile is None:
            path_to_agso_codefile=os.path.join(os.path.abspath('.'),
                                             'pycsamt/geodrill/_geocodes' )
            sampelm= GU.get_agso_properties(
                config_file = os.path.join(path_to_agso_codefile,
                                           'AGSO_STCODES.csv') )
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
                    mess=''.join([
                        'The Sample Name({0}) given in is wrong'.format(elm),
                        'Please provide a right name the right Name.', 
                        'Please consult the AGSO_STCODES.csv file located in ', 
                        '<pycsamt/geodrill/_geocodes> dir. Please keep the'
                        '  header safe and untouchable.'])
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
                              self.wdico['DH_From'].reshape(
                                  (self.wdico['DH_From'].shape[0],1)),
                              self.wdico['DH_To'].reshape(
                                  (self.wdico['DH_To'].shape[0],1)),
                              dhsampseudorock.reshape(
                                  (dhsampseudorock.shape[0],1)),
                              self.wdico['Sample'].reshape(
                                  (self.wdico['Sample'].shape[0],1)),
                              dh_sampeast.reshape(
                                  (dh_sampeast.shape[0],1)),
                              dh_sampnorths.reshape(
                                  (dh_sampnorths.shape[0],1)),
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
        
        # sizep=self.wdico['DH_East'].shape[0]
        if self._f == 0 :
            if add_elevation is None :
                #No topography is added , set to 0 
                add_elevation=np.full((len(self.wdico['DH_East']),1),0,
                                      dtype='<U12')
            elif add_elevation is not None :
                if type(add_elevation ) is list :
                    add_elevation =np.array(add_elevation)
                assert add_elevation.shape[0]==self.wdico[
                    'DH_East'].shape[0],"INDEXERROR:"\
                    " The the current dimention of Elevation data is {0}.It's must be"\
                        " the size {1}.".format(
                            add_elevation.shape[0],self.wdico['DH_East'].shape[0])
            
            self.wdico.__setitem__("Elevation", add_elevation)
                    
        elif self._f == 1 :
            
            if add_elevation is not None:
                
                if type(add_elevation ) is list :
                    add_elevation =np.array(add_elevation)
                try :
                    np.concat((add_elevation,self.wdico['DH_East']))
                except :
                    mess= ''.join([
                        'SIZEERROR! Try to set the elevation dimentional. ', 
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
            self.wdico.__setitem__("DH_RL",np.full(
                (self.daTA[1].shape[0]),0,dtype='<U12'))
        
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
        
        self.wdico.__setitem__("DH_RL",np.full(
            (self.daTA[1].shape[0]),0,dtype='<U12'))
        # add Hole and Depth 
        
        surveyELEV =np.concatenate((self.wdico['DH_Hole'].reshape(
            (self.wdico['DH_Hole'].shape[0],1)),
                                    self.wdico["DH_Bottom"].reshape(
             (self.wdico["DH_Bottom"].shape[0],1))),
                                       axis=1)
        surveyAZIM=np.concatenate((self.wdico['DH_Hole'].reshape(
            (self.wdico['DH_Hole'].shape[0],1)),
                                    self.wdico["DH_Bottom"].reshape(
             (self.wdico["DH_Bottom"].shape[0],1))),
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
                            
        
        self.surveyDHELEV=pd.DataFrame(
            data=surveyELEV, columns=elevazKeys[:5])
        # pop the elevation elm on the list 
        [elevazKeys.pop(ii) for ii, elm in 
         enumerate(elevazKeys) if elm=='Elevation']
        
        self.surveyDHAZIM=pd.DataFrame(data=surveyAZIM, 
                                       columns=elevazKeys)
        
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
    
        savepath =kwargs.pop("savepath",None )
        writeIndex=kwargs.pop('write_index_on_sheet',False)
        writeType =kwargs.pop('writeType', 'xlsx')
        # csvencoding =kwargs.pop('encoding','utf-8')
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
        
        # for df_ in  [df_collar, df_geology, df_sample,
        # df_elevation,df_azimuth]: 
        # df_.set_index(setIndex) # this is unnecessary 
        _dHDico ={'collar': [['1','c'], df_collar],
                 'geology':[['2','g'],df_geology],
                 'sample': [['3','s'],df_sample],
                 'survey_elevation':[['4','elev', 'topo','topography','e'],
                                     df_elevation],
                 'survey_azimuth': [['5','-1','azim','a'],df_azimuth]}
        # skip the sample building  geochemistry doesnt exists
        if self.sampleDHDATA is None :   
            data2write =['1','2','4','5']
          
        if data2write is None or data2write in _all :  # write all 
            with pd.ExcelWriter(''.join([self.daTA[0][:-1],'.xlsx'])) as writer :
                for keys, df_ in _dHDico.items():
                    df_[1].to_excel(writer,sheet_name=keys, index =writeIndex)

                                
        elif data2write is not None :
            
            if type(data2write) is not list:
                data2write=str(data2write)

                try :
                    if writeType in ['xlsx','.xlsx', 'excell',
                                     'Excell','excel','Excel','*.xlsx']:
                        for keys, df in _dHDico.items():
                            if data2write ==keys or data2write.lower(
                                    ) in keys or  data2write in df[0]:
                              df[1].to_excel('.'.join(
                                  [self.daTA[0][:-1],'xlsx']),
                                  sheet_name=keys,index =writeIndex)  

                        
                    elif writeType in ['csv','.csv', 'comma delimited','*.csv',
                                       'comma-separated-value',
                                       'comma seperated value',
                                       'comsepval']:
                        # print('passed')
                        for keys, df_ in _dHDico.items():
                            if data2write == keys or data2write.lower(
                                    ) in keys or data2write in df_[0]:
                              df_[1].to_csv(''.join(
                                  [self.daTA[0][:-1],'.csv']),
                                  header=csvsetHeader,
                                    index =writeIndex,sep=csvsep)  

                except Exception as error :
                    self._logging.error (
                        'The type you provide as WriteType argument is wrong.'
                                ' Support only *.xlsx and *.csv format',error)
                    warnings.warn (
                        'Argument writeType support only [xlsx or csv] format.'
                        ' Must change your *.{0} format'.format(writeType))

                
            elif type(data2write) is list :
                data2write=[str(elm) for elm in data2write] # check the string format
                with pd.ExcelWriter(''.join(
                        [self.daTA[0][:-1],'xlsx'])) as writer :
                    for ii, df in enumerate (data2write):
                        for keys, df__ in _dHDico.items():
                            if df.lower() in keys or df in df__[0] : 
                                df__[1].to_excel(
                                    writer,sheet_name=keys, index =writeIndex)
            else :
                self._logging.error (
                    'The key you provide  as agrument of data2write is wrong. '
                    'the data2write argument should be either [collar, geology,'
                        ' sample, elevation, azimuth] or all (*). ')
                warnings.warn (
                    'Wrong format of input data2write ! Argument dataType is str,'
                    ' or list of string element choosen among [collar, geology,'
                        'sample, elevation, azimuth] or all (*),'
                        ' not {0}'.format(data2write))
 
         # export to savepath 
        if savepath is not None : self.savepath = savepath 
        # create a folder in your current work directory
        if self.savepath is None : 
            try :
                self.savepath  = os.path.join(os.getcwd(), '_outputDH_')
                if not os.path.isdir(self.savepath):
                    os.mkdir(self.savepath)#  mode =0o666)
            except : 
                warnings.warn("It seems the path already exists !")
        
        
        if self.savepath is not None  :
            import shutil
            
            if writeType in ['csv','.csv', 'comma delimited',
                             'comma-separated-value','comma sperated value',
                                       'comsepval']:
                shutil.move ( os.path.join(os.getcwd(),
                                           ''.join(
                                               [self.daTA[0][:-1],'csv'])),
                             self.savepath)
                print('---> Borehole output <{0}> has been written to {1}.'.\
                      format(os.path.basename(
                    ''.join([self.daTA[0][:-1],'.csv'])), self.savepath))
                
            elif writeType in ['xlsx','.xlsx', 'excell','Excell','excel','Excel']:
                try :
                    shutil.move (os.path.join(os.getcwd(),
                                               '.'.join([self.daTA[0][:-1],'xlsx'])),
                                 self.savepath)
                except: 
                    print("--> It seems the destination path "
                          f"{self.savepath} already exists")
                
                print('---> Borehole output <{0}> has been written to {1}.'.\
                      format(os.path.basename(
                      '.'.join([self.daTA[0][:-1],'xlsx'])), self.savepath))
                
                                                                           
@geoplot2d(reason='model', cmap='jet_r', plot_style ='pcolormesh')
def geoModel( **kwargs ):
    """
    Get the type of geoData and plot the model . func will call the decorator 
    `pycsamt.wiewer.plot.geoplot2d`. Decorator can be filled by using matplotlib
    properties. To plot geomodel `misfit `, please 
    change the `reason` argument to `misfit`. 
    
    :see also:: see documentation of decorator `_geoplot2d.__doc__`
    
    :param kwargs: kwargs arguments of "Geodrill obj"  
    :type model_fn: dict  
    
    :param etaCoef: cofficent to divide the model blocks
    :type etaCoef: int 
    
    :param  input_resistivities: Truth values of resistivities 
    :type input_resistivities: array_like or list 
    
    :param plot_misfit: plot misfit of geodata 
    :type plot_misfit: bool

    ==============  =========  ================================================
    Params          Type       Description 
    ==============  =========  ================================================
    geodtype        str         Plot geodata type, can be `rr` for stratigraphy 
                                log, `sd` for roughnness model 
                                and `aver` for rho average model.
                                *default*plot is `rr`.
    iter_fn         str         full path to occam iteration file 
    mesh_fn         str         full path to mesh_fn file 
    data_fn         str         full path to occam_data file 
    doi             str         depth of investigation might 
                                be float or str like "1km" =1000
    depth_scale     str         scale of imaging depth can be 
                                "km" or "m". *Default* is"m"
    step_descent    float       block size for roughning , if not 
                                provided the step will be 20% of D0I
    input_layers    list        true input_layers names : geological 
                                informations of encountered layers             
    ==============  =========  ================================================
                            
    :Example:
        
        >>> from pycsamt.geodrill.geocore import geoModel 
        >>> path=r'F:\ThesisImp\occam2D\invers+files\inver_res\K1'
        >>> inversion_files = {'model_fn':'Occam2DModel', 
        ...                   'mesh_fn': 'Occam2DMesh',
        ...                    "iter_fn":'ITER17.iter',
        ...                   'data_fn':'OccamDataFile.dat'
        ...                    }
        >>> input_resistivity_values =[10, 66, 70, 180, 1000, 2000, 
        ...                               3000, 7000, 15000 ] 
        >>> input_layer_names =['river water', 'fracture zone', 'granite']
        >>> inversion_files = {key:os.path.join(path, vv) for key,
        ...                   vv in inversion_files.items()}
        >>> geoModel(**inversion_files, 
        ...            input_resistivities=input_resistivity_values, 
        ...         input_layers=input_layer_names, geodtype ='rr',
        ...            plot_misfit=True
        ...         )
    
    """
    depth_scale =kwargs.pop('scale', 'm')
    geodtype=kwargs.pop('kind', 'rr')
    plot_misfit = kwargs.pop('plot_misfit', False)

    #create geodrill obj 
    csamtpylog().get_csamtpy_logger().info(
        'Plot Geodatatype = {}'.format(geodtype))
    
    geodrill_obj = Geodrill(**kwargs)
    geodrill_obj.geo_build_strata_logs()

 
    if depth_scale is None : depth_scale= 'm'

    matrix_rhoaver ,matrix_rhorr, \
        matrix_rho_stepdescent, raw_model =[[]for i in range(4)]
    for site in geodrill_obj.station_names : 
        if site ==geodrill_obj.station_names[0]:
            # station id to read is the first index =0
            index =0 
        if site ==geodrill_obj.station_names [-1] :
            index =2 # station id  to read is the 2index 
        else : index =1  # read th middle array 

        matrix_rhorr.append(geodrill_obj.geo_drr[site][:, index])
        matrix_rhoaver.append(geodrill_obj.geo_daver[site][:, index])
        matrix_rho_stepdescent.append(
            geodrill_obj.geo_dstep_descent[site][:, index ])
        raw_model.append(geodrill_obj.geo_d[site][:, index])
        
    # build a matrix of all data collected 
    if geodtype =='sd':
        geo_rho_data= func.concat_array_from_list( 
            matrix_rho_stepdescent, concat_axis=1)
    elif geodtype=='aver':
        geo_rho_data= func.concat_array_from_list(
            matrix_rhoaver, concat_axis=1)
    else : 
         geo_rho_data= func.concat_array_from_list(
             matrix_rhorr, concat_axis=1)
         
    # convert data into log10 values
    # compute model_misfit
    if plot_misfit is True : 
        
        print('You are ploting geoDataType = {} misfit ! '.format(geodtype))
        csamtpylog().get_csamtpy_logger().info(
            'Plot {0} misfit is running !'.format(geodtype))
        raw_model = func.concat_array_from_list(
            raw_model, concat_axis=1)
        
        model_misfit = .01* np.sqrt (
            (raw_model - geo_rho_data)**2 /raw_model**2 ) 
        
        print('{0:-^77}'.format('GeoMisfit info'))
        print('** {0:<37} {1} {2} {3}'.format('Misfit max ','=',
                                              model_misfit.max()*100., '%' ))
        print('** {0:<37} {1} {2} {3}'.format('Misfit min','=',
                                              model_misfit.min()*100., '%' ))
        print('-'*77)
        
        geo_rho_data= model_misfit *100.

    else : 
        geo_rho_data = np.log10(geo_rho_data)


    return (geo_rho_data, geodrill_obj.station_names,
        geodrill_obj.station_location,
        geodrill_obj.geo_depth, geodrill_obj.doi, depth_scale,
        geodrill_obj.model_rms, geodrill_obj.model_roughness , plot_misfit ) 


class GeoStratigraphy(Geodrill):
    """
    Inherit the :class:`pycsamt.geodrill.geoCore.geodrill` to create new model 
    NM using the model get from occam 2D inversion results. 
    
    The challenge of this class  is firstly to delineate with much 
    accuracy the existing layer boundary (top and bottom) and secondly,
    to predict the stratigraphy log before the drilling operations at each 
    station. Moreover, it’s a better way to select the right drilling location
    and also to estimate the thickness of existing layer such as water table 
    layer as well as to figure out the water reservoir rock in the case of 
    groundwater exploration. 
    
    Arguments
    ----------
        **crm** : str,  
                    full path to Occam model file.                             
        **beta** :  int,                
                Value to  divide into the CRM blocks to improve 
                the computation times. default is`5`                               
        **n_epochs** :  int,  
                Number of iterations. default is `100`
        **tres** :  array_like, 
                Truth values of resistivities. Refer to 
                :class:`~.geodrill.Geodrill` for more details
        **ptols**: float,   
                Existing tolerance error between the `tres` values given and 
                the calculated resistivity in `crm` 
        **input_layers** : list or array_like  
                True input_layers names : geological 
                informations of collected in the area.
                
    Hold the attributes from :class:`~.geodrill.Geodrill`
    
    =============  ===========  ===============================================
    Attributes     Type         Explanation 
    =============  ===========  ===============================================
    kind            str         Kind of model function to compute the best fit
                                model to replace the value in `crm` . Can be 
                                'linear' or 'polynomial'. if `polynomial` is 
                                set, specify the `degree. Default is 'linear'. 
    alpha           float       Learning rate for gradient descent computing.  
                                *Default* is ``1e+4`` for linear. If `kind` is 
                                set to `polynomial` the default value should 
                                be `1e-8`. 
    degree          int         Polynomail function degree to implement 
                                gradient descent algorithm. If `kind` is set to 
                                `Polynomial` the default `degree` is 3. 
                                and details sequences 
    nm              ndarray     The NM matrix with the same dimension with 
                                `crm` model blocks. 
    =============  ===========  ===============================================
    
    :Example: 
        
        >>> from pycsamt.geodrill.geocore import Geostratigraphy 
        >>> path=r'F:\ThesisImp\occam2D\invers+files\inver_res\K4'
        >>> inversion_files = {'model_fn':'Occam2DModel', 
                           'mesh_fn': 'Occam2DMesh',
                            "iter_fn":'ITER27.iter',
                           'data_fn':'OccamDataFile.dat'
                            }
        >>> input_resistivity_values =[10, 66, 70, 180, 1000, 2000, 
                                   3000, 50, 7
                                   # 7000,
                                   # 15000 
                                   ] 
        >>> input_resistivity_values =[10, 66, 70, 180, 1000, 2000, 
        ...                               3000, 7000, 15000 ] 
        >>> input_layer_names =['river water', 'fracture zone', 'granite']
        >>> inversion_files = {key:os.path.join(path, vv) for key,
                        vv in inversion_files.items()}
        >>> geosObj = GeoStratigraphy(**inversion_files, 
                             input_resistivities=input_resistivity_values)
        >>> zmodel = geosObj._zmodel
        >>> geosObj._createNM(ptol =0.1)
        >>> geosObj.nm 

    """
    def __init__(self, crm=None, beta=5, ptol=0.1 , n_epochs=100,  **kwargs):
        Geodrill.__init__(self, **kwargs)

        self.crm =crm 
        self._beta =beta 
        self._ptol =ptol 
        self._n_epochs = n_epochs
        
        self._tres = kwargs.pop('tres', None)
        self._eta = kwargs.pop('eta', 1e-4)
        self._kind =kwargs.pop('kind', 'linear')
        self._degree = kwargs.pop('degree', 1)
        self._b = kwargs.pop('build', False)
        
        self.s0 =None 
        self._zmodel =None
        self.nm= None 
        self.z =None
        self.nmSites=None
        self.crmSites=None 
        
        if self.model_res is not None : 
            self.crm = self.model_res 
            self.s0= np.zeros_like(self.model_res)

        for key in list(kwargs.keys()): 
            setattr(self, key, kwargs[key])
            
        if self.input_resistivities is not None: 
            self.tres = self.input_resistivities
     

        if self.crm is not None: 
            self._makeBlock()
            
     
        self.build 

 
     
    @property 
    def n_epochs(self): 
        """ Iteration numbers"""
        return self._n_epochs 
    @n_epochs.setter 
    def n_epochs(self, n_iterations): 
        """ n_epochs must be in integers value and greater than 0"""
        try : 
            self._n_epochs = int(n_iterations)
        except: 
             TypeError('Iteration number must be `integer`') 
        else: 
            if self._n_epochs  <=0: 
                self._logging.debug(
                 " Unaceptable iteration value! Must be a positive value not "
                f"{'a negative value.' if self._n_epochs <0 else 'equal to 0.'}")
                warnings.warn(f" {self._n_epochs} is unaceptable value."
                          " Could be resset to the default value=100.")
                self._n_epochs = 100 
                
    @property 
    def beta (self): 
        """ Block constructor param"""
        return self._beta 
    @beta.setter 
    def beta(self, beta0 ):
        """ Block constructor must be integer value."""
        try : 
            self._beta = int(beta0)
        except Exception: 
            raise TypeError
        else: 
            if self._beta <=0 :
                self._logging.debug(
                    f'{self._beta} is unaceptable. Could resset to 5.')
                warnings.warn(
                    f'`{self._beta}` is unaceptable. Could resset to 5.')
                self._beta= 5
    @property 
    def ptol(self) :
        """ Tolerance parameter """
        return self._ptol 
    @ptol.setter 
    def ptol(self, ptol0): 
        """ Tolerance parameter must be different to zero and includes 
        between 0 and 1"""
        try : 
            self._ptol =float(ptol0)
        except Exception :
            raise TypeError ('Tolerance parameter `ptol` should be '
                             f'a float number not {type (ptol0)!r}.')
        else : 
            if 0 >= self._ptol >1: 
                self._logging.debug(f"Tolerance value `{self._ptol}` is "
                  "{'greater' if self._ptol >1 else 'is unacceptable value'}`."
                    "Could resset to 10%")
                warnings.warn(
                    f'Tolerance value `{self._ptol}` is unacceptable value.'
                    'Could resset to 10%')
                self._ptol = 0.1
    @property 
    def tres(self): 
        """ Input true resistivity"""
        return self._tres 
    @tres.setter 
    def tres(self, ttres):
        """ Convert Tres to log 10 resistivity """
        try : 
            self._tres =[np.log10(t) for t in ttres]
        except : 
            raise ValueError('Unable to convert TRES values') 
        
    @property 
    def build (self): 
        """ Trigger the NM build and return the NM building option """
        
        ntres ='True resistivity values'
        nln ='collected layer names (True)'
        mes =''.join([
            '{0} {1} not defined. Unable to triggered the NM construction. '
            'Please, provide the list/array of {2} of survey area.'])
         
        if self._b:
            if (self.tres and self.input_layers ) is None: 
                warnings.warn(mes.format(
                    'TRES and LN', 'are', ntres +'and'+nln))
                self._b=False 
            elif self.tres is None and self.input_layers is not None: 
                warnings.warn(mes.format('TRES', 'is', ntres))
                self._b=False 
            elif self.input_layers is None and self.tres is not None: 
                warnings.warn(mes.format('LN', 'is', nln))
                self._b=False 
                
            if not self._b:
                self._logging.debug ( "Build is set to TRUE, however,"
                    '{0}'.mes.format(
                        f'{"TRES" if self.tres is None else "LN"}',
                        'is', f'{ntres if self.tres is None else nln}')
                )
                
        if self._b: 
            self._createNM()
        
        #return self._b
   
        
    def _createNM(self, crm =None, beta =5 , ptol= 0.1, **kws): 
        """ Create NM through the differents steps of NM creatings. 
        
        - step 1 : soft minimal computing 
        - step2 : model function computing 
        - step 3: add automatic layers
        - step 4: use ANN to find likehood layers
     
        :param crm: calculated resistivity model blocks 
        :param beta: number of block to build.
        :param ptol: Error tolerance parameters 
  
        """
        def s__auto_rocks (listOfauto_rocks): 
            """ Automatick rocks collected during the step 3
            :param listOfauto_rocks: List of automatic rocks from differents
             subblocks. 
             
            :returns:rocks sanitized and resistivities. 
            """

            listOfauto_rocks= np.concatenate((listOfauto_rocks), axis =1)
            rho_= listOfauto_rocks[1, :]
            rho_=np.array([float(ss) for ss in rho_])
            r_= list(set(listOfauto_rocks[0, :]))
            hres= np.zeros((len(r_), 1))
            h_= []
            for ii, rock  in enumerate(r_): 
                for jj, ro in enumerate(listOfauto_rocks[0, :]): 
                    if rock == ro: 
                        h_.append(rho_[jj])
                m_= np.array(h_)
                hres[ii]= m_.mean()
                h_=[]
            return r_, hres 
        
        iln =kws.pop('input_layers', None)
        tres =kws.pop('tres', None)
        subblocks =kws.pop('subblocks', None)
        disp= kws.pop('display_infos', True)
        n_epochs = kws.pop('n_epochs', None)
        hinfos =kws.pop('headerinfos',
                        ' Layers [auto=automatic]')
        if subblocks is not None: 
            self.subblocks = subblocks
        
        if iln is not None: 
            self.input_layers = iln 
        if tres is not None: 
            self.tres = tres 
        
        if crm is not None:
            self.crm = crm 
        if beta is not None: 
            self.beta = beta 
        if ptol is not None:
            self.ptol = ptol 
        if n_epochs is not None: 
            self.n_epochs = n_epochs 
            
        self.s0 , errors=[], []
        #step1 : SOFMINERROR 
        if itqdm : 
            pbar =tqdm.tqdm(total= 3,
                             ascii=True,unit='B',
                             desc ='WEgeophysics-pycsamt[NM construction]', 
                             ncols =77)
            
        for ii in range(len(self.subblocks)):
            s1, error = self._softMinError(subblocks= self.subblocks[ii])
            self.s0.append(s1)
            errors.append(error)
            if itqdm: pbar.update(1)
        #step2 : MODELFUNCTION USING DESCENT GRADIENT 
        for ii in range(len(self.s0)):
            if 0 in self.s0[ii][:, :]: 
                s2, error = self._hardMinError(subblocks =self.subblocks[ii], 
                                            s0= self.s0[ii])
                self.s0[ii]=s2
                errors[ii]= error 
            if itqdm: pbar.update(2)
        arp_=[]
        #Step 3: USING DATABASE 
        for ii in range(len(self.s0)):
            if 0 in self.s0[ii][:, :]: 
                s3, autorock_properties= self._createAutoLayer(
                    subblocks=self.subblocks[ii], s0=self.s0[ii]  )
                arp_.append(autorock_properties)
                self.s0[ii]=s3       
        # Assembly the blocks 
        self.nm = np.concatenate((self.s0))
        self.z=self.nm[:, 0]
        self.nm = self.nm[:, 1:]
        
        if itqdm: 
            pbar.update(3)
            print(' process completed')
            pbar.close()
            
        # make site blocks 
        self.nmSites= makeBlockSites(x_nodes=self.model_x_nodes, 
                        station_location= self.station_location, 
                             block_model=self.nm )
        self.crmSites = makeBlockSites(x_nodes=self.model_x_nodes, 
                            station_location= self.station_location, 
                             block_model=self.model_res)

        #Update TRES and LN 
        gammaL, gammarho = s__auto_rocks(arp_) 
        
        if self.input_layers is not None: 
            print_layers = self.input_layers  + [ ' {0} (auto)'.format(l) 
                                                 for l in gammaL ]
            self.input_layers = self.input_layers + gammaL
        # keep the auto_layer found     
        self.auto_layers =gammaL
        self.tres = list(np.power(10,self._tres))  + list (np.power(10, 
                  np.array([float(rv) for rv in gammarho])))
        # display infos 
        if disp:
            display_infos(infos=print_layers,
                          header= hinfos)
        #STEP 4: Train ANN: see pycsamt.geodrill.ml.py to predict your 
        #layer: No need to plot the NM 
        
        # copy main attributes for pseudostratigraphic plot purpose 
        import copy 
        for name , attrval in zip(['TRES', 'LNS'], 
                              [self.tres , self.input_layers]):
            setattr(self, name, copy.deepcopy(attrval))
        # memorize data 
        _ps_memory_management(self)
        
        return self.nm

        
    def _softMinError(self, subblocks=None, **kws ): 
        """
        Replace the calculated resistivity by the true resistivity 
        using the soft minimal error (ξ)
        
        :param crm: Is the calculated resistivity model from Occam2D 
        inversion results 
        
        """
        buffer =self.ptol +1  #bufferr error  
        _z = subblocks[:, 0]
        subblocks = subblocks[:, 1:]
        # Hold the columns of depth values 
        s0 = np.zeros_like(subblocks.T)
        error =[]
        for ii in range(subblocks.shape[1]): # hnodes N
            for jj in range(subblocks.shape[0]): # znodes V
                for k in range(len(self.tres)) :
                   sfme_k = (subblocks.T[ii, jj]-self.tres[k])**2\
                       /subblocks.T[ii, jj] **2
                   error.append(sfme_k )
                   if sfme_k  <= self.ptol : 
                       if sfme_k  < buffer : # keep the best minimum 
                           buffer = sfme_k  
                           s0[ii, jj] = self.tres[k]
                           
                buffer = self.ptol +1      # initilize buffer 

        s0= np.concatenate((_z.reshape(_z.shape[0], 1), s0.T), axis =1)
        return s0, error 


    def _hardMinError(self, tres=None, subblocks=None, s0=None, ptol = None,
                         kind='linear', **kwargs ): 
        """The second step introduces the model function F=W∙Z  where W
        contains the weights of parameters number and Z is V×2 matrix 
        that contains a “bias” column. If the parameter number P equal to two, 
        the model function f(z)=∑_(p=1)^P▒〖w_(p-1) z^(p-1) 〗   becomes a
        linear function with 〖f_1〗^((1) ) (z)=  wz+r_0  with w_1=w and w_0=r_0
        he gradient descent algorithm  is used to find the best parameters w
        and r_0  that  minimizes the  MSE loss function  J .
        
        :param subblocks: `crm` block  
        :param s0: blocks from the first step :meth:`~._sofminError`
        :param kind: Type of model function to apply. Can also be 
                a `polynomial` by specifying the `degree` 
                into argument `degree`.
        :Example: 
            
            >>> from pycsamt.geodrill.geocore import GeoStratigraphy
            >>> geosObj = GeoStratigraphy(**inversion_files, 
                              input_resistivities=input_resistivity_values) 
            >>> ss0, error = geosObj._hardMinError(subblocks=geosObj.subblocks[0],
                                     s0=geosObj.s0[0])
        """
        
        if tres is not None:
            self.tres = tres 
        if ptol is not None: 
            self.ptol = ptol 
        
        eta = kwargs.pop('eta', None)
        if eta is not None: 
            self._eta = eta 
        n_epochs =kwargs.pop('n_epochs', None)
        if n_epochs is not None: 
            self.n_epochs = n_epochs 
        kind = kwargs.pop('kind', None)
        if kind is not None:
            self._kind = kind 
        degree = kwargs.pop('degree', None) 
        if degree is not None: 
            self._degree = degree 
        
        buffer =self.ptol +1  #bufferr error 
        _z= s0[:, 0]
        s0 = s0[:, 1:].T

        subblocks=subblocks[:, 1:].T
        error =[]
        for ii in range(s0.shape[0]): # hnodes N
            F, *_= self.gradient_descent(z=_z,s=subblocks[ii,:],
                                         alpha= self._eta,
                                         n_epochs= self.n_epochs, 
                                         kind= self._kind)
            for jj in range(s0.shape[1]): # znodes V
                 if s0[ii, jj] ==0. : 
                    rp =F[jj]
                    for k in range(len(self.tres)) :
                        with np.errstate(all='ignore'): 
                            sfme_k = (rp -self.tres[k])**2\
                                 /rp**2
                            _ermin = abs(rp-subblocks[ii, jj])
                        error.append(sfme_k)
                        if sfme_k <= self.ptol and _ermin<= self.ptol: 
                             if sfme_k  < buffer : # keep the best minimum 
                                buffer = sfme_k  
                                s0[ii, jj]= self.tres[k]
                               
                    buffer = self.ptol +1      # initialize buffer 
                    
        s0= np.concatenate((_z.reshape(_z.shape[0], 1), s0.T), axis =1)    
        return s0, error 
        
    @staticmethod    
    def gradient_descent(z, s, alpha, n_epochs, **kws): 
        """ Gradient descent algorithm to  fit the best model parameter. 
        
        :param z: vertical nodes containing the values of depth V
        :param s: vertical vector containin the resistivity values 
        :param alpha: step descent parameter or learning rate. 
                    *Default* is ``0.01`
        :param n_epochs: number of iterations. *Default* is ``100``
                        Can be changed to other values
        :returns:
            - `F`: New model values with the best `W` parameters found.
            - `W`: vector containing the parameters fits 
            - `cost_history`: Containing the error at each Itiretaions. 
            
        :Example:
            
            >>> z= np.array([0, 6, 13, 20, 29 ,39, 49, 59, 69, 89, 109, 129, 
                             149, 179])
            >>> res= np.array( [1.59268,1.59268,2.64917,3.30592,3.76168,
                                4.09031,4.33606, 4.53951,4.71819,4.90838,
                  5.01096,5.0536,5.0655,5.06767])
            >>> fz, weights, cost_history = gradient_descent(z=z, s=res,
                                                 n_epochs=10,
                                                 alpha=1e-8,
                                                 degree=2)
            >>> import matplotlib.pyplot as plt 
            >>> plt.scatter (z, res)
            >>> plt.plot(z, fz)
        """
        kind_=kws.pop('kind', 'linear')
        kind_degree = kws.pop('degree', 1)
        
        if kind_degree >1 : kind_='poly'
        
        if kind_.lower() =='linear': 
            kind_degree = 1 
        elif kind_.lower().find('poly')>=0 : 
            if kind_degree <=1 :
                _logger.debug(
                    'The model function is set to `Polynomial`. '
                    'The degree must be greater than 1. Degree wil reset to 2.')
                warnings.warn('Polynomial degree must be greater than 1.'
                              'Value is ressetting to `2`.')
                kind_degree = 2
            try : 
                kind_degree= int(kind_degree)
            except Exception :
                raise ValueError(f'Could not `{kind_degree}` convert to integer.')
                
        
        def kindOfModel(degree, x, y) :
            """ Generate kind of model. If degree is``1`` The linear subset 
             function will use. If `degree` is greater than 2,  Matrix will 
             generate using the polynomail function.
             
            :param x: X values must be the vertical nodes values 
            :param y: S values must be the resistivity of subblocks at node x 
            
             """
            c= []
            deg = degree 
            w = np.zeros((degree+1, 1)) # initialize weights 
            
            def init_weights (x, y): 
                """ Init weights by calculating the scope of the function along 
                 the vertical nodes axis for each columns. """
                for j in range(x.shape[1]-1): 
                    a= (y.max()-y.min())/(x[:, j].max()-x[:, j].min())
                    w[j]=a
                w[-1] = y.mean()
                return w   # return weights 
        
            for i in range(degree):
                c.append(x ** deg)
                deg= deg -1 
        
            if len(c)> 1: 
                x= func.concat_array_from_list(c, concat_axis=1)
                x= np.concatenate((x, np.ones((x.shape[0], 1))), axis =1)
        
            else: x= np.vstack((x, np.ones(x.shape))).T # initialize z to V*2
        
            w= init_weights(x=x, y=y)
            return x, w  # Return the matrix x and the weights vector w 
        
        
        def model(Z, W): 
            """ Model function F= Z.W where `Z` id composed of vertical nodes 
            values and `bias` columns and `W` is weights numbers."""
            return Z.dot(W)
        
        # generate function with degree 
        Z, W = kindOfModel(degree=kind_degree,  x=z, y=s)
        
        # Compute the gradient descent 
        cost_history = np.zeros(n_epochs)
        s=s.reshape((s.shape[0], 1))
        
        for ii in range(n_epochs): 
            with np.errstate(all='ignore'): # rather than divide='warn'
                #https://numpy.org/devdocs/reference/generated/numpy.errstate.html
                W= W - (Z.T.dot(Z.dot(W)-s)/ Z.shape[0]) * alpha 
                cost_history[ii]= (1/ 2* Z.shape[0]) * np.sum((Z.dot(W) -s)**2)
            
        F= model(Z=Z, W=W)     # generate the new model with the best weights 
                 
        return F,W, cost_history
     
    @classmethod 
    def geoArgumentsParser(cls, config_file =None): 
        """ Read and parse the `GeoStratigraphy` arguments files from 
        the config [JSON|YAML] file.
        :param config_file: configuration file. Can be [JSON|YAML]
        
        :Example: 
            >>> GeoStratigraphy.geoArgumentsParser(
                'data/saveJSON/cj.data.json')
            >>> GeoStratigraphy.geoArgumentsParser(
                'data/saveYAML/cy.data.yml')
        """
        if config_file.endswith('json'): 
            args = BS.parse_json(config_file)
        elif config_file.endswith('yaml') or config_file.endswith('yml'):
            args =BS.parse_yaml(config_file)
        else: 
            raise ValueError('Can only parse JSON and YAML data.')
        
        return cls(**args)
            
        
    def _makeBlock (self): 
        """ Construct the differnt block  based on `beta` param. Separate blocks 
        from number of vertical nodes generated by the first `beta` value applied 
        to the `crm`."""

        self.zmodel = np.concatenate((self.geo_depth.reshape(
            self.geo_depth.shape[0], 1),  self.model_res), axis =1) 
                                    
        vv = self.zmodel[-1, 0] / self.beta 
        for ii, nodev in enumerate(self.zmodel[:, 0]): 
            if nodev >= vv: 
                npts = ii       # collect number of points got.
                break 
        self._subblocks =[]
        
        bp, jj =npts, 0
        if len(self.zmodel[:, 0]) <= npts: 
            self._subblocks.append(self.zmodel)
        else: 
            for ii , row in enumerate(self.zmodel) : 
                if ii == bp: 
                    _tp = self.zmodel[jj:ii, :]
                    self._subblocks.append(_tp )
                    bp +=npts
                    jj=ii
                    
                if len(self.zmodel[jj:, 0])<= npts: 
                    self._subblocks.append(self.zmodel[jj:, :])
                    break 
                
        return self._subblocks 
 
    @property 
    def subblocks(self): 
        """ Model subblocks divised by `beta`"""
        return self._subblocks 
    
    @subblocks.setter 
    def subblocks(self, subblks):
        """ keep subblocks as :class:`~GeoStratigraphy` property"""
        
        self._subblocks = subblks 
        
    def _createAutoLayer(self, subblocks=None, s0=None,
                          ptol = None,**kws):
        """ 
        The third step of replacement using the geological database. 
        
        The third step consists to find the rock  γ_L in the Γ with the 
         ceiled mean value γ_ρ  in E_props column is close to the calculated 
        resistivity r_11. Once the rock γ_L  is found,the calculated 
        resistivity r_11 is replaced by γ_ρ. Therefore, the rock γ_L is
         considered as an automatic layer. At the same time,the TRES and LN
         is updated by adding   GeoStratigraphy_ρ  and  γ_L respectively to 
         the existing given data. 
         
        """
        db_properties = kws.pop('properties',['electrical_props', 
                                              '__description'] )
        tres = kws.pop('tres', None)
        disp = kws.pop('display_infos', False)
        hinfos = kws.pop('header', 'Automatic layers')
        
        if tres is not None :
            self.tres = tres 
        if ptol is not None: 
            self.ptol = ptol 
        
        def _findGeostructures(_res): 
            """ Find the layer from database and keep the ceiled value of 
            `_res` calculated resistivities"""
            
            structures = self.get_structure(_res)
            if len(structures) !=0 or structures is not None:
                if structures[0].find('/')>=0 : 
                    ln = structures[0].split('/')[0].lower() 
                else: ln = structures[0].lower()
                return ln, _res
            else: 
                valEpropsNames = self._getProperties(db_properties)
                indeprops = db_properties.index('electrical_props')
                for ii, elecp_value  in enumerate(valEpropsNames[indeprops]): 
                    if elecp_value ==0.: continue 
                    elif elecp_value !=0 : 
                        try : 
                            iter(elecp_value)
                        except : pass 
                        else : 
                            if  min(elecp_value)<= _res<= max(elecp_value):
                                ln= valEpropsNames[indeprops][ii]
                                return ln, _res
                    
        def _normalizeAutoresvalues(listOfstructures,listOfvalues):                            
            """ Find the different structures that exist and
            harmonize value. and return an array of originated values and 
            the harmonize values and the number of automatics layer found as 
            well as their harmonized resistivity values. 
            """
            autolayers = list(set(listOfstructures))
            hvalues= np.zeros((len(autolayers,)))
            
            temp=[]
            for ii , autol in enumerate(autolayers): 
                for jj, _alay in enumerate(listOfstructures):
                    if _alay ==autol: 
                        temp.append(listOfvalues[jj])
                hvalues[ii]= np.array(list(set(temp))).mean()
                temp=[]
            
            # build values array containes the res and the harmonize values 
            h= np.zeros((len(listOfvalues),))
            for ii, (name, values) in enumerate(zip (listOfstructures,
                                     listOfvalues)):
                for jj, hnames in enumerate(autolayers) : 
                    if name == hnames: 
                        h[ii]= hvalues[jj]
            
            finalres= np.vstack((np.array(listOfvalues),h) )
            finalln = np.vstack((np.array(autolayers), hvalues))
            return  finalln, finalres 
            
        _z= s0[:, 0]
        s0 = s0[:, 1:].T
        _temptres , _templn =[], []
        subblocks=subblocks[:, 1:].T

        for ii in range(s0.shape[0]): # hnodes N
            for jj in range(s0.shape[1]): # znodes V
                if s0[ii, jj] ==0. : 
                    lnames, lcres =_findGeostructures(
                        np.power(10, subblocks[ii, jj]))
                    _temptres.append(np.log10(lcres))
                    _templn.append(lnames)

        auto_rocks_names_res, automatics_resistivities =\
            _normalizeAutoresvalues(_templn,_temptres )
        
        for ii in range(s0.shape[0]): # hnodes N
           for jj in range(s0.shape[1]): # znodes V
               if s0[ii, jj] ==0. :
                   for k in range(automatics_resistivities.shape[1]): 
                       subblocks[ii, jj] == automatics_resistivities[0,:][k]
                       s0[ii, jj]= automatics_resistivities[1,:][k]
                       break 
        
                                   
        s0= np.concatenate((_z.reshape(_z.shape[0], 1), s0.T), axis =1) 
        
        # display infos 
        if disp:
            display_infos(infos=self.input_layers,
                          header= hinfos)
        
        return  s0, auto_rocks_names_res
        
    @staticmethod
    def _getProperties(properties =['electrical_props', '__description'], 
                       sproperty ='electrical_props'): 
        """ Connect database and retrieve the 'Eprops'columns and 'LayerNames'
        
        :param properties: DataBase columns.
        :param sproperty : property to sanitize. Mainly used for the properties
            in database composed of double parenthesis. Property value 
            should be removed and converted to tuple of float values.
        :returns:
            - `_gammaVal`: the `properties` values put on list. 
                The order of the retrieved values is function of 
                the `properties` disposal.
        """
        def _fs (v): 
            """ Sanitize value and put on list 
            :param v: value 
            :Example:
                
                >>> _fs('(416.9, 100000.0)'))
                ...[416.9, 100000.0]
            """
            try : 
                v = float(v)
            except : 
                v = tuple([float (ss) for ss in 
                         v.replace('(', '').replace(')', '').split(',')])
            return v
        # connect to geodataBase 
        try : 
            _dbObj = GeoDataBase()
        except: 
            _logger.debug('Connection to database failed!')
        else:
            _gammaVal = _dbObj._retreive_databasecolumns(properties)
            if sproperty in properties: 
                indexEprops = properties.index(sproperty )
                try:
                    _gammaVal [indexEprops] = list(map(lambda x:_fs(x),
                                                   _gammaVal[indexEprops]))
                except TypeError:
                    _gammaVal= list(map(lambda x:_fs(x),
                                         _gammaVal))
        return _gammaVal 
       
    @geoplot2d(reason='model',cmap='jet_r', plot_style ='pcolormesh',
               show_grid=False )
    def strataModel(self, kind ='nm', **kwargs): 
        """ 
        Visualize the   `strataModel` after `nm` creating using decorator from 
        :class:'~.geoplot2d'. 
        
        :param kind: can be : 
            - `nm` mean new model plots after inputs the `tres`
            - `crm` means calculated resistivity from occam model blocks 
            *default* is `nm`.
        :param plot_misft:  Set to ``True`` if you want to visualise the error 
            between the `nm` and `crm`. 
        :param scale: Can be ``m`` or ``km`` for plot scale 
        :param in_percent`: Set to `True` to see your plot map scaled in %.
        
        :Example: 
            
            >>> from pycsamt.geodrill.geocore import Geostratigraphy 
            >>> geosObj = GeoStratigraphy(**inversion_files,
                                  input_resistivities=input_resistivity_values, 
                                  input_layers=input_layer_names)
            >>> geosObj.strataModel(kind='nm', misfit_G =False)
        """
        m_='pycsamt.geodrill.geocore.GeoStratigraphy.strataModel'
        def compute_misfit(rawb, newb, percent=True): 
            """ Compute misfit with calculated block and new model block """
            m_misfit = .01* np.sqrt (
                (rawb - newb)**2 /rawb**2 ) 
            if percent is True: 
                m_misfit= m_misfit *100.
            return m_misfit 
        
        if isinstance(kind, bool): 
            kind ='nm' if kind else 'crm'

        depth_scale = kwargs.pop('scale', 'm')
        misfit_G =kwargs.pop('misfit_G', False)
        misfit_percentage = kwargs.pop('in_percent', True)
        
        kind = GU._assert_model_type(kind)
        if self.nm is None: 
            self._createNM()  
                
        if kind =='nm':
            data = self.nmSites 
        if kind =='crm': 
            data = self.crmSites

      # compute model_misfit
        if misfit_G is True : 
            if kind =='crm': 
                warnings.warn("Use `pycsamt.modeling.occam2d.getMisfit` "
                              "decorated function to visualize occam2d misfit"
                              "  model. By default, the plot should be the"
                              " stratigraphic misfit<misfit_G>.")
    
            self._logging.info('Visualize the stratigraphic misfit.')
            data = compute_misfit(rawb=self.crmSites , 
                                  newb= self.nmSites, 
                                  percent = misfit_percentage)
            
            print('{0:-^77}'.format('StrataMisfit info'))
            print('** {0:<37} {1} {2} {3}'.format(
                'Misfit max ','=',data.max()*100., '%' ))                      
            print('** {0:<37} {1} {2} {3}'.format(
                'Misfit min','=',data.min()*100., '%' ))                          
            print('-'*77)
            
        warnings.warn(
            f'Data stored from {m_!r} should be moved on binary drive and'
            ' method arguments should be keywordly only.', FutureWarning)
        
        return (data, self.station_names, self.station_location,
            self.geo_depth, self.doi, depth_scale, self.model_rms, 
            self.model_roughness, misfit_G ) 
    
    def stratigraphyModel (self, kind ='nm', misfit_G= False, **kwargs): 
        """ Make stratigraphy model 
        
        :param kind: 
            - `nms` mean new model plots after inputs the `tres`
            - `crm` means calculated resistivity from occam model blocks 
        :param misfit_G: bool, 
            compute error between CRM and NM if set to `True` 
        """
        # get attribute from Geostratigraph object `
        for file in ['model_fn', 'iter_fn', 'data_fn', 'mesh_fn',
                     'input_resistivities', 'input_layers']: 
            if hasattr(self, file):
                kwargs[file]= getattr(self, file)
                  
        geoModel( kind =kind,
                plot_misfit=misfit_G,
              **kwargs) 
        
        warnings.warn(
            'Please use `pycsamt.geodrill.geocore.GeoStratigraphy.strataModel`'
            'instead of `~stratigraphyModel`.', category =DeprecationWarning) 
    
    @staticmethod
    def _strataPropertiesOfSite(obj, station =None, display_s=True): 
        """ Build all properties of strata under each station.  
        
        Parameters
        ----------
            station: str or int
                Use normal count to identify the number of site to plot or use 
                the name of station preceed of letter `S`. For instance site 
                1 matches the station `S00` litterally
            display_s:bool
                Display the log layer infos as well as layers thicknesses
        
        Examples
        --------
            >>> from pycsamt.geodrill.geocore import GeoStratigraphy 
            >>> import pycsamt.utils.geo_utils as GU 
            >>> geosObj = GeoStratigraphy( input_resistivities=GU.TRES, 
            ...              input_layers=GU.LNS,**GU.INVERS_KWS)
            >>> geosObj._strataPropertiesOfSite(geosObj,station= 'S05')
        """
        
        def stamping_ignored_rocks(fittedrocks, lns): 
            """ Stamping the pseudo rocks and ignored them during plot."""
            ir= set(fittedrocks).difference(set(lns))
            for k in range( len(fittedrocks)): 
                if fittedrocks[k] in ir : 
                    fittedrocks[k] ='$(i)$'
            return fittedrocks
        
        # assert the station, get it appropriate index and take the tres 
        # at that index 
        if station is None: 
            stns = ["S{:02}".format(i) for i in range(obj.nmSites.shape[1])]
            obj._logging.error('None station is found. Please select one station'
                                f' between {func.smart_format(stns)}')
            warnings.warn("NoneType can not be read as station name."
                            " Please provide your station name. list of sites" 
                             f" are {func.smart_format(stns)}")
            raise ValueError("NoneType can not be read as station name."
                             " Please provide your station name.")
        
        if  obj.nmSites is None: 
                obj._createNM()
        
        try : 
            id0= int(station.lower().replace('s', ''))
        except : 
            id_ =GU.assert_station(id= station, nm = obj.nmSites)
            station_ = 'S{0:02}'.format(id_)
        else : 
            id_ =GU.assert_station(id= id0 + 1, nm = obj.nmSites)
            station_ = 'S{0:02}'.format(id0)
        obj.logS = obj.nmSites[:, id_] 
     
        # assert the input  given layers and tres 
        is_the_same_length, msg  = GU.assert_len_lns_tres(
            obj.LNS , obj.TRES)
        
        if is_the_same_length :
            # then finf 
            pslns = obj.LNS
            pstres= obj.TRES 
            ps_lnstres   = [(a, b) for a , b 
                            in zip(obj.LNS, obj.TRES)] 
        if not is_the_same_length:
            # find the pseudoTRES and LNS for unknowrocks or layers
            msg +=  "Unknow layers should be ignored."   
            _logger.debug(msg)
            warnings.warn(msg)
            
            pslns, pstres, ps_lnstres =  fit_tres(
                                            obj.LNS, obj.TRES, 
                                            obj.auto_layers)
        # now build the fitting rocks 
        fitted_rocks =GU.fit_rocks(logS_array= obj.nmSites[:,id_],
                                   lns_=pslns , tres_=pstres)
        # set the raws fitted rocks 
        import copy
        setattr(obj, 'fitted_rocks_r', copy.deepcopy(fitted_rocks) )
        # change the pseudo-rocks  located in fitted rocks par ignored $i$
        # and get  the stamped rocks 
        setattr(obj, 'fitted_rocks',stamping_ignored_rocks(
            fitted_rocks, obj.LNS ) )
        
        # fit stratum property 
        sg, _, zg, _= GU.fit_stratum_property (obj.fitted_rocks,
                                obj.z, obj.logS)
        obj.log_thicknesses, obj.log_layers,\
            obj.coverall = GU.get_s_thicknesses( 
            zg, sg,display_s= display_s, station = station_ )
            
        # set the dfault layers properties hatch and colors from database
        obj.hatch , obj.color =  fit_default_layer_properties(
            obj.log_layers) 
        
        return obj
    
    @staticmethod
    def plotPseudostratigraphic(station, zoom=None, annotate_kws=None, **kws):
        """ Build the pseudostratigraphic log. 
        :param station: station to visualize the plot.
        :param zoom: float  represented as visualization ratio
            ex: 0.25 --> 25% view from top =0.to 0.25* investigation depth 
            or a list composed of list [top, bottom].
        
        :Example: 
            >>> input_resistivity_values =[10, 66,  700, 1000, 1500,  2000, 
                                    3000, 7000, 15000 ] 
            >>> input_layer_names =['river water', 'fracture zone', 'granite']
            # Run it to create your model block alfter that you can only use 
            #  `plotPseudostratigraphic` only
            # >>> obj= quick_read_geomodel(lns = input_layer_names, 
            #                             tres = input_resistivity_values)
            >>> plotPseudostratigraphic(station ='S00')
        
        """
        
        if annotate_kws is None: annotate_kws = {'fontsize':12}
        if not isinstance(annotate_kws, dict):
            annotate_kws=dict()
        obj = _ps_memory_management(option ='get' )  
        obj = GeoStratigraphy._strataPropertiesOfSite (obj,station=station,
                                                       **kws )
        # plot the logs with attributes 
        GU.pseudostratigraphic_log (obj.log_thicknesses, obj.log_layers,station,
                                    hatch =obj.hatch ,zoom=zoom,
                                    color =obj.color, **annotate_kws)
        GU.print_running_line_prop(obj)
        
        return obj 
    
def _ps_memory_management(obj=None, option='set'): 
    """ Manage the running times for stratigraphic model construction.
    
    The script allows to avoid running several times the GeoStratigraphy model
    construction to retrieve the pseudostratigraphic (PS) log at each station.  
    It memorizes the model data for the first run and used it when calling it
    to  visualize the PS log at each station. Be aware to edit this script.
    """
    memory, memorypath='__memory.pkl', 'pycsamt/geodrill/_geocodes'
    mkeys= ('set', 'get', 'recover', 'fetch', set)
    if option not in mkeys: 
        raise ValueError('Wrong `option` argument. Acceptable '
                         f'values are  {func.smart_format(mkeys[:-1])}.')
        
    if option in ('set', set): 
        if obj is None: 
            raise TypeError('NoneType object can not be set.') 
        psobj_token = __build_ps__token(obj)
        data = (psobj_token, list(obj.__dict__.items()))
        BS.serialize_data ( data, memory, savepath= memorypath )

        return 
    
    elif option in ('get', 'recover', 'fetch'): 
        memory_exists =  os.path.isfile(os.path.join(memorypath, memory))
        if not memory_exists: 
            _logger.error('No memory found. Run the GeoStratigraphy class '
                          'beforehand to create your first model.')
            warnings.warn("No memory found. You need to build your "
                          " GeoStratigraphy model by running the class first.")
            raise GeoMemoryError(
                "Memory not found. Use the GeoStratigraphy class to "
                "create your model first.")
        psobj_token, data_ = BS.load_serialized_data(
            os.path.join(memorypath, memory))
        data = dict(data_)
        # create PseudoStratigraphicObj from metaclass and inherits from 
        # dictattributes of GeoStratigraphy class
        psobj = type ('PseudoStratigraphic', (), { 
            k:v for k, v in data.items()})
        psobj.__token = psobj_token
        
        return psobj
                                
        
def makeBlockSites(station_location, x_nodes, block_model): 
    """ Build block that contains only the station locations values
    
    :param station_location: array of stations locations. Must be  
                self contains on the horizontal nodes (x_nodes)
    :param x_nodes: Number of nodes in horizontal 
    :param block_model: Resistivity blocks model 
    
    :return: 
        - `stationblocks`: Block that contains only the
        station location values.
        
    :Example:
        
        >>> from pycsamt.geodrill.geocore import makeBlockSite
        >>> mainblocks= get_location_value(
            station_location=geosObj.makeBlockSite,
             x_nodes=geosObj.model_x_nodes, block_model=geosObj.model_res )
    """
    
    index_array =np.zeros ((len(station_location), ), dtype =np.int32)
    for ii, distance in enumerate(station_location): 
        for jj , nodes in enumerate(x_nodes): 
            if nodes == distance : 
                index_array [ii]= jj
                break 
            elif nodes> distance: 
                min_= np.abs(distance-x_nodes[jj-1])
                max_= np.abs(distance - x_nodes[jj+1])
                if min_<max_: 
                    index_array [ii]= jj-1
                else: index_array [ii]=jj
                break 
    _tema=[]
    for ii in range(len(index_array )):
        a_= block_model[:, int(index_array [ii])]
        _tema.append(a_.reshape((a_.shape[0], 1)))
        
    stationblock = np.concatenate((_tema), axis=1)
    
    return stationblock 
  

def display_infos(infos, **kws):
    """ Display unique element on list of array infos
    
    :param infos: Iterable object to display. 
    :param header: Change the `header` to other names. 
    :Example: 
        
        >>> from pycsamt.geodrill.geocore import display_infos
        >>> ipts= ['river water', 'fracture zone', 'granite', 'gravel',
             'sedimentary rocks', 'massive sulphide', 'igneous rocks', 
             'gravel', 'sedimentary rocks']
        >>> display_infos('infos= ipts,header='TestAutoRocks', 
                          size =77, inline='~')
    """

    inline =kws.pop('inline', '-')
    size =kws.pop('size', 70)
    header =kws.pop('header', 'Automatic rocks')

    if isinstance(infos, str ): 
        infos =[infos]
        
    infos = list(set(infos))
    print(inline * size )
    mes= '{0}({1:02})'.format(header.capitalize(),
                                  len(infos))
    mes = '{0:^70}'.format(mes)
    print(mes)
    print(inline * size )
    am=''
    for ii in range(len(infos)): 
        if (ii+1) %2 ==0: 
            am = am + '{0:>4}.{1:<30}'.format(ii+1, infos[ii].capitalize())
            print(am)
            am=''
        else: 
            am ='{0:>4}.{1:<30}'.format(ii+1, infos[ii].capitalize())
            if ii ==len(infos)-1: 
                print(am)
    print(inline * size )
    
def fit_default_layer_properties(layers, dbproperties_= ['hatch', 'colorMPL']): 
    """ Get the default layers properties  implemented in database. 
     
    For instance get the hatches and colors from given layers implemented in 
    the database by given the database `dbproperties_`.
    
    :param layers: str or list of layers to retrieve its properties
        If specific property is missing , ``'none'`` will be return 
    :param db_properties_: str, list or database properties 
    :return: property items sanitized
    
    :Example: 
        
        >>> import pycsamt.geodrill.geocore as GC
        >>> GC.fit_default_layer_properties(
        ...    ['tuff', 'granite', 'evaporite', 'saprock']))
        ... (['none', 'none', 'none', 'none'],
        ...     [(1.0, 1.0, 0.0), (1.0, 0.0, 1.0), (0.0, 0.0, 1.0),
        ...     (1.0, 0.807843137254902, 1.0)])
    """
    # for consistency check again and keep the DB properties untouchable.
    dbproperties_= ['colorMPL' if g.lower().find('mpl') >=0 else 
                    'FGDC' if g.lower()=='fgdc'else g.lower() 
                    for g in dbproperties_]
    if isinstance(layers , str): layers =[layers]
    assert_gl = ['yes' if isinstance(ll, str) else 'no' for ll in layers]
    if not len(set(assert_gl)) ==1: 
        raise TypeError("Wrong given layers. Names should be a string!")
    if 'name' or'__description' not in  dbproperties_: 
        dbproperties_.insert(0, 'name')
    
    __gammaProps = GeoStratigraphy._getProperties(dbproperties_)
    
    r_props =[['none' for i in layers] for j in range(len(__gammaProps)-1)]
    for k  , l in enumerate(layers): 
        if l  in __gammaProps[0] :
            ix= __gammaProps[0].index(l)
            for kk, gg in enumerate(r_props) : 
                gg[k]= __gammaProps[1:][kk][ix]
                
    r_props = [GU._sanitize_db_items(r_props[k], force=True )
               for k in range(len (r_props))]
    return tuple(r_props)    
 
def __build_ps__token(obj):
    """ Build a special token for each GeoStratigraphic model. Please don't 
    edit anything here. Force editing is your own risk."""
    import random 
    random.seed(42)
    __c =''.join([ i for i in  [''.join([str(c) for c in obj.crmSites.shape]), 
     ''.join([str(n) for n in obj.nmSites.shape]),
    ''.join([l for l in obj.input_layers]) + str(len(obj.input_layers)), 
    str(len(obj.tres))] + [''.join(
        [str(i) for i in [obj._eta, obj.beta, obj.doi,obj.n_epochs,
                          obj.ptol, str(obj.z.max())]])]])
    __c = ''.join(random.sample(__c, len(__c))).replace(' ', '')                                               
    n= ''.join([str(getattr(obj, f'{l}'+'_fn'))
                         for l in ['model', 'iter', 'mesh', 'data']])
    n = ''.join([s.lower() 
                 for s in random.sample(n, len(n))]
                ).replace('/', '').replace('\\', '')
    
    return ''.join([n, __c]).replace('.', '')
        
def fit_tres(lns, tres, autorocks, force=False, **kws): 
    """ Read and get the resistivity values from tres that match the 
     the given layers.
     
    Find the layers and  their corresponding resistivity values from the 
    database especially when values in the TRES and LN are not the same
    length. It's not possible to match each value to its
    correspinding layer name. Therefore the best approach is to read the
    TRES and find the layer name in the database based on the closest value.

    :param lns: list of input layers 
    :param tres: list of input true resistivity values 
    :param autorocks: list of the autorocks found when building the new model.
    :param force: bool, force fitting resistivity value with the rocks in 
            the database whenever the size of rocks match perfectly 
            the number of the rocks. Don't do that if your are sure that the 
            TRES provided fit the  layers in LNS.
    :param kws: is database column property. Default is
        `['electrical_props', '__description']`
        
    :returns: new pseudolist contains the values of rocks retrived from 
        database as well as it closest value in TRES.
    """
    def flip_back_to_tuple(value , substitute_value, index=1): 
        """convert to tuple to list before assign values and 
          reconvert to tuple  after assignment for consistency. 
          `flip_back_to_tuple` in line in this code is the the same like :
                newTRES[ii] = list(newTRES[ii])
                newTRES[ii][1] = val
                newTRES[ii] = tuple (newTRES[ii]) 
          """ 
        value = list(value)
        if index is not None: 
            value[index] = substitute_value
        else : value = substitute_value
        return tuple (value) 

     
    ix = len(autorocks)
    lns0, tres0, rlns, rtres= GU.lns_and_tres_split(ix,  lns, tres)
    if len(lns0) > len(tres0): 
        msg= ''.join(['Number of given layers `{0}` should not be greater ',
                      ' than the number of given resistivity values` {1}`.'])
        msg= msg.format(len(lns0), len(tres0))
        
        n_rock2drop = len(tres0)-len(lns0) 
        msg += f" Layer{'s' if abs(n_rock2drop)>1 else ''} "\
            f"{func.smart_format(lns0[n_rock2drop:])} should be ignored."
    
        lns0 = lns0[: n_rock2drop]
        warnings.warn(msg)
        _logger.debug(msg)
       
    if sorted([n.lower() for n in lns0]
              ) == sorted([n.lower() for n in lns]): 
        if not force: 
            return lns0, tres0, [(a, b) for a , b in zip(lns0, tres0)]
        
    r0 =copy.deepcopy(tres0)
    # for consistency, lowercase the layer name
    # get the properties [name and electrical properties]  
    # from geoDataBase try to build new list with none values 
    # loop for all layer and find their index then 
    # their elctrical values 
    #           if name exist in database then:
    #           loop DB layers names 
    #           if layer is found then get it index 
    lns0 =[ln.lower().replace('_', ' ') for ln in lns0 ]
    _gammaRES, _gammaLN = GeoStratigraphy._getProperties(**kws)

    newTRES =[None for i in tres0]
    temp=list()
    for ii, name in enumerate(lns0) : 
        if name in _gammaLN: 
            ix = _gammaLN.index (name) 
            temp.append((name,_gammaRES[ix])) 
            
    # keep the lns0 rocks that exists in the database 
    # and replace the database value by the one given 
    #in tres0 and remove the tres value with 
    # unknow layer by its corresponding value.
    if len(temp)!=0: 
        for name, value in temp: 
            ix, val = GU.get_closest_gap (value= value, iter_obj=tres0)
            newTRES[ix]= (name, val) 
            tres0.pop(ix) 
    # try to set the values of res of layer found in 
    # the database is not set = 0 by their corresponding
    # auto -layers. if value is in TRES. We consider 
    #that the rocks does not exist and set to None
    for ii, nvalue in enumerate(newTRES):
        try: 
            iter(nvalue[1])
        except:
            if nvalue is not None and nvalue[1]==0. :
                newTRES[ii]= None 
            continue 
        else: 
            # if iterable get the index and value of layers
            # remove this values in the tres 
            ix, val = GU.get_closest_gap (value=nvalue[1], iter_obj=tres0)
            newTRES[ii] = flip_back_to_tuple (newTRES[ii], val, 1) 
            tres0.pop(ix) 
            
    for ii, nvalue in enumerate(tres0):
        ix,_val=  GU.get_closest_gap (value=nvalue,status ='isoff', 
                                   iter_obj=_gammaRES, 
                          condition_status =True, skip_value =0 )
        # get the index of this values in tres
        index = r0.index (_val) 
        newTRES[index] = (_gammaLN[ix], nvalue)
        
    # create for each tres its pseudorock name 
    # and pseudorock value
    pseudo_lns = [a [0] for a in newTRES] + rlns 
    pseudo_tres = [b[1] for b in newTRES] + rtres 
    newTRES += [(a, b) for a , b in zip(rlns, rtres)]
    
    return pseudo_lns , pseudo_tres , newTRES 

def quick_read_geomodel(lns=GU.LNS, tres=GU.TRES):
    """Quick read and build the geostratigraphy model (NM) 
    
    :param lns: list of input layers 
    :param tres: list of input true resistivity values 
    
    :Example: 
        >>> import pycsamt.geodrill.geocore as GC 
        >>> obj= GC.quick_read_geomodel()
        >>> GC.fit_tres(obj.input_layers, obj.tres, obj.auto_layer_names)
    """
    if len(GU.INVERS_KWS) ==0: 
        _logger.error("NoneType can not be read! Need the basics Occam2D"
                         f" inversion {func.smart_format(GU.k_)} files.")

        raise ValueError("NoneType can not be read! Need the basics Occam2D"
                         f" inversion {func.smart_format(GU.k_)} files.")
        
    geosObj = GeoStratigraphy( input_resistivities=tres, 
                      input_layers=lns,**GU.INVERS_KWS)
    geosObj._createNM()
    
    return geosObj 

def assert_len_layers_with_resistivities(
        real_layer_names:str or list, real_layer_resistivities: float or list ): 
    """
    Assert the length of of the real resistivites with their
    corresponding layers. If the length of resistivities is larger than 
    the layer's names list of array, the best the remained resistivities
    should be topped up to match the same length. Otherwise if the length 
    of layers is larger than the resistivities array or list, layer'length
    should be reduced to fit the length of the given resistivities.
    
    Parameters
    ----------
        * real_layer_names: array_like, list 
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
         # get the last value of resistivities  to match the structures
         # names and its resistivities 
        sec_res = real_layer_resistivities[len(real_layer_names):]        
       # fetch the name of structure 
        geos =Geodrill.get_structure(resistivities_range=sec_res) 
        if len(geos)>1 : tm = 's'
        else :tm =''
        print(f"---> Temporar{'ies' if tm=='s' else 'y'} "
              f"{len(geos)} geological struture{tm}."
              f" {'were' if tm =='s' else 'was'} added."
              " Uncertained layers should be ignored.")
        
        real_layer_names.extend(geos)       
        return real_layer_names 
    elif len(real_layer_names) > len(real_layer_resistivities): 
        real_layer_names = real_layer_names[:len(real_layer_resistivities)]        
        return real_layer_names
                    
                    


