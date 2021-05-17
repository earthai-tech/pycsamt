# -*- coding: utf-8 -*-
"""
===========================================================================
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

===========================================================================

.. _module-cs:: `pycsamt.ff.core.cs` 
    
   :synopsis: Super class or CSAMT Far Field implementation 
               Deal with AVG file and EDI files 
    
Created on Wed Dec  2 11:29:32 2020

@author: KouaoLaurent alias @Daniel03
"""

import os , re, warnings, shutil
import time  
 
import numpy as np
import scipy as sp
from pycsamt.etc import infos as inFO

import pycsamt.ff.core.avg as CSAMTavg 
import pycsamt.ff.core.edi as CSAMTedi 
import pycsamt.ff.core.j as CSAMTj 
import pycsamt.ff.core.z as CSAMTz 

from pycsamt.utils import exceptions as CSex
from pycsamt.utils import gis_tools as gis 
from pycsamt.utils import func_utils  as func
from pycsamt.utils._csamtpylog import csamtpylog

from pycsamt.ff.processing import callffunc as cfunc 
import pycsamt.ff.processing.zcalculator as Zcc

_logger = csamtpylog.get_csamtpy_logger(__name__)

try : 
    import scipy as sp 
    import scipy.stats as spSTATS
    scipy_version = [int(vers) for vers in sp.__version__.split('.')] #[1,4,1]
    if scipy_version [0] == 1 : 
        if scipy_version [1] < 4 :
            warnings.warn('Note: need scipy version 1.4.0 or more . It may probably  get a '
                          'trouble when import "stats" attribute'
                          'under such version.It may probably move from ', ImportWarning)
            _logger.warning('Note: need scipy version 0.14.0 or higher or for'
                            ' stats.linearegress. Under such version'
                            'it might not work.')
    # from sp import stats 
    stats_import =True 
            
except :
    warnings.warn('Could not find scipy.stats, cannot use method linearregression'
                  'check installation you can get scipy from scipy.org.')
    _logger.warning('Could not find scipy.stats, cannot use method linearegress'
                    'check installation you can get scipy from scipy.org.')
    
    stats_import =False 


########################## end import #######################################


class CSAMT(object): 
    """
      CSAMT class is super class container of all the information of the other classes,
     J, Avg and Edi. In fact, the CS object collect all the information of the other classes,
     once questioned for specific uses. The purpose of the construction
     of this object to avoid the repetition of scripts throughout the project.
     Objet CSAMT can recognize the typycal input obj of file and set attributes for use . 
     
    Arguments 
    ---------
        **data_fn** : str  
                full path to EDI, J or AVG file or  
                files directory 
        **edipath** : str  
                path to edifiles directory.,
                where edifiles are located 
        **jpath** : str 
                 full path to AG.Jones file ,
                 directory where jfiles are located 
        **profile_fn** : str 
                full path to Zonge station profile file  "*.stn"
        
    .. note:: to call multiple files , better to specify only  the path. 
    
    =================  ===================  ===================================
    Attributes              Type                Description  
    =================  ===================  ===================================
    lat                 float/ndarray,1     sation latitude 
    lon                 float/ndarray,1     station longitude 
    elev                float/ndarray       station elevantion 
                                            in m or ft 
    east                float/ndarray.1     station easting coordinate (m)
    north               float/ndarray,1     station northing coordinate (m)
    azim                float/ndarray,1     station azimuth in meter (m)
    station             ndarray,1           sation id from survey
    utm_zone            str                 UTM location zone 
    resistivity         dict                resistivity value at each
                                            station (ohm.m)
    resistivity_err     dict                res.error at each  station 
    phase               dict                phase value at each station
                                            in degree 
    phase_err           dict                phase error at each station
    zxy                 dict                impedanceTensor at each station 
                                            from xy
    zxy_error           dict                imped. Tensor error  at each  
                                            station from yx
    zyx                 dict                impedanceTensor  at each station 
                                             from xy
    zyx_err             dict                imped. Tensor error  at each
                                            station from yx
    freq                ndarray,1           frequency array from survey   
    =================  ===================  ===================================

    ===========================  ==============================================
    Methods                         Description 
    ===========================  ==============================================
    _read_csamt_objs                read in CSAMT file [EDI|AVG|J]
    _read_edi_obj                   read_edi_obj and set attributes 
    _read_avg_obj                   read Zonge Eng. Avg file and set attributes
    _read_j_obj                     read A.G. Jones file and set attributes.
    ===========================  ==============================================
     
    :Example:
        
        >>> from pycsamt.ff.core.cs import CSAMT
        >>> profile_stn='K1.stn'
        >>> for ii in ['csi000.dat', 'testemap3.edi', 'K1.AVG']:
        >>> path =  os.path.join(os.environ["pyCSAMT"], 
        ...                  'pycsamt','data', ii)
        >>> csamt_obj = CSAMT(fn = path, 
        ...                      profile_fn= os.path.join(
        ...                          os.path.dirname(path), file_stn))
        ... print(csamt_obj.resistivity['S00'])
        ... print(csamt_obj.phase['S00'])
        ... print(csamt_obj.phase_err['S00'])
    """
    
    def __init__(self, data_fn=None , profile_fn = None , **kwargs): 
        
        self._logging =csamtpylog.get_csamtpy_logger(self.__class__.__name__)
        
        
        self.Site =Site()
        self.Location =Location()
        self.Profile =Profile()
   
        self._z_xy = None 
        self._z_yx =None
        self._z_err =None
        
        self._res=None
        self._res_err=None 
        
        self._phs = None 
        self._phs_err =None 
        
        self._freq = None 
        
        self._fn = data_fn 
        
        self._pfn =profile_fn
        self._fpath =None 
        
        self.save_path = kwargs.pop('savepath', None)
        
        for key  in list(kwargs.keys()): 
            setattr(self, key, kwargs[key])
        
        if self._fn is not None or self._fpath  is not None: 
            self._read_csamt_objs()

    @property 
    def lat (self):
        return self.Location.latitude 
        
    @property 
    def lon (self):
        return self.Location.longitude 
    
    @property 
    def elev (self): 
        return self.Location.elevation 

    @property 
    def east(self): 
        return self.Location.easting 
    @property 
    def north (self): 
        return self.Location.northing 
    
    @property 
    def station (self): 
        return self.Site.stn_name 
   
    
    
    #----- SET ATTRIBUTES ---------
    
    @station.setter 
    def station (self, stn):
        self.Site.stn_name = stn 
        
    @elev.setter 
    def elev (self, elev) : 
        self.Location.elevation = elev 
        
    @lat.setter 
    def lat(self, latitude) : 
        self.Location.latitude =latitude 
        
    @lon.setter 
    def lon(self, longitude):
        self.Location.longitude = longitude 
    @north.setter 
    def north (self, northing): 
        self.Location.northing =northing
    @east.setter 
    def east(self, easting): 
        self.Location.easting =easting
        
      
    @property 
    def fpath (self): 
        """
        assert path and redirect to file either single file  ,edifiles or Jfiles. 
        find the specific path [EDI|J] path.
        """
        return CSAMT.find_path(path = self._fn )
    
 

    def _read_csamt_objs (self, fn=None ): 
        """
        read cs object and set attributes , actualy can read J, Avg and Edifile 

        Parameters 
        ----------
            * fn : str  
                   full path to csamtfile either edi , avg or jfile
            * edipath : str  
                    full path to edifiles ,  edi directory
            * jpath : str 
                    full path to AG.Jones files , j directory
        """
        self._logging.info("Reading <%s> file from <%s> " % (self._fn ,
                                                             self._read_csamt_objs.__name__))
        
        if fn is not None : self._fn = fn 
     
        # if path is not None : self._path =path 
        
        if  self.fpath =='isfile' : 
            if inFO._sensitive.which_file(filename= self._fn, deep =False) =='avg':
                    self._read_avg_obj(avgfile = self._fn)
            elif inFO._sensitive.which_file(filename= self._fn, deep =False) =='edi':
                    self._read_edi_obj(edi_fn = self._fn )
            elif inFO._sensitive.which_file(filename= self._fn, deep =False) =='j':  
                    self._read_j_obj(j_fn= self._fn  )
            else : 
                warnings.warn('Currently  pyCSAMT work with "edi", "j" and "avg" file. '\
                              'Please provided one amomg them.')
                self._logging.warn(
                    'It seems file provided <%s> does not match '\
                        '"edi" nor "avg" nor "j" file.' % self._fn)
                raise CSex.pyCSAMTError_file_handling(
                    'Error reading <%s> file. Please check your'\
                        ' file.'%os.path.basename(self._fn) )
                
        elif self.fpath   =='edipath' :
            self._read_edi_obj(edi_fn = self._fn  )
            
        elif self.fpath =='jpath' : 
            self._read_j_obj(j_fn = self._fn  )
            
            
        
    def _read_edi_obj (self, edi_fn =None ):
        """
        read edifiles and populates important attributes 
   
        :param edi_fn: full path to edifile, 
                        can be  edi directory , 
                        where edifiles are located . 
        :type edi_fn: str 
            
        
        1. Read single edifile
        
        :Example:
  
            >>> from pycsamt.ff.core.cs import CSAMT 
            >>> path = os.path.join(os.environ['pyCSAMT'], 
            ...                        'data', 'edi')
            >>> csamt_obj = CSAMT(
                data_fn = os.path.join(path,'csi000.dat')
                )
            ... print(csamt_obj.resistivity['S00'])
            ... print(csamt_obj.lat)
            ... print(csamt_obj.freq)
                
        2. Read read multiple edifiles  
         
        :Example:
      
            >>> from pycsamt.ff.core.cs import CSAMT 
            >>> csamt_obj = CSAMT( data_fn = edipath )
            ... print(csamt_obj.resistivity['S05'])
            ... print(csamt_obj.resistivity_err['S00'])
        """
        
        self._logging.info("Reading <%s> edifile from edi_obj "% edi_fn)
         
        if edi_fn is not None : self._fn =edi_fn
        
        # if self._fn is not None : 
        if self.fpath  =='edipath': #---> set the path then collect all jfiles .
        
            edifiles =[os.path.join(self._fn, eddfile) for
                       eddfile in os.listdir(self._fn) if eddfile.endswith('.edi')]
        
        # if edifile is not None : self._fn =edifile 
        
        elif self.fpath  == 'isfile' : #-- > once only edifile is provided then put on list ...
            edifiles =[self._fn]

        edi_obj = CSAMTedi.Edi_collection(list_of_edifiles =edifiles)

        self._logging.info ('Compute Z impedance Tensor  and '\
                            'set attributes. _ Func <%s>'%self._read_edi_obj.__name__)
        
        self._res = edi_obj.res_xy 
        self._res_err = edi_obj.res_err_xy
        
        self._phs_err =edi_obj.phs_err_xy
        self._phs = edi_obj.phs_xy
        
        self._z_xy =edi_obj.z_xy 
        self._z_yx =edi_obj.z_yx
        self.z_err_xy =edi_obj.z_err_xy
        self.z_err_yx =edi_obj.z_err_yx
        
        self._freq = edi_obj.freq_array
        
    
        self.lat =edi_obj.latitude
        self.lon= edi_obj.longitude
    
        self.Location.convert_location_2_utm(latitude =self.lat, longitude =self.lon)
        self.east = self.Location.easting 
        self.north =self.Location.northing
        self.elev =edi_obj.elevation 
        self.station = edi_obj.id
        
        
        
    def _read_j_obj(self, j_fn =None, jpath =None): 
        """
        Read A. G. Jones files and populate attributes. 
        
        :param j_fn:  full path to AG Jones format 
                    or full path to all jfiles .
        :type j_fn: str
  
        1.  read single jfile 
        
        :Example: 
 
            >>> from pycsamt.ff.core.cs import CSAMT 
            >>> path = os.path.join(os.environ['pyCSAMT'], 
                               'data', 'j')
            >>> csamt_obj = CSAMT(
                fn = os.path.join(data_fn,
                                  'csi000.dat'))
            ... print(csamt_obj.resistivity['S00'])
            ... print(csamt_obj.lat)
            ... print(csamt_obj.freq)
             
        2. Read multiple files 

        :Example:
                     
            >>> from pycsamt.ff.core.cs import CSAMT 
            >>> csamt_obj = CSAMT( data_fn = path )
            ... print(csamt_obj.resistivity['S05'])
            ... print(csamt_obj.lat)
        """
        self._logging.info("Reading <%s> edifile from edi_obj "% j_fn)
        
        if j_fn is not None : self._fn =j_fn 
        
        # if self._fn is not None : 
        if self.fpath =='jpath':        #---> set the path then collect all jfiles .
            jfiles =[os.path.join(self._fn, jjfile) for jjfile in os.listdir(self._fn)]

                # if jfile is not None : self._fn = jfile 
        elif self.fpath =='isfile' :  #-- > once only jfile is provided then put on list ...
                jfiles =[self._fn]
                

        j_obj = CSAMTj.J_collection(list_of_jfiles= jfiles )
        
        self.lon = j_obj.longitude
        self.lat = j_obj.latitude
        self.elev=j_obj.elevation
        
  
        self.station= j_obj.stnames
        
        self.jsites_infos = j_obj.J.jinfo.site_infos
        self.jprogvers= j_obj.J.jinfo.progvers 
  
        self.Location.convert_location_2_utm(latitude =self.lat, longitude =self.lon)
        
        self.station =j_obj.id 
        # self.azim =np.stack ((self.east, self.north), axis= 1)

        
        self._res = {key :rhovalue for key ,
                     rhovalue in zip (self.station , j_obj.app_rho)}

        self._phs = {key : phsvalue  for key, 
                     phsvalue in zip (self.station, j_obj.phase)}
        
        self._res_err ={key: res_err - self._res[key] for key ,
                        res_err in zip (self.station ,  j_obj.rhomax)}
        
        self._phs_err = {key: phs_err -self._phs [key] for key ,
                         phs_err in zip (self.station , j_obj.phamax) }


        #---> get frequency through the first list of period element 
   
        self._freq = 1/j_obj.period[0]
        
        
        # self.z_err_xy , self._z_xy, self.z_err_yx , self._z_yx = [{} for ii in range(4)]
        
        #--- > compute impedance Tensor and set attributes accordingly -----
        
        #--> build Z_obj
        
        self._logging.info ('Compute Z impedance Tensor and set attributes. _ Func <%s>'%self._read_j_obj.__name__)
        
        self._z_err_xy,  self._z_xy,   self._z_yx = [{} for ii in range(3)]
        csamt_z_obj =CSAMTz.Z()
        for stn, values in self._res.items(): 
            # csamt_obj =  CSAMTz.Z()
            phs_array , phs_err_array = np.zeros ((self._freq.size, 2, 2),
                                                  dtype =np.float),  np.zeros (
                                                      (self._freq.size, 2, 2), 
                                                      dtype= np.float)
            res_array , res_err_array = np.zeros ((self._freq.size,2,2),
                                                  dtype =np.float), np.zeros (
                                                      (self._freq.size, 2, 2),
                                                      dtype=np.float)
            phs_array[:, 0, 1]  , phs_err_array [:, 0, 1] = self._phs[stn] , self._phs_err[stn]
            res_array [:, 0, 1] , res_err_array [:, 0, 1]  = values , self._res_err[stn]
  
            
            csamt_z_obj.set_res_phase( res_array =res_array ,phase_array = phs_array,  
                                      freq=self._freq, 
                                      res_err_array= res_err_array, 
                                      phase_err_array= phs_err_array)
                                  
            self._z_err_xy [stn] = csamt_z_obj.z_err [:, 0, 1]
            self._z_xy [stn] = csamt_z_obj.z [:, 0, 1]
            self._z_yx [stn] = csamt_z_obj.z[:, 1, 0]
        
        
      
    def _read_avg_obj (self, avgfile, **kwargs) : 
        """
        Read Zonge Engineering AVg file and and set attributes.
        
        :param avgfile: full path to avg file
        :type avgfile: str 
                
        :Example:  
            
            >>> path =  os.path.join(os.environ["pyCSAMT"], 
                              'pycsamt','data', file_stn)
            >>> from csampty.ff.core.cs import CSAMT
            >>> csamt_obj = CSAMT(fn = os.path.join(
                os.path.dirname(path), file_1), profile_fn= path)
            ... print(csamt_obj.resistivity)
            ... print(csamt_obj.phase)
            ... print(csamt_obj.resistivity_err)
        """
        profile_fn = kwargs.pop('profile_fn', None)
        
        if profile_fn is not None :self._pfn = profile_fn 
            
        
        northing = kwargs.pop('easting', None)
        easting =kwargs.pop('northing', None)
        
        latitude = kwargs.pop('latitude', None)
        longitude=kwargs.pop('longitude', None)
        
        elevation =kwargs.pop('elevation', None)
        utm_zone = kwargs.pop('utm_zone', '49N')
        
        #---- . Build avg_obj and stn profile obj 
        if self._pfn is not None : 
            self.Profile.read_stnprofile(profile_fn = self._pfn, 
                                         easting =easting, 
                                         northing= northing , 
                                         elevation =elevation , 
                                         latitude =latitude , 
                                         longitude =longitude )
            
        # profile_obj =Profile()
        # profile_obj.read_stnprofile(profile_fn =self._pfn , 
        #                             easting =easting , 
        #                             northing =northing ,
        #                             elevation  =elevation, 
        #                             latitude =latitude , 
        #                             longitude =longitude )
        
        avg_obj = CSAMTavg.Avg(data_fn=avgfile)
        
        self._logging.info ('Compute Z impedance Tensor,Rho and phase '\
                            ' and set attributes. _ Func <%s>'%self._read_avg_obj.__name__)
        #--- set attribute --- 
        self._res =avg_obj.Data_section.Resistivity.loc
        self.station = avg_obj.Data_section.Station.names
        self._data_section = avg_obj.Data_section._data_array
        
        emag_obj ,hmag_obj =avg_obj.Data_section.Emag.loc, avg_obj.Data_section.Hmag.loc  
        c_var_emag_obj , c_var_hmag_obj = avg_obj.Data_section.pcEmag.loc , avg_obj.Data_section.pcHmag.loc
        c_var_app_rho_obj , std_phase_obj =avg_obj.Data_section.pcRho.loc, avg_obj.Data_section.sPhz.loc

        
        self._res_err = {key : Zcc.compute_sigmas_e_h_and_sigma_rho(pc_emag= c_var_emag_obj[key],
                                                                      pc_hmag= c_var_hmag_obj[key] ,
                                                                      pc_app_rho=c_var_app_rho_obj[key],
                                                                      app_rho= self._res[key], 
                                                                      emag= emag_obj[key],
                                                                      hmag=hmag_obj[key])[-1] for key , value in self._res.items()}
        # phase value in Zonge AVG is on milirad , may convert into rad before to degree. 
        self._phs = { stn: 180 * phase_value *1e-3 / np.pi for stn, phase_value in avg_obj.Data_section.Phase.loc.items()}
        
        # self._freq = avg_obj.Data_section.Frequency.loc 
        self._freq=avg_obj.Data_section.Frequency.value 
        self._phs_err = {stn : 180 * sphz_value /100 * 1e-3/ np.pi 
                         for stn , sphz_value  in  std_phase_obj.items()}
        
        #--> build Z_obj 
        self._z_err_xy,  self._z_xy,   self._z_yx = [{} for ii in range(3)]
        csamt_z_obj =CSAMTz.Z()
        for stn, values in self._res.items(): 
            # csamt_obj =  CSAMTz.Z()
            phs_array , phs_err_array = np.zeros ((self._freq.size, 2, 2),
                                                  dtype =np.float),  np.zeros (
                                                      (self._freq.size, 2, 2), 
                                                      dtype= np.float)
            res_array , res_err_array = np.zeros ((self._freq.size,2,2),
                                                  dtype =np.float), np.zeros (
                                                      (self._freq.size, 2, 2),
                                                      dtype=np.float)
            phs_array[:, 0, 1]  , phs_err_array [:, 0, 1] = self._phs[stn] , self._phs_err[stn]
            res_array [:, 0, 1] , res_err_array [:, 0, 1]  = values , self._res_err[stn]
  
            
            csamt_z_obj.set_res_phase( res_array =res_array ,phase_array = phs_array,  
                                      freq=self._freq, 
                                      res_err_array= res_err_array, 
                                      phase_err_array= phs_err_array)
                                  
            self._z_err_xy [stn] = csamt_z_obj.z_err [:, 0, 1]
            self._z_xy [stn] = csamt_z_obj.z [:, 0, 1]
            self._z_yx [stn] = csamt_z_obj.z[:, 1, 0]
        
        #---> set coordinates attributes 
        
        # self.north= profile_obj.north
        self.north = self.Profile.north
        # self.east =profile_obj.east 
        self.east = self.Profile.east 
        self.azim =np.stack ((self.east, self.north), axis= 1)

        if self.lat is None: 
            
            self.Location.convert_location_2_latlon(utm_zone=utm_zone )

        self.lat = self.Location.latitude  
        self.lon= self.Location.longitude
        # self.azim =profile_obj.azimuth
        # self.elev =profile_obj.elev
       
        self.elev = self.Profile.elev
 
    #-----------------------------------------------------
    # reseetting specific attributes 
    #----------------------------------------------------
    @property 
    def freq (self): 
        return self._freq 
    @property 
    def resistivity (self): 
        return self._res 
    
    @property 
    def phase (self): 
        return {stn: value %90 for stn , value in self._phs.items()} 
    
    @property 
    def zxy (self): 
        return self._z_xy 
    
    @property 
    def zyx (self): 
        return self._z_yx 
    
    @property 
    def resistivity_err (self): 
        return self._res_err 
        
        
    @property 
    def phase_err (self): 
        return self._phs_err 
    
    @property 
    def z_err (self): 
        return self._z_err_xy 
    
    @property 
    def dipolelength (self): 
        return self.Profile.compute_dipolelength_from_coords(easting =self.east, 
                                                             northing =self.north)[0]
    @property 
    def skindepth (self): 
        return {stn: 503 * np.sqrt( self.resistivity [stn]/ self.freq) 
                for stn in self.resistivity}
    @property 
    def doi (self): 
        return {stn :  2 ** 0.5 * self.skindepth[stn]  for stn in self.skindepth }
    
    @property 
    def station_distance (self): 
        return self.Profile.compute_dipolelength_from_coords(easting =self.east, 
                                                         northing =self.north)[1]
    @property 
    def station_separation(self): 
        return self.Profile.stn_separation(easting=self.north, 
                                           northing =self.north, 
                                           interpolate=True)[0]
 
    
    
    
    @staticmethod    
    def find_path (path =None, ptol =0.7):
        """
        Check path and return filepath , edipath or jpath .
        
        :param path:  full path to file  or directory 
        :type path: str or pathlike 
        
        :param ptol: tolerance given by the program ,
                    less or egal to 1 
                     
        :type ptol: float
        :returns: specific path
        :rtype: str 
        
        .. note :: tolerence param inspects  the number of EDI or J file
                    located on the path and determine the typical
                     path of files either edipath or jpath.
        """
        
        if path is None : return 
        if path is not None : 
            if os.path.isfile (path) is True : return 'isfile'
            elif os.path.isdir(path) is True :
                if os.listdir(path) is not None : 
                    ex = [file for file in os.listdir(path) if 
                          os.path.splitext(file)[1] =='.edi']
                    if len(ex)/len(os.listdir(path))>= ptol :return 'edipath'
                    elif len(ex)/len(os.listdir(path)) < ptol : 
                        m=[]
                        try : m= [file for file in os.listdir(path) 
                                  if inFO._sensitive.which_file(
                                          filename = os.path.join(path, file)) =='j' ]
                        except : pass 
                        if len(m)/len(os.listdir(path)) >= ptol :return 'jpath'
                            
                    return 
    
#=================================================================================
# PROFILE CLASS 
#=================================================================================
    
class Profile (object): 
    """
    Profile class deal with  AVG Zonge station file and statation locations 
    coordinates could be find in *.stn* file or SEG-EDI file. 
     
     
    :param profile_fn:  Path to Zonge *STN file of 
                        SEG-EDI locations or Zonge station file
    :type profile_fn: str 
                
    .. Note :: When  EDI file is called , EDI-collecton
                 auto populated  profile attributes and 
                 coordinates are automatically rescalled. 
       
    ================  ============  ===========================================
    Attributes         Type         Explanation
    ================  ============  ===========================================
    Location           class        Location class for Easting  Northing 
                                    azimuth details. 
    profile_angle      float        If user doesnt Know   the angle profile . 
                                    He can use the method "get_profile_angle"  
                                    to get the value of profile angle. 
    stn_interval       (ndarray,1)  Array of station separation. 
    dipole_length      float        Dipole length value 
                                    computed automatically.          
    lat/lon            (ndarray,1)  latitude/longitude of  stations points . 
    east/north         (ndarray,1)  Easting and northing of stations points.  
    azimuth            (ndarray,1)  Azimuth array stations .    
    ele                (ndarray,1)  Elevation array at each station points.
    stn_position       (ndarray,1)  Station position occupied  
                                    at each stations.
    ================  ============  ===========================================
    
    ==========================  ===============================================
    Methods                     Description
    ==========================  ===============================================
    read_stnprofile             Read the profile file .
    get_profile_angle           Compute the profile  line and strike angle .
    reajust_coordinates_values  Reajustment of coordinates values.  
    straighten_profileline      Redress the profile line . 
    rewrite_station_profile     Rewrite the station profile
                                and generate a new stn file.
    stn_separation              Compute the stations  separations .      
    compute_dipole\
        length_from_coords      compute dipolelength
    ==========================  ===============================================

    More attributes can be added by inputing a key word dictionary
    
    :Example: 
        
        >>> from pycsamt.ff.core.cs import Profile 
        >>> file_stn = 'K1.stn'
        >>> path =  os.path.join(os.environ["pyCSAMT"], 
        ...                      'pycsamt','data', file_stn)
        >>> profile =Profile (profile_fn=path)
        >>> profile.straighten_profileline(
            X=profile.east, Y=profile.north,straight_type='n')
        >>> profile.rewrite_station_profile(
        ...            easting=profile.east, 
        ...            northing=profile.north, 
        ...            elevation =profile.elev, 
        ...            add_azimuth=True)
        >>> separation = profile.stn_separation(
        ...    easting = profile.east, 
        ...    northing =profile.north)
    """
    def __init__(self , profile_fn=None ,  **kwargs):
        
        self._logging = csamtpylog.get_csamtpy_logger(self.__class__.__name__)
        self.profile_fn =profile_fn 
        self.Location =Location ()
        self.Site =Site()
         
        self.profile_angle=None
        self._stn_position =None 

        self.stn_interval=None
        self.dipole_length =None 
        self.azimuth =None 
        
        self.utm_zone =kwargs.pop('utm_zone', '49N')
        
        for keys in list(kwargs.keys()): 
            setattr(self, keys, kwargs[keys])
        
        
        if self.profile_fn is not None : 
            self.read_stnprofile()
            
    #---- set profile properties -----
    @property 
    def lat (self): 
        return self.Location.latitude 
    @property 
    def lon (self): 
        return self.Location.longitude 
    @property 
    def elev (self): 
        return self.Location.elevation  
    @property 
    def north (self): 
        return self.Location.northing
    @property 
    def east (self): 
        return self.Location.easting 
    @property 
    def stn_position (self): 
        return self.Location.stn_pos
  
    # ---- set functions ----
    @lat.setter 
    def lat(self, latitude): 
        self.Location.latitude =latitude 

    @lon.setter 
    def lon (self, longitude): 
        self.Location.longitude = longitude 

    @north.setter 
    def north (self, northing): 
        self.Location.northing =northing
        
    @east.setter 
    def east(self, easting): 
        self.Location.easting = easting
 
    @elev.setter 
    def elev (self, elevation): 
        self.Location.elevation = elevation 
    @stn_position.setter 
    def stn_position (self, position_array): 
        self.Location.stn_pos = position_array

        
    def read_stnprofile (self, profile_fn =None ,
                         easting=None , northing=None ,
                         elevation =None , split_type=None,
                         **kwargs):
        """ 
        Method to read profile station file.
        user can use its special file .user can specify 
        a head of its file. method will read and will parse 
        easting , northing , elevation , or 
        lon, lat, elev or station . User can also provided 
        easting , northing and elevation value . 
        
        Parameters  
        -----------
            * profile_fn :str 
                        path to station profile file 
            * split_type :str 
                        How data is separed . 
                        Default  is "".
            * easting : array_like  
                        easting coordinate (m), 
            * northing : array_like  
                        northing coordinate value (m)
            * lat : array_like 
                    latitude coordinate in degree 
            * lon : array_like  
                    longitude coordinate in degree
            * azim : array_like ,
                    azimuth in degree 
                    If not provided can computed automatically 
            * utm_zone :str      
                    survey utm zone 
                    if not porvided and lat and lon is set ,
                    can compute automatically 
        """
        lat = kwargs.pop('latitude', None)
        lon = kwargs.pop('longitude', None)
        azim =kwargs.pop('azimuth', None)
        utm_zone =kwargs.pop('utm_zone', None)
        
        if utm_zone is not None : self.utm_zone = utm_zone

        _pflag =0
        if profile_fn is not None : 
            self.profile_fn =profile_fn

        if self.profile_fn is not None : _pflag = 1
        elif easting is not None and northing is not None : _pflag = 2
        elif lat is not None and lon is not None : _pflag = 3
        elif self.profile_fn  is None and (easting is not None or 
                                           northing is None) and (lat is None
                                                                  or lon is None) : 
            raise CSex.pyCSAMTError_profile('Provided at least easting or lon or Northing or lat value or profile stn file.')    
        

        if _pflag ==1:
            if inFO._sensitive.which_file(filename =self.profile_fn ) =='stn': 
                ref='none' 
                if os.path.basename(self.profile_fn).find('_reaj') >=0  or\
                    os.path.basename(self.profile_fn).find('_cor') >=0 or \
                        os.path.basename(self.profile_fn).find('_sca') >=0 : 
                    ref ='scalled'
   
                with open(self.profile_fn, 'r', encoding='utf8') as fn :
                    data_lines =fn.readlines ()
                    # if user doesnt provide the split type : program will search automatically 
                    if split_type is None : 
                        split_type = func.stn_check_split_type(data_lines=data_lines)
                    
                    # in the case where user provide the file written by that software ,
                    #ignore the  #healines start by '> or '!
                   
                    temp=[]             # rebuild a new data lines for safety
                    for hh , values  in enumerate(data_lines) : 
                        if re.match(r'^>', values) is not None or re.match('r^!', values) is  not None \
                            or values.find('++++++++') >=0  : pass 
                        else :temp.append(values)
                    if 'len' in temp[0] and 'sta' in temp[0] and 'azim' in temp[0] : 
                        ref ='pyCSAMT'    
                    data_lines=temp
                    
                    decision , stn_headlines_id = _validate_stnprofile(
                        profile_lines= data_lines , spliting=split_type)
                    if decision <2:
                        raise CSex.pyCSAMTError_profile(
                            'Please provide  for profile at least the Easting'
                             ' and northing coordinates values or lat/lon values for parsing. ')
            
                    else :
                        data_list_of_array= [np.array (line.strip().split(split_type)) 
                                             for ii,  line in enumerate(data_lines)]

                        data_array =func.concat_array_from_list(
                            list_of_array=data_list_of_array, concat_axis=0)
                #  catch atributes and set attributes wherever the position it takes on the *.stn files.

                for lab, index in stn_headlines_id: 
                    if ref =='scalled': 
                        warnings.warn("It seems the data <{0}> is already"
                                      "scaled. We don't need to rescale. "
                                      "To force scaling , please use "
                                      "func.{Profile.reajust_coordinates_value}. ")
                        coords_array = data_array[1:, index]
                    else : 
                        self._logging.info('Rescaling stations positions'
                                           ' from file <%s>'% os.path.basename(self.profile_fn))
                        warnings.warn('Scaling stations positions : Usually Zonge'
                            ' Hardware provides stations locations at'
                            '  each electrodes NOT in center of each dipoles. '
                              'Station locations and EM component orientations are taken'
                              ' at each electrodes point.So , we will rescale stations  '
                              ' to dipolecenter position.Distance will stay the same along'
                              '  profile but the number of points will be minus 1.')
                        try :
                            coords_array =cfunc.dipole_center_position(
                                dipole_position = data_array[1:,index])
                        except :pass 
          
                    if lab == inFO.suit.easting[0] or lab.lower().find('east')>=0 : 
                        self.__setattr__('east', coords_array)
                    if lab == inFO.suit.northing[0]  or lab.lower().find('north')>=0:
                        self.__setattr__('north', coords_array)
                    if lab == inFO.suit.elevation[0] or lab.lower().find('elev')>=0 : 
                        self.__setattr__('elev', coords_array)
                    if lab.lower().find('lat')>=0 :
                        self.__setattr__('lat', coords_array)
                    if lab .lower().find('lon')>=0: 
                        self.__setattr__('lon', coords_array)
                    if ref =='pyCSAMT': 
                        if lab=='len' or lab.find('len') >=0 :
                            self.__setattr__('stn_position', coords_array)
                        if lab=='azim' or lab.find('azim')>=0 : 
                            self.__setattr__('azimuth', coords_array)
                    elif ref =='none':
                        if lab in inFO.suit.station :
                            self.__setattr__('stn_position', coords_array)

            #---- compute azimuth , station interval and dipole_length ---------

            if (self.east is not None ) and (self.north is not None ): 
                if self.azimuth is None :
                    self.azimuth = func.compute_azimuth(easting=self.east,
                                                        northing =self.north, 
                                                    extrapolate=True)
  
                self.stn_interval=self.stn_separation(easting= self.east, 
                                                      northing =self.north,
                                                      interpolate=True)[0]
                
                self.dipole_length = round_dipole_length(self.stn_interval.mean())
                
                try : 
                    self.Location.convert_location_2_latlon(utm_zone=self.utm_zone)
                except : 
                    self._logging.debug(
                        'Could not convert easting/northing to lat/lon'
                        ' with utm_zone = {}.'.format(self.utm))
        
            
            elif (self.lon is not None ) and (self.lat is not None): 
                
                #convert lat lon to utm and get the attribute easting and northing 
                self.Location.convert_location_2_utm(longitude=self.lon,
                                                     latitude=self.lat)
                self.azimuth =func.compute_azimuth(easting=self.east,
                                                   northing=self.north, extrapolate=True)
                self.stn_interval=self.stn_separation(easting= self.Location.easting, 
                                                      northing =self.Location.northing, 
                                                      interpolate=True )
                
                self.dipole_length = round_dipole_length(self.stn_interval.mean())
                
    
        
        elif _pflag ==2 or _pflag ==3 : 
            if _pflag == 2 : 
                assert easting.size == northing.size ,\
                    CSex.pyCSAMTError_profile('Easting and Northing must have the same size.'
                                 ' easting|northing size are currentlysize is <{0}|{1}>.'.\
                                     format(easting.size, northing.size))
                self.east =easting 
                self.north =northing
                
                self.Location.convert_location_2_latlon(utm_zone=self.utm_zone)
                
                if self.azimuth is None :
                    self.azimuth = func.compute_azimuth(easting=self.east, 
                                                        northing=self.north)
                else :self.azimuth =azim 
                
            elif _pflag ==3 :
                assert lat.size == lon.size ,\
                    CSex.pyCSAMTError_profile('Easting and Northing must have the same size.'
                          ' lat|lon size are currentlysize is <{0}|{1}>.'.format(lat.size, lon.size))
                self.Location.convert_location_2_utm(latitude=lat ,
                                                     longitude= lon )
                self.easting =self.Location.easting  
                self.northing =self.Location.northing 
                if self.azimuth is None : 
                    self.azimuth = func.compute_azimuth(easting=self.east,
                                                        northing=self.north,
                                                        interpolate=True)
            
            self.dipole_length , self.stn_position =\
                Profile.compute_dipolelength_from_coords(easting =self.east, 
                                                         northing =self.north )
            self.stn_interval=self.stn_separation(easting= self.east, 
                                                      northing =self.north,
                                                      interpolate=True )
            if elevation is None : self.elev= np.full((self.east.size,), .0 ) # set elevation to 0
            else : self.elev = elevation 

                                                        
            
    def get_profile_angle(self, easting =None, northing =None ):
        """
        Method to compute profile angle . 
        
        :param easting:  easting corrdinate of the station point 
        :type easting: array_like 
        
        :param northing: northing coordinate of station point.
        :type northing:array-like 
        
        :returns: profile angle in degrees. 
        :rtype: float
        """
        self._logging.info ('Computing _ profile angle of geoelectric strike.')
        
        if easting is not None : self.east = easting 
        if self.east is None :
            raise CSex.pyCSAMTError_location('Easting coordinates must be provided . ')
        if northing is not None : self.north=northing 
        if self.north is None : 
            raise CSex.pyCSAMTError_location('Northing coordinate must be provided.')
        
        if self.east is not None and self.north is not None :
            try : 
                
              self.east =self.east.astype('float')
              self.north = self.north.astype('float')
            except : 
                raise CSex.pyCSAMTError_float('Could not convert easting/northing to float.')
                   
        # use the one with the lower standard deviation
        if stats_import is True : 
            profile1 =spSTATS.linregress(self.east, self.north)
   
            profile2 =spSTATS.linregress(self.north, self.east)
        else :
            warnings.warn('Could not find scipy.stats, cannot use method linearRegression '
                  'check installation you can get scipy from scipy.org.')
            self._logging.warning('Could not find scipy.stats, cannot use method lineaRegress '
                    'check installation you can get scipy from scipy.org.')
            raise ImportError('Could not find scipy.stats, cannot use method lineaRegress '
                    'check installation you can get scipy from scipy.org.')
            
        
        profile_line = profile1[:2]
        # if the profile is rather E=E(N), the parameters have to converted  into N=N(E) form:
        
        if profile2[4] < profile1[4]:
            profile_line = (1. / profile2[0], -profile2[1] / profile2[0])

        # if self.profile_angle is None:
        self.profile_angle = (90 - (np.arctan(profile_line[0]) * 180 / np.pi)) % 180

        # otherwise: # have 90 degree ambiguity in strike determination# choose strike which offers larger angle with profile
        # if profile azimuth is in [0,90].
    
        self.profile_angle=np.around(self.profile_angle,2)
        
        print('----> profile angle = {} degrees N.E'.format(self.profile_angle))

        # compute pseudo_elect
        # if self.geoelectric_strike is None :
        if  0<= self.profile_angle < 90 :
            self.geoelectric_strike =self.profile_angle + 90 
        elif 90<= self.profile_angle < 180 :
            self.geoelectric_strike =self.profile_angle -90
        elif 180 <= self.profile_angle <270 :
            self.geoelectric_strike = -self.profile_angle +90 
        else :
            self.geoelectric_strike = -self.profile_angle -90 

        self.geoelectric_strike = self.geoelectric_strike % 180 
            

        self.geoelectric_strike =np.floor(self.geoelectric_strike)
        
        print('----> geoelectrike  strike = {} degrees N.E'.format(self.geoelectric_strike))
        
        return self.profile_angle, self.geoelectric_strike
    
    @staticmethod              
    def reajust_coordinates_values ( x=None, y=None,
                                    stn_fn=None,  rewrite=False,
                                    savepath=None,  **kwargs): 
        """
        Simple staticmethod to readjut  coordinates values and
        write new station file. 
        by default , the reajustment substract value. 
        to add value to you old coordiantes , use negative
        X and Y method offer possibility of output new
        file by setting write to True.
        By convention we use X as EASTING correction and Y 
        for NORTHING correction.
        
        Parameters 
        ----------
            * x: float 
                value for ajusting X coordinates _EASTING 
            * y: float 
                value for ajustig Y coordinates ._NORTHING
            * stn_file: str 
                    station profile file . it may be a STN file . 
            * rewrite: bool 
                    rewrite a new station file after reajust coordinates. 
            * savepath : str 
                    outdir pathLike to save your new profile file. 
        
        Returns 
        ---------
            * array_like  
                stations_pk , station profile pka value(m) .
                Electrode fixed point value.
            * array_like 
                easting coordinate value (m)
            * array_like 
                 northing coordinate value (m)
            * elevation : array_like 
                evelation point  at each station (m)
                    
        :Example :  
            
            >>> from pycsamt.ff.core.cs import Profile 
            >>> stn_file =K1.stn
            >>> path =  os.path.join(os.environ["pyCSAMT"],
            ...                         'pycsamt','data',
            ...                         stn_file)
            >>> profile =Profile.reajust_coordinates_values(
            ...           x=-300238.702 ,y=-2369.252  )                                                                        
        """
        
        elev =kwargs.pop('elevation', None)
        dot=kwargs.pop('station_pk', None) 
        east =kwargs.pop('easting', None)
        north=kwargs.pop('northing', None)
        

        if isinstance(x, str) or isinstance(y,str) : 
            try : float(x), float(y)
            except :
                raise CSex.pyCSAMTError_profile('Readjustment not possible with str number.'\
                                                ' Please provide the right values of x and y .')
        if x is None : 
            warnings.warn('In principle, readjustment is not possible with '
                          'NoneType number . However , we gonna set X to 0.')
            x=0. 
        if y is None :
            y=0.
            warnings.warn('In principle, readjustment is not possible with'
                          ' NoneType number . However , we gonna set Y to 0.')
                
        if stn_fn is None :
            if east is None and north is None : 
                warnings.warn ('Not possible to reajust coordinates.'
                               ' Provide at least easting and northing values. ')
                
                raise CSex.pyCSAMTError_station('Could not find station file '
                                                'to read. Can not readjust profile coordinates.')
            
            if elev is None : elev= 0. 
            if dot is None : dot =0. 
        
        if stn_fn is not None :
            if inFO._sensitive().which_file(filename=stn_fn , deep=False )=='stn':
                if inFO._sensitive.which_file(filename=stn_fn ) !='stn': 
                    raise CSex.pyCSAMTError_profile('File provided <{0}>is unacceptable. '
                                                    'Please provide your right stn file.'.\
                                                        format(stn_fn))
            
            with open(stn_fn, 'r', encoding='utf8') as fn :
                fstn= fn.readlines()
            #substract the info line generate by pyCSAMT when create a new STN file 
            for ss , infolines in enumerate(fstn): 
                if re.match(r'^>', infolines) is  None or re.match(r'^!', infolines) is  None : 
                    fstn =fstn[ss:] # that mean no info is on the file
                    break
                
            spT = func.stn_check_split_type(data_lines= fstn)
            if spT is None :spT=' '
            headc = fstn[0].strip().split(spT)
            
            eastindex,dotindex ,northindex,  elevindex = [0 for ii in range (4)]   # initialise index values
            for ii ,item in enumerate(headc): # get index . sometimes stn files are littel messy 
                
                if item.find('"e')>=0 or item.lower().find('eas') >=0: eastindex= ii
                elif item.find('"dot') >=0 or item.lower().find('sta')>=0:dotindex =ii 
                elif item.find('"n') >= 0 or item.lower().find('nor')>=0:northindex = ii 
                elif item.find('"h') >= 0  or item.lower().find('elev')>=0: elevindex = ii
                else :
                    warnings.warn('Prior to provide the right station'
                                  ' profile file. The station file must have as '
                                  'headlines , at least :'
                                  ' ["dot|sta, "e|easting, "n|northing, "h|elev].'
                                  ' Only this head can be parsed.')
                    
                    raise CSex.pyCSAMTError_profile('Your station file profided is wrong. '
                                                    'Only <"dot|sta, "e|easting, '
                                                    '"n|northing, "h|elevation> can be parsed. ')

                
            dot, east, north , elev =[[] for jj in range(4)]
            for ii, item in enumerate(fstn):
                item=item.strip().split(spT)
                for jj, value in enumerate(item): 
                    try : value = float(value)
                    except :pass 
                    else : 
                        if jj == eastindex : east.append(value)
                        if jj == northindex :north.append(value)
                        if jj == elevindex : elev.append(value)
                        if jj == dotindex : dot.append(value)
        
        
        easting, northing ,station_pk , elevation  =np.array(east), \
            np.array(north), np.array(dot), np.array(elev)
        
        if elevation.size  <2 : elevation = np.repeat(0.,east.size)
        if station_pk.size <2 : station_pk = np.repeat(0., east.size )
        
        # easting += x
        # northing +=y
        easting = np.apply_along_axis(lambda xx: xx + x ,0, easting )
        northing =np.apply_along_axis(lambda yy: yy + y ,0, northing) 
        
        if  rewrite is  True : 
            if stn_fn is None : 
                spT = ','
                run_fstn = easting.size
            else :run_fstn = len(fstn)-1
                
            stn_write_lines =[]
            stn_write_lines.append(''.join(['{0:<10}{1}'.format('"""dot"""',spT), 
                                            '{0:>10}{1}'.format('"""e"""', spT), 
                                            '{0:>10}{1}'.format('"""n"""',spT), 
                                            '{0:>10}'.format('"""h"""')]))
            stn_write_lines.append('\n')

            for ii in range (run_fstn ):  # skip the headline.
                stn_write_lines.append(''.join(['{0:<7}{1}'.format(station_pk[ii], spT), 
                                                '{:>12.3f}{}'.format(easting[ii], spT), 
                                                 '{:>12.3f}{}'.format(northing[ii], spT), 
                                                 '{:>8}'.format(elevation[ii])
                                                ]))
                stn_write_lines.append('\n')
            if stn_fn is None : fnew_='STN{0}_reaj'.format(time.timezone)    
            else : 
                fnew_= os.path.basename(stn_fn).split('.')[0] + '_reaj'
                
            with open(''.join([fnew_, '.stn']), 'w', encoding='utf8') as fid : 
                fid.writelines(stn_write_lines)
                
            if savepath is not None :
                shutil.move( os.path.join(os.getcwd(),''.join([fnew_, '.stn'])),
                                          savepath )
            elif savepath is None : savepath =os.getcwd()
    
   
            print('-'*77)
            print('---> New <{0}> station file has been rewritten.'.\
                  format(''.join([fnew_, '.stn'])))
            print('---> savepath : <{0}>'.format(savepath) )
            print('-'*77)
                
            
        return station_pk, easting , northing , elevation      
  
        
    def straighten_profileline (self, X=None, Y=None ,
                                straight_type ='classic', 
                                reajust=(0,0), output =False,**kwargs):
        """
        Method to straighten profile line and/or rescaled 
        coordinates.  User can readjust coordinateq 
        of profile by adding coordinate of readjustation  
        Method provides 3 type of straighten profile.
        Default is "classic", it could be 
        "'natural or distorded', equidistant" . 
        "natural or distorded Type" is not to straight a 
        profile like a straight line but , it keeps the
        equidistant point of the station at normal place
        that the survey must be. sometimes on the field ,
        crew may get around some obstacle and despite the
        line is not straight , the distance  between station 
        is distorded. Using 'distord or natural type ' ,
        it will show the right place station must be.
        
        .. note::  for easier approch we use X
                 as easting and Y as northing. 

        :param X: easting coordinates array. 
        :type X: array_like (ndarry, array,1)
            
        :param Y: northing coordinates array 
        :type Y: array_like (ndarray,1)
        
        :param straight_type: type of straighten ,it could 
                            be "equistant or egal, natural
                            or distord". *default* is "classic"
        :type straight_type: str 

        :param reajust:  coordinates for reajustment (
                            index 0 :x index 1 : y )  
        :type reajust: tuple
        """
        self._logging.info ('Profile coordinate X<{0}>, y<{1}> '\
                            'are readjusting. Method to reajust is <{2}>'.\
                                format(X, Y, straight_type))
        
        savepath =kwargs.pop('savepath', None)
        
        REW=False               # coordinates scaling flag and control new outputs
        
        if X is not None : self.east = X 
        if Y is not None : self.north = Y
        if self.east is None or self.north is None :
            raise CSex.pyCSAMTError_profile('No possible way to straighten up '\
                'the profile line . Please provide the right coordinates. ')
        
        if self.east.size != self.north.size : 
            raise CSex.pyCSAMTError_profile('X and Y must be the same size. '
                       'X has a size <{0}> while Y has the size <{1}>. '
                       'Line cannot be redressed. '
                         'Please provide the same size of both arrays'.\
                             format(self.east.size, self.north.size))
                
        # rescalled the coordinates 
        if reajust is not None :
            if len(reajust) !=2 : 
                raise CSex.pyCSAMTError_profile(
                    "Could not reajust coordinates beyond three values."
                    "Only x =reajust[0] and y=reajust[1] can be accepted. ")
                
            if self.dipole_length is None: 
                print('--> Dipolelength  is computing to straighten out  profile.')
                self._logging.info ('We are computing dipolelength  to straight  profile. ')
                
                dipolLeng , stn_pk = self.compute_dipolelength_from_coords(
                    easting=self.east, northing =self.northing)
                warnings.warn('We are computing dipole length and we assume '
                              'the stations coordinates are in the center of each dipole.'
                              ' Dipole Length is = {0} m and the total '
                              'length is = {1} m.'.format(dipolLeng, (self.east.size -1) *\
                                                          dipolLeng ))
                
            else : 
                stn_pk = np.arange(int(self.dipole_length/2), 
                                    self.east.size * self.dipole_length,
                                    self.dipole_length)
                warnings.warn('Dipole length is = {0} m. Stations are assumed '
                              'to be in center of each dipole.'
                              ' Total length is = {1} m'.\
                                  format(self.dipole_length,
                                      (self.east.size - 1) * self.dipole_length ))
        
            if self.elev is not None : elev =np.around(self.elev,2) 
            else : elev = np.repeat(0., self.east.size)
            # if output is True : REW =False # firtly ,
            #adjust coordinates without output file  set REW to False 
            _, self.east, self.north ,*_= self.reajust_coordinates_values(
                                        easting= self.east, northing =self.north , 
                                        x=reajust[0], y=reajust[1], station_pk=stn_pk, 
                                        rewrite=REW, elevation= elev)
            
            print("---> Locations coordinates are reajusting to straightening out"
                  " profile. Elevation is added.")
            self._logging.info ("Locations coordinates are reajusting to"
                                " straighten out profile and we'll top elevation.")
            REW =True 
        
        #then reascaled 
        #  keep the adjustments values for others purposes.
        self.__setattr__('e_east', self.east)
        self.__setattr__('n_north', self.north)
        
        if re.match(r'^clas+', straight_type) is not None :
            # rj_factor = (self.east [-1] -self.east[0]) /(self.east.size -1) 
            # r_east =np.ones_like(self.east)* np.arange(self.east.size)
            self.east = np.apply_along_axis(lambda rr: rr * (
                (self.east [-1] -self.east[0]) /(self.east.size -1)) + self.east[0], 
                                            0, np.ones_like(self.east)* np.arange(self.east.size))
            
            self.north = np.apply_along_axis(lambda rr: rr * (
                (self.north [-1] -self.north[0]) /(self.north.size -1 )) + self.north[0],
                                             0, np.ones_like(self.north)* np.arange(self.north.size))
            
        if re.match(r'^eq+', straight_type) is not None : 
            self.east =np.linspace(self.east[0], self.east[-1], self.east.size )
            self.north =np.linspace(self.north[0], self.east[-1], self.north.size )
        elif (re.match(r'^nat+', straight_type) is not None) or (
                re.match(r'^dis+', straight_type) is not None) :
            
            #a kind of straighthen with equisdistant value between station but not straight ,
            #king or natural aspect on the site . 
            rf_X = np.array([value - self.east[ii+1] for ii , value in enumerate(self.east) \
                              if ii <= self.east.size - 2 ]).mean()

            rf_Y = np.array([value - self.north[ii+1] for ii , value in enumerate(self.north) \
                              if ii <= self.north.size - 2 ]).mean()
            
            temp_ee , temp_nn = np.zeros_like(self.east), np.zeros_like(self.north) 
            temp_ee [0], temp_nn[0] = self.east[0] , self.north[0]
            for  jj , value in enumerate(temp_ee): 
                if  jj > 0 : temp_ee [jj]=self.east[jj-1]- rf_X
            for  jj , value in enumerate(temp_nn): 
                if  jj > 0 : temp_nn [jj]=self.north[jj-1]- rf_Y
            
            self.east , self.north = temp_ee , temp_nn 
        
        if output is True : # export the file without scaling , set X, Y to 0. 
        # and keep the new easting and Northing coordinates 
            _,self.east, self.north,*_= self.reajust_coordinates_values(
                                        easting= self.east, northing =self.north ,
                                        x=0., y=0., station_pk=stn_pk, 
                                        rewrite=REW, elevation= elev, savepath =savepath)
            
            self._logging.info ("Locations coordinates are reajusting to"
                                " straighten out profile and we'll top elevation.")
            
            print("---> Locations coordinates are reajusting to straightening out"
                  " profile. Elevation is added.")
        
    def rewrite_station_profile (self, easting =None , 
                                 northing=None ,
                                 elevation =None , 
                                 area_name =None,
                                 username =None, 
                                 add_azimuth =False,
                                 **kwargs): 
        """
        Mthod to rewrite station_profile or output new profile 
        by straightening profile throught reajusting location coordinates values.
         User can use this method to create zonge *stn* file if coordinates are known.
         
        :param easting:  easting coordinates (m)
        :type easting: array_like 
        
        :param northing: northing coordinates (m) 
        :type northing: array_like 
        
        :param elevation: elevation values (m)
        :type elevation: array_like 
        
        :param username: name  of user 
        :type username: str 
        
        :param add_azimuth: compute azimuth  
                            positive down(clockwise)
        :type add_azimuth: bool   
        """
        
    
        utm_zone = kwargs.pop('UTM_Zone', None)
        dipole_length =kwargs.pop('dipole_length', None)
        output_name =kwargs.pop('output_name', None)
        savepath =kwargs.pop('savepath', None)
        if output_name is None : output_name='new_profile'
        
        if utm_zone is not None : self.utm_zone = utm_zone

        write_profile_lines =[]
        if area_name is None : area_name ='' 
        if username is None : username =''
        if easting is not None : self.east =np.array(easting )
        if northing is not None : self.north = np.array(northing) 
        if elevation is not None :self.elev = np.array(elevation)
        
        if self.east is None or self.north is None : 
            self._logging.warn("It seems you didnt provide any easting"
                               " values nor northing value.")
            raise CSex.pyCSAMTError_profile("You may provide value before"
                                            " rewriten station profile file.!")
        
        write_profile_lines.append(''.join(['{0:<17}'.format(
            '>LOCATION'),':'," {0:<55}".format(area_name) ])+'\n')
        write_profile_lines.append(''.join(['{0:<17}'.format(
            '>USERNAME'),':'," {0:<55}".format(username) ])+'\n')
        write_profile_lines.append(''.join(['{0:<17}'.format('>DATE'),':'," {0:<70}".\
                                            format(time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime())) ])+'\n')
        write_profile_lines.append(''.join(['{0:<17}'.format(
            '>UTM_ZONE'),':'," {0:<5}-{1:<7}".format('WGS84',self.utm_zone) ])+'\n')    
        write_profile_lines.append(''.join(['{0:<17}'.format(
            '>SOFTWARE'),':'," {0:55}".format('pyCSAMT') ])+'\n')
        
        if self.lat is None  or self.lon is None : 
            self.Location.convert_location_2_latlon(utm_zone=self.utm_zone)
        
        if self.east is None or self.north is None : 
            self.Location.convert_location_2_utm(latitude=self.lat, longitude=self.lon)
        if add_azimuth is True :
            
            if self.azimuth is not None : azim = self.azimuth
            else :
                self._logging.info('Computing azimuth value.')
                azim = func.compute_azimuth(easting=self.east,
                                            northing=self.north)
                if azim.size != self.east.size or azim.size != self.north.size:
                    # recompute azimuth using extrapolation 
                    azim = func.compute_azimuth(easting=self.east,
                                            northing=self.north, extrapolate=True)
                    
        if dipole_length is None :
            self._logging.info('Automatic dipole length calculation and '
                               'stations position  values are set with.')
            
            self.dipole_length, self.stn_position = self.compute_dipolelength_from_coords(
                easting =self.east, northing = self.north)
            warnings.warn("Dipole length is computing automatically and set"
                          " station position with by default. Dipole length is ={0} m"
                          "and length along the line is ={1} m .".\
                              format(self.dipole_length,
                                     (self.east.size -1)* self.dipole_length))
            
        elif dipole_length is not None : 
            self.dipole_length =dipole_length 

            self._logging.infos('We will set compute station position '
                                'considering your dipole_length value provided. ')
            warnings.warn('We will compute station position using your default dipole length .'\
                          'We will the new station step to <{0}>.'.format(self.dipole_length))
            
            stn_position =np.arange(0,self.east.size * self.dipole_length, self.dipole_length )

            self.stn_position= stn_position
            
                
                
        normalize_distance = np.apply_along_axis(
            lambda dd: dd - self.stn_position[0], 0, self.stn_position)
        stnNames =['S{0:02}'.format(ii) for ii in range(self.east.size)]
        
        #set elevation is not provided .set elevation to 0. 
        if self.elev is None : self.elev = np.full((self.east.size,),0.)
        #write_data 

        write_profile_lines.append(''.join(["{0:^5}".format('""sta'), 
                                            "{0:^9}".format('""len(m)')]))
        write_profile_lines.append(''.join([ "{0:^12}".format(jj) for jj 
                                            in ['""east(m)','""north(m)',
                                                '""lat(°)', '""long(°)',
                                                   '""elev(m)']]))
        if add_azimuth  : 
            write_profile_lines.append('{0:^12}'.format('""azim(°)'))
            write_profile_lines.append(''.join(['\n', '+'*(12*6+14) ,'\n']))
        else :write_profile_lines.append(''.join(['\n', '+'*(12*5+14),'\n']))
        
        
        for ii  in range(self.east.size): 
            
            write_profile_lines.append(''.join(['{:<5}'.format(stnNames[ii]), 
                                                   '{:<7}'.format(normalize_distance[ii]),
                                                   '{:<11.2f}'.format(self.east[ii]),
                                                   '{:<13.2f}'.format(self.north[ii]), 
                                                   '{:<15.10f}'.format(self.lat[ii]), 
                                                   '{:<17.10f}'.format(self.lon[ii]), 
                                                   '{:<10.2f}'.format(self.elev[ii])
                                               ]))
            if add_azimuth is True : 
                write_profile_lines.append('{0:<7.3f}'.format(azim[ii]))
            write_profile_lines.append('\n')
        
        with open (''.join([output_name,'.stn']),'w', encoding='utf8') as fid: 
                fid.writelines(write_profile_lines)
        
        if savepath is not None :
            shutil.move (os.path.join(os.getcwd(),
                                     ''.join([output_name,'.stn'])),
                         savepath)
            
    def stn_separation(self, easting  =None , northing =None ,interpolate =False): 
        """
        Compute the station separation 
        Distance between every two stations 
        
        Parameters 
        ----------
            * easting : array_like (ndarray, 1)
                        easting coordinates  
            * northing : array_like (ndarray,1) 
                        northing coordinates 
            * interpolate : bool 
                        if interpolate is True  will extend to N+1 number to much
                        excatly the number  of electrode. If false ,
                        it match the number of dipole N.
                        *Default* is False.
                
        Returns  
        --------
             array_like 
                 separation value array  
             float 
                separation mean  or average separation value  
        """
        if (easting.dtype !='float') or (northing.dtype !='float'): 
            try :
                easting = np.array([float(ss) for ss in easting])
                northing =np.array([float(ss) for ss in northing])
            except : raise CSex.pyCSAMTError_profile("NoneType can not be computed."\
                                                    "Please provide the right coordinates!.")
        
        if easting.size != northing.size :raise CSex.pyCSAMTError_profile(
                "Both coordinates array must have the same size."
                 "The first argument size is <{0}>, while the second is <{1}>.".\
                     format(easting.size, northing.size))
        
        # for kk in range(self.east.size): 
        #     if kk <=self.east.size -2 : 
        #         np.sqrt((self.east[kk+1]-self.east[kk-1])**2 + (self.north[kk+1]-self.north[kk-1])**2)
        stn_sep =np.array([np.sqrt((easting[kk+1]-easting[kk])**2 + (
            northing[kk+1]-northing[kk])**2) for 
                           kk in range (northing.size) if
                           kk <=northing.size -2])

        if interpolate is True : 
            xx , xx_new = np.arange(stn_sep.size),np.arange(stn_sep.size +1) 
            ff = sp.interpolate.interp1d(xx, stn_sep, fill_value='extrapolate')
            stn_sep= ff(xx_new)
  

            
        return stn_sep , stn_sep.mean()
    
    @staticmethod
    def compute_dipolelength_from_coords(easting=None, northing=None,
                                         **kwargs): 
        """
        Fonction to compute dipole length from coordinates easting 
        and northing values.
        
        :param easting: array of easting coordinate in meters
        :type easting: array_like 
                
        :param northing: array of northing coordinate in meters
        :type northing:  array_like 
        
        :param lat:  latitude coordinate  in degree 
        :type lat: array_like (ndarray,1) 
               
        :param lon: longitude coordinate  in degree 
        :type lon: array_like (ndarray,1) 
        
        :param reference_ellipsoid: id, Ellipsoid name,
                                Equatorial Radius, square of 
                                eccentricity ,default is 23 
        :type reference_ellipsoid: int 
        
        :returns: length of dipole during survey approximated . 
        :rtype: float 
         
        :returns:  position of dipole from reference station  
        :rtype: array_like(ndarray,1)
        
        .. note:: the first electrode is located at 0 and second electrode to 
                 dipole length i.e [0, 50 , ..., nn*50] where nn number 
                 of point -1. Data are relocated in center position
                 of dipole.
        """
        latitude =kwargs.pop('latitude', None )
        longitude =kwargs.pop('longitude', None)
        
        reference_ellipsoid =kwargs.pop('reference_ellipsoid', 23)
        
        if easting is not  None  and northing is not None : 
            assert easting.size == northing.size , CSex.pyCSAMTError_profile(
                'Easting and Northing must have the same size.'
                ' Easting|northing size are currentlysize is <{0}|{1}>.'.\
                    format(easting.size, northing.size))

        elif latitude is not None and longitude is not None : 
            assert latitude.size == longitude.size , \
                CSex.pyCSAMTError_profile('Latitude and longitude must have the same size.'
                                           ' Easting|northing size are currentlysize'
                                               ' is <{0}|{1}>.'.format(easting.size, 
                                                                       northing.size))
            
            #-- convert lat lon to easting northing 
            utm_zone , easting, northing = gis.ll_to_utm(reference_ellipsoid,
                                                         latitude, longitude)
        else : 
            raise CSex.pyCSAMTError_profile('May input at least easting and '
                                            'northing coordinates values or'
                                                ' latitude and longitude values.')
            
            
        distance_east = np.array([easting[kk+1]-easting[kk] 
                                  for kk in range (easting.size -1)])
        distance_north = np.array([northing[kk+1]-northing[kk] 
                                   for kk in range (northing.size -1)])
        
        distance = np.array([np.sqrt(distance_east[kk]**2 + distance_north[kk]**2)
                             for kk in range(distance_east.size) ])
        dipole_length =round_dipole_length(distance.mean()) 
        
        #build stn_pk_position : 
            
        stn_pk = np.array([ ii* np.full((easting.size), dipole_length)[ii]
                           for ii in range(np.full((easting.size), dipole_length).size)])

        
        
        return dipole_length , stn_pk    

#=================================================================================
# LOCATION CLASS 
#=================================================================================
                
class Location (object): 
    
    """
    Details of sation location . Classe used to convert 
    cordinnates and check values for lat/lon , east/north 

    ==================  ====================  =================================
    Attributes              Type                Description  
    ==================  ====================  =================================
    latitude            float/ndarray,1         sation latitude 
    longitude           float/ndarray,1         station longitude 
    elevation           float/ndarray           station elevantion in m or ft 
    easting             float/ndarray.1         station easting coordinate (m)
    northing            float/ndarray,1         station northing coordinate (m)
    azimuth             float/ndarray,1         station azimuth  in meter 
    stn_pos             ndarray,1               sation dipoleposition
    utm_zone            str                     UTM location zone 
    ==================  ====================  =================================
    
    ============================  =============================================
    Methods                         Description 
    ============================  =============================================
    convert_location_2_utm          convert position location lon/lat in 
                                    utm easting northing 
    convert_location_2_latlon       convert location postion  from east/north 
                                    to latitude/longitude 
    ============================  =============================================
    """


    def __init__(self, **kwargs) :
        

        self.datum ='WGS84'
        self._latitude =None 
        self._longitude =None 
        self._easting =None 
        self._northing =None
        self._elevation =None 
        self._utm_zone =None
        self._stn_pos = None 
        self._azimuth = None 
        
    
        for key in list(kwargs.keys()): 
            self.__setattr__(key, kwargs[key])
    
    @property 
    def utm_zone (self): 
        return self._utm_zone
        
    @property 
    def latitude(self): 
        return self._latitude
    @property 
    def longitude (self) : 
        return self._longitude
    @property 
    def easting(self ): 
        return self._easting 
    @property 
    def northing (self):
        return self._northing 
    @property 
    def elevation(self): 
        return self._elevation
    @property 
    def stn_pos (self): 
        return self._stn_pos
    @property
    def azimuth(self): 
        if self.east is not None and self.north is not None : 
                return func.compute_azimuth(easting=np.array([float(eas) for eas in self.easting]),
                                            northing=np.array([float(nor) for nor in self.northing]))
        else : return self._azimuth 
    
    #-----------setting ---- 
    @utm_zone.setter 
    def utm_zone (self, utm_zone): 
        if isinstance(utm_zone,(float,int)):
            warnings.warn('Wrong UTM zone input. Must be a str number.')
            raise CSex.pyCSAMTError_location('UTM Zone must be string, '
                                             'not <{0}>type.'.flroat(type(utm_zone)))
        else : 
            try :
                float(utm_zone[:-2])
            except : raise CSex.pyCSAMTError_location(
                    'Error UTM Zone designator. Both first letters'
                    ' provided are not acceotable !')
            else : 
                if utm_zone[-1] not in list(inFO.notion.utm_zone_dict_designator.keys()): 
                    raise CSex.pyCSAMTError_location(
                        'Wrong UTM Zone letter designator. Letter must be among <{0}>'.\
                                format('|'.join(list(
                                    inFO.notion.utm_zone_dict_designator.keys()))))
            self._utm_zone =utm_zone
    
    @latitude.setter 
    def latitude (self, latitude): 
        if isinstance(latitude, (float,int,str)):
            self._latitude= gis.assert_lat_value(latitude)
        else: self._latitude = np.array([gis.assert_lat_value(lat) for lat in latitude ])
        
    @longitude.setter 
    def longitude(self, longitude) : 
        if isinstance(longitude, (float,int,str)):
            self._longitude= gis.assert_lon_value(longitude)
        else :
            self._longitude = np.array([gis.assert_lon_value(lon) for 
                                        lon in longitude ])
        
    @elevation.setter 
    def elevation(self, elevation) : 
        if isinstance(elevation, (float,int,str)): 
            self._elevation= gis.assert_elevation_value(elevation)
        else :
            self._elevation = np.array([gis.assert_elevation_value (elev)
                                        for elev in elevation]) 
        
    @easting.setter 
    def easting (self, easting): 
        if isinstance(easting, (float,int,str)): 
            self._easting= np.array(easting, dtype=float)
        else : 
            try : self._easting = np.array([ float(east) for east in easting])
            except :raise CSex.pyCSAMTError_float(
                    'Easting must be float or an array of float number.')
            
    @northing.setter 
    def northing (self, northing): 
        if isinstance(northing, (float,int,str)): 
            self._northing= np.array(northing, dtype=float)
        else : 
            try : self._northing = np.array([float(north) for north in northing])
        
            except : raise CSex.pyCSAMTError_float(
                    'northing must be float or an array of float number. ')
        
    @stn_pos.setter 
    def stn_pos (self, stn_pk): 
        try : self._stn_pos =np.array([float(stn) 
                                       for stn in stn_pk])
        except:raise CSex.pyCSAMTError_station('Station pka must be float number.! ')


    @azimuth.setter 
    def azimuth (self, easting_northing):
        print(easting_northing)
        if easting_northing.shape[1] !=2 : 
            raise CSex.pyCSAMTError_azimuth('Azimuth expected to get two array_like(ndarray,2)')                                                                                                    
        elif easting_northing.shape[1] ==2 : 
            easting , northing =np.hsplit(easting_northing, 2)
            self._azimuth =func.compute_azimuth(easting=np.array([float(eas)
                                                                  for eas in easting]),
                                           northing=np.array([float(nor) 
                                                              for nor in northing]))

    def convert_location_2_utm (self, latitude =None  , 
                                longitude=None ): 
        """
        Project coordinates to utm if coordinates are in degrees at  
        given reference ellipsoid constrained to WGS84.
         
        :param latitude:  latitude number 
        :type latitude: float 
        
        :param longitude: longitude number 
        :type longitude: float 
        """
        
        if latitude is not None : self.latitude = latitude
        if longitude is not None : self.longitude =longitude
        
        if isinstance(self.latitude, np.ndarray):
            if self.latitude.size >= 1 or self.longitude >=1: # 
                assert  self.latitude.size ==self.longitude.size, CSex.pyCSAMTError_location(
                    'latitude and longitude must be the same size.') 
    
                self.easting = np.array([ gis.ll_to_utm(reference_ellipsoid=23,
                                              lat=self.latitude[ii],
                                              lon=self.longitude[ii]) [1]
                                         for ii in range(self.latitude.size)])
                self.northing =np.array( [ gis.ll_to_utm(reference_ellipsoid=23,
                                              lat=self.latitude[ii],
                                              lon=self.longitude[ii]) [2] 
                                          for ii in range(self.latitude.size)])
            
        else : 
            _data_info_utm =gis.ll_to_utm(reference_ellipsoid=23,
                                          lat=self.latitude,
                                          lon=self.longitude)
            self.utm_zone = _data_info_utm[0]
            self.easting =_data_info_utm[1]
            self.northing=_data_info_utm [2]
        
        
    def convert_location_2_latlon(self, utm_zone =None ): 
        """
        Project coodinate on longitude latitude once  
        data are utm at  given reference ellispoid 
        constrained to WGS-84.
        """
        if utm_zone is not None : self._utm_zone =utm_zone 
        if self._utm_zone is None :
            raise CSex.pyCSAMTError_location(
                'Try to input the utm_zone : e.g.: 49N')
                                                                     
        _data_info_ll = gis.utm_to_ll(reference_ellipsoid=23,
                                      northing= self.northing, 
                                      easting=self.easting, 
                                      zone =self._utm_zone )
        
        self.latitude= _data_info_ll[0]
        self.longitude = _data_info_ll[1]
        
    def get_eastnorth_array_from_latlon(self,arr_lat , arr_lon):
        """
        Method to quicly convert array of latitude and 
        northing into easting northing
        
        :param arr_lat: array of latitude value 
        :type arr_lat: array_like 
        
        :param array_lon: array of longitude value. 
        :type array_lon: array_like 
        
        :returns: easting array
        :rtype : array_like 
        
        :returns: northing array
        :rtype:  array_like 
        """
        array_easting , array_northing =[[] for ii in range(2)]
        for ii in range(arr_lat.size):
            self.convert_location_2_utm(latitude =arr_lat[ii],
                                        longitude=arr_lon[ii])
            array_easting.append(self.easting)
            array_northing.append(self.northing)
        
        return np.array(array_easting), np.array(array_northing)
    
     
#=================================================================================
# SITE CLASS 
#=================================================================================
class Site(object): 
    """
    Specific site object Easy pack data :lat, lon, elev, azim,
    east, north, into dictionnary for easy access .
        
    Arguments 
    ---------
        **data_fn** :str 
                 path to site file , the same file as profile 
                 or X,Y coordinates values 

    :Example: 
        
        >>>  from pycsamt.ff.core.cs import Site 
        >>>  site=Site(data_fn=path)
        >>>  print(site.east['S07'])
        >>>  print(site.north['S09'])
    """
    
    def __init__(self, data_fn =None , **kwargs):
        
        self.sitedata =data_fn
        self.Location=Location()
        self._lat =None
        self._lon =None 
        self._east= None 
        self._north =None 
        self._azim = None 
        self._elev =None 
        self._len_stn=None 
        self._stn_name=None 
        self.utm_zone =kwargs.pop('utm_zone', '49N')
        
        for keys in list(kwargs.keys()):
            setattr(self, keys, kwargs[keys])
            
        if self.sitedata is not None : 
            self.set_site_info()
            
    @property 
    def stn_name(self): 
        return self._stn_name 
    @stn_name.setter 
    def stn_name (self, names_or_numbOfStations ): 
        
        if isinstance(names_or_numbOfStations, int): 
            self._stn_name = ['S{0:02}'.format(ii) 
                              for ii in range(names_or_numbOfStations)]
        else :
            if isinstance(names_or_numbOfStations, np.ndarray):
                self._stn_name = names_or_numbOfStations.tolist()
            elif isinstance(names_or_numbOfStations, list): 
                self._stn_name =names_or_numbOfStations
            
    @property 
    def lat(self): 
        return self._lat 
    @lat.setter 
    def lat (self, latitude):
        
        if isinstance(latitude, list) :latitude =np.array(latitude)
        if self.stn_name is None : 
            self.stn_name = latitude.size
            warnings.warn('You are not provided stations_names : We will defined stations names'
                          ' automatically starting by S-XX [{0},..,{1}]so '
                          'to zip data with latitude. If you dont want this station nomenclature '
                          'please provide station names.'.format(
                              self.stn_name[0], self.stn_name[-1]))
            
        elif self.stn_name is not None : 
            assert len(self.stn_name)== latitude.size,  CSex.pyCSAMTError_site(
                'Stations_names provided '
                 'must have the same size with latitude data. latitude size is <{0}>'.\
                                      format(latitude.size))
        self._lat ={stn:lat for stn, lat in zip (self.stn_name, latitude)}
    
    @property 
    def lon(self): 
        return self._lon
    @lon.setter 
    def lon(self, longitude): 
        if self.stn_name is None : self.stn_name = longitude.size
        if not isinstance(longitude,np.ndarray): longitude =np.array(longitude)
        if len(self.stn_name) != longitude.size :
            raise CSex.pyCSAMTError_site('Station_names|longitude must'
                                         ' have the same size.Longitude size is <{0}>'.\
                                             format(longitude.size))
        self._lon ={stn:lon for stn, lon in zip (self.stn_name, longitude)}
        
    @property 
    def elev(self): 
        return self._elev
    @elev.setter 
    def elev(self, elevation): 
        if self.stn_name is None : self.stn_name = elevation.size
        if not isinstance(elevation,np.ndarray):
            elevation =np.array(elevation)
        if len(self.stn_name) != elevation.size :raise CSex.pyCSAMTError_site(
                'Station_names|Elevation must'
                 ' have the same size.Elevation size is <{0}>'.format(elevation.size))
        self._elev ={stn:elev for stn, elev in zip (self.stn_name, elevation)} 
        
    @property 
    def azim(self): 
        return self._azim
    @azim.setter 
    def azim(self, azimuth): 
        if self.stn_name is None : self.stn_name = azimuth.size
        if not isinstance(azimuth,np.ndarray): azimuth =np.array(azimuth)
        if len(self.stn_name) != azimuth.size :
            raise CSex.pyCSAMTError_site('Station_names|Azimuth must'
                                            ' have the same size.Azimuth size is <{0}>'.\
                                                format(azimuth.size))
        self._azim ={stn:azim for stn, azim in zip (self.stn_name, azimuth)}   
        
    @property 
    def north(self): 
        return self._north
    @north.setter 
    def north(self, northing): 
        if self.stn_name is None :
            self.stn_name = northing.size
            warnings.warn('You are not provided stations_names : We will defenied stations names'
                  ' automatically starting by S-XX [{0},..,{1}]so to zip data with latitude. '
                  'If you dont want this station nomenclature '
                  'please provide station names.'.format(self.stn_name[0], self.stn_name[-1]))
                
        if not isinstance(northing,np.ndarray): northing =np.array(northing)
        if len(self.stn_name) != northing.size :
            raise CSex.pyCSAMTError_site('Station_names|Northing must'
                                            ' have the same size.Northing size is <{0}>'.\
                                                format(northing.size))
        self._north={stn:north for stn, north in zip (self.stn_name, northing)}  
        
    @property 
    def east(self): 
        return self._east
    @east.setter 
    def east(self, easting): 
        if self.stn_name is None :
            self.stn_name = easting.size 
        if not isinstance(easting,np.ndarray): easting =np.array(easting)
        if len(self.stn_name) != easting.size :
            raise CSex.pyCSAMTError_site('Station_names|Easting must'
                                            ' have the same size.Easting size is <{0}>'.\
                                                format(easting.size))
        self._east ={stn:east for stn, east in zip (self.stn_name, easting)} 
        

    def set_site_info(self, data_fn = None, easting =None ,
                      northing =None , utm_zone=None ):
        """
        Set-info from site file, can read zonge *stn* profile fine or set easting 
        and northing coordinates.
 
        :param data_fn:  path to site data file . 
                       may Use Stn  or other file of
                       coordinates infos
        :type data_fn: str 
                
        :param utm_zone: Utm zone WGS84 
        :type utm_zone: str            
        """
        if utm_zone is not None : self.utm_zone =utm_zone 
        if data_fn is not None : self.sitedata = data_fn 
        if self.sitedata is None  and easting is None  and  northing is None :
            
            raise CSex.pyCSAMTError_profile('Could not find the file to read. '
                                            'Please provide the right path.')
        elif self.sitedata is not None : 
            if os.path.isfile (self.sitedata) : 
                profile_obj =Profile(profile_fn=self.sitedata)
            
            else : raise CSex.pyCSAMTError_profile('Unable to read  a file.'
                                                   ' Please provide the right file.')
            
        elif easting is not None and northing is not None : 
            profile_obj =Profile().read_stnprofile(easting =easting, northing =northing)
        
        #call profile object and populate the missing coordinate lat/lon, or east/north 
        
        if profile_obj.lat is None : 
            self.Location.easting=profile_obj.east 
            self.Location.northing =profile_obj.north 
            self.Location.convert_location_2_latlon(utm_zone=self.utm_zone )
            self.lon , self.lat = self.Location.longitude, self.Location.latitude 
        else : 
            self.lat, self.lon=profile_obj.lat, profile_obj.lon 

        if (profile_obj.east is None) or (profile_obj.north is None): 
            self.Location.latitude =profile_obj.lat 
            self.Location.longitude =profile_obj.lon 
            self.Location.convert_location_2_utm(latitude=self.Location.latitude, 
                                                 longitude=self.Location.longitude) 
    
            self.east, self.north =self.Location.easting, self.Location.northing 
            
        else : self.east, self.north =profile_obj.east, profile_obj.north
            
        self.elev, self.azim =profile_obj.elev , profile_obj.azimuth
        

        
#--------------- special functions -----------------

def _validate_stnprofile(profile_lines, spliting=None): 
    """ 
    Validate station profile, cross ccheck the *stn* file. 
    The Number verification value is set to 2. 
    
    Parameters
    ----------
        * profile_lines: list 
            reading lines 
        * spliting: str 
            type t0 split the stn line . 
    Returns
    -------
        int 
            number of info provided 
        list 
            head_stn infos and the index occupied in the file . 
    """
    
    stn_label =['dot', 'e', 'n', 'h', 'easting',
                'northing', 'elevation', 'elev', 'station', 
                'Northing', 'Easting', 'Elevation','sta' , 
                'east_°', 'north_°', 'east_DMS', 'north_DMS',
                'UTM_zone', 'east_utm', 'north_utm',   'azim_m',  
                'elev_m', 'east', 'lon', 'lat', 'long', 'lon','len' ]
    ok=0
    for ss in stn_label : 
        if ss in profile_lines [0]:
            ok +=1 

    head_index =[ (lab, jj) for jj , val in enumerate(profile_lines[0].strip().split(spliting)) \
                  for lab in stn_label if lab in val]
    return  ok, head_index   


def round_dipole_length(value): 
    """ 
    small function to graduate dipole length 5 to 5.
    Goes to be reality and simple computation 
    """ 
    mm = value % 5 
    if mm < 3 :return np.around(value - mm)
    elif mm >= 3 and mm < 7 :return np.around(value -mm +5) 
    else:return np.around(value - mm +10.)
    


if __name__== '__main__': 
    
    file_1='K1.AVG'
    file_2='LCS01.avg'
    file_3='LCS01_2_to_1.avg'
    
    testj= 'csi000.dat'
    testedi='testemap3.edi'
    testavg ='K1.AVG'
 
    file_stn='K1.stn'
    for ii in ['csi000.dat', 'testemap3.edi', 'K1.AVG']:
        path =  os.path.join(os.environ["pyCSAMT"], 
                              'pycsamt','data', ii)
    
        csamt_obj = CSAMT(data_fn=path, profile_fn= os.path.join(os.path.dirname(path), file_stn))
        # print(CSAMT.find_path(path =path))

        print(csamt_obj.resistivity['S00'])
        print(csamt_obj.phase['S00'])
        print(csamt_obj.phase_err['S00'])
        
    
    