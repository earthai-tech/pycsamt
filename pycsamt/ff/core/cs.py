# -*- coding: utf-8 -*-
#       Author: Kouadio K.Laurent<etanoyau@gmail.con>
#       Licence: LGPL
#       Created on Wed Dec  2 11:29:32 2020
"""
.. _module-cs:: `pycsamt.ff.core.cs` 
    
   :synopsis: Super class or CSAMT Far Field implementation 
               Deal with AVG file and EDI files 
"""

import os 
import warnings
import datetime
import numpy as np

import pycsamt.ff.core.avg as CSAMTavg 
import pycsamt.ff.core.edi as CSAMTedi 
import pycsamt.ff.core.z as CSAMTz 
import pycsamt.utils.zcalculator as Zcc
from pycsamt.utils import _p as inFO
from  pycsamt.ff.core import j as CSAMTj 
from pycsamt.ff.site import (Site, Location, Profile) 
from pycsamt.utils import exceptions as CSex
from pycsamt.utils import func_utils as func
from pycsamt.utils._csamtpylog import csamtpylog
_logger = csamtpylog.get_csamtpy_logger(__name__)


class CSAMT(object): 
    """
      CSAMT class is super class container of all the information of the 
      other classes,J, Avg and Edi. In fact, the CS object collect all 
     the information of the other classes, once questioned for specific uses.
     The purpose of the construction of this object to avoid the repetition
     of scripts throughout the project.Objet CSAMT can recognize the typycal
      input obj of file and set attributes for use . 
     
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
        assert path and redirect to file either single file,edifiles or Jfiles. 
        find the specific path [EDI|J] path.
        """
        self._fpath = CSAMT.find_path(path = self._fn )
        return self._fpath
    
    # @fpath.setter 
    # def fpath(self, fp): 
    #     self._fpath = fp
 

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
        
        if fn is not None :
            self._fn = fn 
     
        # if path is not None : self._path =path 

        # if  self.fpath =='isfile' : 
        #     if inFO._sensitive.which_file(filename= self._fn,
        #                                   deep =False) =='avg':
                
        #             self._read_avg_obj(avgfile = self._fn)
        #     elif inFO._sensitive.which_file(filename= self._fn, 
        #                                     deep =False) =='edi':
        #             self._read_edi_obj(edi_fn = self._fn )
        #     elif inFO._sensitive.which_file(filename= self._fn,
        #                                     deep =False) =='j':  
        #             self._read_j_obj(j_fn= self._fn  )
        #     else : 
        #         warnings.warn(
        #             'Currently  pyCSAMT work with "edi", "j" and "avg" file. '
        #                       'Please provided one amomg them.')
        #         self._logging.warn(
        #             'It seems file provided <%s> does not match '\
        #                 '"edi" nor "avg" nor "j" file.' % self._fn)
        #         raise CSex.pyCSAMTError_file_handling(
        #             'Error reading <%s> file. Please check your'\
        #                 ' file.'%os.path.basename(self._fn) )
                
        # elif self.fpath   =='edipath' :
        #     self._read_edi_obj(edi_fn = self._fn  )
            
        # elif self.fpath =='jpath' : 
        #     self._read_j_obj(j_fn = self._fn  )
        fp =None 
        if  self.fpath in ('isfile', 'jpath', 'edipath') : 
            fp = inFO._sensitive.which_file(filename= self._fn, deep =False)
            
        else :
            warnings.warn(
                'Currently  pyCSAMT work with "edi", "j" and "avg" file. '
                          'Please provided one amomg them.')
            self._logging.warn('It seems file provided <%s> does not match '
                               '"edi" nor "avg" nor "j" file.' % self._fn)
            raise CSex.pyCSAMTError_file_handling(
                'Error reading <%s> file. Please check your'
                    ' file.'%os.path.basename(self._fn) )
            
        if fp =='avg': 
            self._read_avg_obj(avgfile = self._fn)
                    
        elif self.fpath =='edipath'  or fp =='edi' :
            self._read_edi_obj(edi_fn = self._fn  )
            
        elif self.fpath =='jpath' or fp =='j': 
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
                       eddfile in os.listdir(self._fn) 
                       if eddfile.endswith('.edi')]
        
        # if edifile is not None : self._fn =edifile 
        #-- > once only edifile is provided then put on list ...
        elif self.fpath  == 'isfile' : 
            edifiles =[self._fn]

        edi_obj = CSAMTedi.Edi_collection(list_of_edifiles =edifiles)

        self._logging.info ('Compute Z impedance Tensor  and '
                            'set attributes. _ Func '
                            '<%s>'%self._read_edi_obj.__name__)
        
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
    
        self.Location.convert_location_2_utm(latitude =self.lat, 
                                             longitude =self.lon)
        self.east = self.Location.easting 
        self.north =self.Location.northing
        self.elev =edi_obj.elevation 
        self.station = edi_obj.id
        
        # inherit all edit object attributes 
        func.make_introspection(self , edi_obj)
            
    
        
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
        if self.fpath =='jpath':   #---> set the path then collect all jfiles 
            jfiles =[os.path.join(self._fn, jjfile) 
                     for jjfile in os.listdir(self._fn)]

                # if jfile is not None : self._fn = jfile 
                #-- > once only jfile is provided then put on list ...
        elif self.fpath =='isfile' :  
                jfiles =[self._fn]
                

        j_obj = CSAMTj.J_collection(list_of_jfiles= jfiles )
        self.jfiles_list= j_obj.jfiles_list
        self.lon = j_obj.longitude
        self.lat = j_obj.latitude
        self.elev=j_obj.elevation
        
  
        self.station= j_obj.stnames
        
        self.jsites_infos = j_obj.J.jinfo.site_infos
        self.jprogvers= j_obj.J.jinfo.progvers 
  
        self.Location.convert_location_2_utm(latitude =self.lat,
                                             longitude =self.lon)
        
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
        
        self._logging.info (
            'Compute Z impedance Tensor ''and set attributes.'
            ' _ Func <%s>'%self._read_j_obj.__name__)
        
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
            phs_array[:, 0, 1]  , phs_err_array [:, 0, 1] =\
                self._phs[stn] , self._phs_err[stn]
            res_array [:, 0, 1] , res_err_array [:, 0, 1]  =\
                values , self._res_err[stn]
  
            
            csamt_z_obj.set_res_phase( res_array =res_array ,
                                      phase_array = phs_array,  
                                      freq=self._freq, 
                                      res_err_array= res_err_array, 
                                      phase_err_array= phs_err_array)
                                  
            self._z_err_xy [stn] = csamt_z_obj.z_err [:, 0, 1]
            self._z_xy [stn] = csamt_z_obj.z [:, 0, 1]
            self._z_yx [stn] = csamt_z_obj.z[:, 1, 0]
        
        
        func.make_introspection(self , j_obj)
      
    def _read_avg_obj (self, avgfile=None, **kwargs) : 
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
        
        if profile_fn is not None :
            self._pfn = profile_fn 
        if avgfile is not None : 
            self._fn =avgfile   
        
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
        
        avg_obj = CSAMTavg.Avg(data_fn=self._fn)
        
        self._logging.info ('Compute Z impedance Tensor,Rho and phase '\
                            ' and set attributes. _ '
                            'Func <%s>'%self._read_avg_obj.__name__)
        #--- set attribute --- 
        self._res =avg_obj.Data_section.Resistivity.loc
        self.station = avg_obj.Data_section.Station.names
        self._data_section = avg_obj.Data_section._data_array
        
        emag_obj ,hmag_obj =avg_obj.Data_section.Emag.loc,\
            avg_obj.Data_section.Hmag.loc  
        c_var_emag_obj , c_var_hmag_obj = avg_obj.Data_section.pcEmag.loc ,\
            avg_obj.Data_section.pcHmag.loc
        c_var_app_rho_obj , std_phase_obj =avg_obj.Data_section.pcRho.loc,\
            avg_obj.Data_section.sPhz.loc

        self._res_err = {key : Zcc.compute_sigmas_e_h_and_sigma_rho(
            pc_emag= c_var_emag_obj[key],
            pc_hmag= c_var_hmag_obj[key] ,
            pc_app_rho=c_var_app_rho_obj[key],
            app_rho= self._res[key], 
            emag= emag_obj[key],
            hmag=hmag_obj[key])[-1] for key , value in self._res.items()}
        
        # phase value in Zonge AVG is on milirad , 
        #may convert into rad before to degree. 
        self._phs = { stn: 180 * phase_value *1e-3 / np.pi 
                     for stn, phase_value
                     in avg_obj.Data_section.Phase.loc.items()}
        
        # self._freq = avg_obj.Data_section.Frequency.loc 
        self._freq=avg_obj.Data_section.Frequency.value 
        self._phs_err = {stn : (180 * sphz_value) /100 * (1e-3/ np.pi) 
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
            phs_array[:, 0, 1]  , phs_err_array [:, 0, 1] = \
                self._phs[stn] , self._phs_err[stn]
            res_array [:, 0, 1] , res_err_array [:, 0, 1]  =\
                values , self._res_err[stn]
  
            
            csamt_z_obj.set_res_phase( res_array =res_array ,
                                      phase_array = phs_array,  
                                      freq=self._freq, 
                                      res_err_array= res_err_array, 
                                      phase_err_array= phs_err_array)
                                  
            self._z_err_xy [stn] = csamt_z_obj.z_err [:, 0, 1]
            self._z_xy [stn] = csamt_z_obj.z [:, 0, 1]
            self._z_yx [stn] = csamt_z_obj.z[:, 1, 0]
        
        #---> set coordinates attributes 
        if self._pfn  is not None: 
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
 
        func.make_introspection(self , avg_obj)
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
        return self.Profile.compute_dipolelength_from_coords(
            easting =self.east, 
            northing =self.north)[0]
    @property 
    def skindepth (self): 
        return {stn: 503 * np.sqrt( self.resistivity [stn]/ self.freq) 
                for stn in self.resistivity}
    @property 
    def doi (self): 
        return {stn :  2 ** 0.5 * self.skindepth[stn] 
                for stn in self.skindepth }
    
    @property 
    def station_distance (self): 
        return self.Profile.compute_dipolelength_from_coords(
            easting =self.east, 
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
                                          filename = os.path.join(path,
                                                        file)) =='j' ]
                        except : pass 
                        if len(m)/len(os.listdir(path)) >= ptol :return 'jpath'
                            
                    return 
                
                
    def j2edi(self, jfn=None, savepath =None, **kwargs): 
        """
        Method to convert j-files to edi files. Method calls CSAMT class object
        and get from this class edi infos 
        
        :param jfn: collection of jfiles or path-like str  
        :type list_of_files: str  

        :rtype: str
        :returns: edifiles from Jobjects
        
        :Example:
            
            >>> from pycsamt.ff.core.j import J_collection as JObjs
            >>> path2j = 'data/j' 
            >>> jObjs= JObjs().j2edi(path2j)
            
    
        """
        prospect =kwargs.pop('contractor_name', None)
        #hardwareInfos = kwargs.pop('hardware_name', None)
        fileby =kwargs.pop('fileby', 'jediSoftware')
        acqby =kwargs.pop('acqby', 'jedi')
        county =kwargs.pop('county', None)
        project_name =kwargs.pop('project_name', None)
        dipole_length =kwargs.pop('dipole_length', 100.)
        rotation_angle =kwargs.pop ('rotation_angle', 0.)
        
        if jfn is not None : 
            self._fn =jfn 
        
        if self._fn is not None:
            self._read_j_obj(j_fn =self._fn, jpath = self.fpath)
            # self.jfiles_list = jfn

        if self.jfiles_list is None : 
            raise CSex.pyCSAMTError_J(
                'No files found !  Please provide A.G. J-files ')
        # export to savepath 
        if savepath is not None: 
            self.save_path =savepath 
        if self.save_path is None : # create a folder in your current work directory
            try :
                self.save_path = os.path.join(os.getcwd(), '_outputJ2EDI_')
                if not os.path.isdir(self.save_path):
                    os.mkdir(self.save_path)#  mode =0o666)
            except : 
                warnings.warn("It seems the path already exists !")
        
        # call CSAMT obj 
        
        # self= cs_obj.CSAMT(data_fn=jfn)
        
        # create edi-obj and set attributes 
        for ii, stn in enumerate (self.station):
            # create an edi_obj and fill attributes 
            edi_obj=CSAMTedi.Edi()
            
            # fill Head info 
            edi_obj.Head.dataid= stn # ii mean for all the list 
            edi_obj.Head.acqby = acqby 
 
            if self.jsites_infos is None  or\
                self.jsites_infos==[None]: 
                self.jsites_infos = datetime.datetime.fromtimestamp(
                    os.stat(self._fn).st_ctime) # return the creation date of file 
            
            try: 
                edi_obj.Head.acqdate = self.jsites_infos[ii]
            except: #TypeError: 'datetime.datetime' object is not subscriptable
                edi_obj.Head.acqdate= self.jsites_infos
                
            edi_obj.Head.fileby = fileby 
            edi_obj.Head.filedate = datetime.datetime.now(
                ).strftime('%m-%d-%Y %H:%M:%S')
            
            # get the name of location
            edi_obj.Head.loc = os.path.basename(self._fn)[:-4]

            if prospect is None : prospect ='MTnet'
            setattr(edi_obj, 'Head.prospect', prospect)
    
            if county is not None : county ='MT'
            setattr(edi_obj,'Head.county', county)
            
    
            edi_obj.Head.lat = self.lat[ii]
            edi_obj.Head.long = self.lon[ii] 
            edi_obj.Head.elev = round(self.elev[ii], 2)
        
            edi_obj.Head.maxsect =1000
            
            if project_name is None :
                project_name= os.path.basename(self._fn) 

            #=====>  set EDI OBJ INFOS
            # edi_obj.Info.maxinfo = 999
  
            edi_obj.Info.Source.__setattr__('project', project_name)
            edi_obj.Info.Source.__setattr__('survey',  edi_obj.Head.dataid)
            edi_obj.Info.Source.__setattr__('sitename', edi_obj.Head.dataid)
            edi_obj.Info.Processing.__setattr__('processedby', 'pyCSAMT' )
            edi_obj.Info.Processing.ProcessingSoftware.__setattr__(
                'name', edi_obj.Head.fileby )
    
           
            #====> definemeas 
            edi_obj.DefineMeasurement.maxchan =4
            edi_obj.DefineMeasurement.maxrun = len(self.station)
            edi_obj.DefineMeasurement.__setattr__('reftype' ,'CARTesien')
            edi_obj.DefineMeasurement.__setattr__('reflat',edi_obj.Head.lat  ) 
            edi_obj.DefineMeasurement.__setattr__('reflong', edi_obj.Head.long) 
            edi_obj.DefineMeasurement.__setattr__('refelev',round(edi_obj.Head.elev,2))
            
            #creating xxmeas object 
            codeID_dec = '{0}'.format((ii+1)/edi_obj.Head.maxsect)
            # codeID=  '{0:04}{1}'.format(ii * 10 + 1 , codeID_dec[1:] )
            edi_obj.DefineMeasurement.__setattr__('meas_ex', 
                                                  CSAMTedi.Emeasurement(**{
                                                      'id':'{0:04}{1}'.format(ii * 10 + 1 , 
                                                       codeID_dec[1:] ), 
                                                        'chtype':'EX', 
                                                        'x': -(dipole_length/2), 
                                                        'y': 0.,
                                                        'x2':dipole_length/2 , 
                                                        'y2':0, 
                                                        }))

            edi_obj.DefineMeasurement.__setattr__('meas_ey',
                                                  CSAMTedi.Emeasurement(**{
                                                      'id':'{0:04}{1}'.format(ii * 10 + 2 , 
                                                         codeID_dec[1:]), 
                                                        'chtype':'EY', 
                                                        'x':0., 
                                                        'y': -(dipole_length/2),
                                                        'x2':0., 
                                                        'y2':dipole_length/2 , 
                                                            }))
                                                                          

            edi_obj.DefineMeasurement.__setattr__('meas_hx', 
                                                  CSAMTedi.Hmeasurement(**{
                                                      'id':'{0:04}{1}'.format(ii * 10 + 3 , 
                                                       codeID_dec[1:] ), 
                                                        'chtype':'HX', 
                                                        'x':0., 
                                                        'y': 0.,
                                                        'x2':0., 
                                                        'y2':0. , 
                                                        }))

            edi_obj.DefineMeasurement.__setattr__('meas_hy',
                                                  CSAMTedi.Hmeasurement(**{
                                                      'id':'{0:04}{1}'.format(ii * 10 + 4 ,
                                                       codeID_dec[1:]), 
                                                    'chtype':'HY', 
                                                    'x':0., 
                                                    'y': 0.,
                                                    'x2':0., 
                                                    'y2':0. , 
                                                          }))
     
                     
            #====> EMAPSECT
            edi_obj.MTEMAP.sectid = stn 
            edi_obj.MTEMAP.__setattr__('nfreq', len(self.freq))
            edi_obj.MTEMAP.__setattr__(
                'ex', '{0:04}{1}'.format(ii * 10 + 1 , codeID_dec[1:] ))
            edi_obj.MTEMAP.__setattr__(
                'ey', '{0:04}{1}'.format(ii * 10 + 2 , codeID_dec[1:] ))
            edi_obj.MTEMAP.__setattr__(
                'hx', '{0:04}{1}'.format(ii * 10 + 3 , codeID_dec[1:] ))
            edi_obj.MTEMAP.__setattr__(
                'hy', '{0:04}{1}'.format(ii * 10 + 4 , codeID_dec[1:] ))
            
            #Frequency blocks , impendance and resistivity blocs 
            edi_obj.Z.freq = self.freq 
            #add rotation angle 
            edi_obj.Z.rotation_angle = rotation_angle
    
            
            # set phase and resistitivity including error propagation 
            # compute error propagation  
            #-->  initialize ndarray(nfreq, 2, 2) 
            res_array = np.zeros((edi_obj.Z.freq.size, 2,2 ), dtype = np.float)
            res_array_err = np.zeros((edi_obj.Z.freq.size, 2,2 ), dtype = np.float)
            phs_array = np.zeros((edi_obj.Z.freq.size, 2,2 ), dtype = np.float)
            phs_array_err = np.zeros((edi_obj.Z.freq.size, 2,2 ), dtype = np.float)
            
            #dictionnary of components . we set only component into XY . 
            res_array [:, 0 , 1 ] = self.resistivity[stn]  
            res_array_err [:, 0 , 1] = self.resistivity[stn] 
            phs_array[: , 0, 1] = self.phase[stn] 
            phs_array_err  [:, 0, 1]  = self.phase_err[stn] 
        
            #---> Recomputing z with resistivities- phase by using propagrations errors 
            edi_obj.Z.set_res_phase(res_array = res_array, phase_array=phs_array, 
                                    freq=  edi_obj.Z.freq, 
                                    res_err_array=res_array_err,
                                    phase_err_array=phs_array_err)
            
            edi_obj.write_edifile(savepath = self.save_path)
               
        if len(self.station) > 1: 
            print('-'*77)    
            print('---> {0} wrote sucessfully from j-files.\n---> see path:<{1}> '.\
                  format(len(self.station), self.save_path))
            print('-'*77)
                    
    def avg2edi (self,
                 data_fn =None ,
                 profile_fn =None , 
                savepath =None ,
                apply_filter =None, 
                reference_frequency =None, 
                number_of_points=7.,
                dipole_length=50.,
                number_of_skin_depth=3., 
                **kwargs):
        """
        Method to write avg file to SEG-EDIfile.Convert both files.Astatic 
        or plainty avg file .if ASTATIC file is provided , will add the filter 
        and filter values .if avg file is not astatic file , user an apply
        filter by setting filter to "tma, ama, or flma".Once apply ,
        edifiles will exported by computing resistivities filtered
        
        Parameters 
        ----------
            * data_fn : str 
                        full path to avgfile 
            * savepath : str 
                        outdir to store edifiles if None ,
                        is your current work directory                
            * profile_fn: str 
                        full path  to station _profile file
            * apply_filter: str 
                        add the name of filter to process the
                        avg file exported in edifiles. 
                        can be  [TMA | AMA| FLMA]
                        TMA - Trimming Moving Average
                        AMA - Adaptative Moving avarage , 
                        FLMA - Fixed dipoleLength moving 
                        average (number of point=7)
                        Dipolelength willbe computed automatically 
        
        Returns
        --------
            obj, edi_obj  
                edi outputfiles or edi filtered file.
            
        Holdings additionals parameters 
            
        ====================  ==========  =====================================
        Optional                Type        Description            
        ====================  ==========  =====================================
        reference_frequency    float        frequency at clean data.Default is 
                                            computed automatically
        number_of_points=7.    int          number of point to for weighted  
                                            window,. *Default* is 7.
        dipole_length          float        length of dipole in meter. For  
                                            CSAMT, *Default* is 50.
        number_of_skin_depth   float        number of skin_depth can be 1 to 10
                                            when using filter AMA. 
                                            *Default* is 3.
        ====================  ==========  =====================================
                                     
        :Example:   
            
            >>> from pycsamt.ff.core import avg 
            >>> avg_obj= avg.Avg()
            >>> avg_obj.avg_to_edifile(data_fn= os.path.join(
            ...    path_to_avgfile, avgfile) , 
            ...           profile_fn = os.path.join(
            ...    path_to_avgfile, station_profile_file), 
            ...           savepath =save_edipath, 
            ...           apply_filter=None ) 
        """
        # import module Hmeasurement and Emeasurement(
        # from pycsamt.ff.core.edi import (Hmeasurement, Emeasurement)
        from pycsamt.ff.processing import corr 
        # import pycsamt.ff.core.avg as CSAMTavg 
        # import pycsamt.ff.core.edi as CSAMTedi 
        
        
        utm_zone = kwargs.pop('utm_zone', None)
        
        
        if data_fn is not None : 
            self._fn = data_fn 
        if self._fn is None : 
            raise CSex.pyCSAMTError_avg_file(
                "Could not find any path to read ."
                "Please provide your right AVG file.")
        
        
        # export to savepath 
        if savepath is None : # create a folder in your current work directory
            try :
                savepath = os.path.join(os.getcwd(), '_outputAVG2EDI_')
                if not os.path.isdir(savepath):
                    os.mkdir(savepath)#  mode =0o666)
            except : 
                warnings.warn("It seems the path already exists !")
        
        if profile_fn is not None :
            self._pfn = profile_fn 
            
        if self._pfn is None : 
            warnings.warn ('No station profile file will detected.'
                           ' Be aware sure , we will set location longitude ,'
                               ' latitude and elevation to <0.>')
        
        # read avg file 
        if self._pfn is not None and self._fn is not None:
            self._read_avg_obj()
        
        if self._pfn is not None : 
            site_obj = Site(data_fn=self._pfn,
                            utm_zone=utm_zone)
                                                                   
        #---> get the list of stations 
        head_dataid_list = sorted(self.Data_section.Station.names)
        

        #-- compute dipole lenght  between station close 
        dipole_length= self.Data_section.Station.loc[head_dataid_list [1]][0]-\
            self.Data_section.Station.loc[head_dataid_list[0]][0]
        
        #--> nfrequency = 
        nfreq = self.Data_section.Frequency.value.size
        
        
        #-------------------------------------------------------------------
        #---> compute error propagation , phase and resistivity 
        app_rho_obj = self.Data_section.Resistivity.loc 

        phase_obj = self.Data_section.Phase.loc
        emag_obj =self.Data_section.Emag.loc
        hmag_obj =self.Data_section.Hmag.loc 
        c_var_emag_obj = self.Data_section.pcEmag.loc 
        c_var_hmag_obj = self.Data_section.pcHmag.loc
        c_var_app_rho_obj =self.Data_section.pcRho.loc
        std_phase_obj =self.Data_section.sPhz.loc
        
        #self.Header.HardwareInfos.astatic_version
         # maker to check whether is plainty avg file or Astatic file 
        AS_flag =0         
        if self.Header.HardwareInfos.numfilterfreq is not None :
            # self.Header.HardwareInfos.numfilterfreq is not None :
            print('---> Reading Zonge `ASTATIC` file !')
            AS_rho =self.Data_section.Resistivity.loc_Sres
            AS_flag =1
            
        # create 2empty  dicts
        error_propag_rho, error_propag_phs =[{} for ii in range(2)] 
        for key, value in self.Data_section.Resistivity.loc.items(): 
            error_propag_rho[key]= app_rho_obj[key] + 1* \
                Zcc.compute_sigmas_e_h_and_sigma_rho(
                            pc_emag= c_var_emag_obj[key],
                              pc_hmag= c_var_hmag_obj[key] ,
                              pc_app_rho=c_var_app_rho_obj[key],
                              app_rho= app_rho_obj[key], 
                              emag= emag_obj[key],
                              hmag=hmag_obj[key])[-1]
                
        for key, value in self.Data_section.Phase.loc.items(): 
            #  std = sPhz / 100 #converted in radians
            error_propag_phs [key]= 180 *( 
                (std_phase_obj[key]/100)/1000)/np.pi 
  
        
        #convert phase millirad value in degree : 
        phase_obj ={key: 180 * phase * 1e-3/np.pi 
                    for key, phase in phase_obj.items()}
    
        #---------------------------------------------------------------
        if apply_filter is not None :
            try :
                apply_filter= apply_filter.lower()
            except:
                print("---> ErrorType: Filters' names must a str of"
                              " `tma` `flma` or `ama`. not {0}".format(
                                  type(apply_filter)))
                
                warnings.warn("TypeError ,Filters' names must be a str of"
                              " `tma` ,`flma` or  `ama`. not {0}".format(
                                  type(apply_filter)))
                apply_filter=None
            else :
                print('{0:-^77}'.format('Filter Infos'))
                print('** {0:<27} {1} {2}'.format(
                    'Filter applied', '=', apply_filter.upper()))
                
        if apply_filter in ['ama', 'tma', 'flma']:
            # create processing_obj , unique obj for all filters  
            corr_obj= corr.shifting(data_fn = self._fn,
                                    profile_fn=self._pfn)
                
            if apply_filter =='tma':
                self._logging.info (
                    'Computing  rho with Trimming moving average (TMA) filter.')
                        
                # corr_obj= corr.shifting()
                res_TMA= corr_obj.TMA(number_of_TMA_points =number_of_points , 
                                      reference_freq=reference_frequency, 
                                      )
                
                reference_frequency =corr_obj.referencefreq
                
                print('** {0:<27} {1} {2}'.format('TMA filter points', 
                                                  '=', number_of_points))
                print('** {0:<27} {1} {2} Hz.'.format('Reference frequency',
                                                      '=', reference_frequency))
                
            elif apply_filter =='flma':
                self._logging.info (
                    'Computing rho with fixed-length moving average (FLMA) filter.')
     
                res_FLMA = corr_obj.FLMA(number_of_dipole=number_of_points, 
                                        reference_freq=reference_frequency)
                
                dipolelength= corr_obj.dipolelength
                reference_frequency =corr_obj.referencefreq

                print('** {0:<27} {1} {2}'.format('FLMA filter points',
                                                  '=', number_of_points))
                print('** {0:<27} {1} {2} Hz.'.format('Reference frequency',
                                                  '=', reference_frequency))
                print('** {0:<27} {1} {2} m.'.format('Dipole length ',
                                                     '=', dipolelength))
            elif apply_filter == 'ama': 
                self._logging.info (
                    'Computing rho with adaptative moving average (AMA) filter.')
    
                res_AMA = corr_obj.AMA(number_of_skin_depth=number_of_skin_depth, 
                                        reference_freq=reference_frequency)
  
                dipolelength= corr_obj.dipolelength
                reference_frequency =corr_obj.referencefreq
                
                print('** {0:<27} {1} {2}'.format('Number of skin depths',
                                                  '=', int(number_of_skin_depth)))
                print('** {0:<27} {1} {2} Hz.'.format('Reference frequency',
                                                  '=', reference_frequency))
                print('** {0:<27} {1} {2} m.'.format('Dipole length ',
                                                     '=', dipolelength))
     
            elif apply_filter is not None and apply_filter  not in ['tma, ama, flma']: 
                warnings.warn('filter %s provided is UNrecognizing.'
                              ' Pogram worsk currently with : TMA or FLMA ,'
                              ' Please provided the right filters'
                              ' for computing.'% apply_filter)
                                  
                raise CSex.pyCSAMTError_processing(
                    'Filters provided is not acceptable.'
                    ' Recognized filters are "TMA","AMA" AND "FLMA"')
        
        if apply_filter is None : 
            if self.Header.HardwareInfos.numfilterfreq is None\
                and self.Header.HardwareInfos.freq_filter is None :
                 avgfilter, descritipfilter = '', ""   
    
            elif self.Header.HardwareInfos.astatic_version is not None : 
                avgfilter ="{0}.{1}.Filter =TMA for 5pts.".format(
                    'Astatic',self.Header.HardwareInfos.astatic_version )
                descritipfilter = self.Header.HardwareInfos.sconfig_dict['']
            
            else : avgfilter, descritipfilter = '', ""
        
        if apply_filter  is not None :
            msf=''
            if apply_filter =='flma':
                descritipfilter = ' {0} dipoles for FLMA filter'\
                    ' at {1} hertz with {2} m dipole length'.format(
                int(number_of_points), reference_frequency, dipolelength)
                msf='dip'
            elif apply_filter =='tma':
                descritipfilter = ' {0} points for TMA filter at {1} hertz'.format(
                    int(number_of_points), reference_frequency)
                msf='pts' 
                
            elif apply_filter =='ama':
                descritipfilter = ' {0} skin depths for AMA filter at {1} hertz'.format(
                    int(number_of_skin_depth), reference_frequency)
                number_of_points=number_of_skin_depth
                msf='sdepths'

            avgfilter ="{0}.{1}.Filter ={2} for {3} {4}.".format(
                   'pycsamt','v1.0.01', apply_filter.upper(),  
                   int(number_of_points) , msf)
            

        if self._pfn  is None : 
            #---> set to 0. lon , lat and elev if profile _fn is not provided
            #--> create dictionnary of zero value of lat and lon 
            import copy 
            #head_edi_lon = np.full((len(head_dataid_list),), 0., dtype=np.float)
            head_edi_lon = {stn: value for stn , value in zip (
                head_dataid_list,
                np.zeros_like(len(head_dataid_list), dtype=np.float))}
            head_edi_lat = copy.deepcopy(head_edi_lon)
            head_edi_elev = copy.deepcopy(head_edi_lon)
            
        elif self._pfn is not None : 
             head_edi_lon = site_obj.lon
             head_edi_lat = site_obj.lat 
             head_edi_elev = site_obj.elev 
        
        #------------------------START SETTING EDI ATTRIBUTE ----------------
        
        # from pycsamt.utils import gis_tools as gis 
        #   for stn in head_dataid_list : #loop for all station or dataid 
        for ii, stn in enumerate(head_dataid_list): 
            #create Ediobj for eac datalist 
            edi_obj =CSAMTedi.Edi()
            
            #====> set edi_obj Header attributes
            edi_obj.Head.dataid= head_dataid_list[ii] # ii mean for all the list 
            edi_obj.Head.acqby = 'Zonge Engineering'
            edi_obj.Head.acqdate = self.Header.HardwareInfos.dated 
            edi_obj.Head.fileby = "AMTAVG_v.{0}".format(
                self.Header.HardwareInfos.version) 
            if self.Header.SurveyAnnotation.project_name is None :
                 self.Header.SurveyAnnotation.__setattr__(
                     'project_name ',os.path.basename(self._fn)) 
    
            edi_obj.Head.loc = os.path.basename(self._fn)[:-4]
            edi_obj.Head.filedate = self.Header.HardwareInfos.processed
            edi_obj.Head.prospect = self.Header.SurveyAnnotation.contractor_name
            edi_obj.Head.county = self.Header.SurveyAnnotation.project_area 
            
            edi_obj.Head.lat = head_edi_lat[stn]
            edi_obj.Head.long = head_edi_lon[stn]
            edi_obj.Head.elev = round(head_edi_elev[stn],2)
            
            edi_obj.Head.maxsect =1000
            
            
            #=====>  set edi_obj_info . 
            # edi_obj.Info.maxinfo = 999
            edi_obj.Info.filter =self.Header.HardwareInfos.freq_filter
            
               
            edi_obj.Info.Source.__setattr__(
                'project', self.Header.SurveyAnnotation.project_name)
            edi_obj.Info.Source.__setattr__('survey',  head_dataid_list[ii])
            edi_obj.Info.Source.__setattr__('sitename', head_dataid_list[ii])
            edi_obj.Info.Processing.__setattr__('processedby', 'pyCSAMT' )
            edi_obj.Info.Processing.ProcessingSoftware.__setattr__(
                'name', edi_obj.Head.fileby )
            
    
            #====> definemeas 
            edi_obj.DefineMeasurement.maxchan =4
            edi_obj.DefineMeasurement.maxrun = len(head_dataid_list)
            edi_obj.DefineMeasurement.__setattr__('reftype' ,'CARTesien')
            edi_obj.DefineMeasurement.__setattr__('reflat',edi_obj.Head.lat  ) 
            edi_obj.DefineMeasurement.__setattr__('reflong', edi_obj.Head.long) 
            edi_obj.DefineMeasurement.__setattr__('refelev',round(edi_obj.Head.elev,2))
            
            # creating xxmeas object 
            codeID_dec = '{0}'.format((ii+1)/edi_obj.Head.maxsect)
            # codeID=  '{0:04}{1}'.format(ii * 10 + 1 , codeID_dec[1:] )
            edi_obj.DefineMeasurement.__setattr__(
                'meas_ex', 
                CSAMTedi.Emeasurement(**{'id':'{0:04}{1}'.format(ii * 10 + 1 , 
                                                        codeID_dec[1:] ), 
                                        'chtype':'EX', 
                                        'x': -(dipole_length/2), 
                                        'y': 0.,
                                        'x2':dipole_length/2 , 
                                        'y2':0, 
                                        'acqchan': '0.00', 
                                        'filter':avgfilter, 
                                        }))

            edi_obj.DefineMeasurement.__setattr__(
                'meas_ey',
                CSAMTedi.Emeasurement(**{'id':'{0:04}{1}'.format(ii * 10 + 2 , 
                                                        codeID_dec[1:]), 
                                        'chtype':'EY', 
                                        'x':0., 
                                        'y': -(dipole_length/2),
                                        'x2':0., 
                                        'y2':dipole_length/2 , 'acqchan': '0.00', 
                                        'filter':avgfilter}))

            edi_obj.DefineMeasurement.__setattr__(
                'meas_hx', 
                CSAMTedi.Hmeasurement(**{'id':'{0:04}{1}'.format(ii * 10 + 3 , 
                                                        codeID_dec[1:] ), 
                                            'chtype':'HX', 
                                            'x':0., 
                                            'y': 0.,
                                            'x2':0., 
                                            'y2':0. , 'acqchan': '0.00', 
                                        'filter':avgfilter}))

            edi_obj.DefineMeasurement.__setattr__(
                'meas_hy',
                CSAMTedi.Hmeasurement(**{'id':'{0:04}{1}'.format(ii * 10 + 4 ,
                                                        codeID_dec[1:]), 
                                    'chtype':'HY', 
                                    'x':0., 
                                    'y': 0.,
                                    'x2':0., 
                                    'y2':0. , 'acqchan': '0.00', 
                                'filter':avgfilter}))
                
                        
            #====> EMAPSECT
            edi_obj.MTEMAP.sectid = stn 
 
            edi_obj.MTEMAP.__setattr__('nfreq', nfreq)
            edi_obj.MTEMAP.__setattr__('maxblks', 64)
            edi_obj.MTEMAP.ndipole=  len(head_dataid_list) - 1 
            edi_obj.MTEMAP.type = descritipfilter
            edi_obj.MTEMAP.__setattr__('hx', '{0:04}{1}'.format(
                ii * 10 + 3 , codeID_dec[1:] ))
            edi_obj.MTEMAP.__setattr__('hy', '{0:04}{1}'.format(
                ii * 10 + 4 , codeID_dec[1:] ))
            edi_obj.MTEMAP.__setattr__('chksum', 64 * nfreq)
            
            # Frequency blocks , impendance and resistivity blocs 
            edi_obj.Z.freq = self.Data_section.Frequency.value 
            
            # set phase and resistitivity including error propagation 
            # compute error propagation  
            #-->  initialize ndarray(nfreq, 2, 2) 
            res_array = np.zeros((nfreq, 2,2 ), dtype = np.float)
            res_array_err = np.zeros((nfreq, 2,2 ), dtype = np.float)
            phs_array = np.zeros((nfreq, 2,2 ), dtype = np.float)
            phs_array_err = np.zeros((nfreq, 2,2 ), dtype = np.float)
   
            #dictionnary of components . we set only component into XY . 
            res_array [:, 0 , 1 ] = app_rho_obj[stn]
            res_array_err [:, 0 , 1] =  error_propag_rho [stn]
            phs_array[: , 0, 1] = phase_obj[stn] 
            phs_array_err  [:, 0, 1]  = error_propag_phs[stn]
            
            if AS_flag ==1 : 
                res_as_array = np.zeros((nfreq, 2,2 ), dtype = np.float)
                res_as_array [:, 0, 1]  = AS_rho [stn]
            
            if apply_filter is not None : 
                fres_array = np.zeros((nfreq, 2,2 ), dtype = np.float)
                if apply_filter.lower() =='tma' :
                    fres_array  [: , 0 , 1]  = res_TMA[stn]
                elif apply_filter.lower() =='flma' :
                    fres_array  [: , 0 , 1]  = res_FLMA[stn]
                elif apply_filter.lower() =='ama' :
                    fres_array  [: , 0 , 1]  = res_AMA[stn]
                    
            #---> computing z with resistivities , phase by using propagrations errors 
            edi_obj.Z.set_res_phase(res_array = res_array, phase_array=phs_array, 
                                    freq=  edi_obj.Z.freq, 
                                    res_err_array=res_array_err,
                                    phase_err_array=phs_array_err)
            # write edifiles ...
            if AS_flag ==1 :
                edi_obj.write_edifile(savepath = savepath, 
                                      add_filter_array =res_as_array )
                
            elif apply_filter is not None:
               edi_obj.write_edifile(savepath = savepath, 
                                     add_filter_array = fres_array )
               
            else : edi_obj.write_edifile(savepath = savepath)
                
        
        print('-'*77) 
        print('---> {0} Edi-files have been rewritten.\n---> see path:<{1}> '.\
              format(len(head_dataid_list), savepath))
        print('-'*77)    

        
# if __name__== '__main__': 
    
#     file_1='K1.AVG'
#     file_2='LCS01.avg'
#     file_3='LCS01_2_to_1.avg'
    
#     testj= 'csi000.dat'
#     testedi='testemap3.edi'
#     testavg ='K1.AVG'
    
#     data = '/data/avg/K1.AVG'
#     data2= r'C:\Users\Administrator\OneDrive\Python\pyCSAMT\data\avg\K1.AVG'
#     csamt_obj = CSAMT(data_fn=data2)
#     print(csamt_obj.resistivity['S00'])

        
    
    