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
import numpy as np


import pycsamt.ff.core.avg as CSAMTavg 
import pycsamt.ff.core.edi as CSAMTedi 
import pycsamt.ff.core.z as CSAMTz 
import pycsamt.utils.zcalculator as Zcc
from pycsamt.utils import _p as inFO
from  pycsamt.ff.core import j as CSAMTj 
from pycsamt.ff.site import (Site, Location, Profile) 
from pycsamt.utils import exceptions as CSex
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
            if inFO._sensitive.which_file(filename= self._fn,
                                          deep =False) =='avg':
                    self._read_avg_obj(avgfile = self._fn)
            elif inFO._sensitive.which_file(filename= self._fn, 
                                            deep =False) =='edi':
                    self._read_edi_obj(edi_fn = self._fn )
            elif inFO._sensitive.which_file(filename= self._fn,
                                            deep =False) =='j':  
                    self._read_j_obj(j_fn= self._fn  )
            else : 
                warnings.warn(
                    'Currently  pyCSAMT work with "edi", "j" and "avg" file. '
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
        if avgfile is not None : self._fn =avgfile   
        
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

        
    
    