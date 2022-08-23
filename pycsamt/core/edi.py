# -*- coding: utf-8 -*-
#       Created on Mon Jan 11 11:37:51 2021
#       Author: Kouadio K.Laurent<etanoyau@gmail.com>
#       Licence: LGPL
"""
.. _module-edi:: `pycsamt.core.edi`
   :synopsis: EDI module can read and write an .edi file as the 'standard '
             formatof magnetotellurics. Each sectionof the .edi file is given 
             its own class, so the elements of each section are attributes for 
             easy access.Edi file will write following the SEG document
             instructions of  EMAP (Electromagnetic  Array Profiling)
             ...
"""
import os
import re
import warnings
import datetime
import shutil
import time
import numpy as np 

import pycsamt
import pycsamt.utils.func_utils as func
import pycsamt.core.z as MTz
from pycsamt.utils. _p import suit 
from pycsamt.site import  Location
from pycsamt.utils. _p import _sensitive as SB 
from pycsamt.utils import gis_tools as gis 
from pycsamt.utils import exceptions as CSex
from pycsamt.utils._csamtpylog import csamtpylog

_logger = csamtpylog.get_csamtpy_logger(__name__)

    
class Edi_collection : 
    """
    Super class to deal with Edifiles .Collect edifiles and set important 
    properties form Controled Source audiofrequency magnetotelluRic  ,
   two(2) components XY and YX will be set and calculated .
    
    Arguments 
    ---------
        **list_of_edifiles** : list 
                    list of edifiles
        **edipath** : str 
                    path to edifiles. If no listis provided ,
                    provided 
        **survey_name** : location name where the date
                    where collected . if surveyname is None  
                    can chech on edifiles. 
                
    =====================  ===============  ===================================
    Attributes              Type             Description 
    =====================  ===============  ===================================
    edifiles                list            List of edifiles ,default is None 
    survey_name             str             location where edidata where 
                                            collected
    Edi                     class           Ediclass for ediproperties 
    Location                class           Location class for location
                                            properties 
    longitude               ndarray, 1      array of longitude collected
                                            from edifiles
    latitude                ndarray,1       array of latitude collected
                                            from edifiles 
    freq_array              nadarry,1       array of frequencies collected
                                            from edifiles 
    elevation               ndarray ,1      array of elevation collected 
                                            from edifiles. 
    station_names           list             list of EDI dataId , Use station
                                            id to work.
    res_xy|res_yx           dict            dict :{stn: res_xy|res_yx} 
                                            ndarray value of resivities
                                            from 2 comps xy|yx 
    phs_xy|phs_yx           dict            dict :{stn: res_phs|phs_yx} 
                                            ndarray value of phase from 
                                            2  comps xy|yx 
    z_xy|res_yx             dict            dict :{stn: z_xy|z_yx} 
                                            (in degree)
                                            ndarray value of impedance
                                            from 2 comps xy|yx 
    XX_err_xy|XX_err_yx     dict            dict  of error values 
                                            {stn: XX_err_xy|XX_err_yx} 
                                            ndarray value of impedance 
                                            from 2 comps xy|yx 
                                            XX : res|phs|z  stn : name of 
                                            site eg stn :S00, S01 , .., Snn
    =====================  ===============  ===================================
    
    :Example: 
        
       >>> from pycsamt.core.edi import Edi_collection 
       >>> edilist = [os.path.join(path, edi)for edi in os.listdir 
                      (path) if edi.endswith('.edi')]
       ... edi_objs = Edi_collection(list_of_edifiles = edilist)
       ... print(edi_objs.res_xy['S01'])
       ... print(edi_objs.phs_err_xy['S00'])
       ... print(edi_objs.z_err_xy['S00'])
       
       
    """
    
    def __init__(self, list_of_edifiles=None, ediObjs =None , 
                 survey_name =None ): 
        self._logging= csamtpylog.get_csamtpy_logger(self.__class__.__name__)
        self.edifiles =list_of_edifiles
        self.ediObjs = ediObjs 
        self.Location =Location ()
        self.survey_name =survey_name

        #---> construction of georeference attributes :
        for _key in ['_latitude', '_longitude', '_azimuth', 
                     '_elevation', 'station_names']: 
            self.__setattr__(_key, None )

        # construction of impedances , resistivities and
        #phasetensor attributes and set all to None value.
        for tensor  in ['_z', '_z_err', '_z_det','_res', '_res_err',
                        '_res_det','_phs', '_phs_det', '_phs_det_err'
                         ] :  self.__setattr__( tensor, None) 
                                 
        if (self.edifiles or self.ediObjs) is not None : 
            self._collect_edifiles()

        
    @property 
    def latitude(self): 
        return self.Location.latitude
    @latitude.setter 
    def latitude(self, latitude): 
        self.Location.latitude = latitude
        
    @property 
    def longitude(self): 
        return self.Location.longitude
    @longitude.setter 
    def longitude(self, longitude): 
        self.Location.longitude= longitude
        
    @property 
    def elevation(self): 
        return self.Location.elevation 
    @elevation.setter 
    def elevation(self, elevation): 
        self.Location.elevation =elevation
        
    @property 
    def stnames(self):
        return self._station_names
    @stnames.setter 
    def stnames (self, edi_stations):
        if isinstance(edi_stations, np.ndarray ) :
            edi_stations = edi_stations.tolist()
        if isinstance(edi_stations, list) :
            if len(self.edifiles) != len(edi_stations): 
                self._station_names =self.id  # use station id to work 
            else :self._station_names =edi_stations
            
      
    def _collect_edifiles(self, list_of_edifiles =None, *,  ediObjs =None):
        """
        collect edifiles and set appropriates attributes for each stations. 

        :param list_of_edifiles: list of edifiles
        :type list_of_edifiles: list 
        
        :param ediobjs:  can provided from  class built.
        :type ediobjs: pycsamt.core.edi.Edi
        
        :Example: 
            
            >>> from pycsamt.core.edi import Edi_collection 
            >>> path =r'F:\__main__csamt__\paper2_data_old\
                data_edifiles - numStations\K1_edi\new_EDI'
            >>> edilist = [os.path.join(path, edi)for
            ...               edi in os.listdir (path) 
            ...               if edi.endswith('.edi')]
            ... edi_objs = Edi_collection._collect_edifiles(
                list_of_edifiles = edilist)
            ... print(edi_objs.phs_xy['S01'])
            ... print(edi_objs.freq_array)
            
        """
        rfiles =[]
        if list_of_edifiles is not None: 
            self.edifiles = list_of_edifiles 
        if ediObjs is not None: 
            self.ediObjs = ediObjs 
            
        #################
        # data can be a single file or a full path  to edifiles
        if self.ediObjs is not None: 
            if not isinstance(self.ediObjs, (list, tuple)): 
                self.ediObjs = [self.ediObjs] 
            self.ediObjs =list(self.ediObjs) 
            if not isinstance(self.ediObjs[0], pycsamt.core.edi.Edi):
                raise CSex.pyCSAMTError_EDI(
                    "Given object does not match the `pycsamt.core.edi.Edi`"
                    " edi objects."
                    )
            rfiles = self.ediObjs
            
        elif isinstance(self.edifiles, str): 
            # single edi and get the path 
            if os.path.isfile (self.edifiles): 
                edipath = os.path.dirname (self.edifiles) 
            elif os.path.dirname (self.edifiles): 
                # path is given and read  
                edipath = self.edifiles 
            else : 
                raise CSex.pyCSAMTError_EDI(f"Wrong path: {self.edifiles}."
                                            " Please provide the rigth path."
                    )
            rfiles = os.listdir (edipath)
            self.edifiles= sorted ([ os.path.join(edipath, edi ) for edi in 
                               rfiles if edi.endswith ('.edi')])  
  
        elif isinstance(self.edifiles, (tuple, list)): 
            try : 
                self.edifiles= sorted(self.edifiles)
            except TypeError: 
                raise CSex.pyCSAMTError_EDI(
                    "It seems argument for `ediObjs`" 
                    " is passed to the `list_of_edifiles`" 
                    " param and vice versa!")
            rfiles = self.edifiles.copy() 
            
        if self.edifiles is not None: 
            self.ediObjs = [Edi(edi_filename= ediobj) 
                            for  ediobj in self.edifiles]
            
        if self.ediObjs is None: 
            raise CSex.pyCSAMTError_EDI("None EDI object detected!")
            
        # OPTIMIZE 
        self.ediObjs , self.edinames = func.fit_by_ll(self.ediObjs) 
        
        # self.ediObjs =list(self.ediObjs) 
        # self.edinames= list(self.edinames)
    
        self._logging.info (
            'Collectiong edilfiles from <%s>'% self.__class__.__name__)
        try:
            func.show_quick_edi_stats(rfiles, self.ediObjs)
        except: pass 
    
        sta , lat,lon, elev , freq , zz, zz_err, rho,\
            rho_err, phs , phs_err =[[] for ii in range (11)]
            
        #--> set gereference attributes 
        for  edi_obj in self.ediObjs :
            sta.append(edi_obj.Head.dataid ) 
            if edi_obj.Head.long is None : 
                edi_obj.Head.long= edi_obj.DefineMeasurement.reflong 
                self._logging.info (
                    'Longitude of Headid <{0}> has'
                    ' been set from Edi.DefineMeasurenement'.
                        format(edi_obj.Head.dataid))
                
            lon.append(edi_obj.Head.long)
            
            if edi_obj.Head.lat is None : 
                edi_obj.Head.lat= edi_obj.DefineMeasurement.reflat 
                self._logging.info ('Latitude of Headid <{0}> has been '
                                    'set from Edi.DefineMeasurenement'.\
                                        format(edi_obj.Head.dataid))
            lat.append(edi_obj.Head.lat)
            if edi_obj.Head.elev is None : 
                edi_obj.Head.elev = edi_obj.DefineMeasurement.refelev 
                self._logging.info('Elevation of Headid <{0}> has been '
                                   'set from Edi.DefineMeasurenement'.\
                                       format(edi_obj.Head.dataid))
            elev.append(edi_obj.Head.elev)
            freq.append(edi_obj.Z.freq) 
        #-----------------------------    
        # correct_lon_lat 
        lon,*_ = func.scale_position(lon)
        lat,*_ = func.scale_position(lat)
        #------------------------------
        # ---> get impednaces , phase tensor and
        # resistivities values form ediobject
        self._logging.debug('Setting impedances tensor , phases tensor and '
                            ' resisvitivities values from ediobjs.')
        
        zz= [edi_obj.Z.z for edi_obj in self.ediObjs]
        zz_err= [edi_obj.Z.z_err for edi_obj in self.ediObjs]

        rho= [edi_obj.Z.resistivity for edi_obj in self.ediObjs]
        rho_err= [edi_obj.Z.resistivity_err for edi_obj in self.ediObjs]

        phs= [edi_obj.Z.phase for edi_obj in self.ediObjs]
        phs_err= [edi_obj.Z.phase_err for edi_obj in self.ediObjs]
        
        # Create the station ids 
        self.id = func.make_ids(self.ediObjs, prefix='S')
        #self.id = ['S{0:03}'.format(ii) for ii in range(len(self.ediObjs))]

        self.longitude, self.latitude, self.elevation = lon , lat, elev
        
        # get frequency array from the first value of edifiles.
        self.freq_array = freq[0] 

        #---> set ino dictionnary the impdance and phase Tensor 
        self._logging.debug('Setting impedances, resistivities and'
                            ' phase tensors into dictionnaries.')
        #---> set attribute on dictionnary
        self._z = {key:value for key , value in zip (self.id, zz)}
        self._z_err ={key:value for key , value in zip (self.id, zz_err)}

        self._res = {key:value for key , value in zip (self.id, rho)} 
        self._res_err ={key:value for key , value in zip (self.id, rho_err)}

        self._phs ={key:value for key , value in zip (self.id, phs)}
        self._phs_err ={key:value for key , value in zip (self.id, phs_err)}
        
    #---> set a few pertinent atttributes thought components xy and yx 
    
    @property 
    def res_xy (self): 
        return {stn :res[: , 0, 1] for stn, res in self._res.items()}
    @property 
    def res_yx (self): 
        return {stn :res[: , 1, 0] for stn, res in self._res.items()}

    @property 
    def res_err_xy (self): 
        return {stn :res_err[: , 0, 1] 
                for stn, res_err in self._res_err.items()}
    @property 
    def res_err_yx (self): 
        return {stn :res_err[: , 1, 0] 
                for stn, res_err in self._res_err.items()}
    
    ###TODO #######################
    @property 
    def res_xx (self): 
        return {stn :res[: , 0, 0] for stn, res in self._res.items()}
    @property 
    def res_yy (self): 
        return {stn :res[: , 1, 1] for stn, res in self._res.items()}

    @property 
    def res_err_xx (self): 
        return {stn :res_err[: , 0, 0] 
                for stn, res_err in self._res_err.items()}
    @property 
    def res_err_yy (self): 
        return {stn :res_err[: , 1, 1] 
                for stn, res_err in self._res_err.items()}
    #################################################################
    @property 
    def z_xy (self): 
        return {stn :z[: , 0, 1] for stn, z in self._z.items()}
    @property 
    def z_yx (self): 
        return {stn :z[: , 1, 0] for stn, z in self._z.items()}

    @property 
    def z_xx (self): 
        return {stn :z[: , 0, 0] for stn, z in self._z.items()}
    @property 
    def z_yy (self): 
        return {stn :z[: , 1, 1] for stn, z in self._z.items()}
    @property 
    def z_err_xx (self): 
        return {stn :z_err[: , 0, 0] for stn, z_err in self._z_err.items()}
    @property 
    def z_err_yy (self): 
        return {stn :z_err[: , 1, 1] for stn, z_err in self._z_err.items()}
    
    @property 
    def z_err_xy (self): 
        return {stn :z_err[: , 0, 1] for stn, z_err in self._z_err.items()}
    @property 
    def z_err_yx (self): 
        return {stn :z_err[: , 1, 0] for stn, z_err in self._z_err.items()}
    
    @property 
    def phs_xy (self): 
        return {stn :phs[: , 0, 1] for stn, phs in self._phs.items()}
    @property 
    def phs_yx (self): 
        return {stn :phs[: , 1, 0] for stn, phs in self._phs.items()}

    
    @property 
    def phs_err_xy (self): 
        return {stn :phs_err[: , 0, 1] 
                for stn, phs_err in self._phs_err.items()}
    @property 
    def phs_err_yx (self): 
        return {stn :phs_err[: , 1, 0] 
                for stn, phs_err in self._phs_err.items()}
    ### TODO #########################
    @property 
    def phs_xx (self): 
        return {stn :phs[: , 0, 0] for stn, phs in self._phs.items()}
    @property 
    def phs_yy (self): 
        return {stn :phs[: , 1, 1] for stn, phs in self._phs.items()}

    
    @property 
    def phs_err_xx (self): 
        return {stn :phs_err[: , 0, 0] 
                for stn, phs_err in self._phs_err.items()}
    @property 
    def phs_err_yy (self): 
        return {stn :phs_err[: , 1, 1] 
                for stn, phs_err in self._phs_err.items()}
   
    
    def rewrite_edis (self, ediObjs=None,  by = 'name' , prefix = None, 
                      dataid =None, savepath = None, how='py', 
                      correct_ll=True, make_coords =False, reflong=None, 
                      reflat=None, step ='1km', edi_prefix =None, **kws): 
        
        """ Rewrite Edis. 
        
        Can rename the dataid,  customize sites and correct the positioning
        latitudes and longitudes. 
        
        Parameters 
        ------------
  
        ediObjs: list 
            Collection of edi object from pycsamt.core.edi.Edi 
       
        dataid: list 
            list of ids to  rename the existing EDI-dataid from  
            :class:`Head.dataid`. If given, it should match the length of 
            the collections of `ediObjs`. A ValueError will occurs if the 
            length of ids provided is out of the range of the number of EDis
            objects 

        by: str 
            Rename according to the inner module Id. Can be ``name``, ``id``, 
            ``number``. Default is ``name``. If :attr:`~.Edi_collection.survey_name`
            is given, the whole survey name should be used. An other argument 
            is ``ix`` out of number formating. 
            
        prefix: str
            Prefix the number of the site. It could be the abbreviation   
            of the survey area. 

        correct_ll: bool,
            Write the scaled positions( longitude and latitude). Default is
            ``True``. 
            
        make_coords: bool 
            Usefful to hide the real coordinates of the sites by generating 
            a 'fake' coordinates for a specific purposes. When setting to ``True``
            be sure to provide the `reflong` and `reflat` values otherwise and 
            error will occurs. 
            
        reflong: float or string 
            Reference longitude  in degree decimal or in DD:MM:SS for the  
            site considered as the origin of the lamdmark.
            
        reflat: float or string 
            Reference latitude in degree decimal or in DD:MM:SS for the reference  
            site considered as the landmark origin.
            
        step: float or str 
            Offset or the distance of seperation between different sites in meters. 
            If the value is given as string type, except the ``km``, it should be 
            considered as a ``m`` value. Only meters and kilometers are accepables.
            Default value of seperation between the site is ``1km``. 
             
        savepath: str 
            Full path of the save directory. If not given, EDIs  should be 
            outputed in the created directory. 
    
        how: str 
            The way to index the stations. Default is the Python indexing
            i.e. the counting starts by 0. Any other value will start counting 
            the site from 1.
            
        kws: dict 
            Additionnal keyword arguments from `~Edi.write_edifile` and 
            :func:`pycsamt.utils.func_utils.make_ll_coordinates`. 
            
        Examples
        ---------
        >>> from pycsamt.core.edi import Edi_Collection
        >>> edipath = r'/Users/Daniel/Desktop/edi'
        >>> savepath =  r'/Users/Daniel/Desktop/ediout'
        >>> cObjs = Edi_collection (edipath)
        >>> cObjs.rewrite_edis(by='id', edi_prefix ='b1',
        ...                       savepath =savepath)
        
        """
        def replace_reflatlon (  olist , nval, kind ='reflat'):
            """ Replace Definemeaseurement Reflat and Reflong by the interpolated
            values.
            
            :param olist: Old list compoing the read EDI measurement infos.
            :type olist: list 
            :param nval: New reflat or reflong list. Mostly is the DD:MM:SS 
                value interpolated. 
            :param kind: Type of measurement to write. 
            :type kind:str 
            
            :return: List of old element replaced. 
            :rtype: list 
            """
            try : 
                for ii, comp in enumerate (olist):
                    if comp.strip().find(kind)>=0: 
                        olist[ii]= f' {kind}={nval}\n'
                        break 
            except:
                pass
            return olist 
        
        regex = re.compile('\d+', re.IGNORECASE)
        by = str(by).lower() 
        if by.find('survey')>=0 :
            by ='name'
        
        prefix = str(prefix) 
        
        if ediObjs is not None:
           self.ediObjs = ediObjs
           self.id = func.make_ids(self.ediObjs, prefix='S', how= how )
           
        if how !='py': 
            self.id = func.make_ids(self.ediObjs, prefix='S',
                                    cmode =None)  
        if dataid is None: 
            if prefix !='None' : 
                dataid = list(map(lambda s: s.replace('S', prefix), self.id))
                
            elif by =='name': 
                # get the first name of dataId of the EDI ediObjs  and filled
                # the rename dataId. remove the trail'_'  
                name = self.survey_name or  regex.sub(
                    '', self.ediObjs[0].Head.dataid).replace('_', '') 
                # remove prefix )'S' and keep only the digit 
                dataid = list(map(lambda n: name + n, regex.findall(
                    ''.join(self.id)) ))
                
            elif by.find('num')>=0: 
               
               dataid = regex.findall(''.join(self.id))  
               
            elif by =='id': 
                dataid = self.id 
                
            elif by =='ix': 
                dataid = list(map(
                    lambda x: str(int(x)), regex.findall(''.join(self.id))))  
            else :
                dataid = list(map(lambda obj: obj.Head.dataid, self.ediObjs))

        elif dataid is not None: 
            if not np.iterable(dataid): 
                raise ValueError('DataId parameter should be an iterable '
                                 f'object, not {type(dataid).__name__!r}')
            if len(dataid) != len(self.ediObjs): 
                raise ValueError (
                    'DataId length must have the same length with the number'
                    ' of collected EDIs({0}). But {1} {2} given.'.format(
                    len(self.ediObjs), len(dataid),
                    f"{'is' if len(dataid)<=1 else 'are'}"))
       
    
        if make_coords: 
            if (reflong or reflat) is None: 
                raise ValueError('Reflong and reflat params must not be None!')
            self.longitude, self.latitude = func.make_ll_coordinates(
               reflong = reflong, reflat= reflat, nsites= len(self.ediObjs),
               step = step , **kws) 
        # clean the old main Edi section info and 
        # and get the new values
        if correct_ll or make_coords:
            londms,*_ = func.scale_position(self.longitude, todms=True)
            latdms,*_ = func.scale_position(self.latitude, todms=True)
  
        for k, (obj, did) in enumerate(zip(self.ediObjs, dataid)): 
            obj.Head.edi_header = None  
            obj.Head.dataid = did 
            
            if correct_ll or make_coords:
                obj.Head.long = float(self.longitude[k])
                obj.Head.lat = float(self.latitude[k])
                oc = obj.DefineMeasurement.define_measurement
                oc= replace_reflatlon(oc, nval= latdms[k])
                oc= replace_reflatlon(oc, nval= londms[k], kind='reflong')
                obj.DefineMeasurement.define_measurement = oc 
            # Empty the previous MTEMAP infos and 
            # fetch the attribute values newly set.
            obj.MTEMAP.mtemapsectinfo =None 
            obj.MTEMAP.sectid= did 
            obj.write_edifile(savepath = savepath ,
                              new_edifilename = edi_prefix,  **kws)
            
        
class Edi : 
    """
    Ediclass  is for especialy dedicated to .edi files, mainly reading 
    and writingwhich are meant to follow the archaic EDI format put forward 
    by SEG.Can read impedance, Tipper but not spectra. To read spectra format
    please consult MTpy documentation https://mtpy2.readthedocs.io/en/develop/ 
    The Edi class contains a class for each major section of the .edi file.
    
    .. note:: Frequency and components are ordered from highest
            to lowest frequency.
    
    
    Arguments 
    ---------
        **edi_filename** :  str 
                    full path to .edi file to be read 
                    *default* is None. 
                   
    ===================   =====================================================
    Attributes            Description                               
    ===================   =====================================================
    edifile               Full path to edifile                          
    MTEMAP                MT section or EMAP DataSection class, contains
                          basic information  on the data  collected . 
                          Can read MT and EMAP section
    DefineMeasurement     DefineMeasurement class, contains information on how
                          the data was collected.
    Head                  Head class, contains metadata  on where, when, and 
                          who collected the data Information class.
    Info                  contains information on how the data was processed 
                          and how thetransfer functions where estimated.
    Z                     Z class, contains the impedance data
    block_size            number of data in one line.  *Default* is 7
    data_header_com       header string for each of   the data section.         
                          Default is '>!****{0}****!'
    bloc_num_format       string format of data, Default is' 16.6e'           
    _t_comps              components for tipper blocks
    _z_comps              components for impedance blocks
    _res_comps            resistivities components 
    _phs_comps            phase components 
    ===================   =====================================================
    
    
    =====================  ===================================================
    Methods                 Description
    =====================  ===================================================
    read_edi                Reads in an edi file and  populates the associated
                            classes and attributes.
    write_edifile           Writes an .edi file following the EDI format given
                            the apporpriate attributes are filled.  Writes out
                            in impedance and Tipper format.
    _fill_data_array        Fill the impedance blocks,Tipper blocks if exists, 
                            if the.edi file is in EMAP section, read_data compute 
                            data to impedance , rho and phase .
    _write_components       Write MT or EMAP components blocks for data blocks
        _blocks_mt          of the .edi file.
    =====================  ===================================================
   
    :Example:  
         
        >>> import pycsamt.core.edi as csedi
        >>> path =  os.path.join(os.environ["pyCSAMT"], 
        ...                        'data','edi', csa000.edi)
        >>> edihead= csedi.Head.get_header_list_from_edi(
        ...                edi_fn=path)
        >>> info_obj =csedi.Info.get_info_list_from_edi(
        ...                    edi_fn=path)
        >>> edi_obj =csedi.Edi(edi_filename=path )
        >>> edi_obj.write_edifile(new_edifilename=None)
        
    """
    
    def __init__(self, edi_filename=None, verbose =0, **kwargs):
        
        self._logging =csamtpylog.get_csamtpy_logger(self.__class__.__name__)
        
        self.edifile =edi_filename 
        self.Head =Head()
        self.Info=Info()
        self.DefineMeasurement =DefineMeasurement()
        self.MTEMAP =MTEMAP()
        self.Z =MTz.Z()
        self.Tip =MTz.Tipper()

        self.verbose =verbose 

        self.block_size = 6 
        self.data_head_com ='>!****{0}****!\n'
        self.bloc_num_format  =' 15.6e'
        
        self. _z_comps = [zz.lower().strip('>') for \
                          zz in SB._edi if re.match("^>Z", zz) is not None]
                
        self._t_comps = [tt.lower().strip('>') for tt in\
                         SB._edi if  (re.match("^>T", tt)is not None \
                                      and tt.find('EXP')>0) ]
            
        self._res_comps = [res.lower().strip('>') for res in\
                         SB._edi if  re.match("^>RHO", res) is not None  ]
            
        self._phs_comps = [phs.lower().strip('>') for phs in\
                         SB._edi if  re.match("^>PHS", phs) is not None  ]
        
        self.edi_data_sectionline =None 
       
        if self.edifile is not None : 
                self.read_edi ()
            
    @property 
    def lon(self): 
        return self.Head.long or self.DefineMeasurement.reflong 
    @lon.setter 
    def lon (self, longitude): 
        self.Head.long =longitude  
    
    @property 
    def lat(self): 
        return self.Head.lat or self.DefineMeasurement.reflat 
    @lat.setter 
    def lat (self, latitude): 
        self.Head.lat=latitude 
        
    @property 
    def elev(self): 
        return self.Head.elev 
    @elev.setter 
    def elev (self, elevation): 
        self.Head.elev =elevation
        

    def read_edi (self, edifile =None): 
        """
        Read edifile and populate attribute to each data section of edi. 

        :param edifile: full path to edifile 
        :type edifile: str 
        
        """
        self._logging.info ("Reading <{0}> edifile.".format(edifile))
        
        if edifile is not None :
            self.edifile = edifile 
        if self.edifile is None : 
            raise CSex.pyCSAMTError_EDI('NoneType can not read. '
                                        'Please provide at least an edifile.')
        if self.edifile is not None :
            if  not os.path.isfile(self.edifile):
                raise CSex.pyCSAMTError_EDI('Can not find edifile to read.'
                                            ' Please check your path.')
                
        #---> read each section and populate attribute 
        if SB.which_file(filename =self.edifile ) =='edi': 
            self.Head = Head.get_header_list_from_edi(edi_fn=self.edifile )
            self.Info = Info.get_info_list_from_edi(edi_fn=self.edifile)   
            self.DefineMeasurement = \
                DefineMeasurement.get_DefineMeasurement_info(
                    edi_fn=self.edifile)
            self.MTEMAP = MTEMAP.get_mtemap_section_list(edi_fn= self.edifile )
            
            self.edi_data_sectionline = MTEMAP.start_data_lines_num
            
            self._get_specific_comp(edifile = self.edifile) #read data section 
           
   
    def _get_specific_comp (self, edifile=None, data_sect_line =None): 
        """
        Method to get a specific components present on the data set .
        put all components into dictionnary  and data 
         like : {'zxxr' :<list of data zxxr>, 
                    'zxyr':<list of data zxyr>, ....,zyyi':<list data of zyyi>}
        Parameters
        -----------
            * edifile : str  
                    full path to edifile 
            * data_sect_line: int  
                    number of line where data section start.
                    default is None': fill automatycally 
        Returns 
        --------
            dict 
                dictionnary of all component values get on edifiles.
            
        .. note:: data_sect_line parameter  is optional.
        
        """
        
        self._logging.info('Read <get_specific_comp> in edi_data'
                           ' section :edifile :{0}.'.format( edifile))
        
        if edifile is not None : self.edifile = edifile 
        else :raise CSex.pyCSAMTError_EDI(
            'Can not find edifile to read.'
             ' NoneType can not be read. Check, your rigth path.')
        if self.edifile is not None :
            if os.path.isfile(self.edifile) is False :
                raise CSex.pyCSAMTError_EDI(
                    'No edifile file .Please check your edipath.')
                
        # we assume that the file will check above , 
        #just to get to microsec running        
        if SB.which_file(filename =self.edifile , deep=False) =='edi': 
            with open (self.edifile, 'r', encoding ='utf8')  as fedi : 
                edilines =fedi.readlines()
            #get all data section of edi     
            edi_data_section =edilines [self.edi_data_sectionline :] 
            
        # read data section and and put data on list of dict. 
        self.comp_dict , _flag={}, 0 
        
        for ii, datalines in enumerate (edi_data_section) : 
            keylines =datalines.strip()
            if '>!' not in keylines and '>' in keylines :
                edi_lines = keylines[1:].strip().split() # out the '>
                if edi_lines==[] or edi_lines[0] ==''or edi_lines[0] ==None :
                    continue # these assement to be sure that we met a none line  
                compkey =edi_lines[0].lower()
                if compkey in self._z_comps or compkey in self._t_comps or\
                    compkey =='freq' or compkey in self._res_comps or\
                        compkey in self._phs_comps: 
                    _flag =1                # turn on to find data 
                    # create dict of key and ready to fill value on list
                    self.comp_dict [compkey]=[]  
                # turn off to find data     
                else : _flag = 0             
            
            elif _flag ==1 and (keylines.find('>!') < 0  or keylines.find('>') < 0 ) : 
                data_lines =keylines.strip().split()
                for jj, data in enumerate(data_lines): 
                    try : 
                        data_lines [jj]= float(data) 
                        # be sure to set all none value of ** value to 0.0
                        if data_lines [jj] == 1.0e32 : data_lines[jj] =.0 
                    except:
                        data_lines [jj] = 0.0 
                self.comp_dict[compkey].extend(data_lines)
                
        self._fill_data_array(data_dict=self.comp_dict)
                
        
    def _fill_data_array (self, data_dict =None  ): 
        """
        Method to fill Impedance Data bloc and Tipper blocks and Resistivities
        and data blocks if provided User dont need to provided data on dictionnay
         provided that your are sure that data providing respect 
        SEG instructions. If not fill automatically by reading edifile .
        
        :param data_dict: dictionnary of data  < get from reading edifiles > 
                *default* is None ,  fill automatically   
        :type data_dict: dict 
                     
        """
        #frequency arranged to highest to lowest . 
        #if 1 mean is lower to highest 
        flag_freqOrder = 0          
        
        if data_dict is not None : 
            self.comp_dict = data_dict
        elif data_dict is None :
            raise CSex.pyCSAMTError_EDI(
                'None value found. Can not read data')
        
        # get frequency array and initialise z_array and Z_error 
        freq_array = np.array(self.comp_dict['freq'], dtype =np.float)
        z_array , z_error_array =  np.zeros(
            (freq_array.size , 2 , 2),dtype =np.complex),\
                            np.zeros ((freq_array.size , 2 , 2), 
                                      dtype =np.float)
        
                               
        if freq_array [-1]  > freq_array [0]   : 
            self._logging.info (
                'Data originally come to lower to highest frequency.'
                     ' We gonna set to Highest to lower frequency')
            
            # set flag to 1 mean value arerange to lowest to hightest                
            flag_freqOrder  = 1                               
                                                
        for key in self.comp_dict.keys(): 

            if key == 'zxxr' :
                z_array[: , 0 , 0] = np.array(
                    self.comp_dict[key]) + 1j * np.array(
                    self.comp_dict[key[:-1]+'i'] )  
                z_error_array [: , 0 , 0]= np.sqrt(
                    np.array(self.comp_dict[key[:-1]+'.var']))
            if key == 'zxyr' :
                z_array[: , 0 , 1] = np.array(
                    self.comp_dict[key]) + 1j * np.array(
                    self.comp_dict[key[:-1]+'i'])  
                z_error_array [: , 0 , 1]= np.sqrt(
                    np.array(self.comp_dict[key[:-1]+'.var']))
            if key == 'zyxr' :
                z_array[: , 1 , 0] = np.array(
                    self.comp_dict[key] ) + 1j * np.array(
                    self.comp_dict[key[:-1]+'i'])  
                z_error_array [: , 1 , 0]= np.sqrt(
                    np.array(self.comp_dict[key[:-1]+'.var']))
            if key == 'zyyr' :
                z_array[: , 1, 1] = np.array(
                    self.comp_dict[key]) + 1j * np.array(
                    self.comp_dict[key[:-1]+'i'])  
                z_error_array [: , 1 , 1]= np.sqrt(
                    np.array(self.comp_dict[key[:-1]+'.var']))

                                  
    
        if flag_freqOrder == 1 : 
            #-->  return matrice to high to low order 
            z_array  , z_error_array, freq_array  = z_array [::-1],\
                z_error_array[::-1], freq_array[::-1]
        
        #---> call Zmodule and populate attributes 
        self.Z._freq, self.Z._z_err , self.Z._z ,=\
            freq_array ,z_error_array,  z_array 
        
        if 'zrot' in self.comp_dict.keys() :
            self.Z.rotation_angle =np.array (self.comp_dict['zrot'])
            
        else : # fill  zeros value as Zrot 
            self.Z.rotation_angle = np.zeros_like (freq_array, dtype =np.float) 
        self.Z.compute_resistivity_phase() # compute resistivity and phase 
        
        
        _flagTip = False  #
        for key in list(self.comp_dict.keys()): 
            if key in self._t_comps : 
                _flagTip =True 
                break 
        if _flagTip ==True : 
            tip_array =np.zeros ((freq_array.size, 1, 2), dtype =np.complex) 
            tip_error_array =np.zeros ((freq_array.size , 1, 2), dtype =np.float)
            
            if 'trot' in self.comp_dict.keys(): 
                self.Tip.rotation_angle=np.array(self.comp_dict['trot'])
            elif 'zrot' in self.comp_dict.keys(): 
                self.Tip.rotation_angle = self.comp_dict ['zrot']
            else : 
                self.Tipper.rotation_angle = np.zeros_like (freq_array, np.float)
            
            for tkey in self.comp_dict.keys ():
                if tkey == 'txr.exp': 
                    tip_array [:, 0, 0 ] =np.array (
                        self.comp_dict[tkey])+ 1j * np.array(
                            self.comp_dict[key.replace('r', 'i')])# txi.exp
                    tip_array  [:, 0, 1] = np.array (
                        self.comp_dict[tkey.replace('txr', 'tyr')])\
                        +1j * np.array(
                            self.comp_dict [tkey.replace('txr', 'tyi')])
                    tip_error_array [:, 0, 0] = np.sqrt(
                        np.array (self.comp_dict [tkey.replace('txr', 'txvar')]))
                    tip_error_array [:, 0, 1] = np.sqrt(
                        np.array (self.comp_dict [tkey.replace('txr', 'tyvar')]))
                    
            if flag_freqOrder ==1 : 
                tip_error_array, tip_array,  = tip_error_array[::-1], tip_array[::-1] 

                self.Tip._freq ,self.Tip._tipper = freq_array , tip_array
                self.Tip._tipper_err  = tip_error_array 
                self.Tip.compute_amp_phase (), self.Tip.compute_mag_direction()
        

    
    def write_edifile (self, edi_fn=None,  new_edifilename=None, datatype =None , 
                       savepath =None, add_filter_array =None, **kwargs  ): 
        """
        Method to write edifiles from data setting oin attribute of Edi 
        or from existing file. Can write also EMAP data are filled attribute 
        of EDI.

        :new_edifilename: new edifile name .If None , will write edi using 
                        station_name plus type of survey (MT of EMAP) plus 
                        year of writing as< S00_emap.2021.edi> or 
                        <S00_mt.2021.edi>
        :type edifile: str                   
                                               
        :param datatype: type of file , "mt" or "emap" if None , program will 
                        detect which file is provided . If datatype is set , 
                        program will be force to rewrite edi into given format.
        :type datatype: str
        
        :param savepath: path to save edifile. 
                        If None  save to your current work directory
        :type savepath: str
        
        :param add_filter_array: ndarray(nfreq, 2, 2),
                                EDI edifile is EMAP section data ,
                                if add filter is provided , will recompute rho.  
        :type add_filter_ array: str
        
        :returns: new_edifile , full path to edifile 
        :rtype: str 
        
        """
        z_rho_phs_labels =[
            [['zxxr', 'zxxi', 'zxx.var'],
              ['zxyr', 'zxyi', 'zxy.var'],
              ['zyxr', 'zyxi', 'zyx.var'],
              ['zyyr', 'zyyi', 'zyy.var']],[['rhoxx', 'rhoxx.var',
                                             'rhoxx.err', 'rhoxx.fit'],
                                            ['rhoxy','rhoxy.var',
                                             'rhoxy.err', 'rhoxy.fit']],
                                            
        [['phxx','phsxx.var', 'phsxx.err', 'phsxx.fit'],
          ['phxy','phsxy.var', 'phsxy.err', 'phsxy.fit']], 
        
                                                 [['frhoxx','frhoxx.var',
                                                   'frhoxx.err', 'frhoxx.fit'],
                                                  ['frhoxy','frhoxy.var', 
                                                   'frhoxy.err', 'frhoxy.fit']],
                                                 
        [['fphsxx','fphsxx.var', 'fphsxx.err', 'fphsxx.fit'],
         ['fphsxy','fphsxy.var', 'fphsxy.err', 'fphsxy.fit']],
        ]
                                                            
        tip_labels =[['txr.exp', 'txi.exp', 'txvar.exp'],
                     ['tyr.exp', 'tyi.exp', 'tyvar.exp']]
                          
        f=0
        verbose = kwargs.pop('verbose', None)
   
        if verbose is not None: 
            self.verbose = verbose 
        if self.verbose is None: self.verbose =0 
        
        if not hasattr(self, 'savepath'):
            self.savepath =None 
            
        if savepath is not None: 
            self.savepath =savepath 
        
        if edi_fn is not None : 
            self.edifile =edi_fn 
        
        if new_edifilename is not None : 
            try: 
                new_edifilename += '.{0}.edi'
                f=2 
            except: 
                new_edifilename =None 
                f=0
            
        if new_edifilename is None : 
            if self.edifile is not None : new_edifilename = '{0}{1}'.format(
                    'new_', os.path.basename(self.edifile)) 
            else : 
                f = 1
                new_edifilename = '{0}_{1}.{2}.edi'
            
        
        if self.Head.dataid is None : 
            self.read_edi()

        # write info, definemeasurement and mtsection or emapsection 
        edi_header_infolines = self.Head.write_head_info()
        edi_info_infolines =self.Info.write_edi_info()
        edi_definemeasurement_infolines =\
            self.DefineMeasurement.write_define_measurement()
        edi_mtemap_infolines =\
            self.MTEMAP.write_mtemap_section(nfreq=self.Z.freq.size)

        edi_mtemap_infolines.append('\n')
        #---> try to force the program to write either emap or mt section format 
        
        if datatype is None : 
        # if self.typefile is None :
            self.typefile = edi_mtemap_infolines[0][2:-5].lower()
        elif datatype is not None : 
            if re.match(datatype.lower(), 'emapsect') is None \
                and re.match(datatype.lower(), 'mtsect') is None :
                warnings.warn (
                    'Currently <pyCSAMT> can write ">=MTSECT" or '
                     ' ">=EMAPSECT". The only acceptables datatype keys '
                     ' are either "mt" or "emap".')
                raise CSex.pyCSAMTError_EDI(
                    'Datatype provided is not acceptable .'
                    'Please try "mt" or "emap".')

            self.typefile= datatype 
        
        if f==1 :
            new_edifilename = new_edifilename.format(
                self.Head.dataid,
                self.typefile, 
                datetime.datetime.now().year) # set the name of new_filename 
        elif f==2:
               new_edifilename=new_edifilename.format(self.Head.dataid)                              
                                                   

        # write frequency >!****FREQUENCIES****!
        edi_freq_infolines = [self.data_head_com.format('frequencies'.upper())]
        edi_freq_infolines =edi_freq_infolines + self._write_components_blocks(
                                                    edi_datacomp= self.Z.freq,
                                                    comp_key='freq')
        
        # print(self.Z.rotation_angle)
        
        if self.typefile == 'mt' or 'mt' in self.typefile: 
            self.typefile = 'mt'
            # write impedance rotation angle : 
                #>!****IMPEDANCE ROTATION ANGLES****!
            #if rotation angle is =0 , for consistency, 
            #fill  zeros value as Zrot 
            try : 
                if isinstance(self.Z.rotation_angle, float) or \
                    len(self.Z.rotation_angle) ==1:
                    self.Z.rotation_angle =np.full_like(self.Z.freq, 
                                                        self.Z.rotation_angle, 
                                                        dtype=np.float)
                    # self.Z.rotation_angle = np.zeros_like (self.Z.freq,
                    #                                        dtype =np.float)
            except : 
                self._logging.error('Error in "edi-file" while setting '
                                    'rotation angle. ')
            else :
                warnings.warn('Rotation angle is set to {0}'.format(
                    self.Z.rotation_angle[0]))
                self._logging.info('Rotation angle is = {0}.'.format(
                    self.Z.rotation_angle[0]))
                
                
            edi_zrot_infolines =[self.data_head_com.format(
                'impedance rotation angles'.upper())]
            edi_zrot_infolines =   edi_zrot_infolines +\
                self._write_components_blocks(edi_datacomp= self.Z.rotation_angle,
                                                           comp_key = 'zrot')
        if  self.typefile =='emap' or 'emap' in self.typefile :
            self.typefile ='emap'
            
            
        #--> Write impedance data and tipper 
        # for consistency , may replace if exist nan number by np.nan

        self.Z.z =np.nan_to_num(self.Z.z)
        self.Z.z_err  = np.nan_to_num (self.Z.z_err)
        
        self.Z.resistivity = np.nan_to_num(self.Z.resistivity)
        self.Z.phase = np.nan_to_num (self.Z.phase)


        #--->  >!****IMPEDANCES****!
        edi_z_data_infolines = [self.data_head_com.format('impedances'.upper())]
        for ii in range (2):
            for jj in range (2): 
                # if np.all(self.Z.z[:, ii, jj].real ==0.0) \
                #     or np.mean (self.Z.z[:, ii, jj].imag)==0.0 or
                # np.mean(self.Z.z_err[:, ii, jj])==0.  :pass  
                # dont write this none value 

                #else :
                zreal_datalines =self._write_components_blocks(
                    edi_datacomp = self.Z.z[:, ii, jj].real , 
                     comp_key =z_rho_phs_labels [0][ii *2 + jj][0])
                zimag_datalines =self._write_components_blocks(
                    edi_datacomp = self.Z.z[:, ii, jj].imag, 
                    comp_key =z_rho_phs_labels [0][ii *2 + jj][1])
                z_error_values  =self.Z.z_err[:, ii, jj] **2
                zvariance_datalines =self._write_components_blocks(
                    edi_datacomp =z_error_values, 
                    comp_key =z_rho_phs_labels [0][ii *2 + jj][2])
                
                zreal_datalines.append('\n'), zimag_datalines.append('\n')
                zvariance_datalines.append('\n')
                edi_z_data_infolines.extend(zreal_datalines)
                edi_z_data_infolines.extend(zimag_datalines)
                edi_z_data_infolines.extend(zvariance_datalines)

        # write EMAP  :
        #---------------------------------------------------------------                      
        if self.typefile =='emap':
            # define rho and phase array 
            edi_rhophs_infolines =[self.data_head_com.format(
                'Resistivities and phases'.upper())]
            rho = np.zeros ((self.Z.freq.size , 2, 4), dtype=np.float)
            phs = np.zeros ((self.Z.freq.size , 2, 4), dtype=np.float)
            if add_filter_array is not None  :
                frho = np.zeros ((self.Z.freq.size , 2, 4), dtype=np.float)
                
                frho [:, 1 , 0] = add_filter_array[:, 0, 1]
                frho [:, 0 , 1] = add_filter_array[:, 1, 0]
            
            rho [:, 0, 0]  = self.Z.res_xx
            rho [:, 0, 1]  = self.Z.res_err_xx **2
            rho  [:, 0, 2] = self.Z.res_err_xx
            
            rho [:, 1, 0]  = self.Z.res_xy
            rho  [:, 1, 1] = self.Z.res_err_xy**2
            rho [:, 1, 2]  = self.Z.res_err_xy  
            
            phs [:, 0, 0]  = self.Z.phase_xx
            phs [:, 0, 1]  = self.Z.phase_err_xx **2
            phs  [:, 0, 2] = self.Z.phase_err_xx
            
            phs [:, 1, 0]  = self.Z.phase_xy
            phs [:, 1, 1]  = self.Z.phase_err_xy**2 
            phs [:, 1, 2]  = self.Z.phase_err_xy        
                
            # convert np.nan to number 
            rho , phs = np.nan_to_num(rho), np.nan_to_num(phs)
            
            for ii in range(2): 
                for jj in range (4) : 
                    if np.all(rho[:, ii , jj]==0 ) or\
                        np.all (phs[:, ii, jj]==0): continue
                    else : 
                        res_datalines = self._write_components_blocks(
                            edi_datacomp = rho[:, ii, jj],
                            comp_key = z_rho_phs_labels[1][ii][jj])
                        phs_datalines = self._write_components_blocks(
                            edi_datacomp = phs[:, ii, jj],
                            comp_key = z_rho_phs_labels[2][ii][jj])
                        res_datalines.append('\n')
                        phs_datalines.append('\n') 
                        
                        edi_rhophs_infolines.extend(res_datalines)
                        edi_rhophs_infolines.extend(phs_datalines)
                        
                        if add_filter_array is not None : 
                            fres_datalines = self._write_components_blocks(
                                edi_datacomp = frho[:, ii, jj],
                                comp_key = z_rho_phs_labels[3][ii][jj])
                        
                            fres_datalines.append('\n')
                            edi_rhophs_infolines.extend(fres_datalines)

            # delete one more return                 
            edi_z_data_infolines = edi_z_data_infolines[:-1] 

        # write EMAP  pass :***TIPPER  *** 
        #---------------------------------------------------------------
        if np.all(self.Tip.tipper ==0 ) and self.Tip.tipper is not None :
            edi_tip_data_infolines , edi_tip_rot_infolines =[''], ['']
        else : 
             #>!****TIPPER ROTATION ANGLES****!
            try :
                edi_tip_rot_infolines =[self.data_head_com.format(
                    'tipper rotation angles'.upper())]
                if  self.Tip.rotation_angle.dtype ==np.float : 
                    #fill tip_rot_value 
                    tip_rotation_value = np.repeat (
                        self.Tip.rotation_angle , self.Tip.freq.size) 
                else :
                    tip_rotation_value = np.repeat (self.Tip.rotation_angle)
                    
                    edi_tip_rot_infolines = edi_tip_rot_infolines+ \
                        self._write_components_blocks(
                            edi_datacomp= np.array(tip_rotation_value),
                            comp_key = 'trot')
                    # write Tipper data 
                    #  >!****TIPPER PARAMETERS****!
                    edi_tip_data_infolines =[self.data_head_com.format(
                        'tipper parameters'.upper())]
                    for ss in range (2) : 
                        tipreal_lines = self._write_components_blocks(
                            edi_datacomp= self.Tip.tipper[:, 0, jj].real,
                                           comp_key=tip_labels[ss][0])
                        tipimag_lines = self._write_components_blocks(
                            edi_datacomp= self.Tip.tipper[:, 0, jj].imag,
                                      comp_key=tip_labels[ss][1])
                        tipvariance_lines = self._write_components_blocks(
                            edi_datacomp= self.Tip.tipper[:, 0, jj]**2,
                                   comp_key=tip_labels[ss][2])
                        edi_tip_data_infolines.extend( tipreal_lines )
                        edi_tip_data_infolines.extend( tipimag_lines )
                        edi_tip_data_infolines.extend( tipvariance_lines)
            except : 
                edi_tip_data_infolines , edi_tip_rot_infolines =[''], ['']
                
        #WRITE EDI
        write_edilines = edi_header_infolines + edi_info_infolines +\
            edi_definemeasurement_infolines + edi_mtemap_infolines  

        if self.typefile =='mt': 
            for ilines in [edi_freq_infolines , 
                           edi_zrot_infolines,
                           edi_z_data_infolines,
                           edi_tip_rot_infolines, 
                           edi_tip_data_infolines ]: 
                            
                if ilines== [''] : continue 
                else : write_edilines += ilines +['\n'] 
            
        if self.typefile =='emap': 
            for ilines in [edi_freq_infolines , edi_z_data_infolines,
                           edi_rhophs_infolines]:
                write_edilines +=ilines +['\n'] 
        
        write_edilines.append('>END')
        # write file : 
        with open (new_edifilename , 'w+', encoding = 'utf8') as fw : 
            fw.writelines(write_edilines)
            
        self.savepath = func.cpath(self.savepath, '_outputEDI_')

        try : 
            shutil.move(new_edifilename, self.savepath )
        except : 
            pass 

        write_edilines=[]
        
        if self.verbose >0 : 
            print('-'*77)
            print('---> Edifile <{0}> has be successfully written to'
                  ' your savepath.\n---> savepath : <{1}>'.\
                  format(os.path.basename(new_edifilename), self.savepath ))

        return new_edifilename
            
    def _write_components_blocks (self, edi_datacomp , 
                                  comp_key, datatype=None ): 
        """
        Method to write blocks  with data components keys . 
        
        :param edi_datacomp: array of data components,
        :type edi_datacomp: ndarray (nfreq, 2, 2)
 
        :param comp_key: component to write in edifile .
        :type comp_key: str 
        
        :datatype: *mt* or *eamap*,   *default* is "mt"
        :type datatype: str 
        
        :returns: data_list, list of line to write in edifile
        :rtype: list 
        
        .. note:: We assume that comp_key provided is found on the edifile
                before using this method.
            
        
        """
        self._logging.info ('Ready to write edi data component blocks info list !')
        
        block_rot =['ROT=ZROT', 'ROT=TROT', 'ROT=NONE']

        if datatype is None : 
            datatype =self.typefile
        
        if comp_key.lower() =='freq':
            comp_block_line=['>FREQ  //{0}\n'.format(edi_datacomp.size)]
        if datatype.lower() in ["mt", '>=mtsect', 'mtsect'] :
 
            if comp_key.lower() =='zrot' or comp_key.lower()=='trot':
                comp_block_line=['>{0}  //{1}\n'.format(comp_key.upper(), 
                                                        edi_datacomp.size)]
            elif comp_key.lower().find('z') >=0 \
                and comp_key.lower() not in ['zrot', 'trot']:
                comp_block_line = ['>{0} {1}  //{2}\n'.format(
                    comp_key.upper(),block_rot[0],edi_datacomp.size)]
 
            elif comp_key.lower().find('t')>=0 \
                and comp_key.lower()  not in ['zrot', 'trot']:
                comp_block_line = ['>{0} {1}  //{2}\n'.format(
                    comp_key.upper(),block_rot[1],edi_datacomp.size)]
                                                                                                    
            else :
                if comp_key.lower()!='freq':
                    raise CSex.pyCSAMTError_EDI(
                        'Could not write the component key <%s>'% comp_key)
            
        elif datatype in ['emap', '>=emapsect','emapsect']: 
            if comp_key.lower()!='freq': 
                comp_block_line = ['>{0} {1}  //{2}\n'.format(
                    comp_key.upper(), block_rot[2],edi_datacomp.size)]
                                         
            # else : raise CSex.pyCSAMTError_EDI(
            #'Could not write the component key <%s>'% comp_key)
            
        else :#datatype.lower() not in ['mt' , '>=mtsect'] 
        #and datatype.lower() not in ['emap', '>=emapsect']: 
            warnings.warn (
                'Two Types of data can be read either  MTsection or EMAPsection.'
                'The dataType provided <{0}> does not match any of them.'
                 ' We suggest to use "mt" or "emap".as datatype key, if not '
                 'please refer to SEG-Edile ''write principles.'.format(datatype)) 
      
            raise CSex.pyCSAMTError_EDI(
                'DataType <{0}> provided is wrong!'
                 ' please use "MT"or "EMAP".'.format(datatype))
        
        #-- > read value form component edi_data_comp 
        for index_data , datacomp in enumerate (edi_datacomp , 1) : 
 
            if datacomp == 0.0 and comp_key.lower() not in ['zrot', 'trot']: 
                datacomp = 1.0E32          
                # if datacomp.mean() == self.Head.empty :
                    #continue  dont write this data when all are empty
                # else : 
            format_value ="{0:{1}}".format(datacomp, self.bloc_num_format) 
            format_value =format_value.upper()          #   1.000000E+32
            if index_data % self.block_size == 0 :
                format_value =format_value + '\n'
            # ---> when block is at the end 
            if index_data == edi_datacomp.size : format_value +='\n' 
            #join block value to its corresponding      
            comp_block_line.append(format_value)  
            
        return comp_block_line
    
    @property 
    def station (self): 
        return self.Head.dataid 

    @property 
    def processingsoftware (self): 
        return self.Info.Processing.ProcessingSoftware.name 
    
    
    @station.setter 
    def station(self, set_station_name): 
        
        if not isinstance(set_station_name,str):
            set_station_name = str(set_station_name)
        self.Head.dataid = set_station_name
        self.MTEMAP.sectid = set_station_name
    
    @processingsoftware.setter 
    def processingsoftware (self, sofware_name): 
        self.Info.Processing.ProcessingSoftware.name = sofware_name
    

    def interpolateZ(self, new_freq_array, interp_type='slinear',
                    bounds_error=True, period_buffer=None):
        """
        Interpolate the impedance tensor onto different frequencies.

        :param new_freq_array: a 1-d array of frequencies to interpolate on
                               to.  Must be with in the bounds of the existing
                               frequency range, anything outside and an error
                               will occur.
        :type new_freq_array: np.ndarray
        :param period_buffer: maximum ratio of a data period and the closest
                              interpolation period. Any points outside this
                              ratio will be excluded from the interpolated
                              impedance array.

        :returns: a new impedance object with the corresponding
                               frequencies and components.
        :rtype: pycsamt.core.z.Z

        :Interpolate: ::

            >>> import pycsamt.core.edi as CSedi
            >>> edi_fn = r"/home/edi_files/cs_01.edi"
            >>> edi_obj = CSedi.Edi(edi_fn)
            >>> # create a new frequency range to interpolate onto
            >>> new_freq = np.logspace(-3, 3, 24)
            >>> new_z_object= edi_obj.interpolate(new_freq)
            >>> edi_obj.write_new_edifile(new_edi_fn=r"/home/edi_files/cs_01_interp.edi",
            >>> ...                   new_Z_obj=new_z_object,
            >>> ...                   )

        """
        
        # if the interpolation module has not been loaded return
        if func.interp_import is False:
            raise ImportError('could not interpolate, need to install scipy')

        # make sure the input is a numpy array
        if not isinstance(new_freq_array, np.ndarray):
            new_freq_array = np.array(new_freq_array)
            
        new_freq_array = np.around (new_freq_array, 2)  
        
        if period_buffer is not None:
            if 0. < period_buffer < 1.:
                period_buffer += 1.
                print("Warning: period buffer must be > 1. Updating to",
                      period_buffer)

        # check the bounds of the new frequency array
        if bounds_error:
            # logger.debug("new freq array %s", new_freq_array)
            if self.Z.freq.min() > new_freq_array.min():
                raise ValueError(
                    'New frequency minimum of {0:.5g}'.format(new_freq_array.min()) + \
                    ' is smaller than old frequency minimum of {0:.5g}'.format(
                        self.Z.freq.min()) + \
                    '.  The new frequency range needs to be within the ' +
                    'bounds of the old one.')
            if self.Z.freq.max() < new_freq_array.max():
                
                raise ValueError(
                    'New frequency maximum of {0:.5g}'.format(new_freq_array.max()) + \
                    ' is larger than old frequency maximum of {0:.5g}'.format(
                        self.Z.freq.max()) + \
                    '.  The new frequency range needs to be within the ' +
                    'bounds of the old one.')

        # make a new Z object
        new_Z = MTz.Z(z_array=np.zeros((new_freq_array.shape[0], 2, 2),
                                       dtype='complex'),
                      z_err_array=np.zeros((new_freq_array.shape[0], 2, 2)),
                      freq=new_freq_array)
        # interpolate the impedance tensor
        for ii in range(2):
            for jj in range(2):
                # need to look out for zeros in the impedance
                # get the indicies of non-zero components
                nz_index = np.nonzero(self.Z.z[:, ii, jj])

                if len(nz_index[0]) == 0:
                    continue
                # get the non-zero components
                z_real = self.Z.z[nz_index, ii, jj].real
                z_imag = self.Z.z[nz_index, ii, jj].imag
                z_err = self.Z.z_err[nz_index, ii, jj]

                # get the frequencies of non-zero components
                f = self.Z.freq[nz_index]

                # get frequencies to interpolate on to, making sure the
                # bounds are with in non-zero components
                new_nz_index = np.where((new_freq_array >= f.min()) & 
                                        (new_freq_array <= f.max()))[0]
                new_f = new_freq_array[new_nz_index]
                
                # apply period buffer
                if type(period_buffer) in [float, int]:
                    new_f_update = []
                    new_nz_index_update = []
                    for ifidx,ifreq in enumerate(new_f):
                        # find nearest data period
                        difference = np.abs(np.log10(ifreq) - np.log10(f))
                        fidx = np.where(difference == np.amin(difference))[0][0]
                        if max(f[fidx] / ifreq, ifreq / f[fidx]) < period_buffer:
                            new_f_update.append(ifreq)
                            new_nz_index_update.append(new_nz_index[ifidx])
                    new_f = np.array(new_f_update)
                    new_nz_index = np.array(new_nz_index_update)

                # create a function that does 1d interpolation
                z_func_real = func.spi.interp1d(f, z_real, kind=interp_type)
                z_func_imag = func.spi.interp1d(f, z_imag, kind=interp_type)
                z_func_err = func.spi.interp1d(f, z_err, kind=interp_type)

                # interpolate onto new frequency range
                new_Z.z[new_nz_index, ii, jj] = z_func_real(
                    new_f) + 1j * z_func_imag(new_f)
                new_Z.z_err[new_nz_index, ii, jj] = z_func_err(new_f)
                
        # compute resistivity and phase for new Z object
        new_Z.compute_resistivity_phase()

        return new_Z 
        
    # --> write new edi file
    def write_new_edifile(self, new_edi_fn=None, new_Z=None,**kws):
        """
        write a new edi file if things have changed.  Note if new_Z is not 
        None, they are not changed in `Edi` object, you  need to change them 
        manually if you want them to be changed. Similarly, the new function 
        name does not change the `Edi` object `edi_filename` attribute but does
        change Edi.edi_object.edi_filename attribute.

        :param edi_fn: full path to new edi file
        :type new_edi_fn: string

        :param new_Z: new Z object
        :type new_Z: pycsamt.core.z.Z

        :param new_Tipper: new Tipper object
        :type new_Tipper: pycsamt.core.z.Tipper

        :returns edi_fn: full path to edi file written
        :rtype edi_fn: string
        """

        # get header information, mostly from site
        edi_obj = Edi(edi_filename=self.edifile) 
        
        if new_Z is not None:
            edi_obj.Z = new_Z
        else:
            edi_obj.Z = self._Z

        # --> write edi file
        edi_fn = edi_obj.write_edifile(new_edifilename= new_edi_fn, **kws)

        return edi_fn
    
     
class Head (object): 
    """
    The edi head block contains a series of options which (1) identity the data  
    set, (2) describe whn , where and by whoom was acquired , and (3) describe   
    when , how  and by whom it was written.
    
    Arguments 
    ----------
        **edi_fn** :str 
            Full path to edi path 
        
    > HEAD
     - DATAID="kap012"
     - ACQBY="Phoenix"
     - FILEBY="EMTF FCU"
     - FILEDATE=01/02/18
     - LAT=-30:52:05.62
     - LONG=21:44:35.00
     - ELEV=1166
     - STDVERS=SEG 1.0
     - PROGVERS="4.0"
     - PROGDATE=06/20/11
     - MAXSECT=999
     - EMPTY=1.0e+32
                
    ==============  ===============================  =============  ===========
    Attribute       Description                      Restriction    Default 
    ==============  ===============================  =============  ===========
    dataid          Identifier for data set             str         Required 
    acqby           Name of contractor or               str 
                    otherparty                                      Required 
    fileby          Name of contractor of other 
                    party                               str         Required 
    acqdate         Date of (start of) data
                    acquisition                         date        Required 
    enddate         Date of end of data acq             date        ""
    filedate        Date EDI was written                Date        Required 
    country         Name of country of acq.             str         ""
    state           state(s) of province(s) of 
                    acquisition                         str         ""
    county          Name of country of acq.             str         ""
    prospect        Name of associated prospect         str         ""
    loc             Description of location             str         ""
    lat             avg.(approx) latitude of acq.       str         ""
                            
    long            avg.(approx)longitude of avq.       str         ""
    elev            avg.(approx)elevation of acq.       str         ""
    units           Units for elevation                 "m"|"ft"    "m"
    stdvers         Version of EDI format for 
                    this file                           str         Required 
    progvers        Version ID for prog. written        str         Required 
    progdate        Last revision of prog writing
                    file                                str         Required 
    coordsys        coordinate system 
                    [geographic|geomagnetic]            str         Geog.North
    declination     geomagnetic declination             float       "10."     
    maxsect         Maximum data section in EDI 
                    file                                int>=1      "16"
    bindata         if not "", tag for binary data      str of ""   ""
                    file                        
    empty           Value which represents"nodata"      float       "1.0E32"
    ==============  ===============================  =============  ===========
    
    """
    head_keys =['dataid', 'acqby',
                'fileby','acqdate',
                'enddate', 'filedate',
                        'country', 'state', 
                        'county','prospect', 
                        'loc', 'lat', 
                                'long', 'elev','declination','datum', 
                                'units', 'stdvers', 'coordsys',
                                        'progvers', 'progdate',
                                        'maxsect', 'bindata','project',  
                                        'survey','empty']
    
    
    
    
    def __init__(self, edi_header_list=None , **kwargs):
        
        self.logging =csamtpylog.get_csamtpy_logger(self.__class__.__name__)
        self.Location =Location ()
        self.dataid =None 
        self.acqby =None 
        self.fileby =None 
        self.acqdate =None
        self.enddate =None 
        self.filedate =datetime.datetime.utcnow().strftime(
                        '%Y/%m/%d %H:%M:%S UTC')
        self.country =None 
        self.state =None 
        self.county =None 
        self.prospect =None 
        self.loc =None 
        self.units ='m'
        self.stdvers='SEG 1.0'
        self.progvers=f'pyCSAMT {pycsamt.__version__}'[:-2]
        self.progdate =datetime.datetime.utcnow().strftime('%Y/%m/%d')
        self.coordsys ='Geomagnetic North'
        self.declination =None 
        self.datum ='WGS84'
        self.maxsect =None 
        self.bindata =None
        self.project =None
        self.survey =None 
        self.empty =1.0E32
        
        self.edi_header = edi_header_list
        
        for key in list(kwargs.keys()):
            setattr(self, key, kwargs[key])
            
        if self.edi_header is not None : 
            self.read_head()
            
    @property 
    def lat (self): 
        return self.Location.latitude 
    @lat.setter 
    def lat (self, lat):
        try :
            float(lat)
        except :
                self.Location.latitude =gis.convert_position_str2float(lat) 
                self.logging.info (
                    'Converted string "dms" input latitude  to decimal degree.')
        else : 
            self.Location.latitude =float(lat)


    @property 
    def long (self): 
        return self.Location.longitude 
    @long.setter 
    def long (self, long):
        try :float (long)
        except : 
            self.Location.longitude= gis.convert_position_str2float(long)
            self.logging.info (
                'Converted string "dms" input longitude  to decimal degree.')
        else : self.Location.longitude = float(long)

    
    @property 
    def elev(self):
        return self.Location.elevation 
    
    @elev.setter 
    def elev(self, elev): 
        self.Location.elevation =elev

        
    @classmethod 
    def get_header_list_from_edi (cls, edi_fn=None):
        """
        Class method to return edi_head_list .
        
        :paramedi_fn: full path to edifile
        :type edi_fn: str 
        
        """
        _logger.info ('Geting <%s> Head info'% edi_fn)
        
        if edi_fn is None : raise CSex.pyCSAMTError_EDI(
                'Edile file not found ! Please check your right path.')
        markhead , ediheadlist=0,[]
        if SB.which_file(filename= edi_fn, deep=False ) =='edi': 
            if SB.which_file(filename= edi_fn, deep=True) !='edi': 
                raise CSex.pyCSAMTError_Header(
                    'Edifile <{0}> - Header is not correct.'
                    'Please check your edifile'.format(
                        os.apth.basename(edi_fn)))
            
            with open(edi_fn , 'r', encoding ='utf8') as fh : 
                head_lines = fh.readlines()
                for hh, headitems in enumerate (head_lines): 
                    
                    if '>HEAD' in headitems or \
                        re.match(r'>HEAD', headitems) is not None:
                        markhead = hh
   
                    if ('>INFO' in headitems or\
                        re.match(r'>INFO', headitems) is not None) or (
                        ###BUG 
                        # sometimes >INFO data are missing so for safety 
                        # we try to use for stop the next EDI section which 
                        # >=DEFINEMEAS'
                        #re.match(r'>=DEFINEMEAS', headitems) is not None
                        headitems.find('>=DEFINEMEAS')>=0 
                        ):
                        
                        ediheadlist.extend(head_lines[markhead:hh])

                        break 
            return cls(edi_header_list = ediheadlist)

    
    def read_head (self, edi_header_list =None):
        """
        read_header_list and set attributes values
        
        :param edi_header_list: list of edifile header infos 
        :type edi_header_list: list 
        
        :Example:
             
            >>> from pycsamt.core.edi import Head
            >>>> file_edi= 'S00_ss.edi'
            >>>> path =  os.path.join(os.environ["pyCSAMT"], 
            ...                      'pycsamt','data', file_edi)
            >>> edihead= Head.get_header_list_from_edi(edi_fn=path)
            >>> print(edihead.lat)
            >>> print(edihead.long)
            >>> print(edihead.elev)
            >>> print(edihead.acqby)
        """
        self.logging.info ('Reading info list.')
        
        if edi_header_list is not None :
            self.edi_header = edi_header_list
        if self.edi_header is None :
            raise CSex.pyCSAMTError_Header('None items found to read.')
        new_header =[]
        if self.edi_header is not None :
            for ii, item in enumerate(self.edi_header) :
                item= item.strip().replace('"','') # replace any blank line 
                for key in self.head_keys :
                    # some Edifile use LON instead of LONG ,
                    # must take into accountfor parsing
                    if key=='long' : key='lon'  
                    # set attribute that exist in  headkeys  
                    if key.upper() in item:         
                        #for more consistency , 
                        #try to clean possible space before setting attribute 
                        # function return a list of clean value
                        keyi =func._strip_item(
                            item_to_clean=item.lower().split('=')[0])[0] 
                        if keyi=='lon': keyi= 'long'
                        if keyi =='coordsys' :keyi =='coordinate_system'

                        try :
                            value = func._strip_item(
                                item_to_clean=item.split('=')[1])[0]  
                        except: value =''
                        self.__setattr__(keyi.lower(), value)
                        new_header.append(''.join([keyi.upper(), '=', value]))

        self.edi_header = new_header
        
        return self 
    
    def write_head_info (self , head_list_infos =None): 
        """
        Write list info . Can read edi and rewrite list or  to provide 
        input as ['key_01=value_01', 'key_02=value_02', ...,'key_nn=value_nn']
        
        .. Note:: If value is  None ,  don't need to write the key . 

        :param head_list_infos: list , list of head info 
        :type head_list_infos: list 
        
        :returns: write_info list , list ready to write to let
                EDI file more visible .
        :rtype: list 
        
        :Example:
            
            >>> from pycsamt.core.edi import Head
            >>> path =  os.path.join(os.environ["pyCSAMT"],
            ...                         'pycsamt','data', S00_ss.edi)
            >>> edihead= Head.get_header_list_from_edi(edi_fn=path)
            >>> print(edihead.write_head_info(
            ...    head_list_infos = edihead.edi_header))
        """
        self.logging.info ('Writing Edifile info .')

        write_header =['>HEAD\n']
        if head_list_infos is not None : 
            self.edi_header = head_list_infos 
        
        if  self.edi_header is None:
            for key in self.head_keys : 
                keyvalue =getattr(self, key)
                if keyvalue == None :  continue#keyvalue=''
                elif keyvalue !=None : 
                    if key =='long' or key=='lat':
                        keyvalue =gis.convert_position_float2str(keyvalue)
                    if key =='dataid' :
                        keyvalue ='"{0}"'.format(keyvalue) # put ""
                    if key =='stdvers' :
                        keyvalue ='"{0}"'.format(keyvalue)
                    if key =='progvers':
                        keyvalue ='"{0}"'.format(keyvalue)
                        write_header.append(''.join(['  {0}'.format(
                            key.upper()),'=',str(keyvalue),'\n']))
                    else : write_header.append(''.join(['  {0}'.format(
                        key.upper()),'=',str(keyvalue).upper(),'\n']))
                    
        elif self.edi_header is not None : #for consistency rewrite .
            for items in self.edi_header : 
                item, itemvalue= items.strip().split('=')

                if item.lower() in self.head_keys :
                    if  itemvalue in ['None' , ' ', '',None]:pass
                    else :
                        # be sure to let item with no space.
                        item =func._strip_item(item_to_clean=item.lower())[0] 

                        if item =='dataid' : 
                            itemvalue ='"{0}"'.format(itemvalue) # put ""
                        if item =='stdvers' : 
                            itemvalue ='"{0}"'.format(itemvalue)
                        if item =='progvers':

                            itemvalue ='"{0}"'.format(
                                getattr(Software(**{'name':'pyCSAMT'}),'name'))
                            write_header.append(''.join(
                                ['  {0}'.format(item.upper()),'=',
                                 itemvalue,'\n']))
                        
                        else : write_header.append(''.join(
                            ['  {0}'.format(item.upper()),'=',
                             itemvalue.upper(),'\n']))
        
        write_header.append('\n')

        return write_header
    
 
class Info :
    """
    Class EDI info class , collect information of the survey. 

    > INFO
     - MAXINFO=999
     - PROJECT=SAMTEX
     - SURVEY=Kaapvaal 2003
     - YEAR=2003
     - PROCESSEDBY=SAMTEX team
     - PROCESSINGSOFTWARE=JONES 2.3
     - PROCESSINGTAG=
     - SITENAME=South Africa
     - RUNLIST=
     - REMOTEREF=
     - REMOTESITE=
     - SIGNCONVENTION=exp(+ i\omega t)

    ================  ================  =======================================
    Attributes         Type             Explanation
    ================  ================  =======================================
    maxrun              int>=1          maximum number of text lines in 
                                        info text(maybe less)
    Source              class obj       Porvenace of data to rewrite 
    Processing          Processing obj  How data where processed   
    Notes               Note class      info additions 
    ================  ================  =======================================
    
    """
    infokeys =['maxinfo', 'project', 'survey', 'creationdate', 
               'processedby','processingsoftware', 'processingtag','sitename',
                'runlist', 'remoteref', 'remotesite','signconvention']
    
    
    def __init__(self,edi_info_list =None , **kwargs ):
        self.logging =csamtpylog.get_csamtpy_logger(self.__class__.__name__)
        
        self.ediinfo =edi_info_list
        self.filter =None #use for EMEAP section
        self.maxinfo =999 
        self.Source=Source ()
        self.Processing=Processing()
        self.Copyright =Copyright()
        self.Head=Head()
        
        for key in list(kwargs.keys()):
            self.__setattar__(key, kwargs[key])
            
        if self.ediinfo is not None : 
            self.read_info()
                          
            
    @classmethod 
    def get_info_list_from_edi (cls, edi_fn =None):
        """
        Class to get edinfo from edifiles 
          
        :param edi_fn: full path to edifile
        :type edi_fn: str
        
        :returns: edi_info_list 
        :rtype: list 
        
        """
        _logger.info (
            'subclass :Reading <%s> Ediinfo and return  class.'% cls.__name__)
        
        info_mark=0
        if edi_fn is None : 
            raise CSex.pyCSAMTError_EDI(
                'None infos to read! Please provide the right path .')
        # we do this simultaneous assement to get a ms of computation .
        if SB.which_file(filename =edi_fn, deep=False)=='edi': 
            if SB.which_file(filename= edi_fn, deep =True )=='edi': 
                with open(edi_fn, 'r', encoding ='utf8') as fi:
                    edi_lines =fi.readlines()
                    
                for ii, info in enumerate(edi_lines): 
                    if '>info'.upper()in info or info.find('>INF')>= 0  or\
                        re.match(r'>INFO', info) is not None : 
         
                        info_mark=ii 
                    elif info == ['\n'] :continue 
                    elif '>=' in info or re.match(r'>=DEFINEMEAS', info) or\
                        info.find('>=DEFINEMEAS')>=0: 
                        #get the info plus one     
                        list_info = edi_lines[info_mark+1:ii] 
                        break 

            return cls (edi_info_list=list_info)
        
    def read_info (self, edi_info_list =None ):
        """
        readinformation and populate attaribute info 
        can set other attributes once read and not present on the file.
        
        :param edi_info_list: list of infos files 
        :type edi_info_list: list 
        
        """
        listinfo=[]
        
        if edi_info_list is not None : self.ediinfo =edi_info_list 
        if self.ediinfo is None :
            raise CSex.pyCSAMTError_EDI('None list found. can not read.')
        
        for ii , iteminfo in enumerate(self.ediinfo ):
            iteminfo=iteminfo.replace('"','')
            if iteminfo =='\n':continue #continue 
            elif iteminfo !='\n':
                try : 
                    item , infovalue = iteminfo.strip().split('=')
                except : #Try to set attribute for consistency  
                #wheter no value is provide as info value .[''] 
                    if len(iteminfo.strip().split('='))==1 : 
                        item , infovalue= iteminfo.strip().split('=')[0], None
                        pass 
            if item ==''  or infovalue == '' or len(item)==0 or\
                item ==[] or infovalue ==None :continue 
            else : 

                item , infovalue = func._strip_item(
                    item_to_clean=item.lower())[0],\
                    func._strip_item(item_to_clean=infovalue)[0]

            if 'project'  in item  or  'survey' in item or\
                'date' in item or  'sitename'in item:
                self.Source.__setattr__(item, infovalue)
                listinfo.append(''.join([' {0}'.format(item.upper()), 
                                         '=', infovalue]))
            elif 'process' in item or 'run' in item or\
                'remo' in item or 'signconv' in item : 
                if item =='processingsoftware':
                    self.Processing.ProcessingSoftware.__setattr__('name', infovalue)
                    listinfo.append(''.join([' {0}'.format(
                        item.upper()), '=', infovalue]))
                else : 
                    self.Processing.__setattr__(item, infovalue)
                    listinfo.append(''.join([' {0}'.format(
                        item.upper()), '=', infovalue]))
            elif 'creating_application' in item : 
                if self.Processing.ProcessingSoftware.name is None : 
                     self.Processing.ProcessingSoftware.__setattr__('name', infovalue)
            else :
                self.__setattr__(item, infovalue)
                listinfo.append(''.join([' {0}'.format(
                    item.upper()), '=', infovalue]))
        
        self.ediinfo = listinfo
        
        return self 
    
            
    def write_edi_info (self, edi_info_list =None ): 
        """
        Write edi information info . Can read edi and rewrite list or to provide 
        input as ['key_01=value_01', 'key_02=value_02', ...,'key_nn=value_nn' ]
        Note : If value is absent i.e None ,  don't write the key . Info write 
        method add somefield notes informations from other softwares if exists.  
 
        :param edi_info_list: list of infos contain in info sections 
        :type edi_info_list: list 
        
        :returns: list of info 
        :rtype: list 
            
        :Example:
            
            >>> from pycsamt.core.edi import Info 
            >>> file_edi_2='SAMTEX.edi_2.edi'
            >>> file_edi= 'S00_ss.edi'
            >>> path =  os.path.join(os.environ["pyCSAMT"], 
            ...                         'pycsamt','data', file_edi_2)
            >>> info_obj =Info.get_info_list_from_edi(edi_fn=path)
            >>> print(info_obj.write_edi_info())
            
        """
        self.logging.info (
            'writing ediinfo from <%s> class.' % self.__class__.__name__)
        
        write_info =['>INFO\n']

        if edi_info_list is not None : self.ediinfo = edi_info_list
        
        if self.ediinfo is None : 
            for key in self.infokeys :
                if key in ['project' , 'survey' ,'creationdate',  'sitename' ]: 
                    valueinfo =getattr(self.Source, key)
                elif 'process' in key or 'run' in key or\
                    'rem'in key or 'signconv' in key : 
                    if key == 'processingsoftware': 
                        valueinfo = getattr(self.Processing.ProcessingSoftware, 'name')
                    elif key =='processedby': 
                        valueinfo = getattr(self.Processing, 'processedby')
                        
                    else :  valueinfo = getattr(self.Processing, key)
                else : valueinfo = getattr(self, key)
                
                if key == 'processedby': 
                    write_info.append(''.join(['  {0}'.format(key.upper()), 
                         '=', '"{0}"'.format(str(valueinfo)), '\n']))
                else :
                    write_info.append(''.join(['  {0}'.format(key.upper()),
                                                  '=', str(valueinfo).upper(), '\n']))
            #add softwares info 
            write_info.append(''.join(['  {0}'.format('CREATINGSOFTWARE'), '=',
                                       getattr(self.Source, 'creatingsoftware'), '\n']))
            # write_info.append(''.join(['  {0}'.format('CREATIONDATE'), '=',
            #                            getattr(self.Source, 'creationdate'), '\n']))
            if self.filter is not None : 
                write_info.append(''.join(['  {0}'.format('FILTER'), '=',
                                           str(getattr(self, 'filter')).upper(), '\n']))
                                       
            
        elif self.ediinfo is not None : #rewrite for consistency
            
            for infok in self.infokeys : 
                if infok in ['project' , 'survey' , 'creationdate' ,  'sitename' ]: 
                    valueinf =getattr(self.Source, infok)
                    if infok == 'creationdate' :
                        valueinf = getattr(self.Source, infok)
                elif 'process' in infok or 'run' in infok or\
                    'rem'in infok or 'signconv' in infok :
                    if infok == 'processedby' :
                        valueinf = '"{0}"'.format(self.Source.creatingsoftware)    
                    elif infok == 'processingsoftware' :
                        valueinf ='"{0}"'.format( 
                            getattr(self.Processing.ProcessingSoftware, 'name'))
                    else : valueinf =getattr(self.Processing, infok)
                elif  self.filter is not None : 
                    write_info.append(''.join(
                        ['  {0}'.format('FILTER'), '=', 
                         getattr(self, 'filter'), '\n']))
                else : valueinf =getattr(self, infok)
                
                if valueinf is None :continue   
                else : write_info.append(''.join(
                    ['  {0}'.format(infok.upper()), '=', str(valueinf), '\n']))
                
            
            for value in self.ediinfo: #add others info from procesing files.
                keyi, valueinf  = value.strip().split('=')
                keyi = func._strip_item(item_to_clean =keyi)[0]
                valueinf = func._strip_item(item_to_clean =valueinf)[0]
                if keyi.lower() in suit.latitude :
                    write_info.append(''.join(
                        ['  {0}'.format(keyi), '=', valueinf, '\n']))
                elif keyi.lower() in suit.azimuth :
                    write_info.append(''.join(
                        ['  {0}'.format(keyi), '=', valueinf, '\n'])) 
                elif keyi.lower() in suit.longitude :
                    write_info.append(''.join(
                        ['  {0}'.format(keyi), '=', valueinf, '\n']))
                elif keyi.lower() in suit.elevation :
                    write_info.append(''.join(
                        ['  {0}'.format(keyi), '=', valueinf, '\n']))
                elif 'notes' in keyi.lower() : #add some other fieldnotes : 
                    write_info.append(''.join(
                        ['  {0}'.format(keyi), '=', valueinf, '\n']))
        
        write_info.append('\n')     
        return write_info
            
class DefineMeasurement: 
    """
    Begins Measurement definition data section . Defines Location of sensors 
    and parameters pertainning to runs for each measurments . 
    
    Arguments
    -----------
        **param defineMeas_list**:  list 
                            list for define measurement infos  
  
    >=DEFINEMEAS
        - MAXCHAN=7
        - MAXRUN=999
        - MAXMEAS=9999
        - UNITS=M
        - REFTYPE=CART
        - REFLAT=-30:52:05.62
        - REFLONG=21:44:35.00
        - REFELEV=1166
         
     >!****CHANNELS USING ORIGINAL SITE LAYOUT. FOR ROTATIONS SEE ZROT****!
        - >HMEAS ID=1001.001 CHTYPE=HX X=      0.0 Y=      0.0 Z=   0.0 AZM=   0.0
        - >HMEAS ID=1002.001 CHTYPE=HY X=      0.0 Y=      0.0 Z=   0.0 AZM=  90.0
        - >HMEAS ID=1003.001 CHTYPE=HZ X=      0.0 Y=      0.0 Z=   0.0 AZM=   0.0
        - >EMEAS ID=1004.001 CHTYPE=EX X=    -50.0 Y=      0.0 Z=   0.0 X2=     50.0 Y2=      0.0 AZM=   0.0
        - >EMEAS ID=1005.001 CHTYPE=EY X=      0.0 Y=    -50.0 Z=   0.0 X2=      0.0 Y2=     50.0 AZM=  90.0
    
   
    ================  ===================================  =======  ===========
    Attributes        Description                          Default  Restriction 
    ================  ===================================  =======  ===========
    defineMeas_list   list of definemeasurment             None     no
    maxchan           Maximum number of channels measured  None     yes
    maxmeas           Maximum number of measurements       9999     yes
    maxrun            Maximum number of measurement runs   999      yes
    meas_####         Hmeasurement or EmEasurment object   None     yes
                      defining the measurement made  
    refelev           Reference elevation (m)              None     yes
    reflat            Reference latitude                   None     yes
    refloc            Reference location                   None     yes
    reflon            Reference longituted                 None     yes
    reftype           Reference coordinate system          'cart'   yes
    units             Units of length                      m        yes
 
    ================  ===================================  =======  ===========
    
     .. note::  To get the list for define measurement
            it's better to call the classmethod <get_define_measurement_info> to 
            full path to .edi file to read in.
            ...
    """
    
    definemeasurementkeys = ['maxchan', 'maxrun', 'maxmeas',
                            'units', 'reftype', 'reflat',
                            'reflong', 'refelev']
    definemeasurement_comment = '>!***channels using original site layout. For rotations see zrot***!'

    nchan =0
    def __init__(self, defineMeas_list=None , **kwargs) :
        
        self.logging=csamtpylog.get_csamtpy_logger(self.__class__.__name__)
        self.define_measurement=defineMeas_list
        self.maxchan =None  
        self.maxrun =999
        self.maxmeas =9999 
        self.units ='m'
        self._reflat =None 
        self._reflong= None 
        self._refelev =None 
        self.reftype ='CART'

   
        for key in list(kwargs.keys()): 
            self.__setattr__(key, kwargs[key])
            
        if self.define_measurement is not None : 
            self.read_define_measurement()
        
    @property 
    def reflat (self): 
        return self._reflat 
    @reflat.setter 
    def reflat (self, reflat): 
        self._reflat = gis.assert_lat_value(reflat)

    
    @property 
    def reflong (self): 
        return self._reflong
    @reflong.setter 
    def reflong (self, reflong): 

        self._reflong = gis.assert_lon_value(reflong)
        
    @property 
    def refelev (self): 
        return self._refelev
    @refelev.setter 
    def refelev (self, elevation): 
        self._refelev = gis.assert_elevation_value(elevation)

            
    @classmethod 
    def get_DefineMeasurement_info (cls, edi_fn=None):
        """
        Class method to get definemeasurement list.  
        
        :param edi_fn: full path to edifiles.
        :type edi_fn: str 
        
        :returns: new class with infos list 
        :rtype: list 
        
        """
        _logger.info(
            'Reading < %s> file from <%s>' % (os.path.basename(edi_fn), 
                                              cls.__name__))
        
        _flagmeas =0
        dfmeasurementlist = [] 
        if edi_fn is None : raise CSex.pyCSAMTError_EDI(
                'Can not read edifile .None value is found.')
        if SB.which_file(filename = edi_fn) =='edi': 
            with open(edi_fn , 'r', encoding= 'utf8') as fdefm : 
                edilines =fdefm.readlines ()
            for dd , dmvalue in enumerate(edilines) : 
                if '>=definemeas'.upper() in dmvalue or dmvalue.find(
                        '>=DEFINEMEAS')>=0 or\
                    re.match(r'>=DEFINEMEAS',dmvalue) is not None: 
                    _flagmeas= dd
                
                if _flagmeas > 0 : 
                    if dmvalue.find('>=mtsect'.upper())>=0  or dmvalue.find(
                            '>=emapsect'.upper())>=0 :
                        if '>=MTSECT' in dmvalue or '>=EMAPSECT' in dmvalue:
 
                            dfmeasurementlist= edilines[_flagmeas+ 1:dd]
                            break 
   
        return cls (defineMeas_list=dfmeasurementlist) 
    
    
    def read_define_measurement (self, define_measurement_list  =None ):
        """
        readmeasurement inedilist and populate attributes .
        
        :param define_measurement_list:  list of measurement data 
                                can be  [key_01=value_01 , ..., key_xx = value_xx]
                                Emeas and Hmeas will be set on dictionnary
                                and call the class to populate attribute
                                *default* is None 
        :type define_measurement_list: list
        
        :Example:
            
            >>> from pycsamt.core.edi import DefineMeasurement    
            >>> file_edi= 'S00_ss.edi'
            >>> path =  os.path.join(os.environ["pyCSAMT"], 
            ...                         'pycsamt','data', file_edi_2)
            >>> definemeas =DefineMeasurement.get_measurement_info(edi_fn=path)
            >>> print(definemeas.define_measurement)
            >>> print(definemeas.meas_ex)
                
        .. note:: to get measurement_hx or measurement_ex  
                    for instance get attribute <id> of Emeasurement, run the 
                    following script
                    
        :Example:
            
          >>> from pycsamt.core.edi import DefineMeasurement 
          >>> definemeas =DefineMeasurement.get_measurement_info(edi_fn=path)
          >>> print(definemeas.meas_ex.id) 
          ... 1004.
          
        """
        self.logging.info (
            'Reading DefineMeasurement info <%s>' % self.__class__.__name__)
        xmeas_list =[]

        if define_measurement_list is not None :
            self.define_measurement= define_measurement_list
        if self.define_measurement is None :
            raise CSex.pyCSAMTError_EDI(
                'Can not find DefineMeasurment list to read.'
                'Please provide the list of definemeasurment.'
                ' e:g:[key_01=value_01 , ..., key_xx = value_xx] ')
        for keys in self.define_measurement :
            if re.match('r^{0}'.format('hmeas'.upper()), keys) is None or\
                re.match('r^{0}'.format('emeas'.upper()),keys) is None :   
                try : 
                    dfmkey, dfmvalue =keys.strip().split('=')
                except : pass 
                else :
                    if 'maxchan' in dfmkey.lower():
                        if len(dfmvalue) ==0 or dfmvalue =='' : 
                            dfmvalue = '7'
                    if len(dfmvalue) == 0 : dfmvalue == 'None'  

                    dfmkey , dfmvalue = func._strip_item(
                        item_to_clean=dfmkey.lower()) [0],\
                    func._strip_item(item_to_clean=str(dfmvalue)) [0]
                
                    if dfmkey  in self.definemeasurementkeys : 
                        #set attribute and convert to degree decimal lon and lat
                        if "reflat"  in dfmkey: 
                            self.reflat= dfmvalue
                            self.__setattr__(dfmkey, self.reflat)
                            xmeas_list.append(''.join(
                                [' {0}'.format(dfmkey), '=',
                                 str(gis.convert_position_float2str(
                                     self.reflat)), '\n' ]))
                        elif  dfmkey == 'reflon' or dfmkey =='reflong':
                            if dfmkey =='reflon':
                                #add "g" to get reflong defalut attribute
                                dfmkey +='g'          
                            self.reflong = dfmvalue 
                            self.__setattr__(dfmkey, self.reflong)
                            xmeas_list.append(''.join(
                                [' {0}'.format(dfmkey), '=', 
                                 str(gis.convert_position_float2str(
                                     self.reflong)), '\n' ]))
                        elif 'refelev' in dfmkey : 
                            self.refelev= dfmvalue
                            self.__setattr__(dfmkey.lower(), self.refelev)
                            xmeas_list.append(''.join([' {0}'.format(dfmkey), 
                                                       '=', str(self.refelev),
                                                       '\n' ]))
                        else : 
                            self.__setattr__(dfmkey, dfmvalue)
                            xmeas_list.append(''.join([' {0}'.format(dfmkey), 
                                                       '=', dfmvalue, '\n' ]))
      
            #elif re.match('r^{0}'.format('hmeas'.upper()),
                #keys) is not None or re.match('r^{0}'.format(
                    #'emeas'.upper()),keys) is not  None : 
            # set the meas into dictionar of Emas and Hmeas             
            if 'hmeas'.upper() in keys :  
                self.nchan +=1 
                keyH=keys.strip().split()
                keyH= gather_measurement_key_value_with_str_parser (
                    old_measurement_list=keyH)
                hdict={}
                for ss in keyH : 
                    if '=' in ss : 
                        hmeask, hmeasv = ss.split('=')
                        hdict [hmeask.lower()]= hmeasv.lower()
                
                xmeas_list.append(hdict)
            elif 'emeas'.upper() in keys :
                self.nchan +=1        #get number of channels
                keyE =keys.strip().split()
                keyE= gather_measurement_key_value_with_str_parser (
                    old_measurement_list=keyE)
                edict ={}
                for ss in keyE : 
                    if '=' in ss : 
                        emeask, emeasv = ss.split('=')   
                        edict[emeask.lower()] =emeasv.lower()
                xmeas_list.append(edict)
                
        #create attribute Emeas and Hmeas fonction to chtype value
        for dict_meas_eh in xmeas_list :  
            if isinstance(dict_meas_eh, dict): 
                newkey ="meas_{0}".format(dict_meas_eh['chtype'].lower())
                if 'h' in newkey :
                    measvalue = Hmeasurement(**dict_meas_eh)# self.meas_hx etc..
                if 'e' in newkey[5:] :
                    measvalue = Emeasurement(**dict_meas_eh) # self.meas_ex etc..
                try :
                    self.__setattr__(newkey, measvalue) 
                except :
                    AttributeError 
                    pass 
         
                
        self.logging.info (
            'Number of channels "emeas" and "hmeas" found is <%s>.' % self.nchan)
        
   
        self.define_measurement = xmeas_list
        
        return self 

    def write_define_measurement(self, define_measurement_list = None): 
        """
        Write definemeasurement method ,intend to write and rewrite measurements 
        infos into list. informations must be on list as possible if not may 
        set attribute manually i.e 
        [key01=value02 , ..., keynn=valuenn + dictXX ]dictXX ={meas_e}, {meas_hx},
            {meas_ey}{meas_hx}, {meas_hy} {meas_hz}

        :param define_measurement_list: list of define measuremenent 
        :type define_measurement_list: list 
        
        :returns: new list of define_measurement
        :rtype: list 
        
        . .notes :: If no edifiles is provided , can write definemeasurement 
                 by creating dict of Eand H measurment. 
            
        :Example:

            >>> ex_dict ={'id':1002, 'chtype':'Ex', 'x':0, 'y':0, 'z':0, 'x2':-50,'y2':0, 'z2':0, 
            ...                   'acqchan':0,'filter':'hanning', 'sensor':'ex', 'gain':None}
            >>> ey_dict ={'id':1003.1, 'chtype':'Ey', 'x':0, 'y':50, 'z':0, 'x2':0,'y2':0, 'z2':0, 
            ...               'acqchan':0}
            >>> hy_dict ={'id':1003, 'chtype':'Hy', 'x':0, 'y':50, 'z':90, 'azm':0.0,'dip':35,
            ...                   'acqchan':'hy','filter':'Hanning', 'sensor':None, 'gain':0, 'measdate':''}
            >>> definemeas =DefineMeasurement()
            >>> ex =Emeasurement(**ex_dict)
            >>> ey =Emeasurement(**ey_dict)
            >>> hy=Hmeasurement(**hy_dict)
            >>> definemeas.__setattr__('meas_ex', ex)
            >>> definemeas.__setattr__('meas_ey', ey)
            >>> definemeas.__setattr__('meas_hy', hy)
            >>> print(definemeas.write_define_measurement())
            
        """
        write_dfmeasurements, dictlist=['>=DEFINEMEAS\n'],[]
        
        self.logging.info (
            'write defineMeasurement infos <%s>.' % self.__class__.__name__)
        
        if define_measurement_list is not None : 
            self.define_measurement =define_measurement_list 

        if self.define_measurement is not None :
            if minimum_parser_to_write_edi(
                    edilines =self.define_measurement) is None : 
                warnings.warn(
                    '.edi_parser provided is wrong.Definemeasurement '
                    'list provided cannot be write.Please provide key values'
                    '  as :[key01=value02 , ..., keynn=valuenn + dictXX ]'
                    'with dictXX ={meas_e}, {meas_hx}, {meas_ey}{meas_hx},'
                    ' {meas_hy} {meas_hz}')
                
                self.logging.warn(
                    'Please check your definemeasurement '
                    'list provided. Can not be written.')
                raise CSex.pyCSAMTError_EDI(
                    "Can not write the items on the definemeasurement "
                     "list provided. It seems no right parser is found.")
            
            for ii, itemkey in enumerate(self.define_measurement):
                if not isinstance(itemkey, dict) :
                   keydfmeas , valuedfmeas = itemkey.strip().split("=")
                   keydfmeas , valuedfmeas =func._strip_item(
                       item_to_clean=keydfmeas)[0],func._strip_item(
                           item_to_clean=valuedfmeas)[0] 
                   write_dfmeasurements.append(''.join(
                       ['  {0}'.format(keydfmeas.upper()), '=',
                        valuedfmeas, '\n']))
                elif isinstance(itemkey, dict):
                    dictlist.append(itemkey)
                
            write_dfmeasurements.append(
                self.definemeasurement_comment.upper() +'\n')
            for itemkey in dictlist:
                if 'h' in itemkey['chtype'] :
                    write_dfmeasurements.append('>hmeas'.upper())
                    # valuemeas= getattr(self, 'meas_{}'.format(itemkey['chtype']))
                elif 'e' in itemkey['chtype'] :
                    write_dfmeasurements.append('>emeas'.upper())

                for keym, valm in itemkey.items():
                    if valm in ['' , 'None' , None] :
                        valm= str("{0:>7.2f}".format(.0))
                    write_dfmeasurements.append(''.join(
                        ['  {0}'.format(keym.upper()), '=',
                         '{0:>7}'.format(valm.upper())]))
                write_dfmeasurements.append('\n')
                
        # for consistency , all options to right file are welcome .
        elif self.define_measurement is  None : 
        
            chtype = ['ex', 'ey', 'hx', 'hy', 'hz']
            for keydefm in self.definemeasurementkeys:
                measvalues = getattr(self, keydefm)
                if keydefm =='reflong':
                    measvalues =gis.convert_position_float2str(measvalues)
                if keydefm =='reflat':
                    measvalues =gis.convert_position_float2str(measvalues)
                
                if measvalues  ==None or measvalues=='': 
                    measvalues='' 
                write_dfmeasurements.append(''.join(
                    ['  {0}'.format(keydefm.upper()), 
                     '=', str(measvalues).upper()+'\n']))
            write_dfmeasurements.append(
                self.definemeasurement_comment.upper() +'\n')
            # then write Hmeasurement and Emeasurements throught their attributes.
            for meastype in chtype : 
                if hasattr(self, 'meas_{0}'.format(meastype)):
                    if 'e' in meastype :
                        write_dfmeasurements.append('>emeas'.upper())
                        for eattr in Emeasurement.emeasurementkey :
                            try :tempattr =getattr(self, 'meas_{0}'.format(meastype))
                            except :pass 
                            else : 
                                
                                valuemeasxx = tempattr.__dict__["{0}".format(eattr)]
    
                                if valuemeasxx ==None : valuemeasxx=''
                                write_dfmeasurements.append(''.join(
                                    ['  {0}'.format(eattr.upper()), '=',
                                     '{0:>10}'.format(str(valuemeasxx).upper())]))
                        write_dfmeasurements.append('\n')
                    elif 'h' in meastype :
                        write_dfmeasurements.append('>hmeas'.upper())
                        for hattr in Hmeasurement.hmeasurementkey : 
                            try:
                                tempattr =getattr(self, 'meas_{0}'.format(meastype))
                            except :
                                AttributeError 
                                pass 
                            else : 
                                valuemeasyy = tempattr.__dict__["{0}".format(hattr)]
                                if valuemeasyy ==None : valuemeasyy =''
                                write_dfmeasurements.append(''.join(
                                    ['  {0}'.format(hattr.upper()),
                                     '=', '{0:>10}'.format(str(valuemeasyy).upper())]))
                        write_dfmeasurements.append('\n')
                        
        write_dfmeasurements.append('\n')
        
        return write_dfmeasurements
                             
class Hmeasurement(object): 
    """
    Define the sensor location and orientation , and run parameters for
    a magnetic field measurement.HMeasurement contains metadata for a magnetic 
    field measurement.

    =====================  ====================================================
    Attributes             Description
    =====================  ====================================================
    id                     Measurement ID , channel number 
    chtype                 Type of Hmeasurment[ HX | HY | HZ | RHX | RHY ]
    x                      x (m) north offset from reference sensor 
    y                      y (m) offest from ref sensor 
    azm                    angle of sensor relative to north = 0
    acqchan                name of the channel acquired usually same as chtype
    dip                    dip angle for sensor ("")
    filter                 description of sensor to run ("")
    gain                   gain used for run ("")
    measdate               date of run ("")
    =====================  ====================================================
     
    To fill Metadata, let get a look of this example 
    
    :Example:
        
        >>> import pycsamt.core.edi as csedi
        >>> hmeas_dict = {'id': '1000.3', 'chtype':'hx', 'x':0, 
        ...                  'y':0, 'azm':0, 'sensor':'None'}
        >>> hmeas = csedi.Hmeasurement(**hmeas_dict )
        >>> print(hmeas.chtype)
        >>> print(hmeas.azm)
        >>> print(hmeas.measdate)
        >>> print(hmeas.sensor)
     """
    hmeasurementkey = ['id', 'chtype', 'x', 'y', 'z', 'azm','dip',
                       'acqchan','filter', 'sensor', 'gain', 'measdate']
    
    def __init__(self, **kws) : 
        
       for  hmeaskey in self.hmeasurementkey: #initialise attribute
            if hmeaskey =='measdate' : 
                self.__setattr__(hmeaskey, time.strftime(
                    '%Y-%m-%d  %H:%M:%S', time.gmtime()))
            elif 'x' in hmeaskey or 'y' in hmeaskey or 'z' in hmeaskey :
                self.__setattr__(hmeaskey, '{:<7.2f}'.format(.0))
            else : self.__setattr__(hmeaskey, None)

       for keys in list (kws.keys()): #set attribute coming from dict
            if keys in self.hmeasurementkey: 
                try : self.__setattr__(keys,'{:<10.3f}'.format( float(kws[keys])))
                except :self.__setattr__(keys, kws[keys])  
                
class Emeasurement(object): 
    """
    Define the electrode location , and run parameters for electric field 
    measurement. EMeasurement contains metadata for an 
    electric field measurement.
    
    =====================  ====================================================
    Attributes             Description (Restriction)
    =====================  ====================================================
    id                     Measurement ID , channel number  ('required ')
    chtype                 Type of Hmeasurement[ Ex | Ey ](required)
    x                      x (m) north offset from first electrode (reauired)
    y                      y (m) offest from ref for first electrode(reauired)
    z                      z offset from the ref for first electrode(reauired) 
    x2                     x offset from the 2nd electrode ('required')
    y2                     y offset from the 2nd electrode (required)
    z2                     z offset from the 2nd electrode ("")
    acqchan                name of the channel acquired usually same as chtype
    filter                 description of sensor to run ("")
    gain                   gain used for run ("")
    measdate               date of run ("")
    =====================  ====================================================

    To Fill Metadata , let take this example 
    
    :Example: 
        
        >>> import pycsamt.core.edi as csedi
        >>> emeas_dict = {'id': '1000.4', 'chtype':'ex', 'x':0,
        ...                  'y':0, 'azm':0, 'acqchan':'ex', }
        >>> emeas = csedi.Hmeasurement(**emeas_dict)
        ... print(emeas.chtype)
        ... print(emeas.azm)
        ... print(emeas.measdate)
     """
    emeasurementkey = ['id', 'chtype', 'x', 'y', 'z', 'x2','y2', 'z2', 
                       'acqchan','filter', 'sensor', 'gain', 'measdate']
    
    def __init__(self, **kws) : 
        
       for  emeaskey in self.emeasurementkey: #initialise attribute
           if emeaskey =='measdate' :
               self.__setattr__(emeaskey, time.strftime(
                   '%Y-%m-%d  %H:%M:%S', time.gmtime()))
           elif 'x' in emeaskey or 'y' in emeaskey or 'z' in emeaskey : 
               self.__setattr__(emeaskey,'{:<7.2f}'.format(.0))
           else : self.__setattr__(emeaskey, None)

       for keys in list (kws.keys()):
            if keys in self.emeasurementkey: 
                try : self.__setattr__(keys,'{:<10.3f}'.format( float(kws[keys])))
                except :self.__setattr__(keys, kws[keys])            

class TSeries(object):
    """
    .. _MTpy::`MTpy <https://github.com/MTgeophysics/mtpy>`

        .. Future plan:: We will call MTpy directy for Tseries 
            Begin a times series data section. Defines the set of measurments  
           for which times series data are presented. 
           refer to MTpy software :ref:`MTpy`
    """
    pass 
class Spectra(object):
    """
    .. Future plan:: We will call MTpy directy for Spectra 
        Begin a spectra  data section. Defines the set of measurments 
       for which spectra  data are presented.  
       Refer to :ref:`MTpy<https://github.com/MTgeophysics/mtpy>`
    """
    pass 

class MTEMAP (object): 
    """
    Begins an MT and EMAP data section .Defines the default measurement for MT 
    sounding and Defines the measurments which makeup an EMAP lines. 
    
    ===========================  ==============================================
         >=MTSECT                        >=EMAPSECT
    ===========================  ==============================================
        - SECTID=""                       - SECTID=S00
        - NFREQ=**                        - NCHAN=4
        - HX= 1001.001                    - MAXBLKS=999
        - HY= 1002.001                    - NDIPOLE =47 
        - HZ= 1003.001                    - NFREQ=17
        - EX= 1004.001                    - HX= 0.
        - EY= 1005.001                    - HY= .0
                                          - HZ= NONE
                                          - CHKSUM =None
    ===========================  ==============================================

    :param mt_or_emap_section_list:  mt and emap section can read ediflies by
                                 calling class method <'get_mtemap_section_list'>
    :type mt_or_emap_section_list: list 
            
    ==============  ==========================  =============  ================
    Attributes        Description                  Default       Restriction 
    ==============  ==========================  =============  ================
    sectid          Name of this section           None         str or""
    nfreq           Number of frequencies          required     int >=1
    maxblks         maximum number of blocks 
                    of this section                None         int>=1
    ndipole         Number of dipoles in 
                    the EMAP line                  required  
    type            Descrip. of spatial filter
                    type used                      None         str or ""
    hx              Meas ID for Hx measurement     None         Def Meas Id or "      
    hy              Meas ID for Hy measurement     None         Def Meas Id or "
    hz              Meas ID for Hz measurement     None         Def Meas Id or "
    ex              Meas ID for Ex measurement     None         Def Meas Id or "
    ey              Meas ID for Ey measurement     None         Def Meas Id or "
    rx              Meas ID for Rx measurement     None         Def Meas Id or "
    ry              Meas ID for Ry measurement     None         Def Meas Id or "
    chksum          checksum total for dvalues     None         Num of ""
    ==============  ==========================  =============  ================
    
    .. note ::
        MTEMAP can recognize function to value provided with type of data acquired 
        either MT or EMAP. More attributes can be added by inputing a
        key word dictionary.
        ...
    
    :Example:
        
        >>> MTEMAP(**{'ex':'0011.001', 'hx':'0013.001','sectid':'Niable','nfreq':47, 
        ...             'ey':'0012.001', 'hy':0014.001 , 'hz': 0015.001 })
        >>> MTEMAP (**{'sectid':'yobkro','nfreq':18 , 'maxblks':100, 'hx':"1003.1", 
        ...           'ex':'1005.4','hz':'1006.3', 'ey':1002.3 , 'hy':'1000.2','chksum':47,
        ...           'ndipole':47, 'type':'hann'})
    """
    mtemapsectkey =['sectid', 'nfreq', 'maxblks', 'hx', 'hy', 'hz', 'ex','ey','rx',
                'ry', 'ndipole', 'type', 'chksum']
    
    start_data_lines_num = None 
    temp_sectid = None 
    def __init__(self, mt_or_emap_section_list =None, **kwargs): 
        
        self._logging =csamtpylog.get_csamtpy_logger(self.__class__.__name__)
        self.mtemapsectinfo = mt_or_emap_section_list
        for keys in self.mtemapsectkey: 
            self.__setattr__(keys, None) 
            
        for keys in list(kwargs.keys()):setattr(self, keys, kwargs[keys])
        if self.mtemapsectinfo is not None:
            self.read_mtemap_section()
            
        
    @classmethod 
    def get_mtemap_section_list (cls, edi_fn =None ): 
        """
        MT or EMAP section_classmethod to get special MT info or EMAP info in edifile

        :param edi_fn: full path to edifile
        :type edi_fn: str 
            
        :returns: newclass contained a list of mtemap infos 
        :rtype: list 
        
        """
        _logger.info (
            "Read MT Section on edifile <%s>: %s" % (
                os.path.basename(edi_fn),cls.__name__))
        if edi_fn is None :CSex.pyCSAMTError_EDI(
                'NoneType can not be read. Please provide your right edipath.')
        mtflag,  mtsection = 0, []
        gmtsectmeasurement=[]
        
        #check whether it is an edifile.
        if SB.which_file(filename =edi_fn) =='edi': 
            with open (edi_fn , 'r',encoding ='utf8') as fmt : 
                edilines = fmt.readlines()
            for ii, mtitems in enumerate(edilines) : 
                if 'dataid'.upper()in mtitems or \
                    re.match(r'DATAID',mtitems ) is not None :
                    mitems = mtitems.replace('"','')
                    try :
                        forsectid = mitems.split('=')
                    except : pass 
                    else : cls.temp_sectid = forsectid[-1] 

                if '>=mtsect'.upper() in mtitems or \
                    '>=emapsect'.upper() in mtitems  or\
                        mtitems.find('>=mtsect')>=0 :
                    mtflag =ii
                
                if  (mtflag > 0) and ('>!' in mtitems or\
                                      re.match(r'>!', mtitems) is not None or\
                                          '>FREQ' in mtitems) : 
                    mtsection = edilines[mtflag+1:ii]
                    # get the number where start data section .
                    cls.start_data_lines_num =ii 
                    break 
        # rebuild the list for consistent to be sure to
        # take all params whatever it's messy

            ss = gather_measurement_key_value_with_str_parser(
                old_measurement_list=mtsection)
            gmtsectmeasurement.extend(ss)
 
        return cls (mt_or_emap_section_list =gmtsectmeasurement)


    def read_mtemap_section (self, mt_or_emap_section_list =None): 
        """
        Read mtsection and set attribute . values can be set as key 
            [key01=value01, ..., keynn=valuenn] 
            
        :param mt_or_emap_section_list: `mt` or `emap` section list 
        :type mt_or_emap_section_list: list 
            
        :Example:
            
            >>> import pycsamt.core.edi as csedi 
            >>> mtsection_obj= csedi.MTEMAP.get_mtemap_section_list(edi_fn =path)
            >>> info = mtsection_obj.read_mtemap_section()
            >>> print(mtsection_obj.sectid) 
        """

        mtemapsectionlist, fmt =[],"{:>12}"
        
        self._logging.info ('Reading MTEMAP section and populate attributes')
        
        if mt_or_emap_section_list is not None :
            self.mtsectinfo =mt_or_emap_section_list
        if self.mtemapsectinfo is None : 
            warnings.warn(
                'Nonelist can not be read. Please '
                'provide the right MT section list info.')
            raise CSex.pyCSAMTError_EDI(
                'Nonelist found. cannot read MTsection info.'
                ' Please provide the riht path.')
        

        for ii , mtinfo in enumerate (self.mtemapsectinfo):
            
            try : 
                mtinfo =mtinfo.replace('"','')
                mtkey, mtvalue = mtinfo.strip().split('=')
                if mtvalue =='' or mtvalue =='**' :
                    mtvalue ='{:.2f}'.format(.0)
                
            except :pass 
            else : 
                mtkey, mtvalue =func._strip_item(item_to_clean=mtkey)[0],\
                func._strip_item(item_to_clean=mtvalue)[0]
                self.__setattr__(mtkey.lower(), mtvalue)
                if mtkey.lower()=='sectid' or mtkey.lower()=='nfreq' :
                    fmt ='{:>7}'
                mtemapsectionlist.append(''.join([' {0}'.format(mtkey.lower()),
                                                  '=', fmt.format(str(mtvalue))]))
        
        #check the mtsection and fill it is None or 0.00
        for nn, olditem  in enumerate (mtemapsectionlist): 
            if 'sectid' in olditem :
                newitem = olditem.split('=')
                try :float(newitem[-1])
                except : pass 
                else :
                    # for consistency try to leave return if exists
                    mtemapsectionlist[nn]=''.join(
                        ['sectid','=',self.temp_sectid]).strip()   
                    
        #--> rewrite the list for consistency 
        self._logging.info (
            'FileType with DataType supposed to be written'
            ' is on <%s>_sect format ' % mtemapsectionlist[0])
        
        self.mtemapsectinfo =mtemapsectionlist
        
        return self 
    
    
    def write_mtemap_section (self, mt_or_emap_section_list =None ,
                              nfreq=None):
        """
        Method to write MT or EMAP section into list by providing list as
        ['key01= value01', ..., keyxx=valuexx ]. Method can recognize whether
         edifile provided is MT or EMAP then can read file according to. 

        :param mt_or_emap_section_list:
            list of mt or eamp sectiionvalidate by egal ('=')
        :type mt_or_emap_section_list: list 
            
        :Example: 
            
            >>> from pycsamt.core.edi  import MTEMAP
            >>> mtemapinfo ={'sectid':'yobkro','nfreq':18 , 'maxblks':100, 'hx':"1003.1",
            ...         'ex':'1005.4','hz':'1006.3', 'ey':1002.3 , 
            ...         'hy':'1000.2'}#'chksum':18, 'ndipole':47, 'type':'hann'
            >>> mtemapsection_obj = MTEMAP(**mtinfo)
            >>> writeinfomtemapsect = mtemapsection_obj.write_mtemap_section()
            >>> print(writeinfomtemapsect)
        """ 
        self._logging.info ('Writing MT or EMAP section')
        fmt='{:>7}'
        write_lines, mtemapList =['>=MTSECT\n'],[]
        infoType =None          #flag to check whether edifile is MT or EMAP
        if mt_or_emap_section_list is not None :
            self.mtemapsectinfo = mt_or_emap_section_list
        if self.mtemapsectinfo  is not None : 
            for keyitem in self.mtemapsectinfo: 
                key, itemval = keyitem.strip().split('=')
                if key.lower() in self.mtemapsectkey : 
                    if  key.lower() in ['type', 'ndipole', 'chksum']:
                        if itemval not in ['', 'None',' ', None] :
                            infoType ='>=EMAPSECT\n'
                            write_lines.append("".join(
                                ['  {0:<7}'.format(key.upper()), 
                                 '=',fmt.format(str(itemval).upper()), '\n' ]))  
                        else :continue
                    elif key.lower() in ['ex', 'ey', 'hx', 'hy', 'hz']:
                        fmt='{:>12}'
                        write_lines.append(''.join(
                            ['  {0}'.format(key.upper()),'=',
                             fmt.format(str(itemval).upper()),'\n']))
                    
                    else :
                        if key.lower()=='nfreq': 
                            if nfreq is not None : itemval= nfreq
                        write_lines.append("".join(
                            ['  {0:<7}'.format(key.upper()), '=',
                             fmt.format(str(itemval).upper()), '\n' ]))
        # the case to convert file , set attribute 
        elif self.mtemapsectinfo is None : 

            for keyi in self.mtemapsectkey : 
                mtvalue =getattr(self, keyi)
                if keyi in ['type', 'ndipole', 'chksum']:
                    if mtvalue not in ['', 'None',' ', None] : 
                        infoType = '>=EMAPSECT\n'
                    else :continue
                if mtvalue ==None : mtvalue = '' 
                if keyi =='nfreq': 
                    if nfreq is not None : mtvalue = nfreq
                if keyi  in ['ex', 'ey', 'hx', 'hy', 'hz']:
                    fmt='{:>12}'
                    write_lines .append("".join(['  {0}'.format(
                        keyi.upper()), '=',fmt.format(
                            str(mtvalue).upper()), '\n' ]))
                else : write_lines .append("".join(
                    ['  {0:<}'.format(keyi.upper()), 
                     '=',fmt.format(str(mtvalue).upper()), '\n' ]))
        
        if infoType is None :return  write_lines
        if infoType is not None : 
            write_lines[0]=infoType
            for newitem in write_lines :# reduce unuseless components for EMAP
                if newitem.find('EX')> 0 or newitem.find('EY')> 0 or\
                    newitem.find('HZ') >0  :pass
                else :
                    mtemapList.append(newitem)
                
            return mtemapList


class References(object):
    """
    References information for a citation.

    Holds the following information:
        
    ================  ==========  ===============================
    Attributes         Type        Explanation
    ================  ==========  ===============================
    author            string      Author names
    title             string      Title of article, or publication
    journal           string      Name of journal
    doi               string      DOI number 
    year              int         year published
    ================  ==========  ===============================

    More attributes can be added by inputing a key word dictionary
    
    :Example: 
        
        >>> from pycsamt.core.edi import References
        >>> refobj = References(**{'volume':18, 'pages':'234--214', 
        ...                           'title':'pyCSAMT :python toolbox for CSAMT' ,
        ...                  'journal':'Computers and Geosciences', 
        ...                  'year':'2021', 'author':'DMaryE'})
        >>> print(refobj.journal)
    """
    def __init__(self, **kwargs):
        
        for ref in ['author', 'title', 'journal', 'volume', 'doi', 'year']: 
            setattr(self, ref, None)

        for key in list(kwargs.keys()):
            setattr(self, key, kwargs[key])



class Copyright(object):
    """
    Information of copyright, mainly about how someone else can use these
    data. Be sure to read over the conditions_of_use.

    Holds the following informations:

    =================  ===========  ===========================================
    Attributes         Type         Explanation
    =================  ===========  ===========================================
    references          References  citation of published work 
                                    using these data
    conditions_of_use   string      conditions of use of data used 
                                    for testing program
    release_status      string      release status [ open | public |proprietary]
    =================  ===========  ===========================================

    More attributes can be added by inputing a key word dictionary
    
    :Example:
        
        >>> from pycsamt.core.edi import Copyright 
        >>> Copyright(**{'owner':'University of CSAMT',
        ...                  'contact':'Cagniard'})
    """

    def __init__(self, **kwargs):
        self.References = References()
        self.conditions_of_use = ''.join([
            'All data sets located in the "data" directory <data/> may be',
            ' used to test the performance of the software "pyCSAMT"', 
            'but cannot be used for commercial and distributive purposes',
            '. They can only be used for understanding the program and cannot',
            ' be distributed to a third party. However, user can go to IRIS ',
            'website: http://ds.iris.edu/ds/tags/magnetotelluric-data/ to download ',
            'metadata from other useful tests. All data and metadata from  this ',
            'site are available free of charge and may be copied freely, duplicated',
            ' and further distributed provided this data set is cited as the',
            ' reference.'])
                                          
        self.release_status = None
        self.additional_info = None
        for key in list(kwargs.keys()):
            setattr(self, key, kwargs[key])


class Source(object):
    """
    Information of the file history, how it was made

    Holds the following information:

    =====================  ==========  ========================================
    Attributes             Type        Explanation
    =====================  ==========  ========================================
    project                string      where the project have been done 
    sitename               string      where the survey have been taken place
    creationdate           string      creation time of file YYYY-MM-DD,hh:mm:ss
    creatingsoftware       string      name of program creating the file
    author                 Person      person whom created the file
    submitter              Person      person whom is submitting file for
                                       archiving
    =====================  ==========  ========================================

    More attributes can be added by inputing a key word dictionary
    
    :Example:
        
        >>> from pycsamt.core.edi import Source
        >>> Source(**{'archive':'IRIS', 
        ...             'reprocessed_by':'grad_student'})
    """

    def __init__(self, **kwargs):
        
        for key in ['project', 'survey', 'sitename']: 
            self.__setattr__(key, None)
        self.creationdate= time.strftime('%Y-%M-%D - %H:%M:%S', time.gmtime())
        self.creatingsoftware = 'pyCSAMT'
        self.Author = Person()
        self.Recipient =Person()

        for key in list(kwargs.keys()):
            setattr(self, key, kwargs[key])



class Person(object):
    """
    Information for a person

    Holds the following information:
        
    ================  ==========  =============================================
    Attributes         Type        Explanation
    ================  ==========  =============================================
    email             string      email of person
    name              string      name of person
    organization      string      name of person's organization
    organization_url  string      organizations web address
    ================  ==========  =============================================

    More attributes can be added by inputing a key word dictionary
    
    :Example:
        
        >>> from pycsamt.core.edi import Person
        >>> person =Person(**{'name':'ABA', 'email':'aba@cagniard.res.org',
        ...                  'phone':'00225-0769980706', 
        ...          'organization':'CagniadRES'})
        >>> person.name 
        ... ABA
        >>> person.organisation
        ... CagniardRES
    """

    def __init__(self, **kwargs):
        for key in ['email', 'name', 'organization', 'organization_url']:
            setattr(self, key, None)

        for key in list(kwargs.keys()):
            setattr(self, key, kwargs[key])



class Processing(object):
    """
    Information for a Edi processing

    Holds the following information:

    ===================  =============  =======================================
    Attributes              Type        Explanation
    ===================  =============  =======================================
    ProcessingSoftware      class       Software obj : Input software info.
    processedby             str         name handler of dataprocessing
    processingtag           str         specifictag
    runlist                 list        ---
    remoteref               str         reference point for remoting 
    remotesite              str         reference site name 
    signconvention          str         convention sign provide default
    ===================  =============  =======================================

    More attributes can be added by inputing a key word dictionary
    """
    
    
    def __init__(self, **kwargs):
        
        for key in ['processedby', 'processingtag', 'runlist',
                    'remoteref', 'remotesite']:self.__setattr__(key, None)
        self.ProcessingSoftware = Software()
        self.signconvention = 'exp(+i \omega t)'

        for key in list(kwargs.keys()):
            setattr(self, key, kwargs[key])


class Software(object):
    """
    software info 
    
    Holds the following information:

    ================= =========== ========================
    Attributes         Type        Explanation
    ================= =========== ========================
    name                string      name of software 
    version             string      version of sotware 
    Author              string      Author of software
    release             string      latest version release
    ================= =========== ========================

    More attributes can be added by inputing a key word dictionary
    
    :Example:
        
        >>> from pycsamt.core.edi import Software
        >>> Software(**{'release':'0.11.23'})
    """

    def __init__(self, **kwargs):
        for key in ['name', 'version', 'release'] :
            self.__setattr__(key, None)
        self.Author = Person()

        for key in kwargs:
            setattr(self, key, kwargs[key])  
            
def minimum_parser_to_write_edi (edilines, parser =None ):
    """
    This function validates edifile for writing , string with egal. We assume 
    that dictionnary in list will be for definemeasurment E and H fieds. 

    :param edilines: list of items to parse
    :type edilines: list 
    
    :param parser:  parser edifile section DefineMeasurement, 
                    can be change. *default* is egal (=)    
    :type parser: str 
                
    :Example:  
        
        >>> from pycsamt.core.edi  import DefineMeasurement 
        >>> file_edi= 'S00_ss.edi'
        >>> path =  os.path.join(os.environ["pyCSAMT"], 
        ...                         'pycsamt','data', file_edi_2)
        >>> definemeas =DefineMeasurement.get_measurement_info(edi_fn=path)
        >>> minimparser = minimum_parser_to_write_edi(
        ...    edilines =definemeas.define_measurement)
        >>> print(minimparser)
        
    """
    if parser is None :parser ='='
    elif parser is not None :
        if type(parser) is not str : 
            raise TypeError('Parser must a string value like <"=", ...>')
        
    if not isinstance(edilines,list):
        if isinstance(edilines , tuple) : 
            edilines =list(edilines)
        else :
            raise TypeError(
                '<edilines> argument  must be on list of string with egal')
    for ii, lines in enumerate(edilines) :
        if isinstance(lines, dict):continue 
        elif lines.find(parser) <0 : 
            warnings.warn (
                'None <"="> found on this item<{0}> of '
                ' the edilines list. list can not be parsed.'
            ' Please put egal between key and value '.format(edilines[ii]))
            return None 
        
    return edilines 

def gather_measurement_key_value_with_str_parser (
        old_measurement_list, parser =None): 
    """
    fonction to rebuild xmeasurement list , to solder list with egal.
    In the case where no value is found at  the last item, we will add "None" . 
   
 
    :param old_measurement_list: measurement list to solder
    :type old_measurement_list: list 
                
    :param parser: can be egal or all you want, 
                    *Default* is None mean parser '='. 
    :type parser: str 

    :returns:  list solded with egal like <key=value> 
    :rtype: list 
    
    :Example: 
        
        >>> from pycsamt.core.edi import gather_measurement_key_value_with_str_parser
        >>> measm = [ ['>HMEAS', 'ID=1001.001', 'CHTYPE=HX', 'X=', '0.0', 'Y=', 
        ...               '0.0', 'Z=', '0.0', 'AZM=', '0.0', 'TS='],
        ...                 ['>HMEAS', 'ID=1002.001', 'CHTYPE=HY', 'X=', '0.0', 
        ...                     'Y=', '0.0', 'Z=', '0.0', 'AZM=', '90.0'],
        ...                 ['>HMEAS', 'ID=1003.001', 'CHTYPE=HZ', 'X=', '0.0', 'Y=',
        ...                     '0.0', 'Z=', '0.0', 'AZM=', '0.0', 'TS=', '']] 
        >>> for item in measm: 
        ...     print(gather_measurement_key_value_with_str_parser(old_measurement_list=item) )
        ... ['ID=1001.001', 'CHTYPE=HX', 'X=0.0', 'Y=0.0', 'Z=0.0', 'AZM=0.0', 'TS=None']
        ... ['ID=1002.001', 'CHTYPE=HY', 'X=0.0', 'Y=0.0', 'Z=0.0', 'AZM=90.0']
        ... ['ID=1003.001', 'CHTYPE=HZ', 'X=0.0', 'Y=0.0', 'Z=0.0', 'AZM=0.0', 'TS=None']
    """
    # print(old_measurement_list)
    
    new_list, temp=[[] for ii in range(2)]
    if parser is None :parser ='=' 
    for ii, item in enumerate(old_measurement_list): 
        if item.count(parser)> 1 :          # where there is many egal =
            temp.extend(item.split())       #split assume that there is space
            old_measurement_list=old_measurement_list[:ii]
            old_measurement_list.extend(temp)
            break 
        
        # chech whether there is eg: [X=,'']
    for ii, item in enumerate (old_measurement_list):  
        if item.count(parser) ==1 : 
            itemf = item.split(parser)
            if itemf[-1] !=''  or len(itemf[-1]) !=0: new_list.append(item)
            else : 
                try : addvalue = float(old_measurement_list[ii+1])
                except :
                    new_list.append(''.join([item, 'None']))
                    pass
                else : new_list.append(''.join([item, str(addvalue)]))

    return new_list
                        

                 

# if __name__== '__main__': 
#     file_14 ='csi000.dat'
#     file_edi_2='SAMTEX.edi_2.edi'
#     file_edi_3 ='csa000j.edi'
    
#     file_edi_1= 'S00_ssmtpy.edi'
#     file_edi_4 ='testemap2.edi'

#     path =  os.path.join(os.environ["pyCSAMT"], 
#                           'pycsamt','data', file_edi_2)
  
#     # Info.
#     edi_obj =Edi(edi_filename=path)
#     edi_obj.write_edifile(new_edifilename= 'TEST.edi', savepath =pathi, )
  
#     # d=edi_obj.DefineMeasurement.define_measurement
#     # e= edi_obj.MTEMAP.mtemapsectinfo
  
    
  
