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

.. _module-j:: `pycsamt.ff.core.j

   :synopsis: Deal with J files.  The J class can read and write an A.G.Jones file
             file, the 'standard format' of magnetotellurics.  Each section
             of the .edi file is given its own class, so the elements of each
             section are attributes for easy access.
             ...

Created on Thu Dec  3 22:31:16 2020

@author: @Daniel03
"""
import os , re, warnings
import datetime, shutil
import webbrowser
import numpy as np 


import pycsamt.ff.core.cs  as cs_obj
from pycsamt.etc.infos import notion

from pycsamt.ff.core.edi import Edi, Hmeasurement, Emeasurement 
from pycsamt.etc.infos import _sensitive as SB 
import pycsamt.utils.func_utils as func
from pycsamt.utils import exceptions as CSex
from pycsamt.utils._csamtpylog import csamtpylog

_logger = csamtpylog.get_csamtpy_logger(__name__)

class J_collection : 
    """
    J collection Class . Collect jfiles
    
    Arguments 
    ---------
        **list_of_jfiles** : list 
                list of Jfiles collection 
        **survey_name** : str 
                path to A.G. Jones Jfiles 
                
    Holds the following information:
        
    ================  ==========  =============================================
    Attributes         Type        Explanation
    ================  ==========  =============================================
    azimuth           ndarray,1     array_of all stations azimuth (deg) 
    latitude          ndarray,1     stations latitudes in degre decimal 
    longitude         ndarray,1     stations longitudes in degree decimal 
    elevation         ndarray,1     stations elevations in meter.
    period            list          list of array of collection periods.
                                     collect  all periods  from all stations . 
    app_rho           list          list of arrays of apparent resistivities
                                    (ohm.m)
    phase             list          list of array of phases in degree from 
                                    allstations.
    rhomax            list          list of  (rho + 1 standard error) arrays
    rhomin            list          list of (rho - 1 standard error) arrays
    phamax            list          list of (pha + 1 standard error )arrays
    phamin            list          list of (pha - 1 standard error) arrays 
    wrho              list          list of (weight for apparent_resistivity) 
                                    arrays  
    wpha              list          list of (weight for phase ) arrays 
    real              list          list of real part of the transfer
                                    function arrays 
    imag              list          list of imaginary part of the transfer
                                    function
    error             list          list of standard error  from alla stations 
    weight            list          list of  weight (Note: some weights are 
                                    in error!)
    ================  ==========  =============================================
    
    ====================  =====================================================
    Method                          Description 
    ====================  =====================================================
    collect_jfiles          read collections jfiles and populates attributes 
    plot_Topo_Stn_Azim      method to plot Topograpy , station seperation and 
                            Azimuth . 
    rewrite_jfiles          rewrite jfiles. 
    ====================  =====================================================
    
    """
    
    def __init__(self, list_of_jfiles =None , survey_name=None,  **kwargs): 
        self._logging =csamtpylog.get_csamtpy_logger(self.__class__.__name__)
        
        self.jfiles_list =list_of_jfiles
        self.survey_name =survey_name
        
        self.J=J()
        self.Location =cs_obj.Location ()
        self.Site =cs_obj.Site()
        
        self._latitude =None 
        self._longitude =None
        self._azimuth =None 
        self._elevation =None 
        self._station_names = None 
        
        self.nperiods = None 
        self.period =None 
        self.app_rho =None
        self.phase =None 
        
        self.rhomax = None 
        self.rhomin =None 
        self.phamax =None 
        self.phamin = None 
        self.wrho =None 
        self.wpha =None 
        
        self.real =None 
        self.imag =None 
        self.error =None 
        self.weight =None
        self.loc =None 
        self.id =None 
        

        if self.jfiles_list is not None :
            self.collect_jfiles()
        
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
    def azimuth(self): 
        return self._azimuth
    @azimuth.setter 
    def azimuth(self, azimuth): 
        self._azimuth =np.array([float(azz) for azz in azimuth])
        
  
    @property 
    def elevation(self): 
        return self._elevation
    @elevation.setter 
    def elevation(self, elevation): 
        self._elevation =np.array([float(azz) for azz in elevation])
        
    @property 
    def stnames(self):
        return self._station_names
    @stnames.setter 
    def stnames(self, jstnames):
        if isinstance(jstnames, np.ndarray): jstnames=jstnames.tolist()
        elif isinstance(jstnames,list) : self._station_names = J.jname(number_of_sites=len(jstnames), 
                                                                survey_name=self.survey_name)
        elif isinstance(jstnames, int) :self._station_names = J.jname(number_of_sites=jstnames, 
                                                              survey_name=self.survey_name)
        else : 
            try : self._station_names = J.jname(number_of_sites= int(jstnames), 
                                          survey_name=self.survey_name)
            except : raise CSex.pyCSAMTError_J('Stations names must be on list or the number of stations not <{0}>.'.format(type(jstnames)))
    

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
       
            
        if jfn is not None:
            self.jfiles_list = jfn

        if self.jfiles_list is None : 
            raise CSex.pyCSAMTError_J('No files found !'
                                      ' Please provide A.G. J-files ')
        # export to savepath 
        if savepath is None : # create a folder in your current work directory
            try :
                savepath = os.path.join(os.getcwd(), '_outputJ2EDI_')
                if not os.path.isdir(savepath):
                    os.mkdir(savepath)#  mode =0o666)
            except : 
                warnings.warn("It seems the path already exists !")
        
        # call CSAMT obj 
        csamt_jobj= cs_obj.CSAMT(data_fn=jfn)
        
        # create edi-obj and set attributes 
        for ii, stn in enumerate (csamt_jobj.station):
            # create an edi_obj and fill attributes 
            edi_obj=Edi()
            
            # fill Head info 
            edi_obj.Head.dataid= stn # ii mean for all the list 
            edi_obj.Head.acqby = acqby 
       
            if csamt_jobj.jsites_infos is None  or\
                csamt_jobj.jsites_infos==[None]: 
                csamt_jobj.jsites_infos = datetime.datetime.fromtimestamp(
                    os.stat(jfn).st_ctime) # return the creation date of file 
                
            edi_obj.Head.acqdate = csamt_jobj.jsites_infos 
            edi_obj.Head.fileby = fileby 
            edi_obj.Head.filedate = datetime.datetime.now(
                ).strftime('%m-%d-%Y %H:%M:%S')
            
            # get the name of location
            edi_obj.Head.loc = os.path.basename(csamt_jobj._fn)[:-4]

            if prospect is None : prospect ='MTnet'
            setattr(edi_obj, 'Head.prospect', prospect)
    
            if county is not None : county ='MT'
            setattr(edi_obj,'Head.county', county)
            
    
            edi_obj.Head.lat = csamt_jobj.lat[ii]
            edi_obj.Head.long = csamt_jobj.lon[ii] 
            edi_obj.Head.elev = round(csamt_jobj.elev[ii], 2)
        
            edi_obj.Head.maxsect =1000
            
            if project_name is None :
                project_name= os.path.basename(csamt_jobj._fn) 

            #=====>  set EDI OBJ INFOS
            # edi_obj.Info.maxinfo = 999
  
            edi_obj.Info.Source.__setattr__('project', project_name)
            edi_obj.Info.Source.__setattr__('survey',  edi_obj.Head.dataid)
            edi_obj.Info.Source.__setattr__('sitename', edi_obj.Head.dataid)
            edi_obj.Info.Processing.__setattr__('processedby', 'pyCSAMT' )
            edi_obj.Info.Processing.ProcessingSoftware.__setattr__('name', edi_obj.Head.fileby )
    
           
            #====> definemeas 
            edi_obj.DefineMeasurement.maxchan =4
            edi_obj.DefineMeasurement.maxrun = len(csamt_jobj.station)
            edi_obj.DefineMeasurement.__setattr__('reftype' ,'CARTesien')
            edi_obj.DefineMeasurement.__setattr__('reflat',edi_obj.Head.lat  ) 
            edi_obj.DefineMeasurement.__setattr__('reflong', edi_obj.Head.long) 
            edi_obj.DefineMeasurement.__setattr__('refelev',round(edi_obj.Head.elev,2))
            
            #creating xxmeas object 
            codeID_dec = '{0}'.format((ii+1)/edi_obj.Head.maxsect)
            # codeID=  '{0:04}{1}'.format(ii * 10 + 1 , codeID_dec[1:] )
            edi_obj.DefineMeasurement.__setattr__('meas_ex', 
                                                  Emeasurement(**{
                                                      'id':'{0:04}{1}'.format(ii * 10 + 1 , 
                                                       codeID_dec[1:] ), 
                                                        'chtype':'EX', 
                                                        'x': -(dipole_length/2), 
                                                        'y': 0.,
                                                        'x2':dipole_length/2 , 
                                                        'y2':0, 
                                                        }))

            edi_obj.DefineMeasurement.__setattr__('meas_ey',
                                                  Emeasurement(**{
                                                      'id':'{0:04}{1}'.format(ii * 10 + 2 , 
                                                         codeID_dec[1:]), 
                                                        'chtype':'EY', 
                                                        'x':0., 
                                                        'y': -(dipole_length/2),
                                                        'x2':0., 
                                                        'y2':dipole_length/2 , 
                                                            }))
                                                                          

            edi_obj.DefineMeasurement.__setattr__('meas_hx', 
                                                  Hmeasurement(**{
                                                      'id':'{0:04}{1}'.format(ii * 10 + 3 , 
                                                       codeID_dec[1:] ), 
                                                        'chtype':'HX', 
                                                        'x':0., 
                                                        'y': 0.,
                                                        'x2':0., 
                                                        'y2':0. , 
                                                        }))

            edi_obj.DefineMeasurement.__setattr__('meas_hy',
                                                  Hmeasurement(**{
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
            edi_obj.MTEMAP.__setattr__('nfreq', len(csamt_jobj.freq))
            edi_obj.MTEMAP.__setattr__('ex', '{0:04}{1}'.format(ii * 10 + 1 , codeID_dec[1:] ))
            edi_obj.MTEMAP.__setattr__('ey', '{0:04}{1}'.format(ii * 10 + 2 , codeID_dec[1:] ))
            edi_obj.MTEMAP.__setattr__('hx', '{0:04}{1}'.format(ii * 10 + 3 , codeID_dec[1:] ))
            edi_obj.MTEMAP.__setattr__('hy', '{0:04}{1}'.format(ii * 10 + 4 , codeID_dec[1:] ))
            
            #Frequency blocks , impendance and resistivity blocs 
            edi_obj.Z.freq = csamt_jobj.freq 
    
            
            # set phase and resistitivity including error propagation 
            # compute error propagation  
            #-->  initialize ndarray(nfreq, 2, 2) 
            res_array = np.zeros((edi_obj.Z.freq.size, 2,2 ), dtype = np.float)
            res_array_err = np.zeros((edi_obj.Z.freq.size, 2,2 ), dtype = np.float)
            phs_array = np.zeros((edi_obj.Z.freq.size, 2,2 ), dtype = np.float)
            phs_array_err = np.zeros((edi_obj.Z.freq.size, 2,2 ), dtype = np.float)
            
            #dictionnary of components . we set only component into XY . 
            res_array [:, 0 , 1 ] = csamt_jobj.resistivity[stn]  
            res_array_err [:, 0 , 1] = csamt_jobj.resistivity[stn] 
            phs_array[: , 0, 1] = csamt_jobj.phase[stn] 
            phs_array_err  [:, 0, 1]  = csamt_jobj.phase_err[stn] 
        
            #---> Recomputing z with resistivities- phase by using propagrations errors 
            edi_obj.Z.set_res_phase(res_array = res_array, phase_array=phs_array, 
                                    freq=  edi_obj.Z.freq, 
                                    res_err_array=res_array_err,
                                    phase_err_array=phs_array_err)
            
            edi_obj.write_edifile(savepath = savepath)
                
        if len(csamt_jobj.station) > 1: 
            print('-'*77)    
            print('---> {0} wrote sucessfully from j-files.\n---> see path:<{1}> '.\
                  format(len(csamt_jobj.station), savepath))
            print('-'*77)
                    
    
        
    def collect_jfiles (self, list_of_jfiles=None, jpath = None):
        """
        Collect the Jfile from jlist from jpath.Read and populate attributes .
        
        :param list_of_jfiles: list of jfiles
        :type list_of_jfiles: list 
        """
        if list_of_jfiles is None and jpath is not None : 
            self.jfiles_list = [os.path.join(jpath, jfile) 
                                for jfile in os.listdir (jpath)
                                if os.path.isfile(jfile) ==True]
        
        if list_of_jfiles is not None : 
            if isinstance( list_of_jfiles, str) : #we assume to read only one  file 
                self.jfiles_list = [list_of_jfiles] 
            self.jfiles_list = sorted(list_of_jfiles)
            
        elif self.jfiles_list is None : 
            raise CSex.pyCSAMTError_J ('Can not find a list of j files.'
                                       ' Please check your path !')
        
        if isinstance(self.jfiles_list, str): # we assume that only file is read than put on list  for looping.
            self.jfiles_list =[self.jfiles_list]
        
        self._logging.info('Number of J_files collected is %s', len(self.jfiles_list))
        
        if self.jfiles_list is not None : 
            jstations , lon,lat, azim, app_rho, pha, rhomax, rhomin, elev, period, \
            phamax, phamin , wrho, wpha, real , imag , error, weight = [[] for ii in range(18)]
            
            # _set info_dic 
            sinfo, exl, eyl, exaz , eyaz , hxaz , cohmin,cohmax, maxt, nty, weig =[[] for ii in range (11)] 
 
            for jfile in self.jfiles_list : 

                self.J.read_j(j_fn =jfile)
                #----> comments blocks 
                if self.J.read_commentblock is True : 
                    sinfo.append(self.J.jinfo.site_infos),
                    exl.append(self.J.jinfo.ex_length)
                    eyl.append(self.J.jinfo.ey_length)
                    exaz.append(self.J.jinfo.ex_azimuth)
                    eyaz.append(self.J.jinfo.ey_azimuth)
                    cohmin.append(self.J.jinfo.coherence_min)
                    cohmax.append(self.J.jinfo.coherence_max) 
                    maxt.append(self.J.jinfo.max_type), 
                    nty.append(self.J.jinfo.ntype)
                    weig.append(self.J.jinfo.weight)
                
                #--- > information blocks 
                jstations.append(self.J._jnames)
                elev.append(self.J.jelev), period.append(self.J.jperiod)
                lon.append(self.J.jlon), lat.append(self.J.jlat), azim.append(self.J.jazim)
                app_rho.append(self.J.japp_rho),pha.append(self.J.jpha)
                rhomax.append(self.J.jrhomax),rhomin.append(self.J.jrhomin), 
                phamax.append(self.J.jphamax),phamin.append(self.J.jphamin)
                wrho.append(self.J.jwrho),wpha.append(self.J.jwpha)
                real.append(self.J.jreal), imag.append(self.J.jimag)
                error.append(self.J.jerror), weight.append(self.J.jweight)
        
        
        self.stnames , self.elevation = jstations , elev
        self.longitude , self.latitude , self.azimuth = lon ,  lat, azim  
        self.period = period
        self.app_rho, self.rhomax, self.rhomin =  app_rho ,rhomax, rhomin 
        self.phase , self.phamax , self.phamin = pha , phamax , phamin 
        self.wrho , self.wpha =wrho , wpha
        self.weight , self.real , self.imag , self.error = weight , real, imag , error
        
        #-----> set the code name of the site by calling site 
        if self.stnames is not None : 
            self.Site.stn_name = len(self.stnames)
            self.id = self.Site.stn_name 
        
        # set_info
        if self.J.read_commentblock: 
            
            self.J.jinfo.site_infos, self.J.jinfo.ex_length , self.J.jinfo.ey_length = sinfo , exl, eyl
            self.J.jinfo.ex_azimuth, self.J.jinfo.ey_azimuth , self.J.jinfo.hx_azimuth = exaz, eyaz, hxaz 
            self.J.jinfo.coherence_max, self.J.jinfo.coherence_min , self.J.jinfo.max_type = cohmax, cohmin, maxt
            self.J.jinfo.ntype, self.J.jinfo.weight = nty, weig
                
            self.J.jinfo.jconfig ={'Site info': self.J.jinfo.site_infos, 
                                       'Ex length': self.J.jinfo.ex_length, 
                                       'Ey length':self.J.jinfo.ey_length, 
                                       'Ex azimuth': self.J.jinfo.ex_azimuth, 
                                       'Ey azimuth': self.J.jinfo.ey_azimuth, 
                                       'Hx azimuth':self.J.jinfo.hx_azimuth, 
                                                   'Coherence minimum':self.J.jinfo.coherence_min, 
                                                   'Coherence maximum': self.J.jinfo.coherence_max, 
                                                   'Maxtype':self.J.jinfo.max_type, 
                                                   'Ntype':self.J.jinfo.ntype , 
                                                   'Weight':self.J.jinfo.weight,
                                                   }
        
    def plot_Topo_Stn_Azim (self, list_of_jfiles =None, plot ='*',
                            show_stations =False , compute_azimuth =True , 
                savefig =None): 
        """
        plot Topography , Stn _Azim from Jfiles

        :param list_of_jfiles: list of jfiles 
        :type list_of_jfiles: list 
        
        :param plot: type of plot '*' means all the 03 plots or 
                      use `topo`, `stn`, `azim` to plot 
                     individually. *default* is '*'
        :type plot: str or int 
        
        :param show_stations: see the stations names as xlabel .
        :type show_stations: bool  
        
        :param compute_azimuth: if add azimuth, set azimuth to False 
                                *Default* is True 
        :type  compute_azimuth: bool  
        
        :param savefig: PathLike - path to save your figure
        :type savefig: str 
        
        .. note:: 
            Work but not stable ...
        """
        
        if list_of_jfiles is not None : self.jfiles_list= list_of_jfiles 
        if self.jfiles_list is None : raise CSex.pyCSAMTError_J(
                'Can not compute NoneType.Check your path')
        elif self.jfiles_list is not None :
            self.collect_jfiles(list_of_jfiles =self.jfiles_list )
            
        from pycsamt.viewer.plot import Plot1d 
        #--> create Plot1d Obj 
        plot_obj= Plot1d()
        #--> get easting northing arrays after converting from lat, lon
        jeasting , jnorthing = self.Location.get_eastnorth_array_from_latlon(arr_lat=self.latitude, 
                                                                             arr_lon = self.longitude)

        # ----> compute station separation using Profile obj 
        jstn_separation = cs_obj.Profile().stn_separation(easting =jeasting , northing =jnorthing)[0]
       
        #---- > interpolate Jstn_separation  so to get the same size as stations.
        from   scipy.interpolate import interp1d 
        xx_old =np.arange(jstn_separation.size)
        ff =interp1d(x=xx_old, y=jstn_separation, fill_value ='extrapolate')
        xx_new=np.arange(jstn_separation.size+1) # --> add one to get the same size 

        jstn_separation_extrap = ff(xx_new)

        
        if compute_azimuth is True :

            self.azimuth =func.compute_azimuth(easting =jeasting, northing=jnorthing , 
                                               extrapolate =True) 
            # self.Location.azimuth = np.concatenate((jeasting, jnorthing), axis=1)
            # self.azimuth = self.Location.azimuth 

        #---> copmute dipole length if not given 
        dipole_length =jstn_separation.mean()

        station_pk = np.apply_along_axis(lambda xx: xx * dipole_length, 0, np.arange(len(self.stnames)))

        #---> call plot func from plot1D
        plot_obj.plot_topo_sep_azim(plot=plot,
                                    elevation = self.elevation,
                                    station_names= self.stnames ,
                                    station_pk = station_pk , 
                                    azimuth=self.azimuth, 
                                    stn_separation = jstn_separation_extrap, 
                                    set_station_names=show_stations , 
                                    savepath =savefig )
                   

    def rewrite_jfiles (self, list_of_jfiles=None , savepath =None, survey_name =None , 
                        j_extension='.dat' ):
        """
        Method to rewrite A.G. johnson file (J-Format file).
        
        .. seealso::  http://mtnet.info/docs/jformat.txt
        
        :param list_of_files:  collection of Jfiles 
        :type list_of_files: list 
        
        :param survey_name:  name of exploration area 
        :type survey_name:   str  
        
        :param j_extension:  output format,
                            *Default* is '.dat'
        :type j_extension: str  
        
        :param savepath: output directory .If None ,
                        file will be store in current work directory
        :type savepath: str   
        """
        
        if list_of_jfiles is not None : self.jfiles_list =list_of_jfiles 
        if survey_name is not None : self.survey_name =survey_name 
        if self.jfiles_list is  None : 
            raise CSex.pyCSAMTError_J('No files found to read . '
                                      'Please check your path <%s>'% os.getcwd())  
        elif self.jfiles_list is not None : 
            self.collect_jfiles(list_of_jfiles =self.jfiles_list )
         
            
        #--- > start writing 
        write_jlines =[]
        if self.survey_name is None : code_name = 'ks{0:02}-{1}'
        elif self.survey_name is not None : code_name = self.survey_name[:-3]+'{0:02}-{1}'
        codespace= 5*' '
        
        for ii in range (len(self.id)) : 
            #-- > write Head j 
            write_jlines.append(''.join([self.J._comment_mark,'{0:} :'.format(self.J.jinfo.progvers) , 
                                        code_name.format(ii, self.id[ii]), codespace, 
                                        datetime.datetime.now().strftime(self.J.jfmtime), codespace, 'RECS'])+'\n')
            
            #---> write comments line if original file provided it.
            try : 
                if self.J.jinfo.site_infos is not None : write_jlines.append(
                        ''.join([self.J_comment_mark,'{0}:'.format('Site info'),
                                 self.J.jinfo.site_infos,'\n'] ))                                                                  
                if self.J.jinfo.ex_length is not None :
                    write_jlines.append(''.join([self.J_comment_mark,
                                                 '{0}:'.format('Ex length'),
                                                 self.J.jinfo.ex_length,'\n'] )) 
                if self.J.jinfo.ey_length is not None : 
                    write_jlines.append(''.join([self.J_comment_mark,
                                                 '{0}:'.format('Ey length'), 
                                                 self.J.jinfo.ey_length,'\n'] ))
                if self.J.jinfo.ex_azimuth is not None : 
                    write_jlines.append(''.join([self.J_comment_mark,
                                                 '{0}:'.format('Ex azimuth'), 
                                                 self.J.jinfo.ex_azimuth,'\n'] ))                                                                                                                                                   
                if self.J.jinfo.ey_azimuth is not None :
                    write_jlines.append(''.join([self.J_comment_mark,
                                                 '{0}:'.format('Ey azimuth'),
                                                 self.J.jinfo.ey_azimuth,'\n'] ) )                                                                    
                if self.J.jinfo.hx_azimuth is not None : 
                    write_jlines.append(''.join([self.J_comment_mark,
                                                 '{0}:'.format('Hx azimuth'),
                                                 self.J.jinfo.hx_azimuth,'\n'] ))                                                                                                                                   
                if self.J.jinfo.coherence_min is not None :
                    write_jlines.append(''.join([self.J_comment_mark,
                                                 '{0}:'.format('Coherence minimum'),
                                                 self.J.jinfo.coherence_min,'\n']))                                                                                                                                
                if self.J.jinfo.coherence_max is not None :
                    write_jlines.append(''.join([self.J_comment_mark,
                                                 '{0}:'.format('Coherence maximum'),
                                                 self.J.jinfo.coherence_max,'\n'] ))                                                                                                                                  
                if self.J.jinfo.max_type is not None : 
                    write_jlines.append(''.join([self.J_comment_mark,
                                                 '{0}:'.format('Maxtype'), 
                                                 self.J.jinfo.max_type,'\n'] ) )                                                            
                if self.J.jinfo.ntpype is not None :
                    write_jlines.append(''.join([self.J_comment_mark,
                                                 '{0}:'.format('Ntype'), 
                                                 self.J.jinfo.ntpype,'\n'] ) )                                                             
                if self.J.jinfo.weight is not None :
                    write_jlines.append(''.join([self.J_comment_mark,
                                                 '{0}:'.format('Ntype'),
                                                 self.J.jinfo.weight,'\n'] ) )  
            except : pass 
            
            #--- > write jlabels 
            write_jlines.append(''.join(['{0:<13}'.format(
                self.J._jlabels[0]),'=', "{0:>13.1f}".format(self.azimuth[ii]), '\n'] ))
            write_jlines.append(''.join(['{0:<13}'.format(
                self.J._jlabels[1]),'=', "{0:>13.5f}".format(self.latitude[ii]), '\n'] ))
            write_jlines.append(''.join(['{0:<13}'.format(
                self.J._jlabels[2]),'=', "{0:>13.5f}".format(self.longitude[ii]), '\n'] ))
            write_jlines.append(''.join(['{0:<13}'.format(
                self.J._jlabels[3]),'=', "{0:>13.1f}".format(self.elevation[ii]), '\n'] ))
            #---> write jnames , components  and nperiods 
            write_jlines.append(''.join(['{0}'.format(
                self.stnames[ii]), '{0:>7.1f}'.format(self.azimuth[ii]),'\n']))
            write_jlines.append('{0}'.format(self.J.jmode)+'\n')
            write_jlines.append(' {0}'.format(self.J.jperiod.size)+'\n')
            # now  write value :

            for jj in range(self.J.jperiod.size): 
                if re.match(r'^Z+', self.J.jmode) is None or re.match(r'^T+', self.J.jmode) is None:
                    write_jlines.append (''.join( ['{0:<12.3e}'.format(self.period[ii][jj]), 
                                                   '{0:^12.3e}'.format(self.app_rho [ii][jj]), 
                                                   '{0:^12.2e}'.format(self.phase [ii][jj]), 
                                                   '{:^12.3e}'.format(self.rhomax[ii][jj]), 
                                                   '{:^12.3e}'.format(self.rhomin [ii][jj]), 
                                                   '{:^12.2e}'.format(self.phamax[ii][jj]), 
                                                   '{:^12.2e}'.format(self.phamin[ii][jj]),
                                                   '{:^7.2f}'.format(self.wrho[ii][jj]), 
                                                   '{:^7.2f}'.format(self.wpha[ii][jj]),
                                                                     ]))
                elif  re.match(r'^Z+', self.J.jmode) is not None or re.match(r'^T+', self.J.jmode) is not None:
                
                    write_jlines.append (''.join( ['{0:<12.4e}'.format(self.period[ii][jj]), 
                                                   '{0:^12.4e}'.format(self.real [ii][jj]), 
                                                   '{0:^12.4e}'.format(self.imag [ii][jj]), 
                                                   '{:^12.3e}'.format(self.error[ii][jj]), 
                                                   '{:^12.2f}'.format(self.weight [ii][jj]), 
                                                                     ]))
                
                write_jlines.append('\n')
                
            with open('{0}{1}'.format(self.id[ii],j_extension), 'w', encoding='utf8') as fj:
                fj.writelines(write_jlines)
            
                
            if savepath is not None : 
                shutil.move ('{0}{1}'.format(self.id[ii],j_extension), savepath)
                                                           
            write_jlines=[]  
            
        print('-'*77) 
        if savepath is None :savepath =os.getcwd()
        print('---> {0} J-files have been rewritten to <{1}> <----'.format(len(self.id), savepath))
        print('-'*77)
                                                            

class J:
    """
    Class deal with A.G. Jones j-file. see : http://mtnet.info/docs/jformat.txt
    
    Arguments 
    ----------
            **j_fn** : str 
                    path to A.G. Jones Jfile . 
            
    Holds the following information:
        
    ==============  ==============  ===========================================
    Attributes         Type         Explanation
    ==============  ==============  ===========================================
    jazim           float           station azimuth in degree . 
    jlat            float           station latitude in degre decimal 
    jlon            float           station longitude in degree decimal 
    jelev           float           station elevation in meter.
    jmode           str             polarization mode.Polarization mode will 
                                    be use to compute jproperties.
    jnperiod        int             number of period . User can provide.If none 
                                    it will be computed automatically                   
    jperiod         (ndarray,1)     array of periods (s-1) or(1/second)1
    japp_rho        (ndarray,1)     array of apparent resistivities (ohm.m)
    jphase          (ndarray,1)     array of phase in degree 
    jrhomax         (ndarray,1)     rho + 1 standard error
    jrhomin         (ndarray,1)     rho - 1 standard error
    jphamax         (ndarray,1)     pha + 1 standard error 
    jphamin         (ndarray,1)     pha - 1 standard error
    jwrho           (ndarray,1)     weight for apparent_resistivity 
    jwpha           (ndarray,1)     weight for phase 
    jreal           (ndarray,1)     real part of the transfer function
    jimag           (ndarray,1)     imaginary part of the transfer function
    jerror          (ndarray ,1)    standard error 
    jweight         (ndarray,1)     weight (Note: some weights are in error!)
    ==============  ==============  ===========================================
    
    ====================  =====================================================
    Method                          Description 
    ====================  =====================================================
    jname                   staticmethod to generate A.G.Jones station codename
    jMode                   method for comformited H-Polarization and E-Polari-
                            zation mode . if mode is provided , il will check 
                            conformities, if not provied, Default is 'RXY'.
    read_j                  read the jfiles. 
    ====================  =====================================================
    
    .. note ::`-999` means missing data marker
            
    More attributes can be added by inputing a key word dictionary
    
    :Example:
        
        >>> from pycsamt.ff.core import J 
        >>> j=J()
        >>> jmode = j.jMode(polarization_type='RXY')
        >>>  print(jmode)
    """

    jdalpha , jdbeta ={'R':'apparent resistivities and phases', 
                          'S':'as "R" but upward-biased estimates', 
                          'Z':'impedances',
                          'Q':'as "Z" but upward biased estimates', 
                          'C':'"Z" expressed as Schmucker C function', 
                          'T':'GDS transfer functions'}, {'XX': 'xx MT impedance element', 
                                                            'XY':'xy MT impedance element', 
                                                            'YX':'yx MT impedance element',
                                                            'YY':' yy MT impedance element',
                                                            'TE':'TE MT impedance element', 
                                                            'TM':'TM MT impedance element', 
                                                            'AV':'Berdichevsky average', 
                                                            'DE':'determinant average', 
                                                            'ZX':'Tzx transfer function', 
                                                            'ZY':'Tzy transfer function'}

    def __init__(self, j_fn =None , **kws):
        
        self._logging=csamtpylog.get_csamtpy_logger(self.__class__.__name__)

        self.jfn =j_fn 
        
        self._jnames= None 
        self.jlon=None 
        self.jlat=None 
        self.jazim=None 
        self.jelev=None 
        self._jmode =None
        self.jnperiod=None
        self.jinfo =J_infos()

        
        self._jperiod=None
        self._japp_rho =None 
        self._jpha =None 
        self._jrhomax=None 
        self._jrhomin=None 
        self._jphamax =None 
        self._jphamin=None 
        self._jwrho =None  
        self._jwpha =None 
        
        
        self._jreal =None
        self._jimag =None 
        self._jerror =None 
        self._jweight =None 
        self._jhead =None 

        self._jnan = -999

        self._comment_mark='#'
        self.jfmtime ='%m-%d-%Y %H:%M:%S'
        self._jlabels =['>AZIMUTH','>LATITUDE', 
                        '>LONGITUDE', '>ELEVATION' ]
        
        self.read_commentblock =False  
        
        for jkey in list (kws.keys()): 
            setattr(self, jkey, kws[jkey])
        
    #----set jproperties and jsetter-functions -----
        
    @property 
    def jperiod (self): 
        return self._jperiod 
    @jperiod.setter 
    def jperiod(self, jperds): 
        for ii, item in enumerate(jperds): 
            if item in ['','*'] or item == str(self._jnan): jperds[ii]=np.nan
        try : self._jperiod=np.array([float(ii) for ii in jperds])
        except : raise CSex.pyCSAMTError_J('Can not convert "str"jperiod value to float.')
    
    @property 
    def japp_rho (self): 
        return self._app_rho 
    @japp_rho.setter 
    def japp_rho (self, japp_rho): 
        for ii, item in enumerate(japp_rho): 
            if item in ['','*'] or item == str(self._jnan): japp_rho[ii]=np.nan
        try : self._app_rho=np.array([float(ii) for ii in japp_rho])
        except : raise CSex.pyCSAMTError_J('Can not convert "str" apparent resistivity value to float.')
    
    @property 
    def jpha(self): 
        return self._jpha 
    @jpha.setter 
    def jpha (self, jphase): 
        for ii, item in enumerate(jphase): 
            if item in ['','*'] or item == str(self._jnan): jphase[ii]=np.nan
        try : self._jpha=np.array([float(ii) for ii in jphase])
        except : raise CSex.pyCSAMTError_J('Could not convert "str" phase value to float.')
        
    @property 
    def jrhomax(self): 
        return self._jrhomax 
    @jrhomax.setter 
    def jrhomax (self, jrhomax): 
        for ii, item in enumerate(jrhomax): 
            if item in ['','*']or item == str(self._jnan): jrhomax[ii]=np.nan
        try : self._jrhomax=np.array([float(ii) for ii in jrhomax])
        except : raise CSex.pyCSAMTError_J('Could not convert "str"jrhomax value to float.')
        
    @property 
    def jphamax(self): 
        return self._jphamax 
    @jphamax.setter 
    def jphamax (self, jphamax): 
        for ii, item in enumerate(jphamax): 
            if item in  ['','*'] or item == str(self._jnan): jphamax[ii]=np.nan
        try : self._jphamax=np.array([float(ii) for ii in jphamax])
        except : raise CSex.pyCSAMTError_J('Could not convert "str"jphamax value to float.')
        
    @property 
    def jrhomin(self): 
        return self._jrhomin 
    @jrhomin.setter 
    def jrhomin (self, jrhomin): 
        for ii, item in enumerate(jrhomin): 
            if item in ['','*']or item == str(self._jnan): jrhomin[ii]=np.nan
        try : self._jrhomin=np.array([float(ii) for ii in jrhomin])
        except : raise CSex.pyCSAMTError_J('Could not convert "str"jrhomin value to float.')
        
    @property 
    def jphamin(self): 
        return self._jphamin 
    @jphamin.setter 
    def jphamin (self, jphamin): 
        for ii, item in enumerate(jphamin): 
            if item in ['','*'] or item == str(self._jnan): jphamin[ii]=np.nan
        try : self._jphamin=np.array([float(ii) for ii in jphamin])
        except : raise CSex.pyCSAMTError_J('Could not convert "str" jphamin value to float.')
        
    @property 
    def jwrho(self): 
        return self._jwrho
    @jwrho.setter 
    def jwrho (self, jwrho): 
        for ii, item in enumerate(jwrho): 
            if item in  ['','*'] or item == str(self._jnan): jwrho[ii]=np.nan
        try : self._jwrho=np.array([float(ii) for ii in jwrho])
        except : raise CSex.pyCSAMTError_J('Could not convert "str" jwrho value to float.')
    
    @property 
    def jwpha(self): 
        return self._jwpha
    @jwpha.setter 
    def jwpha (self, jwpha): 
        for ii, item in enumerate(jwpha): 
            if item in ['','*'] or item == str(self._jnan): jwpha[ii]=np.nan
        try : self._jwpha=np.array([float(ii) for ii in jwpha])
        except : raise CSex.pyCSAMTError_J('Could not convert "str" jwphase value to float.')
    
    @property 
    def jreal(self): 
        return self._jreal
    @jreal.setter 
    def jreal (self, jreal): 
        for ii, item in enumerate(jreal): 
            if item in ['','*'] or item == str(self._jnan): jreal[ii]=np.nan
        try : self._jreal=np.array([float(ii) for ii in jreal])
        except : raise CSex.pyCSAMTError_J('Could not convert "str" jreal value to float.')   


    @property 
    def jimag(self): 
        return self._jimag
    @jimag.setter 
    def jimag (self, jimag): 
        for ii, item in enumerate(jimag): 
            if item in ['','*'] or item == str(self._jnan): jimag[ii]=np.nan
        try : self._jimag=np.array([float(ii) for ii in jimag])
        except : raise CSex.pyCSAMTError_J('Could not convert "str" jimag value to float.')   

    @property 
    def jerror(self): 
        return self._jerror
    @jerror.setter 
    def jerror (self, jerror): 
        for ii, item in enumerate(jerror): 
            if item in ['','*'] or item == str(self._jnan): jerror[ii]=np.nan
        try : self._jerror=np.array([float(ii) for ii in jerror])
        except : raise CSex.pyCSAMTError_J('Could not convert "str" jerror value to float.')  

    @property 
    def jweight(self): 
        return self._jweight
    @jweight.setter 
    def jweight (self, jweight):
        if jweight is None : self._jweight = np.full((self.japp_rho.size,),1.)
        for ii, item in enumerate(jweight): 
            if item in ['','*'] or item == str(self._jnan): jweight[ii]=np.nan
        try : self._jweight=np.array([float(ii) for ii in jweight])
        except : raise CSex.pyCSAMTError_J('Could not convert "str" jweight value to float.')  


    @property 
    def jmode (self):
        return self._jmode 
    
    @jmode.setter 
    def jmode(self, jpolar):
        if not isinstance(jpolar, str):raise CSex.pyCSAMTError_J('jMode polarization must be on str not <{}>'.format(type(jpolar)))
        self._jmode =self.jMode(polarization_type=jpolar)
        

    @staticmethod
    def jname( number_of_sites, survey_name= None,):
        """
        Staticmethod for generate alphanumeric `station name` (case sensitive 
        when reading A.G Jones files) survey XXX, station 001. 
        
        :param  number_of_sites: number of stations . 
        :type number_of_sites: int 
        
        :param survey_name: place location name . 
        :type survey_name: str 
        
        :returns: list of alphanumeric station names 
        :rtype: list
              
        """
        alpha, cunm ='abcdefghijklmnopqrstuvwxyz', []
        if not isinstance(number_of_sites, int): 
            try :number_of_sites = int(number_of_sites)
            except : raise CSex.pyCSAMTError_J('Number of sites must be int not <{0}>.'.\
                                               format(type(number_of_sites)))
        for ss in range(number_of_sites):
            if survey_name is not None :symb=survey_name[:2].lower() 
            else:  symb='cs'
            if ss == 0 :cunm.append('{0}{1}000'.format(symb, alpha[0])) 
            elif (ss % 10  == 0): 
                cunm.append('{0}{1}{2}'.format(symb, alpha[ss//10], "000"))
            else :cunm.append('{0}{1}{2:03}'.format(symb, alpha[ss//10], ss%10))
        
        return cunm
    
   
    def jMode(self, polarization_type='RXY'):
        """
        Jmode is conformited a different set of H-Polarization ans E-Polarization. 
        for more detail :
            
        .. seealso:: http://mtnet.info/docs/jformat.txt
        
        :Convention:The convention used is for RXY to represent
                the E-polarization (TE) mode, and for RYX
                to represent the B-polarization mode.
        """
        comb =[ ('{0}{1}'.format(jakeys, jbkeys), '{0}/{1}'.format(javalues, jbvalues))\
               for jakeys, javalues  in self.jdalpha.items() for jbkeys, jbvalues in self.jdbeta.items()]
            
        jmodeDICT = {jkeym:jvaluem for jkeym, jvaluem in zip \
                     ([pol for pol, info in comb],[info for pol, info in comb] )}

        if polarization_type not in list(jmodeDICT.keys()) :
            warnings.warn('Polarization Type  provided <{0}> for JMode is unacceptable.'
                          'Please consults the Dictionnary of jPolarisation mode above.')           
            [print('{0:<5}:{1:<55}'.format(keys,values)) for keys,
             values in zip(list(jmodeDICT.keys()), list(jmodeDICT.values()))]
                         
            print(notion.j)
            
            raise CSex.pyCSAMTError_J('Value provided is not in polarization mode '
                                      '.Please consult the dict above.') 
        
        return polarization_type
    
    def read_j (self, j_fn=None  ): 
        """
        Method to read Jformat file.
        
        :param j_fn: path to jfile 
        :type j_fn: str
        """ 

        jdata=[]
        if j_fn is not None : self.jfn =j_fn 
        if self.jfn is None : 
            raise CSex.pyCSAMTError_J('Error file. Please Provide the right path!')
  
        self._logging.info('Reading A.G.Jones J-format "%s"'%os.path.basename(
            self.jfn))
        
        if self.jfn is not None : 

            if SB.which_file(filename=self.jfn) =='j': 
                with open(self.jfn , 'r', encoding='utf8') as fj: 
                    j_data_lines =fj.readlines()
            else :
                self._logging.warn("File <%s>is not J-format file." %self.jfn)
                warnings.warn('File <%s> is not J-Format File.'% self.jfn)
                raise CSex.pyCSAMTError_J('File provided doesn no match the J-format. Please consult :"http://mtnet.info/docs/jformat.txt" {0}'.\
                                                  format(webbrowser.open('http://mtnet.info/docs/jformat.txt')))
            
            #set attribute when loop in the file . 
            for ss , jitems in enumerate(j_data_lines):
                # set comments infos from jfile :
                if re.match(r'^#', jitems) is not None : 

                    
                    j_comminfo, jcommval,*ignore =jitems.strip().split(':')
                    
                    if 'Site info'.lower() in j_comminfo.lower():
                        self.J.jinfo.site_infos = jcommval
                    elif 'Ex length'.lower() in j_comminfo.lower():
                        self.J.jinfo.ex_length = jcommval 
                    elif 'Ey length'.lower() in j_comminfo.lower():
                        self.J.jinfo.ey_length = jcommval
                    elif 'Ex azimuth'.lower() in j_comminfo.lower():
                        self.J.jinfo.ex_azimuth = jcommval 
                    elif 'Ey azimuth'.lower() in j_comminfo.lower():
                        self.J.jinfo.ey_azimuth = jcommval
                    
                    elif 'Hx azimuth'.lower() in j_comminfo.lower():
                        self.J.jinfo.hx_azimuth = jcommval
                    elif 'Coherence minimum'.lower() in j_comminfo.lower():
                        self.J.jinfo.coherence_min = jcommval
                    elif 'Coherence maximum'.lower() in j_comminfo.lower():
                        self.J.jinfo.coherence_max = jcommval
                    
                    elif 'Maxtype'.lower() in j_comminfo.lower():
                        self.J.jinfo.max_type = jcommval
                    elif 'Ntype'.lower() in j_comminfo.lower():
                        self.J.jinfo.ntype = jcommval
                    elif 'Weight'.lower() in j_comminfo.lower():
                        self.J.jinfo.weight = jcommval
                    
                    self.read_commentblock = True 
                #set longitude value 
                elif  re.match(r'^>', jitems) is not None :
                    jparam, jparamvalue =jitems.strip().split('=')
                    if 'azimuth'  in jparam.lower() :
                        self.jazim =jparamvalue 
                    elif 'longitude' in jparam.lower():
                        self.jlon=jparamvalue
                    elif 'latitude' in jparam.lower():
                        self.jlat =jparamvalue
                    elif 'elevation' in jparam.lower():
                        self.jelev =jparamvalue
                
                else :
                    #---> Get the station name , jmode and number of periods 
                    jmain = jitems.strip().split()
                    if jmain == None or jmain ==[]: pass 
                    elif len(jmain) < 5 :
                        try : 
                            self.jnperiods= float(jmain[0])
                        except :
                            if jmain[0] in SB._j :self.jmode=jmain[0]
                            else : self._jnames = jmain[0]
                            pass 
                        else :self.jnperiod = jmain[0] 
                    elif len(jmain) >= 5 and jmain !=['']:

                            jdata.append(np.array(jmain))

                
            JDAT= func.concat_array_from_list(list_of_array= jdata, concat_axis=0)

            if (re.match(r'^Z', self.jmode) is None) or (re.match(r'^T', self.jmode)  is None) :
                if JDAT.shape[1] == 9: 
                    jDATA = [np.reshape(jitem,( jitem.shape[0],)) for jitem in np.hsplit(JDAT,9) ]

                    self.jperiod, self.japp_rho, self.jpha, self.jrhomax , self.jrhomin , \
                        self.jphamax, self.jphamin , self.jwrho , self.jwpha = jDATA
                        
                else :
                    warnings.warn ('J-FORMAT expects to get "9" records'
                                   ' values like <{0}>'.format(*list(notion.j_recordR.keys())))
                    raise CSex.pyCSAMTError_J(
                        'For data type=R?? Only  9 Range values are not '
                        'acceptable. You provided <{0}>'.format(JDAT.shape[1]))
            elif (re.match(r'^Z', self.jmode) is True) or (re.match(r'^T', self.jmode)  is True):
                if JDAT.shape[1] == 5 : 
                    jDATA = [np.reshape(ii,( ii.shape[0],)) for ii in np.hsplit(JDAT,5) ]
                    self.jperiod , self.jreal , self.jimag , self.jerror, self.jweight = jDATA
                else :
                    warnings.warn('JFORMAT for GSC responses expects'
                                  ' to get "5" records values like<{0}>'.\
                                      format(*list(notion.j_recordZ.keys())))
                    raise CSex.pyCSAMTError_J(
                        'For data type=R?? Only  9 Range values'
                        ' are not acceptable. You provided <{0}>'.format(JDAT.shape[1]))
        
     
class J_infos (object):
    """
     J_infos class - set the information of J_file
    
    """

    def __init__(self, **kwargs): 
        
        self.progvers = 'WRITTEN BY pyCSAMT'
        self.site_infos = None 
        self.ex_length = None 
        self.ey_length =None 
        self.ex_azimuth =None 
        self.ey_azimuth =None 
        self.hx_azimuth =None
        self.coherence_min = None
        self.coherence_max =None 
        self.max_type =None 
        self.ntype =1
        self.weight ='Y'
        #---> use j config to quick write file .
        self.jconfig ={'Site info': self.site_infos, 
                       'Ex length': self.ex_length, 
                       'Ey length':self.ey_length, 
                       'Ex azimuth': self.ex_azimuth, 
                       'Ey azimuth': self.ey_azimuth, 
                       'Hx azimuth':self.hx_azimuth, 
                       'Coherence minimum':self.coherence_min, 
                       'Coherence maximum': self.coherence_max, 
                       'Maxtype':self.max_type, 
                       'Ntype':self.ntype , 
                       'Weight':self.weight,
                       }
        
        for keys in list(kwargs.keys()): 
            self.__setattr__(keys, kwargs[keys])

 
 
    
    
    
    
    
    
    
    
    
    
    
