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

.. _module-avg:: `csamtpy.ff.core.avg`  

    :synopsis: Module to deal with AVG file. each section of Avg file is a class 
    Avg file from Zonge Engeneering company is two types . module will read 
    the main file type as Type 1 . the othe type of AVG will rewrite on t
    the main type. the module can expor that file for other use .
    
Created on Wed Nov 18 16:45:35 2020

@author: KLaurent K. alias @Daniel03
"""

import os, re 
import warnings,shutil
import numpy as np
import pandas as pd
from csamtpy.etc import infos as infOS
from  datetime import datetime, timezone

from csamtpy.ff.core.cs import Site 

from csamtpy.ff.core.j import J_collection 
from csamtpy.ff.core.edi import Edi

from csamtpy.ff.processing import corr 
from csamtpy.ff.processing import callffunc as cfunc
from csamtpy.ff.processing import zcalculator as Zcc
from csamtpy.utils._csamtpylog import csamtpylog
from csamtpy.utils import func_utils as func
from csamtpy.utils import exceptions as CSex


_logger=csamtpylog.get_csamtpy_logger(__name__)



class Avg (object):
    """
    A super class of all content of AVG Zonge file 
    
        Deal with Zonge Engeering Avg file .
        
    Arguments 
    ----------
        **data_fn**: str 
                    path to AVG file 
        
    ================  ===========  =======================================
    Attributes         Type        Explanation
    ================  ===========  =======================================
    Header             class        container of Header informations 
    Data_section       class        container of all Data informations   
    Skip_flag          str          value to jugde the quality of data 
                                    can filled automatically . If no skp 
                                    value is provide , the program will
                                    judge quality of provided data 
                                    and set to him a score .
    ================  ===========  =======================================

    More attributes can be added by inputing a key word dictionary
    
    :Example:
        
        >>> Avg_object =Data(data_array=avg_data )
        >>> stn=DATA.Station.name 
        >>> freq=DATA.Frequency.value
    """
    encodage ='utf8'
    def __init__(self, data_fn=None , **kwargs):
        self._logging =csamtpylog.get_csamtpy_logger(self.__class__.__name__)
        self.data_fn =data_fn
        self.Header=Header()
        self.Data_section =Data()
        self.Skip_flag =Skip_flag()
        
        

        self._f1_labels , self._f2_labels=['skp','Station',
                                           ' Freq','Comp',
                                           ' Amps','Emag',
                                           'Ephz','Hmag',
                                           'Hphz','Resistivity',
                                           'Phase','%Emag',
                                           'sEphz','%Hmag',
                                           'sHphz','%Rho',
                                           'sPhz'],['skp',
                                                        'Freq','Tx.Amp',
                                                        'E.mag','E.phz',
                                                        'B.mag','B.phz',
                                                        'Z.mag','Z.phz',
                                                        'ARes.mag','SRes',
                                                        'E.wgt', 'B.wgt',
                                                        'E.%err','E.perr',
                                                        'B.%err','B.perr',
                                                        'Z.%err','Z.perr',
                                                         'ARes.%err']
        self._zmark='$'
                                                    
        self.EM_zonge2j = {'ExHx':'RXX', 
                            'ExHy':'RXY',
                                   'EyHx':'RYX', 
                                   'EyHy':'RYY'}
  
            
        if self.data_fn is not None :                       
            self._read_avg_file ()
                                      
    def _read_avg_file (self, data_fn=None ): 
        """
        Method to read avg file and set each objects .
        
        The main source protocole to transfert data into corresponding classes. 
        
        Parameters
        ------------
            * data_fn: str  
                    path to AVG filename.      
        """
        
        if data_fn is not None : 
            self.data_fn =data_fn 
        if self.data_fn is None  :
            raise CSex.pyCSAMTError_avg_file('No AVG file found.Please check '\
                                             'your pathfile :<{0}>'.format(os.getcwd()))
        # print(self.data_fn)       
        if self.data_fn is not None : 
            if os.path.isfile(self.data_fn) is False :
                raise CSex.pyCSAMTError_AVG('could not find %s', self.data_fn)
                
        self._logging.info('reading and checking the AVG file %s', self.data_fn)         
        #-------------------------------------------------------------
        
        if  (self.data_fn.endswith('avg') or self.data_fn.endswith('avg'.upper())) == False :#self.data_fn[-4:].lower() != '.avg'or self.data_fn.endswith('avg'.upper())==False or
            warnings.warn('<Get more infos about AVG> :{0}'.format(infOS.notion.AVG))
            raise CSex.pyCSAMTError_AVG('file provided is not an AVGfile.Please check your file!')
            
            
        elif self.data_fn.endswith('avg') or self.data_fn.endswith('.AVG'):
                with open (self.data_fn, 'r', encoding='utf8') as fid :
                    avg_data_file=fid.readlines()
            
        # ----now chech the avg file --
        avg_file_type_num = cfunc._validate_avg_file(avg_dataLines =avg_data_file)
        #------------------------------------
                
        if avg_file_type_num == 1: #-------------  main AVG file :F1 and keep data infos --------------

            # --> run the headlines thin to find the word 'skp'
            for idx , infos in enumerate (avg_data_file): 
                if 'skp' in infos : 
                    avg_data_file_infos = avg_data_file[:idx]
                    break

            self.Header.HardwareInfos.set_harware_infos([line.strip() for line  in \
                                                        avg_data_file_infos if line.startswith('\\')] )
            
            self.Header.Tx.set_transmitter_properties (Tx_data =[line.strip() for line  in \
                                                        avg_data_file_infos if line.startswith('$')])
            self.Header.Rx.set_receiver_properties (Rx_data =[line.strip() for line  in \
                                                        avg_data_file_infos if line.startswith('$')])
            
            data_to_array =avg_data_file [5:]
            # --- put data on list and array : 
            data_to_array = np.array([ datum.strip().split() for datum in data_to_array ])
            
            # self.skip.skp_flag==data_to_array[0,0] # fill the skip flag 
            self.Skip_flag.setandget_skip_flag(skip_flag=data_to_array[0,0])
            # delete the skp array and comp _array 
            
            data_to_array =np.delete(data_to_array, np.array([0,3]), axis =1 )

            self.Data_section._avg_set_components_from_data(data_array = data_to_array,
                                               data_type =avg_file_type_num) 
            # Skp_flag(skip_flag='2')
            
        elif  avg_file_type_num ==2 : #-------------Read ASTATIC file : F2 and keep data infos --------------
            
            # set zonge Hardwares infos:
            zonge_hardw_infos =[line.strip() for line  in avg_data_file if line.startswith('\\')]
            self.Header.HardwareInfos.set_harware_infos(zonge_hardw_infos=zonge_hardw_infos)
            
            self.Header.Tx.set_transmitter_properties(Tx_data=[line.strip()for \
                                                                      line in avg_data_file  if line.startswith('$Tx')])
            
            for ss, line_avg in enumerate (avg_data_file): 
                if line_avg.startswith('$Rx'):
                    temp_data =avg_data_file[ss:] # keep the remains data + Rx-properties
                    for_surv_configAnnot=avg_data_file[:ss] # truncate the data an keep the infos above thin $Rx
                    break 
                
            self.Header.SurveyAnnotation.set_survey_annotations_infos (survey_annotations_data =for_surv_configAnnot)
            self.Header.SurveyConfiguration.set_survey_configuration_infos (survey_config_data=for_surv_configAnnot)
            
            tem_d, rx_infos, comp_receiv_info=[],[],[]
            for jj, item in enumerate(temp_data): 
               #---> seek the data without the mark of the receiver of component list 
               if item.find('Z.mag')<0 and item.find('$Rx') < 0 :
                   ss__temp=item.strip().split(',')
                   for kk, value in enumerate ( ss__temp ) : 
                       try : 
                           ss__temp[kk]=float(value)
                       except :     #--> in the case we meet '*' mean no value .
                           if value.strip() == '*' :
                               ss__temp[kk]=np.nan
                           pass 
                       # if value=='*': ss__temp[kk]="0"
                   tem_d.append(ss__temp)
               if '$Rx' in item :
                   comp_receiv_info.append(item)
                
               if '$Rx.GdpStn' in item :
                   item = item.strip().split('=')[-1]
                   rx_infos.append(item)
               if item.find('ARes.mag')> 0 and item.find('SRes')> 0 : #- keep the head component infos
                    comp_sect=item.strip().split(',')
            
            #-set receivers infos 
            self.Header.Rx.set_receiver_properties(Rx_data =comp_receiv_info)
            
            #--> delete the none value on data 
            for sk , itemm in enumerate(tem_d): 
                if itemm==['']: tem_d.pop(sk)
            data=np.array(tem_d)
            data_to_array, other_comps = cfunc.straighten_cac2CSfile (data_array=data,
                                                   component_column_section = comp_sect)
            
            freq_counts,array_of_repeat_frequency=np.unique(data_to_array[:,1],
                                                    return_counts=True)
            
            station_list=rx_infos * freq_counts.size
            
            station_array_for_data=[float(ii) for ii in station_list]
            station_array_for_data.sort()
            station_array_for_data=np.array(station_array_for_data)
            if station_array_for_data.shape[0] !=data_to_array.shape[0]: 
                _logger.error('Trouble occured when reading Astatic file. It seems the frequency are not the '\
                              'same length along all stations.' )
                warnings.warn_explicit('it seems something wrong while reading Astatic file.', category=DeprecationWarning,
                                       filename= func.__name__,  lineno=1)
                raise CSex.pyCSAMTError_avg_file('File might be corrupted. number of frequency must'\
                                                 ' be the same along along all stations.  ')
                    
            
            data_to_array =np.delete(data_to_array, 0, axis =1 )
            # print(data_to_array.shape)
            # print( station_array_for_data.shape)
            
            data_to_array =np.concatenate((station_array_for_data.reshape(station_array_for_data.shape[0],1), 
                                              data_to_array), axis=1) 
            # data_to_array=np.insert(data_to_array, 1, station_array_for_data, axis=1 )

            self.Data_section._avg_set_components_from_data(data_array = data_to_array,
                                               data_type =avg_file_type_num, 
                                               add_astatic_data= other_comps)
            # -- we assume that if not problem occured thin now ,
            #  we set skip flag to '2' as good quality data 
            self.Skip_flag.setandget_skip_flag(skip_flag='2' ) 
        else : 
            raise CSex.pyCSAMTError_avg_file('The number provided doesnt not match Avgfile.')
            
            
            
            
    def avg_write_2_to_1 (self, data_fn=None , savepath =None ):
        """ 
        Method to rewrite avg Astatic file (F2) to main file F1 .
        
        :param data_fn: ASTATIC FILE 
        :type data_fn: str 
        
        :param savepath: path to save your rewritten file . 
        :type savepath: str 
        
        """
        if data_fn is not None : self.data_fn =data_fn 
        if self.data_fn is None : 
            raise CSex.pyCSAMTError_avg_file('No file detected. check your right path .')
        
        if self.data_fn is not None :
            if os.path.isfile (self.data_fn) is not True :
                raise CSex.pyCSAMTError_avg_file('No file detected.'\
                                                 'please check your right path. ')
            #-- check whether the file is astatic file : 
                
            checker = cfunc._validate_avg_file(self.data_fn)
            if checker == 1 : raise CSex.pyCSAMTError_avg_file(
                    'Input Avgfile < {0}> is not Zonge ASTATIC file. '\
                        'Please put the right file !'.format(os.path.basename(self.data_fn)))
            
        write_lines =['\\']
        
        self._read_avg_file(data_fn=self.data_fn)
        HW_dict =self.Header.HardwareInfos.sconfig_dict
        for keys, values in HW_dict.items():
            if 'astatic' == (str(keys)).lower() or keys=='' :continue
            elif type(values) is list:values =''.join([ss for ss in values])
            write_lines.append(''.join([keys, ':',values]+[',']))
                
        
        #add copyleft.
        write_lines.append(self.Header.HardwareInfos._writer_copyleft+
                           ': {0}-{1}'.format(datetime.now(), timezone.utc) +'\n')       
                
         # start writing configure part.
        write_lines.append("".join([self._zmark,'{0:<7}'.format(' ASPACE'),
                                    '=','{0:>7}'.format( self.Header.Rx.rxLength),
                                    self.Header.Rx.unit])+'\n')
        write_lines.append("".join([self._zmark,'{0:<7}'.format(' XMTR'),
                                    '=', '{0:>7}'.format(self.Header.Tx.txGdpStn)])+'\n')
        # write_lines.append(''.join(["{0:^9}".format(idx) for idx in self._f1_labels])+'\n')
        
        #------write the interligne labels and markers ---------

        def draw_fc(format_len =['{:^3}','{:^7}','{:^7}',
                                 '{:^7}','{:^6}','{:^12}',
                                 '{:^8}','{:^12}','{:^9}',
                                 '{:^12}','{:^9}','{:^7}',
                                 '{:^7}','{:^7}','{:^7}',
                                 '{:^7}','{:^7}'] , value_to_format=self._f1_labels):
            """ 
            Drawn formatage chain . 
            
            :param format_len: list of formatage indicator , flexible .
            :type format_len: list 
            
            :param value-to-format: the list of item we want to format.
            :type value_to_format: list 
                
            .. note:: user can change the disposal of avg_labels.
    
            """
            for ss, item in enumerate(format_len):value_to_format[ss]=item.format(value_to_format[ss])
            
            return ''.join([new_carf for new_carf in value_to_format])
        
        write_lines.append(draw_fc()+'\n')
        write_lines.append('\\')
        write_lines.append('++'.join([int(ss) * '-' for ss in ['1','5','5','5','4','10',
                                                               '6','10','7','10','7','5',
                                                               '5','5','5','5','5'] ])+'+\n')

        #------------start write block data -------------------
        stnNames =self.Data_section.Station.names

        
        def nan_to_zstar(nan_array, swaper='*' ): # not usefull , we replace nan at the bottom .
            """
            Swapp np.nan number on the data .
            """
            nan_array= np.nan_to_num(nan_array, nan=-999)
            return np.where(nan_array==-999, swaper, nan_array)
        
            #--bring back all data ----
        stn=np.concatenate(([self.Data_section.Station.loc[name] for  name in  stnNames]), axis=0)
        freq=np.concatenate(([self.Data_section.Frequency.loc[name] for  name in  stnNames]), axis=0)
        amps=np.concatenate(([self.Data_section.Amps.loc[name] for name in stnNames]), axis=0)
        emag=np.concatenate(([self.Data_section.Emag.loc[name] for  name in  stnNames]), axis=0)
        ephz=np.concatenate(([self.Data_section.Ephz.loc[name] for  name in  stnNames]), axis=0)
        bmag=np.concatenate(([self.Data_section.Hmag.loc[name] for name in  stnNames]), axis=0)
        bphz=np.concatenate(([self.Data_section.Hphz.loc[name] for  name in  stnNames]), axis=0)
        rho=np.concatenate(([self.Data_section.Resistivity.loc[name] for  name in  stnNames]), axis=0)
        phase=np.concatenate(([self.Data_section.Phase.loc[name] for name in  stnNames]), axis=0)
        pcemag=np.concatenate(([self.Data_section.pcEmag.loc[name] for name in stnNames]), axis=0)
        sephz=np.concatenate(([self.Data_section.sEphz.loc[name] for name in  stnNames]), axis=0)
        pchmag=np.concatenate(([self.Data_section.pcHmag.loc[name] for  name in  stnNames]), axis=0)
        shphz=np.concatenate(([self.Data_section.sHphz.loc[name] for name in stnNames]), axis=0)
        pcrho=np.concatenate(([self.Data_section.pcRho.loc[name] for  name in stnNames]), axis=0)
        sphase=np.concatenate(([self.Data_section.sPhz.loc[name] for name in  stnNames]), axis=0)

        # build skp_flag  and fill measure component type ...
        skpf =np.full((stn.shape[0],), self.Skip_flag.skip_flag)
        
        comps =np.full((stn.shape[0],),self.Header.Rx.rxComps )
        
        #----> write data line by line ------
        for ii in range (skpf.shape[0]):
            write_lines.append(''.join(['{0:^3}'.format(skpf[ii]),
                                                        '{0:^7}'.format((stn[ii])),
                                                        '{0:>7}'.format(freq[ii]), 
                                                        '{0:^7}'.format(comps[ii]),
                                                        '{0:<6}'.format(amps[ii]),
                                                        '{:12.4e}'.format(emag[ii]),
                                                        '{:>8}'.format(ephz[ii]),
                                                        '{:12.4e}'.format(bmag[ii]),
                                                        '{:>9}'.format(bphz[ii]),
                                                        '{:12.4e}'.format(rho[ii]),
                                                        '{:>9}'.format(phase[ii]),
                                                        '{:>7}'.format(pcemag[ii]),
                                                        '{:>7}'.format(sephz[ii]),
                                                        '{:>7}'.format(pchmag[ii]),
                                                        '{:>7}'.format(shphz[ii]),
                                                        '{:>7}'.format(pcrho[ii]),
                                                        '{:>7}'.format(sphase[ii])]))
            write_lines.append('\n')
            
        #---- replace np.nan number by zonge star who mean no value ----
        for jj, lines in enumerate(write_lines ): 
            if 'nan' in lines : lines =lines.replace('nan','{:>3}'.format('*'))
            write_lines[jj]=lines
        #---------------------end nan replacing ----------------        
        with open(''.join([os.path.basename(self.data_fn)[:-4],'_2_to_1.avg']), 
                  'w', encoding='utf8') as fw:
            fw.writelines (write_lines)
        #---move to your savepath---
        if savepath is not None : 
            shutil.move(''.join([os.path.basename(self.data_fn)[:-4],'_2_to_1.avg']),
                        savepath)
        print('--->  Your <{0}> is rewritten to <{1}> successfully. <-----'.\
              format(os.path.basename(self.data_fn), 
                   ''.join([os.path.basename(self.data_fn)[:-4],'_2_to_1.avg'])))
                
    def avg_to_jfile (self, avg_data_fn=None , station_fn =None, j_extension='dat', 
                      utm_zone='49N',**kws):
        """
        Method to write avg file to Jfile, convert both files , Astatic or plainty
         avg file to A.G. Jones format.
        
        Parameters 
        -----------
            * avg_data_fn:str 
                        pathLike , path to your avg file 
            * station_fn:  str
                        pathLike, path to your profile/station file .
            * j_extension: str
                        Extension type you want to export file .
                        *Default* is ".dat"
            * utm_zone: str
                        add  if station_profile are not referenced yet.
                        later , it would be removed .
            * savepath: str
                    pathLike , path to save your outpufile 
            * write_info: bool
                        write the informations of your input file ,  
                        export informations into Jfile,*Default* is False.
            * survey_name:   bool, 
                        survey_area  
        """
        savepath =kws.pop('savepath', None)
        write_info =kws.pop('writeInfos', False)
        survey_name =kws.pop('survey_name', None)
        
        
        if j_extension is not None and '.' not in j_extension: j_extension ='.'+j_extension
        
        
        self._logging.info('write Jformat file')
        
        if avg_data_fn is not None : self.data_fn = avg_data_fn 
        if self.data_fn is None : raise CSex.pyCSAMTError_avg_file("Could not find any path to read .'\
                                                                   Please provide your right AVG file.")
        if station_fn is  None : raise CSex.pyCSAMTError_station('Need  absolutely station file. please provide your ".STN" file.')
        
        #---> read the avg file 
        self._read_avg_file(data_fn=self.data_fn)
        
        #---> read profile_obj from site 
        site_obj = Site(data_fn=station_fn, utm_zone=utm_zone)#.set_site_info(data_fn=, utm_zone = utm_zone)

        #---> get the list of stations 
        station_list = sorted(self.Data_section.Station.names)
        
        #---> Build jnames with survey name 
        if survey_name is None : 
            survey_name = self.Header.SurveyAnnotation.project_area 
            
        j_obj=J_collection()
        jstnnames = j_obj.J.jname(number_of_sites = len(station_list),
                                  survey_name=survey_name)
        #---> get component and check it with j 
        self.Header.Rx.rxComps
        for key, compvalue in self.EM_zonge2j.items(): 
            if self.Header.Rx.rxComps  == key : 
                jcomp = j_obj.J.jMode(polarization_type = compvalue) 
        #---> number of periods 
        jnperiod = self.Data_section.Frequency.numfreq
        jperiod = 1/self.Data_section.Frequency.value # frequency value 
        # call data section 
        app_rho_obj = self.Data_section.Resistivity.loc 
        phase_obj = self.Data_section.Phase.loc
        emag_obj ,hmag_obj =self.Data_section.Emag.loc, self.Data_section.Hmag.loc  
        c_var_emag_obj , c_var_hmag_obj = self.Data_section.pcEmag.loc , self.Data_section.pcHmag.loc
        c_var_app_rho_obj , std_phase_obj =self.Data_section.pcRho.loc, self.Data_section.sPhz.loc
        
        #---> compute rhomax and rhomin 
        japp_rho_max, japp_rho_min , jphase_max, jphase_min,jwapp_rho =[{} for ii in range(5)] # create 4 empy  dicts
        for key, value in self.Data_section.Resistivity.loc.items(): 
            japp_rho_max[key]= app_rho_obj[key] + 1* Zcc.compute_sigmas_e_h_and_sigma_rho(pc_emag= c_var_emag_obj[key],
                                                                  pc_hmag= c_var_hmag_obj[key] ,
                                                                  pc_app_rho=c_var_app_rho_obj[key],
                                                                  app_rho= app_rho_obj[key], 
                                                                  emag= emag_obj[key],
                                                                  hmag=hmag_obj[key])[-1]
            japp_rho_min[key]= app_rho_obj[key] - 1 * Zcc.compute_sigmas_e_h_and_sigma_rho(pc_emag= c_var_emag_obj[key],
                                                                  pc_hmag= c_var_hmag_obj[key] ,
                                                                  pc_app_rho=c_var_app_rho_obj[key],
                                                                  app_rho= app_rho_obj[key], 
                                                                  emag= emag_obj[key],
                                                                  hmag=hmag_obj[key])[-1]
        for key, value in self.Data_section.Phase.loc.items(): 
            jphase_max [key]= 180 * (phase_obj[key] + 1 * std_phase_obj[key]/100)/np.pi
            jphase_min [key] = 180 * (phase_obj[key] - 1 * std_phase_obj[key]/100)/np.pi
        
        #convert phase millirad value in degree : 
        phase_obj ={key: 180 * phase * 1e-3/np.pi for key, phase in phase_obj.items()}
        # compute jwro and jwpha 
        #we assume to take weight to 1.
        jwapp_rho ={ key:value for key ,
                    value in zip(station_list, 
                     [np.ones((self.Data_section.Frequency.numfreq,), 
                               dtype='float64') for ii in range(len(station_list))])}
        
        import copy  # use deep copy to copy the same weight and re  
        jwpha = copy.deepcopy(jwapp_rho)
        
        #-----> starting writing Jfiles <----- 
                #--- > start writing 
        write_jlines =[]
        if survey_name is None : code_name = 'ks{0:02}-{1}'
        elif survey_name is not None : code_name = survey_name[:-3]+'{0:02}-{1}'
        codespace= 5*' '
        
        if savepath  is None : # create a folder in your current work directory
           try :
               savepath = os.path.join(os.getcwd(), '_outputJ_')
               if not os.path.isdir(savepath):
                   os.mkdir(savepath)
           except : 
               warnings.warn("It seems the path already exists !")
        
        for ii , stn in enumerate(station_list) : 
            #-- > write Head j 
            write_jlines.append(''.join([j_obj.J._comment_mark,'{0:} :'.format(j_obj.J.jinfo.progvers) , 
                                        code_name.format(ii, station_list[ii]), codespace, 
                                        datetime.now().strftime(j_obj.J.jfmtime), codespace, 'RECS'])+'\n')
            #--> write info : 
            if write_info : 
                try :
                    counter_break=0
                    lineinfo=[]
                    #--> write zongeHardware infos 
                    for zongeinfo in [self.Header.HardwareInfos.sconfig_dict,
                                      self.Header.SurveyAnnotation.sconfig_dict,
                                      self.Header.SurveyConfiguration.sconfig_dict,
                                      self.Header.Tx.sconfig_dict,
                                      self.Header.Rx.sconfig_dict]:
                        lineinfo.append(j_obj.J._comment_mark)
                        for keyinfo, valueinfo in zongeinfo.items():
                            if isinstance(valueinfo, list): valueinfo = ''.join([ss for ss in valueinfo if ss !=None])
                            if 'None' in valueinfo: valueinfo=valueinfo.replace('None','*')
                            if valueinfo is not None : 
                                lineinfo.append("{0}{1:^2}{2}".format(keyinfo,':',valueinfo))
                                lineinfo.append('{0:3^}'.format(','))
                                counter_break +=1 
                                if counter_break==7 :#--> go to line  to let visible info on jfile more visible.
                                    lineinfo.append('\n#')
                                    counter_break=0
                        counter_break=0           
                        write_jlines.append(''.join(lineinfo)+'\n')  
                        lineinfo=[]

                except :pass 
     
            #--- > write jlabels 
            write_jlines.append(''.join(['{0:<13}'.format(j_obj.J._jlabels[0]),
                                         '=', "{0:>13.1f}".format(site_obj.azim[stn]), '\n'] ))
            write_jlines.append(''.join(['{0:<13}'.format(j_obj.J._jlabels[1]),
                                         '=', "{0:>13.5f}".format(site_obj.lat[stn]), '\n'] ))
            write_jlines.append(''.join(['{0:<13}'.format(j_obj.J._jlabels[2]),
                                         '=', "{0:>13.5f}".format(site_obj.lon[stn]), '\n'] ))
            write_jlines.append(''.join(['{0:<13}'.format(j_obj.J._jlabels[3]),
                                         '=', "{0:>13.1f}".format(site_obj.elev[stn]), '\n'] ))
            #---> write jnames , components  and nperiods 
            write_jlines.append(''.join(['{0}'.format(jstnnames[ii]), 
                                         '{0:>7.1f}'.format(site_obj.azim[stn]),'\n']))
            write_jlines.append('{0}'.format(jcomp)+'\n')
            write_jlines.append(' {0}'.format(jnperiod)+'\n')
            # now  write value :
                
            for jj in range(jnperiod):

                write_jlines.append (''.join( ['{0:<12.3e}'.format(jperiod[jj]), 
                                               '{0:^12.3e}'.format(app_rho_obj [stn][jj]), 
                                               '{0:^12.2e}'.format(phase_obj[stn][jj]), 
                                               '{:^12.3e}'.format(japp_rho_max[stn][jj]), 
                                               '{:^12.3e}'.format(japp_rho_min[stn][jj]), 
                                               '{:^12.2e}'.format(jphase_max[stn][jj]), 
                                               '{:^12.2e}'.format(jphase_min[stn][jj]),
                                               '{:^7.2f}'.format(jwapp_rho[stn][jj]), 
                                               '{:^7.2f}'.format(jwpha[stn][jj]),
                                                                 ]))

                write_jlines.append('\n')
                
            with open('{0}{1}'.format(station_list[ii],j_extension), 'w', encoding='utf8') as fj:
                fj.writelines(write_jlines)  
            
                
            if savepath is not None : 
                shutil.move ('{0}{1}'.format(station_list[ii],j_extension), savepath)
                                                           
            write_jlines=[]  
            
            
        print('-'*77) 
        # if savepath is None :savepath =os.getcwd()
        print('---> {0} J-files have been rewritten.\n---> see path:<{1}> '.\
              format(len(station_list), savepath))
        print('-'*77)
        
        
    def avg_to_edifile (self, data_fn =None , profile_fn =None , 
                        savepath =None , apply_filter =None, 
                        reference_frequency =None, number_of_points=7.,
                        dipole_length=50., **kwargs):
        """
        Method to write avg file to SEG-EDIfile.Convert both files.Astatic or plainty avg file .
        if ASTATIC file is provided , will add the filter and filter values .
        if avg file is not astatic file , user an apply
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
                        can be ; [TMA | AMA| FLMA]
                        TMA - Trimming Moving Average
                        AMA - Adaptative Moving avarage , 
                        FLMA - Fixed dipoleLength moving 
                        average (number of point=7)
                        Dipolelength willbe computed automatically 
        
        Holdings additionals parameters :
            
        ====================  ==========  =====================================
        Optionalparams          Type        Description 
        ====================  ==========  =====================================
        reference_frequency     float       frequency at clean data.Default is 
                                            computed automatically
        number_of_points=7.    int          number of point to for weighted  
                                            window,. *Default* is 7.
        dipole_length          float        length of dipole in meter. 
                                             For CSAMT, Default is 50.m
        ====================  ==========  =====================================
                                     
        :Example:   
            
            >>> from csamtpy.ff.core import avg 
            >>> avg_obj= avg.Avg()
            >>> avg_obj.avg_to_edifile(data_fn= os.path.join(path_to_avgfile, avgfile) , 
            ...           profile_fn = os.path.join(path_to_avgfile, station_profile_file), 
            ...           savepath =save_edipath, 
            ...           apply_filter=None ) 
        """
        
        utm_zone = kwargs.pop('utm_zone', None)
        
        if data_fn is not None : self.data_fn = data_fn 
        if self.data_fn is None : 
            raise CSex.pyCSAMTError_avg_file("Could not find any path to read ."\
                                            "Please provide your right AVG file.")
        
        
        # export to savepath 
        if savepath is None : # create a folder in your current work directory
            try :
                savepath = os.path.join(os.getcwd(), '_outputAVG2EDI_')
                if not os.path.isdir(savepath):
                    os.mkdir(savepath)#  mode =0o666)
            except : 
                warnings.warn("It seems the path already exists !")
        
        
        if profile_fn is None : 
            warnings.warn ('No station profile file will detected.'\
                           ' Be aware sure , we will set location longitude ,'\
                               ' latitude and elevation to <0.>')
        
        # read avg file 
        self._read_avg_file(data_fn=self.data_fn)
        
        if profile_fn is not None : 
             site_obj = Site(data_fn=profile_fn, utm_zone=utm_zone)
                                                                   
        #---> get the list of stations 
        head_dataid_list = sorted(self.Data_section.Station.names)
        

        #-- compute dipole lenght  between station close 
        dipole_length= self.Data_section.Station.loc[head_dataid_list [1]][0]-\
            self.Data_section.Station.loc[head_dataid_list[0]][0]
        
        #--> nfrequency = 
        nfreq = self.Data_section.Frequency.value.size
        
        
        #------------------------------------------------------------------------------------
        #---> compute error propagation , phase and resistivity 
        app_rho_obj = self.Data_section.Resistivity.loc 

        phase_obj = self.Data_section.Phase.loc
        emag_obj ,hmag_obj =self.Data_section.Emag.loc, self.Data_section.Hmag.loc  
        c_var_emag_obj , c_var_hmag_obj = self.Data_section.pcEmag.loc , self.Data_section.pcHmag.loc
        c_var_app_rho_obj , std_phase_obj =self.Data_section.pcRho.loc, self.Data_section.sPhz.loc
        #self.Header.HardwareInfos.astatic_version
        AS_flag =0          # maker to check whether is plainty avg file or Astatic file 
        if self.Header.HardwareInfos.numfilterfreq is not None :
            # self.Header.HardwareInfos.numfilterfreq is not None :
            print('---> Reading Zonge `ASTATIC` file !')
            AS_rho =self.Data_section.Resistivity.loc_Sres
            AS_flag =1
        
        error_propag_rho, error_propag_phs =[{} for ii in range(2)] # create 2empty  dicts
        for key, value in self.Data_section.Resistivity.loc.items(): 
            error_propag_rho[key]= app_rho_obj[key] + 1* Zcc.compute_sigmas_e_h_and_sigma_rho(
                                                                pc_emag= c_var_emag_obj[key],
                                                                  pc_hmag= c_var_hmag_obj[key] ,
                                                                  pc_app_rho=c_var_app_rho_obj[key],
                                                                  app_rho= app_rho_obj[key], 
                                                                  emag= emag_obj[key],
                                                                  hmag=hmag_obj[key])[-1]
        for key, value in self.Data_section.Phase.loc.items(): 
            error_propag_phs [key]= 180 *( (std_phase_obj[key]/100)/1000)/np.pi #  std = sPhz / 100 #converted in radians
  
        
        #convert phase millirad value in degree : 
        phase_obj ={key: 180 * phase * 1e-3/np.pi for key, phase in phase_obj.items()}
    
        #-----------------------------------------------------------------------------------------
        if apply_filter is not None :
            try :
                apply_filter= apply_filter.lower()
            except:
                print("---> ErrorType: Filters' names must a str of"
                              " `tma` or `ama`. not {0}".format(type(apply_filter)))
                
                warnings.warn("TypeError ,Filters' names must a str of"
                              " `tma` or `ama`. not {0}".format(type(apply_filter)))
                apply_filter=None
            else :
                print('{0:-^77}'.format('Filter Infos'))
                print('** {0:<27} {1} {2}'.format('Filter applied', '=', apply_filter.upper()))
                
            if apply_filter =='tma':
                self._logging.info (
                    'Computing static shift Trimming moving average (TMA).')
                        
                corr_obj= corr.shifting()
                res_TMA= corr_obj.TMA(data_fn =data_fn,
                                      number_of_TMA_points =number_of_points , 
                                      reference_freq=reference_frequency, 
                                      profile_fn = profile_fn) # dictionnary of value computed 
                reference_frequency =corr_obj.referencefreq
                
                print('** {0:<27} {1} {2}'.format('TMA filter points', 
                                                  '=', number_of_points))
                print('** {0:<27} {1} {2} Hz.'.format('Reference frequency',
                                                      '=', reference_frequency))
                
            elif apply_filter =='flma':
                self._logging.info (
                    'Computing static shift fixed-length moving average (FLMA).')
                        
                corr_obj = corr.shifting()
                res_FLMA = corr_obj.FLMA(data_fn =data_fn,
                                        profile_fn =profile_fn ,
                                        number_of_points=number_of_points , 
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
                       msg =''.join([
                            '!--> Sorry , Actually the available filter are only',
                                ' `tma` and `flma`.adaptive moving average `AMA`'
                                ' will be available soon!you may use '
                                ' filters `tma`or `flma`.'])
                       
                       print(' >{0}'.format(msg))
                       self._logging.debug(msg)
                       warnings.warn(msg)
                       apply_filter = None   
                       
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

            avgfilter ="{0}.{1}.Filter ={2} for {3} {4}.".format(
                   'pycsamt','v1.0.01', apply_filter.upper(),  
                   int(number_of_points) , msf)
            
        if profile_fn  is None : 
            #---> set to 0. lon , lat and elev if profile _fn is not provided
            #--> create dictionnary of zero value of lat and lon 
            import copy 
            #head_edi_lon = np.full((len(head_dataid_list),), 0., dtype=np.float)
            head_edi_lon = {stn: value for stn , value in zip (
                head_dataid_list, np.zeros_like(len(head_dataid_list), dtype=np.float))}
            head_edi_lat = copy.deepcopy(head_edi_lon)
            head_edi_elev = copy.deepcopy(head_edi_lon)
            
        elif profile_fn is not None : 
             head_edi_lon = site_obj.lon
             head_edi_lat = site_obj.lat 
             head_edi_elev = site_obj.elev 
        
        #------------------------START SETTING EDI ATTRIBUTE -----------------------------
        # import module Hmeasurement and Emeasurement
        from csamtpy.ff.core.edi import Hmeasurement, Emeasurement
        # from csamtpy.utils import gis_tools as gis 
        
        for ii, stn in enumerate(head_dataid_list): #   for stn in head_dataid_list : #loop for all station or dataid 
            #create Ediobj for eac datalist 
            edi_obj =Edi()
            
            #====> set edi_obj Header attributes
            edi_obj.Head.dataid= head_dataid_list[ii] # ii mean for all the list 
            edi_obj.Head.acqby = 'Zonge Engineering'
            edi_obj.Head.acqdate = self.Header.HardwareInfos.dated 
            edi_obj.Head.fileby = "AMTAVG_v.{0}".format(self.Header.HardwareInfos.version) 
            if self.Header.SurveyAnnotation.project_name is None :
                 self.Header.SurveyAnnotation.__setattr__('project_name ',os.path.basename(self.data_fn)) 
    
            edi_obj.Head.loc = os.path.basename(self.data_fn)[:-4]
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
            
               
            edi_obj.Info.Source.__setattr__('project', self.Header.SurveyAnnotation.project_name)
            edi_obj.Info.Source.__setattr__('survey',  head_dataid_list[ii])
            edi_obj.Info.Source.__setattr__('sitename', head_dataid_list[ii])
            edi_obj.Info.Processing.__setattr__('processedby', 'pyCSAMT' )
            edi_obj.Info.Processing.ProcessingSoftware.__setattr__('name', edi_obj.Head.fileby )
            
    
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
            edi_obj.DefineMeasurement.__setattr__('meas_ex', 
                                                  Emeasurement(**{'id':'{0:04}{1}'.format(ii * 10 + 1 , 
                                                                                          codeID_dec[1:] ), 
                                                                          'chtype':'EX', 
                                                                          'x': -(dipole_length/2), 
                                                                          'y': 0.,
                                                                          'x2':dipole_length/2 , 
                                                                          'y2':0, 
                                                                          'acqchan': '0.00', 
                                                                          'filter':avgfilter, 
                                                                          }))

            edi_obj.DefineMeasurement.__setattr__('meas_ey',
                                                  Emeasurement(**{'id':'{0:04}{1}'.format(ii * 10 + 2 , 
                                                                                          codeID_dec[1:]), 
                                                                          'chtype':'EY', 
                                                                          'x':0., 
                                                                          'y': -(dipole_length/2),
                                                                          'x2':0., 
                                                                          'y2':dipole_length/2 , 'acqchan': '0.00', 
                                                                          'filter':avgfilter}))

            edi_obj.DefineMeasurement.__setattr__('meas_hx', 
                                                  Hmeasurement(**{'id':'{0:04}{1}'.format(ii * 10 + 3 , 
                                                                                          codeID_dec[1:] ), 
                                                                              'chtype':'HX', 
                                                                              'x':0., 
                                                                              'y': 0.,
                                                                              'x2':0., 
                                                                              'y2':0. , 'acqchan': '0.00', 
                                                                          'filter':avgfilter}))

            edi_obj.DefineMeasurement.__setattr__('meas_hy',
                                                  Hmeasurement(**{'id':'{0:04}{1}'.format(ii * 10 + 4 ,
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
            edi_obj.MTEMAP.__setattr__('hx', '{0:04}{1}'.format(ii * 10 + 3 , codeID_dec[1:] ))
            edi_obj.MTEMAP.__setattr__('hy', '{0:04}{1}'.format(ii * 10 + 4 , codeID_dec[1:] ))
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
            
    
class Header (object):
    
    """
    Read the info Header of AVG file and rewrite Avgfile (main type, Type1) .
    
    Arguments  
    ---------
        **data_fn** : str  
                path to avgfile
        
    ======================  =========  =====================================
    Attributes              Type        Explanation
    ======================  =========  =====================================
    HardwareInfos           class        zonge Hardware  
    SurveyAnnotation        class        zonge project  informations 
    SurveyConfiguration     class        zonge project/locations settings 
    Tx                      class        special Transmitter  properties 
    Rx                      class        special receiver properties
    ======================  =========  =====================================

    More attributes can be added by inputing a key word dictionary
    
    :Example:
        
        >>> from csamtpy.ff.core.avg import Avg
        >>> path=os.path.join(os.environ["pyCSAMT"], 
        ...          "csamtpy", "data", "LC101.avg")
        >>> avg_obj= Avg(data_fn=path) 
        >>> surv_area= avg_obj.Header.SurveyAnnotation.project_area 
        >>>> print(survey.area)
        
    """
    def __init__(self, header_infos = None , **kwargs):

        self.HardwareInfos =ZongeHardware()
        self.SurveyAnnotation=SurveyAnnotation()
        self.SurveyConfiguration=SurveyConfiguration()
        self.Tx=TransmitterProperties()
        self.Rx=ReceiverProperties()

        
        for keys  in list(kwargs.keys()):       
            self.__setattr__(keys, kwargs[keys])
        
        self.head_block_infos =['Hardware infos', 
                                'Survey annotation', 
                                'Survey configuration', 
                                'Transmitter properties', 
                                'Receiver properties', 
                                'end']

    def write_header_log(self, data_fn=None , savepath=None ): 
        """ 
        Method to write your head log. In fact just to see whether your Zonge Infos 
        was set correctly.
        
        Parameters 
        ----------
            * data_fn : str 
                        path to Avgfile
            * savepath : str 
                        path to your destination file . if None , 
                        path is your current work directory .
            
        :Example:
            
            >>> from csamtpy.ff.core.avg.Avg import Header
            >>> Header().write_header_log(data_fn=path,
            ...               )
            >>> he=Header()
            ... he.write_header_log(data_fn=path,
            ...                   savepath=r'C:/Users\Administrator\Desktop')
            ... avg_obj.Header.write_header_log(data_fn=path,
            ...               savepath=r'C:/Users/Administrator/Desktop') 
        """
         
        if data_fn is not None : self.data_fn =data_fn 
        if savepath is None : savepath=os.getcwd()

        if os.path.isfile(self.data_fn ) is True :
            avg_obj=Avg(data_fn=self.data_fn)
            avg_obj._read_avg_file(data_fn=self.data_fn)
            
        rmlines, spac_b, retc =[], " "*4,'\n'
        for ss , items in enumerate(self.head_block_infos) : 
            rmlines.append('>'+ items + '\n')
            
            if ss==0 :
                for keys, values in avg_obj.Header.HardwareInfos.sconfig_dict.items():
                    rmlines.append(spac_b)
                    rmlines.append(''.join(['{0}:'.format(str(keys)),
                                            '{0:}'.format(str(values))]))
                    rmlines.append(retc)
            elif ss==1 : 
                for keys , values in avg_obj.Header.SurveyAnnotation.sconfig_dict.items():
                    rmlines.append(spac_b)
                    rmlines.append(''.join(['{0}:'.format(str(keys)),
                                            '{0:}'.format(str(values))]))
                    rmlines.append(retc)
            elif ss==2 :
                for keys , values in avg_obj.Header.SurveyConfiguration.sconfig_dict.items():
                    rmlines.append(spac_b)
                    rmlines.append(''.join(['{0}:'.format(str(keys)),
                                            '{0:}'.format(str(values))]))
                    rmlines.append(retc)
            elif ss ==3 :
                for keys , values in avg_obj.Header.Tx.sconfig_dict.items():
                    rmlines.append(spac_b)
                    rmlines.append(''.join(['{0}:'.format(str(keys)),
                                            '{0:}'.format(str(values))]))
                    rmlines.append(retc)
            elif ss == 4 : 
                for keys , values in avg_obj.Header.Rx.sconfig_dict.items():
                    rmlines.append(spac_b)
                    rmlines.append(''.join(['{0}:'.format(str(keys)),
                                            '{0:}'.format(str(values))]))
                    rmlines.append(retc)
                    
        with open (''.join([os.path.basename(self.data_fn),'_pyCS',".log"]),
                   'w', encoding='utf8') as fd:
            fd.writelines(rmlines)
            
        if savepath is not None :
            shutil.move (''.join([os.path.basename(self.data_fn),'_pyCS',".log"]),savepath)
        
        print('---> Your Hlog file <{0}> has been written successfully !<---'.\
              format(''.join([os.path.basename(self.data_fn),'_pyCS',".log"])))
            

class SurveyAnnotation (object) : 
    """
    Class for survey annotations.
         
    Arguments  
    -----------
        **survey_annotations_data**: list , str  
                                    path to your annotation file or  
                                    Filled automatically
        
    ======================  =========  =====================================
    Attributes              Type        Explanation
    ======================  =========  =====================================
    project_name            str         project name 
    project_area            str         project area 
    custumer_name           str         company name 
    contractor_name         str         contractor name 
    projet_label            str         label or number 
    acqdate                 str         data acquisition date 
    sconfig_dict            dict        survey annotation
                                        dictionnary       
    ======================  =========  =====================================
    
    More attributes can be added by inputing a key word dictionary
    
    :Example:
         
        >>> from csamtpy.ff.core.avg import Avg 
        >>> avg_obj= Avg(data_fn=path) 
        >>> surv_area= avg_obj.Header.SurveyAnnotation.project_area 
        ... survey_area
        ... avg_obj.Header.SurveyAnnotation.sconfig_dict['Job.By']
    """
    def __init__(self, survey_annotations_data =None  ,  **kwargs):

        self.zon_serv_annotations = survey_annotations_data
        self.project_name ='CSAMT_survey'
        self.project_area =None
        self.custumer_name ='Zonge Engineering"' 
        self.contractor_name ='Zonge'
        self.project_label = 'pyCSAMT' 
        self.acqdate ='{0}-{1}'.format(datetime.now(), timezone.utc)
        self.sconfig_dict={'Job.Name':self.project_name , 
                      'Job.Area':self.project_area , 
                      'Job.For':self.custumer_name, 
                      'Job.By':self.contractor_name, 
                      'Job.Number':self.project_label , 
                      'Job.Date': self.acqdate, 
                      }
        
        for keys in list(kwargs.keys()):
            if keys not in SurveyAnnotation.__dict__.keys():
                self.__setattr__(keys, kwargs[keys])
                self.sconfig_dict.__setitem__(keys, kwargs[keys])
                
    def set_survey_annotations_infos (self, survey_annotations_data =None): 
        """
        Method to set _survey annotations informations .
        
        Parameters 
        ----------
            * survey_annotation : list 
                                container of survey annotations infos. 
        """
        _logger.info('Reading and setting zonge survey Annotations infos !')
        
        if survey_annotations_data is not None : 
            self.zon_serv_annotations = survey_annotations_data
            
        if self.zon_serv_annotations is None : 
            CSex.pyCSAMTError_inputarguments('No survey_annotations informations found!')
            
            
        if type(self.zon_serv_annotations ) is not list: 
            try :
                if os.path.isfile(self.zon_serv_annotations ) is True :
                    with open (self.zon_serv_annotations , 'r', encoding ='utf8') as fid :
                        self.zon_serv_annotations =fid.readlines()
            except : raise CSex.pyCSAMTError_Header("Please check your annotation data ! "\
                                                    "Must be be a list or file ") 
            
        for zonsurva in self.zon_serv_annotations :

            zonsurva =zonsurva.strip().split('=')
            if 'job.name' in zonsurva[0].lower(): 
                if zonsurva[1]  !='': self.project_name=zonsurva[1]
            if 'job.area' in zonsurva[0].lower() :self.project_area =zonsurva[1]
            if 'job.for' in zonsurva[0].lower() :
                if zonsurva[1]  !='' :self.custumer_name = zonsurva[1]
            if 'job.by'in zonsurva[0].lower() :
                if zonsurva[1]  !='' :self.contractor_name =  zonsurva[1]
            if 'job.number'in zonsurva[0].lower() :
                if zonsurva[1]  !='' :self.project_label=  zonsurva[1] 
            if 'job.date'in zonsurva[0].lower() :
                if zonsurva[1]  !='' : self.acqdate =  zonsurva[1] 
                   
        self.sconfig_dict={'Job.Name':self.project_name , 
                  'Job.Area':self.project_area , 
                  'Job.For':self.custumer_name, 
                  'Job.By':self.contractor_name, 
                  'Job.Number':self.project_label , 
                  'Job.Date': self.acqdate, 
                  }
        


class SurveyConfiguration(object) : 
    """
    Class for survey Survey configuration.
         
    Arguments  
    ---------
        **survey_config_data** : list 
                                configuration list
        
    =================  ==========  ============================================
    Attributes         Type        Explanation
    =================  ==========  ============================================
    surveyType          str         survey type (CSAMT, TEM, CR, TDIP)        
    surveyArray         str         array type(for CSAMT:Scalar, Vector, Tensor)
    lineName            str         line label 
    lineNumber          float       line number 
    azimuth             str         line azimuth ( deg E of N) 
    stnGdpBeg           float       first GDP station number 
    StnGdpInc           float       GDP station number increment
    stnBeg              float       possibly rescaled first station number
    stnInc              float       possibly rescaled station number increment
    stnLeft             float       rescaled station number on  left edge of
                                    pseudosection plot 
    stnRight            float       rescaled station number on right  edge of  
                                    pseudosection plot 
    unitLength          float       length units (m,ft)
    unitEmag            str         field units (nV/m,nV/Am)
    unitHfield          str         field units (pT,pT/A
    unitPhase           str         units (mrad,deg)
    sconfig_dict        dict        surveyconfiguration dict elmts
    =================  ==========  ============================================
        
    More attributes can be added by inputing a key word dictionary
    
    :Example:
        
        >>> from csamtpy.ff.core import Avg
        >>> avg_obj= Avg(data_fn=path) 
        >>> surv_nameLine= avg_obj.Header.SurveyConfiguration.lineName 
        ... surv_nameLine
        ... avg_obj.Header.SurveyConfiguration.sconfig_dict['Stn.Right']
        
    """
    def __init__(self, survey_config_data =None , **kwargs):
        
        self.zon_surv_config = survey_config_data 
        self.surveyType ='CSAMT, TEM, CR, TDIP' 
        self.surveyArray ='CSAMT: Scalar, Vector, Tensor' 
        self.lineName='pyCSAMTline' 
        self.lineNumber='float:00' 
        self.azimuth =' 0 deg E of N'
        self.stnGdpBeg=None 
        self.stnGdpInc =None 
        self.stnBeg =None 
        self.stnInc=None 
        self.stnLeft=None 
        self.stnRight =None 
        self.unitLength ='m,ft' 
        self.unitEmag ='nV/mA or micoV/kmA' 
        self.unitHfield ='pT, picoT/A' 
        self.unitPhase ='mrad,deg'
        
        
        self.sconfig_dict= {'Survey.Type':self.surveyType,
                                    'Survey.Array':self.surveyArray , 
                                    'Line.Name':self.lineName, 
                                    'Line.Number':self.lineNumber , 
                                    'Line.Azimuth':self.azimuth, 
                                    'Stn.GdpBeg':self.stnGdpBeg , 
                                    'Stn.GdpInc':self.stnGdpInc, 
                                    'Stn.Beg':self.stnBeg , 
                                    'Stn.Inc':self.stnInc , 
                                    'Stn.Left':self.stnLeft , 
                                    'Stn.Right':self.stnRight, 
                                    'Unit.Length':self.unitLength, 
                                    'Unit.E':self.unitEmag, 
                                    'Unit.B':self.unitHfield, 
                                    'Unit.Phase':self.unitPhase, 
                                    }
        
        for keys in list(kwargs.keys()):
            # lista=[key for key in Test.__dict__.keys() if not (key.startswith('__') and key.endswith('__'))]
            if keys not in SurveyConfiguration.__dict__.keys():
                self.__setattr__(keys, kwargs[keys])
                self.sconfig_dict.__setitem__(keys, kwargs[keys])
                
    def set_survey_configuration_infos (self, survey_config_data =None): 
        """
        Method to set _survey configurations informations . 
        
        Parameters 
        ----------
            * survey_config_data : list or pathLike str 
                                container of survey configurations  infos. 
        """
        _logger.info('Reading and setting survey configurations informations!')

        
        if survey_config_data is not None : 
            self.zon_surv_config = survey_config_data
        if self.zon_surv_config is None :
            raise CSex.pyCSAMTError_Header("No Survey configuration data found."
                                             " Please check your configuration file.")
        if type(self.zon_surv_config) is not list: 
            try :
                if os.path.isfile(self.zon_surv_config) is True :
                    with open (self.zon_surv_config, 'r', encoding ='utf8') as fid :
                        self.zon_surv_config=fid.readlines()
            except : raise CSex.pyCSAMTError_Header(
                    "your configuration data must be either a list or a file.")
            
        # print(self.zon_surv_config)
        for configInfos in self.zon_surv_config :
            configInfos=configInfos.strip().split('=')

            if 'survey.type' in configInfos[0].lower():
                if configInfos[1] != '': self.surveyType=configInfos[1]
            if 'survey.array' in configInfos[0].lower():
                if configInfos[1] != '': self.surveyArray=configInfos[1]
            if 'line.name' in configInfos[0].lower():
                if configInfos[1] != '': self.lineName=configInfos[1] 
            if 'line.number' in configInfos[0].lower():
                if configInfos[1] != '': self.lineNumber=configInfos[1] 
            if 'line.azimuth' in configInfos[0].lower():
                if configInfos[1] != '': self.azimuth=configInfos[1] 
            if 'stn.gdpbeg' in configInfos[0].lower():self.stnGdpBeg=configInfos[1] 
            if 'stn.gdpinc' in configInfos[0].lower():self.stnGdpInc=configInfos[1] 
            if 'stn.beg' in configInfos[0].lower(): self.stnBeg=configInfos[1] 
            if 'stn.inc' in configInfos[0].lower(): self.stnInc=configInfos[1] 
            if 'stn.left' in configInfos[0].lower():self.stnLeft=configInfos[1] 
            if 'stn.right' in configInfos[0].lower(): self.stnRight=configInfos[1]
            if 'stn.length' in configInfos[0].lower(): self.unitLength=configInfos[1] 
            if 'unit.e' in configInfos[0].lower():
                if configInfos[1] != '': self.unitEmag=configInfos[1] 
            if 'unit.b' in configInfos[0].lower():
                if configInfos[1] != '': self.unitHfield=configInfos[1] 
            if 'unit.phase' in configInfos[0].lower():
                if configInfos[1] != '': self.unitPhase=configInfos[1] 

        self.sconfig_dict= {'Survey.Type':self.surveyType,
                                    'Survey.Array':self.surveyArray , 
                                    'Line.Name':self.lineName, 
                                    'Line.Number':self.lineNumber , 
                                    'Line.Azimuth':self.azimuth, 
                                    'Stn.GdpBeg':self.stnGdpBeg , 
                                    'Stn.GdpInc':self.stnGdpInc, 
                                    'Stn.Beg':self.stnBeg , 
                                    'Stn.Inc':self.stnInc , 
                                    'Stn.Left':self.stnLeft , 
                                    'Stn.Right':self.stnRight, 
                                    'Unit.Length':self.unitLength, 
                                    'Unit.E':self.unitEmag, 
                                    'Unit.B':self.unitHfield, 
                                    'Unit.Phase':self.unitPhase, 
                                    }
        
        
class TransmitterProperties(object):
    
    """
    Class for transmitter  properties.
     
    Arguments  
    -----------
        **data_fn** : str  
                path to avgfile
        
    =================  ==========  ============================================
    Attributes         Type        Explanation
    =================  ==========  ============================================
    txType              str         source type (for NSAMT /CSAMT:Natural, 
                                     Bipole,Loop)
    txGdpStn            str         transmitter ID from GDP Tx field 
                                    (GDP stn #) (alias = XMTR)
    txStn               str         rescaled Tx ID 
                                    (rescaled client stn #)
    txCenter            float       transmitter center-point easting,
                                    northing,elevation
                                    (float, length units)
                                    (aliases = TxCX, TxCY)
    txHPR               int         transmitter orientation 
                                    heading, pitch, 
                                    roll (Tx heading alias = TxBrg)
                                    (heading = bipole azimuth, pitch=0, roll=0 
                                    for z+up or 180 for z positive down)
    txLength            str         transmitter bipole length or square
                                    loop width 
                                    (positive float, length units)
    sconfig_dict        dict        params config dictionnary  
                                    
    =================  ==========  ============================================
    
    More attributes can be added by inputing a key word dictionary
     
    :Example:
        
        >>> from csamtpy.ff.core import Avg
        >>> path= os.path.join(os.environ['pyCSAMT'], 
        ...                         'csamtpy', data, 'LCS01.AVG')
        >>> avg_obj= Avg(data_fn=path) 
        >>> avg_obj.Header.Tx.txGdpStn
        >>> avg_obj.Header.Tx.txStn
        >>> avg_obj.Header.Tx.sconfig_dict['Tx.Type']
        
     """
    def __init__(self, **kwargs):

    
        self.txType ='NSAMT /CSAMT: Natural, Bipole,Loop' 
        self.txGdpStn ='XMTR ' 
        self.txStn='Tx ID' 
        self.txCenter='TxCX, TxCY' 
        self.txHPR='bipole azimuth, pitch=0, roll=0 for z+up or 180 for z positive down' 
        self.txLength='positive float, length units'
        

        
        for keys in list(kwargs.keys()):
            if keys not in TransmitterProperties.__dict__.keys():
                self.__setattr__(keys, kwargs[keys])
                self.sconfig_dict.__setitem__(keys, kwargs[keys]) 
                
    def set_transmitter_properties(self, Tx_data =None):
        """
        Method to set_Tx_properties.
        
        :param tx_data: list of Tx-infos from AVG filename 
        :type tx_data:list 
        
        """
        
        if Tx_data is None : 
            CSex.pyCSAMTError_inputarguments('No Tx-proprerties found !')
        
        if Tx_data is not None : 
            for tx_infos in Tx_data : 
                if '$ XMTR' in tx_infos or 'Tx.GdpStn' in tx_infos: 
                    self.txGdpStn=np.float(tx_infos.strip().split('=')[-1])
                if '$Tx.Stn' in tx_infos : 
                    self.txStn =np.float(tx_infos.strip().split('=')[-1])
                if 'Tx.Type' in tx_infos :
                    self.txType= tx_infos.strip().split('=')[-1]
                if 'Tx.HPR' in tx_infos :self.txHPR=[float(ss) for ss in \
                                                     tx_infos.strip().split('=')[-1].split(',')]
                if 'Tx.Length'in tx_infos : self.txLength=float(
                                                     tx_infos.strip().split('=')[-1].split()[-2])
                
         # --> set configuration -dictionnary               
        self.sconfig_dict={'Tx.Type':self.txType, 
                  'Tx.GdpStn':self.txGdpStn, 
                  'Tx.Stn':self.txStn, 
                  'Tx.Center':self.txCenter, 
                  'Tx.HPR': self.txHPR,
                  'Tx.Length':self.txLength,
              }
        
class ReceiverProperties(object):
    
    """
    Class for receiver  properties.
      
    Arguments 
    ---------- 
        **Rx_data** : list or str 
                    Filled automatically or path to your
                    receiver properties file.

    =================  ==========  ============================================
    Attributes         Type        Explanation
    =================  ==========  ============================================
    rxGdpStn            str         receiver GDP station number (alias = Station) 
    rxStn               str         receiver client station number
    rxHPR               int         EM component heading, pitch, roll 
                                    (floats, Rx heading  alias = RxBrg, ExAzm)
                                    (heading = x component azimuth in degrees 
                                    east of north, pitch = x  component
                                    orientation in degrees.  up from 
                                    horizontal, pitch = z cmp  rotation about 
                                    x cmp,  0=z+up, 180 = z+ down)
    rxLength            str         E-field dipole length or  loop widths(
                                    positive float(s),length units) 
                                    (alias = aspace)
                                    (positive float, length units)
    rxComps             str         EM component/impedance label(ExHx, ExHy, 
                                    EyHx, EyHy, Zxx, Zxy, Zyx, Zyy,  Zvec, Zdet)   
    sconfig_dict        dict        parameters dict config            
    =================  ==========  ============================================
         
    More attributes can be added by inputing a key word dictionary
     
    :Example:
         
        >>> from csamtpy.ff.core import Avg
        >>> path= os.path.join(os.environ['pyCSAMT'], 
                                 'csamtpy', data, 'LCS01.AVG')
        >>> avg_obj= Avg(data_fn=path) 
        ... avg_obj.Header.Rx.txGdpStn
        ... avg_obj.Header.Rx.txStn
        ... avg_obj.Header.Rx.sconfig_dict['Tx.Length']
     """
    def __init__(self, Rx_data=None ,  **kwargs):

        self.rx_data =Rx_data 
        self.rxGdpStn =None 
        self.rxStn=None 
        self.rxHPR=None 
        self.rxLength=None
        self.rxComps='ExHy'
        self.unit=None
        
        self.sconfig_dict={'Rx.GdpStn':self.rxGdpStn, 
                          'Rx.Stn':self.rxStn, 
                          'Rx.HPR':self.rxHPR , 
                          'Rx.Length': self.rxLength,
                          'Rx.Cmp':self.rxComps,
                          }


        for keys in list(kwargs.keys()):
            if keys not in ReceiverProperties.__dict__.keys():
                self.__setattr__(keys, kwargs[keys])
                self.sconfig_dict.__setitem__(keys, kwargs[keys])
                
                
    def set_receiver_properties(self, Rx_data =None ):
        """
        Methods  to set Receivers - properties infos from AVG files. 

        Parameters
        ----------
        * Rx_data : list, optional
                    container infos of receivers. The default is None.
        """
        if Rx_data is not None : 
            self.rx_data = Rx_data
        
        if self.rx_data is None : 
            CSex.pyCSAMTError_inputarguments('No Zonge Rx-proprerties found !')
            
        if type(self.rx_data) is not list :
            try : 
                if os.path.isfile(self.rx_data) ==True : 
                    with open(self.rx_data, 'r', encoding='utf-8') as frx:
                        self.rx_data =frx.readlines()
            except : raise CSex.pyCSAMTError_Header('Argument provided for rx_data are wrong !'\
                                                    ' Must be a list of file.')
                        

        for rx_infos in self.rx_data : 
            
            if '$Rx.Stn' in rx_infos : self.rxStn =np.float(rx_infos.strip().split('=')[-1])
                
            if 'Rx.Cmp' in rx_infos :self.rxComps= rx_infos.strip().split('=')[-1]
                
            if 'Rx.HPR' in rx_infos :self.rxHPR=[float(ss) for ss in \
                                                 rx_infos.strip().split('=')[-1].split(',')]
            if ('$ ASPACE' in rx_infos) or  ('Rx.Length'in rx_infos) : 
                    #-- it seems sometimes file got the unit at the end either 'm' or ft'
                    # should remove it before converting on float
                tx = rx_infos.strip().split('=')[-1]
                if not 46 < ord(rx_infos.strip().split('=')[-1][-1]) < 58 :  
                    tx =func._remove_str_word(char=tx,word_to_remove=rx_infos.strip().split('=')[-1][-1],
                                              deep_remove=False)
                    self.rxLength ,self.unit =np.float(tx), rx_infos.strip().split('=')[-1][-1]
                else :self.rxLength=np.float(tx)
 
            if '$Rx.GdpStn' in rx_infos : self.rxGdpStn =np.float(rx_infos.strip().split('=')[-1])
                
        self.sconfig_dict={'Rx.GdpStn':self.rxGdpStn, 
                          'Rx.Stn':self.rxStn, 
                          'Rx.HPR':self.rxHPR , 
                          'Rx.Length': self.rxLength,
                          'Rx.Cmp':self.rxComps,
                          }

      

class Skip_flag (object) : 
    """
    skip flag values 
    
        - 0. drop data 
        - 1. skip data 
        - 2. good quality of data
    
    Arguments 
    ----------
        **skip_flag** : str or int , 
                    skip_flag value set automatically .
                    some astatic file doesnt mention the number of 
                    skip_flag .Default is '2'.
 
    :Example:
            
         >>> from csamtpy.ff.core import Avg
         >>> path= os.path.join(os.environ['pyCSAMT'], 
         ...                         'csamtpy', data, 'LCS01.AVG')
         >>> avg_obj= Avg(data_fn=path) 
         >>> avg_obj.Skip_flag.skip_flag
         >>> avg_obj.Skip_flag.get_skip_flag
    """
    
    def __init__(self, skip_flag=None , **kwargs):
        
        self.skip_flag= skip_flag
        
        self.get_skip_flag=None 
        
        self.skip_flag_dict ={'2':'Good data quality', 
                            '1': 'Data skipped', 
                            '0' :'Data rejected', 
                            '*': 'No data obtained'}
        
        
        for keys in list(kwargs.keys()): 
            setattr(self, keys , kwargs[keys])
 
        
        
    def setandget_skip_flag (self, skip_flag=None):
        """ 
        simple method to set and get skip_flag.
        
        Parameters
        ----------
             * skip_flag:  str 
        """
        
        if skip_flag is not None : 
            if skip_flag not in list(self.skip_flag_dict.keys()):
                raise CSex.pyCSAMTError_config_file('Wrong Input ! skip flag must be'
                                                    ' among : {0}'.format(list(self.skip_flag_dict.keys())))
            self.skip_flag = str(skip_flag)

        for keys , values  in self.skip_flag_dict.items(): 
            if keys==self.skip_flag : 
                self.get_skip_flag = values
            
  
class ZongeHardware(object):
    """
    Some features of zonge AMTAVG software. 
    
    Arguments
    ---------
        **zonge_hardw_infos** : list or str , 
                            filled automatically or path to your 
                            hardware setting file.
            
    =================  ==========  =========================================
    Attributes         Type        Explanation
    =================  ==========  =========================================
    version             str         version of AMTAVG software
    dated               str         date of building Avg file 
    processed           str         date of data preprocessing
    from_fld_file       str         parent of avg file 
    astactic_version    str         zonge astactic sofware version
    updated_data        str         date were data were updated 
    numfilterfreq       int         number of point that TMA used to filter 
    freq_filter         float       value of TMA filter frequency
    end_astatic_file    str         end file description.
    =================  ==========  =========================================
        
    More can be added by inputing a key word dictionary
    
    :Example:
         
        >>> from csamtpy.ff.core import Avg
        >>> path= os.path.join(os.environ['pyCSAMT'], 
        ...                     'csamtpy', data, 'LCS01.AVG')
        >>> avg_obj= Avg(data_fn=path) 
        >>> avg_obj.Header.HardwareInfos.freq_filter)
        >>> avg_obj.Header.HardwareInfos.sconfig_dict['AMTAVG'])
    """

    def __init__(self,zonge_hardw_infos=None ,  **kwargs):
        
        self.zonge_hardw_infos=zonge_hardw_infos
        self.version ='7.76'
        self.from_fld_file = None 
        self.dated=None 
        self.processed= None 
        self.astatic_version = 'v3.60'
        self.update_data =None 
        self.numfilterfreq=None
        self.freq_filter=None 
        self.sconfig_dict={'AMTAVG': [self.version,self.from_fld_file], 
                           'Dated':self.dated,
                           'Processed':self.processed, 
                           'ASTATIC':[self.astatic_version, self.update_data ],
                           "":'{0} TMA Filter at {1} hertz'.format(self.numfilterfreq,self.freq_filter )}
        self._writer_copyleft ='@by_pyCSAMT'
            
        for key in list (kwargs.keys()):
            if key not in ZongeHardware.__dict__.keys():
                setattr(key, kwargs[key])
        
    def set_harware_infos (self, zonge_hardw_infos =None ):
        """
        Method to set Harwares infos.
        
        Parameters 
        ----------
            * zonge_hardw_infos : list 
                                Hardware informations collected  
        """
        _logger.info ('Reading and setting Harwware informations !')
        
        if zonge_hardw_infos is not None : 
            self.zonge_hardw_infos=zonge_hardw_infos 
        if self.zonge_hardw_infos is None : raise CSex.pyCSAMTError_Header('No informations from hardware found !')
        
        if type (self.zonge_hardw_infos) is not list : 
            try : 
                with open(self.zonge_hardw_infos, 'r', encoding='utf8') as fzh : 
                    self.zonge_hardw_infos= fzh.readlines()
            except : raise CSex.pyCSAMTError_Header("Harware - infos must be either a list or  file !")

        for  infos in self.zonge_hardw_infos : 
            if 'AMTAVG' in infos  : 
                tem_h=infos.strip().split(',') 
                self.processed =func._remove_str_word(char=tem_h[-1], 
                                                   word_to_remove='Processed')
                self.dated=func._remove_str_word(char=tem_h [-2], 
                                                   word_to_remove='Dated')
   
                self.from_fld_file=tem_h[0].split(':')[-1]
                self.version = tem_h[0].split()[2]
            elif 'ASTATIC' in infos : 
                tem_as =infos.strip().split()
                self.astatic_version=tem_as[2]
                self.update_data =tem_as [-1]
            elif 'TMA' in infos :
                self.end_astatic_file =infos.strip().strip('\\')
                tem_tma=infos.strip().split()
                try : 
                  tem_tma= [float(ss) for ss in tem_tma ]
                  
                except : pass 
                self.numfilterfreq =  tem_tma[1]
                self.freq_filter =tem_tma[-2]
                
        self.sconfig_dict={'AMTAVG': [self.version,self.from_fld_file], 
                           'Dated':self.dated,
                           'Processed':self.processed, 
                           'ASTATIC':[self.astatic_version, self.update_data ], 
                           "":'{0} TMA Filter at {1} hertz'.format(self.numfilterfreq,
                                                                   self.freq_filter )}

class Data (object): 
    """
    AVG Data informations ,  Container of all data infos .
     
    Arguments
    ---------
        **data_array**  : np.ndarray  
                       Data collected form the site, It contains all the informations 
                       during survey. It could read a  user own data build 
                       by classifying data according the main AVG file 

    ==================  ==============  =======================================
    Attributes          Type            Explanation
    ==================  ==============  =======================================
    Station             (ndarray, 1)    station length infos . unit : m 
    Frequency           (ndarray, 1)    frequency using during survey ,unit : Hz  
    Amps                (ndarray, 1)    current induction values using during 
                                        survey , unit : Amp  
    Emag                (ndarray, 1)    Electric field magnitude(Unit.nV/m,nV/Am) 
    Ephz                (ndarray, 1)    Electric field phase  (unit. mrad)
    Hmag                (ndarray, 1)    Magnetic field magnitude (Unit.pT,pT/A)
    Hphz                (ndarray, 1)    Magnetic field phase (unit. mrad)
    Resistivity         (ndarray, 1)    Cagnard apparent Resistivity
                                        magnitude (unit.ohm.m)
    Phase               (ndarray, 1)    Impedance phase (unit. mrad)
    Z                   (ndarray, 1)    Impedance magnitude  (unit. km/s)
    pcEmag              (ndarray, 1)    relative |E| error (%)
    pcHmag              (ndarray, 1)    relative |H| error (%)
    sEphz               (ndarray, 1)    phase(E) error (unit.mrad)
    sHphz               (ndarray, 1)    phase(H) error (unit.mrad)
    sPhz                (ndarray, 1)    phase(Z) error (unit.mrad)
    pcRho               (ndarray, 1)    relative apparent resistivity
                                        (rho) error (%)
    added_astatic_data  (ndarray, 1)    static_corrected-apparent
                                        resistivity (unit.ohm.m)
    set_to_degree           bool        phase angle value :if set to True , 
                                        valueof phase will be on degree. make 
                                        sure to converte the data on radian 
                                        anglebefore not in mrad.
    ==================  ==============  =======================================
    
     .. notes :: when data come straightforwardly from AVG file ,
             data was set respecting the unit provided  User must takes account 
             for data processing .A little bit  conversion is needed.
             User can use a special function. Normalized_station_value parameter 
             normalize station value . if set to true station will started name
             from 0 to n
                ...
                
    More attributes can be added by inputing a key word dictionary

    :Example:
        
        >>> from csamtpy.ff.core import Avg
        >>> DATA=Data(data_array=avg_data )
        >>> stn=DATA.Station.name 
        >>> freq=DATA.Frequency.value
    """
    
    def __init__(self, data_array=None , **kwargs):
        self._logging =csamtpylog.get_csamtpy_logger(self.__class__.__name__)
        self._data_array = data_array 
        
        self.Station =Station()
        self.Frequency=Frequency()
        self.Amps =Amps()
        self.Emag=Emag()
        self.Ephz=Ephz()
        self.Hmag=Hmag()
        self.Hphz=Hphz()
        self.Resistivity=Resistivity()
        self.Phase =Phase()
        self.Z =Z_Tensor()
        self.pcEmag=pcEmag()
        self.sEphz=sEphz()
        self.pcHmag=pcHmag()
        self.sHphz=sHphz()  
        self.pcRho=pcRho()
        self.sPhz=sPhz()
        
   
        self.added_astatic_data= None 
        self._number_of_frequencies=None 
        self._num_of_stations=None 
        self.normalized_station_value=kwargs.pop('normalized_station_value', False)
        self.set_to_degree =kwargs.pop('set_phase_angle_to_deg', False)
        
 
        
        self._f = None 
        
        for key in list(kwargs.keys()):
            setattr(self, key, kwargs[key])
            
        if self._data_array is not None :
            self._avg_set_components_from_data()
        
    def _avg_set_components_from_data (self, data_array=None ,data_type=None, 
                                       add_astatic_data =None ):
        """
        Method to set 
        avg_data_components  from Avgfile 
        
        Parameters
        ----------
            * data_array : np.ndarray 
                        data recovered from avg file. 
                        It must be on numpy array. The default is None.
            * data_type : int, optional
                        show the type of Avg file, 
                        must be 1 for plainty file 2 for Astatic file. 
                         *Default* is None.
            * add_astatic_data : pd.coreDataFrame , optional
                                added infos from astatic file. 
                                *Default* is None.
        """
        if data_array is not None : 
            self._data_array =data_array
            
        if self._data_array is None : 
            raise CSex.pyCSAMTError_AvgData('No Data from avgfile to read.')
            
        if data_type is not None : 
            self._f =data_type
        if self._f not in [1,2]: 
            raise CSex.pyCSAMTError_AvgData(' May check your avg Datafile.')

        
        number_of_freq_array , rep =np.unique(self._data_array[:,1], return_counts=True)
        number_of_station_array , rep =np.unique(self._data_array[:,0], return_counts=True)

        self._number_of_frequencies, self._number_of_stations = number_of_freq_array.size, number_of_station_array.size
        # print(self.numfreq, self.numstation)

        if add_astatic_data is not None :
            added_astatic_data=add_astatic_data
            pdf_Z= added_astatic_data.drop('SRes', axis=1)
            pdf_Z.reset_index(drop=True, inplace =True )
            self.Z._read_zAS_data(z_abs_array =pdf_Z , number_of_frequencies=self._number_of_frequencies, 
                                 number_of_stations=self._number_of_stations )
            
            #--> set asticatic
            self.added_astatic_data=added_astatic_data['SRes'].to_numpy()

            
        def zstar_array_to_nan (zstar_array, nan_value=np.nan, keep_str =False): 
            """
            Parameters
            ----------
                * zstar_array : ndarray 
                            array contain unfloat converter value.
                            the unconverter value can be a '*'
                * nan_value : float or np.nan type  
                            nan_value could be any value either int, float or str.  
                            If The default is np.nan.
                * keep_str : bool, optional  
                            keep the str item on your array.
            Returns 
            -------
                array_like 
                    zstrar_array converted 
            """
            
            for kk , value in enumerate(zstar_array): 
                try : zstar_array[kk] =float(value)
                except :zstar_array[kk]=nan_value
            if type(nan_value) is str : keep_str=True
            if keep_str : return zstar_array 
            return np.array([float(kk) for kk in zstar_array])   
            
        # try to change the dtype of data if dtype of the dataset is not float 

        if self._data_array.dtype not in ['float', 'int']:
            data_transitory =[zstar_array_to_nan(zstar_array=rowline) for dd, rowline in enumerate(self._data_array)]
            self._data_array=func.concat_array_from_list(list_of_array=data_transitory)
  
        # populate classes .  
        #cut of special freq_array and stations
        # for more certaintly , any value will come must be check whether its possible to convert , if not value 
        #       will set to np.nan 

        self.Station._set_station_(station_data_array=zstar_array_to_nan (zstar_array=self._data_array[:,0]), 
                                   normalized_station_value=self.normalized_station_value )
        self.Frequency._set_frequency(freq_array=zstar_array_to_nan (zstar_array=self._data_array[:,1]))
        self.Amps._read_amps_data(amps_array=zstar_array_to_nan (zstar_array=self._data_array[:,2]), 
                                 number_of_frequencies=self._number_of_frequencies, 
                                 number_of_stations=self._number_of_stations)
        # Amps(amps_array=self._data_array[:,2])
        self.Emag._read_emag_data(e_mag_array=zstar_array_to_nan (zstar_array=self._data_array[:,3]),
                                 number_of_frequencies=self._number_of_frequencies, 
                                 number_of_stations=self._number_of_stations )
        # Emag(e_mag_array =self._data_array[:,3])
        self.Ephz._read_ephz_data(e_phz_array=zstar_array_to_nan (zstar_array=self._data_array[:, 4]),
                                 number_of_frequencies=self._number_of_frequencies, 
                                 number_of_stations=self._number_of_stations, to_degree =self.set_to_degree  )
        # Ephz(e_phz_array=self._data_array[:, 4])
        self.Hmag._read_hmag_data(h_mag_array=zstar_array_to_nan (zstar_array=self._data_array[:,5]),
                                 number_of_frequencies=self._number_of_frequencies, 
                                 number_of_stations=self._number_of_stations )
        # Hmag (h_mag_array =self._data_array[:,5])
        self.Hphz._read_hphz_data(h_phz_array=zstar_array_to_nan (zstar_array=self._data_array[:, 6]),
                                 number_of_frequencies=self._number_of_frequencies, 
                                 number_of_stations=self._number_of_stations, to_degree =self.set_to_degree  )
        # Hphz(h_phz_array=self._data_array[:,6])
        
        if self._f==1 : 
        #     Resistivity(res_array=self._data[:,7])
            self.Resistivity._read_res_data(res_array =zstar_array_to_nan (zstar_array=self._data_array[:,7]),
                                            number_of_frequencies=self._number_of_frequencies, 
                                 number_of_stations=self._number_of_stations)
        #     Phase (phase_array=self._data[:8])
            self.Phase._read_phase_data(phase_array=zstar_array_to_nan (zstar_array=self._data_array[:,8]),
                                            number_of_frequencies=self._number_of_frequencies, 
                                 number_of_stations=self._number_of_stations )
        #     sPhz(s_phz_data=self._data_array[:-1])
            self.sPhz._read_sPhz_data(sPhase_array=zstar_array_to_nan (zstar_array=self._data_array[:,-1]), 
                                      number_of_frequencies=self._number_of_frequencies, 
                                 number_of_stations=self._number_of_stations)
        #     pcRho(pc_rho_data=self._data[:,-2])
            self.pcRho._read_pcRes_data(pcRes_array =zstar_array_to_nan (zstar_array=self._data_array[:,-2]),
                                        number_of_frequencies=self._number_of_frequencies, 
                                 number_of_stations=self._number_of_stations)
        if self._f==2 : 
        #     Resistivity(res_array=self._data[:,8])
            self.Resistivity._read_res_data(res_array =zstar_array_to_nan (zstar_array=self._data_array[:,8]),
                                            number_of_frequencies=self._number_of_frequencies, 
                     number_of_stations=self._number_of_stations, Sres=self.added_astatic_data)
        #     Phase (phase_array=self._data[:7])
            self.Phase._read_phase_data(phase_array=zstar_array_to_nan (zstar_array=self._data_array[:,7]),
                                        number_of_frequencies=self._number_of_frequencies, 
                                 number_of_stations=self._number_of_stations )
        #     sPhz(s_phz_data=self._data_array[:-2])
            self.sPhz._read_sPhz_data(sPhase_array=zstar_array_to_nan (zstar_array=self._data_array[:,-2]),
                                      number_of_frequencies=self._number_of_frequencies, 
                                 number_of_stations=self._number_of_stations)
        #     pcRho(pc_rho_data=self._data[:,-1])
            self.pcRho._read_pcRes_data(pcRes_array =zstar_array_to_nan (zstar_array=self._data_array[:,-1]), 
                                        number_of_frequencies=self._number_of_frequencies, 
                                 number_of_stations=self._number_of_stations)
            
        # pcEmag(pc_emag_data=self._data[:,-6])
        self.pcEmag._read_pcEmag_data(pc_e_mag_array=zstar_array_to_nan (zstar_array=self._data_array[:,-6]), 
                                      number_of_frequencies=self._number_of_frequencies, 
                                 number_of_stations=self._number_of_stations)
        # sEphz(s_phz_data=self._data[:,-5])
        self.sEphz._read_sEphz_data (sEphz_array=zstar_array_to_nan (zstar_array=self._data_array[:,-5]),
                                     number_of_frequencies=self._number_of_frequencies, 
                                 number_of_stations=self._number_of_stations )
        # pcHmag(pc_hmag_data =self._data[:, -4])
        self.pcHmag._read_pcHmag_data(pc_h_mag_array=zstar_array_to_nan (zstar_array=self._data_array[:,-4]),
                                      number_of_frequencies=self._number_of_frequencies, 
                                 number_of_stations=self._number_of_stations)
        # sHphz(s_hphz_data=self._data[:,-3])
        self.sHphz._read_sHphz_data (sHphz_array=zstar_array_to_nan (zstar_array=self._data_array[:,-3]),
                                     number_of_frequencies=self._number_of_frequencies, 
                                 number_of_stations=self._number_of_stations )   
    
       

class  Station(object):
    """
    Stations informations 
    
    Arguments
    -----------
        **norm_station_value**: bool 
               If True, station will numbered starting by 0.*Default* is True

    ================  ===========  ============================================
    Attributes        Type          Explanation
    ================  ===========  ============================================
    value             nd.array      data of station  columns components 
    max               float         maximum distance of survey lines 
    min               float         minimum distance strated recorded data  
    names              str          names of stations name filled automaticaly
    end               str           last station name 
    lengh             float         lenght of the survey 
    unit              str           lines lenght units  may be on [m, km , ft]
    alldata           nd.array      ndarray of all data value from raw file
    loc               dict          main attribute location of station: eg
                                    loc ['S01'] show the distance of station at
                                    that point. 
    ================  ===========  ============================================

    More attributes can be added by inputing a key word dictionary

    :Example: 
        
        >>> from csamtpy.ff.core.avg import Station
        >>> path=os.path.join(os.environ["pyCSAMT"], 
        ...                      "csamtpy", "data", "K1.AVG")
        >>> station=Station(path)
        >>> print(station.loc['S07'])
    """
    
    def __init__(self, station_data_array =None,  **kwargs):

        self._logging=csamtpylog.get_csamtpy_logger(self.__class__.__name__)
        

        self.station_data_array =station_data_array
        

        self.normalized_station_value = kwargs.pop('normalized_station_value', False)
        self.rename_station =kwargs.pop('rename_sation', None)
        self.unit=kwargs.pop('unit_of_distance_measured', 'm')
        
        self.names =None 
        self.value =None 
        self.max=None 
        self.min =None 
        self.last_station_name =None 
        self.lenght =None 
        self.value =None 
        

    def _set_station_(self,station_data_array =None, rename_station =None , 
                      normalized_station_value=False )  : 
        """
        Method to set station on dictionnary of each distance. in fact to easily export for 
        plotting each station . 

        Parameters
        ----------
            * station_data_array : np.ndarray (ndarray,1), optional
                                station value recorded in the field at each point.
                                The *default* is None.
            * rename_station : np.ndarray(ndarray,1), optional
                                list of station name provided. 
                                *Default* is None.
            * normalized_station_value : bool, optional
                                start the value of station by 0.
                                *Default* is False.
        """
        
        if station_data_array  is not None : 
            self.station_data_array=np.array([np.float(kk) for kk in station_data_array ])
  
        if self.station_data_array is None : 
            raise CSex.pyCSAMTError_station('No stations data to read .')
    

        num_station_counts , repsta = np.unique (self.station_data_array, return_counts=True)
 
        #---> check if stations_data provided are each the same length. 
        if np.all(repsta, axis=0) != True : 
            raise CSex.pyCSAMTError_station('Stations provided must be the same lenght of reccurency.')
        if repsta[0] ==1 : self.value = self.station_data_array 
        elif  repsta[0] !=1  : self.value =num_station_counts
            
        self.value =np.array([np.float(ss) for ss in self.value])
        #--- convert station  unit -- Default is 'm'---
        if self.unit =='km': self.value, self.station_data_array = self.value /1e3, self.station_data_array/1e3
        if self.unit == 'ft': self.value, self.station_data_array=self.value /3.280839895, self.station_data_array/3.280839895
        if self.unit not in ['km', 'ft', 'm']:
            
            self._logging.warn('Station units provided is incorect.Try  "km" or "ft." Default is "m."')
            raise CSex.pyCSAMTError_station('Unit provided <{0}> doesnt not match correct units.'\
                                      'acceptable units are : "m", "km" or "ft"'.format(self.unit))
                    
            
        self.min, self.max, self.num_of_station =self.value.min(), self.value.max(), self.value.size
        self.length=np.abs(self.max-self.min)
        
        if normalized_station_value is not None :
            self.normalized_station_value=normalized_station_value 
        
        if self.normalized_station_value  : 
            self.value =np.apply_along_axis(lambda x : np.abs(x)-self.min, 0 , self.value)
            self.station_data_array=np.apply_along_axis(lambda x : np.abs(x)-self.min, 0, self.station_data_array)
            
        # the case where the user provide its own name
        if rename_station is not None :
            self.rename_station = rename_station 
            
        if self.rename_station is not None : 
            if type (self.rename_station) is list : 
                self.rename_station=np.array(self.rename_station)
            if self.value.size != self.rename_station.size : 
                raise CSex.pyCSAMTError_station('Stations rename array must have the same length as'\
                                                ' the aray_data provided. lenght or stations_data is :<{0}>'
                                                ' '.format(self.value.size))
            self.names =self.rename_station
        elif rename_station is None : 
            self.names =['S{:02}'.format(ii) for ii in range (self.value.size)]
        
        #---> duplicate stations names and put on dict
        if repsta[0] ==1 : 
             self.loc ={keys:values for keys, values in zip(self.names , self.value)}
        if repsta [0] !=1 : 
            # self.names =self.names * repsta [0]
            self.names.sort()
            # print(self.station_data_array)
            # self.station_data_array=self.station_data_array.reshape ((self.station_data_array.shape[0],1))
            
            # print(self.station_data_array)
            #--- build a commun array station-value and truncated on dict.
            # data_to_truncate =np.concatenate((tem_names, self.station_data_array), axis =1 )
            station_list_truncated =cfunc.truncated_data(data=self.station_data_array, 
                                                         number_of_reccurence= repsta [0])
                
            self.loc= {keys:values for keys, values in zip(self.names, station_list_truncated)}
            
            

class Frequency (object): 
    """
    Frequency informations - Frequency at which data 
    was measured (Hertz). Frequency on Hz 
        
    Arguments  
    ---------
        **freq_array**  : arrya_like 
                        frequence data array 
                        
        **normalize_freq_betw** : list  
                                must be on list of integer value 
                                If float value provided 
                                should be convert on list, same as tuple 
 
    ======================  ==========  =======================================
    Attributes              Type        Explanation
    ======================  ==========  =======================================
    value                   nd.array    frequency data columns
                                        components 
    max                     float       Hight frequency 
    min                     float       lower frequency  
    loc                     dict        main attribute : location 
                                        of station loc ['S01'] show the frequency 
                                        of station at that point. 
    normalize_freq_betw     list        must be a list of integer  be sure that
                                        the numberprovided will be power 
                                        to 10. eg .(start, end,  nfreq)   
                                        start - first number will be power by 10 
                                        0=1:10exp(0) end- last value to be 
                                        interpolated : 4 --> 10(exp)4
                                        nfreq- number of frequency  : eg : 17 
                                        *Default* is None 
    ======================  ==========  =======================================

    More attributes can be added by inputing a key word dictionary

    :Example:
        
        >>> from csamtpy.ff.core.avg import Avg
        >>> path=os.path.join(os.environ["pyCSAMT"], 
              "data", "avg", "K1.AVG")
        >>> freq_obj=Avg(path)
        >>> freq_obj =avg_obj.Data_section.Frequency.loc['S07']
        ... print(freq_obj)
    """
    
    def __init__(self, freq_array =None, 
                 normalize_freq_betw =None,  **kwargs):

        self._freq_array =freq_array 
        
        self.normalize_freq_betw=kwargs.pop('normalize_freq_betw', None )
        self.value =None 
        self.min=None 
        self.max=None 
        self.numfreq=None 
        
        
        for keys in list (kwargs.keys()): 
            setattr(self, keys, kwargs[keys])
            
        if self._freq_array is not None : 
            self._set_frequency()
        
            
    def _set_frequency (self, freq_array =None, normalize_freq_betw = None ):
        """
        Method to set frequency according to each station.

        Parameters
        ----------
            * freq_array : np.ndarray(ndarray, 1)
                        frequency data provided . The default is None.
            * normalize_freq_betw : list 
                        list of frequency range to normalize, optional. 
                        The *default* is None.
        """
        
        
        if freq_array is not None : 
            self._freq_array =freq_array 
            
        if self._freq_array is None : 
            raise CSex.pyCSAMTError_frequency('No Frequency Data found.')
            
            
        self._freq_array =np.array([np.float(kk) for kk in self._freq_array])  
        
        #--> read frequency from AVG data file : 
        vacount_freq, freq_repeat =np.unique(self._freq_array, return_counts=True)

        
        self.numfreq = vacount_freq.size
        if np.all(freq_repeat, axis =0) != True : 
            raise CSex.pyCSAMTError_frequency('Problem occured when reading frequency data.'\
                                              'All frequency on Avgfile Must be the same length. ')
        # self.value =vacount_freq
        
        # if normalize_freq_betw is not None : 
        #      self.value =self.normalize_freq_betw()
        
        # self.max, self.min=self.value.max(), self.value.min()
  
        truncated_freq =cfunc.truncated_data( data =self._freq_array , 
                                       number_of_reccurence=self.numfreq )
        
        name, polyname = cfunc._numbering_station(number_of_station=freq_repeat[0], 
                               number_of_freq =self.numfreq)
        
 
        self.loc ={ key:value for key, value in zip(name, truncated_freq )}
        
        self.value = self.loc[name[0]]
        if normalize_freq_betw is not None : 
             self.value =self.normalize_freq_betw()
        
        self.max, self.min=self.value.max(), self.value.min()
        
        
        
    
    def normalize_frequency (self , normalize_freq_betw=None ):
        """
        method to interpolate frequencies.
        
        .. deprecated:: waist_method will replaced on cs module.
 
        Raises
        ------
            CSex
                raise Error if input arguments for frequency interpolating is wrong.

        Returns
        -------
            array_like
                frequency normalized.
        """
        
        if normalize_freq_betw is not None : 
            self.normalize_freq_betw = normalize_freq_betw
        
        if self.normalize_freq_betw is not None : 
            if type(self.normalize_freq_betw) ==tuple :
                self.normalize_freq_betw =list(self.normalize_freq_betw)
                try :
                    self.normalize_freq_betw =[int(ss) for ss in self.normalize_freq_betw]
                except :
                    self._logging.warn ("Input interpolated arguments are wrong! arguments type must"\
                                        " a list of integers. ")
                    raise CSex.pyCSAMTError_inputarguments("Input value of frequency is wrong. Must be"\
                                                  " a list of integers.")
            lengh_interp =np.linspace(self.normalize_freq_betw[0],self.normalize_freq_betw[1],
                                      self.normalize_freq_betw[2])
        freq_new=np.power(10, lengh_interp)
        # range value provide if value not in order .  
        freq_new.sort()
        # freq_new=func.interpol_scipy(x_value=np.array([ii for ii in range (self.numfreq)]),
        #                              y_value =self.vacount_freq, x_new=)
        return freq_new
        

class Comp (object):
    """
    Components measured
    
    Holds the following information:
         
    Arguments 
    ----------
        **comp_name** : str
                       path to avg filename 
        **new_component** :str  
                       component type which data will acquired . make sure 
                       once change it , the method to calculate impedande Z and 
                       rotating angle will also change to
                       
                       - ExHy Zxy (Ex/Hy)
                       - EyHz Zyx (Ey/Hx)
                       - EyHy Zyy (Ey/Hy)
                       - ExHx Zxx (Ex/Hx)

    =================  ==========  ============================================
    Attributes         Type        Explanation
    =================  ==========  ============================================
    name               str          name of components will filled automatically
    length             float        lenght of the survey lines 
    =================  ==========  ============================================

    More attributes can be added by inputing a key word dictionary
    
    :Example:
        
        >>> from csamtpy.ff.core.avg import Comp
        >>> path=os.path.join(os.environ["pyCSAMT"], 
                          "csamtpy", "data", "K1.AVG")
        >>> component=Comp(path, component_name='EyHx')
        >>> print(component.name)
    """
    component_type =['ExHy', "EyHx", 'EyHy', 'ExHx']
    
    def __init__(self, comp_name =None ,  **kwargs): 

        self.name =comp_name
        
        self.new_comp=kwargs.pop('component_name', None)
        if self.new_comp is not None :
            if type (self.new_comp) is str :
                if self.new_comp in self.component_type :
                    self.name =self.new_comp
                else : 
                    self._logging.warn ('Component provided is wrong !')
                    CSex.pyCSAMTError_inputarguments("Component provide as new component is wrong"\
                                             "list of components:{0}".format(self.component_type))
            else : 
                warnings.warn("Components must be a string : a list of arguments below:"\
                              "{0}".format(self.component_type))
                CSex.pyCSAMTError_inputarguments("Component type is wrong. must a string."\
                                             "list of components:{0}".format(self.component_type))
                    

class Amps (object): 
    """
    Square-Wave current (amps)
    
    Arguments :
    ----------
        **amps_array** : np.ndarray (ndarray, 1)
                        evolution of current data on the site 
        **number_of_stations** : int  
                            number of survey_stations . 
        **number_of_frequencies**: int 
                             number of frequencies during
                             the survey for each station location
            
    ================  ===========  ============================================
    Attributes         Type        Explanation
    ================  ===========  ============================================
    value             nd.array      data of Amp  column  on avgfile   
    max               float         maximum current enforce  a that point
    min               float         minimun current unit :amp (A)  
    loc               dict          main attribute location of station with 
                                    the current enforce a that point .eg loc 
                                    ['S07'] show the current Amps-data
                                    at that point
    ================  ===========  ============================================

    More attributes can be added by inputing a key word dictionary
    
    :Example:
         
        >>> from csamtpy.ff.core.avg import Avg
        >>> path=os.path.join(os.environ["pyCSAMT"], 
        ...      "csamtpy", "data", "K1.AVG")
        >>> avg_obj=Avg(path)
        >>> amps_obj =avg_obj.Data_section.Amps.loc['S00']
        >>> print(amp_obj)
    """

    def __init__(self, amps_array= None, number_of_frequencies=None , number_of_stations=None , **kwargs):

        
        self._amps_array =amps_array
        
        self.number_of_stations =kwargs.pop("number_of_stations", None)
        self.nfreq =kwargs.pop('number_of_frequencies', None)
        
        self.value =None 
        self.min = None 
        self.max = None 
    
        for key in list(kwargs.keys()):
            setattr(self, key, kwargs[key])
            
        if self._amps_array is not None : 
            self._read_amps_data()
        
        
    def _read_amps_data (self, amps_array= None, number_of_frequencies=None ,
                         number_of_stations=None ):
        """
        Method to read and arange data according each station /amps values .
        
        Parameters
        ----------
            * amps_array : ndarray, optional
                            amps _current observed sites. 
                            The *default* is None.
            * number_of_frequencies : int , optional
                            number of frequencies for survey.
                            The *default* is None.
            * number_of_stations : int , optional
                            number_of stations. 
                            The *default* is None.
        """
        if amps_array is not None : 
            self._amps_array =amps_array 
        if self._amps_array is None : 
            raise CSex.pyCSAMTError_inputarguments('No Ampers data !')
            
        if number_of_frequencies is not None : self.nfreq =number_of_frequencies 
        else :raise CSex.pyCSAMTError_inputarguments("Please specify the number of frequency !")
        if number_of_stations is not None : self.number_of_stations =number_of_stations 
        else : raise CSex.pyCSAMTError_inputarguments('Please specify the number of stations')
            
        try : 
            self.nfreq , self.number_of_stations =np.int(self.nfreq), np.int(self.number_of_stations)
            self._amps_array=np.array([np.float(cc) for cc in self._amps_array])

        except : raise CSex.pyCSAMTError_value("Values provided for the amp are wrong."\
                                          " must be float or integer .")
            
        vcounts_amps, repeat_amp =np.unique (self._amps_array, return_counts=True)
        self.value =vcounts_amps
        
        self.max, self.min =vcounts_amps.max(), vcounts_amps.min()
        
        
        truncated_amps=cfunc.truncated_data( data =self._amps_array, 
                                       number_of_reccurence=self.nfreq)
        name, polyname = cfunc._numbering_station(number_of_station=self.number_of_stations, 
                               number_of_freq =self.nfreq)
        self.loc ={ key:value for key, value in zip(name, truncated_amps )}
        
        
        
class Emag (object):
    """
    E-field magnitude (microVolts/(kiloMeter*Amp))
    
    Arguments 
    ----------
        **e_mag_array**: ndarray
                E_field magnitude Class  
                  
    ================  ===========  ============================================
    Attributes         Type        Explanation
    ================  ===========  ============================================
    value             nd.array      data of Emag  column on avgfile   
    max               float         maximum E-Field magnitude unit :mV/km*Amp 
    loc               dict          location of E-Field  magnitude value at 
                                    the station lambda.e:g : loc['S00'] show 
                                    the current  E-mag data of station 
                                    at that point.
    ================  ===========  ============================================

    More attributes can be added by inputing a key word dictionary
    
    :Example: 
        
        >>> from csamtpy.ff.core.avg import Avg
        >>> path=os.path.join(os.environ["pyCSAMT"], 
        ...      "csamtpy", "data", "K1.AVG")
        >>> avg_obj=Avg(path)
        >>> emag_obj =avg_obj.Data_section.Emag.loc['S02']
        ... print(emag_obj)
    """
    def __init__(self, e_mag_array=None, number_of_frequencies=None , number_of_stations=None , **kwargs):
        
        self.e_mag_array =e_mag_array
        
        self.number_of_stations =kwargs.pop("number_of_stations", None)
        self.nfreq =kwargs.pop('number_of_frequencies', None)
        
        self.value =None 
        self.min = None 
        self.max = None 
    
        for key in list(kwargs.keys()):
            setattr(self, key, kwargs[key])
            
        if self.e_mag_array  is not None : 
            self._read_emag_data ()
        
        
    def _read_emag_data (self,e_mag_array=None, number_of_frequencies=None , 
                         number_of_stations=None ):
        """
        Method to read and arange data according each 
        station /e_field magnitude values .
            
        Parameters
        ----------
            * e_mag_array : ndarray, optional
                            B-Field magnitude observed at each station.
                            The *default* is None.
            * number_of_frequencies : int , optional
                                number of frequencies for survey. 
                                The default is None.
            * number_of_stations : int , optional
                                number_of stations. The default is None.
        """
        if e_mag_array is not None : 
            self.e_mag_array =e_mag_array 
        if self.e_mag_array is None : 
            raise CSex.pyCSAMTError_inputarguments('No E-Field  data found  !')
            
        if number_of_frequencies is not None : self.nfreq =number_of_frequencies 
        else :raise CSex.pyCSAMTError_inputarguments("Please specify the number of frequency !")
        if number_of_stations is not None : self.number_of_stations =number_of_stations 
        else : raise CSex.pyCSAMTError_inputarguments('Please specify the number of stations')
            
        try : 
            self.nfreq , self.number_of_stations =np.int(self.nfreq), np.int(self.number_of_stations)
            self.e_mag_array=np.array([np.float(cc) for cc in self.e_mag_array])

        except : raise CSex.pyCSAMTError_value("Values provided for the amp are wrong."\
                                          " must be float or integer .")
            
        vcounts_e_mag, repeat_amp =np.unique (self.e_mag_array, return_counts=True)
        self.value =vcounts_e_mag
        
        self.max, self.min =vcounts_e_mag.max(), vcounts_e_mag.min()
        
        
        truncated_e_mag=cfunc.truncated_data( data =self.e_mag_array, 
                                       number_of_reccurence=self.nfreq)
        name, polyname = cfunc._numbering_station(number_of_station=self.number_of_stations, 
                               number_of_freq =self.nfreq)
        
        self.loc ={ key:value for key, value in zip(name, truncated_e_mag )}
        
    
class Ephz (object):
    """
    E-field phase (milliRadians)
    
    Arguments 
    ----------
        **e_phz_array** : ndarray
                   E_phase data field
        **to_deg**  : bool 
                   put the Ephz value on degree . 
                  
    ================  ===========  ============================================
    Attributes         Type        Explanation
    ================  ===========  ============================================
    value             nd.array      data of Emag  column  unit - mrad 
    max               float         maximum E-Field magnitude  unit -mV/km*Amp 
    loc               dict          location of E-phase   value at the
                                    station lambda . e.g - loc['S46] show the 
                                    current  E-phase data  at the station S46
    ================  ===========  ============================================

    More attributes can be added by inputing a key word dictionary
    
    :Example:
        
        >>> from csamtpy.ff.core.avg import Avg 
        >>> path=os.path.join(os.environ["pyCSAMT"], 
        ...      "csamtpy", "data", "K1.AVG")
        >>> avg_obj=Avg(path)
        >>> ephz_obj =avg_obj.Data_section.Ephz.loc['46']
        >>> print(ephz_obj)
    """
    def __init__(self, e_phz_array=None, number_of_frequencies=None , 
                 number_of_stations=None , **kwargs):
        
        self. e_phz_array= e_phz_array
        
        self.number_of_stations =kwargs.pop("number_of_stations", None)
        self.nfreq =kwargs.pop('number_of_frequencies', None)
        self.to_degree =kwargs.pop('to_degree', False)
        
        self.value =None 
        self.min = None 
        self.max = None 
    
        for key in list(kwargs.keys()):
            setattr(self, key, kwargs[key])
            
        if self.e_phz_array  is not None : 
            self._read_ephz_data ()
        
        
    def _read_ephz_data (self,e_phz_array=None, number_of_frequencies=None , 
                         number_of_stations=None, **kwargs  ):
        """
        Method to read and arange data according each
        station /E-phase angles values .
            
        Parameters
        ----------
            * e_phz_array : ndarray, en rad 
                            E-Field phase observed .
                            The *default* is None.
            * number_of_frequencies : int , optional
                                number of frequencies for survey.
                                The *default* is None.
            * number_of_stations : int , optional
                                number_of stations. 
                                The *default* is None.
            * to_degree : bool,
                        compute angle to degree . 
        """
        to_degree =kwargs.pop("to_degree", False)
        if to_degree :
            self.to_degree = to_degree
        
        if e_phz_array is not None : 
            self.e_phz_array =e_phz_array
        if self.e_phz_array is None : 
            raise CSex.pyCSAMTError_inputarguments('No E-phase data found  !')

        if number_of_frequencies is not None : self.nfreq =number_of_frequencies 
        else :raise CSex.pyCSAMTError_inputarguments("Please specify the number of frequency !")
        if number_of_stations is not None : self.number_of_stations =number_of_stations 
        else : raise CSex.pyCSAMTError_inputarguments('Please specify the number of stations')
            
        try : 
            self.nfreq , self.number_of_stations =np.int(self.nfreq), np.int(self.number_of_stations)
            self.e_phz_array=np.array([np.float(cc) for cc in self.e_phz_array])

        except : raise CSex.pyCSAMTError_value("Values provided for the E-phase are wrong."\
                                          " must be float or integer .")
        if self.to_degree : #---> set angle to degree .
            self.e_phz_array =np.apply_along_axis (lambda ephz: ephz * 180/np.pi, 0,self.e_phz_array)        
        
            
        vcounts_e_phz, repeat_ephz=np.unique (self.e_phz_array, return_counts=True)
        self.value =vcounts_e_phz
        
        self.max, self.min =vcounts_e_phz.max(), vcounts_e_phz.min()
        
        
        truncated_e_phz=cfunc.truncated_data( data =self.e_phz_array, 
                                       number_of_reccurence=self.nfreq)
        name, polyname = cfunc._numbering_station(number_of_station=self.number_of_stations, 
                               number_of_freq =self.nfreq)
        
        self.loc ={ key:value for key, value in zip(name, truncated_e_phz )}
        
        
class Hmag (object):
    """
    H-field magnitude (picoTesla/Amp)
    (milliGammas/Amp)  
    
    Arguments 
    ----------
        **h_mag_array** : ndarray
                       H/B_field magnitude Class . 
                  
    ================  ===========  ============================================
    Attributes         Type        Explanation
    ================  ===========  ============================================
    value             nd.array      data of Hmag  column  on avgfile  maximum
    max               float         H-mag magnitude unit -  mV/km*Amp 
    loc               dict          location of H-mag magnitude value at  the 
                                    station lambda . e:g : loc['S05']
                                    show the current  H-mag data
                                    of station at S05.
    ================  ===========  ============================================

    More attributes can be added by inputing a key word dictionary
    
    :Example: 
         
        >>> from csamtpy.ff.core import Avg 
        >>> path=os.path.join(os.environ["pyCSAMT"], 
        ...      "csamtpy", "data", "K1.AVG")
        >>> avg_obj=Avg(path)
        >>> hmag_obj =avg_obj.Data_section.Hmag.loc['S05']
        ... print(hmag_obj)
    """
    
    def __init__(self, h_mag_array=None, number_of_frequencies=None , number_of_stations=None , **kwargs):
        
        self.h_mag_array =h_mag_array
        
        self.number_of_stations =kwargs.pop("number_of_stations", None)
        self.nfreq =kwargs.pop('number_of_frequencies', None)
        
        self.value =None 
        self.min = None 
        self.max = None 
    
        for key in list(kwargs.keys()):
            setattr(self, key, kwargs[key])
            
        if self.h_mag_array  is not None : 
            self._read_hmag_data ()
        
        
    def _read_hmag_data (self, h_mag_array=None, number_of_frequencies=None , 
                         number_of_stations=None ):
        """
        Method to read and arange data according 
        each station /B_field magnitude  values .
           
        Parameters
        ----------
            * h_mag_array : ndarray, optional
                            B Field measured data on the sites.
                            The default is None.
            * number_of_frequencies : int , optional
                                    number of frequencies for survey. 
                                    The default is None.
            * number_of_stations : int , optional
                                number_of stations. The default is None.
        """
        
        if h_mag_array is not None : 
            self.h_mag_array =h_mag_array 
            
        if self.h_mag_array is None : 
            raise CSex.pyCSAMTError_inputarguments('No B-Field  data found  !')
            
        if number_of_frequencies is not None : self.nfreq =number_of_frequencies 
        else :raise CSex.pyCSAMTError_inputarguments("Please specify the number of frequency !")
        if number_of_stations is not None : self.number_of_stations =number_of_stations 
        else : raise CSex.pyCSAMTError_inputarguments('Please specify the number of stations')
            
        try : 
            self.nfreq , self.number_of_stations =np.int(self.nfreq), np.int(self.number_of_stations)
            self.h_mag_array=np.array([np.float(cc) for cc in self.h_mag_array])

        except : raise CSex.pyCSAMTError_value("Values provided for the B-Field are wrong."\
                                          " must be float or integer .")
            
        vcounts_h_mag, repeat_amp =np.unique (self.h_mag_array, return_counts=True)
        self.value =vcounts_h_mag
        
        self.max, self.min =vcounts_h_mag.max(), vcounts_h_mag.min()
        
        
        truncated_h_mag=cfunc.truncated_data( data =self.h_mag_array, 
                                       number_of_reccurence=self.nfreq)
        name, polyname = cfunc._numbering_station(number_of_station=self.number_of_stations, 
                               number_of_freq =self.nfreq)
        
        self.loc ={ key:value for key, value in zip(name, truncated_h_mag )}
            
class Hphz (object):
    """
    H-field phase (milliRadians)
    
    Arguments 
    ---------
        **h_phz_array** : ndarray
                        E_phase data field
         **to_degree**  : bool 
                    put the Hphz value on degree . 
                  
    ================  ===========  ============================================
    Attributes         Type        Explanation
    ================  ===========  ============================================
    value             nd.array      data of Emag  column  unit : mrad   
    max               float         maximum E-Field magnitude unit :mV/km*Amp 
    loc               dict          location of E-phase  value at the station
                                    lambda , e.g loc['S43] show the
                                    current  E-phase data at the stationS43
    ================  ===========  ============================================

    More attributes can be added by inputing a key word dictionary
    
    :Example: 
        
        >>> from csamtpy.ff.core.avg import Data_section 
        >>> path=os.path.join(os.environ["pyCSAMT"], 
              "csamtpy", "data", "K1.AVG")
        >>> avg_obj=Avg(path)
        >>> hphz_obj =avg_obj.Data_section.Ephz.loc['43']
        ... print(hphz_obj)
    """
    def __init__(self, h_phz_array=None, number_of_frequencies=None , 
                 number_of_stations=None , **kwargs):
        
        self.h_phz_array= h_phz_array
        
        self.number_of_stations =kwargs.pop("number_of_stations", None)
        self.nfreq =kwargs.pop('number_of_frequencies', None)
        self.to_degree =kwargs.pop('to_degree', False)
        
        self.value =None 
        self.min = None 
        self.max = None 
    
        for key in list(kwargs.keys()):
            setattr(self, key, kwargs[key])
            
        if self. h_phz_array is not None : 
            self._read_hphz_data ()
        
        
    def _read_hphz_data (self,h_phz_array=None, number_of_frequencies=None , 
                         number_of_stations=None, **kwargs  ):
        """
        Method to read and arange data according each
        station /E-phase angles 
            values .
        
        Parameters
        ----------
            * e_phz_array : ndarray, en rad 
                            E-Field phase observed . The default is None.
            * number_of_frequencies : int , optional
                            number of frequencies for survey. The default is None.
            * number_of_stations : int , optional
                            number_of stations.
                            The default is None.
            * to_degree : bool,
                            compute angle todegree . 
        """
        to_degree =kwargs.pop("to_degree", False)
        if to_degree :
            self.to_degree = to_degree
        
        if h_phz_array is not None : 
            self.h_phz_array=h_phz_array
        if self.h_phz_array is None : 
            raise CSex.pyCSAMTError_inputarguments('No E-phase data found  !')

        if number_of_frequencies is not None : self.nfreq =number_of_frequencies 
        else :raise CSex.pyCSAMTError_inputarguments("Please specify the number of frequency !")
        if number_of_stations is not None : self.number_of_stations =number_of_stations 
        else : raise CSex.pyCSAMTError_inputarguments('Please specify the number of stations')
            
        try : 
            self.nfreq , self.number_of_stations =np.int(self.nfreq), np.int(self.number_of_stations)
            self.h_phz_array=np.array([np.float(cc) for cc in self.h_phz_array])

        except : raise CSex.pyCSAMTError_value("Values provided for the E-phase are wrong."\
                                          " must be float or integer .")
        if self.to_degree : #---> set angle to degree .
            self.h_phz_array =np.apply_along_axis (lambda hphz: hphz * 180/np.pi, 0,self.h_phz_array)        
        
            
        vcounts_h_phz, repeat_hphz=np.unique (self.h_phz_array, return_counts=True)
        
        if not np.all(repeat_hphz): CSex.pyCSAMTError_value('Values of B-phases provided must have the same length'\
                                                            ' for each stations. ')
        self.value =vcounts_h_phz
        
        self.max, self.min =vcounts_h_phz.max(), vcounts_h_phz.min()
        
        
        truncated_h_phz=cfunc.truncated_data( data =self.h_phz_array, 
                                       number_of_reccurence=self.nfreq)
        name, polyname = cfunc._numbering_station(number_of_station=self.number_of_stations, 
                               number_of_freq =self.nfreq)
        
        self.loc ={ key:value for key, value in zip(name, truncated_h_phz )}
        
        
        
class Resistivity (object):
    """
    based on Cagniard Resistivity (Ohm-Meters) calculation 
    
     .. see also::Cagniard, L., 1953, Basic theory of the magnetotelluric
          method of  geophysical prospecting: Geophysics, 18, 605–635,
          *doi*: 10.1190/1.1437915. 
    
    Arguments 
    ----------
        **res_array**  : str
                      data array of apparents resistivy calculation.
                  
    ================  ===========  ============================================
    Attributes         Type        Explanation
    ================  ===========  ============================================
    value             nd.array      data of resistivity column - ohm.m 
    max               float         higher value of Resistivity 
    min               float         minimum value of resitivity on the area
    mean              float         everage value of the rho 
    loc               dict          location of rho  value at the station lambda 
                                    e.g - loc['S13] show the current  rho data 
                                    at the station S13
    Sres              ndarray       Astatic value corrected
    ================  ===========  ============================================

    More attributes can be added by inputing a key word dictionary
    
    :Example:  
         
        >>> from csammtpy.ff.core.avg import Avg 
        >>> path=os.path.join(os.environ["pyCSAMT"], 
        ...                  "csamtpy", "data", "K1.AVG")
        >>> avg_obj=Avg(path)
        >>> rho_obj =avg_obj.Data_section.Resistivity.loc['S00']
        >>> print(rho.loc['S00'])
    """
    
    def __init__(self, res_array =None , number_of_frequencies=None , 
                 number_of_stations=None , **kwargs):

        
        self.res_array=res_array 
        
        self.Sres=None 
        
        self.number_of_stations =kwargs.pop("number_of_stations", None)
        self.nfreq =kwargs.pop('number_of_frequencies', None)
        self.Sres =kwargs.pop('added_astatic_rho_array', None)
        
        self.value =None 
        self.min = None 
        self.max = None 
        self.loc = None 
        for key in list(kwargs.keys()):
            setattr(self, key, kwargs[key])
            
        if self.res_array is not None : 
            self._read_res_data ()
        
        
    def _read_res_data (self,res_array =None, number_of_frequencies=None , 
                         number_of_stations=None, Sres=None,  **kwargs  ):
        """
        Method to read resitivity  according each station /Resistivity angles 
            values .
        
        Parameters
        ----------
            * res_array : array_like 
                         Resistivity calculated  on the field. 
                         The *default* is None.
            * number_of_frequencies : int , optional
                            number of frequencies for survey. 
                            The *default* is None.
            * number_of_stations : int , optional
                            number_of stations. 
                            The *default* is None.
            * Sres : ndarray, 
                        Zonge Astatic rho calculated.
        """
        
        if res_array  is not None : 
            self.res_array= res_array
        if self.res_array is None : 
            raise CSex.pyCSAMTError_inputarguments('No Resistivity data found !')

        if number_of_frequencies is not None : self.nfreq =number_of_frequencies 
        else :raise CSex.pyCSAMTError_inputarguments("Please specify the number of frequency !")
        if number_of_stations is not None : self.number_of_stations =number_of_stations 
        else : raise CSex.pyCSAMTError_inputarguments('Please specify the number of stations')
        
        if Sres is not None :

            self.Sres = np.array([np.float(res) for res in Sres])
            if self.res_array.shape[0] != self.Sres.size : raise CSex.pyCSAMTError_rho('Resistivity calculated & '\
                                                                                 'Astatic_array must get the same length!.')
        try : 
            self.nfreq , self.number_of_stations =np.int(self.nfreq), np.int(self.number_of_stations)
            self.res_array =np.array([np.float(cc) for cc in self.res_array])

        except : raise CSex.pyCSAMTError_value("Values provided for the Resistivities are wrong."\
                                          " must be float or integer .")
        
        vcounts_res, repeat_hphz=np.unique (self.res_array, return_counts=True)
        if not np.all(repeat_hphz): CSex.pyCSAMTError_value('Values of Rho provided must have the same length'\
                                                           ' for each stations. ')
        self.value =vcounts_res
        self.max, self.min =vcounts_res.max(), vcounts_res.min()
        self.mean=self.res_array.mean()
        truncated_res=cfunc.truncated_data( data =self.res_array, 
                                       number_of_reccurence=self.nfreq)
        name, polyname = cfunc._numbering_station(number_of_station=self.number_of_stations, 
                               number_of_freq =self.nfreq)
        
        self.loc ={ key:value for key, value in zip(name, truncated_res )}
        
        if self.Sres is not None : 
            vcounts_sres, repeat_hphz=np.unique (self.Sres, return_counts=True)
            self.vSres , self.mean_Sres = vcounts_sres, self.Sres.mean()
            self.max, self.min =(vcounts_res.max(),vcounts_sres.max()), (vcounts_res.min(),vcounts_sres.min())
            self.mean =(self.mean, self.mean_Sres)
            truncated_sres=cfunc.truncated_data( data =self.Sres, 
                                       number_of_reccurence=self.nfreq)
            self.loc_Sres ={ key:value for key, value in zip(name, truncated_sres )}

       
class Phase (object): 
    """
    Impedance phase  on milliRadians can be calculate using Ephz
     object and Hphz object 
     
         * **Phase** = E-phase - H-phase
        
    Arguments 
    ----------
        **phase_array**  : ndarray
                        array of phase values got on the field
         **to_degree** : bool 
                        put the Phase value on degree . 
                  
    ================  ===========  ============================================
    Attributes         Type        Explanation
    ================  ===========  ============================================
    value             nd.array      data of phase  column  
    max               float         maximum phase value (mrad)
    min               float         minimum  Phase  (mrad)
    loc               dict          location of H-phase  value  at the 
                                    station lambda . e:g : loc['S00'] show the 
                                    current phase data at the station S00
    ================  ===========  ============================================

    More attributes can be added by inputing a key word dictionary

    :Example:
         
        >>> from csamtpy.ff.core.avg import Avg 
        >>> path=os.path.join(os.environ["pyCSAMT"], 
        ...                  "csamtpy", "data", "K1.AVG")
        >>> avg_obj=Avg(path)
        >>> phs_obj =avg_obj.Data_section.Phase.loc['S00']
        ... print(phs.loc['S00'])
    """
    
    
    def __init__(self, phase_array =None,number_of_frequencies=None , 
                         number_of_stations=None,   **kwargs):

        self._phase_array=phase_array 

        
        self.number_of_stations =kwargs.pop("number_of_stations", None)
        self.nfreq =kwargs.pop('number_of_frequencies', None)
        self.to_degree =kwargs.pop('to_degree', False)
        
        self.value =None 
        self.min = None 
        self.max = None 
    
        for key in list(kwargs.keys()):
            setattr(self, key, kwargs[key])
            
        if self._phase_array is not None : 
            self._read_res_data ()
            
   
    def _read_phase_data (self,phase_array =None,number_of_frequencies=None , 
                     number_of_stations=None, **kwargs):
        """
        method to read phase values obtained on the sites.units on mrads.  
        Zonge keep the phase values on mrad.
        
        Parameters
        ------------
            * phase_array : ndarray, optional
                            phase data aray on the field.
                            The *default* is None.
            * number_of_frequencies : TYPE, optional
                            number of fequency used .
                            The *default* is None.
            * number_of_stations : int, optional
                            number of the stations on the site.
                            The *default* is None.
        """
    
 
        to_degree =kwargs.pop("to_degree", False)
        
        if to_degree :
            self.to_degree = to_degree
        
        if phase_array is not None : 
            self._phase_array =phase_array
        if self._phase_array is None : 
            raise CSex.pyCSAMTError_inputarguments('No Phase data found  !')

        if number_of_frequencies is not None : self.nfreq =number_of_frequencies 
        else :raise CSex.pyCSAMTError_inputarguments("Please specify the number of frequency !")
        if number_of_stations is not None : self.number_of_stations =number_of_stations 
        else : raise CSex.pyCSAMTError_inputarguments('Please specify the number of stations')

            
        try : 
            self.nfreq , self.number_of_stations =np.int(self.nfreq), np.int(self.number_of_stations)
            self._phase_array=np.array([np.float(cc) for cc in self._phase_array])

        except : raise CSex.pyCSAMTError_value("Values provided for the Phase are wrong."\
                                          " must be float or integer .")
        if self.to_degree : #---> set angle to degree .
            self._phase_array =np.apply_along_axis (lambda phz: phz * 180/np.pi, 0,self._phase_array)        
        
            
        vcounts_phz, repeat_phz=np.unique (self._phase_array, return_counts=True)
        self.value =vcounts_phz
        
        self.max, self.min =vcounts_phz.max(), vcounts_phz.min()
        
        
        truncated_phz=cfunc.truncated_data( data =self._phase_array, 
                                       number_of_reccurence=self.nfreq)
        name, polyname = cfunc._numbering_station(number_of_station=self.number_of_stations, 
                               number_of_freq =self.nfreq)
        
        self.loc ={ key:value for key, value in zip(name, truncated_phz )}
            
        
class Z_Tensor(object):
    """
    Impedance Tensor Z Calculation :  
       class can recompute the apparent resistity rho base on impedance Tensor Z.
       
    .. seealso:: Zonge, K.L. and Hughes, L.J., 1991, Controlled source audio-frequency 
                magnetotellurics,in Electromagnetic Methods in Applied Geophysics, ed.
                Nabighian, M.N., Vol. 2,Society of Exploration Geophysicists, pp. 713-809.
            
    Arguments 
    ----------
        **data_fn** : str
                path to avg filename. 
                  
    ================  ===========  ===========  ===============================
    Attributes         Type          units      Explanation
    ================  ===========  ===========  ===============================
    Z                 np.ndarray    V/m         Impedance tensor 
    real              dict          V/m         real part of Z 
    imag              dict                      imaginary part of Z  
    complex           dict          complex     Z on complex type
    rho               dict          ohm.m       Resistivity 
    loc               dict                      location of Z value at 
                                                the station lambda .
                                                e:g loc['S00] show the
                                                current Z data at station S00
    ================  ===========  ===========  ===============================
    
    .. note:: all dict-attributes  can be loc with the station name .
        
    More attributes can be added by inputing a key word dictionary
    
    :Example:
         
        >>> from csamtpy.ff.core.avg impor Z_tensor
        >>> path=os.path.join(os.environ["pyCSAMT"], 
        ...                  "csamtpy", "data", "K1.AVG")
        >>> z_obj=Z_Tensor(path)
        >>> z_z=Zxy.Z
        >>> print(z_obj.real["S00"])
        >>> print(z_obj.imag["S00"])
        >>> print(z_obj.complex["S00"])
        >>> print(z_obj.loc["S00"])
        >>> print(z_obj.rho["S00"])
    """
    def __init__(self, z_array =None , phase_array =None, freq=None, z_error=None, 
                 **kwargs): 

        self._logging =csamtpylog.get_csamtpy_logger(self.__class__.__name__)
        self._z  =z_array
        self._phase =phase_array 
        self._z_err =z_error 
        self._freq = None
        
 
        self._rho =None
        self.loc=None 

        self.zAS =None 
        
        
        for key in list (kwargs.keys()):
            setattr(self, key , kwargs[key])
        
        self.z_pwgt =None 
        self.z_megt =None 
            
        self.z_as_labels=['Z.mag', 'Z.mwgt','Z.pwgt',
                                  'E.wgt', 'B.wgt','Z.%err' ]
        
             
    @property 
    def freq (self): 
        return self._freq 
    
    @freq.setter 
    def freq (self, freq_array):
        if freq_array.dtype not in ['float', 'int']:
            try : 
                self._freq =np.array([np.float(ff) for ff in freq_array])
            except :raise CSex.pyCSAMTError_frequency('wrong Input frequencies'\
                                                      ' values.must be a float number. ')
            
    @property 
    def phase (self):
        return self._phase 
    @phase.setter 
    def phase(self, phz_array):
        if phz_array.dtype not in ['float', 'int']:
            try :self._phase =np.array([float(phz) for phz in phz_array])
            except : CSex.pyCSAMTError_Phase('Arguments phase values must be int, or float .')
            
    @property 
    def z_error (self):
        return self._z_err 
    @z_error.setter 
    def z_error(self, z_error_):
        if z_error_.dtype not in ['float', 'int']:
            try : 
                self._z_err =np.array([float(zz) for zz in z_error_])
            except : CSex.pyCSAMTError_Z("Error z input values are incorrects.")
    
    @property 
    def z  (self): 
        return self._z
    
    @z.setter 
    def z (self, zz_array):
        if zz_array.dtype == 'str':
            try : 
                zz_array= np.array([np.float(zz) for zz in zz_array])
            except : CSex.pyCSAMTError_Z('z_impedance value must be a complex_number.')
        if zz_array.dtype in ['float', 'int']: #---> provide freq value and phase . 
            self._z = Zcc.rhophi2z(phase =self._phase , freq =self._freq , 
                                   z_array=zz_array)[-1]
        else : self._z= zz_array
    
                
    
    def z_and_zerr_2rhophi(self, z_array=None , freq=None):
        """
        Method to compute resistivity and phase phase 
        using Z_values and Zerror _values 
        
        Parameters
        ----------
            * z_array : complex
                    Impedance tensor complex_number 
                    The *default* is None.
            * freq : ndarray, 
                    frequency value.
                    The *default* is None.

        Returns
        -------
            ndarray(ndarray,1)  
                resistivity value computed in ohm.m
            ndarray (ndarray, 1), 
                phase angle value  in degree.
        """
        
        
        if z_array is not None : self.z= z_array 
        
        if freq is not None : self.freq=freq
        
        if self.z is None or self.freq is None : 
            raise CSex.pyCSAMTError_Z('None values can not be computed. Check values !')
            
        self.rho =np.apply_along_axis(lambda valz : np.abs(valz)**2 * (0.2 / self.freq), 0, self.z)
        
        #-- z is complex number , can be compute using np.angle

        self.phase= np.rad2deg(np.angle(self.z))
            
        return self.rho , self.phase
    
    
    @property
    def rho (self): 
        return self._rho 
        
    @rho.setter 
    def rho(self, res_array):
        if res_array.dtype not in ['float', 'int']:
            try : self._rho = np.array([float(res) for res in res_array]) 
            except : raise CSex.pyCSAMTError_rho('Resistivities values must be float '\
                                                 'or integer , not a None type !')
                
    
        
    def rhophi2rhoph_errors (self, res_array=None , phase_array = None,
                             z_error =None , freq=None ): 
        """
        compute the phase and resistivities error via res_array , phase_array and _z error.
        
        Parameters
        ----------
            * res_array : ndarray , 
                        resistivity value in ohm.m 
            * phase_array : ndarray , 
                        phase angle value in mradians 
            * z_error : ndarray , 
                        impedance Tensor error 
            * freq : ndarray , 
                        frequency numbers in Hz 
        Returns 
        --------
           ndarrays 
              resistivities error and phase errors values  .
        """
        
        
        self._logging.info('compute rho_phi to resistivities and phases_errors')
        
        if res_array is not None : self.rho =res_array
        if z_error is not None : self.z_error =z_error 
        if phase_array is not None : self.phase =phase_array 
        if freq is not None : self.freq = freq 
        if self.rho is None  or self.phase is None or self.freq is None : 
            self._logging.warn('NoneType can not be computed. please check your data arrays.')
            raise CSex.pyCSAMTError_Z('could note compute a Nonetype number. please check numbers.')
        
        # compute imag part and real part of Z 
             
        if self.rho.size != self.freq.size : 
            raise CSex.pyCSAMTError_Z('Resistivity , freq and phase_array must be the same size. ')
        
        # zz_= np.array([np.sqrt(0.2 * self._freq[ii] * self._rho[ii]) for ii in range (self._freq.size)])
        zz_comp = Zcc.rhophi2z(phase=self.phase/1e3, freq=self.freq, 
                               resistivity=self.rho)
        self.z =zz_comp[-1]

        
        # zz_ =np.sqrt (0.2 * self._freq * self._rho)
        # z_real=(np.apply_along_axis(lambda xxr : np.cos(xxr), 0, self._phase)) * zz_ 
        # z_imag =(np.apply_along_axis(lambda xxi: np.cos(xxi), 0, self._phase))* zz_

        res_err, phase_err = Zcc.z_error2r_phi_error(z_real=zz_comp[1],
                                                     z_imag=zz_comp[2], 
                                                     error=self.z_error)
        
        return res_err , phase_err 
            
    def _read_zAS_data (self, z_abs_array =None , number_of_frequencies=None ,
                      number_of_stations =None ):
        """
        Method to set Z_values straightforwardly from Zonge Astatic file . 
       
        Parameters
        -----------
            * z_array : pd.DataFrame ,
                        added values from astatic files. 
            * number_of_frequency :  int  , 
                        number of frequency used during surveys.
            * number_of_stations :int , 
                        differents survey frequencies.
            * Frequency : array_like
                        frequency array units in Hz  
            
        .. Notes :: 
             appropriate converter  should be 
             - 1Hz -->  6.28140704, 
             - deg-rad --> theta(deg) *np.pi/180 
             - mu  vaccum permitivity with  scipy.constants.epsilon_0
        """
        
        self._logging.info('Transfering Astatic Z data .')
        
        if z_abs_array is not None : 
            self.zAS = z_abs_array 
        if self.zAS is None : 
            raise CSex.pyCSAMTError_inputarguments('No zonge Astatic  data found  !')            
            
        if number_of_frequencies is not None : self.nfreq = number_of_frequencies 
        else :raise CSex.pyCSAMTError_inputarguments("Please specify the number of frequency !")
        
        if number_of_stations is not None : self.number_of_stations= number_of_stations
        else : raise CSex.pyCSAMTError_inputarguments('Please specify the number of stations')

        
        if type(self.zAS) == pd.core.frame.DataFrame : 
            for wgt in list(self.zAS.columns): 
                if wgt == 'Z.mwgt' or wgt == 'E.wgt': 
                    self.z_mwgt =self.zAS[wgt].to_numpy()
                if wgt == 'Z.pwgt' or wgt =='B.wgt':
                    self.z_pwgt =self.zAS['Z.pwgt'].to_numpy()
                if wgt == 'Z.%err': 
                    self._z_err =self.zAS[wgt].to_numpy()
                if wgt == 'Z.mag' : zabs_array =self.zAS[wgt].to_numpy()
                
                
        try : 
            self.nfreq , self.number_of_stations =np.int(self.nfreq), np.int(self.number_of_stations)
            zabs_array=np.array([np.float(cc) for cc in zabs_array])
        except : raise CSex.pyCSAMTError_value("Values provided for Z-astatic are wrong."\
                                          " check your Zonge AVG file .")        
        
  
        self.max, self.min =zabs_array.max(), zabs_array.min()
        
        truncated_zz=cfunc.truncated_data( data =zabs_array, 
                                       number_of_reccurence=self.nfreq)
        name, polyname = cfunc._numbering_station(number_of_station=self.number_of_stations, 
                               number_of_freq =self.nfreq)
        
        self.loc ={ key:value for key, value in zip(name, truncated_zz )}
        

class pcEmag (object):
    """
    Statistical variation of magnitude values from 
    averaged data blocks.
    
    Standard Deviation/Average Emag (percent)
        
    Arguments 
    ----------
        **pcEmag** : ndarray
                    data array of statistical variataion
                    of Emag value on the field.
            
    ================  ===========  ============================================
    Attributes         Type        Explanation
    ================  ===========  ============================================
    value             nd.array      data of pcEmag  column   on avgfile
    max               float         maximum % value of pc Emag
    min               float         minimum % value of pcEmag 
    loc               dict          location of %Ephz magnitude
                                    value at the station lambda e:g .loc['S05']
                                    show the current  H-mag data 
                                    of station at S05.
    ================  ===========  ============================================

    More attributes can be added by inputing a key word dictionary
    
    :Example:
        
        >>> from csamtpy.ff.core.avg import Avg 
        >>> path=os.path.join(os.environ["pyCSAMT"], 
        ...      "csamtpy", "data", "K1.AVG")
        >>> avg_obj=Avg(path)
        >>> pcemag_obj =avg_obj.Data_section.pcEmag.loc['S05']
        >>> print(pcemag_obj)
    """
    def __init__(self, pc_e_mag_array=None, 
                 number_of_frequencies=None , number_of_stations=None , **kwargs):
        
        self._pcEmag =pc_e_mag_array
        
        self.number_of_stations =kwargs.pop("number_of_stations", None)
        self.nfreq =kwargs.pop('number_of_frequencies', None)
        
        self.value =None 
        self.min = None 
        self.max = None 
    
        for key in list(kwargs.keys()):
            setattr(self, key, kwargs[key])
            
        if self._pcEmag  is not None : 
            self._read_pcEmag_data ()
        
        
    def _read_pcEmag_data (self,pc_e_mag_array=None, number_of_frequencies=None , 
                         number_of_stations=None ):
        """
        Method to read and arange data according each station /Statistical variation of 
            E-Field magnitude values .
        
        Parameters
        ----------
            * pc_e_mag_array : ndarray, optional
                                E-Field std variation  observed at each station.
                                The default is None.
            * number_of_frequencies : int , optional
                                    number of frequencies for survey. 
                                    The default is None.
            * number_of_stations : int , optional
                                number_of stations. The default is None.
        """
        if pc_e_mag_array is not None : 
            self._pcEmag =pc_e_mag_array 
        if self._pcEmag is None : 
            raise CSex.pyCSAMTError_inputarguments('No E-mag statistical variation  data found  !')
            
        if number_of_frequencies is not None : self.nfreq =number_of_frequencies 
        else :raise CSex.pyCSAMTError_inputarguments("Please specify the number of frequency !")
        if number_of_stations is not None : self.number_of_stations =number_of_stations 
        else : raise CSex.pyCSAMTError_inputarguments('Please specify the number of stations')
            
        try : 
            self.nfreq , self.number_of_stations =np.int(self.nfreq), np.int(self.number_of_stations)
            self._pcEmag=np.array([np.float(cc) for cc in self._pcEmag])

        except : raise CSex.pyCSAMTError_value("Values provided for the E-mag stat.variation  are wrong."\
                                          " must be float or integer .")
            
        vcounts_Emag, repeat_Emag =np.unique (self._pcEmag, return_counts=True)
        self.value =vcounts_Emag
        
        self.max, self.min =vcounts_Emag.max(), vcounts_Emag.min()
        
        
        truncated_pcEmag=cfunc.truncated_data( data =self._pcEmag, 
                                       number_of_reccurence=self.nfreq)
        name, polyname = cfunc._numbering_station(number_of_station=self.number_of_stations, 
                               number_of_freq =self.nfreq)
        
        self.loc ={ key:value for key, value in zip(name, truncated_pcEmag)}
        
        
        
class sEphz(object) : 
    """
    Statistical variation of the data blocks averaged for this data point.
    100 * Standard Deviation of Ephz values (milliradians)
        
    Arguments 
    ----------
        **sEphz**  : ndarray,
                   data array of statistical variation of 
                   Ephase value on the field.
                  
    ================  ===========  ============================================
    Attributes         Type        Explanation
    ================  ===========  ============================================
    value             nd.array      data of sHphz  column  on avgfile   
    max               float         maximum % value of sHphz
    min               float         minimum  value of sHphz
    loc               dict          location of sHphz magnitude
                                    value at the station lambda
                                    e.g : loc['S44'] show the 
                                    current  sHphz data of station at S44
    ================  ===========  ============================================
    
    More attributes can be added by inputing a key word dictionary
    
    :Example:
         
        >>> from csamtpy.ff.core.avg import Avg 
        >>> path=os.path.join(os.environ["pyCSAMT"], 
              "csamtpy", "data", "K1.AVG")
        >>> avg_obj=Avg(path)
        >>> shphz_obj =avg_obj.Data_section.sHphz.loc['S44']
        ... print(shphz_obj)
    """
    def __init__(self, sEphz_array=None, number_of_frequencies=None , 
                 number_of_stations=None , **kwargs):
        
        self._sEphz= sEphz_array
        
        self.number_of_stations =kwargs.pop("number_of_stations", None)
        self.nfreq =kwargs.pop('number_of_frequencies', None)
        self.to_degree =kwargs.pop('to_degree', False)
        
        self.value =None 
        self.min = None 
        self.max = None 
    
        for key in list(kwargs.keys()):
            setattr(self, key, kwargs[key])
            
        if self._sEphz is not None : 
            self._read_sEphz_data ()
        
        
    def _read_sEphz_data (self, sEphz_array=None, number_of_frequencies=None , 
                         number_of_stations=None, **kwargs  ):
        """
        Method to read and arange data according each station /Statistical variation 
        of E-Field angles values 
 
        Parameters
        ----------
            * shphz_array : ndarray, 
                            en rad E-Field phase observed 
                            The default is None.
            * number_of_frequencies : int , optional
                                    number of frequencies for survey.
                                    The default is None.
            * number_of_stations : int , optional
                                number_of stations. The default is None.
            * to_degree : bool,
                        compute angle todegree . 
        """
        to_degree =kwargs.pop("to_degree", False)
        if to_degree :
            self.to_degree = to_degree
        
        if sEphz_array is not None : self._sEphz =sEphz_array
        if self._sEphz is None :  raise CSex.pyCSAMTError_inputarguments('No stat.variation'\
                                                                         ' B-field-phase data found  !')
        if number_of_frequencies is not None : self.nfreq = number_of_frequencies 
        else :raise CSex.pyCSAMTError_inputarguments("Please specify the number of frequency !")
        if number_of_stations is not None : self.number_of_stations =number_of_stations 
        else : raise CSex.pyCSAMTError_inputarguments('Please specify the number of stations')
        self.nfreq , self.number_of_stations =np.int(self.nfreq), np.int(self.number_of_stations)    
            
           
        def _inspect_input (input_obj , to_degree=False):
            """
            Ascertain that the value is possible for computing. and put on 
            ndarray (ndarray,1)
            """            
            try : 
                input_obj =np.array([float(ss) for ss in input_obj])
            except : raise  CSex.pyCSAMTError_value("Values provided for computing are wrong."\
                                          " must be float or integer .")
            
            if to_degree : input_obj =np.apply_along_axis (lambda x: x * 180/np.pi, 0, input_obj)
            
            return input_obj 

        self._sEphz= _inspect_input(self._sEphz)

        if self.to_degree : #---> set angle to degree .
            self._sEphz =_inspect_input(input_obj =self._sEphz, to_degree =True)        
        
        # --> check the maximum angle value and minimum angle value.
        vcounts_sEphz, repeat_hphz=np.unique (self._sEphz, return_counts=True)
        self.value =vcounts_sEphz
        self.max, self.min =vcounts_sEphz.max(), vcounts_sEphz.min()
        
        #---> tuncated AVG sHphz value on  dictionnary for easy acces . 
        truncated_sEphz=cfunc.truncated_data( data =self._sEphz, 
                                       number_of_reccurence=self.nfreq)
        name, polyname = cfunc._numbering_station(number_of_station=self.number_of_stations, 
                               number_of_freq =self.nfreq)
        
        self.loc ={ key:value for key, value in zip(name, truncated_sEphz)}
        
class pcHmag (object):
    """
    Statistical variation of magnitude values from averaged data blocks.
        Standard Deviation / Average Hmag (percent)
        
     Arguments :
     ----------
        **pcHmag**  : ndarray
                   data array of statistical variataion 
                   of Hmag value on the field.
                  
    ================  ===========  ============================================
    Attributes         Type        Explanation
    ================= =========== =============================================
    value             nd.array      data of pcHmag  column  on avgfile
    max               float         maximum % value of pc Hmag
    min               float         minimum % value of pcHmag 
    loc               dict          location of %Hphz magnitude
                                    value at the station lambda
                                    e:g . loc['S05'] show the current  H-mag 
                                    data  of station at S05.
    ================  ===========  ============================================

    More attributes can be added by inputing a key word dictionary
    
    :Example:
         
        >>> from csamtpy.ff.core.avg import Avg 
        >>> path=os.path.join(os.environ["pyCSAMT"], 
        ...      "csamtpy", "data", "K1.AVG")
        >>> avg_obj=Avg(path)
        >>> pchmag_obj =avg_obj.Data_section.pcHmag.loc['S05']
        >>> print(pchmag_obj)
        
    """
    
    def __init__(self, pc_h_mag_array=None, 
                 number_of_frequencies=None , number_of_stations=None , **kwargs):
        
        self._pcHmag =pc_h_mag_array
        
        self.number_of_stations =kwargs.pop("number_of_stations", None)
        self.nfreq =kwargs.pop('number_of_frequencies', None)
        
        self.value =None 
        self.min = None 
        self.max = None 
    
        for key in list(kwargs.keys()):
            setattr(self, key, kwargs[key])
            
        if self._pcHmag  is not None : 
            self._read_pcHmag_data ()
        
        
    def _read_pcHmag_data (self,pc_h_mag_array=None, number_of_frequencies=None , 
                         number_of_stations=None ):
        """
        Method to read and arange data according each station /Statistical variation of 
            B-Field magnitude values .
        
        Parameters
        ----------
            * pc_h_mag_array : ndarray, optional
                            E-Field std variation  observed at each station. 
                            The default is None.
            * number_of_frequencies : int , optional
                                    number of frequencies for survey. 
                                    The default is None.
            * number_of_stations : int , optional
                                number_of stations. The default is None.
        """
        if pc_h_mag_array is not None : 
            self._pcHmag =pc_h_mag_array 
        if self._pcHmag is None : 
            raise CSex.pyCSAMTError_inputarguments('No B-mag statistical variation  data found  !')
            
        if number_of_frequencies is not None : self.nfreq =number_of_frequencies 
        else :raise CSex.pyCSAMTError_inputarguments("Please specify the number of frequency !")
        if number_of_stations is not None : self.number_of_stations =number_of_stations 
        else : raise CSex.pyCSAMTError_inputarguments('Please specify the number of stations')
            
        try : 
            self.nfreq , self.number_of_stations =np.int(self.nfreq), np.int(self.number_of_stations)
            self._pcHmag=np.array([np.float(cc) for cc in self._pcHmag])

        except : raise CSex.pyCSAMTError_value("Values provided for the B-mag stat.variation  are wrong."\
                                          " must be float or integer .")
            
        vcounts_Hmag, repeat_Hmag =np.unique (self._pcHmag, return_counts=True)
        self.value =vcounts_Hmag
        
        self.max, self.min =vcounts_Hmag.max(), vcounts_Hmag.min()
        
        
        truncated_pcHmag=cfunc.truncated_data( data =self._pcHmag, 
                                       number_of_reccurence=self.nfreq)
        name, polyname = cfunc._numbering_station(number_of_station=self.number_of_stations, 
                               number_of_freq =self.nfreq)
        
        self.loc ={ key:value for key, value in zip(name, truncated_pcHmag)} 


class sHphz(object) : 
    """
    Statistical variation of the data blocks averaged for this data point.
        100 * Standard Deviation of Hphz values (milliradians)
        
    Arguments
    ----------
        **shphz**  : ndarray
                data array of statistical variation of Emag 
                value on the field.
                  
    ================  ===========  ============================================
    Attributes         Type        Explanation
    ================  ===========  ============================================
    value             nd.array      data of sHphz  column on avgfile   
    max               float         maximum % value of sHphz
    min               float         minimum  value of sHphz
    loc               dict          location of sHphz magnitude
                                    value at the station lambda. 
                                    e:g : loc['S43'] show the
                                    current  sHphz data of station at S43.
    ================  ===========  ============================================

    More attributes can be added by inputing a key word dictionary
    
    :Example: 
         
        >>> from csamtpy.ff.core.avg impor Avg 
        >>> path=os.path.join(os.environ["pyCSAMT"], 
        ...      "csamtpy", "data", "K1.AVG")
        >>> avg_obj=Avg(path)
        >>> shphz_obj =avg_obj.Data_section.sHphz.loc['S43']
        >>> print(shphz_obj)
        
    """
    def __init__(self, shphz_array=None, number_of_frequencies=None , 
                 number_of_stations=None , **kwargs):
        
        self._sHphz= shphz_array
        
        self.number_of_stations =kwargs.pop("number_of_stations", None)
        self.nfreq =kwargs.pop('number_of_frequencies', None)
        self.to_degree =kwargs.pop('to_degree', False)
        
        self.value =None 
        self.min = None 
        self.max = None 
    
        for key in list(kwargs.keys()):
            setattr(self, key, kwargs[key])
            
        if self._sHphz is not None : 
            self._read_sHphz_data ()
        
        
    def _read_sHphz_data (self, sHphz_array=None, number_of_frequencies=None , 
                         number_of_stations=None, **kwargs  ):
        """
        Method to read and arange data according each station /Statistical variation 
        of B-Field angles. 
            values .
        
        Parameters
        ----------
            * shphz_array : ndarray, en rad 
                            E-Field phase observed . 
                            The default is None.
            * number_of_frequencies : int , optional
                            number of frequencies for survey.
                            The default is None.
            * number_of_stations : int , optional
                         number_of stations. The default is None.
            * to_degree : bool,
                         compute angle todegree .
                     
        """
        to_degree =kwargs.pop("to_degree", False)
        if to_degree :
            self.to_degree = to_degree
        
        if sHphz_array is not None : self._sHphz =sHphz_array
        if self._sHphz is None :  raise CSex.pyCSAMTError_inputarguments('No stat.variation'\
                                                                         ' B-field-phase data found  !')
        if number_of_frequencies is not None : self.nfreq = number_of_frequencies 
        else :raise CSex.pyCSAMTError_inputarguments("Please specify the number of frequency !")
        if number_of_stations is not None : self.number_of_stations =number_of_stations 
        else : raise CSex.pyCSAMTError_inputarguments('Please specify the number of stations')
        self.nfreq , self.number_of_stations =np.int(self.nfreq), np.int(self.number_of_stations)    
            
           
        def _inspect_input (input_obj , to_degree=False):
            """
            Ascertain that the value is possible for computing. and put on 
            ndarray (ndarray,1)
            """            
            try : 
                input_obj =np.array([float(ss) for ss in input_obj])
            except : raise  CSex.pyCSAMTError_value("Values provided for computing are wrong."\
                                          " must be float or integer .")
            
            if to_degree : 
                input_obj = np.apply_along_axis (lambda x: x * 180/np.pi, 0, input_obj)
            
            return input_obj 

        self._sHphz= _inspect_input(self._sHphz)


        if self.to_degree : #---> set angle to degree .
            self._sHphz =_inspect_input(input_obj =self._sHphz, to_degree =True)        
        
        # --> check the maximum angle value and minimum angle value.
        vcounts_shphz, repeat_hphz=np.unique (self._sHphz, return_counts=True)
        self.value =vcounts_shphz
        self.max, self.min =vcounts_shphz.max(), vcounts_shphz.min()
        
        #---> tuncated AVG sHphz value on  dictionnary for easy acces . 
        truncated_shphz=cfunc.truncated_data( data =self._sHphz, 
                                       number_of_reccurence=self.nfreq)
        name, polyname = cfunc._numbering_station(number_of_station=self.number_of_stations, 
                               number_of_freq =self.nfreq)
        
        self.loc ={ key:value for key, value in zip(name, truncated_shphz)}
        
    
class pcRho (object):
    """
    Statistical variation of magnitude values 
    from averaged data blocks.
    Standard Deviation / Average Rho (percent)
        
    Arguments 
    ----------
        **pcRes** : ndarray
               statistical variation averaged rho.


    ================  ===========  ============================================
    Attributes         Type        Explanation
    ================  ===========  ============================================
    value             ndarray       data of pcRho column  on avgfile 
    max               float         maximum % value of pcRho
    min               float         minimu % value of pcRho 
    loc               dict          location of %Rho magnitude
                                    value at  the station lambda .
                                    e:g : loc['S07'] show  the current  H-mag 
                                    data of station at S07.
    ================  ===========  ============================================
    
     More attributes can be added by inputing a key word dictionary
    
    :Example: 
         
        >>> from csamtpy.ff.core.avg import Avg
        >>> path=os.path.join(os.environ["pyCSAMT"], 
                          "csamtpy", "data", "K1.AVG")
        >>> avg_obj=Avg(path)
        >>> pcrho_obj =avg_obj.Data_section.pcRho.loc['S00']
        ... print(pcrho.loc['S00'])
        
    """
    
    def __init__(self, pcRes_array =None ,
                 number_of_frequencies=None , 
                 number_of_stations=None , **kwargs):

        
        self._pcRes=pcRes_array
        
        self.number_of_stations =kwargs.pop("number_of_stations", None)
        self.nfreq =kwargs.pop('number_of_frequencies', None)
        self.Sres =kwargs.pop('added_astatic_rho_array', None)
        
        self.value =None 
        self.min = None 
        self.max = None 
    
        for key in list(kwargs.keys()):
            setattr(self, key, kwargs[key])
            
        if self._pcRes is not None : 
            self._read_pcRes_data ()
        
        
    def _read_pcRes_data (self, pcRes_array =None, number_of_frequencies=None , 
                         number_of_stations=None, Sres=None,  **kwargs  ):
        """
        Method to read resitivity  according each station /standard deviation *100
        for resistivity  values .
        
        Parameters
        ----------
            * pcRes_array: ndarray, en ohm.m
                            Resistivity calculated  on the field.
                            The default is None.
            * number_of_frequencies : int , optional
                            number of frequencies for survey. 
                            The default is None.
            * number_of_stations : int , optional
                            number_of stations. 
                            The default is None.
        """
        
        if pcRes_array  is not None : 
            self._pcRes= pcRes_array 
        if self._pcRes is None : 
            raise CSex.pyCSAMTError_inputarguments('No Resistivity data found !')

        if number_of_frequencies is not None : self.nfreq =number_of_frequencies 
        else :raise CSex.pyCSAMTError_inputarguments("Please specify the number of frequency !")
        if number_of_stations is not None : self.number_of_stations =number_of_stations 
        else : raise CSex.pyCSAMTError_inputarguments('Please specify the number of stations')
        
        try : 
            self.nfreq , self.number_of_stations =np.int(self.nfreq), np.int(self.number_of_stations)
            self._pcRes =np.array([np.float(cc) for cc in self._pcRes])

        except : raise CSex.pyCSAMTError_value("Values provided for the Resistivities are wrong."\
                                          " must be float or integer .")
        
        vcounts_pcres, repeat_hphz=np.unique (self._pcRes, return_counts=True)
        if not np.all(repeat_hphz): CSex.pyCSAMTError_value('Values of Rho provided must have the same length'\
                                                           ' for each stations. ')
        self.value =vcounts_pcres
        self.max, self.min =vcounts_pcres.max(), vcounts_pcres.min()
        self.mean=self._pcRes.mean()
        truncated_pcres=cfunc.truncated_data( data =self._pcRes, 
                                       number_of_reccurence=self.nfreq)
        name, polyname = cfunc._numbering_station(number_of_station=self.number_of_stations, 
                               number_of_freq =self.nfreq)
        
        self.loc ={ key:value for key, value in zip(name, truncated_pcres )}
        

class sPhz(object) : 
    """
    Statistical variation of the data blocks averaged for this data point.
        100 * Standard Deviation of Phase values (milliradians)
        
    Arguments 
    ----------
        **sphz array** : array_like 
                coefficient of  phase variation arrqy 
        **path to avg**:str
                full path to filename 
                  
    Holds the following informations:
        
    ================  ==========  =============================================
    Attributes         Type        Explanation
    ================  ==========  =============================================
    value             nd.array      data of sPhz  column on avgfile   
    max               float         maximum % value of sHphz
    min               float         minimu  value of sHphz
    loc               dict          location of sPhz magnitude value at 
                                    the station lambda .
                                    e:g - loc['S07'] show the current  sPhz data 
                                    of station at S07.
    ================  ==========  =============================================
    
    More attributes can be added by inputing a key word dictionary
    
    :Example:

        >>> from csamtpy.ff.core.avg impor Avg
        >>> path=os.path.join(os.environ["pyCSAMT"], 
        ...     "csamtpy", "data", "K1.AVG")
        >>> avg_obj=Avg(path)
        >>> shphz_obj =avg_obj.Data_section.sHphz.loc['43']
        ... print(shphz_obj)
    """
    def __init__(self, sPhase_array=None, number_of_frequencies=None , 
                 number_of_stations=None , **kwargs):
        
        self._sPhs= sPhase_array
        
        self.number_of_stations =kwargs.pop("number_of_stations", None)
        self.nfreq =kwargs.pop('number_of_frequencies', None)
        self.to_degree =kwargs.pop('to_degree', False)
        
        self.value =None 
        self.min = None 
        self.max = None 
    
        for key in list(kwargs.keys()):
            setattr(self, key, kwargs[key])
            
        if self._sPhs is not None : 
            self._read_sPhs_data ()
        
        
    def _read_sPhz_data (self, sPhase_array=None, number_of_frequencies=None , 
                         number_of_stations=None, **kwargs  ):
        """
        Method to read and arange data according each station /Statistical variation 
        of zPhase -Field angles. 

        Parameters
        ------------
            * shphz_array : ndarray, en rad 
                        E-Field phase observed . 
                        The default is None.
            * number_of_frequencies : int , optional
                        number of frequencies for survey.
                        The default is None.
            * number_of_stations : int , optional
                        number_of stations. 
                        The default is None.
            * to_degree : bool,
                        compute angle todegree . 
        """
        to_degree =kwargs.pop("to_degree", False)
        if to_degree :self.to_degree = to_degree
            
        
        if sPhase_array is not None : self._sPhs =sPhase_array
        if self._sPhs is None :  raise CSex.pyCSAMTError_inputarguments('No stat.variation'\
                                                                         ' Impedance Z -phase data found  !')
        if number_of_frequencies is not None : self.nfreq = number_of_frequencies 
        else :raise CSex.pyCSAMTError_inputarguments("Please specify the number of frequency !")
        if number_of_stations is not None : self.number_of_stations =number_of_stations 
        else : raise CSex.pyCSAMTError_inputarguments('Please specify the number of stations')
        self.nfreq , self.number_of_stations =np.int(self.nfreq), np.int(self.number_of_stations)    
            
           
        def _inspect_input (input_obj , to_degree=False):
            """
            Ascertain that the value is possible for computing. and put on 
            ndarray (ndarray,1) and convert angle as possible to degree.
            """            
            try : 
                input_obj =np.array([float(ss) for ss in input_obj])
            except : raise  CSex.pyCSAMTError_value("Values provided for computing are wrong."\
                                          " must be float or integer .")
            
            if to_degree : input_obj =np.apply_along_axis (lambda x: x * 180/np.pi, 0,input_obj)
            
            return input_obj 

        self._sPhs=_inspect_input(self._sPhs)


        if self.to_degree : #---> set angle to degree .
            self._sPhs =_inspect_input(input_obj =self._sPhs, to_degree =True)        
        
        # --> check the maximum angle value and minimum angle value.
        vcounts_sphz, repeat_hphz=np.unique (self._sPhs, return_counts=True)
        self.value =vcounts_sphz
        self.max, self.min =vcounts_sphz.max(), vcounts_sphz.min()
        
        #---> tuncated AVG sHphz value on  dictionnary for easy acces . 
        truncated_sphz=cfunc.truncated_data( data =self._sPhs, 
                                       number_of_reccurence=self.nfreq)
        name, polyname = cfunc._numbering_station(number_of_station=self.number_of_stations, 
                               number_of_freq =self.nfreq)
        
        self.loc ={ key:value for key, value in zip(name, truncated_sphz)}
        
        

        


    

    
    
    
    
        
    

