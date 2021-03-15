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
 
.. _module-Shifting::`csamtpy.ff.processing.corr` 
        :synopsis: Deal with all data files. It corrects apparent resistivity
                by shitibg value to rj static factor . Apply correction and filters
                Some filters applied to correct apparent resistivities are 
                TMA , AMA and FLMA. 
        ...
        
Created on Sat Dec 12 13:55:47 2020

@author: @Daniel03
"""

import os, warnings 
import numpy as np 
import scipy as sp
import matplotlib.pyplot as plt
import scipy.interpolate  as spi 

from csamtpy.etc import infos 
from csamtpy.ff.core  import avg as CSAMTavg
from csamtpy.ff.processing import zcalculator as Zcc
from csamtpy.utils._csamtpylog import csamtpylog
from csamtpy.utils import func_utils as func
from csamtpy.utils import exceptions as CSex

#-------------------- end import module ---------------------------
class shifting(object):
    """ 
    processing class : shifting processing workflow 
     coorection  class deal with  AVG Zonge station file "*.stn" or SEG-EDI file. 
     
    Arguments
    ----------
        **data_fn** : str 
                    path to Zonge *AVG file or SEG-EDI files or jfiles
            
        **res_array** : array_like (ndarray,1) 
                    apparent resistivities uncorrected data 
            
        **freq_array** : array_like (ndarray,1)
                   frequency array during survey 
           
        **phase_array** : array_like(ndarray,1)
                   phase array during survey 
            
    More attribute will populate : 
        
    =================  ==========  ============================================
    Attributes         Type        Explanation
    =================  ==========  ============================================
    _rj                array_like   static shift factor  
    rho_static         array_like   corrected data from files
    _reference_freq    float        reference value to start corrected data . 
                                    If provided , function will use for data  
                                    correction.if notprovided, program will 
                                    search the reference frequency automatically
                                    and set it for computation. 
    =================  ==========  ============================================
    
    ========================  =================================================
    Methods                    Description
    ========================  =================================================
    TMA                         Trimming Moving Average computation. 
    FLMA                        Fixed Length Dipole Moving average computation. 
    AMA                         Adaptative Moving average . see Biblio. 
    ========================  =================================================
    
    More attributes can be added by inputing a key word dictionary
    
     .. note:: Reference frequency doesnt need to  be on frequency range
               If value is not in frequency array, program will 
               interpolate value
    
   :Example:
       
       >>> from csamtpy.ff.processing.corr import Shifting
       >>> path =  os.path.join(os.environ["pyCSAMT"], 
       ...                       'csamtpy','data', LCS01.AVG)
       ... static_cor =shifting().TMA (data_fn=path,
       ...                                 reference_freq=1024.,
       ...                                 number_of_TMA_points =5 )
       ... print(static_cor)   
    """
    
    def __init__(self, data_fn=None , freq_array=None, res_array=None, phase_array =None , 
                 **kwargs): 
        
        self._logging =csamtpylog.get_csamtpy_logger(self.__class__.__name__)
        
        self.data_fn =data_fn 
        self._freq_array =freq_array 
        self._res_array =res_array
        
        self._reference_frequency =None 
        self._phase_array =None
        self._norm_freq=None 
        
        self._rj =None 
        self.rho_static=None 
        
        for keys in list(kwargs.keys()): 
            setattr(self, keys, kwargs[keys])
        
        
    @property 
    def frequency (self): 
        return self._freq_array
    @frequency.setter 
    def frequency(self, freq): 
        if freq.dtype not in ['float', 'int']:
            try : freq = np.array([float(ff) for ff in freq])
            except : raise CSex.pyCSAMTError_frequency('Frequency must be float number!')
        else : self._freq_array=freq
        
 
    @property 
    def app_rho (self): 
        return self._res_array 
    @app_rho.setter 
    def app_rho (self, app_res): 
        if app_res.dtype not in ['float', 'int']: 
            try : app_res =np.array([float(res) for res in app_res])
            except : raise CSex.pyCSAMTError_rho('Apparent resistivities values must be float number.!')
        self._res_array=app_res
            
    @property 
    def referencefreq(self):
        return self._reference_frequency 
    @referencefreq.setter
    def referencefreq(self, reffreq):
        try :reffreq =float(reffreq)
        except:raise CSex.pyCSAMTError_frequency('Reference frequency must be a float or int number.')
        if reffreq not in self.frequency:

            self._reference_frequency =Zcc.find_reference_frequency(freq_array=self.frequency, reffreq_value=reffreq ,
                                                                    sharp=True, etching=True)
        else : self._reference_frequency= reffreq

    @property 
    def phase (self): 
        return self._phase_array 
    @phase.setter 
    def phase(self, phz): 
        try : 
            phz=np.array([float(pzz) for pzz in phz])
        except:raise CSex.pyCSAMTError_Phase('Phase input must be on float number and radians.')
            # for phi in phz : 
            #     if not -np.pi < float(phi) < np.pi:raise CSex.pyCSAMTError_Phase('Phase value must be in range '\
            #                                             'of < [-pi; pi]> . Please check the phase array provide!. ')
        self._phase_array=phz
        
    
    
    
    def TMA (self, data_fn=None , freq_array=None, res_array=None, phase_array=None, stnVSrho_loc=None,   
                            reference_freq=None,number_of_TMA_points =5 ): 
        """ 
        Trimmed-moving-average filter to estimate average apparent resistivities at a
        single static-correction-reference frequency. User can compute TMA by inputing only the data file . 
        the program vill find automaticall other parameters . If not may provide all the parameters excepth the 
        data file . 
        
        Parameters
        -----------
            * data_fn : str 
                    path to avg file or edi file .
                
            * freq_array : array_like (ndarray,1) 
                    frequency array of at normalization frequency (reference value)
                    of all stations. station j to n .( units =  Hz )
                
            * res_array :     dict of array_like (ndarra,1) 
                    dict of array of app.resistivity at reffreq. from station j to n.
                    
            * phase_array : dict of array_lie(ndarray,1), dict of array of phase at reffreq.
                     from station j to n. (unit=rad)value of frequency with clean data . (unit=Hz)
                 
            * stnVSrho_loc : dict 
                    set of dictionnary of all app.resistivity data from station j to n . (optional)
                
            * num_of_TMA_point  :int 
                    window to apply filter .
            
        Returns
        -------
            dict 
               rho_corrected , value corrected with TMA filter  from station j to n. 
        
  
        1.  corrected data from [AVG]
        
        :Example:
            
            >>> from csamtpy.ff.processing.corr import Shifting
            >>> path =  os.path.join(os.environ["pyCSAMT"], 
            ...         'csamtpy','data', LCS01.AVG)
            ... static_cor =shifting().TMA (data_fn=path, 
            ...                            reference_freq=1024.,number_of_TMA_points =5 )
                  
        2. corrected from edifiles [EDI] 
        
        :Example:
            
            >>> from csamtpy.ff.core.cs import CSAMT
            >>> from csamtpy.ff.processing.corr import Shifting
            >>> edipath = r'C:/Users\Administrator\Desktop\test\edirewrite'
            >>> csamt_obj =CSAMT(edipath =edipath)
            >>> static_cor =shifting().TMA( reference_freq =256. ,
            ...                               freq_array = csamt_obj.freq ,
            ...                               res_array = csamt_obj.resistivity , 
            ...                            phase_array =csamt_obj.phase ,
            ...                            number_of_TMA_points=5)
            ... print(static_cor)

        """
        flag= 0
        self._logging.info ('Computing Trimming Moving Average of apparent resistivities.!')
        
        if data_fn is not None : 
            self.data_fn =data_fn 
            if  (self.data_fn.endswith('avg') or self.data_fn.endswith('avg'.upper())) == True  : 
                csamt_obj =CSAMTavg.Avg(self.data_fn)
                flag= 1 
   
            elif (self.data_fn.endswith('edi') or self.data_fn.endswith('edi'.upper())) == True:
                flag==2 
                pass 
            else :raise CSex.pyCSAMTError_file_handling('The input file doesnt not match any AVG or EDI file ! '\
                                                        'Please  check your file.')   
        else : 
            if freq_array is not None :self.frequency= freq_array 
            if (self.frequency is None)  or (res_array is None ) or \
                (phase_array is None)  : raise CSex.pyCSAMTError_inputarguments('None values can note be computed . '\
                                                                                      ' Please check your value.')

        if flag==1:
            res_app_obj=csamt_obj.Data_section.Resistivity.loc
            phase_obj =csamt_obj.Data_section.Phase.loc
            self.frequency=csamt_obj.Data_section.Frequency.value
            
        elif flag !=2 or flag != 1 :
            
            res_app_obj=res_array
            phase_obj =phase_array
            self.frequency=freq_array
            
        #---> set reference frequency . if not will detect automatically   as higher frequency with clean data 
        if reference_freq is not None :self.referencefreq=reference_freq
        elif reference_freq is None : 
            if flag ==1 : 
                self.referencefreq= Zcc.perforce_reference_freq(dataset=CSAMTavg.Avg().Data_section._data_array, 
                                                                  frequency_array=self.frequency) 
            else : self.referencefreq= self.frequency.max() # interpolate to highestv frequency value
        
        self.app_rho =Zcc.get_data_from_reference_frequency(array_loc=res_app_obj,
                                                          freq_array=self.frequency,
                                                          reffreq_value=self.referencefreq)
        self.phase=Zcc.get_data_from_reference_frequency(array_loc=phase_obj,
                                                          freq_array=self.frequency,
                                                          reffreq_value=self.referencefreq)
        
        #Zonge Avg file usaually put phase en mrad , we need to convert on rad by dividing by 1000.
        if flag==1 : self.phase =self.phase/1e3
        
        self.stnVSrho_loc=res_app_obj                                       # make a copy of dict_loc of apparent resistivity 
            

            

        slopej =np.arctan((self.phase/(np.pi/4)-1))*(np.pi/2)**-1           # compute the slope 
        
        
        rho_app_jplus =np.log10(self.app_rho)+ np.log10(np.sqrt(2))*slopej    #extrapolate up in frequency 
        
                                # collect a group of five log(rj+), i.e. for station index j, i = j-2 to j+2.
                                # Discard the lowest and highest valued log(rj+) from the group of five and average the remaining
                                # three => avg_log
        log_app_rj_tma = Zcc.compute_TMA (data_array=rho_app_jplus,
                                          number_of_TMApoints=number_of_TMA_points)
                                # The target static-correction apparent resistivity for station j
        rho_static_targetj =self.app_rho * (np.power(10,log_app_rj_tma))/np.power(10,rho_app_jplus)
        
                                # compute the rstatic factor rj : 
        self._rj =self.app_rho*(rho_static_targetj)**-1
                                # shiffting all value to corrected data : we make a copy.
                                # we assume that user can use the argument values already defined. no need to compute again
        
        stnNames =sorted(self.stnVSrho_loc.keys())
        self.TMA={key:value for key, value in zip (stnNames,self._rj.tolist())}
        for jj, stn in enumerate(stnNames) :
            for rhokeys, rhovalues in sorted(self.stnVSrho_loc.items()):
                if rhokeys ==stn : 
                    TMA_values=np.apply_along_axis(lambda rhoS: rhoS /self.TMA[stn],0, rhovalues)
                    self.TMA[stn]=TMA_values
        
        return self.TMA
    
    def FLMA (self):
        """
        ***Future plan ****
        """
        pass 
    def AMA (self):
        """
        ** future plan ****
        """
        pass 



def interp_to_reference_freq(freq_array, rho_array,  reference_freq, plot=False): 
    """
    Interpolate frequencies to the reference frequencies.
    
    :param freq_array:  frequency array
    :type freq_array: array_like
    
    :param reference_freq: frequency at clean data 
    :type reference_freq: float 
    """
    #find number of frequency
    freqObj , _reffreq =freq_array , reference_freq 
    rhoObj=rho_array
    if freqObj.dtype not in ['float', 'int'] :  
        raise TypeError('Frequency values must be float number.')
    if rhoObj.dtype not in ['float', 'int'] :
       raise TypeError('Apparent Resistivities values must be on float number')
      
    freqObj=np.array([float(kk) for kk in freqObj]) # do it in the case dtype value are integer.
    rhoObj=np.array([float(kk) for kk in rhoObj])
    
    if reference_freq not in freqObj : 
        raise CSex.pyCSAMTError_frequency('Reference frequency selected as input '\
                                          'argument is not Found on the frequency array.'\
                                              ' Please select the right frequency ')
    
    def cut_off_array_to_interp(freq_array, reference_freq, kind_interp='linear'):

        """
        Seek the array for interpolation with reference frequency.
        Return size of array to interpolate and the array 
        
        :param freq_array:  frequency array
        :type freq_array: array_like
        
        :param reference_freq: frequency at clean data 
        :type reference_freq: float 
        
        :param kind_interp: kind of interpolation 
                            Default is `linear`
        :type kind_interp: str 
        """
        for ss , value in enumerate(freq_array): 
            if reference_freq == freq_array[-1]:
                array_to_interp  = freq_array
                break
            if value == reference_freq:
                array_to_interp  =freq_array[:ss+1]
                break 
        return array_to_interp.size,array_to_interp  
    
    x_num_of_freq , interparray= cut_off_array_to_interp(freq_array=freqObj, reference_freq=_reffreq)
    
    func_interp =spi.interp1d(interparray, rhoObj[:x_num_of_freq], kind='linear')

    new_rhoObj= func_interp(freqObj[:x_num_of_freq] )
    if plot: 
        fig, ax =plt.subplots()
        ax.loglog(freqObj[:x_num_of_freq], new_rhoObj, c='r', lw=5)
        ax.loglog(freq_array,rhoObj )
        plt.show()
        
    return new_rhoObj


# if __name__=='__main__':
    
#     from csamtpy.pyCS.core.cs import CSAMT
    
#     edipath = r'C:\Users\Administrator\Desktop\test\edirewrite'
#     csamt_obj =CSAMT(edipath =edipath)
#     # freq_array=None, res_array=None, phase_array=None, stnVSrho_loc=None
#     static_cor =shifting().TMA( reference_freq =256. , freq_array = csamt_obj.freq , res_array = csamt_obj.resistivity , 
#                                 phase_array =csamt_obj.phase , number_of_TMA_points=5)
#     print(static_cor)
    
