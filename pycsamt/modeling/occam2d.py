# -*- coding: utf-8 -*-
#       Created on Fri Jan 22 20:31:14 2021
"""
.. _module-Occam2D :: `pycsamt.modeling.occam2d`
    :synopsis:: The Occam 2D MT inversion code (v3.0) presented here is an
                implementation of the general Occam procedure of  Constable, 
                et al. (1987)extended to  2D by deGroot-Hedlin and Constable (1990).
                The 2D MT forward  calculations are carried out with codeprovided
                 by Wannamaker, et al (1987)  using reciprocity to calculate the
                Jacobian (de Lugao and Wannamaker, 1996).
                ...

.. seealso::
        - http://marineemlab.ucsd.edu/Projects/Occam/sharp/index.html 
        - https://marineemlab.ucsd.edu/~kkey/Software.php
        - https://github.com/MTgeophysics/mtpy.git>`
        ...

"""
import os
import warnings
import datetime
import numpy as np 

import pycsamt.utils.exceptions as CSex
from pycsamt.modeling.__init__ import  SUCCESS_IMPORT_MTPY
from pycsamt.utils import func_utils as func
from pycsamt.utils import plot_utils as punc
from pycsamt.utils import  plotdecorator  as mdeco
from pycsamt.utils._p import _sensitive as SB 
from pycsamt.utils._csamtpylog import csamtpylog 

_logger =csamtpylog.get_csamtpy_logger(__name__)

if  SUCCESS_IMPORT_MTPY :

    _logger.info('`MTpy` successfull loaded!')
    try : 
        from mtpy.modeling import occam2d as MToccam2d
        _logger.info('`occam2d` module sucessfully imported  from `MTpy`.')
        
    except :
        _logger.info('loading  `occam2d` module from '
                     '`MTpy`packages  failed')
        warnings.warn('Loading occam2d module from "MTpy" '
                      'package failed! Please try again')
        
        SUCCESS_IMPORT_MTPY  =False 
    
if not SUCCESS_IMPORT_MTPY:
    _logger.info('Unable to import `MTpy`packages. Loading failed!')
    warnings.warn('Try to import `MTpy` module failed! You can get "MTPY"'
                  ' from <https://github.com/MTgeophysics/mtpy.git>`.')
     

class occamLog (object) : 
    """
    Class to deal with occcam 2d logfile . output File  after inverted data.
    most likely called `logFile.logfile`.
    
    Arguments
    ---------
        **fn** : str 
            full path to occam logfile 
            
    ==================  ==============  =====================================
    Attributes          Type            Description 
    ==================  ==============  =====================================
    rms                 float           RootMeansquared computation Error 
    iteration           float           number of iteration 
    fminocc             float           minimum tolerance 
    stepsize            float           cutoff evaluation function 
    roughness           float           deGrootHeldin roughness parameters.
    ==================  ==============  =====================================
    
    """
    occam_params =['rms', 'iteration', 'fminocc', 'stepsize',
                   'roughness', 'and_is']
    
    def __init__(self, fn =None, **kwargs):
        self._logging=csamtpylog.get_csamtpy_logger(self.__class__.__name__)
        self._fn =fn 
        for key in self.occam_params : self.__setattr__(key, [])
        
        if self.fn is not None : 
            
            self.read_occam2d_logfile()
            
            
    @property 
    def fn (self):
        return self._fn 
    
    @fn.setter 
    def fn (self, logpath):
        """
        assert  the occam logfile. 
        
        """
        self._logging.info ('Asserting <%s> occam2d logfile.'% logpath)
        
        if SB.which_file(filename=logpath) =='logfile':
            self._fn = logpath
        else :
            warnings.warn('No occam 2d logfile detected. Please check your path.')
            raise CSex.pyCSAMTError_occam2d(
                'Doesnt not recognize the Logfile <{0}>'\
                ' provided.'.format(os.path.basename(logpath)))
            
    def read_occam2d_logfile(self, fn =None ) : 
        """
        Read occam file and populate attributes. 
        
        :param fn:  full path to occam2d log file 
        :type fn: str

        :Example: 
            
            >>> from pycsamt.modeling import occam2d 
            >>> path =os.path.join(os.environ ['pyCSAMT'], 'pycsamt', 'data',
            ...                       'occam2D', 'logFile.logfile')
            >>> occamlog_obj =occam2d.Log(fn = path)
            >>> print(occamlog_obj.rms.shape)
            >>> print(occamlog_obj.roughness.shape)
        """
        
        if fn is not None : self.fn = fn 
        self._logging.info (
            'Reading and setting OccamLog  attributes'
            ' from logfiles <%s>.'% self.fn )
        if self.fn is None : 
            raise CSex.pyCSAMTError_occam2d(
                'None occam2dlog file detected. Please check your right path .')

        with open(self.fn , 'r', encoding ='utf8') as focc: 
            occamlog_lines =focc.readlines()
        for oc ,  occam  in enumerate(occamlog_lines) : 
            
            if occam.find(
                    'r.m.s'.upper())>=0 or occam.find('starting'.upper())>=0 : 
                self.rms.append(float(occam.strip().split('=')[1]))
        # for safety , we collect the iteration number from logile.
            if occam.find('iteration'.upper())>= 0 :
                 for iocc in occam.strip().split(): 
                     try : int(iocc)
                     except:  pass 
                     else : self.iteration.append(int(iocc))
            if occam.find('fminocc')>=0 or occam.find(
                    'minimum tol from'.upper())>=0 : 
                self.fminocc.append(float(occam.strip().split('=')[1]))
            if occam.find('and is'.upper()) >= 0 : 
                self.and_is.append(float(occam.strip().split('=')[1]))
            if occam.find('stepsize is'.upper()) >=0 :
                self.stepsize.append(float(occam.strip().split('=')[1]))
            if occam.find('roughness is'.upper())>=0 :
                self.roughness.append(float(occam.strip().split('=')[1]))
                
        #---< resetting attributes 
        for key in self.occam_params : 
            if hasattr(self,key) :
                value = getattr(self,key)
                self.__setattr__(key, np.array(value))
                
class Data(object): 
    """
    Read and write Occam data
    
    Arguments
    ---------
        **fn** : str 
            full path to occam Data file 
        
     .. note:: Replace the marker "**" in followig attributes  by "occam_data" 
     
    ===================  ================  ====================================
    Attributes            Type              Description
    ===================  ================  ====================================
    **format              str               occam sofware version name 
    **title               str               title given the building files 
    **sites               array_like        array of stations names 
    **offsets             array_like        array of site locations along
                                            the 2D line in meters  
    **frequencies         array_like        array of frequencies in decreasing 
    **data_nblocks        int               order number of data blocks for  
                                            building models 
    **data                (ndarray,4)       array of {site, freq, type, datum,
                                            error} . datacollected                  
    ===================  ================  ====================================
    
    Reads data files and find list of explanations of model model mode .  

    ===================== =====================================================
    Model Modes           Description                     
    ===================== =====================================================
    1 or log_all          Log resistivity of TE and TM plus Tipper
    2 or log_te_tip       Log resistivity of TE plus Tipper
    3 or log_tm_tip       Log resistivity of TM plus Tipper
    4 or log_te_tm        Log resistivity of TE and TM
    5 or log_te           Log resistivity of TE
    6 or log_tm           Log resistivity of TM
    7 or all              TE, TM and Tipper
    8 or te_tip           TE plus Tipper
    9 or tm_tip           TM plus Tipper
    10 or te_tm           TE and TM mode
    11 or te              TE mode
    12 or tm              TM mode
    13 or tip             Only Tipper
    ===================== =====================================================
    
    """
    occam_data_type= {'te_log10':1, 
                 'te_phase':2, 
                 'real_tip':3, 
                 'imag_tip':4, 
                 'tm_log10':5,
                 'tm_phase':6,
                 'tm_tip':7,
                 'tm_tip':8,
                 'te':9,
                 'tm':10,
                 'real_zxx':11,
                 'imag_zxx':12,
                 'real_zxy':13, 
                 'imag_zyx':14,
                 'real_zyx':15,
                 'imag_zyx':16,
                 'real_zyy':17,
                 'imag_zyy':18 }
    occam_dataType , occam_modes = {
                    1 :'TE apparent resistivity (log10)',
                    2 : 'TE phase',
                    3 : 'real(tipper) Hz/Hy',
                    4 : 'imag(tipper) Hz/Hy',
                    5 : 'TM apparent resistivity (log10)',
                    6 : 'TM phase',
                    7 : 'reserved for TM Ez/Ey tipper',
                    8 : 'reserved for TM Ez/Ey tipper',
                    9 : 'TE apparent resistivity (linear)',
                    10 : 'TM apparent resistivity (linear)',
                    12 : 'imag(Zxx) ! not used until we incorporate site-rotation angle',
                    11 : 'real(Zxx) ! not used until we incorporate site-rotation angle',
                    13 : 'real(Zxy) ! TE',
                    14 : 'imag(Zxy) ! TE',
                    15 : 'real(Zyx) ! TM',
                    16 : 'imag(Zyx) ! TM',
                    17 : 'imag(Zyy) ! not used until we incorporate site-rotation angle',
                    18 : 'imag(Zyy) ! not used until we incorporate site-rotation angle'},{
                                        'log_all': [1, 2, 3, 4, 5, 6],
                                        'log_te_tip': [1, 2, 3, 4],
                                        'log_tm_tip': [5, 6, 3, 4],
                                        'log_te_tm': [1, 2, 5, 6],
                                        'log_te': [1, 2],
                                        'log_tm': [5, 6],
                                        'all': [9, 2, 3, 4, 10, 6],
                                        'te_tip': [9, 2, 3, 4],
                                        'tm_tip': [10, 6, 3, 4],
                                        'te_tm': [9, 2, 10, 6],
                                        'te': [9, 2],
                                        'tm': [10, 6],
                                        'tip': [3, 4],
                                        '1': [1, 2, 3, 4, 5, 6],
                                        '2': [1, 2, 3, 4],
                                        '3': [5, 6, 3, 4],
                                        '4': [1, 2, 5, 6],
                                        '5': [1, 2],
                                        '6': [5, 6],
                                        '7': [9, 2, 3, 4, 10, 6],
                                        '8': [9, 2, 3, 4],
                                        '9': [10, 6, 3, 4],
                                        '10': [9, 2, 10, 6],
                                        '11': [9, 2],
                                        '12': [10, 6],
                                        '13': [3, 4]}


    data_blocks_header=['site', 'freq', 'type', 'datum', 'error']
    
    def __init__(self, data_fn =None , **kwargs): 
        self._logging =csamtpylog.get_csamtpy_logger(self.__class__.__name__)
        self.fn = data_fn
        self.data_format = kwargs.pop('format', 'OCCAM2MTDATA_1.0')
        self.data_title =kwargs.pop('title', 'Forward model response from Occam2D')
        self.data_sites =kwargs.pop('sites', None)
        self.data_offsets =kwargs.pop('offsets', None)
        self.data_frequencies =kwargs.pop('frequencies', None)
        self.data_nblocks = kwargs.pop('data_blocks', None)
        self.verbose = kwargs.pop('verbose', 0)
        
        self.occam_data =None

        for key in list(kwargs.keys()): 
            self.__setattr__(key, kwargs[key])
        
        if self.fn is not None : 
            self.read_occam2d_datafile ()
            
            
    def read_occam2d_datafile (self, data_fn =None ): 
        """
        read_occam_data file. and populates attributes 

        :param fn: full path to occam data file .
        :type fn: str
            
        :Example:
            
            >>> path =os.path.join(os.environ ['pyCSAMT'], 
            ...                       'pycsamt', 'data', 'occam2D', 
            ...                       'OccamDataFile.dat')
            >>> occamdata_obj =Data(fn = path)
            >>> print(occamdata_obj.occam_data_format)
            >>> print(occamdata_obj.occam_data_blocks)
            >>> print(occamdata_obj.occam_data_title)
            >>> print(occamdata_obj.occam_data_sites)
            >>> print(occamdata_obj.occam_data_frequencies)
            >>> print(occamdata_obj.occam_data_offsets)
            >>> occamDATA = occamdata_obj.occam_data
        """
        if data_fn is not None : 
            self.fn =data_fn
            
        self._logging.info ('Read Occam2D data file :<%s>' % self.fn)
        
        if SB.which_file(filename=self.fn)=='occamdat':
            with open(self.fn , 'r', encoding='utf8') as focd : 
                occam_data_lines =focd.readlines()
                
        # get data from occam files 
        indextemp =[]
        for ii, dataparams in enumerate(occam_data_lines) : 
            if dataparams.lower().find('format')>= 0 :
                self.data_format= dataparams.strip().split(':')[-1]
            if dataparams.lower().find('title') >=0 : 
                self.data_title =dataparams.strip().split(':')[-1]
            if dataparams.lower().find('site') >=0 :
                indextemp.append(ii)
            if dataparams.lower().find('offset') >=0 :
                indextemp.append(ii)
            if dataparams.lower().find('frequencies')>=0 :
                indextemp.append(ii)
            if dataparams.lower().find('block')>=0 :
                self.data_nblocks= float(dataparams.strip().split(':')[-1])
    
            if dataparams.lower().find('type') >=0 and dataparams.lower(
                    ).find('error')>=0 : 
                indextemp.append(ii)
                
                break
        #--> now get each specials arrays and eliminate (the first tiles )
        self.data_sites=[site.strip() 
                         for site in occam_data_lines[
                                 indextemp[0]+1:indextemp[1]]]
        self.data_offsets = [float(off.strip()) 
                             for off in occam_data_lines[
                                     indextemp[1]+1 : indextemp[2]]]
        self.data_frequencies =[float(freq.strip())
                                for freq in occam_data_lines[
                                        indextemp[2]+1:indextemp[3]-1]]
        
        #cchek if frequencies are sorted in decrease order 
        #: highest to lower frequencies 
        self.data_flip_freq=False 
        if self.data_frequencies [0] <self.data_frequencies[-1]: 
            self.data_frequencies= self.data_frequencies[::-1]
            self.data_flip_freq =True 
            
        temp=[]
        for data in occam_data_lines[indextemp[-1]+1:]:
            temp.append(np.array(data.strip().split()))
        
        self.data = func.concat_array_from_list(list_of_array=temp)
        
        # set datablocks attributes  
        self.parse_block_data()
        
       #XXXFIXME
        # for ikey, key in enumerate(self.data_titles ) : 
        occam_data_dict ={key :self.data[:, ikey] 
                          for ikey, key in enumerate(
                                  self.data_blocks_header)}
        
        self.__setattr__('occam_data_blocks', DataBlock(**occam_data_dict))
        #straighten out offset  and get the step 
        # graduate to min to max for consistency
        if self.data_offsets[0]> self.data_offsets[-1] :  
            self.data_offsets =self.data_offsets[::-1]
            
        step=  func.round_dipole_length( (
            self.data_offsets[-1] - self.data_offsets[0])/ (len(
                self.data_offsets) -1))
        
        #regraduate offset
        self.data_offsets= np.arange(0, step* len(self.data_offsets), step)
        self.data_frequencies = np.array(self.data_frequencies )
        
        if self.verbose >0: 
            
            print('{0:-^77}'.format('Occam 2D Data infos'))
            print('** {0:<27} {1} {2}'.format('Sites num.', 
                                              '=', len(self.data_offsets)))
            print('** {0:<27} {1} {2}'.format('Frequencies num.',
                                              '=', len(self.data_frequencies)))
            print('** {0:<27} {1} {2}'.format('Highest frequency (Hz)',
                                              '=', self.data_frequencies.max()))
            print('** {0:<27} {1} {2}'.format('Lowest frequency (Hz)', '=',
                                              self.data_frequencies.min()))
            
            print('** {0:<27} {1} {2}'.format('Minimum offset (m)', '=',
                                              self.data_offsets.min()))
            print('** {0:<27} {1} {2}'.format('Maximum offset (m)', '=',
                                              self.data_offsets.max()))
            print('-'*77)
        
        
    def parse_block_data(self, datablocks=None):
        """
        From block data, retreive the value of each blocks 
        each blocks is ndarray(nfreq, nstations)
        Attributes are : 
            - `self.appRho` for apparent resistivity blocks
            - `self.phase` for phase  values 
            - `self.error_appRho` for error in resistivity  (block )
            - `self.error_phase` for error in phase 9blocks
        
        """    
        # got frequency array  and site array which uniques values 
        if datablocks is not None: 
            self.data = datablocks 

        # self.data.astype(np.float)  # for consistency 
        self.data = np.asarray(self.data,dtype = np.float64, order= 'K')
        
        sites, *_= np.unique(self.data[:, 0], 
                             return_counts =True)
        freq , *_=np.unique(self.data[:, 1], 
                            return_counts= True)
        datum , *_= np.unique(self.data[:, 2],
                              return_counts= True)
        sites.astype('int32')
        freq.astype('int32')
    

        # get mode and create model array 
        k_=None 
        m_= [int(s) for s in datum]
    
        self.occam_dtype_value = m_
        for key, values in self.occam_modes.items(): 
            if values ==m_: 
                k_= key 
                self.occam_dtype_ =k_
                self.occam_dtype_value = values
                break 
            
        cf =False 
        if self.occam_dtype_ =='log_te_tm': 
            self.te_appRho =np.zeros((len(freq), len(sites)))
            self.te_phase = np.zeros((len(freq), len(sites)))
            self.tm_appRho =np.zeros((len(freq), len(sites)))
            self.tm_phase = np.zeros((len(freq), len(sites)))
            self.te_error_appRho = np.zeros((len(freq), len(sites)))
            self.te_error_phase = np.zeros((len(freq), len(sites))) 
            self.tm_error_appRho = np.zeros((len(freq), len(sites)))
            self.tm_error_phase = np.zeros((len(freq), len(sites))) 
            cf =True
        elif self.occam_dtype_=='log_te': 
            self.te_appRho =np.zeros((len(freq), len(sites)))
            self.te_phase = np.zeros((len(freq), len(sites)))
            self.te_error_appRho = np.zeros((len(freq), len(sites)))
            self.te_error_phase = np.zeros((len(freq), len(sites))) 
            
        elif self.occam_dtype_=='log_tm': 
            self.tm_appRho =np.zeros((len(freq), len(sites)))
            self.tm_phase = np.zeros((len(freq), len(sites)))
            self.tm_error_appRho = np.zeros((len(freq), len(sites)))
            self.tm_error_phase = np.zeros((len(freq), len(sites))) 
       
        # create for each mode id datablocks 
  
        for kk, rowlines in enumerate(self.data): 
            index_site = np.where(sites ==rowlines[0])[0]
            index_freq = np.where(freq ==rowlines[1])[0]
            if int(rowlines[2]) ==1: 
                self.te_appRho[int(index_freq),
                           int(index_site)] = rowlines[3]
                self.te_error_appRho[int(index_freq),
                           int(index_site)] = rowlines[4]
                
            if int(rowlines[2]) ==2:
                self.te_phase[int(index_freq),
                           int(index_site)] = rowlines[3]
                self.te_error_phase[int(index_freq),
                           int(index_site)] = rowlines[4]
            if int(rowlines[2]) ==5: 
                self.tm_appRho[int(index_freq),
                           int(index_site)] = rowlines[3]
                self.tm_error_appRho[int(index_freq),
                           int(index_site)] = rowlines[4]
            if int(rowlines[2]) ==6: 
                self.tm_phase[int(index_freq),
                           int(index_site)] = rowlines[3]
                self.tm_error_phase[int(index_freq),
                           int(index_site)] = rowlines[4]
                
        # flip data so of frequency in decrease order   
        if self.data_flip_freq is True : 
            
            if cf is True or self.occam_dtype_=='log_te': 
                self.te_appRho=self.te_appRho[::-1]
                self.te_phase = self.te_phase[::-1]
                self.te_error_appRho= self.te_error_appRho[::-1]
                self.te_error_phase= self.te_error_phase[::-1]
            if cf is True or self.occam_dtype_=='log_tm': 
                self.tm_appRho=self.tm_appRho[::-1]
                self.tm_phase = self.tm_phase[::-1]
                self.tm_error_appRho=self.tm_error_appRho[::-1]
                self.tm_error_phase=self.tm_error_phase[::-1]
                    
        
class DataBlock (object): 
    """
    Read OccamDataBlock aand set corresponding attributes 
     
    ======================  ===================================================
    Attributes              Description (Restriction)
    ======================  ===================================================
    site                    number of the site from the site list that this data 
                            belongs to.
    freq                    number of the frequency from the frequency list.
    type                    type of data
    datum                   data value
    error                   size of the error for this measurement  
                            For log10 resistivity this value can look a little 
                            strange. It is derived from the calculation
                            d(log_10(x))/dx = 1/[x ln(10)]. So for 10% error, 
                            dx = 0.10x thus d(log_10(x)) = 0.10x / x ln(10) 
                            =.1/ln(10) = 0.0434
    ======================  ===================================================
    
    :Example:
        
        >>> form pycsamt.modeling import DataBlock 
        >>> np.random.seed(1983)
        >>> dblocks_dict = {'site' : ['S{0:02}'.format(ii)  for ii in range (10)], 
        ...                'freq': np.linspace(1, 8912, 10), 
        ...                'type' :(5,6), 
        ...                'data' : np.power(10, np.random.randn(20*2)), 
        ...                'error' : np.random.randn(20)
        ...                }
        >>> dblock_obj = DataBlock(**dblocks_dict)
        >>> print(dblock_obj.freq )
        >>> print(dblock_obj.error)
        >>> print(dblock_obj.site)
    """
    
    # blocks_params = ['site', 'feq', 'type', 'datum', 'error']

    def __init__(self, **kwargs):
        self._logging =csamtpylog.get_csamtpy_logger(self.__class__.__name__)

        for keyb in list(kwargs.keys()): 
            setattr(self, keyb, kwargs[keyb])
        
        #cehck main attributes []
        if hasattr(self, 'type'): # find the occam2d mode provided 
            if self.type is not None : 
                type_value =np.unique(self.type)
                for key, items in Data.occam_modes.items(): 
                    if items ==[int(ii) for ii in list(type_value)]: 
                        self.__setattr__('data_mode', key)
        
        for ikey in ['site', 'freq', 'error']: # 
            if hasattr(self, ikey):
                if getattr(self, ikey) is not None : 
                    if ikey =='error' : 
                        setattr(self, ikey, np.array(
                            [float(err)for err in getattr(self, ikey)]))
                    else : setattr(self, ikey, 
                                   [int(sf) for sf in getattr(self, ikey)])
  
        # check if all attribute are set from data  
        self._logging.info ('Ckeck dataBlocks and share corresponding'
                            ' data params.')
        
        _dbk= np.array([1 for attr in Data.data_blocks_header 
                        if getattr(self, attr) is not None ])
        
        if np.all(_dbk==1) and _dbk.size ==5 : 
            self._set_dbk_properly() 

                
    def _set_dbk_properly (self) : 
        """
        method to set property data 
        Ckeck dataBlocks and share corresponding data params separately.
        
        for instance:
             site (1=site number 1=S00), freq(1: freq numb1),
                   dataType(5: logTM, 6:phaseTM) 
        like:
            dict_mode ={ '6': }
        """
        def fill_mode_data (d, mode, keep_otype =False):
            """ Build each occam data from logtem res/phase values. 
            :param d: ndarray data of np.ndarray(nsites, 5)
            :param mode: Occam mode from :attar:`Data.occam_modes`
            :param keep_otype: Keep occam type inside the data. if ``True`` 
                type column should be remove in the data. 
            :return: dict of occam_mode values as dictkeys and their 
                corresponding  data. 
            """
            mvalues = Data.occam_modes[mode]
            m = {str(vk):d [d[:, 0] ==int(vk)] for ii, vk in enumerate(mvalues)} 
            return m if keep_otype else {k: v[:, 1:] for k , v in m.items()}
            
        ar = func.concat_array_from_list(
            [self.type, self.site, self.freq, self.datum, self.error], 
            concat_axis=1 )
        
        self.dict_mode= fill_mode_data(ar, mode =self.data_mode)
        
        return self 
    
        
    @staticmethod
    def decode_each_site_data (data_blocks, data_type_index)  : 
        """
        Decode each site data and get differents values array 
        :param datablocks:  datafrom differents blocks
        :type datablocks: ndarray 
        
        :param data_type_index: specify the data type index 
        :type data_type_index: int 
        
        """
        block_type = data_blocks[ :,  data_type_index]
        get_occam_mode  = np.array([int(ff) for ff in np.unique(block_type)])
        

        for keymode, valuemode in Data.occam_modes.items() : 
            if list(get_occam_mode) ==valuemode : 
                occam_data_mode = keymode
                occam_value_mode = valuemode 
  
        if occam_data_mode =='log_tm' or occam_data_mode =='6': 
            temp, decode_dict =[], {}
            # get the first value for occma mmode 
            value_to_keep = Data.occam_modes['log_tm'][0]   
            for ii, respitems in enumerate(data_blocks):
                if int(respitems[int(data_type_index)]) == int(value_to_keep) : 
                    temp.append(respitems)
        decode_dict [occam_data_mode] = func.concat_array_from_list(
            list_of_array=temp) 
        
        if occam_data_mode =='log_te' or occam_data_mode =='5': 
            temp,decode_dict =[], {}
            # get the first value for occma mmode 
            value_to_keep = Data.occam_modes['log_te'][0]   
            for ii, respitems in enumerate(data_blocks):
                if int(respitems[int(data_type_index)]) == int(value_to_keep): 
                    temp.append(respitems)
            
            decode_dict [occam_data_mode] = func.concat_array_from_list(
                list_of_array=temp) 
    
            
        return decode_dict, occam_data_mode,  occam_value_mode
                   

class Model (object): 
    """
    Occam Model , Actually read model file
    
    Arguments
    ----------
        **model_fn** : str 
                full path to occam model file 
        
    ..note:: Replace the marker "**" in followig attributes  by "model_"
    
    ==================  =================  ====================================
    Attributes          Type                Description 
    ==================  =================  ====================================
    **format              str               occam sofware version name 
    **description         str               describe the model 
    **name                str               name of the models 
    **meshfile            str               The name of the mesh file.
    **mesh_type           str               the type of the mesh ,default is PW2D
    **static_file         str               name of static file if given . 
                                            Default  is 'none'.
    **prejudice_file      str               name of the prejudice file if given.
                                            default is 'none'.models 
    **binding_offset      float             param value for a model block can be 
                                            aggregation into  many mesh blocks
    **num_layers          float             number of layers of model blocks          
    ==================  =================  ====================================
    
    :Example:   
        
        >>> from pycsamt.modeling.occam2d import Model
        >>> data='OccamDataFile.dat'
        >>> mesh = 'Occam2DMesh'
        >>> model = 'Occam2DModel'
        >>> iter_='ITER17.iter'
        >>> path =os.path.join(os.environ ['pyCSAMT'], 'pycsamt', 'data', 'occam2D', mesh)
        ...                   #,'OccamDataFile.dat')
        >>> occam_model_obj =Model(mesh_fn=path, 
        ...                    iter_fn = os.path.join(os.path.dirname(path), iter_), 
        ...                    model_fn =os.path.join(os.path.dirname(path), model) )
        ... print(occam_model_obj.model_binding_offset)
        ... print(occam_model_obj.model_resistivity)
        ... print(occam_model_obj.model_values)
    """

    def __init__(self, model_fn =None ,iter_fn =None , mesh_fn =None,  **kwargs): 

        self._logging =csamtpylog.get_csamtpy_logger(self.__class__.__name__)
        
        self.model_fn = model_fn
        self.iter_fn =iter_fn 
        self.mesh_fn =mesh_fn 
        
        self.model_format=kwargs.pop('format', 'OCCAM2MTMOD_1.0')
        self.model_name =kwargs.pop('model_name','Descriptive text here' )
        self.model_description =kwargs.pop('description',
                                           'SIMPLE INVERSION|More descriptive text')
        self.model_mesh_file =kwargs.pop('mesh_file', 'Occam.mesh')
        self.model_mesh_type =kwargs.pop('mesh_type', 'PW2D')
        
        self.model_statics_file =kwargs.pop('prejudice_file', None) 
        self.model_prejudice_file =kwargs.pop('prejudice_file', None )
        self.model_binding_offset =kwargs.pop('binding_offset', None )
        self.model_num_layers =kwargs.pop('num_layers', None )
        self.model_exceptions = kwargs.pop('numbers_exceptions', None)
        self.verbose = kwargs.pop('verbose', 0)
   
        for key in list(kwargs.keys()):
            self.__setattr__(key, kwargs[key])
            
        if self.model_fn  is not None : 
            self.read_occam2d_modelfile ()
            
    def read_occam2d_modelfile(self, model_fn=None, mesh_fn=None , iter_fn=None):
        """
        read   and ascertain modelfile.
        
        :param model_fn: full path to OCCAM2D model file 
        :type model_fn: str 
        
        :param mesh_fn: full path to MESH model file 
        :type mesh_fn: str 
        
        :param iter_fn: full path to ITER model file 
        :type iter_fn: str 
        
        :Example:
            
            >>> path =os.path.join(os.environ ['pyCSAMT'], 'pycsamt', 
            ...                       'data', 'occam2D', 'Occam2DModel')
            >>> occammodel_obj =Model(model_fn = path)
            ... print(occammodel_obj.model_name)
            ... print(occammodel_obj.model_mesh_file)
            ... print(occammodel_obj.model_num_layers)
        """
        if model_fn is not None : self.model_fn= model_fn 
        if iter_fn is not None : self.iter_fn = iter_fn 
        if mesh_fn is not None : self.mesh_fn =mesh_fn 
        
        self._logging.info ('Read Occam2D MODEL file :<%s>' % self.model_fn )
        
        #initialise attrinute as empty list for row and colums appended 
        for key in ['model_columns', 'model_rows']: self.__setattr__(key, [])
            
        # make ascertainement for all files provided 
        for  ifile , fname in zip ([self.model_fn, self.iter_fn, self.mesh_fn], 
                                   ['model', 'iteration', 'mesh']): 
            
            mess='Error Input file ! No {0} file found.'\
                ' Please input {0} file.'.format(fname)
             
            if ifile is None : 
                self._logging.error(mess)
                warnings.warn(mess)
                raise CSex.pyCSAMTError_occam2d(mess)
                
                
        if self.iter_fn is None :
            raise CSex.py
        
        if SB.which_file(filename =self.model_fn )=='model':
            with open(self.model_fn , 'r', encoding='utf8')as fmod : 
                occam_model_datalines =fmod.readlines ()
    
        #initialise empty list for model columns and rows and populates attributes
        for ii, modelitem in enumerate(occam_model_datalines): 
            if modelitem.find(':')> 0: 
                modelitem =modelitem .strip().split(':')
                if modelitem[0].lower().find('exception') >=0 :
                    modelitem [0]='exception' 
                # keep the index of num_layer so to loop for fill model values 
                try : 
                    self.__setattr__('model_' + modelitem[0].replace(
                        ' ', '_').lower(), float(modelitem[1]))
                except : 
                    self.__setattr__('model_' + modelitem[0].replace(
                        ' ', '_').lower(), modelitem[1])
                    pass 
            else : 
                modelitem =modelitem.strip().split()
                if len(modelitem) ==2 : 
                    self.model_rows.append([int(ii) for ii in (modelitem)])
                elif len(modelitem)> 2 : 
                    self.model_columns.append([int(ii) for ii in modelitem])
        
        
        # get model value from iteration object 
        occam_iter_obj =Iter(iter_fn = self.iter_fn)
        self.model_values = occam_iter_obj.iter_data  # flatterned array
        # get iteration info 
        self.model_rms = round(occam_iter_obj.iter_misfit_value,2)
        self.model_roughness= round(occam_iter_obj.iter_roughness_value)

        # get xnodes and znodes from mesh object  # dont need to get grid because 
        occam_mesh_obj = Mesh(mesh_fn = self.mesh_fn )
        occam_mesh_x_nodes = occam_mesh_obj.mesh_x_nodes 
        occam_mesh_z_nodes =occam_mesh_obj.mesh_z_nodes 

        # initilise model array to put resistivity into agragated mesh 
        self.model_resistivity =np.zeros(
            (occam_mesh_z_nodes.size, occam_mesh_x_nodes.size))
      
        self.model_rows=np.array(self.model_rows) 
    
        self._logging.info ('Read Model Blocks and put model'
                            ' resistivity inot aggregated mesh blocks ')
        #---> make sure the size of model column = modelsize 

        if len(self.model_columns) != len(self.model_rows):
            mess =''.join([
                'Model rows and Models columns should be  the same length.',
                           'However lengths model rows =',
                           ' {0} and model colums = {1}.'.format(
                               len(self.model_columns), len(self.model_rows)), 
                           'Please check your <{0}> modelfile.'.format(
                               os.path.basename(self.model_fn) )
                           ])
                                                                          
            warnings.warn (mess)
            self._logging.error(mess)
            raise CSex.pyCSAMTError_occam2d(mess)
        
        #---> loop the model rows suppose the number of models layer 
        occm=0
        for mindex  in range(len(self.model_rows)): 
            #get the first inthe first in the pair of lines 
            #that describes a layer of model blocks.
            # sum to get the number of model layer at the end
            fp_sta = self.model_rows[:mindex, 0].sum()  
            #find next value for slicing 
            fp_end = fp_sta + self.model_rows[mindex][0] 
            # find the second pair of lines that describes
            #  a layer of model blocks. is the columns 
            sp_columns = np.array(self.model_columns[mindex])
            # then loop the collum and put the the number of 
            #mesh blocks that are aggregated into each model block
            for im in range(len(sp_columns)): 
                cp_sta = sp_columns [:im].sum()
                cp_end = cp_sta + sp_columns[im]    # stop add when 
                #Put model resistivities into aggregated mesh 
                self.model_resistivity[
                    fp_sta:fp_end, cp_sta:cp_end] = self.model_values[occm]
                occm +=1 
      
        #-------MAKE A PLOT PURPOSE ----------------------
        # making model z_offset
        self.z_grid , self.x_grid =occam_mesh_z_nodes.copy(
            ), occam_mesh_x_nodes.copy()
        
        # for dep, offs in zip( [self.z_grid, self.x_grid],
        #   [occam_mesh_z_nodes, occam_mesh_x_nodes])  : 
        #     dep = np.array([dep[:ii].sum() for ii in range(len(offs))])
        # self.z_grid = np.array([self.z_grid[:ii].sum()
        #   for ii in range(len(occam_mesh_z_nodes))])
        # self.x_grid =np.array([self.x_grid[:ii].sum() 
        #   for ii in range (len(occam_mesh_x_nodes))])
        
        # make some arrays for plotting the model
        self.model_station_offsets = np.array([occam_mesh_x_nodes[:ii + 1].sum()
                                for ii in range(len(occam_mesh_x_nodes))])
        self.model_depth_offsets = np.array([occam_mesh_z_nodes[:ii + 1].sum(
            )for ii in range(len(occam_mesh_z_nodes))])
                                
        # center the grid onto the station coordinates
        x0 = self.model_binding_offset - self.model_station_offsets[
            self.model_columns[0][0]]
        self.model_station_offsets += x0
        
        # rescal the model offset to to have the top elevation to zero  
        self.model_depth_offsets = self.model_depth_offsets -\
            self.model_depth_offsets[0]
        
        # self.model_plot_z = self.model_plot_z[::-1]  #--> flipping depth 
        if self.verbose > 0: 
            print('{0:-^37}'.format('Boundaries X (Horizontal nodes)'))
            print('**{0:<27} {1} {2}'.format(' Minimum offset (m)','=' ,
                                             self.model_station_offsets.min()))
            print('**{0:<27} {1} {2}'.format(' Maximum offset (m)','=',
                                             self.model_station_offsets.max()))
            print('{0:-^37}'.format('Boundaries Z (Vertical nodes)'))
            print('** {0:<20} {1} {2}'.format('Minimum depth (m)', '=',
                                              self.model_depth_offsets.min()))
            print('** {0:<20} {1} {2}'.format('Maximum depth (m)', '=',
                                              self.model_depth_offsets.max()))
            
            print('{0:-^77}'.format('Occam 2D Models params'))
            print('** {0:<27} {1} {2}'.format('Model layer num.', '=',
                                              len(self.model_rows)))
            print('** {0:<27} {1} {2}'.format('Model param count', '=',
                                              occam_iter_obj.iter_param_count))
    
            print('** {0:<27} {1} {2}'.format('Iteration num.', '=', 
                                              occam_iter_obj.iter_iteration))
            print('** {0:<27} {1} {2}'.format('Occam Misfit value', '=',
                                              occam_iter_obj.iter_misfit_value))
            print('** {0:<27} {1} {2}'.format('Occam Misfit reached', '=',
                                              occam_iter_obj.iter_misfit_reached))
            print('** {0:<27} {1} {2}'.format('Occam Misfit target', '=',
                                              occam_iter_obj.iter_target_misfit))
            print('** {0:<27} {1} {2}'.format('Occam Roughness params', '=',
                                              self.model_roughness))

        


class Startup(object): 
    """
    Occam startup file  Actually read startup file.
    for more detail:
    
    .. seealso:: http://marineemlab.ucsd.edu/Projects/Occam/sharp/index.html
        
    Arguments 
    ----------
        **startup_fn** : str 
                full path to occam startup file 
                
    .. note: to get "startup" attributes , replace the marker"**" by "startup_"

    ======================  ==========  =======================================
    Attributes              Type        Description 
           
    ======================  ==========  =======================================
    **format                str          occam sofware version name 
    **description           str         describe the startup file
    **model_file            str         name of the model file . 
    **data_file             str         The name of the data file.
    **iteration_to_run      float       number of iteration to run 
    **target_misfit         float       misfit target. default is 1.0 
    **roughness_type        int         type of roughness params:
    **diagonal_penalties    int         penalties horizontally & vertically
                                        [0|1].default is 0
    **stepsize_cut_count    int        The parameter limits the number of times
                                        Occam cuts the step size in a
                                        search for a better fitting model.
    **model_limits          tuple       (min, max):params imposes limits
                                        on the values that the model parameters 
                                        may take.min & max are also log10 res.
    **model_value_steps     flaot      stepsize : discretizes model parameter
                                        space into steps “stepsize” large.
                                        the model space is in log10 resistivity, 
                                        so thesesteps arein the exponent (e.g. 
                                        for stepsize=0.1, values allowed will be 
                                        100.1, 100.2,100.3, 100.4, etc…
    **debug_level           int         may be [0, 1, 2] : 
                                        0 minimizes the amount of “chatter” that 
                                        Occam2D displays while running.
                                        1 (the default), displays progress on a 
                                        frequency-by-frequency basis.
                                        2 displays more diagnostic info and 
                                        changes the way Occam2D searches model 
                                        space.
    **iteration             int         Occam2D iteration number represented 
                                        on the file.
    **param_count           int         assumes that the model parameter
    ======================  ==========  =======================================
    
    """

    def __init__(self, startup_fn =None, **kwargs): 
        
        self._logging =csamtpylog.get_csamtpy_logger(self.__class__.__name__)
        self.startup_fn =startup_fn 
        self.startup_format=kwargs.pop('format','OCCAMITER_FLEX' )
        self.startup_description=kwargs.pop('description', 'Some description')
        self.startup_model_file =kwargs.pop('model_file', 'occam.model')
        self.startup_data_file =kwargs.pop('data_file', 'datafilename.dat')
        self.startup_iterations_to_run=kwargs.pop('iterations_to_run', None)
        self.startup_target_misfit = kwargs.pop('target_misfit', 1.)
        self.startup_roughness_type =kwargs.pop('roughness_type', 1.)
        self.startup_diagonal_penalties =kwargs.pop('diagonal_penalties', 0)
        self.startup_stepsize_cut_count =kwargs.pop('stepsize_cut_count', 8)
        self.startup_model_limits=kwargs.pop('model_limits', (None, None))
        self.startup_value_steps = kwargs.pop('stepsize', None)
        self.startup_debug_level =kwargs.pop('debug_level', 1)
        self.startup_iteration =kwargs.pop('iteration', 0)
        
        self.startup_lagrange_value = 5. 
        self.startup_roughness_value =0.1E+11
        self.startup_misfit_value = 1000.
        self.startup_misfit_reached = 0. 
        self.startup_param_count =kwargs.pop('param_count', None)
        
        self.startup_date_time =datetime.datetime.now(
            ).strftime('%m-%d-%Y %H:%M:%S')
        
        for key in list(kwargs.keys()): 
            setattr(self, key, kwargs[key])
        if self.startup_fn is not None : 
            self.read_occam2d_startupfile()
            
    def read_occam2d_startupfile(self, startup_fn=None):
        """
        read occam2d_startupfile
        
        :param startup_fn: full path to startup_file 
        :type startup_fn: str 
        
        :Example:
            
            >>> from pycsamt.modeling.occam2d import Startup
            >>> path =os.path.join(os.environ ['pyCSAMT'],
            ...                       'pycsamt', 'data', 'occam2D', 'Startup')
            >>> startup_obj=Startup(startup_fn = path)
            ... print(startup_obj.occam_startup_data_file)
        """
        if startup_fn is not None : self.startup_fn = startup_fn
        self._logging.info (
            'Read Occam2D STARTUP file :<%s>' % self.startup_fn)
        
        if SB.which_file(filename=self.startup_fn)=='startup': 
            with open(self.startup_fn, 'r', encoding='utf8') as fsup : 
                occam_startup_datalines =fsup.readlines ()
                
        # buid list of all attributs then collapse to 
        #find the name in the data file and loop it 
        attr_list = [ attr for attr in  list(
            self.__dict__.keys()) if 'startup_' in attr ]
        # maybe usefull when writing startup files 
        self.twin_attr = [attr.replace(
            'startup_', '').replace('_', ' ') for attr in attr_list] 

        for ii, startupitems in enumerate(occam_startup_datalines):
            for ii, attr in enumerate(self.twin_attr) : 
                if startupitems.lower().find(attr) >=0 :
                    try : 
                        self.__setattr__(attr_list[ii], float(
                            startupitems.strip().split(':')[-1]) )
                    except : 
                        self.__setattr__(attr_list[ii], startupitems.strip(
                            ).split(':')[-1]) 
                        pass 

class Iter (Startup): 
    """
    Occam iteration file, inherets from startup obj , 
    in fact two object a similar the same 
    
    Arguments
    ----------
        **iter_fn** : str 
            full path to occam iteration file 
    
    .. note:: Most attributes are describe in Startup.__doc__
            Iter is a copy of the startup file with the model parameters 
            set to the values from the endof  iteration “??” (?? means =
            number of iteration output file.) In addition 
            to get "iter" attributes , replace the marker"**" by "ccam_iter_"
            
    ==================  ==========  ==========================================
    Attributes          Type        Description     
    ==================  ==========  ==========================================
    **iteration         int         Occam2D iteration number represented on
                                    the file.
    **param_count       int         assumes that the model parameter
    **iteraion_data     array       iteration data output after inversion
    ==================  ==========  ==========================================
    
    """
    def __init__(self, iter_fn =None, **kwargs): 
        Startup.__init__(self, **kwargs)
        
        self._logging =csamtpylog.get_csamtpy_logger(self.__class__.__name__)
        
        self.iter_fn =iter_fn
        self.iter_data =None 
        
        self.iter_date_time =datetime.datetime.now().strftime('%m-%d-%Y %H:%M:%S')
        
        for key in list(kwargs.keys()): 
            setattr(self, key, kwargs[key])
            
        if self.iter_fn is not None : 
            self.read_occam2d_iterfile()
        
    def read_occam2d_iterfile(self, iter_fn=None ):
        """
        read_occam iteration file  and populate attributes 
        
        :param iter_fn: full path to iteration file 
        :type iter_fn: str 
        
        :Example:
            
            >>> from pycsamt.modeling.occam2d import Iter 
            >>> path =os.path.join(os.environ ['pyCSAMT'], 'pycsamt', 
            ...                       'data', 'occam2D', 'ITER17.iter')
            >>> iter_obj=Iter(iter_fn= path)
            ... print(iter_obj.occam_iter_param_count)
            ... iterDATA= iter_obj.occam_iter_data
        """

        if iter_fn is not None : self.iter_fn = iter_fn
        self._logging.info (
            'Read Occam2D iteration file and set convenient'
            ' attributes :<%s>' % self.iter_fn)

        if SB.which_file(filename=self.iter_fn)=='iter': 
            with open(self.iter_fn, 'r', encoding='utf8') as fsup : 
                occam_iter_datalines =fsup.readlines ()

        for ii, iteritems in enumerate(occam_iter_datalines):
            if iteritems.lower().find('param count') >=0 : 
                iterdata =occam_iter_datalines[ii+1:]
            if iteritems.find(':')>=0 : 
                iteritems = iteritems.strip().split(':')
                iter_attr = iteritems[0].lower().replace(
                    ' ', '_').replace('!','')
                try :
                    self.__setattr__(
                        'iter_' + iter_attr, float(iteritems[-1]) )
                except : 
                    self.__setattr__(
                        'iter' + iter_attr, iteritems[-1]) 
                    pass 

        # temp =[] # get the data from iterfile  
        #and collect model resistivities .
        temp=[]
        for iteritems in iterdata : 
            sfitems = iteritems.strip().split()
            temp.extend(sfitems)
        # temp=[np.array([float(item) for item in 
        # iteritems.strip().split()]) for iteritems in iterdata ]

        self.iter_data =np.array([float(mm) for mm in temp])

      
class Response (Data): 
    """
    Occam response file . Response is paired with the iteration file above.
    One is output at the end of each successful iteration.
    Inherets from Occam 2D data file .
    
    Arguments
    ---------
        **response_fn** : str 
                full path to Occam 2D  response_file 
        **data_fn** : str 
                full path to Occam 2D data file 
            
    .. note:: To get "response attributes" replace "**" by "resp_" and '***'
                means 'resp_+occam_dtype'.An example below to help 
                understanding how to get attributes from this class.
                
    ====================  ===============  ====================================
    Attributes              Type            Description 
    ====================  ===============  ====================================
    occam_mode              str             occam  mode of survey area 
    occam_dtye              str             Type of data provided 
    **receiver_number       int             number of sites names 
    **frequencies           array_like      frequency array 
    **datatype              array_like      datatype from OccamData file (TE, 
                                            TM etc..) 
    ** data_value           array_like      occam2d data value in OccamDatafile .
    ** forward_data         array_like      data value produced by the forward 
                                            code using
                                            the model from ITER??.iter  
    ***residual              array_like     normalized difference i.e
                                            (input data – forward data) / error.
                                            The “error” value is from the input
                                            data file.
    ***forward              ndarray         array of forward mode resistivities 
    ***forward_phase        ndarray         array of computed phase forward mode 
    ***strike_angle         ndarray         array of strike angle , it rotation 
                                            is applied 
    ====================  ===============  ====================================
    
    .. note:: easy way to get all attributes for other purposes  
            "occam_mode" means the data type provided 
            If "occam_dtype"= tm_log10 then 
            resp_obj.resp_{occam_dtype}_residual = resp_obj.resp_tm_log10_residual
            
    .. seealso:: `Data.occam_mode` and ` Data.data_type `
        
    :Example :
        
        >>> from pycsamt.modeling.occam2d import Response
        >>> resp_obj =Response (data_fn =os.path.join(os.environ[pyCSAMT], 
        ...                                 'pycsamt', data, occam2D,
        ...                                 'OccamDataFile.dat'), 
        ...                    response_fn =os.path.join(os.environ[pyCSAMT], 
        ...                                     'pycsamt', data, occam2D,
        ...                                      'RESP17.resp')
        >>> occam_mode =resp_obj.resp_occam_mode 
        >>> occam_data_type =resp_obj.occam_dtype 
        >>> occam_mode_forward =resp_obj.resp_{occam_dtype}_forward
        >>> occam_mode_residual =resp_obj.resp_{occam_dtype}_residual 
        >>> occam_mode_phase_error = resp_obj.resp_{occam_dtype}_residual_phase_error
    """
    response_params =['sites', 'frequencies', 'data_type',
                      'strike_angle', 'data', 'forward', 'residual']
    
    def __init__(self, response_fn =None , data_fn =None,  **kwargs):
        Data.__init__(self, data_fn, **kwargs)
        
        self.resp_fn =response_fn 
        self.data_fn =data_fn 
        
        self._logging =csamtpylog.get_csamtpy_logger(self.__class__.__name__)
 
        for key in self.response_params : 
            self.__setattr__('resp_'+key, None)
        self.resp_data =None 
        
        for keys in list(kwargs.keys()): setattr(self, key, kwargs[keys])
        if self.resp_fn is not None : 
            self.read_occam2d_responsefile()
            
    def read_occam2d_responsefile(self, response_fn=None): 
        """
        Read Occam2D response file.
        
        :param response_fn: full path to response file 
        :type response_fn:str
        
        :Example:
            
            >>> from pycsamt.modeling.occam2d import Response 
            >>> pathresp =os.path.join(os.environ ['pyCSAMT'], 'pycsamt',
            ...                           'data', 'occam2D',RESP17.resp )
            >>> path_data =os.path.join(os.environ ['pyCSAMT'], 'pycsamt',
            ...                            'data', 'occam2D',OccamDataFile.dat )
            >>> resp_obj = Response(response_fn=pathresp, data_fn = path_data )
            >>> respDATA= resp_obj.resp_data_value 
            >>> resp_obj.occam_mode 
            >>> resp_obj.occam_dtype
            >>> forward_data = resp_obj.resp_forward_value
            >>> residual_data = resp_obj.resp_residual_value
            >>> tm_log10= resp_obj.resp_tm_log10
            >>> tm_phase=resp_obj.resp_tm_phase
            >>> tm_forward = resp_obj.resp_tm_log10_forward
            >>> tm_residual = resp_obj.resp_tm_log10_forward
            >>> tm_forward = resp_obj.resp_tm_log10_residual
            >>> tm_phase_error = resp_obj.resp_tm_phase_err
            >>> tm_phase_err = resp_obj.resp_tm_log10_forward_phase
            >>> tm_residual_phase_err = resp_obj.resp_tm_log10_residual_phase_err
        """
        if response_fn is not None : 
            self.resp_fn =response_fn 
        self._logging.info ('Read Occam 2D response file <%s>' % self.resp_fn)
        
        if SB.which_file(filename= self.resp_fn) =='resp': 
            with open (self.resp_fn , 'r', encoding ='utf8' ) as fresp : 
                occam_response_lines = fresp.readlines()
    
        # for respitems in occam_response_lines : 
        temp= [np.array([float(newitem) for newitem 
                         in respitems.strip().split()]) 
               for respitems in occam_response_lines] 
        self.resp_data = func.concat_array_from_list(list_of_array=temp)  
        
        # container of original data 
        for ii, keyattr in enumerate(['receiver_number', 'frequencies', 'data_type',
                                      'strike_angle_value', 'data_value',
                                      'forward_value', 'residual_value']) : 
            self.__setattr__(''.join(['resp_',keyattr]),self.resp_data[:, ii])
        
        #--->  populate specific  attributes 
        
        def truncate_response (data , station_data, station_names ):  
            """
            function to truncate data and put on dictionnay for all stations  and 
            their coresponding values 
            
            :param data: data ndarray(depth , station_location)
            :type data: ndarray 
            
            :param station_data: from response data 
            :type station_data: array_like  
            
            :param station_names: list of sites names 
            :param station_names: list
            
  
            :returns dico: dictionnary of station value 
            :rtype dico: dict 
            """ 
            dico, sta={},0
   
            for ii, item in enumerate(station_data): 
                if int(item) != int(station_data[ii+1]):
                    if int(item) == len(station_names) -1 :
                        # tale before last array and 
                        dico[station_names[int(item)-1]]= data [sta:ii+1] 
                        # add the last array
                        dico[station_names[int(item)]] =data [ii+1:]
                        break
                    else : 
                        dico[station_names[int(item)-1]]= data [sta:ii+1]
                        sta=ii+1 
            return dico 
  
        #----------> build matrix iter2data [resp_{occam mode} ,
        #resp_{occam mode}_forward,
        # -------> resp_{occam mode}_residual ]
        #Note to access attribute , use "resp_+occam_mode _
        #{response data{forward, strike_angle or residual }}
        # set frequencies into corresponding value 
        
        self.resp_frequencies =np.array([
            int(freq) for freq in self.resp_frequencies])
        self.resp_data_type = np.array([
            int(dt) for dt in self.resp_data_type])
        
        # get occam mode from resp_data_type 
        occam_mode , *ignore = np.unique (
            self.resp_data_type, return_counts=True)
        
        # Once occam mode is get then create attribute with occam datatype

        keymode =[] # container of new datatype  eg [tm_log10, tm_phase]
        occam_mode =[int(oo) for oo in occam_mode]  # PUT occam mode on list 
        
        for imode  in  occam_mode : 
            for key, mode in Data.occam_data_type.items(): 
                if imode == mode : 
                    keymode.append(key)
                    #set attrinut as numb freq -1 * num sites 
                    #np.zeros((self.resp_data.shape[0],))) 
                    setattr(self, 'resp_' + key ,
                            np.zeros((
                                len(self.data_frequencies)-1) * len(
                                    self.data_sites),)) 
                    
                    # for each mode , try to figure out each forward,
                    # residual value 
                    #--- > initialise attributes eg : 
                        #resp_+{occam_mode}_{forward|residual}
                        # add forward attribut for all datatype 
                    for dt in ['data', 'forward', 'residual', 'strike_angle']: 
                        setattr(self, 'resp_' + key +'_' + dt,
                                np.zeros((len(
                                    self.data_frequencies)-1, len(
                                        self.data_sites)))) 
                        
                    for dt in ['forward_phase', 
                               'residual_phase_err', 'strike_angle_err']: 
                        setattr(self, 'resp_' + key +'_' + dt,
                                np.zeros((len(self.data_frequencies
                                              )-1, len(self.data_sites)))) 
                    
        #set attribute occam mode and occam _data type get from response file  
        self.__setattr__('occam_mode', occam_mode)
        self.__setattr__('occam_dtype', keymode)

        if self.verbose > 0: 
                
            print('{0:-^77}'.format('Occam 2D Response  infos'))
            print('** {0:<27} {1} {2}'.format('Occam data type', '=',
                                              tuple(self.occam_dtype)))
            print('** {0:<27} {1} {2}'.format('Occam data mode', '=', 
                                              tuple(self.occam_mode)))
        
       
        # 
        if len(self.occam_mode)> 2 :
            mess ="".join(["Your data seems a Magnetotelluric (MT) data.",
                           " pyCSAMT works only with Controlled Source",
                           " Audio-frequency Magnetotelluric (CSAMT).", 
                           "Althrough the both methods are close , we are ",
                           "not able to process these data properly beyond ", 
                           "4 modes simultaneously and more.To avoid ",
                           "misleading and miscomputation on your data, ",
                           ", we strongly recommend you to use  MTpy software;",
                           " an efficiently python toolbox solve ",
                           " magnetotelluric problem. You can easy download "
                           " and install via the following link :", 
                           " https://github.com/MTgeophysics/mtpy/wiki ."])
            warnings.warn(mess)
            self._logging.debug(mess)
                           
        # resp_tm_log10' for mode = 5, and resp_tm_phase for mode =6 ,
        #initialise values  to get their corresponding value 
        # raw value are not truncated , get directly from response file                  
        raw_forward = np.zeros((len(self.data_frequencies)-1
                                ) * len(self.data_sites),)
        raw_forward_phase = np.zeros((len(self.data_frequencies)-1)
                                     * len(self.data_sites),)
        
        raw_residual = np.zeros((len(self.data_frequencies)-1)
                                * len(self.data_sites),)
        raw_residual_phase_err = np.zeros((len(self.data_frequencies)-1)
                                          * len(self.data_sites),)
        
        raw_stike_angle =np.zeros((len(self.data_frequencies)-1)
                                  * len(self.data_sites),)
        raw_stike_angle_err =np.zeros((len(self.data_frequencies)-1)
                                      * len(self.data_sites),)
        
        raw_observedData =np.zeros((len(self.data_frequencies)-1) 
                                   * len(self.data_sites),)
        raw_observedData =np.zeros((len(self.data_frequencies)-1) 
                                   * len(self.data_sites),)
        
        # get an corresponding value for each attributes (looping response file )
        step=0          # gap to find the phase value , if not find at 
                        #the next reading then add one as step
        mm=0            # indice for value (data, forward , residual )
        nn=0            # indice for phase 
        # loop the response file and get the different data -Type value
        #[5,6][tm_log10, tm_phase]
        for ii, rowline in enumerate(self.resp_data) :  
            for jj , (imode , ikey) in enumerate(zip (
                    self.occam_mode, self.occam_dtype)): 
                if int(self.resp_data[ii, 2]) == imode: 
                    
                    if jj%2 ==0 : # mean is rho value len of data type =2
                        # then repeat if find the same value for next reading
                        #to_fill=  getattr(self, 'resp_' + keymode[jj]) 
                        # self.resp_tm_log10 |self.resp_tm_phase
                        getattr(self, 'resp_' + ikey)[mm]=rowline[4]
                        
                        raw_forward[mm] = rowline[-2]
                        raw_residual[mm]=rowline[-1] 
                        raw_stike_angle[mm]=rowline[3] 
                        raw_observedData[mm]= rowline[-3]
                        
                        mm +=1  # incremente to get new index to fill data 
                        step +=1 # step to get the new
                        
                    if jj %2  != 0 : # if phase for mm # when is odd 
                                    #--> mean a phase, need to add step 
                        nn = nn + step 
                        try:# self.resp_tm_phase
                            getattr(self, 'resp_' + ikey)[nn]=rowline[4] 
    
                            raw_forward_phase[nn], raw_residual_phase_err[nn]\
                                , raw_stike_angle_err[nn]  = \
                                rowline[-2], rowline[-1], rowline[3]
                            raw_observedData[nn]= rowline[-3]
                            step =0 # initialise  for next reading 
                        except :
                            self._logging.error(
                                'Index Error ! out of bounds for axis 0.')
                            pass
        
        # build a block data of forward residual with 
        #shape = number frequency , num_sites 
        # set response value

        rr=0
        for dd , dval in zip(['forward', 'residual', 'strike_angle'],
                             [raw_forward, raw_residual, raw_stike_angle]) :
            
            for ii in range(len(self.data_sites)): 
                for kk in range(len(self.data_frequencies)-1):
                    for mode in self.occam_dtype: 
                        getattr(self, 'resp_{0}_{1}'.format(
                            mode, dd))[kk, ii] = dval[rr] 
                        
                    rr +=1
           
            rr=0 # 
        # create matric for observed data : to get observedData: 
            #use mode +observedData   ex: tm_obsData
        rr=0
        # for dd , dval in zip(['data', raw_observedData]) :
        for ii in range(len(self.data_sites)): 
            for kk in range(len(self.data_frequencies)-1):
                for mode in self.occam_dtype: 
                    getattr(self, 'resp_{0}_{1}'.format(mode,'data')
                            )[kk, ii] = raw_observedData[rr] 
                    # self.__setattr__('resp_{0}_{1}'.format(mode, dd), vv[])
                rr +=1
           
            # rr=0 # 
       
        # set phase and phase error # do the dame for errors 
        rr=0   
        for dd , dval in zip([
                'forward_phase', 'residual_phase_err', 'strike_angle_err'],
            [raw_forward_phase, raw_residual_phase_err, raw_stike_angle_err]) : 
            for ii in range(len(self.data_sites)): 
                for kk in range(len(self.data_frequencies)-1):
                    for mode in self.occam_dtype: 
                        getattr(self, 'resp_{0}_{1}'.format(mode, dd)
                                )[kk, ii] = dval[rr] 
                        # self.__setattr__('resp_{0}_{1}'.format(mode, dd), vv[])
                    rr +=1
            rr=0 # 

        # resetting mode value  like [self.resp_tm_log10 | self.resp_tm_phase]
        # let get the value of resp_tm_log10 as example 
        for mode in self.occam_dtype : 
            vv= getattr(self, 'resp_{}'.format(mode)) 
            # let reinitialize our  resp_tm_log10 attribute
            # so to get the matrix length as nfreq-1 , nsites
            ll=0
            self.__setattr__('resp_{}'.format(mode), np.zeros((
                len(self.data_frequencies)-1, len(self.data_sites)))) 
            # let loop to repopulative value 'in 
            #shape len(data frequency, len(sitenames) )
            for ii in range(len(self.data_sites)): 
                for kk in range(len(self.data_frequencies)-1):
                    getattr(self, 'resp_{0}'.format(mode))[kk, ii] = vv[ll] 
                    ll +=1
        # sometimes some processing softwares doesnt respect 
        #th Occam2D criteria to range frequency from Highest to lower . 
        # to be sure that this criteria is respected , let flip
        # all data to highest to lowest frequency
        # especially avoid miscomputation when the frequencies are
        # not sorted according this principle.
        # if freq is range to increase order , let flip the data so
        # to get the frequency as occam 2D range 
        #  Highest frequency to lowest frequency 
        
        if self.data_flip_freq : # get attribute from inheret Data object 
            for occ in self.occam_dtype :
                # flip the mode attribute e.g : tm_log10
                self.__setattr__('resp_{0}'.format(occ), 
                                     np.flipud(getattr(
                                         self, 'resp_{0}'.format(occ))))
                # flip other attribute eg : tm_log10_forward
                for attr in ['forward', 'residual', 'strike_angle','data'] : 
                    self.__setattr__('resp_{0}_{1}'.format(occ, attr), 
                                     np.flipud(getattr(
                                         self, 'resp_{0}_{1}'.format(occ, attr))))
                #flip other attribute eg : tm_log10_forward_phase
                for attr in [
                        'forward_phase', 'residual_phase_err', 'strike_angle_err']: 
                    self.__setattr__('resp_{0}_{1}'.format(occ, attr), 
                                     np.flipud(getattr(
                                         self, 'resp_{0}_{1}'.format(occ, attr))))
                
        # display infos 
        for occ  in self.occam_dtype:
            for fwo in ['Forward', 'Residual']:
                print('** {0:<27} {1} {2}'.format(
                    fwo+ ' shape:' + occ.split('_')[0].upper(
                        ) + " "+ occ.split('_')[1],
                                                  '=', 
                    getattr(self, 'resp_{0}_{1}'.format(
                        occ, fwo.lower())).shape))
    
    
          
class Mesh(object): 
    """
    Read Occam read mesh file 
    
    Arguments
    -----------
        **mesh_fn** : str 
            full path to occam_2D mesh file 
            
    .. note: to get "mesh" attributes , replace the marker"**" by "mesh_"
    
    ====================  ===============  ====================================
    Attributes(**=mesh_)  Type              Description 
    ====================  ===============  ====================================
    **x_nodes              array_like       horizontal “nodes” (i.e. blocks + 1)
    **z_nodes              array_like       vertical “nodes” (i.e. layers + 1)
    **mesh_values          ndarray          array of mesh files parameters 
                                            xnodes,z_nodes
    ====================  ===============  ====================================
    
    """
    def  __init__(self, mesh_fn=None , **kwargs)  : 
        self._logging =csamtpylog.get_csamtpy_logger(self.__class__.__name__)
        self.mesh_fn =mesh_fn 
        for key in ['x_nodes', 'z_nodes']:
            self.__setattr__('mesh_'+key, None )
        
        for key in list (kwargs.keys()): 
            setattr(key, kwargs[key])
            
        if not hasattr(self, 'verbose'): 
            setattr(self, 'verbose', 0)
        if self.mesh_fn is not None : 
            self.read_occam2d_mesh()
    
    def read_occam2d_mesh(self, mesh_fn =None): 
        """
        Read mesh and get mesh values data and populates attributes 
        
        :param mesh_fn:   full path to mesh file 
        :type mesh_fn: str 
        """
        if mesh_fn is not None : 
            self.mesh_fn= mesh_fn 
        if self.mesh_fn is not None :
            self._logging.info('read Occam2d meshfiles <%s>'% self.mesh_fn)
            if SB.which_file(filename =self.mesh_fn )=='mesh': 
                 with open(self.mesh_fn, 'r') as fmesh : 
                    occam_mesh_lines = fmesh.readlines()
            else : 
                mess ='No Mesh file detected. Please provide the'\
                    ' right occam2d mesh files.'
                warnings.warn(mess), self._logging.error(mess)
                raise CSex.pyCSAMTError_occam2d(mess)
         # characteristic of the mesh (nblocks +1)
        mesh_char = occam_mesh_lines[1].strip().split()  
        # horizontal nodes , verticales nodes (nlayers+1)
        num_h_nodes ,num_v_nodes = int(mesh_char[1]), int(mesh_char[2]) 

        #initialise attribute x nodes and z nodes 
        self.mesh_x_nodes =np.zeros(num_h_nodes)
        self.mesh_z_nodes =np.zeros (num_v_nodes)
        self.mesh_values = np.zeros((num_h_nodes, num_v_nodes, 4), dtype=str)
        
        # take off the both descriptive lines 
        starting_lines_counter = 2          
        count_xnodes , count_z_nodes, mesh_index= 0 , 0, 0
        for mii, mitems in enumerate(occam_mesh_lines[starting_lines_counter:]): 
            mitems = mitems.strip().split()
            for mval in mitems: 
                self.mesh_x_nodes[count_xnodes]=float(mval)
                count_xnodes += 1 
            # asssume to get (number of blocks =x_nodes
                # -1 -1 of Python readlines )
            if count_xnodes >= num_h_nodes -1 :  
                starting_lines_counter += mii
                break 
        # print(self.mesh_x_nodes)
        # print(starting_lines_counter)
        #sometimes there is between xnodes and verticales 
        #nodes a lines of spaces , skip it if exists 
        # and go to the line to count the vertical nodes 
        if occam_mesh_lines[starting_lines_counter +1] =='' or\
            occam_mesh_lines[starting_lines_counter +1] =='\n': 
            starting_lines_counter +=1
 
        #--> read vertical nodes part 
        for mii, mitems in enumerate(occam_mesh_lines[
                starting_lines_counter +1:]): 
            mitems =mitems.strip().split() 
            for mj, mvalues in enumerate(mitems): 
                self.mesh_z_nodes [count_z_nodes]= float(mvalues)
                count_z_nodes +=1 
                
            # substract the number of blocks =1 with the Python index 
            if count_z_nodes >= num_v_nodes -1 : 
                starting_lines_counter += mii 
                break 
        # be sure line counter match the 
        
        # Check whether we are at the end of nodes  for consistency , it 
        #There should always be a zero before the next section  
        # so add +1 to index line 
        try :      
            if int(occam_mesh_lines[starting_lines_counter +1]) == 0:                               
                starting_lines_counter +=1
        except :pass

        # --> fill model values with parameter specs {???}
        self._logging.info("Reading of parameters"
                           " specification of mesh and fill model values")
        
        for specs, paramline in enumerate(occam_mesh_lines[
                                    starting_lines_counter+ 1:], 
                                   starting_lines_counter):
            paramline = paramline.strip()
            if mesh_index == num_v_nodes or paramline.lower(
                    ).find('exception') > 0:
                break
            else:
                list_paramspecs = list(paramline)
                if len(list_paramspecs) != num_h_nodes - 1:
                    if self.verbose >0:
                        print('-'*77)
                        mess= ''.join([
                            "---> Params specification for line = {0}"
                            " from <{1}> have too many columns.", 
                            "Should be {2} instead of {3}. "
                            "Please check your mesh file."])
                        print(mess.format(specs, os.path.basename(self.mesh_fn),
                                          num_h_nodes,len(list_paramspecs)  ))
                        print('-'*77)
       
                    list_paramspecs =  list_paramspecs[0:num_h_nodes]
                for kk in range(4):
                    for jj, mvalue in enumerate(list( list_paramspecs )):
                        self.mesh_values[jj, mesh_index, kk] = paramline [jj]
                mesh_index += 1

        #Note : from coocam 2D,  number of hozotonal nodes = nblocks +1 
        #                       number of vertical nodes = nlayer +1 
        # remove me zero stay on array if number of horizontal node 
        #deosnt match the array x_nodes x_shape 
        self.mesh_x_nodes = self.mesh_x_nodes[np.nonzero(self.mesh_x_nodes)]
        # do the same for z_nodes 
        self.mesh_z_nodes = self.mesh_z_nodes[np.nonzero(self.mesh_z_nodes)] 
        
        if self.mesh_x_nodes.shape[0] != num_h_nodes: 
            new_h_nodes = self.mesh_x_nodes.shape[0]  # reshape the mesh vales 
            #self.mesh_values.resize(new_h_nodes , num_v_nodes, 4)
            self.mesh_values =np.resize(self.mesh_values,
                                        (new_h_nodes , num_v_nodes, 4) )
        else: new_h_nodes= num_h_nodes
            
        if self.mesh_z_nodes.shape[0] != num_v_nodes:
            new_v_nodes = self.mesh_z_nodes.shape[0]
            #self.mesh_values.resize(new_h_nodes, num_v_nodes, 4)
            self.mesh_values =np.resize(self.mesh_values,
                                        (new_h_nodes , num_v_nodes, 4))

        if self.verbose> 0: 
            print('{0:-^77}'.format('Occam 2D Mesh params '))
            
            for im , nodes  in zip (['Horizontal', 'Vertical'],
                                    [[new_h_nodes,num_h_nodes ], 
                                     [new_v_nodes,num_v_nodes]]) : 
                mess='*** {0} nodes read = {1} instead of {2} in '\
                    'mesh files.'.format(im, nodes[0], nodes[1])
                    
                print(mess) 
                self._logging.info(mess)
      
        # generate mesh x_grid and z_grid

        self.mesh_x_grid = self.mesh_x_nodes.copy()
        self.mesh_x_grid = np.append(self.mesh_x_grid, self.mesh_x_grid[-1])
       
        self.mesh_x_grid = np.array([self.mesh_x_grid[:ii].sum()
                                for ii in range(self.mesh_x_grid.shape[0])])
        self.mesh_x_grid -= self.mesh_x_grid.mean()
        self.mesh_z_grid = np.array([self.mesh_z_nodes[:ii].sum()
                                for ii in range(self.mesh_z_nodes.shape[0])])
        
        if self.verbose >0 : 
            print('---> {0:<20} {1} {2}'.format('Horizontal nodes',
                                                '=',self.mesh_x_nodes.shape[0]))
            print('---> {0:<20} {1} {2}'.format('Vertical nodes',
                                                '=',self.mesh_z_nodes.shape[0]))
        
  
class Iter2Dat (object): 
    """
    Iter2Dat is a format converter which convert *.iter file and related mesh
    files to so called 'x,y,z' *.data file for post-processing.
    The class was inpired from the Bo Yang matlab script . reading functions come
    from 'plotOccam2DMT.m' routine
    
    .. seealso:: `Occam routine <http://marineemlab.ucsd.edu/Projects/Occam/sharp/index.html >`
                of SIO, UCSD.
            
    Bo Yang matlab-script is on `add.info` sub-packages of pyCSAMT.
    If your are familiar with matlab and you want to rewrite  the script
    please contact the author at: 
        
    Bo Yang, 2011.
    ------------------
        China Univ. of Geosciences, Wuhan, China.
        yangbo.cug@163.com.
        Copyright 2011-2016 Bo Yang.
        $Revision: 1.0 $ $Date: 2012/04/05 21:50:20 $
        Revision log:
        2012/04/04 : Version 1.0 released.
        
    
    Arguments
    ----------
        **path_to_occamfiles** : str 
                path to occamfiles : put on one directory 
        
    ==================  ==============  =======================================
    Attributes              Type         Description
    ==================  ==============  =======================================
    ** Data             obj             Occam Data object 
    **Iter              obj             occam Iteration object
    **Response          obj             OccamResponse object 
    model_fn            str             full path to occam model file 
    iter_fn             str             full path to iteration file 
    mesh_fn             str             full path to occam mesh file 
    data_fn             str             full path to occam data file 
    doi                 str of float    depth of investigation. default is 1km
                                        if your povide  value on float,
    model_x_nodes       array_like      rescalled station offsets according 
                                        to doi.
    model_z_nodes       array_like      rescalled depth offsets accordng to doi
    model_res           ndarray         rescalled model resistivity according
                                        to doi.shape len(model_z_nodes,
                                        model_x_nodes)
    station_names       list            list of sites names 
    station_location    array_like      station offsets  
    elevation           array_like      elevation value  , to rewrite a file 
                                        if elevation is given ,
                                        will take elevation              
    ==================  ==============  =======================================
    
    """
    iter2d_params =['model_x_nodes', 'model_z_nodes', 'model_res', 'station_names', 
                    'station_location', 'elevation']
    
    def __init__(self, iter2dat_fn =None , model_fn =None , data_fn=None , 
                 iter_fn=None,  mesh_fn=None,
                   **kwargs): 
                 
        self._logging =csamtpylog.get_csamtpy_logger(self.__class__.__name__)
        self.iter2dat_fn   =iter2dat_fn 
        self.model_fn =model_fn 
        self.data_fn =data_fn
        self.iter_fn = iter_fn 
        self.mesh_fn =mesh_fn 
        
        self.bln_fn =kwargs.pop('bln_fn', None)
        self.savepath =kwargs.pop('savepath', None)
        self.OccamModel=kwargs.pop('occam_model_obj', None)
        self.verbose  = kwargs.pop('verbose', 0)
        
        for keys in list(kwargs.keys()):
            self.__setattr__(keys, kwargs[keys])
            
        # initialise the main attributes
        for key in self.iter2d_params :self.__setattr__(key, None )  
       

        if self.iter2dat_fn is not None or self.model_fn is not None : 
            self.read_iter2dat()
        

    def read_occamfiles(self , path_to_occamfiles=None,  **kws): 
        """
        getffiles and readfiles to populates speciales attributes 

        :param path_to_occamfiles: full path to occamfiles [ITER|DATA|RESP] files. 
        :type path_to_occamfiles: str
                
        :Example:
            
            >>> path =os.path.join(os.environ ['pyCSAMT'], 'pycsamt', 'data', 'occam2D')
            >>> iter2_obj=Iter2Dat(path_to_occamfiles= path)
            >>> sites= iter2_obj.OccamData.occam_data_sites
            >>> freq =iter2_obj.OccamData.occam_data_frequencies
            >>> iter_roughn =iter2_obj.OccamIter.occam_iter_roughness_value
            >>> iter_misfit=iter2_obj.OccamIter.occam_iter_misfit_reached
            >>> forward_data = iter2_obj.OccamResponse.forward_data
            >>> residual_data = iter2_obj.OccamResponse.residual
        """
        def truncate_response (data , station_data, station_names ):  
            """
            function to truncate data and put on dictionnay for all stations  and 
            their coresponding values 
  
            :param data: arrray of data ndarray(depth , station_location)
            :type data: array_like 
            
            :param station_data: from response data 
            :type station_data: array_like , dtype '<U12'
                
            :param station_names: list , list of sites names
            :type station_names: list 
                
            :returns:  dictionnary  of station value 
            :rtype: dict 
            """ 
            dico, sta={},0
   
            for ii, item in enumerate(station_data): 
                if int(item) != int(station_data[ii+1]):
                    if int(item) == len(station_names) -1 :
                        # tale before last array and 
                        dico[station_names[int(item)-1]]= data [sta:ii+1] 
                        # add the last array
                        dico[station_names[int(item)]] =data [ii+1:]
                        break
                    else : 
                        dico[station_names[int(item)-1]]= data [sta:ii+1]
                        sta=ii+1 
            return dico 

       

    def  write_iter2dat_file(self,  filename=None ,model_fn=None, iter_fn=None ,
                             mesh_fn=None , data_fn=None , 
                              doi ="1km",  savepath =None ,occam_model_obj=None , 
                              negative_depth=True, scale='km', **kws): 
        """
        write  'x,y,z' *.data file for post-processing
        can read and rewrite iter2dat file 
        
        :param x: offset
        :type x: array_like 
        
        :param y: depth 
        :type y: array_like 
        
        :param z: log10 resistivities
        :type z: array_like 
        
        ==============  ==============  =======================================
        params          Type            Description 
        ==============  ==============  =======================================
        model_fn        str             full path to occam model file 
        iter_fn         str             full path to iteration file 
        mesh_fn         str             full path to occam mesh file 
        data_fn         str             fll path to occam data file 
        filename        str             output of iter2dat file 
        doi             str of float    depth of investigation. default is 1km
                                        if your povide  value on float,  
                                        default unit is "meter" eg : 2000=2000m 
        negative_depth  bool            if True , will provide file with 
                                        negative depth value 
        scale           str             scaled the offset value and elevation . 
                                        might be [m|km]
        elevation       ndarray|list    can be provided if usefull 
        ==============  ==============  =======================================
        
        1. Write with Occam 2D files  
        
        :Examples:
 
            >>> from pycsamt.modeling.occam2d import Iter2Dat as i2d
            >>> data='OccamDataFile.dat'
            >>> mesh = 'Occam2DMesh'
            >>> model = 'Occam2DModel'
            >>> path =os.path.join(os.environ ['pyCSAMT'], 'pycsamt', 
            ...                       'data', 'occam2D')
            ...                   #,'OccamDataFile.dat')
            >>> pathi2d =os.path.join(os.environ ['pyCSAMT'], 'pycsamt', 'data', '_iter2dat_')
            >>> occam_iter2dat_obj =i2d(mesh_fn=os.path.join(path, mesh), 
            ...                    iter_fn = os.path.join(path, iter_), 
            ...                    model_fn =os.path.join(path, model), 
            ...                    data_fn =os.path.join(path, data))
            >>> occam_iter2dat_obj.write_iter2dat_file()
                
        2. Rewrite the  with iter2dat file 
        
        :Example:
             
            >>> from pycsamt.modeling.occam2d import Iter2Dat as i2d
            >>> iter_='ITER17.iter'
            >>> idat = K1.iter.20142.dat
            >>> bln = K1.iter.20142.bln
            >>> occam_iter2dat_obj =i2d(iter2dat_fn=os.path.join(pathi2d, idat),
                                        bln_fn = os.path.join(pathi2d, bln)
            >>> occam_iter2dat_obj.write_iter2dat_file()
        """
        
        iter2dat_fn =kws.pop('iter2dat_fn', None)
        bln_fn =kws.pop('bln_fn', None)
        elev=kws.pop('elevation', None)
        
        if elev is not None :
            self.elevation = elev 
        if bln_fn is not None : 
            self.bln_fn = bln_fn 
        if doi is not None :
            self.doi =doi
        
        if savepath is not None : 
            self.savepath =savepath
        # statement to do so to 
        if iter2dat_fn is not None : 
            self.iter2dat_fn= iter2dat_fn
            
        elif self.iter2dat_fn is None :
            if model_fn is not None :
                self.model_fn =model_fn 
            if occam_model_obj is not None :
                self.OccamModel = occam_model_obj
            
            # resseting other main attributes 
            if iter_fn is not None : setattr(self, 'iter_fn', iter_fn)
            if mesh_fn is not None : setattr(self, 'mesh_fn', mesh_fn)
            if data_fn is not None : setattr(self, 'data_fn', data_fn)
      
        # then read file 
        
        self.read_iter2dat()
        
        self._logging.info (
            'Get ready to write iter2dat_file with model file or model_obj.')
        
        
        #--------WRITE station location bln file ---------------------
        #--specs the name of the file is not given 
      
        if filename is None : 
            if self.iter2dat_fn is not None : 
                filename = os.path.basename(
                    self.iter2dat_fn)[:-4].lower() +\
                    '.{0}'.format(datetime.datetime.now().month)
            else : 
                filename ='iter{0}.{1}{2}'.format(
                    int(self.iter_num), self.iter_roughness,
                    datetime.datetime.now().month)
            
        #-
        if self.elevation is None : 
            if self.verbse :
                mess='!Elevation is not provided. Set to 0.'
                print('-->'+mess)
                
            self._logging.debug(mess)
            
            if self.station_location is not None :  # assume station location exist
                self.elevation  = np.repeat(0., len(self.station_location))
            
        else : 
            if isinstance(self.elevation , (list, tuple)) :
                self.elevation  = np.array(self.elevation )
            assert len(self.elevation) == len(self.station_location) ,\
                CSex.pyCSAMTError_occam2d(
                    'Elevation provided must be the same with offsets :'\
                      ' length|size ={0}.'.format(self.station_location.size))
        
        
        self._logging.info('Now write an associate file call *bln include'
                           ' :offset+elev+station name: Can be used'
                           ' straiforwardly for Golden software plot.')
           
        write_bln_lines =[]
        if scale is None :scale ='km'
        if scale =='km' :  # convert to  kilometer elev ,  offset  and  depth.
            dz = 1e3
        elif scale =='m' :
            dz= 1.
            
        self.model_z_nodes /= dz
        self.model_x_nodes /=dz
        
        if self.station_location is not None :
            self.station_location /= dz
            self.elevation  /=dz
            # then can write bln file 
            for  ikey, stn in enumerate(self.station_names) : 
                write_bln_lines.append(''.join([
                    '{0:<7.6f},'.format(self.station_location[ikey]), 
                    '{0:<7.6f},'.format(self.elevation[ikey]),
                    '{0:<4}'.format(stn), '\n']))
        #---- end of writing bln file ---- 
        
        #---WRITE iteration to DATA ----
        # loop all the grid makes :flatern all array s (expensive way)
        # x , y , z  = mesh_x_grid.flattern(),
        # mesh_z_grid.flattern(),new_model_res_matrix.flattern() 
        
        self._logging.info ('Writing a iter2dat file from specifics objects.'
                            ' Can be visualized  by Golden software "Surfer".')
        
        if negative_depth is True : self.model_z_nodes *=-1
        write_iter2data_lines =[]
        for ii in range(len(self.model_z_nodes)) : 
            for jj in range(len(self.model_x_nodes)): 
                write_iter2data_lines.append(''.join([
                    '{0:>15.7f}'.format(self.model_x_nodes[jj]),
                    '{0:>15.7f}'.format(self.model_z_nodes[ii]), 
                    '{0:>15.7f}'.format(self.model_res[ii,jj]),
                    '\n',
                    ]))

        #writes files 
        for ii , wfiles in enumerate([write_iter2data_lines, write_bln_lines ]): 
            if ii == 0 : mm='.dat'
            else :mm='.bln'
            if ii ==1 :
                if self.station_location is  None : break

            with open(''.join([filename,'{}'.format(mm)]), 'w') as fw : 
                fw.writelines(wfiles)
                
        #savefile
        self.savepath = func.cpath(self.savepath, '_savei2dout_')
        if self.savepath is not None : 
            import shutil 
            try :
                for jj, file  in enumerate([filename +'.dat', filename+'.bln']): 
                    if  jj==1 : 
                        if self.station_location is None : break 
         
                    shutil.move(os.path.join(os.getcwd(),file) ,
                        os.path.join(self.savepath , file))
            except : 
                warnings.warn("It seems the files already exists !")
        if self.station_location is None :
            print('---> Only file {0}.dat have been'
                  ' successfully written to  <{1}>.'.format(
                      filename, self.savepath))
            
            mess =''.join(['! Error writing "{0}.bln" file . Could not ', 
                           'write Golden software (*.bln) file without ',
                           'station names and station location.'])
            
            mess=mess.format(filename)
            
            warnings.warn('---> !'+mess)
            self._logging.debug(mess)
 
        else : 
          
            print('---> files {0}.dat & {0}.bln have been '
                     'successfully written to  <{1}>.'.format(
                         filename, self.savepath))
        
    
    def read_iter2dat(self, iter2dat_fn =None, bln_fn =None, scale ='km',
                      model_fn=None, iter_fn=None , mesh_fn=None , doi='1km',
                      data_fn=None , occam_model_obj=None, **kws): 
        """
        Read YangBo iter2data file or provided  Occam 2D specific files .

        :param iter2dat_fn: full path to iter2dat file 
        :type iter2dat_fn: str  
        
        :param bln_file: full path to bln file 
        :type  bln_file: str  
        
        :param scale: str ,  scale of output data . Most of time ,
                    Bo yang iter2Dat file is converted in kilometer . If not turn 
                        the scale to `m`.
        :type param scale: str 
            
        :Example:
            
            >>> from pycsamt.modeling.occam2d import Iter2Dat 
            >>> i2d='iter17.2412.dat'
            >>> i2d_2='K1.iter.dat'
            >>> bln='iter17.2412.bln'
            >>> bln_2 = 'K1.bln'
            >>> pathiter = os.path.join(os.path.dirname(os.environ ['pyCSAMT']),
            ...                            '_iter2dat_', i2d_2)
            >>> pathbln = os.path.join(os.path.dirname(os.environ ['pyCSAMT']),
            ...                           '_iter2dat_', bln_2)
            >>> occam_iter2dat_obj =Iter2Dat(iter2dat_fn=pathiter, 
            ...                                 bln_fn =pathbln ) 
            >>> i2d_data = occam_iter2dat_obj.read_iter2datfile()
            >>> i2d_sta, i2d_depth = occam_iter2dat_obj.model_x_nodes,
            >>> occam_iter2dat_obj.model_z_nodes
            >>> i2d_model_res = occam_iter2dat_obj.model_res
        """
        f=0 # flag to read data from file or attributes populating 
        
        if scale.lower() in ['km', 'kilometer']:
            dz= 1000.
        elif scale.lower()in ['m', 'meter'] : 
            dz=1.
        else : 
            warnings.warn('---> Iter2Dat file is generally converted into km.'
                  ' so provided the specify unit of yourf file.')
                
            raise CSex.pyCSAMTError_occam2d_iter2dat(
                'Wrong scale ! . Might be "m" or "km".')
            
        elev =kws.pop('elevation', None)
        if elev is not None :self.elevation  =elev 

        def validate_iter2dat_file (iter2dat_fn) : 
            """
            validate iter2dat fn and return array of data 
            """
            newf =[]
            if  os.path.isfile(iter2dat_fn) is False:
                mess ='<%s> is not a file. Please check'\
                    ' your right path.' % iter2dat_fn
                self._logging.error(mess)
                warnings.warn(mess)
                
                raise CSex.pyCSAMTError_occam2d_iter2dat(mess)
            else : 
                with open(iter2dat_fn, 'r') as fi2d: 
                    iter2dat_readlines =fi2d.readlines()
                for ii, item in enumerate(iter2dat_readlines) :
                    try : 
                        item =[float(i2d) for i2d in item.strip().split()] 
                    except : 
                        mess ='<{0}> is not Bo Yang data file.'\
                            ' Please provide the right iterartion file.'
                        warnings.warn(mess)
                        self._logging.error(mess)
                        raise CSex.pyCSAMTError_occam2d_iter2dat(mess)
                    else :    
                        newf.append(item)
            return func.concat_array_from_list(list_of_array=newf)
        
         # statement , which file |data should be read 
        if iter2dat_fn is not None : 
            self.iter2dat_fn =iter2dat_fn 
    
            
        if self.iter2dat_fn is not None :  f=2
        else :f= 1      # read module and fill attributes 

        if f == 1 : 
            if model_fn is not None : self.model_fn =model_fn 
            
            self._logging.info (
                'Get ready to write iter2dat_file with model file or model_obj.')
            
            if occam_model_obj is not None : 
                self.OccamModel = occam_model_obj 
       
    
            #-->  get exception if files are poperly given 
            if self.model_fn  is not None :
                if iter_fn is not None : setattr(self, 'iter_fn', iter_fn)
                if mesh_fn is not None : setattr(self, 'mesh_fn', mesh_fn)
                if data_fn is not None : setattr(self, 'data_fn', data_fn)
                
                if self.iter_fn is None and self.mesh_fn is \
                    None and self.data_fn is None : 
                     mess =''.join([
                         'Can not read iter2data files !None object is detected.', 
                    ',  Must supply at least mesh_file, iteration file ',
                    'and data file  or provided model obj.'])
                     warnings.warn(mess)
                     self._logging.error(mess)
                     raise CSex.pyCSAMTError_occam2d_iter2dat(mess)
                    
                self.OccamModel = Model(model_fn = self.model_fn ,
                                        iter_fn =self.iter_fn , 
                                        mesh_fn = self.mesh_fn )
                
            elif self.model_fn is None and self.OccamModel is\
                None  and self.data_fn is None : 
                mess =''.join([
                    'Can not read iter2data files !None object is detected.', 
                    'Please provide modelfile, mesh file,',
                    ' data file or occam model obj .'])
                warnings.warn(mess)
                self._logging.error(mess)
                raise CSex.pyCSAMTError_occam2d_iter2dat(mess)
            elif self.OccamModel is not None :
                # try to get an spacific attribute from Model 
                try : 
                    getattr(self, 'model_station_offsets')
                except : 
                    mess=''.join([
                        'Object provided is not an Python Model class object.',
                    ' Please provide a model object from occam2d module. e.g', 
                    '--> from pycsamt.modeling.occam2d import Model .'])
                    
                    warnings.warn(mess)
                    self._logging.error(mess)
                    raise CSex.pyCSAMTError_occam2d_iter2dat(
                        'Object provided is not a Model class object.'\
                        ' Please provide a right model obj. ')
            # thin now everything is fine then :    
            #----Get model_attribute from model_obj ------- 

            # get data offset _obj 
            occam_data_obj =Data(data_fn = self.data_fn)
            offs = occam_data_obj.data_offsets
            
             # # get the last iteration number
             # and roughness value from Iter obj 
            occam_iter_obj =Iter(iter_fn=self.iter_fn)
            
            self.__setattr__('iter_num', occam_iter_obj.iter_iteration)
            self.__setattr__('iter_roughness',round(
                occam_iter_obj.iter_roughness_value) )
            
            # get largest model attributes (nodes x ,
            # nodes z  and model resistivity )
            model_plot_x = self.OccamModel.model_station_offsets
            model_plot_z = self.OccamModel.model_depth_offsets
            model_res= self.OccamModel.model_resistivity
            
            # Matrix is large enough , let get the part
            #we need from investigation depth 
            #Note : from CSAMT : 1km is enough 
            #them we set 1km as a ddefault value .
            
            self._logging.info(
                'Get a corresponding range of matrix we need ,'
                ' offsets, depth and rho and populate attributes.')

            self.doi =punc.depth_of_investigation(doi=doi)
            # check the doi and the max z_nodes depth 
            #if user change the default doi value 
            
            if self.doi > model_plot_z.max() : 
                mess ='doi provided = {0} m is much larger '\
                    'than maximum investigation depth ={1}.'\
                    'Set new doi = {1} m.'.format(doi, model_plot_z.max())
                
                warnings.warn(mess)
                self._logging.debug(mess)
                print('---> ! Input doi ={0} m was '
                      'set to maxdepth ={1} m.'.format(
                          doi,model_plot_z.max() ))
                
                doi = model_plot_z.max()
            # slice matrix and get range 
            new_station_offsets,  new_depth_offsets, new_model_res_matrix= \
                punc.slice_csamt_matrix(block_matrix =model_res ,
                                   offset_MinMax=(offs[0], offs[-1]),
                                    station_offsets=model_plot_x, 
                                    depth_offsets=model_plot_z , 
                                    doi =self.doi)
            
            # resetting attribute :
                # Defalut values are in  meter as default units
            
            for i2dkey, i2dvalues  in zip(self.iter2d_params,
                                          [new_station_offsets,new_depth_offsets,
                                           new_model_res_matrix,
                                           occam_data_obj.data_sites, 
                                           occam_data_obj.data_offsets,
                                           self.elevation ]) :
                setattr(self, i2dkey, i2dvalues)
                
         # ---> NOW Go to the case where data iter2dat file
         # is provided, read file and populate attributes   
        elif f==2 : 
            
            iter2dat_data = validate_iter2dat_file(self.iter2dat_fn)
           
            # Bo Yang  iteration data  are range into  negative 
            #depth , then we will resset to positive depth 
            iter2dat_data[:,1] = np.apply_along_axis(
                lambda i2d :i2d * -1 , 0, iter2dat_data[:,1])
          
            # get depth offset and station offset  and sorted values
            #according to station range the min to max . 
            # sometimes  Bo yang files  are little missy, reorganize it  
            
            iter2dat_data=iter2dat_data[iter2dat_data[:,1].argsort(
                kind="mergesort")]
            
            # loop and ressetting attributes ---
            # resetting attributes 
            for ii , name  in enumerate(self.iter2d_params[:2]): 
                value,*ign  = np.unique(iter2dat_data[:,ii], return_counts=True) 
                setattr(self,  name, value *dz)
    
            # build resistivitvity model
            self._logging.info('Building Resisistivity'
                               ' model from file <%s>' % os.path.basename(
                                   self.iter2dat_fn))
            # initialise model resistivity array 
            self.model_res = np.zeros ((self.model_z_nodes.shape[0],
                                              self.model_x_nodes.shape[0])) 
            mm=0
            for ii in range(len(self.model_z_nodes)) : 
                for jj in range(len(self.model_x_nodes)): 
                    self.model_res[ii, jj] = iter2dat_data[:, 2][mm] 
                    mm +=1 
                    
            if self.verbose >0: 
                print('{0:-^77}'.format('Iter2Dat *Data* info'))
                
                print('**{0:<37} {1} {2}'.format(
                    ' Mininum Limit vertical nodes (m)','=' ,
                    self.model_z_nodes.min() +1 ))
                print('**{0:<37} {1} {2}'.format(
                    ' Maximum limit vertical nodes (m)','=' , 
                    self.model_z_nodes.max() +1 ))
                print('**{0:<37} {1} {2}'.format(
                    ' Minimum limit Horizontal nodes (m) ',
                    '=' , self.model_x_nodes.min()+1))
                print('**{0:<37} {1} {2}'.format(
                    ' Maximum limit Horizontal nodes(m)',
                    '=' , self.model_x_nodes.max()+1))
                print('**{0:<37} {1} {2}'.format(
                    ' Horizontal nodes get from file ',
                    '=' , self.model_res.shape[1] +1))
                print('**{0:<37} {1} {2}'.format(
                    ' Vertical nodes get from file ',
                    '=' , self.model_res.shape[0] +1))
    
            if bln_fn is not None :
                self.bln_fn = bln_fn 
                
            if self.bln_fn is None or os.path.isfile(self.bln_fn) is False : 
                mess=''.join([
                    ' ! Site file is not provided. Need  a Bo',
                    ' Yang bln file read from Golden software',
                    ' to set station location , elevation and sites names.', 
                    'Please profile a file <*.bln>, if not, station location,',
                    ' elevation and sites names will set to None.'])
                warnings.warn(mess)

            if self.bln_fn is not None : 
                if os.path.isfile (self.bln_fn ) is True : 

                    with open(self.bln_fn , 'r') as fbln :
                        bln_readlines = fbln.readlines()
                    temb=np.array([item.strip().split(',') 
                                   for item in bln_readlines ])
                    
                    for ii , ss in  enumerate([
                            'station_location', 'elevation' , 'station_names']) :
                        if ii ==2 :
                            self.__setattr__(ss,   temb[:, 2])
                        else : 
                            try :
                                # self.__setattr__(ss, np.apply_along_axis
                                #(lambda off : float(off) * dz, 0 , temb[:, ii] ))
                                self.__setattr__(ss, np.array([
                                    float(off) for off in temb[:,ii]]) * dz)
                            except : 
                                mess =" It seem all data from station and elevation "\
                                    "can not be convert into float value."\
                                        " Please check your data."
                                raise CSex.pyCSAMTError_occam2d_iter2dat(mess)
            
            if self.verbose > 0:
                print('{0:-^77}'.format('Iter2Dat *Station* info'))  
                try : 
                    for ffmt, fmtvalue in zip ([
                            ' Stations num.',' Minimum offset (m)',
                            ' Maximum offset (m)' ], 
                            [len(self.station_names),
                             self.station_location.min(),
                             self.station_location.max() ]):
                        print('**{0:<27} {1} {2}'.format(ffmt,'=' , fmtvalue))
                except : 
                    for ffmt in ([' Stations num.',' Minimum offset (m)',
                                  ' Maximum offset (m)' ]):
                        print('**{0:<27} {1} {2}'.format(ffmt,'=' , None))
                    pass 
    
                    
                if np.all(self.elevation == 0.) or self.elevation is None: 
                    print('*** ! Elevation not provided.')
                else : 
                    print('**{0:<27} {1} {2}'.format(' Minimum elevation (m)',
                                                     '=' , self.elevation.min()))
                    print('**{0:<27} {1} {2}'.format(' Maximum elevation (m)',
                                                     '=' , self.elevation.max()))
                
                
class occam2d_write(object):
    """
    Special class to build occam2d imput files with :ref:`MTpy` module .
    
    Arguments
    --------
        **edi_fn**: str 
            full path to edifiles locations
        **freq_num** :float 
            number of frequencies to use in inversion
        **interpolate_freq**: bool,
            frequency interpolation , default is *False*
        **geoelectric_strike**: bool
             geoelectric strike angle assuming N = 0, E = 90.
             If True , provided , losgspace interpolation as tuple value
        
    Others important attributes can be  found in : 
        
    ======================  ===================================================
    Key Words/Attributes    Description
    ======================  ===================================================
    edi_fn                  full path to data file
    n_layers                number of vertical layers in mesh
                            *default* is 31. 
    num_layers              [ int ] number of regularization layers.
    num_z_pad_cells         number of vertical padding cells below 
    iterations_to_run       maximum number of iterations to run
                            *default* is 20
    resistivity_start       starting resistivity value.  If model_values is
                            not given, then all values with in model_values
                            array will be set to resistivity_start
    save_path               directory path to save startup file to
                            *default* is current working directory  
    startup_basename        basename of startup file name. 
                            *default* is Occam2DStartup
    startup_fn              full path to startup file.
                            *default* is save_path/startup_basename  
    target_misfit           target misfit value.
                            *default* is 1.
    x_pad_multiplier        horizontal padding cells will increase by this
                            multiple out to the edge of the grid.
                            *default* is 1.7
    z1_layer                thickness of the first layer in the model.
                            Should be at least 1/4 of the first skin depth
                            *default* is 10
    z_bottom                bottom depth of the model (m).  Needs to be large 
                            enough to be 1D at the edge. 
                            *default* is 200000.0 
    z_target_depth          depth to deepest target of interest.  Below this
                            depth cells will be padded to z_bottom
    cell_width              width of cells with in station area in meters
                            *default* is 100
    phase_te_err            percent error in phase for TE mode. *default* is 5
    phase_tm_err            percent error in phase for TM mode. *default* is 20. 
    res_te_err              percent error in resistivity for TE mode. 
                              *default* is 10
    res_tm_err              percent error in resistivity for TM mode.
    trigger                 [ float ] multiplier to merge model blocks at 
                            depth.  A higher number increases the number of
                            model blocks at depth.  *default* is .1.12
    model_mode              model mode to use for inversion, see module`Data`. 
    ======================  ===================================================
    
    .. note:: We consider occam2D buildingInputs file is focused on CSAMT data in
            TM mode then default configuration as  setting according this feature. 
            If input `edi_fn` are MT data , please resetting configuration
            you SEG `edi-fn` data provided !.
            ... 
            
    :Example:
        
        >>> from pycsamt.modeling.occam2d import occam2d_write 
        >>>  occam2d_write.buildingInputfiles(os.path.join(os.environ['pyCSAMT'], 
        ...                                              'data', 'edi'), 
        ...                        savepath =os.path.join(os.environ['pyCSAMT'],
        ...                                               'data', 'tesocc2' ),
        ...                        geolectricke_strike=34.)
  
    """
    verbose = 0 
    if SUCCESS_IMPORT_MTPY is False :
        
            mess = ''.join([
                'Failed to import MTpy module ! '
                'Could not build inputOccam2D files.', 
                            ' Please try to install :ref:`MTpy` manually !'])
            print('---> '+ mess)
    
   
    @staticmethod 
    def buildingInputfiles(edi_fn, freq_num =None , savepath =None ,
                               interpolate_freq =False, 
                               geoelectric_strike =None , **kwargs): 
        """
        Method to build Occam2D inputfiles. Deal with :ref:`MTpy` module. 
        Try to install :ref:`MTpy` is not installed yet.
        
        :param edi_fn: full path to edifiles 
        :type edi_fn: str 
        
        :param freq_num:  number of frequencies to use in inversion
        :type freq_num: float
        
        :param interpolate_freq: frequency interpolation , default is *False*
        :type interpolate_freq: bool 
        
        :param geoelectric_strike: geoelectric strike angle assuming N = 0, E = 90.
                            If True, provided, losgspace interpolation as tuple value
        :type geoelectric_strike: float
        
        :return: outfiles building occam2d inputfiles 
        :rtype: str , sys.stdout
   
        """
        
        if SUCCESS_IMPORT_MTPY is False : return 
           
        freq_logspace =kwargs.pop('intp_freq_logspace', None)
        startup_basename =kwargs.pop('startup_basename', 'Startup')
        
        res_tm_error_floor =kwargs.pop('res_tm_err', 10.)
        phase_tm_error_floor =kwargs.pop('phase_tm_err', 20.)
        
        occam_model_mode =kwargs.pop('occam_mode', '6')
        
        res_te_error_floor =kwargs.pop('res_te_err', 10.)
        phase_te_error_floor =kwargs.pop('phs_te_err', 20.)
        
        number_of_layers = int(kwargs.pop('n_layers', 31.))
    
        cell_width =kwargs.pop('cell_width', 5.)
        x_pad_multiplier = kwargs.pop('x_pad_multiplier', 1.7)
        trigger =kwargs.pop('trigger', 1.12)
        
        z_bottom =kwargs.pop('z_bottom', 5000.)
        z1_layer=kwargs.pop('z1_layer', 5.)
        z_target_depth =kwargs.pop('z_target', 1100.)
        
        iterations_to_run =kwargs.pop('iteration_to_run', 100.)
        resistivity_start =kwargs.pop('resistivity_start', 1.)
        OccamDataFile = kwargs.pop('occamDataFile','OccamDataFile.dat' )
        
        
        savepath = func.cpath (savepath ,'_saveoc2dInputfiles_' )
        # collected the list of stations from edifile object 
        
        slst=[edi[0:-4] for edi in
              os.listdir(edi_fn) if edi.find('.edi')>0]
        
    
        # create an occam data object
         # if not specified will calculate from the data
        if freq_logspace is not None: 
            #freq_logspace= (-1,4,17 )
            freq_logspace = np.logspace( *freq_logspace)  
        _logger.info('Read occam2d Data,'
                     ' write data and build regularization mesh ')
        
        ocd = MToccam2d.Data(edi_path=edi_fn,
                           station_list=slst,
                           interpolate_freq=interpolate_freq,
                           freq= freq_logspace 
                           )
        ocd.save_path = savepath
        ocd.freq_num = freq_num# number of frequencies to invert
        
        #### make data file
        # geoelectric strike for rotation
       
        ocd.geoelectric_strike = geoelectric_strike
        
        if occam2d_write.verbose >0: 
            print('---> geoelectric strike added !')
        
       
        if occam_model_mode not in ['{}'.format(mm) 
                                    for mm in list(Data.occam_dataType.keys())]: 
            msg=''.join(['Occam model mode provided is wrong !', 
                         ' Please select the convenient mode between',
                        '  the following list modes: {0}'])
                
            _logger.error('occam model value given is out'
                          ' of the model mode range . ')
            print('-- > '+ msg.format([
                '{}'.format(mm) for mm in list(Data.occam_dataType.keys())]))
            return 
        
         # added error floors
        if occam_model_mode in ['5','6']: 
            ocd.model_mode= "6"
            ocd.res_tm_err = res_tm_error_floor
            ocd.phase_tm_err = phase_tm_error_floor
        elif occam_model_mode in ['1', '2']: 
            ocd.res_te_err = res_te_error_floor 
            ocd.phase_te_err = phase_te_error_floor 
            
        elif occam_model_mode =='4': 
            ocd.model_mode ="4"
            ocd.res_tm_err = res_tm_error_floor
            ocd.phase_tm_err = phase_tm_error_floor
            ocd.res_te_err = res_te_error_floor 
            ocd.phase_te_err = phase_te_error_floor 
            
        if occam2d_write.verbose >0:     
            print('---> Errors floors successfully  added !')
        
        # now write occam_data_file 
        ocd.write_data_file()
        if occam2d_write.verbose >0: 
            print('---> Read occam2D data and write occam Data done !')
        # make model and mesh files
  
        ocr = MToccam2d.Regularization(ocd.station_locations, 
                                       n_layers =number_of_layers, 
                                       cell_width = cell_width, 
                                       x_pad_multiplier = x_pad_multiplier, 
                                       trigger= trigger, 
                                       z1_layer = z1_layer,
                                       z_target_depth = z_target_depth , 
                                       save_path=ocd.save_path,
                                       z_bottom =z_bottom ,
                                        )
    
        ocr.build_mesh()
        ocr.build_regularization()
        ocr.write_mesh_file()
        ocr.write_regularization_file()
        ocr.plot_mesh()
        if occam2d_write.verbose >0: 
            print('---> Build occam2d Regularization mesh  done !')
        # make startup file
        ocs=MToccam2d.Startup(iterations_to_run=iterations_to_run, 
                              startup_basename=startup_basename, 
                              data_fn=os.path.join(ocd.save_path,
                                                   OccamDataFile), 
                              resistivity_start=resistivity_start)
        ocr.get_num_free_params()
        ocs.param_count=ocr.num_free_param
        ocs.save_path=ocd.save_path
        ocs.model_fn=ocr.reg_fn
        ocs.write_startup_file()
        
        if occam2d_write.verbose >0: 
            print('---> Build write occam2D startup file done !')
            
            print('{0:-^77}'.format('Summary *occam2d input params* infos'))  
            
            print('**{0:<27} {1} {2}'.format(' Given frequency number',
                                             '=' , freq_num))
            print('**{0:<27} {1} {2}'.format(' Interpolate frequencies range',
                                             '=' , (freq_logspace[0], freq_logspace[1])
                                             ))
            
            print('**{0:<27} {1} {2}'.format(' Occam model mode',
                                             '=' , occam_model_mode))
            
            if occam_model_mode in ['6', '5']:
                print('**{0:<27} {1} {2} %.'.format(' TM rho error floors',
                                                    '=' , res_tm_error_floor))
                print('**{0:<27} {1} {2} %.'.format(' TM phase error floors',
                                                    '=' , phase_tm_error_floor))
            elif occam_model_mode in ['1', '2']:
                print('**{0:<27} {1} {2} %.'.format(' TE rho error floors',
                                                    '=' , res_te_error_floor))
                print('**{0:<27} {1} {2} %.'.format(' TE phase error floors',
                                                    '=' , phase_te_error_floor))
                
            print('**{0:<27} {1} {2}'.format(' Model cell width',
                                             '=' , cell_width ))
            print('**{0:<27} {1} {2}'.format(' model horizontal pad',
                                             '=' , x_pad_multiplier ))
            print('**{0:<27} {1} {2}'.format(' Model bricks trigger',
                                             '=' , trigger ))
            
            print('**{0:<27} {1} {2}'.format(' Number of model layers',
                                             '=' , int(number_of_layers) ))
            print('**{0:<27} {1} {2} m.'.format(' Top layer thickness',
                                                '=' , z1_layer ))
            print('**{0:<27} {1} {2} m.'.format(' Expected image depth',
                                                '=' , z_target_depth))
            print('**{0:<27} {1} {2} m.'.format(' Model bottom',
                                                '=' , z_bottom ))
            
            print('**{0:<27} {1} {2}'.format(' Expected iteration to run',
                                             '=' , iterations_to_run ))
            print('**{0:<27} {1} {2} ohm.m'.format(' starting model resistivity',
                                                   '=' , np.power(10,
                                                        resistivity_start) ))
            
            print('**{0:<27} {1} +{2} degrees E of N'.format(' Geoelectric strike',
                                                        '=' ,geoelectric_strike))
            
            print('-'*77)  

        print('---> Building occamInputfiles  successfully done. ')  
 
         
@mdeco.geoplot2d(reason='misfit', cmap ='bwr', climits=(-2,2), show_grid=False)
def getMisfit(resp_fn =None, data_fn =None, kind='rho', **kwargs) : 
    """
    Calling *getMisfit* to plot misfit value using `geoplot2d` decorator 
    see :ref:`pycsamt.viewer.plot.geoplot2d` to get the kwargs arguments. 
    Most keywords arguments are the same used by :mod:`pycsamt.viewer.plot.Plot2d`
    Please refer to `mod` documentation. 
    For more details , refer to :class:pycsamt.viewer.plot.geoplot2d
    
    :warning:: when using decorator`geoplot2d`, set always *reason* 
                argument to `misfit`.
    
    :param kwargs: use the main argument from getMisfit function
    :type kwargs: dict 
    
    :return:* resp_misfit, 
            * resp_sites_names,
            * resp_sites_offsets,
            * resp_freq, doi , 
            * model_rms, 
            * model_roughness

    :Example:
        >>> from pycsamt.modeling.occamd2d import getMisfit 
        >>> occamPath ='data/occam2d
        >>> response_file ='RESP27.resp'
        >>> data_file = 'OccamDataFile.dat' 
        >>> log_file = 'LogFile.logfile'
        >>> getMisfit(response_fn = os.path.join(occamPath, response_file)
        ...            data_fn = os.path.join(occamPath, data_file), 
        ...            logfile =  os.path.join(occamPath, log_file))

    """
    
    # call reponse Object 
    mode =kwargs.pop('mode', None)
    oclogfile =kwargs.pop('logfile', None)
    
    useResiValue= kwargs.pop('residual', False)
    if kind.lower().find('rho')>=0 or \
        kind.lower().find('res')>=0  or kind==1: 
        kind='rho'
    elif kind.lower().find('ph')>=0 or kind==2: 
        kind = 'phase'
    else :
        raise CSex.pyCSAMTError_occam2d(
            f'Wrong argument kind {kind!r}.'
            'Should be `rho` or `phase`')

    resp_obj  =Response(response_fn =resp_fn , data_fn =data_fn)
    
    # get response values 
    #------------------------STATEMENT RESPONSE OBJECT ----------------------
   
    make_mode = resp_obj.occam_dtype + [str(mm)for mm in resp_obj.occam_mode]
    resp_occam_dtype_obj = resp_obj.occam_dtype
    
    #---------------------------------MANAGE OCCAM PLOT MODE -----------------
    # if mode is not provided , then take the first occam mode 
    if mode is None : 
        mode = resp_occam_dtype_obj [0]
        mode_phase = resp_occam_dtype_obj[1]

    elif mode is not None: 
        mode =str(mode).lower() 
        # check the mode if provided 
        if mode not in  make_mode  : 
            mess =''.join([
                'Occam mode provided ={0} is wrong !. Occam2D data mode is ={1}'.
                           format(mode, make_mode ),
                           'Please select the right mode.'])
            warnings.warn(mess)
            _logger().error (mess)                     
            raise CSex.pyCSAMTError_occam2d(f' Mode `{mode}` is wrong!'
                         'Do you mean {0}{1}?'.format(*resp_occam_dtype_obj))
        elif mode not in resp_occam_dtype_obj: 
            raise CSex.pyCSAMTError_occam2d(f' Mode `{mode}` is wrong!'
                         'Do you mean {0}{1}?'.format(*resp_occam_dtype_obj))

    # check the mode provided , can be str 
    for im , imode in enumerate(resp_obj.occam_mode) :
        if imode ==mode : 
            resp_occam_dtype_obj  =  resp_occam_dtype_obj [im]
            
    for im , imode  in enumerate( [str(mm) for mm in resp_obj.occam_mode]): 
        if imode ==mode : 
            resp_occam_dtype_obj  =  resp_occam_dtype_obj [im]  
        

    resp_forward = getattr(resp_obj, 'resp_{0}_forward'.format(mode))
    occam_data = getattr(resp_obj, 'resp_data_value')
    
    occam_data =getattr(resp_obj, 'resp_{0}_data'.format(mode))
    resp_freq = np.log10(getattr(resp_obj, 'data_frequencies'))

    resp_sites_names = getattr(resp_obj, 'data_sites')
    resp_sites_offsets = getattr(resp_obj, 'data_offsets')
    

    # get rms and roughness value if occamLogfile is provided 
    if oclogfile is not None : 
        log_obj = occamLog(fn =oclogfile )
        model_rms = log_obj.rms.min()
        rms_index = np.where(log_obj.rms == model_rms)
        index= rms_index[0]
        if len(index)>1 : 
            index =int(index[0])-1
        else : index = int(index) -1 
        model_rms = round(model_rms,2)
        model_roughness= round(log_obj.roughness[index],3) 
        
    else : 
        mess ='OccamLogfile is not given.'+\
            'R.M.S & roughness =<?>'
        _logger.info(mess)
        print(mess)
        model_rms ='?'
        model_roughness ='?'
        
    if kind=='rho': 
        resp_misfit = .01 * np.sqrt(
            (occam_data -resp_forward)**2/(occam_data**2) )
        
        if useResiValue is True : 
            print('---|> Plot residual !')
            _logger.info('Plot Occam2D residual data !')
            
            resp_misfit =getattr(resp_obj, 'resp_{0}_residual'.format(mode))
        
        else : resp_misfit *= 100 # we plot in percentage 
        
    elif kind =='phase': 
 
        phase_resi= getattr(resp_obj, 'resp_{0}_residual'.format(mode_phase) )
        if mode_phase.find('tm')>=0 : 
            phase_data= resp_obj.tm_phase 
        else : phase_data = resp_obj.te_phase  
        
        fw_phase = getattr(resp_obj, 'resp_{}_forward_phase'.format(mode))
        # resp_misfit =phase_resi
        resp_misfit = ((phase_data - fw_phase)/phase_resi)%90/100
        
    print('{0:=^77}'.format('Misfit Infos : Kind = {}'.format(kind.upper()))) 

    print('** {0:<37} {1} {2} {3}'.format('Misfit max ','=',
                                          resp_misfit.max(), '%' ))
    print('** {0:<37} {1} {2} {3}'.format('Misfit min','=',
                                          resp_misfit.min(), '%' ))
    print('-'*77)
    
    return (resp_misfit, resp_sites_names,
            resp_sites_offsets, resp_freq, model_rms, 
            model_roughness)  
  

@mdeco.geoplot1d(reason = 'resp', color_mode='color', font_size=3.,
           grid_kws= {'color':'k', 'ls':'-', 'lw':.2, 'alpha':.5, 
                      'which':'major'},
           linebbox_kws={'boxstyle':'sawtooth','facecolor':
                         'whitesmoke', 'color':'white'})
def plotResponse(data_fn =None, resp_fn =None, stations =None, **kws):
    """ 
    Find Responses properties and plot values and could read multiples 
    lines provided that each line as its occam datafile and its response file.
    
    :class:`~modeling.occam2d.Response` inherits of `Data` class then 
    easy to data attributes for plot purposes. Visualize all error from 
    raw to inversion using `error_type` arguments. For instance:: 
        
        - `1` : visualize the misfit compute manually
        - `2`:  visualize the raw error from raw occam data 
        - `3`: visualzie the error data and phase  defined as 
                * error = (input data - forward data)/ RESI 
        - `4` or 'residual`': Visualize only  the residual data . 
        
        Default is `residual.`

    :param data_fn: Occam data file. Can be a string or list of datafiles
    :param resp_fn: Response file and can be string or list of response file. 
    :param stations:  str or list of station to visualize. 
    
    .. note:: To visualize multiple files. Can provide only the path 
    where occam (data and response) files are located. In that case, set
    datafile and response file the name like:: 
        
        `data_fn` =['line01.dat' , 'line01.resp', 'line02.dat',
                  'line02.resp' , ...]
    
    :returns: 
        - resp_lines : coloection of survey lines id 
        - resp_stations : Collection of station to visualize 
        - resp_freq: collection of frequencies
        - resp_appRHO: collection of te/tm data and
                    forward data as a tuple(te|tm , fw'te|tm')
        - resp_phase: collection of the te/tm phase data as tuple (data, fw)
        - resp_appRho_err. collection of  misfit rho response 
        - resp_phase_err: collection of the misfit phase  

    :Example: 
        
        >>> from pysamt.modeling.occam2d import plotResponse 
        >>> resPath =r'F:\ThesisImp\occam2D\old\K1'
        >>> plotResponse(data_fn =resPath,
        ...                 stations = ['S00', 'S04', 's08', 'S12']
        ...                  rms =['1.013', '1.451', '1.00', '1.069'], 
        ...                  error_type ='resi' )
    """
    
    error_type = kws.pop('error_type', 'residual')
    model_rms =kws.pop('rms', None)
    
    try : 
        if isinstance(error_type, str): 
            if error_type.lower().find('resi')>=0: 
                error_type=4
        error_type = int(error_type )
    except :
        error_type =1

    # try automatically to find occam data file and response file 
    # providin only the path 
    try : 
        if os.path.isdir(data_fn): 
            resPath = data_fn
        elif os.path.isdir(resp_fn): 
            resPath = resp_fn 
    except : pass

    else : 
        data_fn = [os.path.join(resPath, file) 
                     for file in os.listdir(resPath) 
                     if file.find('.dat')>=0]
        resp_fn = [os.path.join(resPath, file)  
                     for file in os.listdir(resPath) 
                     if file.find('.resp')>=0]
        
    if data_fn is None: 
        raise CSex.pyCSAMTError_AVG(
            'No `Data` file detected. Please provide an occam data file.')
    elif isinstance(data_fn, str): 
        data_fn = [data_fn]
        
    if resp_fn is None: 
        raise CSex.pyCSAMTError_occam2d(
            'No OccamResponse file detected.'
            ' Please provided an occcam Response file')
        
    elif isinstance(resp_fn, str): 
        resp_fn = [resp_fn]
    
    if len(data_fn) != len(resp_fn):
        a_, b_= len(data_fn) , len(resp_fn) 
        raise CSex.pyCSAMTError_inputarguments('Number of Occam datafiles and'
                        'responses files must be the same length.'
                        f' `{a_}` and `{b_}` are given respectively.')
    
    if stations is None : 
        _logger.debug ( 'None sation found can not be plotted.'
                       'Should be considered the default value =`S00`.')
        stations =['S00']
           
    # check the sation and duplicate until the max gridspec =4 
    if isinstance(data_fn, list): 
        if isinstance(stations, str): 
            # duplicate the 
            station_list =[stations for i in range(len(data_fn))]
        elif isinstance(stations, list): 
            if len(stations)< len(data_fn):# max grid spec = 4 
                station_list = stations + ['S00' 
                                          for i in range(4-len(stations))]
            elif len(stations) > len(data_fn): 
                station_list =stations[:4]
            else : station_list= stations 
            
    if len(data_fn) != len(station_list):
        if len(data_fn)< len(station_list): 
            _logger.debug(
                f"Number of stations ={len(station_list)} are larger."
                " Note that each station must fit a given line. "
                f"{'A single' if len(data_fn) ==1 else len(data_fn)}"
                f" line{'s are' if len(data_fn)>1 else' is'} detected. Only "
                f"{len(data_fn)} can be plotted. Data are shrunked to match"
                f" the number of line ={len(data_fn)}.")
            station_list =station_list [:len(data_fn)]
        else:
            raise CSex.pyCSAMTError_value(
                "Data and stations to visualize must have the same "
                f"length. But {len(data_fn)} and {len(station_list)} "
                "were given respectively.")
            
    def read_singleDataResponse( stn_INDEX, dfn , respfn): 
        """ Read single occam response file and datafile 
        Collect attribute for plots. Note `Response` inherits 
        of `Data` class.  
            
        :param stn_INDEX: station index for plotting 
        :param dfn: Data file to get data properties 
        :param resp_fn: Response file to retrieve main properties 
        
        :returns:
            - line: line name 
            - stn: station to plot 
            - frequencies 
            - raw tm/te data and phase te/tm  on tuple 
            - forward te/tm data and phase te/tm data on tuple 
            - error te/tm 
            -error phase 
        """
         
        # line name id from response file
        ln_id = os.path.basename(os.path.splitext(respfn)[0])
    
        respObj = Response(response_fn = respfn , data_fn =dfn )
        # get the name of site data
        try : 
            stnN = respObj.data_sites[stn_INDEX]
        except: 
            len_stnobj =len(respObj.data_sites) 
            _logger.debug(
            f" Too {'large' if stnN >len_stnobj else 'small'}."
            "values. Default station should be `S00`.")
            stnN  = 'S00'
        
        #--> Get occam mode and main response properties # tm_log10, tm_phase
        oc_mode, oc_phs_mode = respObj.occam_dtype [0], respObj.occam_dtype[1]
        # Get frequency in decrease order and remove the last freq 
        freq = respObj.data_frequencies  # array 
        if freq [0] > freq[-1]:  # frequeny ranged in decrease order
            freq = freq [:-1]
            # 
        else : freq =freq [1:]
        
        #--- for resistivity tm data and forward 

        # flag to read occamtype 'log_te_tm'
        cflag=False
    
        if respObj.occam_dtype_ =='log_te_tm': 
            cflag=True 
            
        if cflag is True or respObj.occam_dtype_ =='log_te': 
            tm_data = respObj.te_appRho
            te_error_appRho =respObj.te_error_appRho 
            te_error_phase =respObj.te_error_phase 
    
        if cflag is True or respObj.occam_dtype_ =='log_tm': 
            # Plot as default TM mode when `occam_dtype` is 'log_te_tm
            tm_error_appRho =respObj.tm_error_appRho 
            tm_error_phase =respObj.tm_error_phase 
            tm_data =respObj.tm_appRho 
        # errors from forward
    
        fw_data = getattr(respObj, 'resp_{}_forward'.format(oc_mode ))
        
        # --for phase data tm_phase and forward phase 
        phase_data = getattr(respObj, 'resp_{}'.format(oc_phs_mode ))
        fw_phase = getattr(respObj, 'resp_{}_forward_phase'.format(oc_mode ))
        # 
         # get residual errors  (input data - foward data)/ RESI 
        fw_residual_rho = getattr(respObj,
                                  'resp_{0}_residual'.format(oc_mode) )
        c_= getattr(respObj, 'resp_{0}_residual'.format(oc_phs_mode) )
        fw_phase_residual =c_
        
        # Raw error from occam data 
        if error_type ==1 : 
            misfit_rho = .01 * np.sqrt((tm_data -fw_data)**2/(tm_data**2) )
            misfit_rho *= 100 # we plot in percentage 
            misfit_phase = .01 * np.sqrt(
                (phase_data -fw_phase)**2/(phase_data**2) )
            misfit_phase = (misfit_phase*100)%90 # we plot in percentage 
        elif error_type ==2: 
            if cflag is True or respObj.occam_dtype_ =='log_tm': 
                misfit_rho = tm_error_appRho
                misfit_phase = tm_error_phase
                
            elif respObj.occam_dtype_ =='log_te': 
                misfit_rho = te_error_appRho
                misfit_phase = te_error_phase

        elif error_type ==3: 
            misfit_rho = (tm_data - fw_data)/fw_residual_rho
            misfit_phase = np.abs((phase_data - fw_phase)/fw_phase_residual)%90
            
        elif error_type ==4: # plot only residual 
        
            misfit_rho = fw_residual_rho
            misfit_phase=fw_phase_residual 
            
        # replace all nan value by 0 if exists. 
        misfit_rho[np.isnan(misfit_rho)]=0.
        misfit_phase[np.isnan(misfit_phase)]=0.
        # put frequency on ohm.m not in log 10 
        tm_data = np.power(10, tm_data)
        fw_data = np.power(10, fw_data)
        
        return (ln_id, stnN, freq , (tm_data[:,stn_INDEX ], 
                                     fw_data[:, stn_INDEX]),
            (phase_data[:, stn_INDEX],fw_phase[:, stn_INDEX] ),\
                misfit_rho[:, stn_INDEX], misfit_phase[:, stn_INDEX])

    # Initialise list # read the data and collect data 
    resp_lines, resp_stations, resp_freq, resp_appRHO, resp_phase,\
        resp_appRho_err, resp_phase_err=[[]for i in range(7)]
        
    # Get the index from stations list 
    # return tuple of index stations 
    stnIndex_lst = punc.station_id(station_list) 

    for ii , iline in enumerate(stnIndex_lst): 
        r_l, r_sta, r_freq, r_appRHO, r_phase,\
        r_appRho_err, r_phase_err = read_singleDataResponse(
            stn_INDEX=iline, 
            dfn = data_fn[ii], 
            respfn = resp_fn[ii], 
            )
        resp_lines.append(r_l)
        resp_stations.append(r_sta)
        resp_freq.append(r_freq)
        resp_appRHO.append(r_appRHO)                # tuple (data , fw)
        resp_phase.append(r_phase)                  # tuple(data, fw)
        resp_appRho_err.append(r_appRho_err)        # misfit rho
        resp_phase_err.append(r_phase_err)          # misfit phase
            
    return (resp_lines, resp_stations, resp_freq, resp_appRHO, resp_phase,\
        resp_appRho_err, resp_phase_err, model_rms)
        





















