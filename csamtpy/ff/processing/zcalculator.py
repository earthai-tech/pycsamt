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

.. _module-ZCalculator::`csamtpy.ff.processing.zcalculator`
    :synopsis:  helper functions special calculator 

Created on Thu Dec  3 16:44:29 2020

@author: @Daniel03
"""
import os,warnings 
import numpy as np 
from csamtpy.etc import infos
from csamtpy.utils import exceptions as CSex
from csamtpy.utils import func_utils as func
from csamtpy.utils.decorator import redirect_cls_or_func
from csamtpy.utils._csamtpylog import csamtpylog 

_logger =csamtpylog.get_csamtpy_logger(__name__)
try:
    import scipy

    scipy_version = [int(ss) for ss in scipy.__version__.split('.')]
    if scipy_version[0] == 0:
        if scipy_version[1] < 14:
            warnings.warn('Note: need scipy version 0.14.0 or higher or interpolation '
                          'might not work.', ImportWarning)
            _logger.warning('Note: need scipy version 0.14.0 or higher or interpolation '
                            'might not work.')
    import scipy.interpolate as spi

    interp_import = True

except ImportError: 
    warnings.warn('Could not find scipy.interpolate, cannot use method interpolate'
                  'check installation you can get scipy from scipy.org.')
    _logger.warning('Could not find scipy.interpolate, cannot use method interpolate'
                    'check installation you can get scipy from scipy.org.')
    interp_import = False

############################# END MODULE IMPORT ############################################

mu0= 4* np.pi * 1e-7 
phase_referencing = 2 * np.pi 

def mag_avg (mag_array , A_spacing, Txcur): 
    """
    RAW, E-, or H-field magnitude values)for each frequency
        units : mV/Km*
        
    :param mag_array:  magnnitude value for each data block
    :type mag_array: np.array (ndarray, 1),

    :param a_spacing: dipole length 
    :type a_spacing: float 
    
    :param txcur:  curv coil, transmitter length 
    :type txcurv: float

    :returns: averaged data of magnitude data Block
    :rtype mag_avg: np.ndarray 
    
    """

    
    if type(mag_array) is list : 
        mag_avg =np.array(mag_array) 
    try : 
        mag_avg =np.array([float(ii) for ii in mag_array])
    except : 
        raise CSex.pyCSAMTError_inputarguments("Input first arguments is wrong. must"\
                                               "be a list of ndarray of numbers.")
        
    mag_avg =mag_array.mean()
    En =mag_avg /(A_spacing *  Txcur)
    
    return mag_avg, En 


def phz_avg (phz_array, to_degree): 
    """
    E-, H-field, or Impedance Phase values unit in  mrad
    
    :param phz_array: array of phase values in mrad 
    :type phz_array: array_like 
    
    :param to_degree: ascertain convertion to degree 
    :type to_degree: bool 

    """
    if type (phz_array) is list : 
        phz_array =np.array(phz_array)
        
    try : 
        phz_array =[float(kk) for kk in phz_array ]
    except : 
        raise CSex.pyCSAMTError_inputarguments("input argument must be a number in array or list.")
    
    if to_degree : return (phz_array.mean()/1e3) *180/np.pi, (phz_array/1000) *180/np.pi
    
    return phz_array.mean()
    
def param_rho (rho_array):
    """
    Parameter Average RHO (RHOp), RHO are from each data block, unit ( ohm.m)
    
    :param rho_array: array of resistivity values 
    :type rho_array: array_like 
    """
    if type(rho_array) is list : rho_array=np.array(rho_array)
    try : 
        phz_array =np.array([float(kk) for kk in rho_array ])
    except : 
        raise CSex.pyCSAMTError_inputarguments('Error inputs arguments . must be list or'\
                                               'ndarray of numbers.')
    return  phz_array.mean()

def param_phz (pphz_array, to_degree =False): 
    """
    PHZn are from each data block ,  units in  mrad 
       
    :param pphz_array: average parameters phase for data blocs.
    :type pphz_array: array_like
    
    :param to_degree: ascertain convertion to degree 
    :type to_degree: bool 
    """
    if type(pphz_array) is list : pphz_array=np.array(pphz_array)
    try : 
        pphz_array =np.array([float(kk) for kk in pphz_array ])
    except : 
        raise CSex.pyCSAMTError_inputarguments('Error inputs arguments. must be list or'\
                                               'ndarray of numbers.')
    if to_degree : return (pphz_array.mean()/1000) *180/np.pi, (pphz_array/1000) *180/np.pi
    
    return  pphz_array.mean()

def comp_rho (mag_E_field, mag_H_field, freq_array, A_spacing, Txcurr ):
    """
    Function to compute component average 
    unit in in ohm.m 

    :param mag_E_field: magnitude of E-filed , averaged 
    :type mag_E_field: np.ndarray(ndarray,1)
  
    :param mag_H_field: magnitude of H-Field ,  averaged 
    :type mag_H_field: np.ndarray(ndarray,1)
    
    :param freq_array:  frequency of station field
    :type freq_array: np.ndarray(ndarray,1)
    
    :param  A_spacing: step_between station
    :type A_spacing: np.float
    
    :param Txcurr:  distance of coil in meter 
    :type Txcurr: np.float,

    :returns: comp_rho , component averaged rho.
    :rtype: np.ndarray 
        
    """
    args =[mag_E_field, mag_H_field, freq_array, A_spacing, Txcurr]
    
    for a_args in args : 
        if type (a_args) is not np.ndarray :
            try :
                a_args =np.array([ float(ss) for ss in a_args ])
            except : 
                raise CSex.pyCSAMTError_inputarguments("Elemts composed of each array must be a number.")
        
    # if comp_rho.__code__.co_argcount ==5 : 
    #     if mag_E_field.ndim !=1 or  mag_H_field.ndim !=1 or freq_array.ndim !=1 : 
    #         raise CSex.pyCSAMTError_inputarguments ("dimension of Input arguments must me equal to 1.")
    
    Emag =mag_E_field /(A_spacing * Txcurr)
    Hmag =mag_H_field /Txcurr 
    comp_rho =(Emag/Hmag)**2 * (1/5*freq_array)*np.power(10,6)
    
    return comp_rho

def comp_phz (comphz_array, units ='deg'):
    
    """
    PHZc are from each data block, units in rad 

    :param comphz_array: average parameters phase for data blocs.
    :type comphz_array: float 
    
    :returns: component phase averaged.
    :rtype: component phase averaged.
    
    :Example:
        
        >>> path =  os.path.join(os.environ["pyCSAMT"], 
        ...              'csamtpy','data', 'K1.AVG')
        >>> from csamtpy.core import avg  
        >>> phs_obj =avg.Phase(path)
        >>> phs_obj.loc['S00']
        >>> value, ss = comp_phz(comphz_array=phs_obj.loc[
        ...            'S00'], to_degree=True)
        ... print(value)
    """
    
    if type(comphz_array) is list : comphz_array=np.array(comphz_array)
    
    if comphz_array.dtype not in ["float", "int"] : 
        raise CSex.pyCSAMTError_inputarguments("Type provided is wrong ! number must be float")  

    try : 
        comphz_array=np.array([float(kk) for kk in comphz_array ])
    except : 
        raise CSex.pyCSAMTError_inputarguments('Error inputs arguments. must be list or'\
                                               'ndarray of numbers.')
    
            
    if units =='deg' : return comphz_array.mean() *180/np.pi, comphz_array *180/np.pi
    
    return  comphz_array.mean()


def compute_components_Z_Phz(magn_E_field , magn_H_field, phz_E_field,
               phz_H_field, freq_value, **kwargs):
    """
    Function to compute all components  derived from Impedance Z. 
    user can  enter specifik units in kwargs arguments . program will compute and converts value 
    automatically.
    
    Parameters
    ----------
        * magn_E_field : np.ndarray 
            E_.field magnitude (ndarray,1) in  microV/KM*A
            
        * magn_H_field : np.ndarray 
             H_.field magnitude (ndarray,1)in  mGammas/A or picoTesla/A
             
        * phz_E_field : np.ndarray 
            E_field phase (ndarray, 1) in  mrad 
            
        * phz_H_field : np.ndarray 
             H_field phase (ndarray,1) in  mrad.
             
        * freq_value : np.ndarray 
             Frequency at which data was measured(ndarray,1)in  Hz
             
        * kwargs : str 
             units conversion.
        
    Raises
    ------
        CSex.pyCSAMTError_z(), 
            Exceptions if units entered by the user doesnt match or are messy.

    Returns
    -------
        rho: ndarray 
            Cagnard resistivity calculation. ohm.m 
        phz: ndarray 
             Impedance phase value.
        Zij: ndarray  
             Impedance Tensor value.
        Zreal: ndarray 
            Value of Real part of impedance Tensor.
        Zimag: float
             Value of Imaginary part of impedance Tensor.
        Zreal_imag: ndarray , complex 
            Complex value of impedance Tensor.
         
    :Example: 
        
        >>> from csamtpy.core import avg 
        >>> path =  os.path.join(os.environ["pyCSAMT"], 
        ...              data', 'avg', 'K1.AVG')
        >>> emag_ob = avg.Emag(path)
        >>> hmag_obj = avg.Hmag(path)
        >>> ephz_obj = avg.Ephz(path)
        >>> hphz_obj = avg.Hphz(path)
        >>> freq_obj =avg.Frequency(path)
        >>> station_name ='S00'
        >>> rho, phz, Z, real, imag, comp =compute_components_Z_Phz( 
        ...    magn_E_field=emag_ob.loc[station_name], 
        ...                            magn_H_field =hmag_obj.loc[station_name], 
        ...                            phz_E_field =ephz_obj.loc[station_name], 
        ...                            phz_H_field=hphz_obj.loc[station_name], 
        ...                            freq_value=freq_obj.loc[station_name])
        ... rho
    """
    
    units_E_field =kwargs.pop('unit_E_field', 'microV/km*A')
    units_H_field =kwargs.pop('unit_E_field', 'mGamma/A')
    unit_phz=kwargs.pop('unit_phase', 'mrad')
    unit_freq=kwargs.pop('unit', 'Hz')
    
    #units converter
    units ={units_E_field: {'microv/km*a':1/np.power(10,6),
                            'nv/ma':1/np.power(10,6),
                           'mmv/km*a':1/np.power(10,3),
                           'v/km*a':np.power(10,0)
                   },               # SI : --> E-Field  V/KM*A
            units_H_field : {'mgamma/a':1/np.power(10,6),
                             'pt/a':1/np.power(10,6),
                            'picotesla/a':1/np.power(10,6),
                            'gamma/km*a':np.power(10,3)
                   }, #SI : --> KGamma/A or microTesla/A 
            unit_phz:{'rad':np.power(10,0), 
                      'mrad': 1/np.power(10,3),
                      'deg': np.pi/180
                  },        # S.I : angle  on radians 
            unit_freq:{"hz":np.power(10,0), 
                         'rad/s': 0.1592,
                         'mrad/s':0.1592/np.power(10,3)} # on Hz 
            }
    #check units 
    for keys ,values in units.items():
        if keys.lower() == units_E_field.lower():
            if keys.lower() in values.keys():
                magn_E_field = magn_E_field * values[keys.lower()]
            else : 
                raise CSex.pyCSAMTError_Emag("Input units is not valid. try : {0}"\
                                                     " ".format(list(values.keys())))
        if keys.lower() ==units_H_field.lower() : 
            if keys.lower() in values.keys():
                magn_H_field= magn_H_field* values[keys.lower()]
            else : 
                raise CSex.pyCSAMTError_Emag("Input units is not valid. try : {0}"\
                                                 .format(list(values.keys()))) 
                    
        if keys.lower() ==unit_phz.lower(): 
            if keys.lower() in values.keys():
                
                phz_E_field = phz_E_field* values[keys.lower()]
                phz_H_field = phz_H_field  * values[keys.lower()]
            else : 
                raise CSex.pyCSAMTError_Emag("Input units is not valid. try : {0}"\
                                                 .format(list(values.keys()))) 
        if keys.lower() ==unit_freq.lower(): 
            if keys.lower() in values.keys():
                freq_value= freq_value * values[keys.lower()]
            else : 
                raise CSex.pyCSAMTError_Emag("Input units is not valid. try : {0}"\
                                                 .format(list(values.keys()))) 
        
    #--> component real part and imag 
    Zij =magn_E_field /magn_H_field
    phz=phz_E_field - phz_H_field
    Zreal, Zimag  =Zij * np.cos(phz), Zij * np.sin(phz)
    # Zreal_imag= np.complex(real=Zreal, imag=Zimag)
    Zreal_imag =np.complex(real=Zreal , imag=Zimag) 
    rho =(1/(5*freq_value))* Zij**2
    
    return rho , phz,  Zij, Zreal, Zimag, Zreal_imag


def z_error2r_phi_error(z_real, z_imag, error):
    """
    Error estimation from rectangular to polar coordinates.
    By standard error propagation, relative error in resistivity is 
    2*relative error in z amplitude. 
    Uncertainty in phase (in degrees) is computed by defining a circle around 
    the z vector in the complex plane. The uncertainty is the absolute angle
    between the vector to (x,y) and the vector between the origin and the
    tangent to the circle.
 
    :param z_real: real component of z (real number or array)
    :type z_real: float 
    
    :param z_imag:  imaginary component of z (real number or array)
    :type z_imag: complex
    
    :param error: absolute error in z (real number or array)
    :type error: float
    
    :returns: containers of relative error in resistivity, absolute error in phase
    :rtupe: tuple
    
    """
        
    z_amp = np.abs(z_real + 1j*z_imag)

    z_rel_err = error/z_amp
    
    res_rel_err = 2.*z_rel_err
    
    #if the relative error of the amplitude is >=100% that means that the relative 
    #error of the resistivity is 200% - that is then equivalent to an uncertainty 
    #in the phase angle of 90 degrees:
    if np.iterable(z_real):
        phi_err = np.degrees(np.arctan(z_rel_err))   
        phi_err[res_rel_err > 1.] = 90.
        
    else:
        if res_rel_err > 1.:
            phi_err = 90
        else:
            phi_err = np.degrees(np.arctan(z_rel_err))    
    
    
    return res_rel_err, phi_err

def rhophi2z ( phase, freq , resistivity=None,  z_array= None):
    """
    Function to compute z , real part and imag part . 

    :param phase: phase angles array in radians
    :type phase: ndarray 
     
    :param freq: frequencies array in Hz
    :type freq: array_like
    
    :param resistivity:  rho array in ohm.m
    :type resistivity: array_like
    
    :param z_array: impedance z array in V/m 
    :type z_array: array_like 

    :returns: z_abs , absolute value of zz
    :rtype: float 
    
    :returns: z_real, real part of complex number
    :rtype: float 
    
    :returns: z_imag, imaginary part of zz
    :rtype: complex 
    
    :returns: ndarray 
    :rtype: zz, array of z_abs, z_imag, z_real 
    """
    if (type(phase) is not np.ndarray) or (type(freq) is not np.ndarray) :
        phase, freq = np.array(phase), np.array(freq)
    
    if phase.dtype not in ['float', 'int'] or freq.dtype not in ['float', 'int']:
        try : 
            phase =np.array([float(ss) for ss in phase ])
            freq =np.array([float(f) for f in freq ])
        except : ValueError('must be flost number of integer . not str !')
    
    try :    
        if resistivity is not None : resistivity =np.array([float(res) for res in resistivity])
        if z_array is not None : z_array =np.array([float(z) for z in z_array])
    except : ValueError('Arguments number must be float or int')
    
    if freq.size != phase.size and (phase.size != resistivity.size or phase.size != z_array.size):
        raise 'Arrays nust get the same size.'
    
    
    if resistivity is not None : z_abs= np.sqrt(0.2 * freq * resistivity)
    if resistivity is  None and z_array is not None : z_abs = z_array 
    
    
    z_real, z_imag  = z_abs * np.cos (phase), z_abs * np.sin(phase)
    # zz =np.complex(z_real , z_imag)
    zz=z_real + z_imag *1j
    
    return z_abs, z_real, z_imag ,zz


def get_reffreq_index(freq_array, reffreq_value): 
    """ 
    Get the index of reference index. From this index ,All array will filter data at this reffreq
    value . 
    :param freq_array: array of frequency values
    :type freq_array: array_like 
    
    :param reffreq_value:  value of frequency at clean data 
    :type reffreq_value: float, int 
    """
    freq_array=np.array([float(ss) for ss in freq_array])
    if float(reffreq_value) not in freq_array : 
        raise CSex.pyCSAMTError_frequency('Reference frequency must be a value of frequency array.')

    for ii, freq in enumerate(freq_array): 
        if freq == reffreq_value:
            index_rf = ii
            break
    return index_rf

def interpolate_sets (array_to, fill_value =None , array_size=None ):
    """
    Function to interpolate data contain of multiple nan values. 
    """
    if array_size is None : array_size == array_to.size 
    
    xx_array_to = np.arange(array_size)
    
    ff = spi.interp1d(x=xx_array_to, y=array_to, fill_value=fill_value)
    # new_array = ff(xx_array_to)
    return ff(xx_array_to)

def  _interpolate_array_fromreffreq(stationNames , freq_array , reffreq_array,
                                    array_dict_loc, x_new=None, rigoureous=True,  order_axis=0): 
    """ 
    Interpolate value from  reference frequency : 
 
    :param stationNames:  list of stations of survey .
    :type stationNames: list
     
    :param freq_array:  array of frequencies 
    :type freq_array: (ndarray,1)
     
    :param reffreq_array: reference frequence for clean Data .
    :type reffreq_array: float 
     
    :param array_dict_loc: location of each value according each station names  
                            keys are stationsNmes and value of array
                            at that station , eg -S00 , valueof RhO  
    :type array_dict_loc: dict 
    """
    
    def _interpolate_array(array_to, kind='linear', fill_value='extrapolate' ):
        """ 
        Interpolate value to according the list of frequencies. 
        
        :param array_to: (ndarray,1), array to interpolate
        :type array_to: array_like
        
        :param fill_value: kind of extrapolation
        :type fill_value: str 
        """ 

        Y_intp_value =[]
        for ss, rowline in enumerate(array_to):
            func_interp = spi.interp1d(x_old, rowline, kind='linear', fill_value=fill_value)
            y_new=func_interp(x_new)

            Y_intp_value.append(y_new)
        new_Y_array =func.concat_array_from_list(list_of_array=Y_intp_value, concat_axis=order_axis)
        
        return new_Y_array
    
    if x_new is None : x_new=freq_array
    x_old =freq_array[:get_reffreq_index(freq_array=freq_array,
                                                    reffreq_value=reffreq_array)+1] # add the value of rffreq_freqq
    if rigoureous is True : 
        x_new = np.linspace(freq_array[0], freq_array[-1], x_old.size)

    temp_list =[values [:get_reffreq_index(freq_array=freq_array,
                                                     reffreq_value=reffreq_array)+1] for  keys , values in array_dict_loc.items() \
                                                       for  stn in stationNames  if keys ==stn  ]
    value_to_intp =func.concat_array_from_list(list_of_array=temp_list, concat_axis=order_axis) 
    return value_to_intp , _interpolate_array(array_to=value_to_intp)
       

def find_reference_frequency(freq_array =None, reffreq_value =None , sharp =False, etching=True): 
    """
    Method to find and interpolate reference value if it is not present on the frequency range. 

    :param freq_array: array_like frequency range 
    :type freq_array: array_like
    
    :param reffreq_value: reference frequency value
    :type reffreq_value: float or int 
    
    :param sharp:  if set to True , it forces the program to find mainly 
                    a value closest inside the  frequency range.
    :type sharp: bool  

    :param etching: bool , if set to True , it will print in your stdout.
    :type etching: bool 
    
    :returns: reference frequency 
    :rtype: float 
            
    """

    if freq_array is None or reffreq_value is None : 
        raise ValueError("None value can not be computed. check your frequency array or reference value.")
        
    if freq_array.dtype not in ['float', 'int']: 
        try : freq_array=np.array([float(ff) for ff in freq_array])
        except :raise TypeError('Frequency data must be on (ndarray,1) of float value not str.')
    if type(reffreq_value) is not int or type(reffreq_value) is not float: 
        try : reffreq_value= np.int(reffreq_value)
        except: raise CSex.pyCSAMTError_frequency("Reference frequency must be either float value or integer.")
        
    if freq_array.min() > reffreq_value > freq_array.max():
        warnings.warn('Reference frequency is out of frequency range. see more infos about reference '\
                                                  ' frequency.|{0}|'.format(infos.notion.reference_frequency))
        raise CSex.pyCSAMTError_frequency ("Input reference frequency is out the frequency range. "\
                                          "Please put value between the frequency range. "\
                                        "Frequency range is [{0} Hz to {1} Hz].".format(freq_array.min(), freq_array.max()))
    
    def force_interpolation (value_to_steep= None): 
        """ 
        Method to force reference value interpolated to find a value in frequency range close to . 

        :param value_to_steep:  reference frequency value to be forced .
        :type value_to_steep: float 
        
        :returns: closet value of interpolated frequency 
        :rtype: float
        """ 
        #find the value close to reference array.
        
        if freq_array[0] > freq_array[-1]: # check whether frequencies are range to higher frequency to lower :  
            temp_freq_array= freq_array[::-1] # sorted to lower to Highest 
        else : temp_freq_array = freq_array
        
        if value_to_steep is not None : 
            for xi in temp_freq_array : 
                if xi == value_to_steep : 
                    index_refreq= np.where(temp_freq_array==xi)[0] #find the index location on the normalization freq 
                    break 
                elif xi > value_to_steep : 
                    # check which frequency is closet one : 
                    index_temp= np.where(temp_freq_array==xi)[0] # find the upper location and take the lower value 
                                                       #assume to be the -1nieme element of the range.
                    # print(xi, index_temp)
                    dx1 = np.abs(temp_freq_array[index_temp-1]- xi )
                    dx2 = np.abs(temp_freq_array[index_temp] -xi)
                    if dx1 < dx2 : index_refreq =index_temp-1 
                    else : index_refreq=index_temp
                    break
        #new_value find  . and 
            new_val = temp_freq_array[index_refreq]
        return freq_array [np.where(freq_array ==new_val)]
        
        # return freq_array[index_refreq]
    
    interp_func =spi.interp1d(x=np.log10(freq_array), y=freq_array)
    new_reference_value = interp_func(np.log10(reffreq_value))

    if sharp ==True :
        if etching:
            print('---> Input reference frequency <{0}> Hz has been interpolated to < {1} > Hz.'.\
                   format(np.float(reffreq_value), np.around(float(new_reference_value ),2)))
        return force_interpolation(value_to_steep=new_reference_value  )
    elif sharp ==False :return new_reference_value 
    
def compute_TMA (data_array=None, number_of_TMApoints=None ):
    """
    function to compute a trimmed-moving-average filter to estimate average apparent resistivities.

    :param data_array:  content of value to be trimmed 
    :type data_array: array_like(ndarray,1) 
    
    :param number_of_TMA points:  number of filter points .
    :type number_of_TMA points: int 
    
    :returns:  value corrected with TMA  
    :rtype: array_like (ndarray, 1)
          
    """
    
    if data_array is None or number_of_TMApoints is None :raise CSex.pyCSAMTError_processing('NoneType arguments can not be computed!.')
    if data_array.dtype not in ['float', 'int']: 
        try :data_array=np.array([float(dd) for dd in data_array])
        except : raise CSex.pyCSAMTError_AvgData('Data must be on array_like float number.!')
    if type(number_of_TMApoints) is not int :
        try : number_of_TMApoints=np.int(number_of_TMApoints)
        except : raise CSex.pyCSAMTError_parameter_number('TMA filter point must be integer.')
    
    roll_TMA,window  =[], np.int(np.trunc(number_of_TMApoints/2)),

    for ii, value in enumerate(data_array): 
        if ii < window:
            get_TMA_value =np.delete(data_array[ii:number_of_TMApoints+ii], 
                                     np.array([data_array[ii:number_of_TMApoints+ii].argmin(), 
                                     data_array[ii:number_of_TMApoints+ii].argmax()])).mean()
            roll_TMA.append(get_TMA_value)
            
        elif ii == window or ii < data_array.size - window :
            if ii == data_array.size - (window+1) :
                get_TMA_value =np.delete(data_array[ii-window:],
                                         np.array([data_array[ii-window:].argmin(),
                                                                                                               data_array[ii-window:].argmax()])).mean() 
            else :get_TMA_value =np.delete(data_array[ii-window:ii+window+1],
                                           np.array([data_array[ii-window:ii+window+1].argmin(), 
                                           data_array[ii-window:ii+window+1].argmax()])).mean()
            roll_TMA.append(get_TMA_value)
            
        elif ii >=data_array.size - window:
            if value == data_array[-1]:
                get_TMA_value =np.delete(data_array[ii-number_of_TMApoints-1:],
                                         np.array([data_array[ii-number_of_TMApoints-1:].argmin(), 
                                         data_array[ii-number_of_TMApoints-1:].argmax()])).mean()
                
            else :get_TMA_value =np.delete(data_array[ii-number_of_TMApoints:],
                                           np.array([data_array[ii-number_of_TMApoints+1:ii+1].argmin(), 
                                           data_array[ii-number_of_TMApoints+1:ii+1].argmax()])).mean()  
            roll_TMA.append(get_TMA_value) 
            
    return np.array(roll_TMA)      
                                                              
def get_data_from_reference_frequency(array_loc, freq_array, reffreq_value):
    """
    Function to get reference frequency  without call especially stations array.
    The function is profitable but . It's less expensive However if something wrong happened
    by using the first step to get a reference array , it will try the traditionnally function 
    to get it. If none value is found , an Error will 
    occurs. 

    :param array_loc: assume to be a dictionnary of stations_data_values. 
    :type  array_loc: dict
    
    :param freq_array:  frequency array 
    :type  freq_array: array_like
    
    :param reffreq_value:  reffrence value, If the reference value is not in frequency array ,
                            function will force to interpolate value and find the correlative array.
                           
    :type reffreq_value: float or int 
            
    :returns:  an array of reference value at specific index 
    :rtype: array_like 
            
    :Example: 
        
        >>> get_data_from_reference_frequency(array_loc=rho,freq_array=freq_array,
        ...                                      reffreq_value=1023.)
        ... Input reference frequency has been interpolated to < 1024.0 > Hz 
    """
    
    if type(reffreq_value) is str : 
        try: reffreq_value=np.float(reffreq_value)
        except:raise CSex.pyCSAMTError_frequency("Reference value must be float or int value, not str!")
    #find the size of array loc : seem to be the size of stationsNames .
    stn_size =np.array(list(array_loc.keys())).size
    # # we sorted data for ensurity
    reffreq_array= np.array([values[np.where(freq_array==reffreq_value)[0]]
                             for keys, values in sorted(array_loc.items())])

    if reffreq_value not in freq_array : 
        reffreq_array= get_data_from_reference_frequency(array_loc=array_loc,
                                                         freq_array=freq_array,
                                                         reffreq_value=find_reference_frequency(
                                                                      freq_array =freq_array, 
                                                                      reffreq_value =reffreq_value ,
                                                                      sharp =True, 
                                                                      etching=True))
    if reffreq_array.size == stn_size :
        return reffreq_array.reshape(reffreq_array.shape[0],)
    else :

        from csamtpy.processing.callffunc import get_array_from_reffreq  as gfreq
        reffreq_array=gfreq(array_loc=array_loc , freq_array=freq_array, 
                            reffreq_value=reffreq_value,
                            stnNames=sorted(array_loc.keys()))
        if reffreq_array is None :
            raise CSex.pyCSAMTError_frequency('Something wrong happened during reference array computing.'\
                                             ' Please check the your inputs data size, and your station size.')
        return reffreq_array

def perforce_reference_freq(dataset, frequency_array=None ):
    """ 
    Function to get automatically the reference frequency. If user doesnt provide the value , 
    function will find automatically value . 
 
    :param data: array of avg DATA,  ndim>1 
    :type data: array_like 
    
    :param frequency_array: array of frequency
    :type frequency_array: array_like   
    
    :returns: reffreq_value float , reference frequency value 
    :rtype: float
    
    :returns: uncover_index,  index of reference value on frequency array 
    :rtype: int 
     
    :returns: nan_ratio , the ratio or the prevalence of nan in the data_set 
    :rtype: float 
    """
    import copy
    new_data_array =copy.deepcopy(dataset)
    if frequency_array is None :
        frequency_array=np.unique(new_data_array[:,1]) #find according to AVGData disposal .
        
    # ifall data are cleaned the take the highst freqeuncy 
    if np.isnan(new_data_array).sum()/new_data_array.size ==0: 
        if frequency_array[0] < frequency_array[-1] :
            reffreq_value = frequency_array[-1] 
        else : reffreq_value = frequency_array[0]
        print('--> Reference frequency is estimated to <{0}> Hz'.format(reffreq_value ))
        
        return  reffreq_value , frequency_array.size -1,0.
    
    elif np.isnan(new_data_array).sum()/new_data_array.size !=0 : 
        nan_ratio = np.isnan(new_data_array).sum()/new_data_array.size 
        new_data_array=np.nan_to_num(new_data_array, nan = -99999)
        # uncover_index = np.where(new_data_array==-99999)   
        reffreq_value = frequency_array[np.where(new_data_array==-99999)[0][0]]
        print('--> Reference frequency is estimated to <{0}> Hz'.format(reffreq_value))
        return reffreq_value, np.where(new_data_array==-99999)[0][0],  nan_ratio 

def hanning (dipole_length=50., number_of_points=7, large_band=False):
    """ 
    Function to compute hanning window .
    
    .. seealso::Torres-Verdin and Bostick, 1992, Principles of spatial surface electric
             field filtering in magnetotellurics: electromagnetic array profiling (EMAP),
             Geophysics, v57, p603-622.
        
    :param dipole_length: the length of dipole , xk is centered between dipole 
    :type dipole_length: float
    
    :param number_of_points: number of filter points 
    :type number_of_points: int  
    
    :returns: windowed hanning 
    :rtype: array_like 
    
    """
    
    if type(dipole_length) !=float or type(number_of_points) !=int: 
        
        try : dipole_length, number_of_points = np.float(dipole_length), np.int(np.floor(number_of_points))
        except :raise CSex.pyCSAMTError_processing('Dipole length and number of'\
                                           ' points must be a float number.')
        
    window_width = number_of_points * dipole_length
    xk , HAN = window_width / 2 , []                           # center point xk 
    X= np.linspace(-xk, xk, number_of_points)
    # X = np.arange(  -xk , xk, dipole_length)
    # xpos= np.arange(number_of_points+1)*xk 
    
    if large_band == True : return np.array([(1/window_width * (1+ np.cos(2* np.pi * xx / window_width))) for xx in X ])

    for xx in X : 
        if xx <= 0: hannx= 1/window_width * (1+ np.cos(2* np.pi * xx / window_width))
        else :hannx = 0 
        HAN.append(hannx)
    # import matplotlib.pyplot as plt 
    # plt.plot(X,np.array(HAN) )
    # plt.show()

    return np.array(HAN)

 
def weight_beta (dipole_length =50. , number_of_points = 7, window_width=None): 
    """
    WeightBeta function is  weight Hanning window . if window width is not provide  , function 
    will compute the width of window. 
    
    .. seealso:: 
            Torres-Verdin and Bostick, 1992, Principles of spatial surface electric field filtering in
            magnetotellurics: electromagnetic array profiling (EMAP), Geophysics, v57, p603-622.
            ...
            
    .. note::SUM(Betaj (j=1..M)=1 .
              
    :param dipole_length: length of dipole in meter (m) 
    :type dipole_length: float
    
    :param number_of_points: number of station points to filter
    :type  number_of_points: int 
    
    :param window_width: the width of window filter 
    :type window_width: float 
      
    :returns: beta_array at each station 
    :rtype: array_like 
    """
    _logger.info( 'Computing weight Hanning window ')
    # beta_infx = - dipole_length/2 
    # beta_posfx = dipole_length /2 
    if window_width is None : window_width = number_of_points *  dipole_length
    if window_width is not None : 
        try : number_of_points = np.int(window_width / dipole_length ) 
        except : 
            warnings.warn('Please see more info about coeff.Beta :<{0}>'.format(infos.notion.weighted_Beta))

            raise CSex.pyCSAMTError_processing('Window width type must be float number.')
    
    xk = window_width / 2 
    num_pk =np.arange(dipole_length/2, window_width , dipole_length)
    beta_range= []
    for xp in num_pk : 
        if np.abs(xp - xk) <= window_width /2 : 
            beta_xp = ((1/window_width) * ((xp + dipole_length/2 - xk) +(window_width/ 2* np.pi)*\
                       (np.sin(2* np.pi * (xp + dipole_length/2 -xk )/window_width))))-\
                ((1/window_width)* ((xp - dipole_length/2 - xk) +(window_width/ 2* np.pi)*\
                 (np.sin(2* np.pi * (xp - dipole_length/2 -xk)/window_width))))
                    
        elif np.abs(xp - xk) > window_width/2 : beta_xp = 0.
        
        beta_range.append(beta_xp) 
    
    return np.array(beta_range)
    # for xj in enra

def compute_AMA ( z_array =None  , weighted_window=None ,
                  dipole_length = None , number_of_points =None ):
    """
    .. note:: 
        We will add this filter later !
    ***future plan ***
    """
     
    if z_array is None : raise CSex.pyCSAMTError_Z('NoneType Impedance can not be computed.')
    if weighted_window is None : 
        if dipole_length is None or number_of_points is None :
            raise CSex.pyCSAMTError_inputarguments('Please check '
                'your Inputs values ! NoneType values of dipole length or number_of_stations'
                ' can not be computed.!')
        try : 
            dipole_length , number_of_points  = np.float(dipole_length) , np.int(number_of_points)
        except : raise CSex.pyCSAMTError_float('Dipole length, number of points must be a float, and '\
                                               ' int  number respectively.!')

        weighted_window =weight_beta(dipole_length=dipole_length, 
                                      number_of_points=number_of_points)

    xc_pk , cFLAM , cut_off = np.int(np.floor(weighted_window.size / 2)), [], -1 

    for kk , zz in enumerate(z_array): 
        if kk < xc_pk : 
            cut_off += 1 
            zz_new = weighted_window[xc_pk - cut_off : weighted_window.size] *\
            z_array[kk - cut_off: weighted_window[xc_pk-cut_off : weighted_window.size].size]
            cFLAM.append(zz_new.sum())
        elif kk == xc_pk or kk <= z_array.size - weighted_window[xc_pk:].size: 
            zz_new, cut_off= weighted_window * z_array[kk- xc_pk : kk + weighted_window[xc_pk :\
                                                                weighted_window.size].size], 0
            cFLAM.append(zz_new.sum())
        elif kk > z_array.size - weighted_window[xc_pk:].size: 
            cut_off -= 1
            zz_new = z_array[kk-xc_pk :] * weighted_window[:(weighted_window.size +cut_off)]
            cFLAM.append(zz_new.sum())
            
    return np.array(cFLAM)
            
def hanning_xk (dipole_length =50. , number_of_points =7):
    """
    compute _hanning window on a wtdth of number of point : 
        integrate value on all the window_bandwidth discrete and continue.
        if value is greater than Hald of the width  value == 0 . 
       
    :param dipole_length: length of dipole 
    :type dipole_length: float 
    
    :param number_of_points:  value of points or survey stations . 
    :type number_of_points: int 
    
    :returns: han_xk , continue value on half bandwidth x0-- xk (center point) 
    :rtype: array_like 
    
    :returns: windowed hanning, discrete _value ,SUM(han(x0, xk))
    :rtype: array_like 
    
    """
    
    window_width = number_of_points *  dipole_length 
    xk , hannx= window_width / 2 , []
    xpk_num =np.arange(dipole_length/2, window_width, dipole_length)

    han_xk = np.array((1/ window_width * ( xk  + (window_width/2*np.pi ) * np.sin(2 * np.pi * (xk) /window_width))) - \
                (1/ window_width * ( xpk_num[0]  + (window_width/2*np.pi ) * np.sin(2 * np.pi * (xpk_num[0]) /window_width))))
    
    # compute Hanning discrete : 
    for  ii, xx in enumerate(xpk_num) : 
        if np.abs(xx) <= window_width /2 :
            han_xx = 1/ window_width * ( (xx + dipole_length/2)  + (window_width/2*np.pi  * np.sin(2 * np.pi * (xx + dipole_length/2) /window_width))) - \
                1/ window_width * ( (xx - dipole_length/2)  + (window_width/2*np.pi  * np.sin(2 * np.pi * (xx - dipole_length/2) /window_width)))
        else : han_xx= 0 
        
        hannx.append(han_xx)

    
    return han_xk , np.array(hannx)
    
def hanning_x (x_point_value, dipole_length =50.,  number_of_points=7, 
               bandwidth_values=False , on_half=True):
    """
    Function to compute point on window width .  Use discrete computing . Function show the value
    at center point of window assume that the point is center locate on the window width . 
    It intergrates  value between dipole length. User can use see_extraband to see the values 
    on the total bandwith. If half is False the value of greater than center point will be
    computed and not be 0 as the normal definition of Hanning window filter. 
    
    :param x_point_value:  value  to intergrate.
    :type x_point_value: float 
    
    :param dipole_length: length of dipole on survey
    :type dipole_length: float 
    
    :param number_of_point: survey point or number point to apply.
    :type number_of_point: int 
    
    :param bandwidth_values: see all value on the bandwith , value greater than x_center
                            point will be computed .
    :type bandwidth_values:bool 
                            
    :param on_half:  value on the bandwith; value greater that x_center point = 0.
    :type on_half: bool 

    :returns: hannx  integrated X_point_value or  array of window bandwidth .
    :rtype: array_like 
    """
    window_width = number_of_points * dipole_length 
    
    xx_points_array = np.arange(dipole_length/2 , window_width , dipole_length )
    
    def see_extraband( apply_on_half =False): 
        """ 
        Function to see data of the large bande integrated hanning. if apply_on_half 
        is set to True < it will show the data on band <= half windowith i.e. thin 
        center point xk.
        
        :param apply_on_half: half band of hanning filter.
        :type apply_on_half: bool  

        :returns: windowed hanning band
        :rtype: array_like(ndarray,1) 
        """
        hv=[]
        
        if not apply_on_half:
            for xv  in  xx_points_array:
                hann_xv= np.array( 1/ window_width * ( (xv + dipole_length/2)  + (window_width/2*np.pi  * np.sin(2 * np.pi * (xv+ dipole_length/2) /window_width))) - \
                        1/ window_width * ( (xv - dipole_length/2)  + (window_width/2*np.pi  * np.sin(2 * np.pi * (xv - dipole_length/2) /window_width))))
                hv.append(hann_xv)
        elif apply_on_half : 
            for xv  in  xx_points_array:
                if xv  <= window_width/2: 
                    hann_xv= np.array( 1/ window_width * ( (xv + dipole_length/2)  + (window_width/2*np.pi  * np.sin(2 * np.pi * (xv+ dipole_length/2) /window_width))) - \
                        1/ window_width * ( (xv - dipole_length/2)  + (window_width/2*np.pi  * np.sin(2 * np.pi * (xv - dipole_length/2) /window_width))))
                else :hann_xv=0 
                hv.append(hann_xv)
                
        return np.array(hv)
        
    if (x_point_value - dipole_length/2 < 0 or x_point_value +dipole_length/2 > window_width):
        raise CSex.pyCSAMTError_processing('Hanning x point value is outside the window.')
    
    if bandwidth_values is True : 
        if on_half :return see_extraband(apply_on_half=True)
        if not on_half :return  see_extraband(apply_on_half=False)
    else :
        if x_point_value <= window_width/2:
            hann_xp= np.array( 1/ window_width * ( (x_point_value + dipole_length/2)  + (window_width/2*np.pi  * np.sin(2 * np.pi * (x_point_value + dipole_length/2) /window_width))) - \
                        1/ window_width * ( (x_point_value - dipole_length/2)  + (window_width/2*np.pi  * np.sin(2 * np.pi * (x_point_value - dipole_length/2) /window_width))))
            
        else : hann_xp = 0
        
    return  hann_xp

def wbetaX2 (Xpos, dipole_length=50., number_of_points=5) :
    """
    weight Beta is  computed following the paper of 
    Torres-verdfn, C., and F. X. Bostick, 1992, Principles of spatial surface electric field 
    filtering in magnetotellurics - Electromagnetic array profiling ( EMAP )- Geophysics, 57(4), 25–34. 
        
    :param Xpos: reference position on the field 
    :type Xpos: str 
    
    :param dipole_length: length of dipole measurement
    :type dipole_length: float 
    
    :param number_of_points:  point to stand filters , window width 
    :type number_of_points: int 
    
    """
    
    window_width = dipole_length * number_of_points 
    xk = window_width/2
    

    # if np.abs(xp - xk) <= window_width /2 : 
    beta_xpos = ((1/window_width) * ((Xpos + dipole_length/2 - xk) +(window_width/ 2* np.pi)*\
               (np.sin(2* np.pi * (Xpos + dipole_length/2 -xk )/window_width))))-\
        ((1/window_width)* ((Xpos - dipole_length/2 - xk) +(window_width/ 2* np.pi)*\
         (np.sin(2* np.pi * (Xpos - dipole_length/2 -xk)/window_width))))
            
        # elif np.abs(xp - xk) > window_width/2 : beta_xpos = 0.
        
        # beta_range.append(beta_xpos) 
    return np.array([Xpos-dipole_length/2, Xpos + dipole_length/2]) , beta_xpos 

def compute_sigmas_e_h_and_sigma_rho(pc_emag , pc_hmag, pc_app_rho, app_rho, emag, hmag):
    """
    function to compute Standard Deviation for E-field (sigma_e), 
    Standard Deviation for H-Field (sigma_h) , 
    & Standard Deviation for Component RHO (sigma_rho)

    Parameters
    -----------
        * pc_emag : float  
                Statistical variation of magnitude values from averaged data blocks.
                Standard Deviation/Average Emag (%)
        * pc_hmag :float 
                Statistical variation of magnitude values from averaged data blocks.
                Standard Deviation / Average Hmag (%)
        * pc_app_rho: float 
                Statistical variation of magnitude values from averaged data blocks.
                Standard Deviation / Average Rho (%)
                 
        * app_rho :float 
                resistivity calculated from averaged component (ohm.m)
        * Emag : float
                average E - field magnitude(microVolt/Km *amp )
        * Hmag : float 
                average H - field magnitude(pTesta/amp) or (milliGammas/Amp)
        
    Returns
    --------
        sigma_rho : float 
             srhoC (Standard Deviation for Component RHO)
        c_var_Rho : float 
            C-varrhoC(  Coefficient of Variation for Component RHO)
    """
    # Standard Deviation for E-field (se) in (%) 1micoV/m*amp= 1mV/Km*amp
    emag /=1e3 # convert microvolt to millivolt 
    # sigma_e, sigma_h, sigma_rho   = emag *(pc_emag/100), hmag * (pc_hmag /100) ,\
    #     pc_app_rho/100 * app_rho
    
    sigma_e , sigma_h, sigma_rho   =  pc_emag/100, pc_hmag /100 , pc_app_rho/100

   
    return  sigma_e , sigma_h, sigma_rho 
    
def compute_weight_factor_betaj(Xpos , dipole_length=50. , number_of_points =5. ):
    """
    weight Beta is  computed following the paper of Torres-verdfn, C., and F. X. Bostick, 1992,
     Principles of spatial surface electric field filtering in magnetotellurics - Electromagnetic 
    array profiling ( EMAP )- Geophysics, 57(4), 25–34. 
        
    :param Xpos: reference position on the field 
    :type Xpos: str 
    
    :param dipole_length: length of dipole measurement
    :type dipole_length: float 
    
    :param number_of_points:  point to stand filters , windowed width 
    :type number_of_points: int 
    
    """
    #figureout the center ot dipole_lenght  and  compute the window width 
    
    window_width = dipole_length * number_of_points  #+  dipole_length # for extrema 
    
    xk = window_width/2
    
    # declare limit of integrating hanning 
    
    xsup_lim = Xpos + dipole_length/2 
    xinf_lim = Xpos - dipole_length /2 
    
    X =np.arange(0, window_width, dipole_length) 

    betapos = (1/ window_width) * ((xsup_lim - xk + (window_width /(2 * np.pi)) *\
                                 np.sin( 2* np.pi *(xsup_lim - xk)/ window_width )) - (
                                xinf_lim - xk + window_width /(2 * np.pi) *\
                                 np.sin( 2* np.pi *(xinf_lim - xk)/ window_width ))) 
            
    return betapos, window_width, X 

def compute_FLMA ( z_array =None  , weighted_window=None ,
                  dipole_length = 50. , number_of_points =5. ):
    """
    Use a fixed-length-moving-average filter to estimate average apparent
    resistivities at a single static-correction-reference frequency
    
    :param z_array:  array of uncorrected impedance values.
    :type z_array: arry_like, complex
    
    :param dipole_length: length of dipole on survey
    :type dipole_length: float 
    
    :param number_of_points: survey point or number point to apply.
    :type number_of_points: int 
    
    :param weighted_window: bandwidth of hanning window,
                            if window is None , will computed automatically once 
                            wwith number of points or dipoles and dipole length
    :type weighted_window: float 
 
    :returns: windowed hanning , filtered window.
    :rtype: array_like 
    
    """
     
    if z_array is None : 
        raise CSex.pyCSAMTError_Z('NoneType Impedance can not be computed.')
    if weighted_window is None : 
        if dipole_length is None or number_of_points is None :
            raise CSex.pyCSAMTError_inputarguments('Please check '
                'your Inputs values ! NoneType values of dipole length '
                ' or number_of_stations can not be computed.!')
        try : 
            dipole_length = np.float(dipole_length) 
            number_of_points= np.int(number_of_points)
            
        except :
            raise CSex.pyCSAMTError_float(
                'Dipole length, number of points must be a float, and '
                    ' int  number respectively.!')

        weighted_window = np.arange(0, dipole_length * number_of_points , 
                                    dipole_length)
        
    # compute weight factor Bj 
    wk = np.array([compute_weight_factor_betaj(Xpos , dipole_length=dipole_length , 
                                number_of_points =number_of_points)[0] 
          for Xpos in weighted_window])
    
    #seek the center point  index xc_pk  
    
    if wk.size %2 ==0 : #mean number is even 
        xc_pk = np.int(wk.size /2) -1  # backward for one index 
    elif wk.size %2 !=0 :
       xc_pk =  np.int(wk.size /2)   # weight array size is odd 
       
    
    cFLMA , cut_off = [], -1  # take a minimum 

    import copy 
    tem_wk = copy.deepcopy(wk)
    
    for kk , zz in enumerate(z_array): 
        if kk < xc_pk : 
            cut_off += 1  # cut off to shrink one index  
            tem_wk [:xc_pk -cut_off] = 0. # where wk 0 outside of the Hanning window
            tem_z= z_array[:kk + xc_pk +1] # keep main array 
            # buid array with number of 0 outside the window 
            add0 = len(wk)- len(tem_z)  
            # build array as the same dim with hanning window centered a k position
            tem_z = np.concatenate((np.repeat(0., add0), tem_z))
            zz_new= tem_z * tem_wk # compute value 
            cFLMA.append(zz_new.sum())
            
            tem_wk =copy.deepcopy(wk) # copy once again for next loop if exist 
            
   
        elif  xc_pk <= kk < z_array.size - xc_pk: 
            if xc_pk %2 ==0 :
                zz_new= z_array[kk - xc_pk: kk + xc_pk +1  ] * wk # when center point is even 
            else : # when odd 
                zz_new= z_array[kk - xc_pk  : kk + xc_pk + 1  ] * wk
                
            cut_off =0
            cFLMA.append(zz_new.sum())
   
        elif kk >= z_array.size - xc_pk:
       
            cut_off -= 1
            tem_z = z_array[kk-xc_pk:] 
            add0 = len(wk) - len(tem_z)  
            tem_wk = np.concatenate((tem_z, np.repeat(0.,add0)))

            zz_new =tem_wk * wk
            
            cFLMA.append(zz_new.sum())
    
            
    return np.array(cFLMA)    

    
if __name__=='__main__': 
    
    z1= np.array([2.46073791 +3.00162006j ,
         9.74193019 +1.82209497j,
         15.68879141 + 12.91164264j ,
         5.84384925 +3.6899018j, 
         2.4430065  +0.57175607j])
  
    X= np.arange(0, 250, 50)
    wk = np.array([compute_weight_factor_betaj(va, dipole_length=50, 
                                                number_of_points=5.)[0] for va in X])

 
    
    