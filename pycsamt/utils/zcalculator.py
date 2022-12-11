# -*- coding: utf-8 -*-
#       Created on Thu Dec  3 16:44:29 2020
#       Author: Kouadio K.Laurent<etanoyau@gmail.com>
#       Licence: LGPL

"""
.. _module-ZCalculator::`pycsamt.utils.zcalculator`
    :synopsis:  helper functions special calculator 
"""
# import os
import warnings 
import numpy as np 

import pycsamt.utils._p as infos
from pycsamt.utils import exceptions as CSex
from pycsamt.utils import func_utils as func
from pycsamt.utils.decorator import (deprecated, deprecated_to)
from pycsamt._csamtpylog import csamtpylog 

_logger =csamtpylog.get_csamtpy_logger(__name__)
try:
    import scipy
    scipy_version = [int(ss) for ss in scipy.__version__.split('.')]
    if scipy_version[0] == 0:
        if scipy_version[1] < 14:
            warnings.warn(func._msg, ImportWarning)
            _logger.warning(func._msg)
    import scipy.interpolate as spi

    interp_import = True

except ImportError: 
    warnings.warn(func._msg0)
    _logger.warning(func._msg0)
    interp_import = False

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
        raise CSex.pyCSAMTError_inputarguments(
            "Input first arguments is wrong. must"
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
        raise CSex.pyCSAMTError_inputarguments(
            "input argument must be a number in array or list.")
    
    if to_degree :
        return (phz_array.mean()/1e3) *180/np.pi, (phz_array/1000) *180/np.pi
    
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
        raise CSex.pyCSAMTError_inputarguments(
            'Error inputs arguments . must be list or ndarray of numbers.')
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
        raise CSex.pyCSAMTError_inputarguments(
            'Error inputs arguments. must be list or ndarray of numbers.')
    if to_degree : 
        return (pphz_array.mean()/1000) *180/np.pi, (pphz_array/1000) *180/np.pi
    
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
                raise CSex.pyCSAMTError_inputarguments(
                    "Elemts composed of each array must be a number.")
        
    # if comp_rho.__code__.co_argcount ==5 : 
    #     if mag_E_field.ndim !=1 or  mag_H_field.ndim !=1 or freq_array.ndim !=1 : 
    #         raise CSex.pyCSAMTError_inputarguments 
    #       ("dimension of Input arguments must me equal to 1.")
    
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
        >>> from pycsamt.core import avg  
        >>> phs_obj =avg.Phase(path)
        >>> phs_obj.loc['S00']
        >>> value, ss = comp_phz(comphz_array=phs_obj.loc[
        ...            'S00'], to_degree=True)
        ... print(value)
    """
    
    if type(comphz_array) is list : comphz_array=np.array(comphz_array)
    
    if comphz_array.dtype not in ["float", "int"] : 
        raise CSex.pyCSAMTError_inputarguments(
            "Type provided is wrong ! number must be float")  

    try : 
        comphz_array=np.array([float(kk) for kk in comphz_array ])
    except : 
        raise CSex.pyCSAMTError_inputarguments(
            'Error inputs arguments. must be list or ndarray of numbers.')
    
            
    if units =='deg' : 
        return comphz_array.mean() *180/np.pi, comphz_array *180/np.pi
    
    return  comphz_array.mean()


def compute_components_Z_Phz(magn_E_field , magn_H_field, phz_E_field,
               phz_H_field, freq_value, **kwargs):
    """
    Function to compute all components  derived from Impedance Z. 
    user can  enter specifik units in kwargs arguments . 
    program will compute and converts value 
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
        
        >>> from pycsamt.utils import zcalculator as Zcc
        >>> from pycsamt.core import avg 
        >>> path =  os.path.join(os.environ["pyCSAMT"], 
        ...              data', 'avg', 'K1.AVG')
        >>> emag_ob = avg.Emag(path)
        >>> hmag_obj = avg.Hmag(path)
        >>> ephz_obj = avg.Ephz(path)
        >>> hphz_obj = avg.Hphz(path)
        >>> freq_obj =avg.Frequency(path)
        >>> station_name ='S00'
        >>> rho, phz, Z, real, imag, comp =Zcc.compute_components_Z_Phz( 
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
                raise CSex.pyCSAMTError_Emag(
                    "Input units is not valid. try : {0}".format(
                        list(values.keys())))
        if keys.lower() ==units_H_field.lower() : 
            if keys.lower() in values.keys():
                magn_H_field= magn_H_field* values[keys.lower()]
            else : 
                raise CSex.pyCSAMTError_Emag(
                    "Input units is not valid. try : {0}".format(
                        list(values.keys()))) 
                    
        if keys.lower() ==unit_phz.lower(): 
            if keys.lower() in values.keys():
                
                phz_E_field = phz_E_field* values[keys.lower()]
                phz_H_field = phz_H_field  * values[keys.lower()]
            else : 
                raise CSex.pyCSAMTError_Emag(
                    "Input units is not valid. try : {0}".format(
                        list(values.keys()))) 
        if keys.lower() ==unit_freq.lower(): 
            if keys.lower() in values.keys():
                freq_value= freq_value * values[keys.lower()]
            else : 
                raise CSex.pyCSAMTError_Emag(
                    "Input units is not valid. try : {0}".format(
                        list(values.keys()))) 
        
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
        if resistivity is not None : 
            resistivity =np.array([float(res) for res in resistivity])
        if z_array is not None : 
            z_array =np.array([float(z) for z in z_array])
    except : ValueError('Arguments number must be float or int')
    
    if freq.size != phase.size and (phase.size != resistivity.size or\
                                    phase.size != z_array.size):
        raise 'Arrays nust get the same size.'
    
    
    if resistivity is not None : z_abs= np.sqrt(0.2 * freq * resistivity)
    if resistivity is  None and z_array is not None : z_abs = z_array 
    
    
    z_real, z_imag  = z_abs * np.cos (phase), z_abs * np.sin(phase)
    # zz =np.complex(z_real , z_imag)
    zz=z_real + z_imag *1j
    
    return z_abs, z_real, z_imag ,zz


def get_reffreq_index(freq_array, reffreq_value): 
    """ 
    Get the index of reference index. From this index, 
    All array will filter data at this reffreq
    value . 
    :param freq_array: array of frequency values
    :type freq_array: array_like 
    
    :param reffreq_value:  value of frequency at clean data 
    :type reffreq_value: float, int 
    """
    freq_array=np.array([float(ss) for ss in freq_array])
    if float(reffreq_value) not in freq_array : 
        raise CSex.pyCSAMTError_frequency(
            'Reference frequency must be a value of frequency array.')

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
                                    array_dict_loc, x_new=None,
                                    rigoureous=True,  order_axis=0): 
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
            func_interp = spi.interp1d(x_old, rowline, kind='linear',
                                       fill_value=fill_value)
            y_new=func_interp(x_new)

            Y_intp_value.append(y_new)
        new_Y_array =func.concat_array_from_list(list_of_array=Y_intp_value,
                                                 concat_axis=order_axis)
        
        return new_Y_array
    
    if x_new is None : x_new=freq_array
    # add the value of rffreq_freqq
    x_old =freq_array[:get_reffreq_index(
        freq_array=freq_array,reffreq_value=reffreq_array)+1] 
    if rigoureous is True : 
        x_new = np.linspace(freq_array[0], freq_array[-1], x_old.size)

    temp_list =[values [:get_reffreq_index(
        freq_array=freq_array,reffreq_value=reffreq_array)+1] 
        for  keys , values in array_dict_loc.items() for  stn in stationNames 
        if keys ==stn  ]
    value_to_intp =func.concat_array_from_list(
        list_of_array=temp_list, concat_axis=order_axis) 
    return value_to_intp , _interpolate_array(array_to=value_to_intp)
       

def find_reference_frequency(freq_array =None, reffreq_value =None ,
                             sharp =False, etching=True): 
    """
    Method to find and interpolate reference value if it is not present
    on the frequency range. 

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
        raise ValueError("None value can not be computed. check "
                         "your frequency array or reference value.")
        
    if freq_array.dtype not in ['float', 'int']: 
        try : freq_array=np.array([float(ff) for ff in freq_array])
        except :
            raise TypeError('Frequency data must be'
                            ' on (ndarray,1) of float value not str.')
    if type(reffreq_value) is not int or type(reffreq_value) is not float: 
        try : reffreq_value= np.int(reffreq_value)
        except: 
            raise CSex.pyCSAMTError_frequency(
                "Reference frequency must be either float value or integer.")
        
    if freq_array.min() > reffreq_value > freq_array.max():
        warnings.warn('Reference frequency is out of frequency '
                      'range. see more infos about reference'
                      ' frequency.|{0}|'.format(
                          infos.notion.reference_frequency))
        raise CSex.pyCSAMTError_frequency (
            "Input reference frequency is out the frequency range. "
             "Please put value between the frequency range. "
             "Frequency range is [{0} Hz to {1} Hz].".format(
                 freq_array.min(), freq_array.max()))
    
    def force_interpolation (value_to_steep= None): 
        """ 
        Method to force reference value interpolated to find a value 
        in frequency range close to . 

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
            print('---> Input reference frequency <{0}> '
                  'Hz has been interpolated to < {1} > Hz.'.
                   format(np.float(reffreq_value), np.around(
                       float(new_reference_value ),2)))
        return force_interpolation(value_to_steep=new_reference_value  )
    elif sharp ==False :return new_reference_value 
                                                           
def get_data_from_reference_frequency(array_loc, freq_array, reffreq_value):
    """
    Function to get reference frequency  without call especially 
    stations array. The function is profitable but  It's less expensive 
    However if something wrong happened by using the first step to get a 
    reference array , it will try the traditionnally function 
    to get it. If none value is found , an Error will 
    occurs. 

    :param array_loc: assume to be a dictionnary of stations_data_values. 
    :type  array_loc: dict
    
    :param freq_array:  frequency array 
    :type  freq_array: array_like
    
    :param reffreq_value:  reffrence value, If the reference value is not 
            in frequency array function will force to interpolate value 
            and find the correlative array.
                           
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
        except:
            raise CSex.pyCSAMTError_frequency(
                "Reference value must be float or int value, not str!")
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

        from pycsamt.processing.callffunc import get_array_from_reffreq  as gfreq
        reffreq_array=gfreq(array_loc=array_loc , freq_array=freq_array, 
                            reffreq_value=reffreq_value,
                            stnNames=sorted(array_loc.keys()))
        if reffreq_array is None :
            raise CSex.pyCSAMTError_frequency(
                'Something wrong happened during reference array computing.'
                 ' Please check the your inputs data size, and your station size.')
        return reffreq_array

def perforce_reference_freq(dataset, frequency_array=None ):
    """ 
    Function to get automatically the reference frequency. If user
    does not provide the value, the function will find automatically value . 
 
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
        frequency_array,*_=np.unique(new_data_array[:,1], return_counts=True) #find according to AVGData disposal .
        
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
        
        try : 
            dipole_length = np.float(dipole_length)
            number_of_points = np.int(np.floor(number_of_points))
        except :
            raise CSex.pyCSAMTError_processing('Dipole length and number of'\
                                           ' points must be a float number.')
        
    window_width = number_of_points * dipole_length
    xk , HAN = window_width / 2 , []                           # center point xk 
    X= np.linspace(-xk, xk, number_of_points)
    # X = np.arange(  -xk , xk, dipole_length)
    # xpos= np.arange(number_of_points+1)*xk 
    
    if large_band == True : return np.array([(1/window_width * \
                                              (1+ np.cos(2* np.pi * xx / window_width))) for xx in X ])

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
    WeightBeta function is  weight Hanning window . 
    if window width is not provide  , function 
    will compute the width of window. 
    
    .. seealso:: 
            Torres-Verdin and Bostick, 1992, Principles of spatial 
            surface electric field filtering in
            magnetotellurics: electromagnetic array profiling (EMAP),
            Geophysics, v57, p603-622.
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
            warnings.warn('Please see more '
                          'info about coeff.Beta :<{0}>'.format(infos.notion.weighted_Beta))

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
    
@deprecated('Function deprecated') 
def compute_adaptative_moving_average ( z_array =None  , weighted_window=None ,
                  dipole_length = None , number_of_points =None ):
    """
    .. note:: 
        see `compute_AMA` !
    """
     
    if z_array is None :
        raise CSex.pyCSAMTError_Z('NoneType Impedance can not be computed.')
    if weighted_window is None : 
        if dipole_length is None or number_of_points is None :
            raise CSex.pyCSAMTError_inputarguments('Please check '
                'your Inputs values ! NoneType values of dipole length '
                ' or number_of_stations can not be computed.!')
        try : 
            dipole_length , number_of_points  = np.float(dipole_length) , np.int(number_of_points)
        except : 
            raise CSex.pyCSAMTError_float(
                'Dipole length, number of points must be a float, and '
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
    Function to compute point on window width .  Use discrete computing .
    Function show the value at center point of window assume that the point
    is center locate on the window width .  It intergrates  value between 
    dipole length. User can use see_extraband to see the values 
    on the total bandwith. If half is False the value of greater than center
    point will be computed and not be 0 as the normal definition 
    of Hanning window filter. 
    
    :param x_point_value:  value  to intergrate.
    :type x_point_value: float 
    
    :param dipole_length: length of dipole on survey
    :type dipole_length: float 
    
    :param number_of_point: survey point or number point to apply.
    :type number_of_point: int 
    
    :param bandwidth_values: see all value on the bandwith, value greater 
                            than x_center point will be computed .
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
    Torres-verdfn, C., and F. X. Bostick, 1992, Principles of spatial 
    surface electric field filtering in magnetotellurics
     - Electromagnetic array profiling ( EMAP )- Geophysics, 57(4), 25–34. 
        
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
    
def compute_weight_factor_betaj(Xpos , dipole_length=50. , number_of_points =5. , 
                                window_width=None):
    """
    weight Beta is  computed following the paper of Torres-verdfn,
    C., and F. X. Bostick, 1992,Principles of spatial surface electric field
      filtering in magnetotellurics - Electromagnetic 
    array profiling ( EMAP )- Geophysics, 57(4), 25–34. 
        
    :param Xpos: reference position on the field 
    :type Xpos: str 
    
    :param dipole_length: length of dipole measurement
    :type dipole_length: float 
    
    :param number_of_points:  point to stand filters , windowed width 
    :type number_of_points: int 
    
    """
    #figureout the center ot dipole_lenght  and  compute the window width 
    if window_width is None : 
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
                  dipole_length = 50. , number_of_points =5., **kwargs ):
    """
    Use a fixed-length-moving-average filter to estimate average apparent
    resistivities at a single static-correction-reference frequency
    
    :param z_array:  array of uncorrected impedance values.
    :type z_array: arry_like, complex
    
    :param dipole_length: length of dipole on survey
    :type dipole_length: float 
    
    :param number_of_points: survey point or number point to apply.
    :type number_of_points: int 
    
    :param weighted_window: bandwidth of hanning window filter lengths,
                           skin_depth arrays at each station.
    :type weighted_window: array_like 
 
    :param app_rho: array of apparent resistivities in ohm.m , if not provided 
                    will use impedance zrho to compute app_rho with 
                    reference_frequency 
    :type app_rho: array_like 
    
    :returns: windowed hanning , filtered window.
    :rtype: array_like 
    
    :Example:
        
        >>> from pycsamt.utils import zcalculator as Zcc
        >>> z1= np.array([2.46073791 +3.00162006j ,
        ...             9.74193019 +1.82209497j,
        ...             15.68879141 + 12.91164264j ,
        ...             5.84384925 +3.6899018j, 
        ...             2.4430065  +0.57175607j])
        >>> flma= Zcc.compute_FLMA(z_array=z1, dipole_length=50. ,
                               number_of_points=4)
        >>> print(flma)
    """
    #reference_freq= kwargs.pop('reference_freq', 8192.)
    
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
            
        elif xc_pk <= kk:
            if len(wk) %2 !=0: # if window weight len is odd 
                if  kk < z_array.size - xc_pk:
                    zz_new= z_array[kk - xc_pk: kk + xc_pk +1  ] * wk # when center point is even 
                elif  kk >= z_array.size - xc_pk:
                    tem_z = z_array[kk-xc_pk:] 
                    add0 = len(wk) - len(tem_z)  
                    tem_wk = np.concatenate((tem_z, np.repeat(0.,add0)))
                    zz_new =tem_wk * wk

            elif len(wk) %2 ==0:# when window length is even.
                if kk <= z_array.size - xc_pk:
                    zz_new= z_array[kk - xc_pk  : kk + xc_pk] * wk
      
                elif  kk > z_array.size - xc_pk:
                    tem_z = z_array[kk-xc_pk:] 
                    add0 = len(wk) - len(tem_z)  
                    tem_wk = np.concatenate((tem_z, np.repeat(0.,add0)))
                    zz_new =tem_wk * wk    
                    
            cFLMA.append(zz_new.sum())    
            
    print('** {0:<27} {1} {2} m.'.format("Filter's width",
                                '=', int(number_of_points)* dipole_length))   
             
    return np.array(cFLMA)    

def compute_TMA (data_array=None, number_of_TMApoints=5. ):
    """
    function to compute a trimmed-moving-average filter 
    to estimate average apparent resistivities.

    :param data_array:  content of value to be trimmed 
    :type data_array: array_like(ndarray,1) 
    
    :param number_of_TMA points:  number of filter points .
    :type number_of_TMA points: int 
    
    :returns:  value corrected with TMA  
    :rtype: array_like (ndarray, 1)
    
    :Example:
        
        >>> from pycsamt.utils import zcalculator as Zcc
        >>> z2 = np.array([2.46073791 , 3.00162006 ,
        ...     9.74193019 ])# 1.82209497,
        ...     # 15.68879141 ,  12.91164264 ,
        ...     # 5.84384925 , 3.6899018, 
        ...     # 2.4430065  ,0.57175607])
        >>> tma= Zcc.compute_TMA(data_array=z2,
        ...                     number_of_TMApoints= 7)
        ... print(tma)

    """
    
    if data_array is None :
        raise CSex.pyCSAMTError_processing(
            'NoneType arguments can not be computed!.')
    
    if data_array.dtype not in ['float', 'int']: 
        try :data_array=np.array([float(dd) for dd in data_array])
        except : 
            raise CSex.pyCSAMTError_float(
                'Data must be on array_like float number.!')
    if type(number_of_TMApoints) is not int :
        try : 
            number_of_TMApoints=np.int(number_of_TMApoints)
        except : 
            raise CSex.pyCSAMTError_parameter_number(
                'TMA filter point must be integer.')
    
    roll_TMA =[]
    # let compute TMA for one points, 1 poinst will return the data array
    # Ascertqin number of TMA point and Data _array
    
    if number_of_TMApoints > len(data_array) : 
        warnings.warn("Number of TMA points = {}  is too large !".
                      format(number_of_TMApoints))
        number_of_TMApoints =5 # reset temporray to 5.
        
    if len(data_array) < 5: #find the TMA point accessible for computation
        # window width 
        number_of_TMApoints = len(data_array)
        
        if len(data_array) <= 2:
            
            warnings.warn(" Number of sites explored are too short !")
            print("--> Investigation sites too short !")
            
            print('** {0:<27} {1} {2} pts.'.format('Filter width  ',
                                '=', int(number_of_TMApoints))) 
            
            if len(data_array) ==1 :  return data_array 
            if len(data_array) ==2 : # take the mean of both value 
                return np.repeat(data_array.mean(), len(data_array))
        
        if len(data_array)==3: # loop the window and take mean value at each time
            print("--> ! Less number of TMA points  for higher frequencies mitigation.!")
                  
            for kk in range(len(data_array)): # add add overall mean to the last value 
                if kk ==len(data_array)-1:
                    tma = np.array([data_array.mean(),
                                    data_array[-1]]).mean()
                    roll_TMA.append(float(tma)) 
                else : 
                    roll_TMA.append(data_array[kk:kk+2].mean())
            
            print('** {0:<27} {1} {2} pts.'.format('Filter width  ',
                                '=', int(number_of_TMApoints)))     
            return np.array(roll_TMA)
        
     
    if len(data_array) >= 5 and number_of_TMApoints < 5: 

        if number_of_TMApoints <= 3 :
            warnings.warn("Number of points provided to compute TMA are not enough "
                          "for higher frequencies reduction.TMA points less than three (3)  "
                          " are not consistent to correct data."
                          )
            print("-->  Number of points are not enough for higher frequencies attenuation."
                  "For consistency, provide TMA points more than three.")
        
        if len(data_array)== 5 : mesf ="equal to"
        else : mesf = "more than"
        
        warnings.warn("Number of sites explored is estimated "
                      "to {}, then  number of TMA points is resetting "
                      "to 5 as default value.".format(len(data_array)))
        
        print("--> Number of explored sites is {0}  five(5).".format(mesf))
        print("--> TMA points was resseting  to 5.")
        
        number_of_TMApoints=5
        
    
    xk_cp  = np.int(np.trunc(number_of_TMApoints/2)) # centered point
    
    cut_off=-1          #starting point 
    for kk, value in enumerate(data_array): 
        if kk < xk_cp : 
            cut_off +=1 
            tma= data_array[kk-cut_off: kk+ xk_cp+1]
            
        elif xk_cp <= kk <= data_array.size - xk_cp:
            tma = data_array[kk- xk_cp:kk + xk_cp +1] # for Python index 
        
        elif kk > data_array.size - xk_cp:
            tma = data_array[kk-xk_cp:]

        get_TMA_value= np.delete(tma, 
                       np.array(tma.argmin(), tma.argmax())
                        )
        #tma = np.delete(tma, np.argmax(rho_tem), axis = 0)
        roll_TMA.append(get_TMA_value.mean())
        
    print('** {0:<27} {1} {2} pts.'.format('Filter width  ',
                                '=', int(number_of_TMApoints)))     
    return np.array(roll_TMA)   

@deprecated('Function too  expensive, redirected to `compute_TMA` more faster.')    
@deprecated_to(compute_TMA, "short and faster than the func below")    
def compute_trimming_moving_average (data_array=None, number_of_TMApoints=None ):
    """
    function to compute a trimmed-moving-average filter
    to estimate average apparent resistivities.

    :param data_array:  content of value to be trimmed 
    :type data_array: array_like(ndarray,1) 
    
    :param number_of_TMA points:  number of filter points .
    :type number_of_TMA points: int 
    
    :returns:  value corrected with TMA  
    :rtype: array_like (ndarray, 1)
          
    :Example:
        
        >>> from pycsamt.utils import zcalculator as Zcc
        >>> z2 = np.array([2.46073791 , 3.00162006 ,
        ...     9.74193019 ])# 1.82209497,
        ...     # 15.68879141 ,  12.91164264 ,
        ...     # 5.84384925 , 3.6899018, 
        ...     # 2.4430065  ,0.57175607])
        >>>  tma= Zcc.compute_trimming_moving_average(data_array=z2,
        ...                                          number_of_TMApoints= 4)
        >>> print(tma)
        >>> print(z2.shape)
        >>> print(tma.shape)
    """
    
    if data_array is None or number_of_TMApoints is None :
        raise CSex.pyCSAMTError_processing(
            'NoneType arguments can not be computed!.')
        
    if data_array.dtype not in ['float', 'int']: 
        try :data_array=np.array([float(dd) for dd in data_array])
        except :
            raise CSex.pyCSAMTError_AvgData(
                'Data must be on array_like float number.!')
    if type(number_of_TMApoints) is not int :
        try : number_of_TMApoints=np.int(number_of_TMApoints)
        except : 
            raise CSex.pyCSAMTError_parameter_number(
                'TMA filter point must be integer.')
    
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
           
def compute_AMA(reference_freq=None, z_array =None, 
                 number_of_skin_depth=1., dipole_length = 50., 
                  **kwargs  ):
    """
    Use an adaptive-moving-average filter to estimate average apparent
    resistivities at a single static-correction-reference frequency. Adaptive
    filtering is based on ideas presented inTorres-Verdin and Bostick, 1992,
    Principles of spatial surface electric field filtering in magnetotellurics
    -electromagnetic array profiling (EMAP), Geophysics, v57, p603-622.
    
    :param reference_frequency: frequency with clean data in Hertz 
    :type reference_frequency: float 
    
    :param z_array:  array of uncorrected impedance values.
    :type z_array: arry_like, complex
    
    :param dipole_length: length of dipole on survey
    :type dipole_length: float 
    
    :param number_of_skin_depth: arbitrary real constant, can be changed 
                to any value between 1 and 10 skin depths, by experiency, 
                value shoulb be 1<= c <= 4
    :type number_of_skin_depth: int , default is 1.
    
    :param weighted_window: bandwidth of hanning window filter lengths,
                           skin_depth arrays at each station.
    :type weighted_window: array_like 
 
    :param app_rho: array of apparent resistivities in ohm.m , if not provided 
                    will use impedance zrho to compute app_rho with reference_frequency 
    :type app_rho: array_like 
    
    :returns: windowed hanning , filtered window.
    :rtype: array_like 
    
    .. note:: If resistivities array is provided, provided also uncorrected phase
            array  to compute `impedance array`.
    
    :Example:
        
        >>> from pycsamt.utils import zcalculator as Zcc
        >>> z1= np.array([2.46073791 +3.00162006j ,
        ...     9.74193019 +1.82209497j,
        ...     15.68879141 + 12.91164264j ,
        ...     5.84384925 +3.6899018j, 
        ...     2.4430065  +0.57175607j])
        ...
        >>> reference_frequency = 8192.
        >>> znew= compute_AMA(reference_frequency=8192.,
        ...                   number_of_skin_depth= 1., dipole_length =50.,
        ...                   z_array=z1)
        ... print(znew)
    
    """
    
    if reference_freq is None : 
        raise CSex.pyCSAMTError_frequency(
            'Please add your reference frequency!')
        
    mu0 = 4* np.pi * 1e-7 
    omega_reffreq  = 2* np.pi *  reference_freq  
    
    resistivities = kwargs.pop('app_rho', None)
    phase =kwargs.pop('phase', None)
    adaptative_value =kwargs.pop('scale', 10.)

    if resistivities is None : 
        if  z_array is None :
            warnings.warn('`Impedance` is None. could not compute'
                          ' `apparent resistivities`')
            raise CSex.pyCSAMTError_processing('`Impedance`  is None.'
                'Could not compute `resistivities`.'
                ' Please provided at least an `impedance array`.')
    
        resistivities = np.abs(z_array) ** 2 / (mu0 * omega_reffreq )
        
    elif resistivities is not None and z_array is None :
        
        if phase is None  and reference_freq is None:
            warnings.warn('Could not compute `impedance` without phase values'
                          ' and reference frequency!')
            raise CSex.pyCSAMTError_processing(
                'None Type could not be computed. '
                'Please provided phase values |reference frequency')
            
        z_array = np.sqrt(resistivities * omega_reffreq * mu0 ) * (
            np.cos(phase)+  1j * np.sin( phase) ) 
    
    # compute skin depths base of uncorrected apparent resistivities
    skin_depths = 503. *  np.sqrt(resistivities/reference_freq)

    # window_width = number_of_skin_depth * skin_depth
    def adapted_skindepth (skin_window, scale_value=10.):
        """
        Adapted windo width and ans scale value to adapt value to be 
        close to dipole length 
        
        :param skin_window: array of skin depth 1D botstick penetration 
                or window width 
        :type skin_window: array_like 
        
        :return: window width scaled and adapted 
        :rtype: array_like 
        
        """
        # Adjusted each filter  iteratively at each station
        window_width = np.array([
            func.round_dipole_length(value, round_value=scale_value) # round dipole 
                                 for value in skin_window ])
        #adjusting filter width to with ten to ten range 
        for ii,  ww in enumerate(window_width): 
            if  ii ==0 or ii == len(window_width)-1:
                if ww < dipole_length:
                    window_width[ii]=func.round_dipole_length(
                        dipole_length, round_value=scale_value)
            elif ww < dipole_length:
                window_width[ii]=func.round_dipole_length(
                    window_width[ii-1:ii+2].mean(), round_value =scale_value)

        return  window_width
    
    def compute_filter_length(skindepth, ref_freq, 
                              dipole_length=dipole_length ):
        """
        :param skindepth: Bostick skin depth calculated with
                    apparent ressistivities at each station 
        :type skin_depth: array_like, 
        
        :param ref_freq: frequency at clean data at each station 
        :type ref_freq: float
        
        :param dipole_length: dipole length of the survey in meters 
        :type dipole_length: float, default is 50.
        
        :return: weight factor beta j at the station 
        :type: array_like 
        
        """
        npoints =np.around(skindepth/dipole_length) 
        # adjusted width  and recompute  and adapthed width 
        npoints_array = np.arange(0, dipole_length * npoints, dipole_length) 
                                    
        
        bj= np.array([compute_weight_factor_betaj(val, dipole_length=dipole_length, 
                                                  number_of_points=npoints)[0]
                      for val in  npoints_array])
  
        return bj, bj.sum()
    
    # collect list of each factor of all  station 
    skin_depths= adapted_skindepth(skin_window=skin_depths,
                                   scale_value=adaptative_value)

    weight_factors= [compute_filter_length(skindepth = skin, ref_freq=reference_freq, 
                              dipole_length=dipole_length )[0] for skin in skin_depths]
    
    # compute new impedance with base impedance and  list of weight factors 
    zxy_impedancef = compute_zxy_xk_omega(z_imp=z_array, coeffs_bj= weight_factors)
    
    # recompute new apparent resistivity 
    resistivities =np.abs( zxy_impedancef) ** 2 / (mu0 * omega_reffreq )
    
    # now recompute new window width that and multiply with arbiraty constant c
    z_bostick_depth = np.abs(zxy_impedancef) / (omega_reffreq * mu0)
    
    z_bostick_depth= adapted_skindepth(skin_window=z_bostick_depth,
                                    scale_value=adaptative_value)
    
    # adaptative window width with number of skin depth 
    w_window_width = number_of_skin_depth * z_bostick_depth 
    
    # recompute again impedance so to get adaptative impedance 
    n_weight_factors= [compute_filter_length(skindepth = skin, ref_freq=reference_freq, 
                              dipole_length=dipole_length )[0] for skin in w_window_width]
    
    #-----------------------text will FLMA computation --------
    # npoints= int(w_window_width.mean()/dipole_length)
    # cAMA= compute_FLMA(z_array=zxy_impedancef, dipole_length=dipole_length, 
    #                     number_of_points= npoints)
    #------------------------------------------------------------
    cAMA = compute_zxy_xk_omega(z_imp=z_array, coeffs_bj= n_weight_factors)
    
    print('** {0:<27} {1} {2} m.'.format("Filter's width",
                                '=', np.around(w_window_width.mean(),2)))    

    return cAMA  

def compute_zxy_xk_omega(z_imp, coeffs_bj):
    """
    Compute adaptative impedance with variable length Hanning window
    :param z_imp: impedance value to be filtered 
    :type z_imp: complex 
    
    :param coeff_bj: list coefficients of weighted window at each stations
    :type coeffs_bj : list 
    
    :return: array of  new impedance filtered 
    :rtype: complex 
    
    :Example: 
        
        >>> from pycsamt.utils import zcalculator as Zcc
        >>> z1= np.array([2.46073791 +3.00162006j ,
        ...     9.74193019 +1.82209497j,
        ...     15.68879141 + 12.91164264j ,
        ...     5.84384925 +3.6899018j, 
        ...     2.4430065  +0.57175607j])
        >>> weight_betaj = [array([0.18169011, 0.81830989]), 
        ...                array([0.02492092, 0.25, 0.47507908, 0.25]),
        ...                array([0.00224272, 0.02771308, 0.09220631, 0.16554531, 0.21341394,
        ...                0.21341394, 0.16554531, 0.09220631, 0.02771308]),
        ...                 array([0.05766889, 0.47116556, 0.47116556]),
        ...                 array([1.])]
        >>> Zcc.compute_zxy_xk_omega(coeffs_bj=weight_betaj , z_imp=z1)
        ... [2.01364616+2.45625537j 9.16556955+4.84395488j 6.83942027+4.21606008j
        ...     4.80923612+2.75254645j 2.4430065 +0.57175607j]   
                      
    """
    znew_f=[]
    for kk, zz in enumerate(z_imp):
        wlen=len(coeffs_bj[kk]) 
        xk_cp=int(wlen/2)
        try :
            if wlen%2 !=0:
                ss=z_imp[kk- xk_cp:kk+ xk_cp+1] * coeffs_bj[kk] # for Python index 
            elif wlen%2 ==0 :
                 ss=z_imp[kk- xk_cp:kk+ xk_cp] * coeffs_bj[kk]
        except:
            # check back lengh and forelength so to fill outside window with wk=0 
            b0 = len(coeffs_bj[kk][:xk_cp]) - len(z_imp[:kk]) # control the pad backward
            f0 = len(coeffs_bj[kk][xk_cp:]) - len(z_imp[kk:]) # control the pad forward 
            if b0 >0:               # if pad doesnt not meet 
                bnew= np.concatenate((np.repeat(0., b0),z_imp[:kk])) # add 0 outside 
            else :bnew = z_imp[kk-xk_cp:kk]
            if f0 >0:
                fnew= np.concatenate((z_imp[kk:], np.repeat(0., f0)))#add 0 outside the window
            else :
                if wlen %2 ==0 : fnew=z_imp[kk: kk + xk_cp] 
                else : fnew=z_imp[kk: kk + xk_cp +1] 
    
            ss = np.concatenate((bnew, fnew))* coeffs_bj[kk]
    
        znew_f.append(ss.sum())
                
    return np.array(znew_f)


def get_array_from_reffreq ( array_loc, freq_array,reffreq_value,
                            stnNames=None):
    """ 
    Get array value at special frequency
    Parameters
    ------------
        * array_loc : dict , 
            dictionnary of stations , array_value e.g: S00:(ndarray,1) 
            rho_values
            
        * freq_array : (ndarray,1) 
            frequency array for CSAMT survey 
            
        * reffreq_value : int or float 
            the value of frequency user want to get the value 
            
        * stnNames : list 
            list of stations names . 
    
    Returns 
    ---------
        array_like
            an array of all station with reffreq_value . 
            e.g reffreq_value =1024. it return all value of the array 
            at 1024Hz frequency . 
    
    """
    if stnNames is None : 
        raise CSex.pyCSAMTError_station('You may at least specify '
                        'the array or list of stations-Names or station id.')
    
    arrayObj,freqObj=array_loc ,freq_array
    rfObj,stnObj = reffreq_value,stnNames
    def get_reffreq_index(freq_array, reffreq_value): 
        """ 
        get the index of reference index. From this index,
        All array will filter data at this reffreq
        value . 
        
        """
    
        for ii, freq in enumerate(freq_array): 
            if freq == reffreq_value:
                index_rf = ii
                break
        return index_rf
    
    return np.array([values[get_reffreq_index(freq_array=freqObj, 
                                              reffreq_value=rfObj)] 
                for stn in stnObj 
                for keys, values in arrayObj.items() if stn==keys])
        
def relocate_on_dict_arrays(data_array, number_of_frequency,
                            station_names =None): 
    """ 
    Put data arrays on dictionnary where keys is each station and
    value the array of that station. if station_names is None ,
    program will create name of station. if station_names is given ,
    function will sorted stations names . please make sure to provide
    correctly station according  the disposal you want . 
   
    
    :param number_of_frequency: array of frequency during survey 
    :type number_of_frequency: array_like
    
    :param station_names:  list of station 
    :type station_names: list of array_like 
        
    :returns: infos at data stations 
    :rtype: dict 
    
    """
        
    if station_names is not None : 
        if type(station_names) is list :
            assert len(station_names) == (data_array.size / number_of_frequency),\
                CSex.pyCSAMTError_station(
                    'Number of Station provided must be <{0}> not'
                    ' <{1}>'.format(np.int(data_array.size / number_of_frequency), 
                                            len(station_names)))
        else :
            assert station_names.size == data_array.size / number_of_frequency, \
            CSex.pyCSAMTError_station(
                'Number of Station provided must be <{0}>'
                ' not <{1}>'.format(np.int(data_array.size / number_of_frequency), 
                                            station_names.size))
            station_names =station_names.tolist()
            
        NUM_STN=sorted(station_names)
    elif station_names is None  :
        station_length =np.int(data_array.size / number_of_frequency )
        
        NUM_STN=sorted((_numbering_station(number_of_station=station_length ,
                                    number_of_freq=number_of_frequency))[0])

    LIST_BROK = truncated_data(data=data_array,
                               number_of_reccurence=number_of_frequency)   
    
    return { key: value for key, value in zip (NUM_STN, LIST_BROK)}

def dipole_center_position (dipole_position=None):
    """
    Generaly positions are taken at each electrode of dipole to that 
    to easy correct data  for ploting and for noise correction ,
    we adjust coordinate by taking the center position that means ,
    the number of points will be substract to one.
    
    :param dipole_position: postion array at each electrodes. 
    :type dipole_position: array_like 
    
    :returns: centered position value  array 
    :rtype: array_like 
    
    """
    if dipole_position is None : 
        raise CSex.pyCSAMTError_inputarguments(
            'NoneType can not be compute. Please provide the right values.')
    try : dipole_position= np.array([float(pp) for pp in dipole_position])
    except : raise CSex.pyCSAMTError_float(
            'Cannot convert position values. Position must be a'
            ' number not str.Please check your data!.')
    # temp_pospp =[(dipole_position[pp]+dipole_position[pp-1])/2 for \
    #              pp, pos in enumerate( dipole_position) if (
    # pp>0 and pp <= dipole_position.size) ] 
    return np.array([(dipole_position[pp]+dipole_position[pp-1])/2 for 
                 pp, pos in enumerate( dipole_position) 
                 if ( pp>0 and pp <= dipole_position.size) ])
        

def get_matrix_from_dict(dict_array, flip_freq =False , 
                         freq_array= None , transpose =False ):
    """
    Function to collect dictionnary  station data with values as array 
    and build matrix along the line, The rowlines is assumed to be frequency 
    and colums as stations names. If `freq_array` is provided , `flip_freq` is
    not usefull.
    
    Parameters
    ----------
    dict_array : dict
        stations dictionnary array, keys are stations names and values are 
        stations values on array_like.
    transpose : bool, optional
        transpose  matrix , if true wwill transpose data and 
        stations should be read as rowlines and frequency as columns.
        The default is False.
    freq_array : array , optional 
        array of frequency, if provided will check if frequency are range from 
        highest to lowest .Defalut is True
    flip_freq :bool, optional 
        set to false if you want frequency to be lowest to highest otherwise
         frequency are sorted Highest to lowest. Default is False

    Returns
    -------
    ndarray(len(stations,), len(frequency))
        matrix of dictionnary values 
    bool, 
        flip_freq, ascertain if frequency array is flipped or not. 
    """
    
    if not isinstance(dict_array, dict): 
        raise CSex.pyCSAMTError_processing('Main argument `dict array`'
                                           ' should be dictionnary'
                                           ', not {0}'.type(dict_array))
        
    if freq_array is not None :
        if isinstance(freq_array, (list, tuple)):
            freq_array =np.array(freq_array)
            try :
                freq_array=freq_array.astype('float')
            except :
                raise CSex.pyCSAMTError_float(
                    'Could not convert array to value !')
        # now check the order 
        if freq_array[0] < freq_array[-1]: 
            freq_array = freq_array[::-1]
            flip_freq = True        # set flip to True to flip data values  
            
    tem= [values for stn , values in dict_array.items()]
    matrix = func.concat_array_from_list(tem, concat_axis=1)
    
    if flip_freq is True : matrix =matrix [::-1]    # reverse matrix 
    if transpose is True : matrix =matrix.T         # Transpose matrix 
    
    return matrix, flip_freq
        

def plot_reference_frequency (reference_array, frequency_array, 
                              data_array,
                              station_names =None ,  **kwargs):
    """
    Function to plot reference frequency 
    
    Parameters
    ----------
    reference_array : array_like, 
        array_of average_impedance Z_avg.
    frequency_array: array_like 
        array of frequency on sites 
    data_array : ndarray,
        nadarray of sounding curve,
        ndarray(len(frequency), len(number_ofstaions).

    Returns
    -------
    viewer
    """
    import matplotlib.pyplot as plt 
    
    ls =kwargs.pop('ls', '-')
    marker =kwargs.pop('marker', 'o')
    ms=kwargs.pop('ms', 0.2)
    lw =kwargs.pop('lw', 2)
    fs =kwargs.pop('fs', 0.8)
    if isinstance(data_array, dict): 
        data_array =sorted(data_array)  # sorted ditionnaries 
        station_names = list(data_array.keys())
        data_array = func.concat_array_from_list(list(data_array.values(), 
                                                      concat_axis=1))
    
    if station_names is None : 
        station_names =['S{0:02}'.format(ii) 
                        for ii in range(data_array.shape[1])]
    
    fig =plt.figure()
    ax =fig.add_subplot(111)
    
    ax.loglog(frequency_array,
                reference_array, 
                c='r', ls =ls, lw=lw, ms=ms *fs , 
                marker =marker)
    for ii in range(data_array.shape[1]):
        
        ax.loglog(frequency_array, data_array[:, ii])
        
    plt.show()
    
def _numbering_station(number_of_station, number_of_freq): 
    """
    small function to numbering stations .
    
    Prameters 
    ---------
        * number_of_station : int  
            number of station found on the site. 
        * number_of freq : int  
            number of frequency found for survey at each station.
        
    Returns
    --------
        numbsta : str 
              station numbered 
        poly_sta : list 
             multplie station for each frequency for each stations.
    """
   
    numbsta =['S{:02}'.format(ii) for ii in range (number_of_station)]
    poly_sta=numbsta*number_of_freq 
    poly_sta.sort()
    
    return numbsta , poly_sta
   
def truncated_data (data, number_of_reccurence, **kwargs):
    """
    Function to truncate all data according to number of frequency.
    
    Parameters
    ----------
        * data : list, or nd.array 
                data must be truncate.
        * number_of_freq : int
                number of frequency imaged.

    Returns
    -------
        list
            loc_list , data truncated on list.

    """
    if type(data) is list : 
        data =np.array(data) 
    value_lengh= data.shape[0]
    # index_loc = [np.int(ss) for ss in np.linspace(0,
    #value_lengh, number_of_station)]
    index_loc = [ss for ss in np.arange(0,value_lengh,
                                        number_of_reccurence)]

    loc_list=[]
    for ss, value  in enumerate(index_loc) : 
        if ss ==len(index_loc)-1 : 
            loc_list.append(data[value:])
            break
        loc_list.append(data[value:index_loc[ss+1]])
            
    return loc_list 
    
                          
if __name__=='__main__': 
    
    z1= np.array([2.46073791 +3.00162006j ,
         9.74193019 +1.82209497j,
         15.68879141 + 12.91164264j ,
         5.84384925 +3.6899018j, 
         2.4430065  +0.57175607j])
    
    reference_frequency = 8192.
    resistivities = np.abs(z1) ** 2 / (4* np.pi*1e-7 * 2* np.pi * reference_frequency )
    dipole_length =50
    znew= compute_AMA(reference_freq=reference_frequency,
                       number_of_skin_depth= 1., dipole_length =50.,
                       z_array=z1)
    print(resistivities)
    new_resistivities = np.abs(znew) ** 2 / (4* np.pi*1e-7 * 2* np.pi * reference_frequency )
    print(new_resistivities)
  
    