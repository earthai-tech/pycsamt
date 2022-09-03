# -*- coding: utf-8 -*-
#       Created on Sat Dec 12 13:55:47 2020
#       Author: Kouadio K.Laurent<etanoyau@gmail.com>
#       Licence: GPL

"""
Processing correction tools 
Created on Mon Aug 15 13:46:33 2022

@author: Daniel

.. _pyCSAMT: https://github.com/WEgeophysics/pycsamt
.. _MTpy: https://github.com/MTgeophysics/mtpy
.. |MT| replace:: Magnetetolluric 
.. |AMT| replace:: Audio-Magnetotellurics 
.. |CSAMT| replace:: Controlled Source |AMT| 
.. |NSAMT| replace:: Natural Source |AMT| 
.. |EM| replace:: electromagnetic
.. |EMAP| replace:: |EM| array profiling

"""
import os 
import re 
import warnings 
from math import factorial
    
import numpy as np
import pandas as pd 
from scipy.signal import fftconvolve
from scipy.integrate import quad 

from pycsamt.core.edi import Edi 
import pycsamt.core.z as MTz
from pycsamt.utils.decorator import donothing  
from pycsamt.utils.func_utils import (
    _assert_all_types, 
    _assert_edi_obj,
    subprocess_module_installation, 
    reshape, 
    scale_values, 
    ismissing, 
    spi,
    smart_format, 
    fillNaN, 
    get_ediObjs 
    )
from pycsamt.__init__ import (
    itqdm, 
    inumba, 
    imtpy 
    )

if itqdm: 
    import tqdm
if inumba: 
    import numba as nb

if imtpy: import mtpy
else:
    warnings.warn("Module 'MTpy' is not detected! Install it mannualy.")
    imtpy = False


try: 
    from sklearn.decomposition import PCA 
except : 
    is_success = subprocess_module_installation('sklearn')
    if not is_success : 
        raise ImportError( 'Could not import module `sklearn`. Please '
                          'install scikit-learn manually.')

mu0 = 4 * np.pi * 1e-7 

""" `d_hanning_window`

Discrete hanning function. For futher details, please refer to
 https://doi.org/10.1190/1.2400625

:param x: variable point along the window width
:param xk: Center of the window `W`. It presumes to host the most weigth.   
:param W: int, window-size; preferably set to odd number. It must be less than
     the dipole length. 
:return: Anonymous function (x,xk, W)
"""
d_hanning_window = lambda x, xk, W: 1/W * (1 + np.cos (2 * np.pi * (x-xk) /W)
                                            ) if np.abs(x-xk) <= W/2 else  0. 

def rhoa2z ( rhoa, phs, freq):
    r""" Convert apparent resistivity to impendance tensor z 
    
    :param rhoa: Apparent resistivity in :math:`\Omega.m` 
    :type rhoa: ndarray, shape (N, M) 
    
    :param phs: Phase in degrees 
    :type phs: ndarray, shape (N, M) 
    :param freq: Frequency in Hertz
    :type freq: array-like , shape (N, )
    :
    :return: Impendance tensor; Tensor is a complex number in :math:`\Omega`.  
    :rtype: ndarray, shape (N, M), dtype = 'complex' 
    
    :example: 
    >>> import numpy as np 
    >>> rhoa = np.array([1623.73691735])
    >>> phz = np.array([45.])
    >>> f = np.array ([1014])
    >>> rhoa2z(rhoa, phz, f)
    ... array([[2.54950976+2.54950976j]])
    
    """
    rhoa = np.array(rhoa); freq = np.array(freq) ; phs = np.array(phs) 
    
    if len(phs) != len(rhoa): 
        raise ValueError ("Phase and rhoa must have the same length."
                          f" {len(phs)} & {len(rhoa)} are given.")

    if len(freq) != len(rhoa): 
        raise ValueError("frequency and rhoa must have the same length."
                         "{len(freq} & {len(rhoa)} are given.")
        
    omega0 = 2 * np.pi * freq[:, None]
    z= np.sqrt(rhoa * omega0 * mu0 ) * (np.cos (
        np.deg2rad(phs)) + 1j * np.sin(np.deg2rad(phs)))
    
    return z 

def z2rhoa (z, freq):
    r""" Convert impendance tensor z  to apparent resistivity
    
    :param z: Impedance tensor  in :math:`\Omega` 
    :type z: ndarray, shape (N, M) 
 
    :param freq: Frequency in Hertz
    :type freq: array-like , shape (N, )
    :
    :return: Apparent resistivity in :math:`\Omega.m`  
    :rtype: ndarray, shape (N, M) 
    
    :example: 
    >>> import numpy as np 
    >>> z = np.array([2 + 1j *3 ])
    >>> f = np.array ([1014])
    >>> z2rhoa(z, f)
    ... array([[1623.73691735]])
        
    """

    z = np.array(z, dtype = 'complex' ) ;  freq = np.array(freq)

    if len(freq) != len(z): 
        raise ValueError("frequency and tensor z must have the same length."
                         "{len(freq} & {len(z)} are given.")
 
    return np.abs(z)**2 / (2 * np.pi * freq[:, None] * mu0 )

def get_reference_frequency (ediObjs= None, arr2d=None, freqs=None,
                             to_log10=False): 
    """ Get the reference frequency from collection Edis objects.
    
    The highest frequency with clean data should be selected as the reference 
    frequency
    
    Parameters 
    ---------- 
    ediObjs: list  of  pycsamt.core.edi.Edi or mtpy.core.edi.Edi objects 
        Collections of EDI-objects from `pyCSAMT`_ and `MTpy`_ packages 
        
    arr2d: np.ndarray (N, M)
        2D dimensional array of number of frequency x number of stations/sites, 
        No need if ediObjs is given.
        
    freqs: array-like, shape (N)
        Frequency array. It should be the complete frequency used during the 
        survey area. It can be get using the :func:`get_full_frequency` 
        No need if ediObjs is provided. 
        
    to_log10: bool, 
        outputs the reference frequency into base 10 logarithm in Hz.
    
    Returns 
    -------
    rf : float 
        the reference frequency at the clean data in Hz 
    Examples 
    ---------
    >>> from pycsamt.processimg import make2d, get_reference_frequency 
    >>> edipath ='data/3edis'
    >>> ediObjs = get_ediObjs (edipath)
    >>> cfs= get_full_frequency(ediObjs)
    >>> res2d = make2d (ediObjs, 'resyx') # make 2d resistivity blocks at yx  
    >>> ref = get_reference_frequency(arr2d= res2d, freqs= cfs, to_log10=True)
    >>> ref 
    ... 4.845098040014257 # in Hz 
    
    References 
    ----------
    http://www.zonge.com/legacy/PDF_DatPro/Astatic.pdf
    
    """

    if isinstance(ediObjs, str) : 
        # presume uses give a path 
        ediObjs = get_ediObjs(ediObjs)
    
    if ediObjs is None and arr2d is None: 
        raise ValueError('NoneType can not be read!')
    if ediObjs is not None: 
        freqs= get_full_frequency(ediObjs)
        # fit z and find all missing data from complete frequency f 
        # we take only the componet xy for fitting.

        zxy = [fit_tensor(freqs, ediObj.Z._freq, ediObj.Z.z[:, 0, 1].real)
              for ediObj in ediObjs
              ]
        # stacked the z values alomx axis=1. 
        arr2d = np.hstack ([ reshape (o, axis=0) for o in zxy])
    if freqs is None:
        raise ValueError("Array of all frequencies is missing.")
        
    freqs = np.array(freqs) # for consistency 
    if len(freqs) != len(arr2d): 
        # shrunk frequency to match the length of 2D array.
        freqs =freqs [: len(arr2d)]
        
    ix_nan = reshape (np.argwhere(np.isnan(arr2d).any(axis =1) ))
        # create bool array and mask the row of NaN 
    mask = np.full_like (freqs, fill_value = True , dtype=bool)
    mask[[*ix_nan] ] = False 
    # get the reference frequency and index 
    return  freqs [mask].max() if not to_log10 else np.log10(
        freqs [mask].max())
       
def interpolate2d (arr2d , method = 'slinear', **kws): 
    """ Interpolate the data in 2D dimensional array. 
    
    If the data contains some missing values. It should be replaced by the 
    interpolated values. 
    
    Parameters 
    -----------
    arr2d : np.ndarray, shape  (N, M)
        2D dimensional data 
        
    method: str, default ``linear``
        Interpolation technique to use. Can be ``nearest``or ``pad``. 
    
    kws: dict 
        Additional keywords. Refer to :func:`~.interpolate1d`. 
        
    Returns 
    -------
    arr2d:  np.ndarray, shape  (N, M)
        2D dimensional data interpolated 
    
    Examples 
    ---------
    >>> from pycsamt.processing import make2d , interpolate2d 
    >>> path2edi ='data/3edis'
    >>> freq2d = make2d (path2edi, 'freq') # make 2d matrix of frequency 
    >>> freq2d_i = interpolate2d(freq2d ) 
    >>> freq2d.shape 
    ...(55, 3)
    >>> freq2d 
    ... array([[7.00000e+04, 7.00000e+04, 7.00000e+04],
           [5.88000e+04, 5.88000e+04, 5.88000e+04],
           ...
            [6.87500e+00, 6.87500e+00, 6.87500e+00],
            [        nan,         nan, 5.62500e+00]])
    >>> freq2d_i
    ... array([[7.000000e+04, 7.000000e+04, 7.000000e+04],
           [5.880000e+04, 5.880000e+04, 5.880000e+04],
           ...
           [6.875000e+00, 6.875000e+00, 6.875000e+00],
           [5.625000e+00, 5.625000e+00, 5.625000e+00]])
    
    References 
    ----------
    
    https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.interpolate.html
    https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.interpolate.interp2d.html        
        
    """ 
    arr2d = np.array(arr2d)
    if len(arr2d.shape) ==1: 
        arr2d = arr2d[:, None] # put on 
    if arr2d.shape[0] ==1: 
        arr2d = reshape (arr2d, axis=0)

    arr2d  = np.hstack ([ reshape (interpolate1d(arr2d[:, ii], kind=method, 
                                        method ='pd', **kws), axis=0)
                         for ii in  range (arr2d.shape[1])])
    return arr2d 


def get_full_frequency (ediObjs, to_log10 =False ): 
    """ Get the frequency with clean data. 
    
    The full or plain frequency is array frequency with no missing  data during 
    the data collection. Note that when using |NSAMT|, some data are missing 
    due to the weak of missing frequency at certain band especially in the 
    attenuation band. 
    
    :param ediObjs: Collections of EDI-objects from `pyCSAMT`_ and `MTpy`_ 
        packages or full path to EDI files. 
    :type ediObjs: str or :mod:`~.core.edi.Edi_collection` object 
    
    :param to_log10: Export frequency to base 10 logarithm 
    :type to_log10: bool, 
    
    :returns: frequency with clean data. Out of `attenuation band` if survey 
        is completed with  |NSAMT|. 
    :rtype: array_like, shape(N, )

    :example: 
        >>> from pycsamt.processing.ctools import get_full_frequency
        >>> edipath = 'data/3edis' 
        >>> ref = get_full_frequency(edipath)  
        >>> ref
        ... array([7.00000e+04, 5.88000e+04, 4.95000e+04, 4.16000e+04, 3.50000e+04,
               2.94000e+04, 2.47000e+04, 2.08000e+04, 1.75000e+04, 1.47000e+04,
               ...
               1.12500e+01, 9.37500e+00, 8.12500e+00, 6.87500e+00, 5.62500e+00])
        >>> len(ref)
        ... 55 
        
    """
    if isinstance(ediObjs, str) : 
        # presume uses give a path 
        ediObjs = get_ediObjs(ediObjs)
    ediObjs = np.array( list(map( lambda o: _assert_edi_obj(o),  ediObjs )),
                       dtype =object) 
    lenfs = np.array([len(ediObj.Z._freq) for ediObj in ediObjs ] ) 
    ix_fm = np.argmax (lenfs) ; f= ediObjs [ix_fm].Z._freq 
    
    return np.log10(f) if to_log10 else f 


def make2d (ediObjs, out = 'resxy', *, kind = 'complex' , **kws ): 
    """ Out 2D resistivity, phase (error) and tensor matrix from EDI-collection
    objects. Matrix for number of frequency x number of sites. 
    
    The function asserts whether all data from all frequencies are available. 
    The missing values should be filled by NaN. 
    
    Parameters 
    ----------- 
    ediObjs: str or :mod:`~.core.edi.Edi_collection` object 
        Collections of EDI-objects from `pyCSAMT`_ and `MTpy`_ packages 
        or full path to EDI files.
    out: str 
        kind of data to output. Be sure to provide the component to retrieve 
        the attribute from the collection object. Except the `error` and 
        frequency attribute, the missing component to the attribute will 
        raise an error. for instance ``resxy`` for xy component. Default is 
        ``zxy``. 
    kind : bool or str 
        focus on the tensor output. Note that the tensor is a complex number 
        of ndarray (nfreq, 2,2 ). If set to``modulus`, the modulus of the complex 
        tensor should be outputted. If ``real`` or``imag``, it returns only
        the specific one. Default is ``complex``.
        
    kws: dict 
        Additional keywords arguments from :func:`~.get_full_frequency`. 
    
    Returns 
    -------- 
    mat2d : np.ndarray(nfreq, nstations) 
        the matrix of number of frequency and number of Edi-collectes which 
        correspond to the number of the stations/sites. 
    
    Examples 
    ---------
    >>> from pycsamt.processing.ctools import make2d, get_ediObjs
    >>> edipath ='data/3edis'
    >>> ediObjs = get_ediObjs (edipath) # get the collection objects 
    >>> # try to export phase at component yx 
    >>> phyx = make2d (ediObjs, 'phaseyx')
    >>> phyx 
    ... array([[ 26.42546593,  32.71066454,  30.9222746 ],
           [ 44.25990541,  40.77911136,  41.0339148 ],
           ...
           [ 37.66594686,  33.03375863,  35.75420802],
           [         nan,          nan,  44.04498791]])
    >>> phyx.shape 
    ... (55, 3)
    >>> # get the real number of the yy componet of tensor z 
    >>> zyy_r = make2d (ediObjs, 'zyx', kind ='real')
    ... array([[ 4165.6   ,  8665.64  ,  5285.47  ],
           [ 7072.81  , 11663.1   ,  6900.33  ],
           ...
           [   90.7099,   119.505 ,   122.343 ],
           [       nan,        nan,    88.0624]])
    >>> # get the resistivity error of component 'xy'
    >>> resxy_err = make2d (ediObjs, 'resxy_err')
    >>> resxy_err 
    ... array([[0.01329037, 0.02942557, 0.0176034 ],
           [0.0335909 , 0.05238863, 0.03111475],
           ...
           [3.33359942, 4.14684926, 4.38562271],
           [       nan,        nan, 4.35605603]])
    >>> phyx.shape ,zyy_r.shape, resxy_err.shape  
    ... ((55, 3), (55, 3), (55, 3))
    
    """
    def fit2dall(objs, attr, comp): 
        """ Read all ediObjs and replace all missing data by NaN value. 
        
        This is useful to let the arrays at each station to  match the length 
        of the complete frequency rather than shrunking  up some data. The 
        missing data should be filled by NaN values. 
        
        """
        zl = [getattr( ediObj.Z, f"{attr}")[tuple (_c.get(comp))]
              for ediObj in objs ]
        
        if name =='z': 
            if kind =='modulus': 
                zl = [ np.abs (v) for v in zl]
                zl = [fit_tensor(ref, ediObj.Z._freq, v)
                      for ediObj ,  v  in zip(objs, zl)]
            if kind in ('real' , 'complex') : 
                zr = [fit_tensor(ref, ediObj.Z._freq, v.real)
                      for ediObj ,  v  in zip(objs, zl)]
                
            if kind in ('imag', 'complex'): 
                zi= [fit_tensor(ref, ediObj.Z._freq, v.imag)
                      for ediObj ,  v  in zip(objs, zl)]
                
            if kind =='complex': 
                zl = [ r + 1j * im for r, im in zip (zr, zi)]
                
                
            zl = zl if kind in ('modulus', 'complex') else (
                zr if kind =='real' else zi )    
        else : 
            zl = [fit_tensor(ref, ediObj.Z._freq, v)
                  for ediObj ,  v  in zip(objs, zl)]
            
        # stacked the z values alomx axis=1. 
        return np.hstack ([ reshape (o, axis=0) for o in zl])
        
    out = str(out).lower().strip () 
    kind = str(kind).lower().strip() 
    if kind.find('imag')>=0 :
        kind ='imag'
    if kind not in ('modulus', 'imag', 'real', 'complex'): 
        raise ValueError(f"Unacceptable argument {kind!r}. Expect 'modulus',"
                         "'imag', 'real', or 'complex'.")
    # get the name for extraction using regex 
    regex1= re.compile(r'res|rho|phase|phs|z|tensor|freq')
    regex2 = re.compile (r'xx|xy|yx|yy')
    regex3 = re.compile (r'err')
    
    m1 = regex1.search(out) 
    m2= regex2.search (out)
    m3 = regex3.search(out)
    
    if m1 is None: 
        raise ValueError (f" {out!r} does not match  any 'resistivity', 'phase'"
                          " 'tensor' nor 'frequency'.")
    m1 = m1.group() 
    
    if m1 in ('res', 'rho'):
        m1 = 'resistivity'
    if m1 in ('phase', 'phs'):
        m1 = 'phase' 
    if m1 in ('z', 'tensor'):
        m1 ='z' 
    if m1  =='freq':
        m1 ='_freq'
        
    if m2 is None or m2 =='': 
        if m1 in ('z', 'resistivity', 'phase'): 
            raise ValueError (f"{'Tensor' if m1=='z' else m1.title()!r} component "
                              f"is missing. Use e.g. '{m1}_xy' for 'xy' component")
    m2 = m2.group() if m2 is not None else m2 
    m3 = m3.group () if m3 is not None else '' 
    
    if m3 =='err':
        m3 ='_err'
    # read/assert edis and get the complete frequency 
    if isinstance(ediObjs, str) : 
       ediObjs = get_ediObjs(ediObjs)

    ediObjs = np.array( list(map( lambda o: _assert_edi_obj(o),  ediObjs )),
                       dtype =object)
    ref = get_full_frequency(ediObjs, **kws)
    
    #=> slice index for component retreiving purpose 
    _c= {
          'xx': [slice (None, len(ref)), 0 , 0] , 
          'xy': [slice (None, len(ref)), 0 , 1], 
          'yx': [slice (None, len(ref)), 1 , 0], 
          'yy': [slice (None, len(ref)), 1,  1] 
    }
    #==> returns mat2d freq 
    if m1 =='_freq': 
        f2d  = [fit_tensor(ref, ediObj.Z._freq, ediObj.Z._freq)
              for ediObj in ediObjs
              ]
        return  np.hstack ([ reshape (o, axis=0) for o in f2d])
    
    # get the value for exportation (attribute name and components)
    name = m1 + m3 if (m3 =='_err' and m1 != ('_freq' or 'z')) else m1 
    #print(name, m1 , m2)
    mat2d  = fit2dall(objs= ediObjs, attr= name, comp= m2)
    
    return mat2d 


def tma (ediObjs=None, res2d=None, phs2d=None, freqs= None,  window_size =5,
         component = 'xy', mode='same', method = 'slinear',  out = 'srho', 
         c= None, **kws): 
    """ A trimmed-moving-average filter to estimate average apparent
    resistivities at a single static-correction-reference frequency. 
    
    The TMA filter option estimates near-surface resistivities by averaging
    apparent resistivities along line at the selected static-correction 
    reference frequency. The highest frequency with clean data should be 
    selected as the reference frequency.
    
    Parameters 
    ----------
    ediObjs: list  of  pycsamt.core.edi.Edi or mtpy.core.edi.Edi objects 
        Collections of EDI-objects from `pyCSAMT`_ and `MTpy`_ packages
        
    res2d : np.ndarray, shape  (N, M)
        2D dimensional resistivity data. N for number of frequency and 
        M for number of stations/sites. No need to provide if ediObjs is given.
        
    phs2d: np.ndarray (N, M)
        2D dimensional phase array of number of frequency x number of stations/sites, 
        No need to provide if ediObjs is given.
        
    freqs: array-like, shape (N)
        Frequency array. It should be the complete frequency used during the 
        survey area. It can be get using the :func:`get_full_frequency` 
        No need if ediObjs is provided. 
        
    window_size : int
        the length of the window. Must be greater than 1 and preferably
        an odd integer number. Default is ``5``
        
    component: str 
       field tensors direction. It can be ``xx``, ``xy``,``yx``, ``yy``. If 
       `arr2d`` is provided, no need to give an argument. It become useful 
       when a collection of EDI-objects is provided. If don't specify, the 
       resistivity and phase value at component `xy` should be fetched for 
       correction by default. Change the component value to get the appropriate 
       data for correction. Default is ``xy``.
       
    mode: str 
        mode of the border trimming. Should be 'valid' or 'same'.'valid' is used 
        for regular trimimg whereas the 'same' is used for appending the first
        and last value of resistivity. Any other argument except 'valid' should 
        be considered as 'same' argument. Default is ``same``.     
       
    method: str, default ``slinear``
        Interpolation technique to use. Can be ``nearest``or ``pad``. Refer to 
        the documentation of :doc:`~.interpolate2d`. 
        
    out : str 
        Value to export. Can be ``sfactor``, ``tensor`` for corrections factor 
        and impedance tensor. Any other values will export the static corrected  
        resistivity. 
        
    c: NoneType, 
        Here, `c` does nothing. It is used for API purpose. 

    Returns 
    -------
    rc or cf: np.ndarray, shape  (N, M)
        EMAP apparent  resistivity static shift corrected or static 
        correction factor or impedance tensor. 
        
    Examples 
    --------
    >>> import matplotlib.pyplot as plt 
    >>> from pycsamt.processing import get_ediObjs, make2d, get_full_frequency 
    >>> from pycsamt.processing import tma 
    >>> edipath = 'data/3edis'
    >>> ediObjs = get_ediObjs(edipath) 
    >>> res2d = make2d (ediObjs, 'resyx')
    >>> phs2d = make2d (ediObjs, 'phaseyx')
    >>> freqs = get_full_frequency(ediObjs)
    >>> rc =tma(res2d =res2d, phs2d =phs2d, freqs= freqs , window_size=5)
    >>> res2d [3, :]  # get the resistivy value of the third frequency  at all stations 
    ... array([ 447.05423001, 1016.54352954, 1415.90992189,  536.54293994,
           1307.84456036,   65.44806698,   86.66817791,  241.76592273,
           ...
            248.29077039,  247.71452712,   17.03888414])
    >>> rc [3, :] # get the resistivity value corrected at the third frequency 
    ... array([ 447.05423001,  763.92416768,  929.33837349,  881.49992091,
            404.93382163,  190.58264151,  160.71917654,  163.30034875,
            394.2727092 ,  679.71542811,  953.2796567 , 1212.42883944,
            ...
            164.58282866,   96.60082159,   17.03888414])
    >>> plt.semilogy (np.arange (res2d.shape[1] ), res2d[3, :], '--',
                      np.arange (res2d.shape[1] ), rc[3, :], 'ok--')
 
    References 
    -----------
    .. [1] http://www.zonge.com/legacy/PDF_DatPro/Astatic.pdf
    
    
    """
    # assert filter arguments 
    res2d , phs2d , freqs, c, window_size, component, out = \
        _assert_emap_filter_args (ediObjs, res2d , phs2d , freqs, c,
                                  window_size,component, out)
    # get the reference frequency 
    rf = get_reference_frequency(arr2d=res2d, freqs= freqs) 
    #  interpolate resistivity and phases 
    phs2d= interpolate2d(phs2d, method =method, **kws)
    res2d= interpolate2d(res2d, method =method, **kws)
    # get the index of the reference frequency  and collect 
    # the resistivity and phase at that frequency 
    ix_rf = np.int(reshape (np.argwhere (freqs==rf)))  
    # normalize log frequency and take the normalize ref freq 
    norm_f = (np.log10(freqs) / np.linalg.norm(np.log10(freqs)))
    # compute the slope at each normalize frequency 
    slope2d = np.arctan( (np.deg2rad(phs2d) / (np.pi /4 )) -1 ) / (np.pi /2 )
    log_rho2d = np.log10 (res2d) + norm_f[:, None] * slope2d 
    # extrapolate up 
    # replace the up frequency thin the index of rf by interpolating up 
    log_rho2d [:ix_rf, :] = np.log10 (
        res2d[:ix_rf, : ]) + np.log10(np.sqrt(2)) * slope2d[:ix_rf, :]
    
    # make an array of weight factor wf 
    wf = np.zeros_like(log_rho2d) # average logj 
    # For each station collect a group of window-size log(rj ), 
    # #i.e. for window size =5 station index j, i = j-2 to j+2. 
    half_window = window_size //2 
    for ii in range(log_rho2d.shape[1]):
        
        if ii ==0 or ii ==log_rho2d.shape[1] -1: 
            w = (log_rho2d[ :, :ii + half_window +1] 
                 if  ii - half_window < 0 else 
                 log_rho2d[:, ii-half_window:] 
                 ) if mode =='valid' else log_rho2d[:, ii][:, None]
  
        elif ii - half_window < 0: 
            w= log_rho2d[ :, :ii + half_window +1]
            
        elif ii + half_window +1 > log_rho2d.shape[1]: 
            w= log_rho2d[:, ii-half_window:]

        else : 
            # Discard the lowest and highest valued log(rj ) 
            # from the group of five and average the remaining
            # three => avg_logj.
            w= log_rho2d[:, ii-half_window : half_window + ii + 1 ]
            try : 
                ls = [ np.delete (w[ii, :] , np.where (
                    np.logical_or(w[ii, :] ==w[ii, :].max(),
                                  w[ii, :] ==w[ii, :].min())
                    )) for ii in range (len(w))]
                
                w = np.vstack (ls)
            except : 
                # in the case the ls has some array with different length 
                # do the average manually and set as an array of axis 1.  
                ls = [np.sum(w[ii, :])/ len(w[ii, :]) for ii in range(len(w))]
                w = np.array(ls)[:, None] # make axis = 1
            
        wf[:, ii] = np.average (w, axis = 1)
        
    # compute the correction factor cf
    cf = np.power(10, wf, dtype =float)/ np. power(10, log_rho2d) 
    
    rc = res2d * cf 
    
    if out =='z': 
        rc = rhoa2z(rc, phs2d, freqs)
        # omega0 = 2 * np.pi * freqs
        # rc = np.sqrt(rc * omega0[:, None] * mu0 ) * (np.cos (
        #     np.deg2rad(phs2d)) + 1j * np.sin(np.deg2rad(phs2d)))
    
    return   cf if out =='sf' else rc   
            
def ama (ediObjs=None, res2d=None, phs2d=None, freqs= None, c=2, window_size =5,
          component = 'xy', method = 'slinear', out='rho', mode ='same', **kws ): 
    """ 
    Use an adaptive-moving-average filter to estimate average apparent 
    resistivities at a single static-correction-reference frequency.. 
    
    The AMA filter estimates static-corrected apparent resistivities at a 
    single reference frequency by calculating a profile of average impedances 
    along the length of the line. Sounding curves are then shifted so that they
    intersect the averaged profile. 
    
    Parameters 
    ----------
    ediObjs: list  of  pycsamt.core.edi.Edi or mtpy.core.edi.Edi objects 
        Collections of EDI-objects from `pyCSAMT`_ and `MTpy`_ packages
        
    res2d : np.ndarray, shape  (N, M)
        2D dimensional resistivity data. N for number of frequency and 
        M for number of stations/sites. No need to provide if ediObjs is given.
        
    phs2d: np.ndarray (N, M)
        2D dimensional phase array of number of frequency x number of stations/sites, 
        No need to provide if ediObjs is given.
        
    freqs: array-like, shape (N)
        Frequency array. It should be the complete frequency used during the 
        survey area. It can be get using the :func:`get_full_frequency` 
        No need if ediObjs is provided. 
        
    c : int, 
        A window-width expansion factor that must be input to the filter 
        adaptation process to control the roll-off characteristics
        of the applied Hanning window. It is recommended to select `c` between 
        ``1``  and ``4``.  Default is ``2``. 
        
    window_size : int
        the length of the window. It is known as the filter width. Must be 
        greater than 1 and preferably be an odd integer number. Default is set 
        to ``5`` dipole lengths. 
        
    component: str 
       field tensors direction. It can be ``xx``, ``xy``,``yx``, ``yy``. If 
       `arr2d`` is provided, no need to give an argument. It become useful 
       when a collection of EDI-objects is provided. If don't specify, the 
       resistivity and phase value at component `xy` should be fetched for 
       correction by default. Change the component value to get the appropriate 
       data for correction mode. Default is ``xy`` i.e TM mode.
       
    method: str, default ``slinear``
        Interpolation technique to use. Can be ``nearest``or ``pad``. Refer to 
        the documentation of :doc:`~.interpolate2d`. 
        
    out : str 
        Value to export. Can be ``tensor`` for corrected impedance tensor. 
        Any other values will export the static corrected resistivity. 
         

    Returns 
    -------
    rc or z: np.ndarray, shape  (N, M)
        EMAP apparent  resistivity static shift corrected  or static 
        correction tensor 
        
    Examples 
    --------
    >>> import matplotlib.pyplot as plt 
    >>> from pycsamt.processing import get_ediObjs, make2d, get_full_frequency 
    >>> from pycsamt.processing import ama 
    >>> edipath = 'data/3edis'
    >>> ediObjs = get_ediObjs(edipath) 
    >>> res2d = make2d (ediObjs, 'resyx')
    >>> phs2d = make2d (ediObjs, 'phaseyx')
    >>> freqs = get_full_frequency(ediObjs)
    >>> rca =ama(res2d =res2d, phs2d =phs2d, freqs= freqs , window_size=5)
    >>> res2d [3, :]  # get the resistivy value of the third frequency  at all stations 
    ... array([ 447.05423001, 1016.54352954, 1415.90992189,  536.54293994,
           1307.84456036,   65.44806698,   86.66817791,  241.76592273,
           ...
            248.29077039,  247.71452712,   17.03888414])
    >>> rca [3, :] # get the resistivity value corrected at the third frequency 
    ... array([ 447.05423001,  707.71978296,  799.30037014,  691.23464805,
            453.48970601,  274.3041412 ,  245.62906983,  422.99829174,
            827.21924566, 1192.63713823, 1150.52690907,  952.57203567,
            ...
            167.21190293,   88.38401361,   26.88590246])
    >>> plt.semilogy (np.arange (res2d.shape[1] ), res2d[3, :], '--',
                      np.arange (res2d.shape[1] ), rca[3, :], 'ok--')
 
    References 
    -----------
    .. [1] http://www.zonge.com/legacy/PDF_DatPro/Astatic.pdf
    .. [2] Torres-Verdin and Bostick, 1992,  Principles of spatial surface 
        electric field filtering in magnetotellurics: electromagnetic array profiling
        (EMAP), Geophysics, v57, p603-622.https://doi.org/10.1190/1.2400625
        
    """

    # assert filter arguments 
    res2d , phs2d , freqs, c, window_size, component, out = \
        _assert_emap_filter_args (ediObjs, res2d , phs2d , freqs, c,
                                  window_size,component, out)
    #  interpolate resistivity and phases 
    phs2d= interpolate2d(phs2d, method =method, **kws)
    res2d= interpolate2d(res2d, method =method, **kws)
    
    # convert app. resistivity and impedance phase  to 
    # impedance values, Zj, for each station
    omega0 = 2 * np.pi * freqs
    zj = np.sqrt(res2d * omega0[:, None] * mu0 ) * (np.cos (
        np.deg2rad(phs2d)) + 1j * np.sin(np.deg2rad(phs2d)))
    
    # compute the weight factor for convoluting 
    # L = dipole length = L : 1 is fixed dipole -length 
    w = np.array([betaj (xj = ii, L= 1 , W= window_size) 
                  for ii in range(window_size)])
    #print(w)
    zjr = np.zeros_like(res2d) 
    zji = zjr.copy() 
    
    for ii in range (len(zj)): 
        w_exp = [ k * window_size for k in range(1, c +1 )]
        zcr=list(); zci = list()
        # compute Zk(xk, w) iteratively
        # with adpatavive W expanded to 1 to c 
        for wind_k  in w_exp : 
            w= np.array([betaj (xj = jj, L= 1, W= wind_k) for jj in range(wind_k)
                         ])
            zcr.append(np.convolve(zj[ii, :].real, w[::-1], 'same'))
            zci.append(np.convolve(zj[ii, :].imag, w[::-1], 'same'))
        # and take the average 
        zjr [ii, :] = np.average (np.vstack (zcr), axis = 0)
        zji[ii, :] = np.average (np.vstack (zci), axis = 0)
           
        # omegaii = 2 * np.pi * freqs[ii]
        # W = c * 1./ (omegaii * mu0) * np.abs (zcr + 1j * zci )
        # W = np.average(W).astype (np.uint32)#//zj.shape[1]  
        #  W.astype (np.uint32)
        # if W < window_size :
        #     W = window_size 
        # if W > zj.shape[1]:
        #     W= zj.shape [1] 
        # wn= np.array([betaj (xj = jj, L= 1, W= W) for jj in range(W)
        #              ])
        # zjr[ii, :] = np.convolve(zj[ii, :].real, wn[::-1], 'same')
        # zji[ii, :] = np.convolve(zj[ii, :].imag, wn[::-1], 'same')
 
    zjc = zjr + 1j * zji 
    rc = z2rhoa(zjc, freqs)  #np.abs(zjc)**2 / (omega0[:, None] * mu0 )
    if mode =='same': 
        rc[:, 0] = res2d[:, 0]
        zjc[:, 0] = zj [:, 0]
    
    return zjc if out =='z' else rc 


def betaj (xj, L, W, **kws): 
    """ Weight factor function for convoluting at station/site j.
    
    The function deals with the discrete hanning window based on ideas presented 
    in Torres-Verdin and Bostick, 1992, https://doi.org/10.1190/1.2400625.
    
    :param xj: int, position of the point to compute its weight
    :param W: int, window size, presumes to be the number of dipole. 
    :param L: int : length of dipole  
    :param kws: dict , additional :func:`scipy.intergate.quad` functions.
    
    :return: Weight value at th eposition `xj`, prefix-`x`is used to specify  
        the direction. Commonly the survey direction is considered as `x`.
        
    :example: 
        >>> from pycsamt.processing import betaj 
        >>> # compute the weight point for window-size = 5 at position j =2
        >>> L= 1 ; W=5 
        >>> betaj (xj = 2 , L=L, W=W )
        ... 0.35136534572813144
    """
    if W < L : 
        raise ValueError("Window-size must be greater than the dipole length.")
        
    xk = W/2 
    # vec_betaj = np.vectorize( betaj ) ; vec_betaj(0, 1, 5)
    return  quad (d_hanning_window, xj - L/2 , xj +L/2, args=( xk, W), 
                  **kws)[0]
    
def flma (ediObjs=None, res2d=None, phs2d=None, freqs= None, c=None, window_size =5,  
          component = 'xy', method = 'slinear', out='rho', mode='same', **kws ): 
    """ Use a fixed-length-moving-average filter to estimate average apparent
    resistivities at a single static-correction-reference frequency. 
    
    The FLMA filter estimates static-corrected apparent resistivities at a 
    single reference frequency by calculating a profile of average impedances 
    along the length of the line. Sounding curves are then shifted so that they
    intersect the averaged profile. 
    
    Parameters 
    ----------
    ediObjs: list  of  pycsamt.core.edi.Edi or mtpy.core.edi.Edi objects 
        Collections of EDI-objects from `pyCSAMT`_ and `MTpy`_ packages
        
    res2d : np.ndarray, shape  (N, M)
        2D dimensional resistivity data. N for number of frequency and 
        M for number of stations/sites. No need to provide if ediObjs is given.
        
    phs2d: np.ndarray (N, M)
        2D dimensional phase array of number of frequency x number of stations/sites, 
        No need to provide if ediObjs is given.
        
    freqs: array-like, shape (N)
        Frequency array. It should be the complete frequency used during the 
        survey area. It can be get using the :func:`get_full_frequency` 
        No need if ediObjs is provided. 
        
    window_size : int
        the length of the window. It is known as the filter width. Must be 
        greater than 1 and preferably be an odd integer number. Default is set 
        to ``5`` dipole lengths. 
        
    component: str 
       field tensors direction. It can be ``xx``, ``xy``,``yx``, ``yy``. If 
       `arr2d`` is provided, no need to give an argument. It become useful 
       when a collection of EDI-objects is provided. If don't specify, the 
       resistivity and phase value at component `xy` should be fetched for 
       correction by default. Change the component value to get the appropriate 
       data for correction. Default is ``xy``.
       
    method: str, default ``slinear``
        Interpolation technique to use. Can be ``nearest``or ``pad``. Refer to 
        the documentation of :doc:`~.interpolate2d`. 
     
    mode: str 
         mode of the border prepending. Should be ``valid`` or ``same``. 
         ``same`` is used for prepending or appending the first t value of
         resistivity for smoothing. Any other argument should keep intact 
         the interpolation. Default is ``same``.     
         
    out : str 
        Value to export. Can be ``tensor`` for corrected impedance tensor. 
        Any other values will export the static corrected resistivity. 
         
    c: NoneType, 
        `c` in this filter does nothing. It is used for API purpose. 
        
    Returns 
    -------
    rc or z : np.ndarray, shape  (N, M)
        EMAP apparent  resistivity static shift corrected  or static 
        correction impedance tensor. 
        
    Examples 
    --------
    >>> import matplotlib.pyplot as plt 
    >>> from pycsamt.processing import get_ediObjs, make2d, get_full_frequency 
    >>> from pycsamt.processing import tma 
    >>> edipath = 'data/3edis'
    >>> ediObjs = get_ediObjs(edipath) 
    >>> res2d = make2d (ediObjs, 'resyx')
    >>> phs2d = make2d (ediObjs, 'phaseyx')
    >>> freqs = get_full_frequency(ediObjs)
    >>> rcf =flma(res2d =res2d, phs2d =phs2d, freqs= freqs , window_size=5)
    >>> res2d [3, :]  # get the resistivy value of the third frequency  at all stations 
    ... array([ 447.05423001, 1016.54352954, 1415.90992189,  536.54293994,
           1307.84456036,   65.44806698,   86.66817791,  241.76592273,
           ...
            248.29077039,  247.71452712,   17.03888414])
    >>> rcf [3, :] # get the resistivity value corrected at the third frequency 
    ... array([ 148.93643044,  581.68186236,  955.93375804,  962.1253775 ,
            738.58658567,  394.19629051,  175.09339708,  146.47552544,
            387.13690304, 1030.53706101, 1562.12079329, 1291.07477836,
            ...
            165.71236274,  165.19885314,   84.46198723])
    >>> plt.semilogy (np.arange (res2d.shape[1] ), res2d[3, :], '--',
                      np.arange (res2d.shape[1] ), rcf[3, :], 'ok--')
 
    References 
    -----------
    .. [1] http://www.zonge.com/legacy/PDF_DatPro/Astatic.pdf
    
    """

    # assert filter arguments 
    res2d , phs2d , freqs, c, window_size, component, out = \
        _assert_emap_filter_args (ediObjs, res2d , phs2d , freqs, c,
                                  window_size,component, out)
    #  interpolate resistivity and phases 
    phs2d= interpolate2d(phs2d, method =method, **kws)
    res2d= interpolate2d(res2d, method =method, **kws)
    
    # convert app. resistivity and impedance phase  to 
    #impedance values, Zj, for each station
    omega0 = 2 * np.pi * freqs
    zj = np.sqrt(res2d * omega0[:, None] * mu0 ) * (np.cos (
        np.deg2rad(phs2d)) + 1j * np.sin(np.deg2rad(phs2d)))
    
    # compute the weight factor for convoluting 
    # L = dipole length = L
    w = np.array([betaj (xj = ii, L= 1 , W= window_size) 
                  for ii in range(window_size)])
    
    zjr = np.zeros_like(res2d) 
    zji = zjr.copy() 
    for ii in range(len(zjr)) :
        zjr[ii, :] = np.convolve(zj[ii, :].real, w[::-1], 'same')
        zji[ii, :] = np.convolve(zj[ii, :].imag, w[::-1], 'same')
    # recover the static apparent resistivity from reference freq 
    zjc = zjr + 1j * zji 
    rc = z2rhoa (zjc, freqs) #np.abs(zjc)**2 / (omega0[:, None] * mu0 )
    
    if mode =='same': 
        rc[:, 0]= res2d[:, 0]
        zjc[:, 0]= zj [:, 0]
    
    return zjc if out =='z' else rc 
    
     
def moving_average (y, /, window_size = 3 , method ='sma',
                    mode ='same', alpha =.5 ): 
    """ A moving average is  used with time series data to smooth out
    short-term fluctuations and highlight longer-term trends or cycles.
    
    Funtion analyzes data points by creating a series of averages of different
    subsets of the full data set. 
    
    Parameters 
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
        
    window_size : int
        the length of the window. Must be greater than 1 and preferably
        an odd integer number.Default is ``3``
        
    method: str 
        variant of moving-average. Can be ``sma``, ``cma``, ``wma`` and ``ema`` 
        for simple, cummulative, weight and exponential moving average. Default 
        is ``wma``. 
        
    mode: str
        returns the convolution at each point of overlap, with an output shape
        of (N+M-1,). At the end-points of the convolution, the signals do not 
        overlap completely, and boundary effects may be seen. Can be ``full``,
        ``same`` and ``valid``. See :doc:`~np.convole` for more details. Default 
        is ``same``. 
        
    alpha: float, 
        smoothing factor. Only uses in exponential moving-average. Default is 
        ``.5``.
    
    Returns 
    --------
    ya: array like, shape (N,) 
        Averaged time history of the signal
    
    Notes 
    -------
    The first element of the moving average is obtained by taking the average 
    of the initial fixed subset of the number series. Then the subset is
    modified by "shifting forward"; that is, excluding the first number of the
    series and including the next value in the subset.
    
    Examples
    --------- 
    >>> import numpy as np ; import matplotlib.pyplot as plt 
    >>> from pycsamt.processing.ctools  import moving_average 
    >>> data = np.random.randn (37) 
    >>> # add gaussion noise to the data 
    >>> data = 2 * np.sin( data)  + np.random.normal (0, 1 , len(data))
    >>> window = 5  # fixed size to 5 
    >>> sma = moving_average(data, window) 
    >>> cma = moving_average(data, window, method ='cma' )
    >>> wma = moving_average(data, window, method ='wma' )
    >>> ema = moving_average(data, window, method ='ema' , alpha =0.6)
    >>> x = np.arange(len(data))
    >>> plt.plot (x, data, 'o', x, sma , 'ok--', x, cma, 'g-.', x, wma, 'b:')
    >>> plt.legend (['data', 'sma', 'cma', 'wma'])
    
    References 
    ----------
    .. * [1] https://en.wikipedia.org/wiki/Moving_average
    .. * [2] https://www.sciencedirect.com/topics/engineering/hanning-window
    .. * [3] https://stackoverflow.com/questions/12816011/weighted-moving-average-with-numpy-convolve
    
    """
    y = np.array(y)
    try:
        window_size = np.abs(_assert_all_types(int(window_size), int))
    except ValueError:
        raise ValueError("window_size has to be of type int")
    if window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if  window_size > len(y):
        raise TypeError("window_size is too large for averaging"
                        f"Window must be greater than 0 and less than {len(y)}")
    
    method =str(method).lower().strip().replace ('-', ' ') 
    
    if method in ('simple moving average',
                  'simple', 'sma'): 
        method = 'sma' 
    elif method  in ('cumulative average', 
                     'cumulative', 'cma'): 
        method ='cma' 
    elif method  in ('weighted moving average',
                     'weight', 'wma'): 
        method = 'wma'
    elif method in('exponential moving average',
                   'exponential', 'ema'):
        method = 'ema'
    else : 
        raise ValueError ("Variant average methods only includes "
                          f" {smart_format(['sma', 'cma', 'wma', 'ema'], 'or')}")
    if  1. <= alpha <= 0 : 
        raise ValueError ('alpha should be less than 1. and greater than 0. ')
        
    if method =='sma': 
        ya = np.convolve(y , np.ones (window_size), mode ) / window_size 
        
    if method =='cma': 
        y = np.cumsum (y) 
        ya = np.array([ y[ii]/ len(y[:ii +1]) for ii in range(len(y))]) 
        
    if method =='wma': 
        w = np.cumsum(np.ones(window_size, dtype = float))
        w /= np.sum(w)
        ya = np.convolve(y, w[::-1], mode ) #/window_size
        
    if method =='ema': 
        ya = np.array ([y[0]]) 
        for ii in range(1, len(y)): 
            v = y[ii] * alpha + ( 1- alpha ) * ya[-1]
            ya = np.append(ya, v)
            
    return ya 

def savitzky_golay1d (y, window_size, order, deriv=0, rate=1, mode='same'):
    r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    
    The Savitzky-Golay filter removes high frequency noise from data. It has the 
    advantage of preserving the original shape and features of the signal better
    than other types of filtering approaches, such as moving averages techniques.
    
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    mode: str 
         mode of the border prepending. Should be ``valid`` or ``same``. 
         ``same`` is used for prepending or appending the first value of
         array for smoothing.Default is ``same``.  
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly suited for 
    smoothing noisy data. The main idea behind this approach is to make for 
    each point a least-square fit with a polynomial of high order over a
    odd-sized window centered at the point.
    
    Examples
    --------
    >>> import numpy as np 
    >>> import matplotlib.pyplot as plt 
    >>> from pycsamt.processing import savitzky_golay1d 
    >>> t = np.linspace(-4, 4, 500)
    >>> y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    >>> ysg = savitzky_golay1d(y, window_size=31, order=4)
    >>> plt.plot(t, y, label='Noisy signal')
    >>> plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    >>> plt.plot(t, ysg, 'r', label='Filtered signal')
    >>> plt.legend()
    >>> plt.show()
    
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    .. [3] https://en.wikipedia.org/wiki/Savitzky%E2%80%93Golay_filter#Moving_average
    
    """

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode=mode)
  
def skew(ediObj, method='swift'): 
    r"""
    The conventional asymmetry parameter based on the Z magnitude. 
    
    Parameters 
    ---------
    ediObj: pycsamt.core.edi.Edi or mtpy.core.edi.Edi 
        EDi object with full impedance tensor Z. 
    
    method: str 
        Kind of correction. Can be ``swift`` for the remove distorsion proposed 
        by Swift in 1967. The value close to 0. assume the 1D and 2D structures 
        and 3D otherwise. Conversly to ``bahr`` for the remove distorsion proposed  
        by Bahr in 1991. The latter threshold is set to 0.3. Above this value 
        the structures is 3D. 
    
    Returns 
    ------- 
    skw, mu : Tuple of array-like , shape (N, )
        - Array of skew at each frequency 
        - rotational invariant ``mu`` at each frequency. 
        
    See also 
    -------- 
    
    The |EM| signal is influenced by several factors such as the dimensionality
    of the propagation medium and the physical anomalies, which can distort the
    |EM| field both locally and regionally. The distortion of Z was determined 
    from the quantification of its asymmetry and the deviation from the conditions 
    that define its dimensionality. The parameters used for this purpose are all 
    rotational invariant because the Z components involved in its definition are
    independent of the orientation system used. The conventional asymmetry
    parameter based on the Z magnitude is the skew defined by Swift (1967) as
    follows:
    
    .. math:: skew_{swift}= |Z_{xx} + Z_{yy} \frac{ Z_{xy} - Z-{yx}}| 
        
    When the :math:`skew_{swift}`  is close to ``0.``, we assume a 1D or 2D model
    when the :math:`skew_{swift}` is greater than ``>=0.2``, we assume 3D local 
    anomaly (Bahr, 1991; Reddy et al., 1977).
    
    Furthermore, Bahr (1988) proposed the phase sensitive skew which calculates
    the skew taking into account the distortions produced in Z over 2D structures
    by shallow conductive anomalies and is defined as follows:
    
    .. math::
        
        skew_{Bahr} & = & \sqrt{ |[D_1, S_2] -[S_1, D_2]|} \frac{|D_2|} \quad \text{where} 
        
        S_1 & = & Z_{xx} + Z_{yy} \quad ; \quad  S_2 = Z_{xy} + Z_{yx} 
        D_1 & = &  Z_{xx} - Z_{yy} \quad ; \quad  D_2 = Z_{xy} - Z_{yx}
        
    Note that The phase differences between two complex numbers :math:`C_1` and 
    :math:`C_2` and the corresponding amplitude  products are now abbreviated 
    by the commutators:
        
    ..math:: 
      
        [C_1, C_2] & = & Im (C_2 C_1 ^{*})
                   & = & R_e(C_1) Im(C_2) - R_e(C_2) Im(C_1)
                    
    Indeed, :math:`skew_{Bahr}` measures the deviation from the symmetry condition
    through the phase differences between each pair of tensor elements,considering
    that phases are less sensitive to surface distortions(i.e. galvanic distortion).
    The :math:`skew_{Bahr}` threshold is set at ``0.3`` and higher values mean 
    3D structures (Bahr, 1991).
    
    
    References 
    ----------
        
    Bahr, K., 1991. Geological noise in magnetotelluric data: a classification 
        of distortion types. Physics of the Earth and Planetary Interiors 66
        (1–2), 24–38.
    Barcelona, H., Favetto, A., Peri, V.G., Pomposiello, C., Ungarelli, C., 2013.
        The potential of audiomagnetotellurics in the study of geothermal fields: 
        A case study from the northern segment of the La Candelaria Range,
        northwestern Argentina. J. Appl. Geophys. 88, 83–93.
        https://doi.org/10.1016/j.jappgeo.2012.10.004   
        
    Swift, C., 1967. A magnetotelluric investigation of an electrical conductivity 
       anomaly in the southwestern United States. Ph.D. Thesis, MIT Press. Cambridge. 
       
       
    """
    method = str(method).lower().strip() 
    
    if method not in ('swift', 'bahr'): 
        raise ValueError(f'Expected argument ``swift`` or ``bahr`` not: {method!r}')
        
    ediObj = _assert_edi_obj(ediObj)
    
    Zxx = ediObj.Z.z[:, 0, 0]; Zyy = ediObj.Z.z[:, 1, 1] 
    Zxy = ediObj.Z.z[:, 0, 1]; Zyx= ediObj.Z.z[:, 1, 0]
    
    S1 =Zxx + Zyy; S2 = Zxy + Zyx; D1 =Zxx-Zyy ;  D2= Zxy-Zyx 
    D1S2 = (S2 * np.conj(D1)).imag ; S1D2 = (D2 * np.conj(S1)).imag 
    
    if method =='swift': 
        skw = np.abs ( S1  / D2 )
    else : 
        skw = np.sqrt(np.abs( D1S2 - S1D2))/np.abs(D2)
        
    mu = np.sqrt(np.abs(D1S2) + np.abs (S1D2))/ np.abs(D2) 
        
    return skw, mu

        
def remove_static_shift( ediObj, ss_x=1., ss_y=1., **kws):
    r"""
    Remove static shift from the apparent resistivity
    
    Assume the original observed tensor Z is built by a static shift S
    and an unperturbated "correct" Z0 such as :math:`Z = S * Z0` therefore the 
    correct Z will be :math:` Z0 = S^(-1) * Z`.
        
    Parameters 
    ---------- 
    ediObj: pycsamt.core.edi.Edi or mtpy.core.edi.Edi 
        EDi object with full impedance tensor Z. 
        
    ss_x: float 
        correction factor for x component


    ss_y: float 
        correction factor for y component

    Returns 
    ------- 
    pycsamt.core.z.Z object -  new Z object with static shift removed


    Notes
    ------
    The factors are in resistivity scale, so the entries of  the matrix "S"
    need to be given by their square-roots!

    Example
    -------

    >>> import pycsamt.processing as Processing 
    >>> edifile = '/Users/Daniel/Desktop/Data/AMT/E1/di_test/new_csa00.edi'
    >>> outputedi = 'rmss_csa00.edi'
    >>> Processing.remove_static_shift(
        edi_fn = edifile, ss_x= .5, ss_y=1.2,new_edi_fn = outputedi
                                       )
    """
    ediObj =_assert_edi_obj(ediObj)
    
    if hasattr(ediObj, 'edifile'): 
        ediObj= mtpy.core.mt.MT(fn =getattr(ediObj, 'edifile'))
    
    new_z_obj = ediObj.remove_static_shift(ss_x=ss_x , ss_y=ss_y)

    return new_z_obj  
    
def remove_distortion( ediObj , num_freq=None, **kws):
    """
    remove distortion following Bibby et al. [2005].

    :param num_freq: number of frequencies to look for distortion from the
                     highest frequency
    :type num_freq: int

    :returns: Distortion matrix
    :rtype: np.ndarray(2, 2, dtype=real)

    :returns: Z with distortion removed
    :rtype: pycsamt.core.z.Z

    :Remove distortion and write new .edi file: ::

        >>> import pycsamt.processing as Processing 
        >>> edifile = '/Users/Daniel/Desktop/Data/AMT/E1/di_test/new_csa00.edi'
        >>> outputedi = 'rmss_csa00.edi'
        >>> Processing.remove_distortion(edi_fn = edifile, new_edi_fn= outputedi
            )

    """
    ediObj =_assert_edi_obj(ediObj)
    
    if hasattr('edifile', ediObj):
        ediObj= mtpy.core.mt.MT(fn =getattr(ediObj, 'edifile'))
        
    D, new_z = ediObj.remove_distortion(**kws)
    
    return D, new_z   
    
def removeNoise(edi_fn = None, kind='sim', **kws): 
    """ Remove multiple noises in EDIs and save to new files. Noise can be 
    either a `staticshift` , `distorsion` noises or other interferences. 
    
    For Electromagnetic Array Profiling (EMAP) data correction, may refer to 
    :meth:`pycsamt.processing.Processing.correct_edi`. The Remove 
    distortion following Bibby et al. [2005]. Remove static shift from the 
    apparent resistivity assume the original observed tensor Z is built by
    a static shift S and an unperturbated 
    
         * Z = S * Z0
 
    therefore the correct Z will be :
        
        * Z0 = S^(-1) * Z
            
    :param edi_fn: Path-Like object. Full path to EDI-files. 
    :type edi_fn: str 
    
    :param kind: Type of noise to remove. Can be a static shift effect for 
        ``ss``, distorsion for ``dist``. The human activites can be removed 
        using the `pca`` and the other interferences can be removed using 
        the simplier filters ``sim``. Default is ``sim``. 
        
    :type kind: str 
    
    :param ss_x: correction factor for x component
    :type ss_x: float

    :param ss_y: correction factor for y component
    :type ss_y: float
    
    :param num_freq: number of frequencies to look for distortion from the
                     highest frequency
    :type num_freq: int

    :returns: Distortion matrix
    :rtype: np.ndarray(2, 2, dtype=real)

    :returns: 
        - new Z object with static shift removed
        - Distortion matrix
        
    :rtype: pycsamt.core.z.Z
    
    .. note:: The factors are in resistivity scale, so the
              entries of  the matrix "S" need to be given by their
              square-roots!
              
    :Example:
        >>> from pycsamt.processing import Processing 
        >>> edipath = '/Users/Desktop/ediout/'
        >>> Processing.removeNoises(edi_fn =edipath)
        
    :See also: 
        
        - MTpy of Alison.Kirkby@ga.gov.au via https://github.com/MTgeophysics/mtpy
        
    """
    kind =str(kind).lower() 
    if kind in ('ss', 'staticshift'): kind =='ss'
    elif kind in ('dt', 'dist', 'distortion'): kind =='dt'
    elif kind in  ('principal component analysis', 'pca', 'sklearn'): 
        kind ='pca'
    elif kind in ('simple', 's0', 'simplier', 'light', 'base'): 
        kind ='sim' 
        
    else : 
        raise ValueError(f"Wrong Param `kind`: {kind}. "
                         "Should be `ss` or `dist`"
                         )
    if os.path.isfile(edi_fn): 
        edi_fn = os.path.dirname (edi_fn) # keep only the path 
    if not os.path.dirname(edi_fn): 
        raise ValueError(f'Wrong given EDI path: {edi_fn}')
        
    D= None 
    # show progress bar 
    if itqdm : 
        pbar =tqdm.tqdm(total= len(os.listdir(edi_fn)),ascii=True,unit='B',
                         desc ='WEgeophysics-pycsamt[---> NoiseRemoval]', 
                         ncols =77)
    for k, edi in enumerate (os.listdir(edi_fn)): 
        edifile = os.path.join(edi_fn, edi)
        with warnings.catch_warnings(): # ignore multiple warnings 
            warnings.simplefilter('ignore')
            if kind =='ss': 
                new_z = remove_static_shift(edifile, **kws)
            elif kind =='dt': 
                D, new_z = remove_distortion (edifile, **kws)
            elif kind =='pca': 
                new_z = pca_filter(edifile, **kws)
            elif kind =='sim': 
                new_z = simpler_filter(edifile, **kws)
        # show the progress bar ;
        if itqdm :
            pbar.update(k)
    # close the progress bar
    if itqdm :
        pbar.close()

    print(' completed' if itqdm else '--- process completed ----')
    
    return new_z, D 
    

def pca_filter (ediObj, var=.95 , **kws): 
    """ Sanitize EDIs data and remove the outliers using the principal 
    component analysis (PCA).
    
    In considering the human activity noises removal using the PCA
    
    :param ediObj: Path-Like object. Full path to EDI-file. 
    :type edi_fn: str 
    
    :param var: component or ratio of moise to remove.  Default value is 
        ``.95``. 
    :type var: float
    
    :param kws: Additional keywords arguments from
        :mod:`sklearn.decomposition.PCA` 
    :type kws: dict 
    
    :returns: New Z impedance object with remove outliers. 
    :rtype: pycsamt.core.z.Z
    
    :Example: 
    >>> from pycsamt.processing import Processing 
    >>> from pycsamt.core.edi import Edi 
    >>> edifile = '/home/edi/m20.E000.edi'
    >>> newZ =Processing.pca_filter ( edifile, var =.80 ) 
    >>> #write a new corrected edifile 
    >>> _= Edi(edifile).write_new_edifile(new_Z= newZ, 
    ...               new_edi_fn ='m20.removeO.E000.edi')
        
    """
 
    def reshape_and_fit_z (z , back =False , fit= 'transform') :
        """ Reshape z for to np.ndarray(nfreq, 1) for scikit-learn API  
        back for filling the z ndarray (nfreq, 2, 2)"""
        z = z.reshape ((1, z.shape [0])) if back else z.reshape (
            (z.shape [1], 1)) 
         
        z = (pca.fit_transform(z) if fit =='transform' else 
             pca.inverse_transform(z) ) if fit is not None else z 
        return z 


    if isinstance (var, str):
        try : var =float(var) 
        except: 
            if '%' in var: 
                var = float(var.replace('%', ''))*1e-2
            else: 
                raise ValueError('Variance could be a float '
                                 f'value not: {type(var).__name__!r}')
    if var > 1.: 
        raise ValueError('Variance ratio should be '
                         f'less than 1 :{str(var)!r}')

    #create an EDI OBJECT 
    ediObj = _assert_edi_obj(ediObj)
    # make a new object 
    new_Z = MTz.Z(z_array=np.zeros((ediObj.Z.freq.shape[0], 2, 2),
                                   dtype='complex'),
                  z_err_array=np.zeros((ediObj.Z.freq.shape[0], 2, 2)),
                  freq=ediObj.Z.freq)
   
    # construct the pca object from sklearn 
    pca = PCA(n_components= var, **kws) 
    # loop to remove outliers on the Z impedance object 
    for ii in range(2):
        for jj in range(2):
            # need to look out for zeros in the impedance
            # get the indicies of non-zero components
            nz_index = np.nonzero(ediObj.Z.z[:, ii, jj])
             
            if len(nz_index[0]) == 0:
                continue
            # get the non_zeros components 
            with np.errstate(all='ignore'):
                z_real = ediObj.Z.z[nz_index, ii, jj].real
                z_imag = ediObj.Z.z[nz_index, ii, jj].imag
                z_err = ediObj.Z.z_err[nz_index, ii, jj]

                z_real_t =reshape_and_fit_z(z_real) 
                z_imag_t = reshape_and_fit_z(z_imag) 
                z_err_t = reshape_and_fit_z(z_err) 
                
                z_real_b =reshape_and_fit_z(z_real_t,back=True,fit =None) 
                z_imag_b = reshape_and_fit_z(z_imag_t,back=True,fit =None) 
                z_err_b = reshape_and_fit_z(z_err_t,back=True,fit =None 
                                            ) 
           # set the new Z object 
            new_Z.z[nz_index, ii, jj] = z_real_b  + 1j * z_imag_b 
            new_Z.z_err[nz_index, ii, jj] = z_err_b
            
    # compute resistivity and phase for new Z object
    new_Z.compute_resistivity_phase()
 
    return new_Z

def simpler_filter (edi_fn, **kws): 
    """ Sanitize EDIs data and remove the outliers using simple 
    fitting function. 
    
    In considering the other interferences noises removal using the base 
    or simpler filter via fitting impedances values. 
    
    
    :param edi_fn: Path-Like object. Full path to EDI-file. 
    :type edi_fn: str 
    
    :param kws: additional keyword arguments for correct values 
    
    :Example: 
    >>> from pycsamt.processing import Processing 
    >>> from pycsamt.core.edi import Edi 
    >>> edifile = '/home/edi/m20.E000.edi'
    >>> newZ =Processing.simple_removal( edifile ) 
    >>> #write a new corrected edifile 
    >>> _= Edi(edifile).write_new_edifile(new_Z= newZ, 
    ...               new_edi_fn ='m20.removeSR.E000.edi')
    
    """
    
    edi_fn = _assert_all_types(edi_fn, str)
    if not os.path.isfile(edi_fn): 
        raise ValueError(f'Wrong given EDI file: {edi_fn}')

    #create an EDI OBJECT 
    ediObj = Edi(edi_fn )
    
    # make a new object 
    new_Z = MTz.Z(z_array=np.zeros((ediObj.Z.freq.shape[0], 2, 2),dtype='complex'),
                  z_err_array=np.zeros((ediObj.Z.freq.shape[0], 2, 2)),
                  freq=ediObj.Z.freq)

    # loop to correct the Z impedance object values 
    for ii in range(2):
        for jj in range(2):
            # need to look out for zeros in the impedance
            # get the indicies of non-zero components
            nz_index = np.nonzero(ediObj.Z.z[:, ii, jj])
    
            if len(nz_index[0]) == 0:
                continue
            # get the non_zeros components 
            with np.errstate(all='ignore'):
                z_real = reshape(
                    ediObj.Z.z[nz_index, ii, jj].real)
                z_imag = reshape(
                    ediObj.Z.z[nz_index, ii, jj].imag)
                z_err = reshape(
                    ediObj.Z.z_err[nz_index, ii, jj]) 
                # correct the values 
                z_real_c, *_ =moving_average (z_real, **kws) 
                z_imag_c, *_ = moving_average (z_imag, **kws) 
                z_err_c, *_ = moving_average (z_err, **kws) 
                
           # set the new Z object 
            new_Z.z[nz_index, ii, jj] = reshape(
                z_real_c, 1)   + 1j * reshape(z_imag_c, 1) 
            new_Z.z_err[nz_index, ii, jj] = reshape(z_err_c, 1)
            
    # compute resistivity and phase for new Z object
    new_Z.compute_resistivity_phase()
 
    return new_Z   

def export2newedis (ediObj, new_Z , savepath =None, **kws):
    """ Export new EDI files from the former object with  a given new impedance 
    tensors. 
    
    The export is assumed a new output EDI resulting from multiples corrections 
    applications. 
    
    Parameters 
    -----------
    ediObj: str or  pycsamt.core.edi.Edi or mtpy.core.edi.Edi 
        Full path to Edi file or object from `pyCSAMT`_ and `MTpy`_ packages 
    
    new_Z: ndarray (nfreq, 2, 2) 
        Ndarray of impendance tensors Z. The tensor Z is 3D array composed of 
        number of frequency `nfreq`and four components (``xx``, ``xy``, ``yx``,
        and ``yy``) in 2X2 matrices. The  tensor Z is a complex number. 
    
    Returns 
    --------
     ediObj from pycsamt.core.edi.Edi or mtpy.core.edi.Edi 
     
    """
    
    ediObj = _assert_edi_obj(ediObj)
    ediObj.write_new_edifile( new_Z=new_Z, **kws)
    return ediObj 


@nb.njit if inumba else donothing("Skip numba when the latter doesn't work "
                                  "with the current numpy.__version__")  
def restoreZ(ediObjs, *, buffer = None, kind='slinear',
                          method ='pd', **kws ): 
    """ Fix the weak and missing signal at the 'dead-band`- and recover the 
    missing impedance tensor values. 
    
    The function uses the complete frequency (frequency with clean data) collected 
    thoughout the survey to recover by inter/extrapolating the missing or weak 
    frequencies thereby restoring the impedance tensors at that 'dead-band'. Note 
    that the 'dead- band' also known as 'attenuation -band' is where the AMT 
    signal is weak or generally abscent. 

    Parameters 
    ---------- 
    ediObjs: list  of  pycsamt.core.edi.Edi or mtpy.core.edi.Edi objects 
        Collections of EDI-objects from `pyCSAMT`_ and `MTpy`_ packages 
        
    buffer: list [max, min] frequency in Hz
        list of maximum and minimum frequencies. It must contain only two values.
        If `None`, the max and min of the clean frequencies are selected. Moreover
        the [min, max] frequency should not compulsory fit the frequency range in 
        the data. The given frequency can be interpolated to match the best 
        closest frequencies in the data. 
    
    kind: str or int, optional
        Specifies the kind of interpolation as a string or as an integer 
        specifying the order of the spline interpolator to use. The string 
        has to be one of ``linear``, ``nearest``, ``nearest-up``, ``zero``, 
        ``slinear``,``quadratic``, ``cubic``, ``previous``, or ``next``. 
        ``zero``, ``slinear``, ``quadratic``and ``cubic`` refer to a spline 
        interpolation of zeroth, first, second or third order; ``previous`` 
        and ``next`` simply return the previous or next value of the point; 
        ``nearest-up`` and ``nearest`` differ when interpolating half-integers 
        (e.g. 0.5, 1.5) in that ``nearest-up`` rounds up and ``nearest`` rounds 
        down. If `method` param is set to ``pd`` which refers to pd.interpolate 
        method , `kind` can be set to ``polynomial`` or ``pad`` interpolation. 
        Note that the polynomial requires you to specify an `order` while 
        ``pad`` requires to specify the `limit`. Default is ``slinear``.
        
    method: str, optional  
        Method of interpolation. Can be ``base`` for `scipy.interpolate.interp1d`
        ``mean`` or ``bff`` for scaling methods and ``pd``for pandas interpolation 
        methods. Note that the first method is fast and efficient when the number 
        of NaN in the array if relatively few. It is less accurate to use the 
        `base` interpolation when the data is composed of many missing values.
        Alternatively, the scaled method(the  second one) is proposed to be the 
        alternative way more efficient. Indeed, when ``mean`` argument is set, 
        function replaces the NaN values by the nonzeros in the raw array and 
        then uses the mean to fit the data. The result of fitting creates a smooth 
        curve where the index of each NaN in the raw array is replaced by its 
        corresponding values in the fit results. The same approach is used for
        ``bff`` method. Conversely, rather than averaging the nonzeros values, 
        it uses the backward and forward strategy  to fill the NaN before scaling.
        ``mean`` and ``bff`` are more efficient when the data are composed of 
        lot of missing values. When the interpolation `method` is set to `pd`, 
        function uses the pandas interpolation but ended the interpolation with 
        forward/backward NaN filling since the interpolation with pandas does
        not deal with all NaN at the begining or at the end of the array. Default 
        is ``pd``.
        
    fill_value: array-like or ``extrapolate``, optional
        If a ndarray (or float), this value will be used to fill in for requested
        points outside of the data range. If not provided, then the default is
        NaN. The array-like must broadcast properly to the dimensions of the 
        non-interpolation axes.
        If a two-element tuple, then the first element is used as a fill value
        for x_new < x[0] and the second element is used for x_new > x[-1]. 
        Anything that is not a 2-element tuple (e.g., list or ndarray,
        regardless of shape) is taken to be a single array-like argument meant 
        to be used for both bounds as below, above = fill_value, fill_value.
        Using a two-element tuple or ndarray requires bounds_error=False.
        Default is ``extrapolate``. 
        
    kws: dict 
        Additional keyword arguments from :func:`~interpolate1d`. 
    
    Returns 
    --------
        Array-like of pycsamt.core.z.Z objects 
        Array collection of new Z impedances objects with dead-band tensor 
        recovered. :class:`pycsamt.core.z.Z` are ndarray (nfreq, 2, 2). 
        2x2 matrices for components xx, xy and yx, yy. 

    See More  
    ---------
    One main problem in collecting |NSAMT| data is the signal level in the 
    'attenuation band'. Compared to the |CSAMT| method (Wang and Tan, 2017; 
    Zonge and Hughes, 1991),the natural signals are not under our control and 
    suffer from frequency  ranges with little or no signal.  Most notably, the 
    |NSAMT| 'dead-band' between approximately 1 kHz and 4 kHz, but also a signal 
    low in the vicinityof 1 Hz where the transition to magnetospheric energy 
    sources occurs (Goldak and Olson, 2015). In this band, natural source signals
    are generally  absent. The EM energy is dissipated and often cultural |EM| 
    noise fills the gap (Zonge, 2000). The response is extrapolated from results 
    observed top frequencies( For instance at 20, 40, 250, and 500 Hz).Experience
    indicates that the natural source signal level at 2000 Hz can be expected 
    to approach 10-6 γ/√Hz (Zheng, 2010; Zonge, 2000).

    References 
    ----------
    Goldak, D.K., Olson, R.W., 2015. New developments in |AMT| exploration :
        Case study from Darnley Bay. CSEG Rec. 22–27.
    Wang, K., Tan, H., 2017. Research on the forward modeling of |CSAMT| in 
        three-dimensional axial anisotropic media. J. Appl. Geophys. 146, 27–36.
        https://doi.org/10.1016/j.jappgeo.2017.08.007
    Zonge, I., 2000. |NSAMT| Imaging. 3322 East Fort Lowell Road, Tucson, AZ 85716 USA. 
    Zonge, L., Hughes, L.J., 1991. |CSAMT|. Soc. Explor. Geophys. 2, 713–809.
       
    Examples 
    --------
    >>> import numpy as np 
    >>> import matplotlib.pyplot as plt 
    >>> from pycsamt.core import Edi_collection 
    >>> from pycsamt.utils import reshape 
    >>> from pycsamt.processing import restoreZ
    >>> from pycsamt.processing import fit_tensor, moving_average
    >>> from pycsamt.processing import control_freq_buffer , export2newedis 
    >>> path2edi = 'data/3edis'
    >>> colObjs = Edi_collection (path2edi)
    >>> ediObjs = colObjs.ediObjs 
    >>> # One can specify the frequency buffer like the example below, However 
    >>> # it is not necessaray at least there is a a specific reason to fix the frequencies 
    >>> buffer = [1.45000e+04,1.11500e+01]
    >>> zobjs_b = restoreZ(ediObjs, buffer = buffer) # with buffer 
    >>> zobjs = restoreZ(ediObjs ) # with no buffer 
    >>> # Now the idea is to compare the first raw z tensor object with its 
    >>> # corresponding restored ( the component |XY| for illustration)
    >>> # export the first restored  |Zxy| object for comparison
    >>> zxy_restored = np.abs (zobjs[0].z[:, 0, 1]) 
    >>> # Export the first raw object with missing Z at the dead dand in ediObjs collection
    >>> z1 = np.abs(ediObjs[0].Z.z) 
    >>> z1freq= ediObjs[0].Z._freq # the frequency of the first z obj 
    >>> # get the frequency of the clean data knonw as reference frequency
    >>> indice_reffreq = np.argmax (list (map(lambda o: len(o.Z._freq), ediObjs)))
    >>> reffreq = ediObjs [indice_reffreq].Z._freq 
    >>> # use the real part of component xy for the test 
    >>> zxy = np.abs (z1[:, 0, 1])  
    >>> # fit zxy to get the missing data in the dead band 
    >>> zfit = fit_tensor(refreq= reffreq, compfreq= z1freq, z=zxy)
    >>> # if one uses a buffer, we can check the interpolated frequency buffer 
    >>> # in the reference frequency 
    >>> control_freq_buffer (reffreq, buffer) 
    ...  array([1.470e+04, 1.125e+01])
    >>> # now we can set buffer to  [1.470e+04, 1.125e+01] and use the value 
    >>> # to find the index in raw zxy for plotting purpose to between in the frequency range 
    >>> # mask the z value out of the buffer frequency range
    >>> ix_buf = reshape(np.argwhere (np.isin (reffreq, buffer ))) 
    >>> # buffered indexes stands for array([ 9, 50], dtype=int64)
    >>> # slice the buffered frequency and mask the useless frequency (out of buff range)
    >>> mask = ~np.isin (reffreq , reffreq[9 : 51] ) # exclude 51e index (ix_buf[1] +1)   
    >>> reffreq_buf = np.ma.masked_where (mask , reffreq) 
    >>> zxy_restored_buf = np.ma.masked_where (mask ,zxy_restored)   
    >>> zfit_buf = np.ma.masked_where (mask ,zfit)  
    >>> # not necessary, one can corrected z to get a smooth resistivity distribution 
    >>> zcorrected , *_ = moving_average (zxy_restored)                     
    >>> # plot the two figures 
    >>> plt.figure(figsize =(10, 5))
    >>> plt.loglog(reffreq, zfit, '^r', reffreq, zxy_restored, 'ok--')
    >>> plt.loglog(reffreq_buf, zfit_buf, 'g.', reffreq_buf, zxy_restored, '<b:')
    >>> plt.loglog( reffreq, zcorrected, '1-.')
    >>> plt.legend (['fit data', 'restored', 'Buffered fit data',
                     'Buffered restored', 'Corrected data' ], loc ='best')
    >>> plt.xlabel ('$Log_{10} frequency [H_z]$') ; plt.ylabel('$ Resistivity [ \Omega.m]$')
    >>> plt.title ('Recovered tensor $|Z_{xy}|$')
    >>> plt.grid (visible =True , alpha =0.8, which ='both', color ='k')
    >>> plt.tight_layout()
    >>> # As the user can see, Zxy is restored at  freq < 10^1 and  >10^4 Hz 
    >>> # write a z restored object in new Edi-files 
    >>> export2newedis (new_Z = zxy_restored, new_edi_fn = '')

    """
    def z_transform (z , rfq, fq,  slice_= None): 
        """ Route to do the same task for real, imaginary and error """
        with np.errstate(all='ignore'):
            z = reshape(z) 
            z = fit_tensor(compfreq= fq, refreq =rfq, z = z  ) 
            z = interpolate1d(arr=z , kind = kind, method = method, **kws )
        return z [slice_] 
        
    #buffer = control_freq_buffer(freq_, buffer =[5.70e7, 2e1])
    if not isinstance( ediObjs, (list, tuple, np.ndarray)): 
        if not hasattr (ediObjs, '__dict__'): 
            raise ValueError ('Object should be an instance or a class'
                              f' not: {type(ediObjs).__name__!r}')
        ediObjs =np.array([ediObjs],dtype =object) 
        
    # get the frequencies obj 
    zObjs = np.array (list(map(lambda o: o.Z, ediObjs)) , dtype =object) 
    #read all frequency length and take the max frequency 
    # known  as the complete frequencies ( That compose all values)
    freqsize = np.array (list (map (lambda o:len(o._freq), zObjs)))
    ix_max  = np.argmax(freqsize)
    # get the complete freq 
    cfreq = zObjs[ix_max]._freq  
    
    # control the buffer and get the the range of frequency 
    buffer = control_freq_buffer(cfreq, buffer)
    ix_buf,  = np.where ( np.isin (cfreq, buffer)) 
    ## index for slice the array in the case the buffer is set 
    ix_s , ix_end = ix_buf ; ix_end += 1 ; slice_= slice (ix_s,  ix_end) 
    s_cfreq = cfreq [slice_] # slice frequency within the buffer 
    
    # make a new Z objects 
    # make a new object 
    new_zObjs =np.zeros_like (zObjs, dtype =object )
    # loop to correct the Z impedance object values 
    for kk, ediObj in enumerate (ediObjs):
        new_Z = MTz.Z(z_array=np.zeros((len(s_cfreq), 2, 2),dtype='complex'),
                    z_err_array=np.zeros((len(s_cfreq), 2, 2)),
                    freq=s_cfreq)
        
        for ii in range(2):
            for jj in range(2):
                # need to look out for zeros in the impedance
                # get the indicies of non-zero components
                nz_index = np.nonzero(ediObj.Z.z[:, ii, jj])
                if len(nz_index[0]) == 0:
                    continue
                # get the non_zeros components and interpolate 
                # frequency to recover the component in dead-band frequencies 
                # Use the whole edi
                with np.errstate(all='ignore'):
                    zfreq = ediObj.Z._freq
                    z_real = reshape(ediObj.Z.z[nz_index, ii, jj].real) 
                    z_real = z_transform (z_real, rfq=cfreq, fq=zfreq, 
                                          slice_=slice_ 
                                          )
                    z_imag = reshape(ediObj.Z.z[nz_index, ii, jj].imag) 
                    z_imag = z_transform (z_imag, rfq=cfreq, fq=zfreq, 
                                          slice_=slice_ 
                                          )
                    z_err = reshape(ediObj.Z.z_err[nz_index, ii, jj]) 
                    z_err = z_transform (z_err, rfq=cfreq, fq=zfreq,
                                         slice_=slice_ 
                                         )
                # Use the new dimension of the z and slice z according 
                # the buffer range. make the new index start at 0. 
                new_nz_index = slice (
                    * np.array([ix_s, ix_end],dtype=np.int32)-ix_s)
       
                new_Z.z[new_nz_index, ii, jj] = reshape(
                    z_real, 1)   + 1j * reshape(z_imag, 1) 
                new_Z.z_err[new_nz_index, ii, jj] = reshape(z_err, 1)
            
        # compute resistivity and phase for new Z object
        new_Z.compute_resistivity_phase()
        new_zObjs[kk] = new_Z 
        
    return new_zObjs 


def control_freq_buffer (freq, buffer = None ) :
    """ Assert the frequency buffer and find the nearest value if the 
    value of the buffer is not in frequency ranges .
    
    :param freq: array-like of frequencies 
    :param buffer: list of maximum and minimum frequency. It should contains 
        only two values. If `None`, the max and min frequencies are selected 
    :returns: Buffer frequency range 
    
    :Example: 
    >>> import numpy as np 
    >>> from pycsamt.processing.ctools import control_freq_buffer
    >>> freq_ = np.linspace(7e7, 1e0, 20) # 20 frequencies as reference
    >>> buffer = control_freq_buffer(freq_, buffer =[5.70e7, 2e1])
    >>> freq_ 
    ... array([7.00000000e+07, 6.63157895e+07, 6.26315791e+07, 5.89473686e+07,
           5.52631581e+07, 5.15789476e+07, 4.78947372e+07, 4.42105267e+07,
           4.05263162e+07, 3.68421057e+07, 3.31578953e+07, 2.94736848e+07,
           2.57894743e+07, 2.21052638e+07, 1.84210534e+07, 1.47368429e+07,
           1.10526324e+07, 7.36842195e+06, 3.68421147e+06, 1.00000000e+00])
    >>> buffer 
    ... array([5.52631581e+07, 1.00000000e+00])
    
    """
    
    if buffer is not None: 
        if np.iterable(buffer): 
            if 1 < len(buffer) > 2 :
                raise ValueError('Frequency buffer expects two values [max, min].'
                                 f' But {"is" if len(buffer)==1 else "are"} given ')
            if len(set(buffer))==1: 
                raise ValueError('Expect distinct values [max, min].'
                                 f'{str(buffer[0])!r} is given in twice.')
                
            for i, v in enumerate(buffer) : 
                
                if str(v).lower() =='min': 
                    buffer[i] = freq.min() 
                elif str(v).lower() =='max': 
                    buffer[i]= freq.max() 
                elif isinstance(v, str): 
                    raise ValueError(f"Expect 'min' or 'max' not: {v!r}")
                # Find the absolute difference with each value   
                # Get the index of the smallest absolute difference
                arr_diff = np.abs(freq - v )
                buffer[i] = freq[arr_diff.argmin()]
            
            buffer = np.array(buffer) 
            buffer.sort() ;  buffer = buffer [::-1]
            
    if buffer is None: 
        buffer = np.array([freq.max(), freq.min()])
        
    if buffer.min() < freq.min(): 
        raise ValueError(
            f'Given value {round(buffer.min(), 4) } is out of frequency range.'
            f' Expect a frequency greater than or equal to {round(freq.min(), 4)}'
            )
        
    if buffer.max() > freq.max() : 
        raise ValueError(
            f'Given value {round(buffer.max(),4)} is out of frequency range.'
            f' Expect a frequency less than or equal {round(freq.max(),4)}')
        
    return buffer 

            
def fit_tensor(refreq, compfreq , z , fill_value = np.nan): 
    """ Fit each tensor component to the complete frequency range. 
    
    The complete frequency is the frequency with clean data. It contain all the 
    frequency range on the site. During the survey since, the missing frequencies 
    lead to missing tensor data. So the function will indicate when the tensor 
    data is missing before any interpolation. 
    
    :param refreq: Reference frequency - Should be the complete frequency 
        collected in the field. 
        
    :param comfreq: array-like, should the frequency of the survey area.
    
    :param z: array-like, should be the  tensor value (real or imaginary part ) 
        at the component  xx, xy, yx, yy. 
        
    :param fill_value: float 
        Value to replace the missing data in tensors. Default is ``NaN``. 
        
    :param return: new Z filled by invalid value `NaN` where the frequency is 
        missing in the data. 
    
    Example::
    >>> import numpy as np 
    >>> from pycsamt.processing.ctools import fit_tensor
    >>> refreq = np.linspace(7e7, 1e0, 20) # 20 frequencies as reference
    >>> freq_ = np.hstack ((refreq.copy()[:7], refreq.copy()[12:] )) 
    >>> z = np.random.randn(len(freq_)) *10 # assume length of  freq as 
    >>>                 # the same like the tensor Z value 
    >>> zn  = fit_tensor (refreq, freq_, z)
    >>> z # some frequency values are missing but not visible. 
    ...array([-23.23448367,   2.93185982,  10.81194723, -12.46326732,
             1.57312908,   7.23926576, -14.65645799,   9.85956253,
             3.96269863, -10.38325124,  -4.29739755,  -8.2591703 ,
            21.7930423 ,   0.21709129,   4.07815217])
    >>> # zn show where the frequencies are missing  
    >>> # the NaN value means in a missing value in  tensor Z at specific frequency  
    >>> zn 
    ... array([-23.23448367,   2.93185982,  10.81194723, -12.46326732,
             1.57312908,   7.23926576, -14.65645799,          nan,
                    nan,          nan,          nan,          nan,
             9.85956253,   3.96269863, -10.38325124,  -4.29739755,
            -8.2591703 ,  21.7930423 ,   0.21709129,   4.07815217])
    >>> # let visualize where the missing frequency value in tensor Z 
    >>> refreq 
    ... array([7.00000000e+07, 6.63157895e+07, 6.26315791e+07, 5.89473686e+07,
           5.52631581e+07, 5.15789476e+07, 4.78947372e+07, 4.42105267e+07*,
           4.05263162e+07*, 3.68421057e+07*, 3.31578953e+07*, 2.94736848e+07*,
           2.57894743e+07, 2.21052638e+07, 1.84210534e+07, 1.47368429e+07,
           1.10526324e+07, 7.36842195e+06, 3.68421147e+06, 1.00000000e+00])
    >>> refreq[np.isnan(zn)] #we can see the missing value between [7:12](*) in refreq 
    ... array([44210526.68421052, 40526316.21052632, 36842105.73684211,
           33157895.2631579 , 29473684.78947368])
    
    """
    
    freqn, mask = ismissing(refarr= refreq , arr =compfreq, return_index='mask',
                            fill_value = fill_value)
    
    #mask_isin = np.isin(refreq, compfreq)
    z_new = np.full_like(freqn,fill_value = fill_value) 
    z_new[mask] = reshape(z) 
    
    return z_new 
    
def interpolate1d (arr, kind = 'slinear', method='mean', order = None, 
                   fill_value ='extrapolate',limit =None, **kws) :
    """ Interpolate array containing invalida values `NaN`
    
    Usefull function to interpolate the missing frequency values in the 
    tensor components. 
    
    Parameters 
    ----------
    arr: array_like 
        Array to interpolate containg invalid values. The invalid value here 
        is `NaN`. 
        
    kind: str or int, optional
        Specifies the kind of interpolation as a string or as an integer 
        specifying the order of the spline interpolator to use. The string 
        has to be one of ``linear``, ``nearest``, ``nearest-up``, ``zero``, 
        ``slinear``,``quadratic``, ``cubic``, ``previous``, or ``next``. 
        ``zero``, ``slinear``, ``quadratic``and ``cubic`` refer to a spline 
        interpolation of zeroth, first, second or third order; ``previous`` 
        and ``next`` simply return the previous or next value of the point; 
        ``nearest-up`` and ``nearest`` differ when interpolating half-integers 
        (e.g. 0.5, 1.5) in that ``nearest-up`` rounds up and ``nearest`` rounds 
        down. If `method` param is set to ``pd`` which refers to pd.interpolate 
        method , `kind` can be set to ``polynomial`` or ``pad`` interpolation. 
        Note that the polynomial requires you to specify an `order` while 
        ``pad`` requires to specify the `limit`. Default is ``slinear``.
        
    method: str, optional  
        Method of interpolation. Can be ``base`` for `scipy.interpolate.interp1d`
        ``mean`` or ``bff`` for scaling methods and ``pd``for pandas interpolation 
        methods. Note that the first method is fast and efficient when the number 
        of NaN in the array if relatively few. It is less accurate to use the 
        `base` interpolation when the data is composed of many missing values.
        Alternatively, the scaled method(the  second one) is proposed to be the 
        alternative way more efficient. Indeed, when ``mean`` argument is set, 
        function replaces the NaN values by the nonzeros in the raw array and 
        then uses the mean to fit the data. The result of fitting creates a smooth 
        curve where the index of each NaN in the raw array is replaced by its 
        corresponding values in the fit results. The same approach is used for
        ``bff`` method. Conversely, rather than averaging the nonzeros values, 
        it uses the backward and forward strategy  to fill the NaN before scaling.
        ``mean`` and ``bff`` are more efficient when the data are composed of 
        lot of missing values. When the interpolation `method` is set to `pd`, 
        function uses the pandas interpolation but ended the interpolation with 
        forward/backward NaN filling since the interpolation with pandas does
        not deal with all NaN at the begining or at the end of the array. Default 
        is ``base``.
        
    fill_value: array-like or (array-like, array_like) or ``extrapolate``, optional
        If a ndarray (or float), this value will be used to fill in for requested
        points outside of the data range. If not provided, then the default is
        NaN. The array-like must broadcast properly to the dimensions of the 
        non-interpolation axes.
        If a two-element tuple, then the first element is used as a fill value
        for x_new < x[0] and the second element is used for x_new > x[-1]. 
        Anything that is not a 2-element tuple (e.g., list or ndarray,
        regardless of shape) is taken to be a single array-like argument meant 
        to be used for both bounds as below, above = fill_value, fill_value.
        Using a two-element tuple or ndarray requires bounds_error=False.
        Default is ``extrapolate``. 
        
    kws: dict 
        Additional keyword arguments from :class:`spi.interp1d`. 
    
    Returns 
    -------
    array like - New interpoolated array. `NaN` values are interpolated. 
    
    Notes 
    ----- 
    When interpolated thoughout the complete frequencies  i.e all the frequency 
    values using the ``base`` method, the missing data in `arr`  can be out of 
    the `arr` range. So, for consistency and keep all values into the range of 
    frequency, the better idea is to set the param `fill_value` in kws argument
    of ``spi.interp1d`` to `extrapolate`. This will avoid an error to raise when 
    the value to  interpolated is extra-bound of `arr`. 
    
    
    References 
    ----------
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.interp1d.html
    https://www.askpython.com/python/examples/interpolation-to-fill-missing-entries
    
    Examples 
    --------
    >>> import numpy as np 
    >>> import matplotlib.pyplot as plt 
    >>> from pycsamt.processing.ctools import interpolate1d,
    >>> z = np.random.randn(17) *10 # assume 17 freq for 17 values of tensor Z 
    >>> z [[7, 10, 16]] =np.nan # replace some indexes by NaN values 
    >>> zit = interpolate1d (z, kind ='linear')
    >>> z 
    ... array([ -1.97732415, -16.5883156 ,   8.44484348,   0.24032979,
              8.30863276,   4.76437029, -15.45780568,          nan,
             -4.11301794, -10.94003412,          nan,   9.22228383,
            -15.40298253,  -7.24575491,  -7.15149205, -20.9592011 ,
                     nan]),
    >>> zn 
    ...array([ -1.97732415, -16.5883156 ,   8.44484348,   0.24032979,
             8.30863276,   4.76437029, -15.45780568,  -4.11301794,
           -10.94003412,   9.22228383, -15.40298253,  -7.24575491,
            -7.15149205, -20.9592011 , -34.76691014, -48.57461918,
           -62.38232823])
    >>> zmean = interpolate1d (z,  method ='mean')
    >>> zbff = interpolate1d (z, method ='bff')
    >>> zpd = interpolate1d (z,  method ='pd')
    >>> plt.plot( np.arange (len(z)),  zit, 'v--', 
              np.arange (len(z)), zmean, 'ok-',
              np.arange (len(z)), zbff, '^g:',
              np.arange (len(z)), zpd,'<b:', 
              np.arange (len(z)), z,'o', 
              )
    >>> plt.legend(['interp1d', 'mean strategy', 'bff strategy',
                    'pandas strategy', 'data'], loc='best')
    
    """
    method =str(method).strip().lower() 
    if method in ('pandas', 'pd', 'series', 'dataframe','df'): 
        method = 'pd' 
    elif method in ('interp1d', 'scipy', 'base', 'simpler', 'i1d'): 
        method ='base' 
    
    # check whether there is nan and masked invalid 
    # and take only the valid values 
    t_arr = arr.copy() 
    
    if method =='base':
        mask = ~np.ma.masked_invalid(arr).mask  
        arr = arr[mask] # keep the valid values
        f = spi.interp1d( x= np.arange(len(arr)), y= arr, kind =kind, 
                         fill_value =fill_value, **kws) 
        arr_new = f(np.arange(len(t_arr)))
        
    if method in ('mean', 'bff'): 
        arr_new = arr.copy()
        
        if method =='mean': 
            # use the mean of the valid value
            # and fill the nan value
            mean = t_arr[~np.isnan(t_arr)].mean()  
            t_arr[np.isnan(t_arr)]= mean  
            
        if method =='bff':
            # fill NaN values back and forward.
            t_arr = fillNaN(t_arr, method = method)
            t_arr= reshape(t_arr)
            
        yc, *_= scale_values(t_arr)
        # replace the at NaN positions value in  t_arr 
        # with their corresponding scaled values 
        arr_new [np.isnan(arr_new)]= yc[np.isnan(arr_new)]
        
    if method =='pd': 
        t_arr= pd.Series (t_arr)
        t_arr = np.array(t_arr.interpolate(
            method =kind, order=order, limit = limit ))
        arr_new = reshape(fillNaN(t_arr, method= 'bff')) # for consistency 
        
    return arr_new 


def make_griddata (Zcol, method ='cubic' , compute_err =False ):
    """ Get the tensor collection from a specific component (xx, xy, yy, or yx) 
    and create the griddata of that component. 
    
    Function can be use for plot purpose to compare the interpolate grid value 
    to the old one with missing data where the frequency values are lost. 
    
    
    Parameters 
    ----------
    Zcol: list of :class:`~.edi.Edi.Z objects  
        A collection of tensor objects. Tensor can be a real or an imaginary  
        part. 
        
    method: str 
        kind of interpolation. Can 'linear', 'nearest', 'cubic', 'auto'. If 
        ``auto`` is set. All the method shoud be computed and find the error 
        difference. Then should select the best minimal mean error. 
        
    compute_err: bool, 
        Compute the error of the interpolated method. Not usefull to set to 
        ``True`` when method is not set  to ``auto``. 
        
    Returns 
    -------- 
        -ndarray(nfrequency, number of tensor object). Number of tensor should 
            correspond to number of ediObject -> number of station. 
        - error : error of the interpolated method 
        - method: The method with the best mean error. Useful when the method 
            is set to ``auto``. 
        
    References 
    -----------
    https://stackoverflow.com/questions/37662180/interpolate-missing-values-2d-python
    answers from @M.T and edited by  
    https://stackoverflow.com/questions/37662180/interpolate-missing-values-2d-python/39596856#39596856
    
    """
    methods = ('linear', 'nearest', 'cubic', 'auto')
    method = str (method).lower() 
    if method not in methods : 
        raise ValueError (f"Unacceptable methods. Can only be "
                          f"{smart_format(methods, 'or')}.")
    if method =='auto': 
        method = ('linear', 'nearest', 'cubic') 
    else: method = [method]
    
    # for consistency reshape array of Zcollections 
    Zcol = list(map (lambda arr: reshape(arr)[:, None], Zcol ))
    # stacked the ndarray (nstatons, nfrequency )
    stackedZ = np.hstack(Zcol ) 

    x0 = np.arange (0, stackedZ.shape [1]) # number of stations 
    y0 = np.arange(0, stackedZ.shape [0] ) # number of frequency 

    xx, yy = np.meshgrid(x0, y0 ) 
    Z_masked = np.ma.masked_invalid(stackedZ) 

    x1 = xx[~Z_masked.mask] 
    y1 = yy[~Z_masked.mask] 
    
    Znew_stacked = Z_masked [~Z_masked.mask] 
    err= None
    if method !='auto': 
        z= spi.griddata ((x1, y1), Znew_stacked.ravel(),(xx, yy), 
                             method =method [0])  
        method = method [0] 
        
    else :
        Zgrid=[spi.griddata ((x1, y1), Znew_stacked.ravel(),(xx, yy), 
                             method =m ) for m in method ] 
        if method =='auto': compute_err =True  
        if compute_err :
            err =np.full_like (len(method), np.nan)
            nan_ix = np.isnan(stackedZ) 
            for i,  z in enumerate (Zgrid)  : 
                # compute the interpolate errors 
                # get the value replaced from index
                # got the interpolated values 
                mean = np.mean(z[nan_ix])
                err[i] = mean  
            min_ix   =  np.argmin(err)
            method = method[ int(min_ix)] 
            z = Zgrid[int (min_ix)]
        else : z = Zgrid [0]
        
    return z, err, method 

def qplot ( *args, leglabels =None, **kwargs): 
     """ Quick plot 
     :param args: list 
         List of  [x, y, c] plot 
    :param leglabels: list 
        List of list of legend labels. The number of labels must have
        the same length to the number of plot 
    :param kwargs: dict 
        Other additionnel `matplotlib.pyplot` plot parameters 
        
    :Example: 
        
       >>> qplot (np.arange(len(z0)), z0, 'o',
                  np.arange(len(z0)), zc, '--',
                  np.arange(len(z0)), ynew, 'ok-')
       >>> qplot(np.arange(len(z0)), znewc , 'g+',
                 np.arange(len(z0)), znew_corr, 'b:')
     """
     import matplotlib.pyplot as plt 
     
     plt.plot(*args, **kwargs)
     plt.xlabel('Log10 Frequency (Hz)')
     plt.ylabel('Impedance tensor Z (ohm.m)')
     plt.legend (leglabels , loc ='best')

def savitzky_golay2d ( z, window_size, order, derivative=None, mode = 'same'):
    """
    Two dimensional data smoothing and least-square gradient estimate. 
    
    Function computes the best fitting interpolating polynomial surface  and 
    its gradient. 
    
    Parameters
    ----------
    z : 2D Ndarray, shape (N,M)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
        
    derivative: str 
        the order of derivative of the 2D data. The maximum order of the
        derivative that can be computed obviously depends on the order of 
        the polynomial used in the fitting. It can be ``row``or ``column``, 
        indicating the direction of the derivative, or ``both``, which returns
        the gradient. Default is ``None``.
    mode: str 
         mode of the border prepending. Should be ``valid`` or ``same``. 
         ``same`` is used for prepending or appending the first value of
         array for smoothing.Default is ``same``.   
    Returns
    -------
    Tuple : ndarray, 
        the Gradient of a two-dimensional function and  the smoothed signal
        (or it's n-th derivative).
                             
    
    Savitsky-Golay filters can also be used to smooth two dimensional data 
    affected by noise. The algorithm is exactly the same as for the one 
    dimensional case, only the math is a bit more tricky. The basic algorithm 
    is as follow::
        
        * [1] for each point of the two dimensional matrix extract a sub-matrix,
            centered at that point and with a size equal to an odd number
            "window_size".
        * [2] for this sub-matrix compute a least-square fit of a polynomial 
            surface, defined as::
            .. math::
                
                p(x,y) = a_0 + a_1*x + a_2*y + a_3*x_2 + a_4*y_2 + a_5*x*y + ... . 
            
            Note that x and y are equal to zero at the central point.
        * [3] replace the initial central point with the value computed with the fit.
    
    Note that because the fit coefficients are linear with respect to the data
    spacing, they can pre-computed for efficiency. Moreover, it is important
    to appropriately pad the borders of the data, with a mirror image of the 
    data itself, so that the evaluation of the fit at the borders of the data
    can happen smoothly.
    
    Examples 
    --------
    >>> import matplotlib.pyplot as plt 
    >>> from pycsamt.processing import  savitzky_golay2d
    >>> # create some sample twoD data
    >>> x = np.linspace(-3,3,100)
    >>> y = np.linspace(-3,3,100)
    >>> X, Y = np.meshgrid(x,y)
    >>> Z = np.exp( -(X**2+Y**2))
    >>> # add noise
    >>> Zn = Z + np.random.normal( 0, 0.2, Z.shape )
    >>> # filter it
    >>> Zf = savitzky_golay2d( Zn, window_size=29, order=4)
    >>> # do some plotting
    >>> plt.matshow(Z)
    >>> plt.matshow(Zn)
    >>> plt.matshow(Zf)
    
    References 
    ----------
    .. [1] https://scipy.github.io/old-wiki/pages/Cookbook/SavitzkyGolay
    
    """
    # number of terms in the polynomial expression
    n_terms = ( order + 1 ) * ( order + 2)  / 2.0
    
    if  window_size % 2 == 0:
        raise ValueError('window_size must be odd')
    
    if window_size**2 < n_terms:
        raise ValueError('order is too high for the window size')

    half_size = window_size // 2
    
    # exponents of the polynomial. 
    # p(x,y) = a0 + a1*x + a2*y + a3*x^2 + a4*y^2 + a5*x*y + ... 
    # this line gives a list of two item tuple. Each tuple contains 
    # the exponents of the k-th term. First element of tuple is for x
    # second element for y.
    # Ex. exps = [(0,0), (1,0), (0,1), (2,0), (1,1), (0,2), ...]
    exps = [ (k-n, n) for k in range(order+1) for n in range(k+1) ]
    
    # coordinates of points
    ind = np.arange(-half_size, half_size+1, dtype=np.float64)
    dx = np.repeat( ind, window_size )
    dy = np.tile( ind, [window_size, 1]).reshape(window_size**2, )

    # build matrix of system of equation
    A = np.empty( (window_size**2, len(exps)) )
    for i, exp in enumerate( exps ):
        A[:,i] = (dx**exp[0]) * (dy**exp[1])
        
    # pad input array with appropriate values at the four borders
    new_shape = z.shape[0] + 2*half_size, z.shape[1] + 2*half_size
    Z = np.zeros( (new_shape) )
    # top band
    band = z[0, :]
    Z[:half_size, half_size:-half_size] =  band -  np.abs( 
        np.flipud( z[1:half_size+1, :] ) - band )
    # bottom band
    band = z[-1, :]
    Z[-half_size:, half_size:-half_size] = band  + np.abs( 
        np.flipud( z[-half_size-1:-1, :] )  -band ) 
    # left band
    band = np.tile( z[:,0].reshape(-1,1), [1,half_size])
    Z[half_size:-half_size, :half_size] = band - np.abs( 
        np.fliplr( z[:, 1:half_size+1] ) - band )
    # right band
    band = np.tile( z[:,-1].reshape(-1,1), [1,half_size] )
    Z[half_size:-half_size, -half_size:] =  band + np.abs(
        np.fliplr( z[:, -half_size-1:-1] ) - band )
    # central band
    Z[half_size:-half_size, half_size:-half_size] = z
    
    # top left corner
    band = z[0,0]
    Z[:half_size,:half_size] = band - np.abs( np.flipud(
        np.fliplr(z[1:half_size+1,1:half_size+1]) ) - band )
    # bottom right corner
    band = z[-1,-1]
    Z[-half_size:,-half_size:] = band + np.abs( np.flipud(
        np.fliplr(z[-half_size-1:-1,-half_size-1:-1]) ) - band ) 
    
    # top right corner
    band = Z[half_size,-half_size:]
    Z[:half_size,-half_size:] = band - np.abs( np.flipud(
        Z[half_size+1:2*half_size+1,-half_size:]) - band ) 
    # bottom left corner
    band = Z[-half_size:,half_size].reshape(-1,1)
    Z[-half_size:,:half_size] = band - np.abs( np.fliplr(
        Z[-half_size:, half_size+1:2*half_size+1]) - band ) 
    
    # solve system and convolve
    if derivative == None:
        m = np.linalg.pinv(A)[0].reshape((window_size, -1))
        return fftconvolve(Z, m, mode=mode)
    
    elif derivative == 'col':
        c = np.linalg.pinv(A)[1].reshape((window_size, -1))
        return fftconvolve(Z, -c, mode=mode)        
    elif derivative == 'row':
        r = np.linalg.pinv(A)[2].reshape((window_size, -1))
        return fftconvolve(Z, -r, mode=mode)        
    elif derivative == 'both':
        c = np.linalg.pinv(A)[1].reshape((window_size, -1))
        r = np.linalg.pinv(A)[2].reshape((window_size, -1))
        return (fftconvolve(Z, -r, mode=mode),
                fftconvolve(Z, -c, mode=mode) )
    

def exportfilterededis (
        ediObjs=None,
        window_size =5,
        c=2,
        *,
        method= 'tma', 
        savepath = None, 
        **kws
) -> object: 
    
    """ Export the filtered EDIs. 
    
    Beforehand, the function interpolates and fixes the missing data in the 
    attenuation band or weak frequencies but it does not bufering the frequencies. 
    To do that, prior use the function :func:`~.restoreZinDeadBand`. 
    
    Parameters 
    ----------
    ediObjs: list  of  pycsamt.core.edi.Edi or mtpy.core.edi.Edi objects 
        Collections of EDI-objects from `pyCSAMT`_ and `MTpy`_ packages
        
    window_size : int
        the length of the window. It is known as the filter width. Must be 
        greater than 1 and preferably be an odd integer number. Default is set 
        to ``5`` dipole lengths. 

    method: str, default ``tma``
        Filtering methods. Can be ``tma``, ``flma`` or ``ama`` for 
        'trimimg-moving-average, 'fixed-dipole-length moving average', and 
        'adaptative moving-average'. Refer to  the documentation of 
        :doc:`~.tma`, :doc:`~.ama` or :doc:`~.flma`.

    c : int, 
        A window-width expansion factor that must be input to the filter 
        adaptation process to control the roll-off characteristics
        of the applied Hanning window. It is recommended to select `c` between 
        ``1``  and ``4``.  Default is ``2``. 
        
    savepath: str 
        Full path to store the Edis outputted.
    
    kws: dict 
        Additional keywords for write new edifiles. Refer to
        :doc:`pycsamt.core.edi.Edi.write_new_edifile`
        
    Returns 
    -------
    ediObjs :  :attr:`pycsamt.core.edi.Edi_collection.ediObjs`
       Collection of  EDI-objects from `pyCSAMT`_ and `MTpy`_ packages 
    
    
    Examples 
    --------
    >>> from pycsamt.processimg import exportfilterededis 
    >>> edipath = 'data/edis' 
    >>> savepath= 'out/edi'
    >>> exportfilterededis (edipath, savepath = savepath3edis,
                            new_edi_fn ='f') # new ediname as prefix 
    
    """    

    ediObjs = get_ediObjs(ediObjs)
    
    if isinstance(ediObjs, str): 
        ediObjs = get_ediObjs(ediObjs)
    else :
        ediObjs = np.array( list(map( lambda o: _assert_edi_obj(o),  ediObjs )),
                            dtype =object) 
        
    method = str(method).lower().strip() 
    if method not in ('tma', 'ama', 'flma'): 
        raise ValueError (f"Expect 'tma', 'ama' or 'flma' filter not {method!r}")
        
    ffunc = tma 
    if method =='ama': 
        ffunc = ama  
    elif method =='flma': 
        ffunc = flma 
    
    # get the frequency 
    freqs = get_full_frequency(ediObjs)
    #-----XX------
    rhoxx = interpolate2d (make2d (ediObjs, 'resxx')) 
    phsxx = interpolate2d (make2d (ediObjs, 'phasexx'))
    rhoxx_err = interpolate2d (make2d (ediObjs, 'resxx_err')) 
    phsxx_err = interpolate2d (make2d (ediObjs, 'phasexx_err'))
    # -----XY-----
    rhoxy = interpolate2d (make2d (ediObjs, 'resxy')) 
    phsxy = interpolate2d (make2d (ediObjs, 'phasexy'))
    rhoxy_err = interpolate2d (make2d (ediObjs, 'resxy_err')) 
    phsxy_err = interpolate2d (make2d (ediObjs, 'phasexy_err'))
    # -----YX-----
    rhoyx = interpolate2d (make2d (ediObjs, 'resyx')) 
    phsyx = interpolate2d (make2d (ediObjs, 'phaseyx'))
    rhoyx_err = interpolate2d (make2d (ediObjs, 'resyx_err')) 
    phsyx_err = interpolate2d (make2d (ediObjs, 'phaseyx_err'))
    # -----YY-----
    rhoyy = interpolate2d (make2d (ediObjs, 'resyy'))
    phsyy = interpolate2d (make2d (ediObjs, 'phaseyy'))
    rhoyy_err = interpolate2d (make2d (ediObjs, 'resyy_err')) 
    phsyy_err = interpolate2d (make2d (ediObjs, 'phaseyy_err'))
        
    #compute all corrected data 
    zxx =ffunc(res2d =rhoxx, phs2d =phsxx, freqs= freqs ,
               window_size=window_size , out='z', c = c )
    zxy =ffunc(res2d =rhoxy , phs2d =phsxy , freqs= freqs ,
               window_size=window_size , out='z', c=c )
    zyx =ffunc(res2d =rhoyx, phs2d =phsyx, freqs= freqs ,
               window_size=window_size , out='z', c= c)
    zyy =ffunc(res2d =rhoyy, phs2d =phsyy, freqs= freqs,
               window_size=window_size , out='z', c= c )

    # -> compute error data
    zxx_err =ffunc(res2d =rhoxx_err, phs2d =phsxx_err, freqs= freqs ,
               window_size=window_size , out='z', c = c )
    zxy_err =ffunc(res2d =rhoxy_err , phs2d =phsxy_err , freqs= freqs ,
               window_size=window_size , out='z', c=c )
    zyx_err =ffunc(res2d =rhoyx_err, phs2d =phsyx_err, freqs= freqs ,
               window_size=window_size , out='z', c= c)
    zyy_err =ffunc(res2d =rhoyy_err, phs2d =phsyy_err, freqs= freqs,
               window_size=window_size , out='z', c= c )
    
    for kk in range  (len(ediObjs)):
        # make a Z object for each Edi
        Z = MTz.Z(z_array=np.zeros((len(freqs), 2, 2),dtype='complex'),
              z_err_array=np.zeros((len(freqs), 2, 2)),
              freq=freqs)

        Z.z[:, 0,  0] = reshape (zxx[:, kk], 1) 
        Z.z[:, 0,  1] = reshape (zxy[:, kk], 1)
        Z.z[:, 1,  0] = reshape (zyx[:, kk], 1) 
        Z.z[:, 0,  1] = reshape (zyy[:, kk], 1)
        
        # compute the absolute error to avoid numpy discarding 
        # the imaginary part of the complex number
        zxx_err = np.abs(zxx_err); zxy_err= np.abs(zxy_err)
        zyx_err = np.abs(zyx_err); zyy_err = np.abs(zyy_err)
        
        Z.z_err[:, 0,  0] = reshape (zxx_err[:, kk], 1) 
        Z.z_err[:, 0,  1] = reshape (zxy_err[:, kk], 1)
        Z.z_err[:, 1,  0] = reshape (zyx_err[:, kk], 1) 
        Z.z_err[:, 0,  1] = reshape (zyy_err[:, kk], 1)
        
        Z.compute_resistivity_phase()
        
        export2newedis(ediObj=ediObjs[kk] , new_Z=Z, 
                      savepath = savepath,  **kws)
        
    return ediObjs     
    
def _assert_emap_filter_args (*args): 
    """ Asserts argument of |EMAP| filter and returns useful arguments.
    
    :param args: Argument of EMAP filter. Refer to functions :func:`~.tma`, 
        :func:`~.flma` and :func:`~.ama` documentation. 
    """
    ediObjs, res2d , phs2d , freqs, c, window_size, component, out = args 
    
    component= str(component).lower().strip() 
    out= str(out).lower().strip() 
    try : 
        c = int (c) 
    except : TypeError(f'Expect an integer value not {type(c).__name__!r}')
    
    if out.find ('factor') >= 0 or out =='sf': 
        out ='sf'
    elif out in ('z', 'impedance', 'tensor'): out ='z'
    
    if component not in ('xx', 'xy', 'yx', 'yy'): 
        raise ValueError(f"Unacceptable component {component!r}. Expect "
                         "'xx', 'xy', 'yx' or 'yy'")
    if ediObjs is not None: 
        ediObjs = np.array( list(map( lambda o: _assert_edi_obj(o),ediObjs)),
                           dtype =object) 
        
        res2d= make2d(ediObjs, out=f'res{component}')
        phs2d = make2d(ediObjs, out=f'phase{component}')
        freqs = get_full_frequency(ediObjs) 
        
    ck= res2d is None  or phs2d is None  or freqs is None  
    if ediObjs is None and ck : 
        ck = f"{'Resistivity' if res2d is None else 'Phase'if phs2d is None else 'Frequency'}"
        raise ValueError(
            f"{ck!r} should not be None."
            )

    if len(res2d) != len(freqs): 
        raise ValueError ("Resistivity and frequency arrays must have a same"
                      f" length. But {len(res2d)} & {len(freqs)} were given")
    if len(res2d) != len(phs2d): 
        raise ValueError ("Resistivity and phase must have the same length."
                          f" But {len(res2d)} & {len(phs2d)} were given.")
    try : 
        window_size = int(window_size)
    except ValueError : 
        raise ValueError ('Could not convert {type(window_size).__name__!r} '
                          ' to integer: {window_size!r}')
 
    res2d = np.array (res2d)
    if window_size > res2d.shape [1]:
        raise ValueError ("window_size might not be greater than"
                          f" {str(res2d.shape [1])!r}")
    
    return res2d , phs2d , freqs, c, window_size, component, out   

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    