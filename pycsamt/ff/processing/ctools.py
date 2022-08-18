# -*- coding: utf-8 -*-
#       Created on Sat Dec 12 13:55:47 2020
#       Author: Kouadio K.Laurent<etanoyau@gmail.com>
#       Licence: LGPL

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
"""
import os 
import warnings 

import numpy as np
import pandas as pd  
import mtpy 
import pycsamt

from pycsamt.ff.core.edi import Edi 
import pycsamt.ff.core.z as MTz
# import  pycsamt.utils.func_utils  as func
from pycsamt.utils.decorator import donothing  
from pycsamt.utils.func_utils import (
    _assert_all_types, 
    subprocess_module_installation, 
    reshape_array, 
    scale_values, 
    ismissing, 
    spi,
    smart_format, 
    fillNaN
    )
try : 
    from pycsamt.__init__ import itqdm 
    if itqdm : 
        import tqdm
except: itqdm = False 

try : 
    from pycsamt.__init__ import inumba  
    if inumba : 
        import numba as nb
        
except: inumba = False 

try: 
    from sklearn.decomposition import PCA 
except : 
    is_success = subprocess_module_installation('sklearn')
    if not is_success : 
        raise ImportError( 'Could not import module `sklearn`. Please '
                          'install scikit-learn manually.')
        

def skew(ediObj, method='swift'): 
    r"""
    The conventional asymmetry parameter based on the Z magnitude. 
    
    Parameters 
    ---------
    ediObj: pycsamt.ff.core.edi.Edi or mtpy.core.edi.Edi 
        EDi object with full impedance tensor Z. 
    
    method: str 
        Kind of correction. Can be ``swift`` for the remove distorsion proposed 
        by Swift in 1967. The value close to 0. assume the 1D and 2D structures 
        and 3D otherwise. Conversly to ``bahr`` for the remove distorsion proposed  
        by Bahr in 1991. The latter threshold is set to 0.3. Above this value 
        the structures is 3D. 
    
    Returns 
    ------- 
        - array-like  (nfreq, ) Array of skew at each frequency 
        - Z object from :class:`~core.edi.Edi.Z.z` 
        
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
    
    if method =='swift': 
        skw = np.abs ( S1  / D2 )
    else : 
        D1S2 = (S2 * np.conj(D1)).imag 
        S1D2 = (D2 * np.conj(S1)).imag 
        skw = np.sqrt(np.abs( D1S2 - S1D2))/np.abs(D2)
        
    return skw 

        
def remove_static_shift( ediObj, ss_x=1., ss_y=1., **kws):
    r"""
    Remove static shift from the apparent resistivity
    
    Assume the original observed tensor Z is built by a static shift S
    and an unperturbated "correct" Z0 such as :math:`Z = S * Z0` therefore the 
    correct Z will be :math:` Z0 = S^(-1) * Z`.
        
    Parameters 
    ---------- 
    ediObj: pycsamt.ff.core.edi.Edi or mtpy.core.edi.Edi 
        EDi object with full impedance tensor Z. 
        
    ss_x: float 
        correction factor for x component


    ss_y: float 
        correction factor for y component

    Returns 
    ------- 
    pycsamt.ff.core.z.Z object -  new Z object with static shift removed


    Notes
    ------
    The factors are in resistivity scale, so the entries of  the matrix "S"
    need to be given by their square-roots!

    Example
    -------

    >>> import pycsamt.ff.processing as Processing 
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

        >>> import pycsamt.ff.processing as Processing 
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
    :meth:`pycsamt.ff.processing.Processing.correct_edi`. The Remove 
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
        
    :rtype: pycsamt.ff.core.z.Z
    
    .. note:: The factors are in resistivity scale, so the
              entries of  the matrix "S" need to be given by their
              square-roots!
              
    :Example:
        >>> from pycsamt.ff.processing import Processing 
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
    :rtype: pycsamt.ff.core.z.Z
    
    :Example: 
    >>> from pycsamt.ff.processing import Processing 
    >>> from pycsamt.ff.core.edi import Edi 
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
    >>> from pycsamt.ff.processing import Processing 
    >>> from pycsamt.ff.core.edi import Edi 
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
                z_real = reshape_array(
                    ediObj.Z.z[nz_index, ii, jj].real)
                z_imag = reshape_array(
                    ediObj.Z.z[nz_index, ii, jj].imag)
                z_err = reshape_array(
                    ediObj.Z.z_err[nz_index, ii, jj]) 
                # correct the values 
                z_real_c, *_ =scale_values (z_real,**kws) 
                z_imag_c, *_ = scale_values(z_imag, **kws) 
                z_err_c, *_ = scale_values(z_err, **kws) 
                
           # set the new Z object 
            new_Z.z[nz_index, ii, jj] = reshape_array(
                z_real_c, 1)   + 1j * reshape_array(z_imag_c, 1) 
            new_Z.z_err[nz_index, ii, jj] = reshape_array(z_err_c, 1)
            
    # compute resistivity and phase for new Z object
    new_Z.compute_resistivity_phase()
 
    return new_Z   
    
        
def _assert_edi_obj (obj)-> object: 
    """Assert that the given argument is an EDI -object from modules 
    EDi of pyCSAMT and MTpy packages. A TypeError will occurs otherwise.
    
    :param obj: Full path EDI file or `pyCSAMT`_ and `MTpy`_ object. 
    :type obj: str or str or  pycsamt.ff.core.edi.Edi or mtpy.core.edi.Edi 
    
    :return: Identical object after asserting.
    
    """
    if isinstance(obj, str): 
        obj = Edi(obj) 
    obj = _assert_all_types (obj, mtpy.core.edi.Edi,
                                  pycsamt.ff.core.edi.Edi)
    return  obj 

def export2newEdi (ediObj, new_Z , savepath =None, **kws):
    """ Export new EDI file from the former bbject with  a given new impedance 
    tensors. 
    
    The export is assumed a new output EDI resulting from multiples corrections 
    applications. 
    
    Parameters 
    -----------
    ediObj: str or  pycsamt.ff.core.edi.Edi or mtpy.core.edi.Edi 
        Full path to Edi file or object from `pyCSAMT`_ and `MTpy`_ packages 
    
    new_Z: ndarray (nfreq, 2, 2) 
        Ndarray of impendance tensors Z. The tensor Z is 3D array composed of 
        number of frequency `nfreq`and four components (``xx``, ``xy``, ``yx``,
        and ``yy``) in 2X2 matrices. The  tensor Z is a complex number. 
    
    Returns 
    --------
     ediObj from pycsamt.ff.core.edi.Edi or mtpy.core.edi.Edi 
     
    """
    
    ediObj = _assert_edi_obj(ediObj)
    ediObj.write_new_edifile( new_Z=new_Z, **kws)
    return ediObj 


@nb.njit if inumba else donothing("Skip numba when the latter doesn't work "
                                  "with the current numpy.__version__")  
def restoreZinDeadBand(ediObjs, *, buffer = None, kind='slinear',
                          method ='pd', **kws ): 
    """ Fix the weak of missing frequencies at the 'dead-band`-  and recover the 
    missing impedance tensors values. 
    
    The function uses the complete frequency (frequency with clean data) collected 
    thoughout the survey to recover by inter/extrapolating the missing or weak 
    frequencies thereby restoring the impedance tensors at that 'dead-band'. Note 
    that the 'dead- band' also known as 'attenuation -band' is where the AMT 
    signal is weak or generally abscent. 

    Parameters 
    ---------- 
    ediObjs: list  of  pycsamt.ff.core.edi.Edi or mtpy.core.edi.Edi objects 
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
        Additional keyword arguments from :func:`~interpolate1d`. 
    
    Returns 
    --------
        Array-like of pycsamt.ff.core.z.Z objects 
        Array collection of new Z impedances objects with dead-band tensor 
        recovered. :class:`pycsamt.ff.core.z.Z` are ndarray (nfreq, 2, 2). 
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
    >>> import numpy as np ; import matplotlib.pyplot as plt 
    >>> from pycsamt.ff.core.edi import Edi_collection 
    >>> from pycsamt.utils.func_utils import reshape_array, scale_values 
    >>> from pycsamt.ff.processing.ctools import restoreZinDeadBand
    >>> from pycsamt.ff.processing.ctools import fit_tensorfromrefreq 
    >>> from pycsamt.ff.processing.ctools import control_freq_buffer 
    >>> path2edi = 'data/3edis'
    >>> colObjs = Edi_collection (path2edi)
    >>> ediObjs = colObjs.ediObjs 
    >>> # One can specify the frequency buffer like the example below, However 
    >>> # it is not necessaray at leastr there is a purpose to fix  the frequencies
    >>> # at a spefic reason.
    >>> buffer = [1.45000e+04,1.11500e+01]
    >>> zobjs_b = restoreZinDeadBand(ediObjs, buffer = buffer) # with buffer 
    >>> zobjs = restoreZinDeadBand(ediObjs ) # with no buffer 
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
    >>> zfit = fit_tensorfromrefreq(refreq= reffreq, compfreq= z1freq, z=zxy)
    >>> # if one uses a buffer, we can check the interpolated frequency buffer 
    >>> # in the reference frequency 
    >>> control_freq_buffer (reffreq, buffer) 
    ...  array([1.470e+04, 1.125e+01])
    >>> # now we can set buffer to  [1.470e+04, 1.125e+01] and use the value 
    >>> # to find the index in raw zxy for plotting purpose to stand in the frequency range 
    >>> # mask the z value out of the buffer frequency range
    >>> ix_buf = reshape_array(np.argwhere (np.isin (reffreq, buffer ))) 
    >>> # buffered indexes stands for array([ 9, 50], dtype=int64)
    >>> # slice the buffered frequency and mask the useless frequency (out of buff range)
    >>> mask = ~np.isin (reffreq , reffreq[9 : 51] ) # exclude 51e index (ix_buf[1] +1)   
    >>> reffreq_buf = np.ma.masked_where (mask , reffreq) 
    >>> zxy_restored_buf = np.ma.masked_where (mask ,zxy_restored)   
    >>> zfit_buf = np.ma.masked_where (mask ,zfit)  
    >>> # not necessary, one can corrected z to get a smooth resistivity distribution 
    >>> zcorrected , *_ = scale_values (zxy_restored)                     
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

    """
    def z_transform (z , rfq, fq,  slice_= None): 
        """ Route to do the same task for real, imaginaral and error """
        with np.errstate(all='ignore'):
            z = reshape_array (z) 
            z = fit_tensorfromrefreq (compfreq= fq, refreq =rfq, z = z  ) 
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
                    z_real = reshape_array (ediObj.Z.z[nz_index, ii, jj].real) 
                    z_real = z_transform (z_real, rfq=cfreq, fq=zfreq, 
                                          slice_=slice_ 
                                          )
                    z_imag = reshape_array (ediObj.Z.z[nz_index, ii, jj].imag) 
                    z_imag = z_transform (z_imag, rfq=cfreq, fq=zfreq, 
                                          slice_=slice_ 
                                          )
                    z_err = reshape_array (ediObj.Z.z_err[nz_index, ii, jj]) 
                    z_err = z_transform (z_err, rfq=cfreq, fq=zfreq,
                                         slice_=slice_ 
                                         )
                # Use the new dimension of the z and slice z according 
                # the buffer range. make the new index start at 0. 
                new_nz_index = slice (
                    * np.array([ix_s, ix_end],dtype=np.int32)-ix_s)
       
                new_Z.z[new_nz_index, ii, jj] = reshape_array(
                    z_real, 1)   + 1j * reshape_array(z_imag, 1) 
                new_Z.z_err[new_nz_index, ii, jj] = reshape_array(z_err, 1)
            
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
    >>> from pycsamt.ff.processing.ctools import control_freq_buffer
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

            
def fit_tensorfromrefreq (refreq, compfreq , z , fill_value = np.nan): 
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
    >>> from pycsamt.ff.processing.ctools import fit_tensorfromrefreq
    >>> refreq = np.linspace(7e7, 1e0, 20) # 20 frequencies as reference
    >>> freq_ = np.hstack ((refreq.copy()[:7], refreq.copy()[12:] )) 
    >>> z = np.random.randn(len(freq_)) *10 # assume length of  freq as 
    >>>                 # the same like the tensor Z value 
    >>> zn  = fit_tensorfromrefreq (refreq, freq_, z)
    >>> z # some frequency values are missing by not visible. 
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
    
    freqn, mask = ismissing(arr =compfreq, refarr= refreq ,return_index='mask',
                            fill_value = fill_value )
    
    #mask_isin = np.isin(refreq, compfreq)
    z_new = np.full_like(freqn,fill_value = fill_value) 
    z_new[mask] = reshape_array(z) 
    
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
    >>> from pycsamt.ff.processing.ctools import interpolate1d,
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
            t_arr= reshape_array(t_arr)
            
        yc, *_= scale_values(t_arr)
        # replace the at NaN positions value in  t_arr 
        # with their corresponding scaled values 
        arr_new [np.isnan(arr_new)]= yc[np.isnan(arr_new)]
        
    if method =='pd': 
        t_arr= pd.Series (t_arr)
        t_arr = np.array(t_arr.interpolate(
            method =kind, order=order, limit = limit ))
        arr_new = reshape_array(fillNaN(t_arr, method= 'bff')) # for consistency 
        
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
    Zcol = list(map (lambda arr: reshape_array(arr)[:, None], Zcol ))
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



    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    