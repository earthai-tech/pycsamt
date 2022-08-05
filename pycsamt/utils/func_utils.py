# -*- coding: utf-8 -*-
#       Created on Sun Sep 13 09:24:00 2020
#       Author: Kouadio K.Laurent<etanoyau@gmail.com>
#       Licence: LGPL

"""
.. _module-Func-utils::`pycsamt.utils.func_utils`  
    :synopsis: helpers functions 
    
    * check_dimensionality
    * subprocess_module_installation
    * cpath
    * smart_format
    * make_introspection
    * show_quick_edi_stats
    * sPath
    * averageData 
    * concat_array_from_list 
    * sort_array_data 
    * transfer_array_ (deprecated)
    * interpol_scipy 
    * _set_depth_to_coeff 
    * broke_array_to_ 
    *  _OlDFUNCNOUSEsearch_fill_data (deprecated)
    * _search_ToFill_Data 
    * straighten_out_list  
    * take_firstValue_offDepth 
    * dump_comma 
    * build_wellData 
    * compute_azimuth 
    * build_geochemistry_sample 
    * _nonelist_checker 
    * _order_well 
    * intell_index 
    * _nonevalue_checker 
    * _clean_space 
    *_cross_eraser 
    * _remove_str_word 
    * stn_check_split_type
    * minimum_parser_to_write_edi
    * round_dipole_length
    * keepmin  
    * get_closest_value
    * geo_length_checker
    * fr_en_parser
    * convert_csvdata_from_fr_to_en
    * make_ids
    * make_ll_coordinates 
    * resize_resphase_values 
    * _assert_all_types 
    * scalePosition 
    *
"""

import os 
import sys 
import subprocess 
import shutil 
import warnings
import inspect
import csv
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt
from copy import deepcopy

import  pycsamt.utils.gis_tools as gis
import pycsamt.utils.exceptions as CSex
from pycsamt.utils.decorator import deprecated 
from pycsamt.utils._csamtpylog import csamtpylog

_logger = csamtpylog.get_csamtpy_logger(__name__)

_msg= ''.join([
    'Note: need scipy version 0.14.0 or higher or interpolation,',
    ' might not work.']
)
_msg0 = ''.join([
    'Could not find scipy.interpolate, cannot use method interpolate'
     'check installation you can get scipy from scipy.org.']
)

try:
    import scipy
    scipy_version = [int(ss) for ss in scipy.__version__.split('.')]
    if scipy_version[0] == 0:
        if scipy_version[1] < 14:
            warnings.warn(_msg, ImportWarning)
            _logger.warning(_msg)
            
    import scipy.interpolate as spi
    from scipy.optimize import curve_fit

    interp_import = True
 # pragma: no cover
except ImportError: 
    
    warnings.warn(_msg0)
    _logger.warning(_msg0)
    
    interp_import = False


    
    
      
def make_ll_coordinates(reflong, reflat, nsites,  *,  r=45.,
                        sep='1km', order= '+', todms=False): 
    """ Generate multiples stations coordinates (longitudes, latitudes)
    from a reference station/sites.
    
    One degree of latitude equals approximately 364,000 feet (69 miles), 
    one minute equals 6,068 feet (1.15 miles), and one-second equals 101 feet.
    One-degree of longitude equals 288,200 feet (54.6 miles), one minute equals
    4,800 feet (0.91 mile), and one second equals 80 feet. Illustration showing
    longitude convergence. (1 feet ~=0.3048 meter)
    
    Parameters 
    ----------
    reflong: float or string 
        Reference longitude  in degree decimal or in DD:MM:SS for the first 
        site considered as the origin of the lamdmark.
        
    reflat: float or string 
        Reference latitude in degree decimal or in DD:MM:SS for the reference  
        site considered as the landmark origin.
        
    nsites: int or float 
        Number of site to generate the coordinates onto. 
        
    r: float or int 
        The rotate angle in degrees. Rotate the angle feature the direction
        of the projection line. Default value is ``45`` degrees. 
        
    sep: float or str 
        Distance of seperation between different sites in meters. If the value 
        is given as string type, except the ``km``, it should be considered as
        a ``m`` value. Only meters and kilometers are accepables.
        
    order: str 
        Direction of the projection line. By default the projected line is 
        in ascending order i.e. from SW to NE with angle `r` set to ``45``
        degrees. Could be ``-`` for descending order. Any other value should 
        be in ascending order. 
    
    todms: bool 
        Convert the degree decimal values into the DD:MM:SS. Default is ``False``. 
        
        
    Returns 
    -------
        Tuple of  generated projected coordinates longitudes and latitudes
        either in degree decimals or DD:MM:SS
        
        
    Notes 
    ------
    The distances vary. A degree, minute, or second of latitude remains 
    fairly constant from the equator to the poles; however a degree, minute,
    or second of longitude can vary greatly as one approaches the poles
    and the meridians converge.
        
    References 
    ----------
    https://math.answers.com/Q/How_do_you_convert_degrees_to_meters
    
    Examples 
    --------
    >>> from pycsamt.utils.func_utils import make_ll_coordinates 
    >>> rlons, rlats = make_ll_coordinates('110:29:09.00', '26:03:05.00', 
    ...                                     nsites = 7, todms=True)
    >>> rlons
    ... array(['110:29:09.00', '110:29:35.77', '110:30:02.54', '110:30:29.30',
           '110:30:56.07', '110:31:22.84', '110:31:49.61'], dtype='<U12')
    >>> rlats 
    ... array(['26:03:05.00', '26:03:38.81', '26:04:12.62', '26:04:46.43',
           '26:05:20.23', '26:05:54.04', '26:06:27.85'], dtype='<U11')
    
    """ 
    def assert_ll(coord):
        """ Assert coordinate when the type of the value is string."""
        try: coord= float(coord)
        except ValueError: 
            if ':' not in coord: 
                raise ValueError(f'Could not convert value to float: {coord!r}')
            else : 
                coord = gis.convert_position_str2float(coord)
        return coord
    
    nsites = int(_assert_all_types(nsites,int, float)) 

    sep=str(sep).lower() 
    if sep.find('km')>=0: # convert to meter 
        sep = float(sep.replace('km', '')) *1e3 
    elif sep.find('m')>=0: sep = float(sep.replace('m', '')) 
    sep = float(sep) # for consistency 
    
    if order in ('descending', 'down', '-'): order = '-'
    else: order ='+'
    # compute length of line using the reflong and reflat
    # the origin of the landmark is x0, y0= reflong, reflat
    x0= assert_ll(reflong) ; y0= assert_ll(reflat) 
    xinf = x0  + (np.sin(np.deg2rad(r)) * sep * nsites) / (364e3 *.3048) 
    yinf = y0 + np.cos(np.deg2rad(r)) * sep * nsites /(2882e2 *.3048)
    
    reflon_ar = np.linspace(x0 , xinf, nsites ) 
    reflat_ar = np.linspace(y0, yinf, nsites)
    #--------------------------------------------------------------------------
    # r0 = np.sqrt(((x0-xinf)*364e3 *.3048)**2 + ((y0 -yinf)*2882e2 *.3048)**2)
    # print('recover distance = ', r0/nsites )
    #--------------------------------------------------------------------------
    if todms:
        reflat_ar = np.array(list(
            map(lambda l: gis.convert_position_float2str(float(l)), reflat_ar)))
        reflon_ar = np.array(list(
            map(lambda l: gis.convert_position_float2str(float(l)), reflon_ar)))
    
    return (reflon_ar , reflat_ar ) if order =='+' else (
        reflon_ar[::-1] , reflat_ar[::-1] )
    


def build_array_from_objattr(obj, attr): 
    """ Quick build array of object attributes value from collections object.
    
    :param obj: Iterable collections objects. 
    :type obj: list, tuple or np.ndarray 
    
    :param attr: attribute value to retrieve from each collection object  
    :type attr: str 
    
    :Example: 
        >>> from pycsamt.ff.core.edi import Edi_Collection
        >>> from pycsamt.utils.func_utils import 
        >>> edipath = r'/Users/Daniel/Desktop/ediout'
        >>> cObjs = Edi_collection (edipath)
        >>> # retrieve the latitude data and make array 
        >>> qar = array_from_obj(cobjb1.ediObjs, 'lat')
        ... array([0.        , 0.        , 0.00027778, 0.00027778, 0.00055556,
               0.00083333, 0.00083333])
    
    """
    return np.array(list(map(lambda r: getattr(r, attr), obj )))


def _assert_all_types (
        obj: object , 
        *expected_objtype: type 
 ) -> object: 
    """ Quick assertion of object type. Raise an `TypeError` if 
    wrong type is given."""
    # if np.issubdtype(a1.dtype, np.integer): 
    if not isinstance( obj, expected_objtype): 
        raise TypeError (
            f'Expected {smart_format(tuple (o.__name__ for o in expected_objtype))}'
            f' type{"s" if len(expected_objtype)>1 else ""} '
            f'but `{type(obj).__name__}` is given.')
            
    return obj 

def scalePosition(ydata , xdata= None, func = None ,c_order= 0,
        show: bool =False, todms=False,
        **kws): 
    """ Correct data location or position and return new corrected location
    or data. 
    
    Parameters 
    ----------
    ydata: array_like, series or dataframe
        The dependent data, a length M array - nominally ``f(xdata, ...)``.
        
    xdata: array_like or object
        The independent variable where the data is measured. Should usually 
        be an M-length sequence or an (k,M)-shaped array for functions with
        k predictors, but can actually be any object. If ``None``, `xdata` is 
        generated by default using the length of the given `ydata`.
        
    func: callable 
        The model function, ``f(x, ...)``. It must take the independent variable 
        as the first argument and the parameters to fit as separate remaining
        arguments. The default `func` is ``linear`` function i.e  for ``f(x)= ax +b``. 
        where `a` is slope and `b` is the intercept value. Setting your own 
        function for better fitting is recommended. 
        
    c_order: int or str
        The index or the column name if ``ydata`` is given as a dataframe to 
        select the right column for scaling.
    todms: bool 
        Convert the decimal long/lat to DD:MM:SS. Default is ``False``. 
        
    show: bool 
        Quick visualization of data distribution.

    kws: dict 
        Additional keyword argument from  `scipy.optimize_curvefit` parameters. 
        Refer to `scipy.optimize.curve_fit`_.  
        
    Returns 
    --------
    - ydata - array -like - Data scaled 
    - popt - array-like Optimal values for the parameters so that the sum of 
        the squared residuals of ``f(xdata, *popt) - ydata`` is minimized.
    - pcov - array like The estimated covariance of popt. The diagonals provide
        the variance of the parameter estimate. To compute one standard deviation 
        errors on the parameters use ``perr = np.sqrt(np.diag(pcov))``. How the
        sigma parameter affects the estimated covariance depends on absolute_sigma 
        argument, as described above. If the Jacobian matrix at the solution
        doesn’t have a full rank, then ‘lm’ method returns a matrix filled with
        np.inf, on the other hand 'trf' and 'dogbox' methods use Moore-Penrose
        pseudoinverse to compute the covariance matrix.
        
    Examples
    --------
    >>> from pycsamt.utils.func_utils  import scalePosition 
    >>> from pycsamt.ff.core.edi import Edi_Collection 
    >>> edipath = r'/Users/Daniel/Desktop/ediout'
    >>> cObjs = Edi_collection (edipath)
    >>> # correcting northing coordinates from latitude data 
    >>> # corrected latitude coordinates using the default x.
    >>> lat_corrected, *_= scalePosition(ydata =cObjs.lat[:12])
    >>> cObjs.lat[:12]
    ... array([0.        , 0.        , 0.00027778, 0.00027778, 0.00055556,
    ...       0.00083333, 0.00083333, 0.00111111, 0.00138889, 0.00138889,
    ...       0.00166667, 0.00194444])
    >>> lat_corrected
    ... array([-1.31766381e-04,  4.79150482e-05,  2.27596478e-04,  4.07277908e-04,
    ...        5.86959337e-04,  7.66640767e-04,  9.46322196e-04,  1.12600363e-03,
    ...        1.30568506e-03,  1.48536649e-03,  1.66504792e-03,  1.84472934e-03])
    >>> lat_corrected_dms, *_= scalePosition(ydata =cObjs.lat[:12], todms=True)
    ... array(['0:00:00.47', '0:00:00.17', '0:00:00.82', '0:00:01.47',
    ...       '0:00:02.11', '0:00:02.76', '0:00:03.41', '0:00:04.05',
    ...       '0:00:04.70', '0:00:05.35', '0:00:05.99', '0:00:06.64'],
          dtype='<U10')
    """
    
    def linfunc (x, a, b): 
        """ Set the simple linear function"""
        return a * x + b 
        
    if str(func).lower() in ('none' , 'linear'): 
        func = linfunc 
    elif not hasattr(func, '__call__') or not inspect.isfunction (func): 
        raise TypeError(
            f'`func` argument is a callable not {type(func).__name__!r}')
        
    ydata = _assert_all_types(ydata, list, tuple, np.ndarray,
                              pd.Series, pd.DataFrame  )
    c_order = _assert_all_types(c_order, int, float, str)
    try : c_order = int(c_order) 
    except: pass 

    if isinstance(ydata, pd.DataFrame): 
        if c_order ==0: 
            warnings.warn("The first column of the data should be considered"
                          " as the `y` target.")
        if c_order is None: 
            raise TypeError('Dataframe is given. The `c_order` argument should '
                            'be defined for column selection. Use column name'
                            ' instead')
        if isinstance(c_order, str): 
            # check whether the value is on the column name
            if c_order.lower() not in list(map( 
                    lambda x :x.lower(), ydata.columns)): 
                raise ValueError (
                    f'c_order {c_order!r} not found in {list(ydata.columns)}'
                    ' Use the index instead.')
                # if c_order exists find the index and get the 
                # right column name 
            ix_c = list(map( lambda x :x.lower(), ydata.columns)
                        ).index(c_order.lower())
            ydata = ydata.iloc [:, ix_c] # series 
        elif isinstance (c_order, (int, float)): 
            c_order =int(c_order) 
            if c_order >= len(ydata.columns): 
                raise ValueError(
                    f"`c_order`'{c_order}' should be less than the number of " 
                    f"given columns '{len(ydata.columns)}'. Use column name instead.")
            ydata= ydata.iloc[:, c_order]
            
    ydata = np.array(ydata)
    if xdata is None: 
        xdata = np.linspace(0, 4, len(ydata))
    if len(xdata) != len(ydata): 
        raise ValueError(" `x` and `y` arrays must have the same length."
                        "'{len(xdata)}' and '{len(ydata)}' are given.")
        
    popt, pcov = curve_fit(func, xdata, ydata, **kws)
    ydata_new = func(xdata, *popt)
    
    if show:
        plt.plot(xdata, ydata, 'b-', label='data')
        plt.plot(xdata, func(xdata, *popt), 'r-',
             label='fit: a=%5.3f, b=%5.3f' % tuple(popt))
        plt.xlabel('x')
        plt.ylabel('y')
        plt.legend()
        plt.show()
        
    if todms: 
        ydata_new= np.array(list(map(lambda l: gis.convert_position_float2str(
            float(l)), ydata_new)))
        
    return ydata_new, popt, pcov 


def make_ids(ediObjs, prefix =None, how ='py'): 
    """ Generate auto Id according to the number of given sites. 
    
    :param ediObjs: list of EDI object , composed of a collection of 
        pycsamt.ff.core.edi.Edi object
    :type ediObjs: pycsamt.ff.core.edi.Edi_Collection 
    
    :param prefix: string value to add as prefix of given id. Prefix can be 
        the site name.
    :type prefix: str 
    
    :param how: Mode to index the station. Default is 'Python indexing' i.e. 
        the counting starts by 0. Any other mode will start the counting by 1.
    :type cmode: str 
    
    :return: ID number formated 
    :rtype: list 
    
    :Example: 
        >>> from pycsamt.utils.func_utils import make_ids 
        >>> values = ['edi1', 'edi2', 'edi3'] 
        >>> make_ids (values, 'ix')
        ... ['ix0', 'ix1', 'ix2']
        
    """ 
    fm='{:0' +'{}'.format(int(np.log10(len(ediObjs))) + 1) +'}'
    id_ =[str(prefix) + fm.format(i if how=='py'else i+ 1 ) if prefix is not 
          None else fm.format(i if how=='py'else i+ 1) 
          for i in range(len(ediObjs))] 
    return id_

def fit_by_ll(ediObjs): 
    """ Fit edi by location reorganize edi according to the site longitude and 
    latitude. 
    
    Edis data are mostly reading an alphabetically order, so the reoganization  
    according to the location(longitude and latitude) is usefull for distance 
    betwen site computing with a right position at each site.  
    
    :param ediObjs: list of EDI object , composed of a collection of 
        pycsamt.ff.core.edi.Edi object 
    :type ediObjs: pycsamt.ff.core.edi.Edi_Collection 
    
    :returns: array splitted into ediObjs and Edifiles basenames 
    :rtyple: tuple 
    
    :Example: 
        >>> import numpy as np 
        >>> from pycsamt.ff.core.edi import Edi_Collection 
        >>> from pycsamt.utils.func_utils import fit_by_ll
        >>> edipath ='data/edi_ss' 
        >>> cediObjs = Edi_Collection (edipath) 
        >>> ediObjs = np.random.permutation(cediObjs.ediObjs) # shuffle the  
        ... # the collection of ediObjs 
        >>> ediObjs, ediObjbname = fit_by_ll(ediObjs) 
        ...
    
    """
    #get the ediObjs+ names in ndarray(len(ediObjs), 2) 
    objnames = np.c_[ediObjs, np.array(
        list(map(lambda obj: os.path.basename(obj.edifile), ediObjs)))]
    lataddlon = np.array (list(map(lambda obj: obj.lat + obj.lon , ediObjs)))
    sort_ix = np.argsort(lataddlon) 
    objnames = objnames[sort_ix ] 
    #ediObjs , objbnames = np.hsplit(objnames, 2) 
    return objnames[:, 0], objnames[:, -1]


def get_interpolate_freqs (ediObjs, to_log10 =False): 
    """ From EDI objects, collected thefrequencies Min, Max and find the 
    frequency for interpolating. 
    
    This method is usefull when a collection of EDI-files have some missing 
    frequency data. 
    
    :param ediObjs: list - Collections of EDI-objects 
    :rtype: pycsamt.ff.core.edi.Edi 
    
    :param to_log10: Put the interpolated min-max frequency to log10 values
    :type to_log10: bool 
    
    :returns: Array-like min max and number of frequency for interpolating. 
    :rtype: tuple
    
    :Example: 
        >>> from pycsamt.ff.core.edi import Edi_collection 
        >>> from pycsamt.utils.func_utils import find_interpolate_freq
        >>> edipath = r'/Users/Daniel/Desktop/ediout'
        >>> cObjs = Edi_collection (edipath)
        >>> ifreqs, nfreq= get_interpolate_freq(cObjs.ediObjs) 
        >>> ifreqs, nfreq
        ... array([5.625e+00, 5.880e+04]), 56 # i.e minimum frequency 
            for interpolation is 5.625e+00 while the max is 5.880e+04
            and the number of frequency that can be interpolated is 56. 
    """
    # create a matrix of ndarray (nedifiles, [nfreq, minfreq, maxfreq])
    if not isinstance(ediObjs, (list, tuple)): 
        ediObjs =[ediObjs]

    freqs= np.zeros((len(ediObjs), 3))
    #edins = np.zeros_like(freqs , dtype ='>U20')
    for ii, ediobj in enumerate (ediObjs) : 
        freqs[ii] = [len(ediobj.Z._freq),ediobj.Z._freq.min(), 
                     ediobj.Z._freq.max() ]
        #edins [ii]= os.path.basename(ediobj.edifile).replace('.edi', '') 
    # frequency for interpolate should be in the bound of frequency ranges 
    # take all minimum in axis 0 and replace the 
    # min ranges frequency by the higher to bound the frequency ranges 
    freqi = freqs.min(0) 
    freqi[1]= freqs[:, 1].max() 

    # find  the index a  check the edifiles 
    #ix = np.where (np.all (freqs ==freqs.max(0), axis = 1)) 
    # get the number of frequency for max freq and in freq columns 
    fmaxs = freqs[:, 0][freqs [:, -1]==freqi[-1]]
    fmins = freqs[:, 0][freqs [:, -2]==freqi[-2]]
    # for consistency if array are composed of repeated max |min values 
    nfreq = max ( fmaxs.max(),  fmins.max()) 
    ifreqs = freqi[1:] # exclude the number of freq of each edi (column 0)

    return (np.log10(ifreqs), np.int(nfreq)) if to_log10 else (
        ifreqs, np.int(nfreq))
                     

def check_dimensionality(obj, data, z, x):
    """ Check dimensionality of data and fix it.
    
    :param obj: Object, can be a class logged or else.
    :param data: 2D grid data of ndarray (z, x) dimensions
    :param z: array-like should be reduced along the row axis
    :param x: arraylike should be reduced along the columns axis.
    """
    def reduce_shape(Xshape, x, axis_name =None): 
        """ Reduce shape to keep the same shape"""
        mess ="`{0}` shape({1}) {2} than the data shape `{0}` = ({3})."
        ox = len(x) 
        dsh = Xshape 
        if len(x) > Xshape : 
            x = x[: int (Xshape)]
            obj._logging.debug(''.join([
                f"Resize {axis_name!r}={ox!r} to {Xshape!r}.", 
                mess.format(axis_name, len(x),'more',Xshape)])) 
                                    
        elif len(x) < Xshape: 
            Xshape = len(x)
            obj._logging.debug(''.join([
                f"Resize {axis_name!r}={dsh!r} to {Xshape!r}.",
                mess.format(axis_name, len(x),'less', Xshape)]))
        return int(Xshape), x 
    
    sz0, z = reduce_shape(data.shape[0], 
                          x=z, axis_name ='Z')
    sx0, x =reduce_shape (data.shape[1],
                          x=x, axis_name ='X')
    data = data [:sz0, :sx0]
    
    return data , z, x 

def subprocess_module_installation (module, upgrade =True , DEVNULL=False,
                                    action=True, verbose =0, **subpkws): 
    """ Install or uninstall a module using the subprocess under the hood.
    
    :param module: str, module name 
    :param upgrade:bool, install the lastest version.
    :param verbose:output a message 
    :param DEVNULL: decline the stdoutput the message in the console 
    :param action: str, install or uninstall a module 
    :param subpkws: additional subprocess keywords arguments.
    
    :Example: 
        >>> from pycsamt.utils.func_utils import subprocess_module_installation
        >>> subprocess_module_installation(
            'tqdm', action ='install', DEVNULL=True, verbose =1)
        >>> subprocess_module_installation(
            'tqdm', action ='uninstall', verbose =1)
    """
    #implement pip as subprocess 
    # refer to https://pythongeeks.org/subprocess-in-python/
    if not action: 
        if verbose > 0 :
            print("---> No action `install`or `uninstall`"
                  f" of the module {module!r} performed.")
        return action  # DO NOTHING 
    
    MOD_IMP=False 

    action_msg ='uninstallation' if action =='uninstall' else 'installation' 

    if action in ('install', 'uninstall', True) and verbose > 0:
        print(f'---> Module {module!r} {action_msg} will take a while,'
              ' please be patient...')
        
    cmdg =f'<pip install {module}> | <python -m pip install {module}>'\
        if action in (True, 'install') else ''.join([
            f'<pip uninstall {module} -y> or <pip3 uninstall {module} -y ',
            f'or <python -m pip uninstall {module} -y>.'])
        
    upgrade ='--upgrade' if upgrade else '' 
    
    if action == 'uninstall':
        upgrade= '-y' # Don't ask for confirmation of uninstall deletions.
    elif action in ('install', True):
        action = 'install'

    cmd = ['-m', 'pip', f'{action}', f'{module}', f'{upgrade}']

    try: 
        STDOUT = subprocess.DEVNULL if DEVNULL else None 
        STDERR= subprocess.STDOUT if DEVNULL else None 
    
        subprocess.check_call(
            [sys.executable] + cmd, stdout= STDOUT, stderr=STDERR,
                              **subpkws)
        if action in (True, 'install'):
            # freeze the dependancies
            reqs = subprocess.check_output(
                [sys.executable,'-m', 'pip','freeze'])
            [r.decode().split('==')[0] for r in reqs.split()]
            _logger.info( f"{action_msg.capitalize()} of `{module}` "
                         "and dependancies was successfully done!") 
        MOD_IMP=True
        
    except: 
        _logger.error(f"Failed to {action} the module =`{module}`.")
        
        if verbose > 0 : 
            print(f'---> Module {module!r} {action_msg} failed. Please use'
                f' the following command: {cmdg} to manually do it.')
    else : 
        if verbose > 0: 
            print(f"{action_msg.capitalize()} of `{module}` "
                      "and dependancies was successfully done!") 
        
    return MOD_IMP 
        
def smart_format(iter_obj): 
    """ Smart format iterable obj 
    :param iter_obj: iterable obj 
    
    :Example: 
        >>> from pycsamt.utils.func_utils import smart_format
        >>> smart_format(['model', 'iter', 'mesh', 'data'])
        ... 'model','iter','mesh' and 'data'
    """
    try: 
        iter(iter_obj) 
    except:  return f"{iter_obj}"
    
    iter_obj = [str(obj) for obj in iter_obj]
    if len(iter_obj) ==1: 
        str_litteral= ','.join([f"{i!r}" for i in iter_obj ])
    elif len(iter_obj)>1: 
        str_litteral = ','.join([f"{i!r}" for i in iter_obj[:-1]])
        str_litteral += f" and {iter_obj[-1]!r}"
    return str_litteral

def make_introspection(Obj , subObj): 
    """ Make introspection by using the attributes of instance created to 
    populate the new classes created.
    :param Obj: callable 
        New object to fully inherits of `subObject` attributes 
    :param subObj: Callable 
        Instance created.
    """
    # make introspection and set the all  attributes to self object.
    # if Obj attribute has the same name with subObj attribute, then 
    # Obj attributes get the priority.
    for key, value in  subObj.__dict__.items(): 
        if not hasattr(Obj, key) and key  != ''.join(['__', str(key), '__']):
            setattr(Obj, key, value)
            
def cpath (savepath=None , dpath= None): 
    """ Control the existing path and create one of it does not exist.
    
    :param savepath: Pathlike obj, str 
    :param dpath: str, default pathlike obj
    """
    if dpath is None:
        file, _= os.path.splitext(os.path.basename(__file__))
        dpath = ''.join(['_', file,
                         '_']) #.replace('.py', '')
    if savepath is None : 
        savepath  = os.path.join(os.getcwd(), dpath)
    if savepath is not None:
        try :
            if not os.path.isdir(savepath):
                os.mkdir(savepath)#  mode =0o666)
        except : pass
    
    return savepath   

def show_quick_edi_stats(nedic , nedir, fmtl='~', lenl=77): 
    """ Format the Edi files and ckeck the number of edifiles
    successfully read. 
    :param nedic: number of input or collected edifiles 
    :param nedir: number of edifiles read sucessfully 
    :param fmt: str to format the stats line 
    :param lenl: length of line denileation."""
    
    def get_obj_len (value):
        """ Control if obj is iterable then take its length """
        try : 
            iter(value)
        except :pass 
        else : value =len(value)
        return value 
    nedic = get_obj_len(nedic)
    nedir = get_obj_len(nedir)
    
    print(fmtl * lenl )
    mesg ='|'.join( ['|{0:<15}{1:^2} {2:<7}',
                     '{3:<15}{4:^2} {5:<7}',
                     '{6:<9}{7:^2} {8:<7}%|'])
    print(mesg.format('Data collected','=',  nedic, 'EDI success. read',
                      '=', nedir, 'Rate','=', round ((nedir/nedic) *100, 2),
                      2))
    print(fmtl * lenl )

def sPath (name_of_path:str):
    """ Savepath func. Create a path  with `name_of_path` 
    if path not exists.
    :param name_of_path: str, Path-like object. If path does not exist,
    `name_of_path` should be created.
    """
    try :
        savepath = os.path.join(os.getcwd(), name_of_path)
        if not os.path.isdir(savepath):
            os.mkdir(name_of_path)#  mode =0o666)
    except :
        warnings.warn("The path seems to be existed!")
        return
    return savepath 


def averageData(np_array, filter_order=0, 
                axis_average=0, astype="float32"): #array_of_average_array=0
    """
    Parameters  
    -----------
        * np_array  : numpy array 
                must be an array data 
        * filter_order  : int 
                must be the index of the column you want to sort 
            
        * axis average  : int  
                axis you want to see data averaged, also
                it is the concatenate axis  
                default is axis=0 
        
        * astype*: str , 
                is the ndarray dtype array .
                change to have an outup arry dtype , you want .
    
    Returns  
    --------
        numpy array 
            Data averaged array
    
    :Example : 
        
        >>> import numpy as np 
        >>> list8=[[4,2,0.1],[8,2,0.7],[10,1,0.18],[4,3,0.1],
        ...               [7,2,1.2],[10,3,0.5],[10,1,0.5],[8.2,0,1.9],
        ...               [10,7,0.5],[10,1,0.5],
        ...               [2,0,1.4],[5,4,0.5],[10,2,0.7],[7,2,1.078],
        ...               [10,2,3.5],[10,8,1.9]]
        >>> np_test=np.array(list8)
        >>> ss=averageData(np_array=np_test,filter_order=1,
        ...               axis_average=0, astype="int")
        >>> ss
    """
    idx,sep_counts=0,0
    global_list=[]
    #Filter the array 
    np_array=np_array[np_array[:,filter_order].argsort(kind="mergesort")]
    #returns the differents value on the filtersort index of array :
    values, counts =np.unique(np_array[:,filter_order], return_counts=True)
    # append array with numpy :(values.shape)
    #             # new_array=np.append(new_array,rowline,axis=0)
    for rowline in np_array :
        
        if rowline[filter_order] ==values.max():
            new_array=np_array[sep_counts:,:]
            temp_array=new_array.mean(axis=axis_average)
            temp_array=temp_array.reshape((1,temp_array.shape[0]))
            global_list.append(temp_array)
            break        
        elif rowline[filter_order] !=np_array[idx+1,filter_order]:
            new_array=np_array[sep_counts:idx+1,:]
            
            temp_array=new_array.mean(axis=axis_average)
            temp_array=temp_array.reshape((1,temp_array.shape[0]))
            global_list.append(temp_array)
            sep_counts=idx+1

            
        idx=idx+1
    
    np_out_put=concat_array_from_list(list_of_array=global_list,
                                      concat_axis=axis_average) 
    np_out_put=np_out_put.astype(astype)       

    return np_out_put

def concat_array_from_list (list_of_array , concat_axis = 0) :
    """ Concat array from list and set the None value in the list as NaN.
    
    :param list_of_array: List of array elements 
    :type list of array: list 
    
    :param concat_axis: axis for concatenation ``0`` or ``1``
    :type concat_axis: int 
    
    :returns: Concatenated array with shape np.ndaarry(
        len(list_of_array[0]), len(list_of_array))
    :rtype: np.ndarray 
    
    :Example: 
    >>> import numpy as np 
    >>> from pycsamt.utils.func_utils import concat_array_from_list 
    >>> np.random.seed(0)
    >>> ass=np.random.randn(10)
    >>> ass = ass2=np.linspace(0,15,10)
    >>> concat_array_from_list ([ass, ass]) 
    
    """
    concat_axis =int(_assert_all_types(concat_axis, int, float))
    if concat_axis not in (0 , 1): 
        raise ValueError(f'Unable to understand axis: {str(concat_axis)!r}')
    
    list_of_array = list(map(lambda e: np.array([np.nan])
                             if e is None else np.array(e), list_of_array))
    # if the list is composed of one element of array, keep it outside
    # reshape accordingly 
    if len(list_of_array)==1:
        ar = (list_of_array[0].reshape ((1,len(list_of_array[0]))
                 ) if concat_axis==0 else list_of_array[0].reshape(
                        (len(list_of_array[0]), 1)
                 )
             ) if list_of_array[0].ndim ==1 else list_of_array[0]
                     
        return ar 

    #if concat_axis ==1: 
    list_of_array = list(map(
            lambda e:e.reshape(e.shape[0], 1) if e.ndim ==1 else e ,
            list_of_array)
        ) if concat_axis ==1 else list(map(
            lambda e:e.reshape(1, e.shape[0]) if e.ndim ==1 else e ,
            list_of_array))
                
    return np.concatenate(list_of_array, axis = concat_axis)


# @deprecated("Should be removed ")
# def concat_array_from_list2(list_of_array, concat_axis=0):
#     """
#     Small function to concatenate a list with array contents 
    
#     Parameters 
#     -----------
#         * list_of_array : list 
#                 contains a list for array data. the concatenation is possible 
#                 if an index array have the same size 
        
#     Returns 
#     -------
#         array_like 
#             numpy concatenated data 
        
#     :Example: 
        
#         >>> import numpy as np 
#         >>> np.random.seed(0)
#         >>> ass=np.random.randn(10)
#         >>> ass2=np.linspace(0,15,12)
#         >>> ass=ass.reshape((ass.shape[0],1))
#         >>>  ass2=ass2.reshape((ass2.shape[0],1))
#         >>> or_list=[ass,ass2]
#         >>> ss_check_error=concat_array_from_list(list_of_array=or_list,
#         ...                                          concat_axis=0)
#         >>>  secont test :
#         >>> ass=np.linspace(0,15,14)
#         >>> ass2=np.random.randn(14)
#         >>> ass=ass.reshape((ass.shape[0],1))
#         >>> ass2=ass2.reshape((ass2.shape[0],1))
#         >>> or_list=[ass,ass2]
#         >>>  ss=concat_array_from_list(list_of_array=or_list, concat_axis=0)
#         >>> ss=concat_array_from_list(list_of_array=or_list, concat_axis=1)
#         >>> ss
#         >>> ss.shape 
#     """
#     #first attemp when the len of list is ==1 :
    
#     if len(list_of_array)==1:
#         if type(list_of_array[0])==np.ndarray:
#             output_array=list_of_array[0]
#             if output_array.ndim==1:
#                 if concat_axis==0 :
#                     output_array=output_array.reshape((1,output_array.shape[0]))
#                 else :
#                     output_array=output_array.reshape((output_array.shape[0],1))
#             return output_array
        
#         elif type(list_of_array[0])==list:
#             output_array=np.array(list_of_array[0])
#             if concat_axis==0 :
#                 output_array=output_array.reshape((1,output_array.shape[0]))
#             else :
#                 output_array=output_array.reshape((output_array.shape[0],1))
#             return output_array
    
#     # check the size of array in the liste when the len of list is >=2
    
#     for ii,elt in enumerate(list_of_array) :
#         if type(elt)==list:
#             elt=np.array(elt)
#         if elt is None :
#             pass
#         elif elt.ndim ==1 :
#             if concat_axis==0 :
#                 elt=elt.reshape((1,elt.shape[0]))
#             else :
#                 elt=elt.reshape((elt.shape[0],1))
#         list_of_array[ii]=elt
 
#     output_array=list_of_array[0]
#     for ii in list_of_array[1:]:
#         output_array=np.concatenate((output_array,ii), axis=concat_axis)
        
#     return output_array


def sort_array_data(data,  sort_order =0,
              concatenate=False, concat_axis_order=0 ):
    """
    Function to sort array data and concatenate 
    numpy.ndarray 
    
    Parameters
    ----------
        * data : numpy.ndarray 
                must be in simple array , list of array and 
                dictionary whom the value is numpy.ndarray 
                
        * sort_order : int, optional 
                index  of colum to sort data. The default is 0.
                
        * concatenate : Boolean , optional
                concatenate all array in the object.
                Must be the same dimentional if concatenate is set to True. 
                The *default* is False.
                
        * concat_axis_order : int, optional
                must the axis of concatenation  . The default is axis=0.

    Returns
    -------
        numpy.ndarray
           data , Either the simple sort data or 
           array sorted and concatenated .
    """
    
    if type(data)==list :
        for ss, val in data :
            val=val[val[:,sort_order].argsort(kind="mergesort")]
            data[ss]=val
            
        if concatenate: 
            data=concat_array_from_list(list_of_array=data,
                        concat_axis=concat_axis_order)
    elif type(data)==dict:
        
        for key, value in data.items(): 
            value=value[value[:,sort_order].argsort(kind="mergesort")]
            data[key]=value
            
        if concatenate==True :
            temp_list=list(data.values())
            data=concat_array_from_list(list_of_array=temp_list,
                                    concat_axis=concat_axis_order)
    else : 
        data=data[data[:,sort_order].argsort(kind="mergesort")]
        
    return data 
        
@deprecated ("Function replaced to another {_search_ToFill_Data}")
def transfer_array_(data, index_key,start_value_depth, end_value_depth, 
                    column_order_selection=0, axis=0): 
    """
    Parameters
    ----------
        * data : dict
                Dictionnary of numpy ndarray .
            
        * index_key : float 
                key of the dictionnary . 
                Must be a number of the first column of offset .
            
        * start_value_depth : float 
                If the depth is not reach must add depth of the closest point.
                give the start value which match to the maxi depth of the data :
                The *default* is -214.
                
        * end_value_depth : float
                Maximum depth of the survey. The default is -904.
            
        * column_order_selection : int,
                the index of depth column. The default is 0.
                
        * axis : int , optional
                numpy.ndarray axis . The default is 0.

    Returns
    -------
        numpy.ndarray
            return the array data we want to top to .
    
    :Example: 
        
        >>> import numpy as np 
        >>> sos=abs(np.random.randn(4,3)*4)
        >>> sos2=abs(np.random.randn(4,3)*10.8)
        >>>  print(sos2)
        >>> sis1=sort_array_data(data=sos,sort_order =1,
        ...                    concatenate=False, concat_axis_order=0)
        >>> sis2=sort_array_data(data=sos2,sort_order =1,
        ...                    concatenate=False, concat_axis_order=0)
        >>> dico={"18.4":sis1,
        ...      "21.4":sis2}
        >>> test=transfer_array_(data=dico, index_key=11.4, 
        ...                      start_value_depth=-14, end_value_depth=23,
        ...                      column_order_selection=1)
        >>> print("sis1:",sis1)
        >>> print("sis2:",sis2)
        >>> print("Finaltest", test)
    """

    start_value_depth=abs(start_value_depth)
    end_value_depth=abs(end_value_depth)
    comp, flag,iter_,translist_=0,-1,0,[]
    
   # chef the depth colum if negative before enter in loop :  

    if type (data)==dict :
        for key, value in data.items(): # scroll the dictionary 
            if float(key)> float(index_key):   # check the key of dictionnary
                                        #before using its value 
                maxi_depth=abs(value[:,column_order_selection]).max() # if yes
                    # calculate the max depth of the value :
                comp=0
                for rowline in value : # scroll the row of the array dict value
                
                    if abs(rowline[column_order_selection])>start_value_depth : # check its
                        # print(abs(rowline[column_order_selection]))

                        # value and compare it to section we must start extract (start)
                        if end_value_depth > maxi_depth:
                            transData_=value[comp:,:] # if yes transfer data :
                            # print(transData_)
                            translist_.append(transData_)
                            
                            start_value_depth=abs(transData_[:,column_order_selection]).max()
                            # start_value_depth=maxi_depth
                            
                            comp,iter_=0,iter_+1
                            flag=1
                            
                        elif end_value_depth<=maxi_depth:
                            indix =abs(value[:,column_order_selection]).tolist().index(end_value_depth)
                            if iter_==0 : 
                                transData_=value[comp:indix,::] 
                                flag=0
                                break 
                            else :
                                transData_=value[comp:indix,::]
                                translist_.append(transData_)
                                flag=1
                                
                        if start_value_depth >= end_value_depth:
                            flag=1
                            break
                    else :
                        comp +=1 
                        
    if flag ==0 : 
        return transData_
    if flag==1 : 
        trans_array=concat_array_from_list(list_of_array=translist_, concat_axis=axis)
        return trans_array 

            
    if type(data)==list: # create a dictionnary of array value 
        dico={}
        for ss, value  in enumerate( data) :
            dico["{0}".format(value[0,0])]=value
        # call recursive function 
        transfer_array_(data=dico,index_key=index_key,start_value_depth=start_value_depth,
                        end_value_depth=end_value_depth,column_order_selection=column_order_selection, axis=0)
        
            
def interpol_scipy (x_value, y_value,x_new,
                    kind="linear",plot=False, fill="extrapolate"):
    
    """
    function to interpolate data 
    
    Parameters 
    ------------
        * x_value : np.ndarray 
                    value on array data : original absciss 
                    
        * y_value : np.ndarray 
                    value on array data : original coordinates (slope)
                    
        * x_new  : np.ndarray 
                    new value of absciss you want to interpolate data 
                    
        * kind  : str 
                projection kind : 
                    maybe : "linear", "cubic"
                    
        * fill : str 
            kind of extraolation, if None , *spi will use constraint interpolation 
            can be "extrapolate" to fill_value.
            
        * plot : Boolean 
            Set to True to see a wiewer graph

    Returns 
    --------
        np.ndarray 
            y_new ,new function interplolate values .
            
    :Example: 
        
        >>> import numpy as np 
        >>>  fill="extrapolate"
        >>>  x=np.linspace(0,15,10)
        >>>  y=np.random.randn(10)
        >>>  x_=np.linspace(0,20,15)
        >>>  ss=interpol_Scipy(x_value=x, y_value=y, x_new=x_, kind="linear")
        >>>  ss
    """
    
    func_=spi.interp1d(x_value, y_value, kind=kind,fill_value=fill)
    y_new=func_(x_new)
    if plot :
        plt.plot(x_value, y_value,"o",x_new, y_new,"--")
        plt.legend(["data", "linear","cubic"],loc="best")
        plt.show()
    
    return y_new


def _set_depth_to_coeff(data, depth_column_index,coeff=1, depth_axis=0):
    
    """
    Parameters
    ----------
        * data : np.ndarray
            must be on array channel .
            
        * depth_column_index : int
            index of depth_column.
            
        * depth_axis : int, optional
            Precise kind of orientation of depth data(axis =0 or axis=1) 
            The *default* is 0.
            
        * coeff : float,
            the value you want to multiplie depth. 
            set depth to negative multiply by one. 
            The *default* is -1.

    Returns
    -------
        data : np.ndarray
            new data after set depth according to it value.
            
    :Example: 

        >>>  import numpy as np 
        >>>  np.random.seed(4)
        >>>  data=np.random.rand(4,3)
        >>>  data=data*(-1)
        >>>  print("data\n",data)
        >>>  data[:,1]=data[:,1]*(-1)
        >>>  data[data<0]
        >>>  print("data2\n",data)
    """
    
    if depth_axis==0:
        data[:,depth_column_index]=data[:,depth_column_index]*coeff
    if depth_axis==1:
        data[depth_column_index,:]=data[depth_column_index,:]*coeff  

    return data
            


def broke_array_to_(arrayData, keyIndex=0, broken_type="dict"):
    """
    broke data array into different value with their same key 

    Parameters
    ----------
        * arrayData :np.array
            data array .
            
        * keyIndex : int 
            index of column to create dict key 

    Returns
    -------
        dict 
           dico_brok ,dictionnary of array.
    """
    
    # find the max_counts
    
    vcounts_temp,counts_temp=np.unique(arrayData[:,keyIndex], return_counts=True)
    vcounts_temp_max=vcounts_temp.max()
    # print(vcounts_temp)
    # print(vcounts_temp_max)
    # print(vcounts_temp.min())

    dico_brok={}
    lis_brok=[]
    index=0
    deb=0
    for rowlines in arrayData :
        if rowlines[0] == vcounts_temp_max:
            value=arrayData[index:,::]
            if broken_type=="dict":
                dico_brok["{0}".format(rowlines[0])]=value
                break
            elif broken_type=="list":
                lis_brok.append(value)
                break
        elif rowlines[0] !=arrayData[index+1,keyIndex]:
            value=arrayData[deb:index+1,::]
            if broken_type=="dict":
                dico_brok["{0}".format(rowlines[0])]=value
            elif broken_type=="list":
                lis_brok.append(value)
            deb=index+1
        index=index+1
    if broken_type=="dict":
        return dico_brok
    elif broken_type=="list":
        return lis_brok
    

@deprecated('this function is replaced by [_search_ToFill_Data] ')
def _OlDFUNCNOUSEsearch_fill_data(dicoReal, arrayTemp , 
                      max_value, index_of_depth,axis=0):
    """
    Deprecated function , very expensive.
    
    Parameters
    ----------
        * data : dict
                Dictionnary of numpy ndarray .
            
        * dataReal : dict
                 dictionnary . must be a dictionnary of real of offset .
            
        * arrayTemp : np.ndarray 
                 must be a numpy array of reserve data , the one , we want to 
                 extract the depth data to fill array 
    
        * max_value : float
                Maximum depth of the survey. 
            
        * index_of_depth : int,
                the index of depth column. The *default* is 0.
                
        * axis : int , optional
            numpy.ndarray axis . The default is 0.
    
    Returns
    -------
        array_like 
            the array data we want to top to .
    
    :Example: 
        
         >>> import numpy as np 
         >>>  np.random.seed(0)
         >>>  sos=abs(np.random.randn(4,3)*4)
         >>>  sos2=abs(np.random.randn(4,3)*10.8)
         >>>  sos3=abs(np.randon.rand(8,3)*12.4)
         >>>  # print(sos2)
         >>>  sis1=sort_array_data(data=sos,sort_order =1,
         ...                    concatenate=False, concat_axis_order=0)
         >>>  sis2=sort_array_data(data=sos2,sort_order =1,
         ...                    concatenate=False, concat_axis_order=0)
         >>>  sis3=sort_array_data(data=sos3,sort_order =1,
         ...                    concatenate=False, concat_axis_order=0)
         >>>  dico={"18.4":sis1,
         ...      "21.4":sis2}
         >>>  test=_search_fill_data(dicoReal=dico, index_key=11.4, 
         ...                      start_value=10, max_value=23,
         ...                      index_of_depth=1)
         >>> print("sis1\n:",sis1)
         >>> print("sis2\n:",sis2)
         >>> print("Finaltest\n", test)
    """
    
    # arange a dictionany : 
    #from keys : for k in sorted(dico.keys()):
    #   print(%s: %s",%(k,names[k]))
    # from values : for k , v  in sorted(dico.items(),key=lambda x:x[1]):
    #   print(%s: %s",%(k,v))
    
    _all_checker,keyToSkip=[],[]
    litemp= broke_array_to_(arrayData=arrayTemp,keyIndex=0,
                            broken_type="list")
    itemp,real,realTem=[],[],[]
    for ii in litemp:
        temp=sort_array_data(data=ii,sort_order=1,
                             concatenate=False)
        itemp.append(temp)
    # print(itemp)
        
    for key , value in sorted(dicoReal.items()):
        # print(sorted(dicoReal.items()))
        # for ii, rowline in enumerate(value)  :
        rowmax=value[:,index_of_depth].max()
        real.append((float(key),rowmax))
        _all_checker.append(rowmax)
        
    
    #chek if all elements are reach the depth max  
    # if all(_all_checker)==True : # if one of the depth is the same 
    #     return dicoReal
    realTem=real.copy()
    for ii , value in enumerate(realTem):
        if value[1] == max_value:
            del realTem[ii]
            #real.pop(ii)
            
    if realTem ==[]:
        return dicoReal
    
    
    print("Real : \n",real)
    idx,flag=0,0
    comp,sp=0,-1
    fin_list,endList=[],[]
    
    while idx < len(realTem):
        
        indexKey=realTem[idx][0] #11,4
        maxKey=realTem[idx][1] #214        
        # if real[idx][1]==max_value : # case where thedepth  value of realdico is get  
        #     flag=3
        # # elif real[idx][1]==max_value : # case where thedepth  value of realdico is get  
        # #     flag=3
        # else :
        if idx==len(realTem)-1:
            indic=real.index(realTem[idx])
            # check whether there is data after a delete the offset with depth value reach 
            if realTem[idx][0]==real[-1][0]:
                nexIndex=itemp[-1][0][0]
            else :
                nexIndex=real[indic+1][0]
        else :
            nexIndex=realTem[idx+1][0]#61.4
        # loop the reserve list :
        for ii, array in enumerate (itemp): #itemp is the  list of reserve broken list
        
            # if array[0][0]> indexKey and array[0][0]> nexIndex:
            #     keyToSkip.append(realTem[idx])
            #     continue 
            if array[0][0]> indexKey and array[0][0]<nexIndex : #check offset and next offset , if True :
                
                for index, rowline in enumerate(array) : # loop the array now 
                    # tem_depth=rowline[1] # take the value of depth of reserve array for the first row 
                    if rowline[1] > maxKey :# maxKey=214, rowline[1]= :# (max_value=904) if True :
                        sp=sp+1
                        if sp==0 :
                            add_array=array[index:,:]
                            # num=index
                            # print("True:\n",add_array)
                            flag=4
                        else :
                            sp=-1
                            pass
                        # print(maxtem)
            if flag==4:
                # maxtem=array[-1][1]       # take the maximum depth of the reserve array , last row
                if array[-1][1]  < max_value : # if the maximum depth not reach 
                    comp=comp+1
                    if array[0][0]<nexIndex :
                        endList.append(add_array)
                    else :
                        fin_list.append(add_array) # keep the array in temporary list 
                    # maxKey = maxtem   # set up new maximum depth 
                    # flag=1
                    maxKey = array[-1][1] 
                    # print("flag4,comp>1")
                    # if maxKey < 
                elif array[-1][1]  ==max_value : # the maximum depth is reached
                    flag=5
    
            if flag==5 :
                
                if comp==0 : # first check is ok 
                    # add_array=array[num:,:] # cut the array 
                    endList.append(add_array)
                    flag=0
                    # print("flag5,comp=0")
                else :
                  add_array=array[index:,:] #  
                  fin_list.append(add_array) # list to create one ar

                  flag=1
                      
        # if flag==3 : # in that case is true , save the value of the offset 
        #     # for key , value in dicoReal.items() :
        #     #     if float(key) == real[idx][0]:
        #     #         endList.append(value)
        #     keyToSkip.append(array[0][0])
            
                    
        if flag==1 :
            arT=concat_array_from_list(list_of_array=fin_list,concat_axis=axis)
            endList.append(arT)
        if flag==0 :
            endList=endList
        
        idx=idx+1
            
    print("keyToSkip\n:",keyToSkip)       
        
    
    #delete the the offset which are full depth on the list 
    # print(keyToSkip)
    # print(real,"\n")
    for ii , value in enumerate(realTem):
        if value[0] in keyToSkip:
            del realTem[ii]
            #real.pop(ii)
    
    # now we are the list of recoverd depth 
    # build dictionnary
    print("Realtem\n",realTem,)

    print("endlist : \n",endList)
    # print(keyToSkip)
    
    for ii , tuple_ in enumerate (realTem) : #take the realkey from dico
        realKey, maxValue=tuple_
        print(realKey,maxValue)
        
        add_array=endList[ii]       # take the add_value generated
        # print(endList)
        for key , value in dicoReal.items(): # search in dictionnary the key
            if float(key)==realKey: # and compare it to the key from tuple ...
                        #... just to have certitude then concatenate 
                val=np.concatenate((value,add_array),axis=0)
                dicoReal[key]=val
        # print(ii)
    return dicoReal 


def _search_ToFill_Data (dicoReal, arrayTemp , 
                      max_value, index_of_depth,axis=0): 
    """
    Parameters 
    ------------
        * data : dict
                Dictionnary of numpy ndarray .
                
        * dicoReal : dict
             dictionnary . must be a dictionnary of real of offset .
            
        * arrayTemp : np.ndarray 
                 must be a numpy array of reserve data , the one , we want to 
                 extract the depth data to fill array 
    
        * max_value : float
                Maximum depth of the survey. 
            
        * index_of_depth : int,
                the index of depth column. The default is 0.
        * axis : int , optional
            numpy.ndarray axis . The default is 0.
        
    Returns
    -------
        dict
         dictionnary of offsets filled the array data we want to top to
        
    :Example: 
    
        >>>  import numpy as np 
        >>>  np.random.seed(0)
        >>>  sos=abs(np.random.randn(4,3)*4)
        >>>  sos2=abs(np.random.randn(4,3)*10.8)
        >>>  sos3=abs(np.random.randn(5,3)*10.8)
        >>>  temp3=abs(np.random.rand(4,3)*12.4)
        >>>  temp2=abs(np.random.rand(4,3)*15.4)
        >>>  temp1=abs(np.random.rand(4,3)*9.4)
        >>>  temp4=abs(np.random.rand(5,3)*9.9)
        >>>  dico,temp={},[]
        >>>  ff=[sos,sos2,sos3,temp1,temp2,temp3,temp4]
        >>>  fin=[sort_array_data(data=ii,sort_order =1,
        ...                concatenate=False, concat_axis_order=0) for ii in  ff ]
        >>>  key=[11.9,61.4,102.7]
        >>>  vat=[214,405,904]
        >>>  for ii in range(3):
        ...        fin[ii][:,0]=key[ii]
        ...        fin[ii][-1][1]=vat[ii]
        >>>  tempi=[(19.4,[11,18,50,120]),(28.4,[12,17,403,904]),
        ...       (78.3,[11,8,202,804]),(124.4,[203,403,604,714,904])]
        >>>  for ss, val in enumerate(tempi) :
        ...        fin[ss+3][:,0]=val[0]
        ...        fin[ss+3][:,1]=np.array(val[1])       
        >>>  for ss, van in enumerate (fin):
        ...        if ss<=3 :
        ...            dico[van[0][0]]=van
        ...        if ss>3 :
        ...            temp.append(van)
        >>>  arrayTemp=concat_array_from_list(list_of_array=temp,
        ...                                 concat_axis=0)
        >>>  sis02=_search_ToFill_Data(dicoReal=dico, arrayTemp=arrayTemp , 
        ...              max_value=904, index_of_depth=1,axis=0)
        >>>  print(sis02)
    """ 
    #Notes :
    # arange a dictionany : 
    #from keys : for k in sorted(dico.keys()):
    #   print(%s: %s",%(k,names[k]))
    # from values : for k , v  in sorted(dico.items(),key=lambda x:x[1]):
    #   print(%s: %s",%(k,v)
    # for ii , value in enumerate(realTem):
    # if value[0] in keyToSkip:
    #     del realTem[ii]
    #     #real.pop(ii)
    #------------------------------------
    # _all_checker,keyToSkip=[],[]
    litemp= broke_array_to_(arrayData=arrayTemp,keyIndex=0,
                            broken_type="list")
    itemp,real=[],[],
    # flag=-1
    for ii in litemp:
        temp=sort_array_data(data=ii,sort_order=1,
                             concatenate=False)
        itemp.append(temp)
    # print(itemp)
    for key , value in sorted(dicoReal.items()):
        # print(sorted(dicoReal.items()))
        # for ii, rowline in enumerate(value)  :
        rowmax=value[:,index_of_depth].max()
        real.append((float(key),rowmax))
        # _all_checker.append(rowmax)
        
    # print(itemp)
    iter_=-1  
    for ss,(offs, maxDepth) in enumerate( real):
        # [(11.4,214),(19.4,120),(61.4,405),(102,904)]
        # print(offs)
        if maxDepth < max_value : # max_value=904 , maxDepth=214
            if ss == len(real)-1:
                for num, arTem in enumerate(itemp):
                    kk=arTem[0][0]
                    if kk > offs :
                        nextIndex=kk
                    else :
                        break
                # nextIndex=itemp[-1][0][0]
            else :
                nextIndex=real[ss+1][0]
                
            for ii, array in enumerate(itemp) : # [array(28),array(78), ....]
            
                if offs < array[0][0]  and array[0][0] < nextIndex: 

                    for index, rowline in enumerate(array) : # scroll array(28)
                        if rowline[1]> maxDepth :  # if rowline > the maxDepth=214
                            iter_=iter_+1 # first iteration on array (28)
                            if iter_==0 : # collect the remain info and add to dico
                                add_array=array[index:,:]
                                maxDepth=array[:,1].max() # caluclate nnex
                                # print(add_array)
                            else :
                                pass
                            if iter_==0 and maxDepth<=max_value :
                                for key, value in dicoReal.items() :
                                    if float(key) ==offs :
                                        new_value=concat_array_from_list(list_of_array=[value,add_array],
                                                                     concat_axis=axis)
                                        dicoReal[key]=new_value
                                        iter_=-1
                                        # break
                else :
                    pass
             
    return dicoReal
        
def straighten_out_list (main_list , list_to_straigh):
    """
    Parameters
    ----------
        * main_list : list
                list of which the data must absolutely appear into 
                the straighen list.in our case , it is the station 
                list : a list of offset 
                
        * list_to_straigh : list
                list contain the data (offset calculated,
                 the depth and the resistivity (log10)), 

    Returns
    -------
        * list
            the straighen list.
            some offset have been replaced by the offsets which are not in the 
            main_list whithout change the lengh of the straighen list. 

    :Example: 
        
        >>>  import numpy as np 
        >>>  np.random.seed(14)
        >>>  ss=np.random.randn(10)*12
        >>>  ss=ss.tolist()
        >>>  ss=[round(float(jj),4) for jj in ss]
        >>>  ss.sort()
        >>>  red=np.random.randn(7)*12
        >>>  red=red.tolist()
        >>>  test=[19, 15.012, 5.5821, 0.7234,3.1, 
        ...      0.7919, 3.445, 4.7398, 5.1, 10.8, 15.51,21]
        >>>  main=[20., 0.7234, 5, 3.445, 15.51,10.7, 3,5.1]
        >>>  test.sort()
        >>>  main.sort()
        >>>  red=[round(float(ss),1) for ss in red]
        >>>  print(test)
        >>>  print(main)
        >>>  sos=straighten_out_list (main_list=main , 
        ...                         list_to_straigh=test)
        >>>  print("sos:\n",sos)   
    """
    
    
    staIter=np.array(list_to_straigh)# ARRAY 
    max_X=staIter.max()
    valuesIter,counts=np.unique(staIter, return_counts=True)
    valuesIter=valuesIter.tolist() # List sorted of straighten value
    # print(valuesIter)
    

    misoffs,vamin,tempi=[],[],[]
    # min_,flag=1,0
    for sta in main_list: # keep the offset not in list to straight 
        if sta not in valuesIter :
            misoffs.append(sta)
    if misoffs ==[]:
        return list_to_straigh
    
    # inject miss offset value in the temporary list and sorted 
    #this to choose the closet value we want to substitude the missoffset value 
    tempi=[ii for ii in valuesIter]
    for ii in misoffs :
        tempi.append(ii)

    tempi=deepcopy(tempi)

    tempi.sort(reverse =False)
    # print("mainlist:\n",main_list)
    # print("tempiSorted:\n",tempi)
    # print("missoff:\n",misoffs)
    # print("lisofStraight_max:\n",max_X)
    # # print(len(misoffs))

    # sp=-1
    for jj, off in enumerate(tempi): #scroll the create temporary list
        for ss, ofi in enumerate(misoffs):  # scroll the miss offset list and compare it 
                            # to the temportary list 
            if off==ofi:
                # print(off,ss)
                indexof =main_list.index(off) # find the index of "off" value in the mainlist
                        # in order to calculate the distance between the value to its next value 
                        # or forward value ex :off=3 , backward =0.724 , forward =3.445
                        
                # if main_list[indexof]==main_list[-1]: # that's mean the "off" value reaches the end of mainlist
                #     if tempi[-1]>main_list[indexof]: # 
                #         deltaMain=main_list[indexof]-tempi[-1] # take the diff between
                #         vamin.append(tempi[-1])
                #     else :
                #         deltaMain=main_list[indexof]-tempi[-2]
                # else :
                #     deltaMain=main_list[indexof]-main_list[indexof+1]
                    
                # if main_list[indexof+1] != tempi[jj+1] :
                if off==tempi[-1]:
                    # deltaOffneg=abs(off-tempi[jj-1])
                    vamin.append(tempi[jj-1])
                    # delta=deltaOffneg
                    
                elif off==tempi[0]:
                    deltaOffpos=abs(off-tempi[jj+1])
                    vamin.append(tempi[jj+1])
                    # delta=deltaOffpos
                else :
                    deltaOffpos=abs(off-tempi[jj+1])
                    deltaOffneg=abs(off-tempi[jj-1])
                    if deltaOffpos >= deltaOffneg :
                        # delta=deltaOffneg
                        vamin.append(tempi[jj-1])
                        # index=jj-1
                    else :
                        # delta=deltaOffpos
                        # index=jj+1
                            vamin.append(tempi[jj+1])

    # print("vamin:\n",vamin)   
    # print(len(vamin))
    for ii, val in enumerate(list_to_straigh):
        for ss, of in enumerate(vamin):
            if val==of :
                list_to_straigh[ii]=misoffs[ss]
                
    # if list_to_straigh[-1] !=max_X :
        
    return list_to_straigh

def take_firstValue_offDepth(data_array,
                             filter_order=1):
    """
    Parameters
    ----------
        * data_array : np.array 
                array of the data .
        * filter_order : int , optional
                the column you want to filter. The default is 1.

    Returns
    -------
        array_like
            return array of the data filtered.
   
    :Example: 
        
        >>>  import numpy as np 
        >>>  list8=[[4,2,0.1],[8,2,0.7],[10,1,0.18],[4,3,0.1],
        ...        [7,2,1.2],[10,3,0.5],[10,1,0.5],[8.2,0,1.9],
        ...        [10,7,0.5],[10,1,0.5],[2,0,1.4],[5,4,0.5],
        ...        [10,2,0.7],[7,2,1.078],[10,2,3.5],[10,8,1.9]]
        >>>  test=np.array(list8)
        >>>   print(np_test)
        >>>  ss=take_firstValue_offDepth(data_array =np_test, filter_order=1)
        >>>   ss=averageData(np_array=np_test,filter_order=1,
        >>>                 axis_average=0, astype="int")
        >>> print(ss)
    """
    
    listofArray=[]#data_array[0,:]]
    data_array=data_array[data_array[:,filter_order].argsort(kind="mergesort")]
    # print(data_array,"\n:")
    # np_array=np_array[np_array[:,filter_order].argsort(kind="mergesort")]
     #returns the differents value on the filtersort index of array :
    values, counts =np.unique(data_array[:,filter_order], return_counts=True)
    
    for ii, rowline in enumerate(data_array ): 
    
        if rowline[filter_order]==values[-1]:
            listofArray.append(data_array[ii])
            break 
        elif rowline[filter_order] !=data_array[ii-1][filter_order]:
            listofArray.append(data_array[ii])
        
        
    array =concat_array_from_list(list_of_array=listofArray, concat_axis=0)
    array=array[array[:,filter_order].argsort(kind="mergesort")]
    listofArray=[]

    return array 

def dump_comma(input_car, max_value=2, carType='mixed'):
    """
    Parameters
    ----------
        * input_car : str,
            Input character.
        * max_value : int, optional
            The default is 2.
            
        * carType: str 
            Type of character , you want to entry
                 
    Returns
    -------
        Tuple of input character
            must be return tuple of float value, or string value
      
    .. note:: carType  may be as arguments parameters like ['value','val',"numeric",
              "num", "num","float","int"] or  for pure character like 
                ["car","character","ch","char","str", "mix", "mixed","merge","mer",
                "both","num&val","val&num&"]
                if not , can  not possible to convert to float or integer.
                the *defaut* is mixed 
                
    :Example: 
        
        >>> import numpy as np
        >>>  ss=dump_comma(input_car=",car,box", max_value=3, 
        ...      carType="str")
        >>>  print(ss)
        ... ('0', 'car', 'box')
    """
    
    
    
        # dump "," at the end of 
    flag=0
    
    if input_car[-1]==",":
        input_car=input_car[:-1]
    if input_car[0]==",":
        input_car="0,"+ input_car[1:]
        
    if carType.lower() in ['value','val',"numeric",
                           "num", "num","float","int"]:
        input_car=eval(input_car)
        
    elif carType.lower() in ["car","character","ch","char","str",
                             "mix", "mixed","merge","mer",
                             "both","num&val","val&num&"]:
        
        input_car=input_car.strip(",").split(",")
        flag=1
        
    # elif carType.lower() in ["mix", "mixed","merge","mer",
    #                          "both","num&val","val&num&"]:
    #     input_car=input_car.split(",")
    
        
    if np.iterable(input_car)==False :
        inputlist=[input_car,0]
        # input_car=tuple(inputlist)
    elif np.iterable(input_car) is True :
        # if flag==1 :
        #     inputlist=input_car
        #     # print(inputlist)
        # else :
        inputlist=list(input_car)
            
        # print("false")
    
    input_car=inputlist[:max_value]
    # print(input_car)
    if flag==1 :
        if len(inputlist)==1 :
            return(inputlist[0])
    
    return tuple(input_car)

                    
def build_wellData (add_azimuth=False, utm_zone="49N",
                    report_path=None,add_geochemistry_sample=False):
    """
    Parameters
    ----------
        * add_azimuth : Bool, optional
                compute azimuth if add_azimut is set to True. 
                The default is False.
             
        *  utm_zone : Str, optional
                 WGS84 utm_projection. set your zone if
                 add_azimuth is turn to True. 
                 The default is "49N".
             
        * report_path : str, optional
                path to save your _well_report. The default is None.
                its match the current work directory 
            
        * add_geochemistry_sample: bool
                add_sample_data.Set to True if you want to add_mannually
                 Geochimistry data. default is False.

    Raises
    ------
        Exception
            manage the dimentionaly of ndarrays .
        OSError
            when report_path is not found in your O.S.

    Returns
    -------
        str
            name of location of well .
        np.ndarray
             WellData , data of build Wells .
        np.ndarray
           GeolData , data of build geology.

    :Example: 
        
        >>>  import numpy as np 
        >>>  import os, shutil
        >>>  import warnings,
        >>>  form _utils.avgpylog import AvgPyLog
        >>>  well=build_wellData (add_azimuth=True, utm_zone="49N")
        >>>  print("nameof locations\n:",well[0])
        >>>  print("CollarData\n:",well[1])
        >>>  print("GeolData\n:", well[2])
        ...  nameof locations
        ...  Shimen
        ...  CollarData
        ...  [['S01' '477205.6935' '2830978.218' '987.25' '-90' '0.0' 'Shi01'
        ...   'Wdanxl0']
        ...   ['S18' '477915.4355' '2830555.927' '974.4' '-90' '2.111' 'Shi18'
        ...   'Wdanxl0']]
        ...  GeolData
        ...  [['S01' '0.0' '240.2' 'granite']
        ...   ['S01' '240.2' '256.4' ' basalte']
        ...   ['S01' '256.4' '580.0' ' granite']
        ...   ['S01' '580.0' '987.25' 'rock']
        ...   ['S18' '0.0' '110.3' 'sand']
        ...   ['S18' '110.3' '520.2' 'agrilite']
        ...    ['S18' '520.2' '631.3' ' granite']
        ...   ['S18' '631.3' '974.4' ' rock']]
        ...   Shimen_wellReports_
    """
    reg_lines=[]
    wellSites,ftgeo,hole_list,Geolist=[],[],[],[]
    
    text=["Enter the name of Location:",
                      "well_name :",
                      "Coordinates (Easting, Northing)_UTM_{0} : ".format(utm_zone),
                      "Hole Buttom  and dip values (Bottom, dip):" ,
                      "Layers-thickness levels in (meters):",
                      "Geology-layers or <stratigraphy> names (Top_To_Buttom):", 
                      "{0:-^70}".format(' Report '),
                      
                      "DH_Hole,DH_Easting, DH_Northing, DH_Buttom,"\
                      " DH_Dip,DH_Azimuth, DH_PlanDepth,DH_Descrip",
                      
                      "GeolData :",
                      "WellData:",
                      "DH_Hole, DH_From, DH_To, Rock",
                      "SampleData",
                      "DH_Hole, DH_From,DH_To, Sample",
                      "{0:-^70}".format(' InputData '),
                      ]
    
    
    name_of_location =input("Enter the name of Location:")
    reg_lines.append(''.join(text[0]+'{0:>18}'.format(name_of_location)+'\n'))
    reg_lines.append('\n')
    reg_lines.append(text[13]+'\n')
    
    comp=-1
    while 1 :
        DH_Hole=input("Enter the well_name :")
        if DH_Hole=="" or DH_Hole=="end":
            Geol=concat_array_from_list(list_of_array=Geolist,
                                            concat_axis=0)
            break
         
        print("Enter the coordinates (Easting, Northing) : ", end="")
        DH_East_North=input()
        DH_East_North=dump_comma(input_car=DH_East_North,
                                 max_value=2, carType='value')
        print("Enter the Hole Bottom value and dip (Bottom, dip):", end='')
        dh_botdip=input()
        dh_botdip=dump_comma(input_car=dh_botdip,
                                 max_value=2, carType='value')
        #check  the dip of the well 
        if float(dh_botdip[1])==0.:
            dh_botdip[1]=(90.)
        elif float(dh_botdip[0])==0.:
            raise Exception(
            "The curent bottom has a value 0.0 . "
            "Must put the bottom of the well as deep as possible !"
                            )
            
        hole_list.append(DH_Hole)

        wellSites.append((DH_Hole, DH_East_North[0],DH_East_North[1],
                          dh_botdip[0],dh_botdip[1] ))

        #DH_Hole (ID)	DH_East	DH_North	DH_RH	DH_Dip	DH_Azimuth	DH_Top	DH_Bottom	DH_PlanDepth	DH_Decr	Mask 
        reg_lines.append("".join(text[1]+'{0:>7}'.format(DH_Hole))+"\n")
        reg_lines.append("".join(text[2])+"\n")
        reg_lines.append(",".join(['{0:>14}'.format(str(ii))
                                   for ii in list(DH_East_North)])+"\n")
        reg_lines.append("".join(text[3])+"\n")
        reg_lines.append(",".join(['{0:>7}'.format(str(ii))
                                   for ii in list(dh_botdip)])+"\n")
        
        comp+=-1
        while DH_Hole :
            # print("Enter the layer thickness (From_, _To, geology):",end='')
            if comp==-1 :
                Geol=concat_array_from_list(list_of_array=ftgeo,
                                            concat_axis=0)
                
                ftgeo=[]        #  initialize temporary list 
                
                break
            # comp=comp+1
            print("Enter the layers-thickness levels in (meters):", end='')
            dh_from_in=input()
            
            if dh_from_in=="" or dh_from_in=="end":
                break

            dh_from=eval(dh_from_in)
            
            dh_from_ar=np.array(dh_from)
            dh_from_ar=dh_from_ar.reshape((dh_from_ar.shape[0],1)) 
            # check the last input bottom : 
            if dh_from_ar[-1] >= float(dh_botdip[0]):
                _logger.info("The input bottom of well {0}, is {1}. "
                             "It's less last layer thickess: {2}."
                             "we add maximum bottom at 1.023km depth.".format(
                                 DH_Hole,dh_botdip[0],dh_from_ar[-1]))
                dh_botdip[0]=(1023.)
                wellSites[-1][3]=dh_botdip[0]  # last append of wellSites 
            #find Dh_to through give dh_from
            
            dh_to=dh_from_ar[1:]
            dh_to=np.append(dh_to,dh_botdip[0])
            dh_to=dh_to.reshape((dh_to.shape[0],1))
            
            print("Enter the geology-layers names (From _To):",end="")
            rock=input()
            rock=rock.strip(",").split(",") # strip in the case where ","appear at the end 
            rock_ar=np.array(rock)
            rock_ar=rock_ar.reshape((rock_ar.shape[0],1))

            try :
                if rock_ar.shape[0]==dh_from_ar.shape[0]:
                    drill_names=np.full((rock_ar.shape[0],1),DH_Hole)
                fromtogeo=np.concatenate((drill_names,dh_from_ar,
                                              dh_to, rock_ar),axis=1)
            except IndexError:
                _logger.warn("np.ndarry sizeError:Check 'geologie', 'Dh_From', and "
                             "'Dh_To' arrays size properly. It seems one size is "
                              " too longeR than another. ")
                warnings.warn(" All the arrays size must match propertly!")
                
            ftgeo.append(fromtogeo)
            comp=-1
            
            reg_lines.append("".join(text[4])+"\n")
            reg_lines.append(",".join(['{0:>7}'.format(str(jj)) 
                                       for jj in list(dh_from)])+"\n")
            reg_lines.append("".join(text[5])+"\n")
            reg_lines.append(",".join(['{0:>12}'.format(jj) 
                                       for jj in rock])+"\n") # rock already in list 
            # reg_lines.append("".join(rock+"\n"))
            
        reg_lines.append("\n")
        
        Geolist.append(Geol)    


    name_of_location=name_of_location.capitalize()
    #set on numpy array
    
    for ii , value in enumerate(wellSites):
        value=np.array(list(value))
        wellSites[ii]=value
    #create a wellsites array 
    
    wellSites=concat_array_from_list(wellSites,concat_axis=0)
    DH_Hole=wellSites[:,0]
    DH_PlanDepth=np.zeros((wellSites.shape[0],),dtype="<U8")
    DH_Decr=np.full((wellSites.shape[0],),"Wdanxl0",dtype="<U9")
    
    
    lenloc=len(name_of_location)
    for ii , row in enumerate(DH_Hole):
         # in order to keep all the well name location
        DH_PlanDepth[ii]=np.array(name_of_location[:-int(lenloc/2)]+row[-2:])

    Geol[:,1]=np.array([np.float(ii) for ii in Geol[:,1]])
    Geol[:,2]=np.array([np.float(ii) for ii in Geol[:,2]])
    
    if add_azimuth==False:
        DH_Azimuth=np.full((wellSites.shape[0]),0)
    elif add_azimuth == True :
        DH_Azimuth=compute_azimuth(easting = np.array(
            [np.float(ii) for ii in wellSites[:,1]]),
                                   northing =np.array(
                                       [np.float(ii) for ii in wellSites[:,2]]),
                                   utm_zone=utm_zone)

        
    DH_Azimuth=DH_Azimuth.reshape((DH_Azimuth.shape[0],1))
    
    
    WellData=np.concatenate((wellSites,DH_Azimuth,
                             DH_PlanDepth.reshape((DH_PlanDepth.shape[0],1)),
                             DH_Decr.reshape((DH_Decr.shape[0],1))), axis=1)
    GeolData=Geol.copy()

    #-----write Report---
    
    reg_lines.append(text[6]+'\n')
    # reg_lines.append(text[7]+"\n")
    
    reg_lines.append(text[9]+'\n')
    reg_lines.append("".join(['{0:>12}'.format(ss) 
                              for ss in text[7].split(",")]) +'\n')
    for rowline in WellData :
        reg_lines.append(''.join(["{0:>12}".format(ss)
                                  for ss in rowline.tolist()])+"\n")
        
    reg_lines.append(text[8]+"\n")
    reg_lines.append("".join(['{0:>12}'.format(ss) 
                              for ss in text[10].split(",")]) +'\n')
    for ii , row in enumerate(GeolData):
        reg_lines.append(''.join(["{0:>12}".format(ss) 
                                  for ss in row.tolist()])+"\n")
        
    if add_geochemistry_sample==True:
        SampleData=build_geochemistry_sample()
        reg_lines.append(text[11]+'\n')
        reg_lines.append("".join(['{0:>12}'.format(ss) 
                                  for ss in text[12].split(",")]) +'\n')
        for ii , row in enumerate(SampleData):
            reg_lines.append(''.join(["{0:>12}".format(ss) 
                                      for ss in row.tolist()])+"\n")
    else :
        SampleData=None        
        

    with open("{0}_wellReport_".format(name_of_location),"w") as fid:
        # for ii in reg_lines :
        fid.writelines(reg_lines)
    fid.close()
    
    #---end write report---
    
    if report_path is None:
        report_path=os.getcwd()
    elif report_path is not None :
        if os.path.exists(report_path):
            shutil.move((os.path.join(os.getcwd(),"{0}_wellReport_".\
                                      format(name_of_location))),report_path)
        else :
            raise OSError (
                "The path does not exit.Try to put the right path")
            warnings.warn (
                "ignore","the report_path doesn't match properly.Try to fix it !")
        
    
    
    return (name_of_location, WellData , GeolData, SampleData)
        
        
def compute_azimuth(easting, northing, utm_zone="49N", extrapolate=False):
    
    """
    Parameters
    ----------
        * easting : np.ndarray
                Easting value of coordinates _UTM_WGS84 
             
        * northing : np.ndarray
                Northing value of coordinates._UTM_WGS84
            
        * utm_zone : str, optional
                the utm_zone . if None try to get is through 
                gis.get_utm_zone(latitude, longitude). 
                latitude and longitude must be on degree decimals. 
                The default is "49N".
        * extrapolate : bool , 
                for other purpose , user can extrapolate azimuth value, 
                in order to get the sizesize as 
                the easting and northing size. The the value will
                repositionate at each point data were collected. 
                    Default is False as originally azimuth computation . 

    Returns
    -------
        np.ndarray
            azimuth.
        
    :Example: 
        
        >>> import numpy as np
        >>> import gis_tools as gis
        >>>  easting=[477205.6935,477261.7258,477336.4355,477373.7903,477448.5,
        ...  477532.5484,477588.5806,477616.5968]
        >>>   northing=[2830978.218, 2830944.879,2830900.427, 2830878.202,2830833.75,
        ...                  2830783.742,2830750.403,2830733.734]
        >>>  test=compute_azimuth(easting=np.array(easting), 
        ...                      northing=np.array(northing), utm_zone="49N")
        >>>  print(test)
    """
    #---**** method to compute azimuth****----
    
    reference_ellipsoid=23
    
    lat,long=gis.utm_to_ll(reference_ellipsoid=reference_ellipsoid,
                                northing=northing, easting=easting, zone=utm_zone)
    
    #i, idx, ic_=0,0,pi/180
    azimuth=0
    
    i,ic_=0,np.pi /180
    
    while i < lat.shape[0]:
        xl=np.cos(lat[i]*ic_)*np.sin(lat[i+1]*ic_) - np.sin(lat[i]*ic_)\
            *np.cos(lat[i+1]*ic_)*np.cos((long[i+1]-long[i])*ic_)
        yl=np.sin((long[i+1]-long[i])*ic_)*np.cos(lat[i+1])
        azim=np.arctan2(yl,xl)
        azimuth=np.append(azimuth, azim)
        i=i+1
        if i==lat.shape[0]-1 :
            # azimuth.append(0)
            break
    
    
    if extrapolate is True : 
        # interpolate azimuth to find the azimuth to first station considered to 0.
    
        ff=spi.interp1d(x=np.arange(1,azimuth.size),
                                   y=azimuth[1:], fill_value='extrapolate')
        y_new , azim = ff(0),np.ones_like(azimuth)
        azim[0], azim[1:] = y_new , azimuth[1:]
    else : 
        azim=azimuth[1:] # substract the zero value added for computation as origin.
    #convert to degree : modulo 45degree
    azim = np.apply_along_axis(lambda zz : zz * 180/np.pi , 0, azim) 
    
    
    return np.around(azim,3)
        
def build_geochemistry_sample():
    """
    Build geochemistry_sample_data
    
    Raises
    ------
        Process to build geochemistry sample data manually .

    Returns
    -------
       np.ndarray
          Sample ,Geochemistry sample Data.

    :Example:
        
        >>> geoch=build_geochemistry_sample()
        >>> print(geoch)
        ... sampleData
        ... [['S0X4' '0' '254.0' 'PUP']
        ...     ['S0X4' '254' '521.0' 'mg']
        ...     ['S0X4' '521' '625.0' 'tut']
        ...     ['S0X4' '625' '984.0' 'suj']
        ...     ['S0X2' '0' '19.0' 'pup']
        ...     ['S0X2' '19' '425.0' 'hut']
        ...     ['S0X2' '425' '510.0' 'mgt']
        ...     ['S0X2' '510' '923.2' 'pyt']]
    """

    tempsamp,SampleList=[],[]
    comp=-1
    while 1:
        print('Enter Hole Name or <Enter/end> to stop:',end='')
        holeName=input()
        if holeName=="" or holeName.lower() in ["stop","end","enter",
                                                "finish","close"]:
            Sample=concat_array_from_list(list_of_array=SampleList,
                                          concat_axis=0)
            break
        comp=comp+1
        samp_buttom=np.float(input("Enter the buttom of the sampling (m):"))
        
        while holeName:
            if comp==-1:
                samP=concat_array_from_list(list_of_array=tempsamp,concat_axis=0)
                tempsamp=[]
                break
            print("Enter the sampling levels:",end='')
            samplevel=input()    
            samplevel=dump_comma(input_car=samplevel, max_value=12, 
                            carType='value')
            samp_ar=np.array(list(samplevel))
            samp_ar=samp_ar.reshape((samp_ar.shape[0],1))
            
            dh_to=samp_ar[1:]
            dh_to=np.append(dh_to,samp_buttom)
            dh_to=dh_to.reshape((dh_to.shape[0],1))
            
            print("Enter the samples' names:",end='')
            sampName=input()
            sampName=dump_comma(input_car=sampName, max_value=samp_ar.shape[0], 
                            carType='mixed')
            sampName_ar=np.array(list(sampName))
            sampName_ar=sampName_ar.reshape((sampName_ar.shape[0],1))
            try :
                holes_names=np.full((sampName_ar.shape[0],1),holeName)
                samfromto=np.concatenate((holes_names,samp_ar,dh_to,sampName_ar),axis=1)
                
            except Exception as e :
                raise ("IndexError!, arrrays sample_DH_From:{0} &DH_To :{1}&"
                       " 'Sample':{2} doesn't not match proprerrly.{3}".\
                           format(samp_ar.shape,dh_to.shape,sampName_ar.shape, e))
                warnings.warn("IndexError !, dimentional problem."
                              " please check np.ndarrays.shape.")
                _logger.warn("IndexError !, dimentional problem."
                              " please check np.ndarrays.shape.")
            tempsamp.append(samfromto)
            comp=-1
        SampleList.append(samP)
        print("\n")
        
    return Sample



def parse_wellData(filename=None, include_azimuth=False,
                   utm_zone="49N"):
    """
    Function to parse well information in*csv file 

    Parameters
    ----------
        * filename : str, optional
                full path to parser file, The default is None.
                
        * include_azimuth: bool , 
            Way to compute azimuth automatically 
            
        * utm_zone : str, 
            set coordinate _utm_WGS84. Defaut is 49N

    Raises
    ------
        FileNotFoundError
            if typical file deoesnt match the *csv file.

    Returns
    -------
       location:  str
            Name of location .
            
       WellData : np.ndarray
              Specificy the collar Data .
              
        GeoData : np.ndarray
             specify the geology data .
             
        SampleData : TYPE
            geochemistry sample Data.


    :Example: 
        
        >>> import numpy as np
        >>> dir_=r"F:\OneDrive\Python\CodesExercices\ex_avgfiles\modules"
        >>> parse_=parse_wellData(filename='Drill&GeologydataT.csv')
        >>> print("NameOflocation:\n",parse_[0])
        >>> print("WellData:\n",parse_[1])
        >>> print("GeoData:\n",parse_[2])
        >>> print("Sample:\n",parse_[3])
    """
    
    
    
    identity=["DH_Hole (ID)","DH_East","DH_North","DH_Dip",
              "Elevation" ,'DH_Azimuth',"DH_Top","DH_Bottom",
              "DH_PlanDepth","DH_Decr","Mask "]
    
    car=",".join([ss for ss in identity])
    # print(car)
    
    #ckeck the if it is the correct file
    _flag=0
    if filename is None :
        filename=[file for file in os.listdir(os.getcwd()) \
                  if (os.path.isfile(file)and file.endswith(".csv"))]
        # print(filename)
        if np.iterable (filename):
            for file in filename :
                with open (file,"r", encoding="utf-8") as f :
                    data=f.readlines()
                    _flag=1
        else :
            _logger.error('The {0} doesnt not match the *csv file'
                          ' You must convert file on *.csv format'.format(filename))
            warnings.error("The input file is wrong ! only *.csv file "
                          " can be parsed.")
    elif filename is not None :
        assert filename.endswith(".csv"), "The input file {0} is not in *.csv format"
        with open (filename,'r', encoding='utf-8') as f:
            data=f.readlines()   
            _flag=1
            
    if _flag==1 :
        try :
            # print(data[0])
            head=data[0].strip("'\ufeff").strip('\n').split(',')[:-1]
            head=_nonevalue_checker(list_of_value=head,
                                    value_to_delete='')
            chk=[1 for ii ,value in enumerate( head) if value ==identity[ii]]
            # data[0]=head
            if not all(chk):
                _logger.error('The {0} doesnt not match the correct file'\
                              ' to parse drill data'.format(filename))
                warnings.warn("The input file is wrong ! must "\
                              "input the correct file to be parsed.")
        except Exception as e :
            raise FileNotFoundError("The *csv file does no match the well file",e)
            
        # process to parse all data 
        # coll,geol,samp=[],[],[]

    for ss, elm in enumerate (data):
        elm=elm.split(',')
        for ii, value in enumerate(elm):
            if value=='' or value=='\n':
                elm=elm[:ii]
        data[ss]=elm
        
    data[0]=head

    [data.pop(jj)for jj, val in enumerate(data) if val==[]]

    if data[-1]==[]:
        data=data[:-1]
        
    ## final check ###
    safeData=_nonelist_checker(data=data, _checker=True ,
                  list_to_delete=['\n'])
    data=safeData[2]
    # identify collar , geology dans sample data 
    comp=0
    for ss , elm in enumerate(data):
        if elm[0].lower()=='geology':
            coll=data[:ss]
            comp=ss
        if elm[0].lower() =="sample":
            geol=data[comp+1:ss]
            samp=data[ss+1:]
            
    # build numpy data array 
    collar_list=coll[1:]
    # print(collar_list)
    collar_ar=np.zeros((len(collar_list), len(identity)),dtype='<U12')
    for ii, row in enumerate(collar_ar):
        collar_ar[ii:,:len(collar_list[ii])]= np.array(collar_list[ii])
    
    bottom_ar=collar_ar[:,7]
    # print(bottom_ar)
        
    geol_list=geol[1:]
    geol_ar=_order_well (data=geol_list,bottom_value=bottom_ar)
    samp_list=samp[1:]
    samp_ar=_order_well (data=samp_list,bottom_value=bottom_ar)
    
    name_of_location =filename[:-3]
    # find Description 
    DH_PlanDepth=collar_ar[:,8]
    DH_Decr=collar_ar[:,9]
    
    lenloc=len(name_of_location)
    
 
    for ss , singleArray in enumerate (DH_PlanDepth):
        # print(singleArray)
        if singleArray == ''or singleArray is None :
            singleArray=name_of_location[:-int(lenloc/2)]+collar_ar[ss,0][-2:]
            DH_PlanDepth[ss]=singleArray
            if DH_Decr[ss] ==''or DH_Decr[ss] is None :
                DH_Decr[ss] = "wdanx"+collar_ar[ss,0][0]+\
                    name_of_location.lower()[0]+collar_ar[ss,0][-1]
    
    collar_ar[:,8]=DH_PlanDepth
    collar_ar[:,9]=DH_Decr
    
    if include_azimuth==False:
        DH_Azimuth=np.full((collar_ar.shape[0]),0)
    elif include_azimuth== True:
        DH_Azimuth=compute_azimuth(easting = np.array([np.float(ii) 
                                                       for ii in collar_ar[:,1]]),
                                  northing = np.array([np.float(ii) 
                                                       for ii in collar_ar[:,2]]),
                                  utm_zone=utm_zone, extrapolate=True)
    
    collar_ar[:,5]=DH_Azimuth
    

    name_of_location=name_of_location.capitalize()
    WellData,GeoData,SampleData=collar_ar,geol_ar, samp_ar 

    return (name_of_location, WellData,GeoData,SampleData)
                   

            
def _nonelist_checker(data, _checker=False ,
                      list_to_delete=['\n']):
    """
    Function to delete a special item on list in data.
    Any item you want to delete is acceptable as long as item is on a list.
    
    Parameters
    ----------
        * data : list
             container of list. Data must contain others list.
             the element to delete should be on list.
             
        *  _checker : bool, optional
              The default is False.
             
        * list_to_delete : TYPE, optional
                The default is ['\n'].

    Returns
    -------
        _occ : int
            number of occurence.
        num_turn : int
            number of turns to elimate the value.
        data : list
           data safeted exempt of the value we want to delete.
        
    :Example: 
        
        >>> import numpy as np 
        >>> listtest =[['DH_Hole', 'Thick01', 'Thick02', 'Thick03',
        ...   'Thick04', 'sample02', 'sample03'], 
        ...   ['sample04'], [], ['\n'], 
        ...  ['S01', '98.62776918', '204.7500461', '420.0266651'], ['prt'],
        ...  ['pup', 'pzs'],[], ['papate04', '\n'], 
        ...  ['S02', '174.4293956'], [], ['400.12', '974.8945704'],
        ...  ['pup', 'prt', 'pup', 'pzs', '', '\n'],
        ...  ['saple07'], [], '',  ['sample04'], ['\n'],
        ...  [''], [313.9043882], [''], [''], ['2'], [''], ['2'], [''], [''], ['\n'], 
        ...  [''], ['968.82'], [], [],[], [''],[ 0.36], [''], ['\n']]
        >>> ss=_nonelist_checker(data=listtest, _checker=True,
        ...                        list_to_delete=[''])
        >>> print(ss)
    """
    
    _occ,num_turn=0,0
    
    
    if _checker is False:
        _occ=0
        return _occ, num_turn, data 
    
   
    while _checker is True :
        if list_to_delete in data :
            for indix , elem_chker in enumerate(data):
                if list_to_delete == elem_chker :
                    _occ=_occ+1
                    del data[indix]
                    if data[-1]==list_to_delete:
                        _occ +=1
                        # data=data[:-1]
                        data.pop()
        elif list_to_delete not in data :
            _checker =False
        num_turn+=1

        
    return _occ,num_turn,data
        
def _order_well (data,**kwargs):
    """
    Function to reorganize value , depth rock and depth-sample 
    the program controls  the input depth value and compare it . 
    with the bottom. It will pay attention that bottom depth must
    be greater or egual of any other value. In the same time ,
    the program check if value entered are sorted on ascending order . 
    well must go deep ( less value to great value). Negative values of
    depths are not acceptable.
    
    Parameters
    ----------
        * data : list,
                data contains list of well thickness and rock description .
    
        * bottom_value  : np.ndarray float
                value of bottom . it may the basement extrapolation.
                default is 1.023 km

    Returns
    -------
         np.ndarray
            data, data aranged to [DH_Hole, DH_From, DH_To, Rock.] arrays
        
    :Example:
        
        >>> import numpy as np
        >>> listtest =[['DH_Hole', 'Thick01', 'Thick02', 'Thick03',
        ...           'Thick04','Rock01', 'Rock02', 'Rock03', 'Rock04'],
        >>> ['S01', '0.0', '98.62776918', '204.7500461','420.0266651',
        ...     'GRT', 'ATRK', 'GRT', 'ROCK'],
        >>> ['S02', '174.4293956', '313.9043882', '400.12', '974.8945704',
        ...     'GRT', 'ATRK', 'GRT', 'ROCK']]
        >>> print(listtest[1:])
        >>> ss=_order_well(listtest[1:])
        >>> print(ss)
    """
    
    this_function_name=inspect.getframeinfo(inspect.currentframe())[2]

    bottomgeo=kwargs.pop("bottom_value", None)
    # bottomsamp=kwargs.pop("sample_bottom_value",None)
    # if type(bottomgeo) is np.ndarray :_flag=1
    # else : _flag=0
    
    temp=[]

    _logger.info ("You pass by {0} function! "
                  "Thin now , everything is ok. ".format(this_function_name))
    
    for jj, value in enumerate (data):
        thickness_len,dial1,dial2=intell_index(datalist=value[1:]) #call intell_index
        value=np.array(value)

        dh_from=value[1:thickness_len+1]
        dh_to =np.zeros((thickness_len,), dtype="<U12")
                # ---- check the last value given        
        max_given_bottom=np.max(np.array([np.float(ii) for ii in dh_from])) 
                # ----it may be less than bottom value                                                                .
        dh_geo=value[thickness_len+1:]
        dh_to=dh_from[1:]
        # if _flag==0 :
        if bottomgeo[jj] =="" or bottomgeo[jj]==None :
            bottomgeoidx=1023.
            if max_given_bottom > bottomgeoidx:
                _logger.warn("value {0} is greater than the "
                             "Bottom depth {1}".format(max_given_bottom ,
                                                       bottomgeoidx))
                warnings.warn (
                    "Given values have a value greater than the depth !")
                
            dh_to=np.append(dh_to,bottomgeoidx)
        else: #elif :_flag==1 :         # check the bottom , any values given 
                                        # must be less or egual to depth not greater.
            if max_given_bottom> np.float(bottomgeo[jj]):
                _logger.warn("value {0} is greater than the Bottom "
                             "depth {1}".format(max_given_bottom ,
                                                bottomgeo[jj]))
                warnings.warn (
                    "Given values have a value greater than the depth !")
            dh_to=np.append(dh_to,bottomgeo[jj])
            
        dh_hole=np.full((thickness_len),value[0])
        
        temp.append(np.concatenate((dh_hole.reshape((dh_hole.shape[0],1)),
                                    dh_from.reshape((dh_from.shape[0],1)),
                                    dh_to.reshape((dh_to.shape[0],1)),
                                    dh_geo.reshape((dh_geo.shape[0],1))),axis=1 ))
        
    data=concat_array_from_list(list_of_array=temp, concat_axis=0)
    
    return data
                                
        
def intell_index (datalist,assembly_dials =False):
    """
    function to search index to differency value to string element like 
    geological rocks and geologicals samples. It check that value are sorted
    in ascending order.

    Parameters
    ----------
        * datalist : list
                list of element : may contain value and rocks or sample .
        * assembly_dials : list, optional
                separate on two list : values and rocks or samples. The default is False.

    Returns
    -------
        index: int
             index of breaking up.
        first_dial: list , 
           first sclice of value part 
        secund_dial: list , 
          second slice of rocks or sample part.
        assembly : list 
             list of first_dial and second_dial
    
    :Example: 
        
        >>> import numpy as np
        >>> listtest =[['DH_Hole', 'Thick01', 'Thick02', 'Thick03',
        ...           'Thick04','Rock01', 'Rock02', 'Rock03', 'Rock04'],
        ...           ['S01', '0.0', '98.62776918', '204.7500461','420.0266651','520', 'GRT', 
        ...            'ATRK', 'GRT', 'ROCK','GRANODIORITE'],
        ...           ['S02', '174.4293956', '313.9043882','974.8945704', 'GRT', 'ATRK', 'GRT']]
        >>> listtest2=[listtest[1][1:],listtest[2][1:]]
        >>> for ii in listtest2 :
        >>> op=intell_index(datalist=ii)
        >>> print("index:\n",op [0])
        >>> print('firstDials :\n',op [1])
        >>> print('secondDials:\n',op [2])
    """
    # assembly_dials=[]
    max_=0              # way to check whether values are in sort (ascending =True) order 
                        # because we go to deep (first value must be less than the next none)
    for ii, value in enumerate (datalist):
        try : 
            thick=float(value)
            if thick >= max_:
                max_=thick 
            else :
                _logger.warning("the input value {0} is less than the previous one."\
                              " Please enter value greater than {1}.".format(thick, max_) )
                warnings.warn("Value {1} must be greater than the previous value {0}."\
                              " Must change on your input data.".format(thick,max_))
        except :
            # pass
            indexi=ii
            break
        
    first_dial=datalist[:indexi]
    second_dial =datalist[indexi:]

    if assembly_dials:
        
        assembly_dials=[first_dial,second_dial]
        
        return indexi, assembly_dials
        
    return indexi, first_dial, second_dial

def _nonevalue_checker (list_of_value, value_to_delete=None):
    """
    Function similar to _nonelist_checker. the function deletes the specific 
    value on the list whatever the number of repetition of value_to_delete.
    The difference with none list checker is value to delete 
    may not necessary be on list.
    
    Parameters
    ----------
        * list_of_value : list
                list to check.
        * value_to_delete : TYPE, optional
            specific value to delete. The default is ''.

    Returns
    -------
        list
            list_of_value , safe list without the deleted value .
    
    :Example: 
        
        >>> import numpy as np 
        >>> test=['DH_Hole (ID)', 'DH_East', 'DH_North',
        ...          'DH_Dip', 'Elevation ', 'DH_Azimuth', 
        ...            'DH_Top', 'DH_Bottom', 'DH_PlanDepth', 'DH_Decr', 
        ...            'Mask', '', '', '', '']
        >>> test0= _nonevalue_checker (list_of_value=test)
        >>> print(test0)
    """
    if value_to_delete is None :
        value_to_delete  =''
    
    if type (list_of_value) is not list :
        list_of_value=list(list_of_value)
        
    start_point=1
    while start_point ==1 :
        if value_to_delete in list_of_value :
            [list_of_value.pop(ii) for  ii, elm in\
             enumerate (list_of_value) if elm==value_to_delete]
        elif value_to_delete not in list_of_value :
            start_point=0 # not necessary , just for secure the loop. 
            break           # be sure one case or onother , it will break
    return list_of_value 

###OPTIMIZE  
def resize_resphase_values (dictobj: dict , fill_value: float = np.nan,
                            c =None,mask_nan:bool =True, return_array =True): 
    """ Get the resistivity and phases values from dictof items 
    and concatenate them. 
    
    Array must have the same length everywhere. Function get the max 
    length and add `np.nan` values to reach the size of max array. 
    if the given data is under the max array. 
    
    :param dictobj: dict - container of array-like values 
    :param fill_value: float - value to append to the max array so 
        to reach the reference array. Default is ``np.nan`.` 
    :param c: array-like - Controler; the array value for reference. If given,  
        its length should be used as the max array length. Default is ``None``.
    :param mask_value: bool - Mask the appender value. Default is ``True``. 
    :param return_array: bool - Return concatenated array otherwise return the 
        dictobject with filled values. 
        
    :return: ndarray - array of `dictobj` values concatenated into two 
        dimensional arrays; otherwise return the `dictobj` resized
        
    """
    max_ = len(c) if c is not None else  max (
        [len(v) for v in dictobj.values()]) 
    for key,  aval in sorted (dictobj.items()) : 
        s = max_-len(aval)
        if s==0: 
            continue 
        elif s > 0: 
            aval = np.append(aval, np.full(
                (s,), fill_value ))
        elif s < 0 : 
            aval = aval [:s]
         
        dictobj[key] = aval
        
    # if not return_array : 
    #     return dictobj 
    
    ar = np.vstack([v for v in dictobj.values()]) 
    
    # masked the fill values
    return dictobj if not return_array else np.ma.masked_array( ar, mask = np.isnan(ar) 
                              ) if mask_nan else ar 
   

def _strip_item(item_to_clean, item=None, multi_space=12):
    """
    Function to strip item around string values.  if the item to clean is None or 
    item-to clean is "''", function will return None value

    Parameters
    ----------
        * item_to_clean : list or np.ndarray of string 
                 List to strip item.
        * cleaner : str , optional
                item to clean , it may change according the use. The default is ''.
        * multi_space : int, optional
                degree of repetition may find around the item. The default is 12.

    Returns
    -------
        list or ndarray
            item_to_clean , cleaned item 
            
    :Example: 
        
     >>> import numpy as np
     >>> new_data=_strip_item (item_to_clean=np.array(['      ss_data','    pati   ']))
     >>>  print(np.array(['      ss_data','    pati   ']))
     ... print(new_data)

    """
    if item==None :item = ' '
    
    cleaner =[(''+ ii*'{0}'.format(item)) for ii in range(multi_space)]
    
    if type(item_to_clean ) != list :#or type(item_to_clean ) !=np.ndarray:
        if type(item_to_clean ) !=np.ndarray:
            item_to_clean=[item_to_clean]
    ###TIP
    if item_to_clean in cleaner or item_to_clean ==['']:
        warnings.warn ('No data found for sanitization; returns None.')
        return None 

  
    try : 
        multi_space=int(multi_space)
    except : 
        raise TypeError('argument <multplier> must be'\
                        ' an integer not {0}'.format(type(multi_space)))
    
    for jj, ss in enumerate(item_to_clean) : 
        for space in cleaner:
            if space in ss :
                new_ss=ss.strip(space)
                item_to_clean[jj]=new_ss
    
    return item_to_clean
    
def _cross_eraser (data , to_del, deep_cleaner =False):
    """
    Function to delete some item present in another list. It may cheCk deeper 

    Parameters
    ----------
        * data : list
                Main data user want to filter.
        * to_del : list
                list of item you want to delete present on the main data.
        * deep_cleaner : bool, optional
                Way to deeply check. Sometimes the values are uncleaned and 
            capitalizeed . this way must not find their safety correspondace 
            then the algorth must clean item and to match all at 
            the same time before erasing. The *default* is False.

    Returns
    -------
        list
         data , list erased.

    :Example: 
        
        >>> data =['Z.mwgt','Z.pwgt','Freq',' Tx.Amp','E.mag','   E.phz',
        ...          '   B.mag','   B.phz','   Z.mag', '   Zphz  ']
        >>> data2=['   B.phz','   Z.mag',]
        ...     remain_data =cross_eraser(data=data, to_del=data2, 
        ...                              deep_cleaner=True)
        >>> print(remain_data)
    """
    
    data , to_del=_strip_item(item_to_clean=data), _strip_item(item_to_clean=to_del)
    if deep_cleaner : 
        data, to_del =[ii.lower() for ii in data], [jj.lower() for jj in to_del]
        
    for index, item in enumerate(data): 
        while item in to_del :
            del data[index]
            if index==len(data)-1 :
                break

    return data 

def _remove_str_word (char, word_to_remove, deep_remove=False):
    """
    Small funnction to remove a word present on  astring character 
    whatever the number of times it will repeated.
    
    Parameters
    ----------
        * char : str
                may the the str phrases or sentences . main items.
        * word_to_remove : str
                specific word to remove.
        * deep_remove : bool, optional
                use the lower case to remove the word even the word is uppercased 
                of capitalized. The default is False.

    Returns
    -------
        str 
            char , new_char without the removed word .
        
    :Example: 
        
        >>> from pycsamt.utils  import func_utils as func
        >>> ch ='AMTAVG 7.76: "K1.fld", Dated 99-01-01,AMTAVG, Processed 11 Jul 17 AMTAVG'
        >>> ss=func._remove_str_word(char=ch, word_to_remove='AMTAVG', deep_remove=False)
        >>> print(ss)
    """
    if type(char) is not str : char =str(char)
    if type(word_to_remove) is not str : word_to_remove=str(word_to_remove)
    
    if deep_remove == True :
        word_to_remove, char =word_to_remove.lower(),char.lower()

    if word_to_remove not in char :
        return char

    while word_to_remove in char : 
        if word_to_remove not in char : 
            break 
        index_wr = char.find(word_to_remove)
        remain_len=index_wr+len(word_to_remove)
        char=char[:index_wr]+char[remain_len:]

    return char

def stn_check_split_type(data_lines): 
    """
    Read data_line and check for data line the presence of 
    split_type < ',' or ' ', or any other marks.>
    Threshold is assume to be third of total data length.
    
    :params data_lines: list of data to parse . 
    :type data_lines: list 
 
    :returns: The split _type
    :rtype: str
    
    :Example: 
        
        >>> from pycsamt.utils  import func_utils as func
        >>> path =  os.path.join(os.environ["pyCSAMT"], 
                          'csamtpy','data', K6.stn)
        >>> with open (path, 'r', encoding='utf8') as f : 
        ...                     data= f.readlines()
        >>>  print(func.stn_check_split_type(data_lines=data))
            
    """

    split_type =[',', ':',' ',';' ]
    data_to_read =[]
    if isinstance(data_lines, np.ndarray): # change the data if data is not dtype string elements.
        if data_lines.dtype in ['float', 'int', 'complex']: 
            data_lines=data_lines.astype('<U12')
        data_lines= data_lines.tolist()
        
    if isinstance(data_lines, list):
        for ii, item in enumerate(data_lines[:int(len(data_lines)/3)]):
             data_to_read.append(item)
             data_to_read=[''.join([str(item) for item in data_to_read])] # be sure the list is str item . 

    elif isinstance(data_lines, str): data_to_read=[str(data_lines)]
    
    for jj, sep  in enumerate(split_type) :
        if data_to_read[0].find(sep) > 0 :
            if data_to_read[0].count(sep) >= 2 * len(data_lines)/3:
                if sep == ' ': return  None  # use None more conventional 
                else : return sep 

def minimum_parser_to_write_edi (edilines, parser = '='):
    """
    This fonction validates edifile for writing , string with egal.
    We assume that dictionnary in list will be for definemeasurment 
    E and H fied. 
    
    :param edilines: list of item to parse 
    :type edilines: list 
    
    :param parser: the egal is use  to parser edifile .
                    can be changed, default is `=`
    :type parser: str 
  
    """
    if isinstance(edilines,list):
        if isinstance(edilines , tuple) : edilines =list(edilines)
        else :raise TypeError('<Edilines> Must be on list')
    for ii, lines in enumerate(edilines) :
        if isinstance(lines, dict):continue 
        elif lines.find('=') <0 : 
            raise 'None <"="> found on this item<{0}> of  the edilines list. list can not'\
            ' be parsed.Please put egal between key and value '.format(edilines[ii])
    
    return edilines 
            

def round_dipole_length(value, round_value =5.): 
    """ 
    small function to graduate dipole length 5 to 5. Goes to be reality and 
    simple computation .
    
    :param value: value of dipole length 
    :type value: float 
    
    :returns: value of dipole length rounded 5 to 5 
    :rtype: float
    """ 
    mm = value % round_value 
    if mm < 3 :return np.around(value - mm)
    elif mm >= 3 and mm < 7 :return np.around(value -mm +round_value) 
    else:return np.around(value - mm +10.)
    
def keepmin(array): 
    """ Keep the minimum in array and array value"""
    os= np.where(array==array.min())[0]
    if len(os)>1 :
        os= os[0]
    return int(os), array[int(os)]


def get_closest_value (values_range, input_value): 
    """
    Function  to select the closest values when input values is
    not in the values range. We assume that value types are arrays.
    If the same value is repeated, should take the first index and 
    the value at that index. 
    
    :param values_range: values to get 
    :type values_range: array_like   
                 
    :param input_value: specific value
    :type input_value: float,

    :returns: the closest value and its index
    :rtye: float 
    
    """

    values_range = np.array(values_range)
    # if all element less than zero than convert 
     #in put value to negative value 
    if np.all(values_range <0) : 
        if input_value >0 : 
            input_value *=-1  # for depth select purpose 
        
    if input_value < values_range.min(): 
        print('--> ! Input value ={0} is out  the range,'
              ' min value = {1}. Value should be reset to ={1}.'.
              format(input_value, values_range.min())
              )
        
        warnings.warn('Input value ={0} is out  '
                      'the range ! min value = {1}'.
                      format(input_value, values_range.min())
                      )
        _logger.debug('Input value ={0} is out '
                      'the range ! min value = {1}'.
                      format(input_value, values_range.min())
                      )  
            
        input_value = values_range.min()
    elif input_value > values_range.max(): 
        
        warnings.warn('Input value ={0} is out '
                      'the range ! max value = {1}'.
                      format(input_value, values_range.max())
                      )
        _logger.debug('Input value ={0} is out '
                      'the range ! max value = {1}'.
                      format(input_value, values_range.max()))
            
        input_value = values_range.max()
        print('!--> Input value ={0} is out '
              'the range , min value = {1}. Value'
              ' should be reset to ={1}.'.format(
                  input_value, values_range.max()))
        
    if input_value in values_range : 
        indexes,*_= np.where(values_range==input_value)
        if len(indexes)==1 : 
            index =int(indexes )
        # mean element is repeated then take the first index    
        elif len(indexes)>1 : 
            index = int(indexes)[0]
            
        value = values_range[index]
        return value , index 
    
    elif values_range.min() < input_value < values_range.max(): 
        # values_range = sorted(values_range)
        for ii, xx in enumerate(values_range): 
            if xx > input_value : 
                #compute distance : 
                # make the diffence between distances 
                d0 = abs(input_value-xx) 
                # need to take the short
                d1= abs(values_range[ii-1]-input_value) 
                if d0 < d1 : 
                    return  xx , ii
                elif d0> d1 : 
                    return   values_range[ii-1], ii-1
                
def geo_length_checker(main_param, 
                       optional_param,
                       force =False, 
                        param_names =('input_resistivities',
                                      'input_layers'),
                        **kws): 
    """
    Geo checker is a function to check the differents length 
    of different geoparams.
    
    The length of optional params  should depend of the length of main params. 
    Therefore if the length of optional params is larger than the length of 
    the main params, the length of optional params will be reduced to
    the length of main params.Otherwise if the length of optional params 
    is shorther than the length of the main params, will filled it either 
    with "None" if dtype param is string or 0. is float or 0 if integer.
    If `force`  is set ``True``, shoud raise errors if the main params and 
    the optional params have are not the same length. 
  
    Parameters 
    ------------
        * main_param : array_like, list 
                 main parameter that must took 
                 its length as reference length 
                 
        * optional params : array_like, list 
                 optional params, whom length depend 
                 to the length of main params
                 
        * param_names : tuple or str 
                 names of main params and optional params 
                 so to generate error if exits.
                 
        * fill_value: str, float, optional  
                Default value to fill thearray in the case where 
                the length of optional param is 
                less than the length of the  main param .If None ,
                will fill according to array dtype
            
    Returns 
    --------
        array_like 
           optional param truncated according to the man params 
    """
    add_v =kws.pop('fill_value', None)
    

    if isinstance(main_param, (str, float, int)):
        main_param = np.array([main_param])
    if isinstance(optional_param, (str, float, int)):
       optional_param = np.array([optional_param])
    if isinstance(main_param, (list, tuple)):
        main_param =np.array(main_param)
    if isinstance(optional_param, (list, tuple)):
        optional_param =np.array(optional_param)
            
    mes=''
    if len(optional_param) > len(main_param): 
        mes ="".join(["---> Note ! {0} will be truncated ",
            "to length = {1}as the same length of {2} ."])
        
        warnings.warn(mes.format(param_names[1],
                                 len(main_param),param_names[0] ))
        
        optional_param= optional_param[:len(main_param)]
        if force is True : 
            mess =''.join(['--> `force` argument is set ``True``,', 
                           ' Can not truncate {0} = {1} to fit the ',
                           'length of {2} = {3}.'])

            raise CSex.pyCSAMTError_parameter_number(
                mess.format(param_names[1], 
                        len(param_names[1]), 
                        param_names[0],
                        len(param_names[0])))

    elif len(optional_param) < len(main_param) : 
        if force is True : 
            mess =''.join([ '--> `force` argument  is set ``True``,',
                           ' Can not fill the value of {0} ',
                            'to match the length of {1} = {2}.'])
            
            raise CSex.pyCSAMTError_parameter_number(
                mess.format(param_names[1], param_names[0],
                            len(param_names[0])))
        if add_v is not None : 
            # repeat the value to add 
            add_v =[add_v for vv in range(
                len(main_param)-len(optional_param))]
            add_v = np.array(add_v)
        if add_v is None :
             if optional_param.dtype not in [ 'float', 'int'] : 
                 add_v =['None' for i in range(
                     len(main_param)-len(optional_param))]
    
             else : 
                 for type_param, fill_value in zip(
                         [ 'float', 'int'],[ 0., 0] ): 
                     if  type_param  == optional_param.dtype :
                         add_v =[fill_value for i in range(
                             len(main_param)-len(optional_param))]

        mes =''.join(["--> Length of {0} is ={1} which ",
                      "length of {2} is ={3}. We'll add {4}", 
                      "to fill {5} value."
                      ])

        warnings.warn(mes.format(param_names[1],
                                 len(optional_param),
                                 param_names[0],
               len(main_param), add_v[0], param_names[1]))
        optional_param= optional_param.tolist()
        optional_param.extend(add_v)

    return np.array(optional_param)

def fr_en_parser (f, delimiter =':'): 
    """ Parse the translated data file. 
    
    :param f: translation file to parse 
    :param delimiter: str, delimiter
    
    :return: generator obj, composed of a list of 
        french  and english Input translation. 
    
    :Example:
        >>> file_to_parse = 'pme.parserf.md'
        >>> path_pme_data = r'C:/Users\Administrator\Desktop\__elodata
        >>> data =list(BS.fr_en_parser(
            os.path.join(path_pme_data, file_to_parse)))
    """
    
    is_file = os.path.isfile (f)
    if not is_file: 
        raise IOError(f'Input {f} is not a file. Please check your file.')
    
    with open(f, 'r', encoding ='utf8') as ft: 
        data = ft.readlines()
        for row in data :
            if row in ( '\n', ' '):
                continue 
            fr, en = row.strip().split(delimiter)
            yield([fr, en])

def convert_csvdata_from_fr_to_en(csv_fn, pf, destfile = 'pme.en.csv',
                                  savepath =None, delimiter =':'): 
    """ Translate variable data from french csva data  to english with 
    varibale parser file. 
    
    :param csv_fn: data collected in csv format 
    :param pf: parser file 
    :param destfile: str,  Destination file, outputfile 
    :param savepath: [Path-Like object, save data to a path 
                      
    :Example: 
        # to execute this script, we need to import the two modules below
        >>> import os 
        >>> import csv 
        >>> path_pme_data = r'C:/Users\Administrator\Desktop\__elodata
        >>> datalist=convert_csvdata_from_fr_to_en(
            os.path.join( path_pme_data, _enuv2.csv') , 
            os.path.join(path_pme_data, pme.parserf.md')
                         savefile = 'pme.en.cv')
    """
    # read the parser file and separed english from french 
    parser_data = list(fr_en_parser (pf,delimiter) )
    
    with open (csv_fn, 'r', encoding ='utf8') as csv_f : 
        csv_reader = csv.reader(csv_f) 
        csv_data =[ row for row in csv_reader]
    # get the index of the last substring row 
    ix = csv_data [0].index ('Industry_type') 
    # separateblock from two 
    csv_1b = [row [:ix +1] for row in csv_data] 
    csv_2b =[row [ix+1:] for row in csv_data ]
    # make a copy of csv_1b
    csv_1bb= deepcopy(csv_1b)
   
    for ii, rowline in enumerate( csv_1bb[3:]) : # skip the first two rows 
        for jj , row in enumerate(rowline): 
            for (fr_v, en_v) in  parser_data: 
                # remove the space from french parser part
                # this could reduce the mistyping error 
                fr_v= fr_v.replace(
                    ' ', '').replace('(', '').replace(
                        ')', '').replace('\\', '').lower()
                 # go  for reading the half of the sentence
                row = row.lower().replace(
                    ' ', '').replace('(', '').replace(
                        ')', '').replace('\\', '')
                if row.find(fr_v[: int(len(fr_v)/2)]) >=0: 
                    csv_1bb[3:][ii][jj] = en_v 
    
    # once translation is done, concatenate list 
    new_csv_list = [r1 + r2 for r1, r2 in zip(csv_1bb,csv_2b )]
    # now write the new scv file 
    if destfile is None: 
        destfile = f'{os.path.basename(csv_fn)}_to.en'
        
    destfile.replace('.csv', '')
    
    with open(f'{destfile}.csv', 'w', newline ='',encoding ='utf8') as csvf: 
        csv_writer = csv.writer(csvf, delimiter=',')
        csv_writer.writerows(new_csv_list)
        # for row in  new_csv_list: 
        #     csv_writer.writerow(row)
    savepath = cpath(savepath , '__pme')
    try :
        shutil.move (f'{destfile}.csv', savepath)
    except:pass 
    
    return new_csv_list

# if __name__=="__main__" :

    # parse_=parse_wellData(filename='shimenDH.csv',
    # include_azimuth=True,utm_zone='49N')
    
    # print("NameOflocation:\n",parse_[0])
    # print("WellData:\n",parse_[1])
    # print("GeoData:\n",parse_[2])
    # print("Sample:\n",parse_[3])
    # data =['Z.mwgt','Z.pwgt','Freq',' Tx.Amp','E.mag','   E.phz',
    #         '   B.mag','   B.phz','   Z.mag', '   Zphz  ']
    
    # data2=[' B.phz',' Z.mag',]
    # # cleaner =[(''+ ii*'*') for ii in range(7)]
    # print(_cross_eraser(data=data2, to_del=data))
    # print(_cross_eraser(data=data, to_del=data2))
    # # print(_strip_item(item_to_clean=data, item=' '))
    # # # print(cleaner)
    # ts ='50.0'
    # ch ='AMTAVG 7.76: "K1.fld", Dated 99-01-01,AMTAVG, Processed 11 Jul 17 AMTAVG'
    # ss=_remove_str_word(char=ts, word_to_remove='m', deep_remove=False)
            


    

    
    
    
    
    
    
    

    
        