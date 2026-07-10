#   License: BSD-3-Clause
#   Author: LKouadio <etanoyau@gmail.com>

from __future__ import annotations

import copy
import csv
import datetime
import inspect
import itertools
import json
import numbers
import os
import pickle
import re
import shutil
import subprocess
import sys
import warnings
from zipfile import ZipFile

import h5py
import joblib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import yaml
from six.moves import urllib

from .._typing import (
    Any,
    ArrayLike,
    DataFrame,
    Dict,
    F,
    Iterable,
    List,
    NDArray,
    Optional,
    Sub,
    T,
    Text,
    Tuple,
)
from .._watexlog import watexlog
from ..exceptions import (
    EDIError,
    ParameterNumberError,
)
from ..property import P
from ._dependency import import_optional_dependency

_logger = watexlog.get_watex_logger(__name__)

_msg= ''.join([
    'Note: need scipy version 0.14.0 or higher or interpolation,',
    ' might not work.']
)
_msg0 = ''.join([
    'Could not find scipy.interpolate, cannot use method interpolate'
     'check installation you can get scipy from scipy.org.']
)

try:
    scipy_version = [int(ss) for ss in scipy.__version__.split('.')]
    if scipy_version[0] == 0:
        if scipy_version[1] < 14:
            warnings.warn(_msg, ImportWarning, stacklevel=2)
            _logger.warning(_msg)

    import scipy.interpolate as spi

    interp_import = True
 # pragma: no cover
except ImportError:

    warnings.warn(_msg0, stacklevel=2)
    _logger.warning(_msg0)

    interp_import = False


def to_numeric_dtypes (
    arr: NDArray|DataFrame, *,
    columns: List[str, ...]=None,
    return_feature_types: bool=... ,
    missing_values: float=np.nan,
    pop_cat_features: bool=...,
    sanitize_columns: bool=...,
    regex: re|str=None,
    fill_pattern: str='_',
    drop_nan_columns: bool=True,
    how: str='all',
    reset_index: bool=...,
    drop_index: bool=True,
    verbose: bool=...,
    )-> DataFrame :
    r""" Convert array to dataframe and coerce arguments to appropriate dtypes.

    Function includes additional tools to manipulate the transformed data
    such as:

    - `pop_cat_features` to remove the categorical attributes,
    - `sanitize_columns` to clean the columns of the dataframe by removing
      the undesirable characters,
    - `drop_nan_columns` to drop all the columns and/or rows that contains
      full NaN, ...

    Parameters
    -----------
    arr: Ndarray or Dataframe, shape (m_samples, n_features)
        Array of dataframe to create, to sanitize or to auto-detect
        feature categories ( numerical or categorical).

    columns: list of str, optional
        Usefull to create a dataframe when array is given. Be aware to fit the
        number of array columns (shape[1])

    return_feature_types: bool, default=False,
        return the list of numerical and categorial features.

    missing_values: float, default='NaN'
        Replace the missing or empty string if exist in the dataframe.

    pop_cat_features:bool, default=False,
        remove the categorial features  from the DataFrame.

    sanitize_columns: bool, default=False,
       remove undesirable character in the data columns using the default
       argument of `regex` parameters.

       .. versionadded:: 0.1.9

    regex: `re` object,
        Regular expresion object used to polish the data columns.
        the default is::

        >>> import re
        >>> re.compile (r'[_#&.)(*@!_,;\s-]\s*', flags=re.IGNORECASE)

       .. versionadded:: 0.1.9

    fill_pattern: str, default=''
        Pattern to replace the non-alphabetic character in each item of
        columns.

    drop_nan_columns: bool, default=True
       Remove all columns filled by NaN values.

       .. versionadded: 0.2.4
          By default, it auto-removes columns with all NaN values. To
          deactive this functionality, set it to ``False``.

    how: str, default='all'
       Drop also the NaN row data. The row data which is composed entirely
       with NaN or Null values.

    reset_index: bool, default=False
       Reset the index of the dataframe.

       .. versionadded: 0.2.4

    drop_index: bool, default=True,
       Drop index in the dataframe after reseting.

       .. versionadded: 0.2.4

    verbose: bool, default=False,
        outputs a message by listing the categorial items dropped from
        the dataframe if exists.

    Returns
    --------
    df or (df, nf, cf): Dataframe of values casted to numeric types
        also return `nf` and `cf`  if `return_feature_types` is set
        to``True``.

    Examples
    ---------
    >>> from watex.datasets.dload import load_bagoue
    >>> from watex.utils.funcutils import to_numeric_dtypes
    >>> X, y = load_bagoue (as_frame =True )
    >>> X0 =X[['shape', 'power', 'magnitude']]
    >>> X0.dtypes
    ... shape        object
        power        object
        magnitude    object
        dtype: object
    >>> df = to_numeric_dtypes(X0)
    >>> df.dtypes
    ... shape         object
        power        float64
        magnitude    float64
        dtype: object

    """
    from .validator import _is_numeric_dtype

    # pass ellipsis argument to False
    ( sanitize_columns, reset_index, verbose, return_feature_types,
     pop_cat_features ) = ellipsis2false(
        sanitize_columns, reset_index, verbose, return_feature_types,
        pop_cat_features
        )

    if not is_iterable (arr, exclude_string=True):
        raise TypeError(f"Expect array. Got {type (arr).__name__!r}")

    if hasattr ( arr, '__array__') and hasattr ( arr, 'columns'):
        df = arr.copy()
        if columns is not None:
            if verbose:
                print("Dataframe is passed. Columns should be replaced.")
            df =pd.DataFrame ( np.array ( arr), columns =columns )

    else: df = pd.DataFrame (arr, columns =columns  )

    # sanitize columns
    if sanitize_columns:
        # Pass in the case columns are all integer values.
        if not _is_numeric_dtype(df.columns , to_array=True):
           # for consistency reconvert to str
           df.columns = df.columns.astype(str)
           df = sanitize_frame_cols(
               df, regex=regex, fill_pattern=fill_pattern )

    #replace empty string by Nan if NaN exist in dataframe
    df= df.replace(r'^\s*$', missing_values, regex=True)

    # check the possibililty to cast all
    # the numerical data
    for serie in df.columns:
        try:
            df= df.astype(
                {serie:np.float64})
        except:continue

    # drop nan  columns if exists
    if drop_nan_columns:
        if verbose:
            nan_columns = df.columns [ df.isna().all()].tolist()
            print("No NaN column found.") if len(
                nan_columns)==0 else listing_items_format (nan_columns,
                    "NaN columns found in the data",
                    " ", inline =True, lstyle='.')
        # drop rows and columns with NaN values everywhere.
        df.dropna ( axis=1, how='all', inplace =True)
        if str(how).lower()=='all':
            df.dropna ( axis=0, how='all', inplace =True)

    # reset_index of the dataframe
    # This is useful after droping rows
    if reset_index:
        df.reset_index (inplace =True, drop = drop_index )
    # collect numeric and non-numeric data
    nf, cf =[], []
    for serie in df.columns:
        if _is_numeric_dtype(df[serie], to_array =True ):
            nf.append(serie)
        else: cf.append(serie)

    if pop_cat_features:
        [ df.pop(item) for item in cf ]
        if verbose:
            msg ="Dataframe does not contain any categorial features."
            b= f"Feature{'s' if len(cf)>1 else ''}"
            e = (f"{'have' if len(cf) >1 else 'has'} been dropped"
                 " from the dataframe.")
            print(msg) if len(cf)==0 else listing_items_format (
                cf , b, e ,lstyle ='.', inline=True)

        return df

    return (df, nf, cf) if return_feature_types else df


def listing_items_format (
        lst, /, begintext ='', endtext='' , bullet='-',
        enum =True , lstyle=None , space =3 , inline =False, verbose=True
        ):
    """ Format list by enumerate them successively with carriage return

    :param lst: list,
        object for listening
    :param begintext: str,
        Text to display at the beginning of listing the items in `lst`.
    :param endtext: str,
        Text to display at the end of the listing items in `lst`.
    :param enum:bool, default=True,
        Count the number of items in `lst` and display it
    :param lstyle: str, default =None
        listing marker.
    :param bullet:str, default='-'
        symbol that is used to introduce item if `enum` is set to False.
    :param space: int,
        number of space to keep before each outputted item in `lst`
    :param inline: bool, default=False,
        Display all element inline rather than carriage return every times.
    :param verbose: bool,
        Always True for print. If set to False, return list of string
        litteral text.
    :returns: None or str
        None or string litteral if verbose is set to ``False``.
    Examples
    ---------
    >>> from watex.utils.funcutils import listing_items_format
    >>> litems = ['hole_number', 'depth_top', 'depth_bottom', 'strata_name',
                'rock_name','thickness', 'resistivity', 'gamma_gamma',
                'natural_gamma', 'sp','short_distance_gamma', 'well_diameter']
    >>> listing_items_format (litems , 'Features' ,
                               'have been successfully drop.' ,
                              lstyle ='.', space=3)
    """
    out =''
    if not is_iterable(lst):
        lst=[lst]

    if hasattr (lst, '__array__'):
        if lst.ndim !=1:
            raise ValueError (" Can not print multidimensional array."
                              " Expect one dimensional array.")
    lst = list(lst)
    begintext = str(begintext); endtext=str(endtext)
    lstyle=  lstyle or bullet
    lstyle = str(lstyle)
    b= f"{begintext +':' } "
    if verbose :
        print(b, end=' ') if inline else (
            print(b)  if  begintext!='' else None)
    out += b +  ('\n' if not inline else ' ')
    for k, item in enumerate (lst):
        sp = ' ' * space
        if ( not enum and inline ): lstyle =''
        o = f"{sp}{str(k+1) if enum else bullet+ ' ' }{lstyle} {item}"
        if verbose:
            print (o , end=' ') if inline else print(o)
        out += o + ('\n' if not inline else ' ')

    en= ' ' + endtext if inline else endtext
    if verbose:
        print(en) if endtext !='' else None
    out +=en

    return None if verbose else out

def parse_attrs (attr, /, regex=None ):
    r""" Parse attributes using the regular expression.

    Remove all string non-alphanumeric and some operator indicators,  and
    fetch attributes names.

    Parameters
    -----------

    attr: str, text litteral containing the attributes
        names

    regex: `re` object, default is
        Regular expresion object. the default is::

            >>> import re
            >>> re.compile (r'per|mod|times|add|sub|[_#&*@!_,;\s-]\s*',
                                flags=re.IGNORECASE)
    Returns
    -------
    attr: List of attributes

    Example
    ---------
    >>> from watex.utils.funcutils import parse_attrs
    >>> parse_attrs('lwi_sub_ohmSmulmagnitude')
    ... ['lwi', 'ohmS', 'magnitude']


    """
    regex = regex or re.compile (r'per|mod|times|add|sub|[_#&*@!_,;\s-]\s*',
                        flags=re.IGNORECASE)
    attr= list(filter (None, regex.split(attr)))
    return attr

def url_checker (url: str , install:bool = False,
                 raises:str ='ignore')-> bool :
    """
    check whether the URL is reachable or not.

    function uses the requests library. If not install, set the `install`
    parameter to ``True`` to subprocess install it.

    Parameters
    ------------
    url: str,
        link to the url for checker whether it is reachable
    install: bool,
        Action to install the 'requests' module if module is not install yet.
    raises: str
        raise errors when url is not recheable rather than returning ``0``.
        if `raises` is ``ignore``, and module 'requests' is not installed, it
        will use the django url validator. However, the latter only assert
        whether url is right but not validate its reachability.

    Returns
    --------
        ``True``{1} for reacheable and ``False``{0} otherwise.

    Example
    ----------
    >>> from watex.utils.funcutils import url_checker
    >>> url_checker ("http://www.example.com")
    ...  0 # not reacheable
    >>> url_checker ("https://watex.readthedocs.io/en/latest/api/watex.html")
    ... 1

    """
    isr =0 ; success = False

    regex = re.compile(
        r'^(?:http|ftp)s?://' # http:// or https://
        #domain...
        r'(?:(?:[A-Z0-9](?:[A-Z0-9-]{0,61}[A-Z0-9])?\.)+(?:[A-Z]{2,6}\.?|[A-Z0-9-]{2,}\.?)|'
        r'localhost|' #localhost...
        r'\d{1,3}\.\d{1,3}\.\d{1,3}\.\d{1,3})' # ...or ip
        r'(?::\d+)?' # optional port
        r'(?:/?|[/?]\S+)$', re.IGNORECASE
        )

    try :
        import requests
    except ImportError:
        if install:
            success  = is_installing('requests', DEVNULL=True)
        if not success:
            if raises=='raises':
                raise ModuleNotFoundError(
                    "auto-installation of 'requests' failed."
                    " Install it mannually.")

    else : success=True

    if success:
        try:
            get = requests.get(url) #Get Url
            if get.status_code == 200: # if the request succeeds
                isr =1 # (f"{url}: is reachable")

            else:
                warnings.warn(
                    f"{url}: is not reachable, status_code: {get.status_code}", stacklevel=2)
                isr =0

        except requests.exceptions.RequestException as e:
            if raises=='raises':
                raise SystemExit(f"{url}: is not reachable \nErr: {e}")
            else: isr =0

    if not success :
        # use django url validation regex
        # https://github.com/django/django/blob/stable/1.3.x/django/core/validators.py#L45
        isr = 1 if re.match(regex, url) is not None else 0

    return isr


def shrunkformat (text: str | Iterable[Any] ,
                  chunksize: int =7 , insert_at: str = None,
                  sep =None,
                 ) :
    """ format class and add ellipsis when classes are greater than maxview

    :param text: str - a text to shrunk and format. Can also be an iterable
        object.
    :param chunksize: int, the size limit to keep in the formatage text. *default*
        is ``7``.
    :param insert_at: str, the place to insert the ellipsis. If ``None``,
        shrunk the text and put the ellipsis, between the text beginning and
        the text endpoint. Can be ``beginning``, or ``end``.
    :param sep: str if the text is delimited by a kind of character, the `sep`
        parameters could be usefull so it would become a starting point for
        word counting. *default*  is `None` which means word is counting from
        the space.

    :example:

    >>> import numpy as np
    >>> from watex.utils.funcutils import shrunkformat
    >>> text=" I'm a long text and I will be shrunked and replaced by ellipsis."
    >>> shrunkformat (text)
    ... 'Im a long ... and replaced by ellipsis.'
    >>> shrunkformat (text, insert_at ='end')
    ...'Im a long ... '
    >>> arr = np.arange(30)
    >>> shrunkformat (arr, chunksize=10 )
    ... '0 1 2 3 4  ...  25 26 27 28 29'
    >>> shrunkformat (arr, insert_at ='begin')
    ... ' ...  26 27 28 29'

    """
    is_str = False
    chunksize = int (_assert_all_types(chunksize, float, int))

    regex = re.compile (r"(begin|start|beg)|(end|close|last)")
    insert_at = str(insert_at).lower().strip()
    gp = regex.search (insert_at)
    if gp is not None:
        if gp.group (1) is not None:
            insert_at ='begin'
        elif gp.group(2) is not None:
            insert_at ='end'
        if insert_at is None:
            warnings.warn(f"Expect ['begining'|'end'], got {insert_at!r}"
                          " Default value is used instead.", stacklevel=2)
    if isinstance(text , str):
        textsplt = text.strip().split(sep) # put text on list
        is_str =True

    elif hasattr (text , '__iter__'):
        textsplt = list(text )

    if len(textsplt) < chunksize :
        return  text

    if is_str :
        rl = textsplt [:len(textsplt)//2][: chunksize//2]
        ll= textsplt [len(textsplt)//2:][-chunksize//2:]

        if sep is None: sep =' '
        spllst = [f'{sep}'.join ( rl), f'{sep}'.join ( ll)]

    else : spllst = [
        textsplt[: chunksize//2 ] ,textsplt[-chunksize//2:]
        ]
    if insert_at =='begin':
        spllst.insert(0, ' ... ') ; spllst.pop(1)
    elif insert_at =='end':
        spllst.pop(-1) ; spllst.extend ([' ... '])

    else :
        spllst.insert (1, ' ... ')

    spllst = spllst if is_str else str(spllst)

    return re.sub(r"[\[,'\]]", '', ''.join(spllst),
                  flags=re.IGNORECASE
                  )


def is_installing (
        module: str ,
        upgrade: bool=True ,
        action: bool=True,
        DEVNULL: bool=False,
        verbose: int=0,
        **subpkws
    )-> bool:
    """ Install or uninstall a module/package using the subprocess
    under the hood.

    Parameters
    ------------
    module: str,
        the module or library name to install using Python Index Package `PIP`

    upgrade: bool,
        install the lastest version of the package. *default* is ``True``.

    DEVNULL:bool,
        decline the stdoutput the message in the console

    action: str,bool
        Action to perform. 'install' or 'uninstall' a package. *default* is
        ``True`` which means 'intall'.

    verbose: int, Optional
        Control the verbosity i.e output a message. High level
        means more messages. *default* is ``0``.

    subpkws: dict,
        additional subprocess keywords arguments
    Returns
    ---------
    success: bool
        whether the package is sucessfully installed or not.

    Example
    --------
    >>> from watex import is_installing
    >>> is_installing(
        'tqdm', action ='install', DEVNULL=True, verbose =1)
    >>> is_installing(
        'tqdm', action ='uninstall', verbose =1)
    """
    #implement pip as subprocess
    # refer to https://pythongeeks.org/subprocess-in-python/
    if not action:
        if verbose > 0 :
            print("---> No action `install`or `uninstall`"
                  f" of the module {module!r} performed.")
        return action  # DO NOTHING

    success=False

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

        success=True

    except:

        if verbose > 0 :
            print(f'---> Module {module!r} {action_msg} failed. Please use'
                f' the following command: {cmdg} to manually do it.')
    else :
        if verbose > 0:
            print(f"{action_msg.capitalize()} of `{module}` "
                      "and dependancies was successfully done!")

    return success

def smart_strobj_recognition(
        name: str  ,
        container: List | Tuple | Dict[Any, Any ],
        stripitems: str | List | Tuple = '_',
        deep: bool = False,
) -> str :
    """ Find the likelihood word in the whole containers and
    returns the value.

    :param name: str - Value of to search. I can not match the exact word in
    the `container`
    :param container: list, tuple, dict- container of the many string words.
    :param stripitems: str - 'str' items values to sanitize the  content
        element of the dummy containers. if different items are provided, they
        can be separated by ``:``, ``,`` and ``;``. The items separators
        aforementioned can not  be used as a component in the `name`. For
        isntance::

            name= 'dipole_'; stripitems='_' -> means remove the '_'
            under the ``dipole_``
            name= '+dipole__'; stripitems ='+;__'-> means remove the '+' and
            '__' under the value `name`.

    :param deep: bool - Kind of research. Go deeper by looping each items
         for find the initials that can fit the name. Note that, if given,
         the first occurence should be consider as the best name...

    :return: Likelihood object from `container`  or Nonetype if none object is
        detected.

    :Example:
        >>> from watex.utils.funcutils import smart_strobj_recognition
        >>> from watex.methods import ResistivityProfiling
        >>> rObj = ResistivityProfiling(AB= 200, MN= 20,)
        >>> smart_strobj_recognition ('dip', robj.__dict__))
        ... None
        >>> smart_strobj_recognition ('dipole_', robj.__dict__))
        ... dipole
        >>> smart_strobj_recognition ('dip', robj.__dict__,deep=True )
        ... dipole
        >>> smart_strobj_recognition (
            '+_dipole___', robj.__dict__,deep=True , stripitems ='+;_')
        ... 'dipole'

    """

    stripitems =_assert_all_types(stripitems , str, list, tuple)
    container = _assert_all_types(container, list, tuple, dict)
    ix , rv = None , None

    if isinstance (stripitems , str):
        for sep in (':', ",", ";"): # when strip ='a,b,c' seperated object
            if sep in stripitems:
                stripitems = stripitems.strip().split(sep) ; break
        if isinstance(stripitems, str):
            stripitems =[stripitems]

    # sanitize the name.
    for s in stripitems :
        name = name.strip(s)

    if isinstance(container, dict) :
        #get only the key values and lower them
        container_ = list(map (lambda x :x.lower(), container.keys()))
    else :
        # for consistency put on list if values are in tuple.
        container_ = list(container)

    # sanitize our dummny container item ...
    #container_ = [it.strip(s) for it in container_ for s in stripitems ]
    if name.lower() in container_:
        try:
            ix = container_.index (name)
        except ValueError:
            raise AttributeError(f"{name!r} attribute is not defined")

    if deep and ix is None:
        # go deeper in the search...
        for ii, n in enumerate (container_) :
            if n.find(name.lower())>=0 :
                ix =ii ; break

    if ix is not None:
        if isinstance(container, dict):
            rv= list(container.keys())[ix]
        else : rv= container[ix]

    return  rv

def repr_callable_obj(obj: F  , skip = None ):
    """ Represent callable objects.

    Format class, function and instances objects.

    :param obj: class, func or instances
        object to format.
    :param skip: str ,
        attribute name that is not end with '_' and whom it needs to be
        skipped.

    :Raises: TypeError - If object is not a callable or instanciated.

    :Examples:

    >>> from watex.utils.funcutils import repr_callable_obj
    >>> from watex.methods.electrical import  ResistivityProfiling
    >>> repr_callable_obj(ResistivityProfiling)
    ... 'ResistivityProfiling(station= None, dipole= 10.0,
            auto_station= False, kws= None)'
    >>> robj= ResistivityProfiling (AB=200, MN=20, station ='S07')
    >>> repr_callable_obj(robj)
    ... 'ResistivityProfiling(AB= 200, MN= 20, arrangememt= schlumberger, ... ,
        dipole= 10.0, station= S07, auto= False)'
    >>> repr_callable_obj(robj.fit)
    ... 'fit(data= None, kws= None)'

    """
    regex = re.compile (r"[{'}]")

    # inspect.formatargspec(*inspect.getfullargspec(cls_or_func))
    if not callable(obj) and not hasattr(obj, '__dict__'):
        raise TypeError (
            f'Format only callabe objects: Got {type (obj).__name__!r}')

    if callable(obj):
        cls_or_func_signature = inspect.signature(obj)
        objname = obj.__name__
        PARAMS_VALUES = {k: None if v.default is (inspect.Parameter.empty
                         or ...) else v.default
                    for k, v in cls_or_func_signature.parameters.items()
                    # if v.default is not inspect.Parameter.empty
                    }
    elif hasattr(obj, '__dict__'):
        objname=obj.__class__.__name__
        PARAMS_VALUES = {k:v  for k, v in obj.__dict__.items()
                         if not (k.endswith('_') or k.startswith('_')
                                  # remove the dict objects
                                  or k.endswith('_kws') or k.endswith('_props')
                                 )
                         }
    if skip is not None :
        # skip some inner params
        # remove them as the main function or class params
        if isinstance(skip, (tuple, list, np.ndarray)):
            skip = list(map(str, skip ))
            exs = [key for key in PARAMS_VALUES.keys() if key in skip]
        else:
            skip =str(skip).strip()
            exs = [key for key in PARAMS_VALUES.keys() if key.find(skip)>=0]

        for d in exs:
            PARAMS_VALUES.pop(d, None)

    # use ellipsis as internal to stdout more than seven params items
    if len(PARAMS_VALUES) >= 7 :
        f = {k:PARAMS_VALUES.get(k) for k in list(PARAMS_VALUES.keys())[:3]}
        e = {k:PARAMS_VALUES.get(k) for k in list(PARAMS_VALUES.keys())[-3:]}

        PARAMS_VALUES= str(f) + ', ... , ' + str(e )

    return str(objname) + '(' + regex.sub('', str (PARAMS_VALUES)
                                          ).replace(':', '=') +')'


def accept_types (
        *objtypes: list ,
        format: bool = False
        ) -> List[str] | str :
    """ List the type format that can be accepted by a function.

    :param objtypes: List of object types.
    :param format: bool - format the list of the name of objects.
    :return: list of object type names or str of object names.

    :Example:
        >>> import numpy as np; import pandas as pd
        >>> from watex.utils.funcutils import accept_types
        >>> accept_types (pd.Series, pd.DataFrame, tuple, list, str)
        ... "'Series','DataFrame','tuple','list' and 'str'"
        >>> atypes= accept_types (
            pd.Series, pd.DataFrame,np.ndarray, format=True )
        ..."'Series','DataFrame' and 'ndarray'"
    """
    return smart_format(
        [f'{o.__name__}' for o in objtypes]
        ) if format else [f'{o.__name__}' for o in objtypes]

def read_from_excelsheets(erp_file: str = None ) -> List[DataFrame]:

    """ Read all Excelsheets and build a list of dataframe of all sheets.

    :param erp_file:
        Excell workbooks containing `erp` profile data.

    :return: A list composed of the name of `erp_file` at index =0 and the
      datataframes.

    """

    allfls:Dict [str, Dict [T, List[T]] ] = pd.read_excel(
        erp_file, sheet_name=None)

    list_of_df =[os.path.basename(os.path.splitext(erp_file)[0])]
    for _sheets , values in allfls.items():
        list_of_df.append(pd.DataFrame(values))

    return list_of_df

def check_dimensionality(obj, data, z, x):
    """ Check dimensionality of data and fix it.

    :param obj: Object, can be a class logged or else.
    :param data: 2D grid data of ndarray (z, x) dimensions.
    :param z: array-like should be reduced along the row axis.
    :param x: arraylike should be reduced along the columns axis.

    """
    def reduce_shape(Xshape, x, axis_name=None):
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


def smart_format(iter_obj, choice ='and'):
    """ Smart format iterable object.

    :param iter_obj: iterable obj
    :param choice: can be 'and' or 'or' for optional.

    :Example:
        >>> from watex.utils.funcutils import smart_format
        >>> smart_format(['model', 'iter', 'mesh', 'data'])
        ... 'model','iter','mesh' and 'data'
    """
    str_litteral =''
    try:
        iter(iter_obj)
    except:  return f"{iter_obj}"

    iter_obj = [str(obj) for obj in iter_obj]
    if len(iter_obj) ==1:
        str_litteral= ','.join([f"{i!r}" for i in iter_obj ])
    elif len(iter_obj)>1:
        str_litteral = ','.join([f"{i!r}" for i in iter_obj[:-1]])
        str_litteral += f" {choice} {iter_obj[-1]!r}"
    return str_litteral

def make_introspection(Obj: object , subObj: Sub[object])->None:
    """ Make introspection by using the attributes of instance created to
    populate the new classes created.

    :param Obj: callable
        New object to fully inherits of `subObject` attributes.

    :param subObj: Callable
        Instance created.
    """
    # make introspection and set the all  attributes to self object.
    # if Obj attribute has the same name with subObj attribute, then
    # Obj attributes get the priority.
    for key, value in  subObj.__dict__.items():
        if not hasattr(Obj, key) and key  != ''.join(['__', str(key), '__']):
            setattr(Obj, key, value)

def cpath (savepath: str=None , dpath: str='_default_path_'):
    """ Control the existing path and create one of it does not exist.

    :param savepath: Pathlike obj, str
    :param dpath: str, default pathlike obj
    """
    # if dpath is None:
    #     #file, _= os.path.splitext(os.path.basename(__file__))
    #     #dpath = ''.join(['_', file,'_']) #.replace('.py', '')
    #     dpath ='_dpath_'
    if savepath is None :
        savepath  = os.path.join(os.getcwd(), dpath)
        try:os.mkdir(savepath)
        except: pass
    if savepath is not None:
        try :
            if not os.path.isdir(savepath):
                os.mkdir(savepath)#  mode =0o666)
        except : pass
    return savepath


def sPath (name_of_path:str):
    """ Savepath func. Create a path  with `name_of_path` if path not exists.

    :param name_of_path: str, Path-like object. If path does not exist,
        `name_of_path` should be created.
    """

    try :
        savepath = os.path.join(os.getcwd(), name_of_path)
        if not os.path.isdir(savepath):
            os.mkdir(name_of_path)#  mode =0o666)
    except :
        warnings.warn("The path seems to be existed!", stacklevel=2)
        return
    return savepath


def format_notes(text:str , cover_str: str ='~', inline=70, **kws):
    """ Format note
    :param text: Text to be formated.

    :param cover_str: type of ``str`` to surround the text.

    :param inline: Nomber of character before going in liine.

    :param margin_space: Must be <1 and expressed in %. The empty distance
        between the first index to the inline text
    :Example:

        >>> from watex.utils import funcutils as func
        >>> text ='Automatic Option is set to ``True``.'\
            ' Composite estimator building is triggered.'
        >>>  func.format_notes(text= text ,
        ...                       inline = 70, margin_space = 0.05)

    """

    headnotes =kws.pop('headernotes', 'notes')
    margin_ratio = kws.pop('margin_space', 0.2 )
    margin = int(margin_ratio * inline)
    init_=0
    new_textList= []
    if len(text) <= (inline - margin):
        new_textList = text
    else :
        for kk, _char in enumerate (text):
            if kk % (inline - margin)==0 and kk !=0:
                new_textList.append(text[init_:kk])
                init_ =kk
            if kk ==  len(text)-1:
                new_textList.append(text[init_:])

    print('!', headnotes.upper(), ':')
    print(f'{cover_str * inline}')
    for _k in new_textList:
        fmtin_str ='{'+ f'0:>{margin}' +'}'
        print('{}{1:>2}{:<51}'.format(fmtin_str.format(cover_str), '', ))

    print('{}{1:>51}'.format(' '* (margin -1), cover_str * (inline -margin+1 )))


def sanitize_fdataset(
    _df: DataFrame
) ->Tuple[DataFrame, int]:
    """ Sanitize the feature dataset.

    Recognize the columns provided  by the users and resset according
    to the features labels disposals :attr:`~Features.featureLabels`."""

    utm_flag =0

    def getandReplace(
            optionsList:List[str],
            params:List[str], df:DataFrame
            ) -> List[str]:
        """
        Function to  get parames and replace to the main features params.

        :param optionsList:
            User options to qualified the features headlines.
        :type optionsList: list

        :param params: Exhaustive parameters names.
        :type params: list

        :param df: pd.DataFrame collected from `features_fn`.

        :return: sanitize columns
        :rtype: list
        """

        columns = [c.lower() for c in df.columns]
        for ii, celemnt in enumerate(columns):
            for listOption, param in zip(optionsList, params):
                 for option in listOption:
                     if param =='lwi':
                        if celemnt.find('eau')>=0 :
                            columns[ii]=param
                            break
                     if re.match(rf'^{option}+', celemnt):
                         columns[ii]=param
                         if columns[ii] =='east':
                             utm_flag =1
                         break

        return columns, utm_flag

    new_df_columns, utm_flag = getandReplace(
        optionsList=P.param_options,
            params=P.param_ids,
            df= _df
                                  )
    df = pd.DataFrame(data=_df.to_numpy(), columns= new_df_columns)
    return df , utm_flag


def interpol_scipy (
        x_value,
        y_value,
        x_new,
        kind="linear",
        plot=False,
        fill="extrapolate"
        ):

    """
    function to interpolate data

    Parameters
    ------------
    * x_value : np.ndarray
        value on array data : original abscissA

    * y_value : np.ndarray
        value on array data : original coordinates (slope)

    * x_new  : np.ndarray
        new value of absciss you want to interpolate data

    * kind  : str
        projection kind maybe : "linear", "cubic"

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

    func_=spi.interp1d(
        x_value,
        y_value,
        kind=kind,
        fill_value=fill
        )
    y_new=func_(x_new)
    if plot :
        plt.plot(
        x_value,
        y_value,
        "o",
        x_new,
        y_new,
        "--"
        )
        plt.legend(["data", "linear","cubic"],loc="best")
        plt.show()

    return y_new


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
        str ; char , new_char without the removed word .

    Examples
    ---------
    >>> from watex.utils import funcutils as func
    >>> ch ='AMTAVG 7.76: "K1.fld", Dated 99-01-01,AMTAVG,
    ...    Processed 11 Jul 17 AMTAVG'
    >>> ss=func._remove_str_word(char=ch, word_to_remove='AMTAVG',
    ...                             deep_remove=False)
    >>> print(ss)

    """
    if type(char) is not str : char =str(char)
    if type(word_to_remove) is not str : word_to_remove=str(word_to_remove)

    if deep_remove :
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
        >>> from watex.utils  import funcutils as func
        >>> path =  data/ K6.stn
        >>> with open (path, 'r', encoding='utf8') as f :
        ...                     data= f.readlines()
        >>>  print(func.stn_check_split_type(data_lines=data))

    """

    split_type =[',', ':',' ',';' ]
    data_to_read =[]
    # change the data if data is not dtype string elements.
    if isinstance(data_lines, np.ndarray):
        if data_lines.dtype in ['float', 'int', 'complex']:
            data_lines=data_lines.astype('<U12')
        data_lines= data_lines.tolist()

    if isinstance(data_lines, list):
        for _ii, item in enumerate(data_lines[:int(len(data_lines)/3)]):
             data_to_read.append(item)
             # be sure the list is str item .
             data_to_read=[''.join([str(item) for item in data_to_read])]

    elif isinstance(data_lines, str): data_to_read=[str(data_lines)]

    for _jj, sep  in enumerate(split_type) :
        if data_to_read[0].find(sep) > 0 :
            if data_to_read[0].count(sep) >= 2 * len(data_lines)/3:
                if sep == ' ': return  None  # use None more conventional
                else : return sep

def minimum_parser_to_write_edi (edilines, parser = '='):
    """
    This fonction validates edifile for writing , string with egal.
    we assume that dictionnary in list will be for definemeasurment
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
            raise EDIError (
             f'None <"="> found on this item<{edilines[ii]}> of '
            ' the edilines list. list can not be parsed. Please'
            ' put egal sign "=" between key and value '
            )

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

def display_infos(infos, **kws):
    """ Display unique element on list of array infos

    :param infos: Iterable object to display.
    :param header: Change the `header` to other names.

    :Example:
    >>> from watex.utils.funcutils import display_infos
    >>> ipts= ['river water', 'fracture zone', 'granite', 'gravel',
         'sedimentary rocks', 'massive sulphide', 'igneous rocks',
         'gravel', 'sedimentary rocks']
    >>> display_infos('infos= ipts,header='TestAutoRocks',
                      size =77, inline='~')
    """

    inline =kws.pop('inline', '-')
    size =kws.pop('size', 70)
    header =kws.pop('header', 'Automatic rocks')

    if isinstance(infos, str ):
        infos =[infos]

    infos = list(set(infos))
    print(inline * size )
    mes= f'{header.capitalize()}({len(infos):02})'
    mes = f'{mes:^70}'
    print(mes)
    print(inline * size )
    am=''
    for ii in range(len(infos)):
        if (ii+1) %2 ==0:
            am = am + f'{ii+1:>4}.{infos[ii].capitalize():<30}'
            print(am)
            am=''
        else:
            am =f'{ii+1:>4}.{infos[ii].capitalize():<30}'
            if ii ==len(infos)-1:
                print(am)
    print(inline * size )

def fr_en_parser (f, delimiter =':'):
    r""" Parse the translated data file.

    :param f: translation file to parse.

    :param delimiter: str, delimiter.

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
        raise OSError(f'Input {f} is not a file. Please check your file.')

    with open(f, encoding ='utf8') as ft:
        data = ft.readlines()
        for row in data :
            if row in ( '\n', ' '):
                continue
            fr, en = row.strip().split(delimiter)
            yield([fr, en])

def convert_csvdata_from_fr_to_en(csv_fn, pf, destfile = 'pme.en.csv',
                                  savepath =None, delimiter =':'):
    r""" Translate variable data from french csv data  to english with
    parser file.

    :param csv_fn: data collected in csv format.

    :param pf: parser file.

    :param destfile: str,  Destination file, outputfile.

    :param savepath: Path-Like object, save data to a path.

    :Example:
        # to execute this script, we need to import the two modules below
        >>> import os
        >>> import csv
        >>> from watex.utils.funcutils import convert_csvdata_from_fr_to_en
        >>> path_pme_data = r'C:/Users\Administrator\Desktop\__elodata
        >>> datalist=convert_csvdata_from_fr_to_en(
            os.path.join( path_pme_data, _enuv2.csv') ,
            os.path.join(path_pme_data, pme.parserf.md')
                         savefile = 'pme.en.cv')
    """
    # read the parser file and separed english from french
    parser_data = list(fr_en_parser (pf,delimiter) )

    with open (csv_fn, encoding ='utf8') as csv_f :
        csv_reader = csv.reader(csv_f)
        csv_data =[ row for row in csv_reader]
    # get the index of the last substring row
    ix = csv_data [0].index ('Industry_type')
    # separateblock from two
    csv_1b = [row [:ix +1] for row in csv_data]
    csv_2b =[row [ix+1:] for row in csv_data ]
    # make a copy of csv_1b
    csv_1bb= copy.deepcopy(csv_1b)

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

def parse_md_data (pf , delimiter =':'):

    if not os.path.isfile (pf):
        raise OSError( " Unable to detect the parser file. "
                      "Need a Path-like object ")

    with open(pf, encoding ='utf8') as f:
        pdata = f.readlines ()
    for row in pdata :
        if row in ('\n', ' '):
            continue
        fr, en = row.strip().split(delimiter)
        fr = sanitize_unicode_string(fr)
        en = en.strip()
        # if capilize, the "I" inside the
        #text should be in lowercase
        # it is better to upper the first
        # character after striping the whole
        # string
        en = list(en)
        en[0] = en[0].upper()
        en = "".join(en)

        yield fr, en

def sanitize_unicode_string (str_) :
    """ Replace all spaces and remove all french accents characters.

    :Example:
    >>> from watex.utils.funcutils import sanitize_unicode_string
    >>> sentence ='Nos clients sont extrêmement satisfaits '
        'de la qualité du service fourni. En outre Nos clients '
            'rachètent frequemment nos "services".'
    >>> sanitize_unicode_string  (sentence)
    ... 'nosclientssontextrmementsatisfaitsdelaqualitduservice'
        'fournienoutrenosclientsrachtentfrequemmentnosservices'
    """
    sp_re = re.compile (r"[.'()-\\/’]")
    e_re = re.compile(r'[éèê]')
    a_re= re.compile(r'[àâ]')

    str_= re.sub(r'\s+', '', str_.strip().lower())

    for cobj , repl  in zip ( (sp_re, e_re, a_re),
                             ("", 'e', 'a')):
        str_ = cobj.sub(repl, str_)

    return str_

def read_main (csv_fn , pf , delimiter =':',
               destfile ='pme.en.csv') :

    parser_data = list(parse_md_data(pf, delimiter) )
    parser_dict =dict(parser_data)

    with open (csv_fn, encoding ='utf8') as csv_f :
        csv_reader = csv.reader(csv_f)
        csv_data =[ row for row in csv_reader]

    # get the index of the last substring row
    # and separate block into two from "Industry_type"
    ix = csv_data [0].index ('Industry_type')

    csv_1b = [row [:ix +1] for row in csv_data]
    csv_2b =[row [ix+1:] for row in csv_data ]
    # make a copy of csv_1b
    csv_1bb= copy.deepcopy(csv_1b)
    copyd = copy.deepcopy(csv_1bb); is_missing =list()

    # skip the first two rows
    for ii, rowline in enumerate( csv_1bb[3:]) :
        for jj , row in enumerate(rowline):
            row = row.strip()
            row = sanitize_unicode_string(row )
            csv_1bb[3:][ii][jj] = row

    #collect the missing values
    for ii, rowline in enumerate( csv_1bb[3:]) :
        for jj , row in enumerate(rowline):
            if row not in parser_dict.keys():
                is_missing.append(copyd[3:][ii][jj])
    is_missing = list(set(is_missing))

    # merge the prior two blocks and build the dataframe
    new_csv_list = [r1 + r2 for r1, r2 in zip(csv_1bb, csv_2b )]
    df = pd.DataFrame (
        np.array(new_csv_list [1:]),
        columns =new_csv_list [0]
                       )
    for key, value in parser_dict.items():
        # perform operation in place and return None
        df.replace (key, value, inplace =True )


    df.to_csv (destfile)
    return  df , is_missing


def _isin (
        arr: ArrayLike | List [float] ,
        subarr: Sub [ArrayLike] |Sub[List[float]] | float,
        return_mask:bool=False,
) -> bool :
    """ Check whether the subset array `subcz` is in  `cz` array.

    :param arr: Array-like - Array of item elements
    :param subarr: Array-like, float - Subset array containing a subset items.
    :param return_mask: bool, return the mask where the element is in `arr`.

    :return: True if items in  test array `subarr` are in array `arr`.

    """
    arr = np.array (arr );  subarr = np.array(subarr )

    return (True if True in np.isin (arr, subarr) else False
            ) if not return_mask else np.isin (arr, subarr)

def _assert_all_types (
    obj: object ,
    *expected_objtype: type,
    objname:str=None,
 ) -> object:
    """ Quick assertion of object type. Raises a `TypeError` if wrong type
    is passed as an argument. For polishing the error message, one can add
    the object name `objname` for specifying the object that raises errors
    for letting the users to be aware of the reason of failure."""
    # if np.issubdtype(a1.dtype, np.integer):
    if not isinstance( obj, expected_objtype):
        n=str(objname) + ' expects' if objname is not None else 'Expects'
        raise TypeError (
            f"{n} type{'s' if len(expected_objtype)>1 else ''} "
            f"{smart_format(tuple (o.__name__ for o in expected_objtype))}"
            f" but {type(obj).__name__!r} is given.")

    return obj


def savepath_ (nameOfPath):
    """
    Shortcut to create a folder
    :param nameOfPath: Path name to save file
    :type nameOfPath: str

    :return:
        New folder created. If the `nameOfPath` exists, will return ``None``
    :rtype:str

    """

    try :
        savepath = os.path.join(os.getcwd(), nameOfPath)
        if not os.path.isdir(savepath):
            os.mkdir(nameOfPath)#  mode =0o666)
    except :
        warnings.warn("The path seems to be existed !", stacklevel=2)
        return
    return savepath


def drawn_boundaries(erp_data, appRes, index):
    """
    Function to drawn anomaly boundary
    and return the anomaly with its boundaries

    :param erp_data: erp profile
    :type erp_data: array_like or list

    :param appRes: resistivity value of minimum pk anomaly
    :type appRes: float

    :param index: index of minimum pk anomaly
    :type index: int

    :return: anomaly boundary
    :rtype: list of array_like

    """
    f = 0 # flag to mention which part must be calculated
    if index ==0 :
        f = 1 # compute only right part
    elif appRes ==erp_data[-1]:
        f=2 # compute left part

    def loop_sideBound(term):
        """
        loop side bar from anomaly and find the term side

        :param term: is array of left or right side of anomaly.
        :type term: array

        :return: side bar
        :type: array_like
        """
        tem_drawn =[]
        maxT=0

        for _ii, tem_rho in enumerate(term) :

            diffRes_betw_2pts= tem_rho - appRes
            if diffRes_betw_2pts > maxT :
                maxT = diffRes_betw_2pts
                tem_drawn.append(tem_rho)
            elif diffRes_betw_2pts < maxT :
                # rho_limit = tem_rho
                break
        return np.array(tem_drawn)
    # first broke erp profile from the anomalies
    if f ==0 or f==2 :
        left_term = erp_data[:index][::-1] # flip left term  for looping
        # flip again to keep the order
        left_limit = loop_sideBound(term=left_term)[::-1]

    if f==0 or f ==1 :
        right_term= erp_data[index :]
        right_limit=loop_sideBound(right_term)
    # concat right and left to get the complete anomaly
    if f==2:
        anomalyBounds = np.append(left_limit,appRes)

    elif f ==1 :
        anomalyBounds = np.array([appRes]+ right_limit.tolist())
    else:
        left_limit = np.append(left_limit, appRes)
        anomalyBounds = np.concatenate((left_limit, right_limit))

    return appRes, index, anomalyBounds

def serialize_data(
        data,
        filename=None,
        force=True,
        savepath=None,
        verbose:int =0
     ):
    """ Store a data into a binary file

    :param data: Object
        Object to store into a binary file.
    :param filename: str
        Name of file to serialize. If 'None', should create automatically.
    :param savepath: str, PathLike object
         Directory to save file. If not exists should automaticallycreate.
    :param force: bool
        If ``True``, remove the old file if it exists, otherwise will
        create a new incremenmted file.
    :param verbose: int, get more message.
    :return: dumped or serialized filename.

    :Example:

        >>> import numpy as np
        >>> import watex.utils.coreutils import serialize_data
        >>> data = np.arange(15)
        >>> file = serialize_data(data, filename=None,  force=True,
        ...                          savepath =None, verbose =3)
        >>> file
    """

    def _cif(filename, force):
        """ Control the file. If `force` is ``True`` then remove the old file,
        Otherwise create a new file with datetime infos."""
        f = copy.deepcopy(filename)
        if force :
            os.remove(filename)
            if verbose >2: print(f" File {os.path.basename(filename)!r} "
                      "has been removed. ")
            return None
        else :
            # that change the name in the realpath
            f= os.path.basename(f).replace('.pkl','') + \
                f'{datetime.datetime.now()}'.replace(':', '_')+'.pkl'
            return f

    if filename is not None:
        file_exist =  os.path.isfile(filename)
        if file_exist:
            filename = _cif (filename, force)
    if filename is None:
        filename =f'__mymemoryfile.{datetime.datetime.now()}__'
        filename =filename.replace(' ', '_').replace(':', '-')
    if not isinstance(filename, str):
        raise TypeError(f"Filename needs to be a string not {type(filename)}")
    if filename.endswith('.pkl'):
        filename = filename.replace('.pkl', '')

    _logger.info (
        f"Save data to {'memory' if filename.find('memo')>=0 else filename}.")
    try :
        joblib.dump(data, f'{filename}.pkl')
        filename +='.pkl'
        if verbose > 2:
            print(f'Data dumped in `{filename} using to `~.externals.joblib`!')
    except :
        # Now try to pickle data Serializing data
        with open(filename, 'wb') as wfile:
            pickle.dump( data, wfile)
        if verbose >2:
            print( 'Data are well serialized using Python pickle module.`')
    # take the real path of the filename
    filename = os.path.realpath(filename)

    if savepath is  None:
        dirname ='_memory_'
        try : savepath = sPath(dirname)
        except :
            # for consistency
            savepath = os.getcwd()
    if savepath is not None:
        try:
            shutil.move(filename, savepath)
        except :
            file = _cif (os.path.join(savepath,
                                      os.path.basename(filename)), force)
            if not force:
                os.rename(filename, os.path.join(savepath, file) )
            if file is None:
                #take the file  in current word
                file = os.path.join(os.getcwd(), filename)
                shutil.move(filename, savepath)
            filename = os.path.join(savepath, file)

    if verbose > 0:
            print(f"Data are well stored in {savepath!r} directory.")

    return os.path.join(savepath, filename)

def load_serialized_data (filename, verbose=0):
    """
    Load data from dumped file.

    :param filename: str or path-like object
        Name of dumped data file.
    :return: Data reloaded from dumped file.

    :Example:

        >>> from watex.utils.functils import load_serialized_data
        >>> data = load_serialized_data(
        ...    filename = '_memory_/__mymemoryfile.2021-10-29_14-49-35.647295__.pkl',
        ...    verbose =3)

    """
    if not isinstance(filename, str):
        raise TypeError(f'filename should be a <str> not <{type(filename)}>')

    if not os.path.isfile(filename):
        raise FileExistsError(f"File {filename!r} does not exist.")

    _filename = os.path.basename(filename)
    _logger.info(
        f"Loading data from {'memory' if _filename.find('memo')>=0 else _filename}.")

    data =None
    try :
        data= joblib.load(filename)
        if verbose >2:
            (f"Data from {_filename !r} are sucessfully"
             " reloaded using ~.externals.joblib`!")
    except :
        if verbose >2:
            print(f"Nothing to reload. It's seems data from {_filename!r}"
                      " are not dumped using ~external.joblib module!")

        with open(filename, 'rb') as tod:
            data= pickle.load (tod)

        if verbose >2: print(f"Data from `{_filename!r} are well"
                      " deserialized using Python pickle module.`!")

    is_none = data is None
    if verbose > 0:
        if is_none :
            print("Unable to deserialize data. Please check your file.")
        else : print(f"Data from {_filename} have been sucessfully reloaded.")

    return data

def savejob(
    job ,
    savefile ,* ,
    protocol =None,
    append_versions=True,
    append_date=True,
    fix_imports= True,
    buffer_callback = None,
    **job_kws
    ):
    """ Quick save your job using 'joblib' or persistent Python pickle module

    Parameters
    -----------
    job: Any
        Anything to save, preferabaly a models in dict

    savefile: str, or path-like object
         name of file to store the model
         The *file* argument must have a write() method that accepts a
         single bytes argument. It can thus be a file object opened for
         binary writing, an io.BytesIO instance, or any other custom
         object that meets this interface.

    append_versions: bool, default =True
        Append the version of Joblib module or Python Pickle module following
        by the scikit-learn, numpy and also pandas versions. This is useful
        to have idea about previous versions for loading file when system or
        modules have been upgraded. This could avoid bottleneck when data
        have been stored for long times and user has forgotten the date and
        versions at the time the file was saved.

    append_date: bool, default=True,
       Append the date  of the day to the filename.

       .. versionadded:: 0.2.3

    protocol: int, optional
        The optional *protocol* argument tells the pickler to use the
        given protocol; supported protocols are 0, 1, 2, 3, 4 and 5.
        The default protocol is 4. It was introduced in Python 3.4, and
        is incompatible with previous versions.

        Specifying a negative protocol version selects the highest
        protocol version supported.  The higher the protocol used, the
        more recent the version of Python needed to read the pickle
        produced.

    fix_imports: bool, default=True,
        If *fix_imports* is True and *protocol* is less than 3, pickle
        will try to map the new Python 3 names to the old module names
        used in Python 2, so that the pickle data stream is readable
        with Python 2.

    buffer_call_back: int, optional
        If *buffer_callback* is None (the default), buffer views are
        serialized into *file* as part of the pickle stream.

        If *buffer_callback* is not None, then it can be called any number
        of times with a buffer view.  If the callback returns a false value
        (such as None), the given buffer is out-of-band; otherwise the
        buffer is serialized in-band, i.e. inside the pickle stream.

        It is an error if *buffer_callback* is not None and *protocol*
        is None or smaller than 5.

    job_kws: dict,
        Additional keywords arguments passed to :func:`joblib.dump`.

    Returns
    --------
    savefile: str,
        returns the filename
    """
    def remove_extension(fn, ex):
        """Remove extension either joblib or pickle """
        return fn.replace (ex, '')

    import sklearn

    versions = f'sklearn_v{sklearn.__version__}.numpy_v{np.__version__}.pandas_v{pd.__version__}'
    date = datetime.datetime.now()

    savefile =str(savefile)
    if (
            '.joblib' in savefile or '.pkl' in savefile
            ):
        ex = '.joblib' if savefile.find('.joblib')>=0 else '.pkl'
        savefile = remove_extension(savefile ,  ex )

    if append_date:
        savefile +=f".{date}"

    if append_versions :
        savefile += ".{}"+ versions
    try :
        if append_versions:
            savefile += f".joblib_v{joblib.__version__}."

        joblib.dump(job, f'{savefile}.joblib', **job_kws)

    except :
        if append_versions:
            savefile +=f".pickle_v{pickle.__version__}.pkl"

        with open(savefile, 'wb') as wfile:
            pickle.dump( job, wfile, protocol= protocol,
                        fix_imports=fix_imports ,
                        buffer_callback=buffer_callback )

    return savefile

def find_position_from_sa(
        an_res_range,
        pos=None,
        selectedPk=None):
    """
    Function to select the main `pk` from both :func:`get_boundaries`.

    :paran an_res_range: anomaly resistivity range on |ERP| line.
    :type an_res_range: array_like

    :param pos: position of anomaly boundaries (inf and sup):
                anBounds = [90, 130]
                - 130 is max boundary and 90 the  min boundary
    :type pos: list

    :param selectedPk:

        User can set its own position of the right anomaly. Be sure that
        the value provided is right position .
        Could not compute again provided that `pos`
        is not `None`.

    :return: anomaly station position.
    :rtype: str 'pk{position}'

    :Example:

        >>> from watex.utils.funcutils import find_positon_from_sa
        >>> resan = np.array([168,130, 93,146,145])
        >>> pk= find_pk_from_selectedAn(
        ...    resan, pos=[90, 13], selectedPk= 'str20')
        >>> pk

    .. |ERP| replace:: Electrical Resistivity Profiling

    """
    #compute dipole length from pos
    if pos is not None :
        if isinstance(pos, list):
            pos =np.array(pos)
    if pos is None and selectedPk is None :
        raise ParameterNumberError(
            'Give at least the anomaly boundaries'
            ' before computing the selected anomaly position.')

    if selectedPk is not None :  # mean is given
        if isinstance(selectedPk, str):
            if selectedPk.isdigit() :
                sPk= int(selectedPk)
            elif selectedPk.isalnum():
                oss = ''.join([s for s in selectedPk
                    if s.isdigit()])
                sPk =int(oss)
        else :
            try :
                sPk = int(selectedPk)

            except : pass

        if pos is not None : # then compare the sPk and ps value
            try :
                if not pos.min()<= sPk<=pos.max():
                    warnings.warn(f'Wrong position given <{selectedPk}>.'
                                  ' Should compute new positions.', stacklevel=2)
                    _logger.debug(f'Wrong position given <{selectedPk}>.'
                                  'Should compute new positions.')

            except UnboundLocalError:
                print("local variable 'sPk' referenced before assignment")
            else :

                return f'pk{sPk}', an_res_range


        else :
            selectedPk=f'pk{sPk}'

            return selectedPk , an_res_range


    if isinstance(pos, list):
        pos =np.array(pos)
    if isinstance(an_res_range, list):
        an_res_range =np.array(an_res_range)
    dipole_length = (pos.max()-pos.min())/(len(an_res_range)-1)

    tem_loc = np.arange(
        pos.min(), pos.max()+dipole_length, dipole_length)

    # find min value of  collected anomalies values
    locmin = np.where (an_res_range==an_res_range.min())[0]
    if len(locmin) >1 : locmin =locmin[0]
    pk_= int(tem_loc[int(locmin)]) # find the min pk

    selectedPk=f'pk{pk_}'

    return selectedPk , an_res_range

def fmt_text(
        anFeatures=None,
        title = None,
        **kwargs) :
    """
    Function format text from anomaly features

    :param anFeatures: Anomaly features
    :type anFeatures: list or dict

    :param title: head lines
    :type title: list

    :Example:

        >>> from watex.utils.funcutils import fmt_text
        >>> fmt_text(anFeatures =[1,130, 93,(146,145, 125)])

    """
    if title is None:
        title = ['Ranking', 'rho(Ω.m)', 'position pk(m)', 'rho range(Ω.m)']
    inline =kwargs.pop('inline', '-')
    mlabel =kwargs.pop('mlabels', 100)
    line = inline * int(mlabel)

    #--------------------header ----------------------------------------
    print(line)
    tem_head ='|'.join([f'{i:^15}' for i in title[:-1]])
    tem_head +=f'|{title[-1]:^45}'
    print(tem_head)
    print(line)
    #-----------------------end header----------------------------------
    newF =[]
    if isinstance(anFeatures, dict):
        for keys, items in anFeatures.items():
            rrpos=keys.replace('_pk', '')
            rank=rrpos[0]
            pos =rrpos[1:]
            newF.append([rank, min(items), pos, items])

    elif isinstance(anFeatures, list):
        newF =[anFeatures]


    for anFeatures in newF:
        strfeatures ='|'.join([f'{str(i):^15}' \
                               for i in anFeatures[:-1]])
        try :
            iter(anFeatures[-1])
        except :
            strfeatures +=f'|{str(anFeatures[-1]):^45}'
        else :
            strfeatures += '|{:^45}'.format(
                ''.join([f'{str(i)} ' for i in anFeatures[-1]]))

        print(strfeatures)
        print(line)


def find_feature_positions (
        anom_infos,
        anom_rank,
        pks_rhoa_index,
        dl):
    """
    Get the pk bound from ranking of computed best points

    :param anom_infos:

        Is a dictionnary of best anomaly points computed from
        :func:`drawn_anomaly_boundaries2` when `pk_bounds` is not given.
        see :func:`find_position_bounds`

    :param anom_rank: Automatic ranking after selecting best points

    :param pk_rhoa_index:

        Is tuple of selected anomaly resistivity value and index in the whole
        |ERP| line. for instance:

            pks_rhoa_index= (80., 17)

        where "80" is the value of selected anomaly in ohm.m and "17" is the
        index of selected points in the |ERP| array.

    :param dl:

        Is the distance between two measurement as `dipole_length`. Provide
        the `dl` if the *default* value is not right.

    :returns:

        Refer to :doc:`.exmath.select_anomaly`

    .. |ERP| replace:: Electrical Resistivity Profiling

    """
    rank_code = f'{anom_rank}_pk'
    for key in anom_infos.keys():
        if rank_code in key:
            pk = float(key.replace(rank_code, ''))

            rhoa = list(pks_rhoa_index[anom_rank-1])[0]
            codec = key
            break

    ind_rhoa =np.where(anom_infos[codec] ==rhoa)[0]
    if len(ind_rhoa) ==0 : ind_rhoa =0
    leninf = len(anom_infos[codec][: int(ind_rhoa)])

    pk_min = pk - leninf * dl
    lensup =len(anom_infos[codec][ int(ind_rhoa):])
    pk_max =  pk + (lensup -1) * dl

    pos_bounds = (pk_min, pk_max)
    rhoa_bounds = (anom_infos[codec][0], anom_infos[codec][-1])

    return pk, rhoa, pos_bounds, rhoa_bounds, anom_infos[codec]


def find_position_bounds(
        pk,
        rhoa,
        rhoa_range,
        dl=10.
        ):
    """
    Find station position boundary indexed in |ERP| line.

    Useful to get the boundaries indexes `pk_boun_indexes` for |ERP|
    normalisation  when computing `anr` or else.

    .. |ERP| replace:: Electrical Resistivity Profiling

    :param pk: Selected anomaly station value
    :type pk: float

    :param rhoa: Selected anomaly value in ohm.m
    :type rhoa: float

    :rhoa_range: Selected anomaly values from `pk_min` to `pk_max`
    :rhoa_range: array_like

    :parm dl: see :func:`find_position_from_sa` docstring.

    :Example:

        >>> from watex.utils.funcutils import find_position_bounds
        >>> find_position_bounds(pk=110, rhoa=137,
                          rhoa_range=np.array([175,132,137,139,170]))


    """

    if isinstance(pk, str):
        pk = float(pk.replace(pk[0], '').replace('_pk', ''))

    index_rhoa = np.where(rhoa_range ==rhoa)[0]
    if len(index_rhoa) ==0 : index_rhoa =0

    leftlen = len(rhoa_range[: int(index_rhoa)])
    rightlen = len(rhoa_range[int(index_rhoa):])

    pk_min = pk - leftlen * dl
    pk_max =  pk + (rightlen  -1) * dl

    return pk_min, pk_max


def wrap_infos (
        phrase ,
        value ='',
        underline ='-',
        unit ='',
        site_number= '',
        **kws) :
    """Display info from anomaly details."""

    repeat =kws.pop('repeat', 77)
    intermediate =kws.pop('inter+', '')
    begin_phrase_mark= kws.pop('begin_phrase', '--|>')
    on = kws.pop('on', False)
    if not on: return ''
    else :
        print(underline * repeat)
        print(f'{begin_phrase_mark} {phrase:<50}',
              f'{value:<10} {unit}',
              f'{intermediate}', f"{site_number}")
        print(underline * repeat )

def drawn_anomaly_boundaries2(
        erp_data,
        appRes,
        index
        ):
    """
    Function to drawn anomaly boundary
    and return the anomaly with its boundaries

    :param erp_data: erp profile
    :type erp_data: array_like or list

    :param appRes: resistivity value of minimum pk anomaly
    :type appRes: float

    :param index: index of minimum pk anomaly
    :type index: int

    :return: anomaly boundary
    :rtype: list of array_like

    """
    f = 0 # flag to mention which part must be calculated
    if index ==0 :
        f = 1 # compute only right part
    elif appRes ==erp_data[-1]:
        f=2 # compute left part

    def loop_sideBound(term):
        """
        loop side bar from anomaly and find the term side

        :param term: is array of left or right side of anomaly.
        :type trem: array

        :return: side bar
        :type: array_like
        """
        tem_drawn =[]
        maxT=0

        for _ii, tem_rho in enumerate(term) :

            diffRes_betw_2pts= tem_rho - appRes
            if diffRes_betw_2pts > maxT :
                maxT = diffRes_betw_2pts
                tem_drawn.append(tem_rho)
            elif diffRes_betw_2pts < maxT :
                # rho_limit = tem_rho
                break
        # print(tem_drawn)
        return np.array(tem_drawn)
    # first broke erp profile from the anomalies
    if f==2 : # compute the left part
        # flip array and start backward counting
        temp_erp_data = erp_data [::-1]
        sbeg = appRes   # initialize value
        for ii, valan in enumerate(temp_erp_data):
            if valan >= sbeg:
                sbeg = valan
            elif valan < sbeg:
                left_term = erp_data[ii:]
                break

        left_term = erp_data[:index][::-1] # flip left term  for looping
        # flip again to keep the order
        left_limit = loop_sideBound(term=left_term)[::-1]

    if f==0 or f ==1 :
        right_term= erp_data[index :]
        right_limit=loop_sideBound(right_term)
    # concat right and left to get the complete anomaly
    if f==2:
        anomalyBounds = np.append(left_limit,appRes)

    elif f ==1 :
        anomalyBounds = np.array([[appRes]+ right_limit.tolist()])
    else:
        left_limit = np.append(left_limit, appRes)
        anomalyBounds = np.concatenate((left_limit, right_limit))

    return appRes, index, anomalyBounds

def get_boundaries(df):
    """
    Define anomaly boundary `upper bound` and `lowerbound` from
    :ref:`define_position_bounds` location.

    :param df: Dataframe pandas contained  the columns
                'pk', 'x', 'y', 'rho', 'dl'.
    returns:
        - `autoOption` triggered the automatic Option if nothing is specified
            into excelsheet.
        - `ves_loc`: Sounding curve location at pk
        - `posMinMax`: Anomaly boundaries composed of ``lower`` and ``upper``
            bounds.
           Specific names can be used  to define lower and upper bounds::

               `lower`: 'lower', 'inf', 'min', 'min', '1' or  'low'
               `upper`: 'upper', 'sup', 'maj', 'max', '2, or 'up'

        To define the sounding location, can use::
            `ves`:'ves', 'se', 'sond','vs', 'loc', '0' or 'dl'

    """
    shape_=[ 'V','W', 'U', 'H', 'M', 'C', 'K' ]
    type__= ['EC', 'NC', 'CP', 'CB2P']
        # - ``EC`` for Extensive conductive.
        # - ``NC`` for narrow conductive.
        # - ``CP`` for conductive PLANE
        # - ``CB2P`` for contact between two planes.
    shape =None
    type_ =None

    def recoverShapeOrTypefromSheet(
            listOfAddedArray,
            param):
        """ Loop the array and get whether an anomaly shape name is provided.

        :param listOfAddedArray: all Added array values except
         'pk', 'x', 'y', 'rho' are composed of list of addedArray.
        :param param: Can be main description of different `shape_` of `type__`

        :returns:
            - `shape` : 'V','W', 'U', 'H', 'M', 'C' or  'K' from sheet or
             `type` : 'EC', 'NC', 'CP', 'CB2P'
            - listOfAddedArray : list of added array
        """
        param_ =None
        for jj, colarray in enumerate(listOfAddedArray[::-1]):
            tem_=[str(ss).upper().strip() for ss in list(colarray)]
            for ix , elem in enumerate(tem_):
                for param_elm in param:
                    if elem ==param_elm :
                        # retrieves the shape and replace by np.nan value
                        listOfAddedArray[::-1][jj][ix]=np.nan
                        return param_elm , listOfAddedArray

        return param_, listOfAddedArray

    def mergeToOne(
            listOfColumns,
            _df
            ):
        """ Get data from other columns annd merge into one array.

        :param listOfColumns: Columns names
        :param _df: dataframe to retrieve data to one
        """
        new_array = np.full((_df.shape[0],), np.nan)
        listOfColumnData = [ _df[name].to_numpy() for name in listOfColumns ]
        # loop from backward so we keep the most important to the first row
        # close the main df that composed `pk`,`x`, `y`, and `rho`.
        # find the shape
        shape, listOfColumnData = recoverShapeOrTypefromSheet(listOfColumnData,
                                                        param =shape_)
        type_, listOfColumnData = recoverShapeOrTypefromSheet(listOfColumnData,
                                                        param =type__)

        for colarray in listOfColumnData[::-1]:
           for ix , val  in enumerate(colarray):
               try:
                   if not np.isnan(val) :
                       new_array[ix]=val
               except :pass

        return shape , type_,  new_array


    def retrieve_ix_val(
            array
            ):
        """ Retrieve value and index  and build `posMinMax boundaries

        :param array: array of main colum contains the anomaly definitions or
                a souding curve location like ::

                sloc = [NaN, 'low', NaN, NaN, NaN, 'ves', NaN,
                        NaN, 'up', NaN, NaN, NaN]
                `low`, `ves` and `up` are the lower boundary, the electric
                sounding  and the upper boundary of the selected anomaly
                respectively.
        For instance, if dipole_length is =`10`m, t he location (`pk`)
            of `low`, `ves` and `up` are 10, 50 and 80 m respectively.
            `posMinMax` =(10, 80)
        """

        lower_ix =None
        upper_ix =None
        ves_ix = None

        array= array.reshape((array.shape[0],) )
        for ix, val in enumerate(array):
            for low, up, vloc in zip(
                    ['lower', 'inf', 'min', 'min', '1', 'low'],
                    ['upper', 'sup', 'maj', 'max', '2', 'up'],
                    ['ves', 'se', 'sond','vs', 'loc', '0', 'dl']
                    ):
                try :
                    floatNaNor123= np.float(val)
                except:
                    if val.lower().find(low)>=0:
                        lower_ix = ix
                        break
                    elif val.lower().find(up) >=0:
                        upper_ix = ix
                        break
                    elif val.lower().find(vloc)>=0:
                        ves_ix = ix
                        break
                else :
                    if floatNaNor123 ==1:
                        lower_ix = ix
                        break
                    elif floatNaNor123 ==2:
                        upper_ix = ix
                        break
                    elif floatNaNor123 ==0:
                        ves_ix = ix
                        break

        return lower_ix, ves_ix, upper_ix

    # set pandas so to consider np.inf as NaN number.

    pd.options.mode.use_inf_as_na = True

    # unecesseray to specify the colum of sounding location.
    # dl =['drill', 'dl', 'loc', 'dh', 'choi']

    _autoOption=False  # set automatic to False one posMinMax
    # not found as well asthe anomaly location `ves`.
    posMinMax =None
    #get df columns from the 4-iem index
    for sl in ['pk', 'sta', 'loc']:
        for val in df.columns:
            if val.lower()==sl:
                pk_series = df[val].to_numpy()
                break

    listOfAddedColumns= df.iloc[:, 4:].columns

    if len(listOfAddedColumns) ==0:
        return True,  shape, type_, None, posMinMax,  df

    df_= df.iloc[:, 4:]
    # check whether all remains dataframe values are `NaN` values
    if len(list(df_.columns[df_.isna().all()])) == len(listOfAddedColumns):
         # If yes , trigger the auto option
        return True,  shape, type_, None, posMinMax,  df.iloc[:, :4]

    # get the colum name with any nan values
    sloc_column=list(df_.columns[df_.isna().any()])
    # if column contains  one np.nan, the sloc colum is found
    sloc_values = df_[sloc_column].to_numpy()

    if len(sloc_column)>1 : #
    # get the value from single array
        shape , type_,  sloc_values =  mergeToOne(sloc_column, df_)


    lower_ix, ves_ix ,upper_ix   = retrieve_ix_val(sloc_values)

    # if `lower` and `upper` bounds are not found then start or end limits of
    # selected anomaly  from the position(pk) of the sounding curve.
    if lower_ix is None :
        lower_ix =ves_ix
    if upper_ix is None:
        upper_ix = ves_ix

    if (lower_ix  and upper_ix ) is None:
        posMinMax =None
    if posMinMax is None and ves_ix is None: _autoOption =True
    else :
        posMinMax =(pk_series[lower_ix], pk_series[upper_ix] )

    if ves_ix is None: ves_loc=None
    else : ves_loc = pk_series[ves_ix]

    return _autoOption, shape, type_, ves_loc , posMinMax,  df.iloc[:, :4]


def reshape(arr , axis = None) :
    """ Detect the array shape and reshape it accordingly, back to the given axis.

    :param array: array_like with number of dimension equals to 1 or 2
    :param axis: axis to reshape back array. If 'axis' is None and
        the number of dimension is greater than 1, it reshapes back array
        to array-like

    :returns: New reshaped array

    :Example:
        >>> import numpy as np
        >>> from watex.utils.funcutils import reshape
        >>> array = np.random.randn(50 )
        >>> array.shape
        ... (50,)
        >>> ar1 = reshape(array, 1)
        >>> ar1.shape
        ... (1, 50)
        >>> ar2 =reshape(ar1 , 0)
        >>> ar2.shape
        ... (50, 1)
        >>> ar3 = reshape(ar2, axis = None)
        >>> ar3.shape # goes back to the original array
        >>> ar3.shape
        ... (50,)

    """
    arr = np.array(arr)
    if arr.ndim > 2 :
        raise ValueError('Expect an array with max dimension equals to 2'
                         f' but {str(arr.ndim)!r} were given.')

    if axis  not in (0 , 1, -1, None):
        raise ValueError(f'Wrong axis value: {str(axis)!r}')

    if axis ==-1:
        axis =None
    if arr.ndim ==1 :
        # ie , axis is None , array is an array-like object
        s0, s1= arr.shape [0], None
    else :
        s0, s1 = arr.shape
    if s1 is None:
        return  arr.reshape ((1, s0)) if axis == 1 else (arr.reshape (
            (s0, 1)) if axis ==0 else arr )
    try :
        arr = arr.reshape ((s0 if s1==1 else s1, )) if axis is None else (
            arr.reshape ((1, s0)) if axis==1  else arr.reshape ((s1, 1 ))
            )
    except ValueError:
        # error raises when user mistakes to input the right axis.
        # (ValueError: cannot reshape array of size 54 into shape (1,1))
        # then return to him the original array
        pass

    return arr


def ismissing(refarr, arr, fill_value = np.nan, return_index =False):
    """ Get the missing values in array-like and fill it  to match the length
    of the reference array.

    The function makes sense especially for frequency interpollation in the
    'attenuation band' when using the audio-frequency magnetotelluric methods.

    :param arr: array-like- Array to be extended with fill value. It should be
        shorter than the `refarr`. Otherwise it returns the same array `arr`
    :param refarr: array-like- the reference array. It should have a greater
        length than the array
    :param fill_value: float - Value to fill the `arr` to match the length of
        the `refarr`.
    :param return_index: bool or str - array-like, index of the elements element
        in `arr`. Default is ``False``. Any other value should returns the
        mask of existing element in reference array

    :returns: array and values missings or indexes in reference array.

    :Example:

    >>> import numpy as np
    >>> from watex.utils.funcutils import ismissing
    >>> refreq = np.linspace(7e7, 1e0, 20) # 20 frequencies as reference
    >>> # remove the value between index 7 to 12 and stack again
    >>> freq = np.hstack ((refreq.copy()[:7], refreq.copy()[12:] ))
    >>> f, m  = ismissing (refreq, freq)
    >>> f, m
    ...array([7.00000000e+07, 6.63157895e+07, 6.26315791e+07, 5.89473686e+07,
           5.52631581e+07, 5.15789476e+07, 4.78947372e+07,            nan,
                      nan,            nan,            nan,            nan,
           2.57894743e+07, 2.21052638e+07, 1.84210534e+07, 1.47368429e+07,
           1.10526324e+07, 7.36842195e+06, 3.68421147e+06, 1.00000000e+00])
    >>> m # missing values
    ... array([44210526.68421052, 40526316.21052632, 36842105.73684211,
           33157895.2631579 , 29473684.78947368])
    >>>  _, m_ix  = ismissing (refreq, freq, return_index =True)
    >>> m_ix
    ... array([ 7,  8,  9, 10, 11], dtype=int64)
    >>> # assert the missing values from reference values
    >>> refreq[m_ix ] # is equal to m
    ... array([44210526.68421052, 40526316.21052632, 36842105.73684211,
           33157895.2631579 , 29473684.78947368])

    """
    return_index = str(return_index).lower()
    fill_value = _assert_all_types(fill_value, float, int)
    if return_index in ('false', 'value', 'val') :
        return_index ='values'
    elif return_index  in ('true', 'index', 'ix') :
        return_index = 'index'
    else :
        return_index = 'mask'

    ref = refarr.copy() ; mask = np.isin(ref, arr)
    miss_values = ref [~np.isin(ref, arr)]
    miss_val_or_ix  = (ref [:, None] == miss_values).argmax(axis=0
                         ) if return_index =='index' else ref [~np.isin(ref, arr)]

    miss_val_or_ix = mask if return_index =='mask' else miss_val_or_ix
    # if return_missing_values:
    ref [~np.isin(ref, arr)] = fill_value
    #arr= np.hstack ((arr , np.repeat(fill_value, 0 if m <=0 else m  )))
    #refarr[refarr ==arr] if return_index else arr
    return  ref , miss_val_or_ix

def make_arr_consistent (
        refarr, arr, fill_value = np.nan, return_index = False,
        method='naive'):
    """
    Make `arr` to be consistent with the reference array `refarr`. Fill the
    missing value with param `fill_value`.

    Note that it does care of the position of the value in the array. Use
    Numpy digitize to compute the bins. The array caveat here is the bins
    must be monotonically decreasing or increasing.

    If the values in `arr` are present in `refarr`, the position of `arr`
    in new consistent array should be located decreasing or increasing order.

    Parameters
    ------------
    arr: array-like 1d,
        Array to extended with fill value. It should be  shorter than the
        `refarr`.

    refarr: array-like- the reference array. It should have a greater
        length than the array `arr`.
    fill_value: float,
        Value to fill the `arr` to match the length of the `refarr`.
    return_index: bool or str, default=True
         index of the position of the  elements in `refarr`.
         Default is ``False``. If ``mask`` should  return the
        mask of existing element in reference array
    method: str, default="naive"
        Is the method used to find the right position of items in `arr`
        based on the reference array.
        - ``naive``, considers the length of ``arr`` must fit the number of
            items that should be visible in the consistent array. This method
            erases the remaining bins values out of length of `arr`.
        - ``strict` did the same but rather than considering the length,
            it considers the maximum values in the `arr`. It assumes that `arr`
            is sorted in ascending order. This methods is usefull for plotting
            a specific stations since the station loactions are sorted in
            ascending order.

    Returns
    ---------
    non_zero_index , mask or t
        index: indices of the position of `arr` items in ``refarr``.
        mask: bool of the position `arr` items in ``refarr``
        t: new consistent array with the same length as ``refarr``

    Examples
    ----------
    >>> import numpy as np
    >>> from watex.utils.funcutils import make_arr_consistent
    >>> refarr = np.arange (12)
    >>> arr = np.arange (7, 10)
    >>> make_arr_consistent (refarr, arr )
    Out[84]: array([nan, nan, nan, nan, nan, nan, nan,  7.,  8.,  9., nan, nan])
    >>> make_arr_consistent (refarr, arr , return_index =True )
    Out[104]: array([7, 8, 9], dtype=int64)
    >>> make_arr_consistent (refarr, arr , return_index ="mask" )
    Out[105]:
    array([False, False, False, False, False, False, False,  True,  True,
            True, False, False])
    >>> a = np.arange ( 12 ); b = np.linspace (7, 10 , 7)
    >>> make_arr_consistent (a, b )
    Out[112]: array([nan, nan, nan, nan, nan, nan, nan,  7.,  8.,  9., 10., 11.])
    >>> make_arr_consistent (a, b ,method='strict')
    Out[114]: array([nan, nan, nan, nan, nan, nan, nan,  7.,  8.,  9., 10., nan])
    """
    try :
        refarr = reshape( refarr).shape[1]
        arr= reshape( arr).shape[1]
    except :pass
    else: raise TypeError ("Expects one-dimensional arrays for both arrays.")

    t = np.full_like( refarr, fill_value = np.nan, dtype =float )
    temp_arr = np.digitize( refarr, arr)
    non_zero_index = reshape (np.argwhere (temp_arr!=0 ) )
    t[non_zero_index] = refarr [non_zero_index]
    # force value to keep only
    # value in array
    if method=='strict':
        index = reshape ( np.argwhere (  (max( arr)  - t) < 0 ) )
        t [index ]= np.nan
    else:
        if len (t[~np.isnan (t)]) > len(arr):
            t [ - (len(t[~np.isnan (t)])-len(arr)):]= np.nan
    # update the non_zeros index
    non_zero_index= reshape ( np.argwhere (~np.isnan (t)))
    # now replace all NaN value by filled value
    t [np.isnan(t)] = fill_value

    return  refarr == t  if return_index =='mask' else (
        non_zero_index if return_index else t )

def find_close_position (refarr, arr):
    """ Get the close item from `arr` in the reference array `refarr`.

    :param arr: array-like 1d,
        Array to extended with fill value. It should be  shorter than the
        `refarr`.

    :param refarr: array-like-
        the reference array. It should have a greater length than the
        array `arr`.
    :return: generator of index of the closest position in  `refarr`.
    """
    for item in arr :
        ix = np.argmin (np.abs (refarr - item))
        yield ix


def fillNaN(arr, method ='ff'):
    """ Most efficient way to back/forward-fill NaN values in numpy array.

    Parameters
    ----------
    arr : ndarray
        Array containing NaN values to be filled
    method: str
        Method for filling. Can be forward fill ``ff`` or backward fill `bf``.
        or ``both`` for the two methods. Default is `ff`.

    Returns
    -------
    new array filled.

    Notes
    -----
    When NaN value is framed between two valid numbers, ``ff`` and `bf` performs
    well the filling operations. However, when the array is ended by multiple
    NaN values, the ``ff`` is recommended. At the opposite the ``bf`` is  the
    method suggested. The ``both``argument does the both tasks at the expense of
    the computation cost.

    Examples
    ---------

    >>> import numpy as np
    >>> from from watex.utils.funcutils import fillNaN
    >>> arr2d = np.random.randn(7, 3)
    >>> # change some value into NaN
    >>> arr2d[[0, 2, 3, 3 ],[0, 2,1, 2]]= np.nan
    >>> arr2d
    ... array([[        nan, -0.74636104,  1.12731613],
           [ 0.48178017, -0.18593812, -0.67673698],
           [ 0.17143421, -2.15184895,         nan],
           [-0.6839212 ,         nan,         nan]])
    >>> fillNaN (arr2d)
    ... array([[        nan, -0.74636104,  1.12731613],
           [ 0.48178017, -0.18593812, -0.67673698],
           [ 0.17143421, -2.15184895, -2.15184895],
           [-0.6839212 , -0.6839212 , -0.6839212 ]])
    >>> fillNaN(arr2d, 'bf')
    ... array([[-0.74636104, -0.74636104,  1.12731613],
           [ 0.48178017, -0.18593812, -0.67673698],
           [ 0.17143421, -2.15184895,         nan],
           [-0.6839212 ,         nan,         nan]])
    >>> fillNaN (arr2d, 'both')
    ... array([[-0.74636104, -0.74636104,  1.12731613],
           [ 0.48178017, -0.18593812, -0.67673698],
           [ 0.17143421, -2.15184895, -2.15184895],
           [-0.6839212 , -0.6839212 , -0.6839212 ]])

    References
    ----------
    Some function below are edited by the authors in pyQuestion.com website.
    There are other way more efficient to perform this task by calling the module
    `Numba` to accelerate the computation time. However, at the time this script
    is writen (August 17th, 2022) , `Numba` works with `Numpy` version 1.21. The
    latter  is older than the one used in for writting this package (1.22.3 ).

    For furher details, one can refer to the following link:
    https://pyquestions.com/most-efficient-way-to-forward-fill-nan-values-in-numpy-array

    """

    if not hasattr(arr, '__array__'):
        arr = np.array(arr)

    def ffill (arr):
        """ Forward fill."""
        idx = np.where (~mask, np.arange(mask.shape[1]), 0)
        np.maximum.accumulate (idx, axis =1 , out =idx )
        return arr[np.arange(idx.shape[0])[:, None], idx ]

    def bfill (arr):
        """ Backward fill """
        idx = np.where (~mask, np.arange(mask.shape[1]) , mask.shape[1]-1)
        idx = np.minimum.accumulate(idx[:, ::-1], axis =1)[:, ::-1]
        return arr [np.arange(idx.shape [0])[:, None], idx ]

    method= str(method).lower().strip()

    if arr.ndim ==1:
        arr = reshape(arr, axis=1)

    if method  in ('backward', 'bf',  'bwd'):
        method = 'bf'
    elif method in ('forward', 'ff', 'fwd'):
        method= 'ff'
    elif method in ('both', 'ffbf', 'fbwf', 'bff', 'full'):
        method ='both'
    if method not in ('bf', 'ff', 'both'):
        raise ValueError ("Expect a backward <'bf'>, forward <'ff'> fill "
                          f" or both <'bff'> not {method!r}")
    mask = np.isnan (arr )
    if method =='both':
        arr = ffill(arr)
        #mask = np.isnan (arr)
        arr = bfill(arr)

    return (ffill(arr) if method =='ff' else bfill(arr)
            ) if method in ('bf', 'ff') else arr


def get_params (obj: object
               ) -> Dict:
    """
    Get object parameters.

    Object can be callable or instances

    :param obj: object , can be callable or instance

    :return: dict of parameters values

    :examples:
        >>> from sklearn.svm import SVC
        >>> from watex.utils.funcutils import get_params
        >>> sigmoid= SVC (
            **{
                'C': 512.0,
                'coef0': 0,
                'degree': 1,
                'gamma': 0.001953125,
                'kernel': 'sigmoid',
                'tol': 1.0
                }
            )
        >>> pvalues = get_params( sigmoid)
        >>> {'decision_function_shape': 'ovr',
             'break_ties': False,
             'kernel': 'sigmoid',
             'degree': 1,
             'gamma': 0.001953125,
             'coef0': 0,
             'tol': 1.0,
             'C': 512.0,
             'nu': 0.0,
             'epsilon': 0.0,
             'shrinking': True,
             'probability': False,
             'cache_size': 200,
             'class_weight': None,
             'verbose': False,
             'max_iter': -1,
             'random_state': None
         }
    """
    if callable(obj):
        cls_or_func_signature = inspect.signature(obj)
        PARAMS_VALUES = {k: None if v.default is (inspect.Parameter.empty
                         or ...) else v.default
                    for k, v in cls_or_func_signature.parameters.items()
                    # if v.default is not inspect.Parameter.empty
                    }
    elif hasattr(obj, '__dict__'):
        PARAMS_VALUES = {k:v  for k, v in obj.__dict__.items()
                         if not (k.endswith('_') or k.startswith('_'))}

    return PARAMS_VALUES


def fit_ll(ediObjs, by ='index', method ='strict', distance='cartesian' ):
    """ Fit EDI by location and reorganize EDI according to the site
    longitude and latitude coordinates.

    EDIs data are mostly reading in an alphabetically order, so the reoganization

    according to the location(longitude and latitude) is usefull for distance
    betwen site computing with a right position at each site.

    :param ediObjs: list of EDI object, composed of a collection of
        watex.edi.Edi or pycsamt.core.edi.Edi or mtpy.core.edi objects
    :type ediObjs: watex.edi.Edi_Collection

    :param by: ['name'|'ll'|'distance'|'index'|'name'|'dataid']
       The kind to sorting EDI files. Default uses the position number
       included in the EDI-files name.
    :type by: str

    :param method:  ['strict|'naive']. Kind of method to sort the
        EDI file from longitude, latitude. Default is ``strict``.
    :type method: str

    :param distance: ['cartesian'|'harvesine']. Use the distance between
       coordinates points to sort EDI files. Default is ``cartesian`` distance.
    :type distance: str

    :returns: array splitted into ediObjs and Edifiles basenames
    :rtyple: tuple

    :Example:
        >>> import numpy as np
        >>> from watex.methods.em import EM
        >>> from watex.utils.funcutils import fit_ll
        >>> edipath ='data/edi_ss'
        >>> cediObjs = EM().fit (edipath)
        >>> ediObjs = np.random.permutation(cediObjs.ediObjs) # shuffle the
        ... # the collection of ediObjs
        >>> ediObjs, ediObjbname = fit_by_ll(ediObjs)
        ...

    """
    method= 'strict' if str(method).lower() =='strict' else "naive"
    if method=='strict':
        return _fit_ll(ediObjs, by = by, distance = distance )

    #get the ediObjs+ names in ndarray(len(ediObjs), 2)
    objnames = np.c_[ediObjs, np.array(
        list(map(lambda obj: os.path.basename(obj.edifile), ediObjs)))]
    lataddlon = np.array (list(map(lambda obj: obj.lat + obj.lon , ediObjs)))
    if len(np.unique ( lataddlon)) < len(ediObjs)//2:
        # then ignore reorganization and used the
        # station names.
        pass
    else:
        sort_ix = np.argsort(lataddlon)
        objnames = objnames[sort_ix ]

    #ediObjs , objbnames = np.hsplit(objnames, 2)
    return objnames[:, 0], objnames[:, -1]

def _fit_ll(ediObjs, distance='cartes', by = 'index'):
    """ Fit ediObjs using the `strict method`.

    An isolated part of :func:`watex.utils.funcutils.fit_by_ll`.
    """
    # get one obj randomnly and compute distance
    obj_init = ediObjs[0]
    ref_lat = 34.0522  # Latitude of Los Angeles
    ref_lon = -118.2437 # Longitude of Los Angeles

    if str(distance).find ('harves')>=0:
        distance='harves'
    else: distance='cartes'

    # create stations list.
    stations = [
        {"name": os.path.basename(obj.edifile),
         "longitude": obj.lon,
         "latitude": obj.lat,
         "obj": obj,
         "dataid": obj.dataid,
         # compute distance using cartesian or harversine
         "distance": _compute_haversine_d (
            ref_lat, ref_lon, obj.lat, obj.lon
            ) if distance =='harves' else np.sqrt (
                ( obj_init.lon -obj.lon)**2 + (obj_init.lat -obj.lat)**2),
         # check wether there is a position number in the data.
         "index": re.search (r'\d+', str(os.path.basename(obj.edifile)),
                            flags=re.IGNORECASE).group() if bool(
                                re.search(r'\d', os.path.basename(obj.edifile)))
                                else float(ii) ,
        }
        for ii, obj in enumerate (ediObjs)
        ]

    ll=( 'longitude', 'latitude')

    by = 'index' or str(by ).lower()
    if ( by.find ('ll')>=0 or by.find ('lonlat')>=0):
        by ='ll'
    elif  by.find ('latlon')>=0:
        ll =ll[::-1] # reverse

    # sorted from key
    sorted_stations = sorted (
        stations , key = lambda o: (o[ll[0]], [ll[-1]])
        if (by =='ll' or by=='latlon')
        else o[by]
             )

    objnames = np.array( list(
        map ( lambda o : o['name'], sorted_stations)))
    ediObjs = np.array ( list(
        map ( lambda o: o['obj'], sorted_stations)),
                        dtype =object )

    return ediObjs, objnames

def _compute_haversine_d(lat1, lon1, lat2, lon2):
    """ Sort coordinates using Haversine distance calculus.
    An isolated part of :func:`watex.utils.funcutils._fit_by_ll"""
    # get reference_lat and reference lon
    # get one obj randomnly and compute distance
    # obj_init = np.random.choice (ediObjs)
    import math
    # Define a function to calculate the distance
    # between two points in kilometers
    # def distance(lat1, lon1, lat2, lon2):
        # Convert degrees to radians
    lat1 = math.radians(lat1)
    lon1 = math.radians(lon1)
    lat2 = math.radians(lat2)
    lon2 = math.radians(lon2)

    # Apply the haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = math.sin(dlat / 2)**2 + math.cos(lat1) * math.cos(
        lat2) * math.sin(dlon / 2)**2
    c = 2 * math.asin(math.sqrt(a))
    r = 6371 # Earth's radius in kilometers

    return c * r


def make_ids(arr, prefix =None, how ='py', skip=False):
    """ Generate auto Id according to the number of given sites.

    :param arr: Iterable object to generate an id site . For instance it can be
        the array-like or list of EDI object that composed a collection of
        watex.edi.Edi object.
    :type ediObjs: array-like, list or tuple

    :param prefix: string value to add as prefix of given id. Prefix can be
        the site name.
    :type prefix: str

    :param how: Mode to index the station. Default is 'Python indexing' i.e.
        the counting starts by 0. Any other mode will start the counting by 1.
    :type cmode: str

    :param skip: skip the long formatage. the formatage acccording to the
        number of collected file.
    :type skip: bool
    :return: ID number formated
    :rtype: list

    :Example:
        >>> import numpy as np
        >>> from watex.utils.func_utils import make_ids
        >>> values = ['edi1', 'edi2', 'edi3']
        >>> make_ids (values, 'ix')
        ... ['ix0', 'ix1', 'ix2']
        >>> data = np.random.randn(20)
        >>>  make_ids (data, prefix ='line', how=None)
        ... ['line01','line02','line03', ... , line20]
        >>> make_ids (data, prefix ='line', how=None, skip =True)
        ... ['line1','line2','line3',..., line20]

    """
    fm='{:0' + ('1' if skip else f'{int(np.log10(len(arr))) + 1}') +'}'
    id_ =[str(prefix) + fm.format(i if how=='py'else i+ 1 ) if prefix is not
          None else fm.format(i if how=='py'else i+ 1)
          for i in range(len(arr))]
    return id_

def show_stats(nedic , nedir, fmtl='~', lenl=77, obj='EDI'):
    """ Estimate the file successfully read reading over the unread files

    :param nedic: number of input or collected files
    :param nedir: number of files read sucessfully
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
    print(mesg.format('Data collected','=',  nedic, f'{obj} success. read',
                      '=', nedir, 'Rate','=', round ((nedir/nedic) *100, 2),
                      2))
    print(fmtl * lenl )

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
    >>> from watex.utils.funcutils import concat_array_from_list
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

def station_id (id_, is_index= 'index', how=None, **kws):
    """
    From id get the station  name as input  and return index `id`.
    Index starts at 0.

    :param id_: str, of list of the name of the station or indexes .

    :param is_index: bool
        considered the given station as a index. so it remove all the letter and
        keep digit as index of each stations.

    :param how: Mode to index the station. Default is
        'Python indexing' i.e.the counting starts by 0. Any other mode will
        start the counting by 1. Note that if `is_index` is ``True`` and the
        param `how` is set to it default value ``py``, the station index should
        be downgraded to 1.

    :param kws: additionnal keywords arguments from :func:`~.make_ids`.

    :return: station index. If the list `id_` is given will return the tuple.

    :Example:

    >>> from watex.utils.funcutils import station_id
    >>> dat1 = ['S13', 's02', 's85', 'pk20', 'posix1256']
    >>> station_id (dat1)
    ... (13, 2, 85, 20, 1256)
    >>> station_id (dat1, how='py')
    ... (12, 1, 84, 19, 1255)
    >>> station_id (dat1, is_index= None, prefix ='site')
    ... ('site1', 'site2', 'site3', 'site4', 'site5')
    >>> dat2 = 1
    >>> station_id (dat2) # return index like it is
    ... 1
    >>> station_id (dat2, how='py') # considering the index starts from 0
    ... 0

    """
    is_iterable =False
    is_index = str(is_index).lower().strip()
    isix=True if  is_index in ('true', 'index', 'yes', 'ix') else False

    regex = re.compile(r'\d+', flags=re.IGNORECASE)
    try :
        iter (id_)
    except :
        id_= [id_]
    else : is_iterable=True

    #remove all the letter
    id_= list(map( lambda o: regex.findall(o), list(map(str, id_))))
    # merge the sequences list and for consistency remove emty list or str
    id_=tuple(filter (None, list(itertools.chain(*id_))))

    # if considering as Python index return value -1 other wise return index

    id_ = tuple (map(int, np.array(id_, dtype = np.int32)-1)
                 ) if how =='py' else tuple ( map(int, id_))

    if (np.array(id_) < 0).any():
        warnings.warn('Index contains negative values. Be aware that you are'
                      " using a Python indexing. Otherwise turn 'how' argumennt"
                      " to 'None'.", stacklevel=2)
    if not isix :
        id_= tuple(make_ids(id_, how= how,  **kws))

    if not is_iterable :
        try: id_ = id_[0]
        except : warnings.warn("The station id is given as a non iterable "
                          "object, but can keep the same format in return.", stacklevel=2)
        if id_==-1: id_= 0 if how=='py' else id_ + 2

    return id_

def assert_doi(doi):
    """
     assert the depth of investigation Depth of investigation converter

    :param doi: depth of investigation in meters.  If value is given as string
        following by yhe index suffix of kilometers 'km', value should be
        converted instead.
    :type doi: str|float

    :returns doi:value in meter
    :rtype: float

    """
    if isinstance (doi, str):
        if doi.find('km')>=0 :
            try: doi= float(doi.replace('km', '000'))
            except :TypeError (" Unrecognized value. Expect value in 'km' "
                           f"or 'm' not: {doi!r}")
    try: doi = float(doi)
    except: TypeError ("Depth of investigation must be a float number "
                       "not: {str(type(doi).__name__!r)}")
    return doi

def strip_item(item_to_clean, item=None, multi_space=12):
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
    if item is None :
        item = ' '

    cleaner =[(''+ ii*f'{item}') for ii in range(multi_space)]

    if isinstance (item_to_clean, str) :
        item_to_clean=[item_to_clean]

    # if type(item_to_clean ) != list :#or type(item_to_clean ) !=np.ndarray:
    #     if type(item_to_clean ) !=np.ndarray:
    #         item_to_clean=[item_to_clean]
    if item_to_clean in cleaner or item_to_clean ==['']:
        #warnings.warn ('No data found for sanitization; returns None.')
        return None
    try :
        multi_space=int(multi_space)
    except :
        raise TypeError('argument <multplier> must be an integer'
                        f'not {type(multi_space)}')

    for jj, ss in enumerate(item_to_clean) :
        for space in cleaner:
            if space in ss :
                new_ss=ss.strip(space)
                item_to_clean[jj]=new_ss

    return item_to_clean

def parse_json(json_fn =None,
               data=None,
               todo='load',
               savepath=None,
               verbose:int =0,
               **jsonkws):
    """ Parse Java Script Object Notation file and collect data from JSON
    config file.

    :param json_fn: Json filename, URL or output JSON name if `data` is
        given and `todo` is set to ``dump``.Otherwise the JSON output filename
        should be the `data` or the given variable name.
    :param data: Data in Python obj to serialize.
    :param todo: Action to perform with JSON:
        - load: Load data from the JSON file
        - dump: serialize data from the Python object and create a JSON file
    :param savepath: If ``default``  should save the `json_fn`
        If path does not exist, should save to the <'_savejson_'>
        default path .
    :param verbose: int, control the verbosity. Output messages

    .. see also:: Read more about JSON doc
            https://docs.python.org/3/library/json.html
         or https://www.w3schools.com/python/python_json.asp
         or https://www.geeksforgeeks.org/json-load-in-python/
         ...

    :Example:
        >>> PATH = 'data/model'
        >>> k_ =['model', 'iter', 'mesh', 'data']
        >>> try :
            INVERS_KWS = {
                s +'_fn':os.path.join(PATH, file)
                for file in os.listdir(PATH)
                          for s in k_ if file.lower().find(s)>=0
                          }
        except :
            INVERS=dict()
        >>> TRES=[10, 66,  70, 100, 1000, 3000]# 7000]
        >>> LNS =['river water','fracture zone', 'MWG', 'LWG',
              'granite', 'igneous rocks', 'basement rocks']
        >>> import watex.utils.funcutils as FU
        >>> geo_kws ={'oc2d': INVERS_KWS,
                      'TRES':TRES, 'LN':LNS}
        # serialize json data and save to  'jsontest.json' file
        >>> FU.parse_json(json_fn = 'jsontest.json',
                          data=geo_kws, todo='dump', indent=3,
                          savepath ='data/saveJSON', sort_keys=True)
        # Load data from 'jsontest.json' file.
        >>> FU.parse_json(json_fn='data/saveJSON/jsontest.json', todo ='load')

    """
    todo, domsg =return_ctask(todo)
    # read urls by default json_fn can hold a url
    try :
        if json_fn.find('http') >=0 :
            todo, json_fn, data = fetch_json_data_from_url(json_fn, todo)
    except:
        #'NoneType' object has no attribute 'find' if data is not given
        pass

    if todo.find('dump')>=0:
        json_fn = get_config_fname_from_varname(
            data, config_fname= json_fn, config='.json')

    JSON = dict(load=json.load,# use loads rather than load
                loads=json.loads,
                dump= json.dump,
                dumps= json.dumps)
    try :
        if todo=='load': # read JSON files
            with open(json_fn) as fj:
                data =  JSON[todo](fj)
        elif todo=='loads': # can be JSON string format
            data = JSON[todo](json_fn)
        elif todo =='dump': # store data in JSON file.
            with open(f'{json_fn}.json', 'w') as fw:
                data = JSON[todo](data, fw, **jsonkws)
        elif todo=='dumps': # store data in JSON format not output file.
            data = JSON[todo](data, **jsonkws)

    except json.JSONDecodeError:
        raise json.JSONDecodeError(f"Unable {domsg} JSON {json_fn!r} file. "
                              "Please check your file.", f'{json_fn!r}', 1)
    except:
        msg =''.join([
        f"{'Unrecognizable file' if todo.find('load')>=0 else'Unable to serialize'}"
        ])

        raise TypeError(f'{msg} {json_fn!r}. Please check your'
                        f" {'file' if todo.find('load')>=0 else 'data'}.")

    cparser_manager(f'{json_fn}.json',savepath, todo=todo, dpath='_savejson_',
                    verbose=verbose , config='JSON' )

    return data

def fetch_json_data_from_url (url:str , todo:str ='load'):
    """ Retrieve JSON data from url
    :param url: Universal Resource Locator .
    :param todo:  Action to perform with JSON:
        - load: Load data from the JSON file
        - dump: serialize data from the Python object and create a JSON file
    """
    with urllib.request.urlopen(url) as jresponse :
        source = jresponse.read()
    data = json.loads(source)
    if todo .find('load')>=0:
        todo , json_fn  ='loads', source

    if todo.find('dump')>=0:  # then collect the data and dump it
        # set json default filename
        todo, json_fn = 'dumps',  '_urlsourcejsonf.json'

    return todo, json_fn, data

def parse_csv(
        csv_fn:str =None,
        data=None,
        todo='reader',
        fieldnames=None,
        savepath=None,
        header: bool=False,
        verbose:int=0,
        **csvkws
   ) :
    """ Parse comma separated file or collect data from CSV.

    :param csv_fn: csv filename,or output CSV name if `data` is
        given and `todo` is set to ``write|dictwriter``.Otherwise the CSV
        output filename should be the `c.data` or the given variable name.
    :param data: Sequence Data in Python obj to write.
    :param todo: Action to perform with JSON:
        - reader|DictReader: Load data from the JSON file
        - writer|DictWriter: Write data from the Python object
        and create a CSV file
    :param savepath: If ``default``  should save the `csv_fn`
        If path does not exist, should save to the <'_savecsv_'>
        default path.
    :param fieldnames: is a sequence of keys that identify the order
        in which values in the dictionary passed to the `writerow()`
            method are written `csv_fn` file.
    :param savepath: If ``default``  should save the `csv_fn`
        If path does not exist, should save to the <'_savecsv_'>
        default path .
    :param verbose: int, control the verbosity. Output messages
    :param csvkws: additional keywords csv class arguments

    .. see also:: Read more about CSV module in:
        https://docs.python.org/3/library/csv.html or find some examples
        here https://www.programcreek.com/python/example/3190/csv.DictWriter
        or find some FAQS here:
    https://stackoverflow.com/questions/10373247/how-do-i-write-a-python-dictionary-to-a-csv-file
        ...
    :Example:
        >>> import watex.utils.funcutils as FU
        >>> PATH = 'data/model'
        >>> k_ =['model', 'iter', 'mesh', 'data']
        >>> try :
            INVERS_KWS = {
                s +'_fn':os.path.join(PATH, file)
                for file in os.listdir(PATH)
                          for s in k_ if file.lower().find(s)>=0
                          }
        except :
            INVERS=dict()
        >>> TRES=[10, 66,  70, 100, 1000, 3000]# 7000]
        >>> LNS =['river water','fracture zone', 'MWG', 'LWG',
              'granite', 'igneous rocks', 'basement rocks']
        >>> geo_kws ={'oc2d': INVERS_KWS,
                      'TRES':TRES, 'LN':LNS}
        >>> # write data and save to  'csvtest.csv' file
        >>> # here the `data` is a sequence of dictionary geo_kws
        >>> FU.parse_csv(csv_fn = 'csvtest.csv',data = [geo_kws],
                         fieldnames = geo_kws.keys(),todo= 'dictwriter',
                         savepath = 'data/saveCSV')
        # collect csv data from the 'csvtest.csv' file
        >>> FU.parse_csv(csv_fn ='data/saveCSV/csvtest.csv',
                         todo='dictreader',fieldnames = geo_kws.keys()
                         )

    """
    todo, domsg =return_ctask(todo)

    if todo.find('write')>=0:
        csv_fn = get_config_fname_from_varname(
            data, config_fname= csv_fn, config='.csv')
    try :
        if todo =='reader':
            with open (csv_fn) as csv_f :
                csv_reader = csv.reader(csv_f) # iterator
                data =[ row for row in csv_reader]

        elif todo=='writer':
            # write without a blank line, --> new_line =''
            with open(f'{csv_fn}.csv', 'w', newline ='',
                      encoding ='utf8') as new_csvf:
                csv_writer = csv.writer(new_csvf, **csvkws)
                csv_writer.writerows(data) if len(
                    data ) > 1 else csv_writer.writerow(data)
                # for row in data:
                #     csv_writer.writerow(row)
        elif todo=='dictreader':
            with open (csv_fn, encoding ='utf8') as csv_f :
                # generate an iterator obj
                csv_reader= csv.DictReader (csv_f, fieldnames= fieldnames)
                # return csvobj as a list of dicts
                data = list(csv_reader)

        elif todo=='dictwriter':
            with open(f'{csv_fn}.csv', 'w') as new_csvf:
                csv_writer = csv.DictWriter(new_csvf, **csvkws)
                if header:
                    csv_writer.writeheader()
                # DictWriter.writerows()expect a list of dicts,
                # while DictWriter.writerow() expect a single row of dict.
                csv_writer.writerow(data) if isinstance(
                    data , dict) else csv_writer.writerows(data)

    except csv.Error:
        raise csv.Error(f"Unable {domsg} CSV {csv_fn!r} file. "
                      "Please check your file.")
    except:

        msg =''.join([
        f"{'Unrecognizable file' if todo.find('read')>=0 else'Unable to write'}"
        ])

        raise TypeError(f'{msg} {csv_fn!r}. Please check your'
                        f" {'file' if todo.find('read')>=0 else 'data'}.")
    cparser_manager(f'{csv_fn}.csv',savepath, todo=todo, dpath='_savecsv_',
                    verbose=verbose , config='CSV' )

    return data

def return_ctask (todo:Optional[str]=None) -> Tuple [str, str]:
    """ Get the convenient task to do if users misinput the `todo` action.

    :param todo: Action to perform:
        - load: Load data from the config [YAML|CSV|JSON] file
        - dump: serialize data from the Python object and
            create a config [YAML|CSV|JSON] file."""

    def p_csv(v, cond='dict', base='reader'):
        """ Read csv instead.
        :param v: str, value to do
        :param cond: str, condition if  found in the value `v`.
        :param base: str, base task to do if condition `cond` is not met.

        :Example:

        >>> todo = 'readingbook'
        >>> p_csv(todo) <=> 'dictreader' if todo.find('dict')>=0 else 'reader'
        """
        return  f'{cond}{base}' if v.find(cond) >=0 else base

    ltags = ('load', 'recover', True, 'fetch')
    dtags = ('serialized', 'dump', 'save', 'write','serialize')
    if todo is None:
        raise ValueError('NoneType action can not be perform. Please '
                         'specify your action: `load` or `dump`?' )

    todo =str(todo).lower()
    ltags = list(ltags) + [todo] if  todo=='loads' else ltags
    dtags= list(dtags) +[todo] if  todo=='dumps' else dtags

    if todo in ltags:
        todo = 'loads' if todo=='loads' else 'load'
        domsg= 'to parse'
    elif todo in dtags:
        todo = 'dumps' if todo=='dumps' else 'dump'
        domsg  ='to serialize'
    elif todo.find('read')>=0:
        todo = p_csv(todo)
        domsg= 'to read'
    elif todo.find('write')>=0:
        todo = p_csv(todo, base ='writer')
        domsg =' to write'

    else :
        raise ValueError(f'Wrong action {todo!r}. Please select'
                         f' the right action to perform: `load` or `dump`?'
                        ' for [JSON|YAML] and `read` or `write`? '
                        'for [CSV].')
    return todo, domsg

def parse_yaml (yml_fn:str =None, data=None,
                todo='load', savepath=None,
                verbose:int =0, **ymlkws) :
    """ Parse yml file and collect data from YAML config file.

    :param yml_fn: yaml filename and can be the output YAML name if `data` is
        given and `todo` is set to ``dump``.Otherwise the YAML output filename
        should be the `data` or the given variable name.
    :param data: Data in Python obj to serialize.
    :param todo: Action to perform with YAML:
        - load: Load data from the YAML file
        - dump: serialize data from the Python object and create a YAML file
    :param savepath: If ``default``  should save the `yml_fn`
        to the default path otherwise should store to the convenient path.
        If path does not exist, should set to the default path.
    :param verbose: int, control the verbosity. Output messages

    .. see also:: Read more about YAML file https://pynative.com/python-yaml/
         or https://python.land/data-processing/python-yaml and download YAML
         at https://pypi.org/project/PyYAML/
         ...

    """

    todo, domsg =return_ctask(todo)
    #in the case user use dumps or loads with 's'at the end
    if todo.find('dump')>= 0:
        todo='dump'
    if todo.find('load')>=0:
        todo='load'
    if todo=='dump':
        yml_fn = get_config_fname_from_varname(data, yml_fn)
    try :
        if todo=='load':
            with open(yml_fn) as fy:
                data =  yaml.load(fy, Loader=yaml.SafeLoader)
                # args =yaml.safe_load(fy)
        elif todo =='dump':

            with open(f'{yml_fn}.yml', 'w') as fw:
                data = yaml.dump(data, fw, **ymlkws)
    except yaml.YAMLError:
        raise yaml.YAMLError(f"Unable {domsg} YAML {yml_fn!r} file. "
                             'Please check your file.')
    except:
        msg =''.join([
        f"{'Unrecognizable file' if todo=='load'else'Unable to serialize'}"
        ])

        raise TypeError(f'{msg} {yml_fn!r}. Please check your'
                        f" {'file' if todo=='load' else 'data'}.")

    cparser_manager(f'{yml_fn}.yml',savepath, todo=todo, dpath='_saveyaml_',
                    verbose=verbose , config='YAML' )

    return data

def cparser_manager (
    cfile,
    savepath =None,
    todo:str ='load',
    dpath=None,
    verbose =0,
    **pkws):
    """ Save and output message according to the action.

    :param cfile: name of the configuration file
    :param savepath: Path-like object
    :param dpath: default path
    :param todo: Action to perform with config file. Can ve
        ``load`` or ``dump``
    :param config: Type of configuration file. Can be ['YAML|CSV|JSON]
    :param verbose: int, control the verbosity. Output messages

    """
    if savepath is not None:
        if savepath =='default':
            savepath = None
        yml_fn,_= move_cfile(cfile, savepath, dpath=dpath)
    if verbose > 0:
        print_cmsg(yml_fn, todo, **pkws)


def get_config_fname_from_varname(data,
                                  config_fname=None,
                                  config='.yml') -> str:
    """ use the variable name given to data as the config file name.

    :param data: Given data to retrieve the variable name
    :param config_fname: Configurate variable filename. If ``None`` , use
        the name of the given varibale data
    :param config: Type of file for configuration. Can be ``json``, ``yml``
        or ``csv`` file. default is ``yml``.
    :return: str, the configuration data.

    """
    try:
        if '.' in config:
            config =config.replace('.','')
    except:pass # in the case None is given

    if config_fname is None: # get the varname
        # try :
        #     from varname.helpers import Wrapper
        # except ImportError:
        #     import_varname=False
        #     import_varname = FU.subprocess_module_installation('varname')
        #     if import_varname:
        #         from varname.helpers import Wrapper
        # else : import_varname=True
        try :
            for c, n in zip(['yml', 'yaml', 'json', 'csv'],
                            ['cy.data', 'cy.data', 'cj.data',
                             'c.data']):
                if config ==c:
                    config_fname= n
                    break
            if config_fname is None:
                raise # and go to except
        except:
            #using fstring
            config_fname= f'{data}'.split('=')[0]

    elif config_fname is not None:
        config_fname= config_fname.replace(
            f'.{config}', '').replace(f'.{config}', '').replace('.yaml', '')

    return config_fname

def pretty_printer(
        clfs: List[F],
        clf_score:List[float]=None,
        scoring: Optional[str] =None,
        **kws
 )->None:
    """ Format and pretty print messages after gridSearch using multiples
    estimators.

    Display for each estimator, its name, it best params with higher score
    and the mean scores.

    Parameters
    ----------
    clfs:Callables
        classifiers or estimators

    clf_scores: array-like
        for single classifier, usefull to provided the
        cross validation score.

    scoring: str
        Scoring used for grid search.
    """
    empty =kws.pop('empty', ' ')
    e_pad =kws.pop('e_pad', 2)
    p=list()

    if not isinstance(clfs, (list,tuple)):
        clfs =(clfs, clf_score)

    for ii, (clf, clf_be, clf_bp, clf_sc) in enumerate(clfs):
        s_=[e_pad* empty + f'{clf.__class__.__name__:<20}:' + '{:<20}:'.format(
                f'Best-estimator <{ii+1}>') +f'{clf_be}',
         e_pad* empty +'{:<20}:'.format(' ')+ '{:<20}:'.format(
            'Best paramaters') + f'{clf_bp}',
         e_pad* empty  +'{:<20}:'.format(' ') + '{:<20}:'.format(
            f'scores<`{scoring}`>') +f'{clf_sc}']
        try :
            s0= [e_pad* empty +'{:<20}:'.format(' ')+ '{:<20}:'.format(
            'scores mean')+ f'{clf_sc.mean()}']
        except AttributeError:
            s0= [e_pad* empty +'{:<20}:'.format(' ')+ '{:<20}:'.format(
            'scores mean')+ 'None']
            s_ +=s0
        else :
            s_ +=s0

        p.extend(s_)

    for i in p:
        print(i)

def move_cfile (
    cfile:str ,
    savepath:Optional[str]=None,
    **ckws
    ) :
    """ Move file to its savepath and output message.

    If path does not exist, should create one to save data.
    :param cfile: name of the configuration file
    :param savepath: Path-like object
    :param dpath: default path

    :returns:
        - configuration file
        - out message
    """
    savepath = cpath(savepath, **ckws)
    try :shutil.move(cfile, savepath)
    except:
        shutil.copy2 (cfile, savepath )
        os.remove (cfile)
        # delete file
        # warnings.warn("It seems the path already exists!")
    cfile = os.path.join(savepath, cfile)

    msg = (f'--> {os.path.basename(cfile)!r} data was successfully'
          f' saved to {os.path.realpath(cfile)!r}.')

    return cfile, msg

def print_cmsg(cfile:str, todo:str='load', config:str='YAML') -> str:
    """ Output configuration message.

    :param cfile: name of the configuration file
    :param todo: Action to perform with config file. Can be
        ``load`` or ``dump``
    :param config: Type of configuration file. Can be [YAML|CSV|JSON]
    """
    if todo=='load':
        msg = ''.join([
        f'--> Data was successfully stored to {os.path.basename(cfile)!r}',
            f' and saved to {os.path.realpath(cfile)!r}.']
            )
    elif todo=='dump':
        msg =''.join([ f"--> {config.upper()} {os.path.basename(cfile)!r}",
                      " data was sucessfully loaded."])
    return msg


def random_state_validator(seed):
    """Turn seed into a np.random.RandomState instance.

    Parameters
    ----------
    seed : None, int or instance of RandomState
        If seed is None, return the RandomState singleton used by np.random.
        If seed is an int, return a new RandomState instance seeded with seed.
        If seed is already a RandomState instance, return it.
        Otherwise raise ValueError.

    Returns
    -------
    :class:`numpy:numpy.random.RandomState`
        The random state object based on `seed` parameter.
    """
    if seed is None or seed is np.random:
        return np.random.mtrand._rand
    if isinstance(seed, numbers.Integral):
        return np.random.RandomState(seed)
    if isinstance(seed, np.random.RandomState):
        return seed
    raise ValueError(
        f"{seed!r} cannot be used to seed a numpy.random.RandomState instance"
    )

def is_iterable (
        y, /, exclude_string= False, transform = False , parse_string =False,
)->bool | list:
    r""" Asserts iterable object and returns boolean or transform object into
     an iterable.

    Function can also transform a non-iterable object to an iterable if
    `transform` is set to ``True``.

    :param y: any, object to be asserted
    :param exclude_string: bool, does not consider string as an iterable
        object if `y` is passed as a string object.
    :param transform: bool, transform  `y` to an iterable objects. But default
        puts `y` in a list object.
    :param parse_string: bool, parse string and convert the list of string
        into iterable object is the `y` is a string object and containg the
        word separator character '[#&.*@!_,;\s-]'. Refer to the function
        :func:`~watex.utils.funcutils.str2columns` documentation.

    :returns:
        - bool, or iterable object if `transform` is set to ``True``.

    .. note::
        Parameter `parse_string` expects `transform` to be ``True``, otherwise
        a ValueError will raise. Note :func:`.is_iterable` is not dedicated
        for string parsing. It parses string using the default behaviour of
        :func:`.str2columns`. Use the latter for string parsing instead.

    :Examples:
    >>> from watex.funcutils.is_iterable
    >>> is_iterable ('iterable', exclude_string= True )
    Out[28]: False
    >>> is_iterable ('iterable', exclude_string= True , transform =True)
    Out[29]: ['iterable']
    >>> is_iterable ('iterable', transform =True)
    Out[30]: 'iterable'
    >>> is_iterable ('iterable', transform =True, parse_string=True)
    Out[31]: ['iterable']
    >>> is_iterable ('iterable', transform =True, exclude_string =True,
                     parse_string=True)
    Out[32]: ['iterable']
    >>> is_iterable ('parse iterable object', parse_string=True,
                     transform =True)
    Out[40]: ['parse', 'iterable', 'object']
    """
    if (parse_string and not transform) and isinstance (y, str):
        raise ValueError ("Cannot parse the given string. Set 'transform' to"
                          " ``True`` otherwise use the 'str2columns' utils"
                          " from 'watex.utils.funcutils' instead.")
    y = str2columns(y) if isinstance(y, str) and parse_string else y

    isiter = False  if exclude_string and isinstance (
        y, str) else hasattr (y, '__iter__')

    return ( y if isiter else [ y ] )  if transform else isiter


def str2columns (text, /, regex=None , pattern = None):
    r"""Split text from the non-alphanumeric markers using regular expression.

    Remove all string non-alphanumeric and some operator indicators,  and
    fetch attributes names.

    Parameters
    -----------
    text: str,
        text litteral containing the columns the names to retrieve

    regex: `re` object,
        Regular expresion object. the default is::

            >>> import re
            >>> re.compile (r'[#&*@!_,;\s-]\s*', flags=re.IGNORECASE)
    pattern: str, default = '[#&*@!_,;\s-]\s*'
        The base pattern to split the text into a columns

    Returns
    -------
    attr: List of attributes

    Examples
    ---------
    >>> from watex.utils.funcutils import str2columns
    >>> text = ('this.is the text to split. It is an: example of; splitting str - to text.')
    >>> str2columns (text )
    ... ['this',
         'is',
         'the',
         'text',
         'to',
         'split',
         'It',
         'is',
         'an:',
         'example',
         'of',
         'splitting',
         'str',
         'to',
         'text']

    """
    pattern = pattern or  r'[#&.*@!_,;\s-]\s*'
    regex = regex or re.compile (pattern, flags=re.IGNORECASE)
    text= list(filter (None, regex.split(str(text))))
    return text

def sanitize_frame_cols(
        d, /, func:F = None , regex=None, pattern:str = None,
        fill_pattern:str =None, inplace:bool =False
        ):
    r""" Remove an indesirable characters and returns new columns

    Use regular expression for columns sanitizing

    Parameters
    -----------

    d: list, columns,
        columns to sanitize. It might contain a list of items to
        to polish. If dataframe or series are given, the dataframe columns
        and the name respectively will be polished and returns the same
        dataframe.

    func: F, callable
       Universal function used to clean the columns

    regex: `re` object,
        Regular expresion object. the default is::

            >>> import re
            >>> re.compile (r'[_#&.)(*@!_,;\s-]\s*', flags=re.IGNORECASE)
    pattern: str, default = '[_#&.)(*@!_,;\s-]\s*'
        The base pattern to sanitize the text in each column names.

    fill_pattern: str, default=''
        pattern to replace the non-alphabetic character in each item of
        columns.
    inplace: bool, default=False,
        transform the dataframe of series in place.

    Returns
    -------
    columns | pd.Series | dataframe.
        return Serie or dataframe if one is given, otherwise it returns a
        sanitized columns.

    Examples
    ---------
    >>> from watex.utils.funcutils import sanitize_frame_cols
    >>> from watex.utils.coreutils import read_data
    >>> h502= read_data ('data/boreholes/H502.xlsx')
    >>> h502 = sanitize_frame_cols (h502, fill_pattern ='_' )
    >>> h502.columns[:3]
    ... Index(['depth_top', 'depth_bottom', 'strata_name'], dtype='object')
    >>> f = lambda r : r.replace ('_', "'s ")
    >>> h502_f= sanitize_frame_cols( h502, func =f )
    >>> h502_f.columns [:3]
    ... Index(['depth's top', 'depth's bottom', 'strata's name'], dtype='object')

    """
    isf , iss= False , False
    pattern = pattern or r'[_#&.)(*@!_,;\s-]\s*'
    fill_pattern = fill_pattern or ''
    fill_pattern = str(fill_pattern)

    regex = regex or re.compile (pattern, flags=re.IGNORECASE)

    if isinstance(d, pd.Series):
        c = [d.name]
        iss =True
    elif isinstance (d, pd.DataFrame ) :
        c = list(d.columns)
        isf = True

    else :
        if not is_iterable(d) : c = [d]
        else : c = d

    if inspect.isfunction(func):
        c = list( map (func , c ) )

    else : c =list(map (
        lambda r : regex.sub(fill_pattern, r.strip() ), c ))

    if isf :
        if inplace : d.columns = c
        else : d =pd.DataFrame(d.values, columns =c )

    elif iss:
        if inplace: d.name = c[0]
        else : d= pd.Series (data =d.values, name =c[0] )

    else : d = c

    return d

def to_hdf5(d, /, fn, objname =None, close =True,  **hdf5_kws):
    """
    Store a frame data in hierachical data format 5 (HDF5)

    Note that is `d` is a dataframe, make sure that the dependency 'pytables'
    is already installed, otherwise and error raises.

    Parameters
    -----------
    d: ndarray,
        data to store in HDF5 format
    fn: str,
        File path to HDF5 file.
    objname: str,
        name of the data to store
    close: bool, default =True
        when data is given as an array, data can still be added if
        close is set to ``False``, otherwise, users need to open again in
        read mode 'r' before pursuing the process of adding.
    hdf5_kws: dict of :class:`pandas.pd.HDFStore`
        Additional keywords arguments passed to pd.HDFStore. they could be:
        *  mode : {'a', 'w', 'r', 'r+'}, default 'a'

             ``'r'``
                 Read-only; no data can be modified.
             ``'w'``
                 Write; a new file is created (an existing file with the same
                 name would be deleted).
             ``'a'``
                 Append; an existing file is opened for reading and writing,
                 and if the file does not exist it is created.
             ``'r+'``
                 It is similar to ``'a'``, but the file must already exist.
         * complevel : int, 0-9, default None
             Specifies a compression level for data.
             A value of 0 or None disables compression.
         * complib : {'zlib', 'lzo', 'bzip2', 'blosc'}, default 'zlib'
             Specifies the compression library to be used.
             As of v0.20.2 these additional compressors for Blosc are supported
             (default if no compressor specified: 'blosc:blosclz'):
             {'blosc:blosclz', 'blosc:lz4', 'blosc:lz4hc', 'blosc:snappy',
              'blosc:zlib', 'blosc:zstd'}.
             Specifying a compression library which is not available issues
             a ValueError.
         * fletcher32 : bool, default False
             If applying compression use the fletcher32 checksum.
    Returns
    -------
    store : Dict-like IO interface for storing pandas objects.

    Examples
    ------------
    >>> import os
    >>> from watex.utils.funcutils import sanitize_frame_cols, to_hdf5
    >>> from watex.utils import read_data
    >>> data = read_data('data/boreholes/H502.xlsx')
    >>> sanitize_frame_cols (data, fill_pattern='_', inplace =True )
    >>> store_path = os.path.join('watex/datasets/data', 'h') # 'h' is the name of the data
    >>> store = to_hdf5 (data, fn =store_path , objname ='h502' )
    >>> store
    ...
    >>> # fetch the data
    >>> h502 = store ['h502']
    >>> h502.columns[:3]
    ... Index(['hole_number', 'depth_top', 'depth_bottom'], dtype='object')

    """
    store =None
    if (
        not hasattr (d, '__array__')
        or not hasattr (d, 'columns')
            ) :
        raise TypeError ("Expect an array or dataframe,"
                         f" not {type (d).__name__!r}")

    if hasattr (d, '__array__') and hasattr (d, "columns"):
        # assert whether pytables is installed
        import_optional_dependency ('tables')
        # remove extension if exist.
        fn = str(fn).replace ('.h5', "").replace(".hdf5", "")
        # then store.
        store = pd.HDFStore(fn +'.h5' ,  **hdf5_kws)
        objname = objname or 'data'
        store[ str(objname) ] = d


    elif not hasattr(d, '__array__'):
        d = np.asarray(d)

        store= h5py.File(f"{fn}.hdf5", "w")
        store.create_dataset("dataset_01", store.shape,
                             dtype=store.dtype,
                             data=store
                             )

    if close : store.close ()

    return store


def find_by_regex (o , /, pattern,  func = re.match, **kws ):
    r""" Find pattern in object whatever an "iterable" or not.

    when we talk about iterable, a string value is not included.

    Parameters
    -------------
    o: str or iterable,
        text litteral or an iterable object containing or not the specific
        object to match.
    pattern: str, default = '[_#&*@!_,;\s-]\s*'
        The base pattern to split the text into a columns

    func: re callable , default=re.match
        regular expression search function. Can be
        [re.match, re.findall, re.search ],or any other regular expression
        function.

        * ``re.match()``:  function  searches the regular expression pattern and
            return the first occurrence. The Python RegEx Match method checks
            for a match only at the beginning of the string. So, if a match is
            found in the first line, it returns the match object. But if a match
            is found in some other line, the Python RegEx Match function returns
            null.
        * ``re.search()``: function will search the regular expression pattern
            and return the first occurrence. Unlike Python re.match(), it will
            check all lines of the input string. The Python re.search() function
            returns a match object when the pattern is found and “null” if
            the pattern is not found
        * ``re.findall()`` module is used to search for 'all' occurrences that
            match a given pattern. In contrast, search() module will only
            return the first occurrence that matches the specified pattern.
            findall() will iterate over all the lines of the file and will
            return all non-overlapping matches of pattern in a single step.
    kws: dict,
        Additional keywords arguments passed to functions :func:`re.match` or
        :func:`re.search` or :func:`re.findall`.

    Returns
    -------
    om: list
        matched object put is the list

    Example
    --------
    >>> from watex.utils.funcutils import find_by_regex
    >>> from watex.datasets import load_hlogs
    >>> X0, _= load_hlogs (as_frame =True )
    >>> columns = X0.columns
    >>> str_columns =','.join (columns)
    >>> find_by_regex (str_columns , pattern='depth', func=re.search)
    ... ['depth']
    >>> find_by_regex(columns, pattern ='depth', func=re.search)
    ... ['depth_top', 'depth_bottom']

    """
    om = []
    if isinstance (o, str):
        om = func ( pattern=pattern , string = o, **kws)
        if om:
            om= om.group()
        om =[om]
    elif is_iterable(o):
        o = list(o)
        for s in o :
            z = func (pattern =pattern , string = s, **kws)
            if z :
                om.append (s)

    if func.__name__=='findall':
        om = list(itertools.chain (*om ))
    # keep None is nothing
    # fit the corresponding pattern
    if len(om) ==0 or om[0] is None:
        om = None
    return  om

def is_in_if (o: iter, /, items: str | iter, error = 'raise',
               return_diff =False, return_intersect = False):
    """ Raise error if item is not  found in the iterable object 'o'

    :param o: unhashable type, iterable object,
        object for checkin. It assumes to be an iterable from which 'items'
        is premused to be in.
    :param items: str or list,
        Items to assert whether it is in `o` or not.
    :param error: str, default='raise'
        raise or ignore error when none item is found in `o`.
    :param return_diff: bool,
        returns the difference items which is/are not included in 'items'
        if `return_diff` is ``True``, will put error to ``ignore``
        systematically.
    :param return_intersect:bool,default=False
        returns items as the intersection between `o` and `items`.
    :raise: ValueError
        raise ValueError if `items` not in `o`.
    :return: list,
        `s` : object found in ``o` or the difference object i.e the object
        that is not in `items` provided that `error` is set to ``ignore``.
        Note that if None object is found  and `error` is ``ignore`` , it
        will return ``None``, otherwise, a `ValueError` raises.

    :example:
        >>> from watex.datasets import load_hlogs
        >>> from watex.utils.funcutils import is_in_if
        >>> X0, _= load_hlogs (as_frame =True )
        >>> is_in_if  (X0 , items= ['depth_top', 'top'])
        ... ValueError: Item 'top' is missing in the object
        >>> is_in_if (X0, ['depth_top', 'top'] , error ='ignore')
        ... ['depth_top']
        >>> is_in_if (X0, ['depth_top', 'top'] , error ='ignore',
                       return_diff= True)
        ... ['sp',
         'well_diameter',
         'layer_thickness',
         'natural_gamma',
         'short_distance_gamma',
         'strata_name',
         'gamma_gamma',
         'depth_bottom',
         'rock_name',
         'resistivity',
         'hole_id']
    """

    if isinstance (items, str):
        items =[items]
    elif not is_iterable(o):
        raise TypeError (f"Expect an iterable object, not {type(o).__name__!r}")
    # find intersect object
    s= set (o).intersection (items)

    miss_items = list(s.difference (o)) if len(s) > len(
        items) else list(set(items).difference (s))

    if return_diff or return_intersect:
        error ='ignore'

    if len(miss_items)!=0 :
        if error =='raise':
            v= smart_format(miss_items)
            verb = f"{ ' '+ v +' is' if len(miss_items)<2 else  's '+ v + 'are'}"
            raise ValueError (
                f"Item{verb} missing in the {type(o).__name__.lower()} {o}.")


    if return_diff :
        # get difference
        s = list(set(o).difference (s))  if len(o) > len(
            s) else list(set(items).difference (s))
        # s = set(o).difference (s)
    elif return_intersect:
        s = list(set(o).intersection(s))  if len(o) > len(
            items) else list(set(items).intersection (s))

    s = None if len(s)==0 else list (s)

    return s

def map_specific_columns (
        X: DataFrame,
        ufunc:F ,
        columns_to_skip:List[str]=None,
        pattern:str=None,
        inplace:bool= False,
        **kws
        ):
    r""" Apply function to a specific columns is the dataframe.

    It is possible to skip some columns that we want operation to not be
    performed.

    Parameters
    -----------
    X: dataframe,
        pandas dataframe with valid columns
    ufunc: callable,
        Universal function that can be applying to the dataframe.
    columns_to_skip: list or str ,
        List of columns to skip. If given as string and separed by the default
        pattern items, it should be converted to a list and make sure the
        columns name exist in the dataframe. Otherwise an error with
        raise.

    pattern: str, default = '[#&*@!,;\s]\s*'
        The base pattern to split the text in `column2skip` into a columns
        For instance, the following string coulb be splitted to::

            'depth_top, thickness, sp, gamma_gamma' ->
            ['depth_top', 'thickness', 'sp', 'gamma_gamma']

        Refer to :func:`~.str2columns` for further details.
    inplace: bool, default=True
        Modified dataframe in place and return None, otherwise return a
        new dataframe
    kws: dict,
        Keywords argument passed to :func: `pandas.DataFrame.apply` function

    Returns
    ---------
    X: Dataframe or None
        Dataframe modified inplace with values computed using the given
        `func`except the skipped columns, or ``None`` if `inplace` is ``True``.

    Examples
    ---------
    >>> from watex.datasets import load_hlogs
    >>> from watex.utils.plotutils import map_specific_columns
    >>> X0, _= load_hlogs (as_frame =True )
    >>> # let visualize the  first3 values of `sp` and `resistivity` keys
    >>> X0['sp'][:3] , X0['resistivity'][:3]
    ... (0   -1.580000
         1   -1.580000
         2   -1.922632
         Name: sp, dtype: float64,
         0    15.919130
         1    16.000000
         2    24.422316
         Name: resistivity, dtype: float64)
    >>> column2skip = ['hole_id','depth_top', 'depth_bottom',
                      'strata_name', 'rock_name', 'well_diameter', 'sp']
    >>> map_specific_columns (X0, ufunc = np.log10, column2skip)
    >>> # now let visualize the same keys values
    >>> X0['sp'][:3] , X0['resistivity'][:3]
    ... (0   -1.580000
         1   -1.580000
         2   -1.922632
         Name: sp, dtype: float64,
         0    1.201919
         1    1.204120
         2    1.387787
         Name: resistivity, dtype: float64)
    >>> # it is obvious the `resistiviy` values is log10
    >>> # while `sp` stil remains the same

    """
    X = _assert_all_types(X, pd.DataFrame)
    if not callable(ufunc):
        raise TypeError ("Expect a function for `ufunc`; "
                         f"got {type(ufunc).__name__!r}")

    pattern = pattern or r'[#&*@!,;\s]\s*'
    if not is_iterable( columns_to_skip):
        raise TypeError ("Columns  to skip expect an iterable object;"
                         f" got {type(columns_to_skip).__name__!r}")

    if isinstance(columns_to_skip, str):
        columns_to_skip = str2columns (columns_to_skip, pattern=pattern  )
    #assert whether column to skip is in
    if columns_to_skip:
        cskip = copy.deepcopy(columns_to_skip)
        columns_to_skip = is_in_if(X.columns, columns_to_skip, return_diff= True)
        if len(columns_to_skip) ==len (X.columns):
            warnings.warn("Value(s) to skip are not detected.", stacklevel=2)
    elif columns_to_skip is None:
        columns_to_skip = list(X.columns)

    if inplace :
        X[columns_to_skip] = X[columns_to_skip].apply (
            ufunc , **kws)
        X.drop (columns = cskip , inplace =True )
        return
    if not inplace:
        X0 = X.copy()
        X0[columns_to_skip] = X0[columns_to_skip].apply (
            ufunc , **kws)

        return  X0

def is_depth_in (X, name, columns = None, error= 'ignore'):
    """ Assert wether depth exists in the data from column attributes.

    If name is an integer value, it assumes to be the index in the columns
    of the dataframe if not exist , a warming will be show to user.

    :param X: dataframe
        dataframe containing the data for plotting

    :param columns: list,
        New labels to replace the columns in the dataframe. If given , it
        should fit the number of colums of `X`.

    :param name: str, int
        depth name in the dataframe or index to retreive the name of the depth
        in dataframe
    :param error: str , default='ignore'
        Raise or ignore when depth is not found in the dataframe. Whe error is
        set to ``ignore``, a pseudo-depth is created using the lenght of the
        the dataframe, otherwise a valueError raises.

    :return: X, depth
        Dataframe without the depth columns and depth values.
    """

    X= _assert_all_types( X, pd.DataFrame )
    if columns is not None:
        columns = list(columns)
        if not is_iterable(columns):
            raise TypeError("columns expects an iterable object."
                            f" got {type (columns).__name__!r}")
        if len(columns ) != len(X.columns):
            warnings.warn("Cannot rename columns with new labels. Expect "
                          "a size to be consistent with the columns X."
                          f" {len(columns)} and {len(X.columns)} are given.", stacklevel=2
                          )
        else :
            X.columns = columns # rename columns

    else:  columns = list(X.columns)

    _assert_all_types(name,str, int, float )

    # if name is given as indices
    # collect the name at that index
    if isinstance (name, (int, float) )  :
        name = int (name )
        if name > len(columns):
            warnings.warn ("Name index {name} is out of the columns range."
                           f" Max index of columns is {len(columns)}", stacklevel=2)
            name = None
        else :
            name = columns.pop (name)

    elif isinstance (name, str):
        # find in columns whether a name can be
        # found. Note that all name does not need
        # to be written completely
        # for instance name =depth can retrieved
        # ['depth_top, 'depth_bottom'] , in that case
        # the first occurence is selected i.e. 'depth_top'
        n = find_by_regex(
            columns, pattern=fr'{name}', func=re.search)

        if n is not None:
            name = n[0]

        # for consistency , recheck all and let
        # a warning to user
        if name not in columns :
            msg = f"Name {name!r} does not match any column names."
            if error =='raise':
                raise ValueError (msg)

            warnings.warn(msg, stacklevel=2)
            name =None

    # now create a pseudo-depth
    # as a range of len X
    if name is None:
        if error =='raise':
            raise ValueError ("Depth column not found in dataframe."
                              )
        depth = pd.Series ( np.arange ( len(X)), name ='depth (m)')
    else :
        # if depth name exists,
        # remove it from X
        depth = X.pop (name )

    return  X , depth


def count_func (path , verbose = 0 ):
    """ Count function and method using 'ast' modules

    Parameters
    -----------
    path: str, Path-like object,
        Path to the python module file
    verbose: int, default=0
        Different to 0 outputs the counting details.

    Returns
    -----------
    cobj or None: Returns the counter object from module `ast` or nothing if
        `verbose` is ``False``.

    """

    cobj ={}
    import_optional_dependency('ast')
    import ast
    class CountFunc (ast.NodeVisitor):
        func_count=0
        # def visit_FunctionDef(self, node):
        #     self.func_count +=1
        # def visit_Lambda(self, node):
        #     self.func_count +=1
        def visit_ClassDef(self, node):
            self.func_count +=1
        # def visit_Module(self, node):
        #     self.func_count +=1
        # def visit_Call(self, node):
        #     self.func_count +=1

    if os.path.isdir (path):
        pyfiles = [ os.path.join (path , f)
                   for f in os.listdir (path) if f.endswith ('.py') ]
    elif os.path.isfile (path) :
        pyfiles = [ path ]
    else :
        raise TypeError (f"Expects a path-like object, got {path!r}")

    val=0

    if verbose :
        print(f"module = {os.path.dirname (pyfiles[0]):^12}")
    for mod in pyfiles :

        p=ast.parse (open(mod, encoding='utf-8').read())
        f= CountFunc()
        f.visit(p)
        cobj[os.path.basename (mod)]= f.func_count
        val += f.func_count
        if verbose:
            print("### {:^7}{:<17} ={:>7}".format (' ', os.path.basename (mod),
                                              f.func_count ))

    print(f">>>Total = {val:>24}") if verbose else print()

    return cobj if not verbose else None


def smart_label_classifier (
        arr: ArrayLike, /, values: float | List[float]= None , labels =None,
        order ='soft', func: F=None, raise_warn=True):
    """ map smartly the numeric array into a class labels from a map function
    or a given fixed values.

    New classes created from the fixed values can be renamed if `labels`
    are supplied.

    Parameters
    -------------
    arr: Arraylike 1d,
        array-like whose items are expected to be categorized.

    values: float, list of float,
        The threshold item values from which the default categorization must
        be fixed.
    labels: int |str| or List of [str, int],
        The labels values that might be correspond to the fixed values. Note
        that the number of `fixed_labels` might be consistent with the fixed
        `values` plus one, otherwise a ValueError shall raise if `order` is
        set to ``strict``.

    order: str, ['soft'|'strict'], default='soft',
        If order is ``strict``, the argument passed to `values` must be self
        contain as item in the `arr`, and raise warning otherwise.

    func: callable, optional
        Function to map the given array. If given, values dont need to be
        supply.

    raise_warn: bool, default='True'
        Raise warning message if `order=soft` and the fixed `values` are not
        found in the `arr`. Also raise warnings, if `labels` arguments does
        not match the number of class from fixed `values`.

    Returns
    ----------
    arr: array-like 1d
        categorized array with the same length as the raw

    Examples
    ----------
    >>> import numpy as np
    >>> from watex.utils.funcutils import smart_label_classifier
    >>> sc = np.arange (0, 7, .5 )
    >>> smart_label_classifier (sc, values = [1, 3.2 ])
    array([0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2], dtype=int64)
    >>> # rename labels <=1 : 'l1', ]1; 3.2]: 'l2' and >3.2 :'l3'
    >>> smart_label_classifier (sc, values = [1, 3.2 ], labels =['l1', 'l2', 'l3'])
    >>> array(['l1', 'l1', 'l1', 'l2', 'l2', 'l2', 'l2', 'l3', 'l3', 'l3', 'l3',
           'l3', 'l3', 'l3'], dtype=object)
    >>> def f (v):
            if v <=1: return 'l1'
            elif 1< v<=3.2: return "l2"
            else : return "l3"
    >>> smart_label_classifier (sc, func= f )
    array(['l1', 'l1', 'l1', 'l2', 'l2', 'l2', 'l2', 'l3', 'l3', 'l3', 'l3',
           'l3', 'l3', 'l3'], dtype=object)
    >>> smart_label_classifier (sc, values = 1.)
    array([0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], dtype=int64)
    >>> smart_label_classifier (sc, values = 1., labels='l1')
    array(['l1', 'l1', 'l1', 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], dtype=object)

    """
    name =None
    from .validator import _is_arraylike_1d
    if hasattr(arr, "name") and isinstance (arr, pd.Series):
        name = arr.name

    arr= np.array (arr)

    if not _is_arraylike_1d(arr):
        raise TypeError ("Expects a one-dimensional array, got array with"
                         f" shape {arr.shape }")

    if isinstance (values, str):
        values = str2columns(values )
    if values is not None:
        values = is_iterable(values, parse_string =True, transform = True )
    # if (values is not None
    #     and not is_iterable( values)):
    #     values =[values ]

    if values is not None:
        approx_vs=list()
        values_ =np.zeros ((len(values), ), dtype =float )
        for i, v in enumerate (values ) :
            try : v= float (v)
            except TypeError as type_error :
                raise TypeError (
                    f"Value {v} must be a valid number." + str(type_error))
            diff_v = np.abs (arr[~np.isnan(arr)] - v )

            ix_v = np.argmin (diff_v)
            if order =='strict' and diff_v [ix_v]!=0. :
                raise ValueError (
                    f" Value {v} is missing the array. {v} must be an item"
                    " existing in the array or turn order to 'soft' for"
                    " approximate values selectors. ")

            # skip NaN in the case array contains NaN values
            values_[i] = arr[~np.isnan(arr)][ix_v]

            if diff_v [ix_v]!=0.:
                approx_vs.append ((v, arr[~np.isnan(arr)][ix_v]))

        if len(approx_vs) !=0 and raise_warn:
            vv, aa = zip (*approx_vs)
            verb ="are" if len(vv)>1 else "is"
            warnings.warn(f"Values {vv} are missing in the array. {aa} {verb}"
                          f" approximatively used for substituting the {vv}.", stacklevel=2)
    arr_ = arr.copy ()

    ####
    if (func is None and values is None ):
        raise TypeError ("'ufunc' cannot be None when the values are not given")

    dfunc =None

    if values is not None:
        def dfunc(k):
            return _smart_mapper (k, kr = values_ )
    func = func or  dfunc

    # func_vectorized  =np.vectorize(func )
    # arr_ = func_vectorized( arr )
    arr_ = pd.Series (arr_, name ='temp').map (func).values

    d={}
    if labels is not None:
        labels = is_iterable(labels, parse_string=True, transform =True )
        # if isinstance (labels, str):
        #     labels = str2columns(labels )
        labels, d = _assert_labels_from_values (
            arr_, values_ , labels , d, raise_warn= raise_warn , order =order
            )

    arr_ = arr_ if labels is None else (
        pd.Series (arr_, name = name or 'tname').map(d))

    # if name is None: # for consisteny if labels is None
    arr_= (arr_.values if labels is not None else arr_
           ) if name is None else pd.Series (arr_, name = name )

    return arr_

def _assert_labels_from_values (ar, values , labels , d=None,
                                raise_warn= True , order ='soft'):
    """ Isolated part of the :func:`~.smart_label_classifier`"""
    from .validator import _check_consistency_size

    if d is None:
        d = {}
    nlabels = list(np.unique (ar))
    if not is_iterable(labels):
        labels =[labels]
    if not _check_consistency_size(nlabels, labels, error ='ignore'):
        if order=='strict':
            verb= "were" if len (labels) > 1 else "was"
            raise TypeError (
                "Expect {len(nlabels)} labels for the {len(values)} values"
                f" renaming. {len(labels)} {verb} given.")

        verb ="s are" if len(values)>1 else " is"
        msg = (f"{len(values)} value{verb} passed. Labels for"
                " renaming values expect to be composed of"
                f" {len(values)+1} items i.e. 'number of values"
                " + 1' for pure categorization.")
        ur_classes = nlabels [len(labels):]
        labels = list(labels ) + ur_classes
        labels = labels [:len(nlabels)]
        msg += (f" Class{'es' if len(ur_classes)>1 else ''}"
                f" {smart_format(ur_classes)} cannot be renamed." )

        if raise_warn:
            warnings.warn (msg, stacklevel=2 )

    d = dict (zip (nlabels , labels ))

    return labels, d

def _smart_mapper (k, /,  kr , return_dict_map =False ) :
    """ Default  mapping using dict to validate the continue  value 'k'
    :param k: float,
        continue value to be framed between `kr`
    :param kr: Tuple,
        range of fixed values  to categorize
    :return: int - new categorical class

    :Example:
    >>> from watex.utils.funcutils import _smart_mapper
    >>> _smart_mapper (10000 , ( 500, 1500, 2000, 3500) )
    Out[158]: 4
    >>> _smart_mapper (10000 , ( 500, 1500, 2000, 3500) , return_dict_map=True)
    Out[159]: {0: False, 1: False, 2: False, 3: False, 4: True}

    """
    import math
    if len(kr )==1 :
        d = {0:  k <=kr[0], 1: k > kr[0]}
    elif len(kr)==2:
        d = {0: k <=kr[0], 1: kr[0] < k <= kr[1],  2: k > kr[1]}
    else :
        d= dict()
        for ii  in range (len(kr) + 1  ):
            if ii ==0:
                d[ii]= k <= kr[ii]
            elif ii == len(kr):
                d[ii] = k > kr [-1]
            else :
                d[ii] = kr[ii-1] < k <= kr[ii]

    if return_dict_map:
        return d

    for v, value in d.items () :
        if value: return v if not math.isnan (v) else np.nan

def hex_to_rgb (c, /):
    """ Convert colors Hexadecimal to RGB """
    c=c.lstrip('#')
    return tuple(int(c[i:i+2], 16) for i in (0, 2, 4))

def zip_extractor(
    zip_file ,
    samples ='*',
    ftype=None,
    savepath = None,
    pwd=None,
):
    """ Extract  ZIP archive objects.

    Can extract all or a sample objects when the number of object is passed
    to the parameter ``samples``.

    .. versionadded:: 0.1.5

    Parameters
    -----------
    zip_file: str
        Full Path to archive Zip file.
    samples: int, str, default ='*'
       Number of data to retrieve from archive files. This is useful when
       the archive file contains many data. ``*`` means extract all.
    savepath: str, optional
       Path to store the decompressed archived files.
    ftype: str,
<<<<<<< HEAD
       Is the extension of a specific file to decompress. Indeed, if the
=======
       Is the extension of a the specific file to decompress. Indeed, if the
>>>>>>> 10707dcecd7d0da55b83bcf73ae48c1e6659f2f8
       archived files contains many different data formats, specifying the
       data type would retrieve this specific files from the whole
       files archieved.
    pwd: int, optional
      Password to pass if the zip file is encrypted.

    Return
    --------
    objnames: list,
     List of decompressed objects.

    Examples
    ----------
    >>> from watex.utils.funcutils import zip_extractor
    >>> zip_extractor ('watex/datasets/data/edis/e.E.zip')

    """
    def raise_msg_when ( objn, ft):
        """ Raise message when None file is detected when the type of
        of file is supplied. Otherwise return the object collected
        from this kind of data-types
        """
        objn = [ o for o  in objn if o.endswith (ft)]
        if len(objn)  ==0:
            get_extension = [s.split('.')[-1] for s in objn if '.'  in s ]
            if len(get_extension)==0 : get_extension=['']
            msg = ( "The available file types are {smart_format(get_extension)}"
                   if len(get_extension)!=0 else ''
                   )
            raise ValueError (f"None objects in the zip collection of matches"
                              f"the {ft!r}. Available file types are {msg}")
        return objn

    if not os.path.isfile (zip_file ):
        raise FileExistsError( f"File {os.path.basename(zip_file)!r} does"
                              " not exist. Expect a Path-like object,"
                              f" got {type(zip_file).__name__!r}")

    if not os.path.basename(zip_file ).lower().endswith ('.zip'):
        raise FileNotFoundError("Unrecognized zip-file.")

    samples = str(samples)
    if samples !='*':
        try :samples = int (samples )
        except:
            raise ValueError ("samples must be an integer value"
                              f" or '*' not {samples}")

    with ZipFile (zip_file, 'r', ) as zip_obj :
        objnames = zip_obj.namelist()
        if samples =='*':
                samples = len(objnames )

        if ftype is not None:
            objnames = raise_msg_when(objn=objnames, ft= ftype)

        if ( samples >= len(objnames)
            and ftype is None
            ) :
            zip_obj.extractall( path = savepath , pwd=pwd)
        else:
            for zf in objnames [:samples ]:
                zip_obj.extract ( zf, path = savepath, pwd = pwd)

    return objnames


def remove_outliers (
    ar,
    method ='IQR',
    threshold = 3.,
    fill_value = None,
    axis = 1,
    interpolate=False,
    kind='linear'
    ):
    """ Efficient strategy to remove outliers in the data.

    Indeed, an outlier is the data point of the given sample,
    observation, or distribution that shall lie outside the overall pattern.
    A commonly used rule says that one will consider a data point an
    outlier if it has more than 1.5 IQR below the first quartile or above
    the third.

    Two approaches are used to remove the outliers.

    - Inter Quartile Range (``IQR``)
      IQR is the most commonly used and most trusted approach used in
      the research field. Said differently, low outliers shall
      lie below Q1-1.5 IQR, and high outliers shall lie Q3+1.5IQR.
      One needs to calculate median, quartiles, including IQR, Q1,
      and Q3.

      .. math::

        Q1 = 1/4(n + 1)

        Q3 = 1/4 (n + 1)

        Q2 = Q3 – Q1

      To define the outlier base value is defined above and below
      datasets normal range namely Upper and Lower bounds, define the
      upper and the lower bound (1.5*IQR value is considered) :

      .. math::

         upper = Q3 +1.5*IQR

         lower = Q1 – 1.5*IQR

      In the above formula as according to statistics, the 0.5
      scale-up of :math:`IQR (new_IQR = IQR + 0.5*IQR)` is taken, to consider
      all the data between 2.7 standard deviations in the Gaussian
      Distribution

    - Z-score
      Is also called a standard score. This value/score helps to understand
      that how far is the data point from the mean. And after setting up
      a threshold value one can utilize z score values of data points
      to define the outliers.

      .. math::

          Zscore = (\text{data_point} -\text{mean}) / \text{std. deviation}

    Now to define an outlier threshold value is chosen which is
    generally 3.0. As 99.7% of the data points lie between +/- 3 standard
    deviation (using Gaussian Distribution approach).

    .. versionadded: 0.1.5

    Parameters
    -----------
    ar: Arraylike, pd.dataframe
       Arraylike  containing outliers to remove.

       .. versionadded:: 0.2.7
          Accepts dataframe and can remove outliers using the `z_score`.

    method: str, default='IQR'
      The selected approach to remove the outliers. It can be
      ['IQR'|'Z-score']. See Above for outlier explanations.  Note that
      when selecting ``"z-score"`` the threshold value greatly influence
      the quality of data considering as ooutliers.

    threshold: float, default=3
      Thershold values is useful for ``"z-score"`` as the value for considering
      data above as outliers.

    fill_value: float, optional
      Value to replace the outliers. If not given, outliers are suppressed
      in the array.

    axis: int, default=1
      axis from which to remove values. This is useful when two dimensional
      array is supplied. Default, delete outlier from the rows.

    interpolate: bool, default=False,
       If ``fill_value='NaN'``, interpolation can be triggered to get the
       closest value in array to replace missing values. Note that
       `fill_value` should be NaN for interpolation to be concise.

    kind: str, default='linear'
      kind of interpolation. It could be ['nearest'|'linear'|'cubic'].

    .. versionadded:: 0.2.8
       Interpolate NaN value after outliers removal.


    Returns
    --------
    arr: Array_like
        New array whith removed outliers.

    Examples
    ---------
    >>> import numpy as np
    >>> np.random.seed (42 )
    >>> from watex.utils.funcutils import remove_outliers
    >>> data = np.random.randn (7, 3 )
    >>> data_r = remove_outliers ( data )
    >>> data.shape , data_r.shape
    (7, 3) (5, 3)
    >>> remove_outliers ( data, fill_value =np.nan )
    array([[ 0.49671415, -0.1382643 ,  0.64768854],
           [ 1.52302986, -0.23415337, -0.23413696],
           [ 1.57921282,  0.76743473, -0.46947439],
           [ 0.54256004, -0.46341769, -0.46572975],
           [ 0.24196227,         nan,         nan],
           [-0.56228753, -1.01283112,  0.31424733],
           [-0.90802408,         nan,  1.46564877]])
    >>> # for one dimensional
    >>> remove_outliers ( data[:, 0] , fill_value =np.nan )
    array([ 0.49671415,  1.52302986,  1.57921282,  0.54256004,  0.24196227,
           -0.56228753,         nan])
    >>> remove_outliers ( data[:, 0] , fill_value =np.nan, interpolate=True  )
    >>> import matplotlib.pyplot as plt
    >>> plt.plot (np.arange (len(data ), data, ))
    """
    method = str(method).lower()
    if (
            hasattr ( ar, "__array__")
            and hasattr (ar, 'columns')
            ):
        return _remove_outliers( ar, n_std= threshold )

    arr =np.array (ar, dtype = float)

    if method =='iqr':
        Q1 = np.percentile(arr[~np.isnan(arr)], 25,)
        Q3 = np.percentile(arr[~np.isnan(arr)], 75)
        IQR = Q3 - Q1

        upper = Q3 + 1.5 * IQR

        upper_arr = np.array (arr >= upper)
        lower = Q3 - 1.5 * IQR
        lower_arr =  np.array ( arr <= lower )
        # replace the oulier by nan
        arr [upper_arr]= fill_value if fill_value else np.nan
        arr[ lower_arr]= fill_value if fill_value else np.nan

    if method =='z-score':
        from scipy import stats
        z = np.abs(stats.zscore(arr[~np.isnan(arr)]))
        zmask  = np.array ( z > threshold )
        arr [zmask]= fill_value if fill_value else np.nan

    if fill_value is None:
        # delete nan if fill value is not provided
        arr = arr[ ~np.isnan (arr ).any(axis =1)
                  ]  if np.ndim (arr) > 1 else arr [~np.isnan(arr)]

    if interpolate:
        arr = interpolate_grid (arr, method = kind )

    return arr

def _remove_outliers(data, n_std=3):
    """Remove outliers from a dataframe."""
    # separate cat feature and numeric features
    # if exists
    df, numf, catf = to_numeric_dtypes(
        data , return_feature_types= True,  drop_nan_columns =True )
    # get on;y the numeric
    df = df[numf]
    for col in df.columns:
        # print('Working on column: {}'.format(col))
        mean = df[col].mean()
        sd = df[col].std()
        df = df[(df[col] <= mean+(n_std*sd))]
    # get the index and select only the index
    index = df.index
    # get the data by index then
    # concatename
    df_cat = data [catf].iloc [ index ]
    df = pd.concat ( [df_cat, df ], axis = 1 )

    return df

def normalizer ( arr, /, method ='naive'):
    """ Normalize values to be between 0 and 1.

    This normlizer handles NaN values translates data individually such
    that it is in the given range on the training set, e.g. between
    zero and one.

    Note that when the transformation is set to the ``method ='MinMax'``,
    The transformation is given by::

        X_std = (X - X.min(axis=0)) / (X.max(axis=0) - X.min(axis=0))
        X_normed = X_std * (max - min) + min

    where min, max = feature_range.

    This transformation is often used as an alternative to zero mean,
    unit variance scaling.

    Parameters
    -----------
    arr: Arraylike,
       Array to normalize, can contain NaN values.
    method: str,
       Can be use 'scikit-learn' :class:`~watex.exlib.MinMaxScaler` for
       normalization. Any other values used the naive normalization.

    Returns
    --------
    arr_norm: Normalized array.

    Examples
    ----------
    >>> import numpy as np
    >>> from watex.utils.funcutils import normalizer
    >>> np.random.seed (42)
    >>> arr = np.random.randn (3, 2 )
    array([[ 0.49671415, -0.1382643 ],
           [ 0.64768854,  1.52302986],
           [-0.23415337, -0.23413696]])
    >>> normalizer (arr )
    array([[4.15931313e-01, 5.45697636e-02],
           [5.01849720e-01, 1.00000000e+00],
           [0.00000000e+00, 9.34323403e-06]])
    >>> normalizer (arr , method ='min-max')  # normalize data along axis=0
    array([[0.82879654, 0.05456093],
           [1.        , 1.        ],
           [0.        , 0.        ]])
    >>> arr [0, 1] = np.nan; arr [1, 0] = np.nan
    >>> normalizer (arr )
    array([[4.15931313e-01,            nan],
           [           nan, 1.00000000e+00],
           [0.00000000e+00, 9.34323403e-06]])
    >>> normalizer (arr , method ='min-max')
    array([[ 1., nan],
           [nan,  1.],
           [ 0.,  0.]])

    """
    method = str(method).lower()
    arr = np.array(arr )

    if method in ( 'sklearn', 'scikit-learn', 'minmax', 'min-max'):
        from ..exlib import MinMaxScaler
        arr = arr.reshape(-1, 1) if arr.ndim ==1 else arr
        return  MinMaxScaler().fit_transform(arr )

    arr_norm  = (arr - np.nanmin(arr))/ (np.nanmax (arr) - np.nanmin(arr))

    return arr_norm

def _validate_name_in (name, /, defaults = '', expect_name= None,
                         exception = None , deep=False ):
    """ Assert name in multiples given default names.

    Parameters
    -----------
    name: str,
      given name to assert
    default: list, str, default =''
      default values used for assertion
    expect_name: str, optional
      name to return in case assertion is verified ( as ``True``)
    deep: bool, default=False
      Find item in a litteral default string. If set  to ``True``,
      `defaults` are joined and check whether an occurence of `name` is in the
      defaults

    exception: Exception
      Error to raise if name is not found in the default values.

    Returns
    -------
    name: str,
      Verified name or boolean if expect name if ``None``.

    Examples
    -------
    >>> from watex.utils.funcutils import _validate_name_in
    >>> dnames = ('NAME', 'FIST NAME', 'SUrname')
    >>> _validate_name_in ('name', defaults=dnames )
    False
    >>> _validate_name_in ('name', defaults= dnames, deep =True )
    True
    >>> _validate_name_in ('name', defaults=dnames , expect_name ='NAM')
    False
    >>> _validate_name_in ('name', defaults=dnames , expect_name ='NAM', deep=True)
    'NAM'
    """

    name = str(name).lower().strip()
    defaults = is_iterable(defaults,
            exclude_string= True, parse_string= True, transform=True )
    if deep :
        defaults = ''.join([ str(i) for i in defaults] )

    # if name in defaults:
    name = ( True if expect_name is None  else expect_name
            ) if name in defaults else False

    #name = True if name in defaults else ( expect_name if expect_name else False )

    if not name and exception:
        raise exception

    return name

def get_confidence_ratio (
        ar, /,
        axis = 0,
        invalid = 'NaN',
        mean=False,
        ):

    """ Get ratio of confidence in array by counting the number of
    invalid values.

    Parameters
    ------------
    ar: arraylike 1D or 2D
      array for checking the ratio of confidence

    axis: int, default=0,
       Compute the ratio of confidence alongside the rows by defaults.

    invalid: int, foat, default='NaN'
      The value to consider as invalid in the data might be listed if
      applicable. The default is ``NaN``.

    mean: bool, default=False,
      Get the mean ratio. Average the percentage of each axis.

      .. versionadded:: 0.2.8
         Average the ratio of confidence of each axis.

    Returns
    ---------
    ratio: arraylike 1D
      The ratio of confidence array alongside the ``axis``.

    Examples
    ----------
    >>> import numpy as np
    >>> np.random.seed (0)
    >>> test = np.random.randint (1, 20 , 10 ).reshape (5, 2 )
    >>> test
    array([[13, 16],
           [ 1,  4],
           [ 4,  8],
           [10, 19],
           [ 5,  7]])
    >>> from watex.utils.funcutils import get_confidence_ratio
    >>> get_confidence_ratio (test)
    >>> array([1., 1.])
    >>> get_confidence_ratio (test, invalid= ( 13, 19) )
    array([0.8, 0.8])
    >>> get_confidence_ratio (test, invalid= ( 13, 19, 4) )
    array([0.6, 0.6])
    >>> get_confidence_ratio (test, invalid= ( 13, 19, 4), axis =1 )
    array([0.5, 0.5, 0.5, 0.5, 1. ])

    """
    def gfc ( ar, inv):
        """ Get ratio in each column or row in the array. """
        inv = is_iterable(inv, exclude_string=True , transform =True,
                              )
        # if inv!='NaN':
        for iv in inv:
            if iv in ('NAN', np.nan, 'NaN', 'nan', None):
                iv=np.nan
            ar [ar ==iv] = np.nan

        return len( ar [ ~np.isnan (ar)])  / len(ar )

    # validate input axis name
    axis = _validate_name_in (axis , ('1', 'rows', 'sites', 'stations') ,
                              expect_name=1 )
    if not axis:
        axis =0

    ar = np.array(ar).astype ( np.float64) # for consistency
    ratio = np.zeros(( (ar.shape[0] if axis ==1 else ar.shape [1] )
                      if ar.ndim ==2 else 1, ), dtype= np.float64)

    for i in range (len(ratio)):
        ratio[i] = gfc ( (ar [:, i] if axis ==0 else ar [i, :])
                        if ar.ndim !=1 else ar , inv= invalid
                        )
    if mean:
        ratio = np.array (ratio).mean()
    return ratio

def assert_ratio(
        v, /, bounds: List[float] = None ,
        exclude_value:float= None,
        in_percent:bool =False , name:str ='rate'
        ):
    """ Assert rate value between a specific range.

    Parameters
    -----------
    v: float,
       ratio value to assert
    bounds: list ( lower, upper)
       The range that value must  be included
    exclude_value: float
       A value that ``v`` must not taken. Exclude it from the ``bounds``.
       Raise error otherwise. Note that  any other value will use the
       lower bound in `bounds` as exlusion.

    in_percent: bool, default=False,
       Convert the value into a percentage.

       .. versionchanged:: 0.2.3
          `as_percent` parameter is changed to `in_percent`.

    name: str, default='rate'
       the name of the value for assertion.

    Returns
    --------
    v: float
       Asserted value.

    Examples
    ---------
    >>> from watex.utils.funcutils import assert_ratio
    >>> assert_ratio('2')
    2.0
    >>> assert_ratio(2 , bounds =(2, 8))
    2.0
    >>> assert_ratio(2 , bounds =(4, 8))
    ValueError:...
    >>> assert_ratio(2 , bounds =(1, 8), exclude_value =2 )
    ValueError: ...
    >>> assert_ratio(2 , bounds =(1, 8), exclude_value ='use bounds' )
    2.0
    >>> assert_ratio(2 , bounds =(0, 1) , in_percent =True )
    0.02
    >>> assert_ratio(2 , bounds =(0, 1) )
    ValueError:
    >>> assert_ratio(2 , bounds =(0, 1), exclude_value ='use lower bound',
                         name ='tolerance', in_percent =True )
    0.02
    """
    msg =("greater than {} and less than {}" )


    if isinstance (v, str):
        if "%" in v: in_percent=True
        v = v.replace('%', '')
    try :
        v = float (v)
    except TypeError :
        raise TypeError (f"Unable to convert {type(v).__name__!r} "
                         f"to float: {v}")
    except ValueError:
        raise ValueError(f"Expects 'float' not {type(v).__name__!r}: "
                         f"{(v)!r}")
    # put value in percentage
    # if greater than 1.
    if in_percent:
        if 1 < v <=100:
            v /= 100.

    bounds = bounds or []
    low, up, *_ = list(bounds) + [ None, None]
    e=(f"Expects a {name} value {msg.format(low, up)}, got: {v}")
    err = ValueError (e)

    if len(bounds)!=0:
        if (
                low is not None  # use is not None since 0. is
                and up is not None # consider as False value
            and  (v < low or v > up)
            ) :
                raise err

    if exclude_value is not None:
        try :
            low = float (str(exclude_value))
        except : # use bounds
            pass
        if low is None:
            warnings.warn("Cannot exclude the lower value in the interval"
                          " while `bounds` argument is not given.", stacklevel=2)
        else:
            if v ==low:
                raise ValueError (e.replace (", got:", ' excluding') + ".")

    if in_percent and v > 100:
         raise ValueError (f"{name.title()} value should be {msg.format(low, up)}, got: {v}")
    return v

def exist_features (df, features, error='raise'):
    """Control whether the features exist or not

    :param df: a dataframe for features selections
    :param features: list of features to select. Lits of features must be in the
        dataframe otherwise an error occurs.
    :param error: str - raise if the features don't exist in the dataframe.
        *default* is ``raise`` and ``ignore`` otherwise.

    :return: bool
        assert whether the features exists
    """
    isf = False

    error= 'raise' if error.lower().strip().find('raise')>= 0  else 'ignore'

    if isinstance(features, str):
        features =[features]

    features = _assert_all_types(features, list, tuple, np.ndarray)
    set_f =  set (features).intersection (set(df.columns))
    if len(set_f)!= len(features):
        nfeat= len(features)
        msg = f"Feature{'s' if nfeat >1 else ''}"
        if len(set_f)==0:
            if error =='raise':
                raise ValueError (f"{msg} {smart_format(features)} "
                                  f"{'does not' if nfeat <2 else 'dont'}"
                                  " exist in the dataframe")
            isf = False
        # get the difference
        diff = set (features).difference(set_f) if len(
            features)> len(set_f) else set_f.difference (set(features))
        nfeat= len(diff)
        if error =='raise':
            raise ValueError(f"{msg} {smart_format(diff)} not found in"
                             " the dataframe.")
        isf = False
    else : isf = True

    return isf

def interpolate_grid (
    arr, / ,
    method ='cubic',
    fill_value='auto',
    view = False,
    ):
    """
    Interpolate data containing missing values.

    Parameters
    -----------
    arr: ArrayLike2D
       Two dimensional array for interpolation
    method: str, default='cubic'
      kind of interpolation. It could be ['nearest'|'linear'|'cubic'].

    fill_value: float, str, default='auto'
       Fill the interpolated grid at the egdes or surrounding NaN with
       a filled value. The ``auto`` uses the forward and backward
       fill strategy.

    view: bool, default=False,
       Quick visualize the interpolated grid.


    .. versionchanged:: 0.2.8
       One-dimensional array is henceforth possible. Error no longer raises.

    Returns
    ---------
    arri: ArrayLike2d
       Interpolated 2D grid.

    See also
    ---------
    spi.griddata:
        Scipy interpolate Grid data
    fillNaN:
        Fill missing data strategy.

    Examples
    ---------
    >>> import numpy as np
    >>> from watex.utils.funcutils import interpolate_grid
    >>> x = [28, np.nan, 50, 60] ; y = [np.nan, 1000, 2000, 3000]
    >>> xy = np.vstack ((x, y)).T
    >>> xyi = interpolate_grid (xy, view=True )
    >>> xyi
    array([[  28.        ,   28.        ],
           [  22.78880663, 1000.        ],
           [  50.        , 2000.        ],
           [  60.        , 3000.        ]])

    """
    is2d = True
    if not hasattr(arr, '__array__'):
        arr = np.array (arr)

    if arr.ndim==1:
        #convert to two dimension array
        arr = np.vstack ((arr, arr ))
        is2d =False
        # raise TypeError(
        #     "Expect two dimensional array for grid interpolation.")

    # make x, y array for mapping
    x = np.arange(0, arr.shape[1])
    y = np.arange(0, arr.shape[0])
    #mask invalid values
    arr= np.ma.masked_invalid(arr)
    xx, yy = np.meshgrid(x, y)
    #get only the valid values
    x1 = xx[~arr.mask]
    y1 = yy[~arr.mask]
    newarr = arr[~arr.mask]

    arri = spi.griddata(
        (x1, y1),
        newarr.ravel(),
        (xx, yy),
        method=method
        )

    if fill_value =='auto':
        arri = fillNaN(arri, method ='both ')
    else:
        arri [np.isnan(arri)] = float( _assert_all_types(
            fill_value, float, int, objname ="'fill_value'" )
            )

    if view :
        fig, ax  = plt.subplots (nrows = 1, ncols = 2 , sharey= True, )
        ax[0].imshow(arr ,interpolation='nearest', label ='Raw Grid')
        ax[1].imshow (arri, interpolation ='nearest',
                      label = 'Interpolate Grid')

        ax[0].set_title ('Raw Grid')
        ax[1].set_title ('Interpolate Grid')

        plt.show ()

    if not is2d:
        arri = arri[0, :]

    return arri


def random_selector (
        arr:ArrayLike, / , value: float | ArrayLike,
        seed: int = None, shuffle =False ):
    """Randomly select the number of values in array.

    Parameters
    ------------
    arr: ArrayLike
       Array of values
    value: float, arraylike
        If ``float`` value is passed, it indicates the number of values to
        select among the length of `arr`. If array (``value``) is passed, it
        should be self contain in the given ``arr`. However if ``string`` is
        given and contain the ``%``, it calculates the ratio of
        number to randomly select.
    seed: int, Optional
       Allow retrieving the identical value randomly selected in the given
       array.

    suffle: bool, False
       If  ``True`` , shuffle the selected values.

    Returns
    --------
    arr: Array containing the selected values

    Examples
    ----------
    >>> import numpy as np
    >>> from watex.utils.funcutils import random_selector
    >>> dat= np.arange (42 )
    >>> random_selector (dat , 7, seed = 42 )
    array([0, 1, 2, 3, 4, 5, 6])
    >>> random_selector ( dat, ( 23, 13 , 7))
    array([ 7, 13, 23])
    >>> random_selector ( dat , "7%", seed =42 )
    array([0, 1])
    >>> random_selector ( dat , "70%", seed =42 , shuffle =True )
    array([ 0,  5, 20, 25, 13,  7, 22, 10, 12, 27, 23, 21, 16,  3,  1, 17,  8,
            6,  4,  2, 19, 11, 18, 24, 14, 15,  9, 28, 26])
    """

    msg = "Non-numerical is not allowed. Got {!r}."

    if seed:
        seed = _assert_all_types(seed , int, float, objname ='Seed')
        np.random.seed (seed )

    v = copy.deepcopy(value )

    if not is_iterable( value, exclude_string= True ):

        value = str(value )

        if '%' in  value:
            try:
               value = float( value.replace ('%', '')) /100
            except :
                raise TypeError(msg.format(v))
            # get the number
            value *= len(arr )


        try :
            value = int(value )

        except :
            raise TypeError (msg.format(v))

        if value > len(arr):
            raise ValueError(f"Number {value} is out of the range."
                             f" Expect value less than {len(arr)}.")

        value = np.random.permutation(value )

    arr = np.array (
        is_iterable( arr, exclude_string=True, transform =True ))

    arr = arr.ravel() if arr.ndim !=1 else arr

    mask = _isin (arr, value , return_mask= True )
    arr = arr [mask ]

    if shuffle : np.random.shuffle (arr )

    return arr

def cleaner (
    data: DataFrame|NDArray,
    / ,
    columns:List[str]= None,
    inplace:bool = False,
    labels: List[int|str] =None,
    func : F= None,
    mode:str ='clean',
    **kws
    )->DataFrame | NDArray | None :
    """ Sanitize data or columns by dropping specified labels
    from rows or columns.

    If data is not a pandas dataframe, should be converted to
    dataframe and uses index to drop the labels.

    Parameters
    -----------
    data: pd.Dataframe or arraylike2D.
       Dataframe pandas or Numpy two dimensional arrays. If 2D array is
       passed, it should prior be converted to a daframe by default and
       drop row index from index parameters

    columns: single label or list-like
        Alternative to specifying axis (
            labels, axis=1 is equivalent to columns=labels).

    labels: single label or list-like
      Index or column labels to drop. A tuple will be used as a single
      label and not treated as a list-like.

    func: F, callable
        Universal function used to clean the columns. If performs only when
        `mode` is on ``clean`` option.

    inplace: bool, default False
        If False, return a copy. Otherwise, do operation
        inplace and return None.

    mode: str, default='clean'
       Options or mode of operation to do on the data. It could
       be ['clean'|'drop']. If ``drop``, it behaves like ``dataframe.drop``
       of pandas.

    Returns
    --------
    DataFrame, array2D  or None
            DataFrame cleaned or without the removed index or column labels
            or None if inplace=True or array is data is passed as an array.

    """
    mode = _validate_name_in(mode , defaults =("drop", 'remove' ),
                      expect_name ='drop')
    if not mode:
        return sanitize_frame_cols(
            data,
            inplace = inplace,
            func = func
            )

    objtype ='ar'
    if not hasattr (data , '__array__'):
        data = np.array (data )

    if hasattr(data , "columns"):
        objtype = "pd"

    if objtype =='ar':
        data = pd.DataFrame(data )
        # block inplace to False and
        # return numpy ar
        inplace = False
    # if isinstance(columns , str):
    #     columns = str2columns(columns )
    if columns is not None:
        columns = is_iterable(
            columns, exclude_string=True ,
            parse_string= True,
            transform =True )

    data = data.drop (labels = labels,
                      columns = columns,
                      inplace =inplace,
                       **kws
                       )

    return np.array ( data ) if objtype =='ar' else data

def rename_files (
    src_files:str | List[str], /,
    dst_files:str | List[str],
    basename:Optional[str]=None,
    extension:Optional[str] =None ,
    how:str ='py',
    prefix:bool =True,
    keep_copy:bool=True,
    trailer:str='_',
    sortby: re |F=None,
    **kws
    ):
    """Rename files in directory.

    Parameters
    -----------
    src_files: str, Path-like object
       Source files to rename

    dst_files: str of PathLike object
       Destination files renamed.

    extension: str, optional
       If a path is given in `src_files`, specifying the `extension` will just
       collect only files with this typical extensions.

    basename: str, optional
       If `dst_files` is passed as Path-object, name should be needed
       for a change, otherwise, the number is incremented using the Python
       index counting defined by the parameter ``how=py`

    how: str, default='py'
       The way to increment files when `dst_files` is given as a Path object.
       For instance, for a  ``name=E_survey`` and ``prefix==True``, the first
       file should be ``E_survey_00`` if ``how='py'`` otherwise it should be
       ``E_survey_01``.

    prefix: bool, default=True
      Prefix is used to position the name before the number incrementation.
      If ``False`` and `name` is given, the number is positionning before the
      name. If ``True`` and not `prefix` for a ``name=E_survey``, it should be
      ``00_E_survey`` and ``01_E_survey``.

    keep_copy: bool, default=True
       Keep a copy of the source files.

    trailer: str, default='_',
       Item used to separate the basename for counter.

    sortby: Regex or Callable,
       Key to sort the collection of the items when `src_files` is passed as
       a path-like object.  This is usefull to keep order as the origin files
       especially  when files includes a specific character.  Furthermore
       [int| float |'num'|'digit'] sorted the files according to the
       number included in the filename if exists.

    kws: dict
       keyword arguments passed to `os.rename`.

    """
    dest_dir =None ; trailer = str(trailer)
    extension = str(extension).lower()

    if os.path.isfile (src_files ):
        src_files = [src_files ]

    elif os.path.isdir (src_files):
        src_path = src_files
        ldir = os.listdir(src_path)

        src_files = ldir if extension =='none' else [
             f for f in ldir  if  f.endswith (extension) ]

        if sortby:
            if sortby in ( int, float, 'num', 'number', 'digit'):
                src_files = sorted(ldir, key=lambda s:int( re.search(
                    r'\d+', s).group()) if re.search(r'\d+', s) else 0 )
            else:
                src_files = sorted(ldir, key=sortby)

        src_files = [  os.path.join(src_path, f )   for f in src_files  ]
        # get only the files
        src_files = [ f for f in src_files if os.path.isfile (f ) ]

    else : raise FileNotFoundError(f"{src_files!r} not found.")

    # Create the directory if it doesn't exist
    if ( dst_files is not None
        and not os.path.exists (dst_files)
        ):
        os.makedirs(dst_files)

    if os.path.isdir(dst_files):
        dest_dir = dst_files

    dst_files = is_iterable(dst_files , exclude_string= True, transform =True )
    # get_extension of the source_files
    _, ex = os.path.splitext (src_files[0])

    if dest_dir:
        if basename is None:
            warnings.warn(
                "Missing basename for renaming file. Should use `None` instead.", stacklevel=2)
            basename =''; trailer =''

        basename= str(basename)
        if prefix:
            dst_files =[ f"{str(basename)}{trailer}" + (
                f"{i:03}" if how=='py' else f"{i+1:03}") + f"{ex}"
                        for i in range (len(src_files))]
        elif not prefix:
            dst_files =[ (f"{i:03}" if how=='py' else f"{i+1:03}"
                        ) +f"{trailer}{str(basename)}" +f"{ex}"
                        for i in range (len(src_files))]

        dst_files = [os.path.join(dest_dir , f) for f in dst_files ]

    for f, nf in zip (src_files , dst_files):
        try:
           if keep_copy : shutil.copy (f, nf , **kws )
           else : os.rename (f, nf , **kws )
        except FileExistsError:
            os.remove(nf)
            if keep_copy : shutil.copy (f, nf , **kws )
            else : os.rename (f, nf , **kws )


def get_xy_coordinates (d, / , as_frame = False, drop_xy = False,
                        raise_exception = True, verbose=0 ):
    """Check whether the coordinate values exist in the data

    Parameters
    ------------
    d: Dataframe
       Frame that is expected to contain the longitude/latitude or
       easting/northing coordinates.  Note if all types of coordinates are
       included in the data frame, the longitude/latitude takes the
       priority.
    as_frame: bool, default= False,
       Returns the coordinates values if included in the data as a frame rather
       than computing the middle points of the line
    drop_xy: bool, default=False,
       Drop the coordinates in the data and return the data transformed inplace

    raise_exception: bool, default=True
       raise error message if data is not a dataframe. If set to ``False``,
       exception is converted to a warning instead. To mute the warning set
       `raise_exception` to ``mute``

    verbose: int, default=0
      Send message whether coordinates are detected.

    returns
    --------
    xy, d, xynames: Tuple
      xy : tuple of float ( longitude, latitude) or (easting/northing )
         if `as_frame` is set to ``True``.
      d: Dataframe transformed (coordinated removed )  or not
      xynames: str, the name of coordinates detected.

    Examples
    ----------
    >>> import watex as wx
    >>> from watex.utils.funcutils import get_xy_coordinates
    >>> testdata = wx.make_erp ( n_stations =7, seed =42 ).frame
    >>> xy, d, xynames = get_xy_coordinates ( testdata,  )
    >>> xy , xynames
    ((110.48627946874444, 26.051952363176344), ('longitude', 'latitude'))
    >>> xy, d, xynames = get_xy_coordinates ( testdata, as_frame =True  )
    >>> xy.head(2)
        longitude   latitude        easting      northing
    0  110.485833  26.051389  448565.380621  2.881476e+06
    1  110.485982  26.051577  448580.339199  2.881497e+06
    >>> # remove longitude and  lat in data
    >>> testdata = testdata.drop (columns =['longitude', 'latitude'])
    >>> xy, d, xynames = get_xy_coordinates ( testdata, as_frame =True  )
    >>> xy.head(2)
             easting      northing
    0  448565.380621  2.881476e+06
    1  448580.339199  2.881497e+06
    >>> # note testdata should be transformed inplace when drop_xy is set to True
    >>> xy, d, xynames = get_xy_coordinates ( testdata, drop_xy =True)
    >>> xy, xynames
    ((448610.25612032827, 2881538.4380570543), ('easting', 'northing'))
    >>> d.head(2)
       station  resistivity
    0      0.0          1.0
    1     20.0        167.5
    >>> testdata.head(2) # coordinates are henceforth been dropped
       station  resistivity
    0      0.0          1.0
    1     20.0        167.5
    >>> xy, d, xynames = get_xy_coordinates ( testdata, drop_xy =True)
    >>> xy, xynames
    (None, ())
    >>> d.head(2)
       station  resistivity
    0      0.0          1.0
    1     20.0        167.5

    """

    def get_value_in ( val, /, col , default):
        """ Get the value in the frame columns if `val` exists in """
        x = list( filter ( lambda x: x.find (val)>=0 , col)
                   )
        if len(x) !=0:
            # now rename col
            d.rename (columns = {x[0]: str(default) }, inplace = True )

        return d

    if not (
            hasattr ( d, 'columns') and hasattr ( d, '__array__')
            ) :
        emsg = ("Expect dataframe containing coordinates longitude/latitude"
                f" or easting/northing. Got {type (d).__name__!r}")

        raise_exception = str(raise_exception).lower().strip()
        if raise_exception=='true':
            raise TypeError ( emsg )

        if raise_exception  not in ('mute', 'silence'):
            warnings.warn( emsg, stacklevel=2 )

        return d

    # check whether coordinates exists in the data columns
    for name, tname in zip ( ('lat', 'lon', 'east', 'north'),
                     ( 'latitude', 'longitude', 'easting', 'northing')
                     ) :
        d = get_value_in(name, col = d.columns , default = tname )

    # get the exist coodinates
    coord_columns  = []
    for x, y in zip ( ( 'longitude', 'easting' ), ( 'latitude', 'northing')):
        if ( x  in d.columns and y in d.columns ):
            coord_columns.extend  ( [x, y] )

    xy  = d[ coord_columns] if len(coord_columns)!=0 else None

    if ( not as_frame
        and xy is not None ) :
        # take the middle of the line and if both types of
        # coordinates are supplied , take longitude and latitude
        # and drop easting and northing
        xy = tuple ( np.nanmean ( np.array ( xy ) , axis =0 )) [:2]

    xynames = tuple ( coord_columns)[:2]
    if (
            drop_xy  and len( coord_columns) !=0
            ):
        # modifie the data inplace
        d.drop ( columns=coord_columns, inplace = True  )

    if verbose:
        print("###", "No" if len(xynames)==0 else (
            tuple (xy.columns) if as_frame else xy), "coordinates found.")

    return  xy , d , xynames


def twinning(
    *d: DataFrame,
    on:str | List[str] = None,
    parse_on:bool=False,
    mode: str='strict',
    coerce:bool = False,
    force:bool =False,
    decimals: int = 7,
    raise_warn:bool =True
)-> DataFrame :
    """ Find indentical object in all data and concatenate them using merge
     intersection (`cross`) strategy.

    Parameters
    ----------
    d: List of DataFrames
       List of pandas DataFrames
    on: str, label or list
       Column or index level names to join on. These must be found in
       all DataFrames. If `on` is ``None`` and not merging on indexes then
       a concatenation along columns axis is performed in all DataFrames.
       Note that `on` works with `parse_on` if its argument is  a list of
       columns names passed into single litteral string. For instance::

        on ='longitude latitude' --[parse_on=True]-> ['longitude' , 'latitude']

    parse_on: bool, default=False
       Parse `on` arguments if given as string and return_iterable objects.

    mode: str, default='strict'
      Mode to the data. Can be ['soft'|'strict']. In ``strict`` mode, all the
      data passed must be a DataFrame, otherwise an error raises. in ``soft``
      mode, ignore the non-DataFrame. Note that any other values should be
      in ``strict`` mode.

    coerce: bool, default=False
       Truncate all DataFrame size to much the shorter one before performing
       the ``merge``.

    force: bool, default=False,
       Force `on` items to be in the all DataFrames, This could be possible
       at least, `on` items should be in one DataFrame. If missing in all
       data, an error occurs.

    decimals: int, default=5
       Decimal is used for comparison between numeric labels in `on` columns
       items. If set, it rounds values of `on` items in all data before
       performing the merge.

     raise_warn: bool, default=False
        Warn user to concatenate data along column axis if `on` is ``None``.

    Returns
    --------
    data: DataFrames
      A DataFrame of the merged objects.

    Examples
    ----------
    >>> import watex as wx
    >>> from watex.utils.funcutils import twinning
    >>> data = wx.make_erp (seed =42 , n_stations =12, as_frame =True )
    >>> table1 = wx.methods.DCProfiling ().fit(data).summary()
    >>> table1
           dipole   longitude  latitude  ...  shape  type       sfi
    line1      10  110.486111  26.05174  ...      C    EC  1.141844
    >>> data_no_xy = wx.make_ves ( seed=0 , as_frame =True)
    >>> data_no_xy.head(2)
        AB   MN  resistivity
    0  1.0  0.4   448.860148
    1  2.0  0.4   449.060335
    >>> data_xy = wx.make_ves ( seed =0 , as_frame =True , add_xy =True )
    >>> data_xy.head(2)
        AB   MN  resistivity   longitude  latitude
    0  1.0  0.4   448.860148  109.332931  28.41193
    1  2.0  0.4   449.060335  109.332931  28.41193
    >>> table = wx.methods.VerticalSounding (
        xycoords = (110.486111,   26.05174)).fit(data_no_xy).summary()
    >>> table.table_
             AB    MN   arrangememt  ... nareas   longitude  latitude
    area                             ...
    None  200.0  20.0  schlumberger  ...      1  110.486111  26.05174
    >>> twinning (table1, table.table_,  )
           dipole   longitude  latitude  ...  nareas   longitude  latitude
    line1    10.0  110.486111  26.05174  ...     NaN         NaN       NaN
    None      NaN         NaN       NaN  ...     1.0  110.486111  26.05174
    >>> twinning (table1, table.table_, on =['longitude', 'latitude'] )
    Empty DataFrame
    >>> # comments: Empty dataframe appears because, decimal is too large
    >>> # then it considers values longitude and latitude differents
    >>> twinning (table1, table.table_, on =['longitude', 'latitude'], decimals =5 )
        dipole  longitude  latitude  ...  max_depth  ohmic_area  nareas
    0      10  110.48611  26.05174  ...      109.0  690.063003       1
    >>> # Now is able to find existing dataframe with identical closer coordinates.

    """
    from .validator import _is_numeric_dtype

    if str(mode).lower()=='soft':
        d = [ o for  o in d if hasattr (o, '__array__') and hasattr (o, 'columns') ]

    is_same = set ( [ hasattr (o, '__array__')
                     and hasattr (o, 'columns') for o in d ] )

    if len(is_same)!=1 or not list (is_same) [0]:
        types = [ type(o).__name__ for o in d ]
        raise TypeError (
            f"Expect DataFrame. Got {smart_format(types)}")

    same_len = [len(o) for o in d]

    if len( set(same_len)) !=1 or not list(set(same_len))[0]:
        if not coerce:
            raise ValueError(
                f"Data must be a consistent size. Got {smart_format(same_len)}"
                " respectively. Set ``coerce=True`` to truncate the data"
                " to match the shorter data length.")
        # get the shorthest len
        min_index = min (same_len)
        d = [ o.iloc [:min_index, :  ]  for o in d ]

    if on is None:
        if raise_warn:
            warnings.warn("'twin_items' are missing in the data. A simple merge"
                          " along the columns axis should be performed.", stacklevel=2)

        return pd.concat ( d, axis = 1 )

    # parse string
    on= is_iterable(on, exclude_string= True ,
                    transform =True, parse_string= parse_on
                    )

    feature_exist = [
        exist_features(o, on, error = 'ignore'
                       ) for o in d ]

    if ( len( set (feature_exist) )!=1
        or not list(set(feature_exist)) [0]
            ):
        if not force:
            raise ValueError(
                f"Unable to fit the data. Items {smart_format(on)} are"
                f" missing in the data columns. {smart_format(on)} must"
                 " include in the data columns. Please check your data.")

        # seek the value twin_items in the data and use if for all
        dtems =[] ;  repl_twd= None
        for o in d:
            # select one valid data that contain the
            # twin items
            try :
                exist_features ( o, on )
            except :
                pass
            else:
                repl_twd = o [ on ]
                break

        if repl_twd is None:
            raise ValueError("To force data that have not consistent items,"
                             f" at least {smart_format(on)} items must"
                             " be included in one DataFrame.")
        # make add twin_data if
        for o in d :
            try : exist_features(o, on)
            except :
                a = o.copy()
                a [ on ] = repl_twd
            else : a = o.copy ()

            dtems.append(a)
        # reinitialize d
        d = dtems

    # check whether values to merge are numerics
    # if True, use decimals to round values to consider
    # that identic
    # round value if value before performed merges
    # test single data with on
    is_num = _is_numeric_dtype (d[0][on] )
    if is_num:
        decimals = int (_assert_all_types(
            decimals, int, float, objname ='Decimals'))
        d_ = []
        for o in d :
            a = o.copy()
            a[on ] = np.around (o[ on ].values, decimals )
            d_.append (a )
        # not a numerick values so stop
        d =d_

    # select both two
    data = pd.merge (* d[:2], on= on )

    if len(d[2:]) !=0:
        for ii, o in enumerate ( d[2:]) :
            data = pd.merge ( *[data, o] , on = on, suffixes= (
                f"_x{ii+1}", f"_y{ii+1}"))

    return data

def read_worksheets(*data):
    """ Read sheets and returns a list of DataFrames and sheet names.

    Parameters
    -----------
    data: list of str
      A collection of excel sheets files. Read only `.xlsx` files. Any other
      files raises an errors.

    Return
    ------
    data, sheet_names: Tuple of DataFrames and sheet_names
       A collection of DataFrame and sheets names.

    Examples
    -----------
    >>> import os
    >>> from watex.utils.funcutils import read_worksheets
    >>> sheet_file= r'F:\repositories\\watex\\data\\erp\\sheets\\gbalo.xlsx'
    >>> data, snames =  read_worksheets (sheet_file )
    >>> snames
    ['l11', 'l10', 'l02']
    >>> data, snames =  read_worksheets (os.path.dirname (sheet_file))
    >>> snames
    ['l11', 'l10', 'l02', 'l12', 'l13']

    """
    dtem = []
    data = [o for o in data if isinstance ( o, str )]

    for o in data:
        if os.path.isdir (o):
            dlist = os.listdir (o)
            # collect only the excell sheets
            p = [ os.path.join(o, f) for f in dlist if f.endswith ('.xlsx') ]
            dtem .extend(p)
        elif os.path.isfile (o):
            _, ex = os.path.splitext( o)
            if ex == '.xlsx':
                dtem.append(o)

    data = copy.deepcopy(dtem)
    # if no excel sheets is found return None
    if len(data) ==0:
        return None, None

    # make d dict to collect data
    ddict = dict()
    regex = re.compile (r'[$& #@%^!]', flags=re.IGNORECASE)

    for d in data :
        try:
            ddict.update ( **pd.read_excel (d , sheet_name =None))
        except : pass

    #collect stations names
    if len(ddict)==0 :
        raise TypeError("Can'find the data to read.")

    sheet_names = list(map(
        lambda o: regex.sub('_', o).lower(), ddict.keys()))

    data = list(ddict.values ())

    return data, sheet_names

def key_checker (
    keys: str , /,
    valid_keys:List[str, ...],
    regex:re = None,
    pattern:str = None ,
    deep_search:bool =...
    ):
    r"""check whether a give key exists in valid_keys and return a list if
    many keys are found.

    Parameters
    -----------
    keys: str, list of str
       Key value to find in the valid_keys

    valid_keys: list
       List of valid keys by default.

    regex: `re` object,
        Regular expresion object. the default is::

            >>> import re
            >>> re.compile (r'[_#&*@!_,;\s-]\s*', flags=re.IGNORECASE)

    pattern: str, default = '[_#&*@!_,;\s-]\s*'
        The base pattern to split the text into a columns

    deep_search: bool, default=False
       If deep-search, the key finder is no sensistive to lower/upper case
       or whether a numeric data is included.

       .. versionadded:: 0.2.5

    Returns
    --------
    keys: str, list ,
      List of keys that exists in the `valid_keys`.

    Examples
    --------

    >>> from watex.utils.funcutils import key_checker
    >>> key_checker('h502', valid_keys= ['h502', 'h253','h2601'])
    Out[68]: 'h502'
    >>> key_checker('h502+h2601', valid_keys= ['h502', 'h253','h2601'])
    Out[69]: ['h502', 'h2601']
    >>> key_checker('h502 h2601', valid_keys= ['h502', 'h253','h2601'])
    Out[70]: ['h502', 'h2601']
    >>> key_checker(['h502',  'h2601'], valid_keys= ['h502', 'h253','h2601'])
    Out[73]: ['h502', 'h2601']
    >>> key_checker(['h502',  'h2602'], valid_keys= ['h502', 'h253','h2601'])
    UserWarning: key 'h2602' is missing in ['h502', 'h2602']
    Out[82]: 'h502'
    >>> key_checker(['502',  'H2601'], valid_keys= ['h502', 'h253','h2601'],
                    deep_search=True )
    Out[57]: ['h502', 'h2601']

    """
    deep_search, =ellipsis2false(deep_search)

    _keys = copy.deepcopy(keys)
    valid_keys = is_iterable(valid_keys , exclude_string =True, transform =True )
    if isinstance ( keys, str):
        pattern = pattern or r'[_#&@!_+,;\s-]\s*'
        keys = str2columns (keys, regex = regex , pattern=pattern )
    # If iterbale object , save obj
    # to improve error
    kkeys = copy.deepcopy(keys)
    if deep_search:
        keys = key_search(
            keys,
            default_keys= valid_keys,
            deep=True,
            raise_exception= True,
            regex =regex,
            pattern=pattern
            )
        return keys[0] if len(keys)==1 else keys
    # for consistency
    keys = [ k for k in keys if ''.join(
        [ str(i) for i in valid_keys] ).find(k)>=0 ]
    # assertion error if key does not exist.
    if len(keys)==0:
        verb1, verb2 = ('', 'es') if len(kkeys)==1 else ('s', '')
        msg = (f"key{verb1} {_keys!r} do{verb2} not exist."
               f" Expect {smart_format(valid_keys, 'or')}")
        raise KeyError ( msg )

    if len(keys) != len(kkeys):
        miss_keys = is_in_if ( kkeys, keys , return_diff= True , error ='ignore')
        miss_keys, verb = (miss_keys[0], 'is') if len( miss_keys) ==1 else (
            miss_keys, 'are')
        warnings.warn(f"key{'' if verb=='is' else 's'} {miss_keys!r} {verb}"
                      f" missing in {_keys}", stacklevel=2)
    keys = keys[0] if len(keys)==1 else keys

    return keys

def random_sampling (
    d, / ,
    samples:int = None  ,
    replace:bool = False ,
    random_state:int = None,
    shuffle=True,
    ):
    """ Sampling data.

    Parameters
    ----------
    d: {array-like, sparse matrix} of shape (n_samples, n_features)
      Data for sampling, where `n_samples` is the number of samples and
      `n_features` is the number of features.
    samples: int,optional
       Ratio or number of items from axis to return.
       Default = 1 if `samples` is ``None``.

    replace: bool, default=False
       Allow or disallow sampling of the same row more than once.

    random_state: int, array-like, BitGenerator, np.random.RandomState, \
        np.random.Generator, optional
       If int, array-like, or BitGenerator, seed for random number generator.
       If np.random.RandomState or np.random.Generator, use as given.

    shuffle:bool, default=True
       Shuffle the data before sampling

    Returns
    ----------
    d: {array-like, sparse matrix} of shape (n_samples, n_features)
    samples data based on the given samples.

    Examples
    ---------
    >>> from watex.utils.funcutils import random_sampling
    >>> from watex.datasets import load_hlogs
    >>> data= load_hlogs().frame
    >>> random_sampling( data, samples = 7 ).shape
    (7, 27)
    """

    n= None ; is_percent = False
    orig= copy.deepcopy(samples )
    if not hasattr(d, "__iter__"):
        d = is_iterable(d, exclude_string= True, transform =True )

    if (
            samples is None
            or str(samples) in ('1', '*')
            ):
        samples =1.

    if "%" in str(samples):
        samples = samples.replace ("%", '')
        is_percent=True
    # assert value for consistency.
    try:
        samples = float( samples)
    except:
        raise TypeError("Wrong value for 'samples'. Expect an integer."
                        f" Got {type (orig).__name__!r}")

    if samples <=1 or is_percent:
        samples  = assert_ratio(
            samples , bounds = (0, 1), exclude_value= 'use lower bound',
            in_percent= True )

        n = int ( samples * ( d.shape[0] if scipy.sparse.issparse(d)
                 else len(d)))
    else:
        # data frame
        n= int(samples)

    # reset samples and use number of samples instead
    samples =None
    # get the total length of d
    dlen = ( d.shape[0] if scipy.sparse.issparse(d) else len(d))
    # if number is greater than the length
    # block to the length so to retrieve all
    # value no matter the arrangement.
    if n > dlen:
        n = dlen
    if hasattr (d, 'columns') or hasattr (d, 'name'):
        # data frame
        return d.sample ( n= n , frac=samples , replace = replace ,
                         random_state = random_state
                     ) if shuffle else d.iloc [ :n , ::]

    np.random.seed ( random_state)
    if scipy.sparse.issparse(d) :
        if scipy.sparse.isspmatrix_coo(d):
            warnings.warn("coo_matrix does not support indexing. Conversion"
                          " should be performed in CSR matrix", stacklevel=2)
            d = d.tocsr()

        return d [ np.random.choice(
            np.arange(d.shape[0]), n, replace=replace )] if shuffle else d [
                [ i for i in range (n)]]

        #d = d[idx ]

    # manage the data
    if not hasattr(d, '__array__'):
        d = np.array (d )

    idx = np.random.randint( len(d), size = n ) if shuffle else [ i for i in range(n)]
    if len(d.shape )==1: d =d[idx ]
    else: d = d[idx , :]

    return d


def make_obj_consistent_if (
        item= ... , default = ..., size =None, from_index: bool =True ):
    """Combine default values to item to create default consistent iterable
    objects.

    This is valid if  the size of item does not fit the number of
    expected iterable objects.

    Parameters
    ------------
    item : Any
       Object to construct it default values

    default: Any
       Value to hold in the case the items does not match the size of given items

    size: int, Optional
      Number of items to return.

    from_index: bool, default=True
       make an item size to match the exact size of given items

    Returns
    -------
       item: Iterable object that contain default values.

    Examples
    ----------
    >>> from watex.utils.funcutils import make_obj_consistent_if
    >>> from watex.exlib import SVC, LogisticRegression, XGBClassifier
    >>> classifiers = ["SVC", "LogisticRegression", "XGBClassifier"]
    >>> classifier_names = ['SVC', 'LR']
    >>> make_obj_consistent_if (classifiers, default = classifier_names )
    ['SVC', 'LogisticRegression', 'XGBClassifier']
    >>> make_obj_consistent_if (classifier_names, from_index =False  )
    ['SVC', 'LR']
    >>> >>> make_obj_consistent_if ( classifier_names,
                                     default= classifiers, size =3 ,
                                     from_index =False  )
    ['SVC', 'LR', 'SVC']

    """
    if default==... or None : default =[]
    # for consistency
    default = list( is_iterable (default, exclude_string =True,
                                 transform =True ) )

    if item not in ( ...,  None) :
         item = list( is_iterable( item , exclude_string =True ,
                                  transform = True ) )
    else: item = []

    item += default[len(item):] if from_index else default

    if size is not None:
        size = int (_assert_all_types(size, int, float,
                                      objname = "Item 'size'") )
        item = item [:size]

    return item


def replace_data(
    X, y =None,
    n_times: int  = 1,
    axis = 0,
    reset_index :bool =...
    ):
    """ Replace items in data :math:`n` times

    Parameters
    ----------
    X: Arraylike 1D or pd.DataFrame
      Data to replace. Note Sparse matrices is not allowed. Use
      :func:`random_sampling` instead.
    y: Arraylike 1d.
       Preferably one dimensional data.

    n_times: int,
      Number of times all items should be replaced in data.

    reset_index: bool, default=False.
      If ``True`` and dataframe,Index is reset and dropped.

    Returns
    --------
    X or (X, y)  : Tuple of data replaced
       Tuple is returned if y is passed.

    Examples
    ---------
    >>> import numpy as np
    >>> from watex.utils.funcutils import replace_data
    >>> X, y = np.random.randn ( 7, 2 ), np.arange(7)
    >>> X.shape, y.shape
    ((7, 2), (7,))
    >>> X_new, y_new = replace_data (X, y, n_times =10 )
    >>> X_new.shape , y_new.shape
    Out[158]: ((70, 2), (70,))
    """

    n = n_times or 1
    n= int (_assert_all_types(n, int, float, objname ='n_times'))

    def concat_data ( ar ,):
        if hasattr ( ar, 'columns') or hasattr ( ar, 'series'):
            d = pd.concat (  [ ar for i in range(n)], axis = axis )
            if reset_index:
                d.reset_index (drop =True, inplace =True  )

        else: d = np.concatenate ([ ar for i in range(n)], axis = axis  )
        return d

    X = is_iterable(X, exclude_string =True, transform = True )

    if y is not None:
        y = is_iterable( y, exclude_string= True, transform= True)

    if not hasattr (X, '__array__'):
        X = np.array(X)

    if not hasattr (y, '__array__') and y is not None :
        y = np.array(y)

    return concat_data ( X) if y is None else (
            concat_data( X) , concat_data(y))

def convert_value_in (v, /, unit ='m'):
    """Convert value based on the reference unit.

    Parameters
    ------------
    v: str, float, int,
      value to convert
    unit: str, default='m'
      Reference unit to convert value in. Default is 'meters'. Could be
      'kg' or else.

    Returns
    -------
    v: float,
       Value converted.

    Examples
    ---------
    >>> from watex.utils.funcutils import convert_value_in
    >>> convert_value_in (20)
    20.0
    >>> convert_value_in ('20mm')
    0.02
    >>> convert_value_in ('20kg', unit='g')
    20000.0
    >>> convert_value_in ('20')
    20.0
    >>> convert_value_in ('20m', unit='g')
    ValueError: Unknwon unit 'm'...
    """
    c= { 'k':1e3 ,
        'h':1e2 ,
        'dc':1e1 ,
        '':1e0 ,
        'd':1e-1,
        'c':1e-2 ,
        'm':1e-3
        }
    c = {k +str(unit).lower(): v for k, v in c.items() }

    v = str(v).lower()

    regex = re.findall(r'[a-zA-Z]', v)

    if len(regex) !=0:
        unit = ''.join( regex )
        v = v.replace (unit, '')

    if unit not in c.keys():
        raise ValueError (
            f"Unknwon unit {unit!r}. Expect {smart_format(c.keys(), 'or' )}."
            f" Or rename the `unit` parameter maybe to {unit[-1]!r}.")

    return float ( v) * (c.get(unit) or 1e0)

def split_list(lst:List[Any, ...],/,  val:int, fill_value:Any=None ):
    """Module to extract a slice of elements from the list

    Parameters
    ------------
    lst: list,
      List composed of item elements
    val: int,
      Number of item to group by default.

    Returns
    --------
    group with slide items

    Examples
    --------
    >>> from watex.utils.funcutils import split_list
    >>> lst = [1, 2, 3, 4, 5, 6, 7, 8]
    >>> val = 3
    >>> print(split_list(lst, val))
    [[1, 2, 3], [4, 5, 6], [7, 8]]

    """

    lst = is_iterable(lst , exclude_string =True , transform =True )
    val = int ( _assert_all_types(val, int, float ))
    try:
        sl= [list(group) for key, group in itertools.groupby(
                lst, lambda x: (x-1)//val)]
    except:
        # when string is given
        sl= list(itertools.zip_longest(
            *(iter(lst),)*val,fillvalue =fill_value),)
    return sl

def key_search (
    keys: str, /,
    default_keys: Text | List[str],
    parse_keys: bool=True,
    regex :re=None,
    pattern :str=None,
    deep: bool =...,
    raise_exception:bool=...,
    ):
    r"""Find key in a list of default keys and select the best match.

    Parameters
    -----------
    keys: str or list
       The string or a list of key. When multiple keys is passed as a string,
       use the space for key separating.

    default_keys: str or list
       The likehood key to find. Can be a litteral text. When a litteral text
       is passed, it is better to provide the regex in order to skip some
       character to parse the text properly.

    parse_keys: bool, default=True
       Parse litteral string using default `pattern` and `regex`.

       .. versionadded:: 0.2.7

    regex: `re` object,
        Regular expresion object. Regex is important to specify the kind
        of data to parse. the default is::

            >>> import re
            >>> re.compile (r'[_#&*@!_,;\s-]\s*', flags=re.IGNORECASE)

    pattern: str, default = '[_#&*@!_,;\s-]\s*'
        The base pattern to split the text into a columns. Pattern is
        important especially when some character are considers as a part of
        word but they are not a separator. For example a data columns with
        a name `'DH_Azimuth'`, if a pattern is not explicitely provided,
        the default pattern will parse as two separated word which is far
        from the expected results.

    deep: bool, default=False
       Not sensistive to uppercase.

    raise_exception: bool, default=False
       raise error when key is not find.

    Return
    -------
    list: list of valid keys or None if not find ( default)

    Examples
    ---------
    >>> from watex.utils.funcutils import key_search
    >>> key_search('h502-hh2601', default_keys= ['h502', 'h253','HH2601'])
    Out[44]: ['h502']
    >>> key_search('h502-hh2601', default_keys= ['h502', 'h253','HH2601'],
                   deep=True)
    Out[46]: ['h502', 'HH2601']
    >>> key_search('253', default_keys= ("I m here to find key among h502,
                                             h253 and HH2601"))
    Out[53]: ['h253']
    >>> key_search ('east', default_keys= ['DH_East', 'DH_North']  , deep =True,)
    Out[37]: ['East']
    key_search ('east', default_keys= ['DH_East', 'DH_North'],
                deep =True,parse_keys= False)
    Out[39]: ['DH_East']
    """
    deep, raise_exception, parse_keys = ellipsis2false(
        deep, raise_exception, parse_keys)
    # make a copy of original keys

    kinit = copy.deepcopy(keys)
    if parse_keys:
        if is_iterable(keys , exclude_string= True ):
            keys = ' '.join ( [str(k) for k in keys ])
             # for consisteny checker
        pattern = pattern or r'[#&@!_+,;\s-]\s*'
        keys = str2columns ( keys , regex = regex , pattern = pattern )

        if is_iterable ( default_keys , exclude_string=True ):
            default_keys = ' '. join ( [ str(k) for k in default_keys ])
            # make a copy
        default_keys =  str2columns(
            default_keys, regex =regex , pattern = pattern )
    else :
        keys = is_iterable(
        keys, exclude_string = True, transform =True )
        default_keys = is_iterable (
            default_keys, exclude_string=True, transform =True )

    dk_init = copy.deepcopy(default_keys )
    # if deep convert all keys to lower
    if deep:
        keys= [str(it).lower() for it in keys  ]
        default_keys = [str(it).lower() for it in default_keys  ]

    valid_keys =[]
    for key in keys :
        for ii, dkey in enumerate (default_keys) :
            vk = re.findall(rf'\w*{key}\w*', dkey)
            # rather than rf'\b\w*{key}\w*\b'
            # if deep take the real values in defaults keys.
            if len(vk) !=0:
                if deep: valid_keys.append( dk_init[ii] )
                else:valid_keys.extend( vk)
                break
    if ( raise_exception
        and len(valid_keys)==0
        ):
        kverb ='s' if len(kinit)> 1 else ''
        raise KeyError (f"key{kverb} {kinit!r} not found."
                       f" Expect {smart_format(dk_init, 'or')}")
    return None if len(valid_keys)==0 else valid_keys

def repeat_item_insertion(text, /, pos, item ='', fill_value=''):
    """ Insert character in  text according from it position.

    Parameters
    -----------
    v: text
       Text
    pos: int
      position where the item must be insert.
    item: str,
      Item to insert at each position.
    fill_value: str,
      Does nothing special; fill the the last position.
    Returns
    --------
    text: str,
      New construct object.

    Examples
    ----------
    >>> from watex.utils.funcutils import repeat_item_insertion
    >>> repeat_item_insertion ( '0125356.45', pos=2, item=':' )
    Out[65]: '01:25:35:6.45'
    >>> repeat_item_insertion ( 'Function inserts car in text.', pos=10, item='TK' )
    Out[69]: 'Function iTKnserts carTK in text.'
    """
    pos= _assert_all_types(pos, int, float,
                           objname=f'Position for {item} insertion')
    # for consistency
    lst = list( str(text))
    # checher whether there is a decimal then remove it
    dec_part=[]
    ent_part= lst
    for i, it in enumerate ( lst)  :
        if it =='.':
            ent_part, dec_part  = lst [:i],  lst[i:]
            break
    # now split list
    value = split_list(ent_part, val= pos , fill_value=fill_value)
    #value = split_list ( ent_part, 2)
    #[[1, 2, 3], [4, 5, 6], [7, 8]]
    join_lst= list (map ( lambda s : ''.join( s), value))
    #[123, 456, 78]
    #join with mark
    return f'{str(item)}'.join(join_lst) +''.join(dec_part)

def numstr2dms (
    sdigit: str, /,
    sanitize: bool=True,
    func: F=None,
    args: tuple=(),
    regex: re=None,
    pattern: str=None,
    return_values: bool=...,
    **kws
    ):
    """ Convert numerical digit string to DD:MM:SS

    Note that the any string digit for Minutes and seconds must be composed
    of two values i.e the function accepts at least six digits, otherwise an
    error occurs. For instance the value between [0-9] must be prefixed by 0
    beforehand. Here is an example for designating 1degree-1min-1seconds::

        sdigit= 1'1'1" --> 01'01'01 or 010101

    where ``010101`` is the right arguments for ``111``.

    Parameters
    -----------
    sdigit: str,
      Digit string composing of unique values.
    func: Callable,
      Function uses to parse digit. Function must return a string values.
      Any other values should be convert to str

    args: tuple
      Function `func` positional arguments

    regex: `re` object,
        Regular expresion object. Regex is important to specify the kind
        of data to parse. the default is::

            >>> import re
            >>> re.compile (r'[_#&@!+,;:"\'\\s-]\\s*', flags=re.IGNORECASE)

    pattern: str, default = '[_#&@!+,;:"\'\\s-]\\s*'
      Specific pattern for sanitizing sdigit. For instance remove undesirable
      non-character.

    sanitize:bool=default=True
       Remove undesirable character using the default argument of `pattern`
       parameter.

    return_values: bool, default=False,
       return the DD:MM:SS into a tuple of (DD,MM,SS)

    Returns
    -------
    sdigit/tuple: str, tuple
      DD:MM:SS or tuple of ( DD, MM, SS)

    Examples
    --------
    >>> from watex.utils.funcutils import numstr2dms
    >>> numstr2dms ("1134132.08")
    Out[17]: '113:41:32.08
    >>> numstr2dms ("13'41'32.08")
    Out[18]: '13:41:32.08'
    >>> numstr2dms ("11:34:13:2.08", return_values=True)
    Out[19]: (113.0, 41.0, 32.08)

    """
    # remove any character from the string digit
    if return_values is ...:return_values=False
    sdigit= str(sdigit)

    if sanitize:
        pattern = pattern or '[_#&@!+,;:"\'\\s-]\\s*'
        sdigit = re.sub(pattern , "", str(sdigit), flags=re.IGNORECASE)

    try : float (sdigit)
    except: raise ValueError ("Wrong value. Expects a string-digit or digit."
                              f" Got {sdigit!r}")
    if callable (func):
        sdigit= func (sdigit, *args, **kws )

    # In the case there is'
    decimal ='0'
    # remove decimal
    sdigit_list  = str(sdigit).split(".")

    if len(sdigit_list)==2:
        sdigit, decimal =sdigit_list

    if len(sdigit) < 6:
        raise ValueError(f"DMS expects at list six digits(DD:MM:SS)."
                         f" Got {sdigit!r}")

    sec , sdigit = sdigit[-2:] , sdigit [:-2]
    mm , sdigit = sdigit[-2:], sdigit [:-2]
    deg = sdigit # the remain part
    # conca second ecimal
    sec +=f".{decimal}"

    return tuple (map ( float, [deg, mm, sec]) ) if return_values \
        else ':'.join([deg, mm, sec])


def storeOrwritehdf5 (
    d, /,
    key:str= None,
    mode:str='a',
    kind: str=None,
    path_or_buf:str= None,
    encoding:str="utf8",
    csv_sep: str=",",
    index: bool=...,
    columns: str |List[Any, ...]=None,
    sanitize_columns:bool=...,
    func: F= None,
    args: tuple=(),
    applyto: str|List[Any, ...]=None,
    **func_kwds,
    )->None|DataFrame:
    """ Store data to hdf5 or write data to csv file.

    Note that by default, the data is not store nor write and
    return data if frame or transform the Path-Like object to data frame.

    Parameters
    -----------
    d: Dataframe, shape (m_samples, n_features)
        data to store or write or sanitize.
    key:str
       Identifier for the group in the store.

    mode: {'a', 'w', 'r+'}, default 'a'
       Mode to open file:

       - 'w': write, a new file is created (an existing file with the
                                          same name would be deleted).
       - 'a': append, an existing file is opened for reading and writing,
         and if the file does not exist it is created.
       - 'r+': similar to 'a', but the file must already exist.

    kind: str, {'store', 'write', None} , default=None
       Type of task to perform:

       - 'store': Store data to hdf5
       - 'write': export data to csv file.
       - None: construct a dataframe if array is passed or sanitize it.

    path_or_buf: str or pandas.HDFStore, or str, path object, file-like \
        object, or None, default=None
       File path or HDFStore object. String, path object
       (implementing os.PathLike[str]), or file-like object implementing
       a write() function. If ``write=True`` and  None, the result is returned
       as a string. If a non-binary file object is passed, it should be
       opened with newline=" ", disabling universal newlines. If a binary
       file object is passed, mode might need to contain a 'b'.

    encoding: str, default='utf8'
       A string representing the encoding to use in the output file,
       Encoding is not supported if path_or_buf is a non-binary file object.

    csv_sep: str, default=',',
       String of length 1. Field delimiter for the output file.

    index: bool, index =False,
       Write data to csv with index or not.

    columns: list of str, optional
        Usefull to create a dataframe when array is passed. Be aware to fit
        the number of array columns (shape[1])

    sanitize_columns: bool, default=False,
       remove undesirable character in the data columns using the default
       argument of `regex` parameters and fill pattern to underscore '_'.
       The default regex implementation is::

           >>> import re
           >>> re.compile (r'[_#&.)(*@!,;\\s-]\\s*', flags=re.IGNORECASE)

    func: callable, Optional
       A custom sanitizing function and apply to each columns of the dataframe.
       If provide, the expected columns must be listed to `applyto` parameter.

    args: tuple, optional
       Positional arguments of the sanitizing columns

    applyto: str or list of str, Optional
       The list of columns to apply the function ``func``. To apply the
       function to all columns, use the ``*`` instead.

    func_kwds: dict,
       Keywords arguments of the sanitizing function ``func``.

    Return
    -------
    None or d: None of dataframe.
      returns None if `kind` is set to ``write`` or ``store`` otherwise
      return the dataframe.

    Examples
    --------
    >>> from watex.utils.funcutils import storeOrwritehdf5
    >>> from watex.datasets import load_bagoue
    >>> data = load_bagoue().frame
    >>> data.geol[:5]
    0    VOLCANO-SEDIM. SCHISTS
    1                  GRANITES
    2                  GRANITES
    3                  GRANITES
    4          GEOSYN. GRANITES
    Name: geol, dtype: object
    >>> data = storeOrwritehdf5 ( data, sanitize_columns = True)
    >>> data[['type', 'geol', 'shape']] # put all to lowercase
      type                    geol shape
    0   cp  volcano-sedim. schists     w
    1   ec                granites     v
    2   ec                granites     v
    >>> # compute using func
    >>> def test_func ( a, times  , to_percent=False ):
            return ( a * times / 100)   if to_percent else ( a *times )
    >>> data.sfi[:5]
    0    0.388909
    1    1.340127
    2    0.446594
    3    0.763676
    4    0.068501
    Name: sfi, dtype: float64
    >>> d = storeOrwritehdf5 ( data,  func = test_func, args =(7,), applyto='sfi')
    >>> d.sfi[:5]
    0    2.722360
    1    9.380889
    2    3.126156
    3    5.345733
    4    0.479507
    Name: sfi, dtype: float64
    >>> storeOrwritehdf5 ( data,  func = test_func, args =(7,),
                          applyto='sfi', to_percent=True).sfi[:5]
    0    0.027224
    1    0.093809
    2    0.031262
    3    0.053457
    4    0.004795
    Name: sfi, dtype: float64
    >>> # write data to hdf5 and outputs to current directory
    >>> storeOrwritehdf5 ( d, key='test0', path_or_buf= 'test_data.h5',
                          kind ='store')
    >>> # export data to csv
    >>> storeOrwritehdf5 ( d, key='test0', path_or_buf= 'test_data',
                          kind ='export')
    """
    kind= key_search (str(kind), default_keys=(
        "none", "store", "write", "export", "tocsv"),
        raise_exception=True , deep=True)[0]

    kind = "export" if kind in ('write', 'tocsv') else kind

    if sanitize_columns is ...:
        sanitize_columns=False
    d = to_numeric_dtypes(d, columns=columns,sanitize_columns=sanitize_columns,
                          fill_pattern='_')

    # get categorical variables
    if ( sanitize_columns
        or func is not None
        ):
        d, _, cf = to_numeric_dtypes(d, return_feature_types= True )
        #( strip then pass to lower case all non-numerical data)
        # for minimum sanitization
        for cat in cf :
            d[cat]= d[cat].str.lower()
            d[cat]= d[cat].str.strip()

    if func is not None:
        if not callable(func):
            raise TypeError(
                f"Expect a callable for `func`. Got {type(func).__name__!r}")

        if applyto is None:
            raise ValueError("Need to specify the data column to apply"
                             f"{func.__name__!r} to.")

        applyto = is_iterable( applyto, exclude_string=True,
                               transform =True ) if applyto !="*" else d.columns
        # check whether the applyto columns are in data columns
        exist_features(d, applyto)

        # map each colum
        for col in applyto:
            d [col]=d[col].apply( func, args=args, **func_kwds )

    # store in h5 file.
    if kind=='store':
        if path_or_buf is None:
            print("Destination file is missing. Use 'data.h5' instead outputs"
                  f" in the current directory {os.getcwd()}")
            path_or_buf= 'data.h5'

        d.to_hdf ( path_or_buf , key =key, mode =mode )
    # export to csv file
    if kind=="export":
        d.to_csv(path_or_buf, encoding = encoding  , sep=csv_sep ,
                 index =False if index is ... else index   )

    return d if kind not in ("store", "export") else None

def ellipsis2false( *parameters , default_value: Any=False ):
    """ Turn all parameter arguments to False if ellipsis.

    Note that the output arguments must be in the same order like the
    positional arguments.

    :param parameters: tuple
       List of parameters
    :param default_value: Any,
       Value by default that might be take the ellipsis.
    :return: tuple, same list of parameters passed ellipsis to
       ``default_value``. By default, it returns ``False``. For a single
       parameters, uses the trailing comma for collecting the parameters

    :example:
        >>> from watex.utils.funcutils import ellipsis2false
        >>> var, = ellipsis2false (...)
        >>> var
        False
        >>> data, sep , verbose = ellipsis2false ([2,3, 4], ',', ...)
        >>> verbose
        False
    """
    return tuple (  default_value  if param is  ... else param
                    for param in parameters )













































