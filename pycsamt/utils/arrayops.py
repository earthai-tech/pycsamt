# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""
Array and list manipulation utilities.
"""
from __future__ import annotations 
import re
from typing import Iterable, Union, Optional, Sequence, Any
import warnings
import matplotlib.pyplot as plt 
import numpy as np
import pandas as pd 
import scipy.interpolate as spi

from .validation import _assert_all_types, _is_numeric_dtype
from .text import ( 
    str2columns, 
    listing_items_format, 
)
from .cleaner import sanitize_frame_cols


__all__ = [
    'concat_array_from_list',
    'is_iterable',
    'reshape',
    'assert_xy_in', 
    'frameify', 
    'interpolate_grid', 
    'drop_nan_in', 
]

def concat_array_from_list(
    list_of_arrays: Iterable[Any],
    concat_axis: int = 1
) -> np.ndarray:
    """
    Concatenate a list of array-like objects into a 2D array,
    padding missing or unequal lengths with NaN.

    Parameters
    ----------
    list_of_arrays : iterable
        Iterable of array-like items (lists, arrays, None).
        Scalars or None will be broadcast to length of others
        when possible.
    concat_axis : {0, 1}, default 1
        Axis along which to concatenate:
        - 1: each input forms a column, result shape (n, m)
        - 0: each input forms a row, result shape (m, n)

    Returns
    -------
    result : ndarray
        2D array of shape (n, m) or (m, n), with NaNs for padding.

    Raises
    ------
    TypeError
        If inputs are not iterable or items are not array-like.
    ValueError
        If `concat_axis` is not 0 or 1.

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.utils.arrayops import concat_array_from_list
    >>> a = np.arange(3)
    >>> b = [10, 20, 30]
    >>> concat_array_from_list([a, b])
    array([[ 0., 10.],
           [ 1., 20.],
           [ 2., 30.]])
    >>> concat_array_from_list([a, None, [5]], concat_axis=1)
    array([[ 0., nan,  5.],
           [ 1., nan, nan],
           [ 2., nan, nan]])
    """
    # Validate axis
    axis = int(_assert_all_types(concat_axis, int, objname='concat_axis'))
    if axis not in (0, 1):
        raise ValueError(f"concat_axis must be 0 or 1, got {axis}")
    # Coerce to list
    try:
        seq = list(list_of_arrays)
    except Exception:
        raise TypeError("`list_of_arrays` must be iterable")
    n_items = len(seq)
    if n_items == 0:
        return np.empty((0, 0))
    # Convert elements to 1D arrays, track lengths
    arrs = []
    lengths = []
    for item in seq:
        if item is None:
            arr = np.array([np.nan])
        else:
            try:
                arr = np.asarray(item)
            except Exception:
                raise TypeError(
                    f"Item {item!r} is not array-like"
                )
        if arr.ndim == 0:
            arr = arr.reshape(1)
        elif arr.ndim > 1:
            # flatten higher dims
            arr = arr.ravel()
        arrs.append(arr.astype(float))
        lengths.append(arr.shape[0])
    # Determine target length (max)
    max_len = max(lengths)
    # Pad arrays to max_len
    padded = []
    for arr, length in zip(arrs, lengths):
        if length < max_len:
            pad = np.full((max_len - length,), np.nan)
            arr = np.concatenate([arr, pad])
        padded.append(arr)
    # Reshape for concatenation
    mats = []
    for arr in padded:
        if axis == 1:
            mats.append(arr.reshape(-1, 1))
        else:
            mats.append(arr.reshape(1, -1))
    # Concatenate
    result = np.concatenate(mats, axis=axis)
    return result


def is_iterable(
    y: Any,
    exclude_string: bool = False,
    transform: bool = False,
    parse_string: bool = False
) -> Union[bool, list]:
    """
    Check if `y` is iterable, with options to transform or parse strings.

    Parameters
    ----------
    y : any
        Object to test.
    exclude_string : bool, default False
        If True, treat strings as non-iterable.
    transform : bool, default False
        If True, return iterable form of `y` (list), else return bool.
    parse_string : bool, default False
        If True and `y` is a string, parse into list via `str2columns`.
        Requires `transform=True`.

    Returns
    -------
    bool or list
        If `transform=False`, returns whether `y` is iterable.
        If `transform=True`, returns list(y) or parsed list.

    Raises
    ------
    ValueError
        If `parse_string=True` but `transform=False`.
    TypeError
        If string parsing fails or `y` is not transformable.
    """
    _assert_all_types(exclude_string, bool, objname='exclude_string')
    _assert_all_types(transform, bool, objname='transform')
    _assert_all_types(parse_string, bool, objname='parse_string')
    if parse_string and not transform:
        raise ValueError("parse_string requires transform=True")
    if isinstance(y, str) and parse_string:
        try:
            y = str2columns(y)
        except Exception as e:
            raise TypeError(f"Error parsing string: {e}")
    is_str = isinstance(y, str)
    is_it = hasattr(y, '__iter__') and not (exclude_string and is_str)
    if not transform:
        return is_it
    try:
        return list(y) if is_it else [y]
    except Exception:
        return [y]


def reshape(
    arr: Any,
    axis: Optional[int] = None
) -> np.ndarray:
    """
    Reshape 1D or singleton 2D arrays to desired orientation.

    Parameters
    ----------
    arr : array-like
        Input array with ndim 1 or 2.
    axis : {None, 0, 1}, optional
        Desired axis for 1D->2D:
        - None: squeeze singleton dims to 1D if possible.
        - 0: return column vector (n,1).
        - 1: return row vector (1,n).

    Returns
    -------
    ndarray
        Reshaped array.

    Raises
    ------
    ValueError
        If `arr` has ndim > 2 or invalid `axis`.
    """
    a = np.asarray(arr)
    if a.ndim > 2:
        raise ValueError(
            f"Input ndim must be <= 2, got {a.ndim}"
        )
    if axis not in (None, 0, 1):
        raise ValueError(
            f"axis must be None, 0, or 1, got {axis!r}"
        )
    # 1D case
    if a.ndim == 1:
        n = a.shape[0]
        if axis is None:
            return a
        if axis == 0:
            return a.reshape(n, 1)
        # axis == 1
        return a.reshape(1, n)
    # 2D case
    n, m = a.shape
    # squeeze when one dim is 1
    if axis is None:
        if n == 1:
            return a.reshape(m,)
        if m == 1:
            return a.reshape(n,)
        return a
    # axis specified for singleton dims
    if axis == 0:
        # output dims (n,1)
        if m != 1:
            return a
        return a.reshape(n, 1)
    # axis == 1
    if axis == 1:
        if n != 1:
            return a
        return a.reshape(1, m)
    return a

def frameify(
    data: Any,
    *,
    columns: Optional[Sequence[str]] = None,
    return_feature_types: bool = False,
    missing_values: Any = np.nan,
    drop_nan_columns: bool = True,
    how: str = 'all',
    sanitize_columns: bool = False,
    regex: Optional[Union[str, re.Pattern]] = None,
    fill_pattern: str = '_',
    reset_index: bool = False,
    drop_index: bool = True,
    pop_cat_features: bool = False,
    verbose: bool = False
) -> Union[pd.DataFrame, tuple]:
    """ Convert array to dataframe and coerce arguments to appropriate dtypes. 
    
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

    # pass ellipsis argument to False 

    if not is_iterable (data, exclude_string=True): 
        raise TypeError(f"Expect array. Got {type (data).__name__!r}")

    if hasattr ( data, '__array__') and hasattr ( data, 'columns'): 
        df = data.copy()
        if columns is not None: 
            if verbose: 
                print("Dataframe is passed. Columns should be replaced.")
            df =pd.DataFrame ( np.array ( data), columns =columns )
            
    else: 
        df = pd.DataFrame (data, columns =columns  ) 
        
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


def assert_xy_in(
    x: Union[str, Iterable[Any], Any],
    y: Union[str, Iterable[Any], Any],
    *,
    data: Optional[pd.DataFrame] = None,
    asarray: bool = True,
    to_frame: bool = False,
    columns: Optional[Sequence[str]] = None,
    xy_numeric: bool = False,
    dropna: bool = False
) -> Union[tuple, pd.DataFrame]:
    """
    Validate and extract paired x, y data from arrays or DataFrame.

    Parameters
    ----------
    x, y : str or array-like
        Column name(s) in `data` or sequence of values.
    data : pd.DataFrame, optional
        DataFrame containing named columns for x and y.
    asarray : bool, default=True
        Return NumPy arrays if True; else pandas Series.
    to_frame : bool, default=False
        Return a new DataFrame with x/y columns.
    columns : sequence of str, optional
        Names for output DataFrame columns (length == 2).
    xy_numeric : bool, default=False
        Coerce x and y values to numeric dtype.
    dropna : bool, default=False
        Drop any pairs with missing values.

    Returns
    -------
    tuple or DataFrame
        x and y as arrays or Series, or a DataFrame.

    Raises
    ------
    TypeError
        If string column name given but `data` is missing.
    KeyError
        If specified column not found in `data`.
    ValueError
        If x and y lengths mismatch or invalid `columns`.
    """

    # convert input to pandas Series
    def _get_series(val):
        # If name given, fetch from DataFrame
        if isinstance(val, str):
            if data is None or not hasattr(data, 'columns'):
                raise TypeError(
                    "`data` DataFrame required when x or y is"
                    " string name"
                )
            if val not in data.columns:
                raise KeyError(f"Column {val!r} not in DataFrame")
            return data[val]
        # If already a Series, return as-is
        if isinstance(val, pd.Series):
            return val
        # Otherwise cast sequence or scalar to Series
        try:
            return pd.Series(val)
        except Exception:
            return pd.Series([val])

    # Retrieve x and y Series
    x_ser = _get_series(x)
    y_ser = _get_series(y)

    # Optionally drop missing pairs
    if dropna:
        df_xy = pd.concat([x_ser, y_ser], axis=1).dropna()
        x_ser, y_ser = df_xy.iloc[:, 0], df_xy.iloc[:, 1]

    # Ensure lengths match
    if len(x_ser) != len(y_ser):
        raise ValueError(
            f"Length mismatch: x has {len(x_ser)}, "
            f"y has {len(y_ser)}"
        )

    # Optionally coerce to numeric dtype
    if xy_numeric:
        x_ser = pd.to_numeric(x_ser, errors='raise')
        y_ser = pd.to_numeric(y_ser, errors='raise')

    # Return as DataFrame if requested
    if to_frame:
        cols = (list(columns)
                if columns is not None else ['x', 'y'])
        if len(cols) != 2:
            raise ValueError(
                "`columns` must be sequence of two names"
            )
        return pd.DataFrame({cols[0]: x_ser.values,
                             cols[1]: y_ser.values})

    # Return numpy arrays
    if asarray:
        return x_ser.values, y_ser.values

    # Return pandas Series
    return x_ser, y_ser

def interpolate_grid (
    arr, 
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
        arri = _fill_nan(arri, method ='both ')
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

def _fill_nan(
    arr: Union[Sequence[Any], np.ndarray],
    *,
    method: str = 'ff'
) -> np.ndarray:
    """
    Efficiently forward/backward fill NaNs in array.

    Parameters
    ----------
    arr : sequence or ndarray
        Input 1D or 2D array with NaNs.
    method : {'ff','bf','both'}, default 'ff'
        'ff': forward-fill; 'bf': backward-fill;
        'both': apply ff then bf.

    Returns
    -------
    filled : ndarray
        Array with NaNs filled. Preserves shape.

    Notes
    -----
    - Works row-wise for 2D arrays.
    - 'ff' fills NaNs after valid values;
      'bf' fills before valid values.

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.utils.arrayops import fill_nan
    >>> a = np.array([np.nan,1,np.nan,2])
    >>> fill_nan(a, method='ff')
    array([nan, 1., 1., 2.])
    >>> fill_nan(a, method='both')
    array([1., 1., 1., 2.])
    """

    # ensure array
    arr2 = np.array(arr, dtype=float)
    # reshape 1D to 2D
    if arr2.ndim == 1:
        arr2 = arr2.reshape(1, -1)
    elif arr2.ndim > 2:
        raise ValueError(
            f"fill_nan supports only 1D or 2D arrays; got ndim={arr2.ndim}"
        )
    # normalize method
    m = method.lower().strip()
    if m in ('forward','ff','fwd'):
        m = 'ff'
    elif m in ('backward','bf','bwd'):
        m = 'bf'
    elif m in ('both','fb','bff','ffbf'):
        m = 'both'
    else:
        raise ValueError(
            f"Unknown method {method!r}; choose 'ff','bf', or 'both'"
        )
    # define forward fill
    def _ffill(a):
        mask = np.isnan(a)
        idx = np.where(~mask, np.arange(a.shape[1]), 0)
        np.maximum.accumulate(idx, axis=1, out=idx)
        return a[np.arange(a.shape[0])[:, None], idx]
    # define backward fill
    def _bfill(a):
        mask = np.isnan(a)
        idx = np.where(~mask, np.arange(a.shape[1]), a.shape[1]-1)
        idx = np.minimum.accumulate(idx[:, ::-1], axis=1)[:, ::-1]
        return a[np.arange(a.shape[0])[:, None], idx]
    # apply methods
    if m == 'both':
        temp = _ffill(arr2)
        filled = _bfill(temp)
    elif m == 'ff':
        filled = _ffill(arr2)
    else:
        filled = _bfill(arr2)
    # return to original shape
    return filled[0] if filled.shape[0] == 1 else filled

def fill_nan(
    arr: Union[Sequence[Any], np.ndarray],
    *,
    method: str = 'ff',
    axis: int = 1
) -> np.ndarray:
    """
    Efficiently forward/backward fill NaNs along a given axis.

    Parameters
    ----------
    arr : sequence or ndarray
        1D or 2D array containing NaNs.
    method : {'ff','bf','both'}, default 'ff'
        - 'ff'   : forward-fill (use last valid value).
        - 'bf'   : backward-fill (use next valid value).
        - 'both' : apply forward then backward fill.
    axis : {0,1}, default 1
        Axis along which to fill (1 → each row; 0 → each column).

    Returns
    -------
    filled : ndarray
        New array with NaNs replaced. Shape is identical to `arr`.

    Notes
    -----
    - For 1D inputs, `axis` is ignored (returns 1D).
    - `ff` fills NaNs _after_ a valid value;
      `bf` fills NaNs _before_ a valid value.
    - If you need multi-dimensional arrays, flatten or apply
      per-slice yourself; this only supports up to 2D.

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.utils.arrayops import fill_nan
    >>> a = np.array([np.nan, 1, np.nan, 2])
    >>> fill_nan(a, method='ff')
    array([nan, 1., 1., 2.])
    >>> fill_nan(a, method='bf')
    array([1., 1., 2., 2.])
    >>> M = np.array([[1, np.nan, 3],
    ...               [np.nan, 5, np.nan]])
    >>> fill_nan(M, method='both', axis=0)
    array([[1., 5., 3.],
           [1., 5., 3.]])
    """


    # Convert to float array
    a = np.array(arr, dtype=float)

    # Ensure 1D or 2D
    if a.ndim == 1:
        orig1d = True
        a = a.reshape(1, -1)
    elif a.ndim == 2:
        orig1d = False
    else:
        raise ValueError(
            f"fill_nan supports 1D/2D arrays; got ndim={a.ndim}"
        )

    # Validate axis
    if axis not in (0, 1):
        raise ValueError(f"axis must be 0 or 1, got {axis!r}")

    # Normalize method name
    m = method.lower().strip()
    if m in ('forward', 'ff', 'fwd'):
        m = 'ff'
    elif m in ('backward', 'bf', 'bwd'):
        m = 'bf'
    elif m in ('both', 'fb', 'bff', 'ffbf'):
        m = 'both'
    else:
        raise ValueError(
            "method must be 'ff', 'bf', or 'both'; got "
            f"{method!r}"
        )

    def _ffill_row(row: np.ndarray) -> np.ndarray:
        """Forward-fill a single 1D array of length N."""
        mask = np.isnan(row)
        # idx = index of last valid value (or 0 before any)
        idx = np.where(~mask, np.arange(row.size), 0)
        np.maximum.accumulate(idx, out=idx)
        return row[idx]

    def _bfill_row(row: np.ndarray) -> np.ndarray:
        """Backward-fill a single 1D array of length N."""
        mask = np.isnan(row)
        # idx = index of next valid value (or last index after any)
        idx = np.where(~mask, np.arange(row.size), row.size - 1)
        # propagate backwards via reversed accumulate
        rev = np.minimum.accumulate(idx[::-1])[::-1]
        return row[rev]

    # Select and apply fill along axis
    if axis == 1:
        filler = lambda mat, fn: np.vstack(fn(row) for row in mat)
    else:
        # for axis=0, transpose, fill, then transpose back
        filler = lambda mat, fn: (
            np.vstack(fn(col) for col in mat.T).T
        )

    if m == 'ff':
        out = filler(a, _ffill_row)
    elif m == 'bf':
        out = filler(a, _bfill_row)
    else:  # both
        tmp = filler(a, _ffill_row)
        out = filler(tmp, _bfill_row)

    # Restore 1D shape if needed
    return out[0] if orig1d else out

def drop_nan_in(
    y_true,
    *y_preds: Any,
    error: str = "raise",
    nan_policy: str | None = None,
):
    r"""
    Drop NaNs from ``y_true`` and align all ``y_preds``.

    Samples where ``y_true`` is NaN are removed from ``y_true``
    and from each array in ``y_preds`` so their lengths match.

    Parameters
    ----------
    y_true : array-like, shape (n,)
        True target values. Must be 1-D and same length as
        each array in ``y_preds``.
    *y_preds : array-like, shape (n,)
        One or more prediction arrays aligned to ``y_true``.
    error : {"raise","warn","ignore"}, default "raise"
        Action when NaNs are seen in ``y_true`` **and**
        ``nan_policy`` is ``None``.
    nan_policy : {"raise","propagate","omit"}, optional
        NaN behavior override:
        - ``"raise"``: error if NaNs exist.
        - ``"propagate"``: return inputs unchanged.
        - ``"omit"``: drop NaN rows (preferred).

    Returns
    -------
    y_true_f : ndarray, shape (m,)
        Filtered true values (``m <= n``).
    *y_preds_f : ndarrays, each shape (m,)
        Filtered predictions in the same order as input.

    Notes
    -----
    If ``nan_policy`` is given, it overrides ``error``.
    With ``propagate``, no filtering is done.

    Examples
    --------
    >>> yt = np.array([1., 2., np.nan, 4.])
    >>> yp = np.array([0.9, 1.8, 3.1, 4.2])
    >>> drop_nan_in(yt, yp, error="warn")
    (array([1., 2., 4.]), array([0.9, 1.8, 4.2]))
    """
    # to arrays (float for NaN detection)
    yt = np.asarray(y_true, dtype=float)
    preds = [np.asarray(p) for p in y_preds]

    # shape checks
    for i, p in enumerate(preds):
        if p.shape != yt.shape:
            raise ValueError(
                f"y_pred #{i} shape {p.shape} != y_true {yt.shape}"
            )

    has_nan = np.isnan(yt)
    any_nan = bool(np.any(has_nan))

    # policy overrides error
    if nan_policy is not None:
        if nan_policy == "raise" and any_nan:
            raise ValueError("NaNs in y_true (nan_policy='raise').")
        if nan_policy == "propagate":
            return yt, *preds
        if nan_policy != "omit" and nan_policy != "raise":
            raise ValueError(
                "nan_policy must be one of "
                "{'raise','propagate','omit'}."
            )
        # fall through to omit

    else:
        # handle error when policy not given
        if any_nan:
            if error == "raise":
                raise ValueError("NaNs in y_true.")
            if error == "warn":
                warnings.warn(
                    "NaNs in y_true; dropping rows.", stacklevel=2
                )
            elif error != "ignore":
                raise ValueError(
                    "error must be one of "
                    "{'raise','warn','ignore'}."
                )

    # build mask (omit or ignore/warn path)
    mask = ~has_nan
    yt_f = yt[mask]
    preds_f = tuple(p[mask] for p in preds)
    return yt_f, *preds_f

