# -*- coding: utf-8 -*-
#   License: BSD-3-Clause
#   Author: LKouadio <etanoyau@gmail.com>
"""
Additional plot utilities. 
"""
from __future__ import annotations 

import copy 
import warnings
import numpy as np
import pandas as pd 

import matplotlib.pyplot as plt

from sklearn.preprocessing import MinMaxScaler 
from sklearn.impute import SimpleImputer 
from sklearn.utils.validation import check_array, check_y 

from ..utils.arrayops import ( 
    is_iterable, 
    frameify, 
    is_depth_in, # revise 
 )

from .validation import _assert_all_types, isin_if 
from ..utils.text import str2columns 

from ..utils.plot import make_plot_colors

def plot_logging ( 
    X, 
    y=None, 
    zname = None, 
    tname = None, 
    labels=None,
    impute_nan=True , 
    normalize = False, 
    log10=False, 
    columns_to_skip =None, 
    pattern = None, 
    strategy='mean',  
    posiy= None, 
    fill_value = None,  
    fig_size = (16, 7),
    fig_dpi = 300, 
    colors = None,
    cs4_colors=False, 
    sns_style =False, 
    savefig = None,
    draw_spines=False, 
    seed=None, 
    verbose=0, 
    **kws
    ): 
    """ Plot logging data  
    
    Plot expects a collection of logging data. Each logging data composes a 
    column of data collected on the field.Note that can also plot anykind of 
    data related that it contains numerical values. The function does not 
    accept categorical data.   If categorical data are given, they should be 
    discarded. 
    
    Parameters 
    -----------
    X : Dataframe of shape (n_samples, n_features)
         where `n_samples` is the number of data, expected to be the data 
         collected at different depths and `n_features` is the number of 
         columns (features) that supposed to be plot. 
         Note that `X` must include the ``depth`` columns. If not given a 
         relative depth should be created according to the number of sample 
         that composes `X`.
 
    y : array-like or series of shape (n_samples,), optional
        Target relative to X for classification or regression; If given, by 
        default the target plot should be located at the last position. 
        However with the argument of `posiy` , target plot can be toggled to  
        the desired position. 

    zname: str, default='depth' or 'None'
        The name of the depth column in `X`. If the name 'depth' is not  
        specified as the main depth columns, an other name in the columns 
        that matches the depth can also be indicated so the function will put 
        aside this columm as depth column for plot purpose. If set to ``None``, 
        `zname` holds the name ``depth`` and assumes that depth exists in 
        `X` columns.
    tname: str, optional, 
        name of the target. This can rename of the target name if given `y`
        as a pandas series  or add the name of target if given as an array-like. 
        If not provided, it should use the name of the target series if `y` is
        not None. 
        
    normalize: bool, default = False
        Normalize all the data to be range between (0, 1) except the `depth`,    
        
    labels: list or str, optional
        If labels are given, they should fit the size of the number of 
        columns. The given labels should replace the old columns in `X` and 
        should figue out in the plot. This is usefull to change the columns 
        labels in the dataframe to a new labels that describe the best the 
        plot ; for instance by inluding the units in the new labels. Note that 
        if the labels do not match the size of the old columns in `X` a warning 
        should be let to the user and none operation will be performed. 
        
    impute_nan: bool, default=True, 
        Replace the NaN values in the dataframe. Note that the default 
        behaviour for replacing NaN is the ``mean``. However if the argument 
        of `fill_value` is provided,the latter should be used to replace 'NaN' 
        in `X`. 
        
    log10: bool, default=False
        Convert values to log10. This can be usefull when using the logarithm 
        data. However, it seems not all the data can be used this operation, 
        for instance, a negative data. In that case, `column_to_skip` argument
        is usefull to provide so to skip that columns when converting values 
        to log10. 
 
    columns_to_skip: list or str, optional, 
        Columns to skip when performing some operation like 'log10'. These 
        columns with not be affected by the 'log10' operations. Note that 
       `columns_to_skip` can also gives as litteral string. In that case, the 
       `pattern` is need to parse the columns into a list of string. 
       
    pattern: str, default = '[#&*@!,;\s]\s*'
        Regex pattern to parse the `columns_to_skip` into a list of string 
        where each item is a column name especially when the latter is given 
        as litteral text string. For instance:: 
            
            columns_to_skip='depth_top, thickness, sp, gamma_gamma'  
            -> ['depth_top', 'thickness', 'sp', 'gamma_gamma']
            
        by using the default pattern. To have full control of columns splitted
        it is recommended to provided your own pattern to avoid wrong parsing 
        and can lead to an error. 
        
    strategy : str, default='mean'
        The imputation strategy.

        - If "mean", then replace missing values using the mean along
          each column. Can only be used with numeric data.
        - If "median", then replace missing values using the median along
          each column. Can only be used with numeric data.
        - If "most_frequent", then replace missing using the most frequent
          value along each column. Can be used with strings or numeric data.
          If there is more than one such value, only the smallest is returned.
        - If "constant", then replace missing values with fill_value. Can be
          used with strings or numeric data.

    fill_value : str or numerical value, optional
        When strategy == "constant", fill_value is used to replace all
        occurrences of missing_values.
        If left to the default, fill_value will be 0 when imputing numerical
        data and "missing_value" for strings or object data types. If not 
        given and `impute_nan` is ``True``, the mean strategy is used instead.

    posiy: int, optional 
        the position to place the target plot `y` . By default the target plot 
        if given is located at the last position behind the logging plots. 
    
    colors: str, list of Matplotlib.colors map, optional 
        The colors for plotting each columns of `X` except the depth. If not
        given, default colors are auto-generated.
        
        If `colors` is string and 'cs4'or 'xkcd' is included. 
        Matplotlib.colors.CS4_COLORS or Matplotlib.colors.XKCD_COLORS 
        should be used instead. In addition if the `'cs4'` or `'xkcd'` is  
        suffixed by colons and integer value like ``cs4:4`` or ``xkcd:4``, the 
        CS4 or XKCD colors should be used from index equals to ``4``. 
        
        .. versionadded:: 0.2.3 
           Matplotlib.colors.CS4_COLORS or Matplotlib.colors.XKCD_COLORS can 
           be used by setting `colors` to ``'cs4'`` or ``'xkcd'``. To reproduce 
           the same CS4 or XKCD colors, set the `seed` parameter to a 
           specific value. 
        
    draw_spines: bool, tuple (-lim, +lim), default= False, 
        Only draw spine between the y-ticks. ``-lim`` and ``+lim`` are lower 
        and upper bound i.e. a range to draw the spines in y-axis. 
        
    fig_size : tuple (width, height), default =(8, 6)
        the matplotlib figure size given as a tuple of width and height
        
    fig_dpi: float or 'figure', default: rcParams["savefig.dpi"] \
        (default: 'figure')
        The resolution in dots per inch. If 'figure', use the figure's dpi value.
        
    savefig: str, default =None , 
        the path to save the figure. Argument is passed to 
        :class:`matplotlib.Figure` class. 

    sns_style: str, optional, 
        the seaborn style.
        
    seed: int, optional 
       Allow to reproduce the Matplotlib.colors.CS4_COLORS if `colors` is 
       set to ``cs4``. 
       
       .. versionadded:: 0.2.3 
       
    verbose: int, default=0 
        Output the number of categorial features dropped in the dataframe.  
        
    kws: dict, 
        Additional keyword arguments passed to :func:`matplotlib.axes.plot`
        
    Examples
    ---------
    >>> from watex.datasets import load_hlogs 
    >>> from watex.utils.plotutils import plot_logging
    >>> X0, y = load_hlogs (as_frame =True) # get the frames rather than object 
    >>> # plot the default logging with Normalize =True 
    >>> plot_logging (X0, normalize =True) 
    >>> # Include the target in the plot 
    >>> plot_logging ( X0,  y = y.kp , posiy = 0, 
                      columns_to_skip=['thickness', 'sp'], 
                      log10 =True, 
                      )
    >>> # draw spines and limit plot from (0, 700) m depth 
    >>> plot_logging (X0 , y= y.kp, draw_spines =(0, 700) )
    """
    
    X = _assert_all_types(X, pd.DataFrame, pd.Series , np.ndarray ) 
    X= check_array (
        X, 
        dtype =object, 
        force_all_finite="allow-nan", 
        input_name ="Logging dataset",
        to_frame =True  
        )
    # Discard all categorical values and 
    # keep only the numerical features.
    # drop the complete Nan columns and rows
    X = frameify(X, pop_cat_features=True, verbose = verbose ) 

    if y is not None: 
       if isinstance (y, (list, tuple)): 
           # in the case a lst is given 
           y = np.array (y) 
       if not is_iterable (y): 
           raise TypeError ("y expects an iterable object."
                              f" got {type(y).__name__!r}")
       y = _assert_all_types(y, pd.Series, pd.DataFrame, np.ndarray)
       
       y=check_y (
            y, 
            to_frame =True, 
            allow_nan= True,
            )
       
       if len(y) !=len(X): 
           raise ValueError ("y and X sizes along axis 0 must be consistent;"
                             f" {len(y)} and {len(X)} are given.")
    # return X and depth 
    X, depth = is_depth_in(X, zname or 'depth', columns = labels 
                            )
    # fetch target if is given  
    X, y   = _is_target_in(X, y = y , tname = tname )
    
    # skip log10 columns if log 10 is set to True 
    if log10: 
        X = _skip_log10_columns (X, column2skip = columns_to_skip , 
                                  pattern= pattern, inplace =False) 
    # if normalize then  
    if normalize:
        msc = MinMaxScaler()
        Xsc = msc.fit_transform (X)
        # set a new dataframe with features
        if hasattr (msc , 'feature_names_in_'): 
            X = pd.DataFrame (Xsc , columns = list(msc.feature_names_in_ )
                                )
        else : X = pd.DataFrame(Xsc, columns =list(X.columns )) 
        # set the x axis and delete the normalize from X 
        # at index 0 supposed to be the x axis 
        # Xsc.iloc [:, 0 ] = x_ser 
        # X= Xsc.copy()  
    # impute_nan 
    if impute_nan: 
        # check whether there is a Nan value  in the data 
        # impute data using mean values
        if X.isnull().values.any(): 
            Xi= SimpleImputer(strategy= strategy if not fill_value else None, 
                             fill_value= fill_value
                             ).fit_transform(X)
            X = pd.DataFrame(Xi, columns= X.columns)
            
    # toggle y 
    if y is not None: 
        X = _toggle_target_in(X, y, pos = posiy)
        
    #manage colors along colors 
    colors = make_plot_colors (
        X, colors = colors , axis = 1, seed = seed , chunk=False )

    fig, ax = plt.subplots (1, ncols = X.shape [1], sharey = True , 
                            figsize = fig_size )
    
    # customize bound and set spines 
    for k in range (X.shape [1]): 
     
        ax[k].plot ( X.iloc[:, k], 
                    depth, 
                    color = colors[k], 
                    **kws
                    )
        ax[k].tick_params(top=True, 
                          labeltop=True, 
                          bottom=False, 
                          labelbottom=False
                       )
        ax[k].set_title (X.columns [k])
        ax[k].spines['right'].set_visible(False)
        ax[k].spines['bottom'].set_visible(False)
        # only show tick on the top and left 
        ax[k].xaxis.set_ticks_position('top')
        if y is not None: 
            # make X axis of the target to red 
            # for differenciation from features. 
            if X.columns [k] ==y.name: 
                ax[k].spines['top'].set_color('red')  
         
        if draw_spines: 
            # Only draw spine between the y-ticks
            if is_iterable(draw_spines): 
                # for consistency check whether values 
                # are numeric
                draw_spines = sorted (
                    list(map (lambda x: float (x) , draw_spines[:2])) 
                    ) 
                if len(draw_spines) <2: 
                    warnings.warn(
                        "Spine bounds is a tuple of (startpoint, endpoint)"
                         " Single limit value is not allowed."
                         )
            else: 
                # in case only True is given 
                # use the default plot
                ytv= ax[0].get_yticks () 
                spacing = (ytv[-1] - ytv[0] )/(len(ytv)-1) 
                # commonly matplotlib axis extrapoled the limit so 
                # start with the first and last index 
                draw_spines=  (ytv[0] + spacing/2 , ytv[-1] - spacing/2 ) 
                
            ax[k].spines['left'].set_bounds(*draw_spines )
     
    # set labels
    ax[0].set_ylabel ("Depth (m)")
        # Tweak spacing between subplots to prevent labels 
        # from overlapping 
        # plt.subplots_adjust(hspace=0.5)-> removed
    plt.gca().invert_yaxis()
    
    if savefig is not None:
        plt.savefig(savefig, dpi = fig_dpi )
        
    plt.close () if savefig is not None else plt.show() 
    

def _is_target_in (X, y=None, tname=None): 
    """ Create new target name for tname if given 
    
    :param X: dataframe 
        dataframe containing the data for plotting 
    :param y: array or series
        target data for plotting. Note that multitarget outpout is not 
        allowed yet. Moroever, it `y` is given as a dataframe, 'tname' must 
        be supplied to retrive y as a pandas series object, otherwise an 
        error will raise. 
    :param tname: str,  
        target name. If given and `y` is ``None``, Will try to find `tname`
        in the `X` columns. If 'tname' does not exist, plot for target is 
        cancelled. 
        
    :return y: Series 
    """
    _assert_all_types(X, pd.DataFrame)
    
    if y is not None: 
        y = _assert_all_types(y , pd.Series, pd.DataFrame, np.ndarray)
        
        if hasattr (y, 'columns'): 
            if tname not in (y.columns): tname = None 
            if tname is None: 
                raise TypeError (
                    "'tname' must be supplied when y is a dataframe.")
            y = y [tname ]
        elif hasattr (y, 'name'): 
            tname = tname or y.name 
            # reformat inplace the name of series 
            y.name = tname 
            
        elif hasattr(y, '__array__'): 
            y = pd.Series (y, name = tname or 'target')
            
    elif y is None: 
        if tname in X.columns :
            y = X.pop(tname)

    return X, y 

def _toggle_target_in  (X , y , pos=None): 
    """ Toggle the target in the convenient position. By default the target 
    plot is the last subplots 
    
    :param X: dataframe 
        dataframe containing the data for plotting 
    :param y: array or series
        the target for  plotting. 
    :param pos: int, the position to insert y in the dataframe X 
        By default , `y` is located at the last position 
        
    :return: Dataframe 
        Dataframe containing the target 'y'
        
    """
    
    pos =  0 if pos ==0  else ( pos or X.shape [1])

    pos= int ( _assert_all_types(pos, int, float ) ) 
    ms= ("The positionning of the target is out of the bound."
         "{} position is used instead.")
    
    if pos > X.shape[1] : 
        warnings.warn(ms.format('The last'))
        pos=X.shape[1]
    elif pos < 0: 
        warnings.warn(ms.format(
            " Negative index is not allowed. The first")
                      )
        pos=0 
 
    X.insert (pos, y.name, y )
    
    return X
    
def _skip_log10_columns ( X, column2skip, pattern =None , inplace =True): 
    """ Skip the columns that dont need to put value in logarithms.
    
    :param X: dataframe 
        pandas dataframe with valid columns 
    :param column2skip: list or str , 
        List of columns to skip. If given as string and separed by the default
        pattern items, it should be converted to a list and make sure the 
        columns name exist in the dataframe. Otherwise an error with 
        raise. 
    :param pattern: str, default = '[#&*@!,;\s]\s*'
        The base pattern to split the text in `column2skip` into a columns
        
    :return X: Dataframe
        Dataframe modified inplace with values computed in log10 
        except the skipped columns. 
        
    :example: 
       >>> from watex.datasets import load_hlogs 
       >>> from watex.utils.plotutils import _skip_log10_columns 
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
       >>> _skip_log10_columns (X0, column2skip)
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
      >>> # while `sp` still remains the same 
      
    """
    X0 = X.copy () 
    if not is_iterable( column2skip): 
        raise TypeError ("Columns  to skip expect an iterable object;"
                         f" got {type(column2skip).__name__!r}")
        
    pattern = pattern or r'[#&*@!,;\s]\s*'
    
    if isinstance(column2skip, str):
        column2skip = str2columns (column2skip, pattern=pattern  )
    #assert whether column to skip is in 
    if column2skip:
        cskip = copy.deepcopy (column2skip) 
        column2skip = isin_if(X.columns, column2skip, return_diff= True)
        if len(column2skip) ==len (X.columns): 
            warnings.warn("Value(s) to skip are not detected.")
        if inplace : 
            X[column2skip] = np.log10 ( X[column2skip] ) 
            X.drop (columns =cskip , inplace =True )
            return  
        else : 
            X0[column2skip] = np.log10 ( X0[column2skip] ) 
            
    return X0