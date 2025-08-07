# -*- coding: utf-8 -*-
#   License: BSD-3-Clause
#   Author: LKouadio <etanoyau@gmail.com>
"""
Additional plot utilities. 
"""
from __future__ import annotations 
import os
import re 
import copy 
import datetime 
import warnings
import itertools 
import numpy as np
import pandas as pd 
import matplotlib as mpl 
from matplotlib.patches import Ellipse
import matplotlib.colors as mcolors
import matplotlib.transforms as transforms 

import seaborn as sns 
import matplotlib.pyplot as plt

from sklearn.preprocessing import MinMaxScaler 
from sklearn.impute import SimpleImputer 
from sklearn.utils.validation import check_array, check_y 

from ..exceptions import  PlotError
from .arrayops import ( 
    assert_xy_in, 
    is_iterable, 
    frameify, 

    is_in_if, # revise 
    is_depth_in, # revise 
    )
from .text import str2columns 

from .validation import ( 
    _assert_all_types, 
)

D_COLORS =[
    'g',
    'gray',
    'y', 
    'blue',
    'orange',
    'purple',
    'lime',
    'k', 
    'cyan', 
    (.6, .6, .6),
    (0, .6, .3), 
    (.9, 0, .8),
    (.8, .2, .8),
    (.0, .9, .4)
]

D_MARKERS =[
    'o',
    '^',
    'x',
    'D',
    '8',
    '*',
    'h',
    'p',
    '>',
    'o',
    'd',
    'H'
]

D_STYLES = [
    '-',
    '-',
    '--',
    '-.',
    ':', 
    'None',
    ' ',
    '',
    'solid', 
    'dashed',
    'dashdot',
    'dotted' 
]

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
    
def make_plot_colors(d , / , colors:str | list[str]=None , axis:int = 0, 
                     seed:int  =None, chunk:bool =... ): 
    """ Select colors according to the data size along axis 
    
    Parameters 
    ----------
    d: Arraylike 
       Array data to select colors according to the axis 
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
           
    axis: int, default=0 
       Axis along with the colors must be generated. By default colors is 
       generated along the row axis 
       
    seed: int, optional 
       Allow to reproduce the Matplotlib.colors.CS4_COLORS if `colors` is 
       set to ``cs4``. 
       
    chunk: bool, default=True 
       Chunk generated colors to fit the exact length of the `d` size 
       
    Returns 
    -------
    colors: list 
       List of new generated colors 
       
    Examples 
    --------
    >>> import numpy as np 
    >>> from watex.utils.plotutils import make_plot_colors
    >>> ar = np.random.randn (7, 2) 
    >>> make_plot_colors (ar )
    ['g', 'gray', 'y', 'blue', 'orange', 'purple', 'lime']
    >>> make_plot_colors (ar , axis =1 ) 
    Out[6]: ['g', 'gray']
    >>> make_plot_colors (ar , axis =1 , colors ='cs4')
    ['#F0F8FF', '#FAEBD7']
    >>> len(make_plot_colors (ar , axis =1 , colors ='cs4', chunk=False))
    150
    >>> make_plot_colors (ar , axis =1 , colors ='cs4:4')
    ['#F0FFFF', '#F5F5DC']
    """
    
    # get the data size where colors must be fitted. 
    # note colors should match either the row axis or colurms axis 
    axis = str(axis).lower() 
    if 'columns1'.find (axis)>=0: 
        axis =1 
    else: axis =0
    
    # manage the array 
    d= is_iterable( d, exclude_string=True, transform=True)
    if not hasattr (d, '__array__'): 
        d = np.array(d, dtype =object ) 
    
    axis_length = len(d) if len(d.shape )==1 else d.shape [axis]
    m_cs = make_mpl_properties(axis_length )
    
     #manage colors 
    # we assume the first columns is dedicated for 
    if colors ==...: colors =None 
    if ( 
            isinstance (colors, str) and 
            ( 
                "cs4" in str(colors).lower() 
                 or 'xkcd' in str(colors).lower() 
                 )
            ): 
        #initilize colors infos
        c = copy.deepcopy(colors)
        if 'cs4' in str(colors).lower() : 
            DCOLORS = mcolors.CSS4_COLORS
        else: 
            # remake the dcolors my removing the xkcd: in the keys: 
            DCOLORS = dict(( (k.replace ('xkcd:', ''), c) 
                            for k, c in mcolors.XKCD_COLORS.items()))  
        
        key_colors = list(DCOLORS.keys ())
        colors = list(DCOLORS.values() )
        
        shuffle_cs4=True 
        
        cs4_start= None
        #------
        if ':' in str(c).lower():
            cs4_start = str(c).lower().split(':')[-1]
        #try to converert into integer 
        try: 
            cs4_start= int (cs4_start)
        except : 
            if str(cs4_start).lower() in key_colors: 
                cs4_start= key_colors.index (cs4_start)
                shuffle_cs4=False
            else: 
                pass 
        
        else: shuffle_cs4=False # keep CS4 and dont shuffle 
        
        cs4_start= cs4_start or 0
        
        if shuffle_cs4: 
            np.random.seed (seed )
            colors = list(np.random.choice(colors  , len(m_cs)))
        else: 
            if cs4_start > len(colors)-1: 
                cs4_start = 0 
    
            colors = colors[ cs4_start:]
    
    if colors is not None: 
        if not is_iterable(colors): 
            colors =[colors]
        colors += m_cs 
    else :
        colors = m_cs 
        
    # shrunk data to map the exact colors 
    chunk =True if chunk is ... else False 
    return colors[:axis_length] if chunk else colors 


def savefigure (fig: object ,
             figname: str = None,
             ext:str ='.png',
             **skws ): 
    """ save figure from the given figure name  
    
    :param fig: Matplotlib figure object 
    :param figname: name of figure to output 
    :param ext: str - extension of the figure 
    :param skws: Matplotlib savefigure keywards additional keywords arguments 
    
    :return: Matplotlib savefigure objects. 
    
    """
    ext = '.' + str(ext).lower().strip().replace('.', '')

    if figname is None: 
        figname =  '_' + os.path.splitext(os.path.basename(__file__)) +\
            datetime.datetime.now().strftime('%m-%d-%Y %H:%M:%S') + ext
        warnings.warn("No name of figure is given. Figure should be renamed as "
                      f"{figname!r}")
        
    file, ex = os.path.splitext(figname)
    if ex in ('', None): 
        ex = ext 
        figname = os.path.join(file, f'{ext}')

    return  fig.savefig(figname, **skws)


def resetting_ticks ( get_xyticks,  number_of_ticks=None ): 
    """
    resetting xyticks  modulo , 100
    
    :param get_xyticks:  xyticks list  , use to ax.get_x|yticks()
    :type  get_xyticks: list 
    
    :param number_of_ticks:  maybe the number of ticks on x or y axis 
    :type number_of_ticks: int
    
    :returns: a new_list or ndarray 
    :rtype: list or array_like 
    """
    if not isinstance(get_xyticks, (list, np.ndarray) ): 
        warnings.warn (
            'Arguments get_xyticks must be a list'
            ' not <{0}>.'.format(type(get_xyticks)))
        raise PlotError (
            '<{0}> found. "get_xyticks" must be a '
            'list or (nd.array,1).'.format(type(get_xyticks)))
    
    if number_of_ticks is None :
        if len(get_xyticks) > 2 : 
            number_of_ticks = int((len(get_xyticks)-1)/2)
        else : number_of_ticks  = len(get_xyticks)
    
    if not(number_of_ticks, (float, int)): 
        try : number_of_ticks=int(number_of_ticks) 
        except : 
            warnings.warn('"Number_of_ticks" arguments is the times to see '
                          'the ticks on x|y axis.'\
                          ' Must be integer not <{0}>.'.
                          format(type(number_of_ticks)))
            raise PlotError(f'<{type(number_of_ticks).__name__}> detected.'
                            ' Must be integer.')
        
    number_of_ticks=int(number_of_ticks)
    
    if len(get_xyticks) > 2 :
        if get_xyticks[1] %10 != 0 : 
            get_xyticks[1] =get_xyticks[1] + (10 - get_xyticks[1] %10)
        if get_xyticks[-2]%10  !=0 : 
            get_xyticks[-2] =get_xyticks[-2] -get_xyticks[-2] %10
    
        new_array = np.linspace(get_xyticks[1], get_xyticks[-2],
                                number_of_ticks )
    elif len(get_xyticks)< 2 : 
        new_array = np.array(get_xyticks)
 
    return  new_array
        
def make_mpl_properties(n ,prop ='color'): 
    """ make matplotlib property ('colors', 'marker', 'line') to fit the 
    numer of samples
    
    :param n: int, 
        Number of property that is needed to create. It generates a group of 
        property items. 
    :param prop: str, default='color', name of property to retrieve. Accepts 
        only 'colors', 'marker' or 'line'.
    :return: list of property items with size equals to `n`.
    :Example: 
        >>> from watex.utils.plotutils import make_mpl_properties
        >>> make_mpl_properties (10 )
        ... ['g',
             'gray',
             'y',
             'blue',
             'orange',
             'purple',
             'lime',
             'k',
             'cyan',
             (0.6, 0.6, 0.6)]
        >>> make_mpl_properties(100 , prop = 'marker')
        ... ['o',
             '^',
             'x',
             'D',
              .
              .
              .
             11,
             'None',
             None,
             ' ',
             '']
        >>> make_mpl_properties(50 , prop = 'line')
        ... ['-',
             '-',
             '--',
             '-.',
               .
               .
               . 
             'solid',
             'dashed',
             'dashdot',
             'dotted']
        
    """ 
    n=int(_assert_all_types(n, int, float, objname ="'n'"))
    prop = str(prop).lower().strip().replace ('s', '') 
    if prop not in ('color', 'marker', 'line'): 
        raise ValueError ("Property {prop!r} is not availabe yet. , Expect"
                          " 'colors', 'marker' or 'line'.")
    # customize plots with colors lines and styles 
    # and create figure obj
    if prop=='color': 
        d_colors =  D_COLORS 
        d_colors = mpl.colors.ListedColormap(d_colors[:n]).colors
        if len(d_colors) == n: 
            props= d_colors 
        else:
            rcolors = list(itertools.repeat(
                d_colors , (n + len(d_colors))//len(d_colors))) 
    
            props  = list(itertools.chain(*rcolors))
        
    if prop=='marker': 
        
        d_markers =  D_MARKERS + list(mpl.lines.Line2D.markers.keys()) 
        rmarkers = list(itertools.repeat(
            d_markers , (n + len(d_markers))//len(d_markers))) 
        
        props  = list(itertools.chain(*rmarkers))
    # repeat the lines to meet the number of cv_size 
    if prop=='line': 
        d_lines =  D_STYLES
        rlines = list(itertools.repeat(
            d_lines , (n + len(d_lines))//len(d_lines))) 
        # combine all repeatlines 
        props  = list(itertools.chain(*rlines))
    
    return props [: n ]
       
def resetting_colorbar_bound(cbmax ,
                             cbmin,
                             number_of_ticks = 5, 
                             logscale=False): 
    """
    Function to reset colorbar ticks more easy to read 
    
    :param cbmax: value maximum of colorbar 
    :type cbmax: float 
    
    :param cbmin: minimum data value 
    :type cbmin: float  minimum data value
    
    :param number_of_ticks:  number of ticks should be 
                            located on the color bar . Default is 5.
    :type number_of_ticks: int 
    
    :param logscale: set to True if your data are lograith data . 
    :type logscale: bool 
    
    :returns: array of color bar ticks value.
    :rtype: array_like 
    """
    def round_modulo10(value): 
        """
        round to modulo 10 or logarithm scale  , 
        """
        if value %mod10  == 0 : return value 
        if value %mod10  !=0 : 
            if value %(mod10 /2) ==0 : return value 
            else : return (value - value %mod10 )
    
    if not(number_of_ticks, (float, int)): 
        try : number_of_ticks=int(number_of_ticks) 
        except : 
            warnings.warn('"Number_of_ticks" arguments '
                          'is the times to see the ticks on x|y axis.'
                          ' Must be integer not <{0}>.'.format(
                              type(number_of_ticks)))
            raise PlotError('<{0}> detected. Must be integer.')
        
    number_of_ticks=int(number_of_ticks)
    
    if logscale is True :  mod10 =np.log10(10)
    else :mod10 = 10 
       
    if cbmax % cbmin == 0 : 
        return np.linspace(cbmin, cbmax , number_of_ticks)
    elif cbmax% cbmin != 0 :
        startpoint = cbmin + (mod10  - cbmin % mod10 )
        endpoint = cbmax - cbmax % mod10  
        return np.array(
            [round_modulo10(ii) for ii in np.linspace(
                             startpoint,endpoint, number_of_ticks)]
            )
    

            
def controle_delineate_curve(res_deline =None , phase_deline =None ): 
    """
    fonction to controle delineate value given  and return value ceilling .
    
    :param  res_deline:  resistivity  value todelineate. unit of Res in `ohm.m`
    :type  res_deline: float|int|list  
    
    :param  phase_deline:   phase value to  delineate , unit of phase in degree
    :type phase_deline: float|int|list  
    
    :returns: delineate resistivity or phase values 
    :rtype: array_like 
    """
    fmt=['resistivity, phase']
 
    for ii, xx_deline in enumerate([res_deline , phase_deline]): 
        if xx_deline is  not None  : 
            if isinstance(xx_deline, (float, int, str)):
                try :xx_deline= float(xx_deline)
                except : raise PlotError(
                        'Value <{0}> to delineate <{1}> is unacceptable.'\
                         ' Please ckeck your value.'.format(xx_deline, fmt[ii]))
                else :
                    if ii ==0 : return [np.ceil(np.log10(xx_deline))]
                    if ii ==1 : return [np.ceil(xx_deline)]
  
            if isinstance(xx_deline , (list, tuple, np.ndarray)):
                xx_deline =list(xx_deline)
                try :
                    if ii == 0 : xx_deline = [
                            np.ceil(np.log10(float(xx))) for xx in xx_deline]
                    elif  ii ==1 : xx_deline = [
                            np.ceil(float(xx)) for xx in xx_deline]
                        
                except : raise PlotError(
                        'Value to delineate <{0}> is unacceptable.'\
                         ' Please ckeck your value.'.format(fmt[ii]))
                else : return xx_deline


def fmt_text (data_text, fmt='~', leftspace = 3, return_to_line =77) : 
    """
    Allow to format report with data text , fm and leftspace 

    :param  data_text: a long text 
    :type  data_text: str  
        
    :param fmt:  type of underline text 
    :type fmt: str

    :param leftspae: How many space do you want before starting wrinting report .
    :type leftspae: int 
    
    :param return_to_line: number of character to return to line
    :type return_to_line: int 
    """

    return_to_line= int(return_to_line)
    begin_text= leftspace *' '
    text= begin_text + fmt*(return_to_line +7) + '\n'+ begin_text

    
    ss=0
    
    for  ii, num in enumerate(data_text) : # loop the text 
        if ii == len(data_text)-1 :          # if find the last character of text 
            #text = text + data_text[ss:] + ' {0}\n'.format(fmt) # take the 
            #remain and add return chariot 
            text = text+ ' {0}\n'.format(fmt) +\
                begin_text +fmt*(return_to_line+7) +'\n' 
      
 
            break 
        if ss == return_to_line :                       
            if data_text[ii+1] !=' ' : 
                text = '{0} {1}- \n {2} '.format(
                    text, fmt, begin_text + fmt ) 
            else : 
                text ='{0} {1} \n {2} '.format(
                    text, fmt, begin_text+fmt ) 
            ss=0
        text += num    # add charatecter  
        ss +=1

    return text 

def plotvec1(u, z, v):
    """
    Plot tips function with  three vectors. 
    
    :param u: vector u - a vector 
    :type u: array like  
    
    :param z: vector z 
    :type z: array_like 
    
    :param v: vector v 
    :type v: array_like 
    
    return: plot 
    
    """
    
    ax = plt.axes()
    ax.arrow(0, 0, *u, head_width=0.05, color='r', head_length=0.1)
    plt.text(*(u + 0.1), 'u')
    
    ax.arrow(0, 0, *v, head_width=0.05, color='b', head_length=0.1)
    plt.text(*(v + 0.1), 'v')
    ax.arrow(0, 0, *z, head_width=0.05, head_length=0.1)
    plt.text(*(z + 0.1), 'z')
    plt.ylim(-2, 2)
    plt.xlim(-2, 2)

def plotvec2(a,b):
    """
    Plot tips function with two vectors
    Just use to get the orthogonality of two vector for other purposes 

    :param a: vector u 
    :type a: array like  - a vector 
    :param b: vector z 
    :type b: array_like 
    
    *  Write your code below and press Shift+Enter to execute
    
    :Example: 
        
        >>> import numpy as np 
        >>> from watex.utils.plotutils import plotvec2
        >>> a=np.array([1,0])
        >>> b=np.array([0,1])
        >>> Plotvec2(a,b)
        >>> print('the product a to b is =', np.dot(a,b))

    """
    ax = plt.axes()
    ax.arrow(0, 0, *a, head_width=0.05, color ='r', head_length=0.1)
    plt.text(*(a + 0.1), 'a')
    ax.arrow(0, 0, *b, head_width=0.05, color ='b', head_length=0.1)
    plt.text(*(b + 0.1), 'b')
    plt.ylim(-2, 2)
    plt.xlim(-2, 2)  

def plot_errorbar0(
        ax,
        x_ar,
        y_ar,
        y_err=None,
        x_err=None,
        color='k',
        marker='x',
        ms=2, 
        ls=':', 
        lw=1, 
        e_capsize=2,
        e_capthick=.5,
        picker=None,
        **kws
 ):
    """
    convinience function to make an error bar instance
    
    Parameters
    ------------
    
    ax: matplotlib.axes 
        instance axes to put error bar plot on

    x_array: np.ndarray(nx)
        array of x values to plot
                  
    y_array: np.ndarray(nx)
        array of y values to plot
                  
    y_error: np.ndarray(nx)
        array of errors in y-direction to plot
    
    x_error: np.ndarray(ns)
        array of error in x-direction to plot
                  
    color: string or (r, g, b)
        color of marker, line and error bar
                
    marker: string
        marker type to plot data as
                 
    ms: float
        size of marker
             
    ls: string
        line style between markers
             
    lw: float
        width of line between markers
    
    e_capsize: float
        size of error bar cap
    
    e_capthick: float
        thickness of error bar cap
    
    picker: float
          radius in points to be able to pick a point. 
        
        
    Returns:
    ---------
    errorbar_object: matplotlib.Axes.errorbar 
           error bar object containing line data, errorbars, etc.
    """
    # this is to make sure error bars 
    #plot in full and not just a dashed line
    eobj = ax.errorbar(
        x_ar,
        y_ar,
        marker=marker,
        ms=ms,
        mfc='None',
        mew=lw,
        mec=color,
        ls=ls,
        xerr=x_err,
        yerr=y_err,
        ecolor=color,
        color=color,
        picker=picker,
        lw=lw,
        elinewidth=lw,
        capsize=e_capsize,
        # capthick=e_capthick
        **kws
         )
    
    return eobj

def plot_errorbar(
    ax,
    x_ar,
    y_ar,
    y_err=None,
    x_err=None,
    color='k',
    marker='x',
    ms=2, 
    ls=':', 
    lw=1, 
    e_capsize=2,
    e_capthick=.5,
    picker=None,
    show_error_bars=True,  
    **kws
 ):
    """
    Convenience function to make an error bar instance.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes instance to put error bar plot on.
    x_ar : np.ndarray
        Array of x values to plot.
    y_ar : np.ndarray
        Array of y values to plot.
    y_err : np.ndarray, optional
        Array of errors in y-direction to plot.
    x_err : np.ndarray, optional
        Array of errors in x-direction to plot.
    color : str or tuple
        Color of marker, line, and error bar.
    marker : str
        Marker type to plot data as.
    ms : float
        Size of marker.
    ls : str
        Line style between markers.
    lw : float
        Width of line between markers.
    e_capsize : float
        Size of error bar cap.
    e_capthick : float
        Thickness of error bar cap.
    picker : float, optional
        Radius in points to be able to pick a point.
    show_error_bars : bool, default True
        If False, skip plotting the error bars.
    **kws : dict
        Additional keyword arguments passed to `ax.errorbar`.

    Returns
    -------
    errorbar_object : matplotlib.container.ErrorbarContainer
        Error bar object containing line data, error bars, etc.
    """
    if show_error_bars:
        yerr = y_err
        xerr = x_err
    else:
        # Skip error bars by setting them to None
        yerr = None
        xerr = None

    eobj = ax.errorbar(
        x_ar, y_ar, marker=marker, ms=ms, mfc='None', mew=lw, mec=color,
        ls=ls, xerr=xerr, yerr=yerr, ecolor=color, color=color,
        picker=picker, lw=lw, elinewidth=lw, capsize=e_capsize,
        **kws
    )

    return eobj

def get_color_palette (RGB_color_palette): 
    """
    Convert RGB color into matplotlib color palette. In the RGB color 
    system two bits of data are used for each color, red, green, and blue. 
    That means that each color runson a scale from 0 to 255. Black  would be
    00,00,00, while white would be 255,255,255. Matplotlib has lots of
    pre-defined colormaps for us . They are all normalized to 255, so they run
    from 0 to 1. So you need only normalize data, then we can manually  select 
    colors from a color map  

    :param RGB_color_palette: str value of RGB value 
    :type RGB_color_palette: str 
        
    :returns: rgba, tuple of (R, G, B)
    :rtype: tuple
     
    :Example: 
        
        >>> from watex.utils.plotutils import get_color_palette 
        >>> get_color_palette (RGB_color_palette ='R128B128')
    """  
    
    def ascertain_cp (cp): 
        if cp >255. : 
            warnings.warn(
                ' !RGB value is range 0 to 255 pixels , '
                'not beyond !. Your input values is = {0}.'.format(cp))
            raise ValueError('Error color RGBA value ! '
                             'RGB value  provided is = {0}.'
                            ' It is larger than 255 pixels.'.format(cp))
        return cp
    if isinstance(RGB_color_palette,(float, int, str)): 
        try : 
            float(RGB_color_palette)
        except : 
              RGB_color_palette= RGB_color_palette.lower()
             
        else : return ascertain_cp(float(RGB_color_palette))/255.
    
    rgba = np.zeros((3,))
    
    if 'r' in RGB_color_palette : 
        knae = RGB_color_palette .replace('r', '').replace(
            'g', '/').replace('b', '/').split('/')
        try :
            _knae = ascertain_cp(float(knae[0]))
        except : 
            rgba[0]=1.
        else : rgba [0] = _knae /255.
        
    if 'g' in RGB_color_palette : 
        knae = RGB_color_palette .replace('g', '/').replace(
            'b', '/').split('/')
        try : 
            _knae =ascertain_cp(float(knae[1]))
        except : 
            rgba [1]=1.
            
        else :rgba[1]= _knae /255.
    if 'b' in RGB_color_palette : 
        knae = knae = RGB_color_palette .replace('g', '/').split('/')
        try : 
            _knae =ascertain_cp(float(knae[1]))
        except :
            rgba[2]=1.
        else :rgba[2]= _knae /255.
        
    return tuple(rgba)       

def _get_xticks_formatage (
        ax,  xtick_range, space= 14 , step=7,
        fmt ='{}',auto = False, ticks ='x', **xlkws):
    """ Skip xticks label at every number of spaces 
    :param ax: matplotlib axes 
    :param xtick_range: list of the xticks values 
    :param space: interval that the label must be shown.
    :param step: the number of label to skip.
    :param fmt: str, formatage type. 
    :param ticks: str, default='x', the ticks axis to format the labels. 
      can be ``'y'``. 
    :param auto: bool , if ``True`` a dynamic tick formatage will start. 
    
    """
    def format_ticks (ind, x):
        """ Format thick parameter with 'FuncFormatter(func)'
        rather than using:: 
            
        axi.xaxis.set_major_locator (plt.MaxNLocator(3))
        
        ax.xaxis.set_major_formatter (plt.FuncFormatter(format_thicks))
        """
        if ind % step ==0: 
            return fmt.format (ind)
        else: None 
        
    # show label every 'space'samples 
    if auto: 
        space = 10.
        step = int (np.ceil ( len(xtick_range)/ space )) 
        
    rotation = xlkws.get('rotation', 90 ) if 'rotation' in xlkws.keys (
        ) else xlkws.get('rotate_xlabel', 90 )
    
    if len(xtick_range) >= space :
        if ticks=='y': 
            ax.yaxis.set_major_formatter (plt.FuncFormatter(format_ticks))
        else: 
            ax.xaxis.set_major_formatter (plt.FuncFormatter(format_ticks))

        plt.setp(ax.get_yticklabels() if ticks=='y' else ax.get_xticklabels(), 
                 rotation = rotation )
    else: 
        
        # ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(3))
        # # ticks_loc = ax.get_xticks().tolist()
        # ax.xaxis.set_major_locator(mpl.ticker.FixedLocator(ticks_loc))
        # ax.set_xticklabels([fmt.format(x) for x in ticks_loc])
        tlst = [fmt.format(item) for item in  xtick_range]
        ax.set_yticklabels(tlst, **xlkws) if ticks=='y' \
            else ax.set_xticklabels(tlst, **xlkws) 
  
def _set_sns_style (s, /): 
    """ Set sns style whether boolean or string is given""" 
    s = str(s).lower()
    s = re.sub(r'true|none', 'darkgrid', s)
    return sns.set_style(s) 

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
        column2skip = is_in_if(X.columns, column2skip, return_diff= True)
        if len(column2skip) ==len (X.columns): 
            warnings.warn("Value(s) to skip are not detected.")
        if inplace : 
            X[column2skip] = np.log10 ( X[column2skip] ) 
            X.drop (columns =cskip , inplace =True )
            return  
        else : 
            X0[column2skip] = np.log10 ( X0[column2skip] ) 
            
    return X0
    
def plot_bar(x, y, wh= .8,  kind ='v', fig_size =(8, 6), savefig=None,
             xlabel =None, ylabel=None, fig_title=None, **bar_kws): 
    """
    Make a vertical or horizontal bar plot.

    The bars are positioned at x or y with the given alignment. Their dimensions 
    are given by width and height. The horizontal baseline is left (default 0)
    while the vertical baseline is bottom (default=0)
    
    Many parameters can take either a single value applying to all bars or a 
    sequence of values, one for each bar.
    
    Parameters 
    -----------
    x: float or array-like
        The x coordinates of the bars. is 'x' for vertical bar plot as `kind` 
        is set to ``v``(default) or `y` for horizontal bar plot as `kind` is 
        set to``h``. 
        See also align for the alignment of the bars to the coordinates.
    y: float or array-like
        The height(s) for vertical and width(s) for horizonatal of the bars.
    
    wh: float or array-like, default: 0.8
        The width(s) for vertical or height(s) for horizaontal of the bars.
        
    kind: str, ['vertical', 'horizontal'], default='vertical'
        The kind of bar plot. Can be the horizontal or vertical bar plots. 
    bar_kws: dict, 
        Additional keywords arguments passed to : 
            :func:`~matplotlib.pyplot.bar` or :func:`~matplotlib.pyplot.barh`. 
    """
    
    assert str(kind).lower().strip() in ("vertical", 'v',"horizontal", "h"), (
        "Support only the horizontal 'h' and vertical 'v' bar plots."
        " Got {kind!r}")
    kind =str(kind).lower().strip()
    
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize =fig_size)
    if kind in ("vertical", "v"): 
        ax.bar (x, height= y, width =  wh , **bar_kws)
    elif kind in ("horizontal", "h"): 
        ax.barh (x , width =y , height =wh, **bar_kws)
        
    ax.set_xlabel (xlabel )
    ax.set_ylabel(ylabel) 
    ax.set_title (fig_title)
    if savefig is not  None: 
        savefigure (fig, savefig, dpi = 300)
        
    plt.close () if savefig is not None else plt.show() 
  

def _format_ticks (value, tick_number, fmt ='S{:02}', nskip =7 ):
    """ Format thick parameter with 'FuncFormatter(func)'
    rather than using `axi.xaxis.set_major_locator (plt.MaxNLocator(3))`
    ax.xaxis.set_major_formatter (plt.FuncFormatter(format_thicks))
    
    :param value: tick range values for formatting 
    :param tick_number: number of ticks to format 
    :param fmt: str, default='S{:02}', kind of tick formatage 
    :param nskip: int, default =7, number of tick to skip 
    
    """
    if value % nskip==0: 
        return fmt.format(int(value)+ 1)
    else: None 
    
#XXX OPTIMIZE 
def plot_confidence (
    data = None, 
    *,  
    y=None, 
    x=None, 
    ci =.95 ,  
    kind ='line', 
    b_samples = 1000, 
    **sns_kws
    ): 
    """ Plot confidence data 
    
    Confidence Interval (CI)  is a type of estimate computed from the statistics 
    of the observed data which gives a range of values that’s likely to 
    contain a population parameter with a particular level of confidence.
    CI as a concept was put forth by Jerzy Neyman in a paper published 
    in 1937. There are various types of the confidence interval, some of 
    the most commonly used ones are: CI for mean, CI for the median, CI for 
    the difference between means, CI for a proportion and CI for the difference 
    in proportions.
    
    Parameters 
    ------------
    data: pandas.DataFrame, numpy.ndarray, mapping, or sequence
       Input data structure. Either a long-form collection of vectors 
       that can be assigned to named variables or a wide-form dataset 
       that will be internally reshaped.

    x, y: vectors or keys in data
       Variables that specify positions on the x and y axes.
       
    ci: float, default=.95 
       Confidence value. 
       
    kind: str, default='line' 
       kind of confidence intervval plot. 
      
    b_samples: int, default=1000
        Number of bootstraps to use for computing the confidence interval.
        
    sns_kws: dict, 
       Keywords arguments passed to the `sns.lineplot` or `sns.regplot`
       
    Returns 
    ----------
    ax: matplotlib.axes.Axes
       The matplotlib axes containing the plot.
       
    """   
    #y = np.array (y) 
    #x= x or ( np.arange (len(y)) if 
    ax=None 
    if 'lin' in str(kind).lower(): 
        ax = sns.lineplot(data= data, x=x, y=y, ci=ci, **sns_kws)
    elif 'reg' in  str(kind).lower(): 
        ax = sns.regplot(data = data, x=x, y=y, ci=ci, **sns_kws ) 
    else: 
        if not y: 
            raise ValueError("y should not be None when using the boostrapping"
                             " for plotting the confidence interval.")
        b_samples = _assert_all_types(
            b_samples, int, float, objname="Bootstrap samples `b_samples`")
        
        from sklearn.metrics import resample 
        # configure bootstrap
        n_iterations = 1000 # here k=no. of bootstrapped samples
        n_size = int(len(y))
          
        # run bootstrap
        medians = list()
        for i in range(n_iterations):
           s = resample(y, n_samples=n_size);
           m = np.median(s);
           medians.append(m)
          
        # plot scores
        plt.hist(medians)
        plt.show()
          
        # confidence intervals
        p = ((1.0-ci)/2.0) * 100
        lower =  np.percentile(medians, p)
        p = (ci+((1.0-ci)/2.0)) * 100
        upper =  np.percentile(medians, p)
  
        print(f"\n{ci*100} confidence interval {lower} and {upper}")
    
    return ax 

def plot_confidence_ellipse (x, y ): 
    """ Plot a confidence ellipse of a two-dimensional dataset 
    
    This function plots the confidence ellipse of the covariance of 
    the given array-like variables x and y. The ellipse is plotted 
    into the given axes-object ax.
    
    The approach that is used to obtain the correct geometry 
    is explained and proved here:
      https://carstenschelp.github.io/2018/09/14/Plot_Confidence_Ellipse_001.html
      
    The method avoids the use of an iterative eigen decomposition 
    algorithm and makes use of the fact that a normalized covariance 
    matrix (composed of pearson correlation coefficients and ones) is 
    particularly easy to handle.
    
    """
    fig, ax_nstd = plt.subplots(figsize=(6, 6))
    # dependency_nstd = [[0.8, 0.75],
    #                    [-0.2, 0.35]]
    mu = 0, 0
    # scale = 8, 5
    ax_nstd.axvline(c='grey', lw=1)
    ax_nstd.axhline(c='grey', lw=1)
    
    #x, y = get_correlated_dataset(500, dependency_nstd, mu, scale)
    ax_nstd.scatter(x, y, s=0.5)
    
    confidence_ellipse(x, y, ax_nstd, n_std=1,
                       label=r'$1\sigma$', edgecolor='firebrick')
    confidence_ellipse(x, y, ax_nstd, n_std=2,
                       label=r'$2\sigma$', edgecolor='fuchsia',
                       linestyle='--')
    confidence_ellipse(x, y, ax_nstd, n_std=3,
                       label=r'$3\sigma$', edgecolor='blue', 
                       linestyle=':')
    
    ax_nstd.scatter(mu[0], mu[1], c='red', s=3)
    ax_nstd.set_title('Different standard deviations')
    ax_nstd.legend()
    plt.show()
    
def confidence_ellipse(
    x, 
    y, 
    ax, 
    n_std=3.0, 
    facecolor='none', 
    **kwargs
    ):
    """
    Create a plot of the covariance confidence ellipse of *x* and *y*.

    Parameters
    ----------
    x, y : array-like, shape (n, )
        Input data.

    ax : matplotlib.axes.Axes
        The axes object to draw the ellipse into.

    n_std : float
        The number of standard deviations to determine the ellipse's radiuses.

    **kwargs
        Forwarded to `~matplotlib.patches.Ellipse`

    Returns
    -------
    mpl.patches.Ellipse
    """
    if x.size != y.size:
        raise ValueError("x and y must be the same size")

    cov = np.cov(x, y)
    pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
    # Using a special case to obtain the eigenvalues of this
    # two-dimensional dataset.
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = Ellipse((0, 0), width=ell_radius_x * 2, height=ell_radius_y * 2,
                      facecolor=facecolor, **kwargs)

    # Calculating the standard deviation of x from
    # the squareroot of the variance and multiplying
    # with the given number of standard deviations.
    scale_x = np.sqrt(cov[0, 0]) * n_std
    mean_x = np.mean(x)

    # calculating the standard deviation of y ...
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_y = np.mean(y)

    transf = transforms.Affine2D() \
        .rotate_deg(45) \
        .scale(scale_x, scale_y) \
        .translate(mean_x, mean_y)

    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)  
   
def plot_strike (
    list_of_edis, /, 
    kind = 2,
    period_tolerance=.05, 
    text_pad =1.65 , 
    rot_z=0. , 
    **kws
    ): 
    
    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    from mtpy.imaging.plotstrike import PlotStrike
    from ..property import IsEdi
    from .context import nullify_output
    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    if isinstance ( list_of_edis, str): 
        if os.path.isdir ( list_of_edis ): 
            list_of_edis = [
                os.path.join(list_of_edis,  f) for f in os.listdir (list_of_edis)
                            if str(f).lower().endswith ('.edi')]
        elif os.path.isfile (list_of_edis): 
            list_of_edis =[list_of_edis ]
          
    # now check whether is valid EDI     
    # list comprehension faster than
    # tuple (map (lambda f:  IsEdi._assert_edi (f ), list_of_edis ) )
    [ IsEdi._assert_edi (f ) for f in list_of_edis ] 
    # suppress third party verbosity 

    with nullify_output():
        PlotStrike(
            fn_list=list_of_edis,
            plot_type=kind,
            **kws
            )
        
plot_strike.__doc__="""\ 
Plot the strike estimated from the invariants and phase tensor. 
in a rose diagram of xy plot. 

Parameters 
------------
list_of_edis: list, 
    full paths to .edi files to plot or list of :term:`EDI` files. 
    
   .. versionchanged:: 0.2.0 
      No need to provide a list of term:`EDI` files. Henceforth `list_of_edis`
      accepts the EDI path-like object of single EDI file then asserts 
      the validity of the EDI files afterward. 

kind: int, default=2 
    Can be [ 1 | 2 ] where:
       
   - *1* to plot individual decades in one plot
   - *2* to plot all period ranges into one polar diagram for each 
     strike angle estimation
         
   One could try also plot_type = 1 to plot by decade
                                                                   
fig_num: int, default=1, 
   figure number to be plotted. *Default* is 1
        
font_size: float, default=10, 
   Figure size 
   
rot_z: float, default=0., 
   angle of rotation clockwise positive. 

period_tolerance: float, default=.05 
    Tolerance level to match periods from different edi files.
    *Default* is 0.05

text_pad: float, default=1.65
   padding of the angle label at the bottom of each
                 polar diagram.  *Default* is 1.65
plot_range:  str, tuple 
   The period range to estimate the strike angle. It can be 
   [ 'data' | (period_min,period_max) ].  Options are:
       
   * *'data'* for estimating the strike for all periods
     in the data.
   * (pmin,pmax) for period min and period max, input as
     (log10(pmin),log10(pmax))

plot_tipper: [ True | False ]
    - True to plot the tipper strike
    - False to not plot tipper strike
   
pt_error_floor: int, optional 
   Maximum error in degrees that is allowed to 
   estimate strike. *Default* is None allowing all 
   estimates to be used.

fold: [ True | False ]
   * True to plot only from 0 to 180
   * False to plot from 0 to 360
            
plot_orthogonal: [ True | False]
    * True to plot the orthogonal strike directions
    * False to not

color: [ True | False ]
    * True to plot shade colors
    * False to plot all in one color
              
color_inv:str, 
   color of invariants plots

color_pt: str, 
   color of phase tensor plots

color_tip: str 
   color of tipper plots

ring_spacing: float, optional 
    spacing of rings in polar plots

ring_limits: tuple of int, 
   plot limits (min count, max count) set each plot have these limits 
                    
plot_orientation: str, [ 'h' | 'v' ] 
   horizontal or vertical plots


See More
--------
Plots the strike angle as determined by invariants of the impedance tensor 
(Weaver et al. [2003] [1]_) and phase tensor azimuth 
(Caldwell et al. [2004] [2]_) 

The data is split into decades where the histogram for each is plotted in
the form of a rose diagram with a range of 0 to 180 degrees.
Where 0 is North and 90 is East.  The median angle of the period band is
set in polar diagram.  The top row is the strike estimated from
the invariants of the impedance tensor.  The bottom row is the azimuth
estimated from the phase tensor.  If tipper is 'y' then the 3rd row is the
strike determined from the tipper, which is orthogonal to the induction
arrow direction.

References 
----------
.. [1] Weaver J.T, Lilley F.E.M.(2003)  Invariants of rotation of axes and indicators of
       dimensionality in magnetotellurics, Australian National University,
       University of Victoria; http://bib.gfz-potsdam.de/emtf/2007/pdf/lilley.pdf
.. [2] T. Grant Caldwell, Hugh M. Bibby, Colin Brown, The magnetotelluric phase tensor, 
       Geophysical Journal International, Volume 158, Issue 2, August 2004, 
       Pages 457–469, https://doi.org/10.1111/j.1365-246X.2004.02281.x
Examples 
----------
>>> import os 
>>> from watex.datasets import fetch_data 
>>> from watex.utils.plotutils import plot_strike 
>>> from watex.datasets._io import get_data # get edidata in cache  
>>> fetch_data ( 'huayuan', samples = 25 ) # store edi in cache 
>>> # get the edi in cache and plotStrike 
>>> edi_fn_lst = [os.path.join(get_data(),ff) for ff in os.listdir(get_data()) 
...         if ff.endswith('.edi')] 
>>> plot_strike(edi_fn_lst )

"""
def plot_text (
    x, y, 
    text=None , 
    data =None, 
    coerce =False, 
    basename ='S', 
    fig_size =( 7, 7 ), 
    show_line =False, 
    step = None , 
    xlabel ='', 
    ylabel ='', 
    color= 'k', 
    mcolor='k', 
    lcolor=None, 
    show_leg =False,
    linelabel='', 
    markerlabel='', 
    ax=None, 
    **text_kws
    ): 
    """ Plot text(s) indicating each position in the line. 
    
    Parameters 
    -----------
    x, y: str, float, Array-like 
        The position to place the text. By default, this is in data 
        coordinates. The coordinate system can be changed using the 
        transform parameter.
        
    text: str, 
        The text
        
    data: pd.DataFrame, 
       Data containing x and y names. Need to be supplied when x and y 
       are given as string names. 
       
    coerce:bool, default=False 
       Force the plot despite the given textes do not match the number of  
       positions `x` and `y`. If ``False``, number of positions must be 
       consistent with x and y, otherwise error raises. 
       
    basename: str, default='S' 
       the text to prefix the position when the text is not given. 
       
    fig_size: tuple, default=(7, 7) 
       Matplotlib figure size.
       
    show_line: bool, default=False 
       Display the line from x, y. 
       
    step: int,Optional 
       The number of intermediate positions to skip in the plotting text. 
       
    xlabel, ylabel: str, Optional, 
       The labels of x and y. 
       
    color: str, default='k', 
       Text color.
       
    mcolor: str, default='k', 
       Marker color. 
       
    lcolor: str, Optional 
       Line color if `show_line` is set to ``True``. 
       
    show_leg: bool, default=False 
       Display the legend of line and marker labels. 
       
    linelabel, markerlabel: str, Optional 
        The labels of the line and marker. 
       
    ax: Matplotlib.Axes, optional 
       Support plot to another axes 
       
       .. versionadded:: 0.2.5 
       
    text_kws: dict, 
       Keyword arguments passed to :meth:`matplotlib.axes.Axes.text`. 

    Return 
    -------
    ax: Matplotlib axes 
    
    Examples 
    --------
    >>> import watex as wx 
    >>> data =wx.make_erp (as_frame =True, n_stations= 7 )
    >>> x , y =[ 0, 1, 3 ], [2, 3, 6] 
    >>> texto = ['AMT-E1147', 'AMT-E1148',  'AMT-E180']
    >>> plot_text (x, y , text = texto)# no need to set  coerce, same length 
    >>> data =wx.make_erp (as_frame =True, n_stations= 20 )
    >>> x , y = data.easting, data.northing
    >>> text1 = ['AMT-E1147', 'AMT-E1148',  'AMT-E180'] 
    >>> plot_text (x, y , coerce =True , text = text1 , show_leg= True, 
                   show_line=True, linelabel='E1-line', markerlabel= 'Site', 
               basename ='AMT-E0' 
               )
    """
    # assume x, y  series are passed 
    if isinstance(x, str) or hasattr ( x, 'name'): 
        xlabel = x  if isinstance(x, str) else x.name 
        
    if isinstance(y, str) or hasattr ( y, 'name'): 
        ylabel = y  if isinstance(y, str) else y.name 
        
    if x is None and  y is None:
        raise TypeError("x and y are needed for text plot. NoneType"
                        " cannot be plotted.")    
        
    x, y = assert_xy_in(x, y, data = data ) 

    if text is None and not coerce: 
       raise TypeError ("Text cannot be plotted. To force plotting text with"
                        " the basename, set ``coerce=True``.")

    text = is_iterable(text , exclude_string= True , transform =True )
    
    if ( len(text) != len(y) 
        and not coerce) : 
        raise ValueError("In principle text array and x/y must be consistent."
                         f" Got {len(text)} and {len(y)}. To plot anyway,"
                         " set ``coerce=True``.")
    if coerce : 
        basename =str(basename)
        text += [f'{basename}{i+len(text):02}' for i in range (len(y) )]

    if step is not None: 
        step = _assert_all_types(step , float, int , objname ='Step') 
        for ii in range(len(text)): 
            if not ii% step ==0: 
                text[ii]=''

    if ax is None: 
        
        fig, ax = plt.subplots(1,1, figsize =fig_size)
    
    # plot = ax.scatter if show_line else ax.plot 
    ax_m = None 
    if show_line: 
        ax.plot (x, y , label = linelabel, color =lcolor 
                 ) 
        
    for ix, iy , name in zip (x, y, text ): 
        ax.text ( ix , iy , name , color = color,  **text_kws)
        if name !='':
           ax_m  = ax.scatter ( [ix], [iy] , marker ='o', color =mcolor, 
                       )
  
    ax.set_xlabel (xlabel)
    ax.set_ylabel (ylabel) 
    
    ax_m.set_label ( markerlabel) if ax_m is not None else None 
    
    if show_leg : 
        ax.legend () 
        
    return ax 


 
def _make_axe_multiple ( n, ncols = 3 , fig_size =None, fig =None, ax= ... ): 
    """ Make multiple subplot axes from number of objects. """
    if is_iterable (n): 
       n = len(n) 
     
    nrows = n // ncols + ( n % ncols ) 
    if nrows ==0: 
       nrows =1 
       
    if ax in ( ... , None) : 
        fig, ax = plt.subplots (nrows, ncols, figsize = fig_size )  
    
    return fig , ax 


  
    
  
    
  
    
  
    
  
    
  
    
  
    
  
    
  
    
  
    
  
    
  
    
  
    
  
    
  
    
  
    
  
    
  
    
    
    
    
