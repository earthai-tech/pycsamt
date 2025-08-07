from __future__ import annotations 
import re
import warnings
import inspect 
from abc import ABCMeta 
import copy 
import numpy as np 
import pandas as pd
import seaborn as sns 
from scipy.cluster.hierarchy import dendrogram 

import matplotlib as mpl 
import matplotlib.pyplot  as plt
import matplotlib.ticker as mticker
from matplotlib import cm 
from matplotlib.colors import BoundaryNorm

from .base import BasePlot 
from ..utils.generic import make_ids 
from ..utils.validation import _check_consistency_size
from ..utils.plot import _get_xticks_formatage

# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# create a shadow class to hold the font and matplotlib properties
# from 'EvalPlot` and giving an option for saving figure

pobj = type ('Plot', (BasePlot, ), {}) 
# setattr(pobj, 'save', _b.save )
# redefine the pobj doc 
pobj.__doc__="""\
Shadow plotting class that holds the :class:`~watex.property.BasePlot`
parameters. 

Each matplotlib properties can be modified as  :class:`~watex.view.pobj`
attributes object. For instance:: 
    
    >>> pobj.ls ='-.' # change the line style 
    >>> pobj.fig_Size = (7, 5) # change the figure size 
    >>> pobj.lw=7. # change the linewidth 
    
.. seealso:: 
    
    Refer to :class:`~watex.property.BasePlot` for parameter details. 
    
"""
  
def save (self, fig): 
    """ savefigure if figure properties are given. """
    if self.savefig is not None: 
        fig.savefig (self.savefig,dpi = self.fig_dpi , 
                     bbox_inches = 'tight', 
                     orientation=self.fig_orientation 
                     )
    plt.show() if self.savefig is None else plt.close () 
        
def plot2d(
    ar, 
    y=None,  
    x =None,  
    distance=50., 
    stnlist =None, 
    prefix ='S', 
    how= 'py',
    to_log10=False, 
    plot_contours=False,
    top_label='', 
    **baseplot_kws
    ): 
    """Two dimensional template for visualization matrices.
    
    It is a wrappers that can plot any matrice by customizing the position 
    X and y. By default X is considering as stations  and y the resistivity 
    log data. 
    
    Parameters 
    -----------
    ar: Array-like 2D, shape (M, N) 
        2D array for plotting. For instance, it can be a 2D resistivity 
        collected at all stations (N) and all frequency (M) 
    y: array-like, default=None
        Y-coordinates. It should have the length N, the same of the ``arr2d``.
        the rows of the ``arr2d``.
    x: array-like, default=None,  
        X-coordinates. It should have the length M, the same of the ``arr2d``; 
        the columns of the 2D dimensional array.  Note that if `x` is 
        given, the `distance is not needed. 

    distance: float 
        The step between two stations. If given, it creates an array of  
        position for plotting purpose. Default value is ``50`` meters. 
        
    stnlist: list of str 
        List of stations names. If given,  it should have the same length of 
        the columns M, of `arr2d`` 
       
    prefix: str 
        string value to add as prefix of given id. Prefix can be the site 
        name. Default is ``S``. 
        
    how: str 
        Mode to index the station. Default is 'Python indexing' i.e. 
        the counting of stations would starts by 0. Any other mode will 
        start the counting by 1.
     
    to_log10: bool, default=False 
       Recompute the `ar`  in logarithm  base 10 values. Note when ``True``, 
       the ``y`` should be also in log10. 
    plot_contours: bool, default=True 
       Plot the contours map. Is available only if the plot_style is set to 
       ``pcolormesh``. 
       
    top_label: str, 
       Name of the top label. 
       
    baseplot_kws: dict, 
       All all  the keywords arguments passed to the property  
       :class:`watex.property.BasePlot` class. 
       
    Returns 
    -------
    axe: <AxesSubplot> object 
    
    Examples 
    -------- 
    >>> import numpy as np
    >>> import watex 
    >>> np.random.seed (42) 
    >>> data = np.random.randn ( 15, 20 )
    >>> data_nan = data.copy() 
    >>> data_nan [2, 1] = np.nan; data_nan[4, 2]= np.nan;  data_nan[6, 3]=np.nan
    >>> watex.view.mlplot.plot2d (data )
    <AxesSubplot:xlabel='Distance(m)', ylabel='log10(Frequency)[Hz]'>
    >>> watex.view.mlplot.plot2d (data_nan ,  plt_style = 'imshow', 
                                  fig_size = (10, 4))
    """
    #xxxxxxxxx update base plot keyword arguments
    for k  in list(baseplot_kws.keys()): 
        setattr (pobj , k, baseplot_kws[k])
        
    if y is not None: 
        if len(y) != ar.shape [0]: 
            raise ValueError ("'y' array must have an identical number " 
                              f" of row of 2D array: {ar.shape[0]}")
            
    if x is not None: 
        if len(x) != ar.shape[1]: 
            raise ValueError (" 'x' array must have the same number " 
                              f" of columns of 2D array: {ar.shape[1]}")

    d= distance or 1.
    try : 
         distance = float(distance) 
    except : 
        raise TypeError (
             f'Expect a float value not {type(distance).__name__!r}')
        
    # put value to log10 if True 
    if to_log10: 
        ar = np.log10 (ar ) # assume the resistivity data 
        y = np.log10(y) if y is not None else y # assume the frequency data 

    y = np.arange(ar.shape [0]) if y is None else y 
    x=  x  or np.arange(ar.shape[1]) * d
         
    stn = stnlist or make_ids ( x , prefix , how = how) 
    #print(stnlis)
    if stn is not None: 
        stn = np.array(stn)
        
    if not _check_consistency_size(stn, x, error ="ignore"): 
        raise ValueError("The list of stations and positions must be"
                         f" consistent. {len(stnlist)} and {len(x)}"
                         " were given respectively")
            
    # make figure 
    fig, axe = plt.subplots(1,figsize = pobj.fig_size, 
                            num = pobj.fig_num,
                            dpi = pobj.fig_dpi
                            )
    
    cmap = plt.get_cmap( pobj.cmap)
    
    if pobj.plt_style not in ('pcolormesh','imshow' ): 
        warnings.warn(f"Unrecognized plot style {pobj.plt_style!r}."
                      " Expect ['pcolormesh'|'imshow']."
                      " 'pcolormesh' ( default) is used instead.")
        pobj.plt_style= 'pcolormesh'
        
    if pobj.plt_style =='pcolormesh': 
        X, Y = np.meshgrid (x, y)
        # ar = np.ma.masked_where(np.isnan(ar), ar)
        #Zm = ma.array(Z,mask=np.isnan(Z))
        pkws = dict (vmax = np.nanmax (ar),
                     vmin = np.nanmin (ar), 
                     ) 
        
        if plot_contours: 
            levels = mticker.MaxNLocator(nbins=15).tick_values(
                    np.nanmin (ar), np.nanmax(ar) )
            # delete vmin and Vmax : not supported 
            # when norm is passed 
            del pkws ['vmin'] ; del pkws ['vmax']
            pkws ['norm'] = BoundaryNorm(
                levels, ncolors=plt.colormaps[pobj.cmap].N, clip=True)
            
        
        ax = axe.pcolormesh ( X, Y, np.flipud (ar),
                    shading= pobj.plt_shading, 
                    cmap =cmap, 
                    **pkws 
            )
        if plot_contours: 
             # contours are *point* based plots, so convert 
             # our bound into point centers
            dx, dy = 0.05, 0.05
            axe.contourf(X+ dx/2.,
                         Y + dy/2., np.flipud (ar) , levels=levels,
                         cmap=plt.colormaps[pobj.cmap]
                         )
    if pobj.plt_style =='imshow': 
        ax = axe.imshow (ar,
                    interpolation = pobj.imshow_interp, 
                    cmap =cmap,
                    aspect = pobj.fig_aspect ,
                    origin= 'lower', 
                    extent=(  np.nanmin(x),
                              np.nanmax (x), 
                              np.nanmin(y), 
                              np.nanmax(y)
                              )
            )
    # set axis limit 
    axe.set_ylim(np.nanmin(y), 
                 np.nanmax(y))
    axe.set_xlim(np.nanmin(x), 
                 np.nanmax (x))

    cbl = 'log_{10}' if to_log10 else ''
    axe.set_xlabel(pobj.xlabel or 'Distance(m)', 
                 fontdict ={
                  'size': 1.5 * pobj.font_size ,
                  'weight': pobj.font_weight}
                 )
      
    axe.set_ylabel(pobj.ylabel or  f"{cbl}Frequency$[Hz]$",
             fontdict ={
                     #'style': pobj.font_style, 
                    'size':  1.5 * pobj.font_size ,
                    'weight': pobj.font_weight})
    if pobj.show_grid is True : 
        axe.minorticks_on()
        axe.grid(color='k', ls=':', lw =0.25, alpha=0.7, 
                     which ='major')
    
   
    labex = pobj.cb_label or f"{cbl}App.Res$[Ω.m]$" 
    
    cb = fig.colorbar(ax , ax= axe)
    cb.ax.yaxis.tick_left()
    cb.ax.tick_params(axis='y', direction='in', pad=2., 
                      labelsize = pobj.font_size )
    
    cb.set_label(labex,fontdict={'size': 1.2 * pobj.font_size ,
                              'style':pobj.font_style})
    #--> set second axis 
    axe2 = axe.twiny() 
    axe2.set_xticks(range(len(x)),minor=False )
    
    # set ticks params to reformat the size 
    axe.tick_params (  labelsize = pobj.font_size )
    axe2.tick_params (  labelsize = pobj.font_size )
    # get xticks and format labels using the auto detection 
    _get_xticks_formatage(axe2, stn, fmt = 'S{:02}',  auto=True, 
                          rotation=pobj.rotate_xlabel )
    
    axe2.set_xlabel(top_label, fontdict ={
        'style': pobj.font_style,
        'size': 1.5 * pobj.font_size ,
        'weight': pobj.font_weight}, )
      
    fig.suptitle(pobj.fig_title,ha='left',
                 fontsize= 15* pobj.fs, 
                 verticalalignment='center', 
                 style =pobj.font_style,
                 bbox =dict(boxstyle='round',
                            facecolor ='moccasin')
                 )
   
    #plt.tight_layout(h_pad =1.8, w_pad =2*1.08)
    plt.tight_layout()  
    if pobj.savefig is not None :
        fig.savefig(pobj.savefig, dpi = pobj.fig_dpi,
                    orientation =pobj.orient)
 
    plt.show() if pobj.savefig is None else plt.close(fig=fig) 
    
    
    return axe        