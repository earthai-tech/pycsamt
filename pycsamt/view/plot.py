# -*- coding: utf-8 -*-
#       Created on Mon Dec 28 14:28:06 2020
#       Author: Kouadio K.Laurent<etanoyau@gmail.com>
#       Licence: LGPL
"""
Module View
============

Some templates for plotting purposes. It gives a quick alternative for uses to 
save their time for writting their own plot scripts. However to have full 
control of the plot, it is recommended to write your own plot scripts. 

Note that the package can not handle all the plots possibilty that offers the
software. In the future, many plots should be removed and replaced by some  
most efficient ones. 
 
"""

import os 
import re
import copy
import warnings 
import numpy as np

import matplotlib as mpl 
import matplotlib.pyplot  as plt
import matplotlib.cm as cm 
import matplotlib.colorbar as mplcb
import matplotlib.gridspec as gspec

import pycsamt.utils.exceptions as CSex
import pycsamt.utils.func_utils as func
import pycsamt.utils.plot_utils as mplotus  
import pycsamt.utils.zcalculator  as Zcc
from pycsamt.core import avg as CSMATavg 
from pycsamt.core.cs import CSAMT
from pycsamt.site import Profile
from pycsamt.modeling import occam2d
from pycsamt.geodrill import geocore  as geoD
from  pycsamt.utils.plotdecorator import geoplot1d 
from pycsamt.utils._p import suit 
from pycsamt.utils._csamtpylog import csamtpylog 
from pycsamt.processing import Processing 

try : 
    from pycsamt.__init__ import imtpy  
    if imtpy : 
        from mtpy.imaging.phase_tensor_pseudosection import (
            PlotPhaseTensorPseudoSection ) 
        from mtpy.modeling.modem import ( 
            plot_response, 
            plot_rms_maps, 
            plot_slices, 
            phase_tensor_maps
            
        )
        
except:
    imtpy = False 

_logger=csamtpylog.get_csamtpy_logger(__name__)

warnings.warn("In the future, many plots of the module 'view' should be"
              " removed and replaced by the most efficient ones.", 
              FutureWarning
              )


class Plot1d :
    """
    Plot One-dimensional. It offers a quick visualisation of one-dimensional plot. 
    
    The methods functions are some templates for plotting puposes. They  are not 
    exhaustive. Users can edit and customize each plot according its will.
    
    Hold other informations of attributes: 

    ==================  =======================================================
    Key Words           Description        
    ==================  =======================================================
    fig_dpi             dots-per-inch resolution of the figure
                        *default* is 300
    fig_num             number of the figure instance
                        *default* is 'Mesh'
    fig_size            size of figure in inches (width, height)
                        *default* is [5, 5]
    fs                  size of font of axis tick labels, axis labels are
                        fs+2. *default* is 6 
    ls                  [ '-' | '.' | ':' ] line style of mesh lines
                        *default* is '-'
    marker              marker. *default* is ``o``.
    ms                  size of marker in points. *default* is 3.
    x_minorticks        minorticks in x-abscissa. default is ``1``. 
    y_minorticks        minortick in y-coordinates . *default* is ``1`` 
    font_size           text font size. Default is ``3``
    font_style          text font style. *default* is ``italic``
    fs                  value for adjusting all `size`in plot. 
                        *default* is ``2``
    mplfont             matplotlib font for rcparameters. *default* is``cursive``
    markerfacecolor     facecolor of the marker.*default* is ``k``. 
    markeredgecolor     edge color of the marker point. *default* is ``k``. 
    fontweight          weight type of the font. *default* is ``italic``.
    lw                  width of the plot line. *default* is ``1.5``.
    alpha               control the transparency of grid or line. 
                        *default* is ``.5``.
    xlim                Bound of x -coordinates (abscissa). Should be a tuple 
                        of e.g.(0, 4 ) where 0 is the start and the end points 
                        respectively. 
    ylim                Bound of y -coordinates (up). Should be a tuple 
                        of e.g.(0, 4 ) where 0 is the start and the end points 
                        respectively.  
    xlabel              label name along the x-coordinates. 
    ylabel              label name along the x-coordinates.
    savefig             figure name or path-like object for saving. 
    ==================  =======================================================
    
    """
    
    def __init__(self,**kwargs):
        
        self._logging =csamtpylog.get_csamtpy_logger(self.__class__.__name__)
  
        self.fig_num= kwargs.pop('fig_num', 1)
        self.fig_size = kwargs.pop('fig_size', [12,8])
        self.fig_dpi =kwargs.pop('fig_dpi', 300)
        
        self.fig_title =kwargs.pop('fig_title', None)
        
        self.x_minorticks=kwargs.pop('xminorticks', 1)
        self.y_minorticks =kwargs.pop('yminorticks', 1)
        
        self.font_size =kwargs.pop('font_size',3.)
        self.font_style=kwargs.pop('font_style', 'italic')
        self.fs =kwargs.pop('fs', 2.)
        
        self.marker =kwargs.pop('maker', 'o')
        self.ms =kwargs.pop('ms', 3)
        self.mplfont =kwargs.pop('font','cursive')
        self.markerfacecolor=kwargs.pop('markerfacecolor', 'k' )
        self.markeredgecolor=kwargs.pop('markeredgecolor', 'k' )
        
        self.fontweight =kwargs.pop('font_weight', 'bold')
        self.ls= kwargs.pop('ls', '-')
        self.lw =kwargs.pop('lw', 1.5)
        self.alpha = kwargs.pop('alpha', .5)
        
        self.xlim =None 
        self.ylim=None 
        self.xlabel = kwargs.pop('xlabel', None) 
        self.ylabel = kwargs.pop('ylabel', None) 
        
        self.subplot_wspace = kwargs.pop('subplot_wspace', .3)
        self.subplot_hspace = kwargs.pop('subplot_hspace', .0)
        self.subplot_right = kwargs.pop('subplot_right', .98)
        self.subplot_left = kwargs.pop('subplot_left', .08)
        self.subplot_top = kwargs.pop('subplot_top', .85)
        self.subplot_bottom = kwargs.pop('subplot_bottom', .1)
        
        self.savefig = kwargs.pop('savefig', None)
        
    
    def plotoverlay1d (self, arr, freqs,  ns= 0, leg = None,
                      show_grid =False, **kws ) : 
        """ Plot the data at each station above to one to another. 
        
        This is useful to see the missing value in the data at each frequency. 
        For |AMT| survey, the plot can show the weak frequency and also the 
        dead band data. 
        
        Parameters 
        -----------
        arr : ndarray , shape (N, M) 
            2D array for plotting. For instance, it can be a 2D resistivity 
            collected at all stations (N) and all frequency (M) 
        freqs: array-like 
            Frequency values to designate the X-coordinates. It should have 
            the length N, the same of the ``arr2d`` i.e. the rows of the 
            ``arr``. Can be a frequency 
        ns: list of int  
            The list of indices of each stations. Start index at 0. 
        leg: list 
            legend labels put into a list. IOt must match the number of given 
            resistivity values. 
        show_grid: bool, 
            Draw a grid space for the plot. User can customize the grid plot 
            by him/herself.
        kws : dict 
            Matplotlib logs functions additional keywords arguments. 

        Returns 
        ------- 
        ax: Matplotlib.pyplot <AxesSubplot> 
        
        Examples 
        --------
        >>> from pycsamt.processing import make2d, get_ediObjs 
        >>> from pycsamt.processing import get_full_frequency 
        >>> from pycsamt.view import Plot1d 
        >>> edipath = 'data/3edis' 
        >>> ediObjs = get_ediObjs(edipath)
        >>> f =get_full_frequency(ediObjs) 
        >>> zxy = make2d (ediObjs, 'zxy', kind ='modulus')
        >>> p1d= Plot1d (fig_size =(5, 4), fig_title ='$Modulus |Z_{xy}|$')
        >>> p1d.xlabel ='$Log_{10}Frequency [Hz]$'
        >>> p1d.ylabel = '$|Z_{xy}|$'
        >>> ns =[2, 16, 28, 40]
        >>> leg =[f'S{int(ii):02}' for ii in ns]
        >>> p1d.plotoverlay1d(zxy, freqs=f)
       
        """
        
        arr = np.array(arr) 
        
        if len(freqs) != len(arr) : 
            raise ValueError ("Frequency and number of rows in array must have"
                              f" the same length. {len(freqs)} and {len(arr)}" 
                              " were given respectively.")
        
        if arr.ndim ==1: 
           arr = arr[: , np.newaxis]
           ns= [0]
            
        mask = np.isnan (arr) 
        # replace each nan if exist by 0 before cumsum 
        ar = arr.copy() 
        ar[ mask] = 0. 
        ar = np.cumsum(ar, axis =1 )
        # now reverse input the NaN in the array 
        ar [mask] = np.nan 

        fig, ax = plt.subplots(1,figsize = self.fig_size, num = self.fig_num,
                                dpi = self.fig_dpi , 
                    )

        if not isinstance(ns, (list, tuple, np.ndarray)):
            ns =[ns]

        try : 
            ns = list(map(int, ns))
        except : 
            raise TypeError ("Could not convert station numbers to integer." 
                             f" Got {type(ns).__name__!r}")
            
        if max(ns) >= len(ar): 
            raise ValueError(f"Too many station numbers {max(ns)}. Expect a "
                             f"maximum station indice equals to {len(ar)-1}")
        for ii in ns: 
            ax.loglog (freqs , ar[:, ii] , 
                        ls = self.ls ,
                        marker = self.marker, 
                        markersize = self.ms ,
                        markeredgecolor = self.markeredgecolor,
                        markerfacecolor = self.markerfacecolor ,
                       linewidth = self.lw, 
                       **kws)
        
        ax.legend( list(map(str, ns)) if leg is None else leg , 
                  loc ='best' , ncol = len(ns))
        
        ax.set_ylabel (self.ylabel  or '$ Rhoa [ \Omega.m]$', 
                        ) 
        ax.set_xlabel(self.xlabel or '$Log10Frequency[Hz]$', 
                      )
        if show_grid :
            ax.grid (visible =True , alpha =self.alpha,
                     which ='both', color ='gray')
            
        plt.tight_layout()   
        
        if self.savefig is not None: 
             plt.savefig(self.savefig , dpi = self.fig_dpi)
             plt.close (fig =fig ) 
            
        return ax
    
    def plotf1d (self, *args,  leg =None, kind ='semilogy', 
                      show_grid =False, distance =50., 
                     stnlist =None, prefix ='S', how='py',  **kws ):
        """ Plot 1D filtered plot template  stations Vs Rhoa filtered. 
        
        User can edit the template to customize its plot. 
        
        Parameters 
        -----------
        *args : args : list 
            Matplotlib logs funtions plot arguments 
            
        x: array-like 
            X-coordinates. It should have the length M, the same of the ``arr2d``; 
            the columns of the 2D dimensional array.  Note that if `x` is 
            given, the `distance is not needed. 
            
        kind: str 
             Type of log plot. Can be ``loglog``, ``semilogx``, semilogy`` or 
             merely a base plot function. Default is ``semilogy``.

        distance: float 
            The step between two stations. If given, it creates an array of  
            position for plotting purpose. Default value is ``50`` meters. 
            
        stnlist: list of str 
            List of stations names. If given,  it should have the same length of 
            the columns M, of `arr2d`` 
       
        prefix: str 
            string value to add as prefix of given id. Prefix can be the site 
            name. Default is ``S``. 
            
        show_grid: bool, 
            Draw a grid space for the plot. User can customize the grid plot 
            by him/herself. 
            
        how: str 
            Mode to index the station. Default is 'Python indexing' i.e. 
            the counting of stations would starts by 0. Any other mode will 
            start the counting by 1.
            
        kws : dict 
            Additional keywords arguments of Matplotlib subsplots function  
            :func:`plt.loglog` or :func:`plt.semilog`
        
        Returns 
        ------- 
        ax: Matplotlib.pyplot <AxesSubplot> 
            
        Examples 
        --------
        >>> import matplotlib.pyplot as plt 
        >>> from pycsamt.view import Plot1d 
        >>> from pycsamt.processing import get_ediObjs, make2d, get_full_frequency 
        >>> from pycsamt.processing import ama , interpolate2d
        >>> edipath = 'data/3edis'
        >>> ediObjs = get_ediObjs(edipath) 
        >>> res2d_raw = make2d (ediObjs, 'resyx')
        >>> res2d = interpolate2d (make2d (ediObjs, 'resyx')) 
        >>> freqs = get_full_frequency(ediObjs)
        >>> filtered =ama(res2d =res2d, phs2d =phs2d, freqs= freqs , window_size=5)
        >>> #Plot the resistivy value of the third frequency  at all stations 
        >>> figtitle = 'Data filtered with Adapative moving average'
        >>> leg= ['dirty data', 'filtered data'] 
        >>> p1d = Plot1d (fig_size =(5, 3), xlabel = 'Stations',
                  ylabel='$ Resistivity [ \Omega.m]$', fig_title= figtitle )
        >>> p1d.plotf1d (x, res2d_raw[3, :] , '+:', x, res2d[3, :], 
                         'ok--',x, filtered[3, :], 'ob--', kind ='semilogy',
                         show_grid =False, leg=leg)
        """
        x, y, *_ =  args 
        
        fig, ax = plt.subplots(1,figsize = self.fig_size, num = self.fig_num,
                                dpi = self.fig_dpi , 
                    )
        
    
        axplot = ax.plot 
        if kind =='loglog' : 
            axplot = ax.loglog 
        elif kind =='semilogx': 
            axplot = ax.semilogx 
        elif kind =='semilogy': 
            axplot = ax.semilogy 
            
        axplot(*args, **kws, 
               markersize = self.ms, 
               lw = self.lw 
               )
        
        ax.legend( [] if leg is None else leg , loc ='best' , ncol = 4 )
        
        d= distance or 1.
        x= x if x is not None else np.arange(len(x))  * d 
        
        if stnlist is not None: 
            stnlist = func._assert_all_types(stnlist, tuple, list, np.ndarray)
            if len(stnlist)!= len(x):
                raise ValueError("Expect  the length of the list of station names"
                                 f" to be the same like to position {len(x)}")
                
        stn = stnlist or func.make_ids ( x , prefix , how = how) 
        
        ax.set_xticks(ticks= x, minor=False )
        ax.set_xticklabels(stn, rotation=90., 
                           fontdict ={'style': self.font_style, 
                          'size':  1.5 * self.font_size ,
                          'weight': self.fontweight},  
                                )
        ax.set_ylabel (self.ylabel  or '$ Rhoa [ \Omega.m]$') 
        
        ax.set_xlabel(self.xlabel or 'Stations', 
                    #   fontdict ={'style': self.font_style, 
                    # 'size': 1.5 * self.font_size ,
                    # 'weight':  self.fontweight}, 
                      )
        
        if show_grid :
            ax.grid (visible =True , alpha =0.8, which ='both', color ='gray')
            
        plt.tight_layout()   
        
        if self.savefig is not None: 
             plt.savefig(self.savefig , dpi = self.fig_dpi)
             plt.close (fig =fig )
             
        return ax

    def plottensors (self,  resargs , phsargs , res_err=None, phs_err =None, 
                     sharex= False,  show_grid =True, errorbar_kw = None, 
                     leg= None, **kws): 
        """ 
        
        Parameters 
        ----------
        resargs : list of array-like 
            List of resistivity values to plot. Each resistivity value is an 
            array-like composed of collected data at each tensor component. 
            
        phsargs: list of array-like 
            List of phases values  in degree to plot. Each phase value is an 
            array-like composed of collected data at each tensor component.
        res_err : list of array-like 
            List of resistivity errors values to plot. Each resistivity error is 
            an  array-like at each tensor component. 
            
        phs_err: list of array-like 
            List of phase errors values in degrees to plot. Each phase error is 
            an  array-like at each tensor component.
            
        kind: str 
             Type of log plot. Can be ``loglog``, ``semilogx``, semilogy`` or 
             merely a base plot function. Default is ``semilogy``.
             
        show_grid: bool, 
            Will show a grid once triggered i.e set to ``True``. Note that some 
            argument are fixed and can not match the purpose of the users. 
            Edit the part of gridding plot to match your custom plot. Default 
            is ``False``.
        sharex, shareyAxes, optional
            Share the x or y axis with sharex and/or sharey. The axis will
            have the same limits, ticks, and scale as the axis of the shared axes. 
            Default is ``False``. 
            
        show_grid: bool, 
            Will show a grid once triggered i.e set to ``True``. Note that some 
            argument are fixed and can not match the purpose of the users. 
            Edit the part of gridding plot to match your custom plot. Default 
            is ``False``.
        errorbar_kw: dict 
            Error bar additional keyword arguments. Refer to
            :func:`pycsamt.utils.plot_utils.plot_errobar` function. 
 
        leg: list 
            legend labels put into a list. IOt must match the number of given 
            resistivity values. 
            
        kws : dict 
            Matplotlib logs functions additional keywords arguments. 
            
        plt.xlabel ('$Log_{10} frequency [H_z]$') ; plt.ylabel('$ Resistivity [ \Omega.m]$')
        
        Returns 
        ------- 
        ax: Matplotlib.pyplot <AxesSubplot> 
        
        Examples 
        ---------
        >>> from pycsamt.view.plot import Plot1d 
        >>> from pycsamt.processing import flma # filter to correct resistivity 
        >>> edipath = 'data/3edis'
        >>> ediObjs = get_ediObjs(edipath) 
        >>> res2d = make2d (ediObjs, 'resyx')
        >>> phs2d = make2d (ediObjs, 'phaseyx')
        >>> freqs = get_full_frequency(ediObjs)
        >>> edipath = 'data/3edis'
        >>> # -----XX-----
        >>> rhoxx = interpolate2d (make2d (ediObjs, 'resxx')) 
        >>> phsxx = interpolate2d (make2d (ediObjs, 'phasexx'))
        >>> rhoxx_err = interpolate2d (make2d (ediObjs, 'resxx_err')) 
        >>> phsxx_err = interpolate2d (make2d (ediObjs, 'phasexx_err'))
        >>> # -----XY-----
        >>> rhoxy = interpolate2d (make2d (ediObjs, 'resxy')) 
        >>> phsxy = interpolate2d (make2d (ediObjs, 'phasexy'))
        >>> rhoxy_err = interpolate2d (make2d (ediObjs, 'resxy_err')) 
        >>> phsxy_err = interpolate2d (make2d (ediObjs, 'phasexy_err'))
        >>> # -----YX-----
        >>> rhoyx = interpolate2d (make2d (ediObjs, 'resyx')) 
        >>> phsyx = interpolate2d (make2d (ediObjs, 'phaseyx'))
        >>> rhoyx_err = interpolate2d (make2d (ediObjs, 'resyx_err')) 
        >>> phsyx_err = interpolate2d (make2d (ediObjs, 'phaseyx_err'))
        >>> # -----YY-----
        >>> rhoyy = interpolate2d (make2d (ediObjs, 'resyy'))
        >>> phsyy = interpolate2d (make2d (ediObjs, 'phaseyy'))
        >>> rhoyy_err = interpolate2d (make2d (ediObjs, 'resyy_err')) 
        >>> phsyy_err = interpolate2d (make2d (ediObjs, 'phaseyy_err'))
        >>> # frame phase to between 0 and 90 degree 
        >>> def frame_phs (arr) : 
                return np.abs(arr)%90 
        >>> phsxy= frame_phs(phsxy);
        >>> phsyx= frame_phs(phsyx);
        >>> phsyy= frame_phs(phsyy); 
        >>> phsxx= frame_phs(phsxx) 
        >>> leg =['XX', 'XY', 'YX', 'YY']
        >>> res_err= [rhoxy_err[:, 3] , rhoyx_err[:, 3] ] 
        >>> phs_err = [phsxy_err[:, 3] , phsyx_err[:, 3] ] 
        >>> res_args = [f, rhoxx [:, 3], 'tab:pink', f, rhoxy [:, 3],
                        'ob-', f, rhoyx [:, 3], 'or-', f, rhoyy [:, 3], 'tab:cyan'] 
        >>> phs_args = [ f, phsxx [:, 3], 'tab:pink', f, phsxy [:, 3], 
                        'ob-', f, phsyx [:, 3], 'or-', f, phsyy [:, 3], 'tab:cyan'] 
        >>> p1d = Plot1d (fig_size =(5, 4), markeredgecolor ='k', markersize = 2 )
        >>> p1d.xlabel ='$Log_{10}Frequency [HZ]$'; p1d.ylabel ='$Rhoa (\Omega.m)$'
        >>> p1d.plottensors(resargs = res_args, phsargs =phs_args, leg=leg, 
                            res_err =res_err, phs_err = phs_err)

        References 
        -----------
        https://matplotlib.org/stable/gallery/color/named_colors.html
        
        """
        errorbar_kw  = errorbar_kw or {}
        res = func._assert_all_types(resargs, list, tuple)
        phs = func._assert_all_types(phsargs, list, tuple)
        
        fig, ax = plt.subplots(2, sharex = sharex , figsize = self.fig_size, 
                               num = self.fig_num, dpi = self.fig_dpi , 
                    )
        
        ax[0].loglog( *res, 
                     markersize = self.ms , 
                     linewidth = self.lw, 
                     **kws)
        ax[1].semilogx( *phs,  
                       markersize = self.ms, 
                       linewidth = self.lw, 
                       **kws)

        if res_err is not None  : 
            res_err = func._assert_all_types(res_err, list, tuple, np.ndarray)
            jj=0
            for ii in  range(len(res_err)) : 
                mplotus.plot_errorbar(ax[0], res[jj], res[jj+1],
                                              y_error = res_err[ii], 
                                              **errorbar_kw 
                                              )
         
                jj +=3  # to jup to next value if 
                
        if phs_err is not None: 
            phs_err = func._assert_all_types(phs_err, list, tuple, np.ndarray)
            jj=0
            for ii in range(len(phs_err)) : 
                mplotus.plot_errorbar(ax[1], phs[jj], phs[jj+1],
                                              y_error=phs_err[ii], 
                                              **errorbar_kw 
                                              )
                jj +=3
                
        ax[0].set_xlabel ( self.xlabel or '$Log_{10} frequency [H_z]$' 
                          if not sharex else None )
        ax[0].set_ylabel (self.ylabel  or '$ Rhoa [ \Omega.m]$') 
        ax[1].set_xlabel ( self.xlabel or '$Log_{10} frequency [H_z]$')
        ax[1].set_ylabel (self.ylabel  or '$ Phase[ deg^o]$') 
        
        ax[0].legend( [] if leg is None else leg , loc ='best' , ncol = 4 )
        
        if show_grid :
            ax[0].grid (visible =True , alpha =0.8, which ='both', color ='gray')
            ax[1].grid (visible =True , alpha =0.8, which ='both', color ='gray')
 
        plt.tight_layout()   
        
        if self.savefig is not None: 
             plt.savefig(self.savefig , dpi = self.fig_dpi)
        
        return ax
    
    def plotlogs (self, *args,  kind ='semilogy', leg=None, show_grid=False, 
                  **kws): 
        """ loglog or semilogs plots. 
        
        Template to quick visualize the 1D-log plots. 
        
        Parameters 
        ----------
        args : list 
            Matplotlib logs funtions plot arguments 
            
        kind: str 
             Type of log plot. Can be ``loglog``, ``semilogx``, semilogy`` or 
             merely a base plot function. Default is ``semilogy``.
             
        show_grid: bool, 
            Will show a grid once triggered i.e set to ``True``. Note that some 
            argument are fixed and can not match the purpose of the users. 
            Edit the part of gridding plot to match your custom plot. Default 
            is ``False``.
            
        leg: list 
            legend labels put into a list.
        kws : dict 
            Matplotlib logs functions additional keywords arguments. 
        
        Examples 
        ---------
        >>> from pycsamt.view.plot import Plot1d 
        >>> from pycsamt.processing import flma # filter to correct resistivity 
        >>> edipath = 'data/3edis'
        >>> ediObjs = get_ediObjs(edipath) 
        >>> res2d = make2d (ediObjs, 'resyx')
        >>> phs2d = make2d (ediObjs, 'phaseyx')
        >>> freqs = get_full_frequency(ediObjs)
        >>> resf =flma(res2d =res2d, phs2d =phs2d, freqs= freqs, window_size=5)
        >>> x= np.arange (0 , res2d.shape[1]*50 , 50) # step between station is 50 m 
        >>> leg= ['uncorrected data', 'filtered data'] # legend labels .
        >>> figtitle = 'Data filtered with fixed length moving average' 
        >>> p1d = Plot1d (fig_size =(10, 5), xlabel = 'Stations',
                          ylabel='$ Resistivity [ \Omega.m]$', fig_title= figtitle )
        >>> p1d.plotlogs (x, res2d[3, :], 'b--',x, resf[3, :], 'ok--',
                          kind ='semilogy', show_grid =True, leg=leg )
      
        """
        
        plt.figure ( figsize = self.fig_size, num = self.fig_num,
                    dpi = self.fig_dpi)
        
        pplot = plt.plot 
        if kind =='loglog' : 
            pplot = plt.loglog 
        elif kind =='semilogx': 
            pplot = plt.semilogx 
        elif kind =='semilogy': 
            pplot = plt.semilogy 
            
        pplot(*args, 
              markersize = self.ms , 
              linewidth = self.lw,
              # if markerface and edge are uncommented, all maker 
              # should hold the same view.
              # markerfacecolor = self.markerfacecolor ,
              # markeredgecolor=self.markeredgecolor, 
                **kws
                )
        plt.xlabel(self.xlabel)
        plt.ylabel(self.ylabel )
        plt.legend ([] if leg is None else leg, ncol = len(args)//3)
        plt.title (self.fig_title)
        if show_grid :
            plt.grid (visible =True , alpha =self.alpha, which ='both', color ='gray')
            
        plt.tight_layout()
        
        
    def plot_topo_sep_azim(self,fn = None , profile_fn=None  , 
                           savefig =None ,  **kwargs):
        """
        Method to plot topographic , stations separation and azimuth profiles .
        User can add station_names and and set _it to let the program to plot 
        on the corresponding  figure. He can alse force the program to plot
        its dipole length otherwise the program will compute it automatically.
         If "set_station_name" is False , No Name of station will be visible.
        User has the possibility to plot one by one figure or all by using a 
        "*"  symbol or 123. To polt one figure , it may use keyword argument 
        "plot" following the king ['topo', azimuth', 'sep'] or integer 1|2|3.
        Method is flexible .User can customize the plot , marker and 
        line as he wants by putting on list the matplotlib labels properties . 
        The program uses the label properties on order to set configuration lines 
        and other properties .Topography plot correspond to index 0 ,
        stations-separation to index 1 and azimuth to index 2.To plot individually 
       , User doesnt need to put properties on list. Programm will
        recognize and set the poperties provided according the figure He wants. 
        
        Parameters
        ------------
            * fn : str
                full path to [EDI|J|AVG] file. 
                
            * profile_fn : str 
                path to file  may Zonge Engineering *.stn  file 
                
            * plot : str     
                type of plot ,  default is '*' mean of three profile.
                
            * Station_Names: list   
                    list of station names , User could provide. Default is None 
                    compute automatically 
                    
            * set_station_names : bool  
                    display the station name on figure axis . Default is False.
                    
            * elevation : (ndarray,1)
                    Array_like of elevation
                    
            * station_pk : array_like, 
                    array_like station dipole center value.
                    
            * savefig : str 
                    path  to save figure. 
                    
        Returns
        ---------
            obj, 
                plot_obj azim- topo and station separation. 
        
        ======================  ===============================================
        Keywords                Description        
        ======================  ===============================================
        lw                      line width . *default* is 1.5
        ls                      line style of lines,[ '-' | '.' | ':' ] 
                                *default* is "['-', ':', '-.']" for 3 profiles.
        marker                  marker of stations  *default* is 'o'
        ms                      size of marker in points. *default* is 6
        color                   color of line .*Default* is 'k'
        alpha                   Marker transparence .*Defaut* is .2 
        markerfacecolor         facecolor or markers .*Default*  is 'k'
        markeredgecolor         eadgecolor of markers .
                                default is "['w','r''gray']"
        xtick_label_rotation    xtick rotation angle . default* is 45.
        ytick_label_rotation    ytick rotatoion angle .default * is 45 
        xtick_labelsize         xtick label size .defalut* is 12.
        ytick_labelsize         ytick label size .defalut* is 12.
        ======================  ===============================================
        
        :Example: 
            
            >>> import os
            >>> file_stn='K6.stn'
            >>> path =  os.path.join(os.environ["pyCSAMT"], 
            ...                          'pycsamt','data', file_stn)
            ...  plot_1d_obj= Plot1d()
            .... plot_1d_obj.plot_topo_sep_azim(profile_fn= path , plot='*', 
                                                set_station_names=True,
                                               dipole_length_curve=False)
        """
        self._logging.info(
            'Plotting Topography-Station_separations and Azimuth.')
        
        mpl.rcParams['font.size']=10
        
        #----------- set arguments ---------------------
        plot_type =kwargs.pop('plot', '*')
        add_stnnames =kwargs.pop('station_names',None)

        set_stnnames=kwargs.pop('set_station_names', False)
        
        elevation_array =kwargs.pop('elevation', None) 
        station_pk =kwargs.pop('station_pk', None )
        azimuth_array =kwargs.pop('azimuth', None)
        stn_separation_array =('stn_separation', None)

        
        DipoleLength =None 

        #----------set mpl properties ----------------- . 
        
        kws = mplotus.share_props_for_each_plot(
            number_of_plot=3,lw=kwargs.pop('lw', 1.5), 
                    ls=kwargs.pop('ls', ['-', ':', '-.']), 
                    alpha=kwargs.pop('alpha', .2), 
                    color=kwargs.pop('c', 'k'),
                    marker= kwargs.pop('marker_style', 'o'),
                    markerfacecolor=kwargs.pop('markerfacecolor', 'k'),
                    markeredgecolor= kwargs.pop('markeredgecolor',
                                                ['w','r', 'gray']),
                    xtick_label_rotation=kwargs.pop('xtick_label_rotation',45.),
                    ytick_label_rotation=kwargs.pop('ytick_label_rotation',45.),
                    xtick_label_size= kwargs.pop('xtick_labelsize', 12), 
                    ytick_label_size= kwargs.pop('ytick_labelsize', 12)
                                          )
        lw ,ls, color, alpha =kws['lw'],kws ['ls'], kws['color'], kws['alpha']
        
        x_ticklabel_rotation, x_ticks_labelsize ,\
            y_ticks_labelsize = kws['xtick_label_rotation'], \
                                    kws['xtick_label_size'],\
                                    kws['ytick_label_size']
            
        marker_style, markerfacecolor, markeredgecolor =\
            kws['marker'], kws['markerfacecolor'], kws['markeredgecolor']

        #--------------- end of setting mpl properties -----------------------
        

        #  define profile object :
        if profile_fn is not None : 
            profile_obj = Profile (profile_fn =profile_fn)
            
            easting , northing ,elevation_array,\
                station_pk  = profile_obj.east, profile_obj.north, \
                              profile_obj.elev, profile_obj.stn_position

            azimuth_array = profile_obj.azimuth
            stn_separation_array =profile_obj.stn_interval
            DipoleLength = profile_obj.dipole_length
            profile_obj.get_profile_angle(easting = easting,
                                          northing =northing )
            lineazimuth =profile_obj.profile_angle 
        # elif profile_fn is None : 
        #     DipoleLength = stn_separation_array.mean()
        elif profile_fn is None and fn is not None : 
            csamt_obj =CSAMT(data_fn =fn )
            easting, northing, elevation_array , station_pk = csamt_obj.east , \
                csamt_obj.north , csamt_obj.elev,csamt_obj.station_distance
            
            azimuth_array = func.compute_azimuth(
                easting = easting, northing =northing, extrapolate=True) 
            stn_separation_array = csamt_obj.station_separation
            DipoleLength = csamt_obj.dipolelength
            lineazimuth, pmess= geoD.geostrike.compute_profile_angle(
                easting= easting , 
                northing =northing)
            print('---> ' + pmess)        # display profile angle message
            
                
        elif profile_fn is None and fn is None  :
            raise CSex.pyCSAMTError_plot(
                'None path is found. Please spceify you data path.')
            
        
        #-----------------------------Build fonction --------------------------
        def plot_topography(axes=None, set_xlabel=True ): 
            """
            plot topography 
            """
            topo_profile = axes.plot(station_pk, elevation_array, 
                                     lw=lw[0], marker=marker_style[0],
                                     markersize=self.ms,
                                     markerfacecolor=markerfacecolor[0], 
                                     markeredgecolor = markeredgecolor[0],
                                     c=color[0])
            
            axes.legend( topo_profile, ['$profile$'])
            # topo_scatter = axes.scatter(station_pk, elevation_array, lw=lw,
            #                                          ls=self.ls[0] ,  
            #                                          alpha=alpha,
            #                                          marker=marker_style[0], 
            #                                          s=self.ms,
            #                                          )
            if set_xlabel : 
                axes.set_xlabel ('Stations', color='black', fontsize = 12)
                # axes.set_xlabel ('Stations',
                # color='k', fontdict={'size': 12, 'weight': 'bold'})

                
            if plot_type != '*': axes.set_title('Topographic profile')
            axes.set_ylabel ('Elevation(m)', color='black', fontsize = 12)


       
        def plot_stn_separation(axes =None, set_xlabel=True): 
            """
            plot station _separations  x_station / h- separtion(m)
             dipole_length_curve show the line of 
             overage value on de site between all stations  
            """

            stn_sep_profile = axes.plot(station_pk, stn_separation_array,
                                        markersize=self.ms, 
                                        lw= lw[1],marker=marker_style[1], 
                                        markerfacecolor=markerfacecolor[1], 
                                        markeredgecolor = markeredgecolor[1],
                                        ls=ls[1], 
                                        c=color[1])
            
            if set_xlabel : 
                axes.set_xlabel ('Stations',  color='black', fontsize = 12)
                
            if plot_type != '*': axes.set_title('Station Separation')
            
            axes.set_ylabel ('Separation(m)', color='black', fontsize = 12)
            
            # axes.grid(color='k', ls='--', lw =0.125, alpha=0.1, which='minor')
            # customize specific grid 

            axes.text(station_pk.min()+ DipoleLength , stn_separation_array.min() , 
                                          '$DipoleLength~:{0}m$'.format(DipoleLength), 
                                          verticalalignment='bottom', 
                                          c='k', fontsize =12 )
            
            # if dipole_length_curve is True : axes.legend( (stn_sep_profile,
            # dip_length_pf  ), labels=[ '$Separation$','$Average Bar$'])
            # else : axes.legend(stn_sep_profile,['$separation$'])
            axes.legend(stn_sep_profile,['$separation$'])
            

        def plot_azimuth(axes =None, set_xlabel=True): 
            """
            plot_azimuth  
            
            :param set_xlabel: show the labels of on x
            :type set_xlabel: bool
            
            """
            # plot the azimuth

            azim_profile = axes.plot(station_pk, azimuth_array , 
                                     lw=lw[2],marker= marker_style[2],
                                    markersize=self.ms ,
                                    markeredgecolor=markeredgecolor[2], 
                                    markerfacecolor=markerfacecolor[2], 
                                     ls=ls[2] ,c=color[2],
                                     )
            
            
            if plot_type != '*':axes.set_title('Azimuth profile ')                                        
            if set_xlabel : 
                axes.set_xlabel ('stations', color='black',
                                 fontsize = x_ticks_labelsize[2])
            axes.minorticks_on()    
            axes.grid(color='k', ls=':', lw =0.25, alpha=0.7,
                      which ='major') # customize specific grid
            #add text on azimuth profile . 
            axes.text(2* station_pk.min(),
                      azimuth_array.min() , 
                      '$Line Azimuth~:{0}°$'.format(
                          np.around(lineazimuth,2)),
                      #np.around(azimuth_array.mean(),2)), 
                      verticalalignment='bottom', 
                      c='k',
                      fontsize =12 )
            
            axes.set_ylabel ('azimuth (°)', color='black', 
                             fontsize =  x_ticks_labelsize[2])
            axes.legend( azim_profile, ['$azimuth$'])
            
        # --- works on exception to let function less expensive ---- 
        
        if isinstance (plot_type, str): 
            if suit.stn_separation[0].find(plot_type.lower()) < 0 :
                if  suit.topography[0].find(plot_type.lower())< 0 : 
                    if suit.azimuth[0].find(plot_type.lower()) < 0 : 
                        if plot_type.lower() not in ['1','2','3','123','*']:
                            raise CSex.pyCSAMTError_plot(
                        'Argument provided for  plot_type is not acceptable.'
                            ' Please use  1,2,3 or the first '
                             'two letter of "Topography, Separation|Station '
                                 'or Azimuth " or "*|123" to plot the 3 profiles.')
            
        elif isinstance(plot_type, int):
            if plot_type not in [1,2,3,123]: 
                raise CSex.pyCSAMTError_plot(
                    'Only 1|2|3 is need to plot. Try again !')
            plot_type=str(plot_type)

        # ----plot only single figures Topo/sep/azim  ------------------------

        if plot_type not in ['*','123']:
            
            mpl.rcParams['figure.figsize']= (12,5)
            fig, ax =  plt.subplots() 
            # ax= fig.add_axes([0,0,1,1])
            if  plot_type == "1" or (re.match(r'^to+',
                                              plot_type.lower()) is not None) :
                plot_topography(axes=ax)
            elif plot_type == "3" or (re.match(r'^az+', 
                                               plot_type.lower()) is not None) :
                plot_azimuth(axes=ax)

            elif plot_type == "2" or (re.match(r'^se+',
                                               plot_type.lower()) is not None)\
                or (re.match(r'^st+', plot_type.lower())is not None):
                # if dipole_length_curve is not True : 
                    #plot_stn_separation(axes=ax, dipole_length_curve=False)
                 plot_stn_separation(axes=ax)

            if set_stnnames ==True : 
                if add_stnnames  is not None : 
                    if len (add_stnnames) != station_pk.size :
                        raise CSex.pyCSAMTError_station(
                            'Stations names provided must'
                              ' be the same length as survey points.'
                              ' Number of survey points are'
                               ' : <{0}>'.format(station_pk.size))
                    xticklabel =add_stnnames 
    
                elif add_stnnames  is  None : 
                    xticklabel = ['S{0:02}'.format(ss)
                                  for ss in range(station_pk.size)]
                    
                ax.set_xlim ([station_pk.min(), station_pk.max()] )   
                ax.set_xticks(ticks= station_pk.tolist(),minor=False )
                ax.set_xticklabels(xticklabel , rotation=x_ticklabel_rotation[0]) 
                # rotate specific station 
            else : 
                ax.set_xlim ([station_pk.min(), station_pk.max()] ) 
            
            #ax.grid(color='k', ls='--', lw =0.125, alpha=0.1) 
            # customize specific grid 
                # custumize xand yticks labelsize
            ax.xaxis.set_tick_params(labelsize=x_ticks_labelsize[0])
            ax.yaxis.set_tick_params(labelsize=x_ticks_labelsize[0])
            
            if savefig is not None : 
                plt.savefig(savefig, dpi=self.fig_dpi)
                
            plt.show()
        #----------------------------- plot all 3  using '*' -----------------   
        elif plot_type =='*' or plot_type=='123':

            mpl.rcParams['figure.figsize']= self.fig_size
            
            fig=plt.figure()
            axe1 =fig.add_subplot(311)
            # fig , ax =plt.subplots(3, figsize =self.fig_size, sharex =True )
            plot_topography(axes=axe1, set_xlabel =False )
            
            axe2 = fig.add_subplot(312, sharex=axe1)
            plot_stn_separation(axes=axe2, set_xlabel=False)
            
            axe3 = fig.add_subplot(313, sharex=axe1)
            # if dipole_length_curve is False :
                #plot_stn_separation(axes=axe3, dipole_length_curve=False)
            # else :
            plot_azimuth(axes=axe3, set_xlabel =False)
            
            if set_stnnames : 
                if add_stnnames  is not None : 
                    if len (add_stnnames) != station_pk.size : 
                        raise CSex.pyCSAMTError_station(
                            'Stations names provided must'
                                ' be the same length as survey points.'
                                ' Number of survey points are'
                                  ' : <{0}>'.format(station_pk.size))
                    xticklabel =add_stnnames 
                elif add_stnnames  is  None : 
                    xticklabel = ['S{0:02}'.format(ss) 
                                  for ss in range(station_pk.size)]
                
            for ii, axx in enumerate([axe1, axe2, axe3]):
                
                if set_stnnames:axx .set_xticks(
                        ticks= station_pk.tolist(), minor=False ) 
                axx .set_xlim([station_pk.min(), station_pk.max()])
                axx .grid(color='k', ls='--', lw =0.25, alpha=0.3)
                if set_stnnames: axx .set_xticklabels(xticklabel , 
                                        rotation=x_ticklabel_rotation[ii]) 
                # ,va='center')
                
            axe3.set_xlabel('stations', )
            fig.tight_layout()
            fig.suptitle(
                'Topography-Stations separation-Station azimuth profiles',
                     fontsize= 4 * self.font_size, 
                     verticalalignment='center', 
                     style ='italic',
                     bbox =dict(boxstyle='round',facecolor ='moccasin'))
        
            plt.show()
        
        
    def plot_static_correction(self, data_fn, profile_fn =None , 
                               dipole_length =50., 
                               frequency_id  =1 , ADD_FILTER ='tma' ,
                               **kwargs): 
        """
        Plot corrected apparent resistivities at different stations to solve 
        the problem of static shift by adding either Trimimg Moving average 
        (TMA)filter or fixed-length-moving-average(FLMA) filter or Adaptative 
        moving-average (AMA). Actually FLMA and TMA filter are available and 
        default filTer is TMA. To plot all filter into one figure add the joker 
        `*` to arguments `ADD_FILER`.
        
            
        :param data_fn: full path to file , can be [AVG|EDI|J] files
        :type data_fn: str 
        
        :param profile_fn: 
            pathLike  full path to Zonge Engeneering *.station file .
        :type profile_fn: str 
        
        :param dipole_length: 
            length of dipole in meters when user applied for FLMA 
        :type dipole_length: float, int
        
        Holding others informations 
            
        ===============  ===========  =========================================
        Params              Default         Description
        ===============  ===========  =========================================
        frequency_id        str,int         plot the filtered frequency,
                                            eg  frequency_id = 1023 means  
                                            plot  uncorrected rho and static rho 
                                            at that frequency . set on list to 
                                            plot multiple frequency [8,1101 ].
        ADD_FILTER          str             name of filter to apply . 
                                            TMA  Trimming moving average
                                            AMA  Adaptative moving average 
                                            FLMA  Fixed Length moving average
        ===============  ===========  =========================================
        
        .. note:: Profile file (*.stn) is compolsory when provide raw Zonge AVG 
                otherwise no need to profile for EDI or J file. In addition 
                when FLMA filter is used m, provided `dipole_length` and 
                `number_of_points` for window width are necessaries.
   
        :Example: 
            
            >>> from from viewer.plot import Plot1d 
            >>> path =  os.path.join(os.environ["pyCSAMT"], 
            ...              'pycsamt','data', file_1)
            >>> plot_1d_obj =Plot1d()
            ... plot_1d_obj.plot_static_correction(data_fn =path , 
            ...                                   profile_fn= os.path.join(
            ...                                    os.path.dirname(path),
            ...                                            'K1.stn'), 
            ...                                   frequency_id =1023)
        """
        
        FILTER =['tma', 'flma','ama' ,'*']
        
            
        #-------------------------------------------Fonction properties -------
        FILTERpoints =kwargs.pop('number_of_points',5)
        number_of_skin_depth =kwargs.pop('number_of_skin_depth', 1.)
        
        savefig =kwargs.pop('savefig', None)
        orient =kwargs.pop('orientation', 'landscape')
        addstn= kwargs.pop('set_station_names', False)
        
        fill_between =kwargs.pop('fill_between', False)
        fillbc=kwargs.pop('fill_between_color', 'thistle')
        
        tma_color = kwargs.pop('tma_color', 'blue')
        flma_color =kwargs.pop('flma_color', 'aqua')
        ama_color =kwargs.pop('ama_color', 'lime')
        

        #control filter provided 
      
        if ADD_FILTER not in FILTER : 
            try :
                ADD_FILTER = int(ADD_FILTER)
            except : 
                if isinstance( ADD_FILTER, str):
                    if ADD_FILTER.lower() not in FILTER:       
                        warnings.warn('Currently , filters availabale are '
                                      'trimming moving-average (TMA).'
                                      ' and  fixed-length-moving-average(FLMA)'
                                      'Please Try to use the availabe filters.')
                        raise CSex.pyCSAMTError_plot(
                            'Input filter <{0}> is not available.'
                              ' Please use `tma` or `flma` filters!')
                        
                ADD_FILTER = ADD_FILTER.lower()  
            else : 
                if ADD_FILTER ==1 or ADD_FILTER ==0: ADD_FILTER='tma'
                elif ADD_FILTER ==2 : ADD_FILTER='flma'
                elif ADD_FILTER ==3 : ADD_FILTER ='ama'
                elif ADD_FILTER ==123: ADD_FILTER ='*'
                else :
                    ADD_FILTER ='tma' # block to default filter

        #if ADD_FILTER =='ama': ADD_FILTER='tma'
            
        #---> corrected data TMA , FLMA and AMA filters   ---------
        # if ADD_FILTER is not None : 
        corr_obj = Processing(data_fn=data_fn , profile_fn=profile_fn,
                        reference_freq=frequency_id )
        csamt_freq_obj=corr_obj.frequency
        csamt_res_obj= corr_obj.res_app_obj
        csamt_stn_num_obj = corr_obj.station_distance 
        stnnames =corr_obj.site_id 
        if frequency_id is None :
            frequency_id =corr_obj.referencefreq

        if ADD_FILTER =='*':
            
            # create correct TMA obj 
            tma_cor_obj= corr_obj.TMA(number_of_TMA_points=FILTERpoints)
            # create correct FLMA obj 
            flma_cor_obj = corr_obj.FLMA(dipole_length =dipole_length , 
                                      number_of_dipole =FILTERpoints)
            # create correct TMA obj 
            ama_cor_obj = corr_obj.AMA(dipole_length =dipole_length , 
                                      number_of_skin_depth =number_of_skin_depth)
            
            ADD_FILTER ='*'
            messf ='TMA & FLMA & AMA'
            print('---> Filters {} are done !'.format(messf) )
        else :
            if ADD_FILTER.lower()=='tma':
                res_cor_obj = corr_obj.TMA(number_of_TMA_points=FILTERpoints)
                
                messf = 'TMA'
                print('---> Filter {} is done !'.format(messf) )
                
            elif ADD_FILTER.lower()=='flma':
                res_cor_obj = corr_obj.FLMA(dipole_length =dipole_length , 
                                      number_of_dipole =FILTERpoints)

                
                messf = 'FLMA'
                print('---> Filter {} is done !'.format(messf) )
                
            elif ADD_FILTER.lower()=='ama':
                res_cor_obj = corr_obj.AMA(dipole_length =dipole_length , 
                                      number_of_skin_depth =number_of_skin_depth)
                messf = 'AMA'
                print('---> Filter {} is done !'.format(messf) )
        # call function to fixe freauency for plotting 
        freqID = mplotus.get_frequency_id(freq_array=csamt_freq_obj,
                                          frequency_id=frequency_id)

        for referfreq in freqID : 
            #-----------------------Build figure -------------------------
            
            mpl.rcParams['figure.figsize']=[12,6]
            mpl.rcParams['font.family'] ='sans-serif'
            fig =plt.figure()
            axe =fig.add_subplot(111)
            #-------------------------------------------Build Plot --------
            #call reference frequency fonction 
            # ---> rho uncorrected data 
            RES_UNCOR = Zcc.get_data_from_reference_frequency(
                                                array_loc=csamt_res_obj, 
                                                 freq_array=csamt_freq_obj, 
                                                  reffreq_value=referfreq )
            
            mark,  = axe.semilogy (csamt_stn_num_obj,RES_UNCOR ,
                                   c= 'white', 
                                   marker =self.marker,
                                  markersize = self.ms*2*self.fs ,
                                  markeredgecolor= self.markeredgecolor 
                                  )
            
            resPlots,= axe.semilogy(csamt_stn_num_obj,
                                    RES_UNCOR ,
                                    lw =self.lw*self.fs,
                                    c='k',
                                ls =self.ls,
                                label='$Uncorrected apparent resistivity$'
                                )
            
   
            
            if ADD_FILTER =='*': 
                
                TMA_RES_COR = Zcc.get_data_from_reference_frequency(
                                            array_loc=tma_cor_obj,
                                            freq_array=csamt_freq_obj, 
                                            reffreq_value=referfreq )
               
                FLMA_RES_COR = Zcc.get_data_from_reference_frequency(
                                           array_loc=flma_cor_obj,
                                           freq_array=csamt_freq_obj, 
                                           reffreq_value=referfreq )
                
                AMA_RES_COR = Zcc.get_data_from_reference_frequency(
                                           array_loc=ama_cor_obj,
                                           freq_array=csamt_freq_obj, 
                                           reffreq_value=referfreq )
                
                marki,  = axe.semilogy (csamt_stn_num_obj,TMA_RES_COR,
                                        c= 'white', 
                                        marker ='D',
                                        markersize = self.ms*2*self.fs ,
                                        markeredgecolor= tma_color ,
                                        alpha =0.8)
                
                marki,  = axe.semilogy (csamt_stn_num_obj,FLMA_RES_COR,
                        c= 'white', marker ='o',
                        markersize = self.ms*2*self.fs ,
                        markeredgecolor= flma_color ,
                        alpha =0.8)
                
                marki,  = axe.semilogy (csamt_stn_num_obj,AMA_RES_COR,
                        c= 'white', marker ='o',
                        markersize = self.ms*2*self.fs ,
                        markeredgecolor= ama_color ,
                        alpha =0.8)
                
                tma_corPlots, = axe.semilogy(csamt_stn_num_obj,TMA_RES_COR ,
                                             lw =self.lw*self.fs, 
                                             c=tma_color,
                                            ls =self.ls, 
                                            label = '$Corrected\resistivity :'\
                                            ' TMA filter at {0} points$'.
                                            format( FILTERpoints))
                    
                flma_corPlots, = axe.semilogy(csamt_stn_num_obj,FLMA_RES_COR ,
                             lw =self.lw*self.fs, 
                            c=flma_color,
                            ls =self.ls, 
                            label = '$Corrected\resistivity :'\
                            ' FLMA filter at {0} dipole$'.
                            format( FILTERpoints))
                    
                ama_corPlots, = axe.semilogy(csamt_stn_num_obj,AMA_RES_COR ,
                         lw =self.lw*self.fs, 
                        c=ama_color,
                        ls =self.ls, 
                        label = '$Corrected\resistivity :'\
                        ' AMA filter at {0} skin depth $'.
                        format( int(number_of_skin_depth)))
                    
                    
                axe.legend( [resPlots, tma_corPlots,flma_corPlots,ama_corPlots ], 
                           ['$Uncorrected\ app. Rho$',
                            '$Rhofiltered: TMA at\ {0}\ points$'.
                                                   format( int(FILTERpoints)), 
                            '$Rhofiltered: FLMA at\ {0}\ dipoles$'.
                                                   format( int(FILTERpoints)), 
                            '$Rhofiltered: AMA at\ {0}\ skin depths$'.
                                        format( int(number_of_skin_depth))
                                                   ])# 
            else:
                if ADD_FILTER=='tma': color =tma_color
                elif ADD_FILTER=='flma': color =flma_color
                elif ADD_FILTER=='ama': color =ama_color
            
                RES_COR = Zcc.get_data_from_reference_frequency(
                                                array_loc=res_cor_obj,
                                                 freq_array=csamt_freq_obj, 
                                                      reffreq_value=referfreq )
            
           
                marki,  = axe.semilogy (csamt_stn_num_obj,RES_COR , 
                                        c= 'white',
                                        marker ='o',
                                     markersize = self.ms*2*self.fs ,
                                     markeredgecolor= color , 
                                     alpha =0.8)
                
                if ADD_FILTER=='ama':
                    ffmt, sfmt= '{0}'.format(int(number_of_skin_depth)),\
                    'skin depth(s).' 
                elif ADD_FILTER=='flma':
                    ffmt,sfmt= '{0}'.format(int(FILTERpoints)),'dipoles'
                else :
                    ffmt,sfmt= '{0}'.format(int(FILTERpoints)),'points'
                    
                corPlots, = axe.semilogy(csamt_stn_num_obj,RES_COR ,
                                         lw =self.lw*self.fs,
                                         c=color ,
                                ls =self.ls, 
                                label = '$Corrected\resistivity :'\
                                ' {0} filter at {1} {2}$'.format( messf, ffmt, sfmt))
            
                # ['$ApparentResistivity$', '$TMA at {0}$'.format(TMApoints)])
                axe.legend( [resPlots, corPlots ], ['$Uncorrected\ app. Rho$',
                                               '$Rhofiltered: {0}\ at\ {1}\ {2}$'.\
                                                   format(messf,ffmt, sfmt) ])
            axe.minorticks_on()
            
            axe.grid(color='k', ls=':', lw =0.25,
                     alpha=0.8, which ='minor') # customize specific grid
       
            axe.set_xlabel('$Station\ distance(m)$' ,
                           fontdict={'color': 'k',
                                    'size': self.x_minorticks*14,
                                    'weight':self.fontweight})
            
            axe.set_ylabel ('$Resistivity(Ω.m)$', color='black',
                            fontsize =  self.y_minorticks*14, 
                            fontweight=self.fontweight)
            
            
            
            axe.text(2* csamt_stn_num_obj.min(), RES_UNCOR.min() ,
                     '$ Filtered\ ref.frequency: <{0}>Hz$'.format(referfreq), 
                      verticalalignment='bottom',
                      fontweight =self.fontweight,
                      c='k', 
                      fontsize =12 )
            
            if ADD_FILTER !='*':
                if fill_between : axe.fill_between(csamt_stn_num_obj, 
                                                   RES_UNCOR, RES_COR,
                                                   facecolor='orange', 
                                                   color=fillbc, 
                                                   alpha =0.7)
            elif ADD_FILTER =='*':
                if fill_between : axe.fill_between(csamt_stn_num_obj,
                                                   TMA_RES_COR,
                                                   FLMA_RES_COR,
                                                   facecolor='orange', 
                                                   color=fillbc, 
                                                   alpha =0.7)
            if addstn ==True : 
                axe.set_xticks(ticks= csamt_stn_num_obj, minor=False ) 

                axe .grid(color='k',
                          ls='--',
                          lw =0.25, 
                          alpha=0.3)
                axe.set_xticklabels(stnnames , rotation=45)
                
          
            plt.tight_layout()
            fig.suptitle('Static correction plot ', 
                          style ='italic', 
                          bbox =dict(boxstyle='round',
                                     facecolor ='lightgrey'))
            if savefig is not None :
                plt.savefig(savefig,
                            dpi=self.fig_dpi,
                            orientation =orient)
        plt.show()
    @geoplot1d(reason = 'staticshift', mtmm='s', ctmm=(.8, 0.2, .9),
                      color_mode='color', ms_r=1.,lw_r=.7, 
                      )                   
    def plot_multiple_corrections(self,  data_fn, profile_fn =None , 
                               dipole_length =50., filter_type ='ama',
                               siteOrfreq='S00', refreq=None, 
                               **kws): 
        """
        Plot multiples static correction for a line using two options like:: 
            - plot rho|phase VS frequency `kind`= linear 
            - plot rho|phase VS station `kind`= site 
            default is ``line` 
        
        :param data_fn: full path to file , can be [AVG|EDI|J] files
        :type data_fn: str 
        
        :param siteOrfreq: 
            list of station or list of frequency number to visualize.
             Mx list tolorate is 4.  can be:: 
      
                    -`siteOrfreq`=['S00', 'S04', 'S08', 'S12']
                    -`siteOrfreq`=[1, 82, 1024, 8192]  in Hz
        :param kind: Type of plot. can be::
            - kind =1 or 'line'
            - kind =2 or 'site'
                    
        :param profile_fn: 
            pathLike  full path to Zonge Engeneering *.station file .If *.edi 
            file or jfile is given, dont need  to provide `profile_fn`
        :type profile_fn: str 
        
        :param dipole_length: length of dipole in meters
        :type dipole_length: float, int
        
        :Example:
            >>> from pycsamt.view.plot import Plot1d 
            >>> data='data/avg/K1.AVG'
            >>> viewfreq=[80, 1024,2000,  8192] # set `kind` to ``1``.
            >>> viewsite= ['S00', 's04', 's08', 'S12'] # set `kind` to ``2``
            >>> plotObj = Plot1d()
            >>> plotObj. plot_multiple_corrections(data_fn= data, kind=2, 
                                      filter_type ='ama',
                                       siteOrfreq=viewfreq
                                      )
        """
        nskinDepth = kws.pop('skin_depth', 1)
        npoints =kws.pop('npoints', 5)
        kindPlot =kws.pop('kind', 'line') # or site
        
        if siteOrfreq is None : 
            _logger.debug ( 'None sation found can not be plotted.'
                           'Should be considered the default value =`S00`.')
            stations =['S00']
            kindPlot =1
         
        # check the sation and duplicate until the max gridspec =4 

        if isinstance(siteOrfreq, (str, float, int)): 
            # duplicate the 
            siteOrfreq_lst =[siteOrfreq]
        elif isinstance(siteOrfreq, list): 
            if len(siteOrfreq) > 4: 
                siteOrfreq_lst =siteOrfreq[:4]
            else : siteOrfreq_lst= siteOrfreq
            
        try : 
            if os.path.isfile(data_fn): 
                line= os.path.basename(os.path.splitext(data_fn)[0])
        except: 
            line = os.path.basename(data_fn)
         # generate lines number multiply by number of stnOrFreq to visualize    
        lines =  [line  for j in range(len(siteOrfreq_lst))]
          
    
        def wf_(w_, lenfreq): 
            """ Build weight factor matrix"""
        
            w = np.zeros((lenfreq, len(w_)))
            for ii,  rowlines in enumerate(w): 
                w[ii, :] =w_
            return w

        dict_filter ={
            1:['ama', 'adaptative moving-average', 'torres verdin'],
            2: ['flma', 'fixed dipole moving-average', 'fixed dipole'],
            3: ['tma', 'trimming moving-average', 'trimming']
                      }
        try : 
            kindPlot = int(kindPlot)
        except: 
            if kindPlot.lower().find('line')>=0: 
                kindPlot =1 
            elif kindPlot.lower().find('site')>=0 or\
                kindPlot.lower().find('sta')>=0: 
                kindPlot =2 
    
        else: 
            if not kindPlot in [1, 2]: 
                kindPlot =1 

                
        #  manage filter 
        try : 
            filter_type = int(filter_type)
        except :
            f_flag = False 
            for key, values in dict_filter.items(): 
                for nval in values : 
                    if nval.find(filter_type.lower())>=0: 
                        filter_type  = int(key)
                        f_flag =True 
                        break
                if f_flag is True: 
                    break 

            if f_flag is False:
                filter_type =1 # set default fiter
       
        else: 
             if not filter_type [1, 2, 3]: 
                 filter_type =1 
        # get the first name of dict_filter value         
        ff = dict_filter[filter_type ][0] 

        # Read original data 
        sshiftObj = Processing(data_fn=data_fn , profile_fn=profile_fn,
                        reference_freq=refreq)
        
        measured_rho = sshiftObj.app_res_
        measured_phase = sshiftObj.phase_
        sites = sshiftObj.site_id 

        if filter_type == 1: 
            sshiftObj .AMA(number_of_skin_depth = nskinDepth)
        elif filter_type ==2: 
            sshiftObj.FLMA(number_of_dipole= npoints)
        elif filter_type ==3: 
            sshiftObj.TMA(num_of_TMA_points = npoints)
            
        corrected_phase = sshiftObj.phase_cor
        corrected_rho = sshiftObj.app_res_cor
        
        freq_array  =sshiftObj.frequency
        W_=wf_(sshiftObj._rj, len(freq_array)) # weigtfactormatrix 
        
        station_distance = sshiftObj.station_distance

        if freq_array [0] < freq_array[-1]: 

            freq_array  =freq_array [::-1]
            measured_rho = measured_rho[::-1]
            measured_phase = measured_phase[::-1]
            corrected_phase=corrected_phase[::-1]
            corrected_rho =corrected_rho[::-1]
   
        def retreive_single_siteOrFreqdata(stnOrFreq):
            """ 
            Get the plot data for one site or survey line after applying the 
            static shift correction. Retrieve essential data for `kind` plot.
            
            :param stnOrFreq: Can be a site name or frequency value 
            
            :returns: 
                - `stn`: station name or frequencu value 
                - `measured_rho`: value of observed rho data in ohm.m  
                - `measured_phase`: value of observed phase in degrees 
                - `corrected_rho`: rho corrected after staticsfit correction 
                - `corrected_phase`:phase corrected after apply 
                                    staticshift cor.
                - `wff`: weight factor 
            """
            
            if kindPlot ==1: # line mean plot rho|phase vs stations
                try: 
                    stnOrFreq=float(stnOrFreq)
                except: 
                    warnings.warn(f'The kind of plot is = `{kindPlot}`.'
                        f' Coud not found `{stnOrFreq}` in '
                                  f'frequency range ={freq_array}.')
                    self._logging.debug(
                        f'{stnOrFreq} not found in frequency range.')
                    raise CSex.pyCSAMTError_frequency(
                        'frequency must be a float '
                        'number not `{}`.'.format(type(stnOrFreq)))
                else: 
                    index_freq , freq_value = mplotus.getcloser_frequency(
                        freq_array, stnOrFreq)
                    
                measured_rho_r = measured_rho[index_freq, :]
                measured_phase_r = measured_phase[index_freq, :]
                corrected_phase_r = corrected_phase[index_freq, :]
                corrected_rho_r = corrected_rho[index_freq, :]
                # collected list of frequency
                stn='{} Hz'.format(round(freq_value, 1)) 
                
                wff = W_[index_freq, :]
                
            elif kindPlot ==2: #plot rho|phase vs frequency 
                
                stn = mplotus.get_stationid(stations=sites,
                                            station_id= stnOrFreq)[0]
                
                stn_index = sites.index(stn)
                # disce_site =station_distance [stn_index] 
                
                measured_rho_r = measured_rho[:, stn_index ]
                measured_phase_r = measured_phase[:, stn_index ]
                corrected_phase_r = corrected_phase[:, stn_index ]
                corrected_rho_r = corrected_rho[:, stn_index ]
                wff = W_[:, stn_index ]
                
            return (stn, (measured_rho_r, corrected_rho_r),\
                    (measured_phase_r, corrected_phase_r ), wff , wff) 
        
  
        # read the data and collect data 
        stations, appRHO, phase,\
            appRho_w, phase_w=[[]for i in range(5)]
       
        if kindPlot==1: # use station as axis plot and start by 0
        
            freq =[station_distance - station_distance.min()
                   for j in range(len(siteOrfreq_lst))] 
        else: freq =[freq_array for j in range(len(siteOrfreq_lst))] 
            
        ssinfos =(kindPlot, ff) 
         
        for ii, lin_ in enumerate(siteOrfreq_lst) : 
            stn_, rho_, phase_ , rho_err_ , phase_err_ =\
                retreive_single_siteOrFreqdata(lin_)
            stations.append(stn_)
            appRHO.append(rho_)
            phase.append(phase_)
            appRho_w.append(rho_err_)
            phase_w.append(phase_err_)
            
        return (ssinfos, lines, stations, freq, 
                appRHO, phase, appRho_w, phase_w)  
               
  
        
        

    def plot_freqVSRhoPhase (self, fn = None , profile_fn =None , 
                             station_id = 1, rename_stations =None , 
                             **kwargs ):
        """
        Method to plot apparent resistivity |phase vs frequency .
        
        :param fn:  full path to [AVG|EDI|J] file  
        :type fn: str 
        
        :param profile_fn:  full path to station profile . 
        :type profile_fn: str  
        
        .. note:: if user use drectly *AVG data 
                 must provide station profile '*.stn'
            
        =================  ==============  ===================================
        Others params       Default             Description 
        =================  ==============  ===================================
        station_id          str or int      plot the name of station if string 
                                            is povided make be sure that 
                                            the station name is on the station  
                                            list eg : station_id = 1 means 
                                            plot  S00 station =[1,13] means 
                                            plot >S00,S12station [S05, 7, 8]
                                            -- [ S05, S06, S07]
        rename_stations     list            bring the station name . Be sure
                                            the length of station name you 
                                            provided match the data station name  
        show_error          bool            if True , see errobar plot.  
                                            *Default* is False.
        =================  ==============  ===================================
        """

        savefig=kwargs.pop('savefig', None)
        orient = kwargs.pop('orientation', 'portrait')
        rename_stations = kwargs.pop('rename_stations', None)
        fontstyle = kwargs.pop('font_style', 'italic')
        show_error_bar= kwargs.pop('show_error', False)
        
        # creat CSAMT Object : 
        csamt_obj = CSAMT(data_fn = fn , profile_fn = profile_fn   )
    
        if csamt_obj.freq is None :
            raise CSex.pyCSAMTError_plot('Need absolutely frequency'
                                         ' value before plotting.'
                                         'Please add your frequency data.')
        if csamt_obj.resistivity  is None :
            raise CSex.pyCSAMTError_plot(
                'Error plotting.Provide your resistivity values.')
        
        #--- assert stations length ------ 
        stations = sorted(list(csamt_obj.resistivity.keys()))
        if rename_stations is not None : 
            if  len(stations )!= len(rename_stations): 
                warnings.warn('Stations provided are the lenght = {0},'
                              ' The defalut stations length ={1}.'
                              ' Please provided new stations '
                              'list with the same length.'.format(
                                  len(rename_stations), len(stations)))
                raise CSex.pyCSAMTError_station(
                    'New stations provided '
                    'must have the same length with default stations. ')
            stations = rename_stations
            
        #-- call function to specifier the user demand ---     
        stations_for_plot =mplotus.get_stationid(stations =stations,
                                                 station_id=station_id)
        #--> set objets 
        csamt_res_obj, csamt_freq_obj, csamt_phase_obj =\
            csamt_obj.resistivity, csamt_obj.freq ,csamt_obj.phase 
        csamt_res_err_obj , csamt_phs_err_obj =\
            csamt_obj.resistivity_err , csamt_obj.phase_err 

        #---> loop stations 

        for stn  in stations_for_plot: 
            mpl.rcParams['figure.figsize']=[8,6]
            
            fig =plt.figure()
            
            #---- > create axis for each plot 
            axe1 =plt.subplot2grid(shape=(3,3), loc=(0,0),
                                   rowspan=2, colspan=3)
            axe2 =plt.subplot2grid(shape=(3,3), loc=(2,0), 
                                   colspan = 3)
            
            if show_error_bar : 
                axe1.errorbar(csamt_freq_obj ,csamt_res_obj[stn],
                              yerr=csamt_res_err_obj[stn], 
                               fmt='none', ecolor = 'r', lolims=True,
                               uplims=True, 
                               xlolims =True, xuplims=True, 
                               lw=0.7, marker = '.', 
                                color='r',)
                
            
            mark,  = axe1.loglog (csamt_freq_obj ,csamt_res_obj[stn], 
                                  c= 'white',
                                  marker =self.marker,
                                  markersize = self.ms* self.fs ,
                                  markeredgecolor= self.markeredgecolor)
            
            
            logfreqres= axe1.loglog(csamt_freq_obj ,
                                    csamt_res_obj[stn] ,
                                    lw =self.lw *self.fs, c='k',
                        ls =self.ls, 
                        label = 'App. resistivity curve')# ,marker ='D', 
                         #alpha = 1)   

            axe1.minorticks_on()    # customize specific grid
            axe1.grid(color='k', ls=':', lw =0.25, alpha=0.8, which ='major') 
       
            axe1.set_xlabel('Frequency (Hz)' ,
                            fontdict={'color': 'k',
                                      'size': self.x_minorticks*14,
                                      'weight':self.fontweight,
                                      'style': fontstyle})
            
            axe1.set_ylabel ('Resistivity(Ω.m)', color='black',
                             fontsize =  self.y_minorticks*14,
                             fontweight=self.fontweight, 
                             fontstyle = fontstyle)
            
            axe1.legend( logfreqres, ['App. resistivity curve'])
            
            # axe1.set_title('log')
            
            #------------------ define axe2 ------------- 
            if show_error_bar : 
                 axe2.errorbar(csamt_freq_obj ,
                               csamt_phase_obj[stn],
                               yerr=csamt_phs_err_obj[stn], 
                                               fmt='none', 
                                               ecolor = 'r', 
                                               lolims=True, 
                                               uplims=True, 
                                               xlolims =True, xuplims=True, 
                                               lw=0.1, marker = '.', 
                                                color='yellow',)
            
            marki, = axe2.semilogx (csamt_freq_obj , csamt_phase_obj[stn],
                                    c= 'white', marker ='o', alpha =1,
                                    markersize = self.ms*self.fs ,
                                        markeredgecolor= 'cyan')
            
                
            
            logfreqphase =axe2.semilogx(csamt_freq_obj , 
                                        csamt_phase_obj[stn] ,
                                        lw=self.lw *self.fs,
                                        ls =':' ,
                                        c='k',
                                        label ='Phase curve') 
            

                                        
            axe2.set_xlabel('Frequency(Hz)', 
                            fontdict={'color': 'k', 
                                      'size':self.x_minorticks *14, 
                                      'weight':'bold', 
                                      'style': fontstyle})
            axe2.set_ylabel('Phase(°)',
                            fontdict={'color': 'k',
                                                  'size':self.y_minorticks *14, 
                                                  'weight': 'bold', 
                                                  'style':fontstyle})
            axe2.minorticks_on()
            axe2.grid(color ='k', ls=':', lw=0.25, alpha= .8 , which ='minor')
            axe2.legend( logfreqphase, ['Phase curve'])
            

            
            
            
            plt.tight_layout()
            fig.suptitle(' ResPhase plot: station <{0}>'.format(stn), 
                         fontsize= self.font_size*4,
                         bbox =dict(boxstyle='round',facecolor ='whitesmoke'),
                         fontstyle =fontstyle)    
            
            if savefig is not None : 
                plt.savefig(savefig, dpi=self.fig_dpi, orientation=orient)
            
            plt.draw()
        
    def penetration1D(self, fn =None , profile_fn =None , 
                     selected_frequency =None,  **kwargs):
        """
        Pentration1D depth : Show skin depth at selected frequencies .
        for multiples frequencies , put argument `selected_frequency` on list.
        If frequency provided is not on the frequency range , it will be interpolated.
   
        :param fn: full path to [AVG|EDI|J] file 
        :type fn: str 
        
        :param profile_fn: full path to *stn* station file . If user used  EDI or
                    J files , Dont need to add profile_file 
        :type profile_fn: str 
            
        :param selected_frequency: list ,  list of freauency want to see the  
                                   penetration depth. must be on a list  . 
                                   i.e [8, 511,1024 ]
        :type selected_frequency: list 
        
        =================  ===============  ==================================
        Params             Default          Description 
        =================  ===============  ==================================
        rename_station      list            Bring the station name . 
                                            Be sure the length of station name 
                                            you provided match  the size of 
                                            the data station name .
        rotate_station      int             rotation station name  . 
                                            *Default* is 90 degree. 
        fs                  float           can change the size of marker. 
                                            *Default* is .7 : eg ms =9*fs
        lw                  float           change the linewdth 
        plot_grid           bool            add grid on your plot . 
                                            Default is False 
        verbose             bool            Control the level of verbosity 
        =================  ===============  ==================================
        
        .. note:: browse to see others plot config.
        
        :Example: 
            
            >>> from viewer.plot import Plot1d 
            >>> path =  os.path.join(os.environ["pyCSAMT"], 
            ...              'pycsamt','data', file_1)
            >>> plot_1d_obj =Plot1d()
            ... plot1d_depth = plot_1d_obj.penetrated1D(fn =path ,
            ...                    profile_fn= os.path.join(
            ...                        os.path.dirname(path), 'K1.stn'), 
            ...                    selected_frequency =511)
        """
        self._logging.info('PlotPenetration1D datapath <%s>'% fn)
        
        rotate_stn =kwargs.pop('rotate_station_names', 90)
        plotgrid=kwargs.pop('plot_grid', False)
        rename_station = kwargs.pop('rename_station', None)
        
        fs =kwargs.pop('fs', 0.7)
        lw = kwargs.pop('lw', 1.5)
        savefig = kwargs.pop('savefigure', None)
        orient =kwargs.pop('orientation', 'landscape')
        verbose = kwargs.pop('verbose', False)

        #------BUILD CSAMTOBJ------------------------------------
        csamt_obj = CSAMT(data_fn =fn , profile_fn =profile_fn )
        csamt_dep1D_obj = csamt_obj.skindepth 
        csamt_freq_obj = csamt_obj.freq 
        csamt_stndis_obj= csamt_obj.station_distance
        
        
        if selected_frequency is None : 
            warnings.warn ('You may selected frequency between frequency range.'
                           ' If you dont know the number of frequency')
            raise ValueError('No frequency detected! Please provide ' 
                             'at least one frequency value in Hz.'
                           )
        
        def depth1D (freq_selected): 
            """
            selected depth penetrate that match the frequency 
            return depth1D array at that frequency.
            
            """
            return Zcc.get_data_from_reference_frequency(
                array_loc=csamt_dep1D_obj,
                freq_array=csamt_freq_obj, 
                reffreq_value=freq_selected)
        
        freqSELECT = mplotus.get_frequency_id(freq_array =csamt_freq_obj, 
                                              frequency_id= selected_frequency)
   
        fig = plt.subplot
        
        mpl.rcParams['figure.figsize']=[12,6]
        # mpl.rcParams['font.family']='sans-serif'
        mpl.rcParams['font.sans-serif']=['Tahoma']
        fig =plt.figure()
        axis =fig.add_subplot(111)
        
        # create a twin y to store the name of station . 
        axis.set_xlabel('Distance(m)', 
                        fontdict={'color': 'k', 
                                    'size':self.x_minorticks *14,
                                    'weight':'bold'})
        axis.set_ylabel('Penetration depth (m)',
                        fontdict={'color': 'k',
                                    'size':self.y_minorticks *14, 
                                    'weight': 'bold'})
        
        # #---set the second axis ---------
        axis2 = axis.twiny()
        axis2.set_xlabel ('stations' ,
                          fontdict={'color': 'k',
                                    'size':self.x_minorticks *14, 
                                    'weight':'bold'})
        
        csamt_stn_num_obj=sorted(list(csamt_obj.skindepth.keys()))
        axis2.set_xticks(ticks= csamt_stndis_obj, minor=False )
        if rename_station is not None : 
            assert len(rename_station)==len(csamt_stn_num_obj),\
                CSex.pyCSAMTError_plot("Error plot !rename_station and station"
                                       " name must have the same lenght.")
        
            csamt_stn_num_obj= rename_station 
            
        axis2.set_xticklabels(csamt_stn_num_obj, rotation=rotate_stn)
        # axis2.minorticks_on()# create legend ob for appending 
        hand_leg, leglabel, depmax, freqmax =[], [] , 0 ,0         
        
        for freqs in freqSELECT:
          
            # mark, = axis.plot(csamt_stndis_obj, dep1D, marker ='*',
            # markersize =self.ms*2fs , markeredgecolor='blue')
            # recover the interpolated frequency 
            if freqs not in  csamt_freq_obj: 
                interpFreq =Zcc.find_reference_frequency(
                    freq_array=csamt_freq_obj, 
                     reffreq_value=freqs,sharp=True, etching=False)
                if verbose:
                    warnings.warn ('Frequency {0} not in frequency range. '
                                   'It will be interpolated to find '
                                   'maximum closest  frequency.'.format(freqs))
                    print('---> Input frequency <{0}> Hz has been interpolated to'
                          ' <{1}Hz>.'.format(freqs,  float(interpFreq)))
                freqs = float(interpFreq)
                
                
            penetration1d,  = axis.plot(csamt_stndis_obj, -1*depth1D(freqs) , lw=lw,
                                        marker ='*', markersize =self.ms*3*fs , 
                                        label = 'freq {0}'.format(freqs)) 
            
            if verbose: 
                print('---> Max depth reached by freq={0} Hz is <{1}> km.'.\
                      format(freqs,
                             np.around(depth1D(freqs).max()*1e-3,2)))
            
            axis.grid(color ='k', ls=':', lw=0.25, alpha= .8 , which ='minor')
    
            axis.set_xlabel('Distance(m)',  fontdict={'color': 'k',
                                                      'size':self.x_minorticks *14,
                                                      'weight':'bold'})
            axis.set_ylabel('Penetration depth (m)',
                            fontdict={'color': 'k', 
                                      'size':self.y_minorticks *14, 
                                      'weight': 'bold'})
            
            axis.set_xlim([csamt_stndis_obj.min(),
                           csamt_stndis_obj.max()])
            # axis.minorticks_on()
            if plotgrid :
                axis.grid(color ='k', ls=':', lw=0.25, alpha= .8 ,
                          which ='major')
            
            # legend 
            hand_leg.append(penetration1d), leglabel.append(
                '$f={0}Hz$'.format(int(freqs)))
            axis.legend( hand_leg, leglabel)
            #axis.legend(bbox_to_anchor=(1.05, 1), loc='upper left',
            # borderaxespad=0.) #place al legend smaller upper left 
            
            if depmax < depth1D(freqs).max() :
                depmax,  freqmax = depth1D(freqs).max(),freqs
            
            #axis.set_title(' Penetraton depth plot : max depth = {0}
            #km max depth'.format(np.around(depmax*1e-3,2)),  
            #fontsize= self.font_size)
                    # #---set the second axis ---------
        
        
        if len(freqSELECT) >1 :fmt ='frequencies'
        else :fmt='frequency'
        if verbose:
            print('---> On the set of {0} {1}, max depth reached '
                  '={2} km at freq ={3} Hz.'.format(len(freqSELECT),
                                                    fmt, np.around(depmax*1e-3,2), 
                                                    int(freqmax)  ))
        # fig.suptitle(' Penetration depth plot at {0} {1} '.format(len(freqSELECT), fmt),
        #              fontsize= self.font_size, verticalalignment='center', 
        #              )
        plt.tight_layout()
        fig.suptitle(
            ' Penetration depth plot at {0} {1} '.format(len(freqSELECT), fmt), 
                    fontsize= self.font_size*4,
                    bbox =dict(boxstyle='round',facecolor ='whitesmoke'),
                    fontstyle ='italic') 
        if savefig is not None :
            plt.savefig(savefig, dpi=self.fig_dpi, orientation =orient )

        plt.show()
        
    def plot_curves (self, fn =None , savefig =None ,
                     selected_stations =1,  ** kws): 
        """
        Plot Zonge Engineering AVG file with different components 
        E and H at differents frequencies.  
        
        :param fn: full path to Zonge Engineering file 
        :type fn: str 
        
        :param profile_fn:   full path to profile file .
        :type profile_fn: str 
        
        :param savefig:  path to figure plot
        :type savefig: str 
        
        =================  ===========  ======================================
        Params             Default      Description 
        =================  ===========  ======================================
        fs                  float       can change the size of marker. 
                                        *Default is .7 : eg ms =9*fs
        lw                  float       change the linewdth 
        error_bar           bool        set to false to let invisible. 
                                        Default is True 
        =================  ===========  ======================================
                
        :Example :
            
            >>> from viewer.plot import Plot1d 
            >>>path =  os.path.join(os.environ["pyCSAMT"], 
            ...          'pycsamt','data', file_1)
            >>> plot_1d_obj =Plot1d()
            >>> plotcurves = plot_1d_obj .plot_curves(fn = path, 
                                                    selected_stations=[1,10, 20], 
                                                      error_bar=True)
        """
        self._logging.info('Plot curves from <%s>'% fn)
        
        #----- call Avg class and Define CSAMTobj ----
        fontstyle =kws.pop ('font_style', 'italic')
        orientation =kws.pop('orientation', 'landscape')
        errBar = kws.pop('error_bar', False)
        xylabelsize =kws.pop('xylabel_size',10)
        
        
        #-----------IMPORT OBJ AND DEFINE COMPONENTS ------------------

        zonge_csamt_obj = CSMATavg.Avg(data_fn= fn )
        zonge_freq_obj = zonge_csamt_obj.Data_section.Frequency.value 
        zonge_res_obj = zonge_csamt_obj.Data_section.Resistivity.loc 
        zonge_emag_obj = zonge_csamt_obj.Data_section.Emag.loc 
        zonge_hmag_obj = zonge_csamt_obj.Data_section.Hmag.loc 
        zonge_stn_obj = zonge_csamt_obj.Data_section.Station.names 
        
        #---------define _error _obj -----------------------------
        zonge_res_err_obj = {stn : value *100 for stn ,
                  value in zonge_csamt_obj.Data_section.pcRho.loc.items()} 
        zonge_phase_err_obj= {stn : (180 * value/1e3) for stn , 
                 value in  zonge_csamt_obj.Data_section.sPhz.loc.items()}
        zonge_emag_err_obj = {stn : value *100 for stn , 
                     value in zonge_csamt_obj.Data_section.pcEmag.loc.items()}
        zonge_hmag_err_obj = {stn : value *100 for stn ,
                     value in zonge_csamt_obj.Data_section.pcHmag.loc.items()}
        
        
        #--convert phase from mrad to degree 
        zonge_phase_obj = {key : (180 * value *1e-3/ np.pi)%90 for key , 
                           value in  zonge_csamt_obj.Data_section.Phase.loc.items()}
        
        # create figure from grid2supplots 
        mpl.rcParams['figure.figsize']=[12,7]
        fig =plt.figure()
        
        #---- > create axis for each plot from subplot2grid-------------
        
        axe_res =plt.subplot2grid(shape=(3,5), loc=(0,0), rowspan=2, colspan=3)
        axe_phase =plt.subplot2grid(shape=(3,5), loc=(2,0),  colspan=3)
        axe_emag = plt.subplot2grid(shape=(3,5), loc=(0,3), colspan =2)
        axe_hmag = plt.subplot2grid(shape=(3,5), loc=(1,3), colspan =2)
        #----> controlled the stations selected for plots --- 
        stnID = mplotus.get_stationid(stations=zonge_stn_obj,
                                      station_id=selected_stations)
        #for legend manager
        handleg , lableg  =[
            {'res':[], 'phs':[], 'emag':[], 'hmag':[]} for ii in range(2)]                 
        
        for stn_id in stnID :
            #-----RESISTIVITY PLOT -----------------------
            mark, =axe_res.loglog(zonge_freq_obj,
                                  zonge_res_obj [stn_id],  c='white', 
                                marker ='o', markersize =self.ms*self.fs ,
                                markeredgecolor = self.markeredgecolor)
            rho,  = axe_res.loglog(zonge_freq_obj,zonge_res_obj [stn_id] , 
                                   lw= self.lw,
                                   label = 'station {0}'.format(stn_id))
            axe_res.set_xlabel('Frequency (Hz)', 
                               fontdict={'color': 'k',
                                        'size':self.x_minorticks *xylabelsize,
                                        'weight':self.fontweight,
                                        'style' :fontstyle})
            axe_res.set_ylabel('Apparent resistivity(Ω.m)', 
                               fontdict={'color': 'k',
                                         'size':self.y_minorticks *xylabelsize,
                                          'weight': self.fontweight, 
                                          'style':fontstyle})
            
            handleg['res'].append(rho), lableg['res'].append(
                '$station\ {0}$'.format(stn_id))
            
            if errBar is True : 
                for xpos , ypos , err in  zip(
                                            zonge_freq_obj,zonge_res_obj[stn_id], 
                                              zonge_res_err_obj[stn_id]): 
                    axe_res.errorbar (xpos, ypos, err ,uplims=True, lolims=True, 
                                      ecolor ='r',lw= 0.7,
                                      #capsize =4 ,
                                      marker = '.', 
                                      #ms=7, 
                                      color='blue',
                                      #color ='magenta'
                                      )
            
            #---------PHASE PLOT --------------------------
            markphs, =axe_phase.semilogx(
                                zonge_freq_obj,zonge_phase_obj [stn_id],
                                c='white', 
                                marker ='o', 
                                markersize =self.ms*self.fs/2 ,
                                markeredgecolor=self.markeredgecolor)
            phase,  = axe_phase .semilogx(zonge_freq_obj,
                                          zonge_phase_obj [stn_id] ,
                                          lw= self.lw, 
                                          label = 'station {0}'.format(stn_id))
            axe_phase .set_xlabel('Frequency (Hz)', 
                                  fontdict={'color': 'k',
                                            'size':self.x_minorticks *xylabelsize,
                                            'weight':self.fontweight,
                                            'style' :fontstyle})
            axe_phase .set_ylabel('Phase(°)',
                                  fontdict={'color': 'k', 
                                            'size':self.y_minorticks *xylabelsize,
                                            'weight': self.fontweight,
                                            'style':fontstyle})
            axe_phase.set_ylim([0, 90])
            handleg['phs'].append(phase), lableg['phs'].append(
                '$station\ {0}$'.format(stn_id))
            if errBar is True : 
                for xpos , ypos , err in  zip(zonge_freq_obj,
                                              zonge_phase_obj [stn_id],
                                              zonge_phase_err_obj[stn_id]): 
                    axe_phase.errorbar (xpos, ypos, err ,
                                        uplims=True, lolims=True, 
                                      ecolor ='r',lw= 0.7,
                                      #capsize =4 ,
                                      marker = '.', 
                                      #ms=7, 
                                      #capthick=4,
                                      color ='magenta')
            
            #---------EMAG PLOT -----------------------------
            
            markemag, =axe_emag.loglog(zonge_freq_obj,
                                       zonge_emag_obj [stn_id], 
                                       c='white', 
                                marker ='o', markersize =self.ms*self.fs ,)
            emag,  = axe_emag.loglog(zonge_freq_obj,zonge_emag_obj [stn_id] ,
                                     lw= self.lw, 
                                     label = 'station {0}'.format(stn_id))
                                        #marker ='*', markersize =self.ms*self.fs , ) 
            axe_emag.set_xlabel('Frequency (Hz)', 
                                fontdict={'color': 'k',
                                          'size':self.x_minorticks *xylabelsize, 
                                          'weight':self.fontweight,
                                          'style' :fontstyle})
            axe_emag.set_ylabel('E-magn.(microV/Km*A)',
                                fontdict={'color': 'k',
                                          'size':self.y_minorticks *xylabelsize,
                                          'weight': self.fontweight, 
                                          'style':fontstyle})
            
            handleg['emag'].append(emag), lableg['emag'].append(
                '$station\ {0}$'.format(stn_id))  
            
            if errBar is True : 
                for xpos , ypos , err in  zip(zonge_freq_obj,
                                              zonge_emag_obj [stn_id],
                                              zonge_emag_err_obj[stn_id]): 
                    axe_emag.errorbar (xpos, ypos, err ,
                                       uplims=True, lolims=True, 
                                      ecolor ='dimgray',lw= 0.7,
                                      #capsize =4 ,
                                      marker = '|', 
                                      #ms=7, 
                                      #capthick=4,
                                      color ='k')
            
            #-----------HMAG PLOT --------------------------------------------
            markhmag, =axe_hmag.loglog(zonge_freq_obj, zonge_hmag_obj [stn_id],
                                c='white', 
                                marker ='o', markersize =self.ms*self.fs ,
                                alpha =0.8 )
            hmag,  = axe_hmag.loglog(zonge_freq_obj,zonge_hmag_obj [stn_id] ,
                                     lw= self.lw,
                                     label = 'station {0}'.format(stn_id))
                                        #marker ='*', markersize =self.ms*self.fs , ) 
            axe_hmag.set_xlabel('Frequency (Hz)', 
                                fontdict={'color': 'k', 
                                          'size':self.x_minorticks *xylabelsize,
                                          'weight':self.fontweight, 
                                          'style' :fontstyle})
            axe_hmag.set_ylabel('H-magn.(picoT/A)', 
                                fontdict={'color': 'k',
                                          'size':self.y_minorticks *xylabelsize, 
                                          'weight': self.fontweight, 
                                          'style':fontstyle})
            
            handleg['hmag'].append(hmag), lableg['hmag'].append(
                '$station\ {0}$'.format(stn_id))  
            
            if errBar is True : 
                for xpos , ypos , err in  zip(zonge_freq_obj,zonge_hmag_obj [stn_id],
                                              zonge_hmag_err_obj[stn_id]): 
                    axe_hmag.errorbar (xpos, ypos, err ,uplims=True,
                                       lolims=True, 
                                      ecolor ='dimgray',lw= 0.7,
                                      #capsize =4 ,
                                      marker = '|', 
                                      color ='k')
        
        # set all legend
        
        for axe , comps in zip ([axe_res, axe_phase, axe_emag, axe_hmag],
                                ['res', 'phs', 'emag', 'hmag']): 
            axe.legend(handleg[comps], lableg[comps])
        #-----------------append necesary infos -------------------------------
        proj_name =zonge_csamt_obj.Header.SurveyAnnotation.project_name
        custumer_info  = zonge_csamt_obj.Header.SurveyAnnotation.custumer_name
        contractor= zonge_csamt_obj.Header.SurveyAnnotation.contractor_name
        acqdate =zonge_csamt_obj.Header.SurveyAnnotation.acqdate
        
        surveyType =zonge_csamt_obj.Header.SurveyConfiguration.surveyType
        surveyline = zonge_csamt_obj.Header.SurveyConfiguration.lineName
        surveyazim =zonge_csamt_obj.Header.SurveyConfiguration.azimuth

        txTypeinfo =zonge_csamt_obj.Header.Tx.txType 
        txcentinfo =zonge_csamt_obj.Header.Tx.txCenter
        rxcomp =zonge_csamt_obj.Header.Rx.rxComps
        # create grid plot from info 
        axe_info = plt.subplot2grid(shape=(3,5), loc=(2,3), colspan =2)
        
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

        for ii, info in enumerate ([proj_name,
                                    custumer_info , 
                                    contractor,acqdate,
                                  surveyType,
                                  surveyType, 
                                  surveyline,
                                  surveyazim ,
                                  txTypeinfo,
                                  txcentinfo,rxcomp   ]): 
            
            axe_info.text(0,1-ii/10, 
                          '${0}$'.format(info), 
                          fontdict={'size':7, 
                                   'color':'k'}, fontstyle ='italic', 
                          verticalalignment ='center' , 
                          bbox = props,
                         )
            axe_info.axis('off') # let not visible axis of info 
            
        #--------------------------figure 
        fig.suptitle('Plot curves rho, phase, E- and Hy-magnitudes',
                     fontsize= 3.5 * self.font_size, 
                     verticalalignment='center', 
                     style =fontstyle,
                     bbox =dict(boxstyle='round',facecolor ='moccasin'))
        
        plt.tight_layout(h_pad =1.8, w_pad =1.08)
        
        if savefig is not None : plt.savefig(savefig,
                                             dpi = self.fig_dpi, 
                                             orientation =orientation)
        plt.show()
    
    def  plotRMS(self, fn=None ,target =1. ,savefig =None,  **kwargs): 
        """
        Plot RMS . If occamlogfile is not available , set rms value , 
        iteration value at each rms and|or roughness.
        
        Parameters
        -----------
           * fn : str 
               full path to occam2D logfile 
           * savefig : str 
               full directory to save fig 
           * target : float 
               target supposed RMS to reach . Default is 1. 
        
    
        Returns 
        -------
            obj , plot_RMS obj 
            
            
        ============  ===============  ========================================
        Params          Type            Description    
        ============  ===============  ========================================
        rms            array_like       RootMeanSquare array . 
        iteration      int              number of interation reached .
                                        iteration starts from "0". so we will
                                        add the number provided plus 1.
        roughness      array_like       deGootHeldlin roughness parameters. 
                                        number of params =Num(RMS)-1 . so we 
                                        excluded the starting R.M.S
        target         float            RMS target weexpected to reach .
                                        Default is   1.0
        ============  ===============  ========================================
        
         .. note:: If occam2d logfile is availbale , dont need other parameters ,
             except the path "fn" and as possible the "target". 
             
        """
       
        orientation =kwargs.pop('fig_orientation', None)  
        plot_target=kwargs.pop ('show_target_line', True)
        show_grid = kwargs.pop('show_grid', False)
        
        if fn is not None : 
            occam2d_obj =occam2d.occamLog(fn =fn )
            csamt_rms_obj = occam2d_obj.rms
            csamt_iteration_obj =occam2d_obj.iteration 
            csamt_roughness_obj =occam2d_obj.roughness
            
        for key in list(kwargs.keys()):
            if key in ['rms', 'iteration', 'roughness']:
                if kwargs['rms'] is not None :
                    csamt_rms_obj = kwargs['rms']
                if kwargs['iteration'] is not None :
                    csamt_iteration_obj= kwargs['rms']
                if kwargs['roughness'] is not None :
                    csamt_roughness_obj =kwargs['roughness']
           
        for comp in [csamt_rms_obj, csamt_iteration_obj ]:
            if comp is None : 
                warnings.warn (
                    'R.M.S and Iteration value are required for plotting.'
                        ' Can not plot "None" value.')
                self._logging.warn(
                    'Error plotting R.M.S VS Iteration. Can not plot "None" value.')
                raise CSex.pyCSAMTError_plot(
                    'Error Plot R.M.S . '
                    'Compulsory need R.M.S value and '
                        'Iteration value. Please check your value. ')
                
        # append None to roughness value for last iteration if not exists
        if csamt_roughness_obj is not None :
            if csamt_roughness_obj.size != csamt_iteration_obj.size :
                addnum = csamt_iteration_obj.size - csamt_roughness_obj.size 
                for kk in range (addnum): csamt_roughness_obj = np.append(
                        csamt_roughness_obj, None)

        #---------------------------------------------------------------------
        mpl.rcParams['figure.figsize']=[10,4]
        mpl.rcParams["image.origin"] ='upper' 
        if orientation is None : orientation = 'landscape'

        fig, axrms =plt.subplots()
        
        leghandles, leglabels = [], []

        #---------RMS 
        RMS_axis, = axrms.plot(csamt_iteration_obj, csamt_rms_obj, c='k',
                   ls =self.ls , 
                   lw =self.lw , marker ='o', 
                    markeredgecolor='k',
                    markerfacecolor='red',)
        
        axrms.set_xlabel ('Iteration number', 
                          fontdict ={'size': 4*self.font_size,
                                        'c': 'k',
                                        'weight':self.fontweight, 
                                        'style':'italic', 
                                        })
        axrms.set_ylabel ('R.M.S', fontdict ={'size': 4* self.font_size,
                                                  'c': 'k',
                                                  'weight':self.fontweight, 
                                                  'style':'italic', 
                                                  })

        leghandles.append(RMS_axis), leglabels.append('$R.M.S$')
        dx= csamt_iteration_obj.min()/2
        axrms.set_xlim(csamt_iteration_obj.min()- dx,
                       csamt_iteration_obj.max()+dx)
        #see only tick on yaxis 
        axrms.minorticks_on()
        axrms.xaxis.set_tick_params(which='minor',
                                    bottom=False)

        

        if show_grid :
            axrms.grid(axis ='y', color ='k', ls=':', lw=0.25, alpha= .8 , 
                       which ='major')
        
        #------plotROUGHNESS 
        if csamt_roughness_obj is not None : 
            axrough = axrms.twinx()
            
            ROUGH_axis, = axrough.plot(csamt_iteration_obj, csamt_roughness_obj,
                                       lw=self.lw , ls='-', c='gray', 
                     marker ='o', markeredgecolor='blue', markerfacecolor='gray')
            axrough.set_ylabel ('Roughness param. value', 
                                fontdict ={'size': 4* self.font_size,
                                                      'c': 'k',
                                                      'weight':self.fontweight, 
                                                      'style':'italic', 
                                                      })
            axrough.minorticks_on()
            axrough.xaxis.set_tick_params(which='minor', bottom=False)
            # axrough.tick_params(axis='x', which='minor', bottom='off')
            leghandles.append(ROUGH_axis), leglabels.append('$Roughness\ value$')
        
        #----PLOT TARGET
        if plot_target :
            
            marki, = axrms.plot(csamt_iteration_obj[0],
                                target,
                                c='white' ,
                                marker ='*', 
                            markeredgecolor='k',
                            markerfacecolor='red',
                            markersize = 3*self.fs)
            
            ytarget_= np.repeat(target,csamt_iteration_obj.size)
            TARG_axis, = axrms.plot(csamt_iteration_obj,  ytarget_,
                                    ls =':' ,
                                    lw =self.lw, 
                                    c='k'
                       )
            
            leghandles.append(TARG_axis), leglabels.append('$R.M.S\ Target$')
        
        # ----set LEGEND 
        axrms.legend(leghandles,leglabels)
        if csamt_roughness_obj is not None : fmt=' and Roughness'
        else: fmt=''

        fig.suptitle('Plot R.M.S{0}'.format(fmt),
                     fontsize= 4 * self.font_size, 
                     verticalalignment='center', 
                     style ='italic',
                     bbox =dict(boxstyle='round',
                                facecolor ='moccasin'))
        
        plt.tight_layout(h_pad =1.8, w_pad =1.08)
        
        if savefig is not None : plt.savefig(savefig,
                                             dpi = self.fig_dpi,
                                             orientation =orientation)
        
    def plot_station_profile (self, fn =None , straighten_type ='classic', 
                              reajust_coordinates=(0,0), savefig =None, 
                              **kwargs ):
        """
        Method to plot original station profile and coordinate reajustment
         profiles. Deal with Zonge AVG file .
        
        Parameters
        -----------
            * fn : str   
                full path to profile station file  of Zonge Engineering station
                profile file . format egal to *.stn
                
            * straighten_type : str 
                type of straingther profile   
                it may be `classic`, `equisistant` or `distord`
                *Default* is 'classic'
                
            * reajust_coordinates : list 
                list of float x, y values  
                                        
        :Example: 
            
            >>> path =  os.path.join(os.environ["pyCSAMT"], 
            ...          'pycsamt','data', 'avg', 'K1.stn')
            >>> plot_1d_obj= Plot1d()
            >>> plot_1d_obj.plot_station_profile(fn = path)
        """
        orientation= kwargs.pop('orientation', 'landscape')
        alpha  =kwargs.pop('alpha', 0.5)
        font_style =kwargs.pop('font_style', None)
        fw =kwargs.pop('font_weight', 'bold')
        output_reajfile =kwargs.pop('outputfile', False)
        
        # call profile _obj abd get special attribute 
        profile_obj =Profile(profile_fn= fn )
        prof_east_obj = profile_obj.east
        prof_north_obj = profile_obj.north 
        # apply readjustment 
        profile_obj.straighten_profileline(X=prof_east_obj, 
                                           Y=prof_north_obj, 
                                           straight_type=straighten_type, 
                                           reajust=reajust_coordinates, 
                                           output=output_reajfile, 
                                           savepath = savefig
                                           )
        
        prof_east_obj = profile_obj.e_east
        prof_north_obj = profile_obj.n_north 
        
        new_prof_east_obj = profile_obj.east 
        new_prof_north_obj = profile_obj.north 

        #get profile name though Site_obj
        #easting =new_prof_east_obj , northing =new_prof_north_obj)
        profile_obj.Site.set_site_info(data_fn =fn)

        profile_station_names =profile_obj.Site.stn_name
        
        # plot station 
        mpl.rcParams['figure.figsize']=[10,6]
        mpl.rcParams["image.origin"] ='upper' 
        
        fig =plt.figure()
        axe =fig.add_subplot(111)
        mark, =  axe.plot(prof_east_obj , prof_north_obj,
                          marker =11, c='white', 
                          markersize= 3* self.ms , 
                          markeredgecolor= 'k', )
        
        unscalled, = axe.plot(prof_east_obj , prof_north_obj, 
                              lw =self.lw , 
                              ls =self.ls , 
                              marker = '1',
                              c='dimgray',
                              markersize = 3* self.ms , 
                              markeredgecolor='k', 
                              markerfacecolor='k')
        mark, =  axe.plot(new_prof_east_obj, new_prof_north_obj,
                          marker ='o', c='white', 
                          markersize= 3* self.ms , 
                          markeredgecolor= 'blue', )

        scalled, = axe.plot(new_prof_east_obj, new_prof_north_obj,
                            lw=2* self.lw, 
                            ls= self.ls , 
                            c='k',
                            marker ='*', 
                            markersize= 2* self.ms , 
                            markeredgecolor= 'r', 
                            markerfacecolor='k')
        axe.minorticks_on()
        axe.grid(color='k', ls=':', lw =0.25, alpha=alpha, which ='major')
        axe.set_xlabel('Easting(m)',  fontdict ={'style': font_style, 
                                                  'size': 4* self.font_size ,
                                                  'weight': fw} )
        axe.set_ylabel('Northing({0})'.format('m'),
                       fontdict ={'style': font_style, 
                                                'size': 4* self.font_size ,
                                                'weight': fw})
        #annotate station
        #d=find dx 
        dx = np.sqrt(profile_obj.dipole_length)*2
        for ii , stn in enumerate(profile_station_names) : 
            axe.annotate(stn, 
                         #(new_prof_east_obj[ii],new_prof_north_obj[ii]), 
                         xy=(new_prof_east_obj[ii]+dx, new_prof_north_obj[ii]+dx), 
                         #xycoords='figure points'
                         )
                
        
        print('----> {0} stations have been plotted and reajusted successfully.'.\
              format(len(profile_station_names )))
        
        fig.suptitle('Plot stations',
                     fontsize= 4 * self.font_size, 
                     verticalalignment='center', 
                     style ='italic',
                     bbox =dict(boxstyle='round',facecolor ='moccasin'))
        
        axe.legend([unscalled, scalled], ['$Unscalled$', '$Scalled$'])
        
        plt.tight_layout(h_pad =1.8, w_pad =1.08)
        
        if savefig is not None : 
            plt.savefig(savefig, dpi = self.fig_dpi, orientation =orientation)
        
        plt.show()
        
    def plot_multiStations(self, X =None , Y=None,  path =None, 
                           profile_lines =None,  **kwargs): 
        """
        Plot the survey area multiple stations   
        
        ..TODO:: For the future plan- added rose diagram from strikes angles   
            computation for all `edi`, `avg` or `j` files. Rose Diagram  
            actually does not work. We keep windrose plot for future version. 
            For the future version because  we will use windrose to 
            add Rose diagram plot. For instance:: 
                
                >>> from pycsamt.utils import func_utils as func 
                >>> IMPWINR= func.subprocess_module_installation('windrose')
                >>> if  IMPWINR:
                ...    import windrose
                ...    gs=gspec.GridSpec(1, 2, figure =self.fig)
                ...    axeProfiles = self.fig.add_subplot(gs[0, 0])
                ...    axeStrikes=self.fig.add_subplot(gs[0,1] ,
                ...                                    projection='polar' )
                
        Parameters
        -----------
            * path : str 
                full path to station profile path . In the case where 
                Zonge avg file is provided , use `stn` profile files. Group 
                all `stn` file on a folder will call automatically
                
            * profile_lines  : list 
                name of profile lines . if profile lines is NOne 
                will tale all `stn` profiles in the path directory 
                
            * X : list 
                list of arrays array of X coordinates values  for each
                 survey line. Can be easting or Northing  
                 
            * Y:  list 
                list of arrays  of Y coordinates values : can be easting 
                     or northing
                     
        .. note::  `X` and `Y` MUST be the same length  
        
        """
        
        savefig =kwargs.pop('savefig', None )
        straigthen_out_profile= kwargs.pop('straigthen_out_profile', True )
        show_station_labels  = kwargs.pop('show_station_labels', True)
        scale =kwargs.pop('scale', 'm')
        
        # scalled plots 
        if scale.lower() is None : scale ='m'
        if scale.lower() =='m' : 
            dz = 1. 
        elif scale.lower() =='km': 
            dz =1000.
        else :
            dz =1.
        
        for attr in ['X', 'Y'] : self.__setattr__(attr, None ) 
        
        if X is not None : self.X = X 
        if Y is not None : self.Y =Y
        if path is not None : self.path =path 
        # initialise containers  to hold eastings and northings values from all
        # survey lines 
        eastings , northings , profile_angles, gstrikes, snames =[
                                            [] for i in range(5)] 
        # assert len  is values are privided 
        
        if self.X is not None and self.Y is not None :
            # controle the length 
    
            if len(self.X) != len(self.Y):
                mess= ''.join(['X and Y must have the same length.', 
                               'The first index has length = {0} and the second has', 
                               ' length = {1}'])
                
                warnings.warn(mess.format(len(self.Y),len(self.Y)))
                self._logging.error(mess.format(len(self.Y),len(self.Y)))
                raise CSex.pyCSAMTError_plot(mess.format(len(self.Y),
                                                         len(self.Y)))
            if isinstance(self.X, np.ndarray) and isinstance(self.Y, np.ndarray): 
                #case where user profide only one line with easting and northings 
                # coordinates. if X and Y not on list , then put on list. 
                
                self.X, self.Y =[self.X] , [self.Y] 
        
        
        if self.path  is not None : 
            self._logging.info (
                'Reading Zonge Engeneering `station profile` file !')
            
            if profile_lines is not None : # create fle fn to read 
                profiles_path =[os.path.join( self.path , file) 
                                for file in profile_lines]
                profile_lines = [file.split('.')[0] for file in profile_lines]
            else : # take all profile in path directory 
                profiles_path = [os.path.join(self.path, file) for
                                 file in os.listdir(self.path) 
                                 if file.lower().endswith ('.stn')] # ckeck when
                profile_lines = [os.path.basename(os.path.splitext(pathfile)[0]) for
                                 pathfile in profiles_path] # keep 
                                        #the name of profile 
   
            # pfofiles are Stn files 
            if len(profiles_path) == 0 : 
                msg =''.join([
                    'No `.stn` profiles found . Please check your right path !', 
                    'or provide  X and Y as lists of coordinates values.'])
                self._logging.error(msg)
                warnings.warn(msg)
                raise CSex.pyCSAMTError_profile(msg)
                
            print('---> {0:02} *station profiles* detected !'.format(
                                            len(profiles_path)))
            for file in profiles_path : 
                profile_obj  = Profile( profile_fn=file) 
                if straigthen_out_profile is True : 
                    profile_obj.straighten_profileline() #straighen out profile 
                    
                profile_obj.get_profile_angle()
                # gstks, pangs, _=
                #geostrike.compute_geoelectric_strike(easting =profile_obj.east,
                #              northing = profile_obj.north)
                # now get value reajusted : 
                eastings.append(profile_obj.east)
                northings.append(profile_obj.north)
                profile_angles.append(profile_obj.profile_angle)
                gstrikes.append(profile_obj.geoelectric_strike)
        
                # get site names
                profile_obj.Site.stn_name= len(profile_obj.east)
                snames.append(profile_obj.Site.stn_name)
                
                
            self.X =eastings 
            self.Y =northings
            
        elif self.path is None and self.X is not None and self.Y is not None : 
            profile_obj =Profile()
            for  ii, slen in enumerate(self.X) : 
                profile_obj.Site.stn_name = len(slen) # build sitenames 
                snames.append(profile_obj.Site.stn_name)
                if straigthen_out_profile is True: # ccompute profile angle and 
                                                    #strike angle 
                    pangles, strikes = profile_obj.get_profile_angle(
                                                    easting =slen, 
                                                   northing =self.Y[ii]) 
                                                                                
                    gstrikes.append(strikes)
                    profile_angles.append(pangles)
                    
            profile_lines =['line_{0:02}' for i in range(len(self.X))]
                    
            print('---> {0:02} * survey* lines read !'.format(
                                            len(self.X)))

        # ------------------------DECLARE FIGURE AND PROPERTIES --------------
        
        # statement of figures 
        self.fig = plt.figure(self.fig_num, self.fig_size, dpi=self.fig_dpi,) 
        # plt.clf()
        # add a subplot to the figure with the specified aspect ratio
        self.fig_aspect ='auto'
                #--- using window rose to plot diagramm ------- 
                # -----**future plan ****

        axeProfiles = self.fig.add_subplot(1,1,1)
 
        # get min max for easting and northing profile
        Xmin, Xmax = self.X[0].min(), self.X[0].max()
        Ymin, Ymax = self.Y[0].min(), self.Y[0].max()
        for ii, (east, north)  in enumerate( zip( self.X , self.Y)):
            
            if east.min() < Xmin : 
                Xmin = east.min()
            if east.max() > Xmax : 
                Xmax = east.max ()
            if north.min() < Ymin : 
                Ymin = north.min()
            if north.max()> Ymax: 
                Ymax = north.max()
        
        self.xlimits=  (Xmin/dz, Xmax/dz ) 
        self.ylimits= (Ymin/dz, Ymax/dz )
        # -------------------------------------------------------------------
        leghandles= []
        for ii in range(len(self.X)): 
            
            mark, =  axeProfiles.plot(self.X[ii]/dz , 
                                      self.Y[ii]/dz ,
                          marker ='o' , 
                          c= 'white' , 
                          markersize= self.ms/2 , 
                          markeredgecolor= 'blue', )

            scalled, =axeProfiles.plot(self.X[ii]/dz ,
                                       self.Y[ii]/dz ,
                                lw= self.lw, 
                                ls= self.ls , 
                                c='k',
                                marker ='*', 
                                markersize=self.ms/3 , 
                                markeredgecolor= self.markeredgecolor , 
                                markerfacecolor= self.markerfacecolor,
                                )
            
            leghandles.append(scalled)
            
            #annotate station
            
            if show_station_labels  is True: 
                dx = np.sqrt(profile_obj.dipole_length)*2 
                # dx = profile_obj.dipole_length*2
    
                for jj , stn in enumerate(snames[ii]) : 
                    axeProfiles.annotate(stn,  
                                 xy=((self.X[ii][jj]+dx)/dz,
                                     (self.Y[ii][jj]+dx)/dz), 
                                 #xycoords='figure points'
                                 fontsize = 1.*self.font_size, 
                                 
                                 )
                    
        axeProfiles.minorticks_on()
        axeProfiles.grid(color='k', ls=':', lw =0.25,
                         alpha=self.alpha, which ='major')
        axeProfiles.set_xlabel('Easting({0})'.format(scale),
                               fontdict ={'style': self.font_style, 
                                          'size':  2*self.font_size ,
                                          'weight': self.fontweight} )
        axeProfiles.set_ylabel('Northing({0})'.format(scale),
                       fontdict ={'style': self.font_style, 
                                  'size': 2*self.font_size ,
                                  'weight': self.fontweight})
        
        #reduce the number of X ad Y  ticks density 
        tick_spacing = int(2000./dz)  
        xticks = axeProfiles.get_xticks()
        yticks =axeProfiles.get_yticks()
        
        axeProfiles.set_xticks(xticks[::tick_spacing])
        axeProfiles.set_yticks(yticks[::tick_spacing])
        axeProfiles.tick_params(axis ='y', 
                                labelsize = 2* self.font_size , 
                                 labelrotation = 90)
        axeProfiles.tick_params(axis ='x', 
                                labelsize = 2* self.font_size , 
                                 )
    
        if show_station_labels: # add profiles labels 
        
            for ii, name in enumerate(profile_lines) : 
     
                axeProfiles.text((self.X[ii][0] +dx)/dz  ,
                            (self.Y[ii][00]+dx)/dz,
                            s= name,
                            horizontalalignment='center',
                            verticalalignment='baseline',
                            fontdict={'size': self.ms*3, 
                                      'color': 'k'},
                            bbox =dict(boxstyle='round',
                                       facecolor ='rosybrown', 
                                       alpha =0.5, pad =.2, 
                                       )
                            )

        
        profile_lines = ['${0}$'.format(''.join([name , 
                                                 ':Azim ={0} °'.format(angle)])) 
                         for name, angle in zip(profile_lines, profile_angles)]
                      
        axeProfiles.legend(leghandles, profile_lines,
                           prop ={'size':2*self.font_size , 
                                'style': self.font_style, 
                                })
        
        print('----> {0} stations plotted successfully.'.
              format(len(snames )))
        
        self.fig.suptitle('Plot {0:02} survey lines'.format(len(self.X)),
                      fontsize= 2. * self.font_size, 
                      verticalalignment='center', 
                      style ='italic',
                      bbox =dict(boxstyle='round',facecolor ='moccasin'), 
                      y=0.95)
        
        
        # plt.tight_layout()
        
        if savefig is not None : plt.savefig(savefig, 
                                             dpi = self.fig_dpi, 
                                                 # orientation =orientation
                                                 )
    
        plt.show()
  
        
class Plot2d (object): 
    """
    class to plot 2D map  
    Deal with all 2D plots
    
    ======================  ===============================================
    keywords                Description
    ======================  ===============================================
    cb_pad                  padding between axes edge and color bar 
    cb_shrink               percentage to shrink the color bar
    climits                 limits of the color scale for resistivity
                            in log scale (min, max)
    cmap                    name of color map for resistivity values
    fig_aspect              aspect ratio between width and height of 
                            resistivity image. 1 for equal axes
    fig_dpi                 resolution of figure in dots-per-inch
    fig_num                 number of figure instance
    fig_size                size of figure in inches (width, height)
    font_size               size of axes tick labels, axes labels is +2
    grid                    [ 'both' | 'major' |'minor' | None ] string 
                            to tell the program to make a grid on the 
                            specified axes.
    ms                      size of station marker 
    plot_yn                 [ 'y' | 'n']
                            'y' --> to plot on instantiation
                            'n' --> to not plot on instantiation
    station_color           color of station marker
    station_font_color      color station label
    station_font_pad        padding between station label and marker
    station_font_rotation   angle of station label in degrees 0 is 
                            horizontal
    station_font_size       font size of station label
    station_font_weight     font weight of station label
    station_id              index to take station label from station name
    station_marker          station marker.  if inputing a LaTex marker
                            be sure to input as r"LaTexMarker" otherwise
                            might not plot properly
    title                   title of plot.  If None then the name of the
                            iteration file and containing folder will be
                            the title with RMS and Roughness.
    xlimits                 limits of plot in x-direction in (km) 
    xminorticks             increment of minor ticks in x direction
    xpad                    padding in x-direction in km
    ylimits                 depth limits of plot positive down (km)
    yminorticks             increment of minor ticks in y-direction
    ypad                    padding in negative y-direction (km)
    yscale                  [ 'km' | 'm' ] scale of plot, if 'm' everything
                            will be scaled accordingly.
    ======================  ===============================================
    
    """
    
    def __init__(self, **kws): 
        
        self._logging= csamtpylog.get_csamtpy_logger(self.__class__.__name__)
        
        self.plot_comp = kws.pop('plot_comp', 'xy')
        
        self.fs =kws.pop('fs', 0.7)
        self.fig_num = kws.pop('fig_num', 1)
        self.fig_size = kws.pop('fig_size', [7,7])
        self.fig_aspect = kws.pop('fig_aspect','auto')
        self.fig_title = kws.pop('fig_title', None)
        
        self.fig_dpi =kws.pop('fig_dpi', 300)
        self.font_size = kws.pop('font_size', 7)
        
        self.aspect = kws.pop('aspect', 'auto')
        self.font_style =kws.pop('font_style', 'italic')
        self.orient=kws.pop('orientation', 'landscape')

        self.plot_style = kws.pop('plot_style', 'imshow')
        self.imshow_interp = kws.pop('imshow_interp', 'bicubic')

        #--> set plot limits
        self.res_limits = kws.pop('res_limits', (0, 4))
        self.phase_limits = kws.pop('phase_limits', (0, 90))
        
        for key in ['fig_title', 'xlimits', 'ylimits']:
            self.__setattr__(key, None )
        
        #--> set colorbar properties
        self.cb_pad = kws.pop('cb_pad', .0375)
        self.cb_orientation = kws.pop('cb_orientation', 'vertical')
        self.cb_shrink = kws.pop('cb_shrink', .75)
        self.cb_position = kws.pop('cb_position', None)
        self.climits = kws.pop('climits', (0, 4))
        self.clabel = kws.pop('clabel', None)

        #--> set text box parameters
        self.text_location = kws.pop('text_location', None)
        self.text_xpad = kws.pop('text_xpad', .95)
        self.text_ypad = kws.pop('text_ypad', .95)
        self.text_size = kws.pop('text_size', self.font_size)
        self.text_weight = kws.pop('text_weight', 'bold')
        self.station_label_rotation = kws.pop('station_label_rotation',45)
        self.show_grid = kws.pop('show_grid',True)
        
        # set makers properties
        self.marker = kws.pop('marker', 'o')
        self.ls =kws.pop('ls', '-')
        self.lc = kws.pop('lc', None)
        self.markeredgecolor= kws.pop('markeredgecolor', 'k')
        self.markerfacecolor=kws.pop('markerfacecolor', 'k')
        self.ms = kws.pop('ms', 2)
        self.lw =kws.pop('lw', 2)
        
        #------ticks parameters ------------------
        self.fw =kws.pop('font_weight', 'bold')
        self.depth_scale=kws.pop('depth_scale', 'm')
        
        # set some figure properties to make plot occam 
        self.station_offsets = None
        self.station_names = kws.pop('station_names', None)
        self.station_id =kws.pop('station_id', None)
  
        self.station_font_size = kws.pop('station_font_size', 8)
        self.station_font_pad = kws.pop('station_font_pad', 1.0)
        self.station_font_weight = kws.pop('station_font_weight', 'bold')
  
        self.station_font_color = kws.pop('station_font_color', 'k')
        self.station_marker = kws.pop('station_marker',
                                         r"$\blacktriangledown$")
        self.station_color = kws.pop('station_color', 'k')

        self.xpad = kws.pop('xpad', 1.0)
        self.ypad = kws.pop('ypad', 1.0)

        self.xminorticks = kws.pop('xminorticks', 5)
        self.yminorticks = kws.pop('yminorticks', 1)

        self.cmap = kws.pop('cmap', 'jet_r')
        
        self.xlabel = kws.pop('xlabel' , None) 
        self.ylabel = kws.pop('ylabel', None) 
        
        self.savefig = kws.pop('savefig', None)
        self.show = kws.pop('show', True) 
        
        for keys in list(kws.keys()): 
            setattr(self, keys, kws[keys])
            
    #XXX TODO     
    def plotphasetensorpseudosection (self, ediObjs, ellip_dict = None, 
                                      stretch = (7000, 20 ), mode ='frequency',
                                      linedir ='ns', **kws): 
        """ Plot phase tensor pseudosection. 
        
        Method plots the phase tensor ellipses in a pseudo section format. Inherits
        from :class:pycsamt.imaging.phase_tensor_pseudosection.PlotPhaseTensorPseudoSection`
        
        Parameters 
        -----------
        ediObjs: list  of  pycsamt.core.edi.Edi or pycsamt.core.edi.Edi objects 
            Collections of EDI-objects from `pyCSAMT`_ and `MTpy`_ packages 

        mode: str ['frequency' | 'period']
              temporal scale of y-axis ('frequency' | 'period').
              Default is ``frequency``. 
              
        stretch : float or tuple (xstretch, ystretch)
            is a factor that scales the distance from one station to the next 
            to make the plot readable. It determines (x,y) aspect ratio of plot.
            *Default* is 200

        linedir: str [ 'ns' | 'ew' ]
            predominant direction of profile line
            * 'ns' -> North-South Line or line is closer to north-south)
            * 'ew' -> East-West line or line is closer to east-west
            *Default* is 'ns'

        station_id: tuple or list
            start and stop of station name indicies. ex: for MT01dr 
            station_id=(0,4) will be MT01.

        ellipse_dict: dictionary
                dictionary of parameters for the phase tensor ellipses with keys:
                * 'size' -> size of ellipse in points
                           *default* is 2

                * 'colorby' : [ 'phimin' | 'phimax' | 'skew' |
                                'skew_seg' | 'phidet' |
                                'ellipticity' ]

                          - 'phimin' -> colors by minimum phase
                          - 'phimax' -> colors by maximum phase
                          - 'skew' -> colors by skew
                          - 'skew_seg' -> colors by skew in
                                         discrete segments
                                         defined by the range
                          - 'normalized_skew' -> colors by skew
                                  see [Booker, 2014]
                          - 'normalized_skew_seg' -> colors by
                                         normalized skew in
                                         discrete segments
                                         defined by the range
                          - 'phidet' -> colors by determinant of
                                       the phase tensor
                          - 'ellipticity' -> colors by ellipticity
                          *default* is 'phimin'

                 * 'range' : tuple (min, max, step)
                       Need to input at least the min and max
                       and if using 'skew_seg' to plot
                       discrete values input step as well
                       *default* depends on 'colorby'

            * 'cmap' : [ 'mt_yl2rd' | 'mt_bl2yl2rd' |
                        'mt_wh2bl' | 'mt_rd2bl' |
                        'mt_bl2wh2rd' | 'mt_seg_bl2wh2rd' |
                        'mt_rd2gr2bl' ]

                     - 'mt_yl2rd' -> yellow to red
                     - 'mt_bl2yl2rd' -> blue to yellow to red
                     - 'mt_wh2bl' -> white to blue
                     - 'mt_rd2bl' -> red to blue
                     - 'mt_bl2wh2rd' -> blue to white to red
                     - 'mt_bl2gr2rd' -> blue to green to red
                     - 'mt_rd2gr2bl' -> red to green to blue
                     - 'mt_seg_bl2wh2rd' -> discrete blue to
                                           white to red
        kws: dict 
            Additional keywords arguments from :class:`~.PlotPhaseTensorPseudoSection` 
            of `MTpy`_ package. 
                              
        Returns 
        --------
        ptsection: `~.PlotPhaseTensorPseudoSection` object 
        
        Examples
        ---------
        >>> from pycsamt.utils import get_ediObjs 
        >>> ediObjs = get_ediObjs('data/edis') 
        >>> p2d= Plot2d (fig_size =(7, 7))
        >>> p2d.plotphasetensorpseudosection(ediObjs) 
        
        """
        if not imtpy: 
            raise ImportError ("Can't plot Phase Tensor Pseudosection. Please "
                               " install `mtpy` manually. ")
             
        ediObjs = func.get_ediObjs(ediObjs)
        edi_list= np.array (list (map (lambda o: o.edifile , ediObjs)))
        ellipse_dict = ellip_dict or  {
            'ellipse_colorby':'phimin', # phimin, phimax, skew,
            # skew_seg
            'ellipse_range':[0,90], # Color limits. If plotting
            'ellip_size': 2, 
            'ellipse_cmap':'mt_bl2wh2rd'
        } 
        # skew_seg need to provide
        # 3 numbers, the 3rd indicates
        # interval, e.g. [-12,12,3]
        # create a plot object
        ptsection = PlotPhaseTensorPseudoSection(
                        fn_list = edi_list,
                        fig_size = self.fig_size, 
                        tscale = mode, #'frequency', 
                        plot_num = self.fig_num, 
                        plot_title = self.fig_title, 
                        xlimits = self.xlimits, 
                        ylimits = self.ylimits,
                        linedir= linedir ,  
                        #(25, 10), # (17,8) determines (x,y) aspect ratio of plot
                        stretch= stretch , #(7000, 20 ), for fig = (5, 5)#20 000, 10
                        station_id=(0,len(ediObjs)), 
                        font_size=self.font_size ,
                        lw=self.lw,
                        **ellipse_dict , 
                        **kws,
        )
        if self.show :
            ptsection.plot()
        
        if self.savefig is not None: 
            ptsection.save_figure(save_fn =self.savefig, 
                                  fig_dpi=self.fig_dpi)
            ptsection.close() 
            
        return ptsection 
    
        
    def plottemp2d (self, arr2d, y=None,  x =None, style ='pcolormesh', 
                distance = 50., stnlist =None, prefix ='S', how= 'py',
                **kws): 
        """ 2D template for quick visualization. 
        
        Parameters 
        -----------
        arr2d : ndarray , shape (N, M) 
            2D array for plotting. For instance, it can be a 2D resistivity 
            collected at all stations (N) and all frequency (M) 
        y: array-like 
            Y-coordinates. It should have the length N, the same of the ``arr2d``.
            the rows of the ``arr2d``.
        x: array-like 
            X-coordinates. It should have the length M, the same of the ``arr2d``; 
            the columns of the 2D dimensional array.  Note that if `x` is 
            given, the `distance is not needed. 
            
        style: str 
            matplotlib plot style.  It could be ``imshow`` or ``pcolormesh``. The 
            default is ``pcolormesh``. 
            
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
            
        kws : dict 
            Additional keywords arguments of Matplotlib subsplots function  
            :func:`plt.subplots`
            
        Examples 
        -------- 
        >>> import numpy as np
        >>> from pycsamt.processing import make2d , get_ediObjs, interpolate2d 
        >>> from pycsamt.processing import get_full_frequency 
        >>> from pycsamt.view import Plot2d 
        >>> # create a 2d grid of resistivity data from edipath 
        >>> edipath = 'data/3edis'
        >>> ediObjs = get_ediObjs (edipath)
        >>> res2d_xy = make2d (ediObjs, 'resxy')
        >>> # can interpolate to fix the missing data 
        >>> # get the full frequency of survey area in the whole collection data 
        >>> freqs = get_full_frequency(ediObjs) # plot station vs frequency 
        >>> res2d_xy = interpolate2d (res2d_xy )
        >>> # create a 2D plot objects 
        >>> p2d = Plot2d(figsize =(6, 3), xlabel= '$Distance(m)$',
                       ylabel = '$Log_{10}Frequency [Hz]$', 
                       clabel = '$Log_{10}Rhoa[\Omega.m$')
        >>> # we want to plotlog10 and log10 resistivity 
        >>> res2d_xy = np.log10(res2d_xy); freqs = np.log10(freqs)
        >>> p2d.plottemp2d (res2d_xy, y = freqs, distance = 50 ) # distance =50 m
        
        """
        
        fig, axe = plt.subplots(1, figsize = self.fig_size, num = self.fig_num,
                                dpi = self.fig_dpi , **kws 
                    )
         
        style = style or 'pcolormesh' 
        sc = copy.deepcopy(style)
        style = str (style).lower().strip() 
        if style not in ('pcolormesh', 'imshow'): 
            raise ValueError("Plot can be either 'pcolormesh' or "
                             f"'imshow' not {sc!r} ")
        
        try : 
            distance = float(distance) 
        except : 
            raise TypeError (f'Expect a float value not {type(distance).__name__!r}')
            
        if y is not None: 
            if len(y) != arr2d.shape [0]: 
                raise ValueError (" 'y' array must have an identical number " 
                                  f" of row of 2D array: {arr2d.shape[0]}")
                
        if x is not None: 
            if len(x) != arr2d.shape[1]: 
                raise ValueError (" 'x' array must have the same number " 
                                  f" of columns of 2D array: {arr2d.shape[1]}")

        d= distance or 1. 
        y = np.arange(arr2d.shape [0]) if y is None else y 
        x= x  or np.arange(arr2d.shape[1])  * d 
        if stnlist is not None: 
            stnlist = func._assert_all_types(stnlist, tuple, list, np.ndarray)
            if len(stnlist)!= len(x):
                raise ValueError("Expect  the length of the list of station names"
                                 f" to be the same like to position {len(x)}")
                
        stn = stnlist or func.make_ids ( x , prefix , how = how) 
        
        cmap = plt.get_cmap( self.cmap)
        
        if style =='pcolormesh': 
            
            X, Y = np.meshgrid (x, y)
            axr = axe.pcolormesh ( X, Y, arr2d,
                            # for consistency check whether array does not 
                            # contain any NaN values 
                            vmax = arr2d[ ~np.isnan(arr2d)].max(), 
                            vmin = arr2d[ ~np.isnan(arr2d)].min(), 
                            shading= 'gouraud', 
                            cmap =cmap, 
                                  )

        if style =='imshow': 

            axr= axe.imshow (arr2d,
                            interpolation = self.imshow_interp, 
                            cmap =cmap,
                            aspect = self.fig_aspect ,
                            origin= 'upper', 
                            extent=( x[~np.isnan(x)].min(),
                                      x[~np.isnan(x)].max(), 
                                      y[~np.isnan(y)].min(), 
                                      y[~np.isnan(y)].max())
                                              )
            axe.set_ylim(y[~np.isnan(y)].min(), y[~np.isnan(y)].max())
        

        axe.set_xlabel(self.xlabel or 'Distance(m)', 
                     fontdict ={
   
                      'size': 1.5* self.font_size ,
                      'weight': self.fw}
                                            )
      
        axe.set_ylabel(self.ylabel or 'log10(Frequency)[Hz]',
                 fontdict ={
                         #'style': self.font_style, 
                                  'size': 1.5* self.font_size ,
                                  'weight': self.fw})
        if self.show_grid is True : 
            axe.minorticks_on()
            axe.grid(color='k', ls=':', lw =0.25, alpha=0.7, 
                         which ='major')
 
        labex , cf = self.clabel or '$log10(App.Res)[Ω.m]$', axr

        cb = fig.colorbar(cf , ax= axe)
        cb.ax.yaxis.tick_left()
        cb.ax.tick_params(axis='y', direction='in', pad=2.)
        
        cb.set_label(labex,fontdict={'size': 1.5* self.font_size ,
                                  'style':self.font_style})
        #--> set second axis 
        axe2 = axe.twiny() 
        axe2.set_xticks(ticks= x,
                    minor=False )
        axe2.set_xticklabels(stn, 
                         rotation=self.station_label_rotation)
     
        axe2.set_xlabel('Stations', 
                        fontdict ={'style': self.font_style, 
                                   'size': 1.5* self.font_size ,
                                   'weight': self.fw}, )
      

        fig.suptitle(self.fig_title,
                     ha='left',
                     fontsize= 15* self.fs, 
                     verticalalignment='center', 
                    style =self.font_style,
                    bbox =dict(boxstyle='round',
                               facecolor ='moccasin'))
   
        #plt.tight_layout(h_pad =1.8, w_pad =2*1.08)

        plt.tight_layout()  
        
        if self.savefig is not None :
            fig.savefig(self.savefig, dpi = self.fig_dpi,
                        orientation =self.orient)
            #plt.close(fig =fig ) 
        plt.show() if self.savefig is None else plt.close(fig=fig) 
        
        
        
        
    def penetration2D(self, fn=None , profile_fn =None, 
                      savefig =None, doi= '2km',  **kwargs): 
        """
        Plot penetration 2D. 
        
        Parameters
        -----------
            * fn : str 
                full path to [EDI|AVG|J] files. 
                
            * doi : float 
                depth assumed to be imaged , default is 2000m For  CSAMT ,
                2km is enought to have more info about near surface.
                 * Default* unit is "m".
                 
            * profile_fn : str 
                full path to profile *stn file 
                
            * savefig : str 
                outdir 
                
         Returns
         --------
             obj ,
                 plot penetration obj. 
             
        ==============  ===================  ==================================
        Params          Type                    Description 
        ==============  ===================  ==================================
        plot_style      str                     pcolormesh or imshow 
                                                Default is **pcolormesh**
        ms              int                     markersize *Default* is .9*fs
        cm              str                     mpl.colormap ,
                                                *Default* is"Purples".
        rename_station  list                    can set new_stationname. 
        ==============  ===================  ==================================
        
        :Example:
            
            >>> path =  os.path.join(os.environ["pyCSAMT"], 
            ...      'pycsamt','data', K1.AVG)
            >>> plot2d_obj = plot2d()
            >>> plot2d_obj.penetration2D(fn = path, 
            ...                        profile_fn=os.path.join(
            ...                            os.path.dirname(path),'K1.stn'), 
            ...                         plot_style='imshow',  doi='10000m')
        """
        
        self._logging.info ('Contructing 2D penetration depth')
        
        stnnames = kwargs.pop('rename_station', None )
        plot_style = kwargs.pop('plot_style', None)
        alpha= kwargs.pop('alpha', 0.5)
        mplcmap =kwargs.pop('cm', 'twilight')
        
        #-----define csmat_obj_ 
        csamt_obj = CSAMT(data_fn= fn , profile_fn=profile_fn )
        csamt_dep1D_obj = csamt_obj.skindepth 
        csamt_freq_obj =csamt_obj.freq 
        csamt_stnDis_obj =csamt_obj.station_distance
        csamt_stn_obj = sorted(list(csamt_obj.skindepth .keys()))      

            
        if stnnames is not None : 
            assert len(stnnames ) ==len(csamt_stn_obj),\
                CSex.pyCSAMTError_station('Station provided must have the '
                                          'same lenght with default stations.')
            csamt_stn_obj = stnnames 
            

        def depth2D (dep_loc,  freq_obj, doi ): 
            """
            Build matrix freq at n-station from 1D skin depth and return 
            use station. 
            
            :param dep_loc: obj to get create matrix  from all station 
            :type dep_loc: dict 
            
            :param freq_obj:  frequency array 
            :type freq_obj:  array_like 
            
            :param doi: limit of investigation depth .
            :type doi: int  
            
            :returns: cutoff matrix fonction to depth investigation. 
            :rtype: ndarray
            """
            dep2D_GRID_obj  =  func.concat_array_from_list(
                list_of_array= [values for keys ,
                                values  in dep_loc.items()], 
                                 concat_axis=1)
             #flip frequency if it is ranged lower to higher
            if freq_obj[0] < csamt_freq_obj[-1]: 
                freq_obj ==freq_obj[::-1]
                dep2D_GRID_obj = dep2D_GRID_obj[::-1]
 
            return mplotus.slice_matrix(base_matrix = dep2D_GRID_obj , 
                                        freq_array=freq_obj, doi=doi)
    
        # select the convenient plages of data  to explore
        dep2D_GRID_obj , csamt_freq_obj,  doi = depth2D(dep_loc =csamt_dep1D_obj,
                                                        doi = doi, 
                                                        freq_obj=csamt_freq_obj)
        

        #---------------------plotfeatures ---------------
        if plot_style is None : plot_style = 'pcolormesh'
        
        mpl.rcParams['figure.figsize']=[12,4]
        mpl.rcParams["image.origin"] ='upper' 
        
        cmap = plt.get_cmap(mplcmap)
        fig, axs =plt.subplots()
        
        #--- create meshgrig obj plot with pcolor mesh  ------------------------
        if plot_style =='pcolormesh': 
            ydepth =np.linspace(dep2D_GRID_obj.min(),
                                dep2D_GRID_obj.max(), 
                                dep2D_GRID_obj.shape[0])
            # ydepth = -1* ydepth 
              
            
            xxDepth,  yyDepth = np.meshgrid(csamt_stnDis_obj , ydepth )
            cf = axs.pcolormesh(xxDepth,
                                yyDepth,
                                dep2D_GRID_obj ,  
                                vmax=dep2D_GRID_obj.max(),
                                vmin= dep2D_GRID_obj.min(),
                                shading ='gouraud', 
                                cmap= cmap,


                                )
            axs.set_ylim(dep2D_GRID_obj.max(), dep2D_GRID_obj.min())
            #ressetting yticks 
            axs.set_yticks(mplotus.resetting_ticks(get_xyticks=axs.get_yticks()))


            # cb.set_label('$Depth\ in\ {0}eters$'.format(depth_scale),
            #               fontdict={'size': 2* self.font_size})
            
            axs.minorticks_on()
            axs.grid(color='k', ls=':', lw =0.25, alpha=alpha, which ='major')
            axs.set_xlabel('Distance(m)',  fontdict ={'style': self.font_style, 
                                                      'size': 2* self.font_size ,
                                                      'weight': self.fw}
                                                                            )
            axs.set_ylabel('Penetration depth({0})'.format('m'),
                           fontdict ={'style': self.font_style, 
                                                    'size': 2* self.font_size ,
                                                    'weight': self.fw})
            # axs.set_yscale ('log',  nonposy='clip')
            

        if plot_style =='imshow': 
            cf = axs.imshow( dep2D_GRID_obj ,
                                cmap=cmap,
                                vmax=dep2D_GRID_obj.max(),
                                vmin= dep2D_GRID_obj.min(),
                                interpolation= self.imshow_interp,
                                aspect =self.aspect, 
                                origin='upper', 
                                extent=(csamt_stnDis_obj.min(),
                                        csamt_stnDis_obj.max(), 
                                        dep2D_GRID_obj.max(), 
                                        dep2D_GRID_obj.min()),
                                )

            #               fontdict={'size': 2* self.font_size,})
            axs.minorticks_on()
            axs.grid(color='k', ls=':', lw =0.25, alpha=alpha, which ='major')
            axs.set_xlabel('Distance(m)',  fontdict ={'style': self.font_style, 
                                                      'size': 2* self.font_size ,
                                                      'weight': self.fw}
                                                                            )
            axs.set_ylabel('Penetration depth({0})'.format('m'),
                           fontdict ={'style': self.font_style, 
                                                    'size': 2* self.font_size ,
                                                    'weight': self.fw})
        #set color bar properties 

        cbarmin = dep2D_GRID_obj.min()
        cbarmax = dep2D_GRID_obj.max()
        
        cbound = mplotus.resetting_colorbar_bound(cbmax =cbarmax, cbmin =cbarmin)
        cb = fig.colorbar(cf , ax= axs)
        # mplcb.ColorbarBase(ax=axs, cmap=cmap , 
        #                     norm =mpl.colors.Normalize(vmin=cbarmin, 
        #                                                 vmax =cbarmax ), 
        #                               orientation =self.cb_orientation, 
        #                               spacing ='uniform', 
        #                         )
        cb.set_ticks(cbound)
        cb.set_ticklabels(['{0}'.format(int(ff)) for ff in cbound], )
        cb.ax.yaxis.set_label_position('right')
        cb.ax.yaxis.set_label_coords(2.,.5)
        cb.ax.yaxis.tick_left()
        cb.ax.tick_params(axis='y', direction='in', pad=3)
        
        cb.set_label('$Depth\ in\ {0}eters$'.format('m'),
                          fontdict={'size': 2* self.font_size,})
        
        
        #--> set second axis 
        axe2 = axs.twiny() 
        axe2.set_xticks(ticks= csamt_stnDis_obj, minor=False  )
        axe2.set_xticklabels(csamt_stn_obj , rotation=self.station_label_rotation)
        axe2.set_xlabel('stations', fontdict ={'style': self.font_style, 
                                                'size': 2* self.font_size ,
                                                'weight': self.fw}, )
        
        fig.tight_layout()

        fig.suptitle('Plot Penetration depth', ha='left',
                     fontsize= 15* self.fs, 
                     verticalalignment='center', 
                             style =self.font_style,
                             bbox =dict(boxstyle='round',facecolor ='moccasin'))
            
        plt.tight_layout(h_pad =1.8, w_pad =2*1.08)
                
        if savefig is not None : plt.savefig(savefig, 
                                             dpi = self.fig_dpi,
                                             orientation =self.orient)
        
        
        plt.show()
        
    def pseudocrossResPhase(self, fn , profile_fn =None ,
                            savefig =None , plot_style=None, **kws ):
        """
        Plot Pseudocrossection of resistivity and phase.
        
        :param fn: full path to ['AVG', 'EDI', 'J'] file .
        :type fn: str 
        
        :param profile_fn:  full path to profile station file  in the case fn 
                            is *AVG. 
        :type profile_fn: str
        
        :param savefig: path to save figure 
        :type savefig: str 
                
        """

        self._logging.info ('Construction of PseudocSection of Resistivity'\
                            ' and Phase from <{0}>'.format(os.path.basename(fn)))
        
        delineate_resistivity_curve =kws.pop('delineate_resistivity', None)
        #tolerance_value =kws.pop('atol', 0.2)
        delineate_phase_curve = kws.pop('delineate_phase', None)
        contourlines =kws.pop('contour_lines_style', '-')
        contourcolors =kws.pop('contour_lines_color', 'white')

        #--create obj ----
        csamt_obj =CSAMT(data_fn=fn , profile_fn=profile_fn, 
                         component =self.plot_comp)
        csamt_phase_obj =csamt_obj.phase 
        csamt_res_obj =csamt_obj.resistivity 
        csamt_freq_obj =csamt_obj.freq
        csamt_stnDis_obj =csamt_obj.station_distance  #*1e3 if \
            #distance =='km' else csamt_obj.station_distance
            
        csamt_stn_obj = sorted(csamt_obj.station)
        if str(self.station_names).lower() =='edi': 
            try : 
                csamt_stn_obj = [
                    e.replace('.edi', '') for e in csamt_obj.edinames]
            except : 
                pass
         
        ### Ignore the fact that arrays should have all the same 
        # length and plot anyway. 
        csamt_RES_obj = func.resize_resphase_values(
            csamt_res_obj,c =csamt_freq_obj )
        csamt_PHS_obj = func.resize_resphase_values(csamt_phase_obj,
                                               c = csamt_freq_obj)
        #convert Res and Phase values onto logarithm scale .
        csamt_RES_obj  = np.log10( csamt_RES_obj)
        csamt_freq_obj = np.log10 (csamt_freq_obj ) 
        
        #--> get delineate curve , if exist .
        if delineate_phase_curve is not None : 
            delineate_phase_curve  = [ss%90 for ss in \
                                      mplotus.control_delineate_curve(
                                          phase_deline=delineate_phase_curve) ] 
        if delineate_resistivity_curve is not None : 
            delineate_resistivity_curve = mplotus.control_delineate_curve(
                res_deline=delineate_resistivity_curve)
                                                      
        #-----------------------PLOT ---------------------------------------
        if plot_style is None : 
            plot_style = 'pcolormesh'
        
        #--------------------figure params -----------
        
        mpl.rcParams['figure.figsize']=self.fig_size 
        fig =plt.figure()
        axe_res =plt.subplot2grid(shape=(2,1), loc=(0,0))
        axe_phase =plt.subplot2grid(shape=(2,1), loc=(1,0))
        
        cmap = plt.get_cmap( self.cmap)
 
        if plot_style.lower() =='pcolormesh': 
            yres_matrix, xres_matrix =np.meshgrid(
                csamt_freq_obj, csamt_stnDis_obj ) 
             
            yphase_matrix , xphase_matrix = np.meshgrid(
                 csamt_freq_obj, csamt_stnDis_obj) 
            #---res map 
            app_rho_axe = axe_res.pcolormesh (xres_matrix,
                                          yres_matrix ,
                                          csamt_RES_obj,
                                        vmax = csamt_RES_obj.max(), 
                                        vmin = csamt_RES_obj.min(), 
                                        shading= 'gouraud', 
                                        cmap =cmap, 
                                              )
            #---phase map 
            phase_axe = axe_phase.pcolormesh (xphase_matrix ,
                                              yphase_matrix ,
                                              csamt_PHS_obj, 
                                               vmax = csamt_PHS_obj.max(), 
                                               vmin = csamt_PHS_obj.min(), 
                                              shading= 'gouraud', 
                                              cmap =cmap, 
                                              )
            
            MAT = [[xres_matrix, yres_matrix ,csamt_RES_obj], 
                   [xphase_matrix , yphase_matrix ,csamt_PHS_obj]]
            for ii, (axe, deline)  in enumerate( zip ([axe_res, axe_phase], 
                                      [delineate_resistivity_curve,
                                       delineate_phase_curve])) : 
                if  deline is not None : 
                    contps = axe.contour(*MAT[ii], 
                                         colors =contourcolors, 
                                         linestyles=contourlines)
                    try :
                        axe.clabel(contps, deline ,
                                        inline=True, fmt='%1.1f',
                                        fontsize =self.font_size,
                                                  )
                    except:
                        # deline =None and set all contours
                        if ii==0 : mesf = 'resistivity'
                        else : mesf ='phase'
                    
                        warnings.warn(
                            'Values {0} given as {1} contours levels does not'
                            ' match ! Contours levels are resseting '
                            'to default levels !'.format(deline, mesf))
                        
                        print('---> Values {0} given can not be set as'
                              '  {1} contours levels.'
                              ' Default levels are {2}.'.format(deline, mesf,
                                                                contps.levels))
    
                        print('.--> ! {0} contours levels = {1} are resseting to '
                              ' default levels!'.format(mesf.capitalize(),
                                                        deline))
                        
                        self._logging.debug (
                            'values {0} given as contours levels does not match ! '
                              'availables contours levels are set '
                              'to default values.')
                        
                        axe.clabel(contps, 
                                    inline=True,
                                    fmt='%1.1f',
                                    fontsize =self.font_size,
                                                  )
                    
        if plot_style.lower() =='imshow': 
      
            #---res map 
            app_rho_axe = axe_res.imshow (csamt_RES_obj.T,
                                          # vmax = csamt_RES_obj.max(), 
                                          # vmin = csamt_RES_obj.min(), 
                                          interpolation = self.imshow_interp, 
                                          cmap =cmap,
                                          aspect = self.fig_aspect ,
                                          origin= 'upper', 
                                          extent=(csamt_stnDis_obj.min(),
                                                    csamt_stnDis_obj.max(), 
                                                    csamt_freq_obj.min(), 
                                                    csamt_freq_obj.max())
                                              )
            axe_res.set_ylim(csamt_freq_obj.min(), csamt_freq_obj.max())
 
            #---phase map 
            phase_axe = axe_phase.imshow ( np.flipud(csamt_PHS_obj.T), 
                                              # vmax = csamt_PHS_obj.max(), 
                                              # vmin = csamt_PHS_obj.min(), 
                                              interpolation = self.imshow_interp, 
                                              aspect =self.fig_aspect ,
                                              cmap =cmap,
                                              origin= 'lower', 
                                               extent=(csamt_stnDis_obj.min(),
                                                       csamt_stnDis_obj.max(), 
                                                       csamt_freq_obj.min(), 
                                                        csamt_freq_obj.max(), 
                                                        
                                                       ),
                                              )
   
            MAT = [csamt_RES_obj, csamt_PHS_obj]
            for ii, (axe, deline)  in enumerate( zip ([axe_res, axe_phase], 
                                      [delineate_resistivity_curve,
                                       delineate_phase_curve])) :
                if  deline is not None :
                    if ii ==0 : origin ='upper'
                    else : origin ='lower'
                    contps = axe.contour(MAT[ii], colors =contourcolors, 
                                          vmax=MAT[ii].max(),
                                          vmin = MAT[ii].min(), 
                                            linestyles=contourlines, 
                                            extent =(csamt_stnDis_obj.min(),
                                                    csamt_stnDis_obj.max(), 
                                                    csamt_freq_obj.min(), 
                                                   csamt_freq_obj.max()),
                                            extend ='both', origin= origin ,  
                                           )
                    axe.clabel(contps, deline ,
                                    inline=True, 
                                    fmt='%1.1f',
                                    fontsize =self.font_size,
                                              )
                    axe.set_ylim (csamt_freq_obj.min(), 
                                  csamt_freq_obj.max())
            
        
        #for twin axes 
        for ii, ax in enumerate([axe_res, axe_phase]):
            
             if ii ==1 :
                 ax.set_xlabel('Distance(m)', 
                               fontdict ={
                                   #'style': self.font_style, 
                                'size': 1.5* self.font_size ,
                                'weight': self.fw}
                                                      )
             #if plot_style =='pcolormesh': 
             ax.set_ylabel('log10(Frequency)[Hz]',
                          fontdict ={
                                  #'style': self.font_style, 
                                           'size': 1.5* self.font_size ,
                                           'weight': self.fw})
             if self.show_grid is True : 
                 axe_res.minorticks_on()
                 axe_res.grid(color='k', ls=':', lw =0.25, alpha=0.7, 
                              which ='major')
                 axe_phase.minorticks_on()
                 axe_phase.grid(color='k', ls=':', lw =0.25, alpha=0.7, 
                                which ='major')
            # congigure color bar 
             if ii == 0 : 
  
                 labex , cf = '$log10(App.Res)[Ω.m]$', app_rho_axe
             elif ii == 1 :
                 labex, cf = '$Phase(deg)$', phase_axe
                 
             cb = fig.colorbar(cf , ax= ax)
             cb.ax.yaxis.tick_left()
             cb.ax.tick_params(axis='y', direction='in', pad=2.)
             
             cb.set_label(labex,fontdict={'size': 1.5* self.font_size ,
                                          'style':self.font_style})
                              
            #--> set second axis 
             axe2 = ax.twiny() 
             axe2.set_xticks(ticks= csamt_stnDis_obj,
                             minor=False )
             axe2.set_xticklabels(csamt_stn_obj , 
                                  rotation=self.station_label_rotation)
             if ii==0:
                 axe2.set_xlabel('Stations', 
                                 fontdict ={'style': self.font_style, 
                                            'size': 1.5* self.font_size ,
                                            'weight': self.fw}, )
             if ii != 0:
                 plt.setp(axe2.get_xticklabels(), visible=False)
                 plt.setp(axe2.get_xlabel(), visible=False)

        fig.tight_layout()

        fig.suptitle(r'PseudocrossResistivity $\rho {0}$ '
                     'and Phase $\phi {0}$'.format(csamt_obj.component),
                     ha='left',
                     fontsize= 15* self.fs, 
                     verticalalignment='center', 
                             style =self.font_style,
                             bbox =dict(boxstyle='round',
                                        facecolor ='moccasin'))
            
        plt.tight_layout(h_pad =1.8, w_pad =2*1.08)
                
        if savefig is not None : plt.savefig(savefig, dpi = self.fig_dpi,
                                             orientation =self.orient)
        
        
        plt.show()

    def plot_occam2dModel(self, model_fn =None, iter_fn=None , 
                          mesh_fn =None , data_fn =None, doi=1000, 
                          **kwargs ):
        
        """
        Plotoccam Model  form Occam Model class 
                   
        :param model_fn: full path to Occam 2Dmodel file 
        :type model_fn: str  
        
        ==============  =========  =======================================
        Params          Type       Description 
        ==============  =========  =======================================
        iter_fn         str         full path to occam iteration file 
        mesh_fn         str         full path to mesh_fn file 
        data_fn         str         full path to occam_data file 
        doi             str         depth of investigation might 
                                    be float or str like "1km" =1000
        depth_scale     str         scale of imaging depth can be 
                                    "km" or "m". *Default* is"m"
        ==============  =========  =======================================
        
        :Example:
            
            >>> data='OccamDataFile.dat'
            >>> mesh = 'Occam2DMesh'
            >>> model = 'Occam2DModel'
            >>> iter_='ITER17.iter'
            >>> path =os.path.join(os.environ ['pyCSAMT'], 
            ...                       'data', 'occam2D', mesh)
            >>> plot2d_obj = plot2d()
            >>> plot2d_obj.plot_occam2dModel(mesh_fn=path, 
            ...                            iter_fn = os.path.join(
            ...                                os.path.dirname(path), iter_), 
            ...                            model_fn =os.path.join(
            ...                                os.path.dirname(path), model) , 
            ...                            data_fn =os.path.join(
            ...                                os.path.dirname(path), data ),
            ...                             doi='1km')
        """
        self._logging.info('Plot occamModel 2D')
        depth_scale =kwargs.pop('depth_scale', 'm')
        
        savefig =kwargs.pop('savefig', None)
        change_station_id =kwargs.pop('new_station_names', None)

        plot_style =kwargs.pop('plot_style', None )
        if plot_style is None : plot_style = 'pcolormesh'
  
        show_contour =kwargs.pop('show_contour', False)
        contourlines =kwargs.pop('contour_lines_styles', '-')
        contourcolors =kwargs.pop('contour_lines_colors', 'white')
        delineate_resistivity_curve =kwargs.pop('delineate_rho', None)
        grid_alpha =kwargs.pop('alpha', 0.5)
             
        set_station_label=kwargs.pop('show_station_id', True)
        
        for file , label in zip ( [iter_fn, mesh_fn , data_fn ], 
                                 ['iteration', 'mesh', 'data']): 
            if file is None : 
                mess= 'No {0}-file found !Please input your {0} file.'.format(label)
                warnings.warn(
                    'Iteration, mesh , data files are essential for plotting.'+ mess)
                self._logging.error 
                raise CSex.pyCSAMTError_occam2d_plot(mess) 
                
        # scaled data and x values plots 
        if depth_scale is not None : self.depth_scale= str(depth_scale).lower() 
        if depth_scale not in ["km", "m"]: 
            self.depth_scale= "m"
            mess ="--> ! Depth scale provided ={} is '\
                unrecognized. We reset to 'm'.".format(self.depth_scale)
            warnings.warn(mess)
            self._logging.debug (mess)
        
        if self.depth_scale == 'km':
            dz  = 1000.

        elif self.depth_scale == 'm': # for CSAMT , we use default as meter"m".
            dz = 1.

        #---> get delineate rho curve --- 
        if delineate_resistivity_curve is not None : 
            delineate_resistivity_curve = mplotus.control_delineate_curve(
                res_deline=delineate_resistivity_curve)
            #assert value to put on float rounded to 1 
            # note that value of resistivity for delineate MUST be
            #on OHM- M not log10 resistivity 
            try : 
                # for consistency 
                if not isinstance (delineate_resistivity_curve, list ): 
                    delineate_resistivity_curve=[delineate_resistivity_curve ]
            except : 
                # be sure to stay on the same output value if something wrong happen 
                
                delineate_resistivity_curve=delineate_resistivity_curve
                pass 
                
        # -------------FIND OBJECTS ----------------------------
        # Read occam 2d model object 
        occam_model_obj = occam2d.Model(iter_fn=iter_fn , 
                                        model_fn = model_fn , mesh_fn =mesh_fn)
        occam_data_obj = occam2d.Data(data_fn = data_fn)
        
        # get station names and station offsetS objs 
        
        occam_data_station_offsets =np.array(occam_data_obj.data_offsets)
        # generally data from occam_station offset are 
        #normalized then get the dipole length to normalize 
        # xpad 
        dl = occam_data_station_offsets.max()/ (len(occam_data_station_offsets)-1)
        self.xpad = (dl/2)/dz 

        occam_data_station_names =occam_data_obj.data_sites
    
        if change_station_id is not None :  # call fonction to build a nu
            occam_data_station_names , mess= mplotus.build_new_station_id(
                station_id =occam_data_station_names ,
                new_station_name =change_station_id )
            self._logging.debug(mess)

                                                
        
        # --> get plot objects for model class 
        plot_x_axis  =  occam_model_obj.model_station_offsets
        plot_z_axis  =  occam_model_obj.model_depth_offsets
        occam_model_resistiviy_obj = occam_model_obj.model_resistivity
        
        #------END Objects statements ------------
        
        #-----handles depth of investigation or depth of imaging {doi}:-----
            
        # --> check doi value provided , and convert to default unit {meters}  
        doi =mplotus.depth_of_investigation(doi=doi)
        
        # set boundaries of stations offsets and depth 
        spec_f = -(doi/5)/dz  
                             
        #-25 +0 : 1300 +25 (xpad = 25 )                              
        self.xlimits=(occam_data_station_offsets.min()/dz -self.xpad  , 
                      occam_data_station_offsets.max()/dz + self.xpad )

        
        #then new_sation offset becomes 
        self.ylimits =(spec_f, doi/dz) 

        # plot ---------------figure and properties  ---------------------

        subplot_right = .99
        subplot_left = .085
        subplot_top = .92
        subplot_bottom = .1
        
        mpl.rcParams['figure.figsize']=[12,6]
        plt.rcParams['figure.subplot.left'] = subplot_left
        plt.rcParams['figure.subplot.right'] = subplot_right
        plt.rcParams['figure.subplot.bottom'] = subplot_bottom
        plt.rcParams['figure.subplot.top'] =subplot_top
        
        
        self.fig = plt.figure(self.fig_num, self.fig_size, dpi=self.fig_dpi)
        # plt.clf()
        
        # add a subplot to the figure with the specified aspect ratio
        self.fig_aspect ='auto'
        axm = self.fig.add_subplot(1, 1, 1, aspect=self.fig_aspect)
        
        #-----PLOTS STATEMENTS -----------------------------
        
                #fist option is "pcolormesh " 
        self._logging.info (
            'Ready to plot Model with matplotlib "{0}" style.'.format(plot_style))     
        if plot_style =='pcolormesh':
            
            # if you keep plot_x_axis and plot_z_axis in meter ,
            #  be sure to divided py dz meshes respectively
            #-----------------------------------
            plot_x_axis = np.linspace(occam_data_station_offsets.min(),
                                      occam_data_station_offsets.max(), 
                                      occam_model_resistiviy_obj.shape[1])
            
            #print(occam_model_resistiviy_obj.shape)
            mesh_x  , mesh_z= np.meshgrid(plot_x_axis,  plot_z_axis )
            axm.pcolormesh (mesh_x/dz  , 
                            mesh_z/dz ,
                            occam_model_resistiviy_obj,
                            # vmin =occam_data_station_offsets.min(),  #self.climits[0],
                            # vmax = occam_data_station_offsets.max(), #self.climits[1],  
                            shading= 'auto', 
                            cmap =self.cmap, 
                            alpha = None, 
                                      
                                      )

            if show_contour is True :
                contps = axm.contour(mesh_x/dz  , mesh_z /dz ,
                                     occam_model_resistiviy_obj,
                                      colors =contourcolors,
                                      linestyles=contourlines)
                if  delineate_resistivity_curve is not None : 
                    axm.clabel(contps,  delineate_resistivity_curve ,
                                    inline=True, fmt='%1.1f',
                                    fontsize =self.font_size,
                                      )
            axm.set_xlim( [int(plot_x_axis.min()),  int(plot_x_axis.max())])
            axm.set_ylim ([self.ylimits[1], self.ylimits[0]]) 
            
        if plot_style.lower() =='imshow': 

            #mesh_x  , mesh_z= np.meshgrid(plot_x_axis  , plot_z_axis )
            # to get the origine =0 of the plot
            axm.imshow (occam_model_resistiviy_obj,
                        vmax = self.climits[1], 
                        vmin =self.climits[0], 
                        interpolation = self.imshow_interp, 
                        cmap =self.cmap,
                        aspect = self.fig_aspect,
                        origin= 'upper', 
                        extent=( self.xlimits[0],
                                self.xlimits[1],
                                self.ylimits[1], 
                                self.ylimits[0] - spec_f),  

                                    )

 
            if delineate_resistivity_curve is not None :
                origin ='upper'
                contps = axm.contour(occam_model_resistiviy_obj,
                                     colors =contourcolors, 
                                      vmax=self.climits[0],
                                      vmin = self.climits[1], 
                                    linestyles=contourlines, 
                                    extent =( self.xlimits[0],
                                              self.xlimits[1],
                                              self.ylimits[1], 
                                              self.ylimits[0]),
                                    extend ='both',
                                    origin= origin ,  
                                       )
                axm.clabel(contps, delineate_resistivity_curve ,
                                inline=True, 
                                fmt='%1.1f',
                                fontsize =self.font_size,
                                          )
        
            axm.set_xlim( [self.xlimits[0],  self.xlimits[1]])
            axm.set_ylim ([self.ylimits[1], self.ylimits[0]]) 
        #-----------------------END PLOTS STATEMENTS---------------
            #set xlimits and y limits for model axes 
        # for making a color bar 
        if type(self.cmap) == str:
            self.cmap = cm.get_cmap(self.cmap)
     
        #-----SET TWIN axes for station ticks and labels ---------
        # create twin axis to set ticks to the top station
        # let keep only the axe lines 
        axe2=axm.twiny()
        axe2.xaxis.set_visible(False) 
  
        # show station maker points : 
        for offset , names in zip (occam_data_station_offsets,
                                   occam_data_station_names):
            # plot the station marker ' black triangle down ' 
            # always plots at the surface.
            axm.text(offset/dz  ,
                    self.ylimits[0] - spec_f,  
                    s= self.station_marker,
                    horizontalalignment='center',
                    verticalalignment='baseline',
                    fontdict={'size': self.ms*5, 
                              'color': self.station_color},
                    )
            
            if set_station_label is True :  # then plot label id 
                axm.text(offset/dz ,
                        # get station name closest to station text.
                        self.ylimits[0]/5,    
                        s= names,
                        horizontalalignment='center',
                        verticalalignment='baseline',
                        fontdict={'size': self.ms*3, 
                                  'color': self.station_color},
                        rotation = self.station_label_rotation,
                            )
     
           
        if set_station_label is True : 
            axm.text ((occam_data_station_offsets.max()/dz)/2, 
                       self.ylimits[0] -(self.ylimits[0]/3), 
                        s= 'Stations',
                        horizontalalignment='center',
                        verticalalignment='baseline',
                        fontdict={'size': self.ms*5, 
                                  'color':'k', 
                                  'style': self.font_style,
                                  'weight': self.fw},
                        )
        
        # put a grid on if set to True 
        if self.show_grid is True:
            # axm.grid(alpha=.3, which='major', lw=.35)
            axm.minorticks_on()
            axm.grid(color='k', ls=':', lw =0.5, alpha=grid_alpha, which ='major')
        

         #set color bar properties 
        cbx = mplcb.make_axes(axm,  shrink=self.cb_shrink,
                              pad=self.cb_pad , location ='right' )
        cb = mplcb.ColorbarBase(cbx[0],
                        cmap=self.cmap,
                        norm=mpl.colors.Normalize(vmin=self.climits[0],
                                        vmax=self.climits[1]))
        
        cb.set_label('Resistivity ($\Omega \cdot$m)',
                  fontdict={'size': self.font_size + 1, 
                            'weight': 'bold'})
        
        cb.set_ticks(np.arange(int(self.climits[0]), int(self.climits[1]) + 1))
        # cb.ax.tick_params(axis='y', direction='in', pad=2.)
        cb.set_ticklabels(['10$^{0}$'.format('{' + str(nn) + '}') for nn in
                            np.arange(int(self.climits[0]),
                                      int(self.climits[1]) + 1)])
       
        # set axes labels
        axm.set_xlabel('Distance ({0})'.format(self.depth_scale),
                      fontdict={'size': self.font_size + 2, 'weight': 'bold'})
        axm.set_ylabel('Depth ({0})'.format(self.depth_scale),
                      fontdict={'size': self.font_size + 2, 'weight': 'bold'})
        
        self.fig.suptitle('Plot Resistivity Model :RMS={0}, Roughness={1}'.\
                          format(occam_model_obj.model_rms, 
                                occam_model_obj.model_roughness),
                  ha='center',
          fontsize= 15* self.fs, 
          verticalalignment='center', 
          style =self.font_style,
          bbox =dict(boxstyle='round',facecolor ='moccasin'), 
          y=0.95)
        

        if savefig is not None : 
            plt.savefig(savefig, dpi = self.fig_dpi)
        
        plt.show()

    def plot_Response(self, data_fn =None ,response_fn =None ,
                      mode =None,   **kws ): 
        """
        Function to plot forward value, and residual value from Occam 2D 
        
        list of params are below :
        
        =================  ===============  ===================================
        Params             Type             Description 
        =================  ===============  ===================================
        response_fn         str             full path to occam iteration file           
        data_fn             str             full path to occam_data file              
        doi                 str             depth of investigation might be
                                            float or str like "1km" =1000
        show_station_id     str             show station names 
        =================  ===============  ===================================
            
        :Example: 
            
            >>> from viewer.plot import Plot2d
            >>> pathresp =os.path.join(os.environ ['pyCSAMT'],
            ...                           'pycsamt', 'data', 
            ...         'occam2D','RESP17.resp')
            >>> path_data =os.path.join(os.environ ['pyCSAMT'],
            ...                            'pycsamt', 'data', 
            ...     'occam2D','OccamDataFile.dat' )
            >>> plot2d_obj = plot2d()
            ... plot2d_obj.plot_Response(data_fn =path_data ,
            ...                             response_fn=  pathresp )
            
        """
        self._logging.info('Plot occam pseudosection of forward, residual value ')
        plot_style =kws.pop('plot_style', None)
        
        contourlines =kws.pop('contour_lines_styles', '-')
        contourcolors =kws.pop('contour_lines_colors', 'white')
        delineate_resistivity_curve =kws.pop('delineate_resistivity', None)
    
        show_report = kws.pop('show_report', True)
        
        grid_alpha =kws.pop('alpha', 0.5)
        savefig =kws.pop('savefig', None)
        set_station_label=kws.pop('show_station_id', True)

        #-----------STATEMENT RESPONSE OBJECT ----------------------
        resp_obj = occam2d.Response (data_fn =data_fn ,
                                     response_fn=response_fn)
        # get occam data type  and build large list of possible mode 
        
        make_mode = resp_obj.occam_dtype + [
            str(mm)for mm in resp_obj.occam_mode]
        
        resp_occam_dtype_obj = resp_obj.occam_dtype
        
        #---------MANAGE OCCAM PLOT MODE -------------------------
        # if mode is not provided , then take the first occam mode 

        if mode is None : 
            mode = resp_occam_dtype_obj [0]
        
        mode =str(mode).lower() 
 
        # check the mode if provided 
        if mode not in  make_mode  : 
            mess =''.join(['Occam mode provided ={0} is wrong !.',
                           ' Occam2D data mode is ={1}', 
                           'Please select the right mode.'])
            warnings.warn(mess.format(mode, make_mode ))
            self._logging.error (mess.format(mode, make_mode ))
            
        # check the mode provided , can be str 
        for im , imode in enumerate(resp_obj.occam_mode) :
            if imode ==mode : 
                resp_occam_dtype_obj  =  resp_occam_dtype_obj [im]
                
        for im , imode  in enumerate( [str(mm) 
                                       for mm in resp_obj.occam_mode]): 
            if imode ==mode : 
                resp_occam_dtype_obj  =  resp_occam_dtype_obj [im]  
            
        #----MANAGE COUNTOUR PLOT ------------------------
         #---> get delineate rho curve --- 
        if delineate_resistivity_curve is not None : 
            delineate_resistivity_curve = mplotus.control_delineate_curve(
                res_deline=delineate_resistivity_curve)
            #assert value to put on float rounded to 1 
            # note that value of resistivity for delineate MUST 
            #be on OHM- M not log10 resistivity 
            try : 
                # for consistency 
                if not isinstance (delineate_resistivity_curve, list ): 
                    delineate_resistivity_curve=[delineate_resistivity_curve ]
            except : 
                # be sure to stay on the same
                # output value if something wrong happen 
                
                delineate_resistivity_curve=delineate_resistivity_curve
                pass 
        #----MANAGE PRINT INFO -----------------------------
        print('{0:=^77}'.format('Occam Response plot infos'))
        print('---> Occam 2D plot Mode  = {}'.format(
            mode.split('_')[0].upper() +' ' + mode.split('_')[1]))
        if delineate_resistivity_curve is not None : 
            print('---> Occam 2D contour delineate value = {}'.format(
                tuple(delineate_resistivity_curve)))
 
        #--- CALL RESPONSE OCJECT USEFULL FOR PLOTTING----
        # Now get attributes of forward and residual values 
        self.resp_forward = getattr(resp_obj,
                                    'resp_{0}_forward'.format(mode))
        self.resp_residual =getattr(resp_obj, 
                                    'resp_{0}_residual'.format(mode))
        self.resp_freq = np.log10(getattr(resp_obj,
                                          'data_frequencies'))
        #let elimitale - frequency the lowest one 
        # self.resp_freq = self.resp_freq [:-1]
        self.resp_sites_names = getattr(resp_obj, 'data_sites')
        self.resp_sites_offsets = getattr(resp_obj, 'data_offsets')
        
        if delineate_resistivity_curve is not None :
            print('---> Occam 2D contour delineate value = {}'.format(
                tuple(delineate_resistivity_curve)))
        
        # --DECLARE FIGURE AND PROPERTIES --------------------

        # statement of figures 
        self.fig = plt.figure(self.fig_num, self.fig_size, 
                              dpi=self.fig_dpi,) # constrained_layout=True)
        # plt.clf()
        # add a subplot to the figure with the specified aspect ratio
        self.fig_aspect ='auto'
        gs=gspec.GridSpec(2, 1, figure =self.fig)
        axeFW = self.fig.add_subplot(gs[0, :])
        # axeresd = self.fig.add_subplot(2, 1, 1, aspect=self.fig_aspect)
        axeRESI=self.fig.add_subplot(gs[1,:]  ,sharex = axeFW )
        #                               
        # axeRESI = self.fig.add_subplot(2, 1, 2, aspect=self.fig_aspect, 
        #                                sharex = axefw )
        #-------- SET XY  LIMITS ---------
        # try to make stations separation the same distance between sites 
        self.resp_sites_offsets =np.linspace(self.resp_sites_offsets[0],
                                             self.resp_sites_offsets[-1], 
                                             len(self.resp_sites_offsets))
        
        self.xlimits=(self.resp_sites_offsets.min(), 
                      self.resp_sites_offsets.max())
        self.ylimits = (self.resp_freq.max(), self.resp_freq.min())
        
        # RESIclimits=(self.resp_residual.min(),
        #              self.resp_residual.max())
        RESIclimits =(-2,2)
        #---PLOT STATEMENT --------------------------
        
        # plotstyle is None take a default as pcolormesh
        if plot_style is not None : plot_style =plot_style.lower()
        if plot_style is None : plot_style ='imshow'
        
        print('---> Occam 2D plot style  = "{}"'.format(plot_style))
        
        #--> check dimensionality 
        self.resp_forward, *rr = func.check_dimensionality(
            self,self.resp_forward, self.resp_freq,self.resp_sites_offsets)
        self.resp_freq,self.resp_sites_offsets=rr
        
        if plot_style =='pcolormesh': 
            
            self._logging.info (
                'Ready to plot Forward with matplotlib "pcolormesh"')
            #make mesh_grid 
            mesh_x , mesh_y =np.meshgrid (self.resp_sites_offsets,
                                          self.resp_freq, 
                                          )
            # if you keep plot_x_axis and plot_z_axis in meter ,
            # be sure to divided py dz  meshes respectively 
            #------plot forward response -------
            axeFW.pcolormesh (mesh_x , 
                            mesh_y ,
                             self.resp_forward,
                                vmin = self.climits[0],
                                vmax = self.climits[1],  
                                shading= 'auto', 
                                cmap =self.cmap, 
                                alpha = None, 
                                              
                                              )

            #------plot residual -------------
            
            axeRESI.pcolormesh (mesh_x , 
                                   mesh_y ,
                                   self.resp_residual,
                                    vmin = RESIclimits[0],
                                    vmax = RESIclimits[1],  
                                    shading= 'auto', 
                                    cmap =self.cmap, 
                                    alpha = None, 
                                    
                                              )

        if plot_style.lower() =='imshow': 
            self._logging.info (
                'Ready to plot forward  with matplotlib "imshow"')
            mesh_x , mesh_y =np.meshgrid (self.resp_sites_offsets,
                                          self.resp_freq, 
                                          )
            axeFW.imshow (self.resp_forward,
                                vmax = self.climits[1], 
                                vmin =self.climits[0], 
                                interpolation = self.imshow_interp, 
                                cmap =self.cmap,
                                aspect = self.fig_aspect,
                                origin= 'upper', 
                                extent=( self.xlimits[0],
                                        self.xlimits[1],
                                        self.ylimits[1], 
                                        self.ylimits[0] ), 

                                    )
            axeRESI.imshow (self.resp_residual,
                                vmax = RESIclimits[1], 
                                vmin =RESIclimits[0] ,
                                interpolation = self.imshow_interp, 
                                cmap =self.cmap,
                                aspect = self.fig_aspect,
                                origin= 'upper', 
                                extent=( self.xlimits[0],
                                        self.xlimits[1],
                                        self.ylimits[1], 
                                        self.ylimits[0] ), 

                                    )

        #-------SET AXIS LIMIT -------------------------
        axeFW.set_xlim( [self.xlimits[0],  self.xlimits[1]])
        axeFW.set_ylim ([self.ylimits[1], self.ylimits[0]]) 
        
        axeRESI.set_xlim( [self.xlimits[0],  self.xlimits[1]])
        axeRESI.set_ylim ([self.ylimits[1], self.ylimits[0]]) 
        
        # set delineate cure 
        if  delineate_resistivity_curve is not None :
            for ii, (axe, grid_response) in enumerate(zip(
                    [axeFW, axeRESI], [self.resp_forward, 
                                    self.resp_residual])): 
                
                if ii==0 : mesf = 'forward model'
                else : mesf ='residual'

                try :
                    contps = axe.contour(mesh_x, 
                                         mesh_y,
                                         grid_response,
                                          colors =contourcolors, 
                                          linestyles=contourlines
                                          )
                except :
                    self._logging.debug(
                        'Shape of {0} x & y does not match that of z: found {1}'
                        ' instead of {2} but we will plot instead.'.
                        format(mesf, mesh_x.shape, grid_response.shape))
                    contps = axe.contour(mesh_x[:-1, ::] , 
                                         mesh_y[:-1, ::],
                                         grid_response,
                                          colors =contourcolors, 
                                          linestyles=contourlines
                                          )

                try : 
                    axe.clabel(contps,  delineate_resistivity_curve ,
                                    inline=True, fmt='%1.1f',
                                    fontsize =self.font_size,
                                      )
                except : 
                    #deline =None and set all contours
                    if ii==0 : mesf = 'Forward model'
                    else : mesf ='residual'
                
                    warnings.warn(
                        'Values {0} given as {1} contours levels does not match !'
                        'Contours levels are resseting to default levels !'.
                        format(delineate_resistivity_curve, mesf))
                    
                    print('.--> ! {0} contours levels = {1} are resseting to '
                          ' default levels!'.format(mesf.capitalize(),
                                                    delineate_resistivity_curve))

                    self._logging.debug (
                        'values {0} given as contours levels does not match ! '
                          'availables contours levels are set to default values.')
                    
                    axe.clabel(contps, 
                                inline=True,
                                fmt='%1.1f',
                                fontsize =self.font_size,
                                              )

        for offs , names in zip (self.resp_sites_offsets,
                                 self.resp_sites_names):
            # plot the station marker ' black triangle down ' 
            # always plots at the surface.
            axeFW.text(offs  ,
                    self.ylimits[0],  
                    s= self.station_marker,
                    horizontalalignment='center',
                    verticalalignment='baseline',
                    fontdict={'size': self.ms*4, 
                              'color': self.station_color},
                    )
            
            if set_station_label is True :  # then plot label id 
                axeFW.text(offs,
                        self.ylimits[0] + np.log10( self.ylimits[0]/2.5),   
                        s= names,
                        horizontalalignment='center',
                        verticalalignment='baseline',
                        fontdict={'size': self.ms*2, 
                                  'color': self.station_color},
                        rotation = self.station_label_rotation,
                            )
        if set_station_label is True : 
            axeFW.text ((self.resp_sites_offsets.max()-\
                         self.resp_sites_offsets.min())/2, 
                       self.ylimits[0] + np.log10( self.ylimits[0]/1.25), 
                        s= 'Stations',
                        horizontalalignment='center',
                        verticalalignment='baseline',
                        fontdict={'size': self.ms*3, 
                                  'color':'k', 
                                  'style': self.font_style,
                                  'weight': self.fw},
                        )
        
   
        # put a grid on if set to True 
        if self.show_grid is True:
            # axm.grid(alpha=.3, which='major', lw=.35)
            axeFW.minorticks_on()
            axeRESI.minorticks_on()
            axeFW.grid(color='k', ls=':', lw =0.5,
                       alpha=grid_alpha, which ='major')
            axeRESI.grid(color='k', ls=':', lw =0.5,
                         alpha=grid_alpha, which ='major')

        #---> set color bar properties 
        if type(self.cmap) == str:
            self.cmap = cm.get_cmap(self.cmap)
        for ii, axe in enumerate([axeFW, axeRESI]):
            
            cbx = mplcb.make_axes(axe,  shrink=self.cb_shrink *1.5,
                                  pad=self.cb_pad , location ='right' )
            
            if ii ==0 : 
                cb = mplcb.ColorbarBase(cbx[0],
                                cmap=self.cmap,
                                norm=mpl.colors.Normalize(vmin=self.climits[0],
                                                vmax=self.climits[1]))
                cb.set_label('Resistivity ($\Omega \cdot$m)',
                      fontdict={'size': self.font_size , 
                                'weight': 'bold'})
            
                cb.set_ticks(np.arange(int(self.climits[0]), 
                                       int(self.climits[1]) + 1))
                # cb.ax.tick_params(axis='y', direction='in', pad=2.)
                cb.set_ticklabels(['10$^{0}$'.format('{' + str(nn) + '}') 
                                    for nn in np.arange(int(self.climits[0]),
                                              int(self.climits[1]) + 1)])
                

            elif ii ==1 : 
                #set new color bar limits ....
                cmin = np.floor(RESIclimits[0])
                cmax = np.ceil (RESIclimits[1])
                
                cb = mplcb.ColorbarBase(cbx[0],
                                cmap=self.cmap,
                                norm=mpl.colors.Normalize(vmin=cmin,
                                                vmax=cmax))

                cb.set_label('Resistivity ($\Omega \cdot$m)',
                          fontdict={'size': self.font_size , 
                                    'weight': 'bold'})
                
                cb.set_ticks(np.arange(int(cmin),
                                       int(cmax) + 1))

                cb.set_ticklabels(['${0}$'.format(str(nn)) for nn in
                                    np.arange(int(cmin),
                                              int(cmax) + 1)])
                
                axeRESI.set_xlabel('Distance ({0})'.format('m'),
                          fontdict={'size': self.font_size, 'weight': 'bold'})
           
        # ---> set axes labels
           
            axeFW.set_ylabel('log10 Frequency ({0})'.format('Hz'),
                          fontdict={'size': self.font_size , 'weight': 'bold'})
            axeRESI.set_ylabel('log10 Frequency ({0})'.format('Hz'),
                          fontdict={'size': self.font_size , 'weight': 'bold'})
        
        # let set the axis of forward invisible 
        plt.setp(axeFW.get_xticklabels(), visible=False)

        if show_report is True : 
            mplotus.get_conductive_and_resistive_zone (
                data = self.resp_forward,
                site_names =self.resp_sites_names)

        self.fig.suptitle('Occam 2D : Forward response & Residual',
                  ha='center',
          fontsize= 7* self.fs, 
          verticalalignment='center', 
          style =self.font_style,
          bbox =dict(boxstyle='round',facecolor ='moccasin'), y=1.)
        
        # tigh_layout   plt.tight_layout

        # ---------------------Print Reprot ----------------------------
       
        if savefig is not None : 
            plt.savefig(savefig, dpi = self.fig_dpi)
        
        plt.show()
        
    def plot_Pseudolog(self, station_id= 'S00', iter_fn=None ,
                       mesh_fn=None , data_fn=None , iter2dat_fn=None , 
                      bln_fn=None , model_fn=None , **kwargs):
        """
        Build pseudodrill from the model resistivity . 
        
        Deal with  True resistivities get on the survey or with others firms.
        Indeed the input  resistivities values into our model,could yield 
        an accuracy underground map.The challenge to build a pseudolog is
        to know how the layers are disposal in underground so to emphasize
        the large conductive zone  especially in the case of groundwater
        exploration. 
        
        Program works in combinaison with geophysic data especially Occam 2D 
        inversion data,  with geological data. Actually the program deals only
        with  Occam 2D inversion files or  Bo Yang (x,y,z) files and outputs some 
        model files for other external softares such as Golder sofware('surfer'). 
        If user has a golder software installed on its computer,  he can use the 
        output files generated to produce 2D map so to compare both maps: 
        the model map and detail-sequences map to see the difference  between
        the "calculated model" and "pseudosequences model" which could match
        better the underground. When `step_descent` parameter is small,
        the detail sequences trends to the model map. So more geological data
        are, better the accuracy of detail sequences logs should be. 
        
        :param station_id:  Number or the site id of the survey area 
                            number starts from 1 to the end .
        :type station_id:  str, int
            
        .. note:: User caneither use Occam 2D inversions files to plot or 
                BoYang (x, y, file)+ station location file (*bln) to plot 
                if the two types of files are provided , program with give 
                priority to Occam 2D inversion files.

        ======================  ==========  ===================================
        Params                   Type       Description 
        ======================  ==========  ===================================
        model_fn                str         full path to Occam model file .
        iter_fn                 str         full path to occam iteration file 
        data_fn                 str         full path to occam_data file 
        doi                     str         depth of investigation might be 
                                            float or str like "1km" =1000
        depth_scale             str         scale of imaging depth can be 
                                            "km" or "m". Default is"m"
        step_descent            float       step to enforce the model 
                                            resistivities to keep truth layers
                                            values as reference data . if step 
                                            descentis egal to doi max, data 
                                            looks like model at 99.99%.Step  
                                            decent is function of depth and rho.
        lc_AD_curves            tuple       customize line color of average
                                            curve and details sequneces 
                                            logs  eg : ((0.5, 0.8, 0.),'blue') 
        default_unknow_lcolor   str         In the case  the name of layer   
                                            is notin our data base ,
                                            customize the layer color .
                                            default is "(1.0, 1.0, 1.0)".
        default_unknow_lpatter  str         In the case  the name of layer  
                                            is not in our data base ,
                                            customize the layer pattern 
                                            default is "+.+.+."
        ======================  ==========  ===================================
        
        .. note:: constrained_electrical _properties_of_rocks param 
             keeps the Truth layers resistivities  as reference resistivities.
             If value is false will check in our    data base to find the 
             resistivities that match better the given resistivities of 
              the layers. *Default* is True.                                   
            
        Customize your plot using matplotlib properties. 
        
        :Example: 
            
            >>> from viewer.plot import Plot2d 
            >>> path =os.path.join(os.environ ['pyCSAMT'],
            ...                       'pycsamt', 'data', 'occam2D')
            >>> plot2d_obj = Plot2d(station_label_rotation=None, 
            ...                    show_grid=True, 
            ...                    font_size =8, 
            ...                    lc='r', 
            ...                    fig_size=[5,8], 
            ...                    markerfacecolor='k', 
            ...                    markeredgecolor='k')
            >>> plot2d_obj.plot_Pseudolog( station_id=[43], 
            ...                          input_resistivities=[300, 500,
            ...                                             1000, 2000,
            ...                                             4000, 6000],
            ...                          input_layers =['alluvium', 
            ...                                         'amphibolite',
            ...                                         'altered rock',
            ...                                          'augen gneiss',
            ...                                          'granite'],
            ...                            mesh_fn=os.path.join(path, 
            ...                                             'Occam2DMesh')
            ...                            iter_fn = os.path.join(path,
            ...                                             'ITER17.iter'), 
            ...                            model_fn =os.path.join(path,
            ...                                          'Occam2DModel') , 
            ...                            data_fn =os.path.join(path, 
            ...                                       'OccamDataFile.dat'),
            ...                            doi='1km', 
            ...                            step_descent=200., 
            ...                            plot_style= 'pcolormesh')
        """
        
        self._logging.info('Building pseudo-drill and pseudostratigraphy .')
        
        #set other important kwargs argumenents 
        
        savefig =kwargs.pop('savefigure', None)
        input_layers =kwargs.pop('input_layers', None)
        input_resistivities =kwargs.pop('input_resistivities', None)
        doi = kwargs.pop('doi', None)           # if None , default is 1km 
        mfontdict =kwargs.pop('font_dict_site', 
                              {'size': self.ms*6 , 
                                                  'color': 'saddlebrown', 
                                                  'weight': self.fw, 
                                                  'style': 'italic'}, 
                                                                       )
        
        pseudo_plot_style =kwargs.pop('plot_style', None)
        depth_scale =kwargs.pop('scale', None)
        step_descent = kwargs.pop('step_descent', None)
        cdict_average_detailsC =kwargs.pop('lc_AD_curves', 
                                           ((0.5, 0.8, 0.), 'blue')
                                           )
        
        constrained_electrical_properties_of_rocks=kwargs.pop(
            'constrained_electrical_properties_of_rocks', True)
        default_unknow_lcolor =kwargs.pop('unknow_layer_color', 'white')
        default_unknow_lpattern = kwargs.pop('unknow_layer_pattern', '+.+.+.')
        
        for key in list(kwargs.keys()):
            setattr(self, key, kwargs[key]) 
            
        # --- Ascertainment of differents files --------------------------
        # indicator with file provided for reading , program withcheck
        f,p=0,0      #the corresponding files # even mayfiles is provied ,
                        # default is read Occam 2D files 
                                
                                
        imf, imp =[[]for i in range(2)]
        
                   # intend to read occam 2D files , let assert files now 
        
        for nn, mm in zip(['model_fn', 'iter_fn', 'mesh_fn', 'data_fn'],
                      [model_fn , iter_fn , mesh_fn, data_fn]): 
            if mm is not None : 
                setattr(self, nn, mm)
                f +=1
            elif hasattr( Plot2d, nn):
                f +=1
            else :  imf.append(nn)
        
        if f !=4 : 
            for nn, mm in zip(['iter2dat_fn', 'bln_fn'], [iter2dat_fn, bln_fn]): 
                if mm is not None : 
                    setattr(self, nn, mm)
                    p +=1
                elif hasattr( Plot2d, mm):
                    p +=1
                else :imp.append(nn)
                    
        if p !=2 and f!=4 :
            if p!=2 and p!=0 : 
                print('--> Expected 02 files , only {0}  is given !'.format(p))
                for ss in imp : 
                    mess =' ! <%s> is not provided ! Could not possible to '\
                        'build pseudodrill and stratigraphy log.'\
                        ' Please provide <%s> file.'%(ss,ss.split('_')[0]) 
                    warnings.warn(mess)
                    self._logging.error (mess)
                    if ss=='bln_fn': 
                         mess =''.join(
                             [' ! No station station locations file '
                        'found like *bln file ! Can not plot pseudodrill and ', 
                          'pseudostratigraphy log. Please provided a station'
                              ' location file. You can use Iter2Dat model ::', 
                          'from pycsamt.modeling.occam2d import Iter2Dat:: '
                              'or Profile module ` from pycsamt.core.cs ',
                              'import Profile` to', 
                          'build a station location file.'])
                         warnings.warn(mess)
                         self._logging.error (mess)
 
                
            if f !=4 and f !=0:
                if f >2 : mt='are'
                else :mt='is'
                print('--> Expected 04 files, '
                      'only {0} one {1} given !'.format(f, mt))
                for oss in imf : 
                    mess =''.join([
                        ' ! {0} is not provided!Could not possible',
                        'to build pseudodrill and stratigraphy log.',
                            ' Please provide {0} file.'.format(
                                oss.split('_')[0])]) 
                    
                    self._logging.error(mess)
            
            if p==0 and f==0 : 
                mess =''.join(['None files are found ! Please provided either',
                    'Occam2D outputfiles {Mesh|Model|Data|iter} ',
                    'or Bo Yang Data File output files {Iter2dat|bln}.'])
                warnings.warn(mess)
                self._logging.error(mess)
        
            raise CSex.pyCSAMTError_plot_geoinputargument(mess)
        #----Finish ascertainement then build object and read ----
        elif p==2 or f==4 :    
            if f==4 :                       # priority to Occam2D data files 
                print('**{0:<37} {1} {2}'.format(
                    ' Occam Input files ','=' ,
                    tuple([os.path.basename(file) for
                        file in [self.model_fn , 
                                 self.mesh_fn , 
                                  self.data_fn, 
                                  self.iter_fn]]))) # show message to user
                
                geo_obj =geoD.Geodrill(model_fn= self.model_fn , 
                              data_fn =self.data_fn , 
                              mesh_fn =self.mesh_fn , 
                              iter_fn =self.iter_fn , 
                              input_layers = input_layers , 
                              input_resistivities =input_resistivities, 
                              doi =doi )
                
                # geo_obj.geo_replace_rho()
                
            elif p==2 : 
                print('**{0:<37} {1} {2}'.format(
                    ' Iter Input files ','=' ,
                    tuple([ os.path.basename(file) for 
                                    file in [self.iter2dat_fn , 
                                     self.bln_fn ]]))) 
                
                geo_obj =geoD.Geodrill(iter2dat_fn =self.iter2dat_fn, 
                                    bln_fn =self.bln_fn , 
                                    input_layers =input_layers , 
                                    input_resistivities =input_resistivities, 
                                    doi =doi )
                #------------------ if input resistivities is None ,
                # let get input fron specail station ID -----------------

        # ----------> Get all other attributes from geo_obj -------------------
        # Note data extract here are all in ohm meter not in log10 resistivities 
        # let get geo_d == dictionnaries of sation and resistivities framed 
        # get geodeth is investigation depth depth 
        self.geo_stations_dict = geo_obj.geo_d 
        self.geo_depth = geo_obj.geo_depth 
    
        self.station_names = geo_obj.station_names 
        self.station_location =geo_obj.station_location 
        
        self.step_descent = step_descent
        
        #--Get the corresponding stationid  to Plot user have possibility for
        #  multiples plot by puting station id into a list 
        # get station names on list 
        self.station_id = mplotus.get_stationid(stations=self.station_names,
                                                station_id=station_id) 
        
            # build em dict from station to get
            #a special  input resistivities layers -------
        
        #for stn , res_values in self.geo_stations_dict.items(): 
        if input_resistivities is None : 
            self.input_resistivities ={stn :
                mplotus.get_station_id_input_resistivities(
                    station_rho_value=res_values) 
                    for stn , res_values in self.geo_stations_dict.items()}
        
        #---Figure statements and properties -------------------------

        self.fig = plt.figure(self.fig_num, self.fig_size, dpi=self.fig_dpi )
  
        self.fig_aspect ='auto'   # blocked to automatic , can not change 

        gs=gspec.GridSpec(nrows =7, ncols= 7, figure =self.fig) 
        
        axe_Titles = self.fig.add_subplot(gs[0, :6])
        axe_Titles.set_visible(False)   # make title axes invisible  
        
        
        axePseudodrill = self.fig.add_subplot(gs[1:6, 0])
        axeLogRho= self.fig.add_subplot(gs[1:6,1:3], sharey= axePseudodrill)
        axePseudosequences = self.fig.add_subplot(gs[1:6,3:], 
                                                  sharey=  axePseudodrill)
  
        axe_Legends_rho_sequences = self.fig.add_subplot(gs[6, :4])
        
        axe_Legends_rho_sequences.set_xticks([])
        axe_Legends_rho_sequences.set_yticks([])

        # ----> Manage XY limits --------------------------------------
        
        # build  the distance separation by 
        # rounding value using function round dipole length 
        
        df = func.round_dipole_length(int(self.station_location.max()
                                          /(len(self.station_location)-1)))
        dx =df/5        # xpad = dx
        dy=30           # ypad =dy      to locate station
   
        # ---> build a smallsites offsets set x 
        #pseudodrill limits  of three stations 
        
        plot_sites_offsets = np.arange(df, 4*df , df)
        # take the indice and find two next statiion and build 
        #corresponding stationnames 
        for ids in self.station_id : 
            ids=int(ids.replace('S', ''))
            if ids ==0 :             # the first station 
                # if station is the first station names 
                plot_sites_names =['S{0:02}'.format(ii) for ii in range(3)]  
            # that means the last stations     
            elif ids == len(self.station_names)-1 :                    
                plot_sites_names= ['S{0:02}'.format(ii) for
                                   ii in range(len(self.station_names)-3, 
                                                len(self.station_names)) ]
            else : # let framed the main station id 
                plot_sites_names= ['S{0:02}'.format(ii) for ii in range(ids-1,  
                                                             ids+2) ]

        #----let get default plot and set scale ----- 
        if pseudo_plot_style==None  :
            pseudo_plot_style = 'imshow'
            
        if depth_scale is None : depth_scale ='m'
       
        if depth_scale.lower()=='m':
            dz= 1.
        elif depth_scale.lower() =='km': 
            dz=1000.
        else :               # in the case a wrong value is provided
            dz=1.
            
         # let graduate to eg [-10 +160] -->if df=50 as step     
        self.xlimits=((plot_sites_offsets[0]-dx)/dz,
                      (plot_sites_offsets[-1]+dx)/dz)           
        self.ylimits = (self.geo_depth.min()/dz,
                        self.geo_depth.max()/dz )    # reversed axis order     
            
            
        for stn in self.station_id : 
            #----READ OBJECT AND PUT THE SPECIAL INPUT RESISTIVITIES -------
            if input_resistivities is None : 
                input_resistivities = self.input_resistivities[stn]
            
            geo_obj.geo_build_strata_logs(step_descent = self.step_descent)
           
            
            self.stations_replaced_rho = geo_obj.geo_drr
            
            # set xlog limits  return to log10 resistivities
 
            maxlog =np.log10(self.geo_stations_dict[stn]).max()
            
            minlog = np.log10(self.geo_stations_dict[stn]).min()
            self.xloglimits= (minlog, maxlog) 
            
                    #--------------------------------BUILD PSEUDO DRILL -----
            self._logging.info('Build Pseudo-drill from Occam 2D models.')
            
            if pseudo_plot_style.lower() == 'pcolormesh': 
                
                self._logging.info (
                    'Ready to plot Pseudodrill  with matplotlib "pcolormesh"')
                
                mesh_x , mesh_y = np.meshgrid(
                    plot_sites_offsets/dz, self.geo_depth/dz)

                axePseudodrill.pcolormesh (mesh_x , 
                                           mesh_y ,
                                           np.log10(self.geo_stations_dict[stn]),
                                    vmin = self.climits[0],
                                    vmax = self.climits[1],  
                                    shading= 'auto', 
                                    cmap =self.cmap, 
                                    alpha = None, 
                                                  )
            
            if pseudo_plot_style.lower() =='imshow': 
                
                self._logging.info (
                    'Ready to plot  Pseudodrill with matplotlib "imshow"')
                
                axePseudodrill.imshow (np.log10(self.geo_stations_dict[stn]),
                                    vmax = self.climits[1], 
                                    vmin =self.climits[0], 
                                    interpolation = self.imshow_interp, 
                                    cmap =self.cmap,
                                    aspect = self.fig_aspect,
                                    origin= 'upper', 
                                    extent=( self.xlimits[0]+dx,
                                            self.xlimits[1]-dx,
                                            self.ylimits[1], 
                                            self.ylimits[0] )
                                    )
                # to get the origine =0 of the plot 
            
            # put a grid on if set to True 
            if self.show_grid is True:
                axePseudodrill.minorticks_on()
                axePseudodrill.grid(color='k',
                                    ls=':',
                                    lw =0.5, 
                                    alpha=0.8,
                                    which ='major')
                

            #---- Locate the axe station id
            for offs , names in zip (plot_sites_offsets,
                                     plot_sites_names):
                # plot the station marker ' black triangle down ' 
                # always plots at the surface.
                if names == stn :

                    axePseudodrill.text(offs/dz  ,
                        self.ylimits[0]-2*dy/dz,  
                        s= names,
                        horizontalalignment='center',
                        verticalalignment='baseline',
                        fontdict=mfontdict
                        )
                    axePseudodrill.text(offs/dz  ,
                        self.ylimits[0],  
                        s= self.station_marker,
                        horizontalalignment='center',
                        verticalalignment='baseline',
                        fontdict={'size': self.ms*6, 
                                  'color': mfontdict['color']},     
                        )
                
                    
                else :
                    fdict ={'size': self.ms*4, 
                                  'color': self.station_color} 
                    
                    axePseudodrill.text(offs/dz  ,
                            self.ylimits[0],  
                            s= self.station_marker,
                            horizontalalignment='center',
                            verticalalignment='baseline',
                            fontdict=fdict
                            )                   
            # create a temporary axe to hold station names 
            tempax= axePseudodrill.twiny()
            # customize frame stations id 
            
            for jj , istn in enumerate(plot_sites_names) : 
                if istn ==stn :
                    new_names = plot_sites_names
                    new_offs = plot_sites_offsets/dz
                    del(new_names[jj])
                    new_offs= new_offs.tolist()
                    del(new_offs[jj])
    
                    tempax.set_xticks(ticks= new_offs, minor=False ) 
                    tempax.set_xticklabels( new_names,               
                               rotation=self.station_label_rotation,
                               fontdict={'size': self.ms*3}
                              )

            # set y labels of pseudodrill plots 
            axePseudodrill.minorticks_on()
            axePseudodrill.set_ylabel('Depth({0})'.format(depth_scale),
                          fontdict={'size': self.font_size , 
                                    'weight': self.fw, 
                                    })
            
            axePseudodrill.tick_params(axis='y', labelsize=self.font_size )
            
   
            
            tempax.set_xlabel('Stations',
                          fontdict={'size': self.font_size , 
                                    'weight': self.fw, 
                                    'style':self.font_style})
            # set the x twins limits the same as the x of pseudodrill 
            tempax.set_xlim([self.xlimits[0], self.xlimits[-1]]) 


             #------------------------------plot LOGRHO--------------------
            self._logging.info('Build  Plot1D resistivities sounding curves.')
            
            # Create twin axe to host log 10 rho values and xlables  
            axisLogcurve=axeLogRho.twiny()
            if self.show_grid is True: 
                axeLogRho.grid(axis ='y', color='gray', ls=':',
                               lw =0.3, alpha=0.8, which ='both')
            #Xloglimits = np.arange(int(minlog)-1, int(maxlog)+2, 1) 
            # No need , let do it automatically
      
            
            # get the index of station id to plot 
            mms =int(stn.replace('S',''))
            if mms ==0 : 
                indexplot =0 
                # mean we are at the last station 
            elif mms == len(self.station_names)-1 :     
                indexplot = -1
            else : indexplot =1  # plot the framed station 
            

            # get the array value at the main station id framed into Two 
            plot_y = self.geo_stations_dict[stn][:, indexplot] #3,1
            # plot the replace rho using the same index stn  
            AverageLogcurve, =  axisLogcurve.semilogx(
                self.stations_replaced_rho[stn][:, indexplot],
                        self.geo_depth/dz, 
                        lw= self.lw*2 , 
                        ls=self.ls,
                        c= cdict_average_detailsC[0], 
                        label ='Average logcurve (Ω.m)/station {0}'.format(stn)
                           
                            )
            # plot the details log curves using the step descent 
            DetailLogcurve, =  axisLogcurve.semilogx( 
                geo_obj.geo_dstep_descent[stn][:, indexplot],
                            self.geo_depth/dz, 
                            lw= self.lw *2 , 
                            ls=self.ls,
                            c= cdict_average_detailsC[1], 
                            alpha=0.5,
                            label='Detailsequence logcurve(Ω.m) :'\
                                'step_decsent ={0}{1}.'.\
                                format(self.step_descent/dz, depth_scale), 
                            
                            )
            
            # make a marker (optional) to better mark the curve 
            mark, =  axisLogcurve.semilogx(plot_y,
                            self.geo_depth/dz,
                          marker ='o',
                          c='white', 
                          markersize= 2* self.ms , 
                          
                          markeredgecolor=self.markeredgecolor, 
                          markerfacecolor=self.markerfacecolor
                          ) 
            # make the sounding cure at  the main station id 
            ResLogcurve , =  axisLogcurve.semilogx(plot_y,
                            self.geo_depth/dz, 
                           lw= self.lw , 
                           ls=self.ls,
                           c= self.lc,
                           label = 'Resistivity logcurve (Ω.m)', 
                           
                          
                           )

            # setlog grid , turn on minorticks 
            if self.show_grid  is True : 
                axisLogcurve.minorticks_on()
                axisLogcurve.grid(axis ='both', 
                                  color='gray', 
                                  ls=':', 
                                  lw =0.3, 
                                  alpha=0.8, 
                                  which ='both')
        
            axisLogcurve.set_xlabel('Log10(Rho)/Ω.m',
                          fontdict={'size': self.font_size , 'weight': 'bold', 
                                    'style':self.font_style})
            
            # annotate step descent 
            axisLogcurve.text(plot_y.min(), 
                              self.ylimits[1], 
                              '$Step\ descent = {0}{1}.$'.\
                                  format(self.step_descent/dz , depth_scale),
                               verticalalignment='bottom', 
                               c='k', fontsize =self.font_size/1.5)
                             
                               
            #                    )
            # ----set LEGEND vers 
            axe_Legends_rho_sequences.legend(
                [ ResLogcurve, AverageLogcurve,DetailLogcurve],
                ['Resistivity soundingcurve (Ω.m)/Station = {0}'.format(stn),
                 'Average resistivity logcurve (Ω.m) based on resistivity calculation',
                ' Pseudo details-sequences logcurve',
                          ], 
                         # ncol= 3,
                         prop={'size':self.ms*3, 
                               'weight': self.fw, 
                               'style': self.font_style}, 
                          mode='expand',
                         edgecolor = 'white', 
                         )

            #----------------------------PLOT PSEUDSEDQUENCES  --------------
            self._logging.info(
                'Build Pseudo sequences with delais logs'
                    ' sequences curve and average curves..')
            # call the speudosecquence object from geo_obj 

            self.geo_dpseudo_sequence_thickness =geo_obj.geo_dpseudo_sequence 
            # resistivities at every layer thickess 
            self.geo_dpseudo_sequence_rho = geo_obj.geo_dpseudo_sequence_rho 
            # resseting layer names , color,
            #and pattern according to their resitivities 

            self.input_layers , self.layer_color , self.layer_pattern = \
                geo_obj.get_geo_formation_properties(
                    structures_resistivities = \
                     self.geo_dpseudo_sequence_rho[stn], 
                     real_layer_names = input_layers, 
                      constrained_electrical_properties_of_rocks=\
                          constrained_electrical_properties_of_rocks, 
                       default_layer_color =default_unknow_lcolor  , 
                       default_layer_pattern = default_unknow_lpattern,
                                                     )
                    
             # create twniny so to get station at the top 
            axePS = axePseudosequences.twiny()         
            
            # prepare a list of cumulative sum for bar plot setting everytimes ,
            # the bottom as the top of the next bar 
            # loop value of pseudo_sequences  and create bar plot 
            
            for ii, pseuds in enumerate(
                    self.geo_dpseudo_sequence_thickness[stn]):
                # sum all the previous bars sequences
                next_bottom_bar = self.geo_dpseudo_sequence_thickness[stn][:ii].sum() 
     
                axePS.bar(plot_sites_offsets[0]/(dz*2), 
                          pseuds, 
                          bottom =next_bottom_bar/dz, 
                          width=df/(dz*2), 
                          linewidth = self.lw/4, 
                          edgecolor ='k', 
                          hatch = self.layer_pattern[ii], 
                          color = self.layer_color [ii], 
                          alpha =1.,
            
                          )
            # prepare a cumul sum  for annotation starting to
            # 0 and terminate at the depth minums 1 
            annotate_cumsum = self.geo_dpseudo_sequence_thickness[stn].cumsum()
            annotate_cumsum = np.concatenate(
                (np.array([0]),annotate_cumsum[:-1]))
              
            _, annotate_lnames = mplotus.annotate_tip(
                layer_thickness= annotate_cumsum ,
                                 layer_names =self.input_layers )

            # get station name closest to station text.  
            axePS.text(plot_sites_offsets[0]/(dz *2),
                      self.ylimits[0] - 2*dy/dz,  
                      s= stn,
                      horizontalalignment='center',
                      verticalalignment='baseline',
                      fontdict=mfontdict,
                      rotation = self.station_label_rotation,
                          )
            # axePS.set_xticks
            # location of the station marker 
            axePS.text(plot_sites_offsets[0]/(dz *2) ,
                        self.ylimits[0],  
                        s= self.station_marker,
                        horizontalalignment='center',
                        verticalalignment='baseline',
                        fontdict={'size': self.ms*6, 
                                  'color': mfontdict['color']} ,
                        )
 
                
            axePS.legend(labels=[ ln.capitalize() 
                                 for ln in annotate_lnames], 
                         loc ='upper right', 
                         prop ={'size':self.font_size , 
                                'style': self.font_style, 
                                }, 
                         edgecolor ='white', 
                         )
  
            # position of site and marker 
            # create a temporary axe t ohost x labels 
            axePS.set_xticks(ticks=[plot_sites_offsets[0]/(dz *2)],
                             minor=False )
            axePS.set_xticklabels([stn] , 
                                  rotation=self.station_label_rotation, 
                                  color ='white',
                                  )

            axePS.set_xlabel('Pseudo-sequences',
                          fontdict={'size': self.font_size , 
                                    'weight': 'bold', 
                                    'style':self.font_style})
   
            
            #---> set color bar properties 
            if type(self.cmap) == str:
                self.cmap = cm.get_cmap(self.cmap)
            
            cbx = mplcb.make_axes(axePseudosequences,  
                                  shrink=self.cb_shrink *1.3,
                                  pad=self.cb_pad , 
                                  location ='right' )

            cb = mplcb.ColorbarBase(cbx[0],
                            cmap=self.cmap,
                            norm=mpl.colors.Normalize(vmin=self.climits[0],
                                            vmax=self.climits[1]),
                            )
            cb.set_label('Resistivity ($\Omega \cdot$m)',
                  fontdict={'size': self.font_size , 
                            'weight': 'bold'})
        
            cb.set_ticks(np.arange(int(self.climits[0]), 
                                   int(self.climits[1]) + 1))

            cb.set_ticklabels(['10$^{0}$'.format('{' + str(nn) + '}') for nn in
                                np.arange(int(self.climits[0]),
                                          int(self.climits[1]) + 1)])
                
            #----SETUP  ALL AXIS LIMIT SALL (3 plots)  -----------------------
            
            axePseudodrill.set_xlim( [self.xlimits[0],  self.xlimits[1]])
            axePseudodrill.set_ylim ([self.ylimits[1], self.ylimits[0]]) 
            #make both axes sites invisibles 
            axePS.set_xlim([0, self.xlimits[-1]])  # start lmits at 0 
            
            # axePS.set_xticks([])
            axePseudosequences.set_xticks([]) 
            
        
            # let set the axis of three plots invisibles 
            plt.setp(axePseudodrill.get_xticklabels(), visible=False)
            plt.setp(axePseudosequences.get_xticklabels(), visible=False)
            
            plt.setp(axeLogRho.get_yticklabels(), visible=False)
            plt.setp(axePseudosequences.get_yticklabels(), visible=False)
            
            # make axis invisible to clearly let the legend 
            axe_Legends_rho_sequences.set_axis_off()
            
            axeLogRho.set_xticks([])

            # axePseudodrill.set_xticks([])  
    
            self.fig.suptitle(
                'Pseudosequence of Station : {0}'.format(stn),
                         fontsize=  self.font_size *1.2, 
                         verticalalignment='center', 
                         style ='italic',
                         bbox =dict(boxstyle='round',facecolor ='moccasin'), 
                         y=0.95)
            
            if savefig is not None :
                plt.savefig(savefig , dpi = self.fig_dpi)

     
@geoplot1d(reason = 'zonge_engineering', color_mode='bw', 
           linebbox_kws={'boxstyle':'square','facecolor':
                         'whitesmoke', 
                          'color':'white'
                         })
def plot_dataAndFits(data_fn =None, stations =None, **kws): 
    """ 
    Plot error Apparent resistivity and phase and error in data  
    for multiples zonge lines(AVG)
    
    :param data_fn:
        Path to zonge *.avg file . Can be a list of avg files 
    :param stations: station id . Can be a list of a staion or 
        the station number. e.g. ::
        
        station ='S00' 
        sttion ='1'
        station =['s00', 13, 46]
    :param kws: keywords arguments. Please refer to 
                :class:`pycsamt.core.cs.CSAMT` for further details
        
    :Returns: A list of:
        - z_lines: List of survey lines
        - z_stations: list of station to visualize 
        - z_freq: list of frequency to plot on x-axis 
        - z_appRHO: list of apparent resistivy in ohm.m to plot 
        - z_phase: list of phase in degrees 
        - z_appRho_err: list of error apparent resistivity
        - z_phase_err: list of error phase in %
        
    :Example:
        
        >>> from pycsamt.view.plot import plot_dataAndFits
        >>> path =r'F:\ThesisImp\avg'
        >>> pathData = [os.path.join(path, file) 
                for file in ['K1.AVG',
                             'K4.AVG', 'K6.AVG', 'K8.AVG'
                             ]]
        >>> z_lines, z_stations, z_freq,  z_appRHO, z_phase,z_appRho_err,\
            z_phase_err=plot_dataAndFits(
                data_fn = pathData, 
                stations=['S00', 'S04', 's06', 's10'])
    """
 
    # manage station and create duplicate list
    if data_fn is None: 
        raise CSex.pyCSAMTError_AVG(
            'No `avg` file detected. Please provide a file.')
    elif isinstance(data_fn, str): 
        data_fn = [data_fn]
    
    if stations is None : 
        _logger.debug ( 'None sation found can not be plotted.'
                       'Should be considered the default value =`S00`.')
        stations=['S00']

    # check the sation and duplicate until the max gridspec =4 
    if isinstance(data_fn, list): 
        if isinstance(stations, str): 
            # duplicate the 
            station_list =[stations for i in range(len(data_fn))]
        elif isinstance(stations, list): 
            if len(stations)< len(data_fn):# max grid spec = 4 
                station_list = stations + ['S00' 
                                          for i in range(4-len(stations))]
            elif len(stations) > len(data_fn): 
                station_list =stations[:4]
            else : station_list= stations 
        
     
    def read_singleLine(stn , fn): 
        """ Read Zonge AVG objet and get the following attributes of : 
            - id : name of the line 
            - station : station to plot 
            - freq : freq to plot 
            - appRho: apperent resistivity in ohm.m 
            - phase: Phase value in degrees 
            - appRho_err: Error in resistivity in % 
            - phase_err: error in phase in % 
            
        :param fn: Aonge engineering filename
        :param_stn: station name: check whether the station provided
        is  able to found in the station list of the survey line.
        """
        line_id = os.path.basename(os.path.splitext(fn)[0])
        
        csamt_obj = CSAMT(data_fn = fn, **kws)
        csamt_station_list = csamt_obj.station 
        
        index_station_lst= mplotus.station_id(stn)
        
        try : 
            stn = csamt_obj.station [index_station_lst]
        except:
            len_stnobj =len(csamt_obj.station) 
            _logger.debug(
            f" Too {'large' if index_station_lst >len_stnobj else 'small'}."
            "values. Default station should be `S00`.")
            stn = csamt_station_list[0]

    
        appRho= csamt_obj.resistivity [stn]
        phase =csamt_obj.phase [stn]
        appRho_err= csamt_obj.resistivity_err[stn]
        phase_err = csamt_obj.phase_err[stn]
        freq = csamt_obj.freq 
        
        return (line_id, stn, freq, appRho, phase,
                appRho_err, phase_err)           
                
    # read the data and collect data 
    z_lines, z_stations, z_freq, z_appRHO, z_phase,\
        z_appRho_err, z_phase_err=[[]for i in range(7)]

    for ii, lines in enumerate(data_fn) : 
        lid_, stn_, freq_ , rho_, phase_ , \
            rho_err_ , phase_err_ = read_singleLine(
                stn= station_list[ii],fn=lines)
        z_lines.append(lid_)
        z_stations.append(stn_)
        z_freq.append(freq_)
        z_appRHO.append((rho_, None))
        z_phase.append((phase_, None))
        z_appRho_err.append(rho_err_)
        z_phase_err.append(phase_err_)
        
    return (z_lines, z_stations, z_freq, 
            z_appRHO, z_phase,z_appRho_err, z_phase_err)
                                                          
            
         
class PlotResponse(plot_response.PlotResponse):
    """
    plot data and response

    Plots the real and imaginary impedance and induction vector if present.

    :Example: ::

        >>> import pycsamt.modeling.modem as modem
        >>> dfn = r"/home/MT/ModEM/Inv1/DataFile.dat"
        >>> rfn = r"/home/MT/ModEM/Inv1/Test_resp_000.dat"
        >>> mrp = modem.PlotResponse(data_fn=dfn, resp_fn=rfn)
        >>> # plot only the TE and TM modes
        >>> mrp.plot_component = 2
        >>> mrp.redraw_plot()

    ======================== ==================================================
    Attributes               Description
    ======================== ==================================================
    color_mode               [ 'color' | 'bw' ] color or black and white plots
    cted                     color for data Z_XX and Z_XY mode
    ctem                     color for model Z_XX and Z_XY mode
    ctmd                     color for data Z_YX and Z_YY mode
    ctmm                     color for model Z_YX and Z_YY mode
    data_fn                  full path to data file
    data_object              WSResponse instance
    e_capsize                cap size of error bars in points (*default* is .5)
    e_capthick               cap thickness of error bars in points (*default*
                             is 1)
    fig_dpi                  resolution of figure in dots-per-inch (300)
    fig_list                 list of matplotlib.figure instances for plots
    fig_size                 size of figure in inches (*default* is [6, 6])
    font_size                size of font for tick labels, axes labels are
                             font_size+2 (*default* is 7)
    legend_border_axes_pad   padding between legend box and axes
    legend_border_pad        padding between border of legend and symbols
    legend_handle_text_pad   padding between text labels and symbols of legend
    legend_label_spacing     padding between labels
    legend_loc               location of legend
    legend_marker_scale      scale of symbols in legend
    lw                       line width data curves (*default* is .5)
    ms                       size of markers (*default* is 1.5)
    lw_r                     line width response curves (*default* is .5)
    ms_r                     size of markers response curves (*default* is 1.5)
    mted                     marker for data Z_XX and Z_XY mode
    mtem                     marker for model Z_XX and Z_XY mode
    mtmd                     marker for data Z_YX and Z_YY mode
    mtmm                     marker for model Z_YX and Z_YY mode
    phase_limits             limits of phase
    plot_component           [ 2 | 4 ] 2 for TE and TM or 4 for all components
    plot_style               [ 1 | 2 ] 1 to plot each mode in a seperate
                             subplot and 2 to plot xx, xy and yx, yy in same
                             plots
    plot_type                [ '1' | list of station name ] '1' to plot all
                             stations in data file or input a list of station
                             names to plot if station_fn is input, otherwise
                             input a list of integers associated with the
                             index with in the data file, ie 2 for 2nd station
    plot_z                   [ True | False ] *default* is True to plot
                             impedance, False for plotting resistivity and
                             phase
    plot_yn                  [ 'n' | 'y' ] to plot on instantiation
    res_limits               limits of resistivity in linear scale
    resp_fn                  full path to response file
    resp_object              WSResponse object for resp_fn, or list of
                             WSResponse objects if resp_fn is a list of
                             response files
    station_fn               full path to station file written by WSStation
    subplot_bottom           space between axes and bottom of figure
    subplot_hspace           space between subplots in vertical direction
    subplot_left             space between axes and left of figure
    subplot_right            space between axes and right of figure
    subplot_top              space between axes and top of figure
    subplot_wspace           space between subplots in horizontal direction
    ======================== ==================================================
    """

    def __init__(self, data_fn=None, resp_fn=None, **kws):
        super ().__init__(data_fn= data_fn , resp_fn= resp_fn , **kws) 
        self._logging=csamtpylog.get_csamtpy_logger(self.__class__.__name__)
        
        for key in list(kws.keys()): 
            setattr(self, key, kws[key])
            
class PlotRMSMaps(plot_rms_maps.PlotRMSMaps):
    """
    plots the RMS as (data-model)/(error) in map view for all components
    of the data file.  Gets this infomration from the .res file output
    by ModEM.

    .. author: `MTpy`_ of Alison Kirkby avd Jared R. Peacock (
            https://github.com/MTgeophysics/pycsamt.git)
      
      
    Arguments:
    ----------
        **residual_fn** : string
                          full path to .res file

    =================== =======================================================
    Attributes                   Description
    =================== =======================================================
    fig                 matplotlib.figure instance for a single plot
    fig_dpi             dots-per-inch resolution of figure *default* is 200
    fig_num             number of fig instance *default* is 1
    fig_size            size of figure in inches [width, height]
                        *default* is [7,6]
    font_size           font size of tick labels, axis labels are +2
                        *default* is 8
    marker              marker style for station rms,
                        see matplotlib.line for options,
                        *default* is 's' --> square
    marker_size         size of marker in points. *default* is 10
    pad_x               padding in map units from edge of the axis to stations
                        at the extremeties in longitude.
                        *default* is 1/2 tick_locator
    pad_y               padding in map units from edge of the axis to stations
                        at the extremeties in latitude.
                        *default* is 1/2 tick_locator
    period_index        index of the period you want to plot according to
                        self.residual.period_list. *default* is 1
    plot_yn             [ 'y' | 'n' ] default is 'y' to plot on instantiation
    plot_z_list         internal variable for plotting
    residual            modem.Data instance that holds all the information
                        from the residual_fn given
    residual_fn         full path to .res file
    rms_cmap            matplotlib.cm object for coloring the markers
    rms_cmap_dict       dictionary of color values for rms_cmap
    rms_max             maximum rms to plot. *default* is 5.0
    rms_min             minimum rms to plot. *default* is 1.0
    save_path           path to save figures to. *default* is directory of
                        residual_fn
    subplot_bottom      spacing from axis to bottom of figure canvas.
                        *default* is .1
    subplot_hspace      horizontal spacing between subplots.
                        *default* is .1
    subplot_left        spacing from axis to left of figure canvas.
                        *default* is .1
    subplot_right       spacing from axis to right of figure canvas.
                        *default* is .9
    subplot_top         spacing from axis to top of figure canvas.
                        *default* is .95
    subplot_vspace      vertical spacing between subplots.
                        *default* is .01
    tick_locator        increment for x and y major ticks. *default* is
                        limits/5
    bimg                path to a geotiff to display as background of
                        plotted maps
    bimg_band           band of bimg to plot. *default* is None, which 
                        will plot all available bands
    bimg_cmap           cmap for bimg. *default* is 'viridis'. Ignored 
                        if bimg is RBG/A
    =================== =======================================================

    =================== =======================================================
    Methods             Description
    =================== =======================================================
    plot                plot rms maps for a single period
    plot_loop           loop over all frequencies and save figures to save_path
    read_residual_fn    read in residual_fn
    redraw_plot         after updating attributes call redraw_plot to
                        well redraw the plot
    save_figure         save the figure to a file
    =================== =======================================================


    :Example: ::

        >>> import pycsamt.modeling.modem as modem
        >>> rms_plot = PlotRMSMaps(r"/home/ModEM/Inv1/mb_NLCG_030.res")
        >>> # change some attributes
        >>> rms_plot.fig_size = [6, 4]
        >>> rms_plot.rms_max = 3
        >>> rms_plot.redraw_plot()
        >>> # happy with the look now loop over all periods
        >>> rms_plot.plot_loop()
    """
    def __init__(self, residual_fn, **kws):
        super ().__init__(residual_fn = residual_fn, **kws)
        self._logging=csamtpylog.get_csamtpy_logger(self.__class__.__name__)
        
        for key in list(kws.keys()): 
            setattr(self, key, kws[key])
 
class PlotPTMaps(phase_tensor_maps.PlotPTMaps):
    """
    Plot phase tensor maps including residual pt if response file is input.

      .. author: `MTpy`_ of Alison Kirkby avd Jared R. Peacock (
              https://github.com/MTgeophysics/pycsamt.git)

    :Plot only data for one period: ::

        >>> import pycsamt.modeling.ws3dinv as ws
        >>> dfn = r"/home/MT/ws3dinv/Inv1/WSDataFile.dat"
        >>> ptm = ws.PlotPTMaps(data_fn=dfn, plot_period_list=[0])

    :Plot data and model response: ::

        >>> import pycsamt.modeling.ws3dinv as ws
        >>> dfn = r"/home/MT/ws3dinv/Inv1/WSDataFile.dat"
        >>> rfn = r"/home/MT/ws3dinv/Inv1/Test_resp.00"
        >>> mfn = r"/home/MT/ws3dinv/Inv1/Test_model.00"
        >>> ptm = ws.PlotPTMaps(data_fn=dfn, resp_fn=rfn, model_fn=mfn,
        >>> ...                 plot_period_list=[0])
        >>> # adjust colorbar
        >>> ptm.cb_res_pad = 1.25
        >>> ptm.redraw_plot()


    ========================== ================================================
    Attributes                 Description
    ========================== ================================================
    cb_pt_pad                  percentage from top of axes to place pt
                               color bar. *default* is 1.2
    cb_res_pad                 percentage from bottom of axes to place
                               resistivity color bar. *default* is 0.5
    cb_residual_tick_step      tick step for residual pt. *default* is 3
    cb_tick_step               tick step for phase tensor color bar,
                               *default* is 45
    data_obj                   data object (read in from ModEM data file)
    data_fn                    full path to data fle
    dscale                     scaling parameter depending on map_scale
    ellipse_cmap               color map for pt ellipses. *default* is
                               mt_bl2gr2rd
    ellipse_colorby            [ 'skew' | 'skew_seg' | 'phimin' | 'phimax'|
                                 'phidet' | 'ellipticity' ] parameter to color
                                 ellipses by. *default* is 'phimin'
    ellipse_range              (min, max, step) min and max of colormap, need
                               to input step if plotting skew_seg
    ellipse_size               relative size of ellipses in map_scale
    ew_limits                  limits of plot in e-w direction in map_scale
                               units.  *default* is None, scales to station
                               area
    fig_aspect                 aspect of figure. *default* is 1
    fig_dpi                    resolution in dots-per-inch. *default* is 300
    fig_list                   list of matplotlib.figure instances for each
                               figure plotted.
    fig_size                   [width, height] in inches of figure window
                               *default* is [6, 6]
    font_size                  font size of ticklabels, axes labels are
                               font_size+2. *default* is 7
    grid_east                  relative location of grid nodes in e-w direction
                               in map_scale units
    grid_north                 relative location of grid nodes in n-s direction
                               in map_scale units
    grid_z                     relative location of grid nodes in z direction
                               in map_scale units
    model_fn                 full path to initial file
    map_scale                  [ 'km' | 'm' ] distance units of map.
                               *default* is km
    mesh_east                  np.meshgrid(grid_east, grid_north, indexing='ij')
    mesh_north                 np.meshgrid(grid_east, grid_north, indexing='ij')
    model_fn                   full path to model file
    nodes_east                 relative distance betwen nodes in e-w direction
                               in map_scale units
    nodes_north                relative distance betwen nodes in n-s direction
                               in map_scale units
    nodes_z                    relative distance betwen nodes in z direction
                               in map_scale units
    ns_limits                  (min, max) limits of plot in n-s direction
                               *default* is None, viewing area is station area
    pad_east                   padding from extreme stations in east direction
    pad_north                  padding from extreme stations in north direction
    period_list                list of periods from data
    plot_grid                  [ 'y' | 'n' ] 'y' to plot grid lines
                               *default* is 'n'
    plot_period_list           list of period index values to plot
                               *default* is None
    plot_yn                    ['y' | 'n' ] 'y' to plot on instantiation
                               *default* is 'y'
    res_cmap                   colormap for resisitivity values.
                               *default* is 'jet_r'
    res_limits                 (min, max) resistivity limits in log scale
                               *default* is (0, 4)
    res_model                  np.ndarray(n_north, n_east, n_vertical) of
                               model resistivity values in linear scale
    residual_cmap              color map for pt residuals.
                               *default* is 'mt_wh2or'
    resp                       np.ndarray(n_stations, n_periods, 2, 2)
                               impedance tensors for model response
    resp_fn                    full path to response file
    save_path                  directory to save figures to
    save_plots                 [ 'y' | 'n' ] 'y' to save plots to save_path
    station_east               location of stations in east direction in
                               map_scale units
    station_fn                 full path to station locations file
    station_names              station names
    station_north              location of station in north direction in
                               map_scale units
    subplot_bottom             distance between axes and bottom of figure window
    subplot_left               distance between axes and left of figure window
    subplot_right              distance between axes and right of figure window
    subplot_top                distance between axes and top of figure window
    title                      titiel of plot *default* is depth of slice
    xminorticks                location of xminorticks
    yminorticks                location of yminorticks
    ========================== ================================================
    """ 
    def __init__(self, data_fn=None, resp_fn=None, model_fn=None, **kws):
        super ().__init__( data_fn=data_fn, resp_fn=resp_fn, model_fn=model_fn,
                          **kws)
        self._logging=csamtpylog.get_csamtpy_logger(self.__class__.__name__)
        
        for key in list(kws.keys()): 
            setattr(self, key, kws[key])
            
class PlotSlices(plot_slices.PlotSlices):
    """
    * Plot all cartesian axis-aligned slices and be able to scroll through the model
    * Extract arbitrary profiles (e.g. along a seismic line) from a model

    .. author: `MTpy`_ of Alison Kirkby avd Jared R. Peacock (
            https://github.com/MTgeophysics/pycsamt.git)
    
    :Example: ::

        >>> import pycsamt.modeling.modem as modem
        >>> mfn = r"/home/modem/Inv1/Modular_NLCG_100.rho"
        >>> dfn = r"/home/modem/Inv1/ModEM_data.dat"
        >>> pds = ws.PlotSlices(model_fn=mfn, data_fn=dfn)

    ======================= ===================================================
    Buttons                  Description
    ======================= ===================================================
    'e'                     moves n-s slice east by one model block
    'w'                     moves n-s slice west by one model block
    'n'                     moves e-w slice north by one model block
    'm'                     moves e-w slice south by one model block
    'd'                     moves depth slice down by one model block
    'u'                     moves depth slice up by one model block
    ======================= ===================================================


    ======================= ===================================================
    Attributes              Description
    ======================= ===================================================
    ax_en                   matplotlib.axes instance for depth slice  map view
    ax_ez                   matplotlib.axes instance for e-w slice
    ax_map                  matplotlib.axes instance for location map
    ax_nz                   matplotlib.axes instance for n-s slice
    climits                 (min , max) color limits on resistivity in log
                            scale. *default* is (0, 4)
    cmap                    name of color map for resisitiviy.
                            *default* is 'jet_r'
    data_fn                 full path to data file name
    draw_colorbar           show colorbar on exported plot; default True
    dscale                  scaling parameter depending on map_scale
    east_line_xlist         list of line nodes of east grid for faster plotting
    east_line_ylist         list of line nodes of east grid for faster plotting
    ew_limits               (min, max) limits of e-w in map_scale units
                            *default* is None and scales to station area
    fig                     matplotlib.figure instance for figure
    fig_aspect              aspect ratio of plots. *default* is 1
    fig_dpi                 resolution of figure in dots-per-inch
                            *default* is 300
    fig_num                 figure instance number
    fig_size                [width, height] of figure window.
                            *default* is [6,6]
    font_dict               dictionary of font keywords, internally created
    font_size               size of ticklables in points, axes labes are
                            font_size+2. *default* is 4
    grid_east               relative location of grid nodes in e-w direction
                            in map_scale units
    grid_north              relative location of grid nodes in n-s direction
                            in map_scale units
    grid_z                  relative location of grid nodes in z direction
                            in map_scale units
    index_east              index value of grid_east being plotted
    index_north             index value of grid_north being plotted
    index_vertical          index value of grid_z being plotted
    initial_fn              full path to initial file
    key_press               matplotlib.canvas.connect instance
    map_scale               [ 'm' | 'km' ] scale of map. *default* is km
    mesh_east               np.meshgrid(grid_east, grid_north)[0]
    mesh_en_east            np.meshgrid(grid_east, grid_north)[0]
    mesh_en_north           np.meshgrid(grid_east, grid_north)[1]
    mesh_ez_east            np.meshgrid(grid_east, grid_z)[0]
    mesh_ez_vertical        np.meshgrid(grid_east, grid_z)[1]
    mesh_north              np.meshgrid(grid_east, grid_north)[1]
    mesh_nz_north           np.meshgrid(grid_north, grid_z)[0]
    mesh_nz_vertical        np.meshgrid(grid_north, grid_z)[1]
    model_fn                full path to model file
    ms                      size of station markers in points. *default* is 2
    nodes_east              relative distance betwen nodes in e-w direction
                            in map_scale units
    nodes_north             relative distance betwen nodes in n-s direction
                            in map_scale units
    nodes_z                 relative distance betwen nodes in z direction
                            in map_scale units
    north_line_xlist        list of line nodes north grid for faster plotting
    north_line_ylist        list of line nodes north grid for faster plotting
    ns_limits               (min, max) limits of plots in n-s direction
                            *default* is None, set veiwing area to station area
    plot_yn                 [ 'y' | 'n' ] 'y' to plot on instantiation
                            *default* is 'y'
    plot_stations           default False
    plot_grid               show grid on exported plot; default False
    res_model               np.ndarray(n_north, n_east, n_vertical) of
                            model resistivity values in linear scale
    save_format             exported format; default png
    save_path               path to save exported plots to; default current working folder
    station_color           color of station marker. *default* is black
    station_dict_east       location of stations for each east grid row
    station_dict_north      location of stations for each north grid row
    station_east            location of stations in east direction
    station_fn              full path to station file
    station_font_color      color of station label
    station_font_pad        padding between station marker and label
    station_font_rotation   angle of station label
    station_font_size       font size of station label
    station_font_weight     weight of font for station label
    station_id              [min, max] index values for station labels
    station_marker          station marker
    station_names           name of stations
    station_north           location of stations in north direction
    subplot_bottom          distance between axes and bottom of figure window
    subplot_hspace          distance between subplots in vertical direction
    subplot_left            distance between axes and left of figure window
    subplot_right           distance between axes and right of figure window
    subplot_top             distance between axes and top of figure window
    subplot_wspace          distance between subplots in horizontal direction
    title                   title of plot
    xminorticks             location of xminorticks
    yminorticks             location of yminorticks
    z_limits                (min, max) limits in vertical direction,
    ======================= ===================================================

    """

    def __init__(self, model_fn, data_fn=None, **kws):     
        super ().__init__( model_fn=model_fn, data_fn=data_fn, **kws)
        self._logging=csamtpylog.get_csamtpy_logger(self.__class__.__name__)
        
        for key in list(kws.keys()): 
            setattr(self, key, kws[key])
                 

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    