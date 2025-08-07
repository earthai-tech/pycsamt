# -*- coding: utf-8 -*-
#   License: BSD-3-Clause
#   Author: LKouadio <etanoyau@gmail.com>

"""
:mod:`~watex.utils.plot` is a set of base plots for :term:`tensor` 
visualization, data exploratory and analyses. 
T-E-Q Plots encompass the tensors plots (:class:`~watex.view.TPlot`) dealing 
with :term:`EM` methods, Exploratory plots ( :class:`~watex.view.ExPlot`) and 
Quick analyses (:class:`~watex.view.QuickPlot`) visualization. 
"""
from __future__ import annotations 

import re 
import copy
import warnings
import itertools
import numpy as np 
import  matplotlib.pyplot  as plt
import matplotlib.ticker as mticker
from matplotlib.gridspec import GridSpec 
import pandas as pd 
from pandas.plotting import ( 
    radviz , 
    parallel_coordinates
    ) 
import seaborn as sns 
from .._docstring import ( 
    DocstringComponents,
    _core_docs,
    _baseplot_params
    )

from ..decorators import temp2d 

from ..exceptions import ( 
    PlotError, 
    NotFittedError, 
    EMError, 
    SiteError, 
    )
from .base import BasePlot
from .._typing import (
    List,
    Optional,
    EDIO
)
from ..utils._dependency import  import_optional_dependency 

from ..utils.exmath import ( 
    moving_average , fittensor)
from ..utils.funcutils import ( 
    _assert_all_types , 
    _validate_name_in, 

    remove_outliers, 
    smart_format,
    reshape, 

    is_iterable, 
    station_id, 
    make_ids,
    )


from ..utils.plot import(
    make_mpl_properties, 
     plot_errorbar
     )

from ..methods.em import EMAP, MT  



class TPlot (BasePlot): 

    _t= (
        "survey_area",
        "distance",
        "prefix",
        "window_size",
        "component",
        "mode",
        "method",
        "out",
        "how",
        "c" 
        )
    
    def __init__ (
        self, 
        survey_area =None , 
        distance = 50., 
        prefix ='S', 
        how= 'py',
        window_size:int =5, 
        component:str ='xy', 
        mode: str ='same', 
        method:str ='slinear', 
        out:str  ='srho', 
        c: int =2,
        **kws
        ): 
        super().__init__(**kws)
        
        self.survey_area=survey_area 
        self.distance=distance 
        self.prefix=prefix
        self.window_size=window_size
        self.component=component 
        self.mode=mode
        self.method=method
        self.out=out
        self.how=how
        self.c=c 
        
    def fit (
        self, 
        data: Optional [str|List[EDIO]]
        ): 
        """
        Fit data and populate attributes. 
        
        Parameters 
        ----------- 
        data : str, or list or :class:`pycsamt.core.edi.Edi` object 
            Full path to EDI files or collection of EDI-objects 
   
        Returns
        -------- 
        ``self``: :class:`watex.view.plot.TPlot` instanciated object
            returns ``self`` for chaining methods.
        
        """

        p = EMAP(
            window_size = self.window_size ,  
            component= self.component, 
            mode= self.mode, 
            method= self.method, 
            out=self.out, 
            c=self.c 
            ) 
        p.fit(data)
        # set EM processing module 
        # as an attr
        setattr (self, "p_", p )
        
        # set component slices into a dict
        self._c_= {
              'xx': [slice (None, len(self.p_.freqs_)), 0 , 0] , 
              'xy': [slice (None, len(self.p_.freqs_)), 0 , 1], 
              'yx': [slice (None, len(self.p_.freqs_)), 1 , 0], 
              'yy': [slice (None, len(self.p_.freqs_)), 1,  1] 
        }
        return self
    
    @property 
    def inspect (self): 
        """ Inspect object whether is fitted or not"""
        msg = ( "{obj.__class__.__name__} instance is not fitted yet."
               " Call 'fit' with appropriate arguments before using"
               " this method"
               )
        
        if not hasattr (self, 'p_'): 
            raise NotFittedError(msg.format(
                obj=self)
            )
        return 1 
    
    
    def plot_multi_recovery(
            self,  
            sites:str |List[str | int], 
            colors: List[str] = None, 
            **kws
            ): 
        
        """ 
        Plots mutiple site/stations with signal recovery. 
        
        Parameters 
        -----------
        sites: list 
            list of sites to visualize. Can also be the index of the sites 
        colors: list of str
            matplotlib colors to customize the raw signal and recovery signal
     
        Returns 
        ----------
        ax: Matplotlib suplot axes 
        
        Examples
        --------
        >>> from watex.view.plot import TPlot 
        >>> from watex.datasets import load_edis 
        >>> # takes the 03 samples of EDIs 
        >>> edi_data = load_edis (return_data= True, samples =3 ) 
        >>> TPlot(fig_size =(5, 3)).fit(edi_data).plot_multi_recovery (
            sites =['S00'], colors =['o', 'ok--'])
        <AxesSubplot:title={'center':'Recovered tensor $|Z_{xy}|$'}, 
        xlabel='$Frequency [H_z]$', ylabel='$ App.resistivity \\quad xy \\quad [ \\Omega.m]$'>
        """
        self.inspect 

        if isinstance (sites, str): 
            sites =[sites ] 
        if not is_iterable(sites): 
            sites =[sites] 
       
        site_index = station_id(sites) 
        for i, s in enumerate (site_index): 
            if s > len(self.p_.ediObjs_): 
                raise PlotError(f"Site {sites[i]!r} is out of the expected"
                                f" sites number: {len(self.p_.ediObjs_)}"
                                )
        # read the component XY 
        res2d_r = self.p_.make2d (out=f'res{self.component}') 
        z_xy_rest = self.p_.zrestore() # no buffered data 
        # extracted the station at index 12, 27 for instance.
        zs = [ z_xy_rest[s].resistivity[
            tuple (self._c_.get(self.component))] for s in site_index ]

        ma = [ moving_average ( res2d_r[:,  s_ix  ]) for s_ix in site_index ]

        f= self.p_.getfullfrequency()
        #>>> # ---> make a plot and color 
        # colors = make_mpl_properties(2*len(ma))
        if colors is None: 
            colors =[]
        if isinstance (colors, str): 
            colors =[colors]
        colors +=  make_mpl_properties(2*len(ma))
        
        fs = [f for i in range(len(ma))] # repeat frequency 
        z_norm_args = list( zip (fs, zs, colors[: len(ma)] )) 
        args  = list(itertools.chain(*z_norm_args))
        # >>> #  make a fitting args 
        colors = ['m-'] + colors[len(ma):]
        z_cor_objs = list( zip (fs, ma, ['m-'] + colors[len(ma):] )) 
        fit_args = list(itertools.chain(*z_cor_objs))
        
        xlim = (f.min() -.5 * f.min(), f.max() +.5 * f.max())
        return self._plot_recovery (
            *args, fit_args= fit_args, xlim=xlim, sites = sites,  **kws )

    def _plot_recovery (
            self,
            *args,  
            fit_args = None, 
            leg= None,  
            xlim=None, 
            sites=None, 
            **kws
            ): 
        """" Template to plot two stations with signal recovery 
        
        Isolated part of :meth:`~.TPlot.plot_multi_recovery`. 
        
        Parameters 
        -----------
        *args : args : list 
            Matplotlib logs funtions plot arguments 
            
        fit_args : list or tuple 
            Matplotlib logs funtions plot arguments put on list. It used to 
            visualize the fitting curve after apply anay correction to the data.
            
            X-coordinates. It should have the length M, the same of the ``arr2d``; 
            the columns of the 2D dimensional array.  Note that if `x` is 
            given, the `distance is not needed. 
            
        leg: list 
            legend labels put into a list. It must fit the number of given 
            plots. 
             
        kws : dict 
            Additional keywords arguments of Matplotlib subsplots function  
            :func:`plt.loglog` or :func:`plt.semilog`
        
        Returns 
        ------- 
        ax: Matplotlib.pyplot <AxesSubplot>  

        """
        fig, ax = plt.subplots(
            1,figsize = self.fig_size, 
            #num = self.fig_num,
             )
        p1= ax.loglog(*args,  
                  markersize = self.ms ,
                  linewidth = self.lw ,
                  **kws 
                  )
        p2 =[]
        if fit_args  is not None: 
            fit_args  = _assert_all_types(
                fit_args , list, tuple, objname="Fit arguments")  
            p2 = ax.loglog(*fit_args,
                      markersize = self.ms ,
                      linewidth = self.lw 
                      )

        ax.set_xlabel (self.xlabel or '$Frequency [H_z]$',
                    fontsize =1.5 * self.font_size ) 
        ax.set_ylabel(self.ylabel or '$ App.resistivity \quad xy \quad [ \Omega.m]$',
                   fontsize =1.5*self.font_size)
        
        p1labels= [f'rec.tensor {i}' for i in sites ]
        p2labels= [f'mov-aver. line {i}' for i in sites 
                   ] if fit_args is not None else []
        
        ax.legend (handles = [*p1 ,*p2], 
                   labels= [*p1labels, *p2labels] #['restored data' , 'recovery trend '] 
                   if leg is  None else leg,
                   loc ='best', 
                   # ncol =len(args)//3 if fit_args  is None else (
                   #      (len(args)+len(fit_args )))//3  ,
                    fontsize =1.5 * self.font_size
                    )
    
        if xlim is not None: 
            ax.set_xlim (xlim)
        ax.tick_params (axis= 'both', labelsize = 1.5 * self.font_size )
        plt.title (self.fig_title or  'Recovered tensor $|Z_{xy}|$',
                   fontsize =1.5*self.font_size)
        
        
        if self.show_grid :
            ax.grid (visible =True , alpha =self.galpha,
                     which =self.gwhich, color =self.gc)
            
        if self.savefig is not None: 
             plt.savefig(self.savefig , dpi = self.fig_dpi)
             plt.close (fig =fig ) 

        return ax
    
    @temp2d("Base template for 2D recovery tensors plot.")
    def plot_tensor2d (
        self,    
        tensor='res', 
        sites =None, 
        to_log10=False, 
        ): 
        """ Plot two dimensional tensor. 
         
        Parameters 
        -----------
        freqs: array-like 
            y-coordinates. It should have the length N, the same of the ``arr2d``.
            the rows of the ``arr2d``.Frequency array. It should be the 
            complete frequency used during the survey area. 

        tensor: str , ['res','phase', 'z'], default='res'
            kind of tensor to plot. Can be resistivity or phase. If `phase`, 
            customize your plot to not fit the default 'res' behaviour. 
        to_log10: bool, defaut=False, 
            Convert the resistivity data and frequeny in log10.  
            
        sites: list of str, optional 
            List of stations/sites names. If given, it must have the same 
            length of the positions in of the EDI data. Must fit the number 
            of 'EDI' succesffully read. 

        Returns 
        -------
        ( arr2d , freqs, positions , sites , base_plot_kws): 
            - arr2d: 2D resistivity array from the tensor `component` 
            - freqs: array-like 1d of frequency in the survey.
            - positions: Sites/stations positions. It is equals to the distance
                between stations times the number of sites 
            - sites: list of the names of the station/sites 
            - base_plot_kws: plot keywords arguments inherits from 
                :class:`watex.property.BasePlot`. It composes the last 
                parameters for customizing plot as decorated return function. 
    
        Examples 
        -------- 
        >>> from watex.view.plot import TPlot 
        >>> from watex.datasets import load_edis 
        >>> # get some 3 samples of EDI for demo 
        >>> edi_data = load_edis (return_data =True, samples =3 )
        >>> # customize plot by adding plot_kws 
        >>> plot_kws = dict( ylabel = '$Log_{10}Frequency [Hz]$', 
                            xlabel = '$Distance(m)$', 
                            cb_label = '$Log_{10}Rhoa[\Omega.m$]', 
                            fig_size =(6, 3), 
                            font_size =7. 
                            ) 
        >>> t= TPlot(**plot_kws ).fit(edi_data)
        >>> # plot recovery2d using the log10 resistivity 

        >>> t.plot_tensor2d (to_log10=True)
        <AxesSubplot:xlabel='$Distance(m)$', ylabel='$Log_{10}Frequency [Hz]$'>
 
        """
        
        self.inspect 
        
        assert  str(tensor).lower() in {"res", 'phase'}, (
            "Expect either a resistivity 'res' or 'phase'. Got {tensor!r}")
        tensor =str(tensor).lower() 
        
        arr2d = self.p_.make2d (out = f'{tensor}{self.component}') 

        return self._make_tensor_utils (arr2d, sites , to_log10, tensor )  
    
    @temp2d("Base template for 2D filtered tensors plot.")
    def plot_ctensor2d  (
            self, 
            tensor ='res' , 
            ffilter ='tma', 
            sites = None, 
            to_log10=False
            ): 
        """ Plot filtered tensors 
        
        Parameters 
        -----------
        tensor: str , ['res','phase', 'z'], default='res'
            kind of tensor to plot. Can be resistivity or phase. If `phase`, 
            customize your plot to not fit the default 'res' behaviour.
            
        ffilter: str ['ama', 'flma', 'tma'], default='tma'
            kind of appropriate filter to corrected tensor data. 
            
        to_log10: bool, defaut=False, 
            Convert the resistivity data and frequeny in log10.  
            
        sites: list of str, optional 
            List of stations/sites names. If given, it must have the same 
            length of the positions in of the EDI data. Must fit the number 
            of 'EDI' succesffully read. 
            
        Returns 
        -------
        ( arr2d , freqs, positions , sites , base_plot_kws): 
            - arr2d: 2D filtered tensor array from the `component` 
            - freqs: array-like 1d of frequency in the survey.
            - positions: Sites/stations positions. It is equals to the distance
                between stations times the number of sites 
            - sites: list of the names of the station/sites 
            - base_plot_kws: plot keywords arguments inherits from 
                :class:`watex.property.BasePlot`. It composes the last 
                parameters for customizing plot as decorated return function.
                
        Examples 
        -------- 
        >>> from watex.view.plot import TPlot 
        >>> from watex.datasets import load_edis 
        >>> # get some 3 samples of EDI for demo 
        >>> edi_data = load_edis (return_data =True, samples =3 )
        >>> # customize plot by adding plot_kws 
        >>> plot_kws = dict( ylabel = '$Log_{10}Frequency [Hz]$', 
                            xlabel = '$Distance(m)$', 
                            cb_label = '$Log_{10}Rhoa[\Omega.m$]', 
                            fig_size =(6, 3), 
                            font_size =7. 
                            ) 
        >>> t= TPlot(**plot_kws ).fit(edi_data)
        >>> # plot filtered tensor using the log10 resistivity 
        >>> t.plot_ctensor2d (to_log10=True)
        <AxesSubplot:xlabel='$Distance(m)$', ylabel='$Log_{10}Frequency [Hz]$'>
  
        """
        self.inspect 
        fd = {"tma": self.p_.tma , "flma":self.p_.flma, "ama":self.p_.ama }

        assert str(ffilter).lower() in fd.keys(), (
            "Supports only base filters {tuple (fd.keys())}. Got {ffilter!r}"
            " To apply a simple filter like 'moving average' to a tensor, refer"
            " to <watex.utils.exmath.moving_average>. For other filters like"
            " 'Savitzky Golay1d/2d', 'remove distorsion' or 'remove outliers'"
            " and else, use the package 'pycsamt' instead. "
            ) 
        ffilter= str (ffilter).lower().strip()         
        arr2d = fd.get(ffilter)()
    
        return self._make_tensor_utils (arr2d, sites, to_log10 , tensor ) 
    
    
    def _make_tensor_utils (
            self, arr2d, sites, to_log10= False, tensor=None ): 
        """ Make utilities for plotting tensors   
        
        Parameters 
        ------------
        arr2d: arraylike of shape (n_freq, n_sites): 
            Array of the tensor composed of number of frequency and number 
            of sites that fit the number of EDI correctly read.
        
        sites: list of str, optional 
            List of stations/sites names. If given, it must have the same 
            length of the positions in of the EDI data. Must fit the number 
            of 'EDI' succesffully read. 
        to_log10: bool, defaut=False, 
            Convert the resistivity data and frequeny in log10.
            
        Returns 
        -------
        ( arr2d , freqs, positions , sites , base_plot_kws): 
            - arr2d: 2D filtered tensor array from the `component` 
            - freqs: array-like 1d of frequency in the survey.
            - positions: Sites/stations positions. It is equals to the distance
                between stations times the number of sites 
            - sites: list of the names of the station/sites 
            - base_plot_kws: plot keywords arguments inherits from 
                :class:`watex.property.BasePlot`. It composes the last 
                parameters for customizing plot as decorated return function.
        """
        try : 
            distance = float(self.distance) 
        except : 
            raise TypeError (
                f'Expect a float value not {type(self.distance).__name__!r}')

        freqs = self.p_.freqs_ 

        positions = np.arange(arr2d.shape[1])  * distance
            
        sites = sites or make_ids (
            positions , self.prefix , how = self.how)  
        
        if isinstance(sites, str): 
            sites =[sites] 
        if not is_iterable(sites): 
            raise TypeError("Sites collection must be an iterable" 
                            f" object. Got {type(sites).__name__!r}"
                    )
        if len(sites)!= len(positions): 
            raise TypeError (f"Sites={len(sites)} length must be consistent."
                             "  Expects positions={len(positions)}.")
            
        if tensor in {'phase', 'phs'}: 
            arr2d %=90
            
        if to_log10: 
            arr2d = arr2d if tensor in ("phase", "phs") else np.log10 (arr2d) 
            freqs = np.log10 (freqs)
            
        base_plot_kws = {
            k: v for k, v in self.__dict__.items () 
            if k not in list(self._t ) +['p_']
            }  
        
        return arr2d, freqs, positions ,sites , base_plot_kws  

    def plot_recovery(self, site = 'S00'): 
        """ visualize the restored tensor per site.
        
        Parameters 
        ------------
        site: str, int, default ="S00"
            Site/station name for 
        
        Returns
        -------- 
        ``self``: :class:`watex.view.plot.TPlot` instanciated object
            returns ``self`` for chaining methods.
            
        Examples 
        --------
        >>> from watex.view import TPlot
        >>> from watex.datasets import load_edis 
        >>> edi_data = load_edis (return_data =True, samples =7) 
        >>> plot_kws = dict( ylabel = '$Log_{10}Frequency [Hz]$', 
                    xlabel = '$Distance(m)$', 
                    cb_label = '$Log_{10}Rhoa[\Omega.m$]', 
                    fig_size =(7, 4), 
                    font_size =7. 
                    ) 
        >>> t= TPlot(**plot_kws ).fit(edi_data)
        >>> # plot recovery of site 'S01'
        >>> t.plot_recovery ('S01')
        
        """
        self.inspect 
        
        if isinstance(site, str): 
            site =[site]
        site_index = station_id(site) 
        
        site_index = site_index [0] if isinstance (
            site_index, tuple ) else site_index 
 
        if site_index  > len(self.p_.ediObjs_): 
            raise PlotError(f"Site {site!r} is out of the expected"
                            f" sites number: {len(self.p_.ediObjs_)}"
                            )
    
        ediObjs = self.p_.ediObjs_ 
        # >>> zobjs_b = restoreZ(ediObjs, buffer = buffer) # with buffer 
        zobjs = self.p_.zrestore() # with no buffer 
    
        zxy_restored = np.abs (zobjs[site_index].z [
            tuple (self._c_.get(self.component))])#[:, 0, 1]) 
        # Export the first raw object with missing Z at 
        # the dead dand in ediObjs collection
        z1 = np.abs(ediObjs[site_index].Z.z) 
        z1freq= ediObjs[site_index].Z._freq # the frequency of the first z obj 
        # get the frequency of the clean data knonw as reference frequency
        indice_reffreq = np.argmax (list (map(lambda o: len(o.Z._freq), ediObjs)))
        reffreq = ediObjs [indice_reffreq].Z._freq 
        # >>> # use the real part of component xy for the test 
        zxy = np.abs (z1[tuple (self._c_.get(self.component))])  #[:, 0, 1])  
        # fit zxy to get the missing data in the dead band 
        zfit = fittensor(refreq= reffreq, compfreq= z1freq, z=zxy)

        # not necessary, one can corrected z to get a 
        # smooth resistivity distribution 
        zcorrected = moving_average (zxy_restored)                     
        # plot the two figures 
        plt.figure(figsize =self.fig_size) #(10, 5)
        plt.loglog(reffreq, zfit, '^r', reffreq, zxy_restored, 'ok--')
        plt.loglog( reffreq, zcorrected, '1-.')
        plt.legend (['raw data', 'tensor $res_{xy}$ restored',
                     'moving-average trend line' ],loc ='best')
        plt.xlabel ('$Frequency [H_z]$') 
        plt.ylabel('$ Resistivity_{xy} [ \Omega.m]$')
        plt.title ('Recovered tensor $|Z_{xy}|$' + f" at site {site[0].upper()}")
        plt.grid (visible =True , alpha =0.8, which ='both', color ='k')
        plt.xlim (reffreq.min() -.5* reffreq.min(), 
                  reffreq.max() + .5 * reffreq.max())
        plt.tight_layout()
        
        return self 
    
    def plot_phase_tensors(
        self,
        mode ='frequency',
        stretch = (7000, 20 ),
        linedir ='ns',
        tensor='phimin',
        ellipse_dict = None,
        **kws
        ): 
        """ Plot phase tensor pseudosection and skew ellipsis 
        visualization. 
        
        Method plots the phase tensor ellipses in a pseudo section format.
        It uses `mtpy` as dependency. 
        
        Parameters 
        -----------
        mode: str, default ='frequency'
            Tempoora scale in y-axis. Can be ['frequency' | 'period']

        stretch : float or tuple (xstretch, ystretch), default=200
            Is a factor that scales the distance from one station to the next 
            to make the plot readable. It determines (x,y) aspect ratio of plot.
    
        linedir: str [ 'ns' | 'ew' ], default='ns'
            The predominant direction of profile line. It can be ['ns' | 'ew']
            where: 
               
            * 'ns' refer to North-South Line or line is closer to north-south)
            * 'ew' refer to  East-West line or line is closer to east-west
            *Default* is 'ns'
        tensor: str, default='phimin' 
            Is the tensor skew or ellipsis visualizations. The color for plot 
            style is referred accordingly. Tensor can be: 
                
                [ 'phimin' | 'phimax' | 'skew' |'skew_seg' | 'phidet' |'ellipticity' ]
           where: 
                  
                - 'phimin' -> colors by minimum phase
                - 'phimax' -> colors by maximum phase
                - 'skew' -> colors by skew
                - 'skew_seg' -> colors by skew indiscrete segments defined 
                   by the range
                - 'normalized_skew' -> colors by skew see [Booker, 2014]
                - 'normalized_skew_seg' -> colors by normalized skew in
                   discrete segments defined by the range
                - 'phidet' -> colors by determinant of the phase tensor
                - 'ellipticity' -> colors by ellipticity *default* is 'phimin'  
                
        ellipse_dict: dict, optional
            Dictionary of parameters for the phase tensor ellipses with keys:
            
            * 'size': float, default =2 , is the size of ellipse in points
            * 'colorby' : str, default='phimin' 
               Is the color for plot style referring either to  tensor, 
               skew or ellipsis visualizations. It can be all the `tensor`
               parameter values. see `tensor` parameter values. 
               [ 'phimin' | 'phimax' | 'skew' |'skew_seg' | 'phidet' |'ellipticity' ]
        
            * 'range' : tuple (min, max, step), default='colorby'
               Need to input at least the min and max  and if using 
               'skew_seg' to plot discrete values input step as well
               
            * 'cmap' : [ 'mt_yl2rd' | 'mt_bl2yl2rd' |'mt_wh2bl' | 'mt_rd2bl' |
                        'mt_bl2wh2rd' | 'mt_seg_bl2wh2rd' |'mt_rd2gr2bl' ]

                     - 'mt_yl2rd' -> yellow to red
                     - 'mt_bl2yl2rd' -> blue to yellow to red
                     - 'mt_wh2bl' -> white to blue
                     - 'mt_rd2bl' -> red to blue
                     - 'mt_bl2wh2rd' -> blue to white to red
                     - 'mt_bl2gr2rd' -> blue to green to red
                     - 'mt_rd2gr2bl' -> red to green to blue
                     - 'mt_seg_bl2wh2rd' -> discrete blue to white to red
        kws: dict 
            Additional keywords arguments passed from |MTpy| pseudosection 
            phase tensor class: :class:`~.PlotPhaseTensorPseudoSection` 

        See Also
        ----------
        mtpy.imaging.phase_tensor_pseudosection.PlotPhaseTensorPseudoSection: 
            PlotPhase pseudo section tensor from |MTpy| package. 
        watex.utils.plot_skew: 
            Phase sensitive skew visualization. 
        
        Examples
        ---------
        >>> import watex as wx 
        >>> edi_data = wx.fetch_data ('edis', key='edi', return_data =True , samples =17 ) 
        >>> tplot = wx.methods.TPlot ().fit(edi_data ) 
        >>> tplot.plot_phase_tensors (tensor ='skew')
        
        """
        extra =("Phase tensor plots or skew ellipsis visualization"
                " uses 'mtpy' as dependency. Alternatively, you may"
                " also use the phase sensitive 'skew' visualization"
                " of plot utilities if plot  only refers to 'skew'."
                )
        import_optional_dependency ('mtpy' , extra = extra )
        from mtpy.imaging.phase_tensor_pseudosection import (
            PlotPhaseTensorPseudoSection ) 
        
        self.inspect 
        
        zobjs = [edi_obj.Z for edi_obj in self.p_.ediObjs_]
        
        elrange =  [-7, 7] if 'skew' in str(tensor).lower()  else [0, 90 ]  
        ellipse_dict = ellipse_dict or  {
            'ellipse_colorby':tensor,
            'ellipse_range':elrange,  # Color limits
            'ellip_size': 2, 
            'ellipse_cmap':'mt_bl2wh2rd'
        } 
        # skew_seg need to provide
        # 3 numbers, the 3rd indicates
        # interval, e.g. [-12,12,3]
        #from contextlib import suppress 
        # suppress as possible the external 
        #lib resources
        #with suppress (Exception): 
        ptsection = PlotPhaseTensorPseudoSection(
                        fn_list = self.p_.edifiles,
                        z_object_list = zobjs, 
                        fig_size = self.fig_size, 
                        tscale = mode, 
                        plot_num = self.fig_num, 
                        plot_title = self.fig_title, 
                        xlimits = self.xlim, 
                        ylimits = self.ylim,
                        linedir= linedir,  
                        stretch= stretch, 
                        station_id=(0, len(self.p_.ediObjs_)), 
                        font_size=self.font_size ,
                        lw=self.lw,
                        **ellipse_dict, 
                        **kws,
            )

        ptsection.save_figure(save_fn =self.savefig, fig_dpi=self.fig_dpi
                              ) if self.savefig else  ptsection.plot()

        return self 
    
    def plotSkew (
        self , 
        method ='Bahr', 
        sensitivity ='skew', 
        mode=None,
        threshold_line=None, 
        show_average_sensistivity=True, 
        suppress_outliers =True, 
        **plot_kws 
        ): 
        """ Plot phase sensistive skew visualization
        
        'Skew' is also knwown as the conventional asymmetry parameter 
        based on the Z magnitude. 
        
        Mosly, the :term:`EM` signal is influenced by several factors such 
        as the dimensionality of the propagation medium and the physical 
        anomalies, which can distort theEM field both locally and regionally. 
        The distortion of Z was determined from the quantification of its 
        asymmetry and the deviation from the conditions that define its 
        dimensionality. The parameters used for this purpose are all rotational 
        invariant because the Z components involved in its definition are
        independent of the orientation system used. The conventional asymmetry
        parameter based on the Z magnitude is the skew defined by Swift (1967)
        [1]_ and Bahr (1991) [2]_.

        Parameters 
        -----------
        method: str, default='Bahr': 
            Kind of correction. Can be:
                
            - ``swift`` for the remove distorsion proposed by Swift in 1967. 
              The value close to 0. assume the 1D and 2D structures, and 3D 
              otherwise.  However, In general case, the  electrical structure 
              of :math:`\eta < 0.4` can be treated as a 2D medium.
            - ``bahr`` for the remove distorsion proposed  by Bahr in 1991. 
              The latter threshold is set to 0.3. Above this value the 
              structures is 3D.
              
        sensitivity: str, default='skew'
           phase sensistive visualization. Can be rotational invariant 
           ``invariant``. In fact, setting to ``mu`` or ``invariant`` does 
           not change any interpretation when since the  distortion of Z 
           are all rotational invariant whether using the ``Bahr`` or ``swift``
           methods. 
           
           .. versionchanged:: 
               Param `view` is deprecated and replaced with `sensistivity`. 
               
        mode:str, optional 
           X-axis coordinates for visualisation. plot either ``'frequency'`` or
           ``'periods'``.  The default is ``'frequency'`` 
           
        threshold_line: float, optional
           Visualize th threshold line. Can be ['bahr', 'swift', 'both']:
               
           - Note that when method is set to ``swift``, the value close 
             to close to :math:`0.` assume the 1D and 2D structures 
             (:math:`\eta <0.4`),  and 3D otherwise( :math:`\eta >0.4`). 
             The threshold line for ``swift`` is set to :math:`0.4`. 
             
           - when method is set to ``Bahr``, :math:`\eta > 0.3``  is 3D 
             structures, between :math:`[0.1 - 0.3]` assumes modified 3D/2D 
             structures whereas :math:`<0.1` 1D, 2D or distorted 2D. 
        show_average_sensistivity: bool, default=True 
           Display the averaged value of skew data at all -frequencies. 
           Value can help a dimensionality interpretation purposes. 
           
        suppress_outliers: bool, default=True
           Remove the outliers in the data if exists. It uses the 
           Inter Quartile Range (``IQR``) approach. See the documentation 
           of :func:`watex.utils.remove_outliers`. This is useful for clear 
           interpretation using the skew threshold value. 
          
        See Also
        ---------
        watex.methods.EMAP.skew: 
            For mathematical skew `Bahr` and `Swift` concept formulations. 
        watex.utils.plot_skew: 
            For phase sensistive skew visualization - naive plot.
  
        Examples
        --------
        >>> import watex 
        >>> test_data = watex.fetch_data ('edis', samples =37, return_data =True )
        >>> watex.TPlot(fig_size =(10,  4), marker ='x').fit(
            test_data).plotSkew(method ='swift', threshold_line=True)
        
        References 
        -----------
        .. [1]  Swift, C., 1967. A magnetotelluric investigation of an 
                electrical conductivity  anomaly in the southwestern United 
                States. Ph.D. Thesis, MIT Press. Cambridge.
        .. [2] Bahr, K., 1991. Geological noise in magnetotelluric data: a 
               classification of distortion types. Physics of the Earth and 
               Planetary Interiors 66 (1–2), 24–38.
        """
        self.inspect 
        
        sensitivity = str(sensitivity).lower() 
        for ix in ('inv', 'rot', 'mu'): 
            if sensitivity.find(ix)>=0: 
                sensitivity ='mu' 
                break 
            
        sensitivity='skew' if sensitivity=='none' else sensitivity 
        assert sensitivity in {"skew", 'mu'}, ("expect 'skew' or 'rotational'"
                                        f" invariant plot, got {sensitivity!r}")
        
        if 'period' in str(mode).lower(): 
            mode ='period'

        skew, mu =self.p_.skew(
            method = method, suppress_outliers = suppress_outliers
            )
        freqs =  1/ self.p_.freqs_ if mode =='period' else self.p_.freqs_ 
        ymat = skew if sensitivity =='skew' else mu 
        
        fig, ax = plt.subplots(figsize = self.fig_size )

        #---manage threshold hline ------
        thr_code = {"bahr": [1] , "swift":[ 2] , 'both':[1, 2] }
        
        if str(threshold_line).lower()=='true': 
            threshold_line = str(method).lower() 
            
        if threshold_line is not None: 
            if str(threshold_line).lower() in ("*", "both" ): 
                threshold_line = 'both'
                
        ct = thr_code.get(str(threshold_line).lower(), None ) 
        
        for i in range (skew.shape[1]): 
            ax.scatter ( freqs, reshape (ymat[:, i]),
                        marker = plot_kws.get ('marker', None) or self.marker, 
                        cmap = plot_kws.get('cmap', None) or self.cmap, 
                        alpha = plot_kws.get('alpha', None) or self.alpha, 
                        s = plot_kws.get('s', None) or self.s , 
                        **plot_kws 
                        )
        if ct: 
            for m in ct: 
                plt.axhline(y=0.4 if m==2 else 0.3 , color="k" if m==1 else "g",
                            linestyle="-",
                            label=f'threshold: $\mu={0.4 if m==2 else 0.3}$'
                            )
                # plt.legend()

        # see phase sensitive trend 
        if show_average_sensistivity: 
            plt.text(x= np.nanmin(freqs) , y= np.nanmax(ymat), 
                     s="aver.-{}:{}={}".format(sensitivity, str(method).capitalize(), 
                    np.around (np.average(ymat[ ~np.isnan(ymat)]), 3)),  
                    fontdict= dict (style ='italic',  bbox =dict(
                         boxstyle='round',facecolor ='#CED9EF'))
                     ) 
        
        ax.set_xscale('log')
        ax.set_xlabel('Period ($s$)' if mode=='period' 
                      else 'Frequency ($H_z$)' or self.xlabel )
        ax.set_ylabel(f"{'Skew' if sensitivity =='skew' else 'Rot.Invariant'}" + "($\mu$)"
                      or self.ylabel )

        plt.xlim ([ freqs.min() , freqs.max()] or self.xlim )
        
        plt.xlim() 
        
        if ct: ax.legend() 
        
        if self.savefig is not  None: 
            plt.savefig (self.savefig, dpi = self.fig_dpi)
            
        plt.close () if self.savefig is not None else plt.show() 
        
        return self 
    
    def _check_component_validity (self, tensor, components ): 
        """Retrieve resistiviy, phase or impedance tensors from 
        EDI objets if component exists. 
        
        Parameters 
        -----------
        tensor: str, 
          Name of tensor. Could be ['resistivity'| 'phase'|'z']
        components: str, list, 
          Name of components. Could be ['xx', 'xy', 'yx', 'yy']
        
        Returns
        --------
        rp: list of valid 2D dimensional tensors and ``None`` if 
          no valid tensors are found. 
        
        """
        rp =[] 
        tensor =str(tensor) 
        components = is_iterable(components, exclude_string =True,
                                transform =True, parse_string =True )
        for c in components : 
            try: 
                mat2d = self.p_.make2d (f'{tensor+c}')
            except :continue 
            else: rp.append(mat2d )
            
        return rp if len(rp)!=0 else None 
    
    def plot_rhoa(
        self, 
        mode ='TE', 
        scale ='period', 
        site =None, 
        seed = None, 
        how ='py', 
        show_site=True,
        survey= None, 
        style=None, 
        errorbar=True, 
        suppress_outliers=False, 
        **kws
        ): 
        """ Plot apparent resistivity and phase curves 
        
        Parameters 
        ----------
        mode: str, default='TE', 
          Electromagnetic mode. Can be ['TM' |'both']. If ``both``, 
          components `xy` and `yx` are expected in the data. 
          
        scale: str, default='period'
          Visualization on axis labell. can be ``'frequency'``. 
          
        site: int,str, optional 
          index of name of the site to plot. `site` must be composed of 
          a position number. For instance ``'S13'``. If not provided, 
          a random station is selected instead. 
          
        seed : int, optional 
           If site is not provided, seed fetches randomly a site. To fetch 
           the same sime everytimes, it is better to set the seed value. 
           
        how: str, default='py'
          The way the site is fetched for plot. For instance, in Python 
          indexing (default), the site is numbered from 0. For instance 
          'site05' will fetch the data at index 4. If this positioning 
          is not wished, set to 'None'.
          
        show_site:bool, default=True, 
          Display the number of site. 
          
        survey: str, optional 
          Method used for the survey. e.g., 'AMT' for |AMT|. 
         
        style:str, default='default'
          Matplotlib style. 
          
        errorbar: bool, default=True 
          display the error bar.  
          
        suppress_outliers: bool, default=False, 
          Remove outliers in the data before plotting 
          
        kws: dict, 
          Addfitional keywords arguments passed to 
          Matplotlib.Axes.Scatter plots. 
         
        Examples
        ---------
        >>> import watex as wx 
        >>> edi_data = wx.fetch_data ('edis', return_data =True, samples =27)
        >>> wx.methods.TPlot(show_grid=True).fit(edi_data).plot_rhoa (
            seed =52, mode ='*')
        
        """
        self.inspect 
        
        m=_validate_name_in(mode,  ('*', 'both', 'tetm'), expect_name='*')
        
        if m!='*':
            m= _validate_name_in(mode, defaults = 'tm transverse-magnetic',
                                     expect_name ='tm' )
        if not m: 
            m='te' 

        scale = _validate_name_in(scale, deep =True, defaults='periods', 
                                 expect_name='period')

        cpm = {'te': ["xy"] , 'tm': ["yx"], '*': ('xy', 'yx') }
        
        components = cpm.get(m)
        
        res, phs, site, *s= self._validate_correction (
                             components = components, 
                             errorbar = errorbar , 
                             how = how, 
                             seed = seed , 
                             sites = site, 
                             style =style , 
                             n_sites = 1.
                             )  
        s,  res_err, phs_err  = s 
        # plot only single data 
        site = site [0] ; s = s[0]
        # get the single site 

        fig = plt.figure(self.fig_num , figsize= self.fig_size,
                         dpi = self.fig_dpi , # layout='constrained'
                         )

        gs = GridSpec(3, 1, figure = fig ) 
        
        ax1 = fig.add_subplot (gs[:-1, 0 ])
        ax2 = fig.add_subplot(gs [-1, 0 ], sharex = ax1 )
        plt.setp(ax1.get_xticklabels(), visible=False)
        
        survey= survey or self.p_.survey_name 
        if not survey: survey=''
        
        colors = [ '#069AF3', '#DC143C']
    
        #==plotlog10 --------
        res= [ np.log10 (r) for r in res] 
        # the complete frequency 
        fp = self.p_.freqs_
        
        fp =  1/ fp if scale =='period' else fp 
        
        fp =  np.log10 ( fp) 
        
        if suppress_outliers: 
            res = remove_outliers(res, fill_value=np.nan) 
            phs = remove_outliers(phs, fill_value=np.nan) 
            if errorbar: 
                res_err = remove_outliers(
                    res_err, fill_value=np.nan) 
                phs_err = remove_outliers(
                    phs_err, fill_value=np.nan) 
                
        min_y =  np.nanmin(res[0][:, site])
        
        # add error bar data to main 
        data = [res, phs ] 
        data +=  [ res_err ,  phs_err ] if errorbar else []
        
        for i, sloop in enumerate (zip (* data )) : 
            r, p, *sl = sloop 
            
            if len(sl) !=0 : 
                e, ep = sl  # mean errorbar is set to True 
            
            y =  reshape (r[:, site])
            if errorbar: 
                plot_errorbar (ax1 , 
                               fp, 
                               y,  
                               y_err = reshape (e[:, site]),
                               )
            ax1.scatter (fp  , y, 
                          marker =self.marker, 
                          color =colors [i],
                          edgecolors='k', 
                          label = fr'{survey}$\rho_a${components[i]}',
                          **kws 
                          ) 
            if errorbar: 
                plot_errorbar (ax2 , 
                               fp, 
                               reshape (p[:, site]),
                               y_err = reshape (ep[:, site]),
                               )
            ax2.scatter( fp, 
                        reshape (p[:, site]),
                        marker =self.marker, 
                        color =colors [i] ,
                        edgecolors='k', 
                        label = f'{survey}$\phi${components[i]}',
                        **kws
                        ) 
            min_y = np.nanmin (y) if np.nanmin (
                y) < min_y else min_y 
            try: 
                ax1.legend(ncols = len(res)) 
                ax2.legend(ncols = len(phs)) 
                
            except: 
                # For consistency in the case matplotlib  is < 3.3. 
                ax1.legend() 
                ax2.legend() 
                
        if show_site:
            ax1.text (np.nanmin(fp),
                      min_y,
                      f'site {s}', 
                      fontdict= dict (style ='italic',  bbox =dict(
                           boxstyle='round',facecolor ='#CED9EF'), 
                          alpha = 0.5 )
                      )
        
        ax2.set_ylim ([0, 90 ])
        xlabel = self.xlabel or ( 'Log$_{10}$Period($s$)' if scale=='period' 
                                 else 'Frequency ($H_z$)') 
        
        ax2.set_xlabel(xlabel ) 
        ax1.set_ylabel(self.ylabel or r'Log$_{10}\rho_a$($\Omega$.m)') 
 
        ax2.set_ylabel('$\phi$($\degree$)')
        
        if self.show_grid :
            for ax in (ax1, ax2 ): 
                ax.grid (visible =True , alpha =self.galpha,
                         which =self.gwhich, color =self.gc)
          
            
        if self.savefig is not  None: 
            plt.savefig (self.savefig, dpi = self.fig_dpi)
            
        plt.close () if self.savefig is not None else plt.show() 
        
        return self 
    
    def plot_rhophi( 
        self, 
        sites =None,  
        mode ='TE', 
        scale ='period', 
        seed = None, 
        how ='py', 
        show_site=True,
        survey= None, 
        style=None, 
        errorbar=True, 
        suppress_outliers=False, 
        kind='2', 
        n_sites= 1, 
        spad=.5, 
        **kws
        ): 
        """ Plot resistivities and phases from multiples stations. 
        
        Parameters 
        ----------
        mode: str, default='TE', 
          Electromagnetic mode. Can be ['TM' |'both']. If ``both``, 
          components `xy` and `yx` are expected in the data. 

        sites: int,str, or list,  optional 
          A collection of index of name of the site . Each `site` must be 
          composed of a position number. For instance ``'S13'``. If not
          provided, a random sites are selected instead using the `n_sites` 
          parameter. 
         
        scale: str, default='period'
          Visualization on axis labell. can be ``'frequency'``. 
        seed : int, optional 
           If site is not provided, seed fetches randomly a site. To fetch 
           the same sime everytimes, it is better to set the seed value. 
           
        how: str, default='py'
          The way the site is fetched for plot. For instance, in Python 
          indexing (default), the site is numbered from 0. For instance 
          'site05' will fetch the data at index 4. If this positioning 
          is not wished, set to 'None'.
          
        show_site:bool, default=True, 
          Display the number of site. 
          
        survey: str, optional 
          Method used for the survey. e.g., 'AMT' for |AMT|. 
         
        style:str, default='default'
          Matplotlib style. 
          
        errorbar: bool, default=True 
          display the error bar.  
          
        suppress_outliers: bool, default=False, 
          Remove outliers in the data before plotting 
           
        n_sites: int, default =1. 
           Number of random sites to select for visualizing. It cannot work 
           if the names of sites are given. 
           
        spad: float, default=.5, 
          pad to display the station in the top of each section plot. 
          
          .. versionadded:: 0.2.1 
             
        kws: dict, 
          Addfitional keywords arguments passed to 
          Matplotlib.Axes.Scatter plots. 
          
        Examples
        ---------
        >>> import watex as wx 
        >>> edi_data = wx.fetch_data ('edis', return_data =True, samples =27)
        >>> wx.methods.TPlot(show_grid=True).fit(edi_data).plot_rhophi (
            seed =52, mode ='*', n_sites =3 )
        """
        
        self.inspect 

        m=_validate_name_in(mode,  ('*', 'both', 'tetm'), 
                            expect_name='*')
        
        if m!='*':
            m= _validate_name_in(mode, defaults = 'tm transverse-magnetic',
                                     expect_name ='tm' )
        if not m: 
            m='te' 

        scale = _validate_name_in(scale, deep =True, defaults='periods', 
                                 expect_name='period')

        cpm = {'te': ["xy"] , 'tm': ["yx"], '*': ('xy', 'yx') }
        
        components = cpm.get(m)
        res, phs, sites, *s= self._validate_correction (
                             components = components, 
                             errorbar = errorbar , 
                             how = how, 
                             seed = seed , 
                             sites = sites , 
                             style =style , 
                             n_sites = n_sites, 
                             )  
        s,  res_err, phs_err  = s 

        survey= survey or self.p_.survey_name 
        if not survey: survey=''
        
        #colors = [ '#069AF3', '#DC143C']
        colors = [ '#0000FF', '#FF00FF']

        #==plotlog10 --------
        #xxxxxxxxxxxxxxxxxxxx
        # res= [ np.log10 (r) for r in res] 
        # the complete frequency 
        fp = self.p_.freqs_
        
        fp =  1/ fp if scale =='period' else fp 
        if suppress_outliers: 
            res = remove_outliers(res, fill_value=np.nan) 
            phs = remove_outliers(phs, fill_value=np.nan) 
            if errorbar: 
                res_err = remove_outliers(
                    res_err, fill_value=np.nan) 
                phs_err = remove_outliers(
                    phs_err, fill_value=np.nan) 
        
        # make sites coordinates to place sites 
        # assert whether the number of sites fit the row values 
        sy =[]
        for ii in sites: 
            exp_sites = (len(res[0][0, :]) -1) if how=='py' else len(res[0][0, :])
            if ii > exp_sites: 
                raise SiteError (
                    f"Expects {exp_sites} sites. Got {ii}. Note that" 
                    f" for how={how!r}, the site numbering starts" 
                    f" at {0 if how=='py' else 1}."
                    )
            sy.append ( (np.nanmax(res[0][:, ii]) - np.nanmin(res[0][:, ii])) /2) 
            
        sy = np.average ( sy )
        # sy=  np.average ( [
        #     ( np.nanmax(res[0][:, ii]) - np.nanmin(res[0][:, ii])) /2  
        #     for ii in sites ] )
        sy += spad 
        sx = np.average (fp)
                                
 
        # add error bar data to main 
        data = [res, phs ] 
        data +=  [ res_err ,  phs_err ] if errorbar else []
        
        # make thoa and phase labels 
        rlabels = [fr'{survey}$\rho_a${components[i]}' 
                   for i in range (len(res))]
        plabels = [f'{survey}$\phi${components[i]}' 
                   for i in range(len(phs))]
        
        self._plot_grid_spec (
                data = data , 
                x= fp,  
                sites =sites, 
                errorbar =errorbar, 
                colors = colors, 
                xysites= ( sx, sy ),
                show_site =show_site, 
                scale =scale, 
                rlabels = rlabels, 
                plabels = plabels,
                kind= kind, 
                **kws
                )
        

        if self.savefig is not  None: 
            plt.savefig (self.savefig, dpi = self.fig_dpi)
            
        plt.close () if self.savefig is not None else plt.show() 
         
        return self 
    
    def _plot_grid_spec (
        self, 
        data, 
        x,  
        sites =None, 
        errorbar =False  , 
        colors = None, 
        show_site =False, 
        scale =None, 
        xysites = None, 
        color_mode='color', 
        kind='2', 
        **kws
        ): 
        """ Plot multiple stations using the SpecGrid  
        
        Parameters 
        -----------
        data: list, 
           A collection of resistivity, errors and phases 
           
        x: arraylike 
          Arraylike one-dimensional for plotting data. It should be the 
          frequency array or periods  

        sites: int,str, optional 
          index of name of the site to plot. `site` must be composed of 
          a position number. For instance ``'S13'``. If not provided, 
          a random station is selected instead. 

        errorbar: bool, default=True 
          display the error bar.
          
        colors: str, list 
          a collection of matplotlib colors 
          
        show_site: bool, default=False, 
          Display the name of the site in each section 

        style:str, default='classic'
          Matplotlib style. 
          
        scale: str, default='period'
          Visualization on axis labell. can be ``'frequency'``.
          
        mode: str, {'1', '2'} , default='2'
          choice of plot style. ``mode='2'`` plots only the errorbar and '1' 
          add scatter plots. 
          
        color_mode: str, {"color", "bw"}, default='color' 
          Plot tensor in different colors by default otherwise plots in 
          black-white. This parameter is triggered only if `mode` is set ``2``.
          
        xysites: tuple , optional 
          The coordinates to locate the text of each station. 
 
        kws: dict, 
           Additional keywords passed to matplotlib.scatter plot. 
           Also to rename the labels of resistivy and phase, pass a list 
           of rho and phase labels in parameters `rlabels` and `plabels` 
           respectively. 
           
        Returns 
        --------
        axr, axp : list of Matplotlib.Axes 
          A collection of Matplotlib axes of each stations 
          
        """

        ncols = len (sites) if sites is not None else  1 
        
        fig = plt.figure(figsize = self.fig_size, dpi=self.fig_dpi)
        h_ratio = [1.5, 1, .5]
        
        gs = GridSpec(2, ncols or 1,
            wspace=0. if kind =='2' else .3, # .3,if 
            left=.08,
            top=.85,
            bottom=.1,
            right=.98,
            hspace=.0,
            height_ratios=h_ratio[:2])
            
        sharey = None
        # make a list of axes 
        # to return axes 
        # for another plots 
        axr, axp =[], []
        #++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        #if kind =='2': 
        # color mode
        x /= 1 # inverse , take a periods 
        if str( color_mode) .lower() == 'color':
            # color for data
            cted =  (0, 0, 1)
            ctmd = (1, 0, 0)
            mted = 's'
            mtmd =  'o'
        # black and white mode
        elif color_mode == 'bw':
            # color for data
            cted = (0, 0, 0)
            ctmd = (0, 0, 0)
            mted = 's'
            mtmd = 'o'
            
        # --> make key word dictionaries for plotting
        ms =  1.5
        # ms_r =  3
        lw = .5
        # lw_r =  1.0
        # ls = ':'
        e_capthick =  .5
        e_capsize =  2
        
        # kw_xx=dict(); kw_yy=dict()

        res_limits =[]; phase_limits=[]
        sharey2 =None
        #np.savetxt ( 'x.txt', x )
        #++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        for j, site in enumerate ( sites ): 
            ax1 = fig.add_subplot (gs [ 0, j] , 
                                   sharey = sharey) 
            if j==0: sharey = ax1 
            
            if errorbar: 
                
                ax2 = fig.add_subplot (gs [1,  j], sharey =sharey2 ) 
                if j==0 and kind =='2': sharey2 = ax2 
                
                
            for i, sloop in enumerate (zip (* data )) : 
                r, p, *sl = sloop 
    
                if len(sl) !=0 : 
                    e, ep = sl  # mean errorbar is set to True 
                
                y =  reshape (r[:, site])
                colors = [cted,ctmd ]
                markers = [mted, mtmd]
                kw_xx = {'color': colors[i],
                         'marker': markers[i],
                         'ms': ms,
                         'ls': ':',
                         'lw': lw,
                         'e_capsize': e_capsize,
                         'e_capthick': e_capthick}

                kw_yy = {'color': colors[i],
                         'marker': markers[i],
                         'ms': ms,
                         'ls': ':',
                         'lw': lw,
                         'e_capsize': e_capsize,
                         'e_capthick': e_capthick}
                    
                #if errorbar: 
                plot_errorbar (ax1 , 
                               x, 
                               y, #if i ==0 else y ,  
                               y_err = reshape (e[:, site]),
                               **kw_xx
                               )
                plot_errorbar (ax2 , 
                               x, 
                               reshape (p[:, site]),
                               y_err = reshape (ep[:, site]),
                               **kw_yy,
                               )
                res_limits.append ((min(y), max(y)))
                phase_limits.append( (min(reshape (p[:, site])), 
                                      max(reshape (p[:, site]))
                                      )
                                    )
            if show_site:
                ax1.set_title( f'site {site}', 
                              fontdict={'size': 8 + 2,
                                       'weight': 'bold'})
            axr.append( ax1);  axp.append (ax2)
             
        # --> set default font size
        self.font_size =  6
        plt.rcParams['font.size'] = self.font_size

        fontdict = {'size': self.font_size + 2, 
                    'weight': 'bold'}
        
        for ax0,  site in zip(axr, sites):
            ax0.set_title(f'S{site}', fontdict={'size': self.font_size + 2,
                                                  'weight': 'bold'})
        #     # set axis properties
        # set ylimit 
        res_limit_max=np.array( list(
            map ( lambda x: x[1], res_limits )) )
        res_limit_min=np.array( list( 
            map ( lambda x: x[0], res_limits )))
        res_limits_d= [10 **np.floor (np.log10(res_limit_min.min())),
                      10 **np.ceil (np.log10(res_limit_max.max())) ]
        # phase limit 
        phase_limit_max=np.array( list(
            map ( lambda x: x[1], phase_limits )) )
        phase_limit_min=np.array( list( 
            map ( lambda x: x[0], phase_limits )))
        phase_limits_d= [np.floor (phase_limit_min.min()),
                      np.ceil (phase_limit_max.max()) ]
        phase_limits_d=None
  
        ax_list = [*axr, *axp ]
        for aa, ax in enumerate(ax_list):
            ax.tick_params(axis='y', pad=1.5)

            ax.set_xlabel('Period (s)', 
                          fontdict=fontdict
                          )
            if aa < len(ax_list)//2 : #4 :
                ylabels = ax.get_yticklabels()
                ylabels[0] = ''
                ax.set_yticklabels(ylabels)
                ax.set_yscale('log', #nonposy='clip'
                              )
                try: 
                    ax.set_ylim(res_limits_d)
                except: 
                    ax.set_ylim(None)
                    res_limits_d=None 

            if aa >= len(ax_list)//2 :
                ax.yaxis.set_major_locator(mticker.MultipleLocator(10.0))
                if phase_limits_d is not None:
                    ax.set_ylim(phase_limits_d)
                    
            # set axes labels
            if aa == 0:
                ax.set_ylabel('App. Res. ($\mathbf{\Omega \cdot m}$)',
                                  fontdict=fontdict)
            elif aa == 0 or aa == len(ax_list)//2:
                ax.set_ylabel('Phase (deg)',
                                  fontdict=fontdict)
            ax.set_xscale('log', # nonposx='clip'
                          )
            # set period limits
            period_limits = (10 ** (np.floor(np.log10(x[0]))) * 1.01,
                                  10 ** (np.ceil(np.log10(x[-1]))) * .99)
            ax.set_xlim(xmin=period_limits[0],
                        xmax=period_limits[1])
            ax.grid(True, alpha=.25)

            
            if kind=='2':
                if aa !=0 or aa != len(ax_list)//2:
                    ax.set_yticklabels('')
            else: 
                ylabels = ax.get_yticks().tolist()
                ylabels[-1] = ''
                ylabels[0] = ''
                ax.set_yticklabels(ylabels)
            
            if aa < len(ax_list)//2: 
                plt.setp(ax.get_xticklabels(), visible=False)

        return axr, axp 
           
    def _axesproperties1 (self, j, ax1, ax2, r, p, sites , scale  ): 
        """ Set properties of plot kind number 1. """
        if j > 0: 
            plt.setp(ax1.get_yticklabels(), visible=False)
            plt.setp(ax2.get_yticklabels(), visible=False)
        # Put the legend in the last image
        if j == len(sites)-1: 
            try: 
                ax1.legend(ncols = len(r)) 
                ax2.legend(ncols = len(p)) 
            except: 
                # For consistency in the case matplotlib  is < 3.3. 
                ax1.legend() 
                ax2.legend() 
             
        ax1.set_xscale ('log') ;  ax1.set_yscale ('log') 
        ax2.set_xscale ('log')
        
        ax2.set_ylim ([0, 90 ])
        
        xlabel = self.xlabel or ( 'Period($s$)' if scale=='period' 
                                 else 'Frequency ($H_z$)') 
        
        ax2.set_xlabel(xlabel) 
        
        if j ==0 : 
            # avoid reapeting this 
            ax1.set_ylabel(self.ylabel or r'$\rho_a$($\Omega$.m)') 
            ax2.set_ylabel('$\phi$($\degree$)')
        
        if self.show_grid :
            for ax in (ax1, ax2 ): 
                ax.grid (visible =True , alpha =self.galpha,
                         which =self.gwhich, color =self.gc)
                
        return ax1 , ax2 
    
    
    def _validate_correction (
        self, 
        components , 
        errorbar , 
        seed , 
        sites , 
        how , 
        style , 
        n_sites, 
        ): 
        """Isolated part to validate the :meth:`plot_corrections` and 
        :meth:`plot_rhoa` arguments. 
        
        Parameters
        ----------
        
        components: str ,
           could be 'xx', 'xy', 'yx' or 'yy' 

        sites: int,str, optional 
          index of name of the site to plot. `site` must be composed of 
          a position number. For instance ``'S13'``. If not provided, 
          a random station is selected instead. 
          
        seed : int, optional 
           Get the same site if site is not provided. `seed` fetches 
           a random number of site. 
           
        how: str, default='py'
          The way the site is fetched for plot. For instance, in Python 
          indexing (default), the site is numbered from 0. For instance 
          'site05' will fetch the data at index 4. If this positioning 
          is not wished, set to 'None'.
        
        style:str, default='default'
          Matplotlib style. 
          
        errorbar: bool, default=True 
          display the error bar.
          
        n_sites: int, 
          Number of sites to randomly diplay when sites is not given. 
          
          
        Returns 
        --------
        ( fp, res, phs, site, s ,  res_err , phs_err) : Tuple 
        
          - fp: frequency array 
          - res:  resistivity tensor collected at a specific components 
          - phs: phase tensor collected at a specific component 
          - site: The site number 
          - s : position of the site 
          - res_err: error in resistivity at a specific component 
          - phs_err: error in phase at a specific components. 
          
        """ 
        
        res = self._check_component_validity('res', components)
        phs = self._check_component_validity('phase', components)
        
        res_err , phs_err =[], []
        if errorbar: 
            res_err = self._check_component_validity(
                'res_err', components)
            phs_err = self._check_component_validity(
                'phase_err', components)
  
        
        terror =("{0!r} does not contain component {}. Provide the"
                 " right component of the valid tensor.")
        if res is None: 
            raise EMError(terror.format('resistivity', components))
        if phs is None: 
            raise EMError(terror.format('phase', components))
            
        # assert sites 
        sites, s = self._validate_sites(res, sites = sites, seed = seed , 
                                        n_sites = n_sites , how = how )
 
        try: 
            plt.style.use ( style or 'default')
        except : 
            warnings.warn(
                f"{style} is not available. Use `plt.style.available`"
                " to get the list of available styles.")
            plt.style.use ('default')

       
        return res, phs, sites, s ,  res_err, phs_err 
 
    def _validate_sites (self, 
            data, /, sites = None, seed = None, n_sites = 1, how ='py' 
                         ): 
        """ validate sites or choose random sites from number of stations 
        in the survey data. 
        
        Parameters 
        -----------
        data: List of resistivity-error and phases 
           A collection of resistivy , errors  and phases from 
           EDI-objects 
        sites: str, list 
          A collection of sites to visualize. 
          
        seed : int,  
          `seed` is used to reproduce the same stations when sites are not 
           given.
         
        n_sites: int, default=1 
          Number of sites to randomly selected for displaying. Note that it 
          only works if sites are ``None``. 
        how: str, default='py' 
          The way to fetch and display sites. By default used the Python 
          Indexing i.e the site starts with 0 
          
        Returns 
        -------
        S, s: Tuple 
          Tuple of collection of sites and sites indexes. 
          
        """
        
        # assert sites 
        if seed: 
            seed = _assert_all_types(seed , int, float, objname ='Seed')
            np.random.seed (seed ) 
           
        if sites is None:
            n_sites = int(n_sites ) if n_sites else n_sites 
            sites = np.random.permutation (range (data[0].shape[1])
                                           )[:int (n_sites)]
            # sites = [ np.random.choice (
            #     range (res[0].shape[1])) for i in range (nsites)] 
            
        # make site as an iterable object 
        sites = is_iterable(sites, exclude_string= True , transform =True )
        
        s= copy.deepcopy(sites)
        sites = [ re.search ('\d+', str(site), flags=re.IGNORECASE).group()  
                 for site in sites ] 
        S = []
        for ii, site in enumerate ( sites)  : 
            try: 
               site= int(site)
            except TypeError: 
                raise TypeError ("Missing position number. Station must prefix"
                                 f" with position, e.g. 'S7', got {s[ii]!r}")
            
            site = abs (site) + 1 if how !='py' else site 
            
            if site > data[0].shape [1] : 
                raise ValueError (
                    f"Site position {site} is out of the range. The total"
                    f" number of sites/stations ={data[0].shape [1]}")
            S.append (site)
            
        return  S, s  
    
    def plot_corrections(
        self, 
        fltr='ama',
        ss_fx =None, 
        ss_fy=None, 
        r=1000., 
        nfreq=21,
        skipfreq=5, 
        tol=.12,
        rotate=0., 
        distortion=None, 
        distortion_err =None, 
        mode ='TE', 
        scale ='period', 
        sites =None, 
        seed = None, 
        how ='py', 
        show_site=True,
        survey= None, 
        style=None, 
        errorbar=True,
        spad =.5, 
        n_sites = 1, 
        mcolors= None,
        markers = None, 
        **kws
        ): 
        """Plot apparent resistivity/phase curves and corrections.  
        
        .. versionchanged:: 0.2.1 
           Can henceforth display multiple sites by providing the 
           sites as a collection. 
           
        Parameters 
        ----------
        fltr: str , default='ama'
           Type of filter to apply. ``ss`` is used to remove the static 
           shift using spatial median filter. Whereas ``dist`` is for 
           distorsion removal. Note that `distortion` might be provided 
           otherwise an error raises. Can also be ['tma'|'ama'|'flma'] for 
           :term:`EMAP` filters. 
           
           - ``tma`` for trimming moving-average 
           - ``ama`` for adaptative moving-average 
           - ``flma`` for fixed-length moving-average 
           
           .. versionadded: 0.2.1 
               Applied EMAP filters for the visualization.
  
        distortion_tensor: np.ndarray(2, 2, dtype=real) 
           Real distortion tensor as a 2x2
   
        error: np.ndarray(2, 2, dtype=real), Optional 
          Propagation of errors/uncertainties included
          
        ss_fx: float, Optional  
           static shift factor to be applied to x components
           (ie z[:, 0, :]).  This is assumed to be in resistivity scale. 
           If None should be automatically computed using  the 
           spatial median filter. 
           
        ss_fy: float, optional 
           static shift factor to be applied to y components 
           (ie z[:, 1, :]).  This is assumed to be in resistivity scale. If
           ``None`` , should be computed using the spatial filter median.
           
        r: float, default=1000. 
           radius to look for nearby stations, in meters.
 
        nfreq: int, default=21 
           number of frequencies calculate the median static shift.  
           This is assuming the first frequency is the highest frequency.  
           Cause usually highest frequencies are sampling a 1D earth.  
    
        skipfreq: int, default=5 
           number of frequencies to skip from the highest frequency.  
           Sometimes the highest frequencies are not reliable due to noise 
           or low signal in the :term:`AMT` deadband.  This allows you to 
           skip those frequencies.
     
        tol: float, default=0.12
           Tolerance on the median static shift correction.  If the data is 
           noisy the correction factor can be biased away from 1.  Therefore 
           the shift_tol is used to stop that bias.  If 
           ``1-tol < correction < 1+tol`` then the correction factor is set 
           to ``1``
           
        rotate: float, default=0.  
            Rotate Z array by angle alpha in degrees.  All angles are referenced
            to geographic North, positive in clockwise direction.
            (Mathematically negative!).
            In non-rotated state, X refs to North and Y to East direction. 
            
        mode: str, default='TE', 
          Electromagnetic mode. Can be ['TM' |'both']. If ``both``, 
          components `xy` and `yx` are expected in the data. 
          
        scale: str, default='period'
          Visualization on axis labell. can be ``'frequency'``.
          
        sites: int,str, optional 
          index of name of the site to plot. `site` must be composed of 
          a position number. For instance ``'S13'``. If not provided, 
          a random station is selected instead. 
          
        seed : int, optional 
           Get the same site if site is not provided. `seed` fetches 
           a random number of site. T
           
        how: str, default='py'
          The way the site is fetched for plot. For instance, in Python 
          indexing (default), the site is numbered from 0. For instance 
          'site05' will fetch the data at index 4. If this positioning 
          is not wished, set to 'None'.
          
        show_site:bool, default=True, 
          Display the number of site. 
          
        survey: str, optional 
          Method used for the survey. e.g., 'AMT' for |AMT|. 
         
        style:str, default='default'
          Matplotlib style. 
          
        errorbar: bool, default=True 
          display the error bar.  
          
        spad: float, default=.5, 
          pad to display the station in the top of each section plot. 
          
          .. versionadded:: 0.2.1 
          
        n_sites: int, default =1. 
           Number of random sites to select for visualizing. It cannot work 
           if the names of sites are given. 

        mcolors: str, list, optional 
           The list of colors for resistivy and phase. 
           
           
        markers : str, list, optional 
           The list of marker for resistivy and phase. 
        markers = None, 
        
        kws: dict, 
          Addfitional keywords arguments passed to 
          Matplotlib.Axes.Scatter plots. 
         
        Examples
        ---------
        >>> import numpy as np 
        >>> import watex as wx 
        >>> edi_data = wx.fetch_data ('edis', return_data =True, samples =27)
        >>> wx.methods.TPlot(show_grid=True).fit(edi_data).plot_corrections (
            seed =52, )
        >>> distortion = np.array([[1.1 , 0.6 ],[0.23, 1.9 ]])
        >>> wx.methods.TPlot(show_grid=True).fit(edi_data).plot_corrections (
             seed =52, mode ='tm', fltr ='dist', distortion =distortion 
             )
        """
        self.inspect 
    
        m=_validate_name_in(mode,  'tm transverse-magnetic', expect_name='tm')
        if not m: 
            m='te' 
        scale = _validate_name_in(scale, deep =True, defaults='periods', 
                                 expect_name='period')

        cpm = {'te': ["xy"] , 'tm': ["yx"]}
       
        components = cpm.get(m)
        res, phs, sites, *s= self._validate_correction (
                             components = components, 
                             errorbar = errorbar , 
                             how = how, 
                             seed = seed , 
                             sites = sites , 
                             style =style ,
                             n_sites = n_sites, 
                             )  
        s,  res_err, phs_err  = s 
        # plot only single correction so 
        # Assert filters 
        mc = _validate_name_in(fltr, defaults =('static shift', 'ss', '1'), 
                               expect_name= 'ss')
        if mc!='ss': 
            mc = _validate_name_in(fltr, defaults=('distortion', 'dist', '2'), 
                                   expect_name ='dist')
        if mc not in ('dist', 'ss') : 
            if str(fltr).lower() not in ( 'tma', 'ama', 'flma'): 
                ff = ('ss', 'dist', 'tma', 'ama', 'flma')
                raise ValueError(f"Wrong filter {fltr!r}. Expect"
                                 f"{smart_format(ff, 'or')} for corrections."
                                 )
            else: mc = str (fltr).lower() 
           
        if mc=='dist' and distortion is None: 
            raise TypeError("Distorsion cannot be None!")
        
        # -> compute the corrected values 
        zo = MT().fit(self.p_.ediObjs_)
        if mc =='ss': 
            zo.remove_static_shift (
                ss_fx = ss_fx , 
                ss_fy = ss_fx, 
                nfreq = nfreq ,         
                r=r, 
                skipfreq=skipfreq , 
                tol=tol, 
                rotate = rotate, 
                )
            
        elif mc =='dist': 
            zo.remove_distortion (
                distortion , 
                error = distortion_err 
                )
        
        else: 
            zo.remove_ss_emap (fltr =mc )
            
        # set zcorrected 
        zc = zo.new_Z_ 
        
        zc_res = [ z.resistivity[tuple (self._c_.get(components[0])) ] 
                  for z in zc  ] 
        # zc_res = [ np.log10(r) for r in zc_res ] # convert to log10 res 
        # --> phase 
        zc_phase = [ z.phase[tuple (self._c_.get(components[0])) ] 
                  for z in zc ] 
        # mofulo the phase to be 0 and 90 degree 
        zc_phase = [ np.abs (p)%90  for p in zc_phase ] 
        
        # ----------------------end ---------------------------------

        survey= survey or self.p_.survey_name 
        if not survey: survey=''
        
        # set defaults colors  and markers 
        
        #colors = [ '#069AF3', '#DC143C']
        colors = [] if mcolors is None else mcolors
        c = is_iterable(colors , exclude_string =True , transform =True )
        colors = list(c ) + [(0, .6, .3), (.9, 0, .8) ] 
        
        markers = [] if markers is None else markers 
        m = is_iterable( markers , exclude_string =True , transform =True ) 
        markers = list(m) + [ 'o', 'D']
        
        #==plotlog10 --------
        # to use frequency for individual site rather than 
        # the complete frequency 
        fp = self.p_.ediObjs_[0].Z._freq
        fp =  1/ fp if scale =='period' else fp 
    
        # min_y =  np.nanmin(res[0][:, site])
        # add error bar data to main 
        data = [res, phs ] 
        data +=  [ res_err ,  phs_err ] if errorbar else []
        
        sy=  np.average ( [
            ( np.nanmax(res[0][:, ii]) - np.nanmin(res[0][:, ii])) /2  
            for ii in sites ] )
        sy += spad 
        sx = np.average (fp)
        
        fig = plt.figure(figsize = self.fig_size, dpi=self.fig_dpi)
        gs = GridSpec(2, len(sites),
                        wspace=0., 
                        left=.08,
                        top=.85,
                        bottom=.1,
                        right=.98,
                        hspace=.0,
                        height_ratios=[ 1.5, 1.]
                        ) 
        sharey = None
        
        #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
        for j , site in enumerate (sites ): 
            
            ax1 = fig.add_subplot (gs [ 0, j] , sharey = sharey) 
            
            if j==0: sharey = ax1 
            
            if errorbar: 
                ax2 = fig.add_subplot (gs [1,  j] )
                
            for i, sloop in enumerate (zip (* data )) : 
                r, p, *sl = sloop 
                
                if len(sl) !=0 : 
                    e, ep = sl  # mean errorbar is set to True 
                    
                y =  reshape (r[:, site])
                if errorbar: 
                    plot_errorbar (ax1 , 
                                   fp, 
                                   y,  
                                   y_err = reshape (e[:, site]),
                                   )
                ax1.scatter (fp  , y, 
                              marker =markers [i], 
                              color =colors [i],
                              edgecolors='k', 
                              label = fr'{survey}$\rho_a${components[i]}',
                              **kws 
                              ) 
                # res_corr 
                ax1.scatter (fp, zc_res [site], 
                              marker ='*', 
                              color="#FF00FF",
                              edgecolors='k', 
                              label = fr'{survey}$\rho_a${components[i]} {mc}',
                              **kws 
                              ) 
                
                if errorbar: 
                    plot_errorbar (ax2 , 
                                   fp, reshape (p[:, site]),
                                   y_err = reshape (ep[:, site]), 
                                   )
                
                ax2.scatter( fp, 
                            reshape (p[:, site]),
                            marker =markers[i], 
                            color =colors [i] ,
                            edgecolors='k', 
                            label = f'{survey}$\phi${components[i]}',
                            **kws
                            ) 
                # ----phase_cor 
                ax2.scatter( fp, 
                            zc_phase [site],
                            marker ='*', 
                            color="#FF00FF" ,
                            edgecolors='k', 
                            label = f'{survey}$\phi${components[i]} {mc}',
                            **kws
                            ) 
                
            # set ticks invisibale 
            if j > 0: 
                plt.setp(ax1.get_yticklabels(), visible=False)
                plt.setp(ax2.get_yticklabels(), visible=False)
                
            # Put the legend in the last images 
            
            if j == len(sites)-1: 
                try: 
                    ax1.legend(ncols = len(r)) 
                    ax2.legend(ncols = len(p)) 
                except: 
                    # For consistency in the case matplotlib  is < 3.3. 
                    ax1.legend() 
                    ax2.legend() 

            ax1.set_xscale ('log') ;  ax1.set_yscale ('log') 
            ax2.set_xscale ('log')
            
            if show_site:
                ax1.text (sx, 
                          sy, 
                          f'S{site}', 
                          fontdict= dict (style ='italic',  bbox =dict(
                               boxstyle='round',facecolor ='#CED9EF'), 
                              alpha = 0.5 )
                          )
            
            ax2.set_ylim ([0, 90 ])
            xlabel = self.xlabel or ( 'Log$_{10}$Period($s$)' if scale=='period' 
                                     else 'Frequency ($H_z$)') 

            # fixing yticks with matplotlib.ticker "FixedLocator"
            # xticks_loc =  ax2.get_xticks()
            # ax2.xaxis.set_major_locator(mticker.FixedLocator(xticks_loc))
            # ax2.set_xticklabels(['{:,.0f}'.format(np.log10(x)) 
            #                      for x in xticks_loc])
            ax2.set_xlabel(xlabel ) 
  
            if j ==0 : 
                
                ax1.set_ylabel(self.ylabel or r'$\rho_a$($\Omega$.m)') 
                ax2.set_ylabel('$\phi$($\degree$)')
                
            if self.show_grid :
                for ax in (ax1, ax2 ): 
                    ax.grid (visible =True , alpha =self.galpha, 
                             linestyle = self.gls, 
                             which =self.gwhich, color =self.gc
                             )
    
            
        if self.savefig is not  None: 
            plt.savefig (self.savefig, dpi = self.fig_dpi)
            
        plt.close () if self.savefig is not None else plt.show() 
        
        return self 
    

    def __repr__(self): 
        """ Represents the output class format """
        outm = ( '<{!r}:' + ', '.join(
            [f"{k}={getattr(self, k)!r}" for k in self._t]) + '>' 
            ) 
        return  outm.format(self.__class__.__name__)
    
    
TPlot.__doc__="""\
Tensor plot from EMAP or AMT processing data.

`TPlot` is a :term:`Tensor` (Impedances , resistivity and phases ) plot class. 
Explore SEG ( Society of Exploration Geophysicist ) class data.  Plot recovery 
tensors. `TPlot` methods returns an instancied object that inherits 
from :class:`watex.property.Baseplots` ABC (Abstract Base Class) for 
visualization.
    
Parameters 
------------

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
   
mode: str , ['valid', 'same'], default='same'
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
    resistivity ``srho``. 
    
c : int, 
    A window-width expansion factor that must be input to the filter 
    adaptation process to control the roll-off characteristics
    of the applied Hanning window. It is recommended to select `c` between 
    ``1``  and ``4``.  Default is ``2``.
    
distance: float 
    The step between two stations/sites. If given, it creates an array of  
    position for plotting purpose. Default value is ``50`` meters. 
 
prefix: str 
    string value to add as prefix of given id. Prefix can be the site 
    name. Default is ``S``. 
    
how: str 
    Mode to index the station. Default is 'Python indexing' i.e. 
    the counting of stations would starts by 0. Any other mode will 
    start the counting by 1.
     
{params.base.savefig}
{params.base.fig_dpi}
{params.base.fig_num}
{params.base.fig_size}
{params.base.fig_orientation}
{params.base.fig_title}
{params.base.fs}
{params.base.ls}
{params.base.lc}
{params.base.lw}
{params.base.alpha}
{params.base.font_weight}
{params.base.font_style}
{params.base.font_size}
{params.base.ms}
{params.base.marker}
{params.base.marker_facecolor}
{params.base.marker_edgecolor}
{params.base.marker_edgewidth}
{params.base.xminorticks}
{params.base.yminorticks}
{params.base.bins}
{params.base.xlim}
{params.base.ylim}
{params.base.xlabel}
{params.base.ylabel}
{params.base.rotate_xlabel}
{params.base.rotate_ylabel}
{params.base.leg_kws}
{params.base.plt_kws}
{params.base.glc}
{params.base.glw}
{params.base.galpha}
{params.base.gaxis}
{params.base.gwhich}
{params.base.tp_axis}
{params.base.tp_labelsize}
{params.base.tp_bottom}
{params.base.tp_labelbottom}
{params.base.tp_labeltop}
{params.base.cb_orientation}
{params.base.cb_aspect}
{params.base.cb_shrink}
{params.base.cb_pad}
{params.base.cb_anchor}
{params.base.cb_panchor}
{params.base.cb_label}
{params.base.cb_spacing}
{params.base.cb_drawedges} 
{params.sns.sns_orient}
{params.sns.sns_style}
{params.sns.sns_palette}
{params.sns.sns_height}
{params.sns.sns_aspect}

Returns
--------
{returns.self}

Examples
---------
>>> from watex.view.plot import TPlot 
>>> from watex.datasets import load_edis 
>>> plot_kws = dict( ylabel = '$Log_{{10}}Frequency [Hz]$', 
                    xlabel = '$Distance(m)$', 
                    cb_label = '$Log_{{10}}Rhoa[\Omega.m$]', 
                    fig_size =(6, 3), 
                    font_size =7., 
                    rotate_xlabel=45, 
                    imshow_interp='bicubic', 
                    ) 
>>> edi_data =load_edis (return_data= True, samples=7 ) 
>>> t= TPlot(**plot_kws ).fit(edi_data)
>>> t.fit(edi_data ).plot_tensor2d (to_log10=True )
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
|Data collected =  7      |EDI success. read=  7      |Rate     =  100.0  %|
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Out[150]: <AxesSubplot:xlabel='$Distance(m)$', ylabel='$Log_{{10}}Frequency [Hz]$'>
""".format(
    params=_param_docs,
    returns= _core_docs["returns"],
)