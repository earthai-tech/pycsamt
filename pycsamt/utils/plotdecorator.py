# -*- coding: utf-8 -*-
# This module is part of the pyCSAMT viewer package
#       Author: Kouadio K.Laurent<etanoyau@gmail.com>
#       Licence: LGPL

"""
Created on Sat Aug 21 15:55:48 2021

.. sypnosis:: Decorator module to plot 2d and 1d inversion of geostrata model. 
            ``reason parameter specify the kind of plot and can be for::
                - `model` : plot model inversion or geomodel 
                - `misfit` :plot error between NM and occam inversion results. 
                - `zonge` or `avg`: plot typical resistivity and phase values 
                - `response`: plot fitting curves of inversion rho and phase
"""
import warnings
import functools 
import numpy as np 

import matplotlib as mpl 
import  matplotlib.pyplot  as plt
import matplotlib.cm as cm 
import matplotlib.colorbar as mplcb
import matplotlib.gridspec as gspec
# from mpl_toolkits.axes_grid1 import make_axes_locatable
# from matplotlib.ticker import MultipleLocator, NullLocator

from pycsamt.utils import plot_utils as mplotus  
from pycsamt.utils._csamtpylog import csamtpylog 

_logger=csamtpylog.get_csamtpy_logger(__name__)


class geoplot1d : 
    """ 
 
    Decorator 1d class to plot response and typicall apparent resistivity 
    and errors misfit. It call 'PlotResponse' objects  from :mod:`~occam2d` 
    and :meth:`~plot_dataAndFits` from typycall apprent resistivity from 
    GDPII-multifunctional receiver or zonge internation company. `Plot_style`
    depend on argument `reason`.

    :param reason: type of plot , can be `avg` or `response`
                    if none , will plot `zonge` is triggered.
    :type reason: str 
    
    :param kws: Matplotlib properties and model properties
    :type kws: dict
    
    ======================== ==================================================
    Attributes               Description
    ======================== ==================================================
    color_mode               [ 'color' | 'bw' ] color or black and white plots
    cted                     color for data `app.res` 
    ctem                     color for model `app.res` 
    ctmd                     color for data `phase`
    ctmm                     color for model `phase`
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
    mted                     marker for data `fitting curve rho`
    mtem                     marker for model `fitting curve phase`
    mtmd                     marker for data `response rho & phase`
    mtmm                     marker for model `response rho & phase`
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
    subplot_bottom           space between axes and bottom of figure
    subplot_hspace           space between subplots in vertical direction
    subplot_left             space between axes and left of figure
    subplot_right            space between axes and right of figure
    subplot_top              space between axes and top of figure
    subplot_wspace           space between subplots in horizontal direction
    ======================== ==================================================
    """
    
  
    def __init__(self, reason =None, **kwargs): 
        
        self._logging =csamtpylog.get_csamtpy_logger(self.__class__.__name__)
        
        self.reason = reason 
        
        self.gridspec = kwargs.pop('gridspec', 4)
        self.xlabel =kwargs.pop('x_label', 'freq')
        
        self.fig_num= kwargs.pop('fig_num', 1)
        self.fig_size = kwargs.pop('fig_size', [6,6])
        self.fig_dpi =kwargs.pop('fig_dpi', 300)
        
        self.fig_title =kwargs.pop('title', None)
        self.savefig = kwargs.pop('savefig', None)
        
        self.x_minorticks=kwargs.pop('xminorticks', 1)
        self.y_minorticks =kwargs.pop('yminorticks', 1)
        
        self.font_size =kwargs.pop('font_size',3.)
        self.font_style=kwargs.pop('font_style', 'italic')
        self.fs =kwargs.pop('fs', 2.)
        
        self.fontw =kwargs.pop('font_weight', 'bold')
 
        self.alpha = kwargs.pop('alpha', 0.5)
        
        self.ms = kwargs.pop('ms', 1.5)
        self.ms_r = kwargs.pop('ms_r', 3)
        self.lw = kwargs.pop('lw', .5)
        self.lw_r = kwargs.pop('lw_r', 1.0)
        self.ls = kwargs.pop('ls',':')
        self.e_capthick = kwargs.pop('e_capthick', .5)
        self.e_capsize = kwargs.pop('e_capsize', 2)
        
        self.subplot_wspace = kwargs.pop('subplot_wspace', .3)
        self.subplot_hspace = kwargs.pop('subplot_hspace', .0)
        self.subplot_right = kwargs.pop('subplot_right', .98)
        self.subplot_left = kwargs.pop('subplot_left', .08)
        self.subplot_top = kwargs.pop('subplot_top', .85)
        self.subplot_bottom = kwargs.pop('subplot_bottom', .1)
    
        self.legend_loc = 'upper center'
        self.legend_pos = (.5, .95) #(.5, 1.18)
        self.legend_pos_tipper = (.5, 1.18)
        self.legend_marker_scale = 1
        self.legend_border_axes_pad = .01
        self.legend_label_spacing = 0.07
        self.legend_handle_text_pad = .2
        self.legend_border_pad = .15
        self.legendProps = kwargs.pop('kw_legendprop', 
                                      {'size':1.5*self.font_size, 
                                     'style': self.font_style})


        self.color_mode = kwargs.pop('color_mode', 'color')
        self.plot_style = kwargs.pop('plot_style', 1)
        
       # color for data
        # color mode
        if self.color_mode == 'color':
            # color for data
            self.cted = kwargs.pop('cted', (0, 0, 1))
            self.ctmd = kwargs.pop('ctmd', (1, 0, 0))
            self.mted = kwargs.pop('mted', 's')
            self.mtmd = kwargs.pop('mtmd', 'o')

            # color for occam2d model
            if self.plot_style == 3:
                # if plot_style is 3, set default color
                #for model response to same as data
                self.ctem = kwargs.pop('ctem',self.cted)
                self.ctmm = kwargs.pop('ctmm',self.ctmd)
            else:
                self.ctem = kwargs.pop('ctem', (0, .6, .3))
                self.ctmm = kwargs.pop('ctmm', (.9, 0, .8))
            self.mtem = kwargs.pop('mtem', '+')
            self.mtmm = kwargs.pop('mtmm', '+')

        # black and white mode
        elif self.color_mode == 'bw':
            # color for data
            self.cted = kwargs.pop('cted', (0, 0, 0))
            self.ctmd = kwargs.pop('ctmd', (0, 0, 0))
            self.mted = kwargs.pop('mted', 's')
            self.mtmd = kwargs.pop('mtmd', 'o')

            # color for occam2d model
            self.ctem = kwargs.pop('ctem', (0.6, 0.6, 0.6))
            self.ctmm = kwargs.pop('ctmm', (0.6, 0.6, 0.6))
            self.mtem = kwargs.pop('mtem', '+')
            self.mtmm = kwargs.pop('mtmm', 'x')
        
        self.h_ratio = kwargs.pop('h_ratio',[1.5, 1, .5])  
        self.show_grid=kwargs.pop('show_grid', True)
        self.grid_alpha =kwargs.pop('alpha', .5)

        self.freqOrperiod_limits=kwargs.pop('freqOrperiod_limits', None)
        
        # kwargs bbbox 
        self.kw_linebbox = kwargs.pop('linebbox_kws', {'boxstyle':'round',
                                               'facecolor': 'whitesmoke'})
        self.kw_figsup_bbox = kwargs.pop('figbbox_kws',
                                               {'boxstyle':'round',
                                               'facecolor': 'moccasin'} )
        self.kw_grid= kwargs.pop('grid_kws', 
                                 {'color':'k', 'ls':':', 'lw':0.5, 
                                  'alpha':self.grid_alpha,
                                  'which':'major'})
    def __call__(self, func ): 
        """ Call function to plot """
        self._func = func 
        
        @functools.wraps(func)
        def fwrap_1d(*args, **kwargs):
            """ Function to be decorated. """

            # --> set default font size
            # plt.rcParams['font.size'] = self.font_size
            fontdict = {'size': self.font_size, 'weight': 'bold'}
            
            # set all legend properties on dictionnary 
            kw_legend_props ={'loc':self.legend_loc , 
                            'bbox_to_anchor':self.legend_pos,
                            'markerscale':self.legend_marker_scale,
                            'borderaxespad':self.legend_border_axes_pad,
                            'labelspacing':self.legend_label_spacing,
                            'handletextpad':self.legend_handle_text_pad,
                            'borderpad':self.legend_border_pad,
                               'prop': self.legendProps
                                   }
            # --> make key word dictionaries for plotting
            kw_xx = {'color': self.cted,
                     'marker': self.mted,
                     'ms': self.ms,
                     'ls': ':',
                     'lw': self.lw,
                     'e_capsize': self.e_capsize,
                     'e_capthick': self.e_capthick}
    
            kw_yy = {'color': self.ctmd,
                     'marker': self.mtmd,
                     'ms': self.ms,
                     'ls': ':',
                     'lw': self.lw,
                     'e_capsize': self.e_capsize,
                     'e_capthick': self.e_capthick}
            
            kw_yp = {'color': self.ctmm,
                    'marker': self.mtmm,
                    'ms': self.ms_r,
                    'ls': ':',
                    'lw': self.lw_r,
                    'e_capsize': self.e_capsize,
                    'e_capthick': self.e_capthick}

            
             # make conditions ----------------------------------------
            station_axis =False 
            if self.reason.lower().find('resp')>=0 : 
                
                lines_id, station, freq, appRho, phase, appRho_err,\
                    phase_err, model_RMS= func(*args, **kwargs)
 
                # arramge the model RMS
                if model_RMS is None: 
                    model_RMS =[' ' for i in range(len(lines_id))]
                elif isinstance(model_RMS, str): 
                    model_RMS=[model_RMS]
                if isinstance(model_RMS, list):
                    if len(model_RMS) < len(lines_id): 
                        model_RMS = model_RMS +[' ' 
                                        for i in range(4-len(lines_id)) ]
                    elif len(model_RMS) > len(lines_id): 
                        model_RMS =model_RMS[:4]

                self.reason = 'resp'
                fmsup ='Fitting curves '
                
            elif  self.reason.lower().find('avg') >=0 or \
                self.reason.lower().find('zonge')>=0 :
                self._logging.info('Plot Zonge `avg` file. ')
                lines_id, station, freq, appRho, phase, appRho_err,\
                    phase_err= func(*args, **kwargs)
                    
                self.reason ='zonge'
                fmsup ='Measured'
            elif self.reason.lower().find('stat') >=0 : 
        
                kind_p , lines_id, station, freq, appRho, phase, appRho_err,\
                    phase_err= func(*args, **kwargs)
                    
                self.reason ='sshift'
                fmsup = 'Static correction'
                # if station axis, then plot the semilog axis rho
                if kind_p[0] ==1: station_axis =True  
                ffilter= kind_p[1] # ama/flma or tma
                
                self._logging.info(
                    'Plot static correction of line {0} using {1} filter.'.
                    format(lines_id[0],ffilter.capitalize() ))
                
            n_stn = len(lines_id) 
    
            if self.reason in ['resp', 'zonge']:
                self._logging.info(
                    f"Plot `{' A single' if n_stn==1 else str(n_stn)+'lines'}` "
                    f"{'Occam Response' if self.reason =='resp' else 'zonge avg'}"
                    f" file{'s' if n_stn >1 else ''}.")
            
            
            if isinstance(station, str): 
                    station =[station]
            #---------------------------------------------
            fig = plt.figure(' '.join([str(stn) for stn in station]),
                             self.fig_size, dpi=self.fig_dpi)
            plt.clf()
            
            # Manage the subtitle 
            #--> Check whether all station name are the same and keep one 
            if len(set([l.lower() for l in station]))==1:
                
                stitle ='{0} data at station {1}'.format(fmsup, station[0])
            else : 
                fmt =''.join(['{'+'{0}'.format(ii) +'}, ' 
                              for ii in range(len(station))])
                # remode the last comma ad put dot
                fmt=fmt[:-2] +'.'
                stitle = '{0} data at stations {1}'.\
                    format(fmsup, fmt).format(*station)
                if self.reason =='resp':
                    stitle = stitle.replace('data', 'of Rho and Phase')
                if station_axis is True: 
                    stitle=stitle.replace('stations', 
                        f"{'frequencies' if n_stn>1  else 'frequency'}")
                if self.reason =='sshift':
                    stitle=stitle.replace('data', '')
                
            fig.suptitle(stitle, verticalalignment='center', 
                                  style =self.font_style,
                                  bbox =self.kw_figsup_bbox, 
                                  y=0.95,
                                   fontdict={'size': self.font_size*2 + 2,
                                             'weight': self.fontw}
                                  )
                              
            if self.reason == 'zonge': 
                self.gridspec = len(lines_id)
                
            gs = gspec.GridSpec(2,self.gridspec,
                                   wspace=self.subplot_wspace,
                                   left=self.subplot_left,
                                   top=self.subplot_top,
                                   bottom=self.subplot_bottom,
                                   right=self.subplot_right,
                                   hspace=self.subplot_hspace,
                                   height_ratios=self.h_ratio[:2]) 
                
            # create empty list hold the axis rho, axis phase and legend  
            rho_axis_col =[]
            resp_rho_col =[]
            phase_axis_col =[]
            legend_ax_list =[]
            # plot_type on xlabel 
            if self.xlabel.lower().find('peri')>=0 or self.xlabel ==2: 
                freqOrPeriod = [1/f for f in freq]
                freqOrperiod_list = 1/freq[0] 
                xlabel_name ='Periods(s)'
                
                
                if self.freqOrperiod_limits is None: 
   
                        self.freqOrperiod_limits = (10 ** (np.floor(np.log10(
                                                  freqOrperiod_list[0]))
                                                  ) * 1.01,
                          10 ** (np.ceil(np.log10(freqOrperiod_list[-1]))) * .99)
            else: 
                
                freqOrperiod_list = freq[0]
                freqOrPeriod = freq 
                xlabel_name ='Frequency (Hz)'
            
                 
                if self.freqOrperiod_limits is None: 
        
                    self.freqOrperiod_limits = (10 ** (np.ceil(np.log10(
                                              freqOrperiod_list[0]))
                                              ) * 1.01,
                      10 ** (np.floor(np.log10(freqOrperiod_list[-1]))) * .99)
                    
            if self.reason =='sshift' and station_axis is True: 
                xlabel_name ='Distance (m)'
  
                self.freqOrperiod_limits=(freq[0].min(), freq[0].max())
                        
     
            # plot error in percentage 
            if self.reason =='zonge': 
                coef = 100. 
                x_line_coef, y_line_coef =1.01 , 1.5
                y_pad =.5 
            elif self.reason =='resp': 
                coef =1.
                x_line_coef, y_line_coef =1.01 , 1.01
                y_pad =.3 
            elif self.reason =='sshift': 
                coef =1.
                x_line_coef, y_line_coef =50. , 1.1
                
                if station_axis is True: 
                    x_line_coef *=10 
                    y_line_coef *=3.5
                    
                y_pad =.5 
                
            if len(lines_id)>=1:
                if station_axis is True: 
                    x_scale =None 
                else: x_scale = 'log'
                
                ax__appRho1 = fig.add_subplot(gs[0, 0],
                                              xscale=x_scale ,
                                              yscale='log',
                                              xlim=self.freqOrperiod_limits)
                
                ax__Phase1= fig.add_subplot(gs[1, 0],xscale =x_scale , 
                                            sharex=ax__appRho1) 
                
                #col1
                err_rho = mplotus.plot_errorbar(ax__appRho1,
                                                freqOrPeriod[0],
                                                appRho[0][0],
                                                appRho_err[0]*coef,
                                                **kw_xx)
                err_phase = mplotus.plot_errorbar(ax__Phase1, 
                                                      freqOrPeriod[0], 
                                                      phase[0][0] , 
                                                      phase_err[0]*coef, 
                                                      **kw_yy)
                
                rho_axis_col.append(ax__appRho1)
                phase_axis_col.append(ax__Phase1 )
                
                legend_ax_list.append((err_rho, err_phase))
                
                

            if len(lines_id)>=2:
                ax__appRho2 = fig.add_subplot(gs[0, 1],
                                              xscale=x_scale ,
                                              yscale='log', 
                                              sharey = ax__appRho1)
                ax__Phase2= fig.add_subplot(gs[1, 1], xscale =x_scale ,
                                            sharex=ax__appRho2, 
                                            sharey=ax__Phase1) 
                #col2 

        
                err_rho = mplotus.plot_errorbar(ax__appRho2,
                                                freqOrPeriod[1],
                                                appRho[1][0],
                                                appRho_err[1]*coef,
                                                **kw_xx)
                err_phase = mplotus.plot_errorbar(ax__Phase2, 
                                                  freqOrPeriod[1], 
                                                      phase[1][0] , 
                                                      phase_err[1]*coef, 
                                                      **kw_yy)   

                rho_axis_col.append(ax__appRho2)
                phase_axis_col.append(ax__Phase2)
                legend_ax_list.append((err_rho, err_phase))
                
            if len(lines_id)>=3:
                ax__appRho3 = fig.add_subplot(gs[0, 2],
                                              xscale=x_scale ,
                                              yscale='log',
                                              sharey = ax__appRho1)
                ax__Phase3= fig.add_subplot(gs[1, 2], xscale =x_scale ,
                                            sharex=ax__appRho3, 
                                            sharey=ax__Phase1) 
                #col3 
 
        
                err_rho = mplotus.plot_errorbar(ax__appRho3,
                                                freqOrPeriod[2],
                                                appRho[2][0],
                                                appRho_err[2]*coef,
                                                **kw_xx)
                err_phase = mplotus.plot_errorbar(ax__Phase3, 
                                                  freqOrPeriod[2], 
                                                      phase[2][0] , 
                                                      phase_err[2]*coef, 
                                                      **kw_yy)
                rho_axis_col.append(ax__appRho3)
                phase_axis_col.append(ax__Phase3)
                legend_ax_list.append((err_rho, err_phase))
                
            if len(lines_id)==4: 
                
                ax__appRho4 = fig.add_subplot(gs[0, 3], 
                                              xscale=x_scale ,
                                              yscale='log',
                                              sharey = ax__appRho1)
                ax__Phase4= fig.add_subplot(gs[1, 3], xscale =x_scale ,
                                            sharex=ax__appRho4, 
                                            sharey=ax__Phase1) 
                
                #col4 
        
                err_rho = mplotus.plot_errorbar(ax__appRho4,
                                                freqOrPeriod[3],
                                                appRho[3][0],
                                                appRho_err[3]*coef,
                                                **kw_xx)
                err_phase = mplotus.plot_errorbar(ax__Phase4, 
                                                  freqOrPeriod[3], 
                                                      phase[3][0] , 
                                                      phase_err[3]*coef, 
                                                      **kw_yy)
            
                rho_axis_col.append(ax__appRho4)
                phase_axis_col.append(ax__Phase4)
                legend_ax_list.append((err_rho, err_phase))
            
            # plor response 
            if self.reason =='resp' or self.reason=='sshift': 
   
                for ii, (axis_rho, axis_phase) in enumerate(
                        zip(rho_axis_col, phase_axis_col )):

                    erro_rho_resp= mplotus.plot_errorbar(axis_rho,
                                            freqOrPeriod[0],
                                            appRho[ii][1],
                                            None,
                                            **kw_yp)
                    mplotus.plot_errorbar(axis_phase, 
                                            freqOrPeriod[0], 
                                            phase[ii][1] , 
                                            None, 
                                            **kw_yp)
                    resp_rho_col.append(erro_rho_resp)
            
            
            # rename line if lines not in name : 
            lines_id =['Line {0}'.format(name) for  name in lines_id 
                       if name.find('line') <0]
            lines_id =[name.replace('K','0') for name in lines_id]
            if self.reason =='sshift': 
                if kind_p[0] ==1: 
                    lines_id =[' Freq:{0}'.format(station[ii])
                            for ii, name in enumerate(lines_id)]
            else: lines_id =[name +' - Site {0}'.format(station[ii])
                        for ii, name in enumerate(lines_id)]
            # --> Take the maximum in heigth y for the appRho plot 
            max_y=-999
            for ii , ap in enumerate(appRho): 
                t_y = np.log10(ap[0]).max()
                if t_y > max_y: 
                    max_y = t_y
                        
            for ii, (ax_rho , ax_phase) in enumerate(
                    zip(rho_axis_col, phase_axis_col)): 
   
                if self.show_grid is True:
                    ax_rho.minorticks_on()
                    ax_rho.grid(**self.kw_grid)

                    ax_phase.minorticks_on()
                    ax_phase.grid(**self.kw_grid)

                
                if self.reason =='zonge': 
                    
                    ax_rho.legend(legend_ax_list[ii][0] ,  ['App.res(Ω.m)'],
                                  **kw_legend_props 
                                       )
                    ax_phase.legend(legend_ax_list[ii][1] , 
                                ['phase(degrees)'],
                                **kw_legend_props 
                               )
                
                elif self.reason =='resp' : 
        
                    ax_rho.legend([legend_ax_list[ii][0], resp_rho_col[ii]],
                                 ['App.res(Ω.m)', 'Rho rms={}'.format(
                                     model_RMS[ii])], **kw_legend_props )
                elif self.reason =='sshift': 
                    ax_rho.legend([legend_ax_list[ii][0], resp_rho_col[ii]],
                                     ['App.res(Ω.m)', 'Rho, filter= {}'.format(
                                         ffilter)], **kw_legend_props ) 
   
                ax_phase.set_xlabel(xlabel_name,  
                          fontdict= fontdict
                          )
                ax_rho.set_xlim(self.freqOrperiod_limits)
                
                if ii ==0 :  
                    if self.reason =='resp' or self.reason =='sshift': 
                        ylabel_name = 'Resistivity (Ω.m)'
                    elif self.reason =='zonge': 
                        ylabel_name = 'Apparent resistivity (Ω.m)'
                        
                    ax_rho.set_ylabel(ylabel_name, 
                              fontdict=fontdict)
                    ax_phase.set_ylabel('Phase (degrees)', 
                                          fontdict=fontdict)
                 
               
                if station_axis is True : 
                    # take the mean value of stations x station-axis 
                    x_text = freq[0].mean()
                else: 
                    # display the text : survey line (.5, 1.18)
                    # taxe the max 7
                    x_text = (10 ** np.floor(np.log10(
                        self.freqOrperiod_limits).mean())) *x_line_coef
                
                y_text = (10 ** (max_y + y_pad ))* y_line_coef
                                 
                ax_rho.text(x_text,
                        y_text,  
                        s= lines_id[ii],
                        horizontalalignment='center',
                        verticalalignment='baseline',
                        fontdict={'size': 2* self.font_size, 
                                  'color': 'k', 
                                  'style':self.font_style},
                        bbox =self.kw_linebbox
                        # rotation = self.station_label_rotation,
                            )
                    
            # savefigure
            if self.savefig is not None : 
                plt.savefig(self.savefig, dpi = self.fig_dpi)
                    
            plt.show()
            
            return func(*args, **kwargs)
        
        return fwrap_1d
  
            
       

class geoplot2d(object): 
    """
    Decorator class to plot geodrill model and geodrill misfit . 
    It call geoModel objects and plot the model or misfit according the *reason`
    argument provided. see `pycsamt.geodrill.geoCore.gedrill.geoModel` to get
    an implementation example. 
    use also matplotlib properties to customize your plot
    
    :param reason: type of plot , can be `misfit` or `model`
                    if none , will plot `model`.
    :type reason: str 
    
    :param kws: Matplotlib properties and model properties
    :type kws: dict
    
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
    def __init__(self, reason =None , **kws): 
        
        self._logging= csamtpylog.get_csamtpy_logger(self.__class__.__name__)

        self.reason =reason 
        
        self.fs =kws.pop('fs', 0.7)
        self.fig_num = kws.pop('fig_num', 1)
        self.fig_size = kws.pop('fig_size', [7,7])
        self.fig_aspect = kws.pop('fig_aspect','auto')
        self.fig_dpi =kws.pop('fig_dpi', 300)
        self.font_size = kws.pop('font_size', 7)
        self.aspect = kws.pop('aspect', 'auto')
        self.font_style =kws.pop('font_style', 'italic')
        self.orient=kws.pop('orientation', 'landscape')

        self.cb_pad = kws.pop('cb_pad', .0375)
        self.cb_orientation = kws.pop('cb_orientation', 'vertical')
        self.cb_shrink = kws.pop('cb_shrink', .75)
        self.cb_position = kws.pop('cb_position', None)
        self.climits = kws.pop('climits', (0, 4))

        self.station_label_rotation = kws.pop('station_label_rotation',45)
        self.imshow_interp = kws.pop('imshow_interp', 'bicubic')
        self.ms = kws.pop('ms', 2)
        self.lw =kws.pop('lw', 2)
  
        self.fw =kws.pop('font_weight', 'bold')
 
        self.station_font_color = kws.pop('station_font_color', 'k')
        self.station_marker = kws.pop('station_marker',
                                         r"$\blacktriangledown$")
        self.station_color = kws.pop('station_color', 'k')

        self.xpad = kws.pop('xpad', 1.0)
        self.ypad = kws.pop('ypad', 1.0)

        self.cmap = kws.pop('cmap', 'jet_r')
    
        self.depth_scale =kws.pop('depth_scale', None)
        self.doi = kws.pop('doi', 1000)
        
        self.savefig =kws.pop('savefig', None)
        self.change_station_id =kws.pop('new_station_names', None)

        self.model_rms =kws.pop('model_rms', None)
        self.model_roughness =kws.pop('model_roughness', None)
        self.plot_style =kws.pop( 'plot_style', 'pcolormesh') 
        self.show_contour =kws.pop('show_contour', False)
        self.contourlines =kws.pop('contour_lines_styles', '-')
        self.contourcolors =kws.pop('contour_lines_colors', 'white')
        self.delineate_resistivity_curve =kws.pop('delineate_rho', None)
        self.grid_alpha =kws.pop('alpha', 0.5)
        self.show_grid = kws.pop('show_grid',True)
   

        self.set_station_label=kws.pop('show_station_id', True)
        
        for keys in list(kws.keys()): 
            setattr(self, keys, kws[keys])
            
    def additional_tools(self):
        """
        Method to add additionnal tools like `computing delineate
        resistivity curve` . Will populate with other tools to 
        customize plots.
        
        """
         #---> get delineate rho curve  : value of resistivities must be on
         #MUST be on OHM- M not log10 resistivity . 
         # check the resistivity value and  round to decimal 1
        if self.delineate_resistivity_curve is not None : 
            self.delineate_resistivity_curve = mplotus.controle_delineate_curve(
                res_deline=self.delineate_resistivity_curve)
            try : 
                # for consistency , check whether value are on list 
                if not isinstance (self.delineate_resistivity_curve, list ): 
                    self.delineate_resistivity_curve=[
                        self.delineate_resistivity_curve ]
            except :
                pass 
        
        return self.delineate_resistivity_curve 
 
    
    def __call__(self, func):  
        """
        Model decorator to hold the input function with arguments 
        :param func: function to be decorated 
        :type func: object 
        """

        return self.plot2DModel(func)
        
    def plot2DModel(self, func):
        @functools.wraps(func)
        def new_func (*args, **kwargs): 
            """
            new decorated function . Plot model data and misfit data 
            
            :args: arguments of  function  to be decorated 
            :type args: list 
        
            :param kwargs: positional arguments of decorated function
            :type kwargs: dict 
            :return: function decorated after visualisation
      
            """
            self._logging.info(
                ' Plot decorated {0}.'.format(func.__name__))
    
            _f=0 # flag to separated strata model misfit and occam model misfit
                #   from occamResponse file 
                
            if self.depth_scale is not None :
                self.depth_scale= str(self.depth_scale).lower() 

            if self.depth_scale not in ["km", "m"]: 
                mess ="Depth scale =`{}` is unacceptable value."\
                    " Should be convert to 'm'.".format(self.depth_scale)
                warnings.warn(mess)
                self.depth_scale= "m"
                self._logging.debug (mess)
            
            if self.depth_scale == 'km':
                dz  = 1000.
            elif self.depth_scale == 'm': # for CSAMT , we use default as meter"m".
                dz = 1.
    
            if self.delineate_resistivity_curve is not None : 
                self.delineate_resistivity_curve= self.additional_tools()
                
            # figure configuration 
            
            self.fig = plt.figure(self.fig_num, self.fig_size, dpi=self.fig_dpi)
            plt.clf()
            self.fig_aspect ='auto'
            axm = self.fig.add_subplot(1, 1, 1, aspect=self.fig_aspect)
    
            # get geomodel data 
            if self.reason is None: 
                self.reason = 'model'# by default
                
            # ----populate special attributes from model or misfit ------------
            if self.reason =='model': 
                occam_model_resistiviy_obj, occam_data_station_names, *m = func(
                    *args, **kwargs)
                occam_data_station_offsets, occam_model_depth_offsets, *ddrms= m
                self.doi, self.depth_scale, self.model_rms, *rmisf = ddrms 
                self.model_roughness, plot_misfit = rmisf
                
                self.doi = occam_model_depth_offsets.max()
                #     self.doi = occam_model_depth_offsets.max()
                # --> check doi value provided, and convert to default unit {meters}  
                self.doi =mplotus.depth_of_investigation(doi=self.doi)
                
                # set boundaries of stations offsets and depth 
                spec_f = -(self.doi/5)/dz  # assume that depth will  start by 
                #0 then substract add value so 
                # to get space for station names text
                if self.climits is None :
                    self.climits =(0,4)  
                if plot_misfit is True : 
                    _f=2
                self.ylimits =(spec_f, self.doi/dz)  
                
            if self.reason =='misfit': 
                occam_model_resistiviy_obj, occam_data_station_names, *m = func(
                    *args, **kwargs)     
                occam_data_station_offsets, occam_model_depth_offsets, *rg=m
                self.model_rms, self.model_roughness= rg
                
                # check if "plotmisfit refers to 'geoStrata model 'geodrill
                # module then keep the doi and set `spec_f
                if 'geodtype' in list(kwargs.keys()): 
                    # means plot `misfit` from geostrata model
                    spec_f = -(self.doi/5)/dz 
                    self.ylimits =(spec_f, self.doi/dz) 
                    
                else :
                    # frequency are in log10 new doi is set according to 
                    #  value of frequency . specf _f =0 
                    # occam misfit from response file 
                    #= occam_model_depth_offsets= resp_frequencies in 
                    # log10 resistivities 
                    # resp_data = occam_model_resistiviy_obj,
                    # occam_data_station_offsets = station location or positions 
                    self.doi =occam_model_depth_offsets.max()
                    #spec_f = (self.doi/5)/dz # o.8
                    spec_f = - 0.
                    _f=1 
       
                    self.ylimits = (self.doi, occam_model_depth_offsets.min())
            
            #------------- manage stations and climits ------------------------  
            occam_data_station_offsets =np.array(occam_data_station_offsets)
            # station separation and get xpad . ex dl=50 then xpad =25 
            dl = occam_data_station_offsets.max()/ (len(
                occam_data_station_offsets)-1)
            self.xpad = (dl/2)/dz 
    
            if self.change_station_id is not None :  # call fonction to build a nu
                occam_data_station_names , mess= mplotus.build_new_station_id(
                    station_id =occam_data_station_names ,
                    new_station_name =self.change_station_id )
                self._logging.debug(mess)
                                   
            self.xlimits=(occam_data_station_offsets.min()/dz -self.xpad  , 
                      occam_data_station_offsets.max()/dz + self.xpad )
            
            # configure climits 
            if self.reason =='misfit':
                if self.climits is None : 
                        self.climits =(-3, 3)
                        
                elif 'min' in self.climits or 'max' in self.climits : 
                            self.climits = (occam_model_resistiviy_obj.min(), 
                                            occam_model_resistiviy_obj.max())
    
             
            if _f==2 : 
                self.reason = 'misfit' 
            self._logging.info ('Ready to plot {0}'
                                ' with matplotlib "{1}" style.'.
                                format(self.reason, self.plot_style))  
            
            # -------------- check dimensionnality ---------------------------
            occam_model_resistiviy_obj, *dm= self._check_dimensionality (
                        occam_model_resistiviy_obj,
                        occam_model_depth_offsets,
                          occam_data_station_offsets
                          )
            occam_model_depth_offsets, occam_data_station_offsets = dm

            
            if self.plot_style.lower() =='pcolormesh':
                mesh_x  , mesh_z= np.meshgrid(occam_data_station_offsets,
                                              occam_model_depth_offsets )
     
                vmin = self.climits[0]
                vmax = self.climits[1] 
    
                axm.pcolormesh (mesh_x/dz  , 
                                mesh_z/dz ,
                                  occam_model_resistiviy_obj,
                                      vmin = vmin,
                                      vmax = vmax,  
                                      shading= 'auto', 
                                      cmap =self.cmap, 
                                      alpha = None, 
                                     
                                      )
         
                if self.show_contour is True :
                    contps = axm.contour(mesh_x/dz  , mesh_z /dz ,
                                          occam_model_resistiviy_obj,
                                          colors =self.contourcolors,
                                          linestyles=self.contourlines)
                    if  self.delineate_resistivity_curve is not None : 
                        axm.clabel(contps,  self.delineate_resistivity_curve ,
                                        inline=True, fmt='%1.1f',
                                        fontsize =self.font_size,
                                          )
       
            if self.plot_style.lower() =='imshow': 
    
                mesh_x  , mesh_z= np.meshgrid(occam_data_station_offsets,
                                              occam_model_depth_offsets  )
    
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
    
                if self.delineate_resistivity_curve is not None :
                    origin ='upper'
                    contps = axm.contour(occam_model_resistiviy_obj,
                                         colors =self.contourcolors, 
                                          vmax=self.climits[0],
                                          vmin = self.climits[1], 
                                            linestyles=self.contourlines, 
                                            extent =( self.xlimits[0],
                                                      self.xlimits[1],
                                                      self.ylimits[1], 
                                                      self.ylimits[0]),
                                            extend ='both',
                                            origin= origin ,  
                                            )
                    axm.clabel(contps, self.delineate_resistivity_curve ,
                                    inline=True, 
                                    fmt='%1.1f',
                                    fontsize =self.font_size,
                                              )
            # get colormap for making a colorbar 
            if type(self.cmap) == str:
                self.cmap = cm.get_cmap(self.cmap)
            
            axm.set_xlim( [self.xlimits[0],  self.xlimits[1]])
            axm.set_ylim ([self.ylimits[1], self.ylimits[0]]) 
    
            # create twin axis to set ticks to the top station
            axe2=axm.twiny()
            axe2.xaxis.set_visible(False) # let keep only the axe lines 
            #set axis and set boundaries 
            if self.reason =='model' or _f==2 : 
                ydown_stiteslbls = self.ylimits[0]/5 
                ydown_stationlbls = self.ylimits[0] -(self.ylimits[0]/3)
                xhorizontal_lbs = (occam_data_station_offsets.max()/dz)/2
    
            elif self.reason =='misfit': 
                ydown_stiteslbls = self.ylimits[0] + 0.1 * self.ylimits[1]
                ydown_stationlbls= self.ylimits[0] +\
                    self.ylimits[1]/self.ylimits[0]
                xhorizontal_lbs = (occam_data_station_offsets.max()- 
                                   occam_data_station_offsets.min())/2
               
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
                
                if self.set_station_label is True :  # then plot label id 
                    axm.text(offset/dz ,
                            ydown_stiteslbls,  
                            s= names,
                            horizontalalignment='center',
                            verticalalignment='baseline',
                            fontdict={'size': self.ms*3, 
                                      'color': self.station_color},
                            rotation = self.station_label_rotation,
                                )
         
               
            if self.set_station_label is True : 
                axm.text (xhorizontal_lbs, 
                            ydown_stationlbls,  
                            s= 'Stations',
                            horizontalalignment='center',
                            verticalalignment='baseline',
                            fontdict={'size': self.ms*5, 
                                      'color':'k', 
                                      'style': self.font_style,
                                      'weight': self.fw},
                            )
    
            #-------------------- manage grid and colorbar -------------------- 
            self.g2dgridandcbManager(axm, _f)
            #------------------------------------------------------------------
            # initialize the reason to keep the default reason  
            self.reason = None
            
            if self.savefig is not None : 
                plt.savefig(self.savefig, dpi = self.fig_dpi)
            
            plt.show()
            
            return func(*args, **kwargs)
        
        return new_func
    
    def _check_dimensionality(self, data, z, x):
        """ Check dimensionality of data and fix it"""

        def reduce_shape(Xshape, x, axis_name =None): 
            """ Reduce shape to keep the same shape"""
            mess ="`{0}` shape({1}) {2} than the data shape `{0}` = ({3})."
            ox = len(x) 
            dsh = Xshape 
            if len(x) > Xshape : 
                x = x[: int (Xshape)]
                self._logging.debug(''.join([
                    f"Resize {axis_name!r}={ox!r} to {Xshape!r}.", 
                    mess.format(axis_name, len(x),'more',Xshape)])) 
                                        
            elif len(x) < Xshape: 
                Xshape = len(x)
                self._logging.debug(''.join([
                    f"Resize {axis_name!r}={dsh!r} to {Xshape!r}.",
                    mess.format(axis_name, len(x),'less', Xshape)]))
                
            return int(Xshape), x 
        
        sz0, z = reduce_shape(data.shape[0],
                              x=z, axis_name ='Z')
        sx0, x =reduce_shape (data.shape[1], 
                              x=x, axis_name ='X')
        # resize theshape 
        # data  = np.resize(data, (sz0, sx0))
        data = data [:sz0, :sx0]
        
        return data , z, x 
                
                
    def g2dgridandcbManager(self, axm, _f=None) :
        """ Plot2d model by configure grid and colorbar. 
        :param axm: 2d axis plot 
        :param _f: resize flag; misfit =2 and model =1 """
        # put a grid on if set to True 
        if self.show_grid is True:
            axm.minorticks_on()
            axm.grid(color='k', ls=':', lw =0.5, 
                      alpha=self.grid_alpha, which ='major')
        

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
        
        if self.reason == 'model' : 
            cb.set_ticks(np.arange(int(self.climits[0]), 
                                   int(self.climits[1]) + 1))
            cb.set_ticklabels(['10$^{0}$'.format('{' + str(nn) + '}') 
                               for nn in np.arange(int(self.climits[0]),
                                          int(self.climits[1]) + 1)])
            
                                
        else : 
            cb.set_ticks(np.linspace(self.climits[0], self.climits[1],5))
            cb.set_ticklabels(['{0}'.format(str(round(nn,2))) for nn in
                                np.linspace(self.climits[0],
                                          self.climits[1],5)])
            cb.set_label('misfitvalue(%)',
                  fontdict={'size': self.font_size + 1, 
                            'weight': 'bold'})
       
        # set axes labels
        axm.set_xlabel('Distance ({0})'.format(self.depth_scale),
                      fontdict={'size': self.font_size + 2,
                                'weight': 'bold'})
        
        if self.reason =='misfit':
            if _f ==2: ylabel = 'Depth ({0})'.format(self.depth_scale)
            else : ylabel= 'Log10Frequency(Hz)'
            mesT ='Plot Misfit'
            
        elif self.reason =='model':
            ylabel = 'Depth ({0})'.format(self.depth_scale)
            mesT = 'Plot strata model' 
        
        axm.set_ylabel(ylabel,fontdict={
            'size': self.font_size + 2, 'weight': 'bold'})
       
       
        self.fig.suptitle('{0}- DataType = {1} :RMS={2}, Roughness={3}'.\
                          format(mesT, self.reason, self.model_rms, 
                                self.model_roughness),
                  ha='center',
          fontsize= 7* self.fs, 
          verticalalignment='center', 
          style =self.font_style,
          bbox =dict(boxstyle='round',facecolor ='moccasin'), 
          y=0.95 if self.reason =='model' else 0.98)
      
        return self 
    
                
              
                
              
                
              
                
              
                
              
                
              
                
              
                
              
                
              
                
        
        
        
        
        
        
        
        
        
        
        
        
        