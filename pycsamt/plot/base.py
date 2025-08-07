# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

from __future__ import annotations

from pathlib import Path
from abc import ABC, abstractmethod 
from collections import OrderedDict
from inspect import signature
import warnings
from typing import Iterable, Literal, overload

import matplotlib.pyplot as plt
import numpy as np 
        

class BasePlot(ABC): 
    r""" Base class  deals with Machine learning and conventional Plots. 
    
    The `BasePlot` can not be instanciated. It is build on the top of other 
    plotting classes  and its attributes are used for external plots.
    
    Hold others optional informations: 
        
    ==================  =======================================================
    Property            Description        
    ==================  =======================================================
    fig_dpi             dots-per-inch resolution of the figure
                        *default* is 300
    fig_num             number of the figure instance. *default* is ``1``
    fig_aspect          ['equal'| 'auto'] or float, figure aspect. Can be 
                        rcParams["image.aspect"]. *default* is ``auto``.
    fig_size            size of figure in inches (width, height)
                        *default* is [5, 5]
    savefig             savefigure's name, *default* is ``None``
    fig_orientation     figure orientation. *default* is ``landscape``
    fig_title           figure title. *default* is ``None``
    fs                  size of font of axis tick labels, axis labels are
                        fs+2. *default* is 6 
    ls                  [ '-' | '.' | ':' ] line style of mesh lines
                        *default* is '-'
    lc                  line color of the plot, *default* is ``k``
    lw                  line weight of the plot, *default* is ``1.5``
    alpha               transparency number, *default* is ``0.5``  
    font_weight         weight of the font , *default* is ``bold``.        
    ms                  size of marker in points. *default* is 5
    marker              style  of marker in points. *default* is ``o``.
    marker_facecolor    facecolor of the marker. *default* is ``yellow``
    marker_edgecolor    edgecolor of the marker. *default* is ``cyan``.
    marker_edgewidth    width of the marker. *default* is ``3``.
    xminorticks         minortick according to x-axis size and *default* is 1.
    yminorticks         minortick according to y-axis size and *default* is 1.
    font_size           size of font in inches (width, height)
                        *default* is 3.
    font_style          style of font. *default* is ``italic``
    bins                histograms element separation between two bar. 
                         *default* is ``10``. 
    xlim                limit of x-axis in plot. *default* is None 
    ylim                limit of y-axis in plot. *default* is None 
    xlabel              label name of x-axis in plot. *default* is None 
    ylabel              label name  of y-axis in plot. *default* is None 
    rotate_xlabel       angle to rotate `xlabel` in plot. *default* is None 
    rotate_ylabel       angle to rotate `ylabel` in plot. *default* is None 
    leg_kws             keyword arguments of legend. *default* is empty dict.
    plt_kws             keyword arguments of plot. *default* is empty dict
    plt_style           keyword argument of 2d style. *default* is ``pcolormesh``
    plt_shading         keyword argument of Axes pycolormesh shading. It can be 
                        ['flat'|'nearest'|'gouraud'|'auto'].*default* is 
                        'auto'
    imshow_interp       ['bicubic'|'nearest'|'bilinear'|'quadractic' ] kind of 
                        interpolation for 'imshow' plot. Click `interpol_imshow`_ 
                        to get furher details about the interpolation method. 
                        *default* is ``None``.
    rs                  [ '-' | '.' | ':' ] line style of `Recall` metric
                        *default* is '--'
    ps                  [ '-' | '.' | ':' ] line style of `Precision `metric
                        *default* is '-'
    rc                  line color of `Recall` metric *default* is ``(.6,.6,.6)``
    pc                  line color of `Precision` metric *default* is ``k``
    s                   size of items in scattering plots. default is ``fs*40.``
    cmap                matplotlib colormap. *default* is `jet_r`
    gls                 [ '-' | '.' | ':' ] line style of grid  
                        *default* is '--'.
    glc                 line color of the grid plot, *default* is ``k``
    glw                 line weight of the grid plot, *default* is ``2``
    galpha              transparency number of grid, *default* is ``0.5``  
    gaxis               axis to plot grid.*default* is ``'both'``
    gwhich              type of grid to plot. *default* is ``major``
    tp_axis             axis  to apply ticks params. default is ``both``
    tp_labelsize        labelsize of ticks params. *default* is ``italic``
    tp_bottom           position at bottom of ticks params. *default*
                        is ``True``.
    tp_top              position at the top  of ticks params. *default*
                        is ``True``.
    tp_labelbottom      see label on the bottom of the ticks. *default* 
                        is ``False``
    tp_labeltop         see the label on the top of ticks. *default* is ``True``
    cb_orientation      orientation of the colorbar. *default* is ``vertical``
    cb_aspect           aspect of the colorbar. *default* is 20.
    cb_shrink           shrink size of the colorbar. *default* is ``1.0``
    cb_pad              pad of the colorbar of plot. *default* is ``.05``
    cb_anchor           anchor of the colorbar. *default* is ``(0.0, 0.5)``
    cb_panchor          proportionality anchor of the colorbar. *default* is 
                        `` (1.0, 0.5)``.
    cb_label            label of the colorbar. *default* is ``None``.      
    cb_spacing          spacing of the colorbar. *default* is ``uniform``
    cb_drawedges        draw edges inside of the colorbar. *default* is ``False``
    cb_format           format of the colorbar values. *default* is ``None``.
    sns_orient          seaborn fig orientation. *default* is ``v`` which refer
                        to vertical 
    sns_style           seaborn style 
    sns_palette         seaborn palette 
    sns_height          seaborn height of figure. *default* is ``4.``. 
    sns_aspect          seaborn aspect of the figure. *default* is ``.7``
    sns_theme_kws       seaborn keywords theme arguments. default is ``{
                        'style':4., 'palette':.7}``
    verbose             control the verbosity. Higher value, more messages.
                        *default* is ``0``.
    ==================  =======================================================
    
    """
    
    @abstractmethod 
    def __init__(
        self,
        savefig: str = None,
        fig_num: int =  1,
        figsize: tuple =  (12, 8),
        fig_size : tuple =None, 
        fig_dpi:int = 300, 
        fig_legend: str =  None,
        fig_orientation: str ='landscape',
        fig_title:str = None,
        fig_aspect:str='auto',
        font_size: float =3.,
        font_style: str ='italic',
        font_weight: str = 'bold',
        fs: float = 5.,
        ms: float =3.,
        marker: str = 'o',
        markerfacecolor: str ='yellow',
        markeredgecolor: str = 'cyan',
        markeredgewidth: float =  3.,
        lc: str =  'k',
        ls: str = '-',
        lw: float = 1.,
        alpha: float =  .5,
        bins: int =  10,
        xlim: list = None, 
        ylim: list= None,
        xminorticks: int=1, 
        yminorticks: int =1,
        xlabel: str  =  None,
        ylabel: str = None,
        rotate_xlabel: int = None,
        rotate_ylabel: int =None ,
        leg_kws: dict = dict(),
        plt_kws: dict = dict(), 
        plt_style:str="pcolormesh",
        plt_shading: str="auto", 
        imshow_interp:str =None,
        s: float=  40.,
        cmap:str='jet_r',
        show_grid: bool = False,
        galpha: float = .5,
        gaxis: str = 'both',
        gc: str = 'k',
        gls: str = '--',
        glw: float = 2.,
        gwhich: str = 'major',               
        tp_axis: str = 'both',
        tp_labelsize: float = 3.,
        tp_bottom: bool =True,
        tp_top: bool = True,
        tp_labelbottom: bool = False,
        tp_labeltop: bool = True,               
        cb_orientation: str = 'vertical',
        cb_aspect: float = 20.,
        cb_shrink: float =  1.,
        cb_pad: float =.05,
        cb_anchor: tuple = (0., .5),
        cb_panchor: tuple=  (1., .5),              
        cb_label: str = None,
        cb_spacing: str = 'uniform' ,
        cb_drawedges: bool = False,
        cb_format: float = None ,   
        sns_orient: str ='v', 
        sns_style: str = None, 
        sns_palette: str= None, 
        sns_height: float=4. , 
        sns_aspect:float =.7, 
        sns_theme_kws: dict = None,
        verbose: int=0, 
        ): 
   
        self.savefig=savefig
        self.fig_num=fig_num
        self.fig_size=fig_size
        self.fig_dpi=fig_dpi
        self.fig_legend=fig_legend
        self.fig_orientation=fig_orientation
        self.fig_title=fig_title
        self.fig_aspect=fig_aspect
        self.font_size=font_size
        self.font_style=font_style
        self.font_weight=font_weight
        self.fs=fs
        self.ms=ms
        self.marker=marker
        self.marker_facecolor=markerfacecolor
        self.marker_edgecolor=markeredgecolor
        self.marker_edgewidth=markeredgewidth
        self.lc=lc
        self.ls=ls
        self.lw=lw
        self.alpha=alpha
        self.bins=bins
        self.xlim=xlim
        self.ylim=ylim
        self.x_minorticks=xminorticks
        self.y_minorticks=yminorticks
        self.xlabel=xlabel
        self.ylabel=ylabel
        self.rotate_xlabel=rotate_xlabel
        self.rotate_ylabel=rotate_ylabel
        self.leg_kws=leg_kws
        self.plt_kws=plt_kws
        self.plt_style=plt_style
        self.plt_shading=plt_shading
        self.imshow_interp=imshow_interp
        self.s=s 
        self.cmap=cmap
        self.show_grid=show_grid
        self.galpha=galpha
        self.gaxis=gaxis
        self.gc=gc
        self.gls=gls
        self.glw=glw
        self.gwhich=gwhich
        self.tp_axis=tp_axis
        self.tp_labelsize=tp_labelsize  
        self.tp_bottom=tp_bottom
        self.tp_top=tp_top
        self.tp_labelbottom=tp_labelbottom
        self.tp_labeltop=tp_labeltop
        self.cb_orientation=cb_orientation
        self.cb_aspect=cb_aspect
        self.cb_shrink=cb_shrink
        self.cb_pad=cb_pad
        self.cb_anchor=cb_anchor
        self.cb_panchor=cb_panchor
        self.cb_label=cb_label
        self.cb_spacing=cb_spacing
        self.cb_drawedges=cb_drawedges
        self.cb_format=cb_format  
        self.sns_orient=sns_orient
        self.sns_style=sns_style
        self.sns_palette=sns_palette
        self.sns_height=sns_height
        self.sns_aspect=sns_aspect
        self.verbose=verbose
        self.sns_theme_kws=sns_theme_kws or {'style':self.sns_style, 
                                         'palette':self.sns_palette, 
                                                      }
        self.cb_props = {
            pname.replace('cb_', '') : pvalues
                         for pname, pvalues in self.__dict__.items() 
                         if pname.startswith('cb_')
                         }

    _PARAMS = OrderedDict(
        (p.name, p)
        for p in signature(__init__).parameters.values()
        if p.name not in ("self", "sns_theme_kws")  # ignore internals
    )

    # human-readable banner 
    def __repr__(self) -> str:  # noqa: D401
        cls = self.__class__.__name__
        items = ", ".join(
            f"{k}={_short(getattr(self, k))}"
            for k in ("fig_num", "fig_size", "plt_style", "cmap", "verbose")
        )
        return f"<{cls} {items}>"

    __str__ = __repr__


    def update(self, **kwargs):
        """
        Override any attribute after construction.

        Only parameters recognised by ``__init__`` are accepted
        (case–insensitive).  Return ``self`` to enable chaining.

        Examples
        --------
        >>> plot = MyFancyPlot()
        >>> plot.update(cmap="viridis", lw=2).draw(...)
        """
        unk = []
        for key, value in kwargs.items():
            k = key.lower()
            hits = [p for p in self._PARAMS if p.lower() == k]
            if hits:
                setattr(self, hits[0], value)
            else:
                unk.append(key)
        if unk:
            raise AttributeError(
                "Unknown parameter(s) for update(): "
                + ", ".join(map(repr, unk))
            )
        return self

    # figure + axis factory 
    def _fig_ax(self, **ax_kws):
        """
        Return a fresh (fig, ax) with *all* global styling applied.

        Sub-classes should rely on this helper instead of calling
        ``plt.subplots`` directly – it guarantees consistency across the
        whole library.
        """
        fsize = self.fig_size or (5, 5)
        fig, ax = plt.subplots(
            num=self.fig_num,
            figsize=fsize,
            dpi=self.fig_dpi,
            subplot_kw=ax_kws or {},
        )

        # title / axis labels
        if self.fig_title:
            ax.set_title(
                self.fig_title,
                fontsize=self.fs + 4,
                fontweight=self.font_weight,
            )
        if self.xlabel:
            ax.set_xlabel(self.xlabel, rotation=self.rotate_xlabel or 0)
        if self.ylabel:
            ax.set_ylabel(self.ylabel, rotation=self.rotate_ylabel or 0)

        # limits
        if self.xlim is not None:
            ax.set_xlim(*self.xlim)
        if self.ylim is not None:
            ax.set_ylim(*self.ylim)

        # minor ticks
        ax.minorticks_on()
        ax.xaxis.set_minor_locator(plt.MultipleLocator(self.x_minorticks))
        ax.yaxis.set_minor_locator(plt.MultipleLocator(self.y_minorticks))

        # grid
        if self.show_grid:
            ax.grid(
                axis=self.gaxis,
                which=self.gwhich,
                linestyle=self.gls,
                color=self.gc,
                linewidth=self.glw,
                alpha=self.galpha,
            )

        # tick params
        ax.tick_params(
            axis=self.tp_axis,
            labelsize=self.tp_labelsize,
            bottom=self.tp_bottom,
            top=self.tp_top,
            labelbottom=self.tp_labelbottom,
            labeltop=self.tp_labeltop,
        )

        return fig, ax


class Plot2D(BasePlot):
    """
    Minimal 2-D image helper built on top of :class:`BasePlot`.

    * Accepts *any* ``BasePlot`` keyword in ``__init__``.
    * Two back-ends are supported: ``"pcolormesh"`` (default) and
      ``"imshow"`` chosen through *style*.
    * Can be used directly **or** sub-classed for specialised
      electromagnetic plots (e.g. sounding sections, resistivity maps).
    """

    def __init__(
        self,
        style: Literal["pcolormesh", "imshow"] | str = "pcolormesh",
        *,
        savefig: str | Path | None = None,
        **bpkws,
    ) -> None:
        # store the back-end then forward everything else to BasePlot
        self.style = style.lower()
        super().__init__(
            savefig=savefig, plt_style = str(style).lower(), **bpkws)

    # public API
    @overload
    def __call__(
        self,
        arr2d: np.ndarray,
        *,
        x: Iterable[float] | None = None,
        y: Iterable[float] | None = None,
        stations: Iterable[str] | None = None,
    ): ...

    def __call__(self, *args, **kwargs):
        """Alias for :meth:`plot` so ``Plot2D(...)()`` works."""
        return self.plot(*args, **kwargs)

    def plot(
        self,
        arr2d: np.ndarray,
        *,
        x: Iterable[float] | None = None,
        y: Iterable[float] | None = None,
        stations: Iterable[str] | None = None,
    ):
        """
        Render *arr2d* with optional X/Y axes or station names.

        Parameters
        ----------
        arr2d
            Two-dimensional array of shape ``(Ny, Nx)``.
        x, y
            Coordinates for the *columns* (X) and *rows* (Y).  When any
            of them is *None* the index is used.
        stations
            Optional list of *Nx* station labels replacing the top
            x-axis ticks.
        """
        arr2d = np.asarray(arr2d)
        if arr2d.ndim != 2:
            raise ValueError("arr2d must be 2-D")

        ny, nx = arr2d.shape
        x = np.arange(nx) if x is None else np.asarray(x, dtype=float)
        y = np.arange(ny) if y is None else np.asarray(y, dtype=float)

        if x.size != nx or y.size != ny:
            raise ValueError("x/y length mismatch with arr2d shape")

        #
        # figure boiler-plate inherited from BasePlot
        fig, ax = self._fig_ax()

        img = self._build_image(ax, x, y, arr2d)
        self._post_process(fig, ax, img, stations)

        if self.savefig is not None:
            fig.savefig(
                self.savefig,
                dpi=self.fig_dpi,
                orientation=self.fig_orientation,
            )
        elif self.verbose:
            print("[plot] figure ready – not saved because savefig=None")
        plt.show(block=False)
        return ax

    # helpers (meant to be overridden in sub-classes)
    def _build_image(self, ax, x, y, z):
        cmap = plt.get_cmap(self.cmap)
        if self.style.startswith("pc"):
            X, Y = np.meshgrid(x, y)
            return ax.pcolormesh(
                X,
                Y,
                z,
                cmap=cmap,
                shading=self.plt_shading,
                vmin=np.nanmin(z),
                vmax=np.nanmax(z),
            )
        elif self.style.startswith("im"):
            extent = (x.min(), x.max(), y.min(), y.max())
            return ax.imshow(
                z,
                aspect=self.fig_aspect,
                interpolation=self.imshow_interp or "nearest",
                cmap=cmap,
                origin="upper",
                extent=extent,
            )
        else:
            raise ValueError(f"unknown style {self.style!r}")

    # add colour-bar, twin-x with station labels, etc.
    def _post_process(self, fig, ax, img, stations):
        cbar = fig.colorbar(img, ax=ax, orientation=self.cb_orientation)
        cbar.set_label(
            self.cb_label or "value",
            size=self.font_size * 1.4,
            style=self.font_style,
        )

        # optional station labels on a twin axis
        if stations is not None:
            if len(stations) != img.get_array().shape[1]:
                warnings.warn("stations length does not match Nx – ignored")
                return
            ax2 = ax.twiny()
            ax2.set_xlim(ax.get_xlim())
            ax2.set_xticks(np.linspace(ax2.get_xlim()[0], ax2.get_xlim(
                )[1], len(stations)))
            ax2.set_xticklabels(
                stations,
                rotation=self.rotate_xlabel or 0,
                fontsize=self.font_size,
            )
            ax2.set_xlabel(
                "Stations",
                fontdict=dict(
                    size=self.font_size * 1.4,
                    style=self.font_style,
                    weight=self.font_weight,
                ),
            )
            

def _short(val, n=40):
    """Helper to truncate long values for display."""
    s = str(val)
    return s if len(s) <= n else s[: n - 3] + "..."