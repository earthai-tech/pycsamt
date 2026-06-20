# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
TDEMController — loads TDEM survey data and renders plots into canvas figures.

Data model
──────────
A single TEMSurvey is the central object: it carries all AVG/Z/LOG files
from a survey folder.  Soundings (list[TEMSounding]) are loaded alongside
it for decay-curve and apparent-resistivity plots.

Each draw_*(fig) method:
  1. Instantiates the appropriate pycsamt.tdem plot class with the right data.
  2. Calls  .plot(ax=ax)  for single-axes plots, or  .plot()  for
     multi-axes / figure-level plots (which returns its own Figure whose
     axes are then transplanted into the canvas Figure).
  3. Applies dark/light theme styling across all axes.
"""
from __future__ import annotations

from typing import List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np

from typing import Optional

from pycsamt.app.desktop.controllers.qc_controller import (
    _style_all_axes,
    _annotate_empty,
)


# ── Plot catalogue ─────────────────────────────────────────────────────────────
# Each entry: (display_label, class_name, has_ax, data_key)
#   data_key:  "survey"   → TEMSurvey
#              "soundings"→ list[TEMSounding]
#              "dashboard"→ special multi-arg call

DECAY_PLOTS: list = [
    ("Decay curves",           "PlotDecayCurve",       True,  "soundings"),
    # PlotTransformedRho.plot() needs 2 axes when show_phase=True; let it
    # create its own figure (has_ax=False) so _call_figure_cls handles it.
    ("Apparent resistivity",   "PlotTransformedRho",   False, "soundings"),
]

SECTION_PLOTS: list = [
    ("TEMAVG section",         "PlotTEMAVGSection",    True,  "survey"),
    ("Z-contour section",      "PlotTEMZSection",      True,  "survey"),
    ("Gate profile",           "PlotGateProfile",      True,  "survey"),
]

MAP_PLOTS: list = [
    ("Station map",            "PlotSurveyMap",        True,  "survey"),
    ("Survey overview",        "PlotSurveyOverview",   False, "survey"),
    ("Elevation profile",      "PlotElevationProfile", True,  "survey"),
]

DASHBOARD_PLOTS: list = [
    ("Full TEM dashboard",     "PlotTEMDashboard",     False, "dashboard"),
]

TDEM_GROUPS: list = [
    ("Decay / Rho",     DECAY_PLOTS),
    ("Survey Section",  SECTION_PLOTS),
    ("Map & Overview",  MAP_PLOTS),
    ("Dashboard",       DASHBOARD_PLOTS),
]


# ── Controller ────────────────────────────────────────────────────────────────

class TDEMController:
    """
    Manages loaded TDEM data and drives all TDEM canvas renders.

    Attributes
    ----------
    dark : bool
        Apply dark-mode figure styling (default True).
    """

    def __init__(self) -> None:
        self._survey   = None          # TEMSurvey
        self._soundings: List = []     # list[TEMSounding]
        self._folder:  str  = ""
        self.dark: bool = True

        # Summary info shown in the params panel
        self.summary: str = ""

    # ── Data loading ──────────────────────────────────────────────────────────

    def load_folder(self, folder: str, progress_cb=None) -> bool:
        """
        Load all TDEM files from *folder* via ``read_temavg_survey``.

        Parameters
        ----------
        folder : str
            Path to a folder containing .AVG / .Z / .LOG files.
        progress_cb : callable(int), optional
            Called with 0–100 during loading.

        Returns
        -------
        bool
            True on success, False if no TDEM files were found.
        """
        import pycsamt.tdem as tdem

        if progress_cb:
            progress_cb(10)

        try:
            self._survey = tdem.read_temavg_survey(folder)
        except Exception as exc:
            self.summary = f"Load error: {exc}"
            return False

        if progress_cb:
            progress_cb(50)

        try:
            self._soundings = tdem.read_temavg_soundings(folder)
        except Exception:
            self._soundings = []

        if progress_cb:
            progress_cb(100)

        self._folder = folder
        self._build_summary()
        return bool(self._survey)

    def is_loaded(self) -> bool:
        return self._survey is not None

    # ── Drawing ───────────────────────────────────────────────────────────────

    def draw(
        self,
        class_name: str,
        has_ax: bool,
        data_key: str,
        fig,
    ) -> Optional[object]:
        """
        Render a TDEM plot.

        Parameters
        ----------
        class_name : str
            Name of the pycsamt.tdem plot class (e.g. ``"PlotDecayCurve"``).
        has_ax : bool
            True  → draw into *fig* in-place; return None.
            False → cls.plot() creates its own figure; return it so the
                    caller can hand it to MplCanvas.show_figure().
        data_key : str
            Which data to pass: ``"survey"``, ``"soundings"``, or
            ``"dashboard"`` (special multi-arg case).
        fig : matplotlib.figure.Figure
            Target figure (mutated when has_ax=True or on error).

        Returns
        -------
        Figure or None
        """
        import pycsamt.tdem as tdem

        fig.clear()

        if not self.is_loaded():
            ax = fig.add_subplot(111)
            _annotate_empty(ax, "Load TDEM data first\n(click Browse…)")
            _style_all_axes(fig, self.dark)
            return None

        cls = getattr(tdem, class_name, None)
        if cls is None:
            ax = fig.add_subplot(111)
            _annotate_empty(ax, f"Class not found:\n{class_name}")
            _style_all_axes(fig, self.dark)
            return None

        try:
            if data_key == "soundings":
                data = self._soundings
            elif data_key == "dashboard":
                data = None
            else:
                data = self._survey

            if has_ax:
                ax = fig.add_subplot(111)
                if data_key == "dashboard":
                    # Dashboard always uses has_ax=False in the catalogue, but
                    # guard here defensively.
                    plot_obj = cls(self._survey, self._survey, self._soundings)
                else:
                    plot_obj = cls(data)
                import inspect
                sig = inspect.signature(plot_obj.plot)
                if "axes" in sig.parameters:
                    plot_obj.plot(axes=ax)
                else:
                    plot_obj.plot(ax=ax)
            else:
                src_fig = self._call_figure_cls(cls, data)
                if src_fig is None:
                    ax = fig.add_subplot(111)
                    _annotate_empty(ax, "No figure produced")
                else:
                    _style_all_axes(src_fig, self.dark)
                    try:
                        src_fig.tight_layout(pad=1.2)
                    except Exception:
                        pass
                    return src_fig

        except Exception as exc:
            fig.clear()
            ax = fig.add_subplot(111)
            _annotate_empty(ax, f"{class_name} error:\n{exc}")

        _style_all_axes(fig, self.dark)
        try:
            fig.tight_layout(pad=1.2)
        except Exception:
            pass
        return None

    # ── Helpers ───────────────────────────────────────────────────────────────

    def _call_figure_cls(self, cls, data) -> Optional[object]:
        """Instantiate cls(data).plot() and return the Figure it creates.

        For the dashboard, `data` is None and both avg/zplot come from the
        survey (TEMSurvey supports both _avg_records and _z_records).
        """
        before = set(plt.get_fignums())
        try:
            if data is not None:
                result = cls(data).plot()
            else:
                # PlotTEMDashboard(avg, zplot, soundings)
                # Pass survey for both avg and zplot — both accept TEMSurvey.
                result = cls(self._survey, self._survey, self._soundings).plot()
        except Exception:
            return None
        after = set(plt.get_fignums())
        if hasattr(result, "get_axes"):
            return result
        new_nums = after - before
        if new_nums:
            return plt.figure(max(new_nums))
        return None

    def _build_summary(self) -> None:
        if self._survey is None:
            self.summary = "No data loaded"
            return
        n_avg = getattr(self._survey, "n_avg_files", 0)
        n_z   = getattr(self._survey, "n_z_files",  0)
        n_snd = len(self._soundings)
        self.summary = (
            f"Folder: {self._folder}\n"
            f"AVG files: {n_avg}   Z files: {n_z}\n"
            f"Soundings: {n_snd}"
        )
