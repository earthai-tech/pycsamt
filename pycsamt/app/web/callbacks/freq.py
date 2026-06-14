# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
callbacks/freq.py — Frequency / period range slider initialisation.

When EDI data loads (STORE_DATA changes), derive the full period range from the
cached Sites object and update the RangeSlider bounds and value.
"""

from __future__ import annotations

import math

from dash import Input, Output, State, no_update

from pycsamt.app.web.cache import cache_get
from pycsamt.app.web.layout import IDs


def _log10_range(sites):
    """Return (log10_T_min, log10_T_max) from the loaded Sites object."""
    try:
        import numpy as np
        from pycsamt.emtools._core import ensure_sites

        sites = ensure_sites(sites)
        freqs = []
        for s in sites:
            try:
                f = getattr(s, "freq", None) or getattr(s, "freqs", None)
                if f is not None and len(f):
                    freqs.extend(float(x) for x in f if x > 0)
            except Exception:
                continue

        if not freqs:
            return -4.0, 4.0

        arr = np.asarray(freqs)
        T_vals = 1.0 / arr[arr > 0]
        log_T  = np.log10(T_vals)
        lo = math.floor(log_T.min() * 4) / 4   # round to nearest 0.25
        hi = math.ceil(log_T.max() * 4)  / 4
        return lo, hi

    except Exception:
        return -4.0, 4.0


def register_freq(app) -> None:
    @app.callback(
        Output(IDs.FREQ_SLIDER, "min"),
        Output(IDs.FREQ_SLIDER, "max"),
        Output(IDs.FREQ_SLIDER, "value"),
        Output(IDs.FREQ_SLIDER, "disabled"),
        Input(IDs.STORE_DATA,  "data"),
        State(IDs.SESSION_ID,  "data"),
        prevent_initial_call=True,
    )
    def init_slider(store_data, session_id):
        if not store_data or not session_id:
            return -4.0, 4.0, [-4.0, 4.0], True

        sites = cache_get(session_id)
        if sites is None:
            return -4.0, 4.0, [-4.0, 4.0], True

        lo, hi = _log10_range(sites)
        return lo, hi, [lo, hi], False
