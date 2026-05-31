# -*- coding: utf-8 -*-
"""
pycsamt.ai.plot
===============

Unified, publication-quality visualisation for AI-based EM results.
"""
from ._style import EMStyle, EM_COLORS, EM_CMAPS, EM_FIGSIZE, em_context, add_colorbar
from .compare import plot_compare, plot_profile_pair
from .convergence import plot_convergence, plot_lr_schedule

__all__ = [
    "EMStyle", "EM_COLORS", "EM_CMAPS", "EM_FIGSIZE",
    "em_context", "add_colorbar",
    "plot_compare", "plot_profile_pair",
    "plot_convergence", "plot_lr_schedule",
]
