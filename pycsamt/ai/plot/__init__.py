# -*- coding: utf-8 -*-
"""
pycsamt.ai.plot
===============

Unified, publication-quality visualisation for AI-based EM results.
"""
from ._style import (
    EMStyle, EM_COLORS, EM_CMAPS, EM_FIGSIZE, em_context, add_colorbar,
    StationTickConfig,
)
from .compare import plot_compare, plot_profile_pair
from .convergence import plot_convergence, plot_lr_schedule
from .section import plot_section, plot_section_pair, plot_pseudo_section
from .diagnostics import (
    plot_confusion_matrix,
    plot_residuals,
    plot_layer_errors,
    plot_uncertainty_bands,
    plot_feature_importance,
)

__all__ = [
    "EMStyle", "EM_COLORS", "EM_CMAPS", "EM_FIGSIZE",
    "em_context", "add_colorbar", "StationTickConfig",
    "plot_compare", "plot_profile_pair",
    "plot_convergence", "plot_lr_schedule",
    "plot_section", "plot_section_pair", "plot_pseudo_section",
    "plot_confusion_matrix", "plot_residuals", "plot_layer_errors",
    "plot_uncertainty_bands", "plot_feature_importance",
]
