# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Callback registration for the map-view platform."""

from __future__ import annotations


def register_all(app) -> None:
    """Register every callback module on *app*."""
    from .chrome import register_chrome
    from .load import register_load
    from .view import register_view
    from .controls import register_controls
    from .lines import register_lines
    from .export import register_export
    from .topo import register_topo
    from .toolbar import register_toolbar

    register_chrome(app)
    register_load(app)
    register_view(app)
    register_controls(app)
    register_lines(app)
    register_export(app)
    register_topo(app)
    register_toolbar(app)
