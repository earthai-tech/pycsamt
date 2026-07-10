# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Callback registration for the map-view platform."""

from __future__ import annotations


def register_all(app) -> None:
    """Register every callback module on *app*."""
    from .chrome import register_chrome
    from .controls import register_controls
    from .export import register_export
    from .inversion_import import register_inversion_import
    from .lines import register_lines
    from .load import register_load
    from .session import register_session
    from .settings import register_settings
    from .toolbar import register_toolbar
    from .topo import register_topo
    from .view import register_view

    register_chrome(app)
    register_load(app)
    register_inversion_import(app)
    register_view(app)
    register_controls(app)
    register_lines(app)
    register_export(app)
    register_topo(app)
    register_toolbar(app)
    register_settings(app)
    register_session(app)
