# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.app.web.callbacks — Dash callback registry.

Each submodule owns one logical slice of the UI:

    _session  → per-tab UUID initialisation
    data      → modal toggle, EDI loading (upload + path), station table
    map       → station map figure, map↔table selection sync
    profile   → MT-curve / pseudosection profile tab renders
    agents    → AI agent execution panel
    theme     → dark / light theme toggle
    freq      → period range slider initialisation
    qc        → QC group/plot dropdowns + background QC render

All submodules export a single ``register_<name>(app)`` function.
"""

from __future__ import annotations

from pycsamt.app.web.callbacks._session    import register_session
from pycsamt.app.web.callbacks.data        import register_data
from pycsamt.app.web.callbacks.map         import register_map
from pycsamt.app.web.callbacks.profile     import register_profile
from pycsamt.app.web.callbacks.agents      import register_agents
from pycsamt.app.web.callbacks.theme       import register_theme
from pycsamt.app.web.callbacks.freq        import register_freq
from pycsamt.app.web.callbacks.qc          import register_qc
from pycsamt.app.web.callbacks.navigation  import register_navigation
from pycsamt.app.web.callbacks.advanced    import register_advanced
from pycsamt.app.web.callbacks.tdem        import register_tdem
from pycsamt.app.web.callbacks.correction  import register_correction
from pycsamt.app.web.callbacks.pipeline    import register_pipeline
from pycsamt.app.web.callbacks.forward     import register_forward
from pycsamt.app.web.callbacks.interp      import register_interp
from pycsamt.app.web.callbacks.inversion   import register_inversion
from pycsamt.app.web.callbacks.tools       import register_tools
from pycsamt.app.web.callbacks.settings    import register_settings
from pycsamt.app.web.callbacks.help        import register_help
from pycsamt.app.web.callbacks.dashboard   import register_dashboard


def register_callbacks(app) -> None:
    """Register every Dash callback on *app* in dependency order."""
    register_session(app)
    register_data(app)
    register_map(app)
    register_profile(app)
    register_agents(app)
    register_theme(app)
    register_freq(app)
    register_qc(app)
    register_navigation(app)
    register_advanced(app)
    register_tdem(app)
    register_correction(app)
    register_pipeline(app)
    register_forward(app)
    register_interp(app)
    register_inversion(app)
    register_tools(app)
    register_settings(app)
    register_help(app)
    register_dashboard(app)
    # Page-level card/selection callbacks
    from pycsamt.app.web.pages.agents_page import register_callbacks as _reg_agents_page
    _reg_agents_page(app)
