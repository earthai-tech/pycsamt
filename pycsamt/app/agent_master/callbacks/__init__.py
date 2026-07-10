# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Register all agent_master callbacks."""

from __future__ import annotations

from .chat import register_chat
from .edi import register_edi
from .help import register_help
from .inv_params import register_inv_params
from .line_sel import register_line_sel
from .outdir import register_outdir
from .params import register_params
from .plus import register_plus
from .postproc import register_postproc
from .preview import register_preview
from .session import register_session
from .settings import register_settings
from .sidebar import register_sidebar
from .splash import register_splash


def register_all(app) -> None:
    register_edi(app)
    register_settings(app)
    register_session(app)
    register_inv_params(app)
    register_chat(app)
    register_params(app)
    register_postproc(app)
    register_line_sel(app)
    register_outdir(app)
    register_preview(app)
    register_splash(app)
    register_plus(app)
    register_sidebar(app)
    register_help(app)
