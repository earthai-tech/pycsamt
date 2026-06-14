# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
web/pages — one module per navigation section.

Each module exposes:
  PAGE_ID : str             — unique page identifier (matches nav-section store)
  layout() -> html.Div     — page content (pre-rendered at app startup)
  register_callbacks(app)  — registers this page's Dash callbacks
"""
from pycsamt.app.web.pages import (
    home, qc_page, correction, advanced, tdem,
    pipeline, forward, inversion, interpretation, agents_page,
)

ALL_PAGES = [
    home, qc_page, correction, advanced, tdem,
    pipeline, forward, inversion, interpretation, agents_page,
]
