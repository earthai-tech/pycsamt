# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Session initialisation callback — fires once per browser tab on first load."""

from __future__ import annotations

import uuid

from dash import Input, Output, no_update

from pycsamt.app.web.layout import IDs


def register_session(app) -> None:
    @app.callback(
        Output(IDs.SESSION_ID, "data"),
        Input(IDs.SESSION_ID,  "data"),
        prevent_initial_call=False,
    )
    def init_session_id(existing_id):
        if existing_id:
            return no_update
        return str(uuid.uuid4())
