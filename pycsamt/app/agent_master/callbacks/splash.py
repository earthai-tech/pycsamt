# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Splash overlay dismiss callback."""

from __future__ import annotations

from dash import Input, Output
from dash.exceptions import PreventUpdate

from .._ids import IDs

_SPLASH_CARDS = [
    IDs.SPLASH_CARD_LOAD,
    IDs.SPLASH_CARD_CHAT,
    IDs.SPLASH_CARD_AI,
    IDs.SPLASH_CARD_REPORT,
]


def register_splash(app) -> None:
    # Any card click also opens the EDI canvas
    @app.callback(
        Output(IDs.CANVAS_EDI, "is_open",
               allow_duplicate=True),
        Input(IDs.SPLASH_CARD_LOAD, "n_clicks"),
        prevent_initial_call=True,
    )
    def _card_load_opens_edi(n):
        if not n:
            raise PreventUpdate
        return True

    # CTA + all cards dismiss the overlay
    # via clientside callback (no round-trip)
    app.clientside_callback(
        """
        function() {
          var el = document.getElementById(
            'am-splash-overlay'
          );
          if (!el) return '';
          el.classList.add('wlc-hiding');
          setTimeout(function() {
            el.classList.add('wlc-gone');
          }, 800);
          return '';
        }
        """,
        Output(IDs.SPLASH_OVERLAY, "data-hidden"),
        Input(IDs.SPLASH_CTA, "n_clicks"),
        Input(IDs.SPLASH_CARD_LOAD, "n_clicks"),
        Input(IDs.SPLASH_CARD_CHAT, "n_clicks"),
        Input(IDs.SPLASH_CARD_AI, "n_clicks"),
        Input(
            IDs.SPLASH_CARD_REPORT, "n_clicks"
        ),
        prevent_initial_call=True,
    )
