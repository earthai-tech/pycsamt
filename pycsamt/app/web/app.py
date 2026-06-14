# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Dash web application for pycsamt (Phase 6).

Mirrors the desktop application's three-zone layout using Dash Bootstrap
Components (DBC 2.x) and Plotly figures.  All scientific logic is delegated
to pycsamt.app.desktop.controllers — only the view layer changes.

Layout zones
------------
Left sidebar  — station DataTable + summary label
Centre top    — Plotly interactive station map (Scattergeo)
Centre bottom — matplotlib profile tabs (ρₐ/φ, pseudosections, tipper, PT)
Right         — AI Agents: picker + log + result image
Footer        — status bar

Entry point::

    python -m pycsamt.app.web
    # → http://localhost:8050
"""

from __future__ import annotations

import os
import sys

# Absolute path to the shared SVG icon set
_ICONS_DIR = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),   # .../app/web/
    "..", "desktop", "resources", "icons",
)
_ICONS_DIR = os.path.normpath(_ICONS_DIR)


def create_app(debug: bool = False):
    """
    Create and return the configured Dash application.

    Parameters
    ----------
    debug : bool
        Enable Dash debug mode (hot-reload, Dev Tools panel).

    Returns
    -------
    dash.Dash
        Fully configured application with layout and callbacks registered.
    """
    try:
        import dash
        import dash_bootstrap_components as dbc
    except ImportError as exc:
        raise ImportError(
            "Dash and dash-bootstrap-components are required for the web app.\n"
            "Install them with:  pip install 'pycsamt[app]'"
        ) from exc

    from pycsamt.app.web.layout import layout
    from pycsamt.app.web.callbacks import register_callbacks
    from pycsamt.app.web.cache import cache, has_diskcache

    # Build the background-callback manager when diskcache is available.
    # Without it, background=True callbacks fall back to synchronous execution.
    bg_manager = None
    if has_diskcache():
        try:
            from dash import DiskcacheManager
            bg_manager = DiskcacheManager(cache)
        except Exception:
            pass

    dash_kwargs: dict = dict(
        external_stylesheets=[
            dbc.themes.CYBORG,
            dbc.icons.BOOTSTRAP,
        ],
        suppress_callback_exceptions=True,
        title="pycsamt",
        update_title="pycsamt — loading…",
        meta_tags=[
            {"name": "viewport", "content": "width=device-width, initial-scale=1"},
        ],
    )
    if bg_manager is not None:
        dash_kwargs["background_callback_manager"] = bg_manager

    app = dash.Dash(__name__, **dash_kwargs)

    # ── Serve the shared icon set as /icons/<filename> ───────────────────────
    from flask import send_from_directory

    @app.server.route("/icons/<path:filename>")
    def serve_icon(filename):
        return send_from_directory(_ICONS_DIR, filename)

    # ── Favicon ───────────────────────────────────────────────────────────────
    @app.server.route("/favicon.ico")
    def favicon():
        return send_from_directory(
            _ICONS_DIR, "pycsamt.ico", mimetype="image/vnd.microsoft.icon"
        )

    app.layout = layout()
    register_callbacks(app)

    return app


def main(
    host: str  = "0.0.0.0",
    port: int  = 8050,
    debug: bool = False,
) -> None:
    """
    Entry point: ``python -m pycsamt.app.web``.

    Parameters
    ----------
    host : str
        Bind address (default ``0.0.0.0`` — accessible on local network).
    port : int
        Port number (default 8050).
    debug : bool
        Enable Dash debug / hot-reload mode.
    """
    app = create_app(debug=debug)
    print(f"\n  pycsamt web app → http://localhost:{port}\n")
    app.run(host=host, port=port, debug=debug)


if __name__ == "__main__":
    main()
