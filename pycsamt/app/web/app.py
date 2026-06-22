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

import argparse
import os
import socket
import sys
import threading
import webbrowser
from typing import Optional, Sequence

# Absolute path to the shared SVG icon set
_ICONS_DIR = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),   # .../app/web/
    "..", "desktop", "resources", "icons",
)
_ICONS_DIR = os.path.normpath(_ICONS_DIR)


def _find_free_port(host: str, preferred: int) -> int:
    """Return *preferred* when available, otherwise ask the OS for a free port."""
    if int(preferred) == 0:
        return _allocate_free_port(host)
    try:
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
            sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
            sock.bind((host, int(preferred)))
    except OSError:
        return _allocate_free_port(host)
    return int(preferred)


def _allocate_free_port(host: str) -> int:
    """Ask the OS for an available TCP port."""
    bind_host = "127.0.0.1" if host == "0.0.0.0" else host
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
        sock.bind((bind_host, 0))
        return int(sock.getsockname()[1])


def _browser_url(host: str, port: int) -> str:
    """Return the user-facing URL for the local web app."""
    display_host = "127.0.0.1" if host in {"0.0.0.0", "::"} else host
    return f"http://{display_host}:{int(port)}"


def _open_browser_later(url: str, delay: float) -> None:
    """Open *url* in the default browser after a short delay."""
    timer = threading.Timer(max(0.0, float(delay)), lambda: webbrowser.open(url))
    timer.daemon = True
    timer.start()


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="pycsamt-web",
        description="Launch the pyCSAMT Dash web application.",
    )
    parser.add_argument(
        "--host",
        default="127.0.0.1",
        help=(
            "Bind address. Use 127.0.0.1 for local-only access or 0.0.0.0 "
            "to allow access from other machines on the network."
        ),
    )
    parser.add_argument(
        "--port",
        type=int,
        default=8050,
        help="Preferred TCP port. Use 0 to always choose a free port.",
    )
    parser.add_argument(
        "--strict-port",
        action="store_true",
        help="Fail if --port is unavailable instead of choosing another port.",
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Enable Dash debug/dev tools.",
    )
    browser = parser.add_mutually_exclusive_group()
    browser.add_argument(
        "--browser",
        dest="open_browser",
        action="store_true",
        default=True,
        help="Open the web app in the default browser. This is the default.",
    )
    browser.add_argument(
        "--no-browser",
        dest="open_browser",
        action="store_false",
        help="Start the server without opening a browser.",
    )
    parser.add_argument(
        "--browser-delay",
        type=float,
        default=1.0,
        help="Seconds to wait before opening the browser.",
    )
    return parser


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
    from pycsamt.app.web.cache import bg_manager

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

    # Favicon is served from assets/favicon.ico — Dash detects it automatically
    # and injects <link rel="icon" href="/assets/favicon.ico"> via {%favicon%}.
    # Using assets/ avoids the browser's aggressive per-origin favicon cache
    # that causes /favicon.ico to stay stale after a Flask-route change.

    app.layout = layout()
    register_callbacks(app)

    return app


def main(
    host: str = "127.0.0.1",
    port: int = 8050,
    debug: bool = False,
    open_browser: bool = True,
    strict_port: bool = False,
    browser_delay: float = 1.0,
    argv: Optional[Sequence[str]] = None,
) -> None:
    """
    Entry point: ``python -m pycsamt.app.web`` or ``pycsamt-web``.

    Parameters
    ----------
    host : str
        Bind address. ``127.0.0.1`` is local-only; ``0.0.0.0`` exposes
        the server on the local network.
    port : int
        Preferred port number. If unavailable, a free port is selected
        unless *strict_port* is True.
    debug : bool
        Enable Dash debug / hot-reload mode.
    open_browser : bool
        Open the app in the default browser after startup.
    strict_port : bool
        Raise immediately if the requested port is unavailable.
    browser_delay : float
        Delay before opening the browser.
    argv : sequence of str, optional
        Command-line arguments. When provided, overrides the keyword values.
    """
    if argv is not None:
        args = _build_parser().parse_args(list(argv))
        host = args.host
        port = args.port
        debug = args.debug
        open_browser = args.open_browser
        strict_port = args.strict_port
        browser_delay = args.browser_delay

    requested_port = int(port)
    run_port = requested_port
    if not strict_port:
        run_port = _find_free_port(host, requested_port)

    url = _browser_url(host, run_port)
    app = create_app(debug=debug)
    if open_browser:
        _open_browser_later(url, browser_delay)

    if run_port != requested_port:
        print(
            f"\n  pycsamt web app → {url}"
            f"\n  requested port {requested_port} was busy; using {run_port}\n"
        )
    else:
        print(f"\n  pycsamt web app → {url}\n")

    app.run(host=host, port=run_port, debug=debug, use_reloader=False)


if __name__ == "__main__":
    main(argv=sys.argv[1:])
