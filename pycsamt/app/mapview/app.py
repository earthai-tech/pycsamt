# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""
pyCSAMT Map View — Dash app factory.

Usage::

    from pycsamt.app.mapview import launch
    launch()                       # http://localhost:8770

    # or, from a script with data already in hand:
    from pycsamt.map import MapView
    MapView.from_folder("data/survey").launch()

    # WSGI server object:
    from pycsamt.app.mapview.app import create_app
    server = create_app().server
"""

from __future__ import annotations

import os

# Force a non-interactive Matplotlib backend before anything can claim
# a display (the static export + contour overlays use it).
import matplotlib

matplotlib.use("Agg", force=True)

# Shared pyCSAMT icon set (logos + view glyphs), reused from the desktop app.
_ICONS_DIR = os.path.normpath(
    os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "..",
        "desktop",
        "resources",
        "icons",
    )
)


def create_app(
    debug: bool = False,
    suppress_exceptions: bool = True,
):
    """Create and return the configured Dash app."""
    import dash
    import dash_bootstrap_components as dbc

    app = dash.Dash(
        __name__,
        external_stylesheets=[
            dbc.themes.BOOTSTRAP,
            dbc.icons.BOOTSTRAP,
        ],
        suppress_callback_exceptions=suppress_exceptions,
        title="pyCSAMT — Map View",
        update_title=None,
    )

    from .layout import create_layout

    app.layout = create_layout()

    from .callbacks import register_all

    register_all(app)

    # serve the shared icon set at /mv-icons/<filename>
    from flask import send_from_directory

    @app.server.route("/mv-icons/<path:filename>")
    def _serve_icon(filename):
        return send_from_directory(_ICONS_DIR, filename)

    # clientside: reflect theme onto <html data-theme> for CSS scoping
    app.clientside_callback(
        """
        function(theme) {
          var t = theme || 'light';
          document.documentElement.setAttribute('data-theme', t);
          return 'mv-root mv-root-' + t;
        }
        """,
        dash.Output("mv-root", "className"),
        dash.Input("mv-store-theme", "data"),
        prevent_initial_call=False,
    )

    return app


def launch(
    host: str = "127.0.0.1",
    port: int = 8770,
    debug: bool = False,
    open_browser: bool = True,
    view: object | None = None,
) -> None:
    """Launch the Map View platform in a local browser.

    Parameters
    ----------
    host, port :
        Bind address and HTTP port (default 8770).
    debug :
        Enable Dash hot-reload.
    open_browser :
        Open a browser tab automatically.
    view :
        Optional :class:`~pycsamt.map.MapView` to seed the first
        session with (used by ``MapView.launch``).
    """
    if view is not None:
        from .cache import set_seed

        set_seed(view)

    app = create_app(debug=debug)

    if open_browser:
        import threading
        import time
        import webbrowser

        def _open():
            time.sleep(1.2)
            webbrowser.open(f"http://{host}:{port}")

        threading.Thread(target=_open, daemon=True).start()

    app.run(
        host=host,
        port=port,
        debug=debug,
        use_reloader=debug,
    )
