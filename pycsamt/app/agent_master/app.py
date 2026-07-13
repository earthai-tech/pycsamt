# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""
pyCSAMT Agent Master — Dash app factory.

Usage::

    from pycsamt.app.agent_master import launch
    launch()               # http://localhost:8765

    # or
    from pycsamt.app.agent_master.app import (
        create_app
    )
    server = create_app().server  # WSGI
"""

from __future__ import annotations

from pathlib import Path

# Force non-interactive backend before any
# other matplotlib import can claim a display.
import matplotlib

matplotlib.use("Agg", force=True)

_ICONS_DIR = Path(__file__).parent / "assets" / "icons"


def create_app(
    debug: bool = False,
    suppress_exceptions: bool = True,
):
    """
    Create and return the configured Dash app.

    Parameters
    ----------
    debug : bool
        Enable Dash debug mode.
    suppress_exceptions : bool
        Suppress callback exceptions
        (needed for pattern-matching IDs).

    Returns
    -------
    dash.Dash
    """
    import dash
    import dash_bootstrap_components as dbc

    # Bootstrap + Bootstrap Icons
    # + Highlight.js (atom-one-dark theme)
    _HJS = "11.9.0"
    _HJS_CDN = f"https://cdnjs.cloudflare.com/ajax/libs/highlight.js/{_HJS}"
    _EXT = [
        dbc.themes.BOOTSTRAP,
        dbc.icons.BOOTSTRAP,
        f"{_HJS_CDN}/styles/atom-one-dark.min.css",
    ]
    _EXT_JS = [
        {
            "src": (f"{_HJS_CDN}/highlight.min.js"),
        },
    ]

    app = dash.Dash(
        __name__,
        external_stylesheets=_EXT,
        external_scripts=_EXT_JS,
        suppress_callback_exceptions=(suppress_exceptions),
        title="pyCSAMT — Agent",
        update_title=None,
    )

    from .layout import create_layout

    app.layout = create_layout()

    from .callbacks import register_all

    register_all(app)

    # serve icons from assets/icons/
    from flask import send_from_directory

    @app.server.route("/am-icons/<path:filename>")
    def _serve_icon(filename):
        return send_from_directory(_ICONS_DIR, filename)

    # clientside: apply data-theme to <html>
    # and className for CSS scoping
    app.clientside_callback(
        """
        function(theme) {
          var t = theme || 'light';
          document.documentElement
            .setAttribute('data-theme', t);
          return 'am-root am-root-' + t;
        }
        """,
        dash.Output("am-root", "className"),
        dash.Input("am-store-theme", "data"),
        prevent_initial_call=False,
    )

    # clientside: Enter key + animated placeholder
    app.clientside_callback(
        """
        function(n) {
          var el = document.getElementById(
            'am-input'
          );
          if (!el) return '';

          // Enter key → send
          el.addEventListener('keydown',
            function(e) {
              if (e.key === 'Enter'
                  && !e.shiftKey) {
                e.preventDefault();
                var btn = document.getElementById(
                  'am-btn-send'
                );
                if (btn) btn.click();
              }
          });

          // Animated placeholder cycle
          var hints = [
            'Load EDIs, run QC and static-shift'
            + ' correction...',
            'Perform AI-assisted 1-D inversion'
            + ' on all lines...',
            'Run phase tensor dimensionality'
            + ' analysis...',
            'Prepare and run Occam2D'
            + ' inversion files...',
            'Evaluate inversion result'
            + ' and generate report...',
            'Detect static shift and'
            + ' apply correction...',
          ];
          var i = 0;
          function _fade(dir, cb) {
            el.style.transition =
              'opacity 0.35s ease';
            el.style.opacity =
              dir === 'out' ? '0' : '1';
            setTimeout(cb, 370);
          }
          function _next() {
            _fade('out', function() {
              i = (i + 1) % hints.length;
              el.placeholder = hints[i];
              _fade('in', function() {});
            });
          }
          el.placeholder = hints[0];
          var _tid = setInterval(_next, 3400);
          el.addEventListener('focus',
            function() {
              clearInterval(_tid);
              el.style.opacity = '1';
              el.placeholder = hints[i];
            }
          );
          el.addEventListener('blur',
            function() {
              if (!el.value) {
                _tid = setInterval(_next, 3400);
              }
            }
          );

          return '';
        }
        """,
        dash.Output(
            "am-root",
            "data-init",
        ),
        dash.Input("am-root", "id"),
        prevent_initial_call=False,
    )

    return app


def launch(
    host: str = "127.0.0.1",
    port: int = 8765,
    debug: bool = False,
    open_browser: bool = True,
) -> None:
    """
    Launch the Agent Master in a local browser.

    Parameters
    ----------
    host : str
        Bind address.
    port : int
        HTTP port (default 8765).
    debug : bool
        Enable Dash hot-reload.
    open_browser : bool
        Open a browser tab automatically.
    """
    app = create_app(debug=debug)

    if open_browser:
        import threading
        import webbrowser

        def _open():
            import time

            time.sleep(1.2)
            webbrowser.open(f"http://{host}:{port}")

        threading.Thread(target=_open, daemon=True).start()

    app.run(
        host=host,
        port=port,
        debug=debug,
        use_reloader=debug,
    )
