# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
PlotlyView — embeds a Plotly figure in a Qt widget via QtWebEngine.

pycsamt's richer map/volume renderers (``pycsamt.map.MapView``,
``pycsamt.map.volume.build_3d_map``, station/profile/pseudosection maps,
...) are backend-neutral: they build a plain Plotly ``Figure`` with no Dash
dependency. Map View (the Dash app) just happens to be the only pycsamt
front-end that has rendered them so far, in a browser. This widget lets the
desktop app reuse those exact figures too, instead of reimplementing their
geometry (fence/block/depth-slice meshes, borehole tubes, pattern-textured
geology, ...) in matplotlib.

Renders by writing the figure to a temporary local HTML file and loading it
via a ``file://`` URL rather than ``QWebEngineView.setHtml()`` -- Qt caps
``setHtml()`` payloads at ~2 MB (percent-encoded into a data: URL
internally), and an inline Plotly bundle alone is already close to that
limit before the figure's own data is added.
"""

from __future__ import annotations

import tempfile
import uuid
from pathlib import Path

from PySide6.QtCore import QTimer, QUrl
from PySide6.QtWidgets import QVBoxLayout, QWidget


def _quiet_page(parent):
    """A web page whose JavaScript console goes to the ``pycsamt`` log
    (debug; errors as warnings) instead of the terminal."""
    import logging

    from PySide6.QtWebEngineCore import QWebEnginePage

    log = logging.getLogger("pycsamt.app.desktop.plotly")

    class _QuietPage(QWebEnginePage):
        def javaScriptConsoleMessage(self, level, message, line, source):  # noqa: N802
            error = level == QWebEnginePage.JavaScriptConsoleMessageLevel.ErrorMessageLevel
            (log.warning if error else log.debug)(
                "js %s:%s %s", source, line, message)

    return _QuietPage(parent)


class PlotlyView(QWidget):
    """A QWidget that renders one Plotly figure at a time."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._figure = None
        self._tmp_dir = Path(tempfile.mkdtemp(prefix="pycsamt_plotly_"))
        self._html_path = self._tmp_dir / f"{uuid.uuid4().hex}.html"

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        from PySide6.QtWebEngineWidgets import QWebEngineView

        self._web = QWebEngineView(self)
        # Chromium's console chatter (e.g. Plotly's "Canvas2D: Multiple
        # readback operations ... willReadFrequently") went to the terminal;
        # route it to the Python log instead.
        self._web.setPage(_quiet_page(self._web))
        self._plot_ready = False  # a plot page (not a placeholder) loaded
        self._web.loadFinished.connect(self._on_load_finished)
        layout.addWidget(self._web)

        self._draw_placeholder("Nothing to display yet.")

    # ── Public API ─────────────────────────────────────────────────────

    @property
    def figure(self):
        """The last figure passed to :meth:`set_figure`, or ``None``."""
        return self._figure

    def set_figure(self, fig, *, post_script: str | None = None) -> None:
        """Render *fig* (a ``plotly.graph_objects.Figure``).

        *post_script* is JavaScript run once the plot exists (Plotly
        replaces ``{plot_id}`` with the plot's div id).
        """
        self._figure = fig
        html = fig.to_html(
            include_plotlyjs=True,
            full_html=True,
            config={"responsive": True, "displaylogo": False},
            post_script=post_script,
        )
        self._load_html(html)

    def show_figure(self, fig, *, post_script: str | None = None) -> None:
        """Show *fig*, keeping the camera when a plot is already on screen.

        The first figure loads a page (:meth:`set_figure`); later ones are
        swapped in with ``Plotly.react`` on that page, so the user's
        orbit/zoom and any running page script (spin) survive an option
        change instead of snapping back on a full reload.
        """
        if not self._plot_ready or self._figure is None:
            self.set_figure(fig, post_script=post_script)
            return
        self._figure = fig
        fig.update_layout(uirevision="pycsamt-keep")
        self.run_js(
            "(function(){var gd=document.querySelector('.js-plotly-plot');"
            "if(!gd){return;}var f=" + fig.to_json() + ";"
            "Plotly.react(gd,f.data,f.layout);})();")

    def run_js(self, code: str) -> None:
        """Run JavaScript in the page (e.g. to toggle an animation)."""
        self._web.page().runJavaScript(code)

    def eval_js(self, code: str, callback) -> None:
        """Evaluate *code* in the page; *callback* gets the result."""
        self._web.page().runJavaScript(code, 0, callback)

    @property
    def plot_ready(self) -> bool:
        return self._plot_ready

    def reload_next(self) -> None:
        """Make the next :meth:`show_figure` rebuild the page from
        scratch (a "hard" render) instead of swapping the figure in."""
        self._plot_ready = False

    def add_overlay_button(self, icon, tooltip: str, callback, *,
                           text: str = ""):
        """A small button floating over the scene's top-left corner
        (like the refresh button on the matplotlib canvases)."""
        from PySide6.QtCore import QSize
        from PySide6.QtWidgets import QToolButton

        b = QToolButton(self)
        if icon is not None and not icon.isNull():
            b.setIcon(icon)
            b.setIconSize(QSize(16, 16))
        else:
            b.setText(text or "↻")
        b.setToolTip(tooltip)
        b.setAutoRaise(False)
        b.setStyleSheet(
            "QToolButton { background: rgba(255,255,255,0.92); border: 1px "
            "solid #c3c9d4; border-radius: 6px; padding: 3px; }"
            "QToolButton:hover { border-color: #1864ab; }")
        b.clicked.connect(callback)
        b.adjustSize()
        self._overlays = getattr(self, "_overlays", []) + [b]
        self._place_overlays()
        b.raise_()
        return b

    def _place_overlays(self) -> None:
        x = 8
        for b in getattr(self, "_overlays", []):
            b.move(x, 8)
            b.raise_()
            x += b.width() + 4

    def resizeEvent(self, event) -> None:  # noqa: N802
        super().resizeEvent(event)
        self._place_overlays()

    def snapshot_png(self, path: str, done=None, *, scale: float = 2.0,
                     timeout_ms: int = 8000) -> None:
        """Save the scene as it is on screen (camera included) to a PNG.

        ``Plotly.toImage`` is asynchronous, so the image is parked on the
        page and collected by polling; *done(path_or_None, error)* is
        called at the end.
        """
        import base64

        self.run_js(
            "window.__pycsamtPng=null;(function(){var gd=document."
            "querySelector('.js-plotly-plot');if(!gd){window.__pycsamtPng="
            "'error:no plot';return;}Plotly.toImage(gd,{format:'png',scale:"
            f"{float(scale)}}}).then(function(u){{window.__pycsamtPng=u;}})"
            ".catch(function(e){window.__pycsamtPng='error:'+e;});})();")
        state = {"left": max(timeout_ms // 200, 1)}

        def poll():
            self.eval_js("window.__pycsamtPng", got)

        def got(value):
            if not value:
                state["left"] -= 1
                if state["left"] > 0:
                    QTimer.singleShot(200, poll)
                elif done:
                    done(None, "timed out")
                return
            if str(value).startswith("error:"):
                if done:
                    done(None, str(value)[6:])
                return
            data = str(value).split(",", 1)[-1]
            Path(path).write_bytes(base64.b64decode(data))
            if done:
                done(path, "")

        QTimer.singleShot(200, poll)

    def clear(self) -> None:
        self._figure = None
        self._plot_ready = False
        self._draw_placeholder("Nothing to display yet.")

    def export_html(self, path: str) -> None:
        """Write the current figure to a standalone HTML file at *path*."""
        if self._figure is None:
            raise ValueError("No figure to export.")
        self._figure.write_html(path, include_plotlyjs=True, full_html=True)

    # ── Internal ──────────────────────────────────────────────────────

    def _draw_placeholder(self, message: str) -> None:
        html = (
            "<html><body style='display:flex;align-items:center;"
            "justify-content:center;height:100vh;margin:0;"
            "font-family:sans-serif;color:#888;'>"
            f"<span>{message}</span></body></html>"
        )
        self._load_html(html)

    def _on_load_finished(self, ok: bool) -> None:
        self._plot_ready = bool(ok) and self._figure is not None

    def _load_html(self, html: str) -> None:
        self._plot_ready = False
        self._html_path.write_text(html, encoding="utf-8")
        self._web.load(QUrl.fromLocalFile(str(self._html_path)))

    def closeEvent(self, event) -> None:  # noqa: N802
        self._cleanup_tmp()
        super().closeEvent(event)

    def _cleanup_tmp(self) -> None:
        import shutil

        try:
            shutil.rmtree(self._tmp_dir, ignore_errors=True)
        except Exception:
            pass
