# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
CanvasResultView — an MplCanvas paired with a QStackedWidget-hosted
"no result yet" placeholder.

A bare MplCanvas shows an empty white subplot with default 0-1 axes,
ticks and spines the moment it is constructed, and again whenever the
caller has nothing to draw. That reads as a broken plot rather than an
intentional empty state. This widget hides the canvas behind a
:class:`~pycsamt.app.desktop.widgets.unavailable_view.UnavailableResultView`
card until a real figure exists, matching the pattern first used by
``QCDashboardWindow``.

Usage::

    view = CanvasResultView(
        parent,
        empty_title="Load survey data to begin",
        empty_reason="No stations are currently available.",
        empty_guidance="Load EDI or EMTF-XML files, then click Refresh.",
    )
    layout.addWidget(view)               # add the container, not view.canvas
    view.canvas.set_refresh_callback(self._on_refresh)  # unchanged API

    # once data/figure is ready:
    view.canvas.show_figure(fig)
    view.show_canvas()

    # if there is nothing to draw:
    view.show_unavailable("Result unavailable", str(exc))
"""

from __future__ import annotations

from PySide6.QtWidgets import QStackedWidget, QVBoxLayout, QWidget

from pycsamt.app.desktop.widgets.mpl_canvas import MplCanvas
from pycsamt.app.desktop.widgets.unavailable_view import UnavailableResultView


class CanvasResultView(QWidget):
    """A plot canvas that starts hidden behind a "nothing to show" card."""

    def __init__(
        self,
        parent: QWidget | None = None,
        toolbar: bool = True,
        empty_title: str = "No plot yet",
        empty_reason: str = "Nothing has been computed for this view yet.",
        empty_guidance: str = "",
    ) -> None:
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)

        self.canvas = MplCanvas(self, toolbar=toolbar)
        self._unavailable = UnavailableResultView(self)
        self._unavailable.set_content(empty_title, empty_reason, empty_guidance)

        self._stack = QStackedWidget(self)
        self._stack.addWidget(self.canvas)
        self._stack.addWidget(self._unavailable)
        self._stack.setCurrentWidget(self._unavailable)
        layout.addWidget(self._stack)

    def show_canvas(self) -> None:
        """Reveal the plot canvas, hiding the placeholder card."""
        self._stack.setCurrentWidget(self.canvas)

    def show_unavailable(
        self, title: str, reason: str = "", guidance: str = ""
    ) -> None:
        """Hide the canvas behind an explanatory placeholder card."""
        if reason or title:
            self._unavailable.set_content(
                title, reason or title, guidance
            )
        self._stack.setCurrentWidget(self._unavailable)

    @property
    def showing_canvas(self) -> bool:
        return self._stack.currentWidget() is self.canvas


__all__ = ["CanvasResultView"]
