# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.app.converter.pages._base
====================================

Shared page base class: every page is a drop zone (or a few), an
options form, a Convert button, and the shared
:class:`~pycsamt.app.converter.widgets.log_panel.LogPanel`. This class
factors out running a :mod:`pycsamt.app.converter.jobs` function on a
background :class:`~pycsamt.app.converter.workers.ConversionWorker` and
reporting the outcome, so each concrete page only wires up its own
inputs and calls :meth:`ConverterPage.run_job`.
"""

from __future__ import annotations

from typing import Any, Callable

from PySide6.QtWidgets import (
    QFrame,
    QGridLayout,
    QLabel,
    QScrollArea,
    QVBoxLayout,
    QWidget,
)

from pycsamt.app.converter.widgets import LogPanel
from pycsamt.app.converter.workers import ConversionWorker


def page_root_layout(widget: QWidget) -> QVBoxLayout:
    """A consistently-margined/spaced, *scrollable* top-level layout for a page.

    Every page uses this instead of a bare ``QVBoxLayout(self)``. Pages with
    many advanced-option rows (inversion, PCBH, ...) can exceed the visible
    window height; without a scroll area Qt has no way to give every row its
    natural size and instead compresses them all toward their minimum size,
    which is what made rows overlap/collide in cramped windows. Wrapping the
    real content in a resizable ``QScrollArea`` fixes that regardless of
    window size, and the returned layout is exactly what callers add widgets
    to, so no page code needs to change.
    """
    outer = QVBoxLayout(widget)
    outer.setContentsMargins(0, 0, 0, 0)
    outer.setSpacing(0)

    scroll = QScrollArea()
    scroll.setObjectName("PageScrollArea")
    scroll.setFrameShape(QFrame.Shape.NoFrame)
    scroll.setWidgetResizable(True)
    outer.addWidget(scroll)

    content = QWidget()
    content.setObjectName("PageContent")
    scroll.setWidget(content)

    layout = QVBoxLayout(content)
    layout.setContentsMargins(18, 16, 18, 16)
    layout.setSpacing(10)
    return layout


def two_col_form(rows: list[tuple[str, QWidget]]) -> QGridLayout:
    """A 2-column grid of (label, widget) pairs for short scalar fields.

    Advanced-option sections tend to accumulate many one-line fields
    (EPSG, UTM zone, thresholds, ...); laying each out on its own
    ``QFormLayout`` row burns roughly twice the vertical space a 2-up grid
    needs. DropZones, checkboxes and anything else that wants the full
    width should stay out of this and be added to the page's own layout
    directly.
    """
    grid = QGridLayout()
    grid.setHorizontalSpacing(12)
    grid.setVerticalSpacing(8)
    for i, (label_text, field) in enumerate(rows):
        r, c = divmod(i, 2)
        label = QLabel(label_text)
        # A long QSpinBox special-value string ("final (default)") or a
        # wide numeric range otherwise forces its *sizeHint* -- and hence
        # the whole grid's minimum width -- well past what a narrow window
        # can offer, which is what was triggering an unwanted horizontal
        # scrollbar. Capping each field's minimum width lowers that floor;
        # the stretch factor below still lets fields grow to fill
        # whatever space is actually available.
        field.setMinimumWidth(90)
        grid.addWidget(label, r, c * 2)
        grid.addWidget(field, r, c * 2 + 1)
    grid.setColumnStretch(1, 1)
    grid.setColumnStretch(3, 1)
    return grid


def mark_primary(button) -> None:
    """Style *button* as the page's main call-to-action (accent-filled)."""
    button.setProperty("class", "primary")


class ConverterPage(QWidget):
    """Base class for every conversion page."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.log = LogPanel()
        self._worker: ConversionWorker | None = None
        self._run_button: Any = None

    def run_job(
        self,
        fn: Callable[..., Any],
        *args: Any,
        run_button: Any = None,
        supports_progress: bool = False,
        on_success: Callable[[Any], None] | None = None,
        success_message: Callable[[Any], str] | str | None = None,
        busy_message: str = "Working…",
        **kwargs: Any,
    ) -> None:
        """Run *fn* on a background thread and report the outcome in ``self.log``."""
        button = run_button or self._run_button
        if button is not None:
            button.setEnabled(False)
        self.log.clear()
        self.log.append(busy_message)
        if not supports_progress:
            self.log.set_indeterminate(True)

        worker = ConversionWorker(fn, *args, supports_progress=supports_progress, **kwargs)

        def _on_progress(cur: int, total: int, name: str) -> None:
            self.log.set_progress(cur, total)
            self.log.append(f"  [{cur}/{total}] {name}")

        def _on_done(result: Any) -> None:
            self.log.set_indeterminate(False)
            self.log.set_progress(1, 1)
            if button is not None:
                button.setEnabled(True)
            if success_message is not None:
                msg = success_message(result) if callable(success_message) else success_message
                self.log.append(msg)
            else:
                self.log.append("Done.")
            if on_success is not None:
                on_success(result)

        def _on_error(message: str) -> None:
            self.log.set_indeterminate(False)
            if button is not None:
                button.setEnabled(True)
            self.log.append(f"Error: {message}")

        worker.progress.connect(_on_progress)
        worker.done.connect(_on_done)
        worker.error.connect(_on_error)
        worker.start()
        self._worker = worker  # keep a reference alive while it runs
