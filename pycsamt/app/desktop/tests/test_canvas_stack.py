"""CanvasResultView — placeholder-until-drawn plot container."""
import pytest

pytest.importorskip("PySide6")
from matplotlib.figure import Figure

from pycsamt.app.desktop.widgets.canvas_stack import CanvasResultView


def test_starts_on_unavailable_card_and_switches_to_canvas(qapp):
    view = CanvasResultView(
        empty_title="Load data",
        empty_reason="Nothing loaded yet.",
    )
    assert view.showing_canvas is False
    assert view._stack.currentWidget() is view._unavailable

    fig = Figure()
    fig.add_subplot().plot([0, 1], [0, 1])
    view.canvas.show_figure(fig)
    view.show_canvas()
    assert view.showing_canvas is True
    assert view._stack.currentWidget() is view.canvas

    view.show_unavailable("Result unavailable", "The plot could not be drawn.")
    assert view.showing_canvas is False
    assert view._unavailable._title.text() == "Result unavailable"
    view.close()


def test_refresh_callback_still_works_through_the_wrapper(qapp):
    calls = []
    view = CanvasResultView()
    view.canvas.set_refresh_callback(lambda: calls.append(1))
    view.canvas.refresh_button.click()
    assert calls == [1]
    view.close()
