"""Regression checks for responsive, high-DPI desktop plots."""
import pytest

pytest.importorskip("PySide6")
from matplotlib.figure import Figure
from pycsamt.api.control import PYCSAMT_CONTROL
from pycsamt.app.desktop.widgets.mpl_canvas import MplCanvas


def test_fit_uses_physical_pixels_and_keeps_window_size(qapp):
    view = MplCanvas(toolbar=False)
    view.resize(900, 600)
    view.show()
    qapp.processEvents()
    view._canvas._set_device_pixel_ratio(1.5)
    view.figure.set_size_inches(3, 2, forward=False)
    view.fit_to_view()
    assert view.figure.bbox.width == pytest.approx(view._canvas.width() * 1.5)
    assert view.figure.bbox.height == pytest.approx(view._canvas.height() * 1.5)
    assert (view.width(), view.height()) == (900, 600)
    view.close()


def test_colorbar_labels_fit_and_background_follows_config(qapp):
    view = MplCanvas(toolbar=False)
    view.resize(900, 600)
    view.show()
    image = view.axes.imshow([[1, 2], [3, 4]], aspect="auto")
    view.figure.colorbar(image, ax=view.axes, label="Resistivity (ohm m)")
    view.axes.set_xlabel("Distance along survey (m)")
    qapp.processEvents()
    view.fit_to_view()
    view._canvas.draw()
    renderer = view._canvas.get_renderer()
    for ax in view.figure.axes:
        bounds = ax.get_tightbbox(renderer)
        assert bounds.x0 >= 0 and bounds.y0 >= 0
        assert bounds.x1 <= view.figure.bbox.width
        assert bounds.y1 <= view.figure.bbox.height
    assert view.figure.get_facecolor() == (1, 1, 1, 1)
    with PYCSAMT_CONTROL.context(panel__background="transparent"):
        view.fit_to_view()
        assert view.figure.get_facecolor()[-1] == 0
    view.close()


def test_replacement_navigation_and_detached_snapshot(qapp):
    view = MplCanvas()
    figure = Figure()
    figure.add_subplot().plot([0, 1], [0, 1])
    view.show_figure(figure)
    assert view._toolbar.canvas is view._canvas
    view.open_detached()
    assert view._detached_dialog is not None
    detached = view._detached_dialog.findChild(MplCanvas)
    assert detached.figure is not view.figure
    detached.axes.set_xlim(10, 20)
    assert view.axes.get_xlim() != detached.axes.get_xlim()
    qapp.processEvents()
    view._detached_dialog.close()
    qapp.processEvents()
    view.close()


def test_open_detached_twice_reraises_instead_of_duplicating(qapp):
    # Regression: clicking "open in separate window" while one is already
    # open used to spawn another dialog every time (an unbounded pile of
    # windows). It must now disable the toolbar action and re-raise the
    # existing dialog instead.
    view = MplCanvas()
    figure = Figure()
    figure.add_subplot().plot([0, 1], [0, 1])
    view.show_figure(figure)

    view.open_detached()
    first = view._detached_dialog
    assert first is not None
    assert view._detach_action.isEnabled() is False

    view.open_detached()
    assert view._detached_dialog is first  # no second dialog created

    qapp.processEvents()
    first.close()
    qapp.processEvents()
    assert view._detached_dialog is None
    assert view._detach_action.isEnabled() is True
    view.close()


def test_detach_action_disabled_state_survives_toolbar_rebuild(qapp):
    # _build_toolbar() reruns on every show_figure() call, recreating the
    # QAction object from scratch — the disabled state must be reapplied
    # from self._detached_dialog rather than living only on the stale
    # action reference.
    view = MplCanvas()
    figure = Figure()
    figure.add_subplot().plot([0, 1], [0, 1])
    view.show_figure(figure)
    view.open_detached()
    assert view._detached_dialog is not None

    other_figure = Figure()
    other_figure.add_subplot().plot([1, 2], [2, 3])
    view.show_figure(other_figure)  # rebuilds the toolbar mid-detach
    assert view._detach_action.isEnabled() is False

    qapp.processEvents()
    view._detached_dialog.close()
    qapp.processEvents()
    view.close()


def test_refresh_overlay_sits_below_toolbar_not_over_it(qapp):
    # Regression: a floating refresh button anchored to the top-right of
    # the whole canvas widget landed directly on top of the built-in
    # toolbar's own icons (including "open in separate window", which
    # sits at the toolbar's right edge). It must sit below the toolbar
    # strip instead, clear of every toolbar icon.
    view = MplCanvas(toolbar=True)
    view.resize(900, 600)
    view.show()
    qapp.processEvents()

    calls = []
    view.set_refresh_callback(lambda: calls.append(1))
    assert view.refresh_button is not None
    assert view.refresh_button.y() >= view._tools.height()

    view.refresh_button.click()
    assert calls == [1]
    view.close()


def test_refresh_overlay_reassigns_callback_without_stacking(qapp):
    view = MplCanvas(toolbar=True)
    calls = []
    view.set_refresh_callback(lambda: calls.append("first"))
    view.set_refresh_callback(lambda: calls.append("second"))
    view.refresh_button.click()
    assert calls == ["second"]  # only the latest callback fires, once
    view.close()


def test_refresh_overlay_repositions_when_toolbar_hidden(qapp):
    from pycsamt.api.control import PYCSAMT_CONTROL

    view = MplCanvas(toolbar=True)
    view.resize(900, 600)
    view.show()
    qapp.processEvents()
    view.set_refresh_callback(lambda: None)
    with PYCSAMT_CONTROL.context(panel__toolbar=False):
        view.draw()  # triggers _apply_panel_style()
        assert view.refresh_button.y() < view._tools.height() + 1
    view.close()


def test_manual_colorbar_and_long_bottom_labels_remain_visible(qapp):
    view = MplCanvas(toolbar=False)
    view.figure.clear()
    ax = view.figure.add_axes([0.10, 0.04, 0.70, 0.85])
    cax = view.figure.add_axes([0.85, 0.04, 0.025, 0.85])
    image = ax.imshow([[1, 2], [3, 4]], aspect="auto")
    view.figure.colorbar(image, cax=cax)
    ax.set_xticks([0, 1], ["Long station identifier 001", "Long station identifier 002"], rotation=45)
    ax.set_xlabel("Distance along survey line\nMeasured from the first station (km)")
    view.resize(850, 450)
    view.show()
    qapp.processEvents()
    for height in (450, 340, 650):
        view.resize(850, height)
        qapp.processEvents()
        view.fit_to_view()
        view._canvas.draw()
        bounds = ax.get_tightbbox(view._canvas.get_renderer())
        assert bounds.y0 >= 8
        assert bounds.y1 <= view.figure.bbox.height
    view.close()


def test_label_boxes_stay_white_after_a_dark_restyle(qapp):
    """The phase-tensor |β| legend / size reference / PT-strip station label
    sit in boxes a dark restyle painted near-black, while the canvas set the
    text dark slate: unreadable.  The publication style resets the boxes."""
    from matplotlib.colors import to_hex

    from pycsamt.app.desktop.widgets.mpl_canvas import MplCanvas

    c = MplCanvas()
    try:
        ax = c.figure.add_subplot(111)
        t = ax.text(0.1, 0.9, "|β| < 3°", transform=ax.transAxes,
                    bbox=dict(fc="white"))
        t.get_bbox_patch().set_facecolor("#1a1a2e")  # dark restyle
        c.apply_theme(True)
        assert to_hex(t.get_bbox_patch().get_facecolor()) == "#ffffff"
    finally:
        c.close()
