# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Unit tests for the small/standalone helpers in :mod:`pycsamt.utils.plot`:
plot2d, set_axis_grid, is_valid_kind, make_plot_colors, savefigure,
resetting_ticks, make_mpl_properties, resetting_colorbar_bound,
controle_delineate_curve, fmt_text, get_color_palette,
_get_xticks_formatage, _set_sns_style, _format_ticks and _make_axe_multiple.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

from pycsamt.exceptions import PlotError
from pycsamt.utils.plot import (
    _format_ticks,
    _get_xticks_formatage,
    _make_axe_multiple,
    _set_sns_style,
    controle_delineate_curve,
    fmt_text,
    get_color_palette,
    is_valid_kind,
    make_mpl_properties,
    make_plot_colors,
    plot2d,
    resetting_colorbar_bound,
    resetting_ticks,
    savefigure,
    set_axis_grid,
)


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


# ------------------------------- plot2d --------------------------------


def test_plot2d_pcolormesh_default():
    ar = np.random.rand(6, 10) * 100 + 1
    ax = plot2d(ar)
    assert isinstance(ax, plt.Axes)
    assert ax.get_figure() is not None


def test_plot2d_imshow_style_and_top_label():
    ar = np.random.rand(4, 6) * 50 + 1
    ax = plot2d(ar, plt_style="imshow", top_label="Line", cb_label="Rho")
    assert isinstance(ax, plt.Axes)


def test_plot2d_contours_and_log10():
    ar = np.random.rand(5, 30) * 100 + 1
    ax = plot2d(ar, plot_contours=True, to_log10=True, y=np.arange(1, 6))
    assert isinstance(ax, plt.Axes)


def test_plot2d_custom_x_y_and_stnlist():
    ar = np.random.rand(3, 4) * 10 + 1
    x = np.array([0.0, 10.0, 20.0, 30.0])
    y = np.array([0.0, 1.0, 2.0])
    ax = plot2d(ar, x=x, y=y, stnlist=["A", "B", "C", "D"], how="m")
    assert isinstance(ax, plt.Axes)


def test_plot2d_external_axes():
    ar = np.random.rand(3, 5) * 10 + 1
    fig, ax_ext = plt.subplots()
    ax = plot2d(ar, ax=ax_ext)
    assert ax is ax_ext


def test_plot2d_invalid_ndim_raises():
    with pytest.raises(ValueError):
        plot2d(np.arange(10))


def test_plot2d_to_log10_without_y():
    ar = np.random.rand(4, 6) * 100 + 1
    ax = plot2d(ar, to_log10=True)
    assert isinstance(ax, plt.Axes)


# ---------------------------- set_axis_grid -----------------------------


def test_set_axis_grid_on_default_props():
    fig, ax = plt.subplots()
    set_axis_grid(ax, True)
    fig.canvas.draw()
    lines = ax.get_xgridlines() + ax.get_ygridlines()
    assert any(gl.get_visible() for gl in lines)


def test_set_axis_grid_off_hides_lines():
    fig, ax = plt.subplots()
    ax.plot([1, 2, 3], [1, 2, 3])
    set_axis_grid(ax, True)
    set_axis_grid(ax, False)
    for gl in ax.get_xgridlines() + ax.get_ygridlines():
        assert gl.get_visible() is False


def test_set_axis_grid_list_of_axes_and_custom_props():
    fig, axes = plt.subplots(1, 2)
    set_axis_grid(
        list(axes),
        True,
        grid_props={"which": "major", "axis": "x", "linestyle": "--"},
    )


# ---------------------------- is_valid_kind -----------------------------


@pytest.mark.parametrize(
    "kind,expected",
    [
        ("boxplot", "box"),
        ("box_plot", "box"),
        ("violinchart", "violin"),
        ("scattergraph", "scatter"),
        ("linechart", "line"),
        ("barplot", "bar"),
        ("histogram", "hist"),
        ("heat_map", "heatmap"),
        ("plot_density", "density"),
        ("areachart", "area"),
    ],
)
def test_is_valid_kind_aliases(kind, expected):
    assert is_valid_kind(kind) == expected


def test_is_valid_kind_pattern_fallback_and_unknown():
    assert is_valid_kind("my-scatter-thing") == "scatter"
    assert is_valid_kind("totally_unknown_kind") == "totallyunknownkind"


def test_is_valid_kind_validates_against_allowed_and_raises():
    assert is_valid_kind("barplot", valid_kinds=["bar", "line"]) == "bar"
    with pytest.raises(ValueError):
        is_valid_kind("scatterplot", valid_kinds=["bar", "line"])


def test_is_valid_kind_error_not_raise_returns_normalized():
    out = is_valid_kind(
        "scatterplot", valid_kinds=["bar", "line"], error="ignore"
    )
    assert out == "scatter"


# -------------------------- make_mpl_properties -------------------------


def test_make_mpl_properties_colors_within_palette():
    colors = make_mpl_properties(5)
    assert len(colors) == 5


def test_make_mpl_properties_colors_wraps_when_n_exceeds_palette():
    colors = make_mpl_properties(20)
    assert len(colors) == 20


def test_make_mpl_properties_markers_and_lines():
    markers = make_mpl_properties(15, prop="marker")
    assert len(markers) == 15
    lines = make_mpl_properties(30, prop="line")
    assert len(lines) == 30


def test_make_mpl_properties_zero_length():
    assert make_mpl_properties(0) == []


def test_make_mpl_properties_invalid_prop_raises():
    with pytest.raises(ValueError):
        make_mpl_properties(3, prop="bogus")


# ---------------------------- make_plot_colors ---------------------------


def test_make_plot_colors_default_axis0():
    ar = np.random.randn(7, 2)
    colors = make_plot_colors(ar)
    assert colors == ["g", "gray", "y", "blue", "orange", "purple", "lime"]


def test_make_plot_colors_axis1():
    ar = np.random.randn(7, 2)
    colors = make_plot_colors(ar, axis=1)
    assert colors == ["g", "gray"]


def test_make_plot_colors_cs4_and_chunking():
    ar = np.random.randn(7, 2)
    cs4 = make_plot_colors(ar, axis=1, colors="cs4")
    assert len(cs4) == 2
    full = make_plot_colors(ar, axis=1, colors="cs4", chunk=False)
    assert len(full) > 2


def test_make_plot_colors_cs4_with_start_index():
    ar = np.random.randn(7, 2)
    colors = make_plot_colors(ar, axis=1, colors="cs4:4")
    assert len(colors) == 2


def test_make_plot_colors_xkcd_with_seed_reproducible():
    ar = np.random.randn(4, 3)
    c1 = make_plot_colors(ar, axis=1, colors="xkcd", seed=42)
    c2 = make_plot_colors(ar, axis=1, colors="xkcd", seed=42)
    assert c1 == c2


def test_make_plot_colors_explicit_colors_list_prepended():
    ar = np.random.randn(3, 2)
    colors = make_plot_colors(ar, axis=1, colors=["red", "blue"])
    assert colors[0] == "red"
    assert colors[1] == "blue"


def test_make_plot_colors_single_string_color():
    ar = np.random.randn(3, 2)
    colors = make_plot_colors(ar, axis=1, colors="red")
    assert colors[0] == "red"


def test_make_plot_colors_1d_input():
    d = np.arange(5)
    colors = make_plot_colors(d)
    assert len(colors) == 5


def test_make_plot_colors_ellipsis_colors_treated_as_none():
    ar = np.random.randn(4, 2)
    colors = make_plot_colors(ar, axis=1, colors=...)
    assert colors == ["g", "gray"]


def test_make_plot_colors_plain_list_input_without_array_protocol():
    colors = make_plot_colors([1, 2, 3])
    assert len(colors) == 3


def test_make_plot_colors_cs4_named_start_key():
    ar = np.random.randn(3, 2)
    colors = make_plot_colors(ar, axis=1, colors="cs4:aliceblue")
    assert len(colors) == 2


def test_make_plot_colors_cs4_out_of_range_start_resets_to_zero():
    ar = np.random.randn(3, 2)
    colors = make_plot_colors(ar, axis=1, colors="cs4:999999")
    assert len(colors) == 2


# ------------------------------ savefigure -------------------------------


def test_savefigure_creates_file(tmp_path):
    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1])
    out = tmp_path / "myfig.png"
    savefigure(fig, str(out))
    assert out.exists()


def test_savefigure_without_extension_creates_directory_style_path(tmp_path):
    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1])
    out = tmp_path / "noext"
    out.mkdir()
    savefigure(fig, str(out), ext="png")
    # figname gets rewritten to <file>/.png when no extension is present;
    # matplotlib then appends its default format since ".png" has no
    # extension of its own (leading dot is not an extension separator).
    assert list(out.glob("*.png*"))


def test_savefigure_warns_when_no_name_given(tmp_path, monkeypatch):
    import pycsamt.utils.plot as plotmod

    class _FakeNow:
        def strftime(self, fmt):
            return "07-20-2026_120000"

    class _FakeDateTime:
        @staticmethod
        def now():
            return _FakeNow()

    monkeypatch.setattr(plotmod.datetime, "datetime", _FakeDateTime)
    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1])
    monkeypatch.chdir(tmp_path)
    with pytest.warns(UserWarning):
        savefigure(fig, None)
    assert list(tmp_path.glob("*.png"))


# ---------------------------- resetting_ticks -----------------------------


def test_resetting_ticks_many_values_default_n():
    out = resetting_ticks(list(range(0, 100, 10)))
    assert isinstance(out, np.ndarray)
    assert len(out) > 0


def test_resetting_ticks_explicit_number_of_ticks():
    out = resetting_ticks(list(range(0, 100, 5)), number_of_ticks=4)
    assert len(out) == 4


def test_resetting_ticks_single_value():
    out = resetting_ticks([5])
    assert list(out) == [5]


def test_resetting_ticks_invalid_type_raises():
    with pytest.raises(PlotError):
        resetting_ticks("not-a-list")


def test_resetting_ticks_non_multiple_of_ten_boundaries():
    out = resetting_ticks([1, 7, 15, 23, 31, 39, 47, 53])
    assert isinstance(out, np.ndarray)
    assert len(out) > 0


# -------------------------- resetting_colorbar_bound ----------------------


def test_resetting_colorbar_bound_exact_multiple():
    out = resetting_colorbar_bound(100, 10)
    assert np.allclose(out, np.linspace(10, 100, 5))


def test_resetting_colorbar_bound_non_multiple():
    out = resetting_colorbar_bound(97, 11, number_of_ticks=5)
    assert len(out) == 5


def test_resetting_colorbar_bound_logscale():
    out = resetting_colorbar_bound(100, 10, logscale=True)
    assert len(out) == 5


# -------------------------- controle_delineate_curve -----------------------


def test_controle_delineate_curve_scalar_resistivity():
    assert controle_delineate_curve(res_deline=100) == [2.0]


def test_controle_delineate_curve_scalar_phase():
    assert controle_delineate_curve(phase_deline=45.4) == [46.0]


def test_controle_delineate_curve_list_resistivity_and_phase():
    assert controle_delineate_curve(res_deline=[100, 1000]) == [2.0, 3.0]
    assert controle_delineate_curve(phase_deline=[10.1, 20.9]) == [11.0, 21.0]


def test_controle_delineate_curve_none_returns_none():
    assert controle_delineate_curve() is None


def test_controle_delineate_curve_invalid_scalar_raises():
    with pytest.raises(PlotError):
        controle_delineate_curve(res_deline="abc")


def test_controle_delineate_curve_invalid_phase_raises():
    # Regression check: a prior bug (fmt list holding a single
    # comma-joined string instead of two entries) made this branch
    # raise IndexError instead of the documented PlotError.
    with pytest.raises(PlotError):
        controle_delineate_curve(phase_deline="abc")


def test_controle_delineate_curve_invalid_list_raises():
    with pytest.raises(PlotError):
        controle_delineate_curve(res_deline=["a", "b"])


# -------------------------------- fmt_text ---------------------------------


def test_fmt_text_empty_text_skips_loop_body():
    out = fmt_text("")
    assert isinstance(out, str)
    assert out.count("~") > 0


def test_fmt_text_short_text():
    out = fmt_text("hello world")
    assert isinstance(out, str)
    assert "hello world" in out


def test_fmt_text_wraps_long_text():
    long_text = "word " * 40
    out = fmt_text(long_text, return_to_line=20)
    assert out.count("\n") > 2


def test_fmt_text_wraps_with_hyphen_when_mid_word():
    # The wrap point falls mid-word: the no-space branch inserts a
    # trailing hyphen before continuing on the next line.
    text = ("a" * 20) + " " + ("b" * 20)
    out = fmt_text(text, return_to_line=20)
    assert isinstance(out, str)
    assert "-" in out


def test_fmt_text_wraps_exactly_at_space_boundary():
    # The wrap point falls right before a space, exercising the
    # no-hyphen wrap branch.
    text = ("a" * 21) + " " + ("b" * 20)
    out = fmt_text(text, return_to_line=20)
    assert isinstance(out, str)
    assert "bb" in out


# ---------------------------- get_color_palette ----------------------------


def test_get_color_palette_numeric_string():
    val = get_color_palette("128")
    assert 0 <= val <= 1


def test_get_color_palette_numeric_float():
    val = get_color_palette(64.0)
    assert 0 <= val <= 1


def test_get_color_palette_rgb_string():
    rgba = get_color_palette("R128G64B32")
    assert isinstance(rgba, tuple)
    assert len(rgba) == 3


def test_get_color_palette_out_of_range_raises():
    with pytest.raises(ValueError):
        get_color_palette(300.0)


def test_get_color_palette_green_only():
    rgba = get_color_palette("g128")
    assert rgba[0] == 0.0 and rgba[2] == 0.0
    assert rgba[1] > 0


def test_get_color_palette_blue_only():
    rgba = get_color_palette("b64")
    assert rgba[0] == 0.0 and rgba[1] == 0.0
    assert abs(rgba[2] - 64 / 255.0) < 1e-9


def test_get_color_palette_full_rgb_string_parses_blue_correctly():
    # Regression test for a fixed bug: the blue-channel parser used to
    # split on 'g' only and read the wrong token, so it always fell
    # back to 1.0 regardless of the requested blue value.
    rgba = get_color_palette("R10G20B30")
    assert abs(rgba[0] - 10 / 255.0) < 1e-9
    assert abs(rgba[1] - 20 / 255.0) < 1e-9
    assert abs(rgba[2] - 30 / 255.0) < 1e-9


def test_get_color_palette_non_numeric_components_fallback_to_one():
    rgba = get_color_palette("rxgxbx")
    assert rgba == (1.0, 1.0, 1.0)


# -------------------------- _get_xticks_formatage ---------------------------


def test_get_xticks_formatage_many_ticks_uses_funcformatter():
    fig, ax = plt.subplots()
    ax.plot(range(20), range(20))
    _get_xticks_formatage(ax, list(range(20)), space=14, step=7)


def test_get_xticks_formatage_few_ticks_sets_labels_directly():
    fig, ax = plt.subplots()
    ax.plot(range(5), range(5))
    _get_xticks_formatage(ax, list(range(5)), space=14, step=7)


def test_get_xticks_formatage_y_axis_and_auto():
    fig, ax = plt.subplots()
    ax.plot(range(20), range(20))
    _get_xticks_formatage(ax, list(range(20)), ticks="y", auto=True)


def test_get_xticks_formatage_renders_formatter_callback():
    fig, ax = plt.subplots()
    ax.plot(range(20), range(20))
    _get_xticks_formatage(ax, list(range(20)), space=14, step=7)
    fig.canvas.draw()


# ------------------------------ _set_sns_style ------------------------------


def test_set_sns_style_darkgrid():
    _set_sns_style("darkgrid")


def test_set_sns_style_true_alias():
    _set_sns_style(True)


def test_set_sns_style_whitegrid():
    _set_sns_style("whitegrid")


# ------------------------------ _format_ticks --------------------------------


def test_format_ticks_multiple_of_nskip():
    assert _format_ticks(7, 0) == "S08"


def test_format_ticks_not_multiple_returns_none():
    assert _format_ticks(3, 0) is None


# ---------------------------- _make_axe_multiple -----------------------------


def test_make_axe_multiple_from_int():
    fig, axes = _make_axe_multiple(5, ncols=2)
    assert isinstance(fig, plt.Figure)
    assert axes.shape == (3, 2)


def test_make_axe_multiple_from_iterable():
    fig, axes = _make_axe_multiple(["a", "b", "c"], ncols=3)
    assert axes.shape == (3,)


def test_make_axe_multiple_zero_count_falls_back_to_one_row():
    fig, axes = _make_axe_multiple(0, ncols=3)
    assert isinstance(fig, plt.Figure)
    assert len(axes) == 3


def test_make_axe_multiple_external_axes_passthrough():
    fig0, ax0 = plt.subplots()
    fig, ax = _make_axe_multiple(3, fig=fig0, ax=ax0)
    assert ax is ax0
    assert fig is fig0
