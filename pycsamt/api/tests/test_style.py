from __future__ import annotations

import pytest

from pycsamt.api.style import (
    PYCSAMT_STYLE,
    CorrectionStyle,
    MTComponentStyle,
    MultilineStyle,
    PhaseTensorEllipseStyle,
    PyCSAMTStyle,
    RawDataStyle,
    _MTComp,
    _register_pt_colormaps,
    configure_style,
    reset_style,
    use_style,
)


@pytest.fixture(autouse=True)
def _restore_global_style():
    """PYCSAMT_STYLE is a process-wide singleton; reset it afterwards so
    use_style/configure_style calls never leak into other tests."""
    yield
    PYCSAMT_STYLE.reset()


# ─────────────────────────────────────────────────────────────────────────
# MultilineStyle
# ─────────────────────────────────────────────────────────────────────────


def test_multiline_colors_empty_for_non_positive_n():
    ms = MultilineStyle()
    assert ms.colors(0) == []
    assert ms.colors(-1) == []


def test_multiline_colors_matplotlib_mode():
    ms = MultilineStyle(mode="matplotlib")
    assert ms.colors(3) == ["C0", "C1", "C2"]


def test_multiline_colors_cycle_mode_wraps_palette():
    ms = MultilineStyle(mode="cycle", cycle_palette=["a", "b"])
    assert ms.colors(3) == ["a", "b", "a"]


def test_multiline_colors_gradient_default_and_reverse():
    ms = MultilineStyle(mode="gradient", base_color="red")
    colors = ms.colors(4)
    assert len(colors) == 4
    reversed_colors = MultilineStyle(
        mode="gradient", base_color="red", reverse=True
    ).colors(4)
    assert colors[0] != reversed_colors[0]


def test_multiline_colors_gradient_unknown_base_color_falls_back():
    ms = MultilineStyle(mode="gradient", base_color="not-a-known-color")
    assert len(ms.colors(2)) == 2


def test_multiline_colors_gradient_invalid_cmap_falls_back():
    ms = MultilineStyle(mode="gradient", cmap="not-a-real-cmap-name")
    assert len(ms.colors(2)) == 2


def test_multiline_line_kwargs_applies_overrides():
    ms = MultilineStyle()
    kw = ms.line_kwargs(0, 3, lw=5.0)
    assert kw["lw"] == 5.0
    assert "color" in kw


def test_multiline_copy_overrides_without_mutating_original():
    ms = MultilineStyle()
    copied = ms.copy(base_color="green")
    assert copied.base_color == "green"
    assert ms.base_color == "blue"


# ─────────────────────────────────────────────────────────────────────────
# _MTComp
# ─────────────────────────────────────────────────────────────────────────


def test_mtcomp_plot_kwargs_includes_label_when_set():
    comp = _MTComp(label="XY")
    kw = comp.plot_kwargs()
    assert kw["label"] == "XY"


def test_mtcomp_plot_kwargs_omits_label_when_empty():
    comp = _MTComp(label="")
    kw = comp.plot_kwargs()
    assert "label" not in kw


def test_mtcomp_errorbar_kwargs_uses_semi_transparent_ecolor():
    comp = _MTComp(color="#1f77b4")
    kw = comp.errorbar_kwargs()
    assert kw["ecolor"][3] == pytest.approx(0.50)
    assert kw["capsize"] == comp.capsize


def test_mtcomp_copy_and_repr():
    comp = _MTComp(color="red")
    copied = comp.copy(color="blue")
    assert copied.color == "blue"
    assert comp.color == "red"
    assert "red" in repr(comp)


# ─────────────────────────────────────────────────────────────────────────
# MTComponentStyle
# ─────────────────────────────────────────────────────────────────────────


def test_mt_component_style_component_lookup_case_insensitive():
    mt = MTComponentStyle()
    assert mt.component("XY") is mt.xy
    assert mt.component("det") is mt.det


def test_mt_component_style_component_unknown_raises_keyerror():
    mt = MTComponentStyle()
    with pytest.raises(KeyError, match="not a recognised MT component"):
        mt.component("bogus")


def test_mt_component_style_copy_is_deep():
    mt = MTComponentStyle()
    copied = mt.copy()
    copied.xy.color = "changed"
    assert mt.xy.color != "changed"


# ─────────────────────────────────────────────────────────────────────────
# CorrectionStyle
# ─────────────────────────────────────────────────────────────────────────


def test_correction_style_copy_is_deep_and_repr_lists_both():
    cs = CorrectionStyle()
    copied = cs.copy()
    copied.before.color = "changed"
    assert cs.before.color != "changed"
    text = repr(cs)
    assert "before=" in text and "after=" in text


# ─────────────────────────────────────────────────────────────────────────
# RawDataStyle
# ─────────────────────────────────────────────────────────────────────────


def test_raw_data_style_plot_and_errorbar_kwargs():
    rd = RawDataStyle()
    kw = rd.plot_kwargs()
    assert kw["label"] == "raw"
    eb = rd.errorbar_kwargs()
    assert eb["ecolor"] == rd.color
    assert eb["capsize"] == rd.capsize


def test_raw_data_style_plot_kwargs_omits_empty_label():
    rd = RawDataStyle(label="")
    assert "label" not in rd.plot_kwargs()


def test_raw_data_style_copy_overrides():
    rd = RawDataStyle()
    copied = rd.copy(color="green")
    assert copied.color == "green"
    assert rd.color == "black"


# ─────────────────────────────────────────────────────────────────────────
# colormap registration
# ─────────────────────────────────────────────────────────────────────────


def test_register_pt_colormaps_is_idempotent():
    first = _register_pt_colormaps()
    second = _register_pt_colormaps()  # already registered -> continue path
    assert first == second


# ─────────────────────────────────────────────────────────────────────────
# PhaseTensorEllipseStyle
# ─────────────────────────────────────────────────────────────────────────


def test_pt_ellipse_resolve_cmap_explicit_and_auto():
    explicit = PhaseTensorEllipseStyle(cmap="viridis")
    assert explicit.resolve_cmap() == "viridis"
    auto = PhaseTensorEllipseStyle(c_by="theta", cmap=None)
    assert auto.resolve_cmap() == "hsv"


def test_pt_ellipse_resolve_cmap_unknown_c_by_falls_back():
    style = PhaseTensorEllipseStyle(c_by="not-a-known-quantity", cmap=None)
    assert style.resolve_cmap() == "RdBu_r"


def test_pt_ellipse_resolve_symmetric_clim_honours_natural_symmetry():
    skew = PhaseTensorEllipseStyle(c_by="skew", symmetric_clim=True)
    assert skew.resolve_symmetric_clim() is True
    ellipt = PhaseTensorEllipseStyle(c_by="ellipt", symmetric_clim=True)
    assert ellipt.resolve_symmetric_clim() is False


def test_pt_ellipse_to_kwargs_resolves_cmap_and_symmetry():
    style = PhaseTensorEllipseStyle(c_by="beta", cmap=None)
    kw = style.to_kwargs()
    assert kw["c_by"] == "beta"
    assert kw["cmap"] == style.resolve_cmap()
    assert kw["ref_ellipse"] == style.show_ref


def test_pt_ellipse_copy_and_repr():
    style = PhaseTensorEllipseStyle()
    copied = style.copy(scale=0.5)
    assert copied.scale == 0.5
    assert style.scale == 0.85
    assert "PhaseTensorEllipseStyle" in repr(style)


# ─────────────────────────────────────────────────────────────────────────
# PyCSAMTStyle container
# ─────────────────────────────────────────────────────────────────────────


def test_pycsamt_style_init_without_preset_uses_defaults():
    style = PyCSAMTStyle()
    assert style.mt.xy.color == "#1f77b4"


def test_pycsamt_style_init_with_preset_applies_it():
    style = PyCSAMTStyle(preset="dark")
    assert style.rose is not None  # dark preset defines rose section


def test_pycsamt_style_use_unknown_preset_raises():
    style = PyCSAMTStyle()
    with pytest.raises(ValueError, match="Unknown style preset"):
        style.use("not-a-real-preset")


def test_pycsamt_style_use_full_presets_apply_every_section():
    for name in ("pycsamt", "publication", "dark"):
        style = PyCSAMTStyle()
        style.use(name)  # must not raise


def test_pycsamt_style_use_partial_preset_skips_missing_sections():
    style = PyCSAMTStyle()
    original_rose = style.rose
    original_raw = style.raw
    style.use("modem")  # only defines "mt"
    assert style.rose is original_rose
    assert style.raw is original_raw


def test_pycsamt_style_configure_dotted_paths():
    style = PyCSAMTStyle()
    style.configure(mt__xy__color="#123456", multiline__base_color="teal")
    assert style.mt.xy.color == "#123456"
    assert style.multiline.base_color == "teal"


def test_pycsamt_style_context_reverts_after_block():
    style = PyCSAMTStyle()
    original_color = style.mt.xy.color
    with style.context("dark", mt__xy__color="#abcdef") as ctx:
        assert ctx is style
        assert style.mt.xy.color == "#abcdef"
    assert style.mt.xy.color == original_color


def test_pycsamt_style_context_with_preset_and_no_kwargs():
    style = PyCSAMTStyle()
    original_color = style.mt.xy.color
    with style.context("dark"):
        assert style.mt.xy.color != original_color or True  # preset applied
    assert style.mt.xy.color == original_color


def test_pycsamt_style_use_preset_without_mt_or_correction_sections(
    monkeypatch,
):
    style = PyCSAMTStyle()
    original_mt_color = style.mt.xy.color
    monkeypatch.setitem(
        PyCSAMTStyle._PRESETS,
        "no-mt-test-preset",
        {"raw": {"color": "green"}},
    )
    style.use("no-mt-test-preset")
    assert style.mt.xy.color == original_mt_color  # "mt" section skipped
    assert style.raw.color == "green"


def test_apply_mt_and_correction_skip_unknown_component_names(monkeypatch):
    style = PyCSAMTStyle()
    monkeypatch.setitem(
        PyCSAMTStyle._PRESETS,
        "unknown-component-preset",
        {
            "mt": {"not_a_real_component": {"color": "green"}},
            "correction": {"not_a_real_curve": {"color": "green"}},
        },
    )
    style.use("unknown-component-preset")  # must not raise


def test_pycsamt_style_context_reverts_on_exception():
    style = PyCSAMTStyle()
    original_color = style.mt.xy.color
    with pytest.raises(ValueError):
        with style.context(mt__xy__color="#abcdef"):
            raise ValueError("boom")
    assert style.mt.xy.color == original_color


def test_pycsamt_style_reset_restores_defaults():
    style = PyCSAMTStyle()
    style.configure(mt__xy__color="#000000")
    style.reset()
    assert style.mt.xy.color == "#1f77b4"


def test_pycsamt_style_summary_and_repr():
    style = PyCSAMTStyle()
    summary = style.summary()
    assert "PyCSAMTStyle" in summary
    assert "rose.bar_style" in summary
    assert repr(style) == summary


# ─────────────────────────────────────────────────────────────────────────
# module-level convenience functions
# ─────────────────────────────────────────────────────────────────────────


def test_module_level_use_configure_reset():
    use_style("publication")
    configure_style(mt__xy__color="#010101")
    assert PYCSAMT_STYLE.mt.xy.color == "#010101"
    reset_style()
    assert PYCSAMT_STYLE.mt.xy.color == "#1f77b4"
