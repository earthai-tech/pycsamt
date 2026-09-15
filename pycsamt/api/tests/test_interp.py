from __future__ import annotations

import pytest

from pycsamt.api.interp import (
    PYCSAMT_INTERP,
    HydroProfileStyle,
    HydroSectionStyle,
    InterpStyle,
    PyCSAMTInterp,
    configure_interp,
    reset_interp,
    resolve_figsize,
    resolve_profile_style,
    resolve_section_style,
    use_interp,
)


@pytest.fixture(autouse=True)
def _restore_global_interp():
    """PYCSAMT_INTERP is a process-wide singleton; snapshot/restore so
    tests that call use_interp/configure_interp/reset_interp never leak
    into other tests."""
    snapshot = PYCSAMT_INTERP._snapshot()
    yield
    PYCSAMT_INTERP._restore(snapshot)


# ─────────────────────────────────────────────────────────────────────────
# HydroSectionStyle
# ─────────────────────────────────────────────────────────────────────────


def test_hydro_section_style_cmap_for_known_quantities():
    style = HydroSectionStyle()
    assert style.cmap_for("K") == "viridis"
    assert style.cmap_for("saturation") == "RdYlBu"


def test_hydro_section_style_cmap_for_unknown_raises():
    style = HydroSectionStyle()
    with pytest.raises(ValueError, match="Unknown quantity"):
        style.cmap_for("bogus")


def test_hydro_section_style_kwargs_helpers_apply_overrides():
    style = HydroSectionStyle()
    assert style.wt_kwargs()["color"] == "deepskyblue"
    assert style.wt_kwargs(color="red")["color"] == "red"
    assert style.station_kwargs()["alpha"] == 0.45
    assert style.cb_kwargs()["fraction"] == 0.025
    assert style.crossplot_scatter_kwargs()["cmap"] == "viridis"
    assert style.hs_fill_kwargs()["alpha"] == 0.12
    assert style.model_curve_kwargs()["color"] == "k"
    assert style.rho_curve_kwargs()["color"] == "0.15"
    zb = style.zone_boundary_kwargs(linewidth=2.0)
    assert zb["linewidth"] == 2.0
    assert zb["color"] == "0.55"


def test_hydro_section_style_copy_overrides_without_mutating_original():
    style = HydroSectionStyle()
    copied = style.copy(cmap_K="plasma")
    assert copied.cmap_K == "plasma"
    assert style.cmap_K == "viridis"


# ─────────────────────────────────────────────────────────────────────────
# HydroProfileStyle
# ─────────────────────────────────────────────────────────────────────────


def test_hydro_profile_style_kwargs_helpers():
    style = HydroProfileStyle()
    assert style.envelope_kwargs("blue")["alpha"] == 0.25
    assert style.line_kwargs("blue")["linewidth"] == 2.0
    assert style.scatter_kwargs("blue")["s"] == 14.0
    assert style.ref_kwargs()["color"] == "tomato"
    assert style.station_kwargs()["color"] == "0.6"
    assert style.grid_kwargs() == {"alpha": 0.30, "axis": "y"}
    assert style.bar_kwargs("orange")["color"] == "orange"
    assert style.hist_kwargs("k")["bins"] == 30
    assert style.kde_kwargs()["color"] == "k"


def test_hydro_profile_style_copy_overrides_without_mutating_original():
    style = HydroProfileStyle()
    copied = style.copy(color_wt="gold")
    assert copied.color_wt == "gold"
    assert style.color_wt == "steelblue"


# ─────────────────────────────────────────────────────────────────────────
# InterpStyle
# ─────────────────────────────────────────────────────────────────────────


def test_interp_style_copy_is_a_deep_copy():
    style = InterpStyle()
    copied = style.copy(figsize_section=(1.0, 1.0))
    copied.section.cmap_K = "plasma"
    assert style.section.cmap_K == "viridis"  # untouched, deep-copied
    assert copied.figsize_section == (1.0, 1.0)
    assert style.figsize_section == (13.0, 5.0)


# ─────────────────────────────────────────────────────────────────────────
# PyCSAMTInterp container
# ─────────────────────────────────────────────────────────────────────────


def test_pycsamt_interp_presets_are_distinct_interp_styles():
    interp = PyCSAMTInterp()
    for name in ("default", "publication", "dark", "accessible"):
        assert isinstance(getattr(interp, name), InterpStyle)


def test_style_for_valid_and_invalid_preset():
    interp = PyCSAMTInterp()
    assert interp.style_for("PUBLICATION ") is interp.publication
    with pytest.raises(ValueError, match="interp preset must be one of"):
        interp.style_for("bogus")


def test_use_copies_preset_into_default():
    interp = PyCSAMTInterp()
    interp.publication.section.cmap_K = "marker-value"
    interp.use("publication")
    assert interp.default.section.cmap_K == "marker-value"
    # deep copy: mutating publication afterwards must not affect default
    interp.publication.section.cmap_K = "changed-again"
    assert interp.default.section.cmap_K == "marker-value"


def test_configure_sets_nested_dotted_paths():
    interp = PyCSAMTInterp()
    interp.configure(section__cmap_K="plasma", figsize_section=(1.0, 2.0))
    assert interp.default.section.cmap_K == "plasma"
    assert interp.default.figsize_section == (1.0, 2.0)


def test_reset_restores_package_defaults():
    interp = PyCSAMTInterp()
    interp.configure(section__cmap_K="plasma")
    interp.reset()
    assert interp.default.section.cmap_K == "viridis"


def test_context_reverts_default_after_block():
    interp = PyCSAMTInterp()
    original = interp.default.section.cmap_K
    with interp.context("dark", section__cmap_K="mid-block") as ctx:
        assert ctx is interp
        assert interp.default.section.cmap_K == "mid-block"
    assert interp.default.section.cmap_K == original


def test_context_without_preset_only_applies_kwargs():
    interp = PyCSAMTInterp()
    with interp.context(section__cmap_K="only-kwargs"):
        assert interp.default.section.cmap_K == "only-kwargs"
    assert interp.default.section.cmap_K == "viridis"


def test_context_reverts_on_exception():
    interp = PyCSAMTInterp()
    with pytest.raises(ValueError):
        with interp.context("dark"):
            raise ValueError("boom")
    assert interp.default.section.cmap_K == "viridis"


def test_summary_and_repr_list_every_preset():
    interp = PyCSAMTInterp()
    summary = interp.summary()
    assert "publication" in summary
    assert "active" in summary
    assert repr(interp) == summary


# ─────────────────────────────────────────────────────────────────────────
# module-level convenience functions
# ─────────────────────────────────────────────────────────────────────────


def test_module_level_configure_use_reset():
    use_interp("dark")
    assert PYCSAMT_INTERP.default.section.cmap_K == PYCSAMT_INTERP.dark.section.cmap_K
    configure_interp(section__cmap_K="tweaked")
    assert PYCSAMT_INTERP.default.section.cmap_K == "tweaked"
    reset_interp()
    assert PYCSAMT_INTERP.default.section.cmap_K == "viridis"


# ─────────────────────────────────────────────────────────────────────────
# resolve_section_style / resolve_profile_style / resolve_figsize
# ─────────────────────────────────────────────────────────────────────────


def test_resolve_section_style_variants():
    assert resolve_section_style(None) is PYCSAMT_INTERP.default.section
    assert resolve_section_style("dark") is PYCSAMT_INTERP.dark.section
    bundle = InterpStyle()
    assert resolve_section_style(bundle) is bundle.section
    leaf = HydroSectionStyle()
    assert resolve_section_style(leaf) is leaf
    with pytest.raises(TypeError, match="style must be"):
        resolve_section_style(42)


def test_resolve_profile_style_variants():
    assert resolve_profile_style(None) is PYCSAMT_INTERP.default.profile
    assert resolve_profile_style("accessible") is PYCSAMT_INTERP.accessible.profile
    bundle = InterpStyle()
    assert resolve_profile_style(bundle) is bundle.profile
    leaf = HydroProfileStyle()
    assert resolve_profile_style(leaf) is leaf
    with pytest.raises(TypeError, match="style must be"):
        resolve_profile_style(42)


def test_resolve_figsize_prefers_user_value():
    assert resolve_figsize((9.0, 9.0), None, "section") == (9.0, 9.0)


def test_resolve_figsize_from_named_preset():
    assert resolve_figsize(None, "publication", "profile") == (
        PYCSAMT_INTERP.publication.figsize_profile
    )


def test_resolve_figsize_from_interp_style_instance():
    bundle = InterpStyle(figsize_uncertainty=(2.0, 3.0))
    assert resolve_figsize(None, bundle, "uncertainty") == (2.0, 3.0)


def test_resolve_figsize_falls_back_to_active_default_for_leaf_style():
    leaf = HydroSectionStyle()
    assert resolve_figsize(None, leaf, "section") == (
        PYCSAMT_INTERP.default.figsize_section
    )
