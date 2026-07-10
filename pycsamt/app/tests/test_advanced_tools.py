# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for AdvancedController, TopoPreviewController, ConversionController
and all advanced-plot catalogue entries.

Real data
---------
data/AMT/TIPPER/        — 2 TIPPER EDIs (HBH03_IMP, ZBR02_IMP)
data/AMT/WILLY_DATA/    — 128 WILLY AMT EDIs across 5 profiles
                          (L18PLT 28, L22PLT 25, L26PLT 25, L30PLT 25, L34PLT 25)

Strategy
--------
* Every catalogue function is exercised against real EDIs.
* AdvancedController.draw() must never *raise*; it writes an error annotation
  in the figure when the underlying function fails or requires extra args.
* No-data paths must never raise.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pytest

# ── Paths ─────────────────────────────────────────────────────────────────────

_ROOT   = Path(__file__).parents[3]          # pycsamt/
_TIPPER = _ROOT / "data" / "AMT" / "TIPPER"
_WILLY  = _ROOT / "data" / "AMT" / "WILLY_DATA"

_HAS_TIPPER = _TIPPER.exists() and any(_TIPPER.glob("*.edi"))
_HAS_WILLY  = _WILLY.exists()  and any(_WILLY.rglob("*.edi"))

# ── Session-scoped fixtures ───────────────────────────────────────────────────

@pytest.fixture(scope="session")
def tipper_sites():
    """2-station TIPPER Sites loaded once for the whole session."""
    pytest.importorskip("pycsamt.emtools")
    if not _HAS_TIPPER:
        pytest.skip("TIPPER data not available")
    from pycsamt.emtools import ensure_sites
    return ensure_sites(str(_TIPPER))


@pytest.fixture(scope="session")
def willy_sites():
    """128-station WILLY Sites (all 5 profiles) loaded once for the whole session."""
    pytest.importorskip("pycsamt.emtools")
    if not _HAS_WILLY:
        pytest.skip("WILLY data not available")
    from pycsamt.emtools import ensure_sites
    return ensure_sites(str(_WILLY))


@pytest.fixture(scope="session")
def adv_ctrl():
    from pycsamt.app.desktop.controllers.advanced_controller import (
        AdvancedController,
    )
    return AdvancedController()


@pytest.fixture(scope="session")
def topo_ctrl():
    from pycsamt.app.desktop.controllers.advanced_controller import (
        TopoPreviewController,
    )
    return TopoPreviewController()


@pytest.fixture(scope="session")
def conv_ctrl():
    from pycsamt.app.desktop.controllers.advanced_controller import (
        ConversionController,
    )
    return ConversionController()


# ── Helpers ───────────────────────────────────────────────────────────────────

def _fig():
    """Return a fresh matplotlib Figure and close it after use."""
    return plt.figure(figsize=(6, 4))


def _close():
    plt.close("all")


# ── AdvancedController construction ──────────────────────────────────────────

class TestAdvancedControllerConstruction:
    def test_creates(self):
        from pycsamt.app.desktop.controllers.advanced_controller import (
            AdvancedController,
        )
        ctrl = AdvancedController()
        assert ctrl is not None

    def test_default_dark(self):
        from pycsamt.app.desktop.controllers.advanced_controller import (
            AdvancedController,
        )
        assert AdvancedController().dark is True

    def test_default_sites_none(self):
        from pycsamt.app.desktop.controllers.advanced_controller import (
            AdvancedController,
        )
        assert AdvancedController()._sites is None

    def test_set_sites(self):
        from pycsamt.app.desktop.controllers.advanced_controller import (
            AdvancedController,
        )
        ctrl = AdvancedController()
        ctrl.set_sites("dummy")
        assert ctrl._sites == "dummy"

    def test_clear(self):
        from pycsamt.app.desktop.controllers.advanced_controller import (
            AdvancedController,
        )
        ctrl = AdvancedController()
        ctrl.set_sites("dummy")
        ctrl.clear()
        assert ctrl._sites is None


# ── draw() — no-data paths ────────────────────────────────────────────────────

class TestDrawNoData:
    """draw() must never raise; it annotates the figure when data is missing."""

    def test_no_sites_does_not_raise(self, adv_ctrl):
        fig = _fig()
        adv_ctrl.draw("plot_strike_rose", False, fig)
        _close()

    def test_no_sites_has_ax_does_not_raise(self, adv_ctrl):
        fig = _fig()
        adv_ctrl.draw("plot_induction_arrows", True, fig)
        _close()

    def test_unknown_function_does_not_raise(self, adv_ctrl):
        fig = _fig()
        adv_ctrl.draw("not_a_real_function", False, fig)
        _close()

    def test_unknown_function_annotates_figure(self):
        # Must have sites set so the controller reaches the "function not found" branch
        from pycsamt.app.desktop.controllers.advanced_controller import (
            AdvancedController,
        )
        ctrl = AdvancedController()
        ctrl.set_sites(object())   # truthy mock — any non-None value is enough
        fig = _fig()
        ctrl.draw("no_such_fn", False, fig)
        texts = [t.get_text() for ax in fig.axes for t in ax.texts]
        assert any("not" in t.lower() or "function" in t.lower() for t in texts)
        _close()

    def test_no_sites_shows_load_message(self, adv_ctrl):
        ctrl_fresh = adv_ctrl.__class__()
        fig = _fig()
        ctrl_fresh.draw("plot_strike_rose", False, fig)
        texts = [t.get_text() for ax in fig.axes for t in ax.texts]
        assert any("load" in t.lower() or "data" in t.lower() for t in texts)
        _close()


# ── Catalogue parametrisation ─────────────────────────────────────────────────
#
# Each catalogue group is split into its own parametrised class so that
# failures are clearly attributed to a group and don't block the others.

def _catalogue_ids(group):
    return [fn for _, fn, _ in group]


def _catalogue_params(group):
    return [(fn, has_ax) for _, fn, has_ax in group]


# Strike ──────────────────────────────────────────────────────────────────────

from pycsamt.app.desktop.controllers.advanced_controller import (
    DEPTH_PLOTS,
    IMPEDANCE_PLOTS,
    INDUCTION_PLOTS,
    PHASE_TENSOR_PLOTS,
    STRIKE_PLOTS,
    SURVEY_PLOTS,
)


class TestStrikePlots:
    @pytest.mark.parametrize("fn_name,has_ax", _catalogue_params(STRIKE_PLOTS),
                              ids=_catalogue_ids(STRIKE_PLOTS))
    def test_draw_does_not_raise(self, fn_name, has_ax, adv_ctrl, willy_sites):
        adv_ctrl.set_sites(willy_sites)
        fig = _fig()
        adv_ctrl.draw(fn_name, has_ax, fig)
        _close()

    @pytest.mark.parametrize("fn_name,has_ax", _catalogue_params(STRIKE_PLOTS),
                              ids=_catalogue_ids(STRIKE_PLOTS))
    def test_draw_produces_axes(self, fn_name, has_ax, adv_ctrl, willy_sites):
        adv_ctrl.set_sites(willy_sites)
        fig = _fig()
        adv_ctrl.draw(fn_name, has_ax, fig)
        assert len(fig.axes) >= 1
        _close()


# Phase Tensor ────────────────────────────────────────────────────────────────

class TestPhaseTensorPlots:
    @pytest.mark.parametrize("fn_name,has_ax", _catalogue_params(PHASE_TENSOR_PLOTS),
                              ids=_catalogue_ids(PHASE_TENSOR_PLOTS))
    def test_draw_does_not_raise(self, fn_name, has_ax, adv_ctrl, willy_sites):
        adv_ctrl.set_sites(willy_sites)
        fig = _fig()
        adv_ctrl.draw(fn_name, has_ax, fig)
        _close()

    @pytest.mark.parametrize("fn_name,has_ax", _catalogue_params(PHASE_TENSOR_PLOTS),
                              ids=_catalogue_ids(PHASE_TENSOR_PLOTS))
    def test_draw_produces_axes(self, fn_name, has_ax, adv_ctrl, willy_sites):
        adv_ctrl.set_sites(willy_sites)
        fig = _fig()
        adv_ctrl.draw(fn_name, has_ax, fig)
        assert len(fig.axes) >= 1
        _close()


# Induction / Tipper (use TIPPER data — has actual tipper channels) ────────────

class TestInductionPlots:
    @pytest.mark.parametrize("fn_name,has_ax", _catalogue_params(INDUCTION_PLOTS),
                              ids=_catalogue_ids(INDUCTION_PLOTS))
    def test_draw_does_not_raise(self, fn_name, has_ax, adv_ctrl, tipper_sites):
        adv_ctrl.set_sites(tipper_sites)
        fig = _fig()
        adv_ctrl.draw(fn_name, has_ax, fig)
        _close()

    @pytest.mark.parametrize("fn_name,has_ax", _catalogue_params(INDUCTION_PLOTS),
                              ids=_catalogue_ids(INDUCTION_PLOTS))
    def test_draw_produces_axes(self, fn_name, has_ax, adv_ctrl, tipper_sites):
        adv_ctrl.set_sites(tipper_sites)
        fig = _fig()
        adv_ctrl.draw(fn_name, has_ax, fig)
        assert len(fig.axes) >= 1
        _close()


# Impedance / Z ───────────────────────────────────────────────────────────────

class TestImpedancePlots:
    @pytest.mark.parametrize("fn_name,has_ax", _catalogue_params(IMPEDANCE_PLOTS),
                              ids=_catalogue_ids(IMPEDANCE_PLOTS))
    def test_draw_does_not_raise(self, fn_name, has_ax, adv_ctrl, willy_sites):
        adv_ctrl.set_sites(willy_sites)
        fig = _fig()
        adv_ctrl.draw(fn_name, has_ax, fig)
        _close()

    @pytest.mark.parametrize("fn_name,has_ax", _catalogue_params(IMPEDANCE_PLOTS),
                              ids=_catalogue_ids(IMPEDANCE_PLOTS))
    def test_draw_produces_axes(self, fn_name, has_ax, adv_ctrl, willy_sites):
        adv_ctrl.set_sites(willy_sites)
        fig = _fig()
        adv_ctrl.draw(fn_name, has_ax, fig)
        assert len(fig.axes) >= 1
        _close()


# Depth Imaging ───────────────────────────────────────────────────────────────

class TestDepthPlots:
    @pytest.mark.parametrize("fn_name,has_ax", _catalogue_params(DEPTH_PLOTS),
                              ids=_catalogue_ids(DEPTH_PLOTS))
    def test_draw_does_not_raise(self, fn_name, has_ax, adv_ctrl, willy_sites):
        adv_ctrl.set_sites(willy_sites)
        fig = _fig()
        adv_ctrl.draw(fn_name, has_ax, fig)   # must never raise
        _close()

    @pytest.mark.parametrize("fn_name,has_ax", _catalogue_params(DEPTH_PLOTS),
                              ids=_catalogue_ids(DEPTH_PLOTS))
    def test_draw_produces_axes(self, fn_name, has_ax, adv_ctrl, willy_sites):
        adv_ctrl.set_sites(willy_sites)
        fig = _fig()
        adv_ctrl.draw(fn_name, has_ax, fig)
        assert len(fig.axes) >= 1
        _close()


# Survey Tools ────────────────────────────────────────────────────────────────

class TestSurveyPlots:
    @pytest.mark.parametrize("fn_name,has_ax", _catalogue_params(SURVEY_PLOTS),
                              ids=_catalogue_ids(SURVEY_PLOTS))
    def test_draw_does_not_raise(self, fn_name, has_ax, adv_ctrl, willy_sites):
        adv_ctrl.set_sites(willy_sites)
        fig = _fig()
        adv_ctrl.draw(fn_name, has_ax, fig)
        _close()

    @pytest.mark.parametrize("fn_name,has_ax", _catalogue_params(SURVEY_PLOTS),
                              ids=_catalogue_ids(SURVEY_PLOTS))
    def test_draw_produces_axes(self, fn_name, has_ax, adv_ctrl, willy_sites):
        adv_ctrl.set_sites(willy_sites)
        fig = _fig()
        adv_ctrl.draw(fn_name, has_ax, fig)
        assert len(fig.axes) >= 1
        _close()


# ── Direct emtools function tests (real data, no controller) ──────────────────
#
# These call the emtools functions directly with real data so that any
# signature mismatch surfaces immediately rather than being silently caught
# by AdvancedController.draw()'s except block.

class TestEmtoolsDirectStrike:
    def test_plot_strike_rose(self, willy_sites):
        import pycsamt.emtools as et
        et.plot_strike_rose(willy_sites, verbose=0)
        _close()

    def test_plot_strike_rose_by_line(self, willy_sites):
        import pycsamt.emtools as et
        et.plot_strike_rose_by_line(willy_sites, verbose=0)
        _close()

    def test_plot_strike_analysis(self, willy_sites):
        import pycsamt.emtools as et
        et.plot_strike_analysis(willy_sites, verbose=0)
        _close()

    def test_plot_theta_rose_grid(self, willy_sites):
        import pycsamt.emtools as et
        et.plot_theta_rose_grid(willy_sites, verbose=0)
        _close()

    def test_plot_theta_vs_period(self, willy_sites):
        import pycsamt.emtools as et
        fig, ax = plt.subplots()
        et.plot_theta_vs_period(willy_sites, ax=ax, verbose=0)
        _close()

    def test_plot_theta_stability_stripe(self, willy_sites):
        import pycsamt.emtools as et
        fig, ax = plt.subplots()
        et.plot_theta_stability_stripe(willy_sites, ax=ax, verbose=0)
        _close()

    def test_plot_strike_ribbon(self, willy_sites):
        import pycsamt.emtools as et
        fig, ax = plt.subplots()
        et.plot_strike_ribbon(willy_sites, ax=ax, verbose=0)
        _close()

    def test_plot_strike_mapsticks(self, willy_sites):
        import pycsamt.emtools as et
        fig, ax = plt.subplots()
        et.plot_strike_mapsticks(willy_sites, ax=ax, verbose=0)
        _close()

    def test_plot_strike_stability_bands(self, willy_sites):
        import pycsamt.emtools as et
        et.plot_strike_stability_bands(willy_sites, verbose=0)
        _close()


class TestEmtoolsDirectPhaseTensor:
    def test_plot_phase_tensor_map(self, willy_sites):
        import pycsamt.emtools as et
        fig, ax = plt.subplots()
        et.plot_phase_tensor_map(willy_sites, ax=ax, verbose=0)
        _close()

    def test_plot_phase_tensor_rose(self, willy_sites):
        import pycsamt.emtools as et
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="polar")
        et.plot_phase_tensor_rose(willy_sites, ax=ax, verbose=0)
        _close()

    def test_plot_phase_tensor_summary(self, willy_sites):
        import pycsamt.emtools as et
        et.plot_phase_tensor_summary(willy_sites, verbose=0)
        _close()

    def test_plot_phasor_wheel(self, willy_sites):
        import pycsamt.emtools as et
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="polar")
        et.plot_phasor_wheel(willy_sites, ax=ax, verbose=0)
        _close()

    def test_plot_pt_period_clock(self, willy_sites):
        import pycsamt.emtools as et
        et.plot_pt_period_clock(willy_sites, verbose=0)
        _close()

    def test_plot_skew_ellipt_density(self, willy_sites):
        import pycsamt.emtools as et
        fig, ax = plt.subplots()
        et.plot_skew_ellipt_density(willy_sites, ax=ax, verbose=0)
        _close()

    def test_plot_dimensionality_ternary(self, willy_sites):
        import pycsamt.emtools as et
        et.plot_dimensionality_ternary(willy_sites, verbose=0)
        _close()


class TestEmtoolsDirectInduction:
    def test_plot_induction_arrows(self, tipper_sites):
        import pycsamt.emtools as et
        fig, ax = plt.subplots()
        et.plot_induction_arrows(tipper_sites, ax=ax, verbose=0)
        _close()

    def test_plot_induction_map(self, tipper_sites):
        import pycsamt.emtools as et
        fig, ax = plt.subplots()
        et.plot_induction_map(tipper_sites, ax=ax, verbose=0)
        _close()

    def test_plot_induction_section(self, tipper_sites):
        import pycsamt.emtools as et
        fig, ax = plt.subplots()
        et.plot_induction_section(tipper_sites, ax=ax, verbose=0)
        _close()

    def test_plot_induction_rose(self, tipper_sites):
        import pycsamt.emtools as et
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="polar")
        et.plot_induction_rose(tipper_sites, ax=ax, verbose=0)
        _close()

    def test_plot_tipper_polar(self, tipper_sites):
        import pycsamt.emtools as et
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="polar")
        et.plot_tipper_polar(tipper_sites, ax=ax, verbose=0)
        _close()

    def test_plot_tipper_hodograms(self, tipper_sites):
        import pycsamt.emtools as et
        et.plot_tipper_hodograms(tipper_sites, verbose=0)
        _close()

    def test_plot_induction_multiperiod_map(self, tipper_sites):
        import pycsamt.emtools as et
        et.plot_induction_multiperiod_map(tipper_sites, verbose=0)
        _close()

    def test_plot_induction_convention(self, tipper_sites):
        import pycsamt.emtools as et
        et.plot_induction_convention(tipper_sites, verbose=0)
        _close()


class TestEmtoolsDirectImpedance:
    def test_plot_impedance_mohr_circles(self, willy_sites):
        import pycsamt.emtools as et
        et.plot_impedance_mohr_circles(willy_sites, verbose=0)
        _close()

    def test_plot_zt_argand(self, willy_sites):
        import pycsamt.emtools as et
        et.plot_zt_argand(willy_sites, verbose=0)
        _close()

    def test_plot_z_invariants_section(self, willy_sites):
        import pycsamt.emtools as et
        et.plot_z_invariants_section(willy_sites, verbose=0)
        _close()

    def test_plot_rho_phase_bode(self, willy_sites):
        import pycsamt.emtools as et
        et.plot_rho_phase_bode(willy_sites, verbose=0)
        _close()

    def test_plot_apparent_resistivity_polar(self, willy_sites):
        import pycsamt.emtools as et
        et.plot_apparent_resistivity_polar(willy_sites, verbose=0)
        _close()

    def test_plot_xyyx_crossover_map(self, willy_sites):
        import pycsamt.emtools as et
        fig, ax = plt.subplots()
        et.plot_xyyx_crossover_map(willy_sites, ax=ax, verbose=0)
        _close()

    def test_plot_offdiag_antisym_residual(self, willy_sites):
        import pycsamt.emtools as et
        fig, ax = plt.subplots()
        et.plot_offdiag_antisym_residual(willy_sites, ax=ax, verbose=0)
        _close()

    def test_plot_anisotropy(self, willy_sites):
        import pycsamt.emtools as et
        fig, ax = plt.subplots()
        et.plot_anisotropy(willy_sites, ax=ax, verbose=0)
        _close()

    def test_plot_ellipticity_psection(self, willy_sites):
        import pycsamt.emtools as et
        fig, ax = plt.subplots()
        et.plot_ellipticity_psection(willy_sites, ax=ax, verbose=0)
        _close()


class TestEmtoolsDirectDepth:
    def test_plot_depth_section(self, willy_sites):
        import pycsamt.emtools as et
        fig, ax = plt.subplots()
        et.plot_depth_section(willy_sites, ax=ax, verbose=0)
        _close()

    def test_plot_apparent_depth_psection(self, willy_sites):
        import pycsamt.emtools as et
        fig, ax = plt.subplots()
        et.plot_apparent_depth_psection(willy_sites, ax=ax, verbose=0)
        _close()

    def test_plot_gradient_section(self, willy_sites):
        import pycsamt.emtools as et
        fig, ax = plt.subplots()
        et.plot_gradient_section(willy_sites, ax=ax, verbose=0)
        _close()

    def test_plot_sensitivity_depth_section(self, willy_sites):
        import pycsamt.emtools as et
        et.plot_sensitivity_depth_section(willy_sites, verbose=0)
        _close()

    def test_plot_mt_composite_section(self, willy_sites):
        import pycsamt.emtools as et
        et.plot_mt_composite_section(willy_sites, verbose=0)
        _close()

    def test_plot_atom_psection_requires_model(self, willy_sites):
        """Verify the function signature: model is a required positional arg."""
        import inspect

        import pycsamt.emtools as et
        sig = inspect.signature(et.plot_atom_psection)
        params = list(sig.parameters.keys())
        assert "model" in params
        assert params.index("model") == 1   # second positional after sites

    def test_plot_atom_psection_no_model_shows_instruction(self, willy_sites):
        """Without a trained model, draw() shows an instruction, not a traceback."""
        from pycsamt.app.desktop.controllers.advanced_controller import (
            AdvancedController,
        )
        ctrl = AdvancedController()
        ctrl.set_sites(willy_sites)
        fig = _fig()
        ctrl.draw("plot_atom_psection", True, fig)
        texts = " ".join(t.get_text() for ax in fig.axes for t in ax.texts).lower()
        assert "model" in texts or "train" in texts or "trained" in texts
        _close()

    def test_train_dim_model_returns_dict(self, willy_sites):
        """train_dim_model must return a dict with keys D, mu, sd, feat."""
        from pycsamt.app.desktop.controllers.advanced_controller import (
            AdvancedController,
        )
        ctrl = AdvancedController()
        ctrl.set_sites(willy_sites)
        model = ctrl.train_dim_model(n_atoms=4, n_iter=10)
        assert isinstance(model, dict)
        for key in ("D", "mu", "sd", "feat"):
            assert key in model, f"model missing key '{key}'"

    def test_train_dim_model_no_sites_raises(self):
        """train_dim_model raises ValueError when no sites are loaded."""
        from pycsamt.app.desktop.controllers.advanced_controller import (
            AdvancedController,
        )
        ctrl = AdvancedController()
        with pytest.raises(ValueError, match="No survey data"):
            ctrl.train_dim_model()

    def test_plot_atom_psection_after_training_renders(self, willy_sites):
        """After training the model, draw() produces real axes (not error text)."""
        from pycsamt.app.desktop.controllers.advanced_controller import (
            AdvancedController,
        )
        ctrl = AdvancedController()
        ctrl.set_sites(willy_sites)
        ctrl.train_dim_model(n_atoms=4, n_iter=10)
        fig = _fig()
        ctrl.draw("plot_atom_psection", True, fig)
        assert len(fig.axes) >= 1
        # Must not be a pure error annotation
        texts = " ".join(t.get_text() for ax in fig.axes for t in ax.texts).lower()
        assert "train" not in texts and "no dictionary" not in texts
        _close()

    def test_clear_resets_dim_model(self, willy_sites):
        """clear() must reset the stored dim model."""
        from pycsamt.app.desktop.controllers.advanced_controller import (
            AdvancedController,
        )
        ctrl = AdvancedController()
        ctrl.set_sites(willy_sites)
        ctrl.train_dim_model(n_atoms=4, n_iter=10)
        assert ctrl._dim_model is not None
        ctrl.clear()
        assert ctrl._dim_model is None


class TestEmtoolsDirectSurvey:
    def test_plot_survey_fingerprint(self, willy_sites):
        import pycsamt.emtools as et
        et.plot_survey_fingerprint(willy_sites, verbose=0)
        _close()

    def test_plot_distortion_radar(self, willy_sites):
        import pycsamt.emtools as et
        et.plot_distortion_radar(willy_sites, verbose=0)
        _close()

    def test_plot_sites_compare(self, willy_sites):
        import pycsamt.emtools as et
        et.plot_sites_compare(willy_sites, verbose=0)
        _close()

    def test_plot_sites_panels(self, willy_sites):
        """Regression: was broken by return_errors / wrong unpack (fixed)."""
        import pycsamt.emtools as et
        et.plot_sites_panels(willy_sites, verbose=0)
        _close()

    def test_plot_normalized_response(self, willy_sites):
        import pycsamt.emtools as et
        et.plot_normalized_response(willy_sites, verbose=0)
        _close()

    def test_plot_tf_coherence_network(self, willy_sites):
        import pycsamt.emtools as et
        et.plot_tf_coherence_network(willy_sites, verbose=0)
        _close()

    def test_plot_dim_map(self, willy_sites):
        import pycsamt.emtools as et
        fig, ax = plt.subplots()
        et.plot_dim_map(willy_sites, ax=ax, verbose=0)
        _close()

    def test_plot_dim_occupancy_area(self, willy_sites):
        import pycsamt.emtools as et
        fig, ax = plt.subplots()
        et.plot_dim_occupancy_area(willy_sites, ax=ax, verbose=0)
        _close()


# ── TopoPreviewController ─────────────────────────────────────────────────────

class TestTopoPreviewController:
    def test_creates(self):
        from pycsamt.app.desktop.controllers.advanced_controller import (
            TopoPreviewController,
        )
        ctrl = TopoPreviewController()
        assert ctrl is not None

    def test_default_dark(self):
        from pycsamt.app.desktop.controllers.advanced_controller import (
            TopoPreviewController,
        )
        assert TopoPreviewController().dark is True

    def test_elevation_profile_no_data_does_not_raise(self, topo_ctrl):
        fig = _fig()
        topo_ctrl.plot_elevation_profile(fig)
        _close()

    def test_fill_preview_no_data_does_not_raise(self, topo_ctrl):
        fig = _fig()
        topo_ctrl.plot_fill_preview(fig)
        _close()

    def test_elevation_histogram_no_data_does_not_raise(self, topo_ctrl):
        fig = _fig()
        topo_ctrl.plot_elevation_histogram(fig)
        _close()

    def test_get_stats_no_data_returns_dict(self, topo_ctrl):
        stats = topo_ctrl.get_stats()
        assert isinstance(stats, dict)
        assert "n_stations" in stats
        assert stats["n_stations"] == 0

    def test_fill_preview_with_sites_does_not_raise(self, topo_ctrl, willy_sites):
        topo_ctrl.set_sites(willy_sites)
        fig = _fig()
        topo_ctrl.plot_fill_preview(fig)
        _close()

    def test_elevation_histogram_with_sites_does_not_raise(self, topo_ctrl, willy_sites):
        topo_ctrl.set_sites(willy_sites)
        fig = _fig()
        topo_ctrl.plot_elevation_histogram(fig)
        _close()

    def test_get_stats_with_sites_returns_n_stations(self, topo_ctrl, willy_sites):
        topo_ctrl.set_sites(willy_sites)
        stats = topo_ctrl.get_stats()
        assert stats["n_stations"] > 0

    def test_fill_preview_creates_axes(self, topo_ctrl):
        """fill_preview uses synthetic data and must always produce a figure."""
        fig = _fig()
        topo_ctrl.plot_fill_preview(fig)
        assert len(fig.axes) >= 1
        _close()

    def test_style_both_modes(self):
        from pycsamt.app.desktop.controllers.advanced_controller import (
            TopoPreviewController,
        )
        ctrl = TopoPreviewController()
        for dark in (True, False):
            ctrl.dark = dark
            style = ctrl._get_style()
            assert "bg" in style and "fg" in style


# ── ConversionController ──────────────────────────────────────────────────────

class TestConversionController:
    def test_creates(self, conv_ctrl):
        assert conv_ctrl is not None

    def test_has_result_false_initially(self, conv_ctrl):
        assert conv_ctrl.has_result is False

    def test_result_none_initially(self, conv_ctrl):
        assert conv_ctrl.result is None

    def test_set_source(self, conv_ctrl):
        conv_ctrl.set_source("AVG -> EDI", "/tmp/test.avg")
        assert conv_ctrl._source_type == "AVG -> EDI"
        assert conv_ctrl._source_path == "/tmp/test.avg"

    def test_build_stats_none_collection(self, conv_ctrl):
        stats = conv_ctrl.build_stats(None, [])
        assert stats["n_total"] == 0
        assert stats["n_failures"] == 0

    def test_build_stats_with_failures(self, conv_ctrl):
        stats = conv_ctrl.build_stats(None, ["fail1", "fail2"])
        assert stats["n_failures"] == 2

    def test_plot_impedance_curves_no_data_does_not_raise(self, conv_ctrl):
        fig = _fig()
        conv_ctrl.plot_impedance_curves(fig)
        _close()

    def test_plot_station_map_no_data_does_not_raise(self, conv_ctrl):
        fig = _fig()
        conv_ctrl.plot_station_map(fig)
        _close()

    def test_plot_impedance_curves_no_data_has_axes(self, conv_ctrl):
        fig = _fig()
        conv_ctrl.plot_impedance_curves(fig)
        assert len(fig.axes) >= 1
        _close()

    def test_plot_station_map_no_data_has_axes(self, conv_ctrl):
        fig = _fig()
        conv_ctrl.plot_station_map(fig)
        assert len(fig.axes) >= 1
        _close()

    def test_invalid_source_raises(self, conv_ctrl):
        conv_ctrl.set_source("UNKNOWN -> EDI", "/tmp/bad.xyz")
        with pytest.raises(ValueError, match="Unknown source type"):
            conv_ctrl.run({})

    def test_style_dark_and_light(self, conv_ctrl):
        for dark in (True, False):
            conv_ctrl.dark = dark
            style = conv_ctrl._get_style()
            assert "bg" in style and "fg" in style

    def test_avg_run_uses_config_freq_order_not_transform_kw(
        self,
        conv_ctrl,
        monkeypatch,
    ):
        calls = {}

        class FakeAVG:
            def transform(self, path, **kw):
                from pycsamt.core.config import get_config
                calls["path"] = path
                calls["kw"] = kw
                calls["freq_order"] = get_config().freq_order
                return ["edi"]

        import pycsamt.transformers as transformers
        monkeypatch.setattr(transformers, "AVGtoEDI", lambda: FakeAVG())

        conv_ctrl.set_source("AVG -> EDI", "/tmp/in.avg")
        collection, failures = conv_ctrl.run({"freq_order": "ascending"})

        assert collection == ["edi"]
        assert failures == []
        assert calls["path"] == "/tmp/in.avg"
        assert calls["kw"] == {}
        assert calls["freq_order"] == "asc"

    def test_avg_run_attaches_station_profile(self, conv_ctrl, monkeypatch):
        calls = {}

        class FakeTopo:
            def convert_coords(self, *, to, inplace):
                calls["convert"] = (to, inplace)

        class FakeAVGObject:
            topo = FakeTopo()

            def add_topography(self, stn_file, *, utm_zone=None, epsg=None):
                calls["topography"] = (stn_file, utm_zone, epsg)
                return self

        class FakeAVG:
            @classmethod
            def from_file(cls, path):
                calls["avg_path"] = path
                return FakeAVGObject()

        class FakeAVGtoEDI:
            def transform(self, source, **kw):
                from pycsamt.core.config import get_config
                calls["source"] = source
                calls["kw"] = kw
                calls["freq_tol"] = get_config().freq_tol
                calls["compute_z"] = get_config().compute_z_from_res
                calls["compute_rho_phi"] = get_config().compute_res_from_z
                return ["edi"]

        import pycsamt.transformers as transformers
        import pycsamt.zonge.avg as avg_mod
        monkeypatch.setattr(transformers, "AVGtoEDI", lambda: FakeAVGtoEDI())
        monkeypatch.setattr(avg_mod, "AVG", FakeAVG)

        conv_ctrl.set_source("AVG -> EDI", "/tmp/in.avg")
        collection, failures = conv_ctrl.run({
            "freq_order": "ascending",
            "freq_tol": 1e-6,
            "compute_z": True,
            "compute_rho_phi": True,
            "stn_path": "/tmp/K1.stn",
            "utm_zone": "50N",
            "epsg": "32650",
            "convert_stn_coords": True,
            "name": "K1",
        })

        assert collection == ["edi"]
        assert failures == []
        assert calls["avg_path"] == "/tmp/in.avg"
        assert calls["topography"] == ("/tmp/K1.stn", "50N", 32650)
        assert calls["convert"] == ("ll", True)
        assert isinstance(calls["source"], FakeAVGObject)
        assert calls["kw"] == {"name": "K1"}
        assert calls["freq_tol"] == 1e-6
        assert calls["compute_z"] is True
        assert calls["compute_rho_phi"] is True

    def test_avg_run_with_real_stn_utm_populates_edi_coordinates(self, conv_ctrl):
        avg_path = _ROOT / "data" / "avg" / "K1.AVG"
        stn_path = _ROOT / "data" / "avg" / "K1.stn"
        if not avg_path.exists() or not stn_path.exists():
            pytest.skip("K1 AVG/STN fixtures are not available")

        conv_ctrl.set_source("AVG -> EDI", str(avg_path))
        collection, failures = conv_ctrl.run({
            "freq_order": "descending",
            "stn_path": str(stn_path),
            "utm_zone": "49N",
            "convert_stn_coords": True,
        })

        assert failures == []
        assert len(collection) >= 1

        first = next(iter(collection))
        head = first.get_section("head")
        definemeas = first.get_section("definemeas")
        assert head is not None
        assert definemeas is not None
        assert head.lat == pytest.approx(26.052401, abs=1e-6)
        assert head.long == pytest.approx(113.487159, abs=1e-6)
        assert head.elev == pytest.approx(574.5, abs=1e-6)
        assert definemeas.reflat == pytest.approx(head.lat, abs=1e-6)
        assert definemeas.reflong == pytest.approx(head.long, abs=1e-6)
        assert definemeas.refelev == pytest.approx(head.elev, abs=1e-6)

        stats = conv_ctrl.build_stats(collection, failures)
        first_row = stats["rows"][0]
        assert first_row["lat"] == pytest.approx(head.lat, abs=1e-6)
        assert first_row["lon"] == pytest.approx(head.long, abs=1e-6)
        assert first_row["elev"] == pytest.approx(head.elev, abs=1e-6)

    def test_j_run_drops_ui_only_options(self, conv_ctrl, monkeypatch):
        calls = {}

        class FakeJ:
            def transform(self, path, **kw):
                calls["path"] = path
                calls["kw"] = kw
                return ["edi"]

        import pycsamt.transformers as transformers
        monkeypatch.setattr(transformers, "JtoEDI", lambda: FakeJ())

        conv_ctrl.set_source("J -> EDI", "/tmp/in.j")
        collection, failures = conv_ctrl.run({
            "freq_order": "descending",
            "station_suffix": "_IMP",
        })

        assert collection == ["edi"]
        assert failures == []
        assert calls["path"] == "/tmp/in.j"
        assert calls["kw"] == {}

    def test_spectra_run_maps_ui_options(self, conv_ctrl, monkeypatch):
        calls = {}

        class FakeResult:
            collection = ["edi"]
            failures = ["bad"]

        class FakeSpectra:
            def __init__(self, **kw):
                calls["init"] = kw

            def transform_batch(self, path, **kw):
                calls["path"] = path
                calls["transform"] = kw
                return FakeResult()

        import pycsamt.transformers as transformers
        monkeypatch.setattr(transformers, "SpectraToEDI", FakeSpectra)

        conv_ctrl.set_source("Spectra -> EDI", "/tmp/spec")
        collection, failures = conv_ctrl.run({
            "e_labels": "EX,EY",
            "h_labels": "HX,HY",
            "estimate_errors": True,
            "use_remote_ref": True,
            "station_suffix": "_IMP",
            "skip_errors": False,
            "output_dir": "/tmp/out",
        })

        assert collection == ["edi"]
        assert failures == ["bad"]
        assert calls["path"] == "/tmp/spec"
        assert calls["init"] == {
            "e_labels": ("EX", "EY"),
            "h_labels": ("HX", "HY"),
            "estimate_error": True,
            "use_remote": True,
            "skip_errors": False,
            "station_suffix": "_IMP",
        }
        assert calls["transform"] == {"output_dir": "/tmp/out"}


# ── AdvancedController dark/light mode ───────────────────────────────────────

class TestAdvancedControllerDarkMode:
    def test_dark_mode_toggle(self, adv_ctrl, willy_sites):
        adv_ctrl.set_sites(willy_sites)
        for dark in (True, False):
            adv_ctrl.dark = dark
            fig = _fig()
            adv_ctrl.draw("plot_strike_rose", False, fig)
            _close()

    def test_draw_light_mode_does_not_raise(self, adv_ctrl, willy_sites):
        adv_ctrl.set_sites(willy_sites)
        adv_ctrl.dark = False
        fig = _fig()
        adv_ctrl.draw("plot_phase_tensor_map", True, fig)
        adv_ctrl.dark = True
        _close()


# ── Regression tests for fixed bugs ──────────────────────────────────────────

class TestRegressions:
    def test_get_z_block_wrong_kwarg_fixed(self):
        """plot_sites_panels was calling _get_z_block(return_errors=True) instead
        of _get_z_block(with_errors=True). Verify the correct kwarg now works."""
        import logging; logging.disable(logging.CRITICAL)
        from pycsamt.emtools import ensure_sites
        from pycsamt.emtools._core import _get_z_block
        if not _HAS_WILLY:
            pytest.skip("WILLY data not available")
        sites = ensure_sites(str(_WILLY))
        ed = sites[0].edi
        result = _get_z_block(ed, with_errors=True)
        assert len(result) == 4

    def test_get_z_block_bad_kwarg_raises(self):
        """The old broken keyword 'return_errors' must raise TypeError."""
        from pycsamt.emtools import ensure_sites
        from pycsamt.emtools._core import _get_z_block
        if not _HAS_WILLY:
            pytest.skip("WILLY data not available")
        sites = ensure_sites(str(_WILLY))
        ed = sites[0].edi
        with pytest.raises(TypeError, match="unexpected keyword argument"):
            _get_z_block(ed, return_errors=True)

    def test_plot_sites_panels_no_longer_raises(self, willy_sites):
        """Regression: plot_sites_panels raised ValueError (too many values to unpack)."""
        import pycsamt.emtools as et
        et.plot_sites_panels(willy_sites, verbose=0)
        _close()
