# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for CorrectionController — the non-destructive correction stack
backing the desktop app's "Data Corrections" window (static-shift removal,
noise removal, frequency filtering, coordinate corrections, tensor rotation,
and the Stratagem-branded convenience pipeline).

Real data
---------
data/AMT/WILLY_DATA/L18PLT/  — 28-station AMT profile, used as the default
                                working set (fast, multi-line-capable).

Strategy
--------
* Items 1-3 of the task (dataclass/properties, undo-stack state machine,
  wrap/coord dispatch table) are exercised through the *public*
  apply()/preview() API rather than by calling the private
  ``_wrap_*``/``_coord_*`` staticmethods directly, so that ``_call_fn``'s
  dispatch table and the stack machinery are covered too.
* Item 4 (Stratagem convenience layer) is driven against a real EDI
  directory via ``load_edi_dir()``.
* Item 5 (plotting) is smoke-tested with matplotlib's Agg backend: a plot
  method must not raise and must leave the figure/axes populated.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pytest

# ── Paths ─────────────────────────────────────────────────────────────────────

_ROOT = Path(__file__).parents[3]  # pycsamt/
_L18 = _ROOT / "data" / "AMT" / "WILLY_DATA" / "L18PLT"
_HAS_L18 = _L18.exists() and any(_L18.glob("*.edi"))


# ── Session-scoped fixtures ───────────────────────────────────────────────────


@pytest.fixture(scope="session")
def l18_dir():
    if not _HAS_L18:
        pytest.skip("WILLY L18PLT data not available")
    return str(_L18)


@pytest.fixture(scope="session")
def l18_files(l18_dir):
    return sorted(str(p) for p in _L18.glob("*.edi"))


@pytest.fixture(scope="session")
def small_files(l18_files):
    """First 10 stations — the default fast working set for parametrized tests."""
    return l18_files[:10]


@pytest.fixture(scope="session")
def small_sites(small_files):
    """10-station Sites, loaded once for the whole session."""
    pytest.importorskip("pycsamt.emtools")
    from pycsamt.emtools._core import ensure_sites

    return ensure_sites(small_files)


@pytest.fixture(scope="session")
def l18_sites(l18_files):
    """Full 28-station L18PLT Sites, loaded once for the whole session."""
    from pycsamt.emtools._core import ensure_sites

    return ensure_sites(l18_files)


# ── Per-test controller ───────────────────────────────────────────────────────


@pytest.fixture
def ctrl():
    from pycsamt.app.desktop.controllers.correction_controller import (
        CorrectionController,
    )

    return CorrectionController()


def _fig():
    return plt.figure(figsize=(6, 4))


def _close():
    plt.close("all")


# ═══════════════════════════════════════════════════════════════════════════
# 1. CorrectionStep dataclass + cheap controller state/property surface
# ═══════════════════════════════════════════════════════════════════════════


class TestCorrectionStep:
    def test_construct_minimal(self):
        from pycsamt.app.desktop.controllers.correction_controller import (
            CorrectionStep,
        )

        step = CorrectionStep(
            label="Rotate 10°",
            fn_name="_wrap_rotate",
            kwargs={"angle": 10.0},
            sites_after="SITES_STUB",
        )
        assert step.label == "Rotate 10°"
        assert step.fn_name == "_wrap_rotate"
        assert step.kwargs == {"angle": 10.0}
        assert step.sites_after == "SITES_STUB"
        assert step.coords_df is None
        assert step.index == 0

    def test_construct_full(self):
        from pycsamt.app.desktop.controllers.correction_controller import (
            CorrectionStep,
        )

        step = CorrectionStep(
            label="Coord shift",
            fn_name="_coord_shift",
            kwargs={"delta_lat": 0.001},
            sites_after="SITES",
            coords_df="DF_STUB",
            index=3,
        )
        assert step.coords_df == "DF_STUB"
        assert step.index == 3

    def test_kwargs_independent_between_instances(self):
        """dataclasses.field default handling: each instance's dict must be
        independent (kwargs has no `default_factory`, but callers always pass
        a fresh dict — verify aliasing doesn't leak across two steps)."""
        from pycsamt.app.desktop.controllers.correction_controller import (
            CorrectionStep,
        )

        s1 = CorrectionStep("a", "fn", {"x": 1}, "S1")
        s2 = CorrectionStep("b", "fn", {"x": 2}, "S2")
        s1.kwargs["x"] = 99
        assert s2.kwargs["x"] == 2


class TestControllerBasicState:
    def test_fresh_controller_has_no_data(self, ctrl):
        assert ctrl.has_data is False
        assert ctrl.has_corrections is False
        assert ctrl.has_coord_corrections is False
        assert ctrl.current_sites is None
        assert ctrl.raw_sites is None
        assert ctrl.raw_coords_df is None
        assert ctrl.current_coords_df() is None
        assert ctrl.stack == []
        assert ctrl.n_steps == 0
        assert ctrl.dark is True

    def test_set_raw_sites_populates_state(self, ctrl, small_sites):
        ctrl.set_raw_sites(small_sites)
        assert ctrl.has_data is True
        assert ctrl.raw_sites is small_sites
        assert ctrl.current_sites is small_sites  # no steps yet
        assert ctrl.raw_coords_df is not None
        assert ctrl.current_coords_df() is ctrl.raw_coords_df
        assert ctrl.n_steps == 0
        assert ctrl.has_corrections is False

    def test_set_raw_sites_clears_previous_stack(self, ctrl, small_sites):
        ctrl.set_raw_sites(small_sites)
        ctrl.apply("_wrap_rotate", {"angle": 5.0}, label="rot")
        assert ctrl.n_steps == 1
        ctrl.set_raw_sites(small_sites)
        assert ctrl.n_steps == 0
        assert ctrl.has_corrections is False

    def test_clear_resets_everything(self, ctrl, small_sites):
        ctrl.set_raw_sites(small_sites)
        ctrl.apply("_wrap_rotate", {"angle": 5.0}, label="rot")
        ctrl.clear()
        assert ctrl.has_data is False
        assert ctrl.raw_sites is None
        assert ctrl.raw_coords_df is None
        assert ctrl.current_sites is None
        assert ctrl.n_steps == 0

    def test_stack_property_returns_copy(self, ctrl, small_sites):
        ctrl.set_raw_sites(small_sites)
        ctrl.apply("_wrap_rotate", {"angle": 5.0}, label="rot")
        s1 = ctrl.stack
        s1.append("intruder")
        assert len(ctrl.stack) == 1  # internal list untouched

    def test_current_sites_reflects_last_step(self, ctrl, small_sites):
        ctrl.set_raw_sites(small_sites)
        step = ctrl.apply("_wrap_rotate", {"angle": 5.0}, label="rot")
        assert ctrl.current_sites is step.sites_after
        assert ctrl.current_sites is not small_sites


# ═══════════════════════════════════════════════════════════════════════════
# 2. Undo-able correction-stack state machine
# ═══════════════════════════════════════════════════════════════════════════


class TestApplyPreviewStateMachine:
    def test_preview_no_data_returns_none(self, ctrl):
        assert ctrl.preview("_wrap_rotate", {"angle": 5.0}) is None

    def test_apply_no_data_returns_none(self, ctrl):
        assert ctrl.apply("_wrap_rotate", {"angle": 5.0}, label="x") is None
        assert ctrl.n_steps == 0

    def test_preview_does_not_touch_stack(self, ctrl, small_sites):
        ctrl.set_raw_sites(small_sites)
        result = ctrl.preview("_wrap_rotate", {"angle": 5.0})
        assert result is not None
        assert ctrl.n_steps == 0  # preview never pushes to the stack
        assert ctrl.current_sites is small_sites

    def test_preview_coord_fn_returns_dataframe(self, ctrl, small_sites):
        import pandas as pd

        ctrl.set_raw_sites(small_sites)
        df = ctrl.preview("_coord_shift", {"delta_lat": 0.001})
        assert isinstance(df, pd.DataFrame)
        assert ctrl._preview_coords_df is not None
        assert ctrl.n_steps == 0

    def test_preview_unknown_fn_raises(self, ctrl, small_sites):
        ctrl.set_raw_sites(small_sites)
        with pytest.raises(Exception):
            ctrl.preview("not_a_real_correction_fn", {})

    def test_apply_unknown_fn_returns_none(self, ctrl, small_sites):
        ctrl.set_raw_sites(small_sites)
        result = ctrl.apply("not_a_real_correction_fn", {}, label="bad")
        assert result is None
        assert ctrl.n_steps == 0

    def test_apply_pushes_step_with_correct_metadata(self, ctrl, small_sites):
        ctrl.set_raw_sites(small_sites)
        step = ctrl.apply("_wrap_rotate", {"angle": 12.5}, label="Rotate 12.5°")
        assert step is not None
        assert step.label == "Rotate 12.5°"
        assert step.fn_name == "_wrap_rotate"
        assert step.kwargs == {"angle": 12.5}
        assert step.index == 0
        assert ctrl.n_steps == 1
        assert ctrl.stack[0] is step

    def test_apply_kwargs_stored_as_independent_copy(self, ctrl, small_sites):
        ctrl.set_raw_sites(small_sites)
        kwargs = {"angle": 7.0}
        step = ctrl.apply("_wrap_rotate", kwargs, label="rot")
        kwargs["angle"] = 999.0
        assert step.kwargs["angle"] == 7.0

    def test_multiple_applies_chain_sequentially(self, ctrl, small_sites):
        ctrl.set_raw_sites(small_sites)
        ctrl.apply("_wrap_rotate", {"angle": 5.0}, label="rot5")
        step2 = ctrl.apply("_wrap_rotate", {"angle": 10.0}, label="rot10")
        assert ctrl.n_steps == 2
        assert step2.index == 1
        assert ctrl.current_sites is step2.sites_after

    def test_undo_last_pops_stack(self, ctrl, small_sites):
        ctrl.set_raw_sites(small_sites)
        ctrl.apply("_wrap_rotate", {"angle": 5.0}, label="rot5")
        ctrl.apply("_wrap_rotate", {"angle": 10.0}, label="rot10")
        assert ctrl.n_steps == 2
        ctrl.undo_last()
        assert ctrl.n_steps == 1
        assert ctrl.stack[0].label == "rot5"

    def test_undo_last_on_empty_stack_is_noop(self, ctrl, small_sites):
        ctrl.set_raw_sites(small_sites)
        ctrl.undo_last()  # must not raise
        assert ctrl.n_steps == 0

    def test_undo_last_clears_preview_coords(self, ctrl, small_sites):
        ctrl.set_raw_sites(small_sites)
        ctrl.preview("_coord_shift", {"delta_lat": 0.001})
        assert ctrl._preview_coords_df is not None
        ctrl.undo_last()
        assert ctrl._preview_coords_df is None

    def test_revert_all_clears_stack_keeps_raw(self, ctrl, small_sites):
        ctrl.set_raw_sites(small_sites)
        ctrl.apply("_wrap_rotate", {"angle": 5.0}, label="rot5")
        ctrl.apply("_wrap_rotate", {"angle": 10.0}, label="rot10")
        ctrl.revert_all()
        assert ctrl.n_steps == 0
        assert ctrl.has_data is True
        assert ctrl.current_sites is small_sites

    def test_remove_step_invalid_index_is_noop(self, ctrl, small_sites):
        ctrl.set_raw_sites(small_sites)
        ctrl.apply("_wrap_rotate", {"angle": 5.0}, label="rot5")
        ctrl.remove_step(-1)
        ctrl.remove_step(100)
        assert ctrl.n_steps == 1

    def test_remove_step_replays_subsequent_impedance_steps(self, ctrl, small_sites):
        ctrl.set_raw_sites(small_sites)
        ctrl.apply("_wrap_rotate", {"angle": 5.0}, label="rot5")
        ctrl.apply("_wrap_rotate", {"angle": 10.0}, label="rot10")
        ctrl.apply("_wrap_rotate", {"angle": 15.0}, label="rot15")
        last_before = ctrl.stack[-1].sites_after

        ctrl.remove_step(1)  # drop "rot10"

        assert ctrl.n_steps == 2
        assert [s.label for s in ctrl.stack] == ["rot5", "rot15"]
        # Replay recomputed the tail step against the new base, so its
        # sites_after object identity changed even though the label didn't.
        assert ctrl.stack[-1].sites_after is not last_before

    def test_remove_step_before_coord_step_replays_correctly(self, ctrl, small_sites):
        """Coord steps are not replayed (DataFrame-only chain); the impedance
        step after them must still be recomputed against the new base."""
        ctrl.set_raw_sites(small_sites)
        ctrl.apply("_wrap_rotate", {"angle": 5.0}, label="rot5")
        ctrl.apply("_coord_shift", {"delta_lat": 0.001}, label="shift")
        ctrl.apply("_wrap_rotate", {"angle": 20.0}, label="rot20")
        assert ctrl.n_steps == 3

        ctrl.remove_step(0)  # drop "rot5" — base becomes raw_sites

        assert ctrl.n_steps == 2
        assert [s.label for s in ctrl.stack] == ["shift", "rot20"]
        # coord step's sites_after should now equal the (new) base = raw sites
        assert ctrl.stack[0].sites_after is small_sites
        assert ctrl.stack[0].coords_df is not None

    def test_remove_step_clears_preview_coords(self, ctrl, small_sites):
        ctrl.set_raw_sites(small_sites)
        ctrl.apply("_wrap_rotate", {"angle": 5.0}, label="rot5")
        ctrl.preview("_coord_shift", {"delta_lat": 0.001})
        assert ctrl._preview_coords_df is not None
        ctrl.remove_step(0)
        assert ctrl._preview_coords_df is None


class TestCoordCorrectionStack:
    def test_apply_coord_fn_leaves_sites_unchanged(self, ctrl, small_sites):
        ctrl.set_raw_sites(small_sites)
        step = ctrl.apply("_coord_shift", {"delta_lat": 0.001}, label="shift")
        assert step is not None
        assert step.sites_after is small_sites  # Sites untouched
        assert step.coords_df is not None
        assert ctrl.has_coord_corrections is True

    def test_current_coords_df_falls_back_to_raw(self, ctrl, small_sites):
        ctrl.set_raw_sites(small_sites)
        assert ctrl.current_coords_df() is ctrl.raw_coords_df
        ctrl.apply("_wrap_rotate", {"angle": 5.0}, label="rot")  # not coord
        assert ctrl.current_coords_df() is ctrl.raw_coords_df

    def test_current_coords_df_uses_latest_coord_step(self, ctrl, small_sites):
        ctrl.set_raw_sites(small_sites)
        step1 = ctrl.apply("_coord_shift", {"delta_lat": 0.001}, label="shift1")
        step2 = ctrl.apply("_coord_shift", {"delta_lat": -0.002}, label="shift2")
        assert ctrl.current_coords_df() is step2.coords_df
        assert ctrl.current_coords_df() is not step1.coords_df

    def test_chained_coord_corrections_build_on_previous_df(self, ctrl, small_sites):
        ctrl.set_raw_sites(small_sites)
        ctrl.apply("_coord_spacing_regularize", {"spacing_m": 150.0}, label="reg")
        step2 = ctrl.apply(
            "_coord_elevation_smooth",
            {"method": "mean", "window": 3},
            label="elev",
        )
        assert step2 is not None
        assert ctrl.n_steps == 2


class TestLoadEdiDir:
    def test_load_edi_dir_returns_station_count(self, ctrl, l18_dir):
        n = ctrl.load_edi_dir(l18_dir)
        assert n == 28
        assert ctrl.has_data is True

    def test_load_edi_dir_resets_stack(self, ctrl, l18_dir, small_sites):
        ctrl.set_raw_sites(small_sites)
        ctrl.apply("_wrap_rotate", {"angle": 5.0}, label="rot")
        assert ctrl.n_steps == 1
        ctrl.load_edi_dir(l18_dir)
        assert ctrl.n_steps == 0

    def test_load_edi_dir_populates_strat_edi_objects(self, ctrl, l18_dir):
        ctrl.load_edi_dir(l18_dir)
        assert len(ctrl._strat_edi_objects) == 28
        assert len(ctrl._strat_raw_edi_objects) == 28
        assert ctrl._strat_edi_dir == l18_dir
        assert ctrl._strat_qc_report is None
        assert ctrl._strat_ss_factors is None

    def test_load_edi_dir_builds_usable_sites(self, ctrl, l18_dir):
        ctrl.load_edi_dir(l18_dir)
        # current_sites must support the normal correction path too
        step = ctrl.apply("_wrap_rotate", {"angle": 3.0}, label="rot")
        assert step is not None

    def test_edi_objects_to_sites_static_helper(self, ctrl, l18_dir):
        from pycsamt.app.desktop.controllers.correction_controller import (
            CorrectionController,
        )

        ctrl.load_edi_dir(l18_dir)
        sites = CorrectionController._edi_objects_to_sites(ctrl._strat_edi_objects)
        assert sites is not None


# ═══════════════════════════════════════════════════════════════════════════
# 3. Wrap / coord correction dispatch sweep — every entry driven through
#    the real apply() path (not the static methods directly).
# ═══════════════════════════════════════════════════════════════════════════

WRAP_COORD_CASES = [
    ("_correct_ss_loess", {"half_window": 2, "poly": 1, "it": 1}),
    ("_correct_ss_bilateral", {"half_window": 2, "max_skew": 6.0}),
    ("_correct_ss_refmedian", {"smooth_sites": 0, "max_skew": 6.0}),
    ("_correct_near_field", {"source_offset": 500.0}),
    ("_wrap_rotate", {"angle": 15.0}),
    ("_wrap_rotate_to_strike", {"method": "swift"}),
    ("_wrap_rotate_pt_strike", {"band_lo": 0.001, "band_hi": 1000.0}),
    ("_wrap_rotate_to_profile", {"azimuth": -1.0}),
    ("_coord_profile_projection", {"azimuth": -1.0, "keep_elevation": True}),
    (
        "_coord_spacing_regularize",
        {"spacing_m": 150.0, "azimuth": -1.0, "preserve_extent": True},
    ),
    ("_coord_outlier_snap", {"threshold_m": 50.0, "azimuth": -1.0}),
    ("_coord_elevation_smooth", {"method": "loess", "window": 3}),
    (
        "_coord_shift",
        {"delta_lat": 0.001, "delta_lon": -0.001, "delta_elev": 5.0},
    ),
    ("_coord_interpolate_missing", {"fill_nan_only": True}),
]


class TestWrapCoordDispatchSweep:
    @pytest.mark.parametrize(
        "fn_name,kwargs",
        WRAP_COORD_CASES,
        ids=[c[0] for c in WRAP_COORD_CASES],
    )
    def test_apply_does_not_raise_and_updates_state(self, fn_name, kwargs, small_sites):
        from pycsamt.app.desktop.controllers.correction_controller import (
            _COORD_FN_NAMES,
            CorrectionController,
        )

        c = CorrectionController()
        c.set_raw_sites(small_sites)
        step = c.apply(fn_name, dict(kwargs), label=fn_name)

        assert step is not None, f"apply() returned None for {fn_name}"
        assert c.n_steps == 1
        assert step.fn_name == fn_name

        if fn_name in _COORD_FN_NAMES:
            assert step.coords_df is not None
            assert step.sites_after is small_sites  # coord: Sites unchanged
            assert c.has_coord_corrections is True
        else:
            assert step.coords_df is None
            assert step.sites_after is not None
            assert c.current_sites is step.sites_after


# Passthrough catalogue entries — not wrapped, dispatched straight to
# emtools.<fn_name>(sites, inplace=False, verbose=0, **kwargs). Exercises
# the final fallback branch of _call_fn.
PASSTHROUGH_CASES = [
    (
        "correct_ss_ama",
        {
            "sort_by": "lon",
            "half_window": 3,
            "weights": "tri",
            "max_skew": 6.0,
        },
    ),
    (
        "notch_powerline",
        {
            "mains_hz": 50.0,
            "n_harm": 5,
            "tol_hz": 0.08,
            "mode": "interp",
            "also": "both",
        },
    ),
    ("smooth_logfreq", {"win": 3, "kind": "tri", "also": "both"}),
    ("antisymmetrize", {"how": "rms"}),
]


class TestPassthroughDispatchSweep:
    @pytest.mark.parametrize(
        "fn_name,kwargs",
        PASSTHROUGH_CASES,
        ids=[c[0] for c in PASSTHROUGH_CASES],
    )
    def test_apply_passthrough_fn(self, fn_name, kwargs, small_sites):
        from pycsamt.app.desktop.controllers.correction_controller import (
            CorrectionController,
        )

        c = CorrectionController()
        c.set_raw_sites(small_sites)
        step = c.apply(fn_name, dict(kwargs), label=fn_name)
        assert step is not None, f"apply() returned None for {fn_name}"
        assert c.current_sites is step.sites_after


class TestPreviewWrapCoordSweep:
    """A lighter pass through preview() (not apply()) for a representative
    subset — preview must compute without mutating the stack."""

    @pytest.mark.parametrize(
        "fn_name,kwargs",
        [
            ("_wrap_rotate", {"angle": 8.0}),
            ("_correct_near_field", {"source_offset": 400.0}),
            ("_coord_shift", {"delta_lat": 0.0005}),
            ("_coord_elevation_smooth", {"method": "mean", "window": 3}),
        ],
        ids=["rotate", "near_field", "coord_shift", "coord_elev_smooth"],
    )
    def test_preview_does_not_raise_or_push_stack(self, fn_name, kwargs, small_sites):
        from pycsamt.app.desktop.controllers.correction_controller import (
            CorrectionController,
        )

        c = CorrectionController()
        c.set_raw_sites(small_sites)
        result = c.preview(fn_name, kwargs)
        assert result is not None
        assert c.n_steps == 0


# ═══════════════════════════════════════════════════════════════════════════
# 4. Stratagem convenience layer
# ═══════════════════════════════════════════════════════════════════════════


class TestStratagemMethods:
    def test_strat_fn_without_edi_dir_raises(self, ctrl):
        with pytest.raises(RuntimeError, match="Load EDI Dir"):
            ctrl._strat_qc()

    def test_strat_qc_populates_report(self, ctrl, l18_dir):
        ctrl.load_edi_dir(l18_dir)
        step = ctrl.apply("_strat_qc", {}, label="QC")
        assert step is not None
        assert ctrl._strat_qc_report is not None
        assert not ctrl._strat_qc_report.empty
        # QC is diagnostic only — does not mutate edi_objects_
        assert len(ctrl._strat_edi_objects) == 28

    def test_strat_static_shift_populates_factors(self, ctrl, l18_dir):
        ctrl.load_edi_dir(l18_dir)
        step = ctrl.apply(
            "_strat_static_shift",
            {"sort_by": "lon", "half_window": 3, "weights": "tri"},
            label="SS",
        )
        assert step is not None
        assert ctrl._strat_ss_factors is not None

    def test_strat_static_shift_with_period_band(self, ctrl, l18_dir):
        """Exercise the pband-set branch (pband_lo>0 and pband_hi>pband_lo)."""
        ctrl.load_edi_dir(l18_dir)
        step = ctrl.apply(
            "_strat_static_shift",
            {"pband_lo": 0.001, "pband_hi": 10.0},
            label="SS-band",
        )
        assert step is not None

    def test_strat_noise_removal(self, ctrl, l18_dir):
        ctrl.load_edi_dir(l18_dir)
        step = ctrl.apply(
            "_strat_noise_removal",
            {"mains_hz": 50.0, "n_harm": 5, "smooth": True},
            label="Noise",
        )
        assert step is not None

    def test_strat_freq_filter(self, ctrl, l18_dir):
        ctrl.load_edi_dir(l18_dir)
        step = ctrl.apply(
            "_strat_freq_filter",
            {"fmin": 0.0, "fmax": 0.0, "snr_thresh": 1.0, "min_frac": 0.1},
            label="Filter",
        )
        assert step is not None

    def test_strat_full_pipeline(self, ctrl, l18_dir):
        ctrl.load_edi_dir(l18_dir)
        step = ctrl.apply("_strat_full_pipeline", {}, label="Full")
        assert step is not None
        assert ctrl._strat_qc_report is not None
        assert ctrl._strat_ss_factors is not None

    def test_stratagem_steps_chain_and_undo(self, ctrl, l18_dir):
        ctrl.load_edi_dir(l18_dir)
        ctrl.apply("_strat_qc", {}, label="qc")
        ctrl.apply("_strat_freq_filter", {}, label="ff")
        ctrl.apply("_strat_static_shift", {}, label="ss")
        assert ctrl.n_steps == 3
        ctrl.undo_last()
        assert ctrl.n_steps == 2
        assert ctrl.stack[-1].label == "ff"

    def test_stratagem_remove_step_replays(self, ctrl, l18_dir):
        ctrl.load_edi_dir(l18_dir)
        ctrl.apply("_strat_qc", {}, label="qc")
        ctrl.apply("_strat_freq_filter", {}, label="ff")
        ctrl.apply("_strat_static_shift", {}, label="ss")
        ctrl.remove_step(0)
        assert ctrl.n_steps == 2
        assert [s.label for s in ctrl.stack] == ["ff", "ss"]

    def test_preview_strat_does_not_persist_state(self, ctrl, l18_dir):
        ctrl.load_edi_dir(l18_dir)
        before_report = ctrl._strat_qc_report
        before_n = len(ctrl._strat_edi_objects)
        result = ctrl.preview("_strat_qc", {})
        assert result is not None
        # preview's `finally` clause restores a *deep copy* of the saved
        # working state — object identity is never preserved (a fresh
        # deepcopy is taken up front), but the report/factors and the
        # working-copy content must be unaffected by the preview call.
        assert ctrl._strat_qc_report is before_report  # still None
        assert len(ctrl._strat_edi_objects) == before_n

    def test_reset_strat_to_raw(self, ctrl, l18_dir):
        ctrl.load_edi_dir(l18_dir)
        ctrl.apply("_strat_static_shift", {}, label="ss")
        mutated = ctrl._strat_edi_objects
        ctrl.reset_strat_to_raw()
        assert ctrl._strat_edi_objects is not mutated
        assert ctrl._strat_qc_report is None
        assert ctrl._strat_ss_factors is None

    def test_export_stratagem_no_data_raises(self, ctrl):
        with pytest.raises(RuntimeError, match="No EDI objects"):
            ctrl.export_stratagem("/tmp/wherever")

    def test_export_stratagem_creates_directory(self, ctrl, l18_dir, tmp_path):
        """
        export_stratagem() must at least create the target directory.

        NOTE (bug, not fixed — see task constraints): the per-EDI write
        inside export_stratagem() calls ``edi.write_edifile(savepath=...,
        new_edifilename=...)``, but ``pycsamt.seg.edi.EDIFile`` has no
        ``write_edifile`` method (it only exposes ``write`` and
        ``write_new_edi``). The AttributeError is swallowed by a bare
        ``except Exception: pass`` inside the per-file loop, so
        export_stratagem() currently always returns 0 and never writes a
        single .edi file, regardless of input. This test documents the
        actual (buggy) behaviour rather than the intended one.
        """
        ctrl.load_edi_dir(l18_dir)
        out_dir = tmp_path / "export_out"
        written = ctrl.export_stratagem(str(out_dir), overwrite=False)
        assert out_dir.is_dir()  # directory creation happens unconditionally
        assert written == 0  # see docstring above — real bug, not fixed here
        assert list(out_dir.glob("*.edi")) == []


# ═══════════════════════════════════════════════════════════════════════════
# 5. Plotting — matplotlib Agg smoke tests
# ═══════════════════════════════════════════════════════════════════════════


@pytest.fixture(scope="session")
def rotated_small_sites(small_sites):
    from pycsamt.app.desktop.controllers.correction_controller import (
        CorrectionController,
    )

    c = CorrectionController()
    c.set_raw_sites(small_sites)
    step = c.apply("_wrap_rotate", {"angle": 12.0}, label="rot12")
    return step.sites_after


class TestModuleHelpers:
    def test_init_ax_sets_background(self):
        from pycsamt.app.desktop.controllers.correction_controller import (
            _DARK,
            _init_ax,
        )

        fig, ax = plt.subplots()
        _init_ax(ax, _DARK)
        assert ax.get_facecolor() is not None
        _close()

    def test_empty_draws_text(self):
        from pycsamt.app.desktop.controllers.correction_controller import (
            _DARK,
            _empty,
        )

        fig, ax = plt.subplots()
        _empty(ax, "No data loaded", _DARK)
        texts = [t.get_text() for t in ax.texts]
        assert "No data loaded" in texts
        _close()

    def test_df_from_passthrough_dataframe(self, ctrl):
        import pandas as pd

        df = pd.DataFrame({"a": [1, 2]})
        assert ctrl._df_from(df) is df

    def test_df_from_none(self, ctrl):
        assert ctrl._df_from(None) is None

    def test_df_from_sites(self, ctrl, small_sites):
        result = ctrl._df_from(small_sites)
        assert result is not None
        assert "lat" in result.columns


class TestRhoPlots:
    def test_plot_rho_curves(self, ctrl, small_sites):
        fig = _fig()
        ax = fig.add_subplot(111)
        ctrl.plot_rho_curves(small_sites, ax)
        assert len(ax.lines) > 0
        _close()

    def test_plot_rho_curves_no_data(self, ctrl):
        fig = _fig()
        ax = fig.add_subplot(111)
        ctrl.plot_rho_curves(None, ax)
        _close()

    def test_plot_rho_pseudosection(self, ctrl, small_sites):
        fig = _fig()
        ctrl.plot_rho_pseudosection(small_sites, fig)
        assert len(fig.axes) >= 1
        _close()

    def test_plot_rho_pseudosection_no_data(self, ctrl):
        fig = _fig()
        ctrl.plot_rho_pseudosection(None, fig)
        assert len(fig.axes) >= 1
        _close()

    def test_plot_rho_pseudosection_with_affected_stations(self, ctrl, small_sites):
        fig = _fig()
        from pycsamt.emtools._core import _iter_items, _name

        names = [_name(ed, i) for i, ed in enumerate(_iter_items(small_sites))]
        ctrl.plot_rho_pseudosection(small_sites, fig, affected_stations=names[:2])
        _close()

    def test_plot_overlay(self, ctrl, small_sites, rotated_small_sites):
        fig = _fig()
        ax = fig.add_subplot(111)
        ctrl.plot_overlay(small_sites, rotated_small_sites, ax)
        _close()

    def test_plot_diff(self, ctrl, small_sites, rotated_small_sites):
        fig = _fig()
        ax = fig.add_subplot(111)
        ctrl.plot_diff(small_sites, rotated_small_sites, ax)
        _close()

    def test_plot_diff_missing_before_or_after(self, ctrl, small_sites):
        fig = _fig()
        ax = fig.add_subplot(111)
        ctrl.plot_diff(None, small_sites, ax)
        ctrl.plot_diff(small_sites, None, ax)
        _close()


class TestStationMapPlots:
    def test_plot_station_map_with_sites(self, ctrl, small_sites):
        fig = _fig()
        ax = fig.add_subplot(111)
        ctrl.plot_station_map(small_sites, ax, title="Stations")
        _close()

    def test_plot_station_map_with_coords_df(self, ctrl, small_sites):
        ctrl.set_raw_sites(small_sites)
        fig = _fig()
        ax = fig.add_subplot(111)
        ctrl.plot_station_map(ctrl.raw_coords_df, ax)
        _close()

    def test_plot_station_map_no_data(self, ctrl):
        fig = _fig()
        ax = fig.add_subplot(111)
        ctrl.plot_station_map(None, ax)
        _close()

    def test_plot_station_map_overlay(self, ctrl, small_sites, rotated_small_sites):
        fig = _fig()
        ax = fig.add_subplot(111)
        ctrl.plot_station_map_overlay(small_sites, rotated_small_sites, ax)
        _close()

    def test_plot_station_map_overlay_no_data(self, ctrl):
        fig = _fig()
        ax = fig.add_subplot(111)
        ctrl.plot_station_map_overlay(None, None, ax)
        _close()

    def test_plot_station_elevation(self, ctrl, small_sites):
        fig = _fig()
        ax = fig.add_subplot(111)
        ctrl.plot_station_elevation(small_sites, ax)
        _close()

    def test_plot_station_elevation_no_data(self, ctrl):
        fig = _fig()
        ax = fig.add_subplot(111)
        ctrl.plot_station_elevation(None, ax)
        _close()

    def test_plot_station_elevation_overlay(
        self, ctrl, small_sites, rotated_small_sites
    ):
        fig = _fig()
        ax = fig.add_subplot(111)
        ctrl.plot_station_elevation_overlay(small_sites, rotated_small_sites, ax)
        _close()

    def test_plot_station_elevation_overlay_no_data(self, ctrl):
        fig = _fig()
        ax = fig.add_subplot(111)
        ctrl.plot_station_elevation_overlay(None, None, ax)
        _close()

    def test_plot_displacement_diff(self, ctrl, small_sites, rotated_small_sites):
        fig = _fig()
        ax = fig.add_subplot(111)
        ctrl.plot_displacement_diff(small_sites, rotated_small_sites, ax)
        _close()

    def test_plot_displacement_diff_no_data(self, ctrl):
        fig = _fig()
        ax = fig.add_subplot(111)
        ctrl.plot_displacement_diff(None, None, ax)
        _close()


class TestRotationRosePlot:
    def test_plot_rotation_rose(self, ctrl, small_sites, rotated_small_sites):
        fig = _fig()
        ctrl.plot_rotation_rose(small_sites, rotated_small_sites, fig)
        assert len(fig.axes) == 2
        _close()

    def test_plot_rotation_rose_no_data(self, ctrl):
        fig = _fig()
        ctrl.plot_rotation_rose(None, None, fig)
        assert len(fig.axes) == 2
        _close()

    def test_draw_strike_rose_direct(self, ctrl, small_sites):
        from pycsamt.app.desktop.controllers.correction_controller import (
            _DARK,
        )

        fig = plt.figure()
        ax = fig.add_subplot(111, projection="polar")
        ctrl._draw_strike_rose(small_sites, ax, _DARK)
        _close()


class TestStratagemPlots:
    def test_plot_strat_qc_no_report(self, ctrl):
        fig = _fig()
        ctrl.plot_strat_qc(fig)
        assert len(fig.axes) >= 1
        _close()

    def test_plot_strat_qc_with_report(self, ctrl, l18_dir):
        ctrl.load_edi_dir(l18_dir)
        ctrl.apply("_strat_qc", {}, label="qc")
        fig = _fig()
        ctrl.plot_strat_qc(fig)
        assert len(fig.axes) >= 1
        _close()

    def test_plot_strat_ss_factors_no_data(self, ctrl):
        fig = _fig()
        ax = fig.add_subplot(111)
        ctrl.plot_strat_ss_factors(ax)
        _close()

    def test_plot_strat_ss_factors_with_data(self, ctrl, l18_dir):
        ctrl.load_edi_dir(l18_dir)
        ctrl.apply("_strat_static_shift", {}, label="ss")
        fig = _fig()
        ax = fig.add_subplot(111)
        ctrl.plot_strat_ss_factors(ax)
        assert len(ax.patches) > 0  # bars drawn
        _close()


class TestDarkLightModePlots:
    def test_plot_rho_curves_light_mode(self, ctrl, small_sites):
        ctrl.dark = False
        fig = _fig()
        ax = fig.add_subplot(111)
        ctrl.plot_rho_curves(small_sites, ax)
        ctrl.dark = True
        _close()

    def test_plot_station_map_light_mode(self, ctrl, small_sites):
        ctrl.dark = False
        fig = _fig()
        ax = fig.add_subplot(111)
        ctrl.plot_station_map(small_sites, ax)
        ctrl.dark = True
        _close()
