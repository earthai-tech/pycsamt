# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for InterpController — Interpretation Studio business logic.

Real data
---------
data/MT/kap03lmt_edis/         — 26 KP TIPPER EDIs
data/AMT/WILLY_DATA/L18PLT/    — 28 WILLY AMT EDIs (one profile, for speed)
data/occam2D/                  — real Occam2D inversion result directory
data/drill_example_files/      — real borehole/geology CSV sample

Strategy
--------
* State/status/setup methods are tested directly against their contract.
* A small synthetic :class:`~pycsamt.interp._base.ResistivityModel` (3
  stations x 5 depths) is used for the geological/hydro/MC/plot sweep —
  built the same way as pycsamt/interp/tests/test_hydromodel.py and
  test_calibrate.py so a water table is actually detectable at 2 of the
  3 stations.
* Real EDI Sites (WILLY / KAP03) back set_sites() and the EM-diagnostics
  plot_* methods that dispatch into pycsamt.emtools.
* generate() and every plot_* method must never raise — they return a
  matplotlib Figure, annotated with an error/placeholder message when
  inputs are missing or the underlying call fails.  Several real bugs
  were found this way (see the end of this file / the task report) and
  are captured as explicit regression tests documenting current, if
  broken, behaviour — production code is intentionally left untouched.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest

from pycsamt.interp._base import ResistivityModel

# ── Paths ─────────────────────────────────────────────────────────────────────

_ROOT = Path(__file__).parents[3]  # pycsamt/
_TIPPER = _ROOT / "data" / "MT" / "kap03lmt_edis"
_WILLY_L18 = _ROOT / "data" / "AMT" / "WILLY_DATA" / "L18PLT"
_OCCAM2D = _ROOT / "data" / "occam2D"
_DRILL_CSV = _ROOT / "data" / "drill_example_files" / "Drill&GeologydataT.csv"

_HAS_TIPPER = _TIPPER.exists() and any(_TIPPER.glob("*.edi"))
_HAS_WILLY = _WILLY_L18.exists() and any(_WILLY_L18.glob("*.edi"))
_HAS_OCCAM2D = _OCCAM2D.exists() and (_OCCAM2D / "Occam2DMesh").exists()


def _close():
    plt.close("all")


# ── Session-scoped real-data fixtures ───────────────────────────────────────


@pytest.fixture(scope="session")
def tipper_sites():
    pytest.importorskip("pycsamt.emtools")
    if not _HAS_TIPPER:
        pytest.skip("TIPPER data not available")
    from pycsamt.emtools import ensure_sites

    return ensure_sites(str(_TIPPER))


@pytest.fixture(scope="session")
def willy_sites():
    """One WILLY profile (28 stations) — enough for smoke coverage, fast."""
    pytest.importorskip("pycsamt.emtools")
    if not _HAS_WILLY:
        pytest.skip("WILLY data not available")
    from pycsamt.emtools import ensure_sites

    return ensure_sites(str(_WILLY_L18))


# ── Controller / model fixtures ─────────────────────────────────────────────


@pytest.fixture()
def ctrl():
    from pycsamt.app.desktop.controllers.interp_controller import (
        InterpController,
    )

    return InterpController()


def _make_model(station_names=("S1", "S2", "S3"), rms=1.5):
    """3 stations x 5 depths.  S1 & S3 show a resistive->conductive
    transition (detectable water table); S2 stays uniformly resistive."""
    rho = np.log10(
        np.array(
            [
                [800.0, 800.0, 800.0],
                [600.0, 800.0, 700.0],
                [15.0, 800.0, 20.0],
                [10.0, 800.0, 12.0],
                [8.0, 800.0, 9.0],
            ]
        )
    )
    return ResistivityModel.from_array(
        rho,
        x_centers=np.array([0.0, 500.0, 1000.0]),
        z_centers=np.array([5.0, 15.0, 30.0, 60.0, 100.0]),
        station_x=np.array([0.0, 500.0, 1000.0]),
        station_names=list(station_names) if station_names else None,
        method="TDEM",
        rms=rms,
    )


@pytest.fixture()
def model():
    return _make_model()


@pytest.fixture()
def loaded_ctrl(ctrl, model):
    """Controller with a synthetic model loaded (no sites)."""
    ctrl.set_model(model)
    return ctrl


def _cfg_kwargs():
    return dict(m=1.8, n=2.0, a=1.0, rho_w=20.0, phi=0.25, d50_m=5e-4)


@pytest.fixture()
def hydro_ctrl(loaded_ctrl):
    """Controller with model + hydro_result populated."""
    loaded_ctrl.set_petro_config(**_cfg_kwargs())
    msg = loaded_ctrl.run_hydro()
    assert msg == "Hydro estimation complete."
    return loaded_ctrl


@pytest.fixture()
def geo_ctrl(loaded_ctrl):
    """Controller with model + strat_logs populated (no boreholes)."""
    msg = loaded_ctrl.run_geological()
    assert msg == "Classified 3 stations."
    return loaded_ctrl


@pytest.fixture()
def mc_ctrl(hydro_ctrl):
    """Controller with model + hydro_result + mc_result (small n_samples)."""
    msg = hydro_ctrl.run_monte_carlo(n_samples=8)
    assert msg == "MC complete (8 samples)."
    return hydro_ctrl


# ═════════════════════════════════════════════════════════════════════════
# InterpState / state-status surface
# ═════════════════════════════════════════════════════════════════════════


class TestInterpState:
    def test_default_fields(self):
        from pycsamt.app.desktop.controllers.interp_controller import (
            InterpState,
        )

        st = InterpState()
        assert st.sites is None
        assert st.model is None
        assert st.boreholes == []
        assert st.db is None
        assert st.petro_cfg is None
        assert st.strat_logs == []
        assert st.hydro_result is None
        assert st.mc_result is None
        assert st.constraints == []
        assert st.fusion_model is None
        assert st.timelapse_surveys == []

    def test_mutable_defaults_are_independent(self):
        from pycsamt.app.desktop.controllers.interp_controller import (
            InterpState,
        )

        a, b = InterpState(), InterpState()
        a.boreholes.append("x")
        assert b.boreholes == []


class TestControllerConstruction:
    def test_creates(self, ctrl):
        assert ctrl is not None

    def test_default_dark(self, ctrl):
        assert ctrl.dark is True

    def test_has_model_false_initially(self, ctrl):
        assert ctrl.has_model is False

    def test_has_sites_false_initially(self, ctrl):
        assert ctrl.has_sites is False

    def test_model_status_not_loaded(self, ctrl):
        assert ctrl.model_status == {"loaded": False}


class TestSetSitesModel:
    def test_set_sites(self, ctrl):
        ctrl.set_sites("dummy-sites")
        assert ctrl.state.sites == "dummy-sites"
        assert ctrl.has_sites is True

    def test_set_model(self, ctrl, model):
        ctrl.set_model(model)
        assert ctrl.state.model is model
        assert ctrl.has_model is True

    def test_model_status_loaded(self, loaded_ctrl):
        st = loaded_ctrl.model_status
        assert st["loaded"] is True
        assert st["n_x"] == 3
        assert st["n_z"] == 5
        assert st["depth_max"] == pytest.approx(100.0)
        assert st["profile_m"] == pytest.approx(1000.0)
        assert st["method"] == "TDEM"
        assert st["rms"] == pytest.approx(1.5)

    def test_model_status_missing_attrs_fallback(self, ctrl):
        class Bare:
            pass

        ctrl.set_model(Bare())
        st = ctrl.model_status
        assert st["loaded"] is True
        assert st["n_x"] == "?"
        assert st["n_z"] == "?"
        assert st["method"] == "—"  # em-dash fallback
        assert st["rms"] is None


class TestSetModelFromOccam2d:
    def test_real_occam2d_dir_raises_due_to_wrong_import_path(self, ctrl):
        """BUG: set_model_from_occam2d() does
        `from pycsamt.models.occam2d.plot import InversionResult`, but
        InversionResult actually lives in `pycsamt.models.occam2d.results`
        (and is re-exported from `pycsamt.models.occam2d`). `plot.py` does
        not define or re-export that name, so this call always raises
        ImportError — even against a fully valid Occam2D result directory
        like data/occam2D/. Documented here rather than fixed."""
        if not _HAS_OCCAM2D:
            pytest.skip("occam2D sample data not available")
        with pytest.raises(ImportError, match="InversionResult"):
            ctrl.set_model_from_occam2d(str(_OCCAM2D))
        assert ctrl.has_model is False

    def test_bad_dir_also_raises_importerror_first(self, ctrl, tmp_path):
        """Even a nonexistent/empty result dir hits the same broken
        import before any directory validation happens."""
        empty = tmp_path / "empty_occam"
        empty.mkdir()
        with pytest.raises(ImportError, match="InversionResult"):
            ctrl.set_model_from_occam2d(str(empty))

    def test_workaround_build_model_directly_for_downstream_tests(self, ctrl):
        """Confirms the documented workaround (used throughout this file):
        since set_model_from_occam2d() is broken, downstream plot_* tests
        build a ResistivityModel directly via set_model(), which is how
        main_window.py itself detects an already-built model."""
        m = _make_model()
        ctrl.set_model(m)
        assert ctrl.has_model is True
        assert ctrl.model_status["loaded"] is True


# ═════════════════════════════════════════════════════════════════════════
# Boreholes / rock DB / petro config
# ═════════════════════════════════════════════════════════════════════════


class TestBoreholes:
    def test_add_borehole_csv(self, ctrl, tmp_path):
        p = tmp_path / "Bo17.csv"
        p.write_text(
            "top,bottom,lithology,resistivity\n" "0,20,Clay,10.0\n" "20,40,Sand,150.0\n"
        )
        name = ctrl.add_borehole_csv(str(p))
        assert name == "Bo17"
        assert len(ctrl.state.boreholes) == 1
        assert ctrl.state.boreholes[0].name == "Bo17"

    def test_add_borehole_csv_real_sample(self, ctrl):
        """Reuses the repo's real drill/geology CSV sample if compatible;
        otherwise this is skipped (its schema is a survey-station log,
        not the top/bottom/lithology/resistivity Borehole.from_csv schema
        — kept as a documented probe rather than forcing a fit)."""
        if not _DRILL_CSV.exists():
            pytest.skip("drill example CSV not available")
        try:
            ctrl.add_borehole_csv(str(_DRILL_CSV))
        except Exception as exc:
            pytest.skip(f"Drill sample CSV schema incompatible: {exc}")

    def test_add_borehole_las(self, ctrl, tmp_path):
        las = """~VERSION INFORMATION
 VERS.                  2.0:   CWLS log ASCII Standard
~WELL INFORMATION
 STRT.M          0.0:   START DEPTH
 STOP.M          4.0:   STOP DEPTH
 STEP.M          1.0:   STEP
 WELL.           TestWell:  WELL
~CURVE INFORMATION
 DEPT.M                 :  DEPTH
 RESD.OHMM              :  RESISTIVITY
 LITH.                  :  LITHOLOGY CODE
~A  DEPT         RESD         LITH
   0.0000      100.00000       1
   1.0000      100.00000       1
   2.0000      100.00000       1
   3.0000      500.00000       2
   4.0000      500.00000       2
"""
        p = tmp_path / "well1.las"
        p.write_text(las)
        name = ctrl.add_borehole_las(str(p))
        assert name == "TestWell"
        assert len(ctrl.state.boreholes) == 1

    def test_remove_borehole(self, ctrl, tmp_path):
        p = tmp_path / "Bo1.csv"
        p.write_text("top,bottom,lithology,resistivity\n0,10,Clay,10.0\n")
        ctrl.add_borehole_csv(str(p))
        p2 = tmp_path / "Bo2.csv"
        p2.write_text("top,bottom,lithology,resistivity\n0,10,Sand,50.0\n")
        ctrl.add_borehole_csv(str(p2))
        assert len(ctrl.state.boreholes) == 2
        ctrl.remove_borehole("Bo1")
        assert [b.name for b in ctrl.state.boreholes] == ["Bo2"]

    def test_remove_nonexistent_borehole_noop(self, ctrl, tmp_path):
        p = tmp_path / "Bo1.csv"
        p.write_text("top,bottom,lithology,resistivity\n0,10,Clay,10.0\n")
        ctrl.add_borehole_csv(str(p))
        ctrl.remove_borehole("DoesNotExist")
        assert len(ctrl.state.boreholes) == 1


class TestRockDb:
    def test_set_rock_db_default(self, ctrl):
        from pycsamt.geology.lithology import RockDatabase

        ctrl.set_rock_db_default()
        assert ctrl.state.db is not None
        assert len(ctrl.state.db) == len(RockDatabase.default())

    def test_set_rock_db_csv(self, ctrl, tmp_path):
        p = tmp_path / "rocks.csv"
        p.write_text(
            "name,rho_min,rho_max,color,description\n"
            "MyClay,1,20,#123456,test clay\n"
            "MySand,20,300,#654321,test sand\n"
        )
        ctrl.set_rock_db_csv(str(p))
        assert len(ctrl.state.db) == 2


class TestPetroConfig:
    def test_defaults(self, ctrl):
        ctrl.set_petro_config()
        cfg = ctrl.state.petro_cfg
        assert cfg.petro.m == pytest.approx(2.0)
        assert cfg.petro.n == pytest.approx(2.0)
        assert cfg.petro.a == pytest.approx(1.0)
        assert cfg.rho_w == pytest.approx(10.0)
        assert cfg.porosity_prior == pytest.approx(0.35)
        assert cfg.d50_m == pytest.approx(5e-4)

    def test_custom_kwargs(self, ctrl):
        ctrl.set_petro_config(m=1.8, n=2.1, a=0.9, rho_w=25.0, phi=0.3, d50_m=1e-3)
        cfg = ctrl.state.petro_cfg
        assert cfg.petro.m == pytest.approx(1.8)
        assert cfg.petro.n == pytest.approx(2.1)
        assert cfg.petro.a == pytest.approx(0.9)
        assert cfg.rho_w == pytest.approx(25.0)
        assert cfg.porosity_prior == pytest.approx(0.3)
        assert cfg.d50_m == pytest.approx(1e-3)


# ═════════════════════════════════════════════════════════════════════════
# Computation: run_geological / run_hydro / run_monte_carlo
# ═════════════════════════════════════════════════════════════════════════


class TestRunGeological:
    def test_no_model(self, ctrl):
        assert ctrl.run_geological() == "No model loaded."

    def test_classifies_stations_default_db(self, loaded_ctrl):
        msg = loaded_ctrl.run_geological()
        assert msg == "Classified 3 stations."
        assert len(loaded_ctrl.state.strat_logs) == 3
        assert loaded_ctrl.state.db is not None  # auto-created

    def test_with_existing_db(self, loaded_ctrl):
        loaded_ctrl.set_rock_db_default()
        db_before = loaded_ctrl.state.db
        loaded_ctrl.run_geological()
        assert loaded_ctrl.state.db is db_before  # reused, not replaced

    def test_with_boreholes_triggers_calibration(self, loaded_ctrl, tmp_path):
        p = tmp_path / "Bo1.csv"
        p.write_text(
            "top,bottom,lithology,resistivity\n"
            "0,20,Clay,700.0\n"
            "20,100,Sand,15.0\n"
        )
        loaded_ctrl.add_borehole_csv(str(p))
        # Borehole x defaults to 0.0 -> coincides with station S1 (x=0)
        assert loaded_ctrl.state.boreholes[0].x == 0.0
        msg = loaded_ctrl.run_geological()
        assert msg.startswith("Classified")
        # Either calibration succeeded (model replaced) or silently
        # fell back (except: pass) — both leave strat_logs populated.
        assert len(loaded_ctrl.state.strat_logs) >= 1


class TestRunHydro:
    def test_no_model(self, ctrl):
        assert ctrl.run_hydro() == "No model loaded."

    def test_default_petro_cfg_auto_created(self, loaded_ctrl):
        assert loaded_ctrl.state.petro_cfg is None
        msg = loaded_ctrl.run_hydro()
        assert msg == "Hydro estimation complete."
        assert loaded_ctrl.state.petro_cfg is not None

    def test_result_has_expected_arrays(self, hydro_ctrl):
        res = hydro_ctrl.state.hydro_result
        assert res.hydraulic_K.shape == (5, 3)
        assert res.saturation.shape == (5, 3)
        assert res.water_table.shape == (3,)
        assert res.transmissivity.shape == (3,)


class TestRunMonteCarlo:
    def test_no_model_or_cfg(self, ctrl):
        assert ctrl.run_monte_carlo() == "Run hydro estimation first."

    def test_model_but_no_cfg(self, loaded_ctrl):
        assert loaded_ctrl.run_monte_carlo() == "Run hydro estimation first."

    def test_small_run_succeeds(self, hydro_ctrl):
        msg = hydro_ctrl.run_monte_carlo(n_samples=8)
        assert msg == "MC complete (8 samples)."
        assert hydro_ctrl.state.mc_result is not None


# ═════════════════════════════════════════════════════════════════════════
# generate() dispatcher
# ═════════════════════════════════════════════════════════════════════════


class TestGenerate:
    def test_unknown_method_not_implemented(self, ctrl):
        fig = ctrl.generate("plot_totally_bogus_method")
        assert fig is not None
        texts = " ".join(t.get_text() for ax in fig.axes for t in ax.texts)
        assert "not yet implemented" in texts
        _close()

    def test_known_method_no_model(self, ctrl):
        fig = ctrl.generate("plot_model_summary")
        assert fig is not None
        _close()

    def test_exception_inside_method_is_caught(self, ctrl, monkeypatch):
        def _boom(**kw):
            raise RuntimeError("kaboom")

        monkeypatch.setattr(ctrl, "plot_model_summary", _boom, raising=False)
        fig = ctrl.generate("plot_model_summary")
        assert fig is not None
        texts = " ".join(t.get_text() for ax in fig.axes for t in ax.texts)
        assert "failed" in texts.lower()
        _close()

    def test_kwargs_forwarded(self, geo_ctrl):
        fig = geo_ctrl.generate("plot_strat_log", station="S2")
        assert fig is not None
        _close()


# ═════════════════════════════════════════════════════════════════════════
# Setup & Model plots
# ═════════════════════════════════════════════════════════════════════════


class TestPlotModelSummary:
    def test_no_model(self, ctrl):
        fig = ctrl.plot_model_summary()
        texts = " ".join(t.get_text() for ax in fig.axes for t in ax.texts)
        assert "no resistivity model" in texts.lower()
        _close()

    def test_with_model(self, loaded_ctrl):
        fig = loaded_ctrl.plot_model_summary()
        assert fig is not None
        assert len(fig.axes) >= 1
        _close()

    def test_dark_and_light(self, loaded_ctrl):
        for dark in (True, False):
            loaded_ctrl.dark = dark
            fig = loaded_ctrl.plot_model_summary()
            assert fig is not None
            _close()


class TestPlotDepthCoverage:
    def test_no_sites(self, ctrl):
        fig = ctrl.plot_depth_coverage()
        texts = " ".join(t.get_text() for ax in fig.axes for t in ax.texts)
        assert "no survey data" in texts.lower()
        _close()

    def test_with_sites(self, ctrl, willy_sites):
        ctrl.set_sites(willy_sites)
        fig = ctrl.plot_depth_coverage()
        assert fig is not None
        _close()


class TestPlotBoreholeMap:
    def test_no_model_no_boreholes(self, ctrl):
        fig = ctrl.plot_borehole_map()
        texts = " ".join(t.get_text() for ax in fig.axes for t in ax.texts)
        assert "no model or boreholes" in texts.lower()
        _close()

    def test_with_model_only(self, loaded_ctrl):
        fig = loaded_ctrl.plot_borehole_map()
        assert fig is not None
        _close()

    def test_with_model_and_boreholes(self, loaded_ctrl, tmp_path):
        p = tmp_path / "Bo1.csv"
        p.write_text("top,bottom,lithology,resistivity\n0,10,Clay,10.0\n")
        loaded_ctrl.add_borehole_csv(str(p))
        fig = loaded_ctrl.plot_borehole_map()
        assert fig is not None
        _close()


# ═════════════════════════════════════════════════════════════════════════
# Geological plots
# ═════════════════════════════════════════════════════════════════════════


class TestPlotStratLog:
    def test_no_logs(self, ctrl):
        fig = ctrl.plot_strat_log()
        texts = " ".join(t.get_text() for ax in fig.axes for t in ax.texts)
        assert "run geological classification" in texts.lower()
        _close()

    def test_with_logs_default_station(self, geo_ctrl):
        fig = geo_ctrl.plot_strat_log()
        assert fig is not None
        _close()

    def test_with_logs_named_station(self, geo_ctrl):
        fig = geo_ctrl.plot_strat_log(station="S2")
        assert fig is not None
        _close()

    def test_unknown_station_falls_back_to_first(self, geo_ctrl):
        fig = geo_ctrl.plot_strat_log(station="NoSuchStation")
        assert fig is not None
        _close()


class TestPlotFenceDiagram:
    def test_no_logs(self, ctrl):
        fig = ctrl.plot_fence_diagram()
        texts = " ".join(t.get_text() for ax in fig.axes for t in ax.texts)
        assert "run geological classification" in texts.lower()
        _close()

    def test_with_logs(self, geo_ctrl):
        fig = geo_ctrl.plot_fence_diagram()
        assert fig is not None
        _close()


class TestPlotCalibratedModel:
    def test_no_model(self, ctrl):
        fig = ctrl.plot_calibrated_model()
        texts = " ".join(t.get_text() for ax in fig.axes for t in ax.texts)
        assert "no resistivity model" in texts.lower()
        _close()

    def test_with_model_shows_error_fig(self, loaded_ctrl):
        """BUG: PlotCalibratedModel(crm, nm, ...) requires *both* crm and nm
        as positional args, but plot_calibrated_model() only ever passes
        `self.state.model` (one positional arg) -> TypeError, caught
        internally and rendered as an error figure. Never fixed here."""
        fig = loaded_ctrl.plot_calibrated_model()
        assert fig is not None
        texts = " ".join(t.get_text() for ax in fig.axes for t in ax.texts)
        assert "✕" in texts or "positional argument" in texts.lower()
        _close()


class TestPlotRockDb:
    def test_via_generate_no_db_does_not_raise(self, ctrl):
        """plot_rock_db() reads `db.entries` (RockDatabase's public,
        read-only accessor), so both the generate() dispatcher and a
        direct call succeed."""
        fig = ctrl.generate("plot_rock_db")
        assert fig is not None
        _close()

    def test_direct_call_renders_one_bar_per_entry(self, ctrl):
        """Regression test for the former RockDatabase-not-iterable bug:
        plot_rock_db() now reads `db.entries` instead of `list(db)`, so a
        direct call renders successfully with one bar per database
        entry."""
        fig = ctrl.plot_rock_db()
        ax = fig.axes[0]
        from pycsamt.geology.lithology import RockDatabase

        assert len(ax.patches) == len(RockDatabase.default())
        _close()


# ═════════════════════════════════════════════════════════════════════════
# Hydrology plots
# ═════════════════════════════════════════════════════════════════════════


class TestHydroSectionPlots:
    def test_no_hydro_result(self, ctrl):
        fig = ctrl.plot_K_map()
        texts = " ".join(t.get_text() for ax in fig.axes for t in ax.texts)
        assert "run hydrology estimation" in texts.lower()
        _close()

    def test_plot_K_map_attribute_name_bug(self, hydro_ctrl):
        """BUG: plot_K_map() reads `hydro_result.hydraulic_conductivity`,
        but EMHydroResult stores the array as `hydraulic_K`. getattr(...,
        default=None) silently returns None, so the section always shows
        the "'hydraulic_conductivity' not in result" placeholder instead
        of the real K section, even though hydro_result IS populated."""
        assert not hasattr(hydro_ctrl.state.hydro_result, "hydraulic_conductivity")
        assert hasattr(hydro_ctrl.state.hydro_result, "hydraulic_K")
        fig = hydro_ctrl.plot_K_map()
        texts = " ".join(t.get_text() for ax in fig.axes for t in ax.texts)
        assert "not in result" in texts
        _close()

    def test_plot_Sw_map_renders_real_data(self, hydro_ctrl):
        """`saturation` IS a real EMHydroResult attribute -> this one
        actually renders (no placeholder text)."""
        fig = hydro_ctrl.plot_Sw_map()
        assert fig is not None
        texts = " ".join(t.get_text() for ax in fig.axes for t in ax.texts)
        assert "not in result" not in texts
        _close()


class TestPlotWaterTable:
    def test_no_hydro_result(self, ctrl):
        fig = ctrl.plot_water_table()
        texts = " ".join(t.get_text() for ax in fig.axes for t in ax.texts)
        assert "run hydrology estimation" in texts.lower()
        _close()

    def test_with_result_shows_error_fig(self, hydro_ctrl):
        """BUG: PlotWaterTableProfile(result, *, ...) takes a single
        positional `result` (EMHydroResult), but plot_water_table() calls
        it as PlotWaterTableProfile(model, hydro_result) — 2 positional
        args -> TypeError, caught and rendered as an error figure."""
        fig = hydro_ctrl.plot_water_table()
        assert fig is not None
        _close()


class TestPlotAquiferZones:
    def test_no_model(self, ctrl):
        fig = ctrl.plot_aquifer_zones()
        texts = " ".join(t.get_text() for ax in fig.axes for t in ax.texts)
        assert "no resistivity model" in texts.lower()
        _close()

    def test_with_model_shows_error_text(self, loaded_ctrl):
        """BUG (double): HydroInterpreter.__init__ takes only keyword-only
        args (no positional `model`), and even if it did,
        HydroInterpreter has no `identify_aquifers()` method at all (the
        real API is `.fit(model).zones`). Both are caught by the local
        try/except inside plot_aquifer_zones() and rendered as text on
        the axes, so the call itself does not raise."""
        fig = loaded_ctrl.plot_aquifer_zones()
        assert fig is not None
        _close()


class TestPlotTransmissivity:
    def test_no_hydro_result(self, ctrl):
        fig = ctrl.plot_transmissivity()
        texts = " ".join(t.get_text() for ax in fig.axes for t in ax.texts)
        assert "run hydrology estimation" in texts.lower()
        _close()

    def test_with_result_renders(self, hydro_ctrl):
        fig = hydro_ctrl.plot_transmissivity()
        assert fig is not None
        assert len(fig.axes) >= 1
        _close()


class TestPlotAquiferChar:
    def test_no_hydro_result(self, ctrl):
        fig = ctrl.plot_aquifer_char()
        texts = " ".join(t.get_text() for ax in fig.axes for t in ax.texts)
        assert "run hydrology estimation" in texts.lower()
        _close()

    def test_with_result_shows_error_fig(self, hydro_ctrl):
        """BUG: same extra-positional-arg pattern as plot_water_table —
        PlotAquiferCharacterization(result, *, ...) vs. the controller's
        PlotAquiferCharacterization(model, hydro_result) call."""
        fig = hydro_ctrl.plot_aquifer_char()
        assert fig is not None
        _close()


class TestPlotPetrophysXplot:
    def test_no_hydro_result(self, ctrl):
        fig = ctrl.plot_petrophys_xplot()
        texts = " ".join(t.get_text() for ax in fig.axes for t in ax.texts)
        assert "run hydrology estimation" in texts.lower()
        _close()

    def test_with_result_shows_error_fig(self, hydro_ctrl):
        """BUG: same pattern — PlotPetrophysicalCrossPlot(result, *,
        petro=None, ...) vs. the controller's 3-positional-arg call
        PlotPetrophysicalCrossPlot(model, hydro_result, petro_cfg)."""
        fig = hydro_ctrl.plot_petrophys_xplot()
        assert fig is not None
        _close()


# ═════════════════════════════════════════════════════════════════════════
# Field Constraints (always placeholders — no compute path implemented)
# ═════════════════════════════════════════════════════════════════════════


class TestFieldConstraintPlots:
    def test_plot_constraint_misfit(self, hydro_ctrl):
        fig = hydro_ctrl.plot_constraint_misfit()
        texts = " ".join(t.get_text() for ax in fig.axes for t in ax.texts)
        assert "constraint calibration" in texts.lower()
        _close()

    def test_plot_calib_history(self, hydro_ctrl):
        fig = hydro_ctrl.plot_calib_history()
        texts = " ".join(t.get_text() for ax in fig.axes for t in ax.texts)
        assert "constraint calibration" in texts.lower()
        _close()


# ═════════════════════════════════════════════════════════════════════════
# EM Diagnostics (_emtools_plot dispatch)
# ═════════════════════════════════════════════════════════════════════════

_EMTOOLS_METHODS_NONPOLAR = [
    "plot_pt_section",
    "plot_strike_profile",
    "plot_dimensionality",
    "plot_bostick_depths",
    "plot_gradient_section",
    "plot_snr_section",
]
_EMTOOLS_METHODS_POLAR = ["plot_strike_rose", "plot_phasor_wheel"]


class TestEmtoolsPlotsNoSites:
    @pytest.mark.parametrize(
        "method", _EMTOOLS_METHODS_NONPOLAR + _EMTOOLS_METHODS_POLAR
    )
    def test_no_sites_shows_message(self, ctrl, method):
        fig = ctrl.generate(method)
        texts = " ".join(t.get_text() for ax in fig.axes for t in ax.texts)
        assert "no survey data" in texts.lower()
        _close()


class TestEmtoolsPlotsWithSites:
    @pytest.mark.parametrize("method", _EMTOOLS_METHODS_NONPOLAR)
    def test_does_not_raise(self, ctrl, willy_sites, method):
        ctrl.set_sites(willy_sites)
        fig = ctrl.generate(method)
        assert fig is not None
        _close()

    @pytest.mark.parametrize("method", _EMTOOLS_METHODS_POLAR)
    def test_polar_does_not_raise(self, ctrl, willy_sites, method):
        ctrl.set_sites(willy_sites)
        fig = ctrl.generate(method)
        assert fig is not None
        _close()

    def test_induction_arrows_with_tipper_sites(self, ctrl, tipper_sites):
        ctrl.set_sites(tipper_sites)
        fig = ctrl.generate("plot_induction_arrows")
        assert fig is not None
        _close()

    def test_unavailable_emtools_fn(self, ctrl, willy_sites, monkeypatch):
        ctrl.set_sites(willy_sites)
        fig = ctrl._emtools_plot("this_fn_does_not_exist_in_emtools", "X")
        texts = " ".join(t.get_text() for ax in fig.axes for t in ax.texts)
        assert "not yet implemented" in texts.lower()
        _close()


# ═════════════════════════════════════════════════════════════════════════
# Uncertainty (MC) plots
# ═════════════════════════════════════════════════════════════════════════


class TestMcPlots:
    def test_no_mc_result_all_three(self, ctrl):
        for name in (
            "plot_mc_K_section",
            "plot_mc_wt_profile",
            "plot_mc_histograms",
        ):
            fig = ctrl.generate(name)
            texts = " ".join(t.get_text() for ax in fig.axes for t in ax.texts)
            assert "run monte carlo" in texts.lower()
            _close()

    def test_plot_mc_K_section_shows_error_fig(self, mc_ctrl):
        """BUG: PlotUncertaintySection(unc, *, quantity=...) takes one
        positional arg; plot_mc_K_section() calls it with 2
        (model, mc_result) -> TypeError, caught -> error figure."""
        fig = mc_ctrl.plot_mc_K_section()
        assert fig is not None
        _close()

    def test_plot_mc_wt_profile_shows_error_fig(self, mc_ctrl):
        """BUG: same extra-positional-arg pattern as plot_mc_K_section."""
        fig = mc_ctrl.plot_mc_wt_profile()
        assert fig is not None
        _close()

    def test_plot_mc_histograms_renders_real_data(self, mc_ctrl):
        """PlotUncertaintyHistogram(mc_result) matches its real signature
        -> this one actually renders."""
        fig = mc_ctrl.plot_mc_histograms()
        assert fig is not None
        assert len(fig.axes) >= 1
        _close()


# ═════════════════════════════════════════════════════════════════════════
# Advanced Plots (_emtools_plot dispatch)
# ═════════════════════════════════════════════════════════════════════════

_ADVANCED_METHODS = [
    "plot_mohr_circles",
    "plot_argand",
    "plot_pt_ternary",
    "plot_distortion_radar",
    "plot_bode",
    "plot_z_invariants",
    "plot_anisotropy",
    "plot_period_clock",
    "plot_composite_section",
]


class TestAdvancedPlots:
    @pytest.mark.parametrize("method", _ADVANCED_METHODS)
    def test_does_not_raise_with_sites(self, ctrl, willy_sites, method):
        ctrl.set_sites(willy_sites)
        fig = ctrl.generate(method)
        assert fig is not None
        _close()

    @pytest.mark.parametrize("method", _ADVANCED_METHODS[:3])
    def test_no_sites_message(self, ctrl, method):
        fig = ctrl.generate(method)
        texts = " ".join(t.get_text() for ax in fig.axes for t in ax.texts)
        assert "no survey data" in texts.lower()
        _close()


# ═════════════════════════════════════════════════════════════════════════
# Fusion & Time-Lapse
# ═════════════════════════════════════════════════════════════════════════


class TestFusedModel:
    def test_no_fusion_model(self, ctrl):
        fig = ctrl.plot_fused_model()
        texts = " ".join(t.get_text() for ax in fig.axes for t in ax.texts)
        assert "run fusion" in texts.lower()
        _close()

    def test_with_fusion_model_does_not_raise(self, ctrl, model):
        ctrl.state.fusion_model = model
        fig = ctrl.plot_fused_model()
        assert fig is not None
        _close()


class TestTimelapseChange:
    def test_no_surveys(self, ctrl):
        fig = ctrl.plot_timelapse_change()
        texts = " ".join(t.get_text() for ax in fig.axes for t in ax.texts)
        assert "2 surveys" in texts.lower() or "time-lapse" in texts.lower()
        _close()

    def test_with_two_compatible_surveys(self, ctrl, model):
        m2 = _make_model()
        m2.rho_2d[:] = m2.rho_2d + 0.3
        ctrl.state.timelapse_surveys = [model, m2]
        fig = ctrl.plot_timelapse_change()
        assert fig is not None
        _close()

    def test_plot_timelapse_sat_always_placeholder(self, ctrl):
        fig = ctrl.plot_timelapse_sat()
        texts = " ".join(t.get_text() for ax in fig.axes for t in ax.texts)
        assert "run time-lapse" in texts.lower() or "Δsw" in texts.lower()
        _close()

    def test_plot_timelapse_wt_always_placeholder(self, ctrl):
        fig = ctrl.plot_timelapse_wt()
        texts = " ".join(t.get_text() for ax in fig.axes for t in ax.texts)
        assert "wt displacement" in texts.lower() or "run analysis" in texts.lower()
        _close()


class TestExportPreview:
    def test_dark_and_light(self, ctrl):
        for dark in (True, False):
            ctrl.dark = dark
            fig = ctrl.plot_export_preview()
            assert fig is not None
            texts = " ".join(t.get_text() for ax in fig.axes for t in ax.texts)
            assert "available exports" in texts.lower()
            _close()


# ═════════════════════════════════════════════════════════════════════════
# Export helpers
# ═════════════════════════════════════════════════════════════════════════


class TestExportHelpers:
    def test_export_xyz_no_logs(self, ctrl, tmp_path):
        msg = ctrl.export_xyz(str(tmp_path / "out.xyz"))
        assert "run geological classification" in msg.lower()

    def test_export_las_no_logs(self, ctrl, tmp_path):
        msg = ctrl.export_las(str(tmp_path / "out.las"))
        assert "run geological classification" in msg.lower()

    def test_export_csv_no_logs(self, ctrl, tmp_path):
        msg = ctrl.export_csv(str(tmp_path / "out.csv"))
        assert "run geological classification" in msg.lower()

    def test_export_vtk_no_model(self, ctrl, tmp_path):
        msg = ctrl.export_vtk(str(tmp_path / "out.vtk"))
        assert "no model loaded" in msg.lower()

    def test_export_xyz_with_logs(self, geo_ctrl, tmp_path):
        out = tmp_path / "out.xyz"
        msg = geo_ctrl.export_xyz(str(out))
        assert "exported" in msg.lower()
        assert out.exists()

    def test_export_csv_with_logs(self, geo_ctrl, tmp_path):
        out = tmp_path / "out.csv"
        msg = geo_ctrl.export_csv(str(out))
        assert "exported" in msg.lower()
        assert out.exists()

    def test_export_las_with_logs_default_station(self, geo_ctrl, tmp_path):
        out = tmp_path / "out.las"
        msg = geo_ctrl.export_las(str(out))
        assert "exported" in msg.lower()
        assert out.exists()

    def test_export_las_with_logs_named_station(self, geo_ctrl, tmp_path):
        out = tmp_path / "out2.las"
        msg = geo_ctrl.export_las(str(out), station="S2")
        assert "exported" in msg.lower()
        assert out.exists()

    def test_export_vtk_with_model(self, loaded_ctrl, tmp_path):
        out = tmp_path / "out.vtk"
        msg = loaded_ctrl.export_vtk(str(out))
        assert "exported" in msg.lower()
        assert out.exists()


# ═════════════════════════════════════════════════════════════════════════
# Internal styling / message helpers
# ═════════════════════════════════════════════════════════════════════════


class TestInternalHelpers:
    def test_apply_fig_style_dark_light(self, ctrl):
        for dark in (True, False):
            ctrl.dark = dark
            fig, ax = plt.subplots()
            ctrl._apply_fig_style(fig, ax)
            _close()

    def test_apply_fig_style_minimal_handles_bad_axes(self, ctrl):
        fig, ax = plt.subplots()
        ctrl._apply_fig_style_minimal(fig)
        _close()

    def test_style_cb(self, ctrl):
        fig, ax = plt.subplots()
        im = ax.pcolormesh(np.random.rand(3, 3))
        cb = fig.colorbar(im, ax=ax)
        ctrl._style_cb(cb)
        _close()

    def test_maybe_draw_topo_noop_without_config(self, ctrl, model):
        fig, ax = plt.subplots()
        # Should never raise even though topo is not configured/enabled.
        ctrl._maybe_draw_topo(ax, model)
        _close()

    def test_maybe_draw_topo_empty_station_x(self, ctrl):
        fig, ax = plt.subplots()

        class M:
            station_x = np.array([])

        ctrl._maybe_draw_topo(ax, M())
        _close()

    def test_msg_fig_variants(self, ctrl):
        for dark in (True, False):
            ctrl.dark = dark
            for fig in (
                ctrl._no_model_fig(),
                ctrl._no_sites_fig(),
                ctrl._needs_run_fig("do something"),
                ctrl._not_implemented("plot_x"),
                ctrl._error_fig("boom"),
            ):
                assert fig is not None
                _close()
