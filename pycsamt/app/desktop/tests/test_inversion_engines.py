# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for the Qt-free Inversion Studio engines and library store.

Real bundled data throughout: the CSAMT line (Occam1D build + run +
reload), the Occam2D / ModEM / MARE2DEM result folders (detection, loading,
views).  External solvers are not launched here (see the window tests and
the solver benchmarks).
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import numpy as np
import pytest
from matplotlib.figure import Figure

from pycsamt.app.desktop.controllers import inversion_engines as ie
from pycsamt.app.desktop.controllers import inversion_runs as store

ROOT = Path(__file__).resolve().parents[4]  # repository root
CSAMT = ROOT / Path("data/CSAMT")
KAP = ROOT / Path("data/MT/kap03lmt_edis")
OCCAM2D = ROOT / Path("data/occam2D")
MODEM3D = ROOT / Path("data/MT/broken-hill/final-models")
M2D_DEMO = ROOT / Path("data/mare2dem/demo_mt_inversion")


def _sites(folder: Path, n: int | None = None):
    if not folder.is_dir():
        pytest.skip(f"bundled data not found: {folder}")
    from pycsamt.io import read_transfer_function
    from pycsamt.site.base import Sites

    paths = sorted(folder.glob("*.edi"))[:n]
    return Sites([read_transfer_function(p) for p in paths])


# ── engine declarations ─────────────────────────────────────────────────


@pytest.mark.parametrize("key", list(ie.ENGINES))
def test_default_settings_make_a_valid_library_config(key):
    """Every declared field must map onto the engine's real config class
    (a typo would only surface when the user clicks Run)."""
    eng = ie.engine(key)
    fields = eng.fields()
    assert fields and len({f.key for f in fields}) == len(fields)
    for f in fields:
        assert f.section in ("Data", "Mesh", "Settings", "Run")
        if f.kind == "choice":
            assert f.default in [v for v, _ in f.choices]
    vals = eng.defaults()
    if key == "occam1d":
        cfg = eng._config(vals, None)
    else:
        cfg = eng._config(vals)
    assert cfg is not None


@pytest.mark.parametrize("key", list(store.BUILTIN_PRESETS))
def test_builtin_presets_only_use_real_fields(key):
    names = {f.key for f in ie.engine(key).fields()}
    for p in store.BUILTIN_PRESETS[key]:
        assert set(p.values) <= names, (p.name, set(p.values) - names)


def test_engine_registry_order_and_binaries():
    assert list(ie.ENGINES) == ["occam1d", "occam2d", "modem2d", "modem3d",
                                "mare2dem"]
    assert ie.engine("occam1d").binary_key is None
    assert ie.engine("modem3d").binary_key == "modem3d"
    assert ie.engine("modem2d").dim == "2D"


def test_check_sites_messages():
    assert "Load EDI" in ie.engine("occam1d").check_sites(None)
    assert "two stations" in ie.engine("occam2d").check_sites(
        _sites(CSAMT, 1))
    assert "three" in ie.engine("modem3d").check_sites(_sites(CSAMT, 2))
    assert ie.engine("mare2dem").check_sites(_sites(CSAMT, 3)) is None


# ── real-data features for the AI 1-D inverter ──────────────────────────


def test_sites_to_1d_features_real_mt():
    sites = _sites(KAP, 3)
    freqs = np.logspace(-3, 1, 12)
    X, names = ie.sites_to_1d_features(sites, freqs)
    assert X.shape == (3, 24) and names == [s.name for s in sites]
    assert np.isfinite(X).all()
    phase = X[:, 12:]
    assert (phase > 0).all() and (phase < 90).all()  # first quadrant


def test_sites_to_1d_features_scalar_csamt_uses_xy():
    X, names = ie.sites_to_1d_features(_sites(CSAMT, 2),
                                       np.logspace(0, 3, 8))
    assert X.shape == (2, 16) and np.isfinite(X).all()


def test_sites_to_1d_features_empty():
    assert ie.sites_to_1d_features(None, [1.0]) == (None, [])


def test_auto_1d_mode():
    assert ie._auto_1d_mode(_sites(CSAMT, 2)) == "xy"
    assert ie._auto_1d_mode(_sites(KAP, 2)) == "determinant"


# ── Occam1D: build → run → reload → views ───────────────────────────────


@pytest.fixture(scope="module")
def o1d_run(tmp_path_factory):
    sites = _sites(CSAMT, 3)
    wd = tmp_path_factory.mktemp("o1d")
    eng = ie.engine("occam1d")
    vals = eng.defaults()
    vals.update(max_iterations=3, n_layers=20)
    build = eng.build(sites, wd, vals)
    rep = ie.RunReporter()
    out = eng.make_task(build, vals, None)(rep)
    return wd, build, rep, out


def test_occam1d_build(o1d_run):
    wd, build, _rep, _out = o1d_run
    assert build.stations == ["csa000", "csa050", "csa100"]
    assert build.mesh["depth"].size == 20
    fig = Figure()
    ie.engine("occam1d").plot_mesh(fig, build)
    assert len(fig.axes) == 2


def test_occam1d_run_reports_every_iteration(o1d_run):
    _wd, _build, rep, out = o1d_run
    assert set(out["done"]) == {"csa000", "csa050", "csa100"}
    per_station = {}
    for label, n, rms in rep.history:
        per_station.setdefault(label, []).append(n)
        assert np.isfinite(rms)
    assert all(v[0] == 0 for v in per_station.values())
    assert any("iter" in line for line in rep.lines)


def test_occam1d_detect_load_and_views(o1d_run):
    wd, *_ = o1d_run
    assert ie.detect_engine(wd) == "occam1d"
    eng = ie.engine("occam1d")
    run = eng.load(wd)
    assert run.stations == ["csa000", "csa050", "csa100"]
    assert run.iterations[0] == 0 and run.final_rms is not None
    for key, _label in eng.views(run):
        fig = eng.render(key, Figure(), run)
        assert fig.axes, key
    early = eng.load(wd, iteration=1, station="csa050")
    assert early.station == "csa050" and early.iteration == 1
    sliced = ie._o1d_slice(early.result["result"], 1)
    assert sliced.final.number == 1


def test_occam1d_uninverted_station_explains(tmp_path):
    eng = ie.engine("occam1d")
    eng.build(_sites(CSAMT, 2), tmp_path, eng.defaults())
    run = eng.load(tmp_path)
    with pytest.raises(ie.PlotUnavailable, match="not been inverted"):
        eng.render("model", Figure(), run)
    with pytest.raises(ie.PlotUnavailable, match="at least two"):
        eng.render("section", Figure(), run)


def test_occam1d_cancel_between_iterations(tmp_path):
    eng = ie.engine("occam1d")
    vals = eng.defaults()
    build = eng.build(_sites(CSAMT, 2), tmp_path, vals)
    rep = ie.RunReporter()
    real = rep.iteration

    def stop_after_first(n, rms, label=""):
        real(n, rms, label)
        rep.cancel()

    rep.iteration = stop_after_first
    with pytest.raises(InterruptedError):
        eng.make_task(build, vals, None)(rep)
    assert len(rep.history) == 1


def test_occam1d_no_usable_station_fails_loudly(tmp_path):
    eng = ie.engine("occam1d")
    vals = eng.defaults()
    vals["mode"] = "determinant"  # CSAMT has no Zyx
    with pytest.raises(ValueError, match="scalar CSAMT"):
        eng.build(_sites(CSAMT, 2), tmp_path, vals)


# ── other engines: builds and bundled results ───────────────────────────


@pytest.mark.parametrize("key", ["modem2d", "modem3d", "mare2dem"])
def test_builds_write_inputs_and_mesh(key, tmp_path):
    eng = ie.engine(key)
    build = eng.build(_sites(CSAMT), tmp_path, eng.defaults())
    assert build.files and all(Path(p).is_file()
                               for p in build.files.values())
    assert build.summary and build.stations
    fig = Figure()
    eng.plot_mesh(fig, build)
    assert len(fig.axes) == 2


def test_mare2dem_warns_about_coincident_receivers(tmp_path):
    eng = ie.engine("mare2dem")
    build = eng.build(_sites(CSAMT), tmp_path, eng.defaults())
    assert any("same position" in w for w in build.warnings)
    assert build.mesh["run_stem"] == "mare2dem.0"


def test_stations_aligned_on_core_cells():
    nodes = np.concatenate([[0, 1000, 1500], 1600 + 50 * np.arange(11),
                            [2200, 2700, 3700]])
    st = ie._align(np.array([0.0, 250.0, 500.0]), nodes)
    assert st.min() >= 1600 and st.max() <= 2100
    assert ie._core_centre(nodes) == pytest.approx(1850.0)


def test_detect_and_load_bundled_results():
    if not OCCAM2D.is_dir():
        pytest.skip("bundled Occam2D run not found")
    assert ie.detect_engine(OCCAM2D) == "occam2d"
    run = ie.engine("occam2d").load(OCCAM2D)
    assert run.rms and run.result.rho_2d is not None
    fig = ie.engine("occam2d").render("convergence", Figure(), run)
    assert fig.axes


def test_modem3d_result_views_and_missing_log():
    # The ~11 MB .rho model is not tracked; only .dat/.res are.
    if not MODEM3D.is_dir() or not list(MODEM3D.glob("*.rho")):
        pytest.skip("bundled ModEM result not found")
    key = ie.detect_engine(MODEM3D)
    assert key == "modem3d"
    eng = ie.engine(key)
    run = eng.load(MODEM3D)
    assert ("section", "Section") in eng.views(run)
    with pytest.raises(ie.PlotUnavailable, match="convergence"):
        eng.render("convergence", Figure(), run)
    fig = eng.render("model", Figure(), run)
    assert len(fig.axes) > 1


def test_mare2dem_demo_views():
    if not M2D_DEMO.is_dir():
        pytest.skip("bundled MARE2DEM demo not found")
    eng = ie.engine("mare2dem")
    run = eng.load(M2D_DEMO)
    assert run.rms and run.iterations[-1] == len(run.rms)
    for key, _ in eng.views(run):
        assert eng.render(key, Figure(), run).axes


def test_detect_engine_rejects_unrelated(tmp_path):
    (tmp_path / "notes.txt").write_text("hello")
    assert ie.detect_engine(tmp_path) is None
    assert ie.detect_engine(tmp_path / "missing") is None


def test_elapsed_eta():
    elapsed, eta = ie.elapsed_eta([1, 2, 3, 4], 0.0, 10, now=40.0)
    assert elapsed == 40.0 and eta == pytest.approx(60.0)
    assert ie.elapsed_eta([1], 0.0, 10, now=5.0)[1] is None


# ── run store & presets ─────────────────────────────────────────────────


def test_record_and_list_recent_runs(tmp_path):
    reg = tmp_path / "runs.json"
    a, b = tmp_path / "a", tmp_path / "b"
    a.mkdir()
    b.mkdir()
    store.record_run("occam1d", a, stations=3, final_rms=1.2, path=reg)
    store.record_run("mare2dem", b, status="error", path=reg)
    store.record_run("occam1d", a, stations=3, final_rms=0.9, path=reg)
    runs = store.recent_runs(path=reg)
    assert [Path(r["workdir"]).name for r in runs] == ["a", "b"]
    assert runs[0]["final_rms"] == 0.9
    b.rmdir()
    assert [Path(r["workdir"]).name for r in store.recent_runs(path=reg)] \
        == ["a"]


def test_user_presets_round_trip(tmp_path):
    p = tmp_path / "presets.json"
    store.save_preset("occam2d", "Mine", {"n_layers": 33}, path=p)
    names = [x.name for x in store.presets("occam2d", path=p)]
    assert names[-1] == "Mine" and "Fast preview" in names
    mine = store.presets("occam2d", path=p)[-1]
    assert not mine.builtin and mine.values == {"n_layers": 33}
    store.delete_preset("occam2d", "Mine", path=p)
    assert "Mine" not in [x.name for x in store.presets("occam2d", path=p)]


def test_modem_nan_divergence_is_stopped(monkeypatch, tmp_path):
    """A NaN misfit makes ModEM's line search retry forever (seen with the
    real Mod2DMT): the task stops it and says why."""
    import pycsamt.models.modem as modem_mod
    from pycsamt.models._process import ProcessCancelled

    class _Runner:
        def __init__(self, workdir, config=None):
            pass

        def run(self, *a, on_output=None, cancel=None, **kw):
            for _ in range(40):
                on_output("  CUBICLS: f=  NaN m2=  NaN rms=        NaN")
                if cancel():
                    raise ProcessCancelled("stopped")
            return None

    monkeypatch.setattr(modem_mod, "ModEmRunner", _Runner)
    eng = ie.engine("modem2d")
    build = ie.BuildInfo("modem2d", tmp_path,
                         files={"model": tmp_path / "m", "data": tmp_path / "d",
                                "control": tmp_path / "c"})
    rep = ie.RunReporter()
    with pytest.raises(RuntimeError, match="misfit is NaN"):
        eng.make_task(build, eng.defaults(), "Mod2DMT")(rep)
    assert len(rep.lines) == 30
