# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.format.convert_engine — the backend-neutral conversion
engine behind ``pycsamt format convert``, the format CLI, and the standalone
converter GUI. No dedicated test file previously existed; heavy solver-based
branches (``build_from_solver``'s occam2d/modem/mare2dem paths) are exercised
by monkeypatching the ``InversionResult`` loaders and adapter functions each
backend calls, following the same technique used elsewhere for optional
dependencies in this repo -- the array-only AI/DL path
(``build_from_ai_arrays``) and the ``mare2dem_mesh``/``resolve_poly`` helpers
are driven end-to-end with real small fixture files instead, since they need
no heavy solver output tree.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import pycsamt.format.convert_engine as ce
from pycsamt.format.detect import SourceKind
from pycsamt.models.mare2dem.iotools.poly import PolyFile, write_poly

# ---------------------------------------------------------------------------
# detect()
# ---------------------------------------------------------------------------


def test_detect_forwards_to_detect_source(monkeypatch, tmp_path):
    calls = {}

    def _fake_detect_source(path, solver_hint=None):
        calls["path"] = path
        calls["solver_hint"] = solver_hint
        return "SENTINEL"

    import pycsamt.format.detect as detect_mod

    monkeypatch.setattr(detect_mod, "detect_source", _fake_detect_source)
    result = ce.detect(tmp_path, solver="modem")
    assert result == "SENTINEL"
    assert calls["path"] == tmp_path
    assert calls["solver_hint"] == "modem"


# ---------------------------------------------------------------------------
# ConvertOptions
# ---------------------------------------------------------------------------


def test_convert_options_defaults():
    opts = ce.ConvertOptions()
    assert opts.iteration is None
    assert opts.created_by == "pycsamt format convert"
    assert opts.station_z_convention == "auto"
    assert opts.air_threshold_ohm_m == 1e8


# ---------------------------------------------------------------------------
# resolve_target
# ---------------------------------------------------------------------------


class TestResolveTarget:
    def test_target_and_to_format_agree(self, tmp_path):
        target = tmp_path / "out.pcsf"
        path, fmt = ce.resolve_target(tmp_path / "src", target, "pcsf", tmp_path)
        assert path == target
        assert fmt == "pcsf"

    def test_target_and_to_format_disagree_raises(self, tmp_path):
        target = tmp_path / "out.pcsm"
        with pytest.raises(ValueError, match="disagree"):
            ce.resolve_target(tmp_path / "src", target, "pcsf", tmp_path)

    def test_target_pcsm_gz_and_to_format_pcsm_agree_on_base(self, tmp_path):
        target = tmp_path / "out.pcsm.gz"
        path, fmt = ce.resolve_target(
            tmp_path / "src", target, "pcsm.gz", tmp_path
        )
        assert fmt == "pcsm.gz"
        assert path == target

    def test_target_without_to_format_infers_pcsm(self, tmp_path):
        target = tmp_path / "out.pcsm"
        path, fmt = ce.resolve_target(tmp_path / "src", target, None, tmp_path)
        assert fmt == "pcsm"
        assert path == target

    def test_target_without_to_format_uninferrable_raises(self, tmp_path):
        target = tmp_path / "out.bin"
        with pytest.raises(ValueError, match="Cannot tell the output format"):
            ce.resolve_target(tmp_path / "src", target, None, tmp_path)

    def test_no_target_defaults_to_pcsf(self, tmp_path):
        src = tmp_path / "model.npz"
        src.write_bytes(b"")
        path, fmt = ce.resolve_target(src, None, None, tmp_path)
        assert fmt == "pcsf"
        assert path == tmp_path / "model.pcsf"

    def test_no_target_explicit_to_format_pcsm(self, tmp_path):
        src = tmp_path / "model.npz"
        src.write_bytes(b"")
        path, fmt = ce.resolve_target(src, None, "pcsm", tmp_path)
        assert fmt == "pcsm"
        assert path == tmp_path / "model.pcsm"


class TestSourceStem:
    def test_strips_pcsm_gz(self):
        assert ce._source_stem(Path("run.pcsm.gz")) == "run"

    def test_strips_pcsf(self):
        assert ce._source_stem(Path("run.pcsf")) == "run"

    def test_strips_npz(self):
        assert ce._source_stem(Path("run.npz")) == "run"

    def test_non_existing_non_matching_uses_name(self):
        assert ce._source_stem(Path("/no/such/run.dat")) == "run.dat"

    def test_existing_file_uses_stem(self, tmp_path):
        p = tmp_path / "run.dat"
        p.write_bytes(b"")
        assert ce._source_stem(p) == "run"


class TestInferFormatFromName:
    @pytest.mark.parametrize(
        "name,expected",
        [
            ("x.pcsm.gz", "pcsm.gz"),
            ("x.pcsm", "pcsm"),
            ("x.pcsf", "pcsf"),
            ("x.dat", None),
        ],
    )
    def test_cases(self, name, expected):
        assert ce._infer_format_from_name(name) == expected


# ---------------------------------------------------------------------------
# build_model dispatch
# ---------------------------------------------------------------------------


class TestBuildModelDispatch:
    def test_pcsf_category_returns_none(self):
        sk = SourceKind(path=Path("dummy"), category="pcsf", is_dir=False)
        assert ce.build_model(sk, ce.ConvertOptions()) is None

    def test_pcsm_category_returns_none(self):
        sk = SourceKind(path=Path("dummy"), category="pcsm", is_dir=False)
        assert ce.build_model(sk, ce.ConvertOptions()) is None

    def test_solver_category_dispatches(self, monkeypatch):
        monkeypatch.setattr(ce, "build_from_solver", lambda sk, opts: "SOLVER")
        sk = SourceKind(path=Path("dummy"), category="solver", is_dir=True, backend="occam2d")
        assert ce.build_model(sk, ce.ConvertOptions()) == "SOLVER"

    def test_ai_arrays_category_dispatches(self, monkeypatch):
        monkeypatch.setattr(ce, "build_from_ai_arrays", lambda sk, opts: "AI")
        sk = SourceKind(path=Path("dummy"), category="ai_arrays", is_dir=False)
        assert ce.build_model(sk, ce.ConvertOptions()) == "AI"

    def test_unknown_category_raises(self):
        sk = SourceKind(path=Path("dummy"), category="mystery", is_dir=False, detail="???")
        with pytest.raises(ValueError, match="Don't know how to convert"):
            ce.build_model(sk, ce.ConvertOptions())


# ---------------------------------------------------------------------------
# build_from_solver
# ---------------------------------------------------------------------------


class TestBuildFromSolver:
    def test_occam2d_backend(self, monkeypatch, tmp_path):
        import pycsamt.format.adapters.occam2d as occam2d_adapter
        import pycsamt.models.occam2d.results as occam2d_results

        captured = {}

        class FakeResult:
            def __init__(self, workdir, iteration=None, verbose=0):
                captured["workdir"] = workdir
                captured["iteration"] = iteration

        def _fake_to_pcsf(result, **kw):
            captured["kwargs"] = kw
            return "OCCAM_MODEL"

        monkeypatch.setattr(occam2d_results, "InversionResult", FakeResult)
        monkeypatch.setattr(occam2d_adapter, "occam2d_to_pcsf", _fake_to_pcsf)

        sk = SourceKind(
            path=Path("dummy"),
            category="solver", is_dir=True, backend="occam2d",
            hints={},
        )
        sk.path = tmp_path
        opts = ce.ConvertOptions(iteration=3)
        result = ce.build_from_solver(sk, opts)
        assert result == "OCCAM_MODEL"
        assert captured["iteration"] == 3
        assert captured["kwargs"]["created_by"] == opts.created_by

    def test_modem_backend(self, monkeypatch, tmp_path):
        import pycsamt.format.adapters.modem3d as modem3d_adapter
        import pycsamt.models.modem.results as modem_results

        captured = {}

        class FakeResult:
            def __init__(self, workdir, load_data=True):
                captured["workdir"] = workdir

        def _fake_to_pcsf(result, **kw):
            captured["kwargs"] = kw
            return "MODEM_MODEL"

        monkeypatch.setattr(modem_results, "InversionResult", FakeResult)
        monkeypatch.setattr(modem3d_adapter, "modem3d_to_pcsf", _fake_to_pcsf)

        sk = SourceKind(path=Path("dummy"), category="solver", is_dir=True, backend="modem")
        sk.path = tmp_path
        opts = ce.ConvertOptions(station_z_convention="positive_down")
        result = ce.build_from_solver(sk, opts)
        assert result == "MODEM_MODEL"
        assert captured["kwargs"]["station_z_convention"] == "positive_down"

    def test_mare2dem_backend(self, monkeypatch, tmp_path):
        import pycsamt.models.mare2dem.results as mare2dem_results

        captured = {}

        class FakeResult:
            def __init__(self, workdir):
                captured["workdir"] = workdir

        def _fake_mesh(poly_path):
            captured["poly_path"] = poly_path
            return "MESH"

        def _fake_to_pcsf(result, mesh, **kw):
            captured["mesh"] = mesh
            captured["kwargs"] = kw
            return "MARE2DEM_MODEL"

        poly_file = tmp_path / "run.poly"
        poly_file.write_text("0\n0\n0\n0\n", encoding="utf-8")

        monkeypatch.setattr(mare2dem_results, "InversionResult", FakeResult)
        monkeypatch.setattr(ce, "mare2dem_mesh", _fake_mesh)
        import pycsamt.format.adapters.mare2dem as mare2dem_adapter

        monkeypatch.setattr(mare2dem_adapter, "mare2dem_to_pcsf", _fake_to_pcsf)

        sk = SourceKind(path=Path("dummy"), category="solver", is_dir=True, backend="mare2dem")
        sk.path = tmp_path
        result = ce.build_from_solver(sk, ce.ConvertOptions())
        assert result == "MARE2DEM_MODEL"
        assert captured["mesh"] == "MESH"
        assert captured["poly_path"] == poly_file

    def test_unsupported_backend_raises(self, tmp_path):
        sk = SourceKind(path=Path("dummy"), category="solver", is_dir=True, backend="unknown_solver")
        sk.path = tmp_path
        with pytest.raises(ValueError, match="Unsupported solver backend"):
            ce.build_from_solver(sk, ce.ConvertOptions())

    def test_uses_parent_dir_when_source_is_file(self, monkeypatch, tmp_path):
        import pycsamt.format.adapters.modem3d as modem3d_adapter
        import pycsamt.models.modem.results as modem_results

        captured = {}

        class FakeResult:
            def __init__(self, workdir, load_data=True):
                captured["workdir"] = workdir

        monkeypatch.setattr(modem_results, "InversionResult", FakeResult)
        monkeypatch.setattr(
            modem3d_adapter, "modem3d_to_pcsf", lambda result, **kw: "M"
        )

        model_file = tmp_path / "m.ws"
        model_file.write_text("", encoding="utf-8")
        sk = SourceKind(path=Path("dummy"), category="solver", is_dir=False, backend="modem")
        sk.path = model_file
        ce.build_from_solver(sk, ce.ConvertOptions())
        assert captured["workdir"] == tmp_path


# ---------------------------------------------------------------------------
# resolve_poly
# ---------------------------------------------------------------------------


class TestResolvePoly:
    def test_opts_poly_wins(self, tmp_path):
        opts = ce.ConvertOptions(poly=Path("/explicit/path.poly"))
        sk = SourceKind(path=Path("dummy"), category="solver", is_dir=True, backend="mare2dem")
        assert ce.resolve_poly(sk, tmp_path, opts) == Path("/explicit/path.poly")

    def test_hint_used_when_no_opts_poly(self, tmp_path):
        sk = SourceKind(
            path=Path("dummy"),
            category="solver", is_dir=True, backend="mare2dem",
            hints={"poly": tmp_path / "hinted.poly"},
        )
        out = ce.resolve_poly(sk, tmp_path, ce.ConvertOptions())
        assert out == tmp_path / "hinted.poly"

    def test_glob_finds_poly_in_workdir(self, tmp_path):
        (tmp_path / "run.poly").write_text("", encoding="utf-8")
        sk = SourceKind(path=Path("dummy"), category="solver", is_dir=True, backend="mare2dem")
        out = ce.resolve_poly(sk, tmp_path, ce.ConvertOptions())
        assert out == tmp_path / "run.poly"

    def test_no_poly_found_raises(self, tmp_path):
        sk = SourceKind(path=Path("dummy"), category="solver", is_dir=True, backend="mare2dem")
        with pytest.raises(ValueError, match="No .poly PSLG found"):
            ce.resolve_poly(sk, tmp_path, ce.ConvertOptions())


# ---------------------------------------------------------------------------
# mare2dem_mesh -- real triangle-backed rebuild
# ---------------------------------------------------------------------------


def test_mare2dem_mesh_rebuilds_trimesh_from_poly(tmp_path):
    pf = PolyFile()
    pf.nodes = np.array(
        [[0.0, 0.0], [10.0, 0.0], [10.0, 10.0], [0.0, 10.0]], dtype=float
    )
    pf.segments = np.array([[1, 2], [2, 3], [3, 4], [4, 1]])
    poly_path = write_poly(pf, tmp_path / "square.poly")

    mesh = ce.mare2dem_mesh(poly_path)
    assert mesh.nodes_m.shape[1] == 2
    assert mesh.triangles.shape[1] == 3
    assert len(mesh.region_ids) == len(mesh.triangles)


def test_mare2dem_mesh_with_regions(tmp_path):
    pf = PolyFile()
    pf.nodes = np.array(
        [[0.0, 0.0], [10.0, 0.0], [10.0, 10.0], [0.0, 10.0]], dtype=float
    )
    pf.segments = np.array([[1, 2], [2, 3], [3, 4], [4, 1]])
    pf.regions = np.array([[5.0, 5.0, 1, 0.01]])
    poly_path = write_poly(pf, tmp_path / "square_regions.poly")

    mesh = ce.mare2dem_mesh(poly_path)
    assert mesh.triangles.shape[0] > 0


# ---------------------------------------------------------------------------
# build_from_ai_arrays
# ---------------------------------------------------------------------------


class TestBuildFromAiArrays:
    def test_grid2d_with_explicit_coordinates(self, tmp_path):
        rho = np.full((3, 4), 100.0)
        x = np.arange(4, dtype=float) * 50.0
        z = np.arange(3, dtype=float) * 20.0
        npz_path = tmp_path / "model.npz"
        np.savez(npz_path, rho=rho, x=x, z=z)

        sk = SourceKind(
            path=Path("dummy"),
            category="ai_arrays", is_dir=False, target_geometry="grid2d",
            hints={"resistivity_key": "rho", "x_key": "x", "z_key": "z"},
        )
        sk.path = npz_path
        model = ce.build_from_ai_arrays(sk, ce.ConvertOptions())
        assert model.kind == "grid2d"
        np.testing.assert_allclose(model.resistivity, rho)

    def test_grid2d_missing_coordinates_falls_back_to_arange(self, tmp_path):
        rho = np.full((3, 4), 50.0)
        npz_path = tmp_path / "model.npz"
        np.savez(npz_path, rho=rho)

        sk = SourceKind(
            path=Path("dummy"),
            category="ai_arrays", is_dir=False, target_geometry="grid2d",
            hints={"resistivity_key": "rho"},
        )
        sk.path = npz_path
        model = ce.build_from_ai_arrays(sk, ce.ConvertOptions())
        assert model.resistivity.shape == (3, 4)

    def test_grid2d_npy_source(self, tmp_path):
        rho = np.full((2, 2), 75.0)
        npy_path = tmp_path / "model.npy"
        np.save(npy_path, rho)

        sk = SourceKind(
            path=Path("dummy"),
            category="ai_arrays", is_dir=False, target_geometry="grid2d",
            hints={},
        )
        sk.path = npy_path
        model = ce.build_from_ai_arrays(sk, ce.ConvertOptions())
        np.testing.assert_allclose(model.resistivity, rho)

    def test_grid3d_with_explicit_coordinates(self, tmp_path):
        rho = np.full((2, 3, 4), 100.0)  # (z, y, x)
        x = np.arange(4, dtype=float)
        y = np.arange(3, dtype=float)
        z = np.arange(2, dtype=float)
        npz_path = tmp_path / "model3d.npz"
        np.savez(npz_path, rho=rho, x=x, y=y, z=z)

        sk = SourceKind(
            path=Path("dummy"),
            category="ai_arrays", is_dir=False, target_geometry="grid3d",
            hints={
                "resistivity_key": "rho", "x_key": "x", "y_key": "y",
                "z_key": "z",
            },
        )
        sk.path = npz_path
        model = ce.build_from_ai_arrays(sk, ce.ConvertOptions())
        assert model.kind == "grid3d"

    def test_grid3d_missing_coordinates_falls_back_to_arange(self, tmp_path):
        rho = np.full((2, 3, 4), 20.0)
        npz_path = tmp_path / "model3d.npz"
        np.savez(npz_path, rho=rho)

        sk = SourceKind(
            path=Path("dummy"),
            category="ai_arrays", is_dir=False, target_geometry="grid3d",
            hints={"resistivity_key": "rho"},
        )
        sk.path = npz_path
        model = ce.build_from_ai_arrays(sk, ce.ConvertOptions())
        assert model.resistivity.shape == (2, 3, 4)

    def test_mesh_unstructured_per_triangle_resistivity(self, tmp_path):
        nodes = np.array(
            [[0.0, 0.0], [10.0, 0.0], [10.0, 10.0], [0.0, 10.0]]
        )
        conn = np.array([[0, 1, 2], [0, 2, 3]])
        rho = np.array([100.0, 200.0])  # per-triangle (n_tri=2, n_node=4)
        npz_path = tmp_path / "mesh_model.npz"
        np.savez(npz_path, nodes=nodes, conn=conn, rho=rho)

        sk = SourceKind(
            path=Path("dummy"),
            category="ai_arrays", is_dir=False,
            target_geometry="mesh_unstructured",
            hints={
                "resistivity_key": "rho", "nodes_key": "nodes",
                "connectivity_key": "conn",
            },
        )
        sk.path = npz_path
        model = ce.build_from_ai_arrays(sk, ce.ConvertOptions())
        assert model.kind == "mesh_unstructured"

    def test_mesh_unstructured_per_node_resistivity(self, tmp_path):
        nodes = np.array(
            [[0.0, 0.0], [10.0, 0.0], [10.0, 10.0], [0.0, 10.0], [5.0, 5.0]]
        )
        conn = np.array([[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]])
        rho = np.array([10.0, 20.0, 30.0, 40.0, 50.0])  # per-node (n_node=5)
        npz_path = tmp_path / "mesh_model_node.npz"
        np.savez(npz_path, nodes=nodes, conn=conn, rho=rho)

        sk = SourceKind(
            path=Path("dummy"),
            category="ai_arrays", is_dir=False,
            target_geometry="mesh_unstructured",
            hints={
                "resistivity_key": "rho", "nodes_key": "nodes",
                "connectivity_key": "conn",
            },
        )
        sk.path = npz_path
        model = ce.build_from_ai_arrays(sk, ce.ConvertOptions())
        assert model.kind == "mesh_unstructured"

    def test_unsupported_geometry_raises(self, tmp_path):
        npz_path = tmp_path / "x.npz"
        np.savez(npz_path, rho=np.zeros((2, 2)))
        sk = SourceKind(
            path=Path("dummy"),
            category="ai_arrays", is_dir=False, target_geometry="unknown_geom",
            hints={"resistivity_key": "rho"},
        )
        sk.path = npz_path
        with pytest.raises(ValueError, match="Unsupported AI geometry"):
            ce.build_from_ai_arrays(sk, ce.ConvertOptions())

    def test_encoding_and_uncertainty_hints_forwarded(self, tmp_path):
        rho = np.full((2, 2), 2.0)  # log10 encoding
        unc = np.full((2, 2), 0.1)
        npz_path = tmp_path / "model_unc.npz"
        np.savez(npz_path, rho=rho, unc=unc)

        sk = SourceKind(
            path=Path("dummy"),
            category="ai_arrays", is_dir=False, target_geometry="grid2d",
            hints={
                "resistivity_key": "rho", "uncertainty_key": "unc",
                "encoding": "log10",
            },
        )
        sk.path = npz_path
        model = ce.build_from_ai_arrays(sk, ce.ConvertOptions())
        np.testing.assert_allclose(model.resistivity, 10.0**rho)
        assert model.uncertainty is not None


# ---------------------------------------------------------------------------
# write_model / model_report -- real PCSFModel round trip
# ---------------------------------------------------------------------------


@pytest.fixture
def real_model():
    from pycsamt.format.adapters.generic import grid2d_to_pcsf

    rho = np.array([[100.0, 110.0], [50.0, 55.0]])
    return grid2d_to_pcsf(
        rho,
        x=np.array([0.0, 100.0]),
        z=np.array([10.0, 50.0]),
        created_by="test-suite",
        description="unit test model",
    )


class TestWriteModel:
    def test_write_pcsf(self, real_model, tmp_path):
        dst = tmp_path / "nested" / "out.pcsf"
        ce.write_model(real_model, dst, "pcsf", log10_view=False)
        assert dst.exists()

    def test_write_pcsm(self, real_model, tmp_path):
        dst = tmp_path / "nested" / "out.pcsm"
        ce.write_model(real_model, dst, "pcsm", log10_view=True)
        assert dst.exists()


class TestModelReport:
    def test_report_contains_expected_fields(self, real_model, tmp_path):
        dst = tmp_path / "out.pcsf"
        ce.write_model(real_model, dst, "pcsf", log10_view=False)
        report = ce.model_report(real_model, dst)
        assert report["kind"] == "grid2d"
        assert report["created_by"] == "test-suite"
        assert report["description"] == "unit test model"
        assert report["resistivity_shape"] == [2, 2]
        assert report["size_bytes"] is not None
        assert "rho_ohm_m" in report
        assert report["rho_ohm_m"]["min"] == pytest.approx(50.0)
        assert report["rho_ohm_m"]["max"] == pytest.approx(110.0)

    def test_report_missing_file_has_none_size(self, real_model, tmp_path):
        report = ce.model_report(real_model, tmp_path / "absent.pcsf")
        assert report["size_bytes"] is None

    def test_report_handles_metadata_dict_exception(self, real_model, tmp_path, monkeypatch):
        monkeypatch.setattr(
            type(real_model),
            "metadata_dict",
            lambda self: (_ for _ in ()).throw(RuntimeError("boom")),
        )
        report = ce.model_report(real_model, tmp_path / "absent.pcsf")
        assert "model_provenance" not in report

    def test_report_with_no_finite_resistivity(self, tmp_path):
        from pycsamt.format.adapters.generic import grid2d_to_pcsf

        rho = np.full((2, 2), np.nan)
        model = grid2d_to_pcsf(
            rho, x=np.array([0.0, 1.0]), z=np.array([0.0, 1.0])
        )
        report = ce.model_report(model, tmp_path / "absent.pcsf")
        assert "rho_ohm_m" not in report
