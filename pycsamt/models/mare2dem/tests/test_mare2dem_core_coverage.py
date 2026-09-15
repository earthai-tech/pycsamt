"""Focused coverage for the high-level MARE2DEM compatibility API."""

from __future__ import annotations

import importlib
from types import SimpleNamespace

import numpy as np
import pytest

from pycsamt.models.mare2dem.base import Mare2DEMBase
from pycsamt.models.mare2dem.builder import InputBuilder
from pycsamt.models.mare2dem.config import Mare2DEMConfig
from pycsamt.models.mare2dem.data import EMData
from pycsamt.models.mare2dem.import_topo import TopoConfig, import_topo
from pycsamt.models.mare2dem.iotools.emdata import EMDataFile
from pycsamt.models.mare2dem.mesh import PolyMesh, ResistivityModel
from pycsamt.models.mare2dem.noise import NoiseConfig, add_synthetic_noise, make_synthetic_data
from pycsamt.models.mare2dem.runner import Mare2DEMRunner, _resolve_binary


def test_base_repr_and_string():
    obj = Mare2DEMBase(verbose=True)
    obj.answer = 42
    assert obj.verbose == 1
    assert "answer=42" in repr(obj)
    assert str(obj) == repr(obj)


def test_emdata_compatibility_wrapper(tmp_path):
    empty = EMData()
    assert empty.header == {} and empty.data is None and empty.n_data == 0
    assert "n_data=0" in repr(empty)
    with pytest.raises(RuntimeError, match="nothing to write"):
        empty.write(tmp_path / "no.emdata")

    src = tmp_path / "source.emdata"
    src.write_text("Format: EMData_2.3\n!note\n# Data: 1\n1 1 1 1 2 0.1\n")
    loaded = EMData(src)
    assert loaded.header["Format"] == "EMData_2.3"
    assert loaded.header["Comment"] == "note"
    assert loaded.data.shape == (1, 6)
    assert loaded.write(tmp_path / "copy.emdata").exists()


def test_mesh_compatibility_wrappers(tmp_path):
    empty_model = ResistivityModel()
    assert empty_model.header == {} and empty_model.n_elements == 0
    assert empty_model.n_nodes == 0
    assert empty_model.write(tmp_path / "empty").exists()

    half = ResistivityModel.halfspace(2.0, n_nodes=5)
    assert half.n_elements == 1
    np.testing.assert_allclose(half._rf.resistivity, [[100.0]])
    model_path = half.write(tmp_path / "half.resistivity")
    loaded = ResistivityModel(model_path)
    assert loaded.header["num_regions"] == 1
    assert "n_elements=1" in repr(loaded)

    empty_mesh = PolyMesh()
    assert empty_mesh.vertices is empty_mesh.segments is empty_mesh.holes is None
    assert empty_mesh.write(tmp_path / "empty.poly").exists()
    poly = tmp_path / "triangle.poly"
    poly.write_text("3 2 0 0\n1 0 0\n2 1 0\n3 0 1\n3 0\n1 1 2\n2 2 3\n3 3 1\n0\n0\n")
    mesh = PolyMesh(poly)
    assert mesh.vertices.shape == (3, 2) and mesh.segments.shape == (3, 2)
    assert mesh.holes.shape == (0, 2)
    assert "n_vertices=3" in repr(mesh)
    assert mesh.write(tmp_path / "copy.poly").exists()


def test_import_topo_distance_modes_and_errors(tmp_path):
    with pytest.raises(FileNotFoundError):
        import_topo(TopoConfig(topo_file=tmp_path / "missing"))
    one = tmp_path / "one.txt"
    one.write_text("2 30\n")
    prof = import_topo(TopoConfig(one, col_distance_km=1, col_depth_m=2))
    np.testing.assert_allclose(prof.y_topo, [2000])
    np.testing.assert_allclose(prof.z_topo, [30])

    many = tmp_path / "many.txt"
    many.write_text("0 10\n100 -5\n")
    prof = import_topo(TopoConfig(many, col_distance_m=1, col_elevation_m=2))
    np.testing.assert_allclose(prof.z_topo, [-10, 5])
    with pytest.raises(ValueError, match="depth_m or col_elevation"):
        import_topo(TopoConfig(many, col_distance_m=1))
    with pytest.raises(ValueError, match="distance column"):
        import_topo(TopoConfig(many, col_depth_m=2))


def test_import_topo_geographic_projection(tmp_path, monkeypatch):
    geo = tmp_path / "geo.txt"
    geo.write_text("-1 5 10\n-2 6 20\n")
    topo_module = importlib.import_module("pycsamt.models.mare2dem.import_topo")
    monkeypatch.setattr(
        topo_module,
        "lonlat_to_utm",
        lambda lon, lat, **kw: (np.array([100, 110]), np.array([200, 220]), 30, "N"),
    )
    monkeypatch.setattr(topo_module, "get_line_orientation", lambda n, e: 90)
    cfg = TopoConfig(geo, col_longitude=1, col_latitude=2, col_depth_m=3)
    prof = import_topo(cfg, utm_north0=200, utm_east0=100, utm_theta=0)
    np.testing.assert_allclose(prof.y_topo, [0, 10])
    np.testing.assert_allclose(prof.northings, [200, 220])
    monkeypatch.setattr(topo_module, "get_line_orientation", lambda n, e: 0)
    with pytest.raises(ValueError, match="orientation"):
        import_topo(cfg, utm_theta=0)


def test_noise_all_code_families_and_file_workflow(tmp_path):
    codes = [123, 103, 104, 133, 134, 27, 37, 21, 31, 22, 32, 999]
    rows = [[c, 1, 1, 1, 0, 0, float(i + 1), 0] for i, c in enumerate(codes)]
    em = EMDataFile(is_response=True, data=np.asarray(rows, dtype=float))
    noisy = add_synthetic_noise(em, NoiseConfig(mt_abs_noise_tipper=0), seed=7)
    assert noisy.data.shape == (len(codes) - 1, 6)
    assert not noisy.is_response and noisy.format == "EMData_2.3"
    with pytest.raises(ValueError, match="8-column"):
        add_synthetic_noise(EMDataFile(), NoiseConfig())

    response = tmp_path / "response.resp"
    from pycsamt.models.mare2dem.iotools.emdata import write_emdata
    write_emdata(em, response)
    result = make_synthetic_data(response, tmp_path / "synthetic.emdata", NoiseConfig(mt_abs_noise_tipper=0), seed=7)
    assert result.n_data == len(codes) - 1


def test_builder_all_source_branches(tmp_path, monkeypatch):
    cfg = Mare2DEMConfig(initial_rho=10, max_iterations=4)
    builder = InputBuilder(cfg, verbose=True)
    assert builder.write_settings(tmp_path, filename="custom.settings", tolerance=2).exists()
    model = builder.write_resistivity(tmp_path, filename="custom.resistivity", poly_file="x.poly")
    assert model.exists()

    em = EMDataFile(data=np.empty((0, 6)))
    built = builder.build(em, tmp_path / "object")
    assert all(path.exists() for path in built.values())

    built = builder.build(None, tmp_path / "none")
    assert not built["data"].exists()
    source = tmp_path / "source.emdata"
    source.write_text("Format: EMData_2.3\n# Data: 0\n")
    built = builder.build(source, tmp_path / "copied", data_filename="renamed.emdata")
    assert built["data"].exists()

    fake = SimpleNamespace()
    monkeypatch.setattr("pycsamt.models.mare2dem.builder.make_data_file", lambda path, topo, **kw: fake)
    builder.build(None, tmp_path / "generated", mt=SimpleNamespace())
    assert builder._em is fake


def test_runner_commands_resolution_and_mocked_run(tmp_path, monkeypatch):
    cfg = Mare2DEMConfig(binary="MARE2DEM", use_mpi=True, n_procs=3)
    runner = Mare2DEMRunner(tmp_path, config=cfg, verbose=True)
    assert "-np 3" in runner.command("model.resistivity")
    assert "-np" not in runner.command("model", use_mpi=False)

    monkeypatch.setattr("pycsamt.models.mare2dem.runner.shutil.which", lambda name: "found")
    assert _resolve_binary("MARE2DEM", runner._source_mgr).name == "MARE2DEM"
    monkeypatch.setattr("pycsamt.models.mare2dem.runner.shutil.which", lambda name: None)
    monkeypatch.setattr(runner._source_mgr, "resolve_binary", lambda: None)
    with pytest.raises(FileNotFoundError, match="not found"):
        runner.run("model")

    monkeypatch.setattr(runner._source_mgr, "resolve_binary", lambda: tmp_path / "bin")
    called = {}
    class Proc:
        def check_returncode(self):
            called["checked"] = True
    def fake_run(cmd, **kwargs):
        called["cmd"] = cmd
        called.update(kwargs)
        return Proc()
    monkeypatch.setattr("pycsamt.models.mare2dem.runner.subprocess.run", fake_run)
    assert runner.run("model.resistivity", extra_args=["--x"], timeout=2, load_result=False) is None
    assert called["checked"] and called["cmd"][-2:] == ["model", "--x"]
