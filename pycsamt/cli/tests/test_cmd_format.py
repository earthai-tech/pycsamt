# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for the ``pycsamt format`` command group."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest
from click.testing import CliRunner

from pycsamt.cli import main

_ROOT = Path(__file__).resolve().parents[3]
_OCCAM = _ROOT / "data" / "occam2D"
_MARE = _ROOT / "data" / "mare2dem" / "demo_mt_inversion"


@pytest.fixture
def runner() -> CliRunner:
    return CliRunner()


# ---------------------------------------------------------------------------
# group + help
# ---------------------------------------------------------------------------


def test_group_help(runner):
    result = runner.invoke(main, ["format", "--help"])
    assert result.exit_code == 0
    for sub in ("convert", "detect", "info", "validate"):
        assert sub in result.output


def test_registered_on_root(runner):
    result = runner.invoke(main, ["--help"])
    assert "format" in result.output


# ---------------------------------------------------------------------------
# detect
# ---------------------------------------------------------------------------


def test_detect_npz(runner, tmp_path):
    p = tmp_path / "pred.npz"
    np.savez(p, resistivity=np.ones((6, 9)), x=np.arange(9.0), z=np.arange(6.0))
    result = runner.invoke(main, ["format", "detect", str(p), "-f", "json"])
    assert result.exit_code == 0
    payload = json.loads(result.output)
    assert payload["category"] == "ai_arrays"
    assert payload["target_geometry"] == "grid2d"


def test_detect_unknown_exits_nonzero(runner, tmp_path):
    p = tmp_path / "x.txt"
    p.write_text("nope")
    result = runner.invoke(main, ["format", "detect", str(p)])
    assert result.exit_code == 1


# ---------------------------------------------------------------------------
# convert — AI arrays (no heavy deps)
# ---------------------------------------------------------------------------


def test_convert_npz_grid2d_to_pcsf(runner, tmp_path):
    src = tmp_path / "unet.npz"
    np.savez(
        src,
        resistivity=10 ** np.random.default_rng(0).uniform(1, 3, size=(12, 30)),
        x=np.linspace(0, 600, 30),
        z=np.linspace(0, 240, 12),
    )
    dst = tmp_path / "unet.pcsf"
    result = runner.invoke(
        main, ["format", "convert", str(src), str(dst), "-f", "json"]
    )
    assert result.exit_code == 0, result.output
    report = json.loads(result.output)
    assert report["kind"] == "grid2d"
    assert dst.exists()

    from pycsamt.format.io import read_pcsf

    model = read_pcsf(dst)
    assert model.resistivity.shape == (12, 30)


def test_convert_npz_log10_encoding_recorded(runner, tmp_path):
    src = tmp_path / "net.npz"
    np.savez(
        src,
        log10_rho=np.random.default_rng(1).uniform(1, 3, size=(5, 6, 7)),
        x=np.arange(7.0),
        y=np.arange(6.0),
        z=np.arange(5.0),
    )
    dst = tmp_path / "net.pcsm"
    result = runner.invoke(main, ["format", "convert", str(src), str(dst)])
    assert result.exit_code == 0, result.output

    from pycsamt.format.text import read_pcsm

    model = read_pcsm(dst)
    assert model.resistivity_native_encoding == "log10"
    # canonical resistivity is always linear ohm-m
    assert np.nanmin(model.resistivity) >= 1.0


def test_convert_dry_run_writes_nothing(runner, tmp_path):
    src = tmp_path / "p.npz"
    np.savez(src, resistivity=np.ones((4, 5)), x=np.arange(5.0), z=np.arange(4.0))
    dst = tmp_path / "p.pcsf"
    result = runner.invoke(
        main, ["format", "convert", str(src), str(dst), "--dry-run"]
    )
    assert result.exit_code == 0
    assert not dst.exists()
    assert "dry run" in result.output.lower()


def test_convert_refuses_overwrite_without_flag(runner, tmp_path):
    src = tmp_path / "p.npz"
    np.savez(src, resistivity=np.ones((4, 5)), x=np.arange(5.0), z=np.arange(4.0))
    dst = tmp_path / "p.pcsf"
    dst.write_text("existing")
    result = runner.invoke(main, ["format", "convert", str(src), str(dst)])
    assert result.exit_code != 0
    assert "overwrite" in result.output.lower()


def test_convert_infers_output_name_and_format(runner, tmp_path):
    src = tmp_path / "myrun.npz"
    np.savez(src, resistivity=np.ones((4, 5)), x=np.arange(5.0), z=np.arange(4.0))
    result = runner.invoke(
        main,
        ["format", "convert", str(src), "--to", "pcsm", "-o", str(tmp_path)],
    )
    assert result.exit_code == 0, result.output
    assert (tmp_path / "myrun.pcsm").exists()


# ---------------------------------------------------------------------------
# convert — PCSF <-> PCSM transcode
# ---------------------------------------------------------------------------


def _make_pcsf(path: Path) -> Path:
    from pycsamt.format.adapters.generic import grid2d_to_pcsf
    from pycsamt.format.io import write_pcsf

    model = grid2d_to_pcsf(
        np.geomspace(10, 1000, 40).reshape(5, 8),
        np.arange(8.0),
        np.arange(5.0),
        source_backend="ai",
    )
    return write_pcsf(model, path)


def _make_rich_pcsf(path: Path) -> Path:
    """PCSF exercising stations/topography/history/provenance/geometry
    extras + a NaN in resistivity and an all-NaN uncertainty field."""
    from pycsamt.format.adapters.generic import grid2d_to_pcsf
    from pycsamt.format.io import write_pcsf
    from pycsamt.format.provenance import ModelProvenance
    from pycsamt.format.schema import TopographyPerStation

    rho = np.geomspace(10, 1000, 40).reshape(5, 8)
    rho[0, 0] = np.nan
    names = [f"S{i:02d}" for i in range(8)]
    model = grid2d_to_pcsf(
        rho,
        np.arange(8.0),
        np.arange(5.0),
        origin=(10.0, 20.0),
        azimuth_deg=15.0,
        uncertainty=np.full((5, 8), np.nan),
        station_names=names,
        station_x=np.arange(8.0),
        topography=TopographyPerStation(
            station_id=names, elevation=np.arange(8.0)
        ),
        history={
            "iteration": np.arange(3),
            "rms": np.array([2.0, 1.5, 1.1]),
        },
        provenance=ModelProvenance(framework="pytorch", architecture="unet"),
        description="a rich test model",
        source_backend="ai",
    )
    return write_pcsf(model, path)


def _make_all_nan_pcsf(path: Path) -> Path:
    from pycsamt.format.adapters.generic import grid2d_to_pcsf
    from pycsamt.format.io import write_pcsf

    model = grid2d_to_pcsf(
        np.full((2, 2), np.nan),
        np.arange(2.0),
        np.arange(2.0),
        source_backend="ai",
    )
    return write_pcsf(model, path)


def test_transcode_pcsf_to_pcsm_and_back(runner, tmp_path):
    pcsf = _make_pcsf(tmp_path / "m.pcsf")
    pcsm = tmp_path / "m.pcsm"
    r1 = runner.invoke(main, ["format", "convert", str(pcsf), str(pcsm)])
    assert r1.exit_code == 0, r1.output
    assert pcsm.exists()

    pcsf2 = tmp_path / "roundtrip.pcsf"
    r2 = runner.invoke(main, ["format", "convert", str(pcsm), str(pcsf2)])
    assert r2.exit_code == 0, r2.output

    from pycsamt.format.io import read_pcsf

    a = read_pcsf(pcsf)
    b = read_pcsf(pcsf2)
    np.testing.assert_allclose(a.resistivity, b.resistivity, rtol=1e-6)


def test_transcode_same_format_refused(runner, tmp_path):
    pcsf = _make_pcsf(tmp_path / "m.pcsf")
    result = runner.invoke(
        main, ["format", "convert", str(pcsf), "--to", "pcsf", "-o", str(tmp_path)]
    )
    assert result.exit_code != 0


# ---------------------------------------------------------------------------
# info + validate
# ---------------------------------------------------------------------------


def test_info_json(runner, tmp_path):
    pcsf = _make_pcsf(tmp_path / "m.pcsf")
    result = runner.invoke(main, ["format", "info", str(pcsf), "-f", "json"])
    assert result.exit_code == 0
    report = json.loads(result.output)
    assert report["geometry_kind"] == "grid2d"
    assert report["resistivity"]["shape"] == [5, 8]


def test_info_rejects_non_pcsf(runner, tmp_path):
    p = tmp_path / "x.txt"
    p.write_text("hi")
    result = runner.invoke(main, ["format", "info", str(p)])
    assert result.exit_code != 0


def test_info_text_rich_fields(runner, tmp_path):
    pcsf = _make_rich_pcsf(tmp_path / "rich.pcsf")
    result = runner.invoke(main, ["format", "info", str(pcsf)])
    assert result.exit_code == 0, result.output
    out = result.output
    assert "resistivity NaNs" in out
    assert "stations" in out
    assert "topography" in out
    assert "history keys" in out
    assert "model provenance" in out
    assert "pytorch/unet" in out
    assert "a rich test model" in out
    assert "origin" in out
    assert "azimuth_deg" in out


def test_info_json_rich_fields(runner, tmp_path):
    pcsf = _make_rich_pcsf(tmp_path / "rich.pcsf")
    result = runner.invoke(main, ["format", "info", str(pcsf), "-f", "json"])
    assert result.exit_code == 0, result.output
    report = json.loads(result.output)
    assert report["resistivity"]["n_nan"] == 1
    assert report["uncertainty"]["n_nan"] == 40
    assert "min" not in report["uncertainty"]
    assert report["geometry"]["origin"] == [10.0, 20.0]
    assert report["geometry"]["azimuth_deg"] == 15.0
    assert report["stations"]["n"] == 8
    assert report["topography"]["kind"] == "per_station"
    assert report["topography"]["n"] == 8
    assert report["history"] == {"iteration": [3], "rms": [3]}
    assert report["model_provenance"]["framework"] == "pytorch"


def test_info_all_nan_resistivity_text(runner, tmp_path):
    pcsf = _make_all_nan_pcsf(tmp_path / "allnan.pcsf")
    result = runner.invoke(main, ["format", "info", str(pcsf)])
    assert result.exit_code == 0, result.output
    assert "all NaN" in result.output


def test_info_header_read_error(runner, tmp_path, monkeypatch):
    pcsf = _make_pcsf(tmp_path / "m.pcsf")
    monkeypatch.setattr(
        "pycsamt.format.text.peek_kind",
        lambda *a, **kw: (_ for _ in ()).throw(RuntimeError("corrupt header")),
    )
    result = runner.invoke(main, ["format", "info", str(pcsf)])
    assert result.exit_code != 0
    assert "Cannot read header" in str(result.output) + str(result.exception)


def test_info_model_read_error(runner, tmp_path, monkeypatch):
    pcsf = _make_pcsf(tmp_path / "m.pcsf")
    monkeypatch.setattr(
        "pycsamt.format.text.read_pcsf_or_pcsm",
        lambda *a, **kw: (_ for _ in ()).throw(RuntimeError("bad payload")),
    )
    result = runner.invoke(main, ["format", "info", str(pcsf)])
    assert result.exit_code != 0
    assert "Cannot read" in str(result.output) + str(result.exception)


def test_validate_ok(runner, tmp_path):
    pcsf = _make_pcsf(tmp_path / "m.pcsf")
    result = runner.invoke(main, ["format", "validate", str(pcsf)])
    assert result.exit_code == 0
    assert "VALID" in result.output


def test_validate_corrupt_file(runner, tmp_path):
    bad = tmp_path / "bad.pcsf"
    bad.write_bytes(b"not really hdf5")
    result = runner.invoke(main, ["format", "validate", str(bad), "-f", "json"])
    assert result.exit_code == 1
    payload = json.loads(result.output)
    assert payload["valid"] is False


# ---------------------------------------------------------------------------
# convert — real bundled solver data
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not _OCCAM.exists(), reason="bundled Occam2D data absent")
def test_convert_real_occam2d(runner, tmp_path):
    dst = tmp_path / "occam.pcsf"
    result = runner.invoke(
        main, ["format", "convert", str(_OCCAM), str(dst), "-f", "json"]
    )
    assert result.exit_code == 0, result.output
    report = json.loads(result.output)
    assert report["kind"] == "grid2d"
    assert report["source_backend"] == "occam2d"
    assert report["n_stations"] > 0


@pytest.mark.skipif(not _MARE.exists(), reason="bundled MARE2DEM data absent")
def test_convert_real_mare2dem(runner, tmp_path):
    pytest.importorskip("triangle")
    dst = tmp_path / "mare.pcsf"
    result = runner.invoke(
        main,
        ["format", "convert", str(_MARE), str(dst), "--solver", "mare2dem"],
    )
    assert result.exit_code == 0, result.output
    assert dst.exists()

    from pycsamt.format.io import read_pcsf

    model = read_pcsf(dst)
    assert model.kind == "mesh_unstructured"
