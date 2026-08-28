# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for :mod:`pycsamt.format.detect`."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pycsamt.format.detect import SourceKind, describe_source, detect_source

_ROOT = Path(__file__).resolve().parents[3]


# ---------------------------------------------------------------------------
# error handling
# ---------------------------------------------------------------------------


def test_missing_path_raises():
    with pytest.raises(FileNotFoundError):
        detect_source("does/not/exist.pcsf")


def test_bad_solver_hint_raises(tmp_path):
    (tmp_path / "x").mkdir()
    with pytest.raises(ValueError):
        detect_source(tmp_path / "x", solver_hint="nope")


def test_unknown_file(tmp_path):
    p = tmp_path / "notes.txt"
    p.write_text("hello")
    sk = detect_source(p)
    assert sk.category == "unknown"
    assert not sk.convertible


# ---------------------------------------------------------------------------
# solver directory fingerprints (synthetic signature files)
# ---------------------------------------------------------------------------


def test_detect_occam_dir_by_iter_suffix(tmp_path):
    (tmp_path / "ITER12.iter").touch()
    (tmp_path / "run.resp").touch()
    sk = detect_source(tmp_path)
    assert sk.category == "solver"
    assert sk.backend == "occam2d"
    assert sk.target_geometry == "grid2d"
    assert Path(sk.hints["iteration_files"][0]).name == "ITER12.iter"


def test_detect_occam_dir_by_name(tmp_path):
    (tmp_path / "Occam2DMesh").touch()
    (tmp_path / "Occam2DModel").touch()
    (tmp_path / "OccamStartup").touch()
    sk = detect_source(tmp_path)
    assert sk.backend == "occam2d"


def test_detect_modem_dir(tmp_path):
    (tmp_path / "Modular_NLCG.log").touch()
    (tmp_path / "run.rho").touch()
    (tmp_path / "run.dat").touch()
    sk = detect_source(tmp_path)
    assert sk.backend == "modem"
    assert sk.target_geometry == "grid3d"


def test_detect_mare_dir(tmp_path):
    (tmp_path / "demo.poly").touch()
    (tmp_path / "demo.0.resistivity").touch()
    (tmp_path / "mare2dem.settings").touch()
    sk = detect_source(tmp_path)
    assert sk.backend == "mare2dem"
    assert sk.target_geometry == "mesh_unstructured"
    assert Path(sk.hints["poly"]).name == "demo.poly"


def test_solver_hint_forces_backend(tmp_path):
    (tmp_path / "whatever.dat").touch()
    sk = detect_source(tmp_path, solver_hint="modem")
    assert sk.backend == "modem"
    assert "Forced" in sk.detail


def test_ambiguous_dir_is_medium_confidence(tmp_path):
    # both a mare poly/resistivity pair and an occam .iter file
    (tmp_path / "demo.poly").touch()
    (tmp_path / "demo.0.resistivity").touch()
    (tmp_path / "ITER03.iter").touch()
    sk = detect_source(tmp_path)
    assert sk.confidence == "medium"
    assert sk.backend in {"mare2dem", "occam2d"}


def test_single_file_routes_to_parent_dir(tmp_path):
    (tmp_path / "ITER12.iter").touch()
    sk = detect_source(tmp_path / "ITER12.iter")
    assert sk.category == "solver"
    assert sk.backend == "occam2d"
    assert sk.is_dir is True


def test_empty_dir_unknown(tmp_path):
    sk = detect_source(tmp_path)
    assert sk.category == "unknown"


# ---------------------------------------------------------------------------
# AI array bundles
# ---------------------------------------------------------------------------


def test_detect_npz_grid2d(tmp_path):
    p = tmp_path / "pred.npz"
    np.savez(p, resistivity=np.ones((10, 20)), x=np.arange(20.0), z=np.arange(10.0))
    sk = detect_source(p)
    assert sk.category == "ai_arrays"
    assert sk.target_geometry == "grid2d"
    assert sk.hints["resistivity_key"] == "resistivity"
    assert sk.hints["encoding"] == "linear"
    assert sk.confidence == "high"


def test_detect_npz_grid3d_log10_alias(tmp_path):
    p = tmp_path / "pred.npz"
    np.savez(
        p,
        log10_rho=np.ones((4, 5, 6)),
        x=np.arange(6.0),
        y=np.arange(5.0),
        z=np.arange(4.0),
    )
    sk = detect_source(p)
    assert sk.target_geometry == "grid3d"
    assert sk.hints["encoding"] == "log10"


def test_detect_npz_mesh(tmp_path):
    p = tmp_path / "gcn.npz"
    np.savez(
        p,
        rho=np.ones(12),
        nodes=np.zeros((9, 2)),
        triangles=np.zeros((12, 3), dtype=int),
    )
    sk = detect_source(p)
    assert sk.target_geometry == "mesh_unstructured"
    assert sk.hints["nodes_key"] == "nodes"
    assert sk.hints["connectivity_key"] == "triangles"


def test_detect_npz_no_resistivity_is_unknown(tmp_path):
    p = tmp_path / "bad.npz"
    np.savez(p, foo=np.ones(3), bar=np.ones(3))
    sk = detect_source(p)
    assert sk.category == "unknown"


def test_detect_npz_no_coords_medium_confidence(tmp_path):
    p = tmp_path / "pred.npz"
    np.savez(p, resistivity=np.ones((10, 20)))
    sk = detect_source(p)
    assert sk.target_geometry == "grid2d"
    assert sk.confidence == "medium"
    assert sk.hints["synthetic_coords"] is True


def test_detect_bare_npy(tmp_path):
    p = tmp_path / "m.npy"
    np.save(p, np.ones((7, 9)))
    sk = detect_source(p)
    assert sk.category == "ai_arrays"
    assert sk.target_geometry == "grid2d"
    assert sk.confidence == "low"


# ---------------------------------------------------------------------------
# real bundled solver data
# ---------------------------------------------------------------------------


@pytest.mark.skipif(
    not (_ROOT / "data" / "occam2D").exists(), reason="bundled data absent"
)
def test_detect_real_occam2d():
    sk = detect_source(_ROOT / "data" / "occam2D")
    assert sk.backend == "occam2d"
    assert sk.convertible


@pytest.mark.skipif(
    not (_ROOT / "data" / "mare2dem" / "demo_mt_inversion").exists(),
    reason="bundled data absent",
)
def test_detect_real_mare2dem():
    sk = detect_source(_ROOT / "data" / "mare2dem" / "demo_mt_inversion")
    assert sk.backend == "mare2dem"
    assert "poly" in sk.hints


# ---------------------------------------------------------------------------
# PCSF / PCSM peek
# ---------------------------------------------------------------------------


def test_detect_pcsf_and_pcsm(tmp_path):
    from pycsamt.format.adapters.generic import grid2d_to_pcsf
    from pycsamt.format.io import write_pcsf
    from pycsamt.format.text import write_pcsm

    model = grid2d_to_pcsf(
        np.ones((5, 8)), np.arange(8.0), np.arange(5.0), source_backend="ai"
    )
    pcsf = write_pcsf(model, tmp_path / "m.pcsf")
    pcsm = write_pcsm(model, tmp_path / "m.pcsm")

    sk_f = detect_source(pcsf)
    assert sk_f.category == "pcsf"
    assert sk_f.geometry == "grid2d"

    sk_m = detect_source(pcsm)
    assert sk_m.category == "pcsm"
    assert sk_m.geometry == "grid2d"


def test_dir_with_single_pcsf(tmp_path):
    from pycsamt.format.adapters.generic import grid2d_to_pcsf
    from pycsamt.format.io import write_pcsf

    model = grid2d_to_pcsf(
        np.ones((5, 8)), np.arange(8.0), np.arange(5.0), source_backend="ai"
    )
    write_pcsf(model, tmp_path / "only.pcsf")
    sk = detect_source(tmp_path)
    assert sk.category == "pcsf"


def test_sourcekind_to_dict_and_describe(tmp_path):
    (tmp_path / "ITER01.iter").touch()
    sk = detect_source(tmp_path)
    d = sk.to_dict()
    assert d["category"] == "solver"
    assert isinstance(d["hints"], dict)
    text = describe_source(sk)
    assert "backend" in text
    assert isinstance(sk, SourceKind)
