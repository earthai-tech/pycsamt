# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Headless tests for :mod:`pycsamt.app.converter.jobs`.

None of these import PySide6 or construct a ``QApplication`` -- every
job function is plain Python, so it's tested directly here the same way
a CLI command's underlying logic would be. Real bundled sample data
(Broken Hill EDIs, the Occam2D demo run) is used where available and
skipped automatically when absent; PCGL/PCGS/PCPT/PCBH-CSV cases build
tiny self-contained fixtures via ``tmp_path``.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from pycsamt.app.converter import jobs

_ROOT = Path(__file__).resolve().parents[4]
_EDI_DIR = _ROOT / "data" / "MT" / "broken-hill" / "edis"
_OCCAM = _ROOT / "data" / "occam2D"


def _edi_dir() -> Path:
    if not _EDI_DIR.exists() or not any(_EDI_DIR.glob("*.edi")):
        pytest.skip(f"No Broken Hill EDI sample data found at {_EDI_DIR}")
    return _EDI_DIR


# ---------------------------------------------------------------------------
# detect_source_job / convert_to_pcsf_job
# ---------------------------------------------------------------------------


def test_detect_source_job_npz(tmp_path):
    p = tmp_path / "pred.npz"
    np.savez(p, resistivity=np.ones((6, 9)), x=np.arange(9.0), z=np.arange(6.0))
    result = jobs.detect_source_job(p)
    assert result["category"] == "ai_arrays"
    assert result["target_geometry"] == "grid2d"


def test_convert_to_pcsf_job_npz_grid2d(tmp_path):
    p = tmp_path / "pred.npz"
    np.savez(p, resistivity=np.full((4, 5), 100.0), x=np.arange(5.0), z=np.arange(4.0))
    dst = tmp_path / "out.pcsf"
    report = jobs.convert_to_pcsf_job(p, dst, "pcsf", tmp_path)
    assert dst.exists()
    assert report["resistivity_shape"] == [4, 5]
    assert report["detected"]["category"] == "ai_arrays"


def test_convert_to_pcsf_job_refuses_overwrite(tmp_path):
    p = tmp_path / "pred.npz"
    np.savez(p, resistivity=np.ones((3, 3)), x=np.arange(3.0), z=np.arange(3.0))
    dst = tmp_path / "out.pcsf"
    jobs.convert_to_pcsf_job(p, dst, "pcsf", tmp_path)
    with pytest.raises(FileExistsError):
        jobs.convert_to_pcsf_job(p, dst, "pcsf", tmp_path)
    # overwrite=True must succeed
    report = jobs.convert_to_pcsf_job(p, dst, "pcsf", tmp_path, overwrite=True)
    assert report["file"].endswith("out.pcsf")


@pytest.mark.skipif(not _OCCAM.exists(), reason="bundled Occam2D data absent")
def test_convert_to_pcsf_job_real_occam2d(tmp_path):
    dst = tmp_path / "occam.pcsf"
    report = jobs.convert_to_pcsf_job(_OCCAM, dst, "pcsf", tmp_path)
    assert dst.exists()
    assert report["detected"]["backend"] == "occam2d"


# ---------------------------------------------------------------------------
# transcode_pcsf_pcsm_job
# ---------------------------------------------------------------------------


def test_transcode_pcsf_to_pcsm_and_back(tmp_path):
    npz = tmp_path / "pred.npz"
    np.savez(npz, resistivity=np.full((3, 4), 50.0), x=np.arange(4.0), z=np.arange(3.0))
    pcsf = tmp_path / "m.pcsf"
    jobs.convert_to_pcsf_job(npz, pcsf, "pcsf", tmp_path)

    pcsm = tmp_path / "m.pcsm"
    report = jobs.transcode_pcsf_pcsm_job(pcsf, pcsm, "pcsm")
    assert pcsm.exists()
    assert report["resistivity_shape"] == [3, 4]

    pcsf2 = tmp_path / "m2.pcsf"
    report2 = jobs.transcode_pcsf_pcsm_job(pcsm, pcsf2, "pcsf")
    assert pcsf2.exists()
    assert report2["resistivity_shape"] == [3, 4]


# ---------------------------------------------------------------------------
# validate_pcsf_job / info_pcsf_job
# ---------------------------------------------------------------------------


def test_validate_and_info_job(tmp_path):
    npz = tmp_path / "pred.npz"
    np.savez(npz, resistivity=np.full((3, 3), 10.0), x=np.arange(3.0), z=np.arange(3.0))
    pcsf = tmp_path / "m.pcsf"
    jobs.convert_to_pcsf_job(npz, pcsf, "pcsf", tmp_path)

    validation = jobs.validate_pcsf_job(pcsf)
    assert validation["valid"] is True
    assert {c["step"] for c in validation["checks"]} >= {"header", "load", "schema", "roundtrip"}

    info = jobs.info_pcsf_job(pcsf)
    assert info["resistivity"]["shape"] == [3, 3]


# ---------------------------------------------------------------------------
# EDI <-> EMTF-XML
# ---------------------------------------------------------------------------


def test_edi_to_xml_and_back_round_trip(tmp_path):
    edi_dir = _edi_dir()
    xml_dir = tmp_path / "xml"
    edi_back = tmp_path / "edi_back"

    seen = []
    written_xml = jobs.edi_to_xml_job(
        edi_dir, xml_dir, progress=lambda i, n, name: seen.append((i, n, name))
    )
    n_edi = len(list(edi_dir.glob("*.edi")))
    assert len(written_xml) == n_edi
    assert len(seen) == n_edi

    written_edi = jobs.xml_to_edi_job(xml_dir, edi_back)
    assert len(written_edi) == n_edi

    original_names = {p.stem for p in edi_dir.glob("*.edi")}
    round_tripped_names = {p.stem for p in edi_back.glob("*.edi")}
    assert original_names == round_tripped_names


def test_edi_to_xml_missing_source_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        jobs.edi_to_xml_job(tmp_path, tmp_path / "out")


# ---------------------------------------------------------------------------
# Batch queue
# ---------------------------------------------------------------------------


def test_classify_batch_item(tmp_path):
    edi_dir = _edi_dir()
    first_edi = sorted(edi_dir.glob("*.edi"))[0]
    assert jobs.classify_batch_item(first_edi) == "edi"

    npz = tmp_path / "pred.npz"
    np.savez(npz, resistivity=np.ones((3, 3)), x=np.arange(3.0), z=np.arange(3.0))
    assert jobs.classify_batch_item(npz) == "inversion"

    junk = tmp_path / "notes.txt"
    junk.write_text("not a source")
    assert jobs.classify_batch_item(junk) == "unknown"


def test_run_batch_job_mixed_queue_isolates_failures(tmp_path):
    edi_dir = _edi_dir()
    first_edi = sorted(edi_dir.glob("*.edi"))[0]

    npz = tmp_path / "pred.npz"
    np.savez(npz, resistivity=np.full((3, 3), 20.0), x=np.arange(3.0), z=np.arange(3.0))

    junk = tmp_path / "notes.txt"
    junk.write_text("not a source")

    items = [
        {"path": first_edi, "kind": "edi"},
        {"path": npz, "kind": "inversion"},
        {"path": junk, "kind": "unknown"},
    ]
    out_dir = tmp_path / "out"
    seen = []
    results = jobs.run_batch_job(
        items, out_dir, progress=lambda i, n, name: seen.append((i, n, name))
    )
    assert len(results) == 3
    assert len(seen) == 3
    statuses = {r["path"]: r["status"] for r in results}
    assert statuses[str(first_edi)] == "done"
    assert statuses[str(npz)] == "done"
    assert statuses[str(junk)] == "skipped"
    assert (out_dir / f"{first_edi.stem}.xml").exists()
    assert (out_dir / "pred.pcsf").exists()


def test_run_batch_job_records_error_without_aborting_queue(tmp_path):
    npz = tmp_path / "pred.npz"
    np.savez(npz, resistivity=np.full((3, 3), 20.0), x=np.arange(3.0), z=np.arange(3.0))
    dst = tmp_path / "out" / "pred.pcsf"
    dst.parent.mkdir(parents=True)
    dst.write_bytes(b"placeholder")

    items = [
        {"path": npz, "kind": "inversion"},  # will fail: output already exists
    ]
    results = jobs.run_batch_job(items, tmp_path / "out", overwrite=False)
    assert results[0]["status"] == "error"


# ---------------------------------------------------------------------------
# PCGL / PCGS / PCPT builders
# ---------------------------------------------------------------------------


def test_build_pcgl_job(tmp_path):
    src = tmp_path / "units.csv"
    src.write_text("name,rho_min,rho_max\nA,1,10\nB,10,100\n")
    dst = tmp_path / "units.pcgl.json"
    result = jobs.build_pcgl_job(src, dst, title="Test")
    assert dst.exists()
    assert result["n_entries"] == 2
    payload = json.loads(dst.read_text())
    assert payload["title"] == "Test"


def test_build_pcgs_job_requires_a_source(tmp_path):
    with pytest.raises(ValueError):
        jobs.build_pcgs_job(tmp_path / "out.pcgs.json")


def test_build_pcgs_job_from_planar(tmp_path):
    planar = tmp_path / "planar.csv"
    planar.write_text(
        "x,kind,strike_deg,dip_deg,dip_direction_deg\n100,bedding,45,60,135\n"
    )
    dst = tmp_path / "s.pcgs.json"
    result = jobs.build_pcgs_job(dst, planar_path=planar)
    assert dst.exists()
    assert result["n_planar"] == 1


def test_build_pcpt_job(tmp_path):
    src = tmp_path / "targets.csv"
    src.write_text("name,x,y,z\nA,1,2,3\n")
    dst = tmp_path / "targets.pcpt.json"
    result = jobs.build_pcpt_job(src, dst)
    assert dst.exists()
    assert result["n_points"] == 1


# ---------------------------------------------------------------------------
# PCBH builder
# ---------------------------------------------------------------------------


def test_build_pcbh_job_from_csv(tmp_path):
    src = tmp_path / "boreholes.csv"
    src.write_text(
        "borehole_id,x,y,z,crs,kind,status,total_depth_md,from_md,to_md,"
        "lithology\nBH-1,100,200,10,EPSG:4326,water,completed,20,0,10,Sand\n"
    )
    dst = tmp_path / "boreholes.pcbh.json"
    result = jobs.build_pcbh_job(src, dst)
    assert dst.exists()
    assert result["n_boreholes"] == 1


def test_build_pcbh_job_las_needs_collar(tmp_path):
    fake_las = tmp_path / "hole.las"
    fake_las.write_text("~V\nVERS. 2.0 :\n")
    dst = tmp_path / "hole.pcbh.json"
    with pytest.raises(ValueError):
        jobs.build_pcbh_job(fake_las, dst)
