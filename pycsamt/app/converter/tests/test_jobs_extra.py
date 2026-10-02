# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Additional headless tests for pycsamt.app.converter.jobs, targeting
branches ``test_jobs.py`` doesn't reach: batch-queue edge cases,
PCSF<->PCSM transcoding via ``convert_to_pcsf_job`` directly, every
``validate_pcsf_job``/``info_pcsf_job`` check-failure path, and the
csv-dir/xlsx/las PCBH source kinds plus the XLSX PCPT path.

Some of these (a corrupted-header file, a resistivity-presence mismatch
between a model and its own round trip) are awkward to produce through
the real format stack, so a few tests monkeypatch
``pycsamt.format.text``/``pycsamt.format.io`` entry points with a
minimal duck-typed fake model -- exercising ``jobs.py``'s own branching
logic in isolation, the same way the rest of the suite isolates page
logic from the worker/thread machinery.
"""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.app.converter import jobs

_ROOT = __import__("pathlib").Path(__file__).resolve().parents[4]
_EDI_DIR = _ROOT / "data" / "MT" / "broken-hill" / "edis"


def _edi_dir():
    if not _EDI_DIR.exists() or not any(_EDI_DIR.glob("*.edi")):
        pytest.skip(f"No Broken Hill EDI sample data found at {_EDI_DIR}")
    return _EDI_DIR


def _npz(tmp_path, name="pred.npz", shape=(3, 3), value=20.0):
    p = tmp_path / name
    np.savez(p, resistivity=np.full(shape, value), x=np.arange(shape[1], dtype=float),
              z=np.arange(shape[0], dtype=float))
    return p


# ---------------------------------------------------------------------------
# classify_batch_item / run_batch_job edge cases
# ---------------------------------------------------------------------------


def test_classify_batch_item_nonexistent_path_is_unknown(tmp_path):
    assert jobs.classify_batch_item(tmp_path / "does_not_exist") == "unknown"


def test_run_batch_job_handles_xml_kind(tmp_path):
    edi_dir = _edi_dir()
    xml_dir = tmp_path / "xml"
    jobs.edi_to_xml_job(edi_dir, xml_dir)
    first_xml = sorted(xml_dir.glob("*.xml"))[0]

    out_dir = tmp_path / "out"
    results = jobs.run_batch_job([{"path": first_xml, "kind": "xml"}], out_dir)
    assert results[0]["status"] == "done"
    assert (out_dir / f"{first_xml.stem}.edi").exists()


def test_run_batch_job_infers_kind_when_not_given(tmp_path):
    npz = _npz(tmp_path)
    out_dir = tmp_path / "out"
    results = jobs.run_batch_job([{"path": npz}], out_dir)
    assert results[0]["kind"] == "inversion"
    assert results[0]["status"] == "done"


# ---------------------------------------------------------------------------
# convert_to_pcsf_job on an existing PCSF/PCSM source (transcode branch)
# ---------------------------------------------------------------------------


def test_convert_to_pcsf_job_transcodes_existing_pcsf_to_pcsm(tmp_path):
    npz = _npz(tmp_path)
    pcsf = tmp_path / "m.pcsf"
    jobs.convert_to_pcsf_job(npz, pcsf, "pcsf", tmp_path)

    dst = tmp_path / "m.pcsm"
    report = jobs.convert_to_pcsf_job(pcsf, dst, "pcsm", tmp_path)
    assert dst.exists()
    assert report["transcoded"] is True
    assert report["detected"]["category"] == "pcsf"


def test_convert_to_pcsf_job_refuses_same_format_without_explicit_target(tmp_path):
    npz = _npz(tmp_path)
    pcsf = tmp_path / "m.pcsf"
    jobs.convert_to_pcsf_job(npz, pcsf, "pcsf", tmp_path)

    with pytest.raises(ValueError, match="already PCSF"):
        jobs.convert_to_pcsf_job(pcsf, None, "pcsf", tmp_path)


# ---------------------------------------------------------------------------
# validate_pcsf_job -- every check branch
# ---------------------------------------------------------------------------


def test_validate_pcsf_job_header_failure_on_corrupt_file(tmp_path):
    bad = tmp_path / "corrupt.pcsf"
    bad.write_bytes(b"not a real pcsf file at all")
    result = jobs.validate_pcsf_job(bad)
    assert result["valid"] is False
    steps = {c["step"]: c["ok"] for c in result["checks"]}
    assert steps["header"] is False
    assert "load" not in steps  # short-circuited: kind stayed None


class _FakeModel:
    def __init__(self, kind="grid2d", resistivity=None, topography=None,
                 raise_on_validate=False):
        self.kind = kind
        self.resistivity = resistivity
        self.uncertainty = None
        self.sensitivity = None
        self.resistivity_by_region = None
        self.resistivity_by_node = None
        self.stations = None
        self.topography = topography
        self.geometry = object()
        self.source_backend = "generic"
        self.created_by = "pytest"
        self.created_at = "2026-01-01T00:00:00Z"
        self.description = ""
        self.crs = None
        self.resistivity_native_encoding = None
        self._raise_on_validate = raise_on_validate

    def validate(self):
        if self._raise_on_validate:
            raise ValueError("schema broken")


def test_validate_pcsf_job_load_failure(tmp_path, monkeypatch):
    import pycsamt.format.text as text_mod

    monkeypatch.setattr(text_mod, "peek_kind", lambda path: "grid2d")

    def _boom(path):
        raise OSError("cannot load")

    monkeypatch.setattr(text_mod, "read_pcsf_or_pcsm", _boom)

    result = jobs.validate_pcsf_job(tmp_path / "whatever.pcsf")
    steps = {c["step"]: c["ok"] for c in result["checks"]}
    assert steps["header"] is True
    assert steps["load"] is False
    assert "schema" not in steps
    assert result["valid"] is False


def test_validate_pcsf_job_schema_failure(tmp_path, monkeypatch):
    import pycsamt.format.text as text_mod

    monkeypatch.setattr(text_mod, "peek_kind", lambda path: "grid2d")
    monkeypatch.setattr(
        text_mod, "read_pcsf_or_pcsm",
        lambda path: _FakeModel(raise_on_validate=True),
    )

    result = jobs.validate_pcsf_job(tmp_path / "whatever.pcsf", roundtrip=False)
    steps = {c["step"]: c["ok"] for c in result["checks"]}
    assert steps["load"] is True
    assert steps["schema"] is False
    assert "roundtrip" not in steps


def test_validate_pcsf_job_roundtrip_no_resistivity_array(tmp_path, monkeypatch):
    import pycsamt.format.io as io_mod
    import pycsamt.format.text as text_mod

    monkeypatch.setattr(text_mod, "peek_kind", lambda path: "grid2d")
    monkeypatch.setattr(text_mod, "read_pcsf_or_pcsm", lambda path: _FakeModel())
    monkeypatch.setattr(io_mod, "write_pcsf", lambda model, dst: dst)

    result = jobs.validate_pcsf_job(tmp_path / "whatever.pcsf")
    steps = {c["step"]: (c["ok"], c["detail"]) for c in result["checks"]}
    assert steps["roundtrip"] == (True, "no resistivity array")
    assert result["valid"] is True


def test_validate_pcsf_job_roundtrip_presence_changed(tmp_path, monkeypatch):
    import pycsamt.format.io as io_mod
    import pycsamt.format.text as text_mod

    calls = {"n": 0}

    def _fake_read(path):
        calls["n"] += 1
        if calls["n"] == 1:
            return _FakeModel(resistivity=None)
        return _FakeModel(resistivity=np.ones((2, 2)))

    monkeypatch.setattr(text_mod, "peek_kind", lambda path: "grid2d")
    monkeypatch.setattr(text_mod, "read_pcsf_or_pcsm", _fake_read)
    monkeypatch.setattr(io_mod, "write_pcsf", lambda model, dst: dst)

    result = jobs.validate_pcsf_job(tmp_path / "whatever.pcsf")
    steps = {c["step"]: (c["ok"], c["detail"]) for c in result["checks"]}
    assert steps["roundtrip"] == (False, "resistivity presence changed")
    assert result["valid"] is False


def test_validate_pcsf_job_roundtrip_exception(tmp_path, monkeypatch):
    import pycsamt.format.io as io_mod
    import pycsamt.format.text as text_mod

    monkeypatch.setattr(text_mod, "peek_kind", lambda path: "grid2d")
    monkeypatch.setattr(text_mod, "read_pcsf_or_pcsm", lambda path: _FakeModel())

    def _boom(model, dst):
        raise OSError("disk full")

    monkeypatch.setattr(io_mod, "write_pcsf", _boom)

    result = jobs.validate_pcsf_job(tmp_path / "whatever.pcsf")
    steps = {c["step"]: c["ok"] for c in result["checks"]}
    assert steps["roundtrip"] is False


def test_validate_pcsf_job_pcsm_roundtrip_uses_write_pcsm(tmp_path, monkeypatch):
    import pycsamt.format.text as text_mod

    monkeypatch.setattr(text_mod, "peek_kind", lambda path: "grid2d")
    monkeypatch.setattr(text_mod, "read_pcsf_or_pcsm", lambda path: _FakeModel())
    seen = []
    monkeypatch.setattr(text_mod, "write_pcsm", lambda model, dst, **k: seen.append(dst) or dst)

    result = jobs.validate_pcsf_job(tmp_path / "whatever.pcsm")
    assert seen  # write_pcsm was used, not write_pcsf
    assert result["valid"] is True


# ---------------------------------------------------------------------------
# info_pcsf_job -- topography flag
# ---------------------------------------------------------------------------


def test_info_pcsf_job_reports_topography(tmp_path, monkeypatch):
    import pycsamt.format.text as text_mod

    monkeypatch.setattr(text_mod, "peek_kind", lambda path: "grid2d")
    monkeypatch.setattr(
        text_mod, "read_pcsf_or_pcsm",
        lambda path: _FakeModel(topography=object()),
    )

    placeholder = tmp_path / "whatever.pcsf"
    placeholder.write_bytes(b"")
    info = jobs.info_pcsf_job(placeholder)
    assert info["has_topography"] is True


def test_info_pcsf_job_rejects_wrong_extension(tmp_path):
    bad = tmp_path / "notes.txt"
    bad.write_text("nope")
    with pytest.raises(ValueError):
        jobs.info_pcsf_job(bad)


# ---------------------------------------------------------------------------
# build_pcbh_job -- csv-dir / xlsx / las kinds
# ---------------------------------------------------------------------------


def test_build_pcbh_job_from_csv_directory(tmp_path):
    from pycsamt.format.borehole import (
        Collar,
        CoordinateReferenceSystem,
        LogInterval,
        PCBHBorehole,
        PCBHDocument,
        SurveyStation,
        Trajectory,
        VocabularyEntry,
        write_csv_directory,
    )

    hole = PCBHBorehole(
        id="BH-1",
        name="BH-1",
        kind="mining_exploration",
        status="completed",
        collar=Collar(100.0, 200.0, 50.0),
        total_depth_md=100.0,
        trajectory=Trajectory(
            method="survey",
            north_reference="grid",
            stations=[SurveyStation(0.0, 0.0, 0.0), SurveyStation(100.0, 20.0, 10.0)],
        ),
        interval_logs={"lithology": [LogInterval(0.0, 100.0, code="GRAN", label="Granite")]},
    )
    document = PCBHDocument(
        document_id="test:relational",
        created_at="2026-08-28T00:00:00Z",
        created_by="pytest",
        crs=CoordinateReferenceSystem("EPSG:32629"),
        boreholes=[hole],
        lithologies=[VocabularyEntry("GRAN", "Granite")],
    )
    project_dir = tmp_path / "project"
    write_csv_directory(document, project_dir)

    dst = tmp_path / "out.pcbh.json"
    result = jobs.build_pcbh_job(project_dir, dst, source_kind="csv-dir")
    assert dst.exists()
    assert result["n_boreholes"] == 1
    assert "rows_read" in result


def test_build_pcbh_job_infers_csv_dir_kind_from_directory(tmp_path):
    from pycsamt.format.borehole import (
        Collar,
        CoordinateReferenceSystem,
        PCBHBorehole,
        PCBHDocument,
        write_csv_directory,
    )

    hole = PCBHBorehole(
        id="BH-2", name="BH-2", kind="mining_exploration", status="completed",
        collar=Collar(0.0, 0.0, 0.0), total_depth_md=10.0,
    )
    document = PCBHDocument(
        document_id="test:auto-dir", created_at="2026-08-28T00:00:00Z",
        created_by="pytest", crs=CoordinateReferenceSystem("EPSG:32629"),
        boreholes=[hole],
    )
    project_dir = tmp_path / "project2"
    write_csv_directory(document, project_dir)

    dst = tmp_path / "out2.pcbh.json"
    result = jobs.build_pcbh_job(project_dir, dst)  # source_kind="auto"
    assert result["n_boreholes"] == 1


def test_build_pcbh_job_from_xlsx(tmp_path):
    openpyxl = pytest.importorskip("openpyxl")
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "log"
    ws.append(["Hole ID", "Rock name", "From_m", "To_m", "Resistivity"])
    ws.append(["BH1", "granodiorite", 0, 10.5, 250])
    ws.append(["BH1", "hornblende", 10.5, 22, 80])
    src = tmp_path / "log.xlsx"
    wb.save(src)

    dst = tmp_path / "out.pcbh.json"
    result = jobs.build_pcbh_job(src, dst)  # auto -> xlsx via suffix
    assert dst.exists()
    assert result["n_boreholes"] == 1


_LAS_TEXT = """~VERSION INFORMATION
 VERS. 2.0: CWLS LOG ASCII STANDARD
 WRAP. NO: ONE LINE PER DEPTH STEP
~WELL INFORMATION
 STRT.M 0: START
 STOP.M 2: STOP
 STEP.M 1: STEP
 NULL. -999.25: NULL
 WELL. TEST-1: WELL
~CURVE INFORMATION
 DEPT.M: MEASURED DEPTH
 RESD.OHMM: DEEP RESISTIVITY
~A DEPT RESD
0 10
1 -999.25
2 100
"""


def test_build_pcbh_job_from_las_with_collar(tmp_path):
    src = tmp_path / "hole.las"
    src.write_text(_LAS_TEXT, encoding="utf-8")
    dst = tmp_path / "out.pcbh.json"
    result = jobs.build_pcbh_job(
        src, dst,
        source_kind="las",
        collar_id="BH-LAS-1",
        x=1.0, y=2.0, z=3.0,
        crs_horizontal="LOCAL:test",
        document_id="doc-las",
    )
    assert dst.exists()
    assert result["document_id"] == "doc-las"
    assert result["n_boreholes"] == 1


# ---------------------------------------------------------------------------
# build_pcpt_job -- xlsx source
# ---------------------------------------------------------------------------


def test_build_pcpt_job_from_xlsx(tmp_path):
    openpyxl = pytest.importorskip("openpyxl")
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.append(["Target", "Easting", "Northing", "Elev"])
    ws.append(["P1", 350000, 3000000, 1200])
    ws.append(["P2", 350120, 3000050, 1210])
    src = tmp_path / "targets.xlsx"
    wb.save(src)

    dst = tmp_path / "targets.pcpt.json"
    result = jobs.build_pcpt_job(src, dst, crs="EPSG:32648")
    assert dst.exists()
    assert result["n_points"] == 2
    assert result["crs"] == "EPSG:32648"
