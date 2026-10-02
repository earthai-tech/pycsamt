# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Targeted coverage tests for :mod:`pycsamt.transformers.jedi`.

Exercises branches that the smoke / true-file / topo test files do not
reach: dead-simple utility classes, the AVG rho/phase-only fallback path,
tipper population, magnetic-only channel sets, ``_apply_topo_coords``
edge cases, ``_enrich_from_avg_info``, and the J-side mirrors of the
same logic.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from pycsamt.core.base import TFBundle
from pycsamt.seg.edi import EDIFile
from pycsamt.transformers import jedi as tr
from pycsamt.zonge.avg import BaseAVG

# --------------------------------------------------------------------- #
# _HeadStub — small, otherwise-dead utility classes
# --------------------------------------------------------------------- #


def test_avgtoedi_headstub_sets_attrs():
    stub = tr.AVGtoEDI._HeadStub(
        "S1", lat=1.0, lon=2.0, elev=3.0, empty=1e32, extra="x"
    )
    assert stub.dataid == "S1"
    assert stub.lat == 1.0
    assert stub.long == 2.0
    assert stub.elev == 3.0
    assert stub.empty == 1e32
    assert stub.extra == "x"


def test_avgtoedi_headstub_empty_none_not_set():
    stub = tr.AVGtoEDI._HeadStub("S2")
    assert not hasattr(stub, "empty")


def test_jtoedi_headstub_sets_attrs():
    stub = tr.JtoEDI._HeadStub("S3", lat=10.0, lon=20.0, elev=30.0, empty=5.0)
    assert stub.dataid == "S3"
    assert stub.lat == 10.0
    assert stub.long == 20.0
    assert stub.elev == 30.0
    assert stub.empty == 5.0


# --------------------------------------------------------------------- #
# AVGtoEDI.compute_z_from_res
# --------------------------------------------------------------------- #


def test_compute_z_from_res_returns_unchanged_when_missing_parts():
    b = TFBundle(freq=None, rho=[1.0], phase=[2.0])
    out = tr.AVGtoEDI().compute_z_from_res(b)
    assert out is b
    assert out.z is None


def test_compute_z_from_res_builds_complex_z_from_tensor_rho():
    # rho/phase are full (n, 2, 2) tensors, matching what
    # AVG.to_tensor(var="rho"/"phase") actually returns (and what the
    # rest of this module assumes, e.g. emit_edi's rho/phase branch).
    f = np.array([1.0, 10.0, 100.0])
    rho = np.ones((3, 2, 2)) * 100.0
    phase = np.ones((3, 2, 2)) * 300.0  # milliradians
    b = TFBundle(freq=f, rho=rho, phase=phase)
    out = tr.AVGtoEDI().compute_z_from_res(b)
    assert out.z is not None
    z = np.asarray(out.z)
    assert z.shape == (3, 2, 2)
    assert np.iscomplexobj(z)
    assert np.all(np.isfinite(z))


def test_compute_z_from_res_builds_complex_z_from_flat_rho():
    f = np.array([1.0, 10.0, 100.0])
    rho = np.array([100.0, 100.0, 100.0])
    phase = np.array([0.0, 500.0, 1000.0])  # milliradians
    b = TFBundle(freq=f, rho=rho, phase=phase)
    out = tr.AVGtoEDI().compute_z_from_res(b)
    assert out.z is not None
    z = np.asarray(out.z)
    assert z.shape == (3,)
    assert np.iscomplexobj(z)
    assert np.all(np.isfinite(z))


# --------------------------------------------------------------------- #
# FakeAVG helpers
# --------------------------------------------------------------------- #


class _Comp:
    def __init__(self, comps):
        self.frame = pd.DataFrame({"comp": comps})


class _CompNoColumn:
    def __init__(self):
        self.frame = pd.DataFrame({"other": [1, 2]})


class _Info:
    def __init__(self, comps):
        self.comp = _Comp(comps)


class _InfoNoCompColumn:
    def __init__(self):
        self.comp = _CompNoColumn()


class _Topo:
    def __init__(self, frame=None):
        self.frame = frame if frame is not None else pd.DataFrame()


class _AlwaysFailAVG(BaseAVG):
    """to_tensor always raises for every var -> no usable TFs."""

    def __init__(self):
        self.info = _Info([])
        self.topo = _Topo()

    def to_tensor(self, var="z", **kw):
        raise RuntimeError("no data")


class _RhoPhaseOnlyAVG(BaseAVG):
    """z / z_err unavailable; rho+phase present -> exercises the
    fallback branch inside _iter_bundles."""

    def __init__(self):
        self.info = _Info(["EX", "EY"])
        self.topo = _Topo()
        self._freq = np.array([1.0, 2.0, 4.0])
        self._rho = np.ones((3, 2, 2)) * 100.0
        self._phase = np.ones((3, 2, 2)) * 300.0  # mrad
        self._station = np.array([200.0])

    def to_tensor(self, var="z", **kw):
        if var == "z":
            raise RuntimeError("no Z in this dataset")
        if var == "z_err":
            raise RuntimeError("no Z_err either")
        if var == "rho":
            return self._rho, self._freq, self._station
        if var == "phase":
            return self._phase, self._freq, self._station
        raise KeyError(var)


class _MagOnlyAVG(BaseAVG):
    """Only magnetic-flavoured components (no EX/EY at all) with all
    three H channels present, to hit the HX/HZ construction branches
    and the want_ex/want_ey False branches in emit_edi."""

    def __init__(self):
        self.info = _Info(["HX", "HY", "HZ"])
        self.topo = _Topo()
        self._freq = np.array([1.0, 2.0, 4.0])
        z = np.zeros((3, 2, 2), complex)
        z[:, 0, 1] = 1.0 + 1.0j
        z[:, 1, 0] = 2.0 + 2.0j
        self._z = z
        self._ze = np.ones_like(z, float) * 0.1
        self._station = np.array([300.0])

    def to_tensor(self, var="z", **kw):
        if var == "z":
            return self._z, self._freq, self._station
        if var == "z_err":
            return self._ze, self._freq, self._station
        raise KeyError(var)


# --------------------------------------------------------------------- #
# extract() / _iter_bundles branches
# --------------------------------------------------------------------- #


def test_extract_raises_when_no_tf_available():
    with pytest.raises(ValueError, match="no transfer functions in AVG"):
        tr.AVGtoEDI().extract(_AlwaysFailAVG())


def test_transform_returns_empty_collection_when_no_tf_available():
    # transform() (unlike extract()) tolerates zero bundles: the loop
    # over an empty _iter_bundles() result just yields an empty
    # collection rather than raising.
    out = tr.AVGtoEDI().transform(_AlwaysFailAVG())
    assert len(out) == 0


def test_extract_returns_first_bundle():
    avg = _MagOnlyAVG()
    b = tr.AVGtoEDI().extract(avg)
    assert b.freq is not None
    assert b.z is not None


def test_iter_bundles_falls_back_to_rho_phase():
    avg = _RhoPhaseOnlyAVG()
    coll = tr.AVGtoEDI().transform(avg)
    assert len(coll) == 1
    ed = next(iter(coll))
    assert ed.Z._resistivity is not None


def test_unique_components_missing_comp_column_returns_empty_set():
    avg = SimpleNamespace(info=_InfoNoCompColumn())
    comps = tr.AVGtoEDI()._unique_components_from_avg(avg)
    assert comps == set()


def test_has_any_magnetic_false_for_electric_only():
    avg = SimpleNamespace(info=_Info(["EX", "EY"]))
    assert tr.AVGtoEDI()._has_any_magnetic(avg) is False


# --------------------------------------------------------------------- #
# emit_edi: magnetic-only channel set + tipper + rho/phase elif branch
# --------------------------------------------------------------------- #


def test_emit_edi_builds_hx_and_hz_channels_when_no_electric():
    avg = _MagOnlyAVG()
    coll = tr.AVGtoEDI().transform(avg)
    ed = next(iter(coll))
    dm = ed.get_section("definemeas")
    ids = {m.chtype for m in dm.hmeas}
    assert "HX" in ids
    assert "HZ" in ids
    # no electric channels were requested
    assert dm.emeas == []


def test_emit_edi_tipper_branch_is_populated():
    f = np.array([1.0, 2.0, 3.0])
    z = np.ones((3, 2, 2), complex)
    tip = np.ones((3, 1, 2), complex) * 0.5
    tip_err = np.ones((3, 1, 2), float) * 0.05
    b = TFBundle(
        freq=f,
        z=z,
        rho=None,
        phase=None,
        tipper=tip,
        tipper_err=tip_err,
        station="TPX",
    )
    ed = tr.AVGtoEDI().emit_edi(b)
    assert ed.Tip._tipper is not None
    assert ed.Tip._tipper_err is not None


def test_emit_edi_rho_phase_only_branch():
    f = np.array([1.0, 2.0, 4.0])
    rho = np.ones((3, 2, 2)) * 50.0
    phase = np.ones((3, 2, 2)) * 400.0  # mrad
    b = TFBundle(freq=f, z=None, rho=rho, phase=phase, station="RPX")
    ed = tr.AVGtoEDI().emit_edi(b)
    assert ed.Z._resistivity is not None


def test_emit_edi_raises_when_neither_z_nor_rho_phase():
    b = TFBundle(freq=np.array([1.0]), z=None, rho=None, phase=None)
    with pytest.raises(ValueError, match="Neither Z nor"):
        tr.AVGtoEDI().emit_edi(b)


# --------------------------------------------------------------------- #
# _apply_topo_coords branches
# --------------------------------------------------------------------- #


def _fresh_edi(station="X1"):
    ed = EDIFile(verbose=0)
    ed.station = station
    return ed


def test_apply_topo_coords_no_matching_column_returns():
    ed = _fresh_edi()
    frame = pd.DataFrame({"latitude": [1.0], "longitude": [2.0]})
    tr.AVGtoEDI()._apply_topo_coords(ed, frame, "X1")
    # nothing crashed and no head was force-created by this path
    assert ed.get_section("head") is None


def test_apply_topo_coords_no_matching_row_returns():
    ed = _fresh_edi()
    frame = pd.DataFrame({"station": ["A", "B"], "latitude": [1.0, 2.0]})
    tr.AVGtoEDI()._apply_topo_coords(ed, frame, "ZZZ")
    assert ed.get_section("head") is None


def test_apply_topo_coords_creates_head_when_absent():
    ed = _fresh_edi()
    frame = pd.DataFrame(
        {
            "station": ["X1"],
            "latitude": [12.5],
            "longitude": [34.5],
            "elevation": [100.0],
        }
    )
    tr.AVGtoEDI()._apply_topo_coords(ed, frame, "X1")
    h = ed.get_section("head")
    assert h is not None
    assert h.lat == pytest.approx(12.5)
    assert h.long == pytest.approx(34.5)


def test_apply_topo_coords_bad_numeric_value_falls_back_to_default():
    ed = _fresh_edi()
    frame = pd.DataFrame(
        {
            "station": ["X1"],
            "latitude": ["not-a-number"],
            "longitude": [34.5],
            "elevation": ["also-bad"],
        }
    )
    tr.AVGtoEDI()._apply_topo_coords(ed, frame, "X1")
    h = ed.get_section("head")
    assert h is not None
    assert h.lat == pytest.approx(0.0)
    assert h.elev == pytest.approx(0.0)


def test_apply_topo_coords_updates_existing_definemeas():
    from pycsamt.seg.meas import DefineMeas

    ed = _fresh_edi()
    ed.add_section("definemeas", DefineMeas())
    frame = pd.DataFrame(
        {
            "station": ["X1"],
            "latitude": [1.0],
            "longitude": [2.0],
            "elevation": [3.0],
        }
    )
    tr.AVGtoEDI()._apply_topo_coords(ed, frame, "X1")
    dm = ed.get_section("definemeas")
    assert dm.reflat == pytest.approx(1.0)
    assert dm.reflong == pytest.approx(2.0)
    assert dm.refelev == pytest.approx(3.0)


def test_apply_topo_coords_string_station_id_match():
    ed = _fresh_edi()
    frame = pd.DataFrame(
        {"station": ["7"], "latitude": [1.0], "longitude": [2.0]}
    )
    # station_id passed as int, only matches after str() coercion
    tr.AVGtoEDI()._apply_topo_coords(ed, frame, 7)
    h = ed.get_section("head")
    assert h is not None
    assert h.lat == pytest.approx(1.0)


# --------------------------------------------------------------------- #
# _enrich_from_avg_info
# --------------------------------------------------------------------- #


def test_enrich_from_avg_info_populates_head_and_info():
    ed = _fresh_edi()
    meta = {
        "stdvers": "SEG 1.0",
        "progvers": "9.9",
        "progdate": "2024-01-01",
        "acqdate": "2024-01-02",
        "filedate": "2024-01-03",
        "acqby": "Alice",
        "fileby": "Bob",
        "prospect": "Prospect1",
        "loc": "Location1",
        "maxsect": "5",
        "empty": "1.0e32",
        "survey_co": "ACME",
        "client_co": "Client1",
        "area": "Area1",
    }
    source = SimpleNamespace(info=meta)
    tr.AVGtoEDI()._enrich_from_avg_info(ed, source)

    h = ed.get_section("head")
    assert h is not None
    assert h.stdvers == "SEG 1.0"
    assert h.progvers == "9.9"
    assert h.acqby == "Alice"
    assert h.fileby == "Bob"
    assert h.prospect == "Prospect1"
    assert h.loc == "Location1"
    assert h.maxsect == 5
    assert h.empty == pytest.approx(1.0e32)

    info = ed.get_section("info")
    assert info is not None
    joined = "\n".join(info.info_text)
    assert "SURVEY CO:ACME" in joined
    assert "CLIENT CO:Client1" in joined
    assert "AREA:Area1" in joined


def test_enrich_from_avg_info_bad_maxsect_and_empty_are_swallowed():
    ed = _fresh_edi()
    meta = {"maxsect": "not-an-int", "empty": "not-a-float"}
    source = SimpleNamespace(info=meta)
    # must not raise
    tr.AVGtoEDI()._enrich_from_avg_info(ed, source)
    h = ed.get_section("head")
    assert h is not None


def test_enrich_from_avg_info_no_meta_is_noop():
    ed = _fresh_edi()
    source = SimpleNamespace()
    tr.AVGtoEDI()._enrich_from_avg_info(ed, source)
    assert ed.get_section("head") is None


def test_enrich_from_avg_info_keeps_existing_survey_id_and_rotation():
    ed = _fresh_edi()
    info = ed.get_section("info")
    if info is None:
        info = tr.AVGtoEDI()._ensure_info(ed, survey_id="X1")
    info.update(info_text=["  SURVEY ID:X1", "  ROTATION=FIX"])
    source = SimpleNamespace(info={"survey_co": "ACME"})
    tr.AVGtoEDI()._enrich_from_avg_info(ed, source)
    joined = "\n".join(info.info_text)
    assert joined.count("SURVEY ID:") == 1
    assert joined.count("ROTATION=") == 1


# --------------------------------------------------------------------- #
# post_emit exception-swallowing branches
# --------------------------------------------------------------------- #


def test_post_emit_swallows_ensure_station_failure(monkeypatch):
    def _boom(*a, **k):
        raise RuntimeError("boom")

    monkeypatch.setattr(tr, "ensure_station", _boom)
    ed = _fresh_edi()
    b = TFBundle(freq=np.array([1.0]), station="S1", station_id=1)
    out = tr.AVGtoEDI().post_emit(ed, SimpleNamespace(), b)
    assert out is ed


def test_post_emit_swallows_topo_coord_failure(monkeypatch):
    def _boom(self, *a, **k):
        raise RuntimeError("boom")

    monkeypatch.setattr(tr.AVGtoEDI, "_apply_topo_coords", _boom)
    ed = _fresh_edi()
    b = TFBundle(freq=np.array([1.0]), station="S1", station_id=1)
    out = tr.AVGtoEDI().post_emit(
        ed, SimpleNamespace(topo=SimpleNamespace(frame=pd.DataFrame())), b
    )
    assert out is ed


def test_post_emit_swallows_enrich_failure(monkeypatch):
    def _boom(self, *a, **k):
        raise RuntimeError("boom")

    monkeypatch.setattr(tr.AVGtoEDI, "_enrich_from_avg_info", _boom)
    ed = _fresh_edi()
    b = TFBundle(freq=np.array([1.0]), station="S1", station_id=1)
    out = tr.AVGtoEDI().post_emit(ed, SimpleNamespace(), b)
    assert out is ed


# ======================================================================
# JtoEDI side
# ======================================================================


class _RP:
    def __init__(self, rho, phase):
        self.rho = rho
        self.phase = phase


def test_bundle_from_j_uses_resphase_object():
    jf = SimpleNamespace(
        Z=None,
        ResPhase=_RP(np.array([10.0, 20.0]), np.array([1.0, 2.0])),
        Tipper=None,
        freq=np.array([1.0, 2.0]),
        station="J1",
    )
    b = tr.JtoEDI()._bundle_from_j(jf)
    assert list(b.rho) == [10.0, 20.0]
    assert list(b.phase) == [1.0, 2.0]


def test_bundle_from_j_freq_falls_back_to_jfile_level():
    jf = SimpleNamespace(
        Z=None,
        ResPhase=None,
        Tipper=None,
        freq=np.array([5.0, 6.0]),
        rho=np.array([1.0, 2.0]),
        phase=np.array([3.0, 4.0]),
        station="J2",
    )
    b = tr.JtoEDI()._bundle_from_j(jf)
    assert list(b.freq) == [5.0, 6.0]


class _StubJFile:
    """Duck-typed stand-in for isinstance checks in ``_as_jfile``.

    ``_as_jfile`` requires ``isinstance(src, JFile)``; the real class
    needs a file on disk to construct meaningfully, so tests that only
    care about ``_bundle_from_j``'s duck-typed attribute access
    monkeypatch the *name* jedi.py binds (``tr.JFile``) rather than
    constructing a real one (see fix_monkeypatch_module_level_import
    convention: patch the consuming module's bound name).
    """


def test_jtoedi_extract_returns_bundle(monkeypatch):
    monkeypatch.setattr(tr, "JFile", _StubJFile)
    jf = _StubJFile()
    jf.Z = SimpleNamespace(z=np.ones((2, 2, 2), complex), freq=np.array([1.0, 2.0]))
    jf.ResPhase = None
    jf.Tipper = None
    jf.station = "J3"
    b = tr.JtoEDI().extract(jf)
    assert b.freq is not None


def test_jtoedi_extract_raises_when_empty(monkeypatch):
    monkeypatch.setattr(tr, "JFile", _StubJFile)
    jf = _StubJFile()
    jf.Z = None
    jf.ResPhase = None
    jf.Tipper = None
    with pytest.raises(ValueError, match="no transfer functions in J file"):
        tr.JtoEDI().extract(jf)


def test_jtoedi_emit_edi_rho_phase_branch():
    f = np.array([1.0, 2.0, 4.0])
    rho = np.ones((3, 2, 2)) * 80.0
    phase = np.ones((3, 2, 2)) * 30.0  # degrees
    b = TFBundle(freq=f, z=None, rho=rho, phase=phase, station="JRP")
    ed = tr.JtoEDI().emit_edi(b)
    assert ed.Z._freq is not None


def test_jtoedi_emit_edi_tipper_branch():
    f = np.array([1.0, 2.0, 3.0])
    z = np.ones((3, 2, 2), complex)
    tip = np.ones((3, 1, 2), complex) * 0.3
    tip_err = np.ones((3, 1, 2), float) * 0.03
    b = TFBundle(
        freq=f, z=z, tipper=tip, tipper_err=tip_err, station="JTIP"
    )
    ed = tr.JtoEDI().emit_edi(b)
    assert ed.Tip._tipper is not None
    assert ed.Tip._tipper_err is not None


def test_jtoedi_apply_j_coords_noop_when_all_none():
    ed = _fresh_edi()
    jf = SimpleNamespace()
    tr.JtoEDI()._apply_j_coords(ed, jf)
    assert ed.get_section("head") is None


def test_jtoedi_apply_j_coords_sets_head_and_definemeas():
    from pycsamt.seg.meas import DefineMeas

    ed = _fresh_edi()
    ed.add_section("definemeas", DefineMeas())
    jf = SimpleNamespace(lat=11.0, lon=22.0, elev=33.0)
    tr.JtoEDI()._apply_j_coords(ed, jf)
    h = ed.get_section("head")
    assert h.lat == pytest.approx(11.0)
    dm = ed.get_section("definemeas")
    assert dm.reflat == pytest.approx(11.0)


def test_jtoedi_seed_meas_and_mtsect_have_xy_only():
    ed = _fresh_edi()
    f = np.array([1.0, 2.0])
    ed.Z._freq = f
    z = np.zeros((2, 2, 2), complex)
    z[:, 0, 1] = 1.0  # Zxy nonzero -> have_xy
    ed.Z._z = z
    jf = SimpleNamespace(azimuth=15.0)
    tr.JtoEDI()._seed_meas_and_mtsect(ed, jf)
    dm = ed.get_section("definemeas")
    hids = {m.chtype for m in dm.hmeas}
    eids = {m.chtype for m in dm.emeas}
    assert "HY" in hids
    assert "EX" in eids
    assert "HX" not in hids


def test_jtoedi_seed_meas_and_mtsect_both_and_tipper_skips_duplicates():
    ed = _fresh_edi()
    f = np.array([1.0, 2.0])
    ed.Z._freq = f
    z = np.zeros((2, 2, 2), complex)
    z[:, 0, 1] = 1.0
    z[:, 1, 0] = 1.0
    ed.Z._z = z
    ed.Tip._tipper = np.ones((2, 1, 2), complex)
    jf = SimpleNamespace()
    tr.JtoEDI()._seed_meas_and_mtsect(ed, jf)
    dm = ed.get_section("definemeas")
    n_h_first = len(dm.hmeas)
    n_e_first = len(dm.emeas)
    # calling again must not duplicate ids
    tr.JtoEDI()._seed_meas_and_mtsect(ed, jf)
    assert len(dm.hmeas) == n_h_first
    assert len(dm.emeas) == n_e_first
    hids = {m.chtype for m in dm.hmeas}
    assert {"HX", "HY", "HZ"} <= hids


def test_jtoedi_post_emit_info_enrichment_from_heads(monkeypatch):
    monkeypatch.setattr(tr, "ensure_station", lambda *a, **k: "JH")
    ed = _fresh_edi()
    jf = SimpleNamespace(
        heads=SimpleNamespace(info={"anything": 1}),
    )
    b = TFBundle(freq=np.array([1.0]), station="JH", station_id=1)
    out = tr.JtoEDI().post_emit(ed, jf, b)
    assert out is ed
    info = out.get_section("info")
    assert info is not None
    # real bug fix: Info has no .setdefault; enrichment must not crash
    # and processedby/signconvention remain populated (already defaulted
    # at construction time).
    assert info.Processing.processedby
    assert info.Processing.signconvention


def test_jtoedi_post_emit_swallows_seed_meas_failure(monkeypatch):
    def _boom(self, *a, **k):
        raise RuntimeError("boom")

    monkeypatch.setattr(tr.JtoEDI, "_seed_meas_and_mtsect", _boom)
    ed = _fresh_edi()
    jf = SimpleNamespace()
    b = TFBundle(freq=np.array([1.0]), station="J1", station_id=1)
    out = tr.JtoEDI().post_emit(ed, jf, b)
    assert out is ed


def test_jtoedi_transform_collection_finalizes_each_member(monkeypatch):
    f = np.array([1.0, 2.0])
    z = np.ones((2, 2, 2), complex)

    class _JF:
        def __init__(self, station):
            self.Z = SimpleNamespace(z=z, freq=f)
            self.ResPhase = None
            self.Tipper = None
            self.station = station

    class _StubJCollection(list):
        pass

    monkeypatch.setattr(tr, "JCollection", _StubJCollection)
    coll = _StubJCollection([_JF("JC1"), _JF("JC2")])
    out = tr.JtoEDI().transform(coll)
    assert len(out) == 2
