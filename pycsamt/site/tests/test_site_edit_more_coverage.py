"""Targeted branch coverage for pycsamt.site.edit.

These tests exercise the many "no-throw by design" skip branches in
rotate/select_freq/rename/fill_missing/set_coords_all plus the private
table/coordinate helpers, using lightweight duck-typed stand-ins for
EDI-like objects (the public functions only ever use getattr/setattr).
"""

from __future__ import annotations

import builtins
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from pycsamt.site import edit as ed
from pycsamt.site.base import Sites


class _Obj:
    """Plain attribute bag, deepcopy-friendly."""


def _mk(**kw):
    o = _Obj()
    for k, v in kw.items():
        setattr(o, k, v)
    return o


def _load_edi(p: Path):
    from pycsamt.seg.edi import EDIFile

    return EDIFile(p)


def _dup_edi(tmp_path: Path, src: Path, stem: str) -> Path:
    dst = tmp_path / f"{stem}.edi"
    dst.write_text(src.read_text(encoding="utf-8"), encoding="utf-8")
    return dst


def _mk_two_edifiles(tmp_path, simulated_edi, s1, s2):
    p1 = _dup_edi(tmp_path, simulated_edi, s1)
    p2 = _dup_edi(tmp_path, simulated_edi, s2)
    return _load_edi(p1), _load_edi(p2)


# --------------------------------------------------------------- rotate ---


def test_rotate_noop_when_no_z_and_no_tipper():
    e = _mk()
    out = ed.rotate(e, 30.0, inplace=True)
    assert out is e


def test_rotate_z_present_but_no_z_array():
    Z = _mk(z_error=np.zeros((2, 2, 2)))
    e = _mk(Z=Z)
    out = ed.rotate(e, 10.0, inplace=True)
    assert out is e
    assert np.all(out.Z.z_error == 0)  # untouched: zz was None


def test_rotate_z_wrong_shape_still_processes_error_with_wrong_shape():
    Z = _mk(z=np.zeros((3, 4)), z_error=np.zeros((3, 4)))
    e = _mk(Z=Z)
    out = ed.rotate(e, 15.0, inplace=True)
    # neither z nor z_error match (n,2,2); both left unrotated
    assert out.Z.z.shape == (3, 4)
    assert out.Z.z_error.shape == (3, 4)


def test_rotate_z_without_error_array():
    Z = _mk(z=np.ones((2, 2, 2), complex))
    e = _mk(Z=Z)
    out = ed.rotate(e, 12.0, inplace=True)
    assert out.Z.z.shape == (2, 2, 2)


def test_rotate_tip_obj_on_t_attribute():
    tip = _mk(tipper=np.array([[1 + 1j, 2 + 2j], [3 + 3j, 4 + 4j]]))
    e = _mk(T=tip)
    out = ed.rotate(e, 25.0, inplace=True)
    assert out.T.tipper.shape == (2, 2)


def test_rotate_tip_obj_missing_tipper_attribute():
    tip = _mk()
    e = _mk(T=tip)
    out = ed.rotate(e, 25.0, inplace=True)
    assert out.T is tip


def test_rotate_z_tipper_else_branch():
    Z = _mk(
        z=np.zeros((1, 2, 2), complex),
        tipper=np.array([[1 + 1j, 2 + 2j]]),
    )
    e = _mk(Z=Z)  # no T/TIP/Tip attribute at all
    out = ed.rotate(e, 5.0, inplace=True)
    assert out.Z.tipper.shape == (1, 2)


def test_rotate_z_tipper_wrong_shape_skipped():
    Z = _mk(z=np.zeros((1, 2, 2), complex), tipper=np.array([1 + 1j, 2 + 2j]))
    e = _mk(Z=Z)
    out = ed.rotate(e, 5.0, inplace=True)
    assert out.Z.tipper.shape == (2,)


# ---------------------------------------------------------- select_freq ---


def test_select_freq_keep_with_no_z_and_no_tipper():
    e = _mk()
    out = ed.select_freq(e, keep=[0], inplace=True)
    assert out is e


def test_select_freq_no_z_returns_unchanged():
    e = _mk()
    out = ed.select_freq(e, fmin=1.0, inplace=True)
    assert out is e


def test_select_freq_z_without_freq_attribute():
    Z = _mk(z=np.zeros((2, 2, 2), complex))
    e = _mk(Z=Z)
    out = ed.select_freq(e, fmin=1.0, inplace=True)
    assert out is e  # no freq -> f empty -> no-op


def test_select_freq_fmax_only(simulated_edi):
    edf = _load_edi(simulated_edi)
    f0 = ed.get_freq(edf)
    if f0.size < 2:
        pytest.skip("simulated EDI has <2 freq")
    fmax = float(np.median(f0))
    out = ed.select_freq(edf, fmax=fmax, inplace=False)
    f1 = ed.get_freq(out)
    assert np.all(f1 <= fmax)


# --------------------------------------------------------------- rename ---


def test_rename_noop_without_name_or_policy(simulated_edi):
    edf = _load_edi(simulated_edi)
    out = ed.rename(edf, inplace=True)
    assert out is edf


def test_rename_policy_raises_falls_back_to_current_name(simulated_edi):
    edf = _load_edi(simulated_edi)
    cur = ed.station_name(edf)

    def bad_policy(_n):
        raise RuntimeError("boom")

    out = ed.rename(edf, policy=bad_policy, inplace=False)
    h = out.get_section("head")
    assert str(h.dataid) == cur


def test_rename_header_setattr_errors_suppressed():
    class StrictHead:
        def __setattr__(self, key, value):
            if key in ("sitename", "name", "STATION"):
                raise RuntimeError("blocked")
            object.__setattr__(self, key, value)

    class StrictEd:
        def __init__(self):
            self.Z = None
            self._head = StrictHead()

        def get_section(self, name):
            if name == "head":
                return self._head
            raise KeyError(name)

        def set_section(self, name, value):
            pass

        def __setattr__(self, key, value):
            if key == "name":
                raise RuntimeError("blocked name")
            object.__setattr__(self, key, value)

    e = StrictEd()
    out = ed.rename(e, name="NEWX", inplace=True)
    assert out is e
    assert e._head.dataid == "NEWX"


def test_rename_ensure_head_dataid_error_suppressed():
    class BadDataidHead:
        def __setattr__(self, key, value):
            if key == "dataid":
                raise RuntimeError("no dataid")
            object.__setattr__(self, key, value)

    class BadDataidEd:
        def __init__(self):
            self.Z = None
            self._head = BadDataidHead()

        def get_section(self, name):
            if name == "head":
                return self._head
            raise KeyError(name)

        def set_section(self, name, value):
            pass

    e = BadDataidEd()
    out = ed.rename(e, name="Y1", inplace=True)
    assert out is e
    assert out.name == "Y1"


# ---------------------------------------------------------- fill_missing --


def test_fill_missing_full_z_and_tip_all_optional_fields():
    Z = _mk(
        z=np.array([[[1.0, np.nan], [2.0, 3.0]]]),
        z_error=np.array([[[0.1, np.inf], [0.1, 0.1]]]),
        rho=np.array([[[1.0, np.nan], [1.0, 1.0]]]),
        phase=np.array([[[1.0, np.nan], [1.0, 1.0]]]),
        phase_err=np.array([[[0.1, np.nan], [0.1, 0.1]]]),
    )
    tp = _mk(
        tipper=np.array([[1.0, np.nan]]),
        tipper_err=np.array([[0.1, np.nan]]),
    )
    e = _mk(Z=Z, T=tp)
    out = ed.fill_missing(e, how="nan", components=("Z", "Tip"), inplace=True)
    assert np.isnan(out.Z.z[0, 0, 1])
    assert np.isnan(out.Z.z_error[0, 0, 1])
    assert np.isnan(out.Z.rho[0, 0, 1])
    assert np.isnan(out.Z.phase[0, 0, 1])
    assert np.isnan(out.Z.phase_err[0, 0, 1])
    assert np.isnan(out.T.tipper[0, 1])
    assert np.isnan(out.T.tipper_err[0, 1])


def test_fill_missing_minimal_z_skips_missing_optional_fields():
    Z = _mk(z=np.array([[[1.0, np.nan], [2.0, 3.0]]]))
    e = _mk(Z=Z)  # no z_error/rho/phase/phase_err, no T/TIP/Tip
    out = ed.fill_missing(e, how="zero", components=("Z", "Tip"), inplace=True)
    assert out.Z.z[0, 0, 1] == 0.0


def test_fill_missing_no_z_section_but_tip_present():
    tp = _mk(tipper=np.array([[1.0, np.nan]]))
    e = _mk(Tip=tp)  # no Z attribute at all
    out = ed.fill_missing(e, how="zero", components=("Z", "Tip"), inplace=True)
    assert out.Tip.tipper[0, 1] == 0.0


def test_fill_missing_invalid_how_raises():
    with pytest.raises(ValueError):
        ed.fill_missing(_mk(), how="bogus")


# ------------------------------------------------------- recompute_res_phase


def test_recompute_res_phase_no_z():
    e = _mk()
    out = ed.recompute_res_phase(e, inplace=True)
    assert out is e


def test_recompute_res_phase_fn_not_callable():
    Z = _mk(compute_resistivity_phase=None)
    e = _mk(Z=Z)
    out = ed.recompute_res_phase(e, inplace=True)
    assert out is e


def test_recompute_res_phase_typeerror_then_success():
    calls = []

    def fn(*args):
        if not args:
            raise TypeError("needs args")
        calls.append(args)

    Z = _mk(compute_resistivity_phase=fn)
    e = _mk(Z=Z)
    ed.recompute_res_phase(e, inplace=True)
    assert calls == [(None, None, None)]


def test_recompute_res_phase_typeerror_then_failure_suppressed():
    def fn(*args):
        if not args:
            raise TypeError("no args")
        raise RuntimeError("still bad")

    Z = _mk(compute_resistivity_phase=fn)
    e = _mk(Z=Z)
    out = ed.recompute_res_phase(e, inplace=True)
    assert out is e


def test_recompute_res_phase_unexpected_error_suppressed():
    class RaisingZ:
        def __getattr__(self, item):
            raise RuntimeError("nope")

    e = _mk(Z=RaisingZ())
    out = ed.recompute_res_phase(e, inplace=True)
    assert out is e


# --------------------------------------------------------- set_coords_all -


def test_set_coords_all_callable_src(tmp_path, simulated_edi):
    e1 = _load_edi(_dup_edi(tmp_path, simulated_edi, "H01"))

    def picker(_edi):
        return (1.0, 2.0, 3.0)

    out = ed.set_coords_all([e1], picker, inplace=False)
    assert tuple(map(float, out.by_index(0).coords)) == (1.0, 2.0, 3.0)


def test_set_coords_all_callable_src_raises_suppressed(tmp_path, simulated_edi):
    e1 = _load_edi(_dup_edi(tmp_path, simulated_edi, "H02"))

    def bad(_edi):
        raise RuntimeError("boom")

    out = ed.set_coords_all([e1], bad, inplace=False)
    assert out is not None


def test_set_coords_all_callable_src_returns_none(tmp_path, simulated_edi):
    e1 = _load_edi(_dup_edi(tmp_path, simulated_edi, "H03"))

    def none_picker(_edi):
        return None

    out = ed.set_coords_all([e1], none_picker, inplace=False)
    assert out is not None


def test_set_coords_all_mapping_get_raises_suppressed(tmp_path, simulated_edi):
    e1 = _load_edi(_dup_edi(tmp_path, simulated_edi, "H04"))

    class BadMapping:
        def get(self, _name):
            raise KeyError("nope")

    out = ed.set_coords_all([e1], BadMapping(), inplace=False)
    assert out is not None


class _FrameHolder:
    def __init__(self, frame=None):
        if frame is not None:
            self.frame = frame


def test_set_coords_all_frame_lookup_variants(tmp_path, simulated_edi):
    e1 = _load_edi(_dup_edi(tmp_path, simulated_edi, "G01"))
    e2 = _load_edi(_dup_edi(tmp_path, simulated_edi, "G02"))
    n1 = ed.station_name(e1)
    n2 = ed.station_name(e2)

    # 1) no frame attribute at all
    out0 = ed.set_coords_all([e1], _FrameHolder(), inplace=False)
    assert out0 is not None

    # 2) frame missing 'station' column
    df_no_station = pd.DataFrame({"name": [n1], "lat": [1.0], "lon": [2.0]})
    out1 = ed.set_coords_all([e1], _FrameHolder(df_no_station), inplace=False)
    assert out1 is not None

    # 3) frame with station col but no matching row
    df_no_match = pd.DataFrame({"station": ["ZZZ"], "lat": [1.0], "lon": [2.0]})
    out2 = ed.set_coords_all([e1], _FrameHolder(df_no_match), inplace=False)
    assert out2 is not None

    # 4) successful match via latitude/longitude aliases, no elev column
    df_alias = pd.DataFrame(
        {"station": [n1], "latitude": [10.0], "longitude": [20.0]}
    )
    out3 = ed.set_coords_all([e1], _FrameHolder(df_alias), inplace=False)
    s3 = out3.by_index(0)
    assert tuple(map(float, s3.coords)) == (10.0, 20.0, 0.0)

    # 5) successful match with lat/lon/elev present
    df_full = pd.DataFrame(
        {"station": [n2], "lat": [5.0], "lon": [6.0], "elev": [7.0]}
    )
    out4 = ed.set_coords_all([e2], _FrameHolder(df_full), inplace=False)
    s4 = out4.by_index(0)
    assert tuple(map(float, s4.coords)) == (5.0, 6.0, 7.0)

    # 6) exception path: non-numeric lat triggers float() failure
    df_bad = pd.DataFrame({"station": [n1], "lat": ["oops"], "lon": [1.0]})
    out5 = ed.set_coords_all([e1], _FrameHolder(df_bad), inplace=False)
    assert out5 is not None


# ---------------------------------------------------- _slice_fields (priv) -


def test_slice_fields_suppresses_setattr_exception():
    class Locked:
        def __init__(self):
            self._locked = False
            self.freq = np.array([1.0, 2.0])

        def __setattr__(self, key, value):
            if key == "freq" and getattr(self, "_locked", False):
                raise RuntimeError("locked")
            object.__setattr__(self, key, value)

    obj = Locked()
    obj._locked = True
    ed._slice_fields(obj, np.array([0]))  # must not raise


def test_slice_fields_recompute_exception_suppressed():
    class Obj2:
        def __init__(self):
            self.freq = np.array([1.0, 2.0])

        def compute_resistivity_phase(self):
            raise RuntimeError("boom")

    obj = Obj2()
    ed._slice_fields(obj, np.array([0]))  # must not raise


# ------------------------------------------------------------- _maybe_df --


def test_maybe_df_missing_elev_column_defaults_none():
    table = [{"station": "S1", "lat": 1.0, "lon": 2.0}]
    _, cols = ed._maybe_df(table)
    assert cols["elev"] is None


def test_maybe_df_csv_read_exception_fallback(tmp_path, monkeypatch):
    p = tmp_path / "coords.csv"
    p.write_text("station,lat,lon,elev\nS1,1.0,2.0,3.0\n", encoding="utf-8")

    orig_read_csv = ed.pd.read_csv
    calls = {"n": 0}

    def flaky_read_csv(path, **kwargs):
        calls["n"] += 1
        if calls["n"] == 1:
            raise ValueError("forced failure")
        return orig_read_csv(path)

    monkeypatch.setattr(ed.pd, "read_csv", flaky_read_csv)
    df, cols = ed._maybe_df(p)
    assert calls["n"] == 2
    assert cols["station"] == "station"


def test_maybe_df_unsupported_table_raises_typeerror():
    with pytest.raises(TypeError):
        ed._maybe_df(object())


def test_maybe_df_explicit_columns_mapping():
    table = [{"nm": "S1", "latx": 1.0, "lonx": 2.0}]
    _, cols = ed._maybe_df(
        table, columns={"station": "nm", "lat": "latx", "lon": "lonx"}
    )
    assert cols["station"] == "nm"
    assert cols["lat"] == "latx"
    assert cols["lon"] == "lonx"


def test_maybe_df_easting_northing_autodetect():
    table = [{"station": "S1", "easting": 400000.0, "northing": 5750000.0}]
    _, cols = ed._maybe_df(table)
    assert cols["easting"] == "easting"
    assert cols["northing"] == "northing"
    assert cols["lat"] is None
    assert cols["lon"] is None


def test_maybe_df_no_station_column_raises():
    table = [{"lat": 1.0, "lon": 2.0}]
    with pytest.raises(ValueError):
        ed._maybe_df(table)


# ------------------------------------------------ _project_en_to_lonlat --


def test_project_en_to_lonlat_importerror_when_pyproj_missing(monkeypatch):
    real_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name == "pyproj":
            raise ImportError("no pyproj")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    with pytest.raises(ImportError):
        ed._project_en_to_lonlat(
            np.array([1.0]), np.array([2.0]), "EPSG:32631"
        )


# ---------------------------------------------------- _frame_to_mapping --


def test_frame_to_mapping_raises_when_no_coords():
    df = pd.DataFrame({"station": ["S1"]})
    cols = {
        "station": "station",
        "lat": None,
        "lon": None,
        "elev": None,
        "easting": None,
        "northing": None,
    }
    with pytest.raises(ValueError):
        ed._frame_to_mapping(df, cols)


def test_frame_to_mapping_without_elev_defaults_nan():
    df = pd.DataFrame({"station": ["S1"], "lat": [1.0], "lon": [2.0]})
    cols = {
        "station": "station",
        "lat": "lat",
        "lon": "lon",
        "elev": None,
        "easting": None,
        "northing": None,
    }
    mp = ed._frame_to_mapping(df, cols)
    assert np.isnan(mp["S1"][2])


# ----------------------------------------------------- _set_attr_first --


def test_set_attr_first_all_targets_fail_suppressed():
    class Locked:
        def __setattr__(self, key, value):
            raise RuntimeError("locked")

    obj = Locked()
    ed._set_attr_first(obj, ("a", "b", "c"), 123)  # must not raise


# ---------------------------------------------------------- _each_site --


def test_each_site_accepts_sites_instance(tmp_path, simulated_edi):
    e1, e2 = _mk_two_edifiles(tmp_path, simulated_edi, "SI1", "SI2")
    sites_obj = Sites([e1, e2])
    out = ed.rotate_all(sites_obj, 10.0, inplace=False)
    assert isinstance(out, Sites)
    assert len(out) == 2


def test_wrap_output_inplace_with_sites_instance(tmp_path, simulated_edi):
    e1, e2 = _mk_two_edifiles(tmp_path, simulated_edi, "SJ1", "SJ2")
    sites_obj = Sites([e1, e2])
    out = ed.rotate_all(sites_obj, 5.0, inplace=True)
    assert out is sites_obj


# ------------------------------------------------------------ _fill_array


def test_fill_array_none_input_zero_and_nan():
    z = ed._fill_array(None, (2, 2), "zero")
    assert np.all(z == 0)
    n = ed._fill_array(None, (2, 2), "nan")
    assert np.all(np.isnan(n))


def test_fill_array_existing_array_nan_mode():
    arr = np.array([1.0, np.inf, np.nan, 2.0])
    out = ed._fill_array(arr, (4,), "nan")
    assert np.isnan(out[1])
    assert np.isnan(out[2])
    assert out[0] == 1.0
