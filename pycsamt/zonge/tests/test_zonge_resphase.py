# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

import numpy as np
import pandas as pd
import pytest

from pycsamt.exceptions import AvgDataError
from pycsamt.zonge.resphase import Phase, Resistivity


def _has_line(lines, pred):
    return any(pred(ln) for ln in lines)


def test_resistivity_read_legacy_and_to_xarray():
    # Legacy-style column name ('Resistivity'); no comp in input.
    df = pd.DataFrame(
        {
            "station": [100.0, 150.0],
            "freq": [1.0, 2.0],
            "Resistivity": [10.0, 20.0],
        }
    )
    r = Resistivity.from_avg((df, {"Unit.Rho": "ohm·m"}))

    # Canonical column present
    assert "rho" in r.frame.columns
    assert np.allclose(r.frame["rho"].to_numpy(), [10.0, 20.0])

    # comp should be injected as 'ExHy' when missing
    assert "comp" in r.frame.columns
    assert set(r.frame["comp"].unique()) == {"ExHy"}

    # xarray dataset has dims/vars/attrs as expected
    pytest.importorskip("xarray")
    ds = r.to_xarray()
    assert set(ds.dims) == {"station", "freq", "comp"}
    assert "rho" in ds.variables
    assert str(ds.attrs.get("Unit.Rho", "")).lower() in {"ohm·m", "ohm*m"}


def test_resistivity_to_tensor_union_and_intersection():
    # Two stations, ragged freqs; only ExHy present (others → NaN)
    rows = [
        # station 100: 1, 2
        dict(station=100.0, freq=1.0, comp="ExHy", ARes_mag=5.0),
        dict(station=100.0, freq=2.0, comp="ExHy", ARes_mag=6.0),
        # station 150: 2, 4
        dict(station=150.0, freq=2.0, comp="ExHy", ARes_mag=7.0),
        dict(station=150.0, freq=4.0, comp="ExHy", ARes_mag=8.0),
    ]
    df = pd.DataFrame.from_records(rows).rename(columns={"ARes_mag": "ARes.mag"})
    r = Resistivity.from_avg((df, {}))

    # UNION → freqs {1,2,4}
    T_u, f_u, st_u = r.to_tensor(var="rho", align="union")
    assert T_u.shape == (2, 3, 2, 2)
    assert np.allclose(np.sort(st_u), [100.0, 150.0])
    assert np.allclose(f_u, [1.0, 2.0, 4.0])

    # value lands at ExHy → (0,1) for station=100, freq=1
    si = int(np.where(st_u == 100.0)[0][0])
    fi = int(np.where(f_u == 1.0)[0][0])
    assert np.isclose(T_u[si, fi, 0, 1], 5.0)
    # XX slot remains NaN
    assert np.isnan(T_u[si, fi, 0, 0])

    # INTERSECTION → only shared freq {2}
    T_i, f_i, _ = r.to_tensor(var="rho", align="intersection")
    assert np.allclose(f_i, [2.0])
    assert T_i.shape == (2, 1, 2, 2)


def test_resistivity_write_block_has_banner_and_unit_meta():
    df = pd.DataFrame(
        {
            "station": [0.0],
            "freq": [1.0],
            "comp": ["ExHy"],
            "ARes.mag": [123.0],
        }
    )
    r = Resistivity.from_avg((df, {"Unit.Rho": "ohm·m"}))
    lines = r.write()
    assert _has_line(lines, lambda s: s.strip().startswith(r"\ $Resistivity Block"))
    assert _has_line(lines, lambda s: s.strip().startswith("$Unit.Rho="))


# --------------------------------- Phase ----------------------------------- #
def test_phase_read_modern_and_unit_convert_roundtrip():
    df = pd.DataFrame(
        {
            "station": [0.0],
            "freq": [1.0],
            "comp": ["ExHy"],
            "Z.phz": [1000.0],  # mrad by default
        }
    )
    p = Phase.from_avg((df, {"Unit.Phase": "mrad"}))

    # Canonical column present
    assert "phase" in p.frame.columns

    # mrad → deg
    p.convert_unit("deg")
    v_deg = float(p.frame["phase"].iloc[0])
    expected_deg = 1000.0 * (180.0 / (np.pi * 1000.0))
    assert np.isclose(v_deg, expected_deg, rtol=0, atol=1e-12)
    assert p.meta["Unit.Phase"].lower() == "deg"

    # deg → mrad
    p.convert_unit("mrad")
    v_mrad = float(p.frame["phase"].iloc[0])
    assert np.isfinite(v_mrad)
    assert p.meta["Unit.Phase"].lower() == "mrad"


@pytest.mark.skipif(pytest.importorskip("xarray") is None, reason="xarray required")
def test_phase_to_xarray_dims_and_var_attrs():
    df = pd.DataFrame(
        {
            "station": [100.0, 100.0, 150.0],
            "freq": [1.0, 2.0, 2.0],
            "comp": ["ExHy", "ExHy", "ExHy"],
            "Phase": [10.0, 20.0, 30.0],  # legacy alias
        }
    )
    p = Phase.from_avg((df, {"Unit.Phase": "mrad"}))
    ds = p.to_xarray()
    assert set(ds.dims) == {"station", "freq", "comp"}
    assert "phase" in ds.variables
    assert str(ds.attrs.get("Unit.Phase", "")).lower() == "mrad"

    # values along station=100 for comp='ExHy'
    vals = ds["phase"].sel(station=100.0, comp="ExHy").values
    vals = np.asarray(vals).reshape(-1)
    # should contain 10 and 20 at the first two freqs
    assert np.isclose(vals[0], 10.0, equal_nan=False)
    assert np.isclose(vals[1], 20.0, equal_nan=False)


def test_phase_to_tensor_single_station_shape_and_slot():
    # Single station; ensure we get (n_freq, 2, 2) and slot mapping
    df = pd.DataFrame(
        {
            "station": [200.0, 200.0],
            "freq": [1.0, 2.0],
            "comp": ["ExHy", "ExHy"],
            "Z.phz": [5.0, 7.5],
        }
    )
    p = Phase.from_avg((df, {}))

    # request single-station explicitly
    T, f, st = p.to_tensor(var="phase", station=200.0)
    assert T.shape == (2, 2, 2)  # (n_freq, 2, 2)
    assert st.size == 0  # single-station path
    assert np.allclose(f, [1.0, 2.0])
    # ExHy → (0,1)
    assert np.isclose(T[0, 0, 1], 5.0)
    assert np.isclose(T[1, 0, 1], 7.5)


def test_phase_write_block_has_banner_and_unit_meta():
    df = pd.DataFrame(
        {
            "station": [0.0],
            "freq": [1.0],
            "comp": ["ExHy"],
            "Phase": [12.3],
        }
    )
    p = Phase.from_avg((df, {}))
    lines = p.write()
    assert _has_line(lines, lambda s: s.strip().startswith(r"\ $Phase Block"))
    assert _has_line(lines, lambda s: s.strip().startswith("$Unit.Phase="))


# --------------------------- Additional coverage ---------------------------- #


def test_resistivity_read_vector_like_list():
    r = Resistivity()
    r.read([1.0, 2.0, 3.0], station=100.0, freq=1.0, comp="ExHy")
    assert list(r.frame["rho"]) == [1.0, 2.0, 3.0]
    assert set(r.frame["station"].unique()) == {100.0}
    assert set(r.frame["comp"].unique()) == {"ExHy"}


def test_resistivity_read_vector_like_numpy_array_defaults():
    r = Resistivity()
    r.read(np.array([1.0, 2.0]))
    assert list(r.frame["rho"]) == [1.0, 2.0]
    assert r.frame["comp"].iloc[0] == "ExHy"


def test_resistivity_read_raises_typeerror_for_invalid_source():
    r = Resistivity()
    with pytest.raises(TypeError):
        r.read({"not": "supported"})


def test_resistivity_read_creates_missing_rho_column():
    df = pd.DataFrame({"station": [1.0], "freq": [1.0], "comp": ["ExHy"]})
    r = Resistivity(verbose=True)
    r.read(df)
    assert "rho" in r.frame.columns
    assert pd.isna(r.frame["rho"].iloc[0])


def test_resistivity_read_injects_missing_station_column():
    df = pd.DataFrame({"freq": [1.0], "ARes.mag": [10.0]})
    r = Resistivity()
    r.read(df)
    assert "station" in r.frame.columns
    assert pd.isna(r.frame["station"].iloc[0])


def test_resistivity_read_raises_when_freq_missing():
    df = pd.DataFrame({"station": [1.0], "ARes.mag": [10.0]})
    r = Resistivity()
    with pytest.raises(AvgDataError, match="freq"):
        r.read(df)


def test_resistivity_to_xarray_merges_extra_attrs():
    df = pd.DataFrame(
        {"station": [0.0], "freq": [1.0], "comp": ["ExHy"], "ARes.mag": [1.0]}
    )
    r = Resistivity.from_avg((df, {}))
    ds = r.to_xarray(attrs={"custom": "value"})
    assert ds.attrs.get("custom") == "value"


def test_phase_read_vector_like_list():
    p = Phase()
    p.read([1.0, 2.0], station=50.0, freq=1.0, comp="ExHy")
    assert list(p.frame["phase"]) == [1.0, 2.0]
    assert set(p.frame["station"].unique()) == {50.0}


def test_phase_read_raises_typeerror_for_invalid_source():
    p = Phase()
    with pytest.raises(TypeError):
        p.read({"not": "supported"})


def test_phase_read_creates_missing_phase_column():
    df = pd.DataFrame({"station": [1.0], "freq": [1.0], "comp": ["ExHy"]})
    p = Phase(verbose=True)
    p.read(df)
    assert "phase" in p.frame.columns
    assert pd.isna(p.frame["phase"].iloc[0])


def test_phase_read_injects_missing_comp_and_station_columns():
    df = pd.DataFrame({"freq": [1.0], "Z.phz": [10.0]})
    p = Phase()
    p.read(df)
    assert "comp" in p.frame.columns
    assert p.frame["comp"].iloc[0] == "ExHy"
    assert "station" in p.frame.columns
    assert pd.isna(p.frame["station"].iloc[0])


def test_phase_read_raises_when_freq_missing():
    df = pd.DataFrame({"station": [1.0], "Z.phz": [10.0]})
    p = Phase()
    with pytest.raises(AvgDataError, match="freq"):
        p.read(df)


def test_phase_to_xarray_merges_extra_attrs():
    df = pd.DataFrame(
        {"station": [0.0], "freq": [1.0], "comp": ["ExHy"], "Z.phz": [1.0]}
    )
    p = Phase.from_avg((df, {}))
    ds = p.to_xarray(attrs={"custom": "value"})
    assert ds.attrs.get("custom") == "value"


def test_phase_convert_unit_noop_when_same_unit():
    df = pd.DataFrame(
        {"station": [0.0], "freq": [1.0], "comp": ["ExHy"], "Z.phz": [10.0]}
    )
    p = Phase.from_avg((df, {"Unit.Phase": "mrad"}))
    before = float(p.frame["phase"].iloc[0])
    p.convert_unit("mrad")
    assert float(p.frame["phase"].iloc[0]) == before


def test_phase_convert_unit_raises_for_invalid_target():
    p = Phase.from_avg(
        (pd.DataFrame({"station": [0.0], "freq": [1.0], "Z.phz": [1.0]}), {})
    )
    with pytest.raises(ValueError, match="mrad.*deg"):
        p.convert_unit("bogus")


def test_phase_convert_unit_noop_when_column_missing():
    p = Phase()
    p._meta = {"Unit.Phase": "mrad"}
    # frame has no "phase" column at all
    p._frame = pd.DataFrame({"station": [0.0], "freq": [1.0]})
    p.convert_unit("deg")  # should not raise
    assert p.meta["Unit.Phase"] == "mrad"  # unchanged: early return


def test_phase_convert_unit_raises_for_unsupported_current_unit():
    df = pd.DataFrame(
        {"station": [0.0], "freq": [1.0], "comp": ["ExHy"], "Z.phz": [1.0]}
    )
    p = Phase.from_avg((df, {"Unit.Phase": "rad"}))
    with pytest.raises(ValueError, match="unsupported conversion"):
        p.convert_unit("deg")


if __name__ == "__main__":  # pragma: no-cover
    pytest.main([__file__])
