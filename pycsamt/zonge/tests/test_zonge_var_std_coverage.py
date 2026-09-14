# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from pycsamt.exceptions import AvgDataError
from pycsamt.zonge.var_std import (
    PhaseStdBase,
    SEphz,
    SHphz,
    SPhz,
    _norm_comp,
    _to_num,
)


def test_read_raises_runtimeerror_when_var_name_unset():
    base = PhaseStdBase()
    with pytest.raises(RuntimeError, match="subclass constants"):
        base.read(pd.DataFrame({"freq": [1.0]}))


def test_read_vector_like_list():
    s = SPhz()
    s.read([1.0, 2.0, 3.0], station=10.0, freq=1.0, comp="ExHy")
    assert list(s.frame["s_phz"]) == [1.0, 2.0, 3.0]
    assert set(s.frame["station"].unique()) == {10.0}
    assert set(s.frame["comp"].unique()) == {"ExHy"}


def test_read_vector_like_numpy_array_defaults():
    s = SPhz()
    s.read(np.array([1.0, 2.0]))
    assert list(s.frame["s_phz"]) == [1.0, 2.0]
    assert s.frame["comp"].iloc[0] == "ExHy"


def test_read_raises_typeerror_for_invalid_source():
    s = SPhz()
    with pytest.raises(TypeError):
        s.read({"not": "supported"})


def test_read_creates_missing_var_column_verbose():
    df = pd.DataFrame({"station": [1.0], "freq": [1.0], "comp": ["ExHy"]})
    s = SPhz(verbose=True)
    s.read(df)
    assert "s_phz" in s.frame.columns
    assert pd.isna(s.frame["s_phz"].iloc[0])


def test_read_injects_missing_comp_station_freq():
    df = pd.DataFrame({"sPhz": [1.0]})
    s = SPhz()
    s.read(df)
    assert s.frame["comp"].iloc[0] == "ExHy"
    assert pd.isna(s.frame["station"].iloc[0])
    assert pd.isna(s.frame["freq"].iloc[0])


def test_to_xarray_raises_when_empty():
    s = SPhz()
    with pytest.raises(AvgDataError, match="empty frame"):
        s.to_xarray()


def test_to_xarray_injects_comp_when_missing():
    s = SPhz()
    s._frame = pd.DataFrame({"station": [1.0], "freq": [1.0], "s_phz": [2.0]})
    s._meta = {}
    ds = s.to_xarray(coords=("station", "freq", "comp"))
    assert "comp" in ds.dims


def test_to_xarray_raises_when_no_coordinate_columns_found():
    s = SPhz()
    s._frame = pd.DataFrame({"s_phz": [1.0]})
    s._frame["comp"] = "ExHy"  # comp gets injected but coords= () means none
    s._meta = {}
    with pytest.raises(AvgDataError, match="no coordinate columns"):
        s.to_xarray(coords=())


def test_to_xarray_raises_when_var_missing_from_frame():
    s = SPhz()
    s._frame = pd.DataFrame({"station": [1.0], "freq": [1.0], "comp": ["ExHy"]})
    s._meta = {}
    with pytest.raises(AvgDataError, match="s_phz"):
        s.to_xarray()


def test_to_xarray_averages_duplicate_coordinate_rows():
    df = pd.DataFrame(
        {
            "station": [100.0, 100.0],
            "freq": [1.0, 1.0],
            "comp": ["ExHy", "ExHy"],
            "sPhz": [2.0, 4.0],
        }
    )
    s = SPhz.from_avg((df, {}))
    ds = s.to_xarray()
    val = ds["s_phz"].sel(station=100.0, freq=1.0, comp="ExHy").item()
    assert np.isclose(val, 3.0)


def test_to_xarray_merges_extra_attrs():
    df = pd.DataFrame(
        {"station": [0.0], "freq": [1.0], "comp": ["ExHy"], "sPhz": [1.0]}
    )
    s = SPhz.from_avg((df, {}))
    ds = s.to_xarray(attrs={"custom": "value"})
    assert ds.attrs.get("custom") == "value"


def test_convert_unit_noop_when_same_unit():
    df = pd.DataFrame(
        {"station": [0.0], "freq": [1.0], "comp": ["ExHy"], "sPhz": [10.0]}
    )
    s = SPhz.from_avg((df, {"Unit.Phase": "mrad"}))
    before = float(s.frame["s_phz"].iloc[0])
    s.convert_unit("mrad")
    assert float(s.frame["s_phz"].iloc[0]) == before


def test_convert_unit_raises_for_invalid_target():
    df = pd.DataFrame({"station": [0.0], "freq": [1.0], "sPhz": [1.0]})
    s = SPhz.from_avg((df, {}))
    with pytest.raises(ValueError, match="mrad.*deg"):
        s.convert_unit("bogus")


def test_convert_unit_noop_when_column_missing():
    s = SPhz()
    s._meta = {"Unit.Phase": "mrad"}
    s._frame = pd.DataFrame({"station": [0.0], "freq": [1.0]})
    s.convert_unit("deg")
    assert s.meta["Unit.Phase"] == "mrad"


def test_convert_unit_raises_for_unsupported_current_unit():
    df = pd.DataFrame(
        {"station": [0.0], "freq": [1.0], "comp": ["ExHy"], "sPhz": [1.0]}
    )
    s = SPhz.from_avg((df, {"Unit.Phase": "rad"}))
    with pytest.raises(ValueError, match="unsupported conversion"):
        s.convert_unit("deg")


def test_write_returns_csv_block_with_banner_and_unit_meta():
    df = pd.DataFrame(
        {"station": [0.0], "freq": [1.0], "comp": ["ExHy"], "sPhz": [1.0]}
    )
    s = SPhz.from_avg((df, {"Unit.Phase": "mrad"}))
    lines = s.write()
    assert isinstance(lines, list) and len(lines) > 0
    assert any("Phase-Stdev" in ln for ln in lines)
    assert any("Unit.Phase" in ln for ln in lines)


def test_to_tensor_like_places_value_at_component_slot():
    df = pd.DataFrame(
        {
            "station": [100.0, 100.0],
            "freq": [1.0, 2.0],
            "comp": ["ExHy", "ExHy"],
            "sPhz": [2.0, 4.0],
        }
    )
    s = SPhz.from_avg((df, {}))
    T, f, st = s.to_tensor_like(align="union")
    assert T.shape == (1, 2, 2, 2)
    assert np.allclose(f, [1.0, 2.0])
    assert np.isclose(T[0, 0, 0, 1], 2.0)
    assert np.isnan(T[0, 0, 0, 0])


def test_convert_unit_mrad_to_deg_and_back_applies_factor():
    df = pd.DataFrame(
        {"station": [0.0], "freq": [1.0], "comp": ["ExHy"], "sPhz": [1000.0]}
    )
    s = SPhz.from_avg((df, {"Unit.Phase": "mrad"}))
    s.convert_unit("deg")
    expected_deg = 1000.0 * (180.0 / (np.pi * 1000.0))
    assert np.isclose(float(s.frame["s_phz"].iloc[0]), expected_deg)
    assert s.meta["Unit.Phase"] == "deg"

    s.convert_unit("mrad")
    assert np.isclose(float(s.frame["s_phz"].iloc[0]), 1000.0)
    assert s.meta["Unit.Phase"] == "mrad"


def test_str_representation():
    df = pd.DataFrame(
        {"station": [0.0], "freq": [1.0], "comp": ["ExHy"], "sPhz": [1.0]}
    )
    s = SPhz.from_avg((df, {}))
    text = str(s)
    assert "SPhz" in text and "s_phz" in text


def test_sephz_and_shphz_basic_read():
    df_e = pd.DataFrame({"station": [1.0], "freq": [1.0], "E.perr": [1.0]})
    e = SEphz.from_avg((df_e, {}))
    assert "s_ephz" in e.frame.columns

    df_h = pd.DataFrame({"station": [1.0], "freq": [1.0], "H.perr": [1.0]})
    h = SHphz.from_avg((df_h, {}))
    assert "s_hphz" in h.frame.columns


def test_to_num_none_returns_nan():
    assert np.isnan(_to_num(None))


def test_to_num_sentinel_strings_return_nan():
    for s in ("", "*", "nan", "NaN", "None", "null"):
        assert np.isnan(_to_num(s))


def test_to_num_valid_string_converts():
    assert _to_num("3.14") == 3.14


def test_to_num_invalid_string_returns_nan():
    assert np.isnan(_to_num("not-a-number"))


def test_norm_comp_none_returns_default():
    assert _norm_comp(None) == "ExHy"


def test_norm_comp_blank_returns_default():
    assert _norm_comp("   ") == "ExHy"


def test_norm_comp_strips_whitespace():
    assert _norm_comp("  ExHx  ") == "ExHx"
