# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later

import numpy as np
import pandas as pd
import pytest

from pycsamt.exceptions import StationError
from pycsamt.zonge.survey import Station


def make_grid_df(stations, reps=3):
    """Helper: build a tidy frame with repeated station rows."""
    s = np.repeat(np.asarray(stations, dtype=float), reps)
    # add a dummy freq column to mimic real AVG tables
    f = np.tile(np.arange(reps, dtype=float), len(stations))
    return pd.DataFrame({"station": s, "freq": f})


def test_read_from_dataframe_grid_defaults_and_keywords():
    df = make_grid_df([0, 50, 100], reps=3)

    st = Station()
    st.read(df, meta={"Unit.Length": "m"})

    # unique values and names
    assert st.n_unique == 3
    assert list(st.values) == [0.0, 50.0, 100.0]
    assert st.names == ["S00", "S01", "S02"]

    # index maps: each station has 3 rows
    assert all(len(ix) == 3 for ix in st.index_by_value.values())
    assert set(st.index_by_name) == set(["S00", "S01", "S02"])

    # increment detection and keyword derivation
    assert st.increment == 50.0
    kw = st.to_keywords()
    assert kw["Stn.Beg"] == 0.0
    assert kw["Stn.Left"] == 0.0
    assert kw["Stn.Right"] == 100.0
    assert kw["Stn.Inc"] == 50.0
    assert kw["Stn.GdpBeg"] == 0.0
    assert kw["Stn.GdpInc"] == 50.0

    # write block contains header and CSV with both columns
    lines = st.write()
    blob = "\n".join(lines)
    assert r"\ Station Geometry" in blob
    assert "$Stn.Beg=0.0" in blob
    assert "station,station_m" in blob  # CSV header


def test_read_from_array_feet_with_normalize_and_station_m():
    arr = [100, 100, 150, 150]  # feet
    st = Station(normalize=True)
    st.read(arr, unit="ft")

    # after normalize: 0 and 50 (still in feet for 'station')
    assert list(np.unique(st._frame["station"])) == [0.0, 50.0]

    # station_m is metres view (from feet) and honors normalization
    mvals = np.unique(st._frame["station_m"])
    assert np.isclose(mvals.min(), 0.0)
    assert np.isclose(mvals.max(), 50.0 / 3.280839895, rtol=0, atol=1e-9)


def test_ragged_enforcement_raises_when_disabled():
    # ragged: station 0 has 2 rows, station 50 has 3 rows
    df = pd.DataFrame({"station": [0, 0, 50, 50, 50]})
    st = Station()
    with pytest.raises(StationError):
        st.read(df, allow_ragged=False)


def test_custom_names_and_index_maps():
    df = make_grid_df([10, 20], reps=2)
    st = Station()
    st.read(df, names=["A", "B"])

    assert st.names == ["A", "B"]
    assert "A" in st.index_by_name and "B" in st.index_by_name
    # each has 2 rows
    assert all(len(ix) == 2 for ix in st.index_by_name.values())


def test_to_keywords_without_increment_for_irregular_spacing():
    # irregular spacing: 0, 30, 80 → no constant increment
    df = make_grid_df([0, 30, 80], reps=1)
    st = Station()
    st.read(df)

    assert st.increment is None
    kw = st.to_keywords()
    assert kw["Stn.Beg"] == 0.0
    assert kw["Stn.Left"] == 0.0
    assert kw["Stn.Right"] == 80.0
    assert "Stn.Inc" not in kw  # not grid-like


def test_label_map_and_str_repr_are_sane():
    df = make_grid_df([5, 25, 55], reps=2)
    st = Station()
    st.read(df)

    lm = st.label_map()
    assert set(map(float, lm.keys())) == set([5.0, 25.0, 55.0])
    assert set(lm.values()) == set(["S00", "S01", "S02"])

    # repr/str should include n and unit keywords
    s = str(st)
    assert "Station(" in s
    assert "n=3" in s
    assert st.unit in s


def test_from_avg_tuple_pathway_works_without_filesystem():
    # Use the AVGComponentBase.from_avg tuple branch: (df, meta)
    df = make_grid_df([0, 25, 50], reps=2)
    meta = {"Unit.Length": "m"}
    st = Station.from_avg((df, meta))  # type: ignore[arg-type]

    assert st.n_unique == 3
    assert st.increment == 25.0
    assert st.to_keywords()["Stn.Inc"] == 25.0

if __name__=='__main__': # pragma: no-cover
   pytest.main( [__file__])
