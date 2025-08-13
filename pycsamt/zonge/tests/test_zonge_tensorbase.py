# -*- coding: utf-8 -*-
# Tests for pycsamt.zonge.tensor.TensorBase

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from pycsamt.zonge.tensor import TensorBase
from pycsamt.exceptions import AvgDataError


class DummyTensor(TensorBase):
    """
    Minimal concrete subclass so we can feed a tidy DataFrame
    directly into `_frame` and use TensorBase helpers.
    """
    def read(self, source, meta=None, **kws):
        if not isinstance(source, pd.DataFrame):
            raise TypeError("DummyTensor.read expects a DataFrame")
        self._frame = source.copy()
        self._meta = dict(meta or {})

    def write(self):
        return []


def _df_single_station_csamt():
    """
    station=100, freqs=[1,2,4], comps=ExHx/ExHy/EyHx/EyHy, var='rho'
    values encode position so we can assert placement:
      ExHx=f+0.10, ExHy=f+0.20, EyHx=f+0.30, EyHy=f+0.40
    """
    rows = []
    for f in [1.0, 2.0, 4.0]:
        rows += [
            {"station": 100.0, "freq": f, "comp": "ExHx", "rho": f + 0.10},
            {"station": 100.0, "freq": f, "comp": "ExHy", "rho": f + 0.20},
            {"station": 100.0, "freq": f, "comp": "EyHx", "rho": f + 0.30},
            {"station": 100.0, "freq": f, "comp": "EyHy", "rho": f + 0.40},
        ]
    return pd.DataFrame.from_records(rows)


def _df_single_station_mt():
    """
    station=100, freqs=[1,2], comps=Zxx/Zxy/Zyx/Zyy, var='zabs'
      Zxx=f+1, Zxy=f+2, Zyx=f+3, Zyy=f+4
    """
    rows = []
    for f in [1.0, 2.0]:
        rows += [
            {"station": 100.0, "freq": f, "comp": "Zxx", "zabs": f + 1},
            {"station": 100.0, "freq": f, "comp": "Zxy", "zabs": f + 2},
            {"station": 100.0, "freq": f, "comp": "Zyx", "zabs": f + 3},
            {"station": 100.0, "freq": f, "comp": "Zyy", "zabs": f + 4},
        ]
    return pd.DataFrame.from_records(rows)


def test_to_tensor_single_station_csamt_positions():
    df = _df_single_station_csamt()
    obj = DummyTensor.from_avg((df, {}))  # uses base.from_avg -> read(df)

    T, freqs, stations = obj.to_tensor(var="rho", station=100.0)
    assert T.shape == (3, 2, 2)   # (n_freq, 2, 2)
    assert stations.size == 0
    assert np.allclose(freqs, [1.0, 2.0, 4.0])

    # freq=1 is index 0
    block = T[0]
    # positions: [0,0]=ExHx, [0,1]=ExHy, [1,0]=EyHx, [1,1]=EyHy
    assert np.isclose(block[0, 0], 1.10)
    assert np.isclose(block[0, 1], 1.20)
    assert np.isclose(block[1, 0], 1.30)
    assert np.isclose(block[1, 1], 1.40)


def test_to_tensor_single_station_mt_positions():
    df = _df_single_station_mt()
    obj = DummyTensor.from_avg((df, {}))

    T, freqs, _ = obj.to_tensor(var="zabs", station=100.0)
    assert T.shape == (2, 2, 2)
    assert np.allclose(freqs, [1.0, 2.0])

    # freq=2 (index 1)
    b = T[1]
    assert np.isclose(b[0, 0], 3.0)  # Zxx: 2+1
    assert np.isclose(b[0, 1], 4.0)  # Zxy: 2+2
    assert np.isclose(b[1, 0], 5.0)  # Zyx: 2+3
    assert np.isclose(b[1, 1], 6.0)  # Zyy: 2+4


def test_to_tensor_multi_station_union_and_intersection():
    # Two stations, ragged frequencies, only ExHy present (others become NaN)
    rows = [
        # station 100: freqs 1,2
        {"station": 100.0, "freq": 1.0, "comp": "ExHy", "rho": 10.0},
        {"station": 100.0, "freq": 2.0, "comp": "ExHy", "rho": 20.0},
        # station 150: freqs 2,4
        {"station": 150.0, "freq": 2.0, "comp": "ExHy", "rho": 200.0},
        {"station": 150.0, "freq": 4.0, "comp": "ExHy", "rho": 400.0},
    ]
    df = pd.DataFrame.from_records(rows)
    obj = DummyTensor.from_avg((df, {}))

    # UNION → freqs = [1,2,4]
    T_u, f_u, st_u = obj.to_tensor(var="rho", align="union")
    assert T_u.shape == (2, 3, 2, 2)
    assert np.allclose(np.sort(st_u), [100.0, 150.0])
    assert np.allclose(f_u, [1.0, 2.0, 4.0])

    # Only [0,1] (ExHy) has values where present
    # station=100 at freq=1 → 10.0
    i_st = int(np.where(st_u == 100.0)[0][0])
    i_f1 = int(np.where(f_u == 1.0)[0][0])
    assert np.isclose(T_u[i_st, i_f1, 0, 1], 10.0)
    # XX position stays NaN
    assert np.isnan(T_u[i_st, i_f1, 0, 0])

    # INTERSECTION → only freq=2
    T_i, f_i, _ = obj.to_tensor(var="rho", align="intersection")
    assert T_i.shape == (2, 1, 2, 2)
    assert np.allclose(f_i, [2.0])
    # Both stations have ExHy at freq=2
    i_f2 = 0
    i_st100 = int(np.where(st_u == 100.0)[0][0])
    i_st150 = int(np.where(st_u == 150.0)[0][0])
    assert np.isclose(T_u[i_st100, int(np.where(f_u == 2.0)[0][0]), 0, 1], 20.0)
    assert np.isclose(T_u[i_st150, int(np.where(f_u == 2.0)[0][0]), 0, 1], 200.0)


def test_to_tensor_invalid_align_raises():
    df = _df_single_station_csamt()
    obj = DummyTensor.from_avg((df, {}))
    with pytest.raises(ValueError):
        obj.to_tensor(var="rho", align="nope")


def test_duplicate_rows_are_averaged():
    # Same (station,freq,comp) appears twice -> average used
    df = pd.DataFrame(
        {
            "station": [100.0, 100.0],
            "freq":    [1.0, 1.0],
            "comp":    ["ExHy", "ExHy"],
            "rho":     [10.0, 14.0],  # mean=12.0
        }
    )
    obj = DummyTensor.from_avg((df, {}))
    T, f, _ = obj.to_tensor(var="rho", station=100.0)
    assert T.shape == (1, 2, 2)
    assert np.allclose(f, [1.0])
    assert np.isclose(T[0, 0, 1], 12.0)  # [0,1] = ExHy


def test_from_tensor_roundtrip_single_station_csamt():
    # Build a simple tensor block and reconstruct a tidy frame
    freqs = np.array([1.0, 2.0])
    T = np.zeros((2, 2, 2), float)
    # freq=1
    T[0, 0, 0] = 1.1  # ExHx
    T[0, 0, 1] = 1.2  # ExHy
    T[0, 1, 0] = 1.3  # EyHx
    T[0, 1, 1] = 1.4  # EyHy
    # freq=2
    T[1, 0, 0] = 2.1
    T[1, 0, 1] = 2.2
    T[1, 1, 0] = 2.3
    T[1, 1, 1] = 2.4

    df_back = TensorBase.from_tensor(T, freqs, var="rho", comp_style="csamt")
    assert set(df_back.columns) >= {"station", "freq", "comp", "rho"}
    # sample checks
    v11 = float(df_back.query("freq==1.0 and comp=='ExHy'")["rho"].iloc[0])
    v22 = float(df_back.query("freq==2.0 and comp=='EyHy'")["rho"].iloc[0])
    assert np.isclose(v11, 1.2)
    assert np.isclose(v22, 2.4)


@pytest.mark.skipif(pytest.importorskip("xarray") is None, reason="xarray required")
def test_to_xarray_tensor_dims_and_coords_single_and_multi():
    df = _df_single_station_csamt()
    obj = DummyTensor.from_avg((df, {}))

    # single-station DA
    da1 = obj.to_xarray_tensor(var="rho", station=100.0)
    assert da1.dims == ("freq", "e", "h")
    assert np.all(da1.coords["e"].values == np.array(["Ex", "Ey"]))
    assert np.all(da1.coords["h"].values == np.array(["Hx", "Hy"]))
    # value spot check
    assert float(da1.sel(freq=1.0, e="Ex", h="Hy")) == pytest.approx(1.2)

    # multi-station DA
    da2 = obj.to_xarray_tensor(var="rho")
    assert da2.dims == ("station", "freq", "e", "h")
    assert np.all(np.unique(da2.coords["station"]) == np.array([100.0]))

if __name__=='__main__': # pragma: no-cover 
   pytest.main( [__file__])