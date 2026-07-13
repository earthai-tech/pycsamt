# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

import numpy as np
import pandas as pd
import pytest

from pycsamt.zonge.var_pc import (
    PcEmag,
    PcHmag,
    PcRho,
)
from pycsamt.zonge.var_std import (
    SEphz,
    # SHphz # test this also
    SPhz,
)


def _has_line(lines, pred):
    return any(pred(ln) for ln in lines)


# ------------------------- PHASE STDEV: SPhz/SEphz/SHphz --------------------
def test_sphz_read_legacy_and_write_block_has_unit_meta():
    df = pd.DataFrame(
        {
            "station": [100.0, 100.0, 150.0],
            "freq": [1.0, 2.0, 2.0],
            "comp": ["ExHy", "ExHy", "ExHy"],
            # legacy name
            "sPhz": [12.0, 9.5, 7.0],
        }
    )
    s = SPhz.from_avg((df, {}))
    # internal, canonical column name exists
    assert "s_phz" in s.frame.columns
    assert np.allclose(s.frame["s_phz"].to_numpy(), [12.0, 9.5, 7.0])

    lines = s.write()
    # banner + Unit.Phase meta present
    assert _has_line(
        lines, lambda ln: ln.strip().startswith("\\ $Phase-Stdev")
    )
    assert _has_line(lines, lambda ln: ln.strip().startswith("$Unit.Phase="))


def test_phase_convert_unit_roundtrip_mrad_deg():
    df = pd.DataFrame(
        {
            "station": [0.0],
            "freq": [1.0],
            "comp": ["ExHy"],
            "Z.perr": [1000.0],
        }
    )
    s = SPhz.from_avg((df, {"Unit.Phase": "mrad"}))
    # mrad -> deg
    s.convert_unit("deg")
    v_deg = float(s.frame["s_phz"].iloc[0])
    expected_deg = 1000.0 * (180.0 / (np.pi * 1000.0))
    assert np.isclose(v_deg, expected_deg, rtol=0, atol=1e-8)
    assert s.meta["Unit.Phase"].lower() == "deg"

    # deg -> mrad
    s.convert_unit("mrad")
    v_mrad = float(s.frame["s_phz"].iloc[0])
    assert np.isfinite(v_mrad)
    assert s.meta["Unit.Phase"].lower() == "mrad"


@pytest.mark.skipif(
    pytest.importorskip("xarray") is None, reason="xarray required"
)
def test_phase_to_xarray_dims_and_var():
    df = pd.DataFrame(
        {
            "station": [100.0, 100.0, 150.0],
            "freq": [1.0, 2.0, 2.0],
            "comp": ["ExHy", "ExHy", "ExHy"],
            "E.perr": [3.0, 2.5, 1.0],
        }
    )
    e = SEphz.from_avg((df, {"Unit.Phase": "mrad"}))
    ds = e.to_xarray()
    assert set(ds.dims) == {"station", "freq", "comp"}
    assert "s_ephz" in ds.variables
    # attrs carry Unit.Phase
    assert str(ds.attrs.get("Unit.Phase", "")).lower() == "mrad"


def test_phase_to_tensor_like_union_and_intersection():
    rows = [
        # station 100: 1, 2 (ExHy)
        dict(station=100.0, freq=1.0, comp="ExHy", sPhz=2.0),
        dict(station=100.0, freq=2.0, comp="ExHy", sPhz=4.0),
        # station 150: 2, 4 (ExHy)
        dict(station=150.0, freq=2.0, comp="ExHy", sPhz=40.0),
        dict(station=150.0, freq=4.0, comp="ExHy", sPhz=80.0),
    ]
    df = pd.DataFrame.from_records(rows)
    z = SPhz.from_avg((df, {}))

    T_u, f_u, st_u = z.to_tensor_like(align="union")
    # stations x freqs x 2 x 2
    assert T_u.shape == (2, 3, 2, 2)
    assert np.allclose(np.sort(st_u), [100.0, 150.0])
    assert np.allclose(f_u, [1.0, 2.0, 4.0])

    # value lands at (ExHy) → (0,1)
    si = int(np.where(st_u == 100.0)[0][0])
    fi = int(np.where(f_u == 1.0)[0][0])
    assert np.isclose(T_u[si, fi, 0, 1], 2.0)
    # XX remains NaN
    assert np.isnan(T_u[si, fi, 0, 0])

    T_i, f_i, _ = z.to_tensor_like(align="intersection")
    # only shared frequency (2.0)
    assert np.allclose(f_i, [2.0])
    assert T_i.shape == (2, 1, 2, 2)


# ------------------------- PERCENT-ERROR: PcEmag/PcHmag/PcRho ---------------
def test_pc_emag_read_modern_and_write_block_banner():
    df = pd.DataFrame(
        {
            "station": [100.0, 150.0],
            "freq": [1.0, 2.0],
            "comp": ["ExHy", "ExHy"],
            "E.%err": [5.0, 10.0],
        }
    )
    p = PcEmag.from_avg((df, {}))
    assert "pc_emag" in p.frame.columns
    assert np.allclose(p.frame["pc_emag"].to_numpy(), [5.0, 10.0])

    lines = p.write()
    assert _has_line(
        lines, lambda ln: ln.strip().startswith("\\ $ Percent |E| Variation")
    )


def test_pc_hmag_read_legacy_and_to_tensor_like_position():
    df = pd.DataFrame(
        {
            "station": [100.0, 100.0],
            "freq": [1.0, 2.0],
            "comp": ["ExHy", "ExHy"],
            "%Hmag": [3.0, 6.0],
        }
    )
    p = PcHmag.from_avg((df, {}))
    T, f, st = p.to_tensor_like()
    # single station path → (n_freq, 2, 2); because no station col?
    # we *do* have station; multi-station with 1 station => (1, n_freq, 2, 2)
    assert T.shape == (1, 2, 2, 2)
    assert np.allclose(f, [1.0, 2.0])
    # Only the ExHy slot is filled
    assert np.isclose(T[0, 0, 0, 1], 3.0)
    assert np.isclose(T[0, 1, 0, 1], 6.0)
    assert np.isnan(T[0, 0, 0, 0])


@pytest.mark.skipif(
    pytest.importorskip("xarray") is None, reason="xarray required"
)
def test_pc_rho_to_xarray_dims_and_var():
    df = pd.DataFrame(
        {
            "station": [100.0, 100.0, 150.0],
            "freq": [1.0, 2.0, 2.0],
            "comp": ["ExHy", "ExHy", "ExHy"],
            "%Rho": [1.1, 2.2, 3.3],
        }
    )
    r = PcRho.from_avg((df, {}))
    ds = r.to_xarray()
    assert set(ds.dims) == {"station", "freq", "comp"}
    assert "pc_rho" in ds.variables
    # grid contains the expected percentage values at known coords
    v = ds["pc_rho"].sel(station=100.0, comp="ExHy").values
    v = np.asarray(v).reshape(-1)
    assert np.allclose(v[:2], [1.1, 2.2])


if __name__ == "__main__":  # pragma: no-cover
    pytest.main([__file__])
