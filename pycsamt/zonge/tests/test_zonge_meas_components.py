# -*- coding: utf-8 -*-
# Tests for pycsamt.zonge.meas
from __future__ import annotations

import math
import numpy as np
import pandas as pd
import pytest

from pycsamt.exceptions import InputError  # type: ignore
from pycsamt.zonge.meas import CompMeas, Amps

# ---------------------------------------------------------------------
# CompMeas
# ---------------------------------------------------------------------

def test_compmeas_defaults_when_missing_column():
    df = pd.DataFrame(
        {
            "station": [100, 100, 150, 150],
            "freq": [1.0, 2.0, 1.0, 2.0],
            # no 'comp' column -> should default to ExHy
        }
    )
    cm = CompMeas.from_avg((df, {}))
    assert cm.frame.shape == (4, 3)  # station,freq,comp
    assert cm.unique == ["ExHy"]
    # write should emit a compact CSV block
    lines = cm.write()
    assert isinstance(lines, list) and len(lines) > 0
    # the CSV header line should include our fields
    csv_header = next((ln for ln in lines if "station" in ln and "comp" in ln), "")
    assert "station" in csv_header and "comp" in csv_header


def test_compmeas_normalises_and_validates_labels():
    df = pd.DataFrame(
        {
            "station": [1, 1, 2, 2],
            "freq": [1.0, 2.0, 1.0, 2.0],
            "comp": ["exhy", "EXHX", "EyHy", "EyHx"],  # mixed case
        }
    )
    cm = CompMeas.from_avg((df, {}))
    assert set(cm.unique) == {"ExHy", "ExHx", "EyHy", "EyHx"}

    # now feed a bad label and expect a validation error
    df_bad = df.copy()
    df_bad.loc[0, "comp"] = "ExZZ"
    with pytest.raises(InputError):
        CompMeas.from_avg((df_bad, {}))

# ---------------------------------------------------------------------
# Amps
# ---------------------------------------------------------------------
def test_amps_reads_modern_and_legacy_columns_and_stats():
    # values include numeric strings, '*' (missing), and blank
    df = pd.DataFrame(
        {
            "station": [10, 10, 20, 20],
            "freq": [1.0, 2.0, 1.0, 2.0],
            "comp": ["ExHy"] * 4,
            "Tx.Amp": ["13", "*", " 8.25 ", ""],  # legacy name
        }
    )

    amps = Amps.from_avg((df, {}))
    out = amps.to_frame()

    # coercion to floats with NaNs
    assert out["amps"].dtype.kind in ("f", "d")
    assert math.isclose(out["amps"].iloc[0], 13.0)
    assert np.isnan(out["amps"].iloc[1])
    assert math.isclose(out["amps"].iloc[2], 8.25)
    assert np.isnan(out["amps"].iloc[3])

    # stats computed only on finite values
    st = amps.stats
    assert st.count == 2
    assert math.isclose(st.vmin or 0.0, 8.25)
    assert math.isclose(st.vmax or 0.0, 13.0)
    assert math.isclose(st.mean or 0.0, (13.0 + 8.25) / 2.0)
    assert math.isclose(st.median or 0.0, (13.0 + 8.25) / 2.0)


@pytest.mark.skipif(
    pytest.importorskip("xarray", reason="xarray not installed") is None,
    reason="xarray not installed",
)
def test_amps_to_xarray_shapes_and_values():
    # full grid: 2 stations × 3 freqs × 1 comp
    df = pd.DataFrame(
        {
            "station": [100, 100, 100, 150, 150, 150],
            "freq": [1.0, 2.0, 4.0, 1.0, 2.0, 4.0],
            "comp": ["ExHy"] * 6,
            "amps": [10.0, 11.0, 12.0, 9.0, 8.0, 7.0],
        }
    )
    amps = Amps.from_avg((df, {}))
    ds = amps.to_xarray()
    # dims
    assert set(ds.dims) == {"station", "freq", "comp"}
    assert ds.sizes["station"] == 2
    assert ds.sizes["freq"] == 3
    assert ds.sizes["comp"] == 1
    # data var present and values line up
    arr = ds["amps"].sel(comp="ExHy").values
    # station axis is usually sorted ascending [100, 150]
    # verify a couple of positions
    f_idx = list(ds.coords["freq"].values).index(2.0)
    s_idx = list(ds.coords["station"].values).index(100.0)
    assert math.isclose(arr[s_idx, f_idx], 11.0)
    s_idx2 = list(ds.coords["station"].values).index(150.0)
    f_idx0 = list(ds.coords["freq"].values).index(1.0)
    assert math.isclose(arr[s_idx2, f_idx0], 9.0)


def test_amps_write_emits_csv_with_context():
    df = pd.DataFrame(
        {
            "station": [1, 2],
            "freq": [1.0, 1.0],
            "comp": ["ExHy", "ExHy"],
            "amps": [5.0, 7.5],
        }
    )
    amps = Amps.from_avg((df, {}))
    lines = amps.write()
    assert isinstance(lines, list) and len(lines) > 0
    # should contain a CSV header with our columns
    header = next((ln for ln in lines if "station" in ln and "amps" in ln), "")
    assert "station" in header and "amps" in header


if __name__=='__main__': # pragma: no-cover 
   pytest.main( [__file__])