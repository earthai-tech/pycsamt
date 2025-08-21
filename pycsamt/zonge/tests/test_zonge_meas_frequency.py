# -*- coding: utf-8 -*-

from __future__ import annotations

import math
import re

import numpy as np
import pandas as pd
import pytest

from pycsamt.zonge.meas import Frequency
from pycsamt.exceptions import FrequencyError


def _has_line(lines, predicate):
    return any(predicate(ln) for ln in lines)


def test_frequency_reads_vector_and_enforces_positive():
    f = Frequency()
    f.read([1.0, 2.0, 4.0])
    assert f.n_unique == 3
    s = str(f)
    assert "unique=3" in s
    assert "Hz" in s or "span=" in s

    # Non-positive must raise
    with pytest.raises(FrequencyError):
        Frequency().read([1.0, 0.0, 2.0])
    with pytest.raises(FrequencyError):
        Frequency().read([-1.0, 1.0])


def test_frequency_reads_dataframe_with_aliases_and_missing_markers():
    df = pd.DataFrame(
        {
            "Station": [100, 100, 100],
            "Freq.": ["1", ".5", "*"],      # legacy quirks: '.5' and '*'
            "Comp": ["ExHy", "ExHy", "ExHy"],
        }
    )
    f = Frequency()
    f.read(df)
    # Parsed numeric: 1.0, 0.5, NaN
    got = pd.to_numeric(f.frame["freq"], errors="coerce").tolist()
    assert math.isclose(got[0], 1.0, rel_tol=0, abs_tol=0)
    assert math.isclose(got[1], 0.5, rel_tol=0, abs_tol=0)
    assert np.isnan(got[2])


def test_frequency_unique_and_by_station_tolerant_dedup():
    # full grid 2 stations × 3 freqs with tiny jitter for the second station
    df = pd.DataFrame(
        {
            "station": [100.0, 100.0, 100.0, 150.0, 150.0, 150.0],
            "freq":    [1.0, 2.0, 4.0, 1.0*(1+1e-12), 2.0, 4.0],
            "comp":    ["ExHy"] * 6,
        }
    )
    f = Frequency.from_avg((df, {}))
    uniq = f.unique()
    assert uniq.size == 3
    assert np.allclose(uniq, np.array([1.0, 2.0, 4.0]))

    per_stn = f.by_station()
    assert set(per_stn.keys()) == {100.0, 150.0}
    assert np.allclose(per_stn[100.0], [1.0, 2.0, 4.0])
    assert np.allclose(per_stn[150.0], [1.0, 2.0, 4.0])


@pytest.mark.skipif(
    pytest.importorskip("xarray") is None, reason="xarray required")
def test_frequency_to_xarray_shapes_and_values():
    # full grid: 2 stations × 3 freqs × 1 comp
    df = pd.DataFrame(
        {
            "station": [100.0, 100.0, 100.0, 150.0, 150.0, 150.0],
            "freq":    [1.0, 2.0, 4.0, 1.0, 2.0, 4.0],
            "comp":    ["ExHy"] * 6,
        }
    )
    f = Frequency.from_avg((df, {}))
    ds = f.to_xarray()

    # dims exist with expected sizes
    assert set(ds.dims) == {"station", "freq", "comp"}
    assert ds.dims["station"] == 2
    assert ds.dims["freq"] == 3
    assert ds.dims["comp"] == 1

    # the 'freq' variable along station=100, comp='ExHy' matches [1,2,4]
    vals = ds["freq_grid"].sel(station=100.0, comp="ExHy").values
    # squeeze to 1-D for comparison
    vals = np.asarray(vals).reshape(-1)
    assert np.allclose(vals, [1.0, 2.0, 4.0])


def test_frequency_write_emits_csv_with_context_and_meta():
    df = pd.DataFrame(
        {
            "station": [100.0, 100.0, 150.0],
            "freq":    [1.0, 2.0, 4.0],
            "comp":    ["ExHy", "ExHy", "ExHy"],
        }
    )
    f = Frequency.from_avg((df, {}))
    lines = f.write()

    # Title banner present
    assert _has_line(lines, lambda ln: ln.strip(
        ).startswith("\\ $Frequency Block"))

    # Default meta exported
    assert _has_line(lines, lambda ln: ln.strip(
        ).startswith("$Unit.Freq=Hz"))

    # Written timestamp is present
    assert _has_line(lines, lambda ln: ln.strip(
        ).startswith("$Written="))

    # CSV header present
    assert _has_line(lines, lambda ln: ln.strip(
        ) == "station,freq,comp")

    # And we have as many CSV rows as data rows
    # find header index and count subsequent rows
    try:
        header_idx = next(i for i, ln in enumerate(
            lines) if ln.strip() == "station,freq,comp")
    except StopIteration:
        pytest.fail("CSV header not found in Frequency.write() output")

    csv_rows = [ln for ln in lines[header_idx + 1 :] if ln.strip()]
    assert len(csv_rows) == len(df)

    # spot-check one row content (tolerant on float formatting)
    assert any(re.match(r"^100(\.0)?,1(\.0)?,ExHy$", ln.strip()) for ln in csv_rows)


if __name__=='__main__': # pragma: no-cover 
   pytest.main( [__file__])