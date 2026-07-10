# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for :class:`pycsamt.transformers.TStoEDI`."""
from __future__ import annotations

import numpy as np
import pytest

from pycsamt.seg.collection import EDICollection
from pycsamt.transformers import TransformResult, TStoEDI
from pycsamt.ts import TSData

# ------------------------------------------------------------- fixtures
_LIMS_HEADER = """\
# time series file from mp2ts
#
# Ex line length (m):      50.00
# Ey line length (m):      60.00
#
>INFO_START:
>STATION   :tst01
>LATITUDE  : -30.5
>LONGITUDE :  20.25
>ELEVATION :  100.
>COORD_SYS :MAGNETIC NORTH
>NCHAN     : 5
>CHAN_1    :HX
>AZIM_1    :  0.
>UNITS_1   :nT
>CHAN_2    :HY
>AZIM_2    :  90.
>UNITS_2   :nT
>CHAN_3    :HZ
>AZIM_3    :  0.
>UNITS_3   :nT
>CHAN_4    :EX
>AZIM_4    :  0.
>UNITS_4   :mV/km
>CHAN_5    :EY
>AZIM_5    :  90.
>UNITS_5   :mV/km
>STARTTIME :2003-11-08 16:30:00
>DELTA_T   :  5.
>MIS_DATA  :  1.00000003E+32
>INFO_END  :
"""


def _synthetic_ts(n=2**14, dt=5.0, seed=7):
    rng = np.random.default_rng(seed)
    z0 = np.array([[2.0, 12.0], [-10.0, -1.5]])
    hx = rng.standard_normal(n)
    hy = rng.standard_normal(n)
    return TSData(
        data={
            "HX": hx,
            "HY": hy,
            "HZ": 0.3 * hx - 0.2 * hy,
            "EX": z0[0, 0] * hx + z0[0, 1] * hy,
            "EY": z0[1, 0] * hx + z0[1, 1] * hy,
        },
        dt=dt,
        station="SYN01",
        lat=-30.5,
        lon=20.25,
        azim={"HX": 0.0, "HY": 90.0, "HZ": 0.0,
              "EX": 0.0, "EY": 90.0},
        dipole={"EX": 50.0, "EY": 60.0},
    ), z0


@pytest.fixture()
def lims_file(tmp_path):
    rng = np.random.default_rng(5)
    rows = rng.standard_normal((2048, 5))
    body = "\n".join(
        "  ".join(f"{v:.6f}" for v in row)
        for row in rows
    )
    p = tmp_path / "tst01.ts"
    p.write_text(
        _LIMS_HEADER + body + "\n", encoding="utf-8"
    )
    return p


# ----------------------------------------------------------------- tests
def test_transform_single_tsdata():
    ts, z0 = _synthetic_ts()
    col = TStoEDI(nfft=2048, per_decade=4).transform(ts)
    assert isinstance(col, EDICollection)
    assert len(col) == 1
    ed = col[0]
    assert (ed.station or "").upper() == "SYN01"
    assert ed.Z.n_freq > 3
    assert ed.has_tipper
    assert np.allclose(
        ed.Z.z,
        np.broadcast_to(z0[None], ed.Z.z.shape),
        rtol=1e-6, atol=1e-8,
    )


def test_station_suffix_and_override():
    ts, _ = _synthetic_ts()
    tr = TStoEDI(
        nfft=2048, per_decade=4, station_suffix="_TS"
    )
    ed = tr.transform(ts)[0]
    assert ed.station == "SYN01_TS"
    # explicit override wins over the suffix
    ed2 = tr.transform(ts, station_name="CUSTOM01")[0]
    assert ed2.station == "CUSTOM01"


def test_transform_lims_file_writes_edi(
    lims_file, tmp_path
):
    out_dir = tmp_path / "edis"
    col = TStoEDI(
        nfft=512, per_decade=4, verbose=1
    ).transform(lims_file, output_dir=out_dir)
    assert len(col) == 1
    ed = col[0]
    assert ed.station == "tst01"
    written = list(out_dir.glob("*.edi"))
    assert len(written) == 1
    assert written[0].name == "tst01_ts.edi"
    # round-trip the written file
    from pycsamt.seg.edi import EDIFile

    ed2 = EDIFile(str(written[0]))
    assert ed2.Z.n_freq == ed.Z.n_freq
    dm = ed2.get_section("definemeas")
    assert dm is not None and len(dm.emeas) == 2


def test_batch_directory_with_failure(
    lims_file, tmp_path
):
    # a junk record that no reader can process
    bad = lims_file.parent / "bad.ts"
    bad.write_text("hello\nworld\n", encoding="utf-8")

    tr = TStoEDI(nfft=512, per_decade=4)
    result = tr.transform_batch(lims_file.parent)
    assert isinstance(result, TransformResult)
    assert result.n_ok == 1
    assert result.n_fail == 1
    assert "bad.ts" in result.failures[0].source

    # strict mode raises on the first failure
    with pytest.raises(RuntimeError):
        TStoEDI(
            nfft=512, per_decade=4, skip_errors=False
        ).transform(lims_file.parent)


def test_resolve_sources_errors(tmp_path):
    tr = TStoEDI()
    with pytest.raises(TypeError):
        tr._resolve_sources(123)
    with pytest.raises(ValueError):
        tr._resolve_sources(tmp_path / "nope.ts")
    # empty directory -> ValueError from transform_batch
    empty = tmp_path / "empty"
    empty.mkdir()
    with pytest.raises(ValueError):
        tr.transform_batch(empty)


def test_mixed_source_list(lims_file):
    ts, _ = _synthetic_ts()
    col = TStoEDI(nfft=512, per_decade=4).transform(
        [ts, lims_file]
    )
    stations = sorted(
        (e.station or "") for e in col
    )
    assert stations == ["SYN01", "tst01"]
