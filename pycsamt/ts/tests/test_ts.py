# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for :mod:`pycsamt.ts` (readers, spectra, ts→EDI)."""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.ts import (
    TSData,
    read_ts,
    sniff_format,
    target_frequencies,
    ts_to_edi,
    ts_to_spectra,
    ts_to_z,
)

# --------------------------------------------------------------- fixtures
_LIMS_HEADER = """\
# time series file from mp2ts
# date: Sun Nov  6 13:24:57 2016
#
# Ex line length (m):      50.00
# Ey line length (m):      60.00
#
>INFO_START:
>STATION   :tst01
>INSTRUMENT: 66
>LATITUDE  : -30.5
>LONGITUDE :  20.25
>ELEVATION :  100.
>COORD_SYS :MAGNETIC NORTH
>DECLIN    :  0.
>NCHAN     : 5
>CHAN_1    :HX
>AZIM_1    :  0.
>UNITS_1   :nT
>GAIN_1    :  1.
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
>ENDTIME   :2003-11-08 16:31:00
>T_UNITS   :s
>DELTA_T   :  5.
>MIS_DATA  :  1.00000003E+32
>INFO_END  :
"""


@pytest.fixture()
def lims_file(tmp_path):
    rows = []
    for k in range(12):
        if k == 6:  # one missing row
            rows.append("  ".join(["1.00000003E+32"] * 5))
            continue
        rows.append("  ".join(f"{0.1 * (k + j):.6f}" for j in range(5)))
    p = tmp_path / "tst01.ts"
    p.write_text(
        _LIMS_HEADER + "\n".join(rows) + "\n",
        encoding="utf-8",
    )
    return p


@pytest.fixture()
def emslab_file(tmp_path):
    def hour_block(stamp, base):
        lines = [f"TEST01{stamp}"]
        vals = []
        for s in range(180):  # 180 samples x 5 chan
            for c in range(5):
                v = base + 10 * c + (s % 3)
                vals.append(v)
        # one explicit missing value
        vals[7] = -9999
        for i in range(0, 900, 15):
            lines.append("".join(f"{v:5d}" for v in vals[i : i + 15]))
        return lines

    lines = hour_block("8507180100", 10)
    # hour 02 missing entirely -> must be NaN-padded
    lines += hour_block("8507180300", 40)
    p = tmp_path / "test01.asc"
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return p


def _synthetic_ts(n=2**15, dt=5.0, z0=None, tipper=(0.3, -0.2), seed=7):
    rng = np.random.default_rng(seed)
    if z0 is None:  # real (instantaneous) tensor
        z0 = np.array([[2.0, 12.0], [-10.0, -1.5]])
    hx = rng.standard_normal(n)
    hy = rng.standard_normal(n)
    ex = z0[0, 0] * hx + z0[0, 1] * hy
    ey = z0[1, 0] * hx + z0[1, 1] * hy
    hz = tipper[0] * hx + tipper[1] * hy
    ts = TSData(
        data={
            "HX": hx,
            "HY": hy,
            "HZ": hz,
            "EX": ex,
            "EY": ey,
        },
        dt=dt,
        station="SYN01",
        lat=-30.5,
        lon=20.25,
        elev=100.0,
        azim={"HX": 0.0, "HY": 90.0, "HZ": 0.0, "EX": 0.0, "EY": 90.0},
        dipole={"EX": 50.0, "EY": 60.0},
    )
    return ts, z0


# ---------------------------------------------------------------- readers
def test_sniff_and_read_lims(lims_file):
    assert sniff_format(lims_file) == "lims"
    ts = read_ts(lims_file)
    assert ts.station == "tst01"
    assert ts.channels() == ["HX", "HY", "HZ", "EX", "EY"]
    assert ts.dt == 5.0
    assert ts.lat == pytest.approx(-30.5)
    assert ts.lon == pytest.approx(20.25)
    assert ts.dipole == {"EX": 50.0, "EY": 60.0}
    assert ts.azim["HY"] == 90.0
    assert ts.units["EX"].lower() == "mv/km"
    # sentinel row -> NaN
    assert ts.n_samples == 12
    for cid in ts.channels():
        assert np.isnan(ts.get(cid)[6])
        assert np.isfinite(np.delete(ts.get(cid), 6)).all()
    # data content preserved
    assert ts.get("HX")[0] == pytest.approx(0.0)
    assert ts.get("EY")[1] == pytest.approx(0.5)


def test_read_emslab(emslab_file):
    assert sniff_format(emslab_file) == "emslab"
    ts = read_ts(emslab_file)
    assert ts.station == "TEST01"
    assert ts.dt == 20.0
    # 3 hours span (hour 02 missing) -> 540 samples
    assert ts.n_samples == 3 * 180
    hx = ts.get("HX")
    # hour-1 first sample: base 10 + (0 % 3) -> 10 * 0.1
    assert hx[0] == pytest.approx(1.0)
    # explicit -9999 (raw index 7 = sample 1, chan HZ) -> NaN
    assert np.isnan(ts.get("HZ")[1])
    # padded missing hour -> NaN everywhere
    assert np.isnan(hx[180:360]).all()
    # hour-3 resumes with base 40
    assert hx[360] == pytest.approx(4.0)
    # units scaled by 0.1, start one dt after hour mark
    assert ts.start == "1985-07-18 01:00:20"


def test_read_ascii(tmp_path):
    p = tmp_path / "plain.txt"
    rows = np.arange(20.0).reshape(10, 2)
    p.write_text(
        "\n".join(f"{a} {b}" for a, b in rows),
        encoding="utf-8",
    )
    ts = read_ts(p, format="ascii", dt=1.0, chan=("HX", "HY"))
    assert ts.channels() == ["HX", "HY"]
    assert ts.get("HY")[3] == pytest.approx(7.0)


# ----------------------------------------------------------------- grids
def test_target_frequencies_descending():
    f = target_frequencies(5.0, 4096, per_decade=6)
    assert (np.diff(f) < 0).all()
    assert f.max() <= 0.05 + 1e-12  # fs/4
    assert f.min() >= 4.0 / (4096 * 5.0) - 1e-12


# ------------------------------------------------------- spectra & Z
def test_synthetic_recovery_exact():
    ts, z0 = _synthetic_ts()
    sp = ts_to_spectra(ts, nfft=2048, per_decade=4)
    assert sp.n_freq > 3
    assert sp.chan_ids[:5] == ["HX", "HY", "HZ", "EX", "EY"]
    zhat, tip = sp.to_Z(estimate_error=True)
    # instantaneous (real) tensor recovers to precision
    assert np.abs(zhat.z - z0[None]).max() < 1e-8
    # tipper recovered
    assert tip is not None
    assert np.abs(tip.tipper[:, 0, 0] - 0.3).max() < 1e-8
    assert np.abs(tip.tipper[:, 0, 1] + 0.2).max() < 1e-8
    # DoF-based errors attached and finite
    assert zhat.z_err is not None
    assert np.isfinite(zhat.z_err).all()


def test_noise_robustness_and_gaps():
    ts, z0 = _synthetic_ts(seed=11)
    rng = np.random.default_rng(3)
    # contaminate E with noise and punch a long gap
    for cid in ("EX", "EY"):
        x = ts.get(cid)
        x += 0.05 * np.std(x) * rng.standard_normal(x.size)
        x[5000:5600] = np.nan
    zhat, tip, sp = ts_to_z(
        ts,
        nfft=2048,
        per_decade=4,
        estimate_error=False,
    )
    # 5% E-noise, gap-dropped segments: still within a few %
    rel = np.abs(zhat.z - z0[None]) / np.abs(z0[None])
    assert np.median(rel) < 0.05


def test_remote_reference_estimator():
    ts, z0 = _synthetic_ts(seed=23)
    # noise-free remote copy of the magnetic channels
    remote = TSData(
        data={"HX": ts.get("HX").copy(), "HY": ts.get("HY").copy()},
        dt=ts.dt,
        station="REM01",
    )
    zrr, tip, sp = ts_to_z(
        ts,
        estimator="rr",
        remote=remote,
        nfft=2048,
        per_decade=4,
        estimate_error=False,
    )
    assert "RHX" in sp.chan_ids and "RHY" in sp.chan_ids
    assert np.abs(zrr.z - z0[None]).max() < 1e-8


# ---------------------------------------------------------------- ts→EDI
def test_ts_to_edi_roundtrip(tmp_path):
    ts, z0 = _synthetic_ts()
    out = ts_to_edi(
        ts,
        out="syn01.edi",
        savepath=tmp_path,
        nfft=2048,
        per_decade=4,
        include_spectra=True,
        include_tseries=True,
        tseries_max_samples=64,
    )
    from pycsamt.seg.edi import EDIFile

    ed = EDIFile(out)
    assert (ed.station or "").upper() == "SYN01"
    assert ed.Z.n_freq > 3
    assert ed.has_tipper
    # impedance round-trips through the EDI text format
    assert np.allclose(
        ed.Z.z,
        np.broadcast_to(z0[None], ed.Z.z.shape),
        rtol=1e-4,
        atol=1e-4,
    )
    # geometry made it into DEFINEMEAS
    dm = ed.get_section("definemeas")
    assert dm is not None
    exm = [m for m in dm.emeas if str(m.chtype).upper() == "EX"]
    assert exm and float(exm[0].x2) == pytest.approx(25.0)
    # optional sections present
    assert ed.get_section("spectra") is not None
    ts_sec = ed.get_section("timeseries")
    assert ts_sec is not None
    assert set(ts_sec.channels()) == {"HX", "HY", "HZ", "EX", "EY"}
    assert ts_sec.get("HX").size == 64


def test_ts_to_edi_from_lims_file(lims_file, tmp_path):
    # too short for spectra -> readers still work end to end
    ts = read_ts(lims_file)
    with pytest.raises(Exception):
        # 12 samples cannot support nfft=256 segmentation
        ts_to_edi(ts, savepath=tmp_path, nfft=256)
