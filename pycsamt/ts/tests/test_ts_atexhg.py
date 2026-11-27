# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for the ATEXHG reader and sensor-response correction."""

from __future__ import annotations

import struct
from pathlib import Path

import numpy as np
import pytest

from pycsamt.exceptions import FileHandlingError
from pycsamt.ts import (
    SensorResponse,
    apply_sensor_response,
    read_atexhg_header,
    read_sensor_response,
    read_ts,
    sniff_format,
)
from pycsamt.ts._base import TSData

_CH = [  # kind, index, dipole, gain code, azimuth
    (1, 1, 30.0, 5, 138.0),
    (1, 2, 30.0, 25, 228.0),
    (2, 1, 0.0, -2, 138.0),
    (2, 2, 0.0, 1, 228.0),
]


def _write_atexhg(path: Path, data: np.ndarray, fs: float = 4000.0,
                  channels=_CH, header_size: int = 2048,
                  magic: bytes = b"<ATEXHG>", bits: int = 32,
                  total: int | None = None) -> Path:
    """Write a synthetic file in the inferred ATEXHG layout."""
    n_ch = len(channels)
    head = bytearray(header_size)
    head[0:8] = magic
    struct.pack_into("<HH", head, 0x08, n_ch, bits)
    struct.pack_into("<d", head, 0x10, fs)
    tot = data.size if total is None else total
    struct.pack_into("<III", head, 0x20, 599, 999750, tot)
    for i, (kind, idx, dip, code, azm) in enumerate(channels):
        o = 0x58 + 54 * i
        head[o:o + 2] = b"\xff\xff"
        struct.pack_into("<HHf", head, o + 2, kind, idx, dip)
        head[o + 10:o + 16] = f"{i + 1:<4d}".encode() + b"\0\0"
        struct.pack_into("<i", head, o + 30, code)
        struct.pack_into("<f", head, o + 38, azm)
    path.write_bytes(bytes(head)
                     + data.astype("<f4").tobytes())
    return path


@pytest.fixture()
def atts(tmp_path):
    rng = np.random.default_rng(0)
    data = rng.normal(size=(1000, 4)).astype("<f4")
    return _write_atexhg(tmp_path / "S-30.atts", data), data


def test_sniff_and_read(atts):
    path, data = atts
    assert sniff_format(path) == "atexhg"
    ts = read_ts(path)
    assert ts.channels() == ["EX", "EY", "HX", "HY"]
    assert ts.dt == pytest.approx(1 / 4000.0)
    assert ts.n_samples == 1000
    assert ts.station == "S-30"
    np.testing.assert_allclose(ts.get("EY"), data[:, 1])
    assert ts.azim == {"EX": 138.0, "EY": 228.0, "HX": 138.0,
                       "HY": 228.0}
    assert ts.dipole == {"EX": 30.0, "EY": 30.0}
    assert ts.meta["gain_code"]["HX"] == -2
    assert ts.meta["calibrated"] is False
    assert ts.units["HX"] == "mV"


def test_header_offset_from_sample_count(tmp_path):
    data = np.ones((10, 4), dtype="<f4")
    p = _write_atexhg(tmp_path / "x.atts", data, header_size=3000)
    h = read_atexhg_header(p)
    assert h.data_offset == 3000 and h.n_samples == 10
    assert h.raw == {"word_0x20": 599, "word_0x24": 999750}


@pytest.mark.parametrize(
    "kw, msg",
    [({"magic": b"<NOTATS>"}, "not an <ATEXHG>"),
     ({"bits": 16}, "float32 only"),
     ({"total": 10 ** 9}, "more than the file"),
     ({"total": 4001}, "not a multiple"),
     ({"channels": [_CH[0], _CH[0], _CH[2], _CH[3]]}, "duplicate")],
)
def test_header_validation(tmp_path, kw, msg):
    data = np.ones((1000, 4), dtype="<f4")
    p = _write_atexhg(tmp_path / "bad.atts", data, **kw)
    with pytest.raises(FileHandlingError, match=msg):
        read_atexhg_header(p)


def test_text_files_still_sniffed(tmp_path):
    p = tmp_path / "t.txt"
    p.write_text("1 2 3\n4 5 6\n")
    assert sniff_format(p) == "ascii"


# ------------------------------------------------ sensor response
def _coil():
    f = np.logspace(-2, 4, 61)
    fc = 5.0   # first-order high-pass coil, corner fc
    s = 0.2 * (1j * f / fc) / (1 + 1j * f / fc)
    return SensorResponse(f, np.abs(s), np.degrees(np.angle(s)))


def _forward(b, dt, resp):
    spec = np.fft.rfft(b)
    f = np.fft.rfftfreq(b.size, dt)
    s = np.zeros(f.size, complex)
    s[1:] = resp.at(f[1:])
    return np.fft.irfft(spec * s, n=b.size)


def test_phase_convention_cosine():
    """A coil output A cos(wt + phi) for a unit cosine input."""
    resp = SensorResponse([0.1, 100.0], [3.0, 3.0], [30.0, 30.0])
    dt, n, f0 = 1 / 128.0, 1024, 4.0
    t = np.arange(n) * dt
    v = _forward(np.cos(2 * np.pi * f0 * t), dt, resp)
    expect = 3.0 * np.cos(2 * np.pi * f0 * t + np.radians(30.0))
    np.testing.assert_allclose(v, expect, atol=1e-10)


def test_deconvolution_recovers_field():
    rng = np.random.default_rng(1)
    dt, n = 1 / 256.0, 8192
    f = np.fft.rfftfreq(n, dt)
    spec = np.zeros(f.size, complex)
    band = (f > 1.0) & (f < 100.0)        # inside coil passband
    spec[band] = rng.normal(size=band.sum()) + 1j * rng.normal(
        size=band.sum())
    b = np.fft.irfft(spec, n=n)
    resp = _coil()
    ts = TSData({"HX": _forward(b, dt, resp)}, dt=dt)
    out = apply_sensor_response(ts, {"HX": resp},
                                water_level=1e-6)
    np.testing.assert_allclose(out.get("HX"), b,
                               atol=1e-8 * np.abs(b).max())
    assert out.units["HX"] == "nT"
    assert out.meta["calibrated"] is True
    # input untouched (inplace=False)
    assert ts.units.get("HX") is None


def test_water_level_bounds_low_frequency_gain():
    dt, n = 1 / 256.0, 4096
    ts = TSData({"HX": np.random.default_rng(2).normal(size=n)},
                dt=dt)
    resp = _coil()
    out = apply_sensor_response(ts, {"HX": resp}, water_level=0.1)
    gain = np.abs(np.fft.rfft(out.get("HX"))[1:]) / np.abs(
        np.fft.rfft(ts.get("HX"))[1:])
    # floor is relative to max |S| inside the recorded band
    smax = np.abs(resp.at(np.fft.rfftfreq(n, dt)[1:])).max()
    assert gain.max() <= 1 / (0.1 * smax) * (1 + 1e-9)
    assert gain.max() == pytest.approx(1 / (0.1 * smax), rel=1e-6)


def test_response_validation_and_reader(tmp_path):
    with pytest.raises(ValueError, match="same length"):
        SensorResponse([1.0], [1.0], [0.0])
    with pytest.raises(ValueError, match="> 0"):
        SensorResponse([0.0, 1.0], [1.0, 1.0], [0.0, 0.0])
    resp = SensorResponse([1.0, 10.0], [1.0, 1.0], [0.0, 0.0])
    with pytest.raises(ValueError, match="outside"):
        resp.at(np.array([100.0]), extrapolate="raise")
    p = tmp_path / "coil.cal"
    p.write_text("# coil 123\nfreq amp phase\n"
                 "0.1, 0.01, 89\n1, 0.1, 80\n10, 0.5, 30\n")
    r = read_sensor_response(p, delimiter=",")
    np.testing.assert_allclose(r.freq, [0.1, 1.0, 10.0])
    assert r.name == "coil"


def test_apply_errors():
    ts = TSData({"HX": np.array([1.0, np.nan, 2.0])}, dt=0.1)
    resp = SensorResponse([1.0, 10.0], [1.0, 1.0], [0.0, 0.0])
    with pytest.raises(ValueError, match="NaN"):
        apply_sensor_response(ts, {"HX": resp})
    with pytest.raises(ValueError, match="not in time series"):
        apply_sensor_response(ts, {"HY": resp})
    with pytest.raises(ValueError, match="water_level"):
        apply_sensor_response(ts, {"HX": resp}, water_level=0)


def test_apply_gain(tmp_path):
    from pycsamt.ts.atexhg import gain_from_code

    data = np.ones((20, 4), dtype="<f4")
    p = _write_atexhg(tmp_path / "g.atts", data)
    ts = read_ts(p, apply_gain=True)
    # codes from _CH: EX 5, EY 25, HX -2 (x1/2), HY 1
    np.testing.assert_allclose(ts.get("EX"), 1 / 5)
    np.testing.assert_allclose(ts.get("EY"), 1 / 25)
    np.testing.assert_allclose(ts.get("HX"), 2.0)
    np.testing.assert_allclose(ts.get("HY"), 1.0)
    assert ts.meta["gain_applied"] is True
    assert read_ts(p).meta["gain_applied"] is False
    with pytest.raises(ValueError):
        gain_from_code(0)
