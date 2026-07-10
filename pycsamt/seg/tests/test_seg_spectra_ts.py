
from __future__ import annotations

import numpy as np
import pytest

from pycsamt.seg.spectra import (
    Spectra,
    SpectraIO,
    SpectraSECT,
    _SpectraBlock,
)
from pycsamt.seg.time_series import (
    TSIO,
    TimeSeries,
    TSect,
    _TSBlock,
)

# ----------------------------- helpers

def _hermitian_example() -> np.ndarray:
    # 3x3 Hermitian with easy-to-check entries
    H = np.zeros((3, 3), complex)
    # diag purely real
    H[0, 0] = 1.1
    H[1, 1] = 2.2
    H[2, 2] = 3.3
    # off-diagonals:
    # (0,1) = 0.5 + 0.7j
    # (0,2) = -0.2 + 0.1j
    # (1,2) = 1.0 - 0.4j
    H[0, 1] = 0.5 + 0.7j
    H[1, 0] = 0.5 - 0.7j
    H[0, 2] = -0.2 + 0.1j
    H[2, 0] = -0.2 - 0.1j
    H[1, 2] = 1.0 - 0.4j
    H[2, 1] = 1.0 + 0.4j
    return H


def _pack_like_spectra(H: np.ndarray) -> np.ndarray:
    # mirror Spectra._pack without importing it
    n = H.shape[0]
    M = np.zeros((n, n), float)
    for i in range(n):
        M[i, i] = float(H[i, i].real)
        for j in range(i + 1, n):
            z = H[i, j]
            M[j, i] = float(z.real)
            M[i, j] = float(z.imag)
    return M.ravel()


# ----------------------------- Spectra

def test_spectra_from_io_roundtrip() -> None:
    # make header
    sect = SpectraSECT(
        sectid="SPECA",
        nchan=3,
        nfreq=2,
    )
    sect.meas_ids = ["HX", "HY", "HZ"]

    # build two blocks (ascending freq -> will reverse)
    H1 = _hermitian_example()
    H2 = _hermitian_example() * (1.0 + 0.0j)
    H2[0, 0] = 9.9  # small tweak

    io = SpectraIO()
    b0 = _SpectraBlock()
    b0.freq = 0.100
    b0.rotspec = 1
    b0.bw = 0.01
    b0.avgt = 32.0
    b0.options["avgf"] = 0.095
    b0.options["segnum"] = 4
    b0.options["band"] = "lf"
    b0.nvals_hint = 9
    b0.values = _pack_like_spectra(H1).tolist()

    b1 = _SpectraBlock()
    b1.freq = 1.000
    b1.rotspec = 0
    b1.bw = 0.10
    b1.avgt = 16.0
    b1.options["avgf"] = 0.98
    b1.options["segnum"] = 7
    b1.options["band"] = "hf"
    b1.nvals_hint = 9
    b1.values = _pack_like_spectra(H2).tolist()

    io.blocks = [b0, b1]

    # construct Spectra
    sp = Spectra.from_io(sect, io, empty=1.0e32)

    # reversed to high->low: [1.0, 0.1]
    assert sp.n_chan == 3
    assert sp.n_freq == 2
    assert np.allclose(sp.freq, [1.0, 0.1])

    # first slice should match b1 (freq 1.0)
    assert np.allclose(sp.S[0], H2)
    # second slice should match b0 (freq 0.1)
    assert np.allclose(sp.S[1], H1)

    # metadata preserved in order
    assert sp.chan_ids == ["HX", "HY", "HZ"]
    assert sp.rotspec[0] == pytest.approx(0.0)
    assert sp.bw[0] == pytest.approx(0.10)
    assert sp.band[0] == "HF"
    assert sp.segnum[1] == 4

    # round-trip to IO
    sect2, io2 = sp.to_io()
    assert sect2.nchan == 3
    assert sect2.nfreq == 2
    assert sect2.sectid == "SPECA"
    assert sect2.meas_ids == ["HX", "HY", "HZ"]

    assert len(io2.blocks) == 2
    # compare packed values block 0
    v0 = np.array(io2.blocks[0].values, float)
    assert np.allclose(v0, _pack_like_spectra(sp.S[0]))


# ----------------------------- TimeSeries

def test_timeseries_from_io_and_to_io() -> None:
    # header with dt fallback
    sect = TSect(
        sectid="TSA",
        nchan=2,
        npts=0,
        dt=0.5,
    )

    # three blocks: HX twice, HY once
    io = TSIO()

    b1 = _TSBlock()
    b1.options["id"] = "HX"
    b1.options["npts"] = 3
    b1.options["dt"] = 0.1
    b1.nvals_hint = 3
    b1.values = [1.0, 2.0, 3.0]

    b2 = _TSBlock()
    b2.options["id"] = "HX"
    b2.options["npts"] = 2
    # no DT -> fallback to sect.dt (0.5)
    b2.nvals_hint = 2
    b2.values = [10.0, 20.0]

    b3 = _TSBlock()
    b3.options["id"] = "HY"
    b3.options["npts"] = 4
    b3.options["dt"] = 0.2
    b3.nvals_hint = 4
    b3.values = [0.0, -1.0, -2.0, -3.0]

    io.blocks = [b1, b2, b3]

    ts = TimeSeries.from_io(sect, io)

    # channels and concatenation
    assert ts.channels() == ["HX", "HY"]
    hx = ts.get("HX")
    hy = ts.get("HY")
    assert np.allclose(hx, [1.0, 2.0, 3.0, 10.0, 20.0])
    assert np.allclose(hy, [0.0, -1.0, -2.0, -3.0])

    # dt rules: HX final dt from second block -> 0.5
    assert ts.dt_map["HX"] == pytest.approx(0.5)
    # HY dt from its own block
    assert ts.dt_map["HY"] == pytest.approx(0.2)

    # time vector uses per-channel dt
    t_hx = ts.time("HX")
    assert np.allclose(t_hx[:4], [0.0, 0.5, 1.0, 1.5])

    # to_io should emit one block per channel
    sect2, io2 = ts.to_io()
    assert sect2.nchan == 2
    assert sect2.nmeas == 2
    assert set(sect2.meas_ids) == {"HX", "HY"}

    # HX block
    b_hx = next(b for b in io2 if b.options["id"] == "HX")
    assert b_hx.options["npts"] == hx.size
    assert b_hx.options["dt"] == pytest.approx(0.5)
    assert np.allclose(b_hx.values, hx)

    # HY block
    b_hy = next(b for b in io2 if b.options["id"] == "HY")
    assert b_hy.options["npts"] == hy.size
    assert b_hy.options["dt"] == pytest.approx(0.2)
    assert np.allclose(b_hy.values, hy)
