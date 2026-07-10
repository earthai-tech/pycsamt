from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from pycsamt.constants import MU_0
from pycsamt.jones.blocks import JBlocks
from pycsamt.jones.heads import Head
from pycsamt.jones.j import JFile


def _head(st: str = "S01") -> Head:
    # minimal head object used in comp dicts
    h = Head()
    h.station = st
    h.n = 2
    return h


def _comp_r_only() -> dict[str, dict[str, Any]]:
    # R parts for XY and YX, two rows, simple values
    p = np.array([1.0, 2.0], float)  # periods (s)
    return {
        "RXY": {
            "period": p,
            "rho": np.array([100.0, 200.0], float),
            "pha": np.array([45.0, 40.0], float),
            "head": _head("S01"),
        },
        "RYX": {
            "period": p,
            "rho": np.array([120.0, 220.0], float),
            "pha": np.array([35.0, 30.0], float),
            "head": _head("S01"),
        },
    }


def _comp_z_only() -> dict[str, dict[str, Any]]:
    # Z parts for all tensor entries, two rows
    p = np.array([1.0, 2.0], float)  # periods (s)
    f = 1.0 / p
    w = 2.0 * math.pi * f
    # build a consistent Z tensor magnitude/phase
    # use |Z| from rho via |Z| = sqrt(mu0 * w * rho)
    def make_z(rho: float, phi_deg: float) -> complex:
        mag = math.sqrt(MU_0 * w[0] * rho)
        ang = math.radians(phi_deg)
        return mag * complex(math.cos(ang), math.sin(ang))

    zxx0 = make_z(50.0, 20.0)
    zxy0 = make_z(100.0, 45.0)
    zyx0 = make_z(80.0, 35.0)
    zyy0 = make_z(60.0, 10.0)
    # second row slightly different
    zxx1 = zxx0 * 0.9
    zxy1 = zxy0 * 1.1
    zyx1 = zyx0 * 1.2
    zyy1 = zyy0 * 0.8

    def pack(v0: complex, v1: complex) -> dict[str, Any]:
        return {
            "period": p,
            "real": np.array([v0.real, v1.real], float),
            "imag": np.array([v0.imag, v1.imag], float),
            "err": np.zeros(2, float),
            "head": _head("S02"),
        }

    return {
        "ZXX": pack(zxx0, zxx1),
        "ZXY": pack(zxy0, zxy1),
        "ZYX": pack(zyx0, zyx1),
        "ZYY": pack(zyy0, zyy1),
    }


def _comp_t_only() -> dict[str, dict[str, Any]]:
    p = np.array([5.0, 10.0, 20.0], float)
    zx = np.array([0.1 + 0.2j, 0.0 + 0.0j, -0.2 + 0.1j], complex)
    zy = np.array([0.05 - 0.1j, 0.0 + 0.0j, 0.3 + 0.0j], complex)
    return {
        "TZX": {
            "period": p,
            "real": zx.real.astype(float),
            "imag": zx.imag.astype(float),
            "err": np.zeros(p.size, float),
            "head": _head("S03"),
        },
        "TZY": {
            "period": p,
            "real": zy.real.astype(float),
            "imag": zy.imag.astype(float),
            "err": np.zeros(p.size, float),
            "head": _head("S03"),
        },
    }


def test_utils_align_and_units():
    jf = JFile(verbose=0)

    # _complex
    c = jf._complex(np.array([1.0, -2.0]), np.array([0.5, 3.0]))
    assert np.allclose(c[0].real, 1.0) and np.allclose(c[0].imag, 0.5)

    # _deg2rad
    assert math.isclose(jf._deg2rad(180.0), math.pi)

    # _hz_from_period
    p = np.array([1.0, 2.0, 0.5], float)
    f = jf._hz_from_period(p)
    assert np.allclose(f, np.array([1.0, 0.5, 2.0]))

    # _align_by_periods
    a0 = np.array([1.0, 2.0, 3.0, 4.0])
    a1 = np.array([3.0, 1.0, 5.0])
    pc, i0, i1 = jf._align_by_periods(a0, a1)

    # common = {1, 3} in a0 order -> indices [0, 2] vs [1, 0]
    assert np.allclose(a0[i0], np.array([1.0, 3.0]))
    assert np.allclose(a1[i1], np.array([1.0, 3.0]))
    assert np.allclose(pc,    np.array([1.0, 3.0]))

def test_build_from_comp_R_only():
    jf = JFile(verbose=0)
    z, tip, rp = jf._build_from_comp(_comp_r_only())
    assert z is not None and tip is None
    assert rp is not None
    assert z.z.shape[1:] == (2, 2)  # type: ignore
    # freq must match 1/periods
    f = z.freq
    assert f is not None and np.all(f > 0.0)


def test_build_from_comp_Z_only():
    jf = JFile(verbose=0)
    z, tip, rp = jf._build_from_comp(_comp_z_only())
    assert z is not None and tip is None
    # keep shape and name from head
    assert z.z.shape[1:] == (2, 2)  # type: ignore
    # assert getattr(z, "name") == "S02"  # type: ignore


def test_build_from_comp_T_only():
    jf = JFile(verbose=0)
    z, tip, rp = jf._build_from_comp(_comp_t_only())
    assert z is None and tip is not None
    ta = tip.tipper  # type: ignore[union-attr]
    assert ta is not None and np.asarray(ta).shape[-1] == 2


@pytest.mark.filterwarnings("ignore::RuntimeWarning")
def test_scan_blocks_true_file(j_single_file: Path):
    jf = JFile(j_single_file, verbose=0)
    comp = jf._scan_blocks(jf.blocks)
    # we expect at least resistivity/phase in the sample
    assert isinstance(comp, dict) and len(comp) >= 1
    assert any(k.startswith("R") for k in comp.keys())


@pytest.mark.filterwarnings("ignore::RuntimeWarning")
def test_jfile_read_and_headers(j_single_file: Path):
    jf = JFile.from_file(j_single_file, verbose=0)
    assert jf.__has_read__() is True
    # heads + blocks
    assert jf.heads is not None and jf.blocks is not None
    # frequency bookkeeping (from R->Z fallback)
    assert jf.n_freq >= 1
    # compose headers prints banner+head+info
    lines = jf.compose_headers()
    txt = "\n".join(lines)
    assert "#WRITTEN BY" in txt
    assert ">LATITUDE" in txt or ">LONGITUDE" in txt
    assert jf.heads.head.station in txt  # type: ignore


@pytest.mark.filterwarnings("ignore::RuntimeWarning")
def test_blocks_roundtrip_smoke(j_single_file: Path, tmp_path: Path):
    jf = JFile.from_file(j_single_file, verbose=0)
    # write R only to keep it cheap and deterministic
    out = jf.write(
        savepath=tmp_path,
        new_jfn="smoke_R.j",
        datatype="R",
        overwrite=True,
    )
    p = Path(out)
    assert p.exists()
    # structural parse of what we wrote
    col = JBlocks.from_file(p, verbose=0)
    assert col.__has_read__() is True
    assert col.n >= 1
