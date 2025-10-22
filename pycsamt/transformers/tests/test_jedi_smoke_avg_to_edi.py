# -*- coding: utf-8 -*-
# tiny smoke test for AVG -> EDI

from __future__ import annotations

import numpy as np
import pandas as pd

from pycsamt.zonge.avg import BaseAVG
from pycsamt.transformers.jedi import AVGtoEDI


class _Comp:
    def __init__(self, comps: list[str]) -> None:
        self.frame = pd.DataFrame({"comp": comps})


class _Info:
    def __init__(self, comps: list[str]) -> None:
        self.comp = _Comp(comps)


class _Topo:
    def __init__(self) -> None:
        # empty frame => leaves LAT/LON at defaults
        self.frame = pd.DataFrame()


class FakeAVG(BaseAVG):
    """
    Minimal AVG-like object exposing:
      - to_tensor(var=...) returning (arr, freq, station)
      - info.comp.frame with 'comp' column
      - topo.frame (optional)
    """

    def __init__(self) -> None:
        self.info = _Info(["EX", "HY"])  # Ex/Hy present
        self.topo = _Topo()

        # 3 frequencies; simple Z with non-zero Zxy
        self._freq = np.array([1.0, 2.0, 4.0], float)
        mag = np.array([1.0, 2.0, 4.0], float)
        pha = np.deg2rad([10.0, 20.0, 30.0])  # radians

        z = np.zeros((3, 2, 2), complex)
        z[:, 0, 1] = mag * (np.cos(pha) + 1j * np.sin(pha))
        self._z = z
        self._ze = np.zeros_like(z, float)
        self._station = np.array([1500.0], float)

    # matches AVG API used by the transformer
    def to_tensor(
        self,
        var: str = "z",
        *,
        station=None,
        agg: str | None = "mean",
        fill_value: float = np.nan,
        sort_freq: bool = True,
        align: str = "union",
    ):
        if var == "z":
            # shape (nfreq, 2, 2)
            return self._z, self._freq, self._station
        if var == "z_err":
            return self._ze, self._freq, self._station
        if var == "rho" or var == "phase":
            # not needed since Z is provided
            raise RuntimeError("rho/phase not used in this test")
        raise KeyError(var)


def _first_edi(coll):
    # EDICollection is iterable by item
    return next(iter(coll))


def _mtsect(ed):
    # normalize any spelling the writer kept
    for k in ("mtsect", "MTSECT", ">=MTSECT", ">=EMAPSECT"):
        sec = ed.get_section(k)
        if sec is not None:
            return sec
    raise AssertionError("No MT/EMAP section found.")


def _nfreq_from(sec) -> int:
    # accept dict-like or object with attr
    if hasattr(sec, "NFREQ"):
        return int(getattr(sec, "NFREQ"))
    if hasattr(sec, "nfreq"):
        return int(getattr(sec, "nfreq"))
    if isinstance(sec, dict):
        if "NFREQ" in sec:
            return int(sec["NFREQ"])
        if "nfreq" in sec:
            return int(sec["nfreq"])
    raise AssertionError("NFREQ missing in MTSECT.")


def test_avgtoedi_minimal_fixture():
    avg = FakeAVG()
    coll = AVGtoEDI().transform(avg)
    assert coll is not None

    ed = _first_edi(coll)

    # invariants: sections exist
    assert ed.get_section("head") is not None
    assert ed.get_section("info") is not None
    assert ed.get_section("definemeas") is not None
    sec = _mtsect(ed)
    assert sec is not None

    # NFREQ matches len(FREQ)
    nfreq = _nfreq_from(sec)
    freq = getattr(ed.Z, "freq", None)
    assert freq is not None
    assert int(nfreq) == int(len(freq))

    # Zxy is finite (not empty sentinel)
    z = getattr(ed.Z, "z", None)
    if z is None:
        z = getattr(ed.Z, "_z", None)
    assert z is not None
    zxy = z[:, 0, 1]
    assert np.isfinite(zxy.real).all()
    assert np.isfinite(zxy.imag).all()

    # If resistivity/phase available, they align with freq
    # (do not force present; just sanity when computed)
    rho = getattr(ed.Z, "res_xy", None)
    phi = getattr(ed.Z, "phase_xy", None)
    if rho is not None:
        assert len(rho) == len(freq)
        assert np.isfinite(rho).all()
    if phi is not None:
        assert len(phi) == len(freq)
        assert np.isfinite(phi).all()
