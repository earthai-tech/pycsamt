# -*- coding: utf-8 -*-
from __future__ import annotations

import sys
import types
import importlib

import numpy as np
import pytest


def _install_stubs(monkeypatch):
    # Ensure parent packages exist
    for name in [
        "pycsamt.seg",
        "pycsamt.seg.edi",
        "pycsamt.seg.collection",
        "pycsamt.zonge",
        "pycsamt.zonge.avg",
        "pycsamt.jones",
        "pycsamt.jones.j",
        "pycsamt.jones.collection",
    ]:
        if name not in sys.modules:
            monkeypatch.setitem(sys.modules, name, types.ModuleType(name))

    # ---- seg stubs -------------------------------------------------------
    seg_edi = sys.modules["pycsamt.seg.edi"]

    class ZBlock:
        def __init__(self):
            self._freq = None
            self._z = None
            self._z_err = None

        def compute_resistivity_phase(self):
            return None

        def set_res_phase(self, rho, phi, freq):
            self._freq = freq
            self._z = np.zeros((len(rho), 2, 2), complex)

    class TipBlock:
        def __init__(self):
            self._freq = None
            self._tipper = None
            self._tipper_err = None

        def compute_amp_phase(self):
            return None

        def compute_mag_direction(self):
            return None

    class EDIFile:
        def __init__(self, verbose: int = 0):
            self.Z = ZBlock()
            self.Tip = TipBlock()
            self.sections = {}
            self.station = None

        def add_section(self, name, obj):
            self.sections[name] = obj

    seg_edi.EDIFile = EDIFile

    seg_coll = sys.modules["pycsamt.seg.collection"]

    class EDICollection:
        def __init__(self, items, verbose: int = 0):
            self.items = list(items)

        def __iter__(self):
            return iter(self.items)

        def __len__(self):
            return len(self.items)

    seg_coll.EDICollection = EDICollection

    # ---- zonge AVG stub --------------------------------------------------
    zonge_avg = sys.modules["pycsamt.zonge.avg"]

    class AVG:
        def __init__(self, z=None, z_err=None, f=None, sites=None):
            self._z = z
            self._z_err = z_err
            self._f = f
            self._sites = sites or [1]

        @classmethod
        def from_file(cls, path):
            return cls()

        def to_tensor(self, var, station=None, sort_freq=True, align="union"):
            if var == "z":
                return self._z, self._f, self._sites
            if var == "z_err":
                return self._z_err, None, None
            if var == "rho":
                return None, self._f, self._sites
            if var == "phase":
                return None, None, None
            raise KeyError(var)

    zonge_avg.AVG = AVG

    # ---- Jones stubs -----------------------------------------------------
    j_mod = sys.modules["pycsamt.jones.j"]

    class JFile:
        def __init__(self, z=None, f=None, rho=None, phi=None, tip=None):
            class Z:
                def __init__(self, arr, freq):
                    self.z = arr
                    self.freq = freq

            class RP:
                def __init__(self, r, p):
                    self.rho = r
                    self.phase = p

            class T:
                def __init__(self, t):
                    self.tipper = t

            self.Z = Z(z, f) if z is not None else None
            self.ResPhase = RP(rho, phi) if rho is not None else None
            self.Tipper = T(tip) if tip is not None else None

        @classmethod
        def from_file(cls, path):
            return cls()

    j_mod.JFile = JFile

    jc_mod = sys.modules["pycsamt.jones.collection"]

    class JCollection(list):
        pass

    jc_mod.JCollection = JCollection

    # Import transformers now that stubs are in place
    from pycsamt.core import transformers as tr

    importlib.reload(tr)
    return tr, seg_coll, seg_edi, zonge_avg, j_mod, jc_mod


@pytest.fixture()
def tr_env(monkeypatch):
    return _install_stubs(monkeypatch)


def test_public_api_all(tr_env):
    tr, *_ = tr_env
    assert set(tr.__all__) == {"AVGtoEDI", "JtoEDI"}


def test_avg_to_edi_builds_collection(tr_env):
    tr, seg_coll, seg_edi, zonge_avg, *_ = tr_env
    n_site, n_freq = 2, 5
    f = np.logspace(1, 3, n_freq)
    z = np.ones((n_site, n_freq, 2, 2), complex)
    ze = np.full((n_site, n_freq, 2, 2), 0.1)
    avg = zonge_avg.AVG(z=z, z_err=ze, f=f, sites=[101, 102])

    out = tr.AVGtoEDI().transform(avg)
    assert isinstance(out, seg_coll.EDICollection)
    assert len(out) == 2
    for ed in out:
        assert isinstance(ed, seg_edi.EDIFile)
        assert ed.Z._freq is not None
        assert ed.Z._z is not None
        assert ed.Z._z.shape == (n_freq, 2, 2)


def test_j_to_edi_single_and_collection(tr_env):
    tr, seg_coll, seg_edi, zonge_avg, j_mod, jc_mod = tr_env

    # single file
    f = np.array([1.0, 2.0, 3.0])
    z = np.ones((f.size, 2, 2), complex)
    jf = j_mod.JFile(z=z, f=f)
    ed = tr.JtoEDI().transform(jf)
    assert isinstance(ed, seg_edi.EDIFile)
    assert ed.Z._freq is not None and ed.Z._z is not None

    # collection
    coll = jc_mod.JCollection([
        j_mod.JFile(z=z, f=f),
        j_mod.JFile(z=z, f=f),
    ])
    out = tr.JtoEDI().transform(coll)
    assert isinstance(out, seg_coll.EDICollection)
    assert len(out) == 2
