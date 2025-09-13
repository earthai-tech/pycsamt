# -*- coding: utf-8 -*-
from __future__ import annotations

import sys
import types
from types import SimpleNamespace
import importlib

import pytest


def _install_stubs(monkeypatch):
    mods = [
        "pycsamt.seg",
        "pycsamt.seg.edi",
        "pycsamt.seg.collection",
        "pycsamt.zonge",
        "pycsamt.zonge.avg",
        "pycsamt.jones",
        "pycsamt.jones.j",
    ]
    for name in mods:
        if name not in sys.modules:
            monkeypatch.setitem(
                sys.modules, name, types.ModuleType(name)
            )

    # seg stubs
    seg_edi = sys.modules["pycsamt.seg.edi"]

    class EDIFile:
        def __init__(self):
            self.Z = SimpleNamespace()

        @classmethod
        def from_file(cls, path):
            return cls()

    seg_edi.EDIFile = EDIFile

    seg_coll = sys.modules["pycsamt.seg.collection"]

    class EDICollection(list):
        def __init__(self, items, verbose: int = 0):
            super().__init__(items)

    seg_coll.EDICollection = EDICollection

    # Zonge AVG stub
    zavg = sys.modules["pycsamt.zonge.avg"]

    class AVG:
        def __init__(self):
            self.topo = None
            self._added = False

        def add_topography(self, t):
            self._added = True
            self.topo = t

        @classmethod
        def from_file(cls, path):
            return cls()

    zavg.AVG = AVG

    # Jones JFile stub
    jj = sys.modules["pycsamt.jones.j"]

    class JFile:
        @classmethod
        def from_file(cls, path):
            return cls()

    jj.JFile = JFile

    return seg_edi, seg_coll, zavg, jj


@pytest.fixture()
def env(monkeypatch):
    return _install_stubs(monkeypatch)


def _import_ctx():
    # Import facade after stubs are in place
    mod = importlib.import_module("pycsamt.session")
    return mod


def test_public_api_all(env):
    ctx = _import_ctx()
    assert set(ctx.__all__) == {
        "Session",
        "work_session",
        "Normalize",
        "normalize_session",
    }


def test_session_wraps_to_edi_and_transformers(
    env, tmp_path, monkeypatch
):
    ctx = _import_ctx()

    # Patch the symbol captured by _session wrapper
    smod = importlib.import_module("pycsamt._session")

    called = {"to_edi": 0, "x": 0, "j": 0}

    def fake_to_edi(src, *a, **k):  # noqa: ANN001
        called["to_edi"] += 1
        return {"ok": True}

    monkeypatch.setattr(smod, "to_edi", fake_to_edi)

    # Also patch transformer methods
    from pycsamt.core import transformers as tr

    def fake_x(self, *a, **k):  # noqa: ANN001
        called["x"] += 1
        return [1]

    def fake_j(self, *a, **k):  # noqa: ANN001
        called["j"] += 1
        return {"edi": 1}

    monkeypatch.setattr(tr.AVGtoEDI, "transform", fake_x, raising=False)
    monkeypatch.setattr(tr.JtoEDI, "transform", fake_j, raising=False)

    with ctx.work_session(tmp_path) as s:
        # call through patched base.to_edi wrapper
        from pycsamt.core import base as b
        b.to_edi("src")
        # transformer paths
        tr.AVGtoEDI().transform("a")
        tr.JtoEDI().transform("b")

    assert called["to_edi"] == 1
    assert called["x"] == 1
    assert called["j"] == 1

    rec1 = s.reg.find(tag="to_edi")
    rec2 = s.reg.find(tag="AVGtoEDI.transform")
    rec3 = s.reg.find(tag="JtoEDI.transform")
    assert len(rec1) == 1 and len(rec2) == 1 and len(rec3) == 1


def test_normalize_routes_avg_path_and_edi_wrap(
    env, tmp_path, monkeypatch
):
    ctx = _import_ctx()
    seg_edi = sys.modules["pycsamt.seg.edi"]
    seg_coll = sys.modules["pycsamt.seg.collection"]

    from pycsamt.core import transformers as tr

    def tx(self, a):  # noqa: ANN001
        return seg_coll.EDICollection([seg_edi.EDIFile()])

    monkeypatch.setattr(tr.AVGtoEDI, "transform", tx, raising=False)

    avg_path = tmp_path / "a.avg"
    avg_path.write_text("", encoding="utf-8")

    with ctx.normalize_session(tmp_path) as nx:
        out = nx.load(str(avg_path))

    from pycsamt.seg.collection import EDICollection

    assert isinstance(out, EDICollection)

    # EDI file path → wrapped to collection
    def fe(p):  # noqa: ANN001
        return seg_edi.EDIFile()

    monkeypatch.setattr(seg_edi.EDIFile, "from_file", classmethod(fe))
    edi_path = tmp_path / "a.edi"
    edi_path.write_text("", encoding="utf-8")

    with ctx.normalize_session(tmp_path) as nx:
        out2 = nx.load(str(edi_path))
    assert isinstance(out2, EDICollection)
    assert len(out2) == 1


def test_normalize_routes_j_path(env, tmp_path, monkeypatch):
    ctx = _import_ctx()
    seg_edi = sys.modules["pycsamt.seg.edi"]
    from pycsamt.core import transformers as tr

    def tj(self, j):  # noqa: ANN001
        return seg_edi.EDIFile()

    monkeypatch.setattr(tr.JtoEDI, "transform", tj, raising=False)

    j_path = tmp_path / "a.j"
    j_path.write_text("", encoding="utf-8")

    with ctx.normalize_session(tmp_path) as nx:
        out = nx.load(str(j_path))

    from pycsamt.seg.edi import EDIFile

    assert isinstance(out, EDIFile)


def test_normalize_with_avg_object_and_topo_injection(
    env, tmp_path, monkeypatch
):
    ctx = _import_ctx()
    seg_edi = sys.modules["pycsamt.seg.edi"]
    seg_coll = sys.modules["pycsamt.seg.collection"]
    zavg = sys.modules["pycsamt.zonge.avg"]

    from pycsamt.core import transformers as tr

    def tx(self, a):  # noqa: ANN001
        return seg_coll.EDICollection([seg_edi.EDIFile()])

    monkeypatch.setattr(tr.AVGtoEDI, "transform", tx, raising=False)

    avg = zavg.AVG()
    topo = SimpleNamespace(frame="df")

    with ctx.normalize_session(tmp_path, topo_src=topo) as nx:
        out = nx.load(avg)
    assert avg.topo is topo

    recs = nx.reg.find(tag="normalized")
    assert len(recs) >= 1


def test_normalize_passthrough_collection(env, tmp_path):
    ctx = _import_ctx()
    seg_edi = sys.modules["pycsamt.seg.edi"]
    seg_coll = sys.modules["pycsamt.seg.collection"]
    coll = seg_coll.EDICollection([seg_edi.EDIFile()])

    with ctx.normalize_session(tmp_path) as nx:
        out = nx.load(coll)
    assert out is coll
