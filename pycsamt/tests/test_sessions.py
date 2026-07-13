from __future__ import annotations

import importlib
import sys
import types
from types import SimpleNamespace

import pytest

from pycsamt.core import base as b
from pycsamt.transformers import jedi as tr


def _install_stubs(monkeypatch):
    mods = [
        "pycsamt.seg",
        "pycsamt.seg.edi",
        "pycsamt.seg.collection",
        "pycsamt.zonge",
        "pycsamt.zonge.avg",
        "pycsamt.jones",
        "pycsamt.jones.j",
        "pycsamt.jones.collection",
    ]
    for name in mods:
        if name not in sys.modules:
            m = types.ModuleType(name)
            if name in {"pycsamt.seg", "pycsamt.zonge", "pycsamt.jones"}:
                m.__path__ = []  # make it package-like
            monkeypatch.setitem(sys.modules, name, m)

    # Jones: JFile + JCollection stubs
    jj = sys.modules["pycsamt.jones.j"]

    class JFile:
        @classmethod
        def from_file(cls, path):
            return cls()

    monkeypatch.setattr(jj, "JFile", JFile, raising=False)

    jcoll = sys.modules["pycsamt.jones.collection"]

    class JCollection(list):
        pass

    monkeypatch.setattr(jcoll, "JCollection", JCollection, raising=False)

    # SEG: EDIFile + EDICollection stubs
    seg_edi = sys.modules["pycsamt.seg.edi"]

    class EDIFile:
        def __init__(self):
            self.Z = SimpleNamespace()

        @classmethod
        def from_file(cls, path):
            return cls()

    monkeypatch.setattr(seg_edi, "EDIFile", EDIFile, raising=False)

    seg_coll = sys.modules["pycsamt.seg.collection"]

    class EDICollection(list):
        def __init__(self, items, verbose: int = 0):
            super().__init__(items)

    monkeypatch.setattr(
        seg_coll, "EDICollection", EDICollection, raising=False
    )

    # Zonge: AVG stub
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

    monkeypatch.setattr(zavg, "AVG", AVG, raising=False)

    # pycsamt._session imports all of these names directly at import
    # time, so patching the provider modules alone has no effect when
    # the real modules were imported first (always true in a full-suite
    # run). Rebind the session module's own references to the stubs.
    smod = importlib.import_module("pycsamt._session")
    monkeypatch.setattr(smod, "JFile", JFile, raising=False)
    monkeypatch.setattr(smod, "JCollection", JCollection, raising=False)
    monkeypatch.setattr(smod, "EDIFile", EDIFile, raising=False)
    monkeypatch.setattr(smod, "EDICollection", EDICollection, raising=False)
    monkeypatch.setattr(smod, "AVG", AVG, raising=False)

    return seg_edi, seg_coll, zavg, jj


@pytest.fixture()
def env(monkeypatch):
    return _install_stubs(monkeypatch)


def _import_ctx():
    # Import after stubs so session binds stub classes
    return importlib.import_module("pycsamt.session")


def test_public_api_all(env):
    ctx = _import_ctx()
    assert set(ctx.__all__) == {
        "Session",
        "work_session",
        "Normalize",
        "normalize_session",
    }


def test_session_wraps_to_edi_and_transformers(env, tmp_path, monkeypatch):
    ctx = _import_ctx()
    smod = importlib.import_module("pycsamt._session")

    called = {"to_edi": 0, "x": 0, "j": 0}

    def fake_to_edi(src, *a, **k):  # noqa: ANN001
        called["to_edi"] += 1
        return {"ok": True}

    monkeypatch.setattr(smod, "to_edi", fake_to_edi)

    def fake_x(self, *a, **k):  # noqa: ANN001
        called["x"] += 1
        return [1]

    def fake_j(self, *a, **k):  # noqa: ANN001
        called["j"] += 1
        return {"edi": 1}

    monkeypatch.setattr(tr.AVGtoEDI, "transform", fake_x, raising=False)
    monkeypatch.setattr(tr.JtoEDI, "transform", fake_j, raising=False)

    with ctx.work_session(tmp_path) as s:
        # hit the wrapper symbol captured in _session
        b.to_edi("src")
        # transformer paths
        tr.AVGtoEDI().transform("a")
        tr.JtoEDI().transform("b")

    assert called == {"to_edi": 1, "x": 1, "j": 1}

    rec1 = s.reg.find(tag="to_edi")
    rec2 = s.reg.find(tag="AVGtoEDI.transform")
    rec3 = s.reg.find(tag="JtoEDI.transform")
    assert len(rec1) == 1
    assert len(rec2) == 1
    assert len(rec3) == 1


def test_normalize_routes_avg_path_and_edi_wrap(env, tmp_path, monkeypatch):
    ctx = _import_ctx()
    seg_edi = sys.modules["pycsamt.seg.edi"]
    seg_coll = sys.modules["pycsamt.seg.collection"]

    def tx(self, a):  # noqa: ANN001
        return seg_coll.EDICollection([seg_edi.EDIFile()])

    monkeypatch.setattr(tr.AVGtoEDI, "transform", tx, raising=False)

    avg_path = tmp_path / "a.avg"
    avg_path.write_text("", encoding="utf-8")

    with ctx.normalize_session(tmp_path) as nx:
        out = nx.load(str(avg_path))

    assert isinstance(out, seg_coll.EDICollection)

    # EDI file path → wrapped to collection
    def fe(cls, p):  # noqa: ANN001
        return cls()

    monkeypatch.setattr(seg_edi.EDIFile, "from_file", classmethod(fe))

    edi_path = tmp_path / "a.edi"
    edi_path.write_text("", encoding="utf-8")

    with ctx.normalize_session(tmp_path) as nx:
        out2 = nx.load(str(edi_path))

    # Accept either an EDICollection or a list-like of EDIFile
    assert hasattr(out2, "__iter__")
    out2 = list(out2)
    assert len(out2) == 1
    assert hasattr(out2[0], "Z")  # EDI-like (duck-typed)


def test_normalize_routes_j_path(env, tmp_path, monkeypatch):
    ctx = _import_ctx()
    seg_edi = sys.modules["pycsamt.seg.edi"]
    # seg_coll = sys.modules["pycsamt.seg.collection"]

    def tj(self, j):  # noqa: ANN001
        return seg_edi.EDIFile()

    monkeypatch.setattr(tr.JtoEDI, "transform", tj, raising=False)

    j_path = tmp_path / "a.j"
    j_path.write_text("", encoding="utf-8")

    with ctx.normalize_session(tmp_path) as nx:
        out = nx.load(str(j_path))

    assert hasattr(out, "__iter__")
    out = list(out)
    assert len(out) == 1
    assert hasattr(out[0], "Z")  # EDI-like


def test_normalize_with_avg_object_and_topo_injection(
    env, tmp_path, monkeypatch
):
    ctx = _import_ctx()
    seg_edi = sys.modules["pycsamt.seg.edi"]
    seg_coll = sys.modules["pycsamt.seg.collection"]
    zavg = sys.modules["pycsamt.zonge.avg"]

    def tx(self, a):  # noqa: ANN001
        return seg_coll.EDICollection([seg_edi.EDIFile()])

    monkeypatch.setattr(tr.AVGtoEDI, "transform", tx, raising=False)

    avg = zavg.AVG()
    topo = SimpleNamespace(frame="df")

    with ctx.normalize_session(tmp_path, topo_src=topo) as nx:
        out = nx.load(avg)

    assert avg.topo is topo
    assert isinstance(out, seg_coll.EDICollection)

    recs = nx.reg.find(tag="normalized")
    assert len(recs) >= 1


def test_normalize_passthrough_collection(env, tmp_path):
    ctx = _import_ctx()
    seg_edi = sys.modules["pycsamt.seg.edi"]
    seg_coll = sys.modules["pycsamt.seg.collection"]
    coll = seg_coll.EDICollection(items=[seg_edi.EDIFile()], verbose=0)

    with ctx.normalize_session(tmp_path) as nx:
        out = nx.load(coll)
    assert out is coll

    assert out is coll
