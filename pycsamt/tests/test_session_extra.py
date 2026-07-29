from __future__ import annotations

import importlib
import sys
from types import SimpleNamespace

import pytest

from pycsamt.tests.test_sessions import _install_stubs


@pytest.fixture()
def env(monkeypatch):
    return _install_stubs(monkeypatch)


def _smod():
    return importlib.import_module("pycsamt._session")


# --------------------------------------------------------------------------
# Session._record
# --------------------------------------------------------------------------


def test_record_add_object_exception_is_ignored(env, tmp_path):
    smod = _smod()
    ses = smod.Session(tmp_path)
    ses.reg.add_object = lambda obj, tags=None: (_ for _ in ()).throw(
        RuntimeError("no")
    )
    ses._record(object(), tag="x")  # must not raise


def test_record_capture_children_respects_max(env, tmp_path):
    smod = _smod()
    ses = smod.Session(tmp_path, capture_children=True, max_children=2)
    calls = []
    ses.reg.add_object = lambda obj, tags=None: calls.append(obj)
    ses._record([1, 2, 3, 4], tag="t")
    assert calls == [[1, 2, 3, 4], 1, 2]


def test_record_capture_children_iteration_exception_ignored(env, tmp_path):
    smod = _smod()

    class BadIter:
        def __iter__(self):
            raise RuntimeError("boom")

    ses = smod.Session(tmp_path, capture_children=True)
    ses._record(BadIter(), tag="t")  # must not raise


def test_record_capture_children_non_iterable_object(env, tmp_path):
    smod = _smod()
    ses = smod.Session(tmp_path, capture_children=True)
    ses._record(object(), tag="t")  # not iterable -> loop body skipped


def test_record_capture_children_empty_iterable(env, tmp_path):
    smod = _smod()
    ses = smod.Session(tmp_path, capture_children=True)
    calls = []
    ses.reg.add_object = lambda obj, tags=None: calls.append(obj)
    ses._record([], tag="t")
    assert calls == [[]]


# --------------------------------------------------------------------------
# Session wrapping
# --------------------------------------------------------------------------


def test_wrap_to_edi_noop_when_auto_capture_false(env, tmp_path):
    smod = _smod()
    from pycsamt.core import base as b

    orig = b.to_edi
    ses = smod.Session(tmp_path, auto_capture=False)
    ses._wrap_to_edi()
    assert b.to_edi is orig
    assert ses._orig_to_edi is None


def test_wrap_transformers_noop_when_auto_capture_false(env, tmp_path):
    smod = _smod()
    from pycsamt.transformers import jedi as tr

    orig_x, orig_j = tr.AVGtoEDI.transform, tr.JtoEDI.transform
    ses = smod.Session(tmp_path, auto_capture=False)
    ses._wrap_transformers()
    assert tr.AVGtoEDI.transform is orig_x
    assert tr.JtoEDI.transform is orig_j


def test_wrap_to_edi_idempotent(env, tmp_path):
    smod = _smod()
    ses = smod.Session(tmp_path)
    ses._wrap_to_edi()
    first_orig = ses._orig_to_edi
    ses._wrap_to_edi()  # second call must not re-capture the original
    assert ses._orig_to_edi is first_orig
    ses._restore()


def test_wrap_transformers_idempotent(env, tmp_path):
    smod = _smod()
    ses = smod.Session(tmp_path)
    ses._wrap_transformers()
    ox, oj = ses._orig_t_x, ses._orig_t_j
    ses._wrap_transformers()  # both already wrapped -> skip re-wrap
    assert ses._orig_t_x is ox
    assert ses._orig_t_j is oj
    ses._restore()


def test_wrap_transformers_tag_fallback_on_class_access_error(
    env, tmp_path, monkeypatch
):
    smod = _smod()
    from pycsamt.transformers import jedi as tr

    def fake_transform(self, *a, **k):
        return "result"

    monkeypatch.setattr(tr.AVGtoEDI, "transform", fake_transform, raising=False)

    ses = smod.Session(tmp_path)
    ses._wrap_transformers()

    class Weird:
        def __getattribute__(self, name):
            if name == "__class__":
                raise RuntimeError("blocked")
            return object.__getattribute__(self, name)

    out = tr.AVGtoEDI.transform(Weird())
    assert out == "result"
    recs = ses.reg.find(tag="transform")
    assert len(recs) == 1
    ses._restore()


def test_restore_noop_when_never_wrapped(env, tmp_path):
    smod = _smod()
    ses = smod.Session(tmp_path)
    ses._restore()  # no-op, must not raise


def test_session_exit_reg_save_exception_ignored(env, tmp_path):
    smod = _smod()
    ses = smod.Session(tmp_path)
    ses.reg.save = lambda: (_ for _ in ()).throw(RuntimeError("no"))
    with ses:
        pass  # __exit__ must swallow the save() exception


# --------------------------------------------------------------------------
# Normalize._as_edi_coll
# --------------------------------------------------------------------------


def test_as_edi_coll_duck_typed_single_item(env, tmp_path):
    smod = _smod()
    seg_coll = sys.modules["pycsamt.seg.collection"]
    nz = smod.Normalize(tmp_path)
    duck = SimpleNamespace(Z=SimpleNamespace())
    out = nz._as_edi_coll(duck)
    assert isinstance(out, seg_coll.EDICollection)
    assert list(out) == [duck]


def test_as_edi_coll_list_passthrough(env, tmp_path):
    smod = _smod()
    seg_edi = sys.modules["pycsamt.seg.edi"]
    nz = smod.Normalize(tmp_path)
    items = [seg_edi.EDIFile(), seg_edi.EDIFile()]
    out = nz._as_edi_coll(items)
    assert out is items


def test_as_edi_coll_list_with_non_edi_item_returns_none(env, tmp_path):
    smod = _smod()
    seg_edi = sys.modules["pycsamt.seg.edi"]
    nz = smod.Normalize(tmp_path)
    items = [seg_edi.EDIFile(), object()]
    assert nz._as_edi_coll(items) is None


# --------------------------------------------------------------------------
# Normalize._try_topo
# --------------------------------------------------------------------------


def test_try_topo_none_source_noop(env, tmp_path):
    smod = _smod()
    nz = smod.Normalize(tmp_path)
    avg = SimpleNamespace()
    nz._try_topo(avg)
    assert not hasattr(avg, "topo")


def test_try_topo_uses_add_topography(env, tmp_path):
    smod = _smod()
    topo = SimpleNamespace(frame="df")
    nz = smod.Normalize(tmp_path, topo_src=topo)
    calls = []
    avg = SimpleNamespace(add_topography=lambda t: calls.append(t))
    nz._try_topo(avg)
    assert calls == [topo]


def test_try_topo_add_topography_raises_falls_back_to_frame(env, tmp_path):
    smod = _smod()
    topo = SimpleNamespace(frame="df")
    nz = smod.Normalize(tmp_path, topo_src=topo)

    def boom(t):
        raise RuntimeError("no")

    avg = SimpleNamespace(add_topography=boom)
    nz._try_topo(avg)
    assert avg.topo is topo


def test_try_topo_no_add_topography_falls_to_frame_assignment(env, tmp_path):
    smod = _smod()
    topo = SimpleNamespace(frame="df")
    nz = smod.Normalize(tmp_path, topo_src=topo)
    avg = SimpleNamespace()
    nz._try_topo(avg)
    assert avg.topo is topo


def test_try_topo_source_without_frame_is_noop(env, tmp_path):
    smod = _smod()
    topo = SimpleNamespace()
    nz = smod.Normalize(tmp_path, topo_src=topo)
    avg = SimpleNamespace()
    nz._try_topo(avg)
    assert not hasattr(avg, "topo")


def test_try_topo_second_try_exception_ignored(env, tmp_path):
    smod = _smod()

    class NoTopoAssign:
        def __setattr__(self, name, value):
            if name == "topo":
                raise RuntimeError("cannot set")
            object.__setattr__(self, name, value)

    topo = SimpleNamespace(frame="df")
    nz = smod.Normalize(tmp_path, topo_src=topo)
    avg = NoTopoAssign()
    nz._try_topo(avg)  # must not raise despite the failed assignment


# --------------------------------------------------------------------------
# Normalize._to_avg
# --------------------------------------------------------------------------


def test_to_avg_instance_passthrough(env, tmp_path):
    smod = _smod()
    zavg = sys.modules["pycsamt.zonge.avg"]
    nz = smod.Normalize(tmp_path)
    avg = zavg.AVG()
    assert nz._to_avg(avg) is avg


def test_to_avg_duck_typed_passthrough(env, tmp_path):
    smod = _smod()
    nz = smod.Normalize(tmp_path)
    duck = SimpleNamespace(add_topography=lambda t: None)
    assert nz._to_avg(duck) is duck


def test_to_avg_path_with_avg_suffix(env, tmp_path):
    smod = _smod()
    zavg = sys.modules["pycsamt.zonge.avg"]
    nz = smod.Normalize(tmp_path)
    p = tmp_path / "x.avg"
    p.write_text("", encoding="utf-8")
    out = nz._to_avg(str(p))
    assert isinstance(out, zavg.AVG)


def test_to_avg_returns_none_for_unrelated_input(env, tmp_path):
    smod = _smod()
    nz = smod.Normalize(tmp_path)
    assert nz._to_avg("not/a/path.xyz") is None
    assert nz._to_avg(12345) is None


# --------------------------------------------------------------------------
# Normalize._to_j
# --------------------------------------------------------------------------


def test_to_j_jfile_instance_passthrough(env, tmp_path):
    smod = _smod()
    jj = sys.modules["pycsamt.jones.j"]
    nz = smod.Normalize(tmp_path)
    j = jj.JFile()
    assert nz._to_j(j) is j


def test_to_j_jcollection_passthrough_nonempty(env, tmp_path):
    smod = _smod()
    jcoll = sys.modules["pycsamt.jones.collection"]
    nz = smod.Normalize(tmp_path)
    coll = jcoll.JCollection([1, 2])
    assert nz._to_j(coll) is coll


def test_to_j_jcollection_empty_returns_none(env, tmp_path):
    smod = _smod()
    jcoll = sys.modules["pycsamt.jones.collection"]
    nz = smod.Normalize(tmp_path)
    coll = jcoll.JCollection([])
    assert nz._to_j(coll) is None


def test_to_j_path_jones_suffix_success(env, tmp_path):
    smod = _smod()
    jj = sys.modules["pycsamt.jones.j"]
    nz = smod.Normalize(tmp_path)
    p = tmp_path / "a.j"
    p.write_text("", encoding="utf-8")
    out = nz._to_j(str(p))
    assert isinstance(out, jj.JFile)


def test_to_j_path_jones_suffix_from_file_raises_returns_none(
    env, tmp_path, monkeypatch
):
    smod = _smod()
    jj = sys.modules["pycsamt.jones.j"]
    nz = smod.Normalize(tmp_path)

    def boom(cls, path):
        raise RuntimeError("no")

    monkeypatch.setattr(jj.JFile, "from_file", classmethod(boom))

    p = tmp_path / "a.dat"
    p.write_text("", encoding="utf-8")
    assert nz._to_j(str(p)) is None


def test_to_j_folder_glob_uses_jcollection_from_sources(env, tmp_path, monkeypatch):
    smod = _smod()
    jcoll = sys.modules["pycsamt.jones.collection"]
    nz = smod.Normalize(tmp_path)

    def fake_from_sources(sources, verbose=0):
        return jcoll.JCollection([1])

    monkeypatch.setattr(
        jcoll.JCollection,
        "from_sources",
        staticmethod(fake_from_sources),
        raising=False,
    )

    out = nz._to_j(str(tmp_path))  # a directory -> empty suffix
    assert isinstance(out, jcoll.JCollection)
    assert len(out) == 1


def test_to_j_folder_glob_from_sources_raises_returns_none(env, tmp_path, monkeypatch):
    smod = _smod()
    jcoll = sys.modules["pycsamt.jones.collection"]
    nz = smod.Normalize(tmp_path)

    def boom(sources, verbose=0):
        raise RuntimeError("no")

    monkeypatch.setattr(
        jcoll.JCollection, "from_sources", staticmethod(boom), raising=False
    )
    assert nz._to_j(str(tmp_path)) is None


def test_to_j_raw_list_of_jfile_normalizes(env, tmp_path, monkeypatch):
    smod = _smod()
    jj = sys.modules["pycsamt.jones.j"]
    jcoll = sys.modules["pycsamt.jones.collection"]

    class WorkingJCollection(list):
        def __init__(self, items=None, verbose=0):
            super().__init__(items or [])

    monkeypatch.setattr(jcoll, "JCollection", WorkingJCollection, raising=False)
    monkeypatch.setattr(smod, "JCollection", WorkingJCollection, raising=False)

    nz = smod.Normalize(tmp_path)
    items = [jj.JFile(), jj.JFile()]
    out = nz._to_j(items)
    assert isinstance(out, WorkingJCollection)
    assert len(out) == 2


def test_to_j_raw_list_construction_type_error_returns_none(env, tmp_path, monkeypatch):
    smod = _smod()
    jj = sys.modules["pycsamt.jones.j"]
    jcoll = sys.modules["pycsamt.jones.collection"]

    class StrictJCollection(list):
        def __init__(self, seq=()):
            super().__init__(seq)

    monkeypatch.setattr(jcoll, "JCollection", StrictJCollection, raising=False)
    monkeypatch.setattr(smod, "JCollection", StrictJCollection, raising=False)

    nz = smod.Normalize(tmp_path)
    items = [jj.JFile(), jj.JFile()]
    # StrictJCollection's __init__ doesn't accept items=/verbose= kwargs,
    # so construction raises TypeError, which is swallowed -> None.
    assert nz._to_j(items) is None


def test_to_j_raw_list_construction_empty_result_returns_none(env, tmp_path):
    smod = _smod()
    jj = sys.modules["pycsamt.jones.j"]
    nz = smod.Normalize(tmp_path)
    items = [jj.JFile(), jj.JFile()]
    # Default stub JCollection(list) silently ignores items=/verbose=
    # kwargs and yields an empty collection -> treated as None.
    assert nz._to_j(items) is None


def test_to_j_raw_list_without_jfile_returns_none(env, tmp_path):
    smod = _smod()
    nz = smod.Normalize(tmp_path)
    assert nz._to_j([1, 2, 3]) is None


# --------------------------------------------------------------------------
# Normalize._normalize
# --------------------------------------------------------------------------


def test_normalize_ultimate_fallback_uses_to_edi(env, tmp_path, monkeypatch):
    smod = _smod()
    nz = smod.Normalize(tmp_path)

    sentinel = object()
    monkeypatch.setattr(smod, "to_edi", lambda src: sentinel)

    out = nz._normalize(12345)  # matches no route
    assert out is sentinel


def test_normalize_j_suffix_returns_none_falls_through(env, tmp_path, monkeypatch):
    smod = _smod()
    nz = smod.Normalize(tmp_path)
    monkeypatch.setattr(nz, "_to_j", lambda src: None)

    sentinel = object()
    monkeypatch.setattr(smod, "to_edi", lambda src: sentinel)

    p = tmp_path / "x.j"
    p.write_text("", encoding="utf-8")
    out = nz._normalize(str(p))
    assert out is sentinel


def test_normalize_avg_suffix_to_avg_none_falls_through(env, tmp_path, monkeypatch):
    smod = _smod()
    nz = smod.Normalize(tmp_path)
    monkeypatch.setattr(nz, "_to_avg", lambda src: None)

    sentinel = object()
    monkeypatch.setattr(smod, "to_edi", lambda src: sentinel)

    p = tmp_path / "x.avg"
    p.write_text("", encoding="utf-8")
    out = nz._normalize(str(p))
    assert out is sentinel


def test_normalize_edi_suffix_from_file_exception_falls_through(
    env, tmp_path, monkeypatch
):
    smod = _smod()
    seg_edi = sys.modules["pycsamt.seg.edi"]
    nz = smod.Normalize(tmp_path)

    def boom(cls, path):
        raise RuntimeError("bad edi")

    monkeypatch.setattr(seg_edi.EDIFile, "from_file", classmethod(boom))
    monkeypatch.setattr(nz, "_to_j", lambda src: None)
    monkeypatch.setattr(nz, "_to_avg", lambda src: None)

    sentinel = object()
    monkeypatch.setattr(smod, "to_edi", lambda src: sentinel)

    p = tmp_path / "x.edi"
    p.write_text("", encoding="utf-8")
    out = nz._normalize(str(p))
    assert out is sentinel


def test_normalize_empty_suffix_uses_edicollection_from_sources(
    env, tmp_path, monkeypatch
):
    smod = _smod()
    seg_coll = sys.modules["pycsamt.seg.collection"]
    seg_edi = sys.modules["pycsamt.seg.edi"]
    nz = smod.Normalize(tmp_path)

    def fake_from_sources(sources, verbose=0):
        return seg_coll.EDICollection([seg_edi.EDIFile()])

    monkeypatch.setattr(
        seg_coll.EDICollection,
        "from_sources",
        staticmethod(fake_from_sources),
        raising=False,
    )

    out = nz._normalize(str(tmp_path))  # directory -> empty suffix
    assert isinstance(out, seg_coll.EDICollection)
    assert len(out) == 1


def test_normalize_empty_suffix_from_sources_raises_falls_through(
    env, tmp_path, monkeypatch
):
    smod = _smod()
    seg_coll = sys.modules["pycsamt.seg.collection"]
    nz = smod.Normalize(tmp_path)

    def boom(sources, verbose=0):
        raise RuntimeError("no")

    monkeypatch.setattr(
        seg_coll.EDICollection,
        "from_sources",
        staticmethod(boom),
        raising=False,
    )
    monkeypatch.setattr(nz, "_to_j", lambda src: None)
    monkeypatch.setattr(nz, "_to_avg", lambda src: None)

    sentinel = object()
    monkeypatch.setattr(smod, "to_edi", lambda src: sentinel)

    out = nz._normalize(str(tmp_path))
    assert out is sentinel


def test_normalize_empty_suffix_from_sources_empty_falls_through(
    env, tmp_path, monkeypatch
):
    smod = _smod()
    seg_coll = sys.modules["pycsamt.seg.collection"]
    nz = smod.Normalize(tmp_path)

    def fake_from_sources(sources, verbose=0):
        return seg_coll.EDICollection([])

    monkeypatch.setattr(
        seg_coll.EDICollection,
        "from_sources",
        staticmethod(fake_from_sources),
        raising=False,
    )
    monkeypatch.setattr(nz, "_to_j", lambda src: None)
    monkeypatch.setattr(nz, "_to_avg", lambda src: None)

    sentinel = object()
    monkeypatch.setattr(smod, "to_edi", lambda src: sentinel)

    out = nz._normalize(str(tmp_path))
    assert out is sentinel


def test_normalize_step3_jones_instance_returns_collection_directly(
    env, tmp_path, monkeypatch
):
    smod = _smod()
    from pycsamt.transformers import jedi as tr

    jj = sys.modules["pycsamt.jones.j"]
    seg_coll = sys.modules["pycsamt.seg.collection"]

    coll = seg_coll.EDICollection([])
    monkeypatch.setattr(tr.JtoEDI, "transform", lambda self, j: coll, raising=False)

    nz = smod.Normalize(tmp_path)
    j = jj.JFile()
    out = nz._normalize(j)
    assert out is coll


# --------------------------------------------------------------------------
# Normalize.__exit__ / Normalize.load
# --------------------------------------------------------------------------


def test_normalize_exit_reg_save_exception_ignored(env, tmp_path):
    smod = _smod()
    nz = smod.Normalize(tmp_path)
    nz.reg.save = lambda: (_ for _ in ()).throw(RuntimeError("no"))
    with nz:
        pass  # __exit__ must swallow the save() exception


def test_normalize_load_auto_register_exception_ignored(env, tmp_path, monkeypatch):
    smod = _smod()
    nz = smod.Normalize(tmp_path, auto_register=True)
    monkeypatch.setattr(nz, "_normalize", lambda source: "result")
    nz.reg.add_object = lambda obj, tags=None: (_ for _ in ()).throw(RuntimeError("no"))
    out = nz.load("anything")
    assert out == "result"


def test_normalize_load_auto_register_false_skips_registration(
    env, tmp_path, monkeypatch
):
    smod = _smod()
    nz = smod.Normalize(tmp_path, auto_register=False)
    monkeypatch.setattr(nz, "_normalize", lambda source: "result")
    calls = []
    nz.reg.add_object = lambda obj, tags=None: calls.append(obj)
    out = nz.load("anything")
    assert out == "result"
    assert calls == []
