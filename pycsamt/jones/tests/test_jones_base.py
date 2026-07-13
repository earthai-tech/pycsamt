from pathlib import Path

import pytest

from pycsamt.jones.base import BaseJones, JComponentBase


class DummyComponent(JComponentBase):
    _repr_keys = ["flag"]

    def __init__(self, *, flag: bool | None = None, verbose: int = 0):
        super().__init__(verbose=verbose)
        self.flag = flag

    @classmethod
    def from_lines(cls, lines, *, verbose: int = 0, **kws):
        inst = cls(flag=True, verbose=verbose)
        # mark read to simulate a successful parse
        inst._mark_read(True)
        return inst

    def read(self, *args, **kws):
        self._mark_read(True)
        return self

    def write(self, obj=None, *, encoding=None, **kws):
        return ["DUMMY"]


class DummyJones(BaseJones):
    @classmethod
    def from_lines(cls, lines, *, verbose: int = 0, **kws):
        inst = cls(verbose=verbose)
        inst._mark_read(True)
        return inst

    def read(self, *args, **kws):
        self._mark_read(True)
        return self

    def write(self, obj=None, *, encoding=None, **kws):
        return ["JONES"]


def test_from_file_sets_path_and_encoding(tmp_path: Path):
    p = tmp_path / "dummy.j"
    p.write_text("# header\n>AZIMUTH=0\nS01\nZXY SI\n1\n", encoding="utf-8")

    inst = DummyComponent.from_file(p, verbose=1)
    assert inst.path == p
    assert inst.encoding == "utf-8"
    assert inst.__has_read__() is True


def test_require_read_helper(tmp_path: Path):
    inst = DummyJones(verbose=1)
    with pytest.raises(Exception):
        inst.require_read()
    inst._mark_read(True)
    # should not raise now
    inst.require_read()


def test_to_dict_and_repr():
    d = DummyComponent(flag=True, verbose=0)
    m = d.to_dict()
    assert "flag" in m and m["flag"] is True
    s = repr(d)
    assert "DummyComponent(" in s
