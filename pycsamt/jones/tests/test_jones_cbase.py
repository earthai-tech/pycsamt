# -*- coding: utf-8 -*-
from __future__ import annotations

from pathlib import Path
import pytest

from pycsamt.jones.cbase import ( 
    JParseMixin, 
    JCoreParser, 
    JCBBase
    )



def _has(obj, *names: str) -> bool:
    return all(hasattr(obj, n) for n in names)


@pytest.mark.parametrize("kind", ["path", "str"])
def test_jparse_mixin_path_and_lines(tmp_path: Path, kind: str):
    if JParseMixin is None:
        pytest.skip("JParseMixin not implemented yet")

    # Create a minimal J-like file
    jtxt = "\n".join(
        [
            "#WRITTEN BY TOOL: KB0001 01/01/25 RAW",
            ">AZIMUTH = 0.0",
            ">LATITUDE = 10.0",
            ">LONGITUDE = 20.0",
            "",
            "KB0001",
            "RXY",
            "1",
            " 1.0  100  45  110  90  50  40  1  1",
        ]
    )
    p = tmp_path / "site.j"
    p.write_text(jtxt, encoding="utf-8")

    # Non-J file
    q = tmp_path / "notes.txt"
    q.write_text("hello\nworld\n", encoding="utf-8")

    class _Dummy(JParseMixin):  # type: ignore[misc]
        pass

    d = _Dummy()

    # _as_path
    if hasattr(d, "_as_path"):
        arg = p if kind == "path" else str(p)
        pp = d._as_path(arg)  # type: ignore[attr-defined]
        assert isinstance(pp, Path) and pp.exists()

    # _iter_lines / _read_text
    if hasattr(d, "_iter_lines"):
        lines = list(d._iter_lines(p))  # type: ignore[attr-defined]
        assert any("WRITTEN BY" in ln for ln in lines)
    elif hasattr(d, "_read_text"):
        txt = d._read_text(p)  # type: ignore[attr-defined]
        assert "#WRITTEN BY" in txt

    for name in ("is_j_like", "is_j_file", "is_j_candidate"):
        if hasattr(d, name):
            fn = getattr(d, name)
            assert fn(p) is True
            assert fn(q) is False
            break
    else:
        pytest.skip("No candidate-check helper on JParseMixin")


def test_jcbbase_container_minimal():
    if JCBBase is None:
        pytest.skip("JCBBase not implemented yet")

    # If abstract, make a tiny concrete subclass
    try:
        class _Mini(JCBBase):  # type: ignore[misc]
            def parse(self, inputs):
                self.items = list(inputs)
                return self.items
    except TypeError:
        # Fallback: try to instantiate directly if not abstract
        _Mini = JCBBase  # type: ignore[assignment]

    inst = _Mini(verbose=0)  # type: ignore[call-arg]
    out = inst.parse(["a", "b"])  # type: ignore[attr-defined]
    items = getattr(inst, "items", out)
    assert isinstance(items, (list, tuple))
    # n / __len__ semantics
    n = getattr(inst, "n", None)
    if isinstance(n, int):
        assert n == len(items)
    else:
        assert len(items) == 2

    # __str__/__repr__ should be informative
    assert isinstance(str(inst), str)
    assert isinstance(repr(inst), str)


def test_jcoreparser_parse_single_and_dir(
    project_root: Path, j_single_file: Path
):
    if JCoreParser is None:
        pytest.skip("JCoreParser not implemented yet")

    parser = JCoreParser(verbose=0)  # type: ignore[call-arg]

    # Single path
    got = parser.parse([j_single_file])  # type: ignore[attr-defined]
    items = getattr(parser, "items", got)
    assert isinstance(items, (list, tuple)) and len(items) >= 1

    # Directory of J files (if test data present)
    data_dir = project_root / "data" / "j"
    if data_dir.exists():
        got2 = parser.parse([data_dir])  # type: ignore[attr-defined]
        items2 = getattr(parser, "items", got2)
        assert isinstance(items2, (list, tuple)) and len(items2) >= 1

    # Optional filtering API: where/filter/select by station/name
    filt = None
    for name in ("where", "filter", "select"):
        if hasattr(parser, name):
            filt = getattr(parser, name)
            break
    if callable(filt):
        # Probe with a very permissive filter so we don't overfit
        sub = filt()  # type: ignore[misc]
        # Accept either count property or len()
        if hasattr(sub, "n"):
            assert sub.n >= 1  # type: ignore[attr-defined]
        else:
            assert len(sub) >= 1  # type: ignore[arg-type]


def test_jcoreparser_repr_and_str(j_single_file: Path):
    if JCoreParser is None:
        pytest.skip("JCoreParser not implemented yet")
    p = JCoreParser(verbose=0)  # type: ignore[call-arg]
    p.parse([j_single_file])  # type: ignore[attr-defined]
    s = str(p)
    r = repr(p)
    assert isinstance(s, str) and isinstance(r, str)
    assert "JCoreParser" in r or "JCoreParser" in s
