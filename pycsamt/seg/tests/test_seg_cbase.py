# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pycsamt.seg.cbase import CBBase, CoreParser, ParseMixin
from pycsamt.seg.edi import EDIFile


def _mk_edi(tmp: Path, station: str, nf: int = 2) -> Path:
    fvals = np.geomspace(1.0, 10.0, nf)
    fstr = "  " + "  ".join(f"{v: .6E}" for v in fvals)
    lines: list[str] = [
        ">HEAD",
        f"  DATAID={station}",
        "  LAT=26:00:00N",
        "  LONG=010:00:00E",
        "  ELEV=1000",
        "",
        ">INFO",
        "  PROJECT=SIM",
        "  PROCESSEDBY=pyCSAMT",
        "  PROCESSINGSOFTWARE=pyCSAMT",
        "",
        ">=MTSECT",
        f"  SECTID={station}",
        f"  NFREQ={nf}",
        "",
        ">!****FREQUENCIES****!",
        f">FREQ  //{nf}",
        fstr,
        "",
        ">END",
    ]
    p = tmp / f"{station}.edi"
    p.write_text("\n".join(lines), encoding="utf-8")
    return p


def test_coreparser_parse_file_and_dir(tmp_path: Path) -> None:
    f1 = _mk_edi(tmp_path, "A1", nf=3)
    _mk_edi(tmp_path, "A2", nf=4)

    pr = CoreParser(recursive=True, strict=False, on_dup="replace")
    items = pr.parse(tmp_path)

    # returns EDIFile objects
    assert isinstance(items, list)
    assert all(isinstance(x, EDIFile) for x in items)
    assert {ed.station for ed in items} == {"A1", "A2"}

    # parse a single file path too
    items2 = pr.parse([str(f1)])
    assert len(items2) == 1
    assert items2[0].station == "A1"

    # no unexpected hard errors collected
    assert pr.errors() == []


def test_coreparser_nullish_dataid_falls_back_to_stem(tmp_path: Path) -> None:
    # EMpower/WinGLink write a literal "None" into DATAID when unset; a
    # folder of such files must not collapse onto one station key.
    p1 = _mk_edi(tmp_path, "None", nf=2)
    p1 = p1.rename(tmp_path / "BH_1_imp.edi")
    p2 = _mk_edi(tmp_path, "None", nf=2)
    p2 = p2.rename(tmp_path / "BH_2_imp_rev.edi")

    assert EDIFile(p1).station == "BH_1_imp"
    assert EDIFile(p2).station == "BH_2_imp_rev"

    pr = CoreParser(recursive=True, strict=False, on_dup="replace")
    items = pr.parse(tmp_path)
    assert {ed.station for ed in items} == {"BH_1_imp", "BH_2_imp_rev"}


def test_coreparser_glob_and_errors(tmp_path: Path) -> None:
    _mk_edi(tmp_path, "B1", nf=2)
    _mk_edi(tmp_path, "B2", nf=2)

    # include a bogus path to exercise error capture
    pr = CoreParser(recursive=False, strict=False, on_dup="replace")
    items = pr.parse([tmp_path / "*.edi", tmp_path / "NOPE.edi"])
    assert len(items) == 2
    assert len(pr.errors()) >= 1  # NOPE.edi recorded


def test_cbbase_add_iter_summary(tmp_path: Path) -> None:
    f1 = _mk_edi(tmp_path, "C1", nf=2)
    f2 = _mk_edi(tmp_path, "C2", nf=5)

    ed1 = EDIFile(f1)
    ed2 = EDIFile(f2)

    col = CBBase(items=[ed1])
    assert len(col) == 1
    col.add(ed2)
    assert len(col) == 2

    st = col.stations()
    assert set(st) == {"C1", "C2"}

    sm = col.summary()
    by = {r["station"]: r for r in sm}
    assert by["C1"]["n_freq"] == 2
    assert by["C2"]["n_freq"] == 5

    # iteration yields EDIFile
    assert all(isinstance(x, EDIFile) for x in col)


def test_cbbase_repr_str(tmp_path: Path) -> None:
    f = _mk_edi(tmp_path, "D1", nf=3)
    ed = EDIFile(f)
    col = CBBase(items=[ed])
    s = str(col)
    r = repr(col)
    # light sanity: contains station and count
    assert "D1" in s
    assert "1" in r


# ─────────────────────────────────────────────────────────────────────────
# ParseMixin: low-level helpers
# ─────────────────────────────────────────────────────────────────────────


class _Finder(ParseMixin):
    recursive = True


def test_as_path_resolves_normal_path(tmp_path: Path) -> None:
    finder = _Finder()
    resolved = finder._as_path(str(tmp_path / "x.edi"))
    assert resolved.is_absolute()
    assert resolved.name == "x.edi"


def test_as_path_falls_back_to_absolute_on_resolve_oserror(
    tmp_path: Path, monkeypatch,
) -> None:
    finder = _Finder()

    def _raise(self, strict=False):
        raise OSError("simulated WinError 123")

    monkeypatch.setattr(Path, "resolve", _raise)
    result = finder._as_path(tmp_path / "*.edi")
    assert result.name == "*.edi"
    assert result.is_absolute()


def test_is_edi_path(tmp_path: Path) -> None:
    finder = _Finder()
    edi = _mk_edi(tmp_path, "P1")
    txt = tmp_path / "notes.txt"
    txt.write_text("x", encoding="utf-8")
    assert finder._is_edi_path(edi) is True
    assert finder._is_edi_path(txt) is False


def test_iter_paths_single_and_sequence(tmp_path: Path) -> None:
    finder = _Finder()
    single = list(finder._iter_paths(str(tmp_path / "a.edi")))
    assert len(single) == 1

    multi = list(finder._iter_paths([tmp_path / "a.edi", tmp_path / "b.edi"]))
    assert len(multi) == 2


def test_push_error_creates_errors_store_when_missing() -> None:
    finder = _Finder()
    assert not hasattr(finder, "_errors")
    finder._push_error("bogus.edi", "not found")
    assert len(finder._errors) == 1
    assert isinstance(finder._errors[0][1], FileNotFoundError)


def test_iter_edi_files_relative_glob_uses_root(tmp_path: Path) -> None:
    _mk_edi(tmp_path, "R1")
    _mk_edi(tmp_path, "R2")
    finder = _Finder()
    found = list(
        finder._iter_edi_files(["*.edi"], root=tmp_path)
    )
    assert len(found) == 2


def test_iter_edi_files_double_star_recursive_glob(tmp_path: Path) -> None:
    sub = tmp_path / "nested"
    sub.mkdir()
    _mk_edi(sub, "N1")
    finder = _Finder()
    found = list(finder._iter_edi_files(["**/*.edi"], root=tmp_path))
    assert len(found) == 1


def test_iter_edi_files_absolute_glob_recursive_dispatches_to_rglob(
    tmp_path: Path,
) -> None:
    sub = tmp_path / "nested2"
    sub.mkdir()
    edi_path = _mk_edi(sub, "N2")
    finder = _Finder()
    # An absolute pattern whose final component is "**" triggers the
    # ``rglob`` (recursive) branch -- the glob-pattern path only checks
    # ``**`` on that final path segment, not the whole pattern string.
    # Whether ``Path.rglob("**")`` enumerates files (not just
    # directories) is pathlib-version-dependent: Python < 3.13 yields
    # only directories for a bare ``"**"`` pattern (so nothing matches
    # here and a "no match" error is recorded), while Python >= 3.13
    # changed ``**`` to also match files. Either way exercises the
    # rglob dispatch; only the resulting outcome differs.
    pattern = str(tmp_path / "**")
    found = list(finder._iter_edi_files([pattern]))
    if any(p.is_file() for p in tmp_path.rglob("**")):
        assert found == [edi_path]
    else:
        assert found == []
        assert finder._errors  # "no match" was recorded


def test_iter_edi_files_absolute_glob_non_recursive(tmp_path: Path) -> None:
    _mk_edi(tmp_path, "N3")
    finder = _Finder()
    pattern = str(tmp_path / "*.edi")
    found = list(finder._iter_edi_files([pattern]))
    assert len(found) == 1
    assert found[0].suffix == ".edi"


def test_iter_edi_files_unsupported_source_records_error() -> None:
    finder = _Finder()
    list(finder._iter_edi_files(["\x00not-a-real-anything"]))
    assert finder._errors


def test_iter_edi_files_non_recursive_directory(tmp_path: Path) -> None:
    _mk_edi(tmp_path, "T1")
    sub = tmp_path / "nested3"
    sub.mkdir()
    _mk_edi(sub, "T2")

    class NonRecursiveFinder(ParseMixin):
        recursive = False

    finder = NonRecursiveFinder()
    found = list(finder._iter_edi_files([tmp_path]))
    assert len(found) == 1


def test_fast_station_reads_dataid_from_head(tmp_path: Path) -> None:
    edi = _mk_edi(tmp_path, "FASTSTATION")
    finder = _Finder()
    assert finder._fast_station(edi) == "FASTSTATION"


def test_fast_station_returns_none_for_missing_dataid(tmp_path: Path) -> None:
    p = tmp_path / "no_dataid.edi"
    p.write_text(">HEAD\n  LAT=0:0:0\n\n>END\n", encoding="utf-8")
    finder = _Finder()
    assert finder._fast_station(p) is None


def test_fast_station_returns_none_for_unreadable_file(tmp_path: Path) -> None:
    finder = _Finder()
    assert finder._fast_station(tmp_path / "does_not_exist.edi") is None


# ─────────────────────────────────────────────────────────────────────────
# CoreParser
# ─────────────────────────────────────────────────────────────────────────


def test_coreparser_rejects_invalid_on_dup() -> None:
    with pytest.raises(ValueError):
        CoreParser(on_dup="bogus")


def test_coreparser_strict_raises_on_bad_file(tmp_path: Path) -> None:
    bad = tmp_path / "bad.edi"
    bad.write_text("not a real edi file", encoding="utf-8")
    pr = CoreParser(strict=True)
    with pytest.raises(BaseException):  # noqa: PT011 - error type is IsEdi's
        pr.parse([bad])


def test_coreparser_non_strict_records_read_error(tmp_path: Path) -> None:
    bad = tmp_path / "bad.edi"
    bad.write_text("not a real edi file", encoding="utf-8")
    pr = CoreParser(strict=False)
    items = pr.parse([bad])
    assert items == []
    assert len(pr.errors()) >= 1


def test_coreparser_on_dup_keep_preserves_first(tmp_path: Path) -> None:
    p1 = _mk_edi(tmp_path, "DUP")
    p1 = p1.rename(tmp_path / "z_second.edi")
    p2 = tmp_path / "a_first.edi"
    p2.write_text(p1.read_text(encoding="utf-8"), encoding="utf-8")

    pr = CoreParser(recursive=True, on_dup="keep")
    items = pr.parse(tmp_path)
    assert len(items) == 1
    assert Path(items[0].path).name == "a_first.edi"


def test_coreparser_dataid_missing_falls_back_to_fast_station(
    tmp_path: Path, monkeypatch,
) -> None:
    edi = _mk_edi(tmp_path, "FS1")
    pr = CoreParser()
    monkeypatch.setattr(EDIFile, "station", property(lambda self: None))
    items = pr.parse([edi])
    assert len(items) == 1


# ─────────────────────────────────────────────────────────────────────────
# CBBase: __getitem__ / load / map / interpolate / write
# ─────────────────────────────────────────────────────────────────────────


def test_cbbase_getitem_by_index_and_station(tmp_path: Path) -> None:
    f1 = _mk_edi(tmp_path, "G1")
    col = CBBase(items=[EDIFile(f1)])
    assert col[0].station == "G1"
    assert col["G1"].station == "G1"


def test_cbbase_getitem_unknown_station_raises_keyerror(tmp_path: Path) -> None:
    f1 = _mk_edi(tmp_path, "G2")
    col = CBBase(items=[EDIFile(f1)])
    with pytest.raises(KeyError):
        col["does-not-exist"]


def test_cbbase_load_classmethod(tmp_path: Path) -> None:
    _mk_edi(tmp_path, "L1")
    _mk_edi(tmp_path, "L2")
    col = CBBase.load(tmp_path, recursive=True)
    assert isinstance(col, CBBase)
    assert len(col) == 2


def test_cbbase_load_reports_errors(tmp_path: Path, caplog) -> None:
    _mk_edi(tmp_path, "L3")
    bad = tmp_path / "bad.edi"
    bad.write_text("not a real edi", encoding="utf-8")
    col = CBBase.load(tmp_path, recursive=True, strict=False)
    assert len(col) == 1


def test_cbbase_map_applies_function(tmp_path: Path) -> None:
    f1 = _mk_edi(tmp_path, "M1")
    f2 = _mk_edi(tmp_path, "M2")
    col = CBBase(items=[EDIFile(f1), EDIFile(f2)])
    result = col.map(lambda ed: ed.station)
    assert result == ["M1", "M2"]


def test_cbbase_interpolate_builds_new_collection(tmp_path: Path) -> None:
    f1 = _mk_edi(tmp_path, "I1", nf=4)
    col = CBBase(items=[EDIFile(f1)])
    new_freq = np.geomspace(1.5, 8.0, 3)
    out = col.interpolate(new_freq, bounds_error=False)
    assert isinstance(out, CBBase)
    assert len(out) == 1
    assert out[0].Z.n_freq == 3


def test_cbbase_write_creates_files(tmp_path: Path) -> None:
    f1 = _mk_edi(tmp_path, "W1")
    col = CBBase(items=[EDIFile(f1)])
    out_dir = tmp_path / "out"
    paths = col.write(out_dir)
    assert len(paths) == 1
    assert Path(paths[0]).exists()


def test_cbbase_items_property(tmp_path: Path) -> None:
    f1 = _mk_edi(tmp_path, "IT1")
    col = CBBase(items=[EDIFile(f1)])
    assert col.items == col._items


def test_cbbase_summary_reports_tipper_present(tmp_path: Path) -> None:
    n = 3
    fvals = np.geomspace(1.0, 10.0, n)
    fstr = "  " + "  ".join(f"{v: .6E}" for v in fvals)
    txr = "  " + "  ".join(f"{v: .6E}" for v in np.linspace(0.1, 0.2, n))
    zeros = "  " + "  ".join(f"{0.0: .6E}" for _ in range(n))
    lines = [
        ">HEAD",
        "  DATAID=TIP1",
        "",
        ">=MTSECT",
        "  SECTID=TIP1",
        f"  NFREQ={n}",
        "",
        ">!****FREQUENCIES****!",
        f">FREQ  //{n}",
        fstr,
        "",
        ">!****TIPPER PARAMETERS****!",
        f">TXR.EXP ROT=TROT  //{n}",
        txr,
        f">TXI.EXP ROT=TROT  //{n}",
        zeros,
        f">TXVAR.EXP ROT=TROT  //{n}",
        zeros,
        f">TYR.EXP ROT=TROT  //{n}",
        txr,
        f">TYI.EXP ROT=TROT  //{n}",
        zeros,
        f">TYVAR.EXP ROT=TROT  //{n}",
        zeros,
        "",
        ">END",
    ]
    p = tmp_path / "TIP1.edi"
    p.write_text("\n".join(lines), encoding="utf-8")

    col = CBBase(items=[EDIFile(p)])
    rows = col.summary()
    assert rows[0]["tipper"] is True


def test_cbbase_add_replaces_existing_station(tmp_path: Path) -> None:
    f1 = _mk_edi(tmp_path, "REPL", nf=2)
    ed1 = EDIFile(f1)
    f2 = _mk_edi(tmp_path, "REPL", nf=5)  # overwrites the same file
    ed2 = EDIFile(f2)
    col = CBBase()
    col.add(ed1)
    col.add(ed2)
    assert len(col) == 1
    assert col["REPL"].Z.n_freq == 5
