from __future__ import annotations

import csv
import zipfile
from pathlib import Path

import pytest

from pycsamt.site import export as xp

# ---------------------------- test doubles --------------------------------- #

class _Head:
    def __init__(self, station: str, lat=1.0, lon=2.0, elev=3.0):
        self.station = station
        self.dataid = station
        self.lat = lat
        # be robust to utils picking either "lon" or "long"
        self.lon = lon
        self.long = lon
        self.elev = elev


class EdiWriteKwarg:
    """Backend exposes write(new_edifn=..., savepath=...)."""

    def __init__(self, station: str):
        self.HEAD = _Head(station)
        self.name = station  # allow station_name fallbacks

    def write(self, *_, new_edifn: str | None = None,
              savepath: str | Path | None = None, **__):
        assert new_edifn is not None
        # Honour savepath like the real backend — ignoring it would write
        # a bare filename relative to the CWD (i.e. into the repo root).
        p = Path(savepath) / new_edifn if savepath else Path(new_edifn)
        p.write_text(f"# {self.HEAD.station}\n", encoding="utf-8")


class EdiWritePos:
    """Backend exposes write(path)."""

    def __init__(self, station: str):
        self.HEAD = _Head(station)
        self.name = station

    def write(self, path: str):
        Path(path).write_text(
            f"# {self.HEAD.station}\n", encoding="utf-8"
        )


class EdiToFile:
    """Backend exposes to_file(path)."""

    def __init__(self, station: str):
        self.HEAD = _Head(station)
        self.name = station

    def to_file(self, path: str):
        Path(path).write_text(
            f"# {self.HEAD.station}\n", encoding="utf-8"
        )


class EdiSave:
    """Backend exposes save(path)."""

    def __init__(self, station: str):
        self.HEAD = _Head(station)
        self.name = station

    def save(self, path: str):
        Path(path).write_text(
            f"# {self.HEAD.station}\n", encoding="utf-8"
        )


# ------------------------------- tests ------------------------------------- #

def test_write_site_supports_all_backends(tmp_path: Path) -> None:
    eds = [
        (EdiWriteKwarg("A01"), tmp_path / "a01.edi"),
        (EdiWritePos("A02"), tmp_path / "a02.edi"),
        (EdiToFile("A03"), tmp_path / "a03.edi"),
        (EdiSave("A04"), tmp_path / "a04.edi"),
    ]
    for ed, dest in eds:
        out = xp.write_site(ed, dest)
        assert out == dest
        assert dest.exists()
        txt = dest.read_text(encoding="utf-8")
        assert ed.HEAD.station in txt


def test_write_sites_template_and_manifest(tmp_path: Path) -> None:
    eds = [EdiToFile("S01"), EdiToFile("S02")]
    outdir = tmp_path / "eds"
    mpath = tmp_path / "manifest.csv"

    written = xp.write_sites(
        eds,
        outdir,
        template="{index:03d}_{station}",
        exist_ok=False,
        manifest_csv=mpath,
    )
    assert len(written) == 2
    # template should auto-append .edi
    names = sorted(p.name for p in written)
    assert names == ["000_S01.edi", "001_S02.edi"]
    for p in written:
        assert p.exists()
        assert p.read_text(encoding="utf-8").startswith("# ")

    # manifest checks
    assert mpath.exists()
    with mpath.open("r", encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 2
    assert set(rows[0].keys()) == {
        "index",
        "station",
        "lat",
        "lon",
        "elev",
        "chainage",
        "filename",
        "path",
    }
    # row content sanity
    stations = {r["station"] for r in rows}
    assert stations == {"S01", "S02"}


def test_write_sites_exist_ok_behavior(tmp_path: Path) -> None:
    eds = [EdiWritePos("D01")]
    outdir = tmp_path / "dup"
    first = xp.write_sites(eds, outdir)
    assert len(first) == 1
    # second attempt without exist_ok should raise
    with pytest.raises(FileExistsError):
        xp.write_sites(eds, outdir, exist_ok=False)
    # but with exist_ok=True it should overwrite
    second = xp.write_sites(eds, outdir, exist_ok=True)
    assert len(second) == 1
    assert second[0].exists()


def test_pack_zip_creates_archive_and_optional_manifest(
    tmp_path: Path,
) -> None:
    eds = [EdiSave("Z01"), EdiSave("Z02")]
    zpath = tmp_path / "bundle" / "sites.zip"
    mpath = tmp_path / "bundle" / "manifest.csv"

    out = xp.pack_zip(
        eds,
        zpath,
        template="{station}.edi",
        manifest_csv=mpath,
    )
    assert out == zpath
    assert zpath.exists()

    # zip content
    with zipfile.ZipFile(zpath, "r") as zf:
        names = sorted(zf.namelist())
    assert names == ["Z01.edi", "Z02.edi"]

    # manifest exists and lists the archive path
    assert mpath.exists()
    with mpath.open("r", encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 2
    assert all(Path(r["path"]) == zpath for r in rows)


def test_template_without_extension_appends_dot_edi(tmp_path: Path) -> None:
    eds = [EdiToFile("N01")]
    outdir = tmp_path / "noext"
    paths = xp.write_sites(eds, outdir, template="X_{station}")
    assert len(paths) == 1
    assert paths[0].name == "X_N01.edi"


def test_write_sites_disambiguates_duplicate_station_names(
    tmp_path: Path,
) -> None:
    """Sites sharing a station name must not overwrite each other.

    Regression: a collision used to silently write a single file while
    returning one path per input, destroying all but the last station —
    the re-loaded export then failed with "no valid impedance".
    """
    eds = [EdiToFile("S"), EdiToFile("S"), EdiToFile("S")]
    outdir = tmp_path / "dups"
    paths = xp.write_sites(eds, outdir, exist_ok=True)
    assert len(paths) == 3
    # every input must produce a distinct file actually on disk
    on_disk = sorted(p.name for p in outdir.glob("*.edi"))
    assert len(on_disk) == 3
    assert len(set(paths)) == 3


def test_write_sites_empty_station_name_is_not_stem_less(
    tmp_path: Path,
) -> None:
    """An empty station name must not render to a stem-less ".edi".

    Regression: such files are hidden/ignored by ``ensure_sites``, so the
    export loaded back as zero stations and corrupted the session.
    """
    eds = [EdiToFile(""), EdiToFile("")]
    outdir = tmp_path / "empties"
    paths = xp.write_sites(eds, outdir, exist_ok=True)
    assert len(paths) == 2
    names = sorted(p.name for p in paths)
    assert ".edi" not in names          # no stem-less file
    assert all(len(n) > len(".edi") for n in names)
    assert len(sorted(outdir.glob("*.edi"))) == 2
