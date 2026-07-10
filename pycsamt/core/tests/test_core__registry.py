from __future__ import annotations

from pathlib import Path

from pycsamt.core import _registry as reg


def test_public_api_all():
    expected = {
        "RegistryError",
        "Record",
        "Manifest",
        "ManifestStore",
        "FileManifestStore",
        "Registry",
        "guess_kind",
    }
    assert set(reg.__all__) == expected


def test_manifest_roundtrip(tmp_path: Path):
    store = reg.FileManifestStore()
    p = tmp_path / "manifest.json"
    man = reg.Manifest(root=str(tmp_path))
    assert man.records == {}
    store.save(man, p)
    back = store.load(p)
    assert back.root == str(tmp_path)
    assert back.records == {}


def test_registry_add_path_and_reload(tmp_path: Path):
    root = tmp_path / "session"
    root.mkdir()
    data = root / "file.dat"
    data.write_bytes(b"hello world\n")

    r = reg.Registry(root)
    rec = r.add_path(
        data,
        kind="edi",
        fmt="edi",
        dataid="S001",
        station_id="1",
        tags=["raw", "import"],
        meta={"a": 1},
        with_hash=True,
    )
    assert rec.rid in r.manifest.records
    assert rec.path and Path(rec.path).exists()
    assert rec.checksum and len(rec.checksum) == 64

    # reload from disk
    r2 = reg.Registry(root)
    assert rec.rid in r2.manifest.records
    got = r2.get(rec.rid)
    assert got.tags == ["raw", "import"]
    assert got.meta["a"] == 1


def test_guess_kind_heuristics():
    # craft dynamic classes with module hints
    Avg = type("Avg", (), {})
    Avg.__module__ = "pycsamt.zonge.avg"
    Jf = type("Jf", (), {})
    Jf.__module__ = "pycsamt.jones.j"
    Jc = type("Jc", (), {})
    Jc.__module__ = "pycsamt.jones.collection"
    Edi = type("Edi", (), {})
    Edi.__module__ = "pycsamt.seg.edi"

    assert reg.guess_kind(Avg()) == "avg"
    assert reg.guess_kind(Jf()) == "j"
    assert reg.guess_kind(Jc()) in {"j_col", "j"}
    assert reg.guess_kind(Edi()) == "edi"
    assert reg.guess_kind([1, 2]) == "list"
    assert reg.guess_kind(object()) == "meta"


def test_registry_add_obj_and_query(tmp_path: Path):
    root = tmp_path / "work"
    root.mkdir()
    r = reg.Registry(root)

    class EdiLike:
        __module__ = "pycsamt.seg.edi"
        def __init__(self):
            self.station = "S001"
            self.station_id = 1

    rec = r.add_obj(EdiLike(), tags=["prod"], meta={"a": 2})
    assert rec.kind == "edi"
    assert rec.dataid == "S001"
    found = r.find(tag="prod")
    assert len(found) == 1 and found[0].rid == rec.rid

    # update and remove
    updated = r.update_meta(rec.rid, note="ok")
    assert updated.meta["note"] == "ok"
    r.remove(rec.rid)
    assert r.find(tag="prod") == []
