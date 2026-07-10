from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pycsamt.core import registry as reg
from pycsamt.core.base import TFBundle


def test_public_api_all():
    expected = {
        "Packer",
        "register_packer",
        "get_packer",
        "list_packers",
        "pack_to_file",
        "unpack_from_file",
        "RegistryAPI",
    }
    assert set(reg.__all__) == expected


def test_default_bundle_packer_roundtrip(tmp_path: Path):
    b = TFBundle(
        freq=np.array([1.0, 2.0, 3.0]),
        z=np.ones((3, 2, 2), complex),
        station="S001",
        attrs={"note": "ok"},
    )
    out = tmp_path / "bundle.npz"
    reg.pack_to_file(b, out)
    assert out.exists()
    back = reg.unpack_from_file(out)
    assert isinstance(back, TFBundle)
    assert back.station == "S001"
    assert np.asarray(back.freq).size == 3
    assert np.asarray(back.z).shape == (3, 2, 2)
    assert back.attrs.get("note") == "ok"


def test_register_custom_packer_and_list():
    calls = {"pack": 0, "unpack": 0}

    def pk(x):
        calls["pack"] += 1
        return {"kind": "echo", "v": 42}

    def up(d):
        calls["unpack"] += 1
        return d.get("v")

    reg.register_packer("echo", (pk, up))
    got = reg.get_packer("echo")
    assert isinstance(got, tuple) and len(got) == 2

    listing = reg.list_packers()
    assert "echo" in listing and "bundle" in listing


def test_pack_to_file_rejects_unknown_kind(tmp_path: Path):
    with pytest.raises(ValueError):
        reg.pack_to_file(TFBundle(), tmp_path / "x.npz", kind="nope")


def test_registry_api_add_object_pack_and_materialize(tmp_path: Path):
    api = reg.RegistryAPI(tmp_path)

    b = TFBundle(
        freq=np.array([1.0]),
        z=np.zeros((1, 2, 2)),
        station="S002",
    )
    rec = api.add_object(b, pack=True)
    assert rec.kind == "bundle"
    assert rec.path and Path(rec.path).exists()

    back = api.materialize(rec.rid)
    assert isinstance(back, TFBundle)
    assert back.station == "S002"


def test_registry_api_add_file_and_reload(tmp_path: Path):
    root = tmp_path
    data = root / "hello.txt"
    data.write_text("hi\n", encoding="utf-8")

    api = reg.RegistryAPI(root)
    rec = api.add_file(data, kind="meta", tags=["t1"])
    assert rec.rid in api.low.manifest.records
    rid = rec.rid

    api2 = reg.RegistryAPI(root)
    got = api2.get(rid)
    assert got.tags == ["t1"]
    assert Path(got.path).exists()


def test_to_edi_calls_underlying_dispatch(tmp_path: Path, monkeypatch):
    api = reg.RegistryAPI(tmp_path)
    b = TFBundle(freq=np.array([1.0]), z=np.zeros((1, 2, 2)))
    rec = api.add_object(b, pack=True)

    called = {"n": 0}

    def fake_to_edi(obj, **kw):  # noqa: ANN001
        called["n"] += 1
        assert isinstance(obj, TFBundle)
        return {"ok": True}

    monkeypatch.setattr(reg, "to_edi", fake_to_edi)

    out = api.to_edi(rec.rid)
    assert called["n"] == 1
    assert out["ok"] is True
