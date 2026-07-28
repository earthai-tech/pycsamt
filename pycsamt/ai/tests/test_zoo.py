# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.ai._zoo — model registry and download utilities."""

from __future__ import annotations

import urllib.request

import pytest

from pycsamt.ai._zoo import (
    _MODEL_ZOO,
    _get_cache_dir,
    _md5_file,
    download_checkpoint,
    get_pretrained_info,
    list_pretrained,
)


def test_list_pretrained_returns_name_to_description_mapping():
    models = list_pretrained()
    assert isinstance(models, dict)
    assert set(models) == set(_MODEL_ZOO)
    for name, desc in models.items():
        assert desc == _MODEL_ZOO[name]["description"]


def test_get_pretrained_info_known_model():
    info = get_pretrained_info("mt1d-resnet-5layer-v1")
    assert info["arch"] == "resnet"
    assert info["n_layers"] == 5
    assert info["solver"] == "mt1d"
    # returned dict must be a copy, not the registry entry itself
    info["arch"] = "mutated"
    assert _MODEL_ZOO["mt1d-resnet-5layer-v1"]["arch"] == "resnet"


def test_get_pretrained_info_unknown_model_raises_key_error():
    with pytest.raises(KeyError, match="Unknown model"):
        get_pretrained_info("does-not-exist-v99")


def test_get_cache_dir_explicit_override(tmp_path):
    d = _get_cache_dir(str(tmp_path / "explicit"))
    assert d == tmp_path / "explicit"


def test_get_cache_dir_env_var(monkeypatch, tmp_path):
    monkeypatch.delenv("PYCSAMT_MODEL_CACHE", raising=False)
    env_dir = tmp_path / "from_env"
    monkeypatch.setenv("PYCSAMT_MODEL_CACHE", str(env_dir))
    d = _get_cache_dir(None)
    assert d == env_dir


def test_get_cache_dir_default_home(monkeypatch):
    monkeypatch.delenv("PYCSAMT_MODEL_CACHE", raising=False)
    d = _get_cache_dir(None)
    assert str(d).endswith(".pycsamt/model_zoo") or str(d).endswith(
        ".pycsamt\\model_zoo"
    )


def test_md5_file(tmp_path):
    p = tmp_path / "sample.bin"
    p.write_bytes(b"hello world")
    import hashlib

    expected = hashlib.md5(b"hello world").hexdigest()
    assert _md5_file(p) == expected


def test_download_checkpoint_unknown_model_raises_key_error(tmp_path):
    with pytest.raises(KeyError):
        download_checkpoint("nope-v0", cache_dir=str(tmp_path))


def test_download_checkpoint_uses_cache_when_present(tmp_path, capsys):
    info = get_pretrained_info("mt1d-cnn-5layer-v1")
    fname = info["url"].split("/")[-1]
    cache = tmp_path / "cache"
    cache.mkdir()
    (cache / fname).write_bytes(b"fake-weights")

    path = download_checkpoint("mt1d-cnn-5layer-v1", cache_dir=str(cache), verbose=True)

    assert path == cache / fname
    assert path.read_bytes() == b"fake-weights"
    out = capsys.readouterr().out
    assert "Using cached" in out


def test_download_checkpoint_force_redownload(monkeypatch, tmp_path):
    info = get_pretrained_info("mt1d-cnn-5layer-v1")
    fname = info["url"].split("/")[-1]
    cache = tmp_path / "cache"
    cache.mkdir()
    (cache / fname).write_bytes(b"stale")

    def _fake_urlretrieve(url, fpath):
        with open(fpath, "wb") as f:
            f.write(b"fresh-weights")

    monkeypatch.setattr(urllib.request, "urlretrieve", _fake_urlretrieve)

    path = download_checkpoint(
        "mt1d-cnn-5layer-v1", cache_dir=str(cache), force=True, verbose=False
    )
    assert path.read_bytes() == b"fresh-weights"


def test_download_checkpoint_success_no_md5(monkeypatch, tmp_path):
    def _fake_urlretrieve(url, fpath):
        with open(fpath, "wb") as f:
            f.write(b"weights-blob")

    monkeypatch.setattr(urllib.request, "urlretrieve", _fake_urlretrieve)

    path = download_checkpoint(
        "tem1d-fcn-5layer-v1", cache_dir=str(tmp_path), verbose=True
    )
    assert path.exists()
    assert path.read_bytes() == b"weights-blob"


def test_download_checkpoint_download_failure_raises_runtime_error(
    monkeypatch, tmp_path
):
    def _fake_urlretrieve(url, fpath):
        raise OSError("network unreachable")

    monkeypatch.setattr(urllib.request, "urlretrieve", _fake_urlretrieve)

    with pytest.raises(RuntimeError, match="Failed to download"):
        download_checkpoint(
            "mt1d-resnet-7layer-v1", cache_dir=str(tmp_path), verbose=False
        )
    # the partially-downloaded file must not linger
    info = get_pretrained_info("mt1d-resnet-7layer-v1")
    fname = info["url"].split("/")[-1]
    assert not (tmp_path / fname).exists()


def test_download_checkpoint_md5_mismatch_raises_runtime_error(monkeypatch, tmp_path):
    import pycsamt.ai._zoo as zoo_mod

    monkeypatch.setitem(
        zoo_mod._MODEL_ZOO,
        "mt1d-resnet-7layer-v1",
        {
            **zoo_mod._MODEL_ZOO["mt1d-resnet-7layer-v1"],
            "md5": "deadbeefdeadbeefdeadbeefdeadbeef",
        },
    )

    def _fake_urlretrieve(url, fpath):
        with open(fpath, "wb") as f:
            f.write(b"wrong-content")

    monkeypatch.setattr(urllib.request, "urlretrieve", _fake_urlretrieve)

    with pytest.raises(RuntimeError, match="MD5 mismatch"):
        download_checkpoint(
            "mt1d-resnet-7layer-v1", cache_dir=str(tmp_path), verbose=False
        )
    info = get_pretrained_info("mt1d-resnet-7layer-v1")
    fname = info["url"].split("/")[-1]
    assert not (tmp_path / fname).exists()


if __name__ == "__main__":
    pytest.main([__file__])
