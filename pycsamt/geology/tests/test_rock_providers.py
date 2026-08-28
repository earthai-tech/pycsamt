# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.geology.rock_providers.

Uses real ``file://`` URLs rather than mocking ``urlopen`` -- these are
genuine end-to-end fetch/cache/fallback runs against the local
filesystem, not simulated network calls.
"""

from __future__ import annotations

import json

import pytest

import pycsamt.geology.rock_providers as providers
from pycsamt.geology.lithology import RockDatabase, RockEntry
from pycsamt.geology.rock_providers import (
    LocalRockPropertyProvider,
    RemoteRockPropertyProvider,
    RockPropertyProvider,
    RockProviderFetchError,
)

# ─────────────────────────────────────────────────────────────────────────────
# LocalRockPropertyProvider
# ─────────────────────────────────────────────────────────────────────────────


def test_local_provider_default_matches_rockdatabase_default():
    entries, metadata = LocalRockPropertyProvider().fetch()
    default_db = RockDatabase.default()
    assert len(entries) == len(default_db)
    assert metadata["origin"] == "default"


def test_local_provider_csv(tmp_path):
    p = tmp_path / "rocks.csv"
    p.write_text("name,rho_min,rho_max\nPeat,1,10\n")
    entries, metadata = LocalRockPropertyProvider(p).fetch()
    assert len(entries) == 1
    assert entries[0].name == "Peat"
    assert metadata["origin"] == "csv"


# ─────────────────────────────────────────────────────────────────────────────
# RemoteRockPropertyProvider — success path
# ─────────────────────────────────────────────────────────────────────────────


def _write_rocks_json(path, rows):
    path.write_text(json.dumps(rows))
    return path.as_uri()


def test_remote_provider_fetches_and_parses(tmp_path):
    src = tmp_path / "rocks.json"
    url = _write_rocks_json(
        src,
        [
            {
                "name": "Test Rock",
                "rho_min": 10,
                "rho_max": 100,
                "source": "unit-test",
            }
        ],
    )
    provider = RemoteRockPropertyProvider(url, cache_dir=tmp_path / "cache")
    entries, metadata = provider.fetch()
    assert len(entries) == 1
    assert isinstance(entries[0], RockEntry)
    assert entries[0].name == "Test Rock"
    assert entries[0].source == "unit-test"
    assert metadata["origin"] == "url"
    assert metadata["cache_hit"] is False


def test_remote_provider_defaults_optional_fields(tmp_path):
    src = tmp_path / "rocks.json"
    url = _write_rocks_json(src, [{"name": "Bare Rock", "rho_min": 1, "rho_max": 2}])
    provider = RemoteRockPropertyProvider(url, cache_dir=tmp_path / "cache")
    entries, _ = provider.fetch()
    assert entries[0].color == "#AAAAAA"
    assert entries[0].description == ""
    assert entries[0].source == ""


def test_remote_provider_missing_required_field_raises_without_fallback(
    tmp_path,
):
    src = tmp_path / "bad.json"
    url = _write_rocks_json(src, [{"name": "Incomplete"}])  # no rho_min/max
    provider = RemoteRockPropertyProvider(
        url, cache_dir=tmp_path / "cache", fallback=False
    )
    with pytest.raises(RockProviderFetchError):
        provider.fetch()


def test_remote_provider_non_list_payload_raises_without_fallback(tmp_path):
    src = tmp_path / "bad.json"
    src.write_text(json.dumps({"not": "a list"}))
    provider = RemoteRockPropertyProvider(
        src.as_uri(), cache_dir=tmp_path / "cache", fallback=False
    )
    with pytest.raises(RockProviderFetchError):
        provider.fetch()


# ─────────────────────────────────────────────────────────────────────────────
# Caching
# ─────────────────────────────────────────────────────────────────────────────


def test_second_fetch_is_a_cache_hit(tmp_path):
    src = tmp_path / "rocks.json"
    url = _write_rocks_json(src, [{"name": "R", "rho_min": 1, "rho_max": 2}])
    cache_dir = tmp_path / "cache"

    _, meta1 = RemoteRockPropertyProvider(url, cache_dir=cache_dir).fetch()
    assert meta1["cache_hit"] is False

    _, meta2 = RemoteRockPropertyProvider(url, cache_dir=cache_dir).fetch()
    assert meta2["cache_hit"] is True
    assert meta2["fetched_at"] == meta1["fetched_at"]


def test_force_bypasses_a_fresh_cache(tmp_path):
    src = tmp_path / "rocks.json"
    url = _write_rocks_json(src, [{"name": "R1", "rho_min": 1, "rho_max": 2}])
    cache_dir = tmp_path / "cache"

    RemoteRockPropertyProvider(url, cache_dir=cache_dir).fetch()

    # Change the source; without force the stale-but-fresh cache would win.
    _write_rocks_json(src, [{"name": "R2", "rho_min": 1, "rho_max": 2}])
    entries, meta = RemoteRockPropertyProvider(
        url, cache_dir=cache_dir, force=True
    ).fetch()
    assert entries[0].name == "R2"
    assert meta["cache_hit"] is False


def test_ttl_expiry_triggers_a_refetch(tmp_path):
    src = tmp_path / "rocks.json"
    url = _write_rocks_json(src, [{"name": "R1", "rho_min": 1, "rho_max": 2}])
    cache_dir = tmp_path / "cache"

    RemoteRockPropertyProvider(
        url, cache_dir=cache_dir, ttl_seconds=100.0
    ).fetch()

    _write_rocks_json(src, [{"name": "R2", "rho_min": 1, "rho_max": 2}])
    # ttl_seconds=0 means "always stale" without needing a real sleep.
    entries, meta = RemoteRockPropertyProvider(
        url, cache_dir=cache_dir, ttl_seconds=0.0
    ).fetch()
    assert entries[0].name == "R2"
    assert meta["cache_hit"] is False


# ─────────────────────────────────────────────────────────────────────────────
# Fallback behaviour
# ─────────────────────────────────────────────────────────────────────────────


def test_unreachable_url_falls_back_to_default(tmp_path):
    provider = RemoteRockPropertyProvider(
        "file:///no/such/path/rocks.json", cache_dir=tmp_path / "cache"
    )
    entries, metadata = provider.fetch()
    assert len(entries) == len(RockDatabase.default())
    assert metadata["origin"] == "default-fallback"
    assert "error" in metadata


def test_unreachable_url_without_fallback_raises(tmp_path):
    provider = RemoteRockPropertyProvider(
        "file:///no/such/path/rocks.json",
        cache_dir=tmp_path / "cache",
        fallback=False,
    )
    with pytest.raises(RockProviderFetchError):
        provider.fetch()


def test_failed_refetch_falls_back_to_stale_cache_before_default(tmp_path):
    src = tmp_path / "rocks.json"
    url = _write_rocks_json(src, [{"name": "Cached", "rho_min": 1, "rho_max": 2}])
    cache_dir = tmp_path / "cache"

    RemoteRockPropertyProvider(url, cache_dir=cache_dir).fetch()
    src.unlink()  # source now unreachable

    entries, metadata = RemoteRockPropertyProvider(
        url, cache_dir=cache_dir, force=True
    ).fetch()
    assert entries[0].name == "Cached"
    assert metadata["origin"] == "url-stale-cache"


# ─────────────────────────────────────────────────────────────────────────────
# RockDatabase.from_url / from_provider integration
# ─────────────────────────────────────────────────────────────────────────────


def test_rockdatabase_from_url(tmp_path):
    src = tmp_path / "rocks.json"
    url = _write_rocks_json(src, [{"name": "Only", "rho_min": 1, "rho_max": 2}])
    db = RockDatabase.from_url(url, cache_dir=tmp_path / "cache")
    assert len(db) == 1
    assert db.entries[0].name == "Only"
    assert db.metadata["origin"] == "url"


def test_rockdatabase_from_provider():
    db = RockDatabase.from_provider(LocalRockPropertyProvider())
    assert len(db) == len(RockDatabase.default())
    assert db.metadata["origin"] == "default"


def test_provider_protocol_is_runtime_checkable():
    assert isinstance(LocalRockPropertyProvider(), RockPropertyProvider)
    assert not isinstance(object(), RockPropertyProvider)
    assert RockPropertyProvider.fetch(LocalRockPropertyProvider()) is None


def test_cache_directory_resolution_precedence(tmp_path, monkeypatch):
    explicit = tmp_path / "explicit"
    monkeypatch.setenv("PYCSAMT_ROCKDB_CACHE", str(tmp_path / "env"))
    configured = RemoteRockPropertyProvider("file:///x", cache_dir=explicit)
    assert configured._resolve_cache_dir() == explicit
    assert (
        RemoteRockPropertyProvider("file:///x")._resolve_cache_dir()
        == tmp_path / "env"
    )

    monkeypatch.delenv("PYCSAMT_ROCKDB_CACHE")
    monkeypatch.setattr(
        providers.Path, "home", classmethod(lambda cls: tmp_path)
    )
    assert (
        RemoteRockPropertyProvider("file:///x")._resolve_cache_dir()
        == tmp_path / ".pycsamt" / "rock_db"
    )


@pytest.mark.parametrize(
    "contents",
    ["not-json", json.dumps({"fetched_at": "yesterday", "entries": []})],
)
def test_malformed_cache_is_ignored(tmp_path, contents):
    cache = tmp_path / "cache.json"
    cache.write_text(contents)
    assert providers._read_cache(cache, max_age_seconds=10) is None


def test_cache_read_error_is_ignored(tmp_path, monkeypatch):
    cache = tmp_path / "cache.json"
    cache.write_text("{}")

    def fail_read(*args, **kwargs):
        raise OSError("unreadable")

    monkeypatch.setattr(providers.Path, "read_text", fail_read)
    assert providers._read_cache(cache, max_age_seconds=10) is None


def test_invalid_json_response_raises_without_fallback(tmp_path):
    src = tmp_path / "invalid.json"
    src.write_text("not-json")
    provider = RemoteRockPropertyProvider(
        src.as_uri(), cache_dir=tmp_path / "cache", fallback=False
    )
    with pytest.raises(RockProviderFetchError, match="not valid JSON"):
        provider.fetch()


def test_non_object_entry_raises_without_fallback(tmp_path):
    src = tmp_path / "invalid-row.json"
    url = _write_rocks_json(src, ["granite"])
    provider = RemoteRockPropertyProvider(
        url, cache_dir=tmp_path / "cache", fallback=False
    )
    with pytest.raises(RockProviderFetchError, match="not a JSON object"):
        provider.fetch()


def test_remote_provider_preserves_all_optional_fields(tmp_path):
    src = tmp_path / "complete.json"
    url = _write_rocks_json(
        src,
        [
            {
                "name": "R",
                "rho_min": "1",
                "rho_max": "2",
                "color": "#123456",
                "description": "desc",
                "code": "7",
                "source": "field",
            }
        ],
    )
    entries, _ = RemoteRockPropertyProvider(
        url, cache_dir=tmp_path / "cache"
    ).fetch()
    assert entries[0] == RockEntry(
        name="R",
        rho_min=1.0,
        rho_max=2.0,
        color="#123456",
        description="desc",
        code=7,
        source="field",
    )
