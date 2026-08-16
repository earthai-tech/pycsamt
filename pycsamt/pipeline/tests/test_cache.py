"""Tests for pycsamt.pipeline._cache — fingerprint_sites, chain_key, StepCache.

Public entry point: :mod:`pycsamt.pipeline` (``StepCache``, ``fingerprint_sites``)
Implementation:     :mod:`pycsamt.pipeline._cache`
"""

from __future__ import annotations

from pathlib import Path

import pytest

from pycsamt.pipeline._cache import (
    _MISS,
    DIAGNOSTIC_OK,
    StepCache,
    chain_key,
    fingerprint_sites,
)

# ─────────────────────────────────────────────────────────────────────────────
# fingerprint_sites
# ─────────────────────────────────────────────────────────────────────────────


class _CountingSites:
    """Same test double used throughout pycsamt/pipeline/tests/test_pipeline.py.

    Deliberately NOT iterable in the Site-yielding sense — this is exactly
    the shape that must fall back to the pickle-based fingerprint path.
    """

    def __init__(self, n: int = 3, count: int = 0):
        self._n = n
        self.count = count

    def __len__(self) -> int:
        return self._n


class TestFingerprintSitesFallback:
    def test_identical_content_same_fingerprint(self):
        a = _CountingSites(n=5, count=2)
        b = _CountingSites(n=5, count=2)
        assert fingerprint_sites(a) == fingerprint_sites(b)

    def test_different_content_different_fingerprint(self):
        a = _CountingSites(n=5, count=2)
        b = _CountingSites(n=5, count=3)
        assert fingerprint_sites(a) != fingerprint_sites(b)

    def test_is_a_hex_digest(self):
        fp = fingerprint_sites(_CountingSites())
        assert len(fp) == 64
        assert all(c in "0123456789abcdef" for c in fp)

    def test_empty_countingsites_does_not_raise(self):
        # len(sites) == 0 and not iterable-as-Sites: must still fall back
        # cleanly rather than raising.
        fp = fingerprint_sites(_CountingSites(n=0, count=0))
        assert len(fp) == 64


class _FakeSite:
    def __init__(self, name, freq, z):
        self.name = name
        self.freq = freq
        self.z = z


class _FakeSites:
    """Minimal stand-in for the real Sites shape: iterable of Site-like objects."""

    def __init__(self, stations):
        self._stations = list(stations)

    def __iter__(self):
        return iter(self._stations)

    def __len__(self):
        return len(self._stations)


class TestFingerprintSitesRealShape:
    def test_identical_arrays_same_fingerprint(self):
        import numpy as np

        a = _FakeSites([_FakeSite("S01", np.array([1.0, 2.0]), np.array([1 + 1j]))])
        b = _FakeSites([_FakeSite("S01", np.array([1.0, 2.0]), np.array([1 + 1j]))])
        assert fingerprint_sites(a) == fingerprint_sites(b)

    def test_different_freq_different_fingerprint(self):
        import numpy as np

        a = _FakeSites([_FakeSite("S01", np.array([1.0, 2.0]), np.array([1 + 1j]))])
        b = _FakeSites([_FakeSite("S01", np.array([1.0, 3.0]), np.array([1 + 1j]))])
        assert fingerprint_sites(a) != fingerprint_sites(b)

    def test_station_order_is_significant(self):
        import numpy as np

        s1 = _FakeSite("S01", np.array([1.0]), np.array([1j]))
        s2 = _FakeSite("S02", np.array([2.0]), np.array([2j]))
        forward = _FakeSites([s1, s2])
        backward = _FakeSites([s2, s1])
        assert fingerprint_sites(forward) != fingerprint_sites(backward)

    def test_real_edi_data_fingerprint_matches_reload(self):
        """Real WILLY EDI data, loaded twice, must fingerprint identically —
        this is the property the whole cache correctness story depends on."""
        from pycsamt.emtools._core import ensure_sites

        data_dir = (
            Path(__file__).resolve().parents[3]
            / "data"
            / "AMT"
            / "WILLY_DATA"
            / "L22PLT"
        )
        if not data_dir.exists():
            pytest.skip(f"WILLY EDI data not found at {data_dir}")
        paths = sorted(data_dir.glob("*.edi"))[:3]
        sites_a = ensure_sites([str(p) for p in paths])
        sites_b = ensure_sites([str(p) for p in paths])
        assert fingerprint_sites(sites_a) == fingerprint_sites(sites_b)


# ─────────────────────────────────────────────────────────────────────────────
# chain_key
# ─────────────────────────────────────────────────────────────────────────────


class TestChainKey:
    def test_deterministic(self):
        k1 = chain_key("fp1", "NR001", {"mains_hz": 50})
        k2 = chain_key("fp1", "NR001", {"mains_hz": 50})
        assert k1 == k2

    def test_sensitive_to_upstream_fingerprint(self):
        k1 = chain_key("fp1", "NR001", {"mains_hz": 50})
        k2 = chain_key("fp2", "NR001", {"mains_hz": 50})
        assert k1 != k2

    def test_sensitive_to_code(self):
        k1 = chain_key("fp1", "NR001", {})
        k2 = chain_key("fp1", "NR002", {})
        assert k1 != k2

    def test_sensitive_to_params(self):
        k1 = chain_key("fp1", "NR001", {"mains_hz": 50})
        k2 = chain_key("fp1", "NR001", {"mains_hz": 60})
        assert k1 != k2

    def test_param_key_order_does_not_matter(self):
        k1 = chain_key("fp1", "NR001", {"a": 1, "b": 2})
        k2 = chain_key("fp1", "NR001", {"b": 2, "a": 1})
        assert k1 == k2

    def test_unhashable_param_value_does_not_raise(self):
        # default=str in the JSON dump: never crashes on odd param types.
        chain_key("fp1", "NR001", {"fn": lambda x: x})  # must not raise


# ─────────────────────────────────────────────────────────────────────────────
# StepCache
# ─────────────────────────────────────────────────────────────────────────────


class TestStepCache:
    def test_miss_on_empty_cache(self, tmp_path):
        cache = StepCache(tmp_path)
        assert cache.get("deadbeef") is _MISS

    def test_put_then_get_roundtrips(self, tmp_path):
        cache = StepCache(tmp_path)
        sites = _CountingSites(n=5, count=3)
        cache.put("abc123", sites)
        loaded = cache.get("abc123")
        assert loaded is not _MISS
        assert loaded.count == 3
        assert loaded._n == 5

    def test_diagnostic_ok_marker_roundtrips(self, tmp_path):
        cache = StepCache(tmp_path)
        cache.put("diagkey", DIAGNOSTIC_OK)
        assert cache.get("diagkey") is DIAGNOSTIC_OK

    def test_sharded_layout(self, tmp_path):
        cache = StepCache(tmp_path)
        key = "ab" + "0" * 62
        cache.put(key, "value")
        assert (tmp_path / "ab" / f"{key}.joblib").exists()

    def test_no_leftover_temp_files_after_put(self, tmp_path):
        cache = StepCache(tmp_path)
        key = "cd" + "1" * 62
        cache.put(key, "value")
        shard = tmp_path / "cd"
        leftovers = [p for p in shard.iterdir() if p.name.startswith(".")]
        assert leftovers == []

    def test_corrupt_entry_is_a_miss_not_a_crash(self, tmp_path):
        cache = StepCache(tmp_path)
        key = "ef" + "2" * 62
        shard = tmp_path / "ef"
        shard.mkdir(parents=True)
        (shard / f"{key}.joblib").write_bytes(b"not a valid joblib file")

        with pytest.warns(UserWarning, match="unreadable"):
            result = cache.get(key)
        assert result is _MISS

    def test_clear_removes_everything(self, tmp_path):
        cache = StepCache(tmp_path)
        cache.put("abc123", "value")
        assert tmp_path.exists()
        cache.clear()
        assert not tmp_path.exists()

    def test_default_root_is_under_home_pycsamt(self):
        cache = StepCache()
        assert cache.root == Path.home() / ".pycsamt" / "pipeline_cache"

    def test_overwrite_existing_key(self, tmp_path):
        cache = StepCache(tmp_path)
        cache.put("abc123", "first")
        cache.put("abc123", "second")
        assert cache.get("abc123") == "second"
