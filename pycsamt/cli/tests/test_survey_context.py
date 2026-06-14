# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Unit tests for SurveyContext and resolve_survey cache logic."""

from __future__ import annotations

import json
import pickle
from pathlib import Path
from unittest.mock import patch

import pytest

from pycsamt.cli.survey import (
    SurveyContext,
    _cache_key,
    _cache_is_valid,
    _cache_meta,
    _cache_pkl,
    _write_cache,
    _read_cache,
    resolve_survey,
    set_survey,
    survey_summary,
)


# ---------------------------------------------------------------------------
# Picklable fake Sites (MagicMock is not picklable across subprocesses)
# ---------------------------------------------------------------------------

class _FakeSite:
    def __init__(self, name: str) -> None:
        self.name = name


class _FakeSites:
    """Minimal picklable Sites-like object for cache tests."""

    def __init__(self, names: list[str]) -> None:
        self._items = [_FakeSite(n) for n in names]

    def __len__(self) -> int:
        return len(self._items)

    def __iter__(self):
        return iter(self._items)


def _make_sites(n: int = 3, prefix: str = "S") -> _FakeSites:
    return _FakeSites([f"{prefix}{i:02d}" for i in range(1, n + 1)])


# ---------------------------------------------------------------------------
# _cache_key
# ---------------------------------------------------------------------------

class TestCacheKey:
    def test_same_path_same_key(self, tmp_path: Path) -> None:
        p = tmp_path / "survey"
        p.mkdir()
        assert _cache_key(p) == _cache_key(p)

    def test_different_paths_different_keys(self, tmp_path: Path) -> None:
        a = tmp_path / "a"
        b = tmp_path / "b"
        a.mkdir()
        b.mkdir()
        assert _cache_key(a) != _cache_key(b)

    def test_key_is_12_chars(self, tmp_path: Path) -> None:
        assert len(_cache_key(tmp_path)) == 12


# ---------------------------------------------------------------------------
# SurveyContext
# ---------------------------------------------------------------------------

class TestSurveyContext:
    def test_load_returns_none_when_missing(
        self, isolated_home: Path
    ) -> None:
        assert SurveyContext.load() is None

    def test_save_and_load_roundtrip(
        self, isolated_home: Path, tmp_path: Path
    ) -> None:
        survey_path = tmp_path / "edis"
        survey_path.mkdir()
        sites = _make_sites(2)

        ctx = SurveyContext.save(survey_path, sites)
        assert ctx.n_stations == 2
        assert ctx.station_names == ["S01", "S02"]

        loaded = SurveyContext.load()
        assert loaded is not None
        assert loaded.survey_path == survey_path.resolve()
        assert loaded.n_stations == 2

    def test_clear_removes_context(
        self, isolated_home: Path, tmp_path: Path
    ) -> None:
        survey_path = tmp_path / "edis"
        survey_path.mkdir()
        SurveyContext.save(survey_path, _make_sites(1))

        assert SurveyContext.load() is not None
        SurveyContext.clear()
        assert SurveyContext.load() is None


# ---------------------------------------------------------------------------
# Cache helpers
# ---------------------------------------------------------------------------

class TestCacheHelpers:
    def test_write_and_read_roundtrip(
        self, isolated_home: Path, tmp_path: Path
    ) -> None:
        path = tmp_path / "survey"
        path.mkdir()
        (path / "dummy.edi").write_text("dummy")

        sites = _make_sites(5)
        key = _cache_key(path)
        _write_cache(key, path, sites)

        assert _cache_pkl(key).exists()
        assert _cache_meta(key).exists()
        loaded = _read_cache(key)
        assert loaded is not None
        assert len(loaded) == 5

    def test_cache_is_valid_after_write(
        self, isolated_home: Path, tmp_path: Path
    ) -> None:
        path = tmp_path / "survey"
        path.mkdir()
        (path / "dummy.edi").write_text("dummy")
        key = _cache_key(path)
        _write_cache(key, path, _make_sites(2))
        assert _cache_is_valid(key, path)

    def test_cache_invalid_when_source_newer(
        self, isolated_home: Path, tmp_path: Path
    ) -> None:
        path = tmp_path / "survey"
        path.mkdir()
        (path / "dummy.edi").write_text("dummy")

        key = _cache_key(path)
        _write_cache(key, path, _make_sites(2))

        meta_path = _cache_meta(key)
        meta = json.loads(meta_path.read_text())
        meta["cached_at_mtime"] = meta["cached_at_mtime"] - 10.0
        meta_path.write_text(json.dumps(meta))

        assert not _cache_is_valid(key, path)

    def test_read_cache_returns_none_on_corrupt_pickle(
        self, isolated_home: Path, tmp_path: Path
    ) -> None:
        path = tmp_path / "survey"
        path.mkdir()
        (path / "dummy.edi").write_text("dummy")
        key = _cache_key(path)
        _write_cache(key, path, _make_sites(2))
        _cache_pkl(key).write_bytes(b"not a pickle")
        assert _read_cache(key) is None


# ---------------------------------------------------------------------------
# resolve_survey
# ---------------------------------------------------------------------------

class TestResolveSurvey:
    def test_raises_usage_error_with_no_context(
        self, isolated_home: Path
    ) -> None:
        import click
        with pytest.raises(click.UsageError, match="No active survey"):
            resolve_survey(None)

    def test_resolves_explicit_path_bypassing_context(
        self, isolated_home: Path, tmp_path: Path
    ) -> None:
        path = tmp_path / "survey"
        path.mkdir()
        (path / "dummy.edi").write_text("dummy")
        fake_sites = _make_sites(3)

        with patch("pycsamt.cli.survey._build_sites", return_value=fake_sites):
            sites = resolve_survey(path)
        assert len(sites) == 3

    def test_uses_active_context_when_no_explicit_path(
        self, isolated_home: Path, tmp_path: Path
    ) -> None:
        path = tmp_path / "survey"
        path.mkdir()
        (path / "dummy.edi").write_text("dummy")
        fake_sites = _make_sites(4)

        with patch("pycsamt.cli.survey._build_sites", return_value=fake_sites):
            set_survey(path)
            sites = resolve_survey(None)
        assert len(sites) == 4

    def test_fresh_flag_bypasses_valid_cache(
        self, isolated_home: Path, tmp_path: Path
    ) -> None:
        path = tmp_path / "survey"
        path.mkdir()
        (path / "dummy.edi").write_text("dummy")

        call_count = {"n": 0}

        def counting_build(p, verbose=0):
            call_count["n"] += 1
            return _make_sites(2)

        with patch("pycsamt.cli.survey._build_sites", side_effect=counting_build):
            resolve_survey(path)              # build #1
            resolve_survey(path)              # from cache — no build
            resolve_survey(path, fresh=True)  # forced build #2

        assert call_count["n"] == 2


# ---------------------------------------------------------------------------
# survey_summary
# ---------------------------------------------------------------------------

class TestSurveySummary:
    def test_returns_none_when_no_context(
        self, isolated_home: Path
    ) -> None:
        assert survey_summary() is None

    def test_returns_dict_after_set(
        self, isolated_home: Path, tmp_path: Path
    ) -> None:
        path = tmp_path / "survey"
        path.mkdir()
        (path / "dummy.edi").write_text("dummy")
        fake_sites = _make_sites(7)

        with patch("pycsamt.cli.survey._build_sites", return_value=fake_sites):
            set_survey(path)

        s = survey_summary()
        assert s is not None
        assert "survey_path" in s
        assert s["n_stations"] == 7
        assert "cache_key" in s
        assert "cache_valid" in s
