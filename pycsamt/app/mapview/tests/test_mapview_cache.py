# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Unit tests for pycsamt.app.mapview.cache — session-scoped view storage."""

from __future__ import annotations


class TestViewStorage:
    def test_roundtrip(self):
        from pycsamt.app.mapview.cache import get_view, set_view

        set_view("sess-1", {"stations": 3})
        assert get_view("sess-1") == {"stations": 3}

    def test_missing_session_returns_none(self):
        from pycsamt.app.mapview.cache import get_view

        assert get_view("no-such-session") is None

    def test_empty_session_id_returns_none(self):
        from pycsamt.app.mapview.cache import get_view

        assert get_view("") is None
        assert get_view(None) is None

    def test_set_view_ignores_empty_session_id(self):
        from pycsamt.app.mapview.cache import get_view, set_view

        set_view("", {"stations": 99})
        assert get_view("") is None

    def test_sessions_are_namespaced_from_web_cache(self):
        """mapview and web share a backend but must not collide on keys."""
        from pycsamt.app.mapview.cache import get_view, set_view

        set_view("shared-id", {"who": "mapview"})
        try:
            from pycsamt.app.web.cache import cache_get

            assert cache_get("shared-id") != {"who": "mapview"}
        except ImportError:
            pass
        assert get_view("shared-id") == {"who": "mapview"}


class TestSeedHandoff:
    def test_take_seed_returns_and_clears(self):
        from pycsamt.app.mapview.cache import set_seed, take_seed

        set_seed({"seeded": True})
        assert take_seed() == {"seeded": True}
        assert take_seed() is None

    def test_take_seed_default_is_none(self):
        from pycsamt.app.mapview.cache import take_seed

        assert take_seed() is None
