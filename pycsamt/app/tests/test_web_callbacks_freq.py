# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for pycsamt.app.web.callbacks.freq (frequency/period slider init).

Real data
---------
data/AMT/WILLY_DATA/L18PLT/ — real EDI-derived Sites via
``pycsamt.app.web.cache.cache_set`` so ``_log10_range`` runs against a
real frequency array, not a mock.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from pycsamt.app.web.callbacks.freq import _log10_range
from pycsamt.app.web.layout import IDs

_ROOT = Path(__file__).parents[3]  # pycsamt/
_WILLY_L18 = _ROOT / "data" / "AMT" / "WILLY_DATA" / "L18PLT"
_HAS_WILLY = _WILLY_L18.exists() and any(_WILLY_L18.glob("*.edi"))


@pytest.fixture(scope="session")
def willy_sites():
    pytest.importorskip("pycsamt.emtools")
    if not _HAS_WILLY:
        pytest.skip("WILLY L18PLT data not available")
    from pycsamt.emtools import ensure_sites

    return ensure_sites(str(_WILLY_L18))


def _unwrap(entry):
    fn = entry["callback"]
    return getattr(fn, "__wrapped__", fn)


def _cb(web_app, output_id_prop):
    return _unwrap(web_app.callback_map[output_id_prop])


def _init_slider(web_app):
    key = next(
        k
        for k in web_app.callback_map
        if f"{IDs.FREQ_SLIDER}.min" in k and f"{IDs.FREQ_SLIDER}.disabled" in k
    )
    return _cb(web_app, key)


# ── _log10_range ─────────────────────────────────────────────────────────────


@pytest.fixture
def bypass_ensure_sites(monkeypatch):
    """
    ``_log10_range`` calls ``ensure_sites(sites)`` before iterating, and
    the real ``ensure_sites`` silently drops anything that doesn't look
    like a genuine EDI/Site object -- a plain fake ``_Site`` stand-in
    passed straight through it comes back as an *empty* Sites collection,
    which would make every "synthetic site" test below tautologically
    hit the same "no freqs found" fallback as an empty list, regardless
    of what the fake site's ``.freq``/``.freqs`` actually contain. Patch
    it to identity so the fakes reach the per-site loop unmodified.
    """
    import pycsamt.emtools._core as core_mod

    monkeypatch.setattr(core_mod, "ensure_sites", lambda sites: sites)


class TestLog10Range:
    def test_real_sites_hits_the_getattr_or_crash_bug(self, willy_sites):
        """
        Real bug: ``getattr(s, "freq", None) or getattr(s, "freqs", None)``
        crashes with numpy's "truth value of an array is ambiguous"
        ValueError whenever ``s.freq`` is already a populated (2+ element)
        array -- which is every real EDI site, since ``or`` evaluates
        ``bool()`` on a truthy left operand before short-circuiting. The
        per-site ``except Exception: continue`` swallows it silently, so
        *every* station is skipped and the function always falls back to
        the hardcoded (-4.0, 4.0) default -- the frequency/period slider
        never reflects a real survey's actual range. Confirmed directly:
        ``np.array([1.,2.,3.]) or None`` raises the same ValueError; the
        only reason a synthetic ``freq=None`` site doesn't trip it is that
        ``or`` never calls ``bool()`` on the right operand once the left
        is falsy. Not fixed here, per instructions. (Real ``ensure_sites``
        does accept real Sites through untouched, so this test needs no
        bypass -- unlike the synthetic-object tests below.)
        """
        lo, hi = _log10_range(willy_sites)
        assert (lo, hi) == (-4.0, 4.0)

    def test_synthetic_site_with_populated_freq_array_also_crashes(
        self, bypass_ensure_sites
    ):
        import numpy as np

        class _Site:
            freq = np.array([1.0, 10.0, 100.0])

        assert _log10_range([_Site()]) == (-4.0, 4.0)

    def test_plain_list_freq_avoids_the_bug_and_computes_real_range(
        self, bypass_ensure_sites
    ):
        """
        The crash is specific to numpy arrays: ``bool()`` on a plain
        Python list/tuple with 2+ elements is well-defined, so a site
        whose ``.freq`` happens to be a *list* rather than an ndarray
        sidesteps the ``or``-truthiness bug and reaches the real
        min/max computation -- this is the only way to exercise that
        code path given real EDI Site objects always expose ``.freq``
        as an ndarray.
        """

        class _Site:
            freq = [1.0, 10.0, 100.0]

        lo, hi = _log10_range([_Site()])
        assert (lo, hi) != (-4.0, 4.0)
        assert lo < hi

    def test_empty_sites_returns_default(self):
        assert _log10_range([]) == (-4.0, 4.0)

    def test_site_missing_freq_attr_skipped(self, bypass_ensure_sites):
        class _NoFreq:
            pass

        assert _log10_range([_NoFreq()]) == (-4.0, 4.0)

    def test_site_freq_iteration_exception_skipped(self, bypass_ensure_sites):
        class _BadFreq:
            @property
            def freq(self):
                raise RuntimeError("boom")

        assert _log10_range([_BadFreq()]) == (-4.0, 4.0)

    def test_ensure_sites_exception_returns_default(self, monkeypatch):
        import pycsamt.emtools._core as core_mod

        def _boom(sites):
            raise RuntimeError("cannot ensure")

        monkeypatch.setattr(core_mod, "ensure_sites", _boom)
        assert _log10_range(object()) == (-4.0, 4.0)

    def test_uses_freqs_alias_when_freq_missing(self, bypass_ensure_sites):
        """``.freq`` is None (falsy) so ``or`` short-circuits to ``.freqs``
        *without* calling bool() on it -- the one way a populated ndarray
        can flow through this line without tripping the bug above."""
        import numpy as np

        class _Site:
            freq = None
            freqs = np.array([1.0, 10.0, 100.0])

        lo, hi = _log10_range([_Site()])
        assert (lo, hi) != (-4.0, 4.0)
        assert lo < hi


# ── init_slider callback ─────────────────────────────────────────────────────


class TestInitSlider:
    def test_no_store_data_disables_slider(self, web_app):
        fn = _init_slider(web_app)
        result = fn(None, "session-1")
        assert result == (-4.0, 4.0, [-4.0, 4.0], True)

    def test_no_session_id_disables_slider(self, web_app):
        fn = _init_slider(web_app)
        result = fn({"n_stations": 1}, None)
        assert result == (-4.0, 4.0, [-4.0, 4.0], True)

    def test_cache_miss_disables_slider(self, web_app):
        fn = _init_slider(web_app)
        result = fn({"n_stations": 1}, "session-not-cached")
        assert result == (-4.0, 4.0, [-4.0, 4.0], True)

    def test_cache_hit_enables_slider_but_range_is_the_broken_default(
        self, web_app, willy_sites
    ):
        """Slider does get enabled on a real cache hit, but -- per the
        _log10_range bug above -- always with the hardcoded default
        range rather than the real survey's period range."""
        from pycsamt.app.web.cache import cache_set

        session_id = "test-freq-session"
        cache_set(session_id, willy_sites)
        fn = _init_slider(web_app)
        lo, hi, value, disabled = fn({"n_stations": 5}, session_id)
        assert disabled is False
        assert (lo, hi) == (-4.0, 4.0)
        assert value == [lo, hi]
