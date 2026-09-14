# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Supplementary coverage tests for :mod:`pycsamt.agents.plotting`.

Complements ``test_plot_agent.py`` with the ``_period_range`` non-numeric
exception guard and ``_has_tipper``'s per-station exception-continue and
success-return branches.
"""

from __future__ import annotations

import pytest

from pycsamt.agents.plotting import _has_tipper, _period_range

pytestmark = pytest.mark.filterwarnings("ignore::RuntimeWarning")


def test_period_range_non_numeric_returns_none():
    assert _period_range({"period_min": "abc", "period_max": "10"}) is None


def test_has_tipper_skips_exception_and_finds_later_station(monkeypatch):
    import pycsamt.emtools._core as core

    class _BoomSite:
        pass

    class _GoodSite:
        pass

    def _fake_get_t_block(ed):
        if isinstance(ed, _BoomSite):
            raise RuntimeError("boom")
        return None, [1, 2, 3], [0.1, 0.2, 0.3]

    monkeypatch.setattr(core, "_get_t_block", _fake_get_t_block)
    assert _has_tipper([_BoomSite(), _GoodSite()]) is True


def test_has_tipper_false_when_no_valid_tipper(monkeypatch):
    import pycsamt.emtools._core as core

    monkeypatch.setattr(core, "_get_t_block", lambda ed: (None, None, None))
    assert _has_tipper([object()]) is False
