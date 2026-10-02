# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Station Statistics card: measured quality (not n_freq // 10), response
preview and error-coloured frequency coverage — on real bundled data."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pycsamt.app.desktop.controllers.station_stats import (
    grade_of,
    station_stats,
    survey_medians,
)

ROOT = Path(__file__).resolve().parents[4]
KAP = ROOT / "data" / "MT" / "kap03lmt_edis"
CSAMT = ROOT / "data" / "CSAMT"


def _sites(folder):
    if not folder.is_dir():
        pytest.skip(f"{folder} not present")
    from pycsamt.emtools import ensure_sites

    return ensure_sites(str(folder))


def test_grades():
    assert [grade_of(v) for v in (0.9, 0.75, 0.6, 0.45, 0.1)] == \
        list("ABCDE")
    assert grade_of(float("nan")) == "–"


def test_station_stats_match_the_library_response():
    s = next(iter(_sites(KAP)))
    st = station_stats(s)
    order = np.argsort(np.asarray(s.freq))
    np.testing.assert_allclose(st.rho_xy,
                               np.asarray(s.rho)[order, 0, 1], rtol=1e-6)
    np.testing.assert_allclose(st.phi_xy,
                               np.asarray(s.phase)[order, 0, 1], atol=1e-6)
    assert 0 <= st.completeness <= 1 and st.snr > 0
    assert 0 <= st.skew_ok <= 1 and st.grade in "ABCDE"
    assert np.all(np.diff(st.freq) > 0)


def test_scalar_csamt_has_no_skew_test():
    st = station_stats(next(iter(_sites(CSAMT))))
    assert st.scalar and np.isnan(st.skew_ok)
    assert st.grade in "ABCDE"


def test_quality_is_not_a_frequency_count():
    """A station with many noisy frequencies must not outrank a clean
    one just by having more of them (the old n_freq // 10 dots)."""
    from types import SimpleNamespace

    f = np.logspace(-2, 3, 60)
    z = np.tile(np.array([[0.0, 10 + 10j], [-10 - 10j, 0.0]]), (60, 1, 1))
    z[:, 0, 0] = z[:, 1, 1] = 0.01
    noisy = SimpleNamespace(name="n", freq=f, z=z, z_err=np.abs(z) * 2)
    clean = SimpleNamespace(name="c", freq=f[:10], z=z[:10],
                            z_err=np.abs(z[:10]) * 0.01)
    assert station_stats(clean).score > station_stats(noisy).score


def test_survey_medians():
    m = survey_medians(_sites(CSAMT))
    assert m["n_stations"] == 10 and m["fmin"] < m["fmax"]


def test_card_populates_and_animates(qapp):
    pytest.importorskip("PySide6")
    from pycsamt.app.desktop.widgets.station_detail import StationDetailCard

    sites = _sites(KAP)
    card = StationDetailCard()
    name = next(iter(sites)).name
    card.update_station(name, sites)
    assert card._badge.text() in "ABCDE"
    assert card._metric_rows[0][0].text().endswith("%")
    assert card._response._anim.state() == card._response._anim.State.Running
    assert len(card._frequency_graph._values) == card.stats.freq.size
    card.set_animations(False)
    card.update_station(name, sites)
    assert card._response.progress == 1.0
    card.clear()
    assert card._badge.text() == "–" and not card._btn_map.isEnabled()
