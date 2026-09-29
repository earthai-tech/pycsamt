"""Survey metadata correctness and responsive dashboard regression tests."""
import numpy as np
import pandas as pd
import pytest

pytest.importorskip("PySide6")
from pycsamt.app.desktop.widgets.survey_overview import (
    SurveyOverviewWidget, _SurveySummary,
)


def test_coordinates_stay_paired_and_invalid_counts_are_unknown():
    summary = _SurveySummary.from_frame(pd.DataFrame({
        "Latitude": [10, None, 91, 20, 0],
        "Longitude": [None, 11, 12, 181, 0],
        "N_freq": [0, 12, -1, np.inf, 2.5],
        "Tipper": ["false", "true", None, False, "unreported"],
    }, index=[7, 7, 9, 10, 11]))
    assert summary.positioned.tolist() == [False, False, False, False, True]
    assert summary.frequency.notna().sum() == 2
    assert (summary.tipper == True).sum() == 1
    assert (summary.tipper == False).sum() == 2
    assert summary.tipper.isna().sum() == 2


def test_unknown_values_are_not_reported_as_zero(qapp):
    view = SurveyOverviewWidget()
    view.update_survey(pd.DataFrame({"ID": ["A", "B"]}))
    assert view._v_nfreq.text() == "—"
    assert view._v_elev.text() == "—"
    assert view._kpi_cards["tipper"].isHidden()
    assert "tipper" not in view._charts._kinds
    view.close()


def test_compact_row_animates_and_hides_tipper_on_reload(qapp):
    view = SurveyOverviewWidget()
    frame = pd.DataFrame({
        "ID": [f"S{i}" for i in range(18)],
        "Line": [f"Line {i % 3}" for i in range(18)],
        "Latitude": np.linspace(5, 5.1, 18),
        "Longitude": np.linspace(-4, -3.9, 18),
        "N_freq": [32] * 18, "Tipper": [True, False] * 9,
        "Elevation": np.arange(18) - 30,
    })
    view.resize(800, 520)
    view.show()
    view.update_survey(frame)
    qapp.processEvents()
    animation = view._charts._animation
    animation.setCurrentTime(300)
    assert 0 < view._charts._progress < 1
    animation.setCurrentTime(animation.duration())
    for name in ("stations", "nfreq", "elev", "tipper"):
        value = getattr(view, "_v_" + name)
        value._animation.setCurrentTime(value._animation.duration())
    qapp.processEvents()
    assert len(view._charts._cards) == 3
    assert sum(view._charts._hist) == 18
    assert view._v_tipper.text() == "50"
    assert view._v_elev.text() == "17"
    assert len({rect.top() for _, rect in view._charts._cards}) == 1
    assert view._scroll.horizontalScrollBar().maximum() == 0
    frame["Tipper"] = "false"
    view.update_survey(frame)
    qapp.processEvents()
    assert view._kpi_cards["tipper"].isHidden()
    assert view._charts._kinds == ["frequency", "footprint"]
    assert len(view._charts._cards) == 2
    view.clear()
    assert view._scroll.isHidden()
    assert animation.state() == animation.State.Stopped
    view.close()
