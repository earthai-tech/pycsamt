from __future__ import annotations

import pytest

from pycsamt.api.labels import (
    FREQUENCY_LABEL,
    LOG10_PERIOD_LABEL,
    PERIOD_LABEL,
    period_axis_label,
)


@pytest.mark.parametrize(
    "kind,expected",
    [
        ("logperiod", LOG10_PERIOD_LABEL),
        ("log10period", LOG10_PERIOD_LABEL),
        ("log10t", LOG10_PERIOD_LABEL),
        ("logt", LOG10_PERIOD_LABEL),
        ("log_period", LOG10_PERIOD_LABEL),
        ("period", PERIOD_LABEL),
        ("t", PERIOD_LABEL),
        ("frequency", FREQUENCY_LABEL),
        ("freq", FREQUENCY_LABEL),
        ("f", FREQUENCY_LABEL),
        ("FREQUENCY", FREQUENCY_LABEL),
        ("  Period  ", PERIOD_LABEL),
    ],
)
def test_period_axis_label_known_kinds(kind, expected):
    assert period_axis_label(kind) == expected


def test_period_axis_label_default_is_logperiod():
    assert period_axis_label() == LOG10_PERIOD_LABEL


def test_period_axis_label_unknown_kind_returned_verbatim():
    assert period_axis_label("depth") == "depth"
