from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pycsamt.api.view.config import reset_api_view
from pycsamt.api.view.frame import APIFrame
from pycsamt.api.view.tables import (
    geology_dataframe,
    geology_table,
    quality_dataframe,
    quality_table,
    sites_summary,
)

_DATA_ROOT = Path(__file__).resolve().parents[4] / "data"
_THREE_EDIS = _DATA_ROOT / "3edis"


@pytest.fixture(autouse=True)
def _reset_config():
    reset_api_view()
    yield
    reset_api_view()


class _FakeSite:
    name = "S1"
    freq = np.asarray([1.0, 10.0])
    z = np.ones((2, 2, 2), dtype=complex)
    tipper = np.ones(2, dtype=complex)


def test_quality_dataframe_returns_api_frame():
    result = quality_dataframe([_FakeSite()])
    assert isinstance(result, APIFrame)
    assert result.kind == "metadata.quality"


def test_quality_table_is_alias_for_quality_dataframe():
    result = quality_table([_FakeSite()])
    assert isinstance(result, APIFrame)
    assert result.kind == "metadata.quality"


def test_geology_dataframe_default_catalog():
    result = geology_dataframe()
    assert isinstance(result, APIFrame)
    assert result.kind == "metadata.geology"
    assert result.name == "geology_catalog"
    assert len(result) > 0


def test_geology_dataframe_custom_name_and_catalog():
    class FakeCatalog:
        def to_dataframe(self):
            import pandas as pd

            return pd.DataFrame({"formation": ["granite"]})

    result = geology_dataframe(FakeCatalog(), name="custom")
    assert result.name == "custom"
    assert result.source == "FakeCatalog"
    assert result.df["formation"].tolist() == ["granite"]


def test_geology_table_is_alias_for_geology_dataframe():
    result = geology_table()
    assert isinstance(result, APIFrame)
    assert result.kind == "metadata.geology"


@pytest.mark.skipif(not _THREE_EDIS.exists(), reason="sample EDI data not found")
def test_sites_summary_returns_api_frame_for_real_data():
    result = sites_summary(_THREE_EDIS, strict=False)
    assert isinstance(result, APIFrame)
    assert result.kind == "edi.summary"
    assert "station" in result.columns
    assert len(result) == 3


@pytest.mark.skipif(not _THREE_EDIS.exists(), reason="sample EDI data not found")
def test_sites_summary_custom_fields():
    result = sites_summary(_THREE_EDIS, fields=("station", "n_freq"))
    assert list(result.columns) == ["station", "n_freq"]
