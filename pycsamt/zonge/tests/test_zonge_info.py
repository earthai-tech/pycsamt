# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later

from __future__ import annotations

import warnings
from pathlib import Path

import pandas as pd
import pytest

from pycsamt.exceptions import AvgDataError
from pycsamt.zonge.base import AVGFrame
from pycsamt.zonge.info import DataInfo


def _basic_df():
    return pd.DataFrame(
        {
            "station": [0.0, 50.0],
            "freq": [1024.0, 1024.0],
            "comp": ["ExHy", "ExHy"],
            "ARes.mag": [10.0, 20.0],
            "Z.phz": [45.0, 46.0],
        }
    )


def test_from_avg_with_path(data_path: Path):
    avg_file = data_path / "K1.AVG"
    if not avg_file.exists():
        pytest.skip(f"missing fixture: {avg_file}")
    info = DataInfo.from_avg(avg_file)
    assert isinstance(info, DataInfo)
    assert info.df is not None


def test_from_avg_with_avgframe_instance():
    frame = AVGFrame(_basic_df(), {})
    info = DataInfo.from_avg(frame)
    assert isinstance(info, DataInfo)
    assert info.df is not None


def test_from_avg_with_tuple():
    info = DataInfo.from_avg((_basic_df(), {"Unit.Rho": "ohm·m"}))
    assert isinstance(info, DataInfo)
    assert info.df is not None


def test_from_avg_with_dataframe_and_meta_kwarg():
    info = DataInfo.from_avg(_basic_df(), meta={"Unit.Rho": "ohm·m"})
    assert isinstance(info, DataInfo)
    assert info.df is not None


def test_from_avg_raises_typeerror_for_unsupported_source():
    with pytest.raises(TypeError):
        DataInfo.from_avg(12345)


def test_read_injects_comp_column_when_missing():
    df = _basic_df().drop(columns=["comp"])
    info = DataInfo()
    info.read(df, {})
    assert "comp" in info.df.columns
    assert set(info.df["comp"].unique()) == {"ExHy"}


def test_read_warns_and_skips_on_avg_data_error(monkeypatch):
    info = DataInfo()

    def _boom(self, df, meta):
        raise AvgDataError("simulated missing column")

    monkeypatch.setattr(type(info.z), "read", _boom)

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        info.read(_basic_df(), {})
    assert any("Could not load component" in str(w.message) for w in caught)


def test_read_logs_via_logger_when_verbose_on_avg_data_error(monkeypatch):
    info = DataInfo(verbose=True)

    def _boom(self, df, meta):
        raise AvgDataError("simulated missing column")

    monkeypatch.setattr(type(info.z), "read", _boom)

    logged = []
    monkeypatch.setattr(info._logger, "warning", lambda msg: logged.append(msg))
    info.read(_basic_df(), {})
    assert any("Could not load component" in m for m in logged)


def test_read_warns_on_unexpected_exception(monkeypatch):
    info = DataInfo()

    def _boom(self, df, meta):
        raise RuntimeError("unexpected failure")

    monkeypatch.setattr(type(info.z), "read", _boom)

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        info.read(_basic_df(), {})
    assert any("Unexpected error loading" in str(w.message) for w in caught)


def test_read_logs_via_logger_when_verbose_on_unexpected_exception(monkeypatch):
    info = DataInfo(verbose=True)

    def _boom(self, df, meta):
        raise RuntimeError("unexpected failure")

    monkeypatch.setattr(type(info.z), "read", _boom)

    logged = []
    monkeypatch.setattr(info._logger, "error", lambda msg: logged.append(msg))
    info.read(_basic_df(), {})
    assert any("Unexpected error loading" in m for m in logged)


def test_str_empty():
    info = DataInfo()
    assert str(info) == "DataInfo(empty)"
    assert repr(info) == str(info)


def test_str_nonempty():
    info = DataInfo()
    info.read(_basic_df(), {})
    text = str(info)
    assert "stations=2" in text
    assert "freqs=1" in text
    assert repr(info) == text
