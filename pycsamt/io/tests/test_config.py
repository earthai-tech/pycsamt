from __future__ import annotations

import pandas as pd

from pycsamt.io.config import Config


def test_parsers_maps_extensions_to_pandas_readers():
    parsers = Config().parsers
    assert parsers[".csv"] is pd.read_csv
    assert parsers[".xlsx"] is pd.read_excel
    assert parsers[".json"] is pd.read_json
    assert callable(parsers[".sql"])


def test_writers_maps_extensions_to_bound_dataframe_methods():
    df = pd.DataFrame({"a": [1, 2]})
    writers = Config.writers(df)
    assert writers[".csv"] == df.to_csv
    assert writers[".json"] == df.to_json
    assert callable(writers[".pkl"])
