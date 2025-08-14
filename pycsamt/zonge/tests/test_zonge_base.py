# -*- coding: utf-8 -*-
import json
import re
from pathlib import Path

import pandas as pd
import pytest

from pycsamt.zonge._transfer import LegacyAVGBase  
from pycsamt.exceptions import AvgDataError
from pycsamt.zonge.base import (
    FieldAliases,
    AvgRow,
    AVGFrame,
    AVGComponentBase,
)

def test_field_aliases_have_expected_members():
    # A few representative spot checks (not exhaustive)
    assert "station" in FieldAliases.station
    assert "stn" in FieldAliases.station
    assert "hmag" in FieldAliases.hmag
    assert "b.mag" in FieldAliases.hmag
    assert "rho" in FieldAliases.rho
    assert "ares.mag" in FieldAliases.rho
    assert "phase" in FieldAliases.phase
    assert "z.phz" in FieldAliases.phase


def test_avgrow_defaults_and_str_repr_json_roundtrip():
    r = AvgRow(station=10, freq=256.0, comp="", rho=123.4)
    # Blank comp should default to "ExHy"
    assert r.comp == "ExHy"

    s = str(r)
    assert "AvgRow" in s and "f=256" in s and "rho=123.4" in s

    # JSON roundtrip sanity
    dct = r.asdict()
    js = json.dumps(dct)
    back = json.loads(js)
    assert back["station"] == 10
    assert pytest.approx(back["freq"], rel=0, abs=1e-12) == 256.0


def test_avgframe_core_helpers_and_reprs(tmp_path: Path):
    df = pd.DataFrame(
        {
            "station": [25, 75],
            "freq": [1.0, 2.0],
            "rho": [100.0, 200.0],
            "phase": [-45.0, -30.0],
        }
    )
    meta = {"Survey.Type": "CSAMT", "Unit.Length": "m"}
    frame = AVGFrame(df, meta, source=tmp_path / "K2.avg")

    assert frame.nrows == 2
    assert set(frame.columns) == {"station", "freq", "rho", "phase"}

    # copy() must be deep (mutations do not leak)
    copy = frame.copy()
    copy.data.loc[0, "rho"] = 999.0
    assert frame.data.loc[0, "rho"] == 100.0

    # reprs should be informative, not crash
    s, r = str(frame), repr(frame)
    assert "AVGFrame" in s and "cols=" in s
    assert "meta_keys" in r

    # JSON helpers produce strings
    assert isinstance(frame.to_json(), str)
    assert isinstance(frame.meta_as_json(), str)


class DummyComp(AVGComponentBase):
    """Minimal concrete component for testing the base API."""
    required = {"station", "freq"}
    provides = {"station", "freq"}

    def read(self, source: pd.DataFrame, meta=None) -> None:  # noqa: D401
        self._frame = source.copy()
        self._meta.update(dict(meta or {}))
        # Validate required columns exist in our working frame
        self._require(*self.required)
        # Keep a stable order for writing
        cols = ["station", "freq"]
        self._frame = self._frame.loc[:, cols]

    def write(self):
        return self._write_csv_block(
            cols=["station", "freq"],
            title="$Dummy Component",
            include_meta=True,
            stamp=True,
        )


def test_component_read_write_and_validation():
    df_ok = pd.DataFrame({"station": [25, 75], "freq": [1.0, 2.0]})
    meta = {"Survey.Type": "CSAMT"}

    comp = DummyComp.from_avg((df_ok, meta))
    assert comp.shape == (2, 2)
    assert comp.meta.get("Survey.Type") == "CSAMT"

    lines = comp.write()
    text = "\n".join(lines)
    # Expect banner, meta lines, stamp, and a CSV header
    assert "$Dummy Component" in text
    assert "$Survey.Type=CSAMT" in text
    assert "$Written=" in text
    assert "station,freq" in text
    assert "25,1" in text

    # Missing required column must raise a clear error
    df_bad = pd.DataFrame({"freq": [1.0, 2.0]})
    comp2 = DummyComp()
    with pytest.raises(Exception):
        comp2.read(df_bad, meta)


@pytest.mark.xfail(
    LegacyAVGBase is None,
    reason="LegacyAVGBase not present yet",
    strict=False,
)
def test_legacyavgbase_minimal_contract():
    # If the class exists, assert a minimal usable contract.
    assert LegacyAVGBase is not None

    # The class should be instantiable without exploding.
    obj = LegacyAVGBase()  # type: ignore[operator]
    # It should expose at least one of common transformation hooks.
    has_api = any(
        hasattr(obj, name) for name in ("to_xarray", "transform", "__call__")
    )
    assert has_api

    # And if it advertises to_xarray/transform, calling them should
    # either return something xarray-like or raise NotImplementedError.
    if hasattr(obj, "to_xarray"):
        with pytest.raises (
                AvgDataError, match=re.escape ("Empty legacy table.")
            ): 
            # try:
            out = obj.to_xarray(pd.DataFrame())
            # duck-type check for xarray.Dataset (no hard import)
            assert hasattr(out, "dims") and hasattr(out, "data_vars")
            # except AvgDataError:
            #     pass

    if hasattr(obj, "transform"):
        with pytest.raises (
                AvgDataError, match=re.escape ("Empty legacy table.")
            ): 
        # try:
            out = obj.transform(pd.DataFrame(), meta={})
            assert isinstance(out, (pd.DataFrame, dict))
        # except NotImplementedError:
        #     pass


if __name__=='__main__': # pragma: no-cover 
   pytest.main( [__file__])