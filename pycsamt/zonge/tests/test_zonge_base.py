import json
import re
from pathlib import Path

import pandas as pd
import pytest

from pycsamt.exceptions import AvgDataError
from pycsamt.zonge._transfer import LegacyAVGBase
from pycsamt.zonge.base import (
    AVGComponentBase,
    AVGFrame,
    AvgRow,
    FieldAliases,
    guess_kind_from_df,
)


def _norm_field(field_value):
    return (str(v).lower() for v in field_value)


def test_field_aliases_have_expected_members():
    # A few representative spot checks (not exhaustive)
    assert "station" in _norm_field(FieldAliases.station)
    assert "stn" in _norm_field(FieldAliases.station)
    assert "hmag" in _norm_field(FieldAliases.hmag)
    assert "b.mag" in _norm_field(FieldAliases.hmag)
    assert "ares.mag" in _norm_field(FieldAliases.rho)
    assert "phase" in _norm_field(FieldAliases.phase)
    assert "z.phz" in _norm_field(FieldAliases.phase)


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

    # asdict() combines data + meta + source into plain types
    d = frame.asdict()
    assert d["meta"] == meta
    assert d["source"] == str(tmp_path / "K2.avg")
    assert isinstance(d["data"], list)


def test_avgframe_asdict_no_source():
    df = pd.DataFrame({"station": [1], "freq": [1.0]})
    frame = AVGFrame(df, {})
    assert frame.asdict()["source"] is None


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


def test_component_from_avg_accepts_avgframe_instance():
    frame = AVGFrame(pd.DataFrame({"station": [1], "freq": [1.0]}), {})
    comp = DummyComp.from_avg(frame)
    assert comp.shape == (1, 2)


def test_component_from_avg_accepts_bare_dataframe_and_meta_kwarg():
    df = pd.DataFrame({"station": [1], "freq": [1.0]})
    comp = DummyComp.from_avg(df, meta={"Survey.Type": "CSAMT"})
    assert comp.meta.get("Survey.Type") == "CSAMT"


def test_component_from_avg_falls_back_on_single_arg_read():
    class SingleArgComp(AVGComponentBase):
        def read(self, source: pd.DataFrame) -> None:  # only one arg
            self._frame = source.copy()

        def write(self):
            return []

    df = pd.DataFrame({"station": [1], "freq": [1.0]})
    comp = SingleArgComp.from_avg((df, {"x": 1}))
    assert comp.shape == (1, 2)


def test_component_from_avg_raises_typeerror_for_unsupported_source():
    with pytest.raises(TypeError):
        DummyComp.from_avg(12345)


def test_component_from_avg_with_legacy_path_transforms(data_path: Path):
    legacy_file = data_path / "K1.AVG"
    if not legacy_file.exists():
        pytest.skip(f"missing fixture: {legacy_file}")

    class PassthroughComp(AVGComponentBase):
        def read(self, source: pd.DataFrame, meta=None) -> None:
            self._frame = source.copy()
            self._meta.update(dict(meta or {}))

        def write(self):
            return []

    comp = PassthroughComp.from_avg(legacy_file)
    assert not comp.frame.empty


def test_component_from_avg_with_modern_path(data_path: Path):
    modern_file = data_path / "K2.AVG"
    if not modern_file.exists():
        pytest.skip(f"missing fixture: {modern_file}")

    class PassthroughComp(AVGComponentBase):
        def read(self, source: pd.DataFrame, meta=None) -> None:
            self._frame = source.copy()
            self._meta.update(dict(meta or {}))

        def write(self):
            return []

    comp = PassthroughComp.from_avg(modern_file)
    assert not comp.frame.empty


def test_component_name_property():
    comp = DummyComp(name="Custom")
    assert comp.name == "Custom"


def test_component_asdict_with_and_without_meta():
    df = pd.DataFrame({"station": [1], "freq": [1.0]})
    comp = DummyComp.from_avg((df, {"k": "v"}))
    with_meta = comp.asdict()
    assert with_meta["meta"] == {"k": "v"}
    without_meta = comp.asdict(include_meta=False)
    assert "meta" not in without_meta


def test_component_to_json():
    df = pd.DataFrame({"station": [1], "freq": [1.0]})
    comp = DummyComp.from_avg((df, {}))
    text = comp.to_json()
    assert isinstance(text, str) and "station" in text


def test_write_csv_block_empty_frame_returns_blank_line():
    comp = DummyComp()
    lines = comp.write()
    assert lines[-1] == ""


def test_component_str_and_repr():
    df = pd.DataFrame({"station": [1], "freq": [1.0]})
    comp = DummyComp.from_avg((df, {}))
    s = str(comp)
    assert "DummyComp" in s and "station" in s
    assert repr(comp) == s


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
    has_api = any(hasattr(obj, name) for name in ("to_xarray", "transform", "__call__"))
    assert has_api

    # And if it advertises to_xarray/transform, calling them should
    # either return something xarray-like or raise NotImplementedError.
    if hasattr(obj, "to_xarray"):
        with pytest.raises(AvgDataError, match=re.escape("Empty legacy table.")):
            # try:
            out = obj.to_xarray(pd.DataFrame())
            # duck-type check for xarray.Dataset (no hard import)
            assert hasattr(out, "dims") and hasattr(out, "data_vars")
            # except AvgDataError:
            #     pass

    if hasattr(obj, "transform"):
        with pytest.raises(AvgDataError, match=re.escape("Empty legacy table.")):
            # try:
            out = obj.transform(pd.DataFrame(), meta={})
            assert isinstance(out, (pd.DataFrame, dict))
        # except NotImplementedError:
        #     pass


# --------------------------- guess_kind_from_df ---------------------------- #


def test_guess_kind_raises_for_unsupported_input():
    with pytest.raises(AvgDataError, match="pandas DataFrame or an AVGFrame"):
        guess_kind_from_df(12345)


def test_guess_kind_from_plain_dataframe_modern_dot_notation():
    df = pd.DataFrame({"ARes.mag": [1.0], "Freq": [1.0]})
    assert guess_kind_from_df(df) == 2


def test_guess_kind_from_avgframe_instance():
    frame = AVGFrame(pd.DataFrame({"station": [1], "freq": [1.0]}), {})
    # No modern dot-notation or legacy/canonical indicators -> falls
    # through to the "no clear indicators" default (kind=2, soft mode).
    assert guess_kind_from_df(frame) == 2


def test_guess_kind_detects_legacy_indicators():
    df = pd.DataFrame({"station": [1], "freq": [1.0], "sPhz": [1.0]})
    assert guess_kind_from_df(df) == 1


def test_guess_kind_detects_canonical_modern_names():
    df = pd.DataFrame({"station": [1], "freq": [1.0], "pc_emag": [1.0]})
    assert guess_kind_from_df(df) == 2


def test_guess_kind_strict_raise_when_no_indicators():
    df = pd.DataFrame({"station": [1], "freq": [1.0]})
    with pytest.raises(AvgDataError, match="Could not determine AVG kind"):
        guess_kind_from_df(df, mode="strict", error="raise")


def test_guess_kind_strict_warn_when_no_indicators():
    df = pd.DataFrame({"station": [1], "freq": [1.0]})
    with pytest.warns(UserWarning, match="Defaulting to modern"):
        out = guess_kind_from_df(df, mode="strict", error="warn")
    assert out == 2


def test_guess_kind_strict_ignore_when_no_indicators():
    df = pd.DataFrame({"station": [1], "freq": [1.0]})
    out = guess_kind_from_df(df, mode="strict", error="ignore")
    assert out == 2


def test_guess_kind_transform_legacy_returns_tuple():
    # Raw legacy column names (pre-standardization) so guess_kind_from_df
    # classifies this as kind=1 and exercises the transform branch.
    df = pd.DataFrame(
        {"station": [0.0], "freq": [1.0], "comp": ["ExHy"], "sPhz": [1.0]}
    )
    out = guess_kind_from_df(df, meta={}, transform=True, verbose=True)
    assert isinstance(out, tuple) and len(out) == 3
    out_df, out_meta, kind = out
    assert isinstance(out_df, pd.DataFrame)
    assert kind == 2


def test_guess_kind_transform_true_but_already_modern_returns_asis():
    df = pd.DataFrame({"station": [1], "freq": [1.0], "ARes.mag": [10.0]})
    out = guess_kind_from_df(df, meta={"a": 1}, transform=True)
    assert isinstance(out, tuple) and len(out) == 3
    out_df, out_meta, kind = out
    assert out_df is df
    assert out_meta == {"a": 1}
    assert kind == 2


if __name__ == "__main__":  # pragma: no-cover
    pytest.main([__file__])
