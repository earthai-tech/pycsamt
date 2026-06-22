from __future__ import annotations

import numpy as np
import pandas as pd

from pycsamt.api import (
    APIFrame,
    APIResult,
    APISurvey,
    api_frame,
    configure_api_view,
    geology_dataframe,
    quality_dataframe,
    read_edis,
    reset_api_view,
    wrap_frame,
    wrap_result,
)
from pycsamt.api.view.progress import progress_enabled


def test_api_frame_keeps_pandas_behavior_and_display():
    reset_api_view()
    df = pd.DataFrame(
        {
            "station": ["S1", "S2"],
            "n_freq": [3, 4],
            "rho": [100.0, None],
        }
    )

    obj = APIFrame(df, name="sites", kind="edi.summary", source="demo")

    assert obj.df is df
    assert obj.data.shape == (2, 3)
    assert obj.station.tolist() == ["S1", "S2"]
    assert obj["n_freq"].sum() == 7
    assert obj.query("n_freq > 3")["station"].tolist() == ["S2"]
    assert obj.stats.rows == 2
    assert obj.stats.columns == 3
    assert obj.missing()["rho"] == 1

    text = str(obj)
    assert "APIFrame: sites" in text
    assert "shape: 2 rows x 3 columns" in text
    assert "source: demo" in text


def test_api_frame_is_editable_and_metadata_can_change():
    reset_api_view()
    obj = wrap_frame({"a": [1, 2]}, name="editable")

    obj["b"] = [3, 4]
    obj.df["c"] = [5, 6]
    obj.set_units(a="ohm.m").update_meta(owner="test")

    assert obj.b.tolist() == [3, 4]
    assert obj.c.tolist() == [5, 6]
    assert obj.units["a"] == "ohm.m"
    assert obj.meta["owner"] == "test"


def test_api_frame_decorator_wraps_dataframe_returns():
    reset_api_view()

    @api_frame(name="quality", kind="metadata.quality")
    def build():
        return pd.DataFrame({"score": [0.9, 0.5]})

    out = build()

    assert isinstance(out, APIFrame)
    assert out.name == "quality"
    assert out.kind == "metadata.quality"
    assert out.score.tolist() == [0.9, 0.5]


def test_api_result_wraps_tables():
    reset_api_view()
    res = wrap_result(
        {
            "summary": pd.DataFrame({"station": ["S1"]}),
            "note": "ok",
        },
        name="run",
        kind="qc",
    )

    assert isinstance(res, APIResult)
    assert isinstance(res.summary, APIFrame)
    assert res.note == "ok"
    assert "summary=1x1" in str(res)


def test_api_result_uses_configured_table_wrapper():
    try:
        configure_api_view(backend=False)
        res = wrap_result(
            {
                "summary": pd.DataFrame({"station": ["S1"]}),
            },
            name="run",
        )

        assert isinstance(res, APIResult)
        assert isinstance(res.summary, pd.DataFrame)
        assert not isinstance(res.summary, APIFrame)
    finally:
        reset_api_view()


def test_read_edis_returns_api_survey_for_empty_match(tmp_path):
    survey = read_edis(
        tmp_path / "*.edi",
        progress=False,
        strict=False,
    )

    assert isinstance(survey, APISurvey)
    assert len(survey) == 0
    assert isinstance(survey.summary(), APIFrame)
    assert survey.errors()


def test_progress_auto_is_boolean():
    assert progress_enabled(True) is True
    assert progress_enabled(False) is False


def test_quality_dataframe_follows_global_config_by_default():
    reset_api_view()  # global backend = "pycsamt"
    from pycsamt.metadata.quality import quality_dataframe as internal_quality_df

    class Site:
        name = "S1"
        freq = np.asarray([1.0, 10.0])
        z = np.ones((2, 2, 2), dtype=complex)
        tipper = np.ones(2, dtype=complex)

    sites = [Site()]

    # Public API uses @api_frame decorator — always wraps via global config.
    public = quality_dataframe(sites)
    assert isinstance(public, APIFrame)
    assert public.kind == "metadata.quality"

    # Internal function with api=None (default) also defers to global → wraps.
    internal_wrapped = internal_quality_df(sites)
    assert isinstance(internal_wrapped, APIFrame)

    # Internal function with api=False always returns raw pandas.
    raw = internal_quality_df(sites, api=False)
    assert isinstance(raw, pd.DataFrame)
    assert not isinstance(raw, APIFrame)
    assert internal_wrapped.df.equals(raw)


def test_public_geology_wrapper_returns_api_frame():
    reset_api_view()
    table = geology_dataframe()

    assert isinstance(table, APIFrame)
    assert table.kind == "metadata.geology"
    assert "name" in table.columns
    assert len(table) > 0


def test_metadata_api_flag_overrides_global():
    """Explicit api=True/False override the global config in either direction."""
    from pycsamt.metadata.geology import CATALOG
    from pycsamt.metadata.quality import quality_dataframe as quality_df

    class Site:
        name = "S1"
        freq = np.asarray([1.0, 10.0])
        z = np.ones((2, 2, 2), dtype=complex)
        tipper = np.ones(2, dtype=complex)

    try:
        # With global disabled, api=True still forces wrapping.
        configure_api_view(backend=False)
        q_api = quality_df([Site()], api=True)
        g_api = CATALOG.to_dataframe(api=True)
        assert isinstance(q_api, APIFrame)
        assert isinstance(g_api, APIFrame)

        # With global enabled, api=False still returns raw pandas.
        configure_api_view(backend="pycsamt")
        q_raw = quality_df([Site()], api=False)
        g_raw = CATALOG.to_dataframe(api=False)
        assert isinstance(q_raw, pd.DataFrame)
        assert not isinstance(q_raw, APIFrame)
        assert isinstance(g_raw, pd.DataFrame)
        assert not isinstance(g_raw, APIFrame)
    finally:
        reset_api_view()


def test_api_view_config_can_disable_wrapping():
    """configure_api_view(backend=False) disables wrapping for api=None calls;
    an explicit api=True per-call override still forces an APIFrame."""
    from pycsamt.metadata.quality import quality_dataframe as quality_df

    class Site:
        name = "S1"
        freq = np.asarray([1.0, 10.0])
        z = np.ones((2, 2, 2), dtype=complex)
        tipper = np.ones(2, dtype=complex)

    try:
        configure_api_view(backend=False)

        # Default call (api=None) respects global disabled setting.
        out_default = quality_df([Site()])
        assert isinstance(out_default, pd.DataFrame)
        assert not isinstance(out_default, APIFrame)

        # Explicit api=True overrides the global disabled setting.
        out_forced = quality_df([Site()], api=True)
        assert isinstance(out_forced, APIFrame)
    finally:
        reset_api_view()


def test_api_view_config_accepts_custom_wrapper():
    calls = []

    def custom(data, **metadata):
        calls.append(metadata)
        return {
            "data": data,
            "name": metadata.get("name"),
            "kind": metadata.get("kind"),
        }

    try:
        configure_api_view(wrapper=custom)
        out = wrap_frame(pd.DataFrame({"a": [1]}), name="custom", kind="demo")

        assert out["name"] == "custom"
        assert out["kind"] == "demo"
        assert out["data"].equals(pd.DataFrame({"a": [1]}))
        assert calls and calls[0]["name"] == "custom"
    finally:
        reset_api_view()
