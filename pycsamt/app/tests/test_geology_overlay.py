# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Shared PCGL (Interpretation legend) app-adapter tests."""

from __future__ import annotations

import base64
import json

import numpy as np
import pytest

from pycsamt.app._geology import (
    TABLE_COLUMNS,
    auto_suggest_legend,
    decode_geology_upload,
    geology_bands_from_store,
    legend_from_store,
    legend_from_table_rows,
    store_from_legend,
    table_rows_from_legend,
)
from pycsamt.format.geology import (
    GeologyLegend,
    GeologyLegendValidationError,
    legend_to_dict,
)
from pycsamt.geology import RockDatabase, RockEntry


def _legend() -> GeologyLegend:
    db = RockDatabase(
        [
            RockEntry("Sand (saturated)", 10.0, 200.0, "#E9C46A"),
            RockEntry(
                "Granodiorite", 1000.0, 5000.0, "#8D99AE",
                pattern_id="cross-hatch-01", pattern_source="builtin",
            ),
        ]
    )
    return GeologyLegend.from_rock_database(db, document_id="pcgl:t")


def _upload_json(legend: GeologyLegend) -> str:
    raw = json.dumps(legend_to_dict(legend)).encode()
    return "data:application/json;base64," + base64.b64encode(raw).decode()


def _upload_csv(text: str) -> str:
    return "data:text/csv;base64," + base64.b64encode(text.encode()).decode()


def test_decode_upload_round_trips_json():
    store = decode_geology_upload(_upload_json(_legend()), "legend.pcgl.json")
    assert store["n_entries"] == 2
    restored = legend_from_store(store)
    assert [e.name for e in restored.entries] == [
        "Sand (saturated)", "Granodiorite",
    ]


def test_decode_upload_accepts_plain_csv():
    csv_text = "name,rho_min,rho_max,color\nClay,1,20,#111111\n"
    store = decode_geology_upload(_upload_csv(csv_text), "legend.csv")
    assert store["n_entries"] == 1
    restored = legend_from_store(store)
    assert restored.entries[0].name == "Clay"


def test_decode_upload_rejects_bad_payload():
    with pytest.raises(ValueError):
        decode_geology_upload("not a data url")


def test_geology_bands_from_store_shape_and_order():
    store = store_from_legend(_legend())
    bands = geology_bands_from_store(store)
    assert bands == (
        (10.0, 200.0, "#E9C46A", "Sand (saturated)"),
        (1000.0, 5000.0, "#8D99AE", "Granodiorite"),
    )


def test_geology_bands_from_store_none_when_empty():
    assert geology_bands_from_store(None) is None
    assert geology_bands_from_store({}) is None


def test_table_round_trip_preserves_pattern_columns():
    legend = _legend()
    rows = table_rows_from_legend(legend)
    assert set(rows[0]) == set(TABLE_COLUMNS)
    restored = legend_from_table_rows(rows, title="from rows")
    assert restored.entries[1].pattern_id == "cross-hatch-01"
    assert restored.issues() == []


def test_legend_from_table_rows_skips_blank_and_invalid_rows():
    rows = [
        {"name": "Clay", "rho_min": "1", "rho_max": "20"},
        {"name": "", "rho_min": "5", "rho_max": "50"},  # blank name
        {"name": "Bad", "rho_min": "x", "rho_max": "50"},  # not numeric
    ]
    legend = legend_from_table_rows(rows)
    assert [e.name for e in legend.entries] == ["Clay"]


def test_legend_from_table_rows_raises_on_invalid_ranges():
    rows = [{"name": "Bad", "rho_min": "100", "rho_max": "10"}]
    with pytest.raises(GeologyLegendValidationError):
        legend_from_table_rows(rows)


def test_auto_suggest_covers_the_actual_data_range_not_the_full_db():
    rng = np.random.default_rng(0)
    rho = np.concatenate(
        [rng.uniform(10, 200, 60), rng.uniform(1000, 5000, 60)]
    )
    legend = auto_suggest_legend(rho, n_bins=10)
    assert legend.entries[0].rho_min == pytest.approx(rho.min(), rel=0.05)
    assert legend.entries[-1].rho_max == pytest.approx(rho.max(), rel=0.05)
    assert legend.issues() == []


def test_auto_suggest_merges_adjacent_bins_with_the_same_rock():
    # A narrow range that classifies to one rock entry throughout should
    # collapse to very few merged rows, not one row per raw bin.
    rng = np.random.default_rng(1)
    rho = rng.uniform(140.0, 160.0, 50)
    legend = auto_suggest_legend(rho, n_bins=10)
    assert len(legend.entries) <= 2


def test_auto_suggest_rejects_degenerate_input():
    with pytest.raises(GeologyLegendValidationError):
        auto_suggest_legend(np.array([100.0, 100.0, 100.0]))
    with pytest.raises(GeologyLegendValidationError):
        auto_suggest_legend(np.array([]))


# ---------------------------------------------------------------------------
# colour consistency with loaded boreholes
# ---------------------------------------------------------------------------


def _pcbh_document():
    from pycsamt.format.borehole import (
        Collar,
        CoordinateReferenceSystem,
        LogInterval,
        PCBHBorehole,
        PCBHDocument,
        VocabularyEntry,
    )

    return PCBHDocument(
        document_id="test:colors",
        created_at="2026-09-04T00:00:00Z",
        created_by="pytest",
        crs=CoordinateReferenceSystem("EPSG:4326"),
        lithologies=[
            VocabularyEntry("A", "Granodiorite", color="#C7AA54"),
            VocabularyEntry("B", "Diorite porphyry"),  # no explicit colour
        ],
        boreholes=[
            PCBHBorehole(
                id="BH1", name="Hole 1", kind="mining_exploration",
                status="completed",
                collar=Collar(0.0, 0.0, 0.0),
                total_depth_md=60.0,
                interval_logs={
                    "lithology": [
                        LogInterval(0.0, 30.0, code="A"),
                        LogInterval(30.0, 60.0, code="B"),
                    ]
                },
            )
        ],
    )


def test_borehole_lithology_colors_uses_explicit_color():
    from pycsamt.app._geology import borehole_lithology_colors

    colors = borehole_lithology_colors(_pcbh_document())
    assert colors["granodiorite"] == "#C7AA54"


def test_borehole_lithology_colors_falls_back_to_the_render_pipelines_own_deterministic_color():
    from pycsamt.format.borehole import deterministic_color

    from pycsamt.app._geology import borehole_lithology_colors

    colors = borehole_lithology_colors(_pcbh_document())
    assert colors["diorite porphyry"] == deterministic_color("B")


def test_borehole_lithology_colors_handles_no_document():
    from pycsamt.app._geology import borehole_lithology_colors

    assert borehole_lithology_colors(None) == {}


def test_sync_legend_colors_matches_exact_names_case_and_whitespace_insensitively():
    from pycsamt.app._geology import (
        borehole_lithology_colors,
        sync_legend_colors_with_boreholes,
    )

    colors = borehole_lithology_colors(_pcbh_document())
    rows = [
        {"name": "  GRANODIORITE  ", "color": "#000000"},
        {"name": "Not drilled here", "color": "#111111"},
    ]
    updated, n = sync_legend_colors_with_boreholes(rows, colors)
    assert n == 1
    assert updated[0]["color"] == "#C7AA54"
    assert updated[1]["color"] == "#111111"


def test_sync_legend_colors_does_not_conflate_similar_but_different_names():
    from pycsamt.app._geology import sync_legend_colors_with_boreholes

    colors = {"pyrite-ized granodiorite porphyry": "#C754AD"}
    rows = [{"name": "Pyrite-mineralized granodiorite porphyry", "color": "#000000"}]
    updated, n = sync_legend_colors_with_boreholes(rows, colors)
    assert n == 0
    assert updated[0]["color"] == "#000000"


def test_sync_legend_colors_handles_empty_inputs():
    from pycsamt.app._geology import sync_legend_colors_with_boreholes

    assert sync_legend_colors_with_boreholes(None, {}) == ([], 0)
    assert sync_legend_colors_with_boreholes([], {"a": "#fff"}) == ([], 0)


# ---------------------------------------------------------------------------
# pattern-texture fill: geology_pattern_stencils_from_store
# ---------------------------------------------------------------------------


@pytest.fixture(autouse=False)
def _isolated_pattern_cache(tmp_path, monkeypatch):
    monkeypatch.setenv("PYCSAMT_PATTERN_CACHE", str(tmp_path / "cache"))
    yield


def test_pattern_stencils_from_store_uses_the_builtin_pack(
    _isolated_pattern_cache,
):
    from pycsamt.app._geology import geology_pattern_stencils_from_store

    db = RockDatabase(
        [
            RockEntry(
                "Sand", 10.0, 200.0, "#E9C46A",
                pattern_id="builtin-01", pattern_source="builtin",
            ),
            RockEntry("Granite", 1000.0, 5000.0, "#8D99AE"),
        ]
    )
    legend = GeologyLegend.from_rock_database(db)
    stencils = geology_pattern_stencils_from_store(legend)
    assert list(stencils) == ["Sand"]
    arr = stencils["Sand"]
    assert arr.ndim == 2
    assert 0.0 <= arr.min() and arr.max() <= 1.0


def test_pattern_stencils_from_store_skips_a_stale_pattern_reference(
    _isolated_pattern_cache,
):
    from pycsamt.app._geology import geology_pattern_stencils_from_store

    db = RockDatabase(
        [
            RockEntry(
                "Sand", 10.0, 200.0, "#E9C46A",
                pattern_id="does-not-exist", pattern_source="builtin",
            ),
        ]
    )
    legend = GeologyLegend.from_rock_database(db)
    assert geology_pattern_stencils_from_store(legend) == {}


def test_pattern_stencils_from_store_handles_no_legend():
    from pycsamt.app._geology import geology_pattern_stencils_from_store

    assert geology_pattern_stencils_from_store(None) == {}
