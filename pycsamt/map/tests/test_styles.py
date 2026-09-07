# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for :mod:`pycsamt.map.styles`."""

from __future__ import annotations

import math

import pytest

from pycsamt.map.styles import (
    GEOLOGY_GAP_COLOR,
    geology_colorbar_ticks,
    geology_colorscale,
    geology_crange,
    geology_legend_shapes_annotations,
    to_plotly_cmap,
)


def test_to_plotly_cmap_falls_back_and_remaps():
    assert to_plotly_cmap(None) == "plasma"
    assert to_plotly_cmap("coolwarm") == "balance"
    assert to_plotly_cmap("plasma") == "plasma"


def test_geology_crange_spans_all_bands_in_log10_space():
    bands = [(10.0, 200.0, "#111"), (1000.0, 20000.0, "#222")]
    lo, hi = geology_crange(bands)
    assert lo == pytest.approx(math.log10(10.0))
    assert hi == pytest.approx(math.log10(20000.0))


def test_geology_crange_ignores_degenerate_bands():
    bands = [(0.0, 0.0, "#111"), (5.0, 5.0, "#222"), (10.0, 100.0, "#333")]
    lo, hi = geology_crange(bands)
    assert lo == pytest.approx(math.log10(10.0))
    assert hi == pytest.approx(math.log10(100.0))


def test_geology_crange_empty_falls_back_to_unit_range():
    assert geology_crange([]) == (0.0, 1.0)


def test_geology_colorscale_bounds_are_exactly_0_and_1():
    bands = [(10.0, 200.0, "#111"), (1000.0, 20000.0, "#222")]
    scale = geology_colorscale(bands)
    assert scale[0][0] == 0.0
    assert scale[-1][0] == 1.0


def test_geology_colorscale_is_a_hard_step_not_a_gradient():
    bands = [(10.0, 200.0, "#111"), (1000.0, 20000.0, "#222")]
    scale = geology_colorscale(bands)
    values = [s[0] for s in scale]
    # every internal boundary appears twice (end of one flat run, start
    # of the next) -- that duplication is what makes Plotly render a
    # cliff instead of interpolating.
    assert len(values) != len(set(values))


def test_geology_colorscale_fills_gaps_with_the_neutral_gap_color():
    # A gap between the two bands, and margin at both ends.
    bands = [(100.0, 200.0, "#111"), (1000.0, 2000.0, "#222")]
    scale = geology_colorscale(bands, log_lo=0.0, log_hi=5.0)
    colors = [s[1] for s in scale]
    assert GEOLOGY_GAP_COLOR in colors
    assert "#111" in colors and "#222" in colors


def test_geology_colorscale_single_band_covering_full_range():
    bands = [(1.0, 100.0, "#abcdef")]
    scale = geology_colorscale(bands)
    assert scale == [[0.0, "#abcdef"], [1.0, "#abcdef"]]


def test_geology_colorscale_empty_bands_is_flat_gap_color():
    scale = geology_colorscale([])
    assert scale == [[0.0, GEOLOGY_GAP_COLOR], [1.0, GEOLOGY_GAP_COLOR]]


def test_geology_colorscale_respects_explicit_log_range_override():
    bands = [(10.0, 100.0, "#111")]
    # Explicit range wider than the band's own span.
    scale = geology_colorscale(bands, log_lo=0.0, log_hi=4.0)
    # band spans log10(10)=1 .. log10(100)=2, normalized into [0,4] ->
    # [0.25, 0.5]; t=0.5 also starts the trailing gap, so both colours
    # legitimately appear at that boundary.
    band_stops = [(round(t, 4), c) for t, c in scale if c == "#111"]
    assert (0.25, "#111") in band_stops
    assert (0.5, "#111") in band_stops


def test_pattern_band_stops_goes_from_light_tint_to_full_color():
    from pycsamt.map.styles import pattern_band_stops

    stops = pattern_band_stops("#8D99AE", 4)
    assert len(stops) == 4
    assert stops[-1] == "#8D99AE"
    assert stops[0] != "#8D99AE"


def test_pattern_band_stops_minimum_two_stops():
    from pycsamt.map.styles import pattern_band_stops

    assert len(pattern_band_stops("#000000", 1)) == 2


def test_geology_colorscale_textured_band_gets_a_gradient_sub_range():
    bands = [
        (10.0, 200.0, "#111111", "Sand"),
        (1000.0, 20000.0, "#222222", "Granite"),
    ]
    plain = geology_colorscale(bands)
    textured = geology_colorscale(bands, textured={"Granite": 8})
    assert len(textured) > len(plain)
    # the untextured "Sand" band is unaffected -- still a flat 2-stop run.
    sand_stops = [s for s in textured if s[1] == "#111111"]
    assert len(sand_stops) == 2
    # the textured band's stops span from a light tint up to its own
    # full colour, strictly increasing in t.
    granite_stops = [s for s in textured if s[0] >= sand_stops[-1][0]]
    ts = [s[0] for s in granite_stops]
    assert ts == sorted(ts)
    assert granite_stops[-1][1] == "#222222"


def test_geology_colorscale_textured_ignores_unknown_band_names():
    bands = [(10.0, 200.0, "#111111", "Sand")]
    scale = geology_colorscale(bands, textured={"Nonexistent": 8})
    assert scale == geology_colorscale(bands)


def test_geology_colorbar_ticks_skips_bands_with_no_name():
    bands = [(10.0, 200.0, "#111"), (1000.0, 5000.0, "#222")]
    tickvals, ticktext = geology_colorbar_ticks(bands)
    assert tickvals == [] and ticktext == []


def test_geology_colorbar_ticks_labels_each_named_band_at_its_geometric_mean():
    bands = [
        (10.0, 200.0, "#111", "Sand"),
        (1000.0, 5000.0, "#222", "Granite"),
    ]
    tickvals, ticktext = geology_colorbar_ticks(bands)
    assert ticktext == ["Sand", "Granite"]
    assert tickvals[0] == pytest.approx(math.log10(math.sqrt(10.0 * 200.0)))
    # every tick lands strictly inside its own band, not on a boundary.
    assert math.log10(10.0) < tickvals[0] < math.log10(200.0)


def test_geology_colorbar_ticks_mixed_named_and_unnamed_bands():
    bands = [
        (10.0, 200.0, "#111", "Sand"),
        (1000.0, 5000.0, "#222"),  # unnamed -- skipped
    ]
    tickvals, ticktext = geology_colorbar_ticks(bands)
    assert ticktext == ["Sand"]
    assert len(tickvals) == 1


# ---------------------------------------------------------------------------
# geology_legend_shapes_annotations
# ---------------------------------------------------------------------------


def test_geology_legend_empty_bands_returns_nothing():
    shapes, annotations = geology_legend_shapes_annotations([])
    assert shapes == [] and annotations == []


def test_geology_legend_one_chip_per_band_named_by_rock():
    bands = [
        (10.0, 200.0, "#111", "Sand"),
        (1000.0, 5000.0, "#222", "Granite"),
    ]
    _, annotations = geology_legend_shapes_annotations(bands)
    assert [a["text"].strip() for a in annotations] == ["Sand", "Granite"]
    assert annotations[0]["bgcolor"] == "#111"
    assert annotations[1]["bgcolor"] == "#222"


def test_geology_legend_unnamed_band_falls_back_to_range_label():
    bands = [(10.0, 200.0, "#111")]
    _, annotations = geology_legend_shapes_annotations(bands)
    assert "10" in annotations[0]["text"] and "200" in annotations[0]["text"]


def test_geology_legend_rows_are_evenly_spaced_by_index_not_value():
    # Two bands with wildly different resistivity separations must still
    # sit exactly one row height apart -- position is index-based.
    bands = [
        (10.0, 11.0, "#111", "A"),
        (1.0, 1e6, "#222", "B"),
    ]
    _, annotations = geology_legend_shapes_annotations(bands, row_height=0.05)
    y0, y1 = annotations[0]["y"], annotations[1]["y"]
    assert y0 - y1 == pytest.approx(0.05)


def test_geology_legend_never_overlaps_with_many_bands():
    bands = [
        (float(10**i), float(10 ** (i + 1)), "#334455", f"Unit {i}")
        for i in range(40)
    ]
    _, annotations = geology_legend_shapes_annotations(bands)
    ys = [a["y"] for a in annotations]
    gaps = [ys[i] - ys[i + 1] for i in range(len(ys) - 1)]
    assert all(gap > 0 for gap in gaps)
    assert min(gaps) >= 0.026 - 1e-9  # min_row_height default
    # too many to fit -- truncated with a summary row instead of
    # shrinking into overlap or running past the canvas.
    assert annotations[-1]["text"].startswith("+")
    assert len(annotations) < 40


def test_geology_legend_respects_custom_min_row_height_and_row_count():
    bands = [
        (float(10**i), float(10 ** (i + 1)), "#334455", f"Unit {i}")
        for i in range(10)
    ]
    _, annotations = geology_legend_shapes_annotations(
        bands, row_height=0.05, min_row_height=0.05, y_top=0.3,
    )
    # Only room for a handful of full-height rows before truncating.
    assert len(annotations) < 10
    assert annotations[-1]["text"].startswith("+")


def test_contrast_text_color_picks_readable_text():
    from pycsamt.map.styles import _contrast_text_color

    assert _contrast_text_color("#ffffff") == "#111111"
    assert _contrast_text_color("#000000") == "#f5f5f5"
