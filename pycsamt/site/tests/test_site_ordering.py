"""Regression tests for spatial ordering in Sites and ensure_sites."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

from pycsamt.emtools._core import ensure_sites
from pycsamt.emtools.inspect import pseudosection
from pycsamt.emtools.ss import _order_sites
from pycsamt.seg.edi import EDIFile
from pycsamt.site.base import Site, Sites


def _make_sites(
    tmp_path: Path,
    simulated_edi: Path,
    records: list[tuple[str, float, float]],
) -> Sites:
    """Build real EDI-backed sites in exactly the supplied order."""

    edis = []
    for serial, (name, lat, lon) in enumerate(records):
        path = tmp_path / f"{serial:02d}_{name}.edi"
        path.write_text(simulated_edi.read_text(encoding="utf-8"), encoding="utf-8")
        edi = EDIFile(path)  # type: ignore[arg-type]
        site = Site(edi).rename(name, inplace=True)
        site.set_coords(lat, lon, float(serial), inplace=True)
        edis.append(site.edi)
    return Sites(edis)


def _line_records(
    azimuth: float,
    order: tuple[int, ...] = (4, 1, 5, 2, 3),
    *,
    cross_track_m: tuple[float, ...] | None = None,
) -> list[tuple[str, float, float]]:
    """Return shuffled geographic points along a metric profile."""

    lat0, lon0 = 35.0, -117.0
    az = math.radians(azimuth)
    along = np.arange(5, dtype=float) * 125.0
    cross = np.zeros(5) if cross_track_m is None else np.asarray(cross_track_m)
    east = along * math.sin(az) + cross * math.cos(az)
    north = along * math.cos(az) - cross * math.sin(az)
    lat = lat0 + north / 111_000.0
    lon = lon0 + east / (111_000.0 * math.cos(math.radians(lat0)))
    return [(f"S{i:02d}", float(lat[i - 1]), float(lon[i - 1])) for i in order]


@pytest.mark.parametrize("azimuth", [0.0, 25.0, 90.0, 135.0, 225.0, 315.0])
def test_auto_orders_shuffled_straight_lines_at_any_azimuth(
    tmp_path: Path, simulated_edi: Path, azimuth: float
) -> None:
    sites = _make_sites(tmp_path, simulated_edi, _line_records(azimuth))

    ordered = sites.ordered("auto")

    assert [s.name for s in ordered] == ["S01", "S02", "S03", "S04", "S05"]
    assert ordered.ordering["applied"] == "chainage"
    assert ordered.ordering["linearity"] > 0.999
    assert ordered.ordering["cross_track_ratio"] < 1e-6
    assert ordered.ordering["span_m"] == pytest.approx(500.0, rel=0.01)


def test_auto_tolerates_small_realistic_gps_cross_track_noise(
    tmp_path: Path, simulated_edi: Path
) -> None:
    records = _line_records(55.0, cross_track_m=(2.0, -3.0, 1.0, 4.0, -2.0))
    sites = _make_sites(tmp_path, simulated_edi, records)

    ordered = sites.ordered("auto")

    assert [s.name for s in ordered] == ["S01", "S02", "S03", "S04", "S05"]
    assert ordered.ordering["applied"] == "chainage"
    assert ordered.ordering["cross_track_ratio"] < 0.03


def test_auto_rejects_non_line_square_geometry(
    tmp_path: Path, simulated_edi: Path
) -> None:
    records = [
        ("S03", 35.0, -117.0),
        ("S01", 35.0, -116.99),
        ("S04", 35.01, -116.99),
        ("S02", 35.01, -117.0),
    ]
    sites = _make_sites(tmp_path, simulated_edi, records)

    ordered = sites.ordered("auto")

    assert [s.name for s in ordered] == [r[0] for r in records]
    assert ordered.ordering["applied"] == "input"
    assert "credible straight line" in ordered.ordering["reason"]


def test_auto_rejects_duplicate_coordinate_cloud(
    tmp_path: Path, simulated_edi: Path
) -> None:
    records = [(f"S{i:02d}", 35.0, -117.0) for i in (3, 1, 2)]
    sites = _make_sites(tmp_path, simulated_edi, records)

    ordered = sites.ordered("auto")

    assert [s.name for s in ordered] == ["S03", "S01", "S02"]
    assert ordered.ordering["applied"] == "input"


def test_auto_requires_minimum_coordinate_fraction(
    tmp_path: Path, simulated_edi: Path
) -> None:
    records = _line_records(90.0)[:2] + [
        ("M01", float("nan"), float("nan")),
        ("M02", float("nan"), float("nan")),
        ("M03", float("nan"), float("nan")),
    ]
    sites = _make_sites(tmp_path, simulated_edi, records)

    ordered = sites.ordered("auto")

    assert [s.name for s in ordered] == [r[0] for r in records]
    assert ordered.ordering["applied"] == "input"
    assert ordered.ordering["coordinate_fraction"] == pytest.approx(0.4)


def test_custom_coordinate_fraction_can_enable_partial_line_order(
    tmp_path: Path, simulated_edi: Path
) -> None:
    records = _line_records(90.0)[:2] + [
        ("M01", float("nan"), float("nan")),
        ("M02", float("nan"), float("nan")),
        ("M03", float("nan"), float("nan")),
    ]
    sites = _make_sites(tmp_path, simulated_edi, records)

    ordered = sites.ordered("auto", min_coordinate_fraction=0.4)

    assert ordered.ordering["applied"] == "chainage"
    assert [s.name for s in ordered][-3:] == ["M01", "M02", "M03"]


def test_forced_chainage_retains_missing_sites_stably(
    tmp_path: Path, simulated_edi: Path
) -> None:
    records = [
        ("S03", 35.0, -116.97),
        ("MISS_B", float("nan"), float("nan")),
        ("S01", 35.0, -116.99),
        ("MISS_A", float("nan"), float("nan")),
        ("S02", 35.0, -116.98),
    ]
    sites = _make_sites(tmp_path, simulated_edi, records)

    ordered = sites.ordered("chainage")

    assert [s.name for s in ordered] == ["S01", "S02", "S03", "MISS_B", "MISS_A"]
    assert ordered.ordering["n_coordinates"] == 3


def test_explicit_natural_station_order_handles_mixed_zero_padding(
    tmp_path: Path, simulated_edi: Path
) -> None:
    records = [(n, 35.0, -117.0) for n in ("L10", "L2", "L001", "L11", "L03")]
    sites = _make_sites(tmp_path, simulated_edi, records)

    ordered = sites.ordered("station")

    assert [s.name for s in ordered] == ["L001", "L2", "L03", "L10", "L11"]


def test_explicit_latitude_and_longitude_modes(
    tmp_path: Path, simulated_edi: Path
) -> None:
    records = [
        ("A", 3.0, 20.0),
        ("B", 1.0, 30.0),
        ("C", 2.0, 10.0),
    ]
    sites = _make_sites(tmp_path, simulated_edi, records)

    assert [s.name for s in sites.ordered("lat")] == ["B", "C", "A"]
    assert [s.name for s in sites.ordered("longitude")] == ["C", "A", "B"]


@pytest.mark.parametrize("alias", ["profile", "spatial", "chainage"])
def test_chainage_aliases_are_equivalent(
    tmp_path: Path, simulated_edi: Path, alias: str
) -> None:
    sites = _make_sites(tmp_path, simulated_edi, _line_records(40.0))
    assert [s.name for s in sites.ordered(alias)] == [
        "S01",
        "S02",
        "S03",
        "S04",
        "S05",
    ]


@pytest.mark.parametrize("alias", ["input", "preserve", "none"])
def test_input_aliases_preserve_exact_order(
    tmp_path: Path, simulated_edi: Path, alias: str
) -> None:
    sites = _make_sites(tmp_path, simulated_edi, _line_records(20.0))
    expected = [s.name for s in sites]
    assert [s.name for s in sites.ordered(alias)] == expected


def test_non_inplace_order_does_not_mutate_source(
    tmp_path: Path, simulated_edi: Path
) -> None:
    sites = _make_sites(tmp_path, simulated_edi, _line_records(75.0))
    original_names = [s.name for s in sites]

    ordered = sites.ordered("auto", inplace=False)

    assert ordered is not sites
    assert [s.name for s in sites] == original_names
    assert [s.name for s in ordered] != original_names


def test_inplace_order_returns_and_mutates_same_container(
    tmp_path: Path, simulated_edi: Path
) -> None:
    sites = _make_sites(tmp_path, simulated_edi, _line_records(75.0))

    result = sites.ordered("auto", inplace=True)

    assert result is sites
    assert [s.name for s in sites] == ["S01", "S02", "S03", "S04", "S05"]


def test_invalid_order_mode_has_actionable_error(
    tmp_path: Path, simulated_edi: Path
) -> None:
    sites = _make_sites(tmp_path, simulated_edi, _line_records(90.0))
    with pytest.raises(ValueError, match="by must be one of"):
        sites.ordered("random_guess")


def test_empty_and_single_site_are_safe(
    tmp_path: Path, simulated_edi: Path
) -> None:
    assert Sites([]).ordered("auto").ordering["applied"] == "input"
    one = _make_sites(tmp_path, simulated_edi, [("ONLY", 35.0, -117.0)])
    ordered = one.ordered("auto")
    assert [s.name for s in ordered] == ["ONLY"]
    assert ordered.ordering["applied"] == "input"


def test_ordering_report_is_returned_as_defensive_copy(
    tmp_path: Path, simulated_edi: Path
) -> None:
    ordered = _make_sites(tmp_path, simulated_edi, _line_records(10.0)).ordered("auto")
    report = ordered.ordering
    report["applied"] = "tampered"
    assert ordered.ordering["applied"] == "chainage"


def test_ensure_sites_defaults_to_auto_and_supports_opt_out(
    tmp_path: Path, simulated_edi: Path
) -> None:
    sites = _make_sites(tmp_path, simulated_edi, _line_records(65.0))
    original = [s.name for s in sites]

    automatic = ensure_sites(sites)
    preserved = ensure_sites(sites, order_by="input")

    assert [s.name for s in automatic] == ["S01", "S02", "S03", "S04", "S05"]
    assert automatic.ordering["applied"] == "chainage"
    assert [s.name for s in preserved] == original


def test_static_shift_shared_order_helper_respects_sites_policy(
    tmp_path: Path, simulated_edi: Path
) -> None:
    sites = _make_sites(tmp_path, simulated_edi, _line_records(115.0))

    automatic = _order_sites(sites, "auto")
    preserved = _order_sites(sites, "input")

    assert [s.name for s in automatic] == ["S01", "S02", "S03", "S04", "S05"]
    assert [s.name for s in preserved] == [s.name for s in sites]


def test_real_l22plt_directory_orders_by_chainage() -> None:
    data = Path(__file__).parents[3] / "data" / "AMT" / "WILLY_DATA" / "L22PLT"
    if not data.exists():
        pytest.skip("L22PLT integration data is not available")

    sites = ensure_sites(data)
    names = [s.name for s in sites]

    assert len(names) == 25
    assert names[:3] == ["22-1BF", "22-2VF", "22-3U"]
    assert names[-3:] == ["22-23VF", "22-24BF", "22-025AF"]
    assert sites.ordering["applied"] == "chainage"
    assert sites.ordering["linearity"] > 0.999
    assert sites.ordering["cross_track_ratio"] < 0.02


def test_real_parallel_l18_and_l22_are_kept_as_separate_ordered_lines() -> None:
    """A combined load must not interleave two nearby parallel profiles."""

    data = Path(__file__).parents[3] / "data" / "AMT" / "WILLY_DATA"
    line18 = data / "L18PLT"
    line22 = data / "L22PLT"
    if not line18.exists() or not line22.exists():
        pytest.skip("L18PLT/L22PLT integration data is not available")

    expected18 = [s.name for s in ensure_sites(line18)]
    expected22 = [s.name for s in ensure_sites(line22)]
    combined = ensure_sites([line18, line22])
    names = [s.name for s in combined]

    assert len(expected18) == 28
    assert len(expected22) == 25
    assert names == expected18 + expected22
    assert all(name.startswith("18-") for name in names[: len(expected18)])
    assert all(name.startswith("22-") for name in names[len(expected18) :])
    assert combined.ordering["applied"] == "chainage_by_line"
    assert combined.ordering["n_profiles"] == 2
    assert combined.ordering["cross_track_gap_ratio"] > 0.25


def test_real_l22_pseudosection_keeps_sites_chainage_order() -> None:
    """Dataframe pivoting must not restore lexical station ordering."""

    data = Path(__file__).parents[3] / "data" / "AMT" / "WILLY_DATA" / "L22PLT"
    if not data.exists():
        pytest.skip("L22PLT integration data is not available")
    sites = ensure_sites(data)

    ax = pseudosection(sites, quantity="rho_xy", topo=False, dark=False)
    labels = [tick.get_text() for tick in ax.get_xticklabels()]
    ax.figure.clear()

    assert labels == [site.name for site in sites]
    assert labels[:3] == ["22-1BF", "22-2VF", "22-3U"]
    assert labels[-3:] == ["22-23VF", "22-24BF", "22-025AF"]
