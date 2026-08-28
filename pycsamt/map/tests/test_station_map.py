# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for station maps."""

from __future__ import annotations

from pycsamt.map import StationMapOptions, build_station_map
from pycsamt.map._core import MapData, StationRecord


class _Z:
    freq = [10.0, 1.0]

    def __init__(self) -> None:
        import numpy as np

        self.resistivity = np.ones((2, 2, 2)) * 100.0
        self.phase = np.ones((2, 2, 2)) * 45.0


class _Edi:
    def __init__(self, station: str) -> None:
        self.station = station
        self.Z = _Z()


class _Sites:
    def as_list(self):
        return [_Edi("S00"), _Edi("S01"), _Edi("S02")]


def test_station_map_options_defaults() -> None:
    opts = StationMapOptions()
    assert opts.overlay == "index"
    assert opts.backend == "plotly"
    assert opts.show_profiles is True
    assert opts.elevation_mode == "markers"


def test_station_map_builds_plotly_figure() -> None:
    data = MapData(
        sites=_Sites(),
        stations=(
            StationRecord("S00", 1.0, 2.0, 10.0, "L1", 0),
            StationRecord("S01", 1.1, 2.1, 20.0, "L1", 1),
            StationRecord("S02", 1.0, 2.2, 30.0, "L1", 2),
        ),
    )
    fig = build_station_map(data, StationMapOptions())
    assert fig.data


def test_station_map_depth_and_density_overlay() -> None:
    data = MapData(
        sites=_Sites(),
        stations=(
            StationRecord("S00", 1.0, 2.0, 10.0, "L1", 0),
            StationRecord("S01", 1.1, 2.1, 20.0, "L1", 1),
            StationRecord("S02", 1.0, 2.2, 30.0, "L1", 2),
        ),
    )
    opts = StationMapOptions(
        overlay="skin_depth",
        frequency=10.0,
        show_contours=True,
        log_color=True,
    )
    fig = build_station_map(data, opts)
    assert len(fig.data) >= 2


def test_station_map_matplotlib_backend() -> None:
    data = MapData(sites=None)
    opts = StationMapOptions(backend="matplotlib")
    fig = build_station_map(data, opts)
    assert hasattr(fig, "savefig")


def test_station_elevation_contour_mode() -> None:
    data = MapData(
        sites=_Sites(),
        stations=(
            StationRecord("S00", 1.0, 2.0, 10.0, "L1", 0),
            StationRecord("S01", 1.1, 2.1, 20.0, "L1", 1),
            StationRecord("S02", 1.0, 2.2, 30.0, "L1", 2),
        ),
    )
    opts = StationMapOptions(
        overlay="elevation",
        elevation_mode="contours",
        backend="matplotlib",
        contour_interp="linear",
    )
    fig = build_station_map(data, opts)
    assert len(fig.axes[0].collections) > 1

    plotly = build_station_map(
        data,
        StationMapOptions(
            overlay="elevation",
            elevation_mode="contours",
            contour_interp="linear",
        ),
    )
    assert len(plotly.layout.map.layers) == 1


def test_station_elevation_mode_rejects_unknown_value() -> None:
    import pytest

    with pytest.raises(ValueError, match="elevation_mode"):
        build_station_map(
            MapData(sites=None),
            StationMapOptions(elevation_mode="raster"),
        )


# ── inversion depth-slice map ──────────────────────────────


def _inv_map_data():
    import numpy as np

    from pycsamt.map._core import MapData, StationRecord

    z = np.linspace(20.0, 600.0, 16)
    sections, recs = {}, []
    for li in range(2):
        n = 6
        sta = np.array([f"L{li}S{i}" for i in range(n)], dtype=object)
        zz, xx = np.meshgrid(z, np.arange(n), indexing="ij")
        # conductive 25 over resistive 800; a 5x block down the right flank
        rho = np.where(zz < 250, 25.0, 800.0)
        rho = np.where(xx > n * 0.6, rho * 5.0, rho)
        for i, s in enumerate(sta):
            recs.append(
                StationRecord(
                    str(s), 6.10 + li * 0.004, 1.20 + i * 0.003,
                    100.0, f"L{li}", len(recs)
                )
            )
        sections[f"L{li}"] = {"z": z, "rho": rho, "stations": sta}
    return MapData(
        sites=None, stations=tuple(recs), metadata={"sections": sections}
    )


def test_depth_rho_overlay_slices_the_inversion_and_draws_a_contour_image():
    import numpy as np

    from pycsamt.map._core import resistivity_at_depth

    data = _inv_map_data()
    vals_shallow = resistivity_at_depth(data, 100.0)
    vals_deep = resistivity_at_depth(data, 400.0)
    # left flank: 25 ohm.m shallow, 800 ohm.m deep -- a real slice, not uniform
    assert np.isclose(vals_shallow["L0S0"], 25.0)
    assert np.isclose(vals_deep["L0S0"], 800.0)
    # right flank is 5x
    assert vals_shallow["L0S5"] > vals_shallow["L0S0"]

    fig = build_station_map(
        data,
        StationMapOptions(overlay="depth_rho", depth=100.0, show_markers=False),
    )
    # the depth slice is rendered as a filled-contour basemap image layer
    layers = fig.layout.map.layers
    assert any(getattr(L, "sourcetype", "") == "image" for L in layers)
    # markers off -> no scattermap line/marker traces except the hidden
    # colour-scale carrier
    assert [t.name for t in fig.data] == ["L0 profile", "L1 profile", "scale"]
    assert fig.layout.paper_bgcolor == "rgba(0,0,0,0)"


def test_depth_rho_colourbar_is_fixed_across_slice_depths():
    import numpy as np

    data = _inv_map_data()

    def scale(depth):
        fig = build_station_map(
            data,
            StationMapOptions(
                overlay="depth_rho", depth=depth, show_markers=False
            ),
        )
        carrier = next(t for t in fig.data if t.name == "scale")
        return carrier.marker.cmin, carrier.marker.cmax

    lo1, hi1 = scale(80.0)
    lo2, hi2 = scale(450.0)
    assert np.isclose(lo1, lo2) and np.isclose(hi1, hi2)
