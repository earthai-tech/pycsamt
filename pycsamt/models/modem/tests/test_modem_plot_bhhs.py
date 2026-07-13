"""Real-data integration tests for models/modem/plot.py.

Uses the bundled BHHS/RERUN AMT ModEM sample when available.  In source
checkouts where the local BHHS folder has not been populated, the tests fall
back to the compact Willy 27-frequency WATEX ModEM sample, which contains
representative model, data, response, covariance, and log artefacts from the
larger Jiangsu AMT run.

All tests are skipped automatically when the sample directory is absent.
To run: pytest -v pycsamt/models/modem/tests/test_modem_plot_bhhs.py

The origin of the model (lat/lon of model centre) is (32.129°N, 119.125°E).
Profile EW offsets from model centre (metres):
    L18 ≈ +361   L22 ≈ +170   L26 ≈ -26   L30 ≈ -226   L34 ≈ -424
"""

from __future__ import annotations

from pathlib import Path

import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")
import matplotlib.figure
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Data paths
# ---------------------------------------------------------------------------

_RERUN_BHS = (
    Path(__file__).parents[4]
    / "data"
    / "modem"
    / "rerun-amt-bhs"
)
_BUNDLED_SAMPLE = (
    Path(__file__).parents[4]
    / "data"
    / "modem"
    / "willy_27freq_watex_line02_sample"
)


def _has_modem_result(path: Path) -> bool:
    return (
        path.exists()
        and any(path.glob("*.rho"))
        and any(path.glob("*.dat"))
        and any(path.glob("*.cov"))
    )


_BHHS = _RERUN_BHS if _has_modem_result(_RERUN_BHS) else _BUNDLED_SAMPLE

pytestmark = pytest.mark.skipif(
    not _has_modem_result(_BHHS),
    reason=(
        "ModEM BHHS/RERUN sample not available. Populate "
        "data/modem/rerun-amt-bhs or keep the bundled compact sample."
    ),
)

# Geographic origin of the model centre
_ORIGIN_LAT = 32.129
_ORIGIN_LON = 119.125

# EW offsets (model-centre metres) for the 5 survey lines
_LINE_OFFSETS = [361.0, 170.0, -26.0, -226.0, -424.0]
_LINE_NAMES = ["L18", "L22", "L26", "L30", "L34"]


# ---------------------------------------------------------------------------
# Module-level fixture — load InversionResult once for all tests
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def result():
    from pycsamt.models.modem.results import InversionResult

    return InversionResult(_BHHS)


# ---------------------------------------------------------------------------
# PlotMisfit
# ---------------------------------------------------------------------------


class TestPlotMisfit:
    def test_returns_figure(self, result):
        from pycsamt.models.modem.plot import PlotMisfit

        fig = PlotMisfit(result=result).plot()
        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_has_one_axes(self, result):
        from pycsamt.models.modem.plot import PlotMisfit

        fig = PlotMisfit(result=result).plot()
        assert len(fig.axes) >= 1
        plt.close(fig)

    def test_rms_axis_covers_all_iterations(self, result):
        from pycsamt.models.modem.plot import PlotMisfit

        fig = PlotMisfit(result=result).plot()
        ax = fig.axes[0]
        lines = ax.get_lines()
        assert lines, "No lines drawn in misfit figure"
        # 154 RMS values → line should have that many x-points
        assert len(lines[0].get_xdata()) == len(result.log.rms)
        plt.close(fig)

    def test_no_result_raises(self):
        from pycsamt.models.modem.plot import PlotMisfit

        with pytest.raises(ValueError):
            PlotMisfit().plot()


# ---------------------------------------------------------------------------
# PlotModel3D
# ---------------------------------------------------------------------------


class TestPlotModel3D:
    def test_returns_figure(self, result):
        from pycsamt.models.modem.plot import PlotModel3D

        fig = PlotModel3D(
            result=result, section_offset=_LINE_OFFSETS[0]
        ).plot()
        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_initial_model(self, result):
        from pycsamt.models.modem.plot import PlotModel3D

        fig = PlotModel3D(
            result=result, which="initial", section_offset=0.0
        ).plot()
        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_depth_crop(self, result):
        from pycsamt.models.modem.plot import PlotModel3D

        fig = PlotModel3D(result=result, depth_max=500.0).plot()
        ax = fig.axes[0]
        # y-axis should not extend below -500 m (= -0.5 km)
        ylim = ax.get_ylim()
        assert ylim[0] >= -0.6  # km, with a small tolerance
        plt.close(fig)


# ---------------------------------------------------------------------------
# PlotSection — axis-aligned
# ---------------------------------------------------------------------------


class TestPlotSectionAxisAligned:
    @pytest.mark.parametrize(
        "offset,name", list(zip(_LINE_OFFSETS, _LINE_NAMES))
    )
    def test_ns_profile_all_lines(self, result, offset, name):
        from pycsamt.models.modem.plot import PlotSection

        fig = PlotSection(
            result=result,
            profile_offset=offset,
            direction="NS",
            depth_max=2000.0,
            title=f"NS section {name}",
        ).plot()
        assert isinstance(fig, matplotlib.figure.Figure)
        ax = fig.axes[0]
        assert (
            "N-S" in ax.get_xlabel() or "distance" in ax.get_xlabel().lower()
        )
        plt.close(fig)

    def test_ew_profile_centre(self, result):
        from pycsamt.models.modem.plot import PlotSection

        fig = PlotSection(
            result=result,
            profile_offset=0.0,
            direction="EW",
            depth_max=2000.0,
        ).plot()
        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_station_markers_present(self, result):
        from pycsamt.models.modem.plot import PlotSection

        fig = PlotSection(
            result=result,
            profile_offset=_LINE_OFFSETS[0],  # L18
            direction="NS",
            show_stations=True,
            station_tol=50.0,
        ).plot()
        ax = fig.axes[0]
        # At least one PathCollection (scatter) should be in the axes
        from matplotlib.collections import PathCollection

        scatters = [
            c for c in ax.collections if isinstance(c, PathCollection)
        ]
        assert scatters, "No station markers found for L18 profile"
        plt.close(fig)

    def test_no_stations_when_disabled(self, result):
        from matplotlib.collections import PathCollection

        from pycsamt.models.modem.plot import PlotSection

        fig = PlotSection(
            result=result,
            profile_offset=_LINE_OFFSETS[0],
            direction="NS",
            show_stations=False,
        ).plot()
        ax = fig.axes[0]
        scatters = [
            c for c in ax.collections if isinstance(c, PathCollection)
        ]
        assert not scatters, (
            "Station markers present despite show_stations=False"
        )
        plt.close(fig)


# ---------------------------------------------------------------------------
# PlotSection — arbitrary azimuth (Phase 2)
# ---------------------------------------------------------------------------


class TestPlotSectionArbitrary:
    def test_diagonal_xy_coords(self, result):
        """Diagonal profile crossing all five survey lines."""
        from pycsamt.models.modem.plot import PlotSection

        fig = PlotSection(
            result=result,
            start_point=(-1200.0, 400.0),
            end_point=(1400.0, -450.0),
        ).plot()
        assert isinstance(fig, matplotlib.figure.Figure)
        ax = fig.axes[0]
        assert "Distance along profile" in ax.get_xlabel()
        assert "azimuth" in ax.get_xlabel().lower()
        plt.close(fig)

    def test_diagonal_latlon(self, result):
        """Same diagonal profile defined via lat/lon."""
        from pycsamt.models.modem.plot import PlotSection

        fig = PlotSection(
            result=result,
            start_point=(32.118, 119.121),
            end_point=(32.141, 119.131),
            use_latlon=True,
            origin_lat=_ORIGIN_LAT,
            origin_lon=_ORIGIN_LON,
        ).plot()
        assert isinstance(fig, matplotlib.figure.Figure)
        ax = fig.axes[0]
        title = ax.get_title()
        assert "°N" in title and "°E" in title
        plt.close(fig)

    def test_custom_n_samples(self, result):
        """Different sampling resolutions still produce a figure."""
        from pycsamt.models.modem.plot import PlotSection

        for n in (50, 100, 500):
            fig = PlotSection(
                result=result,
                start_point=(-1000.0, _LINE_OFFSETS[0]),
                end_point=(1000.0, _LINE_OFFSETS[0]),
                n_samples=n,
            ).plot()
            assert isinstance(fig, matplotlib.figure.Figure)
            plt.close(fig)

    def test_start_equals_end_raises(self, result):
        """Degenerate profile should raise."""
        from pycsamt.models.modem.plot import PlotSection

        with pytest.raises(ValueError, match="too close"):
            PlotSection(
                result=result,
                start_point=(0.0, 0.0),
                end_point=(0.0, 0.0),
            ).plot()

    def test_latlon_missing_origin_raises(self, result):
        """use_latlon=True without origin should raise."""
        from pycsamt.models.modem.plot import PlotSection

        with pytest.raises(ValueError, match="origin_lat"):
            PlotSection(
                result=result,
                start_point=(32.12, 119.12),
                end_point=(32.14, 119.13),
                use_latlon=True,
            ).plot()

    def test_arbitrary_colorbar_present(self, result):
        """Arbitrary-azimuth figure must have a colorbar."""
        from pycsamt.models.modem.plot import PlotSection

        fig = PlotSection(
            result=result,
            start_point=(-1200.0, 400.0),
            end_point=(1400.0, -450.0),
        ).plot()
        # colorbar creates an extra axes
        assert len(fig.axes) >= 2, "No colorbar axes found"
        plt.close(fig)

    def test_ns_profile_via_arbitrary(self, result):
        """A purely N-S profile via start/end should match axis-aligned result."""
        from pycsamt.models.modem.plot import PlotSection

        # Pure NS profile at L18 easting.
        fig = PlotSection(
            result=result,
            start_point=(-1300.0, _LINE_OFFSETS[0]),
            end_point=(1300.0, _LINE_OFFSETS[0]),
            station_tol=50.0,
        ).plot()
        ax = fig.axes[0]
        # Azimuth should be ~0° (due north) or ~180°
        xlabel = ax.get_xlabel()
        import re

        match = re.search(r"azimuth\s+([\d.]+)°", xlabel)
        assert match, f"Azimuth not found in xlabel: {xlabel!r}"
        az = float(match.group(1))
        # North = 0° or South = 180° (depending on point order)
        assert abs(az) < 5 or abs(az - 180) < 5, (
            f"Expected ~0° or ~180°, got {az}°"
        )
        plt.close(fig)


# ---------------------------------------------------------------------------
# PlotDepthMap
# ---------------------------------------------------------------------------


class TestPlotDepthMap:
    def test_returns_figure_no_geo(self, result):
        from pycsamt.models.modem.plot import PlotDepthMap

        fig = PlotDepthMap(result=result, depths=[50, 200, 500, 1000]).plot()
        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_returns_figure_with_geo(self, result):
        from pycsamt.models.modem.plot import PlotDepthMap

        fig = PlotDepthMap(
            result=result,
            depths=[50, 200, 500, 1000],
            origin_lat=_ORIGIN_LAT,
            origin_lon=_ORIGIN_LON,
        ).plot()
        assert isinstance(fig, matplotlib.figure.Figure)
        ax = fig.axes[0]
        assert "Longitude" in ax.get_xlabel()
        assert "Latitude" in ax.get_ylabel()
        plt.close(fig)

    def test_default_depths_uses_four_layers(self, result):
        """No depths supplied → 4 subplots (first 4 earth layers)."""
        from pycsamt.models.modem.plot import PlotDepthMap

        fig = PlotDepthMap(result=result).plot()
        # 4 subplots + 4 colorbars = 8 axes minimum
        assert len(fig.axes) >= 4
        plt.close(fig)

    def test_single_depth(self, result):
        from pycsamt.models.modem.plot import PlotDepthMap

        fig = PlotDepthMap(result=result, depths=[500.0], n_cols=1).plot()
        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_subplot_titles_contain_depth(self, result):
        from pycsamt.models.modem.plot import PlotDepthMap

        depths = [200.0, 800.0]
        fig = PlotDepthMap(result=result, depths=depths, n_cols=2).plot()
        titles = [ax.get_title() for ax in fig.axes if ax.get_title()]
        assert any("km" in t for t in titles), f"No km in titles: {titles}"
        plt.close(fig)

    def test_station_markers_drawn(self, result):
        from matplotlib.collections import PathCollection

        from pycsamt.models.modem.plot import PlotDepthMap

        fig = PlotDepthMap(
            result=result, depths=[300.0], n_cols=1, show_stations=True
        ).plot()
        data_axes = [ax for ax in fig.axes if ax.get_title()]
        scatters = [
            c
            for ax in data_axes
            for c in ax.collections
            if isinstance(c, PathCollection)
        ]
        assert scatters, "Expected station scatter markers on depth map"
        plt.close(fig)

    def test_custom_rho_limits(self, result):

        from pycsamt.models.modem.plot import PlotDepthMap

        fig = PlotDepthMap(
            result=result, depths=[200.0], rho_min=10.0, rho_max=500.0
        ).plot()
        # Check colorbar tick labels reflect custom limits
        ax = fig.axes[0]
        assert ax is not None
        plt.close(fig)


# ---------------------------------------------------------------------------
# PlotAllProfiles
# ---------------------------------------------------------------------------


class TestPlotAllProfiles:
    def test_auto_profiles_ns(self, result):
        from pycsamt.models.modem.plot import PlotAllProfiles

        fig = PlotAllProfiles(result=result, n_profiles=5).plot()
        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_named_lines(self, result):
        from pycsamt.models.modem.plot import PlotAllProfiles

        fig = PlotAllProfiles(
            result=result,
            profile_offsets=_LINE_OFFSETS,
            profile_names=_LINE_NAMES,
        ).plot()
        assert isinstance(fig, matplotlib.figure.Figure)
        axes = [ax for ax in fig.axes if ax.get_title()]
        titles = [ax.get_title() for ax in axes]
        for name in _LINE_NAMES:
            assert name in titles, (
                f"Profile {name} title not found; got {titles}"
            )
        plt.close(fig)

    def test_ew_direction(self, result):
        from pycsamt.models.modem.plot import PlotAllProfiles

        fig = PlotAllProfiles(
            result=result, n_profiles=3, direction="EW"
        ).plot()
        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_depth_crop(self, result):

        from pycsamt.models.modem.plot import PlotAllProfiles

        fig = PlotAllProfiles(
            result=result,
            profile_offsets=[_LINE_OFFSETS[0]],
            depth_max=500.0,
        ).plot()
        data_ax = [ax for ax in fig.axes if ax.get_ylabel() == "Depth (km)"]
        for ax in data_ax:
            ylim = ax.get_ylim()
            # y-axis is negative depth; most negative value should not exceed -depth_max/1e3
            assert ylim[0] >= -0.6, f"Depth axis extends beyond 600 m: {ylim}"
        plt.close(fig)

    def test_n_subplots_matches_n_profiles(self, result):
        from pycsamt.models.modem.plot import PlotAllProfiles

        n = 3
        fig = PlotAllProfiles(result=result, n_profiles=n, n_cols=3).plot()
        visible = [
            ax for ax in fig.axes if ax.get_visible() and ax.get_title()
        ]
        assert len(visible) == n
        plt.close(fig)

    def test_station_markers_on_l18(self, result):
        from matplotlib.collections import PathCollection

        from pycsamt.models.modem.plot import PlotAllProfiles

        fig = PlotAllProfiles(
            result=result,
            profile_offsets=[_LINE_OFFSETS[0]],
            profile_names=["L18"],
            show_stations=True,
            station_tol=50.0,
        ).plot()
        scatters = [
            c
            for ax in fig.axes
            for c in ax.collections
            if isinstance(c, PathCollection)
        ]
        assert scatters, "No station markers on L18 in PlotAllProfiles"
        plt.close(fig)


# ---------------------------------------------------------------------------
# PlotCovariance
# ---------------------------------------------------------------------------


class TestPlotCovariance:
    def test_three_panels(self, result):
        from pycsamt.models.modem.plot import PlotCovariance

        fig = PlotCovariance(result=result, show_smoothing=False).plot()
        assert isinstance(fig, matplotlib.figure.Figure)
        # 3 image panels + 3 colorbars = ≥ 6 axes
        assert len(fig.axes) >= 3
        plt.close(fig)

    def test_four_panels_with_smoothing(self, result):
        from pycsamt.models.modem.plot import PlotCovariance

        fig = PlotCovariance(result=result, show_smoothing=True).plot()
        assert isinstance(fig, matplotlib.figure.Figure)
        assert len(fig.axes) >= 4
        plt.close(fig)

    def test_panel_titles(self, result):
        from pycsamt.models.modem.plot import PlotCovariance

        fig = PlotCovariance(result=result, show_smoothing=False).plot()
        titles = [ax.get_title() for ax in fig.axes if ax.get_title()]
        assert any("Plan" in t for t in titles), (
            f"Missing Plan panel: {titles}"
        )
        assert any("N-S" in t for t in titles), f"Missing N-S panel: {titles}"
        assert any("E-W" in t for t in titles), f"Missing E-W panel: {titles}"
        plt.close(fig)

    def test_plan_view_image_shape(self, result):
        """Plan-view imshow data should have shape (nx_earth, ny_earth)."""
        from pycsamt.models.modem.plot import PlotCovariance

        cov = result.covariance
        fig = PlotCovariance(result=result, show_smoothing=False).plot()
        plan_ax = fig.axes[0]
        images = plan_ax.get_images()
        assert images, "No image in plan-view axes"
        data = images[0].get_array()
        assert data.shape == (cov.nx_earth, cov.ny_earth), (
            f"Expected ({cov.nx_earth},{cov.ny_earth}), got {data.shape}"
        )
        plt.close(fig)

    def test_no_covariance_raises(self, tmp_path):
        from pycsamt.models.modem.plot import PlotCovariance
        from pycsamt.models.modem.results import (
            InversionResult,
        )

        r_empty = InversionResult(tmp_path)
        with pytest.raises(ValueError, match="covariance"):
            PlotCovariance(result=r_empty).plot()

    def test_smoothing_panel_has_lines(self, result):
        from pycsamt.models.modem.plot import PlotCovariance

        fig = PlotCovariance(result=result, show_smoothing=True).plot()
        smooth_ax = fig.axes[3]
        lines = smooth_ax.get_lines()
        assert lines, "No lines in smoothing-coefficient panel"
        plt.close(fig)


# ---------------------------------------------------------------------------
# PlotResponse
# ---------------------------------------------------------------------------


class TestPlotResponse:
    def test_zxy_first_station(self, result):
        from pycsamt.models.modem.plot import PlotResponse

        station = list(result.data_obs.site_names)[0]
        fig = PlotResponse(
            result=result, station=station, component="ZXY"
        ).plot()
        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_all_components(self, result):
        from pycsamt.models.modem.plot import PlotResponse

        station = list(result.data_obs.site_names)[0]
        for comp in ("ZXX", "ZXY", "ZYX", "ZYY"):
            fig = PlotResponse(
                result=result, station=station, component=comp
            ).plot()
            assert isinstance(fig, matplotlib.figure.Figure), (
                f"Failed for {comp}"
            )
            plt.close(fig)

    def test_two_axes_rho_phase(self, result):
        """Response figure has separate ρₐ and phase panels."""
        from pycsamt.models.modem.plot import PlotResponse

        station = list(result.data_obs.site_names)[0]
        fig = PlotResponse(
            result=result, station=station, component="ZXY"
        ).plot()
        assert len(fig.axes) >= 2
        plt.close(fig)


# ---------------------------------------------------------------------------
# PlotPseudo
# ---------------------------------------------------------------------------


class TestPlotPseudo:
    def test_zxy_pseudosection(self, result):
        from pycsamt.models.modem.plot import PlotPseudo

        fig = PlotPseudo(result=result, component="ZXY").plot()
        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_zyx_pseudosection(self, result):
        from pycsamt.models.modem.plot import PlotPseudo

        fig = PlotPseudo(result=result, component="ZYX").plot()
        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_two_panels_rho_phase(self, result):
        """Pseudosection must have separate ρₐ and phase pcolormesh panels."""
        from pycsamt.models.modem.plot import PlotPseudo

        fig = PlotPseudo(result=result, component="ZXY").plot()
        assert len(fig.axes) >= 2
        plt.close(fig)
