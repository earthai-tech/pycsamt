"""Tests for pycsamt.emtools.resphase_psection.PlotResPhasePseudoSection.

Uses the bundled 28-station WILLY L18PLT AMT line, which carries all four
impedance components.  Skipped when that data is absent from the checkout.
"""

from __future__ import annotations

from pathlib import Path

import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")
import matplotlib.figure  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

from pycsamt.emtools import ensure_sites  # noqa: E402
from pycsamt.emtools.resphase_psection import (  # noqa: E402
    PlotResPhasePseudoSection,
    plot_res_phase_pseudosection,
)

_L18 = Path(__file__).parents[3] / "data" / "AMT" / "WILLY_DATA" / "L18PLT"

pytestmark = pytest.mark.skipif(
    not (_L18.exists() and any(_L18.glob("*.edi"))),
    reason="WILLY L18PLT sample not available",
)


@pytest.fixture(scope="module")
def sites():
    return ensure_sites(str(_L18), recursive=True)


def _panel_axes(fig):
    return [ax for ax in fig.axes if ax.get_images() or ax.collections]


class TestPlotResPhasePseudoSection:
    def test_returns_figure(self, sites):
        fig = PlotResPhasePseudoSection(sites).plot()
        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_auto_detects_four_components(self, sites):
        built = PlotResPhasePseudoSection(sites)._build()
        assert built["components"] == ["xy", "yx", "xx", "yy"]
        plt.close("all")

    def test_component_order_is_preserved(self, sites):
        fig = PlotResPhasePseudoSection(
            sites, components=["yx", "xy"]
        ).plot()
        titles = [
            ax.get_title() for ax in fig.axes if ax.get_title()
        ]
        assert titles[0] == r"$Z_{yx}$"
        assert titles[1] == r"$Z_{xy}$"
        plt.close(fig)

    def test_two_components_two_column_grid(self, sites):
        fig = PlotResPhasePseudoSection(
            sites, components=["xy", "yx"]
        ).plot()
        # 2 res + 2 phase panels + 2 colorbars
        meshes = [ax for ax in fig.axes if ax.collections]
        assert len(meshes) >= 4
        plt.close(fig)

    def test_phase_range_preset_sets_limits(self, sites):
        fig = PlotResPhasePseudoSection(
            sites, components=["xy"], phase_range="-90-90"
        ).plot()
        qm = [
            c
            for ax in fig.axes
            for c in ax.collections
            if hasattr(c, "get_clim")
        ]
        clims = [c.get_clim() for c in qm]
        assert (-90.0, 90.0) in clims
        plt.close(fig)

    def test_phase_range_tuple(self, sites):
        fig = PlotResPhasePseudoSection(
            sites, components=["xy"], phase_range=(0.0, 60.0)
        ).plot()
        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_bad_phase_range_raises(self, sites):
        with pytest.raises(ValueError, match="phase_range"):
            PlotResPhasePseudoSection(
                sites, components=["xy"], phase_range="weird"
            ).plot()

    def test_res_phase_ratio_changes_panel_height(self, sites):
        thin = PlotResPhasePseudoSection(
            sites, components=["xy"], res_phase_ratio=0.5
        ).plot()
        fat = PlotResPhasePseudoSection(
            sites, components=["xy"], res_phase_ratio=2.0
        ).plot()

        def _res_height(fig):
            axes = sorted(
                (ax for ax in fig.axes if ax.collections and ax.get_title()),
                key=lambda a: -a.get_position().y0,
            )
            return axes[0].get_position().height

        assert _res_height(fat) > 1.6 * _res_height(thin)
        plt.close(thin)
        plt.close(fat)

    def test_log_period_axis_label(self, sites):
        fig = PlotResPhasePseudoSection(
            sites, components=["xy"], log_period=True
        ).plot()
        ylabels = {ax.get_ylabel() for ax in fig.axes}
        assert any("log" in lbl for lbl in ylabels)
        plt.close(fig)

    def test_res_range_and_linear_res(self, sites):
        fig = PlotResPhasePseudoSection(
            sites,
            components=["xy"],
            log_res=False,
            res_range=(10.0, 500.0),
        ).plot()
        labels = {ax.get_ylabel() for ax in fig.axes} | {
            cb.get_label() for ax in fig.axes for cb in [ax.yaxis]
        }
        # the resistivity colorbar drops the log10 prefix
        cbar_labels = {
            ax.get_ylabel() for ax in fig.axes if ax.get_ylabel()
        }
        assert any("rho" in lbl and "log" not in lbl for lbl in cbar_labels)
        plt.close(fig)

    def test_multiple_groups_stack(self, sites):
        from pycsamt.site.base import Sites

        items = list(sites)
        groups = {
            "A": Sites(items[: len(items) // 2]),
            "B": Sites(items[len(items) // 2:]),
        }
        fig = PlotResPhasePseudoSection(
            groups, components=["xy", "yx"]
        ).plot()
        texts = {t.get_text() for t in fig.texts}
        assert "A" in texts and "B" in texts
        plt.close(fig)

    def test_period_range_filter(self, sites):
        fig = PlotResPhasePseudoSection(
            sites, components=["xy"], period_range=(1e-3, 1e-1)
        ).plot()
        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_functional_wrapper(self, sites):
        fig = plot_res_phase_pseudosection(sites, components=["xy", "yx"])
        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_draw_into_supplied_axes(self, sites):
        fig, axgrid = plt.subplots(2, 2, figsize=(8, 6))
        out = PlotResPhasePseudoSection(
            sites, components=["xy", "yx"], axes=axgrid
        ).plot()
        assert out is fig
        # the four supplied axes each carry a pcolormesh
        assert all(ax.collections for ax in axgrid.ravel())
        plt.close(fig)

    def test_too_few_supplied_axes_raises(self, sites):
        fig, axgrid = plt.subplots(1, 2)
        with pytest.raises(ValueError, match="axes must supply"):
            PlotResPhasePseudoSection(
                sites, components=["xy", "yx"], axes=axgrid
            ).plot()
        plt.close(fig)

    def test_station_side_none_drops_markers(self, sites):
        with_markers = PlotResPhasePseudoSection(
            sites, components=["xy"], station_side="top"
        ).plot()
        without = PlotResPhasePseudoSection(
            sites, components=["xy"], station_side="none"
        ).plot()

        def _scatter_count(fig):
            return sum(
                len(
                    [
                        c
                        for c in ax.collections
                        if type(c).__name__ == "PathCollection"
                    ]
                )
                for ax in fig.axes
            )

        assert _scatter_count(with_markers) > _scatter_count(without)
        plt.close(with_markers)
        plt.close(without)

    def test_bad_station_side_raises(self, sites):
        with pytest.raises(ValueError, match="station_side"):
            PlotResPhasePseudoSection(sites, station_side="left")

    def test_panel_labels(self, sites):
        from pycsamt.site.base import Sites

        items = list(sites)
        groups = {
            "north": Sites(items[: len(items) // 2]),
            "south": Sites(items[len(items) // 2:]),
        }
        fig = PlotResPhasePseudoSection(
            groups, components=["xy"], panel_labels=True
        ).plot()
        texts = {
            t.get_text() for ax in fig.axes for t in ax.texts
        }
        assert "(a)" in texts and "(b)" in texts
        plt.close(fig)

    def test_share_period_false_gives_per_group_grids(self, sites):
        from pycsamt.site.base import Sites

        items = list(sites)
        groups = {
            "a": Sites(items[: len(items) // 2]),
            "b": Sites(items[len(items) // 2:]),
        }
        built = PlotResPhasePseudoSection(
            groups, components=["xy"], share_period=False
        )._build()
        grids = [g["period_grid"] for g in built["groups"]]
        assert len(grids) == 2
        # identical station physics here, so the windows match, but the
        # per-group grid path is exercised and returns finite grids
        assert all(g.size > 4 and (g > 0).all() for g in grids)
        plt.close("all")
