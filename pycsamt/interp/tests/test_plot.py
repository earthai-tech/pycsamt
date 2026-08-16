# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.interp.plot — all Plot* classes.

Uses the non-interactive Agg backend.  Each test builds minimal but
physically-consistent fixtures (reusing the same patterns as the other
pycsamt.interp test modules), calls .plot(), asserts a matplotlib Figure
comes back, and closes it.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

from unittest import mock  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pytest  # noqa: E402

from pycsamt.interp import plot as ip  # noqa: E402
from pycsamt.interp._base import ResistivityModel  # noqa: E402
from pycsamt.geology.borehole import Borehole, Interval  # noqa: E402
from pycsamt.interp.hydromodel import (  # noqa: E402
    EMHydroModel,
    PetrophysicalConfig,
)
from pycsamt.geology.lithology import (  # noqa: E402
    RockDatabase,
    StratigraphicLog,
)
from pycsamt.interp.petrophysics import (  # noqa: E402
    ArchieModel,
    WaxmanSmitsModel,
)
from pycsamt.interp.timelapse import TimeLapseEM  # noqa: E402
from pycsamt.interp.uncertainty import (  # noqa: E402
    MonteCarloHydro,
    UncertaintyBounds,
)

# ─────────────────────────────────────────────────────────────────────────────
# Shared fixtures
# ─────────────────────────────────────────────────────────────────────────────


def _model(rms=1.5, station_names=("S1", "S2")):
    rho = np.log10(
        np.array(
            [
                [800.0, 800.0],
                [600.0, 800.0],
                [15.0, 800.0],
                [10.0, 800.0],
                [8.0, 800.0],
            ]
        )
    )
    return ResistivityModel.from_array(
        rho,
        x_centers=np.array([0.0, 500.0]),
        z_centers=np.array([5.0, 15.0, 30.0, 60.0, 100.0]),
        station_names=list(station_names) if station_names else None,
        method="TDEM",
        rms=rms,
    )


def _cfg():
    return PetrophysicalConfig(
        petro=ArchieModel(m=1.8, n=2.0), rho_w=20.0, porosity_prior=0.25
    )


def _hydro_result():
    return EMHydroModel(_model(), _cfg(), method_tag="TDEM").fit()


def _log(name="S1", x=0.0):
    db = RockDatabase.default()
    z = np.array([5.0, 15.0, 30.0, 60.0, 100.0])
    rho = np.log10(np.array([800.0, 600.0, 15.0, 10.0, 8.0]))
    return StratigraphicLog.from_column(name, x, z, rho, db=db)


def _uncertainty_result(n_samples=6):
    bounds = UncertaintyBounds(rho_w_range=(15.0, 25.0), m_range=(1.6, 2.0))
    mc = MonteCarloHydro(_model(), _cfg(), bounds, n_samples=n_samples, seed=0)
    return mc.run()


def _model_wt_detected_both_stations():
    """Both columns show a resistive->conductive transition, so the water
    table is detected everywhere -- needed for histogram tests whose
    Gaussian-synthesis path breaks on an all-NaN (never detected) column."""
    rho = np.log10(
        np.array(
            [
                [800.0, 900.0],
                [600.0, 700.0],
                [15.0, 20.0],
                [10.0, 12.0],
                [8.0, 9.0],
            ]
        )
    )
    return ResistivityModel.from_array(
        rho,
        x_centers=np.array([0.0, 500.0]),
        z_centers=np.array([5.0, 15.0, 30.0, 60.0, 100.0]),
        station_names=["S1", "S2"],
        method="TDEM",
        rms=1.5,
    )


def _uncertainty_result_wt_detected(n_samples=6):
    bounds = UncertaintyBounds(rho_w_range=(15.0, 25.0), m_range=(1.6, 2.0))
    mc = MonteCarloHydro(
        _model_wt_detected_both_stations(),
        _cfg(),
        bounds,
        n_samples=n_samples,
        seed=0,
    )
    return mc.run()


def _timelapse():
    def _survey(rho_linear):
        rho = np.log10(np.array(rho_linear))
        return ResistivityModel.from_array(
            rho,
            x_centers=np.array([0.0, 500.0]),
            z_centers=np.array([5.0, 15.0, 30.0, 60.0, 100.0]),
            method="TDEM",
        )

    dry = _survey([[800, 800], [700, 800], [600, 800], [500, 800], [400, 800]])
    wet = _survey([[500, 800], [200, 800], [50, 800], [20, 800], [10, 800]])
    recharge = _survey([[300, 800], [100, 800], [30, 800], [10, 800], [5, 800]])
    return TimeLapseEM([dry, wet, recharge], labels=["dry", "wet", "recharge"])


def _borehole():
    return Borehole(
        "Bo1",
        x=0.0,
        intervals=[
            Interval(top=0.0, bottom=20.0, lithology="Clay", resistivity=10.0),
            Interval(
                top=20.0,
                bottom=100.0,
                lithology="Sandstone",
                resistivity=800.0,
            ),
        ],
    )


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


# ─────────────────────────────────────────────────────────────────────────────
# _hatch_for
# ─────────────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "lithology,expected",
    [
        ("Granite (fresh)", "+"),
        ("Igneous (basement)", "+"),
        ("Basement complex", "+"),
        ("Fractured zone", "x"),
        ("Fault gouge", "x"),
        ("Aquifer", "o"),
        ("Water-bearing sand", "o"),
        ("Clay", "---"),
        ("Shale", "---"),
        ("Sand (wet)", "..."),
        ("Alluvium (wet)", "..."),
        ("Basalt (fresh)", "///"),
        ("Gabbro", "///"),
        ("Limestone", r"\\\\"),
        ("Dolomite", r"\\\\"),
        ("Marble", r"\\\\"),
        ("Schist", "|||"),
        ("Gneiss", "|||"),
        ("Quartzite", "|||"),
        ("Sulfide ore body", "**"),
        ("Graphite / coal", "**"),
        ("Unknown rock", ""),
    ],
)
def test_hatch_for(lithology, expected):
    assert ip._hatch_for(lithology) == expected


def test_require_mpl_returns_modules():
    m, p = ip._require_mpl()
    assert m is matplotlib
    assert p is plt


def test_require_mpl_raises_when_matplotlib_missing():
    import sys

    with mock.patch.dict(sys.modules, {"matplotlib": None}):
        with pytest.raises(ImportError, match="matplotlib is required"):
            ip._require_mpl()


def test_cmap_with_bad_sets_nan_color():
    cmap = ip._cmap_with_bad("viridis", "red")
    rgba = cmap(np.nan)
    np.testing.assert_array_equal(rgba, cmap.get_bad())


# ─────────────────────────────────────────────────────────────────────────────
# PlotStratigraphicLog
# ─────────────────────────────────────────────────────────────────────────────


def test_plot_stratigraphic_log():
    fig = ip.PlotStratigraphicLog(_log()).plot()
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_stratigraphic_log_custom_title_and_kws():
    fig = ip.PlotStratigraphicLog(
        _log(), title="Custom", annotation_kws={"fontsize": 10}
    ).plot()
    assert fig._suptitle.get_text() == "Custom"


# ─────────────────────────────────────────────────────────────────────────────
# PlotFenceDiagram
# ─────────────────────────────────────────────────────────────────────────────


def test_plot_fence_diagram_multi():
    logs = [_log("S1", 0.0), _log("S2", 500.0)]
    fig = ip.PlotFenceDiagram(logs).plot()
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_fence_diagram_single_log():
    fig = ip.PlotFenceDiagram([_log("S1", 0.0)]).plot()
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_fence_diagram_empty_raises():
    with pytest.raises(ValueError, match="No logs"):
        ip.PlotFenceDiagram([]).plot()


def test_plot_fence_diagram_max_depth_truncates():
    logs = [_log("S1", 0.0), _log("S2", 500.0)]
    fig = ip.PlotFenceDiagram(logs, max_depth=20.0).plot()
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_fence_diagram_legend_matches_distinct_lithologies():
    logs = [_log("S1", 0.0), _log("S2", 500.0)]
    fig = ip.PlotFenceDiagram(logs).plot()
    legend = fig.legends[0]
    lithologies = {log.lithology for log in logs[0].layers + logs[1].layers}
    assert {t.get_text() for t in legend.get_texts()} == lithologies


def test_plot_fence_diagram_legend_can_be_disabled():
    logs = [_log("S1", 0.0)]
    fig = ip.PlotFenceDiagram(logs, show_legend=False).plot()
    assert fig.legends == []


def test_plot_fence_diagram_elevation_strip_adds_axes():
    logs = [_log("S1", 0.0), _log("S2", 500.0), _log("S3", 1000.0)]
    elev = np.array([100.0, 120.0, 90.0])
    fig = ip.PlotFenceDiagram(logs, elevation_m=elev).plot()
    # 3 station panels + 1 elevation strip
    assert len(fig.axes) == 4


def test_plot_fence_diagram_without_elevation_has_no_strip():
    logs = [_log("S1", 0.0), _log("S2", 500.0)]
    fig = ip.PlotFenceDiagram(logs).plot()
    assert len(fig.axes) == 2


def test_plot_fence_diagram_elevation_length_mismatch_raises():
    logs = [_log("S1", 0.0), _log("S2", 500.0)]
    with pytest.raises(ValueError, match="elevation_m has"):
        ip.PlotFenceDiagram(logs, elevation_m=np.array([100.0])).plot()


# ─────────────────────────────────────────────────────────────────────────────
# PlotBoreholeFence
# ─────────────────────────────────────────────────────────────────────────────


def _borehole(name="A1", x=300.0):
    return Borehole(
        name,
        x=x,
        intervals=[
            Interval(top=0.0, bottom=65.0, lithology="Overburden", resistivity=80.0),
            Interval(top=65.0, bottom=321.0, lithology="Weathered basement", resistivity=300.0),
            Interval(top=321.0, bottom=450.0, lithology="Fresh basement", resistivity=3000.0),
        ],
    )


def _regional_db():
    from pycsamt.geology.lithology import RockEntry

    return RockDatabase([
        RockEntry(name="Overburden", rho_min=30, rho_max=180, color="#D4AC0D"),
        RockEntry(name="Weathered basement", rho_min=180, rho_max=900, color="#A9780C"),
        RockEntry(name="Fresh basement", rho_min=900, rho_max=8000, color="#4A4A4A"),
    ])


def test_plot_borehole_fence_multi():
    boreholes = [_borehole("A1", 300.0), _borehole("A2", 2100.0)]
    fig = ip.PlotBoreholeFence(boreholes, db=_regional_db()).plot()
    assert isinstance(fig, matplotlib.figure.Figure)
    assert len(fig.axes) >= 2


def test_plot_borehole_fence_empty_raises():
    with pytest.raises(ValueError, match="No boreholes"):
        ip.PlotBoreholeFence([]).plot()


def test_plot_borehole_fence_legend_matches_db_lithologies():
    boreholes = [_borehole("A1", 300.0)]
    fig = ip.PlotBoreholeFence(boreholes, db=_regional_db()).plot()
    legend = fig.legends[0]
    assert {t.get_text() for t in legend.get_texts()} == {
        "Overburden", "Weathered basement", "Fresh basement",
    }


def test_plot_borehole_fence_without_db_uses_lithology_fallback_color():
    boreholes = [_borehole("A1", 300.0)]
    fig = ip.PlotBoreholeFence(boreholes, db=None, show_legend=False).plot()
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_borehole_fence_max_depth_truncates():
    boreholes = [_borehole("A1", 300.0)]
    fig = ip.PlotBoreholeFence(boreholes, max_depth=100.0).plot()
    assert isinstance(fig, matplotlib.figure.Figure)


# ─────────────────────────────────────────────────────────────────────────────
# PlotCalibratedModel
# ─────────────────────────────────────────────────────────────────────────────


def test_plot_calibrated_model_auto_misfit():
    crm = _model()
    nm = _model()
    nm.rho_2d = nm.rho_2d + 0.1
    fig = ip.PlotCalibratedModel(crm, nm).plot()
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_calibrated_model_explicit_misfit():
    crm = _model()
    nm = _model()
    misfit = np.full(crm.rho_2d.shape, 5.0)
    fig = ip.PlotCalibratedModel(crm, nm, misfit_map=misfit).plot()
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_calibrated_model_no_stations():
    crm = _model(station_names=None)
    crm.station_x = np.array([])
    nm = _model(station_names=None)
    nm.station_x = np.array([])
    fig = ip.PlotCalibratedModel(crm, nm).plot()
    assert isinstance(fig, matplotlib.figure.Figure)


# ─────────────────────────────────────────────────────────────────────────────
# PlotHydroSection
# ─────────────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("quantity", ["K", "saturation", "porosity"])
def test_plot_hydro_section_quantities(quantity):
    fig = ip.PlotHydroSection(_hydro_result(), quantity=quantity).plot()
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_hydro_section_invalid_quantity_raises():
    with pytest.raises(ValueError, match="quantity must be"):
        ip.PlotHydroSection(_hydro_result(), quantity="bogus")


def test_plot_hydro_section_no_water_table_and_depth_bounds():
    fig = ip.PlotHydroSection(
        _hydro_result(),
        show_water_table=False,
        depth_min=10.0,
        depth_max=60.0,
        vmin=-8.0,
        vmax=-3.0,
        cmap="plasma",
        style="publication",
        title="Custom K section",
    ).plot()
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_hydro_section_no_stations():
    result = _hydro_result()
    result.resistivity_model.station_x = np.array([])
    fig = ip.PlotHydroSection(result).plot()
    assert isinstance(fig, matplotlib.figure.Figure)


# ─────────────────────────────────────────────────────────────────────────────
# PlotWaterTableProfile
# ─────────────────────────────────────────────────────────────────────────────


def test_plot_water_table_profile_basic():
    fig = ip.PlotWaterTableProfile(_hydro_result()).plot()
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_water_table_profile_with_reference_and_colors():
    fig = ip.PlotWaterTableProfile(
        _hydro_result(),
        reference_depth=25.0,
        color_wt="red",
        color_T="green",
        title="Custom",
        style="dark",
    ).plot()
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_water_table_profile_no_stations():
    result = _hydro_result()
    result.resistivity_model.station_x = np.array([])
    fig = ip.PlotWaterTableProfile(result).plot()
    assert isinstance(fig, matplotlib.figure.Figure)


# ─────────────────────────────────────────────────────────────────────────────
# PlotTimeLapseSection
# ─────────────────────────────────────────────────────────────────────────────


def test_plot_timelapse_section_rho():
    fig = ip.PlotTimeLapseSection(_timelapse(), quantity="rho").plot()
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_timelapse_section_saturation():
    fig = ip.PlotTimeLapseSection(
        _timelapse(), quantity="saturation", petro=ArchieModel()
    ).plot()
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_timelapse_section_invalid_quantity_raises():
    with pytest.raises(ValueError, match="quantity must be"):
        ip.PlotTimeLapseSection(_timelapse(), quantity="bogus")


def test_plot_timelapse_section_saturation_requires_petro():
    with pytest.raises(ValueError, match="petro"):
        ip.PlotTimeLapseSection(_timelapse(), quantity="saturation")


def test_plot_timelapse_section_survey_idx_out_of_range():
    with pytest.raises(IndexError):
        ip.PlotTimeLapseSection(_timelapse(), quantity="rho", survey_idx=99).plot()


def test_plot_timelapse_section_no_water_table_and_custom_vmax():
    fig = ip.PlotTimeLapseSection(
        _timelapse(),
        quantity="rho",
        survey_idx=1,
        show_water_table=False,
        vmax=2.0,
        depth_max=60.0,
        cmap="PiYG",
        title="Custom",
    ).plot()
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_timelapse_section_no_stations():
    tl = _timelapse()
    for s in tl.surveys:
        s.station_x = np.array([])
    fig = ip.PlotTimeLapseSection(tl).plot()
    assert isinstance(fig, matplotlib.figure.Figure)


# ─────────────────────────────────────────────────────────────────────────────
# PlotUncertaintySection
# ─────────────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("quantity", ["K", "saturation", "porosity"])
def test_plot_uncertainty_section_quantities(quantity):
    fig = ip.PlotUncertaintySection(_uncertainty_result(), quantity=quantity).plot()
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_uncertainty_section_invalid_quantity_raises():
    with pytest.raises(ValueError, match="quantity must be"):
        ip.PlotUncertaintySection(_uncertainty_result(), quantity="bogus")


def test_plot_uncertainty_section_custom_limits_and_no_wt():
    fig = ip.PlotUncertaintySection(
        _uncertainty_result(),
        show_water_table=False,
        vmin_p50=-8.0,
        vmax_p50=-3.0,
        vmax_spread=2.0,
        depth_min=10.0,
        depth_max=60.0,
        cmap_p50="plasma",
        cmap_spread="Greys",
        title="Custom",
    ).plot()
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_uncertainty_section_no_stations():
    unc = _uncertainty_result()
    unc.resistivity_model.station_x = np.array([])
    fig = ip.PlotUncertaintySection(unc).plot()
    assert isinstance(fig, matplotlib.figure.Figure)


# ─────────────────────────────────────────────────────────────────────────────
# PlotUncertaintyProfile
# ─────────────────────────────────────────────────────────────────────────────


def test_plot_uncertainty_profile_basic():
    fig = ip.PlotUncertaintyProfile(_uncertainty_result()).plot()
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_uncertainty_profile_with_reference_and_colors():
    fig = ip.PlotUncertaintyProfile(
        _uncertainty_result(),
        reference_depth=25.0,
        color_wt="red",
        color_T="green",
        title="Custom",
    ).plot()
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_uncertainty_profile_no_stations():
    unc = _uncertainty_result()
    unc.resistivity_model.station_x = np.array([])
    fig = ip.PlotUncertaintyProfile(unc).plot()
    assert isinstance(fig, matplotlib.figure.Figure)


# ─────────────────────────────────────────────────────────────────────────────
# PlotPetrophysicalCrossPlot
# ─────────────────────────────────────────────────────────────────────────────


def test_plot_crossplot_default():
    fig = ip.PlotPetrophysicalCrossPlot(_hydro_result()).plot()
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_crossplot_color_by_depth_no_hs_bounds_linear_rho():
    fig = ip.PlotPetrophysicalCrossPlot(
        _hydro_result(),
        color_by="depth",
        show_hs_bounds=False,
        log_rho=False,
        title="Custom",
    ).plot()
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_crossplot_waxman_smits_petro():
    result = _hydro_result()
    fig = ip.PlotPetrophysicalCrossPlot(
        result, petro=WaxmanSmitsModel(sigma_s=0.01), Sw_for_curve=0.8
    ).plot()
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_crossplot_empty_depth_range_shows_placeholder():
    fig = ip.PlotPetrophysicalCrossPlot(
        _hydro_result(), depth_range=(9000.0, 9999.0)
    ).plot()
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_crossplot_restricted_depth_range():
    fig = ip.PlotPetrophysicalCrossPlot(_hydro_result(), depth_range=(0.0, 40.0)).plot()
    assert isinstance(fig, matplotlib.figure.Figure)


# ─────────────────────────────────────────────────────────────────────────────
# PlotAquiferCharacterization
# ─────────────────────────────────────────────────────────────────────────────


def test_plot_aquifer_characterization_default():
    fig = ip.PlotAquiferCharacterization(_hydro_result()).plot()
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_aquifer_characterization_no_transmissivity_linear_tr():
    fig = ip.PlotAquiferCharacterization(
        _hydro_result(),
        show_transmissivity=False,
        log_TR=False,
        reference_depth=25.0,
        title="Custom",
    ).plot()
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_aquifer_characterization_single_station_no_stations():
    result = _hydro_result()
    result.resistivity_model.station_x = np.array([])
    fig = ip.PlotAquiferCharacterization(result).plot()
    assert isinstance(fig, matplotlib.figure.Figure)


# ─────────────────────────────────────────────────────────────────────────────
# PlotMultiTimeLapseGrid
# ─────────────────────────────────────────────────────────────────────────────


def test_plot_multi_timelapse_grid_rho_all_surveys():
    fig = ip.PlotMultiTimeLapseGrid(_timelapse(), quantity="rho").plot()
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_multi_timelapse_grid_delta_rho_subset():
    fig = ip.PlotMultiTimeLapseGrid(
        _timelapse(), quantity="delta_rho", surveys=[0, 2]
    ).plot()
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_multi_timelapse_grid_delta_saturation():
    fig = ip.PlotMultiTimeLapseGrid(
        _timelapse(), quantity="delta_saturation", petro=ArchieModel()
    ).plot()
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_multi_timelapse_grid_single_panel():
    fig = ip.PlotMultiTimeLapseGrid(_timelapse(), quantity="rho", surveys=[1]).plot()
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_multi_timelapse_grid_invalid_quantity_raises():
    with pytest.raises(ValueError, match="quantity must be"):
        ip.PlotMultiTimeLapseGrid(_timelapse(), quantity="bogus")


def test_plot_multi_timelapse_grid_delta_saturation_requires_petro():
    with pytest.raises(ValueError, match="petro"):
        ip.PlotMultiTimeLapseGrid(_timelapse(), quantity="delta_saturation")


def test_plot_multi_timelapse_grid_empty_surveys_raises():
    with pytest.raises(ValueError, match="empty"):
        ip.PlotMultiTimeLapseGrid(_timelapse(), surveys=[]).plot()


def test_plot_multi_timelapse_grid_no_stations():
    tl = _timelapse()
    for s in tl.surveys:
        s.station_x = np.array([])
    fig = ip.PlotMultiTimeLapseGrid(tl).plot()
    assert isinstance(fig, matplotlib.figure.Figure)


# ─────────────────────────────────────────────────────────────────────────────
# PlotResistivityDepthProfile
# ─────────────────────────────────────────────────────────────────────────────


def test_plot_resistivity_depth_profile_from_model_int_station():
    fig = ip.PlotResistivityDepthProfile(_model(), station=0).plot()
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_resistivity_depth_profile_from_hydro_result_named_station():
    fig = ip.PlotResistivityDepthProfile(
        _hydro_result(), station="S1", show_zones=True
    ).plot()
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_resistivity_depth_profile_with_borehole_panel():
    fig = ip.PlotResistivityDepthProfile(
        _hydro_result(), station="S2", borehole=_borehole(), log_rho=False
    ).plot()
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_resistivity_depth_profile_resistive_basement_zone_color():
    """Sw < 0.25 in the vadose zone selects the 'resistive basement' shade."""
    rho = np.log10(
        np.array(
            [
                [9000.0, 800.0],
                [3000.0, 800.0],
                [15.0, 800.0],
                [10.0, 800.0],
                [8.0, 800.0],
            ]
        )
    )
    model = ResistivityModel.from_array(
        rho,
        x_centers=np.array([0.0, 500.0]),
        z_centers=np.array([5.0, 15.0, 30.0, 60.0, 100.0]),
        station_names=["S1", "S2"],
        method="TDEM",
        rms=1.5,
    )
    result = EMHydroModel(model, _cfg(), method_tag="TDEM").fit()
    fig = ip.PlotResistivityDepthProfile(result, station=0).plot()
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_resistivity_depth_profile_borehole_with_missing_bounds():
    """An interval-like object missing top/bottom is skipped without error."""

    class _FakeInterval:
        top = None
        bottom = None
        lithology = "unknown"
        color = "0.7"

    class _FakeBorehole:
        intervals = [_FakeInterval()]

    fig = ip.PlotResistivityDepthProfile(
        _model(), station=0, borehole=_FakeBorehole()
    ).plot()
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_resistivity_depth_profile_no_zones_depth_max():
    fig = ip.PlotResistivityDepthProfile(
        _hydro_result(), station=0, show_zones=False, depth_max=60.0
    ).plot()
    assert isinstance(fig, matplotlib.figure.Figure)


# ─────────────────────────────────────────────────────────────────────────────
# PlotUncertaintyHistogram
# ─────────────────────────────────────────────────────────────────────────────


def test_plot_uncertainty_histogram_water_table_gaussian_synth():
    fig = ip.PlotUncertaintyHistogram(
        _uncertainty_result_wt_detected(), quantity="water_table"
    ).plot()
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_uncertainty_histogram_transmissivity_with_ensembles():
    unc, wt_ens, T_ens = MonteCarloHydro(
        _model(),
        _cfg(),
        UncertaintyBounds(rho_w_range=(15.0, 25.0)),
        n_samples=8,
        seed=0,
    ).run_ensemble()
    # station S2 (index 1) never detects a water table -> wt_ens[:,1] is
    # all-NaN, exercising the "len(raw) < 2 -> hide axis" branch.
    fig = ip.PlotUncertaintyHistogram(
        unc,
        quantity="transmissivity",
        wt_ensemble=wt_ens,
        T_ensemble=T_ens,
    ).plot()
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_uncertainty_histogram_water_table_ensemble_all_nan_hides_axis():
    unc, wt_ens, T_ens = MonteCarloHydro(
        _model(),
        _cfg(),
        UncertaintyBounds(rho_w_range=(15.0, 25.0)),
        n_samples=8,
        seed=0,
    ).run_ensemble()
    fig = ip.PlotUncertaintyHistogram(
        unc, quantity="water_table", wt_ensemble=wt_ens, T_ensemble=T_ens
    ).plot()
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_uncertainty_histogram_invalid_quantity_raises():
    with pytest.raises(ValueError, match="quantity must be"):
        ip.PlotUncertaintyHistogram(_uncertainty_result_wt_detected(), quantity="bogus")


def test_plot_uncertainty_histogram_explicit_stations_and_no_kde_no_pct():
    fig = ip.PlotUncertaintyHistogram(
        _uncertainty_result_wt_detected(),
        stations=["S1", 1],
        show_kde=False,
        show_percentiles=False,
        log_x=False,
        ncols=2,
        title="Custom",
    ).plot()
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_uncertainty_histogram_kde_exception_is_swallowed():
    with mock.patch("scipy.stats.gaussian_kde", side_effect=RuntimeError("boom")):
        fig = ip.PlotUncertaintyHistogram(
            _uncertainty_result_wt_detected(),
            quantity="water_table",
            show_kde=True,
        ).plot()
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_uncertainty_histogram_leftover_axes_hidden():
    """3 stations with ncols=2 leaves a 4th (empty) axis to be hidden."""
    fig = ip.PlotUncertaintyHistogram(
        _uncertainty_result_wt_detected(),
        stations=[0, 1, 0],
        ncols=2,
    ).plot()
    assert isinstance(fig, matplotlib.figure.Figure)


def test_plot_uncertainty_histogram_no_station_names():
    unc = _uncertainty_result_wt_detected()
    unc.resistivity_model.station_names = []
    fig = ip.PlotUncertaintyHistogram(unc).plot()
    assert isinstance(fig, matplotlib.figure.Figure)
