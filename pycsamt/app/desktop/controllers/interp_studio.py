# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
What the Interpretation Studio offers, declared once (Qt-free).

The studio has six tabs.  Each tab has an optional **step** (a computation
on :class:`~pycsamt.app.desktop.controllers.interp_controller.
InterpController`), its settings (:class:`~pycsamt.app.desktop.controllers.
inversion_engines.Field`, so the window builds forms the same way the
Inversion Studio does) and its **views** (controller ``plot_*`` methods).

Every view names what it needs (a model, EDI/XML sites, the step's
result, ...).  :func:`view_status` checks that against the controller, so
the window can mark each view ready or say what is missing *before* the
user asks for it.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from pycsamt.app.desktop.controllers.inversion_engines import Field

__all__ = [
    "NEEDS",
    "TABS",
    "StudioTab",
    "StudioView",
    "tab",
    "view_status",
]


@dataclass(frozen=True)
class StudioView:
    method: str  # InterpController.plot_* method
    label: str
    help: str = ""
    needs: tuple[str, ...] = ()
    per_station: bool = False  # takes the header's station


@dataclass(frozen=True)
class StudioTab:
    key: str
    label: str
    views: tuple[StudioView, ...]
    step: str = ""  # InterpController method run by the tab's button
    step_label: str = ""
    step_help: str = ""
    fields: tuple[Field, ...] = field(default_factory=tuple)


# requirement -> (test on the controller, what to tell the user)
NEEDS = {
    "model": (lambda c: c.state.model is not None,
              "Load a resistivity model (Model ▸ Load)."),
    "sites": (lambda c: c.state.sites is not None,
              "Load EDI/XML data in the main window."),
    "logs": (lambda c: bool(c.state.strat_logs),
             "Run Classify & calibrate (Geology)."),
    "calibrated": (lambda c: c.state.original_model is not None,
                   "Add boreholes, then run Classify & calibrate."),
    "boreholes": (lambda c: bool(c.state.boreholes),
                  "Add boreholes (Evidence ▸ Boreholes)."),
    "structure": (
        lambda c: c.state.structural_model is not None and bool(
            c.state.structural_model.faults
            or c.state.structural_model.planar
            or c.state.structural_model.linear),
        "Add structural measurements (Evidence ▸ Structure)."),
    "hydro": (lambda c: c.state.hydro_result is not None,
              "Run the hydrology estimate (Hydrology)."),
    "mc": (lambda c: c.state.mc_result is not None,
           "Run Monte Carlo (Uncertainty)."),
    "timelapse": (lambda c: len(c.state.timelapse_surveys) >= 2,
                  "Add a repeat survey (Evidence ▸ Monitoring)."),
    "fusion": (lambda c: c.state.fusion_model is not None,
               "Load a second model and run Fuse models (Monitoring)."),
}


def view_status(ctrl, view: StudioView) -> tuple[bool, str]:
    """``(ready, reason)``: whether *view* can draw now, and if not why."""
    for need in view.needs:
        ok, why = NEEDS[need]
        try:
            if not ok(ctrl):
                return False, why
        except Exception:
            return False, why
    return True, ""


V = StudioView

_PETRO = (
    Field("petro_model", "Petrophysics", "choice", "archie",
          choices=(("archie", "Archie (clean sands)"),
                   ("waxman_smits", "Waxman–Smits (shaly sands)")),
          help="Resistivity → saturation law."),
    Field("rho_w", "Pore-water ρw", "float", 10.0, 0.01, 1e4, 1.0, 2,
          unit="Ω·m", help="Formation-water resistivity."),
    Field("phi", "Porosity φ", "float", 0.35, 0.01, 0.8, 0.01, 2,
          help="Prior porosity (fraction)."),
    Field("m", "Cementation m", "float", 2.0, 1.0, 3.5, 0.1, 2),
    Field("n", "Saturation n", "float", 2.0, 1.0, 3.5, 0.1, 2),
    Field("a", "Tortuosity a", "float", 1.0, 0.3, 2.0, 0.05, 2,
          advanced=True),
    Field("sigma_s", "Clay conductance σs", "float", 0.0, 0.0, 10.0, 0.01,
          3, unit="S/m", advanced=True,
          help="Waxman–Smits surface conduction (ignored by Archie)."),
    Field("d50_m", "Grain size d50", "float", 5e-4, 1e-6, 0.1, 1e-4, 5,
          unit="m", advanced=True,
          help="Median grain size for the Kozeny–Carman K estimate."),
)

_MC = (
    Field("n_samples", "Samples", "int", 200, 20, 5000, 50,
          help="Monte Carlo realisations."),
    Field("rho_w_lo", "ρw min", "float", 5.0, 0.01, 1e4, 1.0, 2,
          unit="Ω·m"),
    Field("rho_w_hi", "ρw max", "float", 50.0, 0.01, 1e4, 1.0, 2,
          unit="Ω·m"),
    Field("phi_lo", "φ min", "float", 0.15, 0.01, 0.8, 0.01, 2),
    Field("phi_hi", "φ max", "float", 0.45, 0.01, 0.8, 0.01, 2),
    Field("m_lo", "m min", "float", 1.5, 1.0, 3.5, 0.1, 2, advanced=True),
    Field("m_hi", "m max", "float", 2.5, 1.0, 3.5, 0.1, 2, advanced=True),
    Field("n_lo", "n min", "float", 1.5, 1.0, 3.5, 0.1, 2, advanced=True),
    Field("n_hi", "n max", "float", 2.5, 1.0, 3.5, 0.1, 2, advanced=True),
)

_FUSION = (
    Field("primary_max_depth", "Primary to", "float", 0.0, 0.0, 1e6, 50.0,
          0, unit="m", auto_zero=True,
          help="Deepest depth the (shallow) loaded model contributes."),
    Field("secondary_min_depth", "Second from", "float", 0.0, 0.0, 1e6,
          50.0, 0, unit="m", auto_zero=True,
          help="Shallowest depth the second (deep) model contributes."),
    Field("blend", "Blend", "choice", "linear",
          choices=(("linear", "Linear ramp"), ("sigmoid", "Sigmoid"),
                   ("rms_weighted", "RMS-weighted"))),
)

TABS: tuple[StudioTab, ...] = (
    StudioTab(
        "geology", "Geology",
        step="run_geological", step_label="Classify && calibrate",
        step_help="Classify the model with the rock database into "
                  "stratigraphic logs; with boreholes, calibrate the model "
                  "against them.",
        views=(
            V("plot_model_summary", "Resistivity model",
              "The loaded 2-D ρ section.", ("model",)),
            V("plot_depth_profile", "Depth profile",
              "ρ against depth under the station (with its borehole).",
              ("model",), per_station=True),
            V("plot_strat_log", "Stratigraphic log",
              "Lithology column under the station.", ("logs",),
              per_station=True),
            V("plot_fence_diagram", "Fence diagram",
              "Pseudo-stratigraphic section from every station's log.",
              ("logs",)),
            V("plot_calibrated_model", "Calibrated model",
              "Model before / after borehole calibration + misfit.",
              ("calibrated",)),
            V("plot_borehole_map", "Borehole positions",
              "Stations and boreholes along the profile.",
              ("model", "boreholes")),
            V("plot_borehole_fence", "Borehole fence",
              "The boreholes' logs side by side.", ("boreholes",)),
            V("plot_rock_db", "Rock database",
              "Resistivity ranges used for classification."),
        )),
    StudioTab(
        "structure", "Structure",
        views=(
            V("plot_structural_section", "Structural section",
              "Faults, planar and linear measurements on the profile.",
              ("structure",)),
            V("plot_strike_rose", "Strike rose",
              "Geoelectric strike directions.", ("sites",)),
            V("plot_strike_profile", "Strike profile",
              "Strike angle along the profile.", ("sites",)),
            V("plot_dimensionality", "Dimensionality",
              "1-D / 2-D / 3-D classification.", ("sites",)),
            V("plot_pt_section", "Phase-tensor section",
              "Phase-tensor ellipses pseudo-section.", ("sites",)),
            V("plot_induction_arrows", "Induction arrows",
              "Tipper arrows along the profile.", ("sites",)),
        )),
    StudioTab(
        "hydrology", "Hydrology",
        step="run_hydro", step_label="Estimate hydrology",
        step_help="Resistivity → saturation, hydraulic conductivity K, "
                  "water table and aquifers with the petrophysics below.",
        fields=_PETRO,
        views=(
            V("plot_K_map", "Hydraulic conductivity",
              "2-D K section.", ("hydro",)),
            V("plot_Sw_map", "Water saturation",
              "2-D Sw section.", ("hydro",)),
            V("plot_water_table", "Water table",
              "Water-table depth along the profile.", ("hydro",)),
            V("plot_aquifer_zones", "Aquifer zones",
              "Aquifer intervals with confidence.", ("model",)),
            V("plot_transmissivity", "Transmissivity",
              "Transmissivity profile.", ("hydro",)),
            V("plot_aquifer_char", "Aquifer characterisation",
              "K, T and Dar-Zarrouk panels.", ("hydro",)),
            V("plot_petrophys_xplot", "Petrophysical cross-plot",
              "ρ against porosity / saturation.", ("hydro",)),
        )),
    StudioTab(
        "monitoring", "Monitoring",
        step="run_fusion", step_label="Fuse models",
        step_help="Merge the loaded (shallow) model with a second (deep) "
                  "model over their overlap.",
        fields=_FUSION,
        views=(
            V("plot_timelapse_change", "Δρ between surveys",
              "Resistivity change of every repeat survey.",
              ("timelapse",)),
            V("plot_timelapse_sat", "ΔSaturation",
              "Saturation change implied by the ρ change.",
              ("timelapse",)),
            V("plot_timelapse_wt", "Water-table change",
              "Water-table rise / fall from the baseline.",
              ("timelapse",)),
            V("plot_fused_model", "Fused model",
              "The two models merged over depth.", ("fusion",)),
        )),
    StudioTab(
        "uncertainty", "Uncertainty",
        step="run_monte_carlo", step_label="Run Monte Carlo",
        step_help="Propagate petrophysical uncertainty (ranges below) "
                  "through the hydrology estimate.",
        fields=_MC,
        views=(
            V("plot_mc_K_section", "K uncertainty",
              "P10 / P50 / P90 hydraulic conductivity.", ("mc",)),
            V("plot_mc_wt_profile", "Water-table uncertainty",
              "Water-table depth bands.", ("mc",)),
            V("plot_mc_histograms", "Parameter histograms",
              "Ensemble distributions.", ("mc",)),
        )),
    StudioTab(
        "diagnostics", "Diagnostics",
        views=(
            V("plot_depth_coverage", "Depth coverage",
              "Bostick depth reach per station.", ("sites",)),
            V("plot_bostick_depths", "Bostick depths",
              "Penetration depth per station and frequency.", ("sites",)),
            V("plot_gradient_section", "Gradient imaging",
              "Joint spatial-frequency gradient.", ("sites",)),
            V("plot_snr_section", "SNR section",
              "Signal-to-noise pseudo-section.", ("sites",)),
            V("plot_z_invariants", "Z invariants",
              "Determinant, skew and SSQ invariants.", ("sites",)),
            V("plot_composite_section", "Composite section",
              "Several metrics in one pseudo-section.", ("sites",)),
            V("plot_anisotropy", "Anisotropy",
              "Apparent anisotropy.", ("sites",)),
            V("plot_bode", "ρa & phase (Bode)",
              "Apparent resistivity and phase against period.", ("sites",)),
            V("plot_phasor_wheel", "Phasor wheel",
              "Impedance hodogram.", ("sites",)),
            V("plot_argand", "Argand diagram",
              "Re Z against Im Z.", ("sites",)),
            V("plot_mohr_circles", "Mohr circles",
              "Impedance Mohr circles.", ("sites",)),
            V("plot_pt_ternary", "Phase-tensor ternary",
              "1-D / 2-D / 3-D ternary.", ("sites",)),
            V("plot_distortion_radar", "Distortion radar",
              "Static shift, twist and shear.", ("sites",)),
            V("plot_period_clock", "Period clock",
              "Phase-tensor clock face.", ("sites",)),
        )),
)


def tab(key: str) -> StudioTab:
    return next(t for t in TABS if t.key == key)
