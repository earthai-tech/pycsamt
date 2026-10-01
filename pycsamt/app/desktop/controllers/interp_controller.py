# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
InterpController — drives all computations inside InterpretationWindow.

Holds the shared interpretation state:
    _sites           — AMT/CSAMT Sites (from MainWindow)
    _model           — ResistivityModel (from inversion or file)
    _boreholes       — list[Borehole] ground-truth data
    _db              — RockDatabase (default or custom CSV)
    _petro_cfg       — PetrophysicalConfig for hydro calculations
    _strat_logs      — list[StratigraphicLog] after classification
    _hydro_result    — EMHydroResult after quantitative run
    _mc_result       — UncertaintyResult after Monte Carlo run

Each generate_*() method returns a matplotlib Figure ready to embed in
an MplCanvas tab.  All methods are pure (no Qt imports).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

# ── Workflow catalogue ─────────────────────────────────────────────────────────
# Maps category name → list of (label, method_name, description)

WORKFLOW_CATALOGUE: dict[str, list[tuple]] = {
    "Setup & Model": [
        ("Model summary", "plot_model_summary", "2-D rho section overview"),
        (
            "Borehole locations",
            "plot_borehole_map",
            "Station + borehole location map",
        ),
        (
            "Depth coverage",
            "plot_depth_coverage",
            "Bostick depth reach per station",
        ),
    ],
    "Geological": [
        (
            "Stratigraphic log",
            "plot_strat_log",
            "Lithology column for selected station",
        ),
        (
            "Fence diagram",
            "plot_fence_diagram",
            "Multi-station pseudo-stratigraphic section",
        ),
        (
            "Calibrated model",
            "plot_calibrated_model",
            "NM vs original model + misfit map",
        ),
        (
            "Rock DB summary",
            "plot_rock_db",
            "Colour-coded resistivity ranges from database",
        ),
    ],
    "Structural": [
        (
            "Structural section",
            "plot_structural_section",
            "Fault traces, planar (strike/dip) and linear (trend/plunge) "
            "measurements on the survey profile",
        ),
    ],
    "Hydrology": [
        (
            "Hydraulic K map",
            "plot_K_map",
            "2-D hydraulic conductivity section",
        ),
        ("Water saturation", "plot_Sw_map", "2-D water saturation section"),
        (
            "Water table profile",
            "plot_water_table",
            "Water-table depth along profile ± confidence",
        ),
        (
            "Aquifer zones",
            "plot_aquifer_zones",
            "Identified aquifer intervals with confidence",
        ),
        (
            "Transmissivity",
            "plot_transmissivity",
            "Transmissivity (m²/s) profile",
        ),
        (
            "Aquifer characterisation",
            "plot_aquifer_char",
            "Multi-panel: K, T, Dar-Zarrouk",
        ),
        (
            "Petrophysical cross-plot",
            "plot_petrophys_xplot",
            "Archie-space: porosity vs Sw vs K",
        ),
    ],
    "Field Constraints": [
        (
            "Constraint misfit",
            "plot_constraint_misfit",
            "Observed vs modelled field data",
        ),
        (
            "Calibration history",
            "plot_calib_history",
            "Optimizer convergence trace",
        ),
    ],
    "EM Diagnostics": [
        (
            "Phase tensor section",
            "plot_pt_section",
            "Phase tensor ellipses pseudo-section",
        ),
        (
            "Strike rose",
            "plot_strike_rose",
            "Geoelectric strike rose diagram",
        ),
        (
            "Strike profile",
            "plot_strike_profile",
            "Strike angle vs profile distance",
        ),
        (
            "Dimensionality",
            "plot_dimensionality",
            "1D / 2D / 3D classification section",
        ),
        (
            "Bostick depth table",
            "plot_bostick_depths",
            "Penetration depth per station/frequency",
        ),
        (
            "Gradient imaging",
            "plot_gradient_section",
            "Joint spatial-frequency gradient (Zhang 2021)",
        ),
        (
            "Phasor wheel",
            "plot_phasor_wheel",
            "Impedance hodogram in complex plane",
        ),
        (
            "Induction arrows",
            "plot_induction_arrows",
            "Induction arrows on profile",
        ),
        (
            "SNR section",
            "plot_snr_section",
            "Signal-to-noise ratio pseudo-section",
        ),
    ],
    "Uncertainty (MC)": [
        (
            "K uncertainty section",
            "plot_mc_K_section",
            "Hydraulic K P10 / P50 / P90 bounds",
        ),
        (
            "WT uncertainty",
            "plot_mc_wt_profile",
            "Water-table depth uncertainty bands",
        ),
        (
            "Parameter histograms",
            "plot_mc_histograms",
            "MC ensemble distributions (K, Sw, WT)",
        ),
    ],
    "Advanced Plots": [
        ("Mohr circles", "plot_mohr_circles", "Impedance Mohr diagram"),
        ("Argand diagram", "plot_argand", "Argand plane (Re Z vs Im Z)"),
        (
            "Phase tensor ternary",
            "plot_pt_ternary",
            "1D / 2D / 3D ternary plot",
        ),
        (
            "Distortion radar",
            "plot_distortion_radar",
            "Static shift · twist · shear",
        ),
        ("Rho-phase Bode", "plot_bode", "Log10(ρ_a) and phase vs period"),
        (
            "Z invariants section",
            "plot_z_invariants",
            "Det, skew, SSQ invariant sections",
        ),
        (
            "Anisotropy section",
            "plot_anisotropy",
            "Apparent anisotropy magnitude",
        ),
        (
            "Period clock",
            "plot_period_clock",
            "Phase tensor clock-face display",
        ),
        (
            "Composite section",
            "plot_composite_section",
            "Multi-metric combined pseudo-section",
        ),
    ],
    "Fusion & Time-Lapse": [
        ("Fused model", "plot_fused_model", "TDEM + AMT merged 2-D section"),
        (
            "Time-lapse change",
            "plot_timelapse_change",
            "ΔΔρ_a between two surveys",
        ),
        (
            "ΔSaturation section",
            "plot_timelapse_sat",
            "ΔSw change from resistivity change",
        ),
        (
            "WT displacement",
            "plot_timelapse_wt",
            "Water-table rise / fall per survey pair",
        ),
    ],
    "Export": [
        (
            "Export summary",
            "plot_export_preview",
            "Preview of exportable data tables",
        ),
    ],
}

CATEGORIES = list(WORKFLOW_CATALOGUE.keys())


# ── State container ────────────────────────────────────────────────────────────


@dataclass
class InterpState:
    sites: Any = None  # pycsamt Sites
    model: Any = None  # ResistivityModel
    boreholes: list = field(default_factory=list)
    db: Any = None  # RockDatabase
    petro_cfg: Any = None  # PetrophysicalConfig
    strat_logs: list = field(default_factory=list)
    hydro_result: Any = None  # EMHydroResult
    mc_result: Any = None  # UncertaintyResult
    constraints: list = field(default_factory=list)
    fusion_model: Any = None  # ResistivityModel (fused)
    timelapse_surveys: list = field(default_factory=list)
    structural_model: Any = None  # pycsamt.geology.structural.StructuralModel
    model_info: Any = None  # interp_sources.ModelInfo of the loaded model
    secondary_model: Any = None  # ResistivityModel to fuse with the model
    timelapse_labels: list = field(default_factory=list)
    original_model: Any = None  # the model before borehole calibration
    misfit_map: Any = None  # calibration misfit (ModelCalibrator)
    calibration_error: str = ""


# ── Controller ─────────────────────────────────────────────────────────────────


class InterpController:
    """
    Pure-Python interpretation controller.

    Holds :class:`InterpState` and exposes ``generate(plot_name, **kwargs)``
    which dispatches to the correct plot routine and returns a
    ``matplotlib.figure.Figure``.
    """

    def __init__(self) -> None:
        self.state = InterpState()
        self.dark: bool = True

    # ── Data setters ──────────────────────────────────────────────────────────

    def set_sites(self, sites) -> None:
        self.state.sites = sites

    def set_model(self, model) -> None:
        self.state.model = model

    def set_model_from_occam2d(self, result_dir: str) -> None:
        from pycsamt.interp._base import ResistivityModel
        # InversionResult lives in occam2d.results (it was imported from
        # occam2d.plot, which never defined it: "Load Occam2D…" always
        # failed with ImportError).
        from pycsamt.models.occam2d import InversionResult

        res = InversionResult(result_dir)
        self.state.model = ResistivityModel.from_occam2d(res)

    def load_model_source(self, path: str, line: str | None = None):
        """Any PCSF/PCSM file or inversion run folder -> the model
        (see :mod:`pycsamt.app.desktop.controllers.interp_sources`)."""
        from pycsamt.app.desktop.controllers.interp_sources import load_model

        model, info = load_model(path, line=line)
        self.state.model = model
        self.state.model_info = info
        # results computed on the previous model no longer apply
        self.state.strat_logs = []
        self.state.hydro_result = None
        self.state.mc_result = None
        self.state.fusion_model = None
        return info

    def add_timelapse_survey(self, path: str, label: str = "") -> int:
        """Append a repeat-survey model (same grid as the others)."""
        from pycsamt.app.desktop.controllers.interp_sources import load_model

        from pycsamt.interp.timelapse import assert_compatible_grids

        model, _info = load_model(path)
        surveys = list(self.state.timelapse_surveys)
        labels = list(self.state.timelapse_labels)
        if not surveys and self.state.model is not None:
            # the baseline is the model as inverted, not its calibration
            surveys.append(self.state.original_model or self.state.model)
            labels.append("baseline")
        surveys.append(model)
        labels.append(label or f"survey {len(surveys)}")
        assert_compatible_grids(surveys)  # raises before anything changes
        self.state.timelapse_surveys = surveys
        self.state.timelapse_labels = labels
        return len(surveys)

    def clear_timelapse(self) -> None:
        self.state.timelapse_surveys = []
        self.state.timelapse_labels = []

    def set_secondary_model(self, path: str) -> str:
        from pycsamt.app.desktop.controllers.interp_sources import load_model

        model, info = load_model(path)
        self.state.secondary_model = model
        return info.source

    def run_fusion(self, primary_max_depth: float = 0.0,
                   secondary_min_depth: float = 0.0,
                   blend: str = "linear") -> str:
        if self.state.model is None or self.state.secondary_model is None:
            return "Load a model and a second (deeper) model first."
        from pycsamt.interp.fusion import MultiMethodEMModel

        mm = MultiMethodEMModel(
            self.state.model, self.state.secondary_model,
            primary_max_depth=primary_max_depth or None,
            secondary_min_depth=secondary_min_depth or None, blend=blend)
        self.state.fusion_model = mm.merge()
        return "Models fused."

    def add_borehole_csv(self, path: str) -> str:
        from pycsamt.geology.borehole import Borehole

        bh = Borehole.from_csv(path)
        self._check_borehole(bh, path)
        self.state.boreholes.append(bh)
        return bh.name

    def add_borehole_las(self, path: str) -> str:
        from pycsamt.geology.borehole import Borehole

        bh = Borehole.from_las(path)
        self._check_borehole(bh, path)
        self.state.boreholes.append(bh)
        return bh.name

    @staticmethod
    def _check_borehole(bh, path) -> None:
        # a file that is not a log reads as a borehole with no intervals,
        # which then silently contributes nothing to calibration
        if not getattr(bh, "intervals", None):
            from pathlib import Path

            raise ValueError(f"no log intervals could be read from "
                             f"{Path(path).name}")

    def pcbh_positions(self, path: str) -> tuple[list[str], dict[str, float]]:
        """Borehole ids in a PCBH file and the profile distance of those
        whose id, name or alias is a model station."""
        from pycsamt.format.borehole import read_pcbh

        doc = read_pcbh(path)
        model = self.state.model
        known: dict[str, float] = {}
        if model is not None and getattr(model, "station_names", None):
            known = {str(n).lower(): float(x) for n, x in
                     zip(model.station_names, model.station_x)}
        ids, pos = [], {}
        for bh in doc.boreholes:
            ids.append(bh.id)
            for key in (bh.id, bh.name, *bh.aliases):
                if key and str(key).lower() in known:
                    pos[bh.id] = known[str(key).lower()]
                    break
        return ids, pos

    def add_borehole_pcbh(self, path: str,
                          profile_x: dict[str, float]) -> list[str]:
        """Add the boreholes of a PCBH file; *profile_x* gives each one's
        distance along the model profile (PCBH collars are map
        coordinates, never silently taken as profile distance)."""
        from pycsamt.format.borehole import read_pcbh
        from pycsamt.format.borehole.adapters import legacy_borehole_views

        doc = read_pcbh(path)
        if self.state.db is None:
            self.set_rock_db_default()
        views = legacy_borehole_views(doc, profile_x=dict(profile_x),
                                      rock_db=self.state.db)
        self.state.boreholes.extend(views)
        return [b.name for b in views]

    def remove_borehole(self, name: str) -> None:
        self.state.boreholes = [
            b for b in self.state.boreholes if b.name != name
        ]

    # ── Structural geology ───────────────────────────────────────────────────
    # Same three-evidence-type model (planar/linear/faults) Map View's Geology
    # rail edits and previews (pycsamt.app._structure); the desktop loads it
    # from CSV rather than a Studio table, matching how boreholes are loaded
    # here (CSV/LAS) instead of interactively.

    def _ensure_structural_model(self):
        from pycsamt.geology.structural import StructuralModel

        if self.state.structural_model is None:
            self.state.structural_model = StructuralModel()
        return self.state.structural_model

    def add_structural_planar_csv(self, path: str) -> int:
        from pycsamt.geology.structural import StructuralModel

        loaded = StructuralModel.from_csv(planar_path=path)
        model = self._ensure_structural_model()
        for m in loaded.planar:
            model.add_planar(m)
        return len(loaded.planar)

    def add_structural_linear_csv(self, path: str) -> int:
        from pycsamt.geology.structural import StructuralModel

        loaded = StructuralModel.from_csv(linear_path=path)
        model = self._ensure_structural_model()
        for m in loaded.linear:
            model.add_linear(m)
        return len(loaded.linear)

    def add_structural_faults_csv(self, path: str) -> int:
        from pycsamt.geology.structural import StructuralModel

        loaded = StructuralModel.from_csv(faults_path=path)
        model = self._ensure_structural_model()
        for f in loaded.faults:
            model.add_fault(f)
        return len(loaded.faults)

    def clear_structural_model(self) -> None:
        self.state.structural_model = None

    def set_rock_db_default(self) -> None:
        from pycsamt.geology.lithology import RockDatabase

        self.state.db = RockDatabase.default()

    def set_rock_db_csv(self, path: str) -> None:
        from pycsamt.geology.lithology import RockDatabase

        self.state.db = RockDatabase.from_csv(path)

    def set_petro_config(self, **kwargs) -> None:
        from pycsamt.interp.hydromodel import (
            PetrophysicalConfig,
        )
        from pycsamt.interp.petrophysics import ArchieModel

        from pycsamt.interp.petrophysics import WaxmanSmitsModel

        m = float(kwargs.get("m", 2.0))
        n = float(kwargs.get("n", 2.0))
        a = float(kwargs.get("a", 1.0))
        rho_w = float(kwargs.get("rho_w", 10.0))
        phi = float(kwargs.get("phi", 0.35))
        d50 = float(kwargs.get("d50_m", 5e-4))
        if kwargs.get("petro_model", "archie") == "waxman_smits":
            # clay surface conduction (shaly sands)
            petro = WaxmanSmitsModel(m=m, n=n, a=a,
                                     sigma_s=float(kwargs.get("sigma_s",
                                                              0.0)))
        else:
            petro = ArchieModel(m=m, n=n, a=a)
        self.state.petro_cfg = PetrophysicalConfig(
            petro=petro,
            rho_w=rho_w,
            porosity_prior=phi,
            d50_m=d50,
        )

    # ── Model status summary ───────────────────────────────────────────────────

    @property
    def model_status(self) -> dict:
        m = self.state.model
        if m is None:
            return {"loaded": False}
        return {
            "loaded": True,
            "n_x": getattr(m, "n_x", "?"),
            "n_z": getattr(m, "n_z", "?"),
            "depth_max": getattr(m, "depth_max", None),
            "profile_m": getattr(m, "profile_length", None),
            "method": getattr(m, "method", "—"),
            "rms": getattr(m, "rms", None),
        }

    @property
    def has_model(self) -> bool:
        return self.state.model is not None

    @property
    def has_sites(self) -> bool:
        return self.state.sites is not None

    # ── Computation ───────────────────────────────────────────────────────────

    def run_geological(self, station: str = "") -> str:
        """Classify model → strat logs + calibrated model if boreholes present."""
        if self.state.model is None:
            return "No model loaded."
        db = self.state.db
        if db is None:
            self.set_rock_db_default()
            db = self.state.db
        from pycsamt.geology.lithology import StratigraphicLog

        model = self.state.model
        logs = []
        try:
            # One log per station, at the model column under it.  (The loop
            # used to run over every grid column -- 576 for an Occam2D mesh
            # -- indexing the 47 station names: "list index out of range".)
            xs = np.asarray(model.x_centers, dtype=float)
            names = list(getattr(model, "station_names", None) or [])
            sx = np.asarray(getattr(model, "station_x", xs), dtype=float)
            if not names:
                names = [f"S{i + 1}" for i in range(sx.size)]
            for st, x in zip(names, sx):
                i = int(np.argmin(np.abs(xs - x)))
                log = StratigraphicLog.from_column(
                    st, float(xs[i]), model.z_centers, model.rho_2d[:, i], db)
                logs.append(log)
        except Exception as exc:
            return f"Geological classification failed: {exc}"
        self.state.strat_logs = logs
        self.state.calibration_error = ""
        if self.state.boreholes:
            try:
                from pycsamt.interp.calibrate import (
                    ModelCalibrator,
                )

                cal = ModelCalibrator(db=db)
                cal.fit(model, self.state.boreholes)
                self.state.original_model = model
                self.state.model = cal.calibrated_model()
                self.state.misfit_map = cal.misfit_map()
                self.state.strat_logs = cal.stratigraphic_logs()
            except Exception as exc:
                # was swallowed silently: say why calibration did not apply
                self.state.calibration_error = str(exc)
                return (f"Classified {len(logs)} stations; borehole "
                        f"calibration failed: {exc}")
            return (f"Classified {len(logs)} stations and calibrated "
                    f"against {len(self.state.boreholes)} borehole(s).")
        return f"Classified {len(logs)} stations."

    def run_hydro(self) -> str:
        """Run quantitative hydro estimation → EMHydroResult."""
        if self.state.model is None:
            return "No model loaded."
        if self.state.petro_cfg is None:
            self.set_petro_config()
        try:
            from pycsamt.interp.hydromodel import EMHydroModel

            hm = EMHydroModel(self.state.model, self.state.petro_cfg)
            self.state.hydro_result = hm.fit()
            return "Hydro estimation complete."
        except Exception as exc:
            return f"Hydro estimation failed: {exc}"

    def run_monte_carlo(
        self,
        n_samples: int = 200,
        rho_w_range=(5.0, 50.0),
        m_range=(1.5, 2.5),
        n_range=(1.5, 2.5),
        phi_range=(0.15, 0.45),
    ) -> str:
        if self.state.model is None or self.state.petro_cfg is None:
            return "Run hydro estimation first."
        try:
            from pycsamt.interp.uncertainty import (
                MonteCarloHydro,
                UncertaintyBounds,
            )

            bounds = UncertaintyBounds(
                rho_w_range=rho_w_range,
                m_range=m_range,
                n_range=n_range,
                phi_prior_range=phi_range,
            )
            mc = MonteCarloHydro(
                self.state.model,
                self.state.petro_cfg,
                bounds,
                n_samples=n_samples,
            )
            self.state.mc_result = mc.run()
            return f"MC complete ({n_samples} samples)."
        except Exception as exc:
            return f"MC failed: {exc}"

    # ── Plot dispatcher ────────────────────────────────────────────────────────

    def generate(self, method_name: str, **kwargs) -> matplotlib.figure.Figure:
        """Dispatch to the correct plot method. Returns a Figure."""
        fn = getattr(self, method_name, None)
        if fn is None:
            return self._not_implemented(method_name)
        try:
            return fn(**kwargs)
        except Exception as exc:
            return self._error_fig(f"{method_name} failed:\n{exc}")

    # ── Plot implementations ───────────────────────────────────────────────────

    # ── Setup & Model ──────────────────────────────────────────────────────

    def plot_model_summary(self, **kw) -> Figure:
        if self.state.model is None:
            return self._no_model_fig()
        import matplotlib.pyplot as plt

        model = self.state.model
        fig, ax = plt.subplots(figsize=(12, 5))
        self._apply_fig_style(fig, ax)
        try:
            rho = model.rho_2d
            x = (
                model.x_centers
                if hasattr(model, "x_centers")
                else np.arange(rho.shape[1])
            )
            z = (
                model.z_centers
                if hasattr(model, "z_centers")
                else np.arange(rho.shape[0])
            )
            cmap = "RdYlBu_r"
            im = ax.pcolormesh(x / 1e3, z, rho, cmap=cmap, shading="auto")
            ax.invert_yaxis()
            ax.set_xlabel("Profile distance (km)", fontsize=9)
            ax.set_ylabel("Depth (m)", fontsize=9)
            ax.set_title(
                f"Resistivity Model  ·  {getattr(model, 'method', '')}  "
                f"RMS {getattr(model, 'rms', '—')}",
                fontsize=10,
            )
            cb = fig.colorbar(im, ax=ax, pad=0.01, aspect=30)
            cb.set_label("log₁₀ ρ (Ω·m)", fontsize=8)
            self._style_cb(cb)
            # Station markers (flat datum)
            if hasattr(model, "station_x") and model.station_x is not None:
                ax.scatter(
                    model.station_x / 1e3,
                    np.zeros(len(model.station_x)) - 5,
                    marker="v",
                    color="#f38ba8",
                    s=40,
                    zorder=5,
                )
            self._maybe_draw_topo(ax, model)
        except Exception as exc:
            ax.text(
                0.5,
                0.5,
                f"Render error: {exc}",
                transform=ax.transAxes,
                ha="center",
                color="#f38ba8",
            )
        return fig

    def plot_depth_coverage(self, **kw) -> Figure:
        if self.state.sites is None:
            return self._no_sites_fig()
        import matplotlib.pyplot as plt

        import pycsamt.emtools as et

        fig, ax = plt.subplots(figsize=(12, 5))
        self._apply_fig_style(fig, ax)
        try:
            fn = getattr(et, "plot_depth_section", None)
            if fn:
                fn(self.state.sites, ax=ax, verbose=0)
            else:
                ax.text(
                    0.5,
                    0.5,
                    "plot_depth_section not available",
                    transform=ax.transAxes,
                    ha="center",
                )
        except Exception as exc:
            ax.text(0.5, 0.5, f"{exc}", transform=ax.transAxes, ha="center")
        return fig

    def plot_borehole_map(self, **kw) -> Figure:
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(10, 4))
        self._apply_fig_style(fig, ax)
        model = self.state.model
        bhs = self.state.boreholes
        if model is None and not bhs:
            ax.text(
                0.5,
                0.5,
                "No model or boreholes loaded",
                transform=ax.transAxes,
                ha="center",
            )
            return fig
        try:
            if model is not None and hasattr(model, "station_x"):
                ax.scatter(
                    model.station_x / 1e3,
                    np.zeros(len(model.station_x)),
                    marker="v",
                    color="#89b4fa",
                    s=60,
                    label="Stations",
                    zorder=3,
                )
            for bh in bhs:
                x = getattr(bh, "x", 0)
                ax.axvline(
                    x / 1e3, color="#a6e3a1", lw=1.5, ls="--", alpha=0.8
                )
                ax.text(
                    x / 1e3,
                    0.5,
                    bh.name,
                    fontsize=7,
                    ha="center",
                    color="#a6e3a1",
                    va="bottom",
                )
            ax.set_xlabel("Profile distance (km)", fontsize=9)
            ax.set_title("Station and Borehole Locations", fontsize=10)
            ax.legend(fontsize=7)
        except Exception as exc:
            ax.text(0.5, 0.5, str(exc), transform=ax.transAxes, ha="center")
        return fig

    # ── Geological ────────────────────────────────────────────────────────

    def plot_strat_log(self, station: str = "", **kw) -> Figure:
        if not self.state.strat_logs:
            return self._needs_run_fig("Run Geological Classification first")
        try:
            from pycsamt.interp.plot import (
                PlotStratigraphicLog,
            )

            log = next(
                (
                    l
                    for l in self.state.strat_logs
                    if l.station_name == station
                ),
                self.state.strat_logs[0],
            )
            plotter = PlotStratigraphicLog(log)
            fig = plotter.plot()
            self._apply_fig_style_minimal(fig)
            return fig
        except Exception as exc:
            return self._error_fig(str(exc))

    def plot_fence_diagram(self, **kw) -> Figure:
        if not self.state.strat_logs:
            return self._needs_run_fig("Run Geological Classification first")
        try:
            from pycsamt.interp.plot import PlotFenceDiagram

            plotter = PlotFenceDiagram(self.state.strat_logs)
            fig = plotter.plot()
            self._apply_fig_style_minimal(fig)
            return fig
        except Exception as exc:
            return self._error_fig(str(exc))

    # ── Structural ────────────────────────────────────────────────────────

    _STRUCTURAL_SENSE_COLOR = {
        "normal": "#3b82f6",
        "reverse": "#ef4444",
        "strike_slip": "#a855f7",
        "unknown": "#6b7280",
    }

    def plot_structural_section(
        self, *, depth_extent: float = 300.0, **kw
    ) -> Figure:
        """Profile-position/depth section of faults + planar/linear picks.

        Matplotlib port of the exact same reading
        :func:`pycsamt.app._structure.structure_section_figure` uses for
        Map View's Geology rail preview: a fault trace is a straight line
        through ``(x, z_top)`` tilted by its apparent dip, a planar
        measurement is a short tilted tick through ``(x, z)`` labelled
        ``strike/dip``, and a linear measurement is a diamond marker
        labelled ``trend/plunge``.
        """
        import matplotlib.pyplot as plt

        from pycsamt.geology.structural import StructuralModel

        model = self.state.structural_model or StructuralModel()
        fig, ax = plt.subplots(figsize=(9, 5))
        self._apply_fig_style(fig, ax)

        if not (model.faults or model.planar or model.linear):
            ax.text(
                0.5,
                0.5,
                "No structural data yet — load planar/linear/fault CSVs",
                transform=ax.transAxes,
                ha="center",
                va="center",
            )
            return fig

        half = max(float(depth_extent) * 0.06, 5.0)

        for fault in model.faults:
            z_top = float(fault.z_top or 0.0)
            z_bot = z_top + float(depth_extent)
            dip = max(1.0, min(89.0, float(fault.dip_deg)))
            direction = 1.0 if fault.downthrown_side == "right" else -1.0
            dx = direction * float(depth_extent) / np.tan(np.radians(dip))
            color = self._STRUCTURAL_SENSE_COLOR.get(
                fault.sense, self._STRUCTURAL_SENSE_COLOR["unknown"]
            )
            ax.plot(
                [fault.x, fault.x + dx],
                [z_top, z_bot],
                color=color,
                lw=2.5,
                solid_capstyle="round",
            )

        for m in model.planar:
            z = float(m.z or 0.0)
            dip = max(1.0, min(89.0, float(m.dip_deg)))
            dx = half / np.tan(np.radians(dip))
            ax.plot(
                [m.x - dx, m.x + dx],
                [z - half, z + half],
                color="#22c55e",
                lw=1.5,
                marker="o",
                ms=3,
            )
            ax.annotate(
                f"{m.strike_deg:.0f}/{m.dip_deg:.0f}",
                (m.x, z),
                fontsize=6,
                color="#22c55e",
            )

        for m in model.linear:
            z = float(m.z or 0.0)
            ax.plot(m.x, z, marker="D", ms=8, color="#f59e0b")
            ax.annotate(
                f"{m.trend_deg:.0f}/{m.plunge_deg:.0f}",
                (m.x, z),
                fontsize=6,
                color="#f59e0b",
            )

        ax.invert_yaxis()
        ax.set_xlabel("Profile position (m)", fontsize=9)
        ax.set_ylabel("Depth (m)", fontsize=9)
        ax.set_title("Structural section", fontsize=10)
        from matplotlib.lines import Line2D

        legend_items = [
            Line2D([0], [0], color=c, lw=2.5, label=sense)
            for sense, c in self._STRUCTURAL_SENSE_COLOR.items()
            if any(f.sense == sense for f in model.faults)
        ]
        if model.planar:
            legend_items.append(
                Line2D(
                    [0], [0], color="#22c55e", lw=1.5, marker="o",
                    ms=3, label="planar",
                )
            )
        if model.linear:
            legend_items.append(
                Line2D(
                    [0], [0], color="#f59e0b", marker="D", ms=6,
                    ls="none", label="linear",
                )
            )
        if legend_items:
            ax.legend(handles=legend_items, fontsize=7)
        return fig

    def plot_calibrated_model(self, **kw) -> Figure:
        if self.state.model is None:
            return self._no_model_fig()
        try:
            from pycsamt.interp.plot import (
                PlotCalibratedModel,
            )

            # (was called with the model only -> TypeError every time)
            if self.state.original_model is None:
                reason = self.state.calibration_error or (
                    "Add boreholes and run Classify & calibrate")
                return self._needs_run_fig(reason)
            plotter = PlotCalibratedModel(self.state.original_model,
                                          self.state.model,
                                          self.state.misfit_map)
            fig = plotter.plot()
            self._apply_fig_style_minimal(fig)
            return fig
        except Exception as exc:
            return self._error_fig(str(exc))

    def plot_rock_db(self, **kw) -> Figure:
        db = self.state.db
        if db is None:
            self.set_rock_db_default()
            db = self.state.db
        import matplotlib.pyplot as plt

        rocks = list(db.entries)
        fig, ax = plt.subplots(figsize=(10, max(4, len(rocks) * 0.4 + 1)))
        self._apply_fig_style(fig, ax)
        try:
            for i, rock in enumerate(rocks):
                color = getattr(rock, "color", "#888888")
                lo = getattr(rock, "rho_min", 0.1)
                hi = getattr(rock, "rho_max", 10000.0)
                ax.barh(
                    i,
                    np.log10(hi) - np.log10(lo),
                    left=np.log10(lo),
                    color=color,
                    alpha=0.85,
                    edgecolor="#313244",
                    linewidth=0.5,
                )
                ax.text(
                    np.log10(lo) - 0.05,
                    i,
                    rock.name,
                    ha="right",
                    va="center",
                    fontsize=7,
                )
            ax.set_xlabel("log₁₀ ρ (Ω·m)", fontsize=9)
            ax.set_title("Rock Database — Resistivity Ranges", fontsize=10)
            ax.set_yticks([])
        except Exception as exc:
            ax.text(0.5, 0.5, str(exc), transform=ax.transAxes, ha="center")
        return fig

    # ── Hydrology ─────────────────────────────────────────────────────────

    def _hydro_section_plot(
        self, attr: str, title: str, cmap: str = "Blues_r", label: str = ""
    ) -> Figure:
        if self.state.hydro_result is None:
            return self._needs_run_fig("Run Hydrology Estimation first")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(12, 5))
        self._apply_fig_style(fig, ax)
        try:
            model = self.state.model
            data = getattr(self.state.hydro_result, attr, None)
            if data is None:
                ax.text(
                    0.5,
                    0.5,
                    f"'{attr}' not in result",
                    transform=ax.transAxes,
                    ha="center",
                )
                return fig
            x = (
                model.x_centers / 1e3
                if hasattr(model, "x_centers")
                else np.arange(data.shape[1])
            )
            z = (
                model.z_centers
                if hasattr(model, "z_centers")
                else np.arange(data.shape[0])
            )
            im = ax.pcolormesh(
                x,
                z,
                np.log10(np.maximum(data, 1e-20)),
                cmap=cmap,
                shading="auto",
            )
            ax.invert_yaxis()
            ax.set_xlabel("Profile distance (km)", fontsize=9)
            ax.set_ylabel("Depth (m)", fontsize=9)
            ax.set_title(title, fontsize=10)
            cb = fig.colorbar(im, ax=ax, pad=0.01, aspect=30)
            cb.set_label(label, fontsize=8)
            self._style_cb(cb)
            self._maybe_draw_topo(ax, model)
        except Exception as exc:
            ax.text(0.5, 0.5, str(exc), transform=ax.transAxes, ha="center")
        return fig

    def plot_K_map(self, **kw):
        return self._hydro_section_plot(
            "hydraulic_K",
            "Hydraulic Conductivity K  (m/s)",
            "viridis_r",
            "log₁₀ K (m/s)",
        )

    def plot_Sw_map(self, **kw):
        return self._hydro_section_plot(
            "saturation", "Water Saturation Sw", "Blues", "Sw"
        )

    def plot_water_table(self, **kw) -> Figure:
        if self.state.hydro_result is None:
            return self._needs_run_fig("Run Hydrology Estimation first")
        try:
            from pycsamt.interp.plot import (
                PlotWaterTableProfile,
            )

            fig = PlotWaterTableProfile(self.state.hydro_result).plot()
            self._apply_fig_style_minimal(fig)
            return fig
        except Exception as exc:
            return self._error_fig(str(exc))

    def plot_aquifer_zones(self, **kw) -> Figure:
        if self.state.model is None:
            return self._no_model_fig()
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(12, 5))
        self._apply_fig_style(fig, ax)
        try:
            from pycsamt.interp.hydro import HydroInterpreter

            db = self.state.db
            if db is None:
                self.set_rock_db_default()
                db = self.state.db
            # (was HydroInterpreter(model, db=) + identify_aquifers(): the
            # constructor takes no model and that method does not exist)
            hi = HydroInterpreter(db=db)
            hi.fit(self.state.original_model or self.state.model,
                   boreholes=self.state.boreholes or None)
            zones = hi.aquifer_zones()
            model = self.state.model
            x = (
                model.x_centers / 1e3
                if hasattr(model, "x_centers")
                else np.arange(model.n_x)
            )
            z = (
                model.z_centers
                if hasattr(model, "z_centers")
                else np.arange(model.n_z)
            )
            rho = model.rho_2d
            ax.pcolormesh(x, z, rho, cmap="Greys_r", shading="auto", alpha=0.4)
            for zone in zones:
                xi = float(getattr(zone, "x", 0)) / 1e3
                top, bot = float(zone.top), float(zone.bottom)
                conf = float(getattr(zone, "confidence", 0.8))
                rect = plt.Rectangle(
                    (xi - 0.05, top),
                    0.1,
                    bot - top,
                    fc="#89b4fa",
                    alpha=conf * 0.7,
                    ec="#89b4fa",
                    lw=1,
                )
                ax.add_patch(rect)
            ax.invert_yaxis()
            ax.set_xlabel("Profile distance (km)", fontsize=9)
            ax.set_ylabel("Depth (m)", fontsize=9)
            ax.set_title("Identified Aquifer Zones", fontsize=10)
        except Exception as exc:
            ax.text(0.5, 0.5, str(exc), transform=ax.transAxes, ha="center")
        return fig

    def plot_transmissivity(self, **kw) -> Figure:
        if self.state.hydro_result is None:
            return self._needs_run_fig("Run Hydrology Estimation first")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(10, 4))
        self._apply_fig_style(fig, ax)
        try:
            T = self.state.hydro_result.transmissivity
            model = self.state.model
            x = (
                model.x_centers / 1e3
                if hasattr(model, "x_centers")
                else np.arange(len(T))
            )
            ax.semilogy(x, T, color="#89b4fa", lw=1.8, marker="o", ms=4)
            ax.set_xlabel("Profile distance (km)", fontsize=9)
            ax.set_ylabel("Transmissivity (m²/s)", fontsize=9)
            ax.set_title("Transmissivity Profile", fontsize=10)
            ax.grid(True, which="both", alpha=0.2, ls="--")
        except Exception as exc:
            ax.text(0.5, 0.5, str(exc), transform=ax.transAxes, ha="center")
        return fig

    def plot_aquifer_char(self, **kw) -> Figure:
        if self.state.hydro_result is None:
            return self._needs_run_fig("Run Hydrology Estimation first")
        try:
            from pycsamt.interp.plot import (
                PlotAquiferCharacterization,
            )

            fig = PlotAquiferCharacterization(self.state.hydro_result).plot()
            self._apply_fig_style_minimal(fig)
            return fig
        except Exception as exc:
            return self._error_fig(str(exc))

    def plot_petrophys_xplot(self, **kw) -> Figure:
        if self.state.hydro_result is None:
            return self._needs_run_fig("Run Hydrology Estimation first")
        try:
            from pycsamt.interp.plot import (
                PlotPetrophysicalCrossPlot,
            )

            # (was called with 3 positionals -> TypeError every time)
            cfg = self.state.petro_cfg
            fig = PlotPetrophysicalCrossPlot(
                self.state.hydro_result,
                petro=getattr(cfg, "petro", None),
                show_hs_bounds=bool(kw.get("show_hs_bounds", True)),
                rho_matrix=float(kw.get("rho_matrix", 5000.0)),
            ).plot()
            self._apply_fig_style_minimal(fig)
            return fig
        except Exception as exc:
            return self._error_fig(str(exc))

    # ── Field Constraints ─────────────────────────────────────────────────

    def plot_constraint_misfit(self, **kw) -> Figure:
        return self._needs_run_fig("Run Constraint Calibration first")

    def plot_calib_history(self, **kw) -> Figure:
        return self._needs_run_fig("Run Constraint Calibration first")

    # ── EM Diagnostics ────────────────────────────────────────────────────

    def _emtools_plot(
        self, fn_name: str, title: str, polar: bool = False, **kw
    ) -> Figure:
        if self.state.sites is None:
            return self._no_sites_fig()
        import inspect

        import matplotlib.pyplot as plt

        import pycsamt.emtools as et

        fn = getattr(et, fn_name, None)
        if fn is None:
            return self._not_implemented(fn_name)
        params = inspect.signature(fn).parameters
        extra = {"verbose": 0} if "verbose" in params else {}
        result = None
        if "ax" in params and not polar:
            # single-panel function: draw into our axes
            fig, ax = plt.subplots(figsize=(12, 5))
            self._apply_fig_style(fig, ax)
            try:
                result = fn(self.state.sites, ax=ax, **extra, **kw)
            except Exception:
                plt.close(fig)
                result = None
            else:
                return fig if not isinstance(result, plt.Figure) else result
        # Multi-panel ("axes=") or polar functions build their own figure:
        # a single placeholder axes made them fail ("axes must provide at
        # least 2 axes", "no attribute set_theta_offset").
        try:
            result = fn(self.state.sites, **extra, **kw)
        except Exception as exc:
            return self._error_fig(f"{fn_name}: {exc}")
        fig = _figure_from(result)
        if fig is None:
            return self._error_fig(f"{fn_name} returned no figure")
        self._apply_fig_style_minimal(fig)
        return fig

    def plot_pt_section(self, **kw):
        return self._emtools_plot(
            "plot_phase_tensor_psection", "Phase Tensor Pseudo-Section"
        )

    def plot_strike_rose(self, **kw):
        return self._emtools_plot(
            "plot_strike_rose", "Strike Rose Diagram", polar=True
        )

    def plot_strike_profile(self, **kw):
        return self._emtools_plot("plot_strike_profile", "Strike Profile")

    def plot_dimensionality(self, **kw):
        return self._emtools_plot(
            "plot_dimensionality_psection", "Dimensionality Classification"
        )

    def plot_bostick_depths(self, **kw):
        return self._emtools_plot(
            "plot_depth_section", "Bostick Depth Coverage"
        )

    def plot_gradient_section(self, **kw):
        return self._emtools_plot(
            "plot_gradient_section", "Gradient Imaging (Zhang 2021)"
        )

    def plot_phasor_wheel(self, **kw):
        return self._emtools_plot(
            "plot_phasor_wheel", "Impedance Phasor Wheel", polar=True
        )

    def plot_induction_arrows(self, **kw):
        return self._emtools_plot("plot_induction_section", "Induction Arrows")

    def plot_snr_section(self, **kw):
        return self._emtools_plot("plot_snr_section", "SNR Pseudo-Section")

    # ── Uncertainty ───────────────────────────────────────────────────────

    def plot_mc_K_section(self, **kw) -> Figure:
        if self.state.mc_result is None:
            return self._needs_run_fig("Run Monte Carlo first")
        try:
            from pycsamt.interp.plot import (
                PlotUncertaintySection,
            )

            fig = PlotUncertaintySection(
                self.state.mc_result, quantity=kw.get("quantity", "K")
            ).plot()
            self._apply_fig_style_minimal(fig)
            return fig
        except Exception as exc:
            return self._error_fig(str(exc))

    def plot_mc_wt_profile(self, **kw) -> Figure:
        if self.state.mc_result is None:
            return self._needs_run_fig("Run Monte Carlo first")
        try:
            from pycsamt.interp.plot import (
                PlotUncertaintyProfile,
            )

            fig = PlotUncertaintyProfile(self.state.mc_result).plot()
            self._apply_fig_style_minimal(fig)
            return fig
        except Exception as exc:
            return self._error_fig(str(exc))

    def plot_mc_histograms(self, **kw) -> Figure:
        if self.state.mc_result is None:
            return self._needs_run_fig("Run Monte Carlo first")
        try:
            from pycsamt.interp.plot import (
                PlotUncertaintyHistogram,
            )

            fig = PlotUncertaintyHistogram(self.state.mc_result).plot()
            self._apply_fig_style_minimal(fig)
            return fig
        except Exception as exc:
            return self._error_fig(str(exc))

    # ── Advanced Plots ────────────────────────────────────────────────────

    def plot_mohr_circles(self, **kw):
        return self._emtools_plot(
            "plot_impedance_mohr_circles", "Impedance Mohr Circles"
        )

    def plot_argand(self, **kw):
        return self._emtools_plot("plot_zt_argand", "Impedance Argand Diagram")

    def plot_pt_ternary(self, **kw):
        return self._emtools_plot(
            "plot_dimensionality_ternary",
            "Dimensionality Ternary (1D / 2D / 3D)",
        )

    def plot_distortion_radar(self, **kw):
        return self._emtools_plot(
            "plot_distortion_radar",
            "Distortion Radar: Static Shift · Twist · Shear",
        )

    def plot_bode(self, **kw):
        return self._emtools_plot(
            "plot_rho_phase_bode", "Bode Plot: log₁₀(ρ_a) and Phase"
        )

    def plot_z_invariants(self, **kw):
        return self._emtools_plot(
            "plot_z_invariants_section", "Z Tensor Invariants Section"
        )

    def plot_anisotropy(self, **kw):
        return self._emtools_plot(
            "plot_apparent_anisotropy_section", "Apparent Anisotropy Section"
        )

    def plot_period_clock(self, **kw):
        return self._emtools_plot(
            "plot_pt_period_clock", "Phase Tensor Period Clock"
        )

    def plot_composite_section(self, **kw):
        return self._emtools_plot(
            "plot_mt_composite_section", "MT Composite Section (multi-metric)"
        )

    # ── Fusion & Time-Lapse ───────────────────────────────────────────────

    def plot_fused_model(self, **kw) -> Figure:
        if self.state.fusion_model is None:
            return self._needs_run_fig("Run Fusion first")
        self_copy = type(
            "_",
            (),
            {"state": type("_", (), {"model": self.state.fusion_model})()},
        )()
        self_copy.dark = self.dark
        self_copy._apply_fig_style = self._apply_fig_style
        self_copy._apply_fig_style_minimal = self._apply_fig_style_minimal
        self_copy._no_model_fig = self._no_model_fig
        self_copy._error_fig = self._error_fig
        self_copy._style_cb = self._style_cb
        return self.plot_model_summary.__func__(self_copy)

    def _timelapse(self):
        from pycsamt.interp.timelapse import TimeLapseEM

        labels = self.state.timelapse_labels or None
        return TimeLapseEM(surveys=self.state.timelapse_surveys,
                           labels=labels)

    def _timelapse_grid(self, quantity: str, **kw) -> Figure:
        if len(self.state.timelapse_surveys) < 2:
            return self._needs_run_fig(
                "Add at least one repeat survey (Monitoring ▸ Surveys)")
        from pycsamt.interp.plot import PlotMultiTimeLapseGrid

        cfg = self.state.petro_cfg
        if quantity == "delta_saturation" and cfg is None:
            self.set_petro_config()
            cfg = self.state.petro_cfg
        fig = PlotMultiTimeLapseGrid(
            self._timelapse(), quantity=quantity,
            petro=getattr(cfg, "petro", None),
            rho_w=float(getattr(cfg, "rho_w", 20.0)),
            phi=float(getattr(cfg, "porosity_prior", 0.25)),
        ).plot()
        self._apply_fig_style_minimal(fig)
        return fig

    def plot_timelapse_change(self, **kw) -> Figure:
        return self._timelapse_grid("delta_rho", **kw)

    def plot_timelapse_sat(self, **kw) -> Figure:
        return self._timelapse_grid("delta_saturation", **kw)

    def plot_timelapse_wt(self, **kw) -> Figure:
        if len(self.state.timelapse_surveys) < 2:
            return self._needs_run_fig(
                "Add at least one repeat survey (Monitoring ▸ Surveys)")
        import matplotlib.pyplot as plt

        if self.state.petro_cfg is None:
            self.set_petro_config()
        cfg = self.state.petro_cfg
        disp = np.asarray(self._timelapse().water_table_displacement(
            cfg.petro, rho_w=float(cfg.rho_w)), dtype=float)
        x = self.state.timelapse_surveys[0].x_centers
        fig, ax = plt.subplots(figsize=(9, 4))
        rows = disp if disp.ndim == 2 else disp[None, :]
        labels = self.state.timelapse_labels[1:] or [
            f"survey {i + 2}" for i in range(len(rows))]
        for row, lab in zip(rows, labels):
            ax.plot(x[: row.size], row, lw=1.6, label=lab)
        ax.axhline(0, color="0.5", lw=0.8)
        ax.set_xlabel("Distance (m)")
        ax.set_ylabel("Water-table change (m, + = rise)")
        ax.set_title("Water-table displacement from the baseline",
                     fontsize=10)
        ax.legend(fontsize=8)
        self._apply_fig_style(fig, ax)
        return fig

    def plot_depth_profile(self, station: str = "", **kw) -> Figure:
        if self.state.model is None:
            return self._no_model_fig()
        from pycsamt.interp.plot import PlotResistivityDepthProfile

        names = list(getattr(self.state.model, "station_names", []) or [])
        st = station if station in names else 0
        bh = next((b for b in self.state.boreholes
                   if getattr(b, "name", None) == station), None)
        fig = PlotResistivityDepthProfile(self.state.model, st,
                                          borehole=bh).plot()
        self._apply_fig_style_minimal(fig)
        return fig

    def plot_borehole_fence(self, **kw) -> Figure:
        if not self.state.boreholes:
            return self._needs_run_fig("Add boreholes first (Evidence)")
        from pycsamt.interp.plot import PlotBoreholeFence

        if self.state.db is None:
            self.set_rock_db_default()
        fig = PlotBoreholeFence(self.state.boreholes, db=self.state.db).plot()
        self._apply_fig_style_minimal(fig)
        return fig

    def plot_export_preview(self, **kw) -> Figure:
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(8, 5))
        self._apply_fig_style(fig, ax)
        lines = [
            "Available exports:",
            "  • Oasis Montaj XYZ (interp.export.to_oasis_montaj_xyz)",
            "  • LAS 2.0 well logs  (interp.export.to_las)",
            "  • Flat CSV table     (interp.export.to_csv)",
            "  • VTK RectilinearGrid (interp.export.to_vtk)",
            "",
            "Use the Export buttons in the panel.",
        ]
        s = _DARK if self.dark else _LIGHT
        for i, ln in enumerate(lines):
            ax.text(
                0.05,
                0.9 - i * 0.12,
                ln,
                transform=ax.transAxes,
                fontsize=9,
                color=s["fg"],
                family="monospace",
            )
        ax.set_axis_off()
        return fig

    # ── Export helpers ─────────────────────────────────────────────────────

    def export_xyz(self, path: str) -> str:
        if not self.state.strat_logs:
            return "Run geological classification first."
        from pycsamt.interp.export import to_oasis_montaj_xyz

        to_oasis_montaj_xyz(self.state.strat_logs, path)
        return f"Exported → {path}"

    def export_las(self, path: str, station: str = "") -> str:
        if not self.state.strat_logs:
            return "Run geological classification first."
        from pycsamt.interp.export import to_las

        log = next(
            (l for l in self.state.strat_logs if l.station_name == station),
            self.state.strat_logs[0],
        )
        to_las(log, path)
        return f"Exported → {path}"

    def export_csv(self, path: str) -> str:
        if not self.state.strat_logs:
            return "Run geological classification first."
        from pycsamt.interp.export import to_csv

        to_csv(self.state.strat_logs, path)
        return f"Exported → {path}"

    def export_vtk(self, path: str) -> str:
        if self.state.model is None:
            return "No model loaded."
        from pycsamt.interp.export import to_vtk

        to_vtk(self.state.model, path)
        return f"Exported → {path}"

    # ── Internal helpers ───────────────────────────────────────────────────

    def _apply_fig_style(self, fig, ax) -> None:
        s = _DARK if self.dark else _LIGHT
        fig.patch.set_facecolor(s["fig_bg"])
        ax.set_facecolor(s["bg"])
        ax.tick_params(colors=s["tick"], labelsize=7)
        ax.xaxis.label.set_color(s["fg"])
        ax.yaxis.label.set_color(s["fg"])
        ax.title.set_color(s["title"])
        for sp in ax.spines.values():
            sp.set_edgecolor(s["spine"])
        ax.grid(True, color=s["grid"], alpha=0.25, ls="--", lw=0.5)

    def _apply_fig_style_minimal(self, fig) -> None:
        s = _DARK if self.dark else _LIGHT
        fig.patch.set_facecolor(s["fig_bg"])
        for ax in fig.get_axes():
            try:
                ax.set_facecolor(s["bg"])
                ax.tick_params(colors=s["tick"], labelsize=7)
                ax.title.set_color(s["title"])
            except Exception:
                pass

    def _style_cb(self, cb) -> None:
        s = _DARK if self.dark else _LIGHT
        cb.ax.yaxis.set_tick_params(color=s["tick"], labelcolor=s["tick"])
        cb.set_label(cb.ax.get_ylabel(), fontsize=8, color=s["fg"])

    def _maybe_draw_topo(self, ax, model) -> None:
        """Overlay terrain on a depth-section axes when PYCSAMT_TOPO is enabled."""
        try:
            from pycsamt.topo.config import PYCSAMT_TOPO

            if not PYCSAMT_TOPO.enabled:
                return
            from pycsamt.topo.drape import interp_elev
            from pycsamt.topo.extract import (
                extract_chainage,
                extract_elevation,
            )
            from pycsamt.topo.overlay import draw_topo_section

            # Try to get station positions and elevations from the model
            sx_m = np.asarray(getattr(model, "station_x", []))
            if len(sx_m) == 0:
                return
            sx_km = sx_m / 1e3
            # Elevation: prefer model attribute, else try from sites state
            elev_m = getattr(model, "station_elev", None)
            if elev_m is None and self.state.sites is not None:
                elev_m = extract_elevation(self.state.sites, warn=False)
                chain = extract_chainage(self.state.sites)
                if len(chain) > 0 and len(sx_km) > 0:
                    elev_m = (
                        interp_elev(chain, elev_m / 1000.0, sx_km) * 1000.0
                    )
            if elev_m is None or len(elev_m) == 0:
                return
            elev_m = np.asarray(elev_m)
            draw_topo_section(
                ax,
                sx_km,
                elev_m,
                dark=self.dark,
            )
        except Exception:
            pass  # topo is always optional; never break the main figure

    def _no_model_fig(self) -> Figure:
        return self._msg_fig(
            "No resistivity model loaded.\n"
            "Load from inversion results or from file."
        )

    def _no_sites_fig(self) -> Figure:
        return self._msg_fig(
            "No survey data loaded.\n"
            "Open EDI files from the main window first."
        )

    def _needs_run_fig(self, msg: str) -> Figure:
        return self._msg_fig(f"⚠  {msg}")

    def _not_implemented(self, name: str) -> Figure:
        return self._msg_fig(f"Plot '{name}' is not yet implemented.")

    def _error_fig(self, msg: str) -> Figure:
        return self._msg_fig(f"✕  {msg}", color="#f38ba8")

    def _msg_fig(self, msg: str, color: str | None = None) -> Figure:
        import matplotlib.pyplot as plt

        s = _DARK if self.dark else _LIGHT
        fig, ax = plt.subplots(figsize=(8, 4))
        fig.patch.set_facecolor(s["fig_bg"])
        ax.set_facecolor(s["bg"])
        ax.set_axis_off()
        ax.text(
            0.5,
            0.5,
            msg,
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=11,
            color=color or s["muted"],
            wrap=True,
            multialignment="center",
        )
        return fig


# ── Theme dicts ────────────────────────────────────────────────────────────────

def _figure_from(obj):
    """The Figure behind a plot function's return value."""
    from matplotlib.figure import Figure

    if isinstance(obj, Figure):
        return obj
    if hasattr(obj, "figure") and isinstance(getattr(obj, "figure"), Figure):
        return obj.figure
    try:
        items = list(np.ravel(obj)) if not isinstance(obj, dict) else \
            list(obj.values())
    except Exception:
        items = []
    for it in items:
        f = _figure_from(it) if it is not obj else None
        if f is not None:
            return f
    import matplotlib.pyplot as plt

    return plt.gcf() if plt.get_fignums() else None


_DARK = dict(
    bg="#1e1e2e",
    fig_bg="#181825",
    fg="#cdd6f4",
    title="#cdd6f4",
    tick="#a6adc8",
    spine="#45475a",
    grid="#313244",
    muted="#585b70",
)
# Publication white (the theme greys #eff1f5/#e6e9ef made every exported
# interpretation figure grey, like the pipeline plots).
_LIGHT = dict(
    bg="#ffffff",
    fig_bg="#ffffff",
    fg="#1f2937",
    title="#111827",
    tick="#374151",
    spine="#6b7280",
    grid="#d1d5db",
    muted="#6b7280",
)
