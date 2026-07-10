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
        ("Model summary",        "plot_model_summary",       "2-D rho section overview"),
        ("Borehole locations",   "plot_borehole_map",        "Station + borehole location map"),
        ("Depth coverage",       "plot_depth_coverage",      "Bostick depth reach per station"),
    ],
    "Geological": [
        ("Stratigraphic log",    "plot_strat_log",           "Lithology column for selected station"),
        ("Fence diagram",        "plot_fence_diagram",       "Multi-station pseudo-stratigraphic section"),
        ("Calibrated model",     "plot_calibrated_model",    "NM vs original model + misfit map"),
        ("Rock DB summary",      "plot_rock_db",             "Colour-coded resistivity ranges from database"),
    ],
    "Hydrology": [
        ("Hydraulic K map",      "plot_K_map",               "2-D hydraulic conductivity section"),
        ("Water saturation",     "plot_Sw_map",              "2-D water saturation section"),
        ("Water table profile",  "plot_water_table",         "Water-table depth along profile ± confidence"),
        ("Aquifer zones",        "plot_aquifer_zones",       "Identified aquifer intervals with confidence"),
        ("Transmissivity",       "plot_transmissivity",      "Transmissivity (m²/s) profile"),
        ("Aquifer characterisation","plot_aquifer_char",     "Multi-panel: K, T, Dar-Zarrouk"),
        ("Petrophysical cross-plot","plot_petrophys_xplot",  "Archie-space: porosity vs Sw vs K"),
    ],
    "Field Constraints": [
        ("Constraint misfit",    "plot_constraint_misfit",   "Observed vs modelled field data"),
        ("Calibration history",  "plot_calib_history",       "Optimizer convergence trace"),
    ],
    "EM Diagnostics": [
        ("Phase tensor section", "plot_pt_section",          "Phase tensor ellipses pseudo-section"),
        ("Strike rose",          "plot_strike_rose",         "Geoelectric strike rose diagram"),
        ("Strike profile",       "plot_strike_profile",      "Strike angle vs profile distance"),
        ("Dimensionality",       "plot_dimensionality",      "1D / 2D / 3D classification section"),
        ("Bostick depth table",  "plot_bostick_depths",      "Penetration depth per station/frequency"),
        ("Gradient imaging",     "plot_gradient_section",    "Joint spatial-frequency gradient (Zhang 2021)"),
        ("Phasor wheel",         "plot_phasor_wheel",        "Impedance hodogram in complex plane"),
        ("Induction arrows",     "plot_induction_arrows",    "Induction arrows on profile"),
        ("SNR section",          "plot_snr_section",         "Signal-to-noise ratio pseudo-section"),
    ],
    "Uncertainty (MC)": [
        ("K uncertainty section","plot_mc_K_section",        "Hydraulic K P10 / P50 / P90 bounds"),
        ("WT uncertainty",       "plot_mc_wt_profile",       "Water-table depth uncertainty bands"),
        ("Parameter histograms", "plot_mc_histograms",       "MC ensemble distributions (K, Sw, WT)"),
    ],
    "Advanced Plots": [
        ("Mohr circles",         "plot_mohr_circles",        "Impedance Mohr diagram"),
        ("Argand diagram",       "plot_argand",              "Argand plane (Re Z vs Im Z)"),
        ("Phase tensor ternary", "plot_pt_ternary",          "1D / 2D / 3D ternary plot"),
        ("Distortion radar",     "plot_distortion_radar",    "Static shift · twist · shear"),
        ("Rho-phase Bode",       "plot_bode",                "Log10(ρ_a) and phase vs period"),
        ("Z invariants section", "plot_z_invariants",        "Det, skew, SSQ invariant sections"),
        ("Anisotropy section",   "plot_anisotropy",          "Apparent anisotropy magnitude"),
        ("Period clock",         "plot_period_clock",        "Phase tensor clock-face display"),
        ("Composite section",    "plot_composite_section",   "Multi-metric combined pseudo-section"),
    ],
    "Fusion & Time-Lapse": [
        ("Fused model",          "plot_fused_model",         "TDEM + AMT merged 2-D section"),
        ("Time-lapse change",    "plot_timelapse_change",    "ΔΔρ_a between two surveys"),
        ("ΔSaturation section",  "plot_timelapse_sat",       "ΔSw change from resistivity change"),
        ("WT displacement",      "plot_timelapse_wt",        "Water-table rise / fall per survey pair"),
    ],
    "Export": [
        ("Export summary",       "plot_export_preview",      "Preview of exportable data tables"),
    ],
}

CATEGORIES = list(WORKFLOW_CATALOGUE.keys())


# ── State container ────────────────────────────────────────────────────────────

@dataclass
class InterpState:
    sites:          Any  = None   # pycsamt Sites
    model:          Any  = None   # ResistivityModel
    boreholes:      list = field(default_factory=list)
    db:             Any  = None   # RockDatabase
    petro_cfg:      Any  = None   # PetrophysicalConfig
    strat_logs:     list = field(default_factory=list)
    hydro_result:   Any  = None   # EMHydroResult
    mc_result:      Any  = None   # UncertaintyResult
    constraints:    list = field(default_factory=list)
    fusion_model:   Any  = None   # ResistivityModel (fused)
    timelapse_surveys: list = field(default_factory=list)


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
        from pycsamt.models.occam2d.plot import (
            InversionResult,
        )
        res = InversionResult(result_dir)
        self.state.model = ResistivityModel.from_occam2d(res)

    def add_borehole_csv(self, path: str) -> str:
        from pycsamt.interp.borehole import Borehole
        bh = Borehole.from_csv(path)
        self.state.boreholes.append(bh)
        return bh.name

    def add_borehole_las(self, path: str) -> str:
        from pycsamt.interp.borehole import Borehole
        bh = Borehole.from_las(path)
        self.state.boreholes.append(bh)
        return bh.name

    def remove_borehole(self, name: str) -> None:
        self.state.boreholes = [b for b in self.state.boreholes if b.name != name]

    def set_rock_db_default(self) -> None:
        from pycsamt.interp.lithology import RockDatabase
        self.state.db = RockDatabase.default()

    def set_rock_db_csv(self, path: str) -> None:
        from pycsamt.interp.lithology import RockDatabase
        self.state.db = RockDatabase.from_csv(path)

    def set_petro_config(self, **kwargs) -> None:
        from pycsamt.interp.hydromodel import (
            PetrophysicalConfig,
        )
        from pycsamt.interp.petrophysics import ArchieModel
        m   = float(kwargs.get("m", 2.0))
        n   = float(kwargs.get("n", 2.0))
        a   = float(kwargs.get("a", 1.0))
        rho_w = float(kwargs.get("rho_w", 10.0))
        phi   = float(kwargs.get("phi", 0.35))
        d50   = float(kwargs.get("d50_m", 5e-4))
        self.state.petro_cfg = PetrophysicalConfig(
            petro=ArchieModel(m=m, n=n, a=a),
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
            "loaded":    True,
            "n_x":       getattr(m, "n_x", "?"),
            "n_z":       getattr(m, "n_z", "?"),
            "depth_max": getattr(m, "depth_max", None),
            "profile_m": getattr(m, "profile_length", None),
            "method":    getattr(m, "method", "—"),
            "rms":       getattr(m, "rms", None),
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
        from pycsamt.interp.lithology import StratigraphicLog
        model = self.state.model
        logs = []
        try:
            for i in range(model.n_x):
                col = model.x_centers[i] if hasattr(model, "x_centers") else i
                st  = (model.station_names[i]
                       if hasattr(model, "station_names") and model.station_names
                       else f"S{i+1}")
                z   = model.z_centers
                rho = model.rho_2d[:, i]
                log = StratigraphicLog.from_column(st, float(col), z, rho, db)
                logs.append(log)
        except Exception as exc:
            return f"Geological classification failed: {exc}"
        self.state.strat_logs = logs
        if self.state.boreholes:
            try:
                from pycsamt.interp.calibrate import (
                    ModelCalibrator,
                )
                cal = ModelCalibrator(db=db)
                cal.fit(model, self.state.boreholes)
                self.state.model = cal.calibrated_model()
                self.state.strat_logs = cal.stratigraphic_logs()
            except Exception:
                pass
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

    def run_monte_carlo(self, n_samples: int = 200,
                        rho_w_range=(5.0, 50.0),
                        m_range=(1.5, 2.5),
                        n_range=(1.5, 2.5),
                        phi_range=(0.15, 0.45)) -> str:
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
                self.state.model, self.state.petro_cfg,
                bounds, n_samples=n_samples,
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
            x   = model.x_centers if hasattr(model, "x_centers") else np.arange(rho.shape[1])
            z   = model.z_centers  if hasattr(model, "z_centers") else np.arange(rho.shape[0])
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
                ax.scatter(model.station_x / 1e3, np.zeros(len(model.station_x)) - 5,
                           marker="v", color="#f38ba8", s=40, zorder=5)
            self._maybe_draw_topo(ax, model)
        except Exception as exc:
            ax.text(0.5, 0.5, f"Render error: {exc}", transform=ax.transAxes,
                    ha="center", color="#f38ba8")
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
                ax.text(0.5, 0.5, "plot_depth_section not available",
                        transform=ax.transAxes, ha="center")
        except Exception as exc:
            ax.text(0.5, 0.5, f"{exc}", transform=ax.transAxes, ha="center")
        return fig

    def plot_borehole_map(self, **kw) -> Figure:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(10, 4))
        self._apply_fig_style(fig, ax)
        model = self.state.model
        bhs   = self.state.boreholes
        if model is None and not bhs:
            ax.text(0.5, 0.5, "No model or boreholes loaded",
                    transform=ax.transAxes, ha="center")
            return fig
        try:
            if model is not None and hasattr(model, "station_x"):
                ax.scatter(model.station_x / 1e3,
                           np.zeros(len(model.station_x)),
                           marker="v", color="#89b4fa", s=60, label="Stations", zorder=3)
            for bh in bhs:
                x = getattr(bh, "x", 0)
                ax.axvline(x / 1e3, color="#a6e3a1", lw=1.5, ls="--", alpha=0.8)
                ax.text(x / 1e3, 0.5, bh.name, fontsize=7, ha="center",
                        color="#a6e3a1", va="bottom")
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
            log = (next((l for l in self.state.strat_logs if l.station_name == station),
                        self.state.strat_logs[0]))
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

    def plot_calibrated_model(self, **kw) -> Figure:
        if self.state.model is None:
            return self._no_model_fig()
        try:
            from pycsamt.interp.plot import (
                PlotCalibratedModel,
            )
            plotter = PlotCalibratedModel(self.state.model)
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
        rocks = list(db)
        fig, ax = plt.subplots(figsize=(10, max(4, len(rocks) * 0.4 + 1)))
        self._apply_fig_style(fig, ax)
        try:
            for i, rock in enumerate(rocks):
                color = getattr(rock, "color", "#888888")
                lo = getattr(rock, "rho_min", 0.1)
                hi = getattr(rock, "rho_max", 10000.0)
                ax.barh(i, np.log10(hi) - np.log10(lo),
                        left=np.log10(lo), color=color, alpha=0.85,
                        edgecolor="#313244", linewidth=0.5)
                ax.text(np.log10(lo) - 0.05, i, rock.name,
                        ha="right", va="center", fontsize=7)
            ax.set_xlabel("log₁₀ ρ (Ω·m)", fontsize=9)
            ax.set_title("Rock Database — Resistivity Ranges", fontsize=10)
            ax.set_yticks([])
        except Exception as exc:
            ax.text(0.5, 0.5, str(exc), transform=ax.transAxes, ha="center")
        return fig

    # ── Hydrology ─────────────────────────────────────────────────────────

    def _hydro_section_plot(self, attr: str, title: str,
                             cmap: str = "Blues_r",
                             label: str = "") -> Figure:
        if self.state.hydro_result is None:
            return self._needs_run_fig("Run Hydrology Estimation first")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(12, 5))
        self._apply_fig_style(fig, ax)
        try:
            model = self.state.model
            data  = getattr(self.state.hydro_result, attr, None)
            if data is None:
                ax.text(0.5, 0.5, f"'{attr}' not in result",
                        transform=ax.transAxes, ha="center")
                return fig
            x = model.x_centers / 1e3 if hasattr(model, "x_centers") else np.arange(data.shape[1])
            z = model.z_centers       if hasattr(model, "z_centers") else np.arange(data.shape[0])
            im = ax.pcolormesh(x, z, np.log10(np.maximum(data, 1e-20)),
                               cmap=cmap, shading="auto")
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
        return self._hydro_section_plot("hydraulic_conductivity",
                                         "Hydraulic Conductivity K  (m/s)",
                                         "viridis_r", "log₁₀ K (m/s)")

    def plot_Sw_map(self, **kw):
        return self._hydro_section_plot("saturation", "Water Saturation Sw",
                                         "Blues", "Sw")

    def plot_water_table(self, **kw) -> Figure:
        if self.state.hydro_result is None:
            return self._needs_run_fig("Run Hydrology Estimation first")
        try:
            from pycsamt.interp.plot import (
                PlotWaterTableProfile,
            )
            fig = PlotWaterTableProfile(
                self.state.model, self.state.hydro_result
            ).plot()
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
            hi   = HydroInterpreter(self.state.model, db=db)
            zones = hi.identify_aquifers()
            model = self.state.model
            x = model.x_centers / 1e3 if hasattr(model, "x_centers") else np.arange(model.n_x)
            z = model.z_centers       if hasattr(model, "z_centers") else np.arange(model.n_z)
            rho = model.rho_2d
            ax.pcolormesh(x, z, rho, cmap="Greys_r", shading="auto", alpha=0.4)
            for zone in zones:
                xi = float(getattr(zone, "x", 0)) / 1e3
                top, bot = float(zone.top), float(zone.bottom)
                conf = float(getattr(zone, "confidence", 0.8))
                rect = plt.Rectangle(
                    (xi - 0.05, top), 0.1, bot - top,
                    fc="#89b4fa", alpha=conf * 0.7, ec="#89b4fa", lw=1,
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
            T  = self.state.hydro_result.transmissivity
            model = self.state.model
            x  = model.x_centers / 1e3 if hasattr(model, "x_centers") else np.arange(len(T))
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
            fig = PlotAquiferCharacterization(
                self.state.model, self.state.hydro_result
            ).plot()
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
            fig = PlotPetrophysicalCrossPlot(
                self.state.model, self.state.hydro_result,
                self.state.petro_cfg
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

    def _emtools_plot(self, fn_name: str, title: str,
                      polar: bool = False, **kw) -> Figure:
        if self.state.sites is None:
            return self._no_sites_fig()
        import inspect

        import matplotlib.pyplot as plt

        import pycsamt.emtools as et
        fn = getattr(et, fn_name, None)
        if fn is None:
            return self._not_implemented(fn_name)
        try:
            params = inspect.signature(fn).parameters
            subplot_kw = {"projection": "polar"} if polar else {}
            fig, ax = plt.subplots(figsize=(12, 5), subplot_kw=subplot_kw)
            self._apply_fig_style(fig, ax)
            if "ax" in params:
                result = fn(self.state.sites, ax=ax, verbose=0, **kw)
            elif "axes" in params:
                result = fn(self.state.sites, axes=ax, verbose=0, **kw)
            else:
                # Function manages its own figure; close the placeholder
                plt.close(fig)
                result = fn(self.state.sites, verbose=0, **kw)
                fig = None
            if isinstance(result, plt.Figure):
                if fig is not None and result is not fig:
                    plt.close(fig)
                fig = result
                self._apply_fig_style_minimal(fig)
            return fig
        except Exception as exc:
            return self._error_fig(f"{fn_name}: {exc}")

    def plot_pt_section(self, **kw):
        return self._emtools_plot("plot_phase_tensor_psection",
                                   "Phase Tensor Pseudo-Section")

    def plot_strike_rose(self, **kw):
        return self._emtools_plot("plot_strike_rose", "Strike Rose Diagram",
                                   polar=True)

    def plot_strike_profile(self, **kw):
        return self._emtools_plot("plot_strike_profile", "Strike Profile")

    def plot_dimensionality(self, **kw):
        return self._emtools_plot("plot_dimensionality_psection",
                                   "Dimensionality Classification")

    def plot_bostick_depths(self, **kw):
        return self._emtools_plot("plot_depth_section", "Bostick Depth Coverage")

    def plot_gradient_section(self, **kw):
        return self._emtools_plot("plot_gradient_section", "Gradient Imaging (Zhang 2021)")

    def plot_phasor_wheel(self, **kw):
        return self._emtools_plot("plot_phasor_wheel", "Impedance Phasor Wheel",
                                   polar=True)

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
                self.state.model, self.state.mc_result
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
            fig = PlotUncertaintyProfile(
                self.state.model, self.state.mc_result
            ).plot()
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
        return self._emtools_plot("plot_impedance_mohr_circles", "Impedance Mohr Circles")

    def plot_argand(self, **kw):
        return self._emtools_plot("plot_zt_argand", "Impedance Argand Diagram")

    def plot_pt_ternary(self, **kw):
        return self._emtools_plot("plot_dimensionality_ternary",
                                   "Dimensionality Ternary (1D / 2D / 3D)")

    def plot_distortion_radar(self, **kw):
        return self._emtools_plot("plot_distortion_radar",
                                   "Distortion Radar: Static Shift · Twist · Shear")

    def plot_bode(self, **kw):
        return self._emtools_plot("plot_rho_phase_bode",
                                   "Bode Plot: log₁₀(ρ_a) and Phase")

    def plot_z_invariants(self, **kw):
        return self._emtools_plot("plot_z_invariants_section",
                                   "Z Tensor Invariants Section")

    def plot_anisotropy(self, **kw):
        return self._emtools_plot("plot_apparent_anisotropy_section",
                                   "Apparent Anisotropy Section")

    def plot_period_clock(self, **kw):
        return self._emtools_plot("plot_pt_period_clock", "Phase Tensor Period Clock")

    def plot_composite_section(self, **kw):
        return self._emtools_plot("plot_mt_composite_section",
                                   "MT Composite Section (multi-metric)")

    # ── Fusion & Time-Lapse ───────────────────────────────────────────────

    def plot_fused_model(self, **kw) -> Figure:
        if self.state.fusion_model is None:
            return self._needs_run_fig("Run Fusion first")
        self_copy = type("_", (), {"state": type("_", (), {"model": self.state.fusion_model})()})()
        self_copy.dark = self.dark
        self_copy._apply_fig_style = self._apply_fig_style
        self_copy._apply_fig_style_minimal = self._apply_fig_style_minimal
        self_copy._no_model_fig = self._no_model_fig
        self_copy._error_fig = self._error_fig
        self_copy._style_cb = self._style_cb
        return self.plot_model_summary.__func__(self_copy)

    def plot_timelapse_change(self, **kw) -> Figure:
        if not self.state.timelapse_surveys:
            return self._needs_run_fig("Load at least 2 surveys for time-lapse")
        try:
            from pycsamt.interp.timelapse import TimeLapseEM
            tl = TimeLapseEM(surveys=self.state.timelapse_surveys)
            changes = tl.resistivity_change()
            import matplotlib.pyplot as plt
            fig, axes = plt.subplots(1, len(changes), figsize=(14, 5))
            if len(changes) == 1:
                axes = [axes]
            for ax, dc in zip(axes, changes):
                im = ax.pcolormesh(dc, cmap="RdBu_r", shading="auto",
                                   vmin=-50, vmax=50)
                ax.invert_yaxis()
                self._apply_fig_style(fig, ax)
            fig.colorbar(im, ax=axes[-1], label="Δρ (%)").set_label(
                "Δρ (%)", fontsize=8
            )
            fig.suptitle("Time-Lapse Resistivity Change", fontsize=11)
            return fig
        except Exception as exc:
            return self._error_fig(str(exc))

    def plot_timelapse_sat(self, **kw) -> Figure:
        return self._needs_run_fig("Time-lapse ΔSw — run time-lapse analysis first")

    def plot_timelapse_wt(self, **kw) -> Figure:
        return self._needs_run_fig("Time-lapse WT displacement — run analysis first")

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
            ax.text(0.05, 0.9 - i * 0.12, ln, transform=ax.transAxes,
                    fontsize=9, color=s["fg"], family="monospace")
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
        log = (next((l for l in self.state.strat_logs if l.station_name == station),
                    self.state.strat_logs[0]))
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
                elev_m = extract_elevation(self.state.sites)
                chain  = extract_chainage(self.state.sites)
                if len(chain) > 0 and len(sx_km) > 0:
                    elev_m = interp_elev(chain, elev_m / 1000.0, sx_km) * 1000.0
            if elev_m is None or len(elev_m) == 0:
                return
            elev_m = np.asarray(elev_m)
            draw_topo_section(
                ax, sx_km, elev_m,
                dark=self.dark,
            )
        except Exception:
            pass   # topo is always optional; never break the main figure

    def _no_model_fig(self) -> Figure:
        return self._msg_fig("No resistivity model loaded.\n"
                             "Load from inversion results or from file.")

    def _no_sites_fig(self) -> Figure:
        return self._msg_fig("No survey data loaded.\n"
                             "Open EDI files from the main window first.")

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
        ax.text(0.5, 0.5, msg, transform=ax.transAxes,
                ha="center", va="center", fontsize=11,
                color=color or s["muted"], wrap=True,
                multialignment="center")
        return fig


# ── Theme dicts ────────────────────────────────────────────────────────────────

_DARK = dict(
    bg="#1e1e2e", fig_bg="#181825",
    fg="#cdd6f4", title="#cdd6f4", tick="#a6adc8",
    spine="#45475a", grid="#313244", muted="#585b70",
)
_LIGHT = dict(
    bg="#eff1f5", fig_bg="#e6e9ef",
    fg="#4c4f69", title="#4c4f69", tick="#6c6f85",
    spine="#bcc0cc", grid="#ccd0da", muted="#9ca0b0",
)
