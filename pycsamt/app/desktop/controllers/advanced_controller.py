# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
AdvancedController — renders advanced emtools plots into matplotlib figures.

Handles two special cases transparently:
  * Polar-projection functions (rose diagrams, phasor wheel, etc.) — creates
    a polar subplot automatically.
  * Figure-creating functions (no ax= parameter) — captures the returned
    figure and transplants its axes into the target canvas figure.
"""
from __future__ import annotations

from typing import Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np

from pycsamt.app.desktop.controllers.qc_controller import (
    _style_all_axes,
    _annotate_empty,
)

# Functions that internally call set_theta_zero_location() or other
# polar-specific methods — must receive a polar subplot.
_POLAR_FNS = {
    "plot_phase_tensor_rose",
    "plot_phasor_wheel",
    "plot_induction_rose",
    "plot_induction_rose_from_spectra",
    "plot_strike_rose",
    "plot_strike_rose_by_line",
    "plot_theta_rose_grid",
    "plot_tipper_polar",
    "plot_tipper_polar_from_spectra",
    "plot_radiation_pattern",
    "plot_polar_errors",
}

# Functions that require extra positional arguments beyond (sites, …, verbose=)
# and are NOT handled by special dispatch logic below.  Extend this dict for
# any future catalogue additions in the same situation.
_EXTRA_ARGS_FUNS: dict[str, str] = {}

# ── Plot catalogue ─────────────────────────────────────────────────────────────
# Each entry: (display_label, fn_name, has_ax_param)

STRIKE_PLOTS: list = [
    ("Strike rose diagram",        "plot_strike_rose",             False),
    ("Strike rose by profile line","plot_strike_rose_by_line",     False),
    ("Strike analysis (full)",     "plot_strike_analysis",         False),
    ("θ rose grid",                "plot_theta_rose_grid",         False),
    ("θ vs period",                "plot_theta_vs_period",         True),
    ("θ stability stripe",         "plot_theta_stability_stripe",  True),
    ("Strike ribbon",              "plot_strike_ribbon",           True),
    ("Strike mapsticks",           "plot_strike_mapsticks",        True),
    ("Strike stability bands",     "plot_strike_stability_bands",  False),
]

PHASE_TENSOR_PLOTS: list = [
    ("Phase tensor map",           "plot_phase_tensor_map",        True),
    ("Phase tensor rose",          "plot_phase_tensor_rose",       True),
    ("Phase tensor summary",       "plot_phase_tensor_summary",    False),
    ("Phasor wheel",               "plot_phasor_wheel",            True),
    ("PT period clock",            "plot_pt_period_clock",         False),
    ("PT skew density",            "plot_skew_ellipt_density",     True),
    ("Dimensionality ternary",     "plot_dimensionality_ternary",  False),
]

INDUCTION_PLOTS: list = [
    ("Induction arrows",           "plot_induction_arrows",        True),
    ("Induction map",              "plot_induction_map",           True),
    ("Induction section",          "plot_induction_section",       True),
    ("Induction rose",             "plot_induction_rose",          True),
    ("Tipper polar",               "plot_tipper_polar",            True),
    ("Tipper hodograms",           "plot_tipper_hodograms",        False),
    ("Induction multiperiod map",  "plot_induction_multiperiod_map", False),
    ("Induction convention",       "plot_induction_convention",    False),
]

IMPEDANCE_PLOTS: list = [
    ("Mohr circles",               "plot_impedance_mohr_circles",  False),
    ("Argand diagram",             "plot_zt_argand",               False),
    ("Z invariants section",       "plot_z_invariants_section",    False),
    ("ρₐ/φ Bode plot",            "plot_rho_phase_bode",          False),
    ("ρₐ polar",                  "plot_apparent_resistivity_polar", False),
    ("XY/YX crossover map",        "plot_xyyx_crossover_map",      True),
    ("Off-diagonal antisym",       "plot_offdiag_antisym_residual", True),
    ("Anisotropy section",         "plot_anisotropy",              True),
    ("Ellipticity psection",       "plot_ellipticity_psection",    True),
]

DEPTH_PLOTS: list = [
    ("Depth section",              "plot_depth_section",           True),
    ("Apparent depth psection",    "plot_apparent_depth_psection", True),
    ("Gradient section",           "plot_gradient_section",        True),
    ("Sensitivity depth section",  "plot_sensitivity_depth_section", False),
    ("MT composite section",       "plot_mt_composite_section",    False),
    ("ATOM psection",              "plot_atom_psection",           True),
]

SURVEY_PLOTS: list = [
    ("Survey fingerprint",         "plot_survey_fingerprint",      False),
    ("Distortion radar",           "plot_distortion_radar",        False),
    ("Sites comparison",           "plot_sites_compare",           False),
    ("Station response panels",    "plot_sites_panels",            False),
    ("Normalised response",        "plot_normalized_response",     False),
    ("TF coherence network",       "plot_tf_coherence_network",    False),
    ("Dimensionality map",         "plot_dim_map",                 True),
    ("Occupancy area",             "plot_dim_occupancy_area",      True),
]

ADVANCED_GROUPS: list = [
    ("Strike Analysis",    STRIKE_PLOTS),
    ("Phase Tensor",       PHASE_TENSOR_PLOTS),
    ("Induction / Tipper", INDUCTION_PLOTS),
    ("Impedance / Z",      IMPEDANCE_PLOTS),
    ("Depth Imaging",      DEPTH_PLOTS),
    ("Survey Tools",       SURVEY_PLOTS),
    ("Topography",         []),   # empty — specialized UI
    ("Conversion",         []),   # empty — specialized UI
]

# Indices for the specialized panels
TOPO_INDEX = 6
CONV_INDEX = 7


# ── Controller ────────────────────────────────────────────────────────────────

class AdvancedController:
    """
    Pure-Python controller that renders advanced emtools plots into
    matplotlib Figure objects provided by MplCanvas.

    Attributes
    ----------
    dark : bool
        Apply dark-mode styling (default True).
    """

    def __init__(self) -> None:
        self._sites = None
        self.dark: bool = True
        self._dim_model: dict | None = None  # produced by train_dim_model()

    def set_sites(self, sites) -> None:
        self._sites = sites

    def clear(self) -> None:
        self._sites = None
        self._dim_model = None

    def train_dim_model(self, n_atoms: int = 6, n_iter: int = 40) -> dict:
        """
        Train a sparse-coding dictionary from the loaded sites and store it.

        Returns the model dict (keys: D, A, mu, sd, feat, meta).
        Raises ValueError if no sites are loaded.
        """
        import pycsamt.emtools as et
        if self._sites is None:
            raise ValueError("No survey data loaded — load EDIs first.")
        model = et.learn_dim_dictionary(
            self._sites,
            n_atoms=n_atoms,
            n_iter=n_iter,
            verbose=0,
        )
        self._dim_model = model
        return model

    # ── Main entry-point ──────────────────────────────────────────────────────

    def draw(self, fn_name: str, has_ax: bool, fig) -> Optional[object]:
        """
        Render *fn_name*.

        Parameters
        ----------
        fn_name : str
            Name of the emtools function.
        has_ax : bool
            True  → draw into *fig* in-place; return None.
            False → function creates its own figure; return it so the caller
                    can hand it to MplCanvas.show_figure().
        fig : matplotlib.figure.Figure
            Target figure (mutated when has_ax=True or on error).

        Returns
        -------
        Figure or None
        """
        import pycsamt.emtools as et

        fig.clear()

        if self._sites is None:
            ax = fig.add_subplot(111)
            _annotate_empty(ax, "Load survey data first")
            _style_all_axes(fig, self.dark)
            return None

        fn = getattr(et, fn_name, None)
        if fn is None:
            ax = fig.add_subplot(111)
            _annotate_empty(ax, f"Function not found:\n{fn_name}")
            _style_all_axes(fig, self.dark)
            return None

        if fn_name in _EXTRA_ARGS_FUNS:
            ax = fig.add_subplot(111)
            _annotate_empty(ax, f"Cannot render from catalogue:\n{_EXTRA_ARGS_FUNS[fn_name]}")
            _style_all_axes(fig, self.dark)
            return None

        # plot_atom_psection needs a trained dictionary model — dispatch specially.
        if fn_name == "plot_atom_psection":
            self._draw_atom_psection(fn, fig)
            return None

        try:
            if has_ax:
                if fn_name in _POLAR_FNS:
                    ax = fig.add_subplot(111, projection="polar")
                else:
                    ax = fig.add_subplot(111)
                fn(self._sites, ax=ax, verbose=0)
            else:
                src_fig = self._call_figure_fn(fn)
                if src_fig is None:
                    ax = fig.add_subplot(111)
                    _annotate_empty(ax, "No figure produced")
                else:
                    _style_all_axes(src_fig, self.dark)
                    try:
                        src_fig.tight_layout(pad=1.2)
                    except Exception:
                        pass
                    return src_fig
        except Exception as exc:
            fig.clear()
            ax = fig.add_subplot(111)
            _annotate_empty(ax, f"{fn_name} error:\n{exc}")

        _style_all_axes(fig, self.dark)
        try:
            fig.tight_layout(pad=1.2)
        except Exception:
            pass
        return None

    # ── Helpers ───────────────────────────────────────────────────────────────

    def _draw_atom_psection(self, fn, fig) -> None:
        """Render plot_atom_psection using the stored dim model."""
        if self._dim_model is None:
            ax = fig.add_subplot(111)
            _annotate_empty(
                ax,
                "No dictionary model trained yet.\n"
                "Set  n_atoms / n_iter  in the params panel\n"
                "and click  'Train Model', then Run again.",
            )
            _style_all_axes(fig, self.dark)
            return
        try:
            ax = fig.add_subplot(111)
            fn(self._sites, self._dim_model, ax=ax, verbose=0)
        except Exception as exc:
            fig.clear()
            ax = fig.add_subplot(111)
            _annotate_empty(ax, f"plot_atom_psection error:\n{exc}")
        _style_all_axes(fig, self.dark)
        try:
            fig.tight_layout(pad=1.2)
        except Exception:
            pass

    def _call_figure_fn(self, fn) -> Optional[object]:
        """Call a multi-axes function and return the Figure it creates."""
        before = set(plt.get_fignums())
        result = fn(self._sites, verbose=0)
        after  = set(plt.get_fignums())
        if hasattr(result, "get_axes"):
            return result
        new_nums = after - before
        if new_nums:
            return plt.figure(max(new_nums))
        return None


# ── Topography preview controller ─────────────────────────────────────────────

class TopoPreviewController:
    """
    Controller for the Topography panel previews.
    Renders elevation profile, terrain fill preview, and elevation histogram
    into matplotlib figures.
    """

    def __init__(self) -> None:
        self._sites = None
        self.dark: bool = True

    def set_sites(self, sites) -> None:
        self._sites = sites

    def _get_style(self) -> dict:
        if self.dark:
            return {
                "bg": "#181825",
                "fig_bg": "#1e1e2e",
                "fg": "#cdd6f4",
                "title": "#cdd6f4",
                "tick": "#a6adc8",
                "spine": "#45475a",
                "grid": "#313244",
            }
        return {
            "bg": "#eff1f5",
            "fig_bg": "#e6e9ef",
            "fg": "#4c4f69",
            "title": "#4c4f69",
            "tick": "#6c6f85",
            "spine": "#bcc0cc",
            "grid": "#ccd0da",
        }

    def plot_elevation_profile(self, fig) -> None:
        """Chainage vs elevation filled-area profile with station markers."""
        s = self._get_style()
        fig.clear()
        fig.patch.set_facecolor(s["fig_bg"])
        ax = fig.add_subplot(111)
        ax.set_facecolor(s["bg"])

        try:
            from pycsamt.topo.extract import (
                extract_chainage,
                extract_elevation,
                extract_station_names,
            )
            from pycsamt.topo.config import PYCSAMT_TOPO

            if self._sites is None:
                ax.text(
                    0.5, 0.5, "Load survey data first",
                    ha="center", va="center",
                    transform=ax.transAxes,
                    color=s["fg"], fontsize=11,
                )
                self._style_ax(ax, s)
                return

            chain = extract_chainage(self._sites)   # km
            elev  = extract_elevation(self._sites)  # m
            names = extract_station_names(self._sites)

            if chain.size == 0 or elev.size == 0:
                ax.text(
                    0.5, 0.5, "No elevation data available",
                    ha="center", va="center",
                    transform=ax.transAxes,
                    color=s["fg"], fontsize=11,
                )
                self._style_ax(ax, s)
                return

            fill_color = getattr(PYCSAMT_TOPO, "fill_color", "#a89070")
            fill_alpha = getattr(PYCSAMT_TOPO, "fill_alpha", 0.4)
            line_color = getattr(PYCSAMT_TOPO, "line_color", "#6b4e2a")
            line_width = getattr(PYCSAMT_TOPO, "line_width", 1.2)

            # Filled area from min elevation
            elev_min = float(np.nanmin(elev))
            ax.fill_between(
                chain, elev_min, elev,
                color=fill_color, alpha=fill_alpha,
                label="Terrain",
            )
            if getattr(PYCSAMT_TOPO, "show_surface_line", True):
                ax.plot(chain, elev, color=line_color, lw=line_width)

            # Station markers as ▼ at the top of the profile
            elev_top = (float(np.nanmax(elev))
                        + 0.02 * (float(np.nanmax(elev)) - elev_min + 1))
            ax.scatter(
                chain, np.full_like(chain, elev_top),
                marker="v", color=line_color, s=30, zorder=5,
                label="Stations",
            )

            ax.set_xlabel("Chainage (km)", fontsize=8, color=s["fg"])
            ax.set_ylabel("Elevation (m a.s.l.)", fontsize=8, color=s["fg"])
            ax.set_title("Elevation Profile", fontsize=9,
                         color=s["title"], pad=5)
            ax.legend(fontsize=7, facecolor=s["bg"], labelcolor=s["fg"],
                      edgecolor=s["spine"])

        except Exception as exc:
            ax.text(
                0.5, 0.5, f"Error: {exc}",
                ha="center", va="center",
                transform=ax.transAxes,
                color=s["fg"], fontsize=9,
            )

        self._style_ax(ax, s)
        try:
            fig.tight_layout(pad=1.2)
        except Exception:
            pass

    def plot_fill_preview(self, fig) -> None:
        """Synthetic section with terrain overlay — preview of current config."""
        s = self._get_style()
        fig.clear()
        fig.patch.set_facecolor(s["fig_bg"])
        ax = fig.add_subplot(111)
        ax.set_facecolor(s["bg"])

        try:
            from pycsamt.topo.config import PYCSAMT_TOPO

            # Determine number of stations
            if self._sites is not None:
                try:
                    from pycsamt.topo.extract import extract_chainage
                    chain = extract_chainage(self._sites)
                    n_st = max(int(chain.size), 4)
                except Exception:
                    n_st = 8
            else:
                n_st = 8

            # Synthetic resistivity grid: 64 depth layers × n_st stations
            n_depth = 64
            rng = np.random.default_rng(42)
            rho = np.exp(
                rng.uniform(np.log(1), np.log(1000), size=(n_depth, n_st))
            )

            # Synthetic chainage and terrain
            chain_km = np.linspace(0, max(n_st * 0.5, 1), n_st)
            elev_m   = 200.0 + 50.0 * np.sin(np.linspace(0, 2 * np.pi, n_st))

            # x and z node grids for pcolormesh
            x_edges  = np.linspace(chain_km[0], chain_km[-1], n_st + 1)
            depth_km = np.linspace(0, 2.0, n_depth + 1)
            X, Z = np.meshgrid(x_edges, depth_km)

            pcm = ax.pcolormesh(
                X, -Z,
                np.log10(np.vstack([rho, rho[-1:]])),
                cmap="jet_r", shading="flat",
                vmin=np.log10(1), vmax=np.log10(1000),
            )
            fig.colorbar(pcm, ax=ax, label="log₁₀(ρ) [Ω·m]",
                         fraction=0.03, pad=0.02)

            fill_color = getattr(PYCSAMT_TOPO, "fill_color", "#a89070")
            fill_alpha = getattr(PYCSAMT_TOPO, "fill_alpha", 0.4)
            line_color = getattr(PYCSAMT_TOPO, "line_color", "#6b4e2a")
            line_width = getattr(PYCSAMT_TOPO, "line_width", 1.2)

            # Overlay terrain
            elev_km = elev_m / 1000.0
            ax.fill_between(
                chain_km, 0, elev_km,
                color=fill_color, alpha=fill_alpha, zorder=3,
            )
            if getattr(PYCSAMT_TOPO, "show_surface_line", True):
                ax.plot(chain_km, elev_km, color=line_color,
                        lw=line_width, zorder=4)

            ax.set_xlabel("Chainage (km)", fontsize=8, color=s["fg"])
            ax.set_ylabel("Depth (km)", fontsize=8, color=s["fg"])
            ax.set_title(
                f"Terrain Fill Preview  (synthetic data, "
                f"config: {PYCSAMT_TOPO.summary()})",
                fontsize=8, color=s["title"], pad=5,
            )

        except Exception as exc:
            ax.text(
                0.5, 0.5, f"Error: {exc}",
                ha="center", va="center",
                transform=ax.transAxes,
                color=s["fg"], fontsize=9,
            )

        self._style_ax(ax, s)
        try:
            fig.tight_layout(pad=1.2)
        except Exception:
            pass

    def plot_elevation_histogram(self, fig) -> None:
        """Bar histogram of station elevations."""
        s = self._get_style()
        fig.clear()
        fig.patch.set_facecolor(s["fig_bg"])
        ax = fig.add_subplot(111)
        ax.set_facecolor(s["bg"])

        try:
            from pycsamt.topo.extract import extract_elevation
            from pycsamt.topo.config import PYCSAMT_TOPO

            if self._sites is None:
                ax.text(
                    0.5, 0.5, "Load survey data first",
                    ha="center", va="center",
                    transform=ax.transAxes,
                    color=s["fg"], fontsize=11,
                )
                self._style_ax(ax, s)
                return

            elev = extract_elevation(self._sites)
            if elev.size == 0:
                ax.text(
                    0.5, 0.5, "No elevation data available",
                    ha="center", va="center",
                    transform=ax.transAxes,
                    color=s["fg"], fontsize=11,
                )
                self._style_ax(ax, s)
                return

            fill_color = getattr(PYCSAMT_TOPO, "fill_color", "#a89070")
            line_color = getattr(PYCSAMT_TOPO, "line_color", "#6b4e2a")

            ax.bar(
                np.arange(len(elev)), elev,
                color=fill_color, edgecolor=line_color, linewidth=0.7,
            )
            ax.set_xlabel("Station index", fontsize=8, color=s["fg"])
            ax.set_ylabel("Elevation (m a.s.l.)", fontsize=8, color=s["fg"])
            ax.set_title("Station Elevation Histogram", fontsize=9,
                         color=s["title"], pad=5)

        except Exception as exc:
            ax.text(
                0.5, 0.5, f"Error: {exc}",
                ha="center", va="center",
                transform=ax.transAxes,
                color=s["fg"], fontsize=9,
            )

        self._style_ax(ax, s)
        try:
            fig.tight_layout(pad=1.2)
        except Exception:
            pass

    def get_stats(self) -> dict:
        """Return a dict of topo/site statistics."""
        from pycsamt.topo.config import PYCSAMT_TOPO

        result = {
            "n_stations":    0,
            "has_elev":      False,
            "elev_min":      float("nan"),
            "elev_max":      float("nan"),
            "elev_mean":     float("nan"),
            "elev_std":      float("nan"),
            "chainage_km":   float("nan"),
            "config_summary": PYCSAMT_TOPO.summary(),
        }
        if self._sites is None:
            return result

        try:
            from pycsamt.topo.extract import (
                extract_elevation,
                extract_chainage,
                has_elevation,
            )
            elev  = extract_elevation(self._sites)
            chain = extract_chainage(self._sites)
            result["n_stations"] = int(elev.size)
            result["has_elev"]   = bool(has_elevation(self._sites))
            if elev.size > 0:
                result["elev_min"]  = float(np.nanmin(elev))
                result["elev_max"]  = float(np.nanmax(elev))
                result["elev_mean"] = float(np.nanmean(elev))
                result["elev_std"]  = float(np.nanstd(elev))
            if chain.size > 0:
                result["chainage_km"] = float(chain[-1] - chain[0])
        except Exception:
            pass

        return result

    # ── Internal helper ───────────────────────────────────────────────────────

    def _style_ax(self, ax, s: dict) -> None:
        ax.tick_params(colors=s["tick"], labelsize=7)
        for sp in ax.spines.values():
            sp.set_edgecolor(s["spine"])
        ax.xaxis.label.set_color(s["fg"])
        ax.yaxis.label.set_color(s["fg"])
        ax.grid(True, color=s["grid"], alpha=0.3, ls="--", lw=0.5)


# ── Qt guard for QThread ──────────────────────────────────────────────────────

try:
    from PySide6.QtCore import QThread
    from PySide6.QtCore import Signal as QtSignal
    _HAS_QT = True
except ImportError:
    QThread = object  # type: ignore[misc,assignment]
    def QtSignal(*args):  # type: ignore[misc]
        return None
    _HAS_QT = False


# ── Dictionary-model training worker ─────────────────────────────────────────

class DimModelWorker(QThread):  # type: ignore[misc]
    """QThread worker that trains the ATOM dictionary model off the GUI thread."""

    if _HAS_QT:
        finished = QtSignal(dict)
        error    = QtSignal(str)
    else:
        finished = None  # type: ignore[assignment]
        error    = None  # type: ignore[assignment]

    def __init__(self,
                 ctrl: AdvancedController,
                 n_atoms: int,
                 n_iter: int) -> None:
        super().__init__()
        self._ctrl    = ctrl
        self._n_atoms = n_atoms
        self._n_iter  = n_iter

    def run(self) -> None:
        try:
            model = self._ctrl.train_dim_model(self._n_atoms, self._n_iter)
            if self.finished is not None:
                self.finished.emit(model)
        except Exception as exc:
            if self.error is not None:
                self.error.emit(str(exc))


# ── Conversion controller ─────────────────────────────────────────────────────

class ConversionController:
    """
    Controller for running file-format conversions (AVG/J/Spectra -> EDI).
    """

    SOURCES = ["AVG -> EDI", "J -> EDI", "Spectra -> EDI"]

    def __init__(self) -> None:
        self._source_type: str = "AVG -> EDI"
        self._source_path: str = ""
        self._result = None          # EDICollection after successful transform
        self._stats: dict = {}
        self._failures: list = []
        self.dark: bool = True

    @property
    def result(self):
        return self._result

    @property
    def has_result(self) -> bool:
        return self._result is not None and len(self._result) > 0

    def set_source(self, type_str: str, path: str) -> None:
        self._source_type = type_str
        self._source_path = path

    def run(self, options: dict):
        """
        Run transformation synchronously.

        Returns
        -------
        (collection, failures) : tuple
        """
        path = self._source_path
        src  = self._source_type
        failures: list = []
        collection = None

        if src in ("AVG -> EDI", "AVG → EDI"):
            from pycsamt.transformers import AVGtoEDI
            collection = AVGtoEDI().transform(path, **options)

        elif src in ("J -> EDI", "J → EDI"):
            from pycsamt.transformers import JtoEDI
            collection = JtoEDI().transform(path, **options)

        elif src in ("Spectra -> EDI", "Spectra → EDI"):
            from pycsamt.transformers import SpectraToEDI
            result     = SpectraToEDI(**options).transform_batch(path)
            collection = result.collection
            failures   = getattr(result, "failures", [])

        else:
            raise ValueError(f"Unknown source type: {src!r}")

        self._result   = collection
        self._failures = failures
        return collection, failures

    def build_stats(self, collection, failures: list) -> dict:
        """Build per-station stats dict from collection."""
        stats: dict = {
            "rows":       [],
            "n_total":    0,
            "n_failures": len(failures) if failures else 0,
        }
        if collection is None:
            return stats

        try:
            from pycsamt.emtools._core import _iter_items, _get_z_block, _name

            rows = []
            for i, ed in enumerate(_iter_items(collection)):
                station = _name(ed, i)
                try:
                    _, z, freqs = _get_z_block(ed)
                    has_z   = z is not None and z.size > 0
                    n_freqs = int(freqs.size) if freqs is not None else 0
                    f_min   = float(np.nanmin(freqs)) if n_freqs > 0 else float("nan")
                    f_max   = float(np.nanmax(freqs)) if n_freqs > 0 else float("nan")
                except Exception:
                    has_z   = False
                    n_freqs = 0
                    f_min   = float("nan")
                    f_max   = float("nan")

                # Check for tipper
                has_tipper = False
                try:
                    t = (getattr(ed, "Tipper", None)
                         or getattr(ed, "tipper", None))
                    if t is not None:
                        tm = getattr(t, "tipper", None)
                        has_tipper = (tm is not None
                                      and np.any(np.isfinite(tm)))
                except Exception:
                    pass

                rows.append({
                    "station":    station,
                    "n_freqs":    n_freqs,
                    "f_min":      f_min,
                    "f_max":      f_max,
                    "has_Z":      has_z,
                    "has_tipper": has_tipper,
                })
            stats["rows"]    = rows
            stats["n_total"] = len(rows)
        except Exception:
            pass

        return stats

    def _get_style(self) -> dict:
        if self.dark:
            return {
                "bg": "#181825", "fig_bg": "#1e1e2e",
                "fg": "#cdd6f4", "title": "#cdd6f4",
                "tick": "#a6adc8", "spine": "#45475a", "grid": "#313244",
            }
        return {
            "bg": "#eff1f5", "fig_bg": "#e6e9ef",
            "fg": "#4c4f69", "title": "#4c4f69",
            "tick": "#6c6f85", "spine": "#bcc0cc", "grid": "#ccd0da",
        }

    def plot_impedance_curves(self, fig) -> None:
        """Multi-station rho_a XY vs period curves."""
        s = self._get_style()
        fig.clear()
        fig.patch.set_facecolor(s["fig_bg"])
        ax = fig.add_subplot(111)
        ax.set_facecolor(s["bg"])

        if not self.has_result:
            ax.text(
                0.5, 0.5, "No conversion result available",
                ha="center", va="center",
                transform=ax.transAxes,
                color=s["fg"], fontsize=11,
            )
            self._style_ax(ax, s)
            return

        try:
            from pycsamt.emtools._core import _iter_items, _get_z_block
            _MU0 = 4.0 * np.pi * 1e-7
            items = list(_iter_items(self._result))
            n = len(items)
            if n == 0:
                ax.text(0.5, 0.5, "Empty collection",
                        ha="center", va="center",
                        transform=ax.transAxes,
                        color=s["fg"], fontsize=11)
                self._style_ax(ax, s)
                return

            cmap = plt.cm.turbo
            plotted = 0
            for idx, ed in enumerate(items):
                _, z, freqs = _get_z_block(ed)
                if z is None or freqs is None or freqs.size == 0:
                    continue
                Z_xy  = z[:, 0, 1]
                omega = 2.0 * np.pi * np.maximum(freqs, 1e-30)
                rho   = np.abs(Z_xy) ** 2 / (omega * _MU0)
                T     = 1.0 / np.maximum(freqs, 1e-30)
                ok    = np.isfinite(rho) & (rho > 0) & np.isfinite(T)
                if ok.sum() < 2:
                    continue
                ax.loglog(T[ok], rho[ok],
                          color=cmap(idx / max(n - 1, 1)),
                          alpha=0.55, lw=0.9)
                plotted += 1

            if plotted == 0:
                ax.text(0.5, 0.5, "No valid Z data",
                        ha="center", va="center",
                        transform=ax.transAxes,
                        color=s["fg"], fontsize=11)
                self._style_ax(ax, s)
                return

        except Exception as exc:
            ax.text(0.5, 0.5, f"Error: {exc}",
                    ha="center", va="center",
                    transform=ax.transAxes,
                    color=s["fg"], fontsize=9)
            self._style_ax(ax, s)
            return

        ax.set_xlabel("Period (s)", fontsize=8, color=s["fg"])
        ax.set_ylabel("ρ_a  [Ω·m]", fontsize=8, color=s["fg"])
        ax.set_title("Apparent Resistivity XY", fontsize=9,
                     color=s["title"], pad=5)
        self._style_ax(ax, s)
        try:
            fig.tight_layout(pad=1.2)
        except Exception:
            pass

    def plot_station_map(self, fig) -> None:
        """Lat/lon scatter of converted stations."""
        s = self._get_style()
        fig.clear()
        fig.patch.set_facecolor(s["fig_bg"])
        ax = fig.add_subplot(111)
        ax.set_facecolor(s["bg"])

        if not self.has_result:
            ax.text(
                0.5, 0.5, "No conversion result available",
                ha="center", va="center",
                transform=ax.transAxes,
                color=s["fg"], fontsize=11,
            )
            self._style_ax(ax, s)
            return

        try:
            from pycsamt.emtools._core import _iter_items, _name

            lats, lons, names_list = [], [], []
            for i, ed in enumerate(_iter_items(self._result)):
                try:
                    head = (getattr(ed, "Head", None)
                            or getattr(ed, "head", None))
                    if head is None:
                        head = ed
                    lat = float(
                        getattr(head, "lat", None)
                        or getattr(head, "latitude", None)
                        or 0.0
                    )
                    lon = float(
                        getattr(head, "lon", None)
                        or getattr(head, "longitude", None)
                        or 0.0
                    )
                    lats.append(lat)
                    lons.append(lon)
                    names_list.append(_name(ed, i))
                except Exception:
                    pass

            if not lats:
                ax.text(0.5, 0.5, "No coordinate data",
                        ha="center", va="center",
                        transform=ax.transAxes,
                        color=s["fg"], fontsize=11)
                self._style_ax(ax, s)
                return

            lats_arr = np.array(lats)
            lons_arr = np.array(lons)
            sc = ax.scatter(
                lons_arr, lats_arr,
                c=np.arange(len(lats)), cmap="turbo",
                s=40, zorder=5,
                edgecolors=s["spine"], linewidths=0.5,
            )
            fig.colorbar(sc, ax=ax, label="Station index",
                         fraction=0.03, pad=0.02)
            for i, nm in enumerate(names_list):
                ax.annotate(nm, (lons_arr[i], lats_arr[i]),
                            fontsize=5, color=s["fg"],
                            xytext=(3, 3), textcoords="offset points")

        except Exception as exc:
            ax.text(0.5, 0.5, f"Error: {exc}",
                    ha="center", va="center",
                    transform=ax.transAxes,
                    color=s["fg"], fontsize=9)
            self._style_ax(ax, s)
            return

        ax.set_xlabel("Longitude (°)", fontsize=8, color=s["fg"])
        ax.set_ylabel("Latitude (°)", fontsize=8, color=s["fg"])
        ax.set_title("Station Map", fontsize=9, color=s["title"], pad=5)
        self._style_ax(ax, s)
        try:
            fig.tight_layout(pad=1.2)
        except Exception:
            pass

    def _style_ax(self, ax, s: dict) -> None:
        ax.tick_params(colors=s["tick"], labelsize=7)
        for sp in ax.spines.values():
            sp.set_edgecolor(s["spine"])
        ax.xaxis.label.set_color(s["fg"])
        ax.yaxis.label.set_color(s["fg"])
        ax.grid(True, color=s["grid"], alpha=0.3, ls="--", lw=0.5)


# ── ConversionWorker ──────────────────────────────────────────────────────────

class ConversionWorker(QThread):  # type: ignore[misc]
    """QThread worker that runs ConversionController.run() off the GUI thread."""

    if _HAS_QT:
        finished = QtSignal(object, list)   # (collection, failures)
        error    = QtSignal(str)
    else:
        finished = None  # type: ignore[assignment]
        error    = None  # type: ignore[assignment]

    def __init__(self, ctrl: ConversionController, options: dict) -> None:
        super().__init__()
        self._ctrl    = ctrl
        self._options = options

    def run(self) -> None:
        try:
            col, failures = self._ctrl.run(self._options)
            if self.finished is not None:
                self.finished.emit(col, failures)
        except Exception as exc:
            if self.error is not None:
                self.error.emit(str(exc))
