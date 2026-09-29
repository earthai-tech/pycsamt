# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
AirborneController — drives AirborneWindow's diagnostics/map plotting.

Loads an :class:`~pycsamt.airborne.site.AirborneSites` via
:func:`~pycsamt.airborne.site.ensure_asites` (EMTF-XML only — see the
MobileMT note below) and dispatches to the relevant
``pycsamt.emtools.{ztem,afmag,mobilemt}`` plotting function. Every
function in the catalogue accepts either a single-axes ``ax=`` (most)
or builds its own multi-panel figure when passed ``axes=None``/
``axis=None`` (the two ``*_grid``/``*_band_mask*`` entries) — both
paths return a ready-to-embed ``matplotlib.figure.Figure``.

MobileMT note
-------------
``pycsamt.airborne``/``ensure_mobilemt_dataset`` only ever reads
EMTF-XML (a generic, already-decoded adapter format) — raw vendor
MobileMT files are a permanent, deliberate restriction (EGS does not
license raw-format access; see ``MEMORY.md``'s
``project_mobilemt_vendor_data_blocked.md``), not a temporary gap this
window works around. The "MobileMT" category's description and the
Load button say so explicitly.
"""

from __future__ import annotations

from collections import namedtuple
from dataclasses import dataclass, field
from typing import Any

# ── Parameter specification (same shape as CorrectionController's) ────────────
ParamSpec = namedtuple(
    "ParamSpec",
    ["name", "label", "kind", "default", "opts", "tip"],
    defaults=[""],
)

# ── Catalogue ──────────────────────────────────────────────────────────────────
# {category: [(label, fn_name, desc, [ParamSpec, ...], multi_axes), ...]}
# multi_axes=True -> function takes axes=/None and returns its own Figure;
# multi_axes=False (default) -> function takes ax=None and returns an Axes.

CATALOGUE: dict[str, list[tuple]] = {
    "ZTEM": [
        (
            "Tipper profile",
            "plot_ztem_tipper_profile",
            "Legault et al. (2012) style raw in-phase/quadrature tipper "
            "profile at one reference frequency.",
            [
                ParamSpec(
                    "component", "Component", "combo", "tzx",
                    ["tzx", "tzy"],
                ),
                ParamSpec(
                    "as_percent", "As percent", "check", True, None,
                ),
            ],
            False,
        ),
        (
            "Divergence profile",
            "plot_ztem_divergence_profile",
            "Along-line total-divergence (Peaker) profile at one "
            "reference frequency.",
            [
                ParamSpec(
                    "component", "Component", "combo", "tzx",
                    ["tzx", "tzy"],
                ),
                ParamSpec(
                    "part", "Part", "combo", "real", ["real", "imag"],
                ),
                ParamSpec(
                    "spacing_m", "Spacing (m)", "dspin", 200.0,
                    (10.0, 5000.0, 10.0),
                ),
            ],
            False,
        ),
        (
            "Divergence pseudosection",
            "plot_ztem_divergence_psection",
            "Total-divergence pseudosection (station x frequency) for "
            "one flight line.",
            [
                ParamSpec(
                    "component", "Component", "combo", "tzx",
                    ["tzx", "tzy"],
                ),
                ParamSpec(
                    "part", "Part", "combo", "real", ["real", "imag"],
                ),
                ParamSpec(
                    "spacing_m", "Spacing (m)", "dspin", 200.0,
                    (10.0, 5000.0, 10.0),
                ),
            ],
            False,
        ),
        (
            "Divergence pseudosection grid (multi-line)",
            "plot_ztem_divergence_psection_grid",
            "Every detected flight line's divergence pseudosection on "
            "one shared colour scale, side by side.",
            [
                ParamSpec(
                    "component", "Component", "combo", "tzx",
                    ["tzx", "tzy"],
                ),
                ParamSpec(
                    "part", "Part", "combo", "real", ["real", "imag"],
                ),
                ParamSpec(
                    "max_lines", "Max lines", "spin", 6, (1, 20, 1),
                ),
            ],
            True,
        ),
        (
            "Phase rotation profile",
            "plot_ztem_phase_rotation_profile",
            "Raw vs. Hilbert-transform phase-rotated response — turns "
            "a crossover anomaly into a peak anomaly (Sattel & "
            "Witherly 2012).",
            [
                ParamSpec(
                    "component", "Component", "combo", "tzx",
                    ["tzx", "tzy"],
                ),
                ParamSpec(
                    "part", "Part", "combo", "real", ["real", "imag"],
                ),
            ],
            False,
        ),
        (
            "Band-mask pseudosection",
            "plot_ztem_band_mask_psection",
            "Before/after |T| pseudosections around the ZTEM usable "
            "bandwidth (default 22-720 Hz).",
            [
                ParamSpec(
                    "component", "Component", "combo", "abs",
                    ["abs", "tzx", "tzy"],
                ),
            ],
            True,
        ),
        (
            "Flight lines map",
            "plot_ztem_flight_lines",
            "Plan-view navigation map, one coloured trace per detected "
            "flight line (Sattel & Witherly 2012, Fig. 7).",
            [],
            False,
        ),
        (
            "Map (tipper / quantity)",
            "plot_ztem_map",
            "Plan-view gridded map of tipper or a derived quantity at "
            "one reference frequency (Legault et al. 2012, Fig. 7).",
            [
                ParamSpec(
                    "quantity", "Quantity", "combo", "tipper",
                    ["tipper", "divergence"],
                ),
                ParamSpec(
                    "component", "Component", "combo", "tzx",
                    ["tzx", "tzy"],
                ),
                ParamSpec(
                    "part", "Part", "combo", "real", ["real", "imag"],
                ),
            ],
            False,
        ),
    ],
    "AFMAG": [
        (
            "Tilt profile",
            "plot_afmag_tilt_profile",
            "Classic AFMAG flight-line tilt-angle profile at one "
            "reference frequency.",
            [
                ParamSpec(
                    "component", "Component", "combo", "real",
                    ["real", "imag"],
                ),
            ],
            False,
        ),
        (
            "Tilt pseudosection",
            "plot_afmag_tilt_psection",
            "AFMAG tilt-angle pseudosection (station x log-period).",
            [
                ParamSpec(
                    "component", "Component", "combo", "resultant",
                    ["resultant", "real", "imag"],
                ),
            ],
            False,
        ),
        (
            "Tilt polar",
            "plot_afmag_tilt_polar",
            "Polar view of AFMAG tilt (azimuth and magnitude vs. "
            "period) for one station.",
            [
                ParamSpec(
                    "component", "Component", "combo", "real",
                    ["real", "imag"],
                ),
            ],
            False,
        ),
        (
            "Motion susceptibility map",
            "plot_motion_susceptibility_map",
            "Plan-view map of each station's motion-noise "
            "susceptibility score, given the survey's aircraft "
            "attitude envelope.",
            [
                ParamSpec(
                    "inclination", "Inclination (°)", "dspin", 60.0,
                    (-90.0, 90.0, 1.0),
                    "Geomagnetic field inclination at the survey.",
                ),
                ParamSpec(
                    "declination", "Declination (°)", "dspin", 0.0,
                    (-180.0, 180.0, 1.0),
                    "Geomagnetic field declination at the survey.",
                ),
                ParamSpec(
                    "roll_amplitude_deg", "Roll amplitude (°)", "dspin",
                    5.0, (0.0, 45.0, 0.5),
                    "Typical peak-to-peak aircraft roll during survey.",
                ),
                ParamSpec(
                    "pitch_amplitude_deg", "Pitch amplitude (°)", "dspin",
                    5.0, (0.0, 45.0, 0.5),
                    "Typical peak-to-peak aircraft pitch during survey.",
                ),
            ],
            False,
        ),
    ],
    "MobileMT": [
        (
            "Admittance profile",
            "plot_mobilemt_admittance_profile",
            "One admittance component along one flight line "
            "(EMTF-XML input only — see the Load button note).",
            [
                ParamSpec(
                    "component", "Component", "combo", "det",
                    ["xx", "xy", "yx", "yy", "hzx", "hzy", "det"],
                ),
                ParamSpec(
                    "part", "Part", "combo", "abs",
                    ["real", "imag", "abs"],
                ),
            ],
            False,
        ),
        (
            "Conductivity pseudosection",
            "plot_mobilemt_conductivity_psection",
            "Apparent-conductivity pseudosection for one flight line "
            "(EMTF-XML input only).",
            [
                ParamSpec(
                    "source", "Source", "combo", "theoretical",
                    ["theoretical", "native"],
                ),
            ],
            False,
        ),
        (
            "Skew profile",
            "plot_mobilemt_skew_profile",
            "Admittance skew (asymmetry) profile along one flight "
            "line (EMTF-XML input only).",
            [],
            False,
        ),
    ],
}

CATEGORIES = list(CATALOGUE.keys())

# emtools module each category dispatches into.
_MODULE_FOR_CATEGORY = {
    "ZTEM": "pycsamt.emtools.ztem",
    "AFMAG": "pycsamt.emtools.afmag",
    "MobileMT": "pycsamt.emtools.mobilemt",
}


@dataclass
class AirborneState:
    asites: Any = None  # pycsamt.airborne.site.AirborneSites
    source_path: str = ""


class AirborneController:
    """Pure-Python controller: load airborne data, dispatch plots.

    No Qt imports — every ``generate()`` call returns a plain
    ``matplotlib.figure.Figure``, ready for ``MplCanvas.show_figure()``.
    """

    def __init__(self) -> None:
        self.state = AirborneState()
        self.dark: bool = True

    # ── Loading ───────────────────────────────────────────────────────────────

    def load(self, path: str) -> int:
        """Load *path* (an EMTF-XML file or a directory of them)."""
        from pycsamt.airborne.site import ensure_asites

        asites = ensure_asites(path, recursive=True, strict=False, verbose=0)
        self.state.asites = asites
        self.state.source_path = path
        return len(asites)

    def clear(self) -> None:
        self.state = AirborneState()

    @property
    def has_data(self) -> bool:
        return self.state.asites is not None and len(self.state.asites) > 0

    # ── Dispatch ──────────────────────────────────────────────────────────────

    def generate(self, category: str, fn_name: str, **kwargs) -> Any:
        """Run *fn_name* from *category*'s emtools module; return a Figure."""
        import matplotlib.pyplot as plt

        if not self.has_data:
            return self._msg_fig("No airborne data loaded.\nUse Load EMTF-XML.")

        module_name = _MODULE_FOR_CATEGORY[category]
        multi_axes = self._is_multi_axes(category, fn_name)

        try:
            import importlib

            mod = importlib.import_module(module_name)
            fn = getattr(mod, fn_name)
            if multi_axes:
                fig = fn(self.state.asites, **kwargs)
            else:
                fig_, ax = plt.subplots(figsize=kwargs.pop("figsize", (9, 5)))
                fn(self.state.asites, ax=ax, **kwargs)
                fig = fig_
            self._apply_fig_style(fig)
            return fig
        except Exception as exc:
            return self._error_fig(str(exc))

    def _is_multi_axes(self, category: str, fn_name: str) -> bool:
        for label, name, desc, params, multi_axes in CATALOGUE[category]:
            if name == fn_name:
                return multi_axes
        return False

    # ── Styling / placeholders ───────────────────────────────────────────────

    def _apply_fig_style(self, fig) -> None:
        s = _DARK if self.dark else _LIGHT
        fig.patch.set_facecolor(s["fig_bg"])
        for ax in fig.get_axes():
            try:
                ax.set_facecolor(s["bg"])
                ax.tick_params(colors=s["tick"], labelsize=7)
                ax.xaxis.label.set_color(s["fg"])
                ax.yaxis.label.set_color(s["fg"])
                ax.title.set_color(s["title"])
                for sp in ax.spines.values():
                    sp.set_edgecolor(s["spine"])
            except Exception:
                pass

    def _msg_fig(self, msg: str, color: str | None = None):
        import matplotlib.pyplot as plt

        s = _DARK if self.dark else _LIGHT
        fig, ax = plt.subplots(figsize=(8, 4))
        fig.patch.set_facecolor(s["fig_bg"])
        ax.set_facecolor(s["bg"])
        ax.set_axis_off()
        ax.text(
            0.5, 0.5, msg,
            transform=ax.transAxes, ha="center", va="center",
            color=color or s["muted"], fontsize=11,
        )
        return fig

    def _error_fig(self, msg: str):
        return self._msg_fig(f"✕  {msg}", color="#f38ba8")


_DARK = dict(
    fig_bg="#1e1e2e", bg="#181825", fg="#cdd6f4", tick="#a6adc8",
    spine="#45475a", title="#cdd6f4", muted="#6c7086",
)
_LIGHT = dict(
    fig_bg="#ffffff", bg="#ffffff", fg="#4c4f69", tick="#6c6f85",
    spine="#ccd0da", title="#4c4f69", muted="#9ca0b0",
)
