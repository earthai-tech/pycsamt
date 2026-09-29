# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
SettingsController — pure-Python controller for all PYCSAMT_* API singletons.

Responsibilities
----------------
* Per-tab apply / reset methods called by each SettingsPage.
* Snapshot and restore all singletons (Cancel support in APIConfigDialog).
* Serialise / deserialise settings to JSON (Save Profile / Load Profile).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any


class SettingsController:
    """
    Mediates between the settings UI and the PYCSAMT_* API singletons.

    Each ``apply_<tab>(**kw)`` method writes only the keys present in *kw*,
    so pages can pass a partial dict.  ``snapshot()`` / ``restore()`` use the
    same key-space, making JSON persistence trivial.
    """

    SETTINGS_PATH = Path.home() / ".pycsamt" / "settings.json"

    # ── Snapshot / restore ────────────────────────────────────────────────────

    def snapshot(self) -> dict[str, Any]:
        """Return a JSON-serialisable snapshot of all configurable singletons."""
        snap: dict[str, Any] = {}

        # View controls
        try:
            from pycsamt.api.control import (
                PYCSAMT_CONTROL as C,
            )

            snap["view_controls"] = {
                "panel_background": C.panel.background,
                "panel_toolbar": C.panel.toolbar,
                "panel_fit_layout": C.panel.fit_layout,
                "rho_view": C.rho.view,
                "phase_range": list(C.phase.range),
                "phase_unit": C.phase.unit,
                "phase_wrap": C.phase.wrap,
                "x_view": C.x.view,
            }
        except Exception:
            pass

        # Plot/export defaults (v2.6)
        try:
            from pycsamt.api.plot import PLOT_CONFIG as P

            snap["plot"] = {
                "fmt": P.fmt,
                "dpi": P.dpi,
                "bbox_inches": P.bbox_inches,
                "savedir": None if P.savedir is None else str(P.savedir),
                "transparent": P.transparent,
                "close_after_save": P.close_after_save,
            }
        except Exception:
            pass

        # Shared contour and mesh rendering (v2.6)
        try:
            from pycsamt.api.contour import PYCSAMT_CONTOUR as C

            snap["contour"] = {
                "enabled": C.default.enabled,
                "levels": C.default.levels,
                "linewidths": C.default.linewidths,
                "alpha": C.default.alpha,
                "labels": C.default.labels,
            }
        except Exception:
            pass
        try:
            from pycsamt.api.mesh import PYCSAMT_MESH as M

            edge = M.review.edge
            snap["mesh"] = {
                "review__edge__show": edge.show,
                "review__edge__linewidth": edge.linewidth,
                "review__edge__alpha": edge.alpha,
                "review__edge__linestyle": edge.linestyle,
            }
        except Exception:
            pass

        # Processing-pipeline engine (v2.3+)
        try:
            from pycsamt.api.pipe import PYCSAMT_PIPE as PP

            snap["pipe"] = {
                "output_root": PP.output_root,
                "processed_subdir": PP.processed_subdir,
                "plots_subdir": PP.plots_subdir,
                "on_step_error": PP.on_step_error,
                "save_intermediate": PP.save_intermediate,
                "plot_dpi": PP.plot_dpi,
                "plot_fmt": PP.plot_fmt,
                "report_formats": list(PP.report_formats),
                "cache_root": PP.cache_root,
                "history_path": PP.history_path,
            }
        except Exception:
            pass

        # Package-wide site ordering (v2.6)
        try:
            from pycsamt.api.ordering import PYCSAMT_ORDERING as O

            snap["ordering"] = {
                "mode": O.mode,
                "min_linearity": O.min_linearity,
                "max_cross_track_ratio": O.max_cross_track_ratio,
                "min_coordinate_fraction": O.min_coordinate_fraction,
            }
        except Exception:
            pass

        # Station rendering (pseudosection preset)
        try:
            from pycsamt.api.station import (
                PYCSAMT_STATION_RENDERING as SR,
            )

            ps = SR.pseudosection
            snap["station"] = {
                "side": ps.side,
                "show_markers": ps.show_markers,
                "marker_symbol": ps.marker.marker,
                "marker_size": ps.marker.size,
                "marker_offset": ps.marker.offset,
                "max_labels": ps.max_labels,
            }
        except Exception:
            pass

        # Section style (pseudosection preset axis)
        try:
            from pycsamt.api.section import (
                PYCSAMT_SECTION as SEC,
            )

            ax = SEC.pseudosection.axis
            snap["section"] = {
                "y_direction": ax.y_direction,
                "station_side": ax.station_side,
            }
        except Exception:
            pass

        # Topography
        try:
            from pycsamt.topo import PYCSAMT_TOPO as T

            snap["topography"] = {
                "enabled": T.enabled,
                "exaggeration": T.exaggeration,
                "marker_pad": T.marker_pad_fraction,
            }
        except Exception:
            pass

        # Existing visual pages must also participate in Cancel and profiles.
        try:
            from pycsamt.api.style import PYCSAMT_STYLE as S

            snap["style"] = {
                **{
                    f"{key}_color": getattr(S.mt, key).color
                    for key in ("xy", "yx", "xx", "yy", "te", "tm")
                },
                **{
                    f"{key}_lw": getattr(S.mt, key).lw
                    for key in ("xy", "yx", "xx", "yy", "te", "tm")
                },
                "correction_before": S.correction.before.color,
                "correction_after": S.correction.after.color,
            }
        except Exception:
            pass
        try:
            from pycsamt.api.interp import PYCSAMT_INTERP as I

            style = getattr(I, "default", I)
            section = getattr(
                style, "section", getattr(I, "pseudosection", None)
            )
            profile = getattr(style, "profile", getattr(I, "profile", None))
            snap["interpretation"] = {
                "section_cmap": getattr(
                    section, "cmap_K", getattr(section, "cmap", "viridis")
                ),
                "water_table_linestyle": getattr(
                    section, "wt_ls", getattr(section, "wt_linestyle", "--")
                ),
                "section_alpha": getattr(
                    section, "station_alpha", getattr(section, "alpha", 0.45)
                ),
                "profile_cmap": getattr(
                    section, "cmap_Sw", getattr(profile, "cmap", "viridis")
                ),
            }
        except Exception:
            pass

        return snap

    def restore(self, snap: dict[str, Any]) -> None:
        """Restore all singletons from a snapshot produced by :meth:`snapshot`."""
        if "view_controls" in snap:
            self.apply_view_controls(**snap["view_controls"])
        if "station" in snap:
            self.apply_station(**snap["station"])
        if "section" in snap:
            self.apply_section(**snap["section"])
        if "topography" in snap:
            self.apply_topography(**snap["topography"])
        for key in (
            "plot", "contour", "mesh", "ordering", "style", "interpretation",
            "pipe",
        ):
            if key in snap:
                getattr(self, f"apply_{key}")(**snap[key])

    # ── Per-tab apply ─────────────────────────────────────────────────────────

    def apply_view_controls(self, **kw: Any) -> None:
        """Write view-control fields to :data:`PYCSAMT_CONTROL`."""
        try:
            from pycsamt.api.control import (
                PYCSAMT_CONTROL as C,
            )

            for field in ("background", "toolbar", "fit_layout"):
                if "panel_" + field in kw:
                    setattr(C.panel, field, kw["panel_" + field])
            if "rho_view" in kw:
                C.rho.view = kw["rho_view"]
            if "phase_range" in kw:
                C.phase.range = tuple(kw["phase_range"])
            if "phase_unit" in kw:
                C.phase.unit = kw["phase_unit"]
            if "phase_wrap" in kw:
                C.phase.wrap = bool(kw["phase_wrap"])
            if "x_view" in kw:
                C.x.view = kw["x_view"]
        except Exception:
            pass

    def apply_station(self, **kw: Any) -> None:
        """Write station-marker fields to the pseudosection preset of
        :data:`PYCSAMT_STATION_RENDERING`."""
        try:
            from pycsamt.api.station import (
                PYCSAMT_STATION_RENDERING as SR,
            )

            ps = SR.pseudosection
            if "side" in kw:
                ps.side = kw["side"]
            if "show_markers" in kw:
                ps.show_markers = bool(kw["show_markers"])
            if "marker_symbol" in kw:
                ps.marker.marker = kw["marker_symbol"]
            if "marker_size" in kw:
                ps.marker.size = float(kw["marker_size"])
            if "marker_offset" in kw:
                ps.marker.offset = float(kw["marker_offset"])
            if "max_labels" in kw:
                ps.max_labels = int(kw["max_labels"])
        except Exception:
            pass

    def apply_section(self, **kw: Any) -> None:
        """Write section-axis fields to all pseudosection-family presets of
        :data:`PYCSAMT_SECTION`."""
        try:
            from pycsamt.api.section import (
                PYCSAMT_SECTION as SEC,
            )

            for preset in (
                "pseudosection",
                "dashboard",
                "compact",
                "publication",
                "dynamic",
            ):
                try:
                    ax = SEC.style_for(preset).axis
                    if "y_direction" in kw:
                        ax.y_direction = kw["y_direction"]
                    if "station_side" in kw:
                        ax.station_side = kw["station_side"]
                except Exception:
                    pass
        except Exception:
            pass

    def apply_topography(self, **kw: Any) -> None:
        """Write fields to :data:`PYCSAMT_TOPO`."""
        try:
            from pycsamt.topo import PYCSAMT_TOPO as T

            if "enabled" in kw:
                T.enabled = bool(kw["enabled"])
            if "exaggeration" in kw:
                T.exaggeration = float(kw["exaggeration"])
            if "marker_pad" in kw:
                T.marker_pad_fraction = float(kw["marker_pad"])
        except Exception:
            pass

    def apply_plot(self, **kw: Any) -> None:
        """Write package-wide plot and export defaults."""
        try:
            from pycsamt.api.plot import PLOT_CONFIG

            with PLOT_CONFIG.context(**kw):
                PLOT_CONFIG.resolve_formats()
            PLOT_CONFIG.configure(**kw)
        except Exception:
            pass

    def apply_contour(self, **kw: Any) -> None:
        """Write the live default contour style."""
        try:
            from pycsamt.api.contour import PYCSAMT_CONTOUR

            PYCSAMT_CONTOUR.configure(**kw)
        except Exception:
            pass

    def apply_mesh(self, **kw: Any) -> None:
        """Write shared mesh rendering fields."""
        try:
            from pycsamt.api.mesh import PYCSAMT_MESH

            PYCSAMT_MESH.configure(**kw)
        except Exception:
            pass

    def apply_ordering(self, **kw: Any) -> None:
        """Write package-wide site-ordering policy."""
        try:
            from pycsamt.api.ordering import PYCSAMT_ORDERING

            PYCSAMT_ORDERING.configure(**kw)
        except Exception:
            pass

    def apply_pipe(self, **kw: Any) -> None:
        """Write pipeline-engine settings to :data:`PYCSAMT_PIPE`."""
        try:
            from pycsamt.api.pipe import PYCSAMT_PIPE

            if "report_formats" in kw:  # JSON round-trips tuples as lists
                kw["report_formats"] = tuple(kw["report_formats"])
            PYCSAMT_PIPE.configure(**kw)
        except Exception:
            pass

    def apply_style(self, **kw: Any) -> None:
        """Write MT component and correction styles used by desktop plots."""
        try:
            from pycsamt.api.style import PYCSAMT_STYLE as S

            for key in ("xy", "yx", "xx", "yy", "te", "tm"):
                comp = getattr(S.mt, key)
                if f"{key}_color" in kw:
                    comp.color = kw[f"{key}_color"]
                if f"{key}_lw" in kw:
                    comp.lw = float(kw[f"{key}_lw"])
            if "correction_before" in kw:
                S.correction.before.color = kw["correction_before"]
            if "correction_after" in kw:
                S.correction.after.color = kw["correction_after"]
        except Exception:
            pass

    def apply_interpretation(self, **kw: Any) -> None:
        """Write the interpretation fields exposed by the desktop dialog."""
        try:
            from pycsamt.api.interp import PYCSAMT_INTERP as I

            style = getattr(I, "default", I)
            section = getattr(
                style, "section", getattr(I, "pseudosection", None)
            )
            profile = getattr(style, "profile", getattr(I, "profile", None))
            if "section_cmap" in kw:
                setattr(
                    section,
                    "cmap_K" if hasattr(section, "cmap_K") else "cmap",
                    kw["section_cmap"],
                )
            if "water_table_linestyle" in kw:
                setattr(
                    section,
                    "wt_ls" if hasattr(section, "wt_ls") else "wt_linestyle",
                    kw["water_table_linestyle"],
                )
            if "section_alpha" in kw:
                setattr(
                    section,
                    (
                        "station_alpha"
                        if hasattr(section, "station_alpha")
                        else "alpha"
                    ),
                    float(kw["section_alpha"]),
                )
            if "profile_cmap" in kw:
                if hasattr(section, "cmap_Sw"):
                    section.cmap_Sw = kw["profile_cmap"]
                elif hasattr(profile, "cmap"):
                    profile.cmap = kw["profile_cmap"]
        except Exception:
            pass

    # ── Reset ─────────────────────────────────────────────────────────────────

    def reset_tab(self, tab: str) -> None:
        """Reset the singleton(s) used by *tab* to package defaults."""
        try:
            if tab == "view_controls":
                from pycsamt.api.control import (
                    PYCSAMT_CONTROL,
                )

                PYCSAMT_CONTROL.reset()
            elif tab == "pseudosections":
                from pycsamt.api.section import (
                    PYCSAMT_SECTION,
                )
                from pycsamt.api.station import (
                    PYCSAMT_STATION_RENDERING,
                )

                PYCSAMT_STATION_RENDERING.reset()
                PYCSAMT_SECTION.reset()
            elif tab == "topography":
                from pycsamt.topo import PYCSAMT_TOPO

                PYCSAMT_TOPO.reset()
            elif tab == "display":
                from pycsamt.api.style import PYCSAMT_STYLE

                PYCSAMT_STYLE.reset()
            elif tab == "interpretation":
                from pycsamt.api.interp import PYCSAMT_INTERP

                PYCSAMT_INTERP.reset()
            elif tab == "output":
                from pycsamt.api.plot import PLOT_CONFIG

                PLOT_CONFIG.reset()
            elif tab == "rendering":
                from pycsamt.api.contour import PYCSAMT_CONTOUR
                from pycsamt.api.mesh import PYCSAMT_MESH

                PYCSAMT_CONTOUR.reset()
                PYCSAMT_MESH.reset()
            elif tab == "ordering":
                from pycsamt.api.ordering import PYCSAMT_ORDERING

                PYCSAMT_ORDERING.reset()
            elif tab == "pipeline":
                from pycsamt.api.pipe import PYCSAMT_PIPE

                PYCSAMT_PIPE.reset()
        except Exception:
            pass

    def reset_all(self) -> None:
        """Reset all API singletons to package defaults."""
        for tab in (
            "view_controls",
            "pseudosections",
            "topography",
            "display",
            "interpretation",
            "output",
            "rendering",
            "ordering",
            "pipeline",
        ):
            self.reset_tab(tab)

    # ── Persistence ───────────────────────────────────────────────────────────

    def save(self, path: Path | str | None = None) -> None:
        """Save current settings to *path* (default: ``~/.pycsamt/settings.json``)."""
        path = Path(path or self.SETTINGS_PATH)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(self.snapshot(), indent=2))

    def load(self, path: Path | str | None = None) -> bool:
        """Load settings from *path*.  Returns True if successful."""
        path = Path(path or self.SETTINGS_PATH)
        if not path.exists():
            return False
        try:
            snap = json.loads(path.read_text())
            self.restore(snap)
            return True
        except Exception:
            return False
