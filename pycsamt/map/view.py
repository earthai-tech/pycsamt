# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
# ruff: noqa: E501
"""High-level :class:`MapView` session façade for pyCSAMT maps.

``MapView`` is the single object scripts, notebooks, and the map-view
platform share.  It holds one combined :class:`~pycsamt.map.MapData`
(possibly spanning many lines) and renders every view through the
existing builders, so the GUI carries no plotting logic of its own.

Examples
--------
>>> from pycsamt.map import MapView
>>> mv = MapView.from_folder("data/survey", detect="folder")
>>> mv.lines                                    # doctest: +SKIP
('L18PLT', 'L22PLT')
>>> fig = mv.station(overlay="rho", frequency=1024)   # doctest: +SKIP
>>> mv.export("out/station.html")                     # doctest: +SKIP
"""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from typing import Any

from ._core import MapData, ensure_map_data, station_dataframe
from .config import (
    ExportOptions,
    ProfileMapOptions,
    StationMapOptions,
    VolumeMapOptions,
)
from .export import export_figure
from .loader import load_lines
from .profile import build_profile_map, build_pseudosection
from .station import build_station_map
from .volume import build_3d_map

__all__ = ["MapView", "open_app"]

_VIEW_BUILDERS = {
    "station": "station",
    "profile": "profile",
    "pseudosection": "pseudosection",
    "map3d": "map3d",
}


class MapView:
    """A code-first session bound to one (multi-line) survey."""

    def __init__(
        self,
        data: Any,
        *,
        theme: str = "light",
        backend: str = "plotly",
    ) -> None:
        self.data: MapData = ensure_map_data(data)
        self.theme = theme
        self.backend = backend

    # ── constructors ───────────────────────────────────

    @classmethod
    def from_folder(
        cls,
        path: str | Path,
        *,
        detect: str = "folder",
        recursive: bool = True,
        theme: str = "light",
        backend: str = "plotly",
        **load_kwargs: Any,
    ) -> MapView:
        """Build a view from a directory of one or more lines."""
        data = load_lines(
            path,
            detect=detect,
            recursive=recursive,
            **load_kwargs,
        )
        return cls(data, theme=theme, backend=backend)

    @classmethod
    def from_lines(
        cls,
        lines: Any,
        *,
        theme: str = "light",
        backend: str = "plotly",
        **load_kwargs: Any,
    ) -> MapView:
        """Build a view from an explicit ``{line: source}`` mapping."""
        data = load_lines(lines, **load_kwargs)
        return cls(data, theme=theme, backend=backend)

    # ── survey introspection ───────────────────────────

    @property
    def lines(self) -> tuple[str, ...]:
        """Return the profile/line names in this survey."""
        return self.data.lines

    @property
    def stations(self) -> tuple[str, ...]:
        """Return the station IDs in this survey."""
        return self.data.station_ids

    @property
    def n_stations(self) -> int:
        """Return the number of stations across all lines."""
        return len(self.data.stations)

    @property
    def has_geo(self) -> bool:
        """Whether stations carry finite map coordinates."""
        return self.data.has_geo

    def table(self):
        """Return survey stations as a :class:`pandas.DataFrame`."""
        return station_dataframe(self.data)

    # ── topography / elevation ─────────────────────────

    def with_elevations(self, elev_map: dict[str, float]) -> MapView:
        """Return a copy with station elevations overridden by *elev_map*."""
        from .topo import apply_elevations

        new = object.__new__(MapView)
        new.data = apply_elevations(self.data, elev_map)
        new.theme = self.theme
        new.backend = self.backend
        return new

    def fetch_elevations(
        self,
        *,
        api_name: str = "open_meteo",
    ) -> dict[str, float]:
        """Fetch station elevations online (see :func:`pycsamt.map.fetch_elevations`)."""
        from .topo import fetch_elevations

        return fetch_elevations(self.data, api_name=api_name)

    # ── view renderers ─────────────────────────────────

    def station(
        self,
        *,
        options: StationMapOptions | None = None,
        **overrides: Any,
    ) -> Any:
        """Build a 2-D station map figure."""
        opts = self._station_options(options, overrides)
        return build_station_map(self.data, opts)

    def profile(
        self,
        *,
        options: ProfileMapOptions | None = None,
        **overrides: Any,
    ) -> Any:
        """Build a profile-view figure."""
        opts = self._profile_options(options, overrides)
        return build_profile_map(self.data, opts)

    def pseudosection(
        self,
        *,
        options: ProfileMapOptions | None = None,
        **overrides: Any,
    ) -> Any:
        """Build a resistivity/phase pseudosection figure."""
        opts = self._profile_options(options, overrides)
        return build_pseudosection(self.data, opts)

    def map3d(
        self,
        *,
        mode: str = "fence",
        options: VolumeMapOptions | None = None,
        **overrides: Any,
    ) -> Any:
        """Build a 3-D fence/block/depth/surface figure."""
        opts = self._volume_options(mode, options, overrides)
        return build_3d_map(self.data, opts)

    def figure(self, view: str = "station", **overrides: Any) -> Any:
        """Build a figure for *view* by name.

        Parameters
        ----------
        view : {"station", "profile", "pseudosection", "map3d"}
            Which renderer to call.
        **overrides :
            Forwarded to the matching renderer.
        """
        key = str(view).lower()
        if key not in _VIEW_BUILDERS:
            raise ValueError(
                f"Unknown view {view!r}. Expected one of "
                f"{sorted(_VIEW_BUILDERS)}."
            )
        return getattr(self, _VIEW_BUILDERS[key])(**overrides)

    # ── export ─────────────────────────────────────────

    def export(
        self,
        path: str | Path,
        *,
        view: str = "station",
        fmt: str | None = None,
        scale: float = 2.0,
        width: int | None = None,
        height: int | None = None,
        **overrides: Any,
    ) -> Path:
        """Render *view* and write it to *path*.

        The format is inferred from the suffix unless *fmt* is given;
        ``.html``, ``.png``, ``.json`` and Plotly static images are
        supported (see :func:`pycsamt.map.export_figure`).
        """
        fig = self.figure(view, **overrides)
        options = ExportOptions(
            path=path,
            format=fmt,
            width=width,
            height=height,
            scale=scale,
        )
        return export_figure(fig, options)

    def export_all(
        self,
        directory: str | Path,
        *,
        fmt: str = "html",
        views: tuple[str, ...] = (
            "station",
            "pseudosection",
            "map3d",
        ),
        **overrides: Any,
    ) -> dict[str, Path]:
        """Export several views into *directory*, keyed by view name."""
        out_dir = Path(directory)
        out_dir.mkdir(parents=True, exist_ok=True)
        written: dict[str, Path] = {}
        for view in views:
            target = out_dir / f"{view}.{fmt.lstrip('.')}"
            written[view] = self.export(
                target,
                view=view,
                fmt=fmt,
                **overrides,
            )
        return written

    # ── platform launch ────────────────────────────────

    def launch(self, **kwargs: Any) -> None:
        """Open this survey in the map-view platform.

        Requires the optional Dash GUI extra
        (``pycsamt.app.mapview``).
        """
        try:
            from pycsamt.app.mapview import launch as _launch
        except ImportError as exc:  # pragma: no cover - GUI extra
            raise ImportError(
                "The map-view platform requires the GUI "
                "dependencies (dash, dash-bootstrap-components). "
                "Install pycsamt with its app extras."
            ) from exc
        _launch(view=self, **kwargs)

    # ── option assembly ────────────────────────────────

    def _station_options(
        self,
        options: StationMapOptions | None,
        overrides: dict[str, Any],
    ) -> StationMapOptions:
        base = options or StationMapOptions(
            theme=self.theme,
            backend=self.backend,
        )
        return replace(base, **overrides) if overrides else base

    def _profile_options(
        self,
        options: ProfileMapOptions | None,
        overrides: dict[str, Any],
    ) -> ProfileMapOptions:
        base = options or ProfileMapOptions(
            theme=self.theme,
            backend=self.backend,
        )
        return replace(base, **overrides) if overrides else base

    def _volume_options(
        self,
        mode: str,
        options: VolumeMapOptions | None,
        overrides: dict[str, Any],
    ) -> VolumeMapOptions:
        base = options or VolumeMapOptions(
            mode=mode,  # type: ignore[arg-type]
            theme=self.theme,
        )
        if options is not None and "mode" not in overrides:
            overrides = {**overrides, "mode": mode}
        return replace(base, **overrides) if overrides else base

    def __repr__(self) -> str:
        return (
            f"MapView(lines={len(self.lines)}, "
            f"stations={self.n_stations}, "
            f"geo={self.has_geo})"
        )


def open_app(source: Any = None, **kwargs: Any) -> None:
    """Launch the map-view platform, optionally seeded with *source*.

    Parameters
    ----------
    source :
        Any input accepted by :meth:`MapView.from_folder` /
        :func:`pycsamt.map.load_lines`.  When ``None`` the platform
        starts empty and the user loads data in the GUI.
    **kwargs :
        Forwarded to the platform ``launch`` (e.g. ``host``, ``port``,
        ``open_browser``).
    """
    if source is None:
        try:
            from pycsamt.app.mapview import launch as _launch
        except ImportError as exc:  # pragma: no cover - GUI extra
            raise ImportError(
                "The map-view platform requires the GUI "
                "dependencies (dash, dash-bootstrap-components)."
            ) from exc
        _launch(**kwargs)
        return
    if isinstance(source, (str, Path)):
        view = MapView.from_folder(source)
    else:
        view = MapView(load_lines(source))
    view.launch(**kwargs)
