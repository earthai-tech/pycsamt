"""Shared topography resolution for inversion-agent outputs."""

from __future__ import annotations

from typing import Any

import numpy as np


def resolve_agent_topography(
    config: Any,
    *,
    sites: Any,
    station_names: list[str],
    coords_m: np.ndarray | None,
    warnings_list: list[str],
) -> dict[str, Any] | None:
    """Resolve terrain and align it to stations retained by preprocessing.

    Terrain changes output geometry only.  It never modifies an inversion
    prediction or implies that its forward solver used a topographic mesh.
    """
    if config is None or config is False:
        return None
    if config is True:
        cfg: dict[str, Any] = {}
    elif isinstance(config, dict):
        cfg = dict(config)
        if not cfg.get("enabled", True):
            return None
    else:
        warnings_list.append("'topography' must be a bool or mapping; ignored.")
        return None

    exaggeration = float(cfg.get("exaggeration", 1.0))
    if not np.isfinite(exaggeration) or exaggeration <= 0:
        warnings_list.append(
            "Topography exaggeration must be finite and positive; ignored."
        )
        return None
    interp_method = str(cfg.get("interp_method", "linear")).lower()
    if interp_method not in {"linear", "nearest", "cubic"}:
        warnings_list.append(
            "Topography interp_method must be linear, nearest, or cubic; ignored."
        )
        return None

    explicit_elev = cfg.get("elevation_m")
    explicit_chain = cfg.get("chainage_km")
    source = "array" if explicit_elev is not None else "sites"
    if explicit_elev is None:
        from ..topo.extract import (
            extract_chainage,
            extract_elevation,
            extract_station_names,
        )

        all_names = extract_station_names(sites)
        all_elev = extract_elevation(sites)
        all_chain = extract_chainage(sites)
        lookup: dict[str, int] = {}
        for i, name in enumerate(all_names):
            lookup.setdefault(str(name).strip().casefold(), i)
        missing = [n for n in station_names if str(n).strip().casefold() not in lookup]
        if missing:
            warnings_list.append(
                "Topography station-name alignment failed for: "
                + ", ".join(missing[:5])
                + (" ..." if len(missing) > 5 else "")
            )
            return _metadata(
                applied=False,
                source=source,
                elevation=None,
                chainage=None,
                exaggeration=exaggeration,
                interp_method=interp_method,
            )
        idx = np.asarray([lookup[str(n).strip().casefold()] for n in station_names])
        elevation = np.asarray(all_elev, dtype=float)[idx]
        chainage = np.asarray(all_chain, dtype=float)[idx]
        chainage = chainage - chainage[0]
    else:
        elevation = np.asarray(explicit_elev, dtype=float).reshape(-1)
        if elevation.size != len(station_names):
            warnings_list.append(
                "Explicit elevation_m length must equal the usable station "
                "count; topography ignored."
            )
            return None
        if explicit_chain is None:
            if coords_m is None:
                warnings_list.append(
                    "Explicit elevation_m also requires chainage_km when "
                    "station coordinates are unavailable; topography ignored."
                )
                return None
            xy = np.asarray(coords_m, dtype=float)
            seg = np.sqrt((np.diff(xy, axis=0) ** 2).sum(axis=1))
            chainage = np.concatenate([[0.0], np.cumsum(seg)]) / 1000.0
        else:
            chainage = np.asarray(explicit_chain, dtype=float).reshape(-1)

    if chainage.size != elevation.size or not np.all(np.isfinite(chainage)):
        warnings_list.append("Topography chainage is invalid; topography ignored.")
        return None
    if not np.all(np.isfinite(elevation)) or not np.any(elevation != 0.0):
        warnings_list.append(
            "No finite, non-zero station elevations were available; "
            "topography was not applied."
        )
        return _metadata(
            applied=False,
            source=source,
            elevation=elevation,
            chainage=chainage,
            exaggeration=exaggeration,
            interp_method=interp_method,
        )
    if np.any(np.diff(chainage) <= 0):
        warnings_list.append(
            "Topography chainage must increase in station order; " "topography ignored."
        )
        return None
    return _metadata(
        applied=True,
        source=source,
        elevation=elevation,
        chainage=chainage,
        exaggeration=exaggeration,
        interp_method=interp_method,
        station_names=station_names,
    )


def _metadata(
    *,
    applied: bool,
    source: str,
    elevation: Any,
    chainage: Any,
    exaggeration: float,
    interp_method: str,
    station_names: list[str] | None = None,
) -> dict[str, Any]:
    data = {
        "requested": True,
        "applied": applied,
        "rendered": False,
        "source": source,
        "affects_forward_physics": False,
        "vertical_datum": "metres above sea level",
        "elevation_m": elevation,
        "chainage_km": chainage,
        "exaggeration": exaggeration,
        "interp_method": interp_method,
    }
    if station_names is not None:
        data["station_names"] = list(station_names)
    return data


__all__ = ["resolve_agent_topography"]
