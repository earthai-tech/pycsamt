# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Hydrogeophysical interpretation from EM resistivity models.

This module turns an EM-derived :class:`~pycsamt.interp.ResistivityModel`
into hydrogeological deliverables: aquifer targets, clay/saline warnings,
fractured/weathered zones, basement indicators, confidence maps, and
station summaries. It deliberately builds on the existing interpretation
objects instead of duplicating them.
"""

from __future__ import annotations

import csv
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Union

import numpy as np

from ..api.property import MetadataMixin, PyCSAMTObject
from ._base import ResistivityModel
from .borehole import Borehole
from .calibrate import ModelCalibrator
from .lithology import RockDatabase, StratigraphicLog

PathLike = Union[str, Path]

__all__ = [
    "AquiferZone",
    "HydroGeophysicalModel",
    "HydroInterpreter",
    "HydroUnit",
]


@dataclass
class HydroUnit(PyCSAMTObject, MetadataMixin):
    """One classified hydrogeophysical cell."""

    x: float
    z: float
    rho_ohm_m: float
    rho_log10: float
    unit: str
    lithology: str
    confidence: float = 1.0
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass
class AquiferZone(PyCSAMTObject, MetadataMixin):
    """Contiguous aquifer-favourable interval at one station/profile column."""

    station: str
    x: float
    top: float
    bottom: float
    mean_rho_ohm_m: float
    confidence: float
    zone_type: str = "aquifer"
    metadata: dict[str, Any] = field(default_factory=dict)

    @property
    def thickness(self) -> float:
        """Interval thickness in metres."""
        return float(self.bottom - self.top)


@dataclass
class HydroGeophysicalModel(PyCSAMTObject, MetadataMixin):
    """Hydrogeological interpretation product for one resistivity model."""

    resistivity_model: ResistivityModel
    unit_map: np.ndarray
    confidence: np.ndarray
    zones: list[AquiferZone] = field(default_factory=list)
    logs: list[StratigraphicLog] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)

    def aquifer_zones(self, *, min_confidence: float = 0.0) -> list[AquiferZone]:
        """Return aquifer-favourable zones above a confidence threshold."""
        return [
            zone
            for zone in self.zones
            if zone.zone_type == "aquifer" and zone.confidence >= min_confidence
        ]

    def station_summary(self) -> list[dict[str, Any]]:
        """Return compact per-station hydrogeological summaries."""
        rows: list[dict[str, Any]] = []
        model = self.resistivity_model
        for ix, x in enumerate(model.x_centers):
            name = _station_name(model, ix)
            col_units = self.unit_map[:, ix]
            col_conf = self.confidence[:, ix]
            aquifer_cells = col_units == "aquifer"
            clay_cells = np.isin(col_units, ["clay", "saline"])
            rows.append(
                {
                    "station": name,
                    "x_m": float(x),
                    "aquifer_cells": int(np.sum(aquifer_cells)),
                    "clay_or_saline_cells": int(np.sum(clay_cells)),
                    "mean_confidence": float(np.nanmean(col_conf)),
                    "n_zones": sum(1 for z in self.zones if z.station == name),
                }
            )
        return rows

    def to_csv(self, path: PathLike) -> Path:
        """Write cell-level hydro units and aquifer zones to CSV."""
        out = Path(path)
        out.parent.mkdir(parents=True, exist_ok=True)
        model = self.resistivity_model
        with out.open("w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(
                [
                    "station",
                    "x_m",
                    "z_m",
                    "rho_log10",
                    "rho_ohm_m",
                    "hydro_unit",
                    "confidence",
                ]
            )
            for ix, x in enumerate(model.x_centers):
                station = _station_name(model, ix)
                for iz, z in enumerate(model.z_centers):
                    rho_log = float(model.rho_2d[iz, ix])
                    writer.writerow(
                        [
                            station,
                            float(x),
                            float(z),
                            rho_log,
                            float(10.0**rho_log),
                            str(self.unit_map[iz, ix]),
                            float(self.confidence[iz, ix]),
                        ]
                    )
        return out

    def zones_to_csv(self, path: PathLike) -> Path:
        """Write interpreted aquifer/fracture target intervals to CSV."""
        out = Path(path)
        out.parent.mkdir(parents=True, exist_ok=True)
        with out.open("w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(
                [
                    "station",
                    "x_m",
                    "top_m",
                    "bottom_m",
                    "thickness_m",
                    "mean_rho_ohm_m",
                    "confidence",
                    "zone_type",
                ]
            )
            for zone in self.zones:
                writer.writerow(
                    [
                        zone.station,
                        zone.x,
                        zone.top,
                        zone.bottom,
                        zone.thickness,
                        zone.mean_rho_ohm_m,
                        zone.confidence,
                        zone.zone_type,
                    ]
                )
        return out


class HydroInterpreter(PyCSAMTObject):
    """Rule-based hydrogeophysical interpreter for EM resistivity sections.

    The default thresholds are conservative and intended as a first-pass
    screening model. Local calibration can be introduced by passing boreholes,
    a custom :class:`RockDatabase`, or explicit threshold overrides.
    """

    def __init__(
        self,
        *,
        context: str = "",
        db: RockDatabase | None = None,
        water_table_depth: float | None = None,
        aquifer_range: tuple[float, float] = (5.0, 300.0),
        fracture_range: tuple[float, float] = (50.0, 800.0),
        clay_max: float = 20.0,
        saline_max: float = 3.0,
        basement_min: float = 1000.0,
        min_zone_thickness: float = 1.0,
        calibration_ptol: float = 0.10,
        max_borehole_distance: float = 500.0,
    ) -> None:
        self.context = context
        self.db = db if db is not None else RockDatabase.default()
        self.water_table_depth = water_table_depth
        self.aquifer_range = tuple(float(v) for v in aquifer_range)
        self.fracture_range = tuple(float(v) for v in fracture_range)
        self.clay_max = float(clay_max)
        self.saline_max = float(saline_max)
        self.basement_min = float(basement_min)
        self.min_zone_thickness = float(min_zone_thickness)
        self.calibration_ptol = float(calibration_ptol)
        self.max_borehole_distance = float(max_borehole_distance)
        self._result: HydroGeophysicalModel | None = None

    def fit(
        self,
        model: Any,
        *,
        boreholes: Sequence[Borehole] | None = None,
    ) -> HydroGeophysicalModel:
        """Interpret a resistivity model or inversion result."""
        res_model = _coerce_resistivity_model(model)
        source_method = res_model.method
        if boreholes:
            calibrator = ModelCalibrator(
                ptol=self.calibration_ptol,
                max_borehole_distance=self.max_borehole_distance,
                db=self.db,
                verbose=False,
            ).fit(res_model, boreholes)
            res_model = calibrator.calibrated_model()

        unit_map = np.empty(res_model.rho_2d.shape, dtype=object)
        confidence = np.full(res_model.rho_2d.shape, np.nan, dtype=float)

        for ix in range(res_model.n_x):
            for iz in range(res_model.n_z):
                rho_log = float(res_model.rho_2d[iz, ix])
                z = float(res_model.z_centers[iz])
                unit, conf = self._classify_cell(10.0**rho_log, z)
                unit_map[iz, ix] = unit
                confidence[iz, ix] = conf

        logs = _logs_from_model(res_model, self.db)
        zones = self._extract_zones(res_model, unit_map, confidence)
        out = HydroGeophysicalModel(
            resistivity_model=res_model,
            unit_map=unit_map,
            confidence=confidence,
            zones=zones,
            logs=logs,
            metadata={
                "context": self.context,
                "source_method": source_method,
                "calibrated": bool(boreholes),
                "thresholds": {
                    "aquifer_range": self.aquifer_range,
                    "fracture_range": self.fracture_range,
                    "clay_max": self.clay_max,
                    "saline_max": self.saline_max,
                    "basement_min": self.basement_min,
                    "water_table_depth": self.water_table_depth,
                },
            },
        )
        self._result = out
        return out

    def aquifer_zones(self, *, min_confidence: float = 0.0) -> list[AquiferZone]:
        """Return aquifer zones from the last fitted model."""
        if self._result is None:
            raise RuntimeError("HydroInterpreter.fit must be called first.")
        return self._result.aquifer_zones(min_confidence=min_confidence)

    def _classify_cell(self, rho: float, z: float) -> tuple[str, float]:
        if not np.isfinite(rho) or rho <= 0:
            return "unknown", 0.0
        if rho <= self.saline_max:
            return "saline", _low_range_confidence(rho, self.saline_max)
        if rho <= self.clay_max:
            return "clay", _low_range_confidence(rho, self.clay_max)
        if self.water_table_depth is not None and z < self.water_table_depth:
            if rho < self.basement_min:
                return "vadose/weathered", 0.45
        aq_lo, aq_hi = self.aquifer_range
        if aq_lo <= rho <= aq_hi:
            return "aquifer", _range_confidence(rho, aq_lo, aq_hi)
        fr_lo, fr_hi = self.fracture_range
        if fr_lo <= rho <= fr_hi:
            return "fractured/weathered", 0.65 * _range_confidence(rho, fr_lo, fr_hi)
        if rho >= self.basement_min:
            return "resistive basement", min(
                1.0, np.log10(rho / self.basement_min + 1.0)
            )
        return "transition", 0.35

    def _extract_zones(
        self,
        model: ResistivityModel,
        unit_map: np.ndarray,
        confidence: np.ndarray,
    ) -> list[AquiferZone]:
        zones: list[AquiferZone] = []
        z_edges = _depth_edges(model.z_centers)
        for ix, x in enumerate(model.x_centers):
            station = _station_name(model, ix)
            mask = unit_map[:, ix] == "aquifer"
            start: int | None = None
            for iz, is_aquifer in enumerate(np.r_[mask, False]):
                if is_aquifer and start is None:
                    start = iz
                elif not is_aquifer and start is not None:
                    stop = iz
                    top = float(z_edges[start])
                    bottom = float(z_edges[stop])
                    if bottom - top >= self.min_zone_thickness:
                        rho_vals = 10.0 ** model.rho_2d[start:stop, ix]
                        zones.append(
                            AquiferZone(
                                station=station,
                                x=float(x),
                                top=top,
                                bottom=bottom,
                                mean_rho_ohm_m=float(np.nanmean(rho_vals)),
                                confidence=float(
                                    np.nanmean(confidence[start:stop, ix])
                                ),
                                zone_type="aquifer",
                            )
                        )
                    start = None
        return zones


def _coerce_resistivity_model(model: Any) -> ResistivityModel:
    if isinstance(model, ResistivityModel):
        return model
    if hasattr(model, "to_resistivity_model"):
        return model.to_resistivity_model()
    if isinstance(model, dict) and "rho_2d" in model:
        return ResistivityModel.from_array(
            model["rho_2d"],
            model.get("x_centers"),
            model.get("z_centers"),
            station_x=model.get("station_x"),
            station_names=model.get("station_names"),
            method=str(model.get("method", "generic")),
            rms=float(model.get("rms", float("nan"))),
        )
    raise TypeError(
        "model must be a ResistivityModel or expose to_resistivity_model()."
    )


def _logs_from_model(
    model: ResistivityModel, db: RockDatabase
) -> list[StratigraphicLog]:
    logs: list[StratigraphicLog] = []
    station_x = model.station_x if len(model.station_x) else model.x_centers
    names = model.station_names or [f"S{i:03d}" for i in range(len(station_x))]
    for ix, x in enumerate(station_x):
        col_idx = int(np.argmin(np.abs(model.x_centers - x)))
        logs.append(
            StratigraphicLog.from_column(
                station_name=names[ix] if ix < len(names) else f"S{ix:03d}",
                x=float(x),
                z_centers=model.z_centers,
                rho_log10=model.rho_2d[:, col_idx],
                db=db,
            )
        )
    return logs


def _station_name(model: ResistivityModel, ix: int) -> str:
    if model.station_names and ix < len(model.station_names):
        return str(model.station_names[ix])
    return f"S{ix:03d}"


def _depth_edges(z_centers: np.ndarray) -> np.ndarray:
    z = np.asarray(z_centers, dtype=float)
    if z.size == 1:
        dz = max(float(z[0]), 1.0)
        return np.array([max(0.0, z[0] - 0.5 * dz), z[0] + 0.5 * dz])
    mids = 0.5 * (z[:-1] + z[1:])
    first = max(0.0, z[0] - (mids[0] - z[0]))
    last = z[-1] + (z[-1] - mids[-1])
    return np.r_[first, mids, last]


def _range_confidence(value: float, low: float, high: float) -> float:
    if high <= low:
        return 0.5
    center = np.sqrt(low * high)
    half_width = max(np.log10(high) - np.log10(center), 1e-6)
    dist = abs(np.log10(value) - np.log10(center))
    return float(np.clip(1.0 - 0.5 * dist / half_width, 0.35, 1.0))


def _low_range_confidence(value: float, high: float) -> float:
    return float(np.clip(1.0 - 0.25 * value / max(high, 1e-12), 0.5, 1.0))
