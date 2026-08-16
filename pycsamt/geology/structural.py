# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Structural geology — field measurements and fault traces.

Independent structural evidence to check an interpretation against, the
same role :class:`~pycsamt.geology.borehole.Borehole` plays for lithology.
None of this is electromagnetic: positions are profile-relative metres
(``x``), exactly like :class:`~pycsamt.geology.borehole.Borehole` and
:class:`~pycsamt.geology.lithology.StratigraphicLog` -- real-world
placement (lat/lon) is handled separately by :mod:`pycsamt.site` and
:mod:`pycsamt.gis`.

Two kinds of field measurement are distinguished, matching how a
stereonet actually plots them:

* :class:`StructuralMeasurement` -- a **planar** feature (bedding,
  foliation, joint, cleavage, contact, fault plane, unconformity
  surface), recorded as strike/dip/dip-direction.
* :class:`LinearMeasurement` -- a **linear** feature (fold axis,
  lineation, slickenline), recorded as trend/plunge.

:class:`FaultTrace` is a distinct, coarser entity: where a fault crosses
the 2-D profile itself, which side is downthrown, and (if known) the
throw -- the piece that plugs directly into the structural-continuity
review questions in :doc:`/user_guide/interpretation/workflow`
("do apparent boundary offsets align with known structures?").
:class:`StructuralModel` collects all three per profile.

Angle conventions
------------------
This is the one place the rest of :mod:`pycsamt` does not already agree
with itself: MT geoelectric strike (:mod:`pycsamt.emtools.strike`) is
axial and wrapped to ``(-90, 90]``; :mod:`pycsamt.gis`/:mod:`pycsamt.site`
azimuth is a directed compass bearing, ``[0, 360)``; the synthetic
geometry in :mod:`pycsamt.ai.geology.lenses` is axial mod-180 for an
unrelated reason (ellipse symmetry). None of those fit real field data,
so this module uses the convention modern digital field-mapping tools
(FieldMove, StraboSpot) record directly from a compass-clinometer:

* ``strike_deg``, ``trend_deg``, ``dip_direction_deg`` -- compass
  bearings, ``[0, 360)``, degrees clockwise from north.
* ``dip_deg``, ``plunge_deg`` -- angle below horizontal, ``[0, 90]``.

Storing ``dip_direction_deg`` alongside ``strike_deg`` (rather than only
a right-hand-rule pair) is deliberate: the two are cross-validated against
each other in :meth:`StructuralMeasurement.validate` (dip direction must
be within tolerance of strike +/- 90 degrees), which catches a transposed
field-notebook entry that a bare right-hand-rule number would not.
:meth:`StructuralMeasurement.from_right_hand_rule` builds one from just
dip direction and dip, for anyone who prefers recording that way.
"""

from __future__ import annotations

import csv
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Union

from ..api.property import MetadataMixin, PyCSAMTObject

__all__ = [
    "StructuralMeasurement",
    "LinearMeasurement",
    "FaultTrace",
    "StructuralModel",
]

PathLike = Union[str, Path]

_DIP_DIRECTION_TOLERANCE_DEG = 20.0


def _norm360(angle_deg: float) -> float:
    """Normalise an angle to ``[0, 360)``."""
    return float(angle_deg) % 360.0


def _angular_diff(a_deg: float, b_deg: float) -> float:
    """Smallest absolute difference between two compass bearings, in [0, 180]."""
    d = abs(_norm360(a_deg) - _norm360(b_deg)) % 360.0
    return d if d <= 180.0 else 360.0 - d


# ---------------------------------------------------------------------------
# StructuralMeasurement (planar)
# ---------------------------------------------------------------------------


@dataclass
class StructuralMeasurement(PyCSAMTObject):
    """A planar structural field measurement.

    Parameters
    ----------
    x : float
        Position along the survey profile, metres.
    kind : str
        Feature type, e.g. ``'bedding'``, ``'foliation'``, ``'joint'``,
        ``'cleavage'``, ``'contact'``, ``'fault_plane'``,
        ``'unconformity'``. Free text -- not an enforced enumeration,
        matching :attr:`~pycsamt.geology.borehole.Interval.lithology`.
    strike_deg : float
        Compass strike, degrees clockwise from north, ``[0, 360)`` as
        measured. Normalised on construction; not reduced modulo 180, so
        the raw field reading is preserved.
    dip_deg : float
        Dip angle below horizontal, degrees, ``[0, 90]``.
    dip_direction_deg : float
        Compass direction the surface dips toward, ``[0, 360)``. Must be
        within *dip_direction_tolerance_deg* of ``strike_deg + 90`` or
        ``strike_deg - 90`` (mod 360); raises :class:`ValueError`
        otherwise, since a wider mismatch usually means one of the two
        readings was transposed in the field notebook.
    z : float, optional
        Depth (positive downward) or elevation of the measurement,
        metres. ``None`` for a surface outcrop reading with no
        associated depth.
    station : str, optional
        Field station or outcrop label.
    confidence : float
        Subjective reading confidence, ``[0, 1]`` (default 1.0).
    notes : str
        Free-text field note.

    Examples
    --------
    >>> m = StructuralMeasurement(
    ...     x=500.0, kind="bedding", strike_deg=45.0, dip_deg=30.0,
    ...     dip_direction_deg=135.0,
    ... )
    >>> m.dip_azimuth_ok
    True
    """

    x: float
    kind: str
    strike_deg: float
    dip_deg: float
    dip_direction_deg: float
    z: float | None = None
    station: str | None = None
    confidence: float = 1.0
    notes: str = ""
    dip_direction_tolerance_deg: float = _DIP_DIRECTION_TOLERANCE_DEG

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        """Re-check and re-normalise this measurement's fields.

        Called automatically by ``__post_init__``, and by
        :meth:`~pycsamt.api.property.PyCSAMTObject.update`/:meth:`~pycsamt.api.property.PyCSAMTObject.clone`
        after they set new attribute values -- both go through this
        method rather than ``__post_init__`` (which only runs once, at
        construction), so a ``clone(dip_direction_deg=...)`` that breaks
        the strike/dip-direction consistency check is caught rather than
        silently accepted.
        """
        if not (0.0 <= self.dip_deg <= 90.0):
            raise ValueError(f"dip_deg ({self.dip_deg}) must be in [0, 90].")
        if not (0.0 <= self.confidence <= 1.0):
            raise ValueError(
                f"confidence ({self.confidence}) must be in [0, 1]."
            )
        self.strike_deg = _norm360(self.strike_deg)
        self.dip_direction_deg = _norm360(self.dip_direction_deg)
        if not self.dip_azimuth_ok:
            raise ValueError(
                f"dip_direction_deg ({self.dip_direction_deg}) is not "
                f"within {self.dip_direction_tolerance_deg} deg of "
                f"strike_deg ({self.strike_deg}) +/- 90 -- check for a "
                "transposed strike/dip-direction reading."
            )

    @property
    def dip_azimuth_ok(self) -> bool:
        """Whether ``dip_direction_deg`` is consistent with ``strike_deg``."""
        tol = self.dip_direction_tolerance_deg
        return (
            _angular_diff(self.dip_direction_deg, self.strike_deg + 90.0)
            <= tol
            or _angular_diff(self.dip_direction_deg, self.strike_deg - 90.0)
            <= tol
        )

    @classmethod
    def from_right_hand_rule(
        cls,
        x: float,
        kind: str,
        dip_direction_deg: float,
        dip_deg: float,
        **kwargs: Any,
    ) -> StructuralMeasurement:
        """Build from a dip-direction/dip pair, deriving strike.

        Strike is set to ``dip_direction_deg - 90`` (mod 360), the
        right-hand-rule convention: facing along strike with the dip
        direction to your right.
        """
        strike_deg = _norm360(float(dip_direction_deg) - 90.0)
        return cls(
            x=x,
            kind=kind,
            strike_deg=strike_deg,
            dip_deg=dip_deg,
            dip_direction_deg=dip_direction_deg,
            **kwargs,
        )

    def __repr__(self) -> str:
        return (
            f"StructuralMeasurement(x={self.x:.1f} m, {self.kind!r}, "
            f"{self.strike_deg:.0f}/{self.dip_deg:.0f}"
            f"->{self.dip_direction_deg:.0f})"
        )


# ---------------------------------------------------------------------------
# LinearMeasurement
# ---------------------------------------------------------------------------


@dataclass
class LinearMeasurement(PyCSAMTObject):
    """A linear structural field measurement.

    Parameters
    ----------
    x : float
        Position along the survey profile, metres.
    kind : str
        Feature type, e.g. ``'fold_axis'``, ``'lineation'``,
        ``'slickenline'``, ``'fold_hinge'``, ``'intersection_lineation'``.
        Free text, as with :class:`StructuralMeasurement`.
    trend_deg : float
        Compass direction the line plunges toward, degrees clockwise
        from north, ``[0, 360)``.
    plunge_deg : float
        Angle below horizontal, degrees, ``[0, 90]``.
    z : float, optional
    station : str, optional
    confidence : float
    notes : str

    Examples
    --------
    >>> LinearMeasurement(x=500.0, kind="fold_axis", trend_deg=210.0, plunge_deg=15.0)
    LinearMeasurement(x=500.0 m, 'fold_axis', 210/15)
    """

    x: float
    kind: str
    trend_deg: float
    plunge_deg: float
    z: float | None = None
    station: str | None = None
    confidence: float = 1.0
    notes: str = ""

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        """Re-check and re-normalise this measurement's fields.

        Called by ``__post_init__`` and by ``update``/``clone``; see
        :meth:`StructuralMeasurement.validate`.
        """
        if not (0.0 <= self.plunge_deg <= 90.0):
            raise ValueError(
                f"plunge_deg ({self.plunge_deg}) must be in [0, 90]."
            )
        if not (0.0 <= self.confidence <= 1.0):
            raise ValueError(
                f"confidence ({self.confidence}) must be in [0, 1]."
            )
        self.trend_deg = _norm360(self.trend_deg)

    def __repr__(self) -> str:
        return (
            f"LinearMeasurement(x={self.x:.1f} m, {self.kind!r}, "
            f"{self.trend_deg:.0f}/{self.plunge_deg:.0f})"
        )


# ---------------------------------------------------------------------------
# FaultTrace
# ---------------------------------------------------------------------------

_SENSES = ("normal", "reverse", "strike_slip", "unknown")
_SIDES = ("left", "right")


@dataclass
class FaultTrace(PyCSAMTObject):
    """Where a fault crosses the 2-D profile.

    Parameters
    ----------
    x : float
        Profile position where the fault is picked, metres.
    dip_deg : float
        Apparent dip of the fault plane *in the section*, degrees,
        ``[0, 90]`` (0 = flat detachment, 90 = vertical). This is the
        angle a 2-D EM section can actually constrain; the true 3-D dip
        differs unless the profile happens to run perpendicular to
        strike. Pass *strike_deg* separately when the true attitude is
        independently known (surface mapping, borehole).
    downthrown_side : {'left', 'right'}
        Which side of ``x`` -- toward decreasing or increasing profile
        distance -- is downthrown.
    sense : {'normal', 'reverse', 'strike_slip', 'unknown'}
        Kinematic sense, where known (default ``'unknown'``).
    throw_m : float, optional
        Vertical displacement, metres (magnitude; direction is carried
        by *downthrown_side*). ``None`` when unknown or unmeasured.
    strike_deg : float, optional
        True compass strike, when independently known.
    z_top : float, optional
        Depth to the top of the picked trace, metres. ``None`` for a
        surface trace or when unconstrained.
    confidence : float
    evidence : str
        Free-text source, e.g. ``'resistivity offset'``, ``'borehole'``,
        ``'surface mapping'``.
    notes : str

    Examples
    --------
    >>> FaultTrace(x=500.0, dip_deg=70.0, downthrown_side="right", throw_m=12.0)
    FaultTrace(x=500.0 m, dip=70 deg, down=right, throw=12.0 m)
    """

    x: float
    dip_deg: float
    downthrown_side: str
    sense: str = "unknown"
    throw_m: float | None = None
    strike_deg: float | None = None
    z_top: float | None = None
    confidence: float = 1.0
    evidence: str = ""
    notes: str = ""

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        """Re-check and re-normalise this trace's fields.

        Called by ``__post_init__`` and by ``update``/``clone``; see
        :meth:`StructuralMeasurement.validate`.
        """
        if not (0.0 <= self.dip_deg <= 90.0):
            raise ValueError(f"dip_deg ({self.dip_deg}) must be in [0, 90].")
        if self.downthrown_side not in _SIDES:
            raise ValueError(
                f"downthrown_side must be one of {_SIDES}, "
                f"got {self.downthrown_side!r}."
            )
        if self.sense not in _SENSES:
            raise ValueError(
                f"sense must be one of {_SENSES}, got {self.sense!r}."
            )
        if self.throw_m is not None and self.throw_m < 0.0:
            raise ValueError(
                f"throw_m ({self.throw_m}) must be >= 0; direction is "
                "carried by downthrown_side, not the sign of throw_m."
            )
        if not (0.0 <= self.confidence <= 1.0):
            raise ValueError(
                f"confidence ({self.confidence}) must be in [0, 1]."
            )
        if self.strike_deg is not None:
            self.strike_deg = _norm360(self.strike_deg)

    def __repr__(self) -> str:
        throw = f"{self.throw_m:.1f} m" if self.throw_m is not None else "?"
        return (
            f"FaultTrace(x={self.x:.1f} m, dip={self.dip_deg:.0f} deg, "
            f"down={self.downthrown_side}, throw={throw})"
        )


# ---------------------------------------------------------------------------
# StructuralModel
# ---------------------------------------------------------------------------


class StructuralModel(PyCSAMTObject, MetadataMixin):
    """Collection of structural evidence along one survey profile.

    Parameters
    ----------
    planar : list of StructuralMeasurement, optional
    linear : list of LinearMeasurement, optional
    faults : list of FaultTrace, optional
    metadata : dict, optional
        Free-form provenance, e.g. survey name or source file paths.

    Examples
    --------
    >>> model = StructuralModel(
    ...     faults=[FaultTrace(x=500.0, dip_deg=70.0, downthrown_side="right")],
    ... )
    >>> len(model.faults)
    1
    """

    def __init__(
        self,
        *,
        planar: Sequence[StructuralMeasurement] | None = None,
        linear: Sequence[LinearMeasurement] | None = None,
        faults: Sequence[FaultTrace] | None = None,
        metadata: dict | None = None,
    ) -> None:
        self.planar: list[StructuralMeasurement] = list(planar or [])
        self.linear: list[LinearMeasurement] = list(linear or [])
        self.faults: list[FaultTrace] = list(faults or [])
        self.metadata: dict = dict(metadata) if metadata else {}

    # ------------------------------------------------------------------
    # Mutation
    # ------------------------------------------------------------------

    def add_planar(self, measurement: StructuralMeasurement) -> None:
        self.planar.append(measurement)

    def add_linear(self, measurement: LinearMeasurement) -> None:
        self.linear.append(measurement)

    def add_fault(self, fault: FaultTrace) -> None:
        self.faults.append(fault)

    # ------------------------------------------------------------------
    # Queries
    # ------------------------------------------------------------------

    def within(self, x_min: float, x_max: float) -> StructuralModel:
        """Return a new model restricted to ``x_min <= x <= x_max``."""
        return StructuralModel(
            planar=[m for m in self.planar if x_min <= m.x <= x_max],
            linear=[m for m in self.linear if x_min <= m.x <= x_max],
            faults=[f for f in self.faults if x_min <= f.x <= x_max],
            metadata=dict(self.metadata),
        )

    def nearest(
        self,
        x: float,
        *,
        kind: str = "faults",
        max_distance: float | None = None,
    ) -> StructuralMeasurement | LinearMeasurement | FaultTrace | None:
        """Return the item of *kind* nearest to profile position *x*.

        Parameters
        ----------
        x : float
        kind : {'faults', 'planar', 'linear'}
        max_distance : float, optional
            Return ``None`` if the nearest item is farther than this
            (metres), instead of returning a distant match silently.
        """
        items = getattr(self, kind, None)
        if items is None:
            raise ValueError(
                "kind must be one of 'faults', 'planar', 'linear', "
                f"got {kind!r}."
            )
        if not items:
            return None
        best = min(items, key=lambda item: abs(item.x - x))
        if max_distance is not None and abs(best.x - x) > max_distance:
            return None
        return best

    # ------------------------------------------------------------------
    # I/O
    # ------------------------------------------------------------------

    @classmethod
    def from_csv(
        cls,
        *,
        planar_path: PathLike | None = None,
        linear_path: PathLike | None = None,
        faults_path: PathLike | None = None,
        delimiter: str = ",",
    ) -> StructuralModel:
        """Load a model from up to three CSV files, one per evidence type.

        Expected columns (case-insensitive header, optional columns may
        be omitted):

        * *planar_path* -- ``x, kind, strike_deg, dip_deg,
          dip_direction_deg[, z, station, confidence, notes]``
        * *linear_path* -- ``x, kind, trend_deg, plunge_deg[, z, station,
          confidence, notes]``
        * *faults_path* -- ``x, dip_deg, downthrown_side[, sense,
          throw_m, strike_deg, z_top, confidence, evidence, notes]``

        Any path left as ``None`` yields an empty list for that
        evidence type.
        """
        planar = (
            _read_planar_csv(planar_path, delimiter)
            if planar_path is not None
            else []
        )
        linear = (
            _read_linear_csv(linear_path, delimiter)
            if linear_path is not None
            else []
        )
        faults = (
            _read_faults_csv(faults_path, delimiter)
            if faults_path is not None
            else []
        )
        return cls(planar=planar, linear=linear, faults=faults)

    def to_dict(self) -> dict:
        return {
            "planar": [vars(m).copy() for m in self.planar],
            "linear": [vars(m).copy() for m in self.linear],
            "faults": [vars(f).copy() for f in self.faults],
        }

    def __repr__(self) -> str:
        return (
            f"StructuralModel({len(self.planar)} planar, "
            f"{len(self.linear)} linear, {len(self.faults)} faults)"
        )


# ---------------------------------------------------------------------------
# CSV helpers
# ---------------------------------------------------------------------------


def _header_lookup(fieldnames: Sequence[str]) -> dict:
    return {h.lower(): h for h in fieldnames}


def _optional_float(row: dict, key: str | None) -> float | None:
    if key is None:
        return None
    raw = row.get(key, "")
    if raw is None or raw.strip() in ("", "nan", "None", "NA"):
        return None
    try:
        return float(raw)
    except ValueError:
        return None


def _optional_str(row: dict, key: str | None) -> str | None:
    if key is None:
        return None
    raw = row.get(key)
    return raw.strip() if raw is not None and raw.strip() else None


def _read_planar_csv(
    path: PathLike, delimiter: str
) -> list[StructuralMeasurement]:
    p = Path(path)
    out: list[StructuralMeasurement] = []
    with p.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter=delimiter)
        if reader.fieldnames is None:
            raise ValueError(f"CSV file has no header: {p}")
        h = _header_lookup(reader.fieldnames)
        for row in reader:
            kwargs: dict[str, Any] = {}
            z = _optional_float(row, h.get("z"))
            if z is not None:
                kwargs["z"] = z
            station = _optional_str(row, h.get("station"))
            if station is not None:
                kwargs["station"] = station
            conf = _optional_float(row, h.get("confidence"))
            if conf is not None:
                kwargs["confidence"] = conf
            notes = _optional_str(row, h.get("notes"))
            if notes is not None:
                kwargs["notes"] = notes
            out.append(
                StructuralMeasurement(
                    x=float(row[h["x"]]),
                    kind=row[h["kind"]].strip(),
                    strike_deg=float(row[h["strike_deg"]]),
                    dip_deg=float(row[h["dip_deg"]]),
                    dip_direction_deg=float(row[h["dip_direction_deg"]]),
                    **kwargs,
                )
            )
    return out


def _read_linear_csv(
    path: PathLike, delimiter: str
) -> list[LinearMeasurement]:
    p = Path(path)
    out: list[LinearMeasurement] = []
    with p.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter=delimiter)
        if reader.fieldnames is None:
            raise ValueError(f"CSV file has no header: {p}")
        h = _header_lookup(reader.fieldnames)
        for row in reader:
            kwargs: dict[str, Any] = {}
            z = _optional_float(row, h.get("z"))
            if z is not None:
                kwargs["z"] = z
            station = _optional_str(row, h.get("station"))
            if station is not None:
                kwargs["station"] = station
            conf = _optional_float(row, h.get("confidence"))
            if conf is not None:
                kwargs["confidence"] = conf
            notes = _optional_str(row, h.get("notes"))
            if notes is not None:
                kwargs["notes"] = notes
            out.append(
                LinearMeasurement(
                    x=float(row[h["x"]]),
                    kind=row[h["kind"]].strip(),
                    trend_deg=float(row[h["trend_deg"]]),
                    plunge_deg=float(row[h["plunge_deg"]]),
                    **kwargs,
                )
            )
    return out


def _read_faults_csv(path: PathLike, delimiter: str) -> list[FaultTrace]:
    p = Path(path)
    out: list[FaultTrace] = []
    with p.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter=delimiter)
        if reader.fieldnames is None:
            raise ValueError(f"CSV file has no header: {p}")
        h = _header_lookup(reader.fieldnames)
        for row in reader:
            kwargs: dict[str, Any] = {}
            sense = _optional_str(row, h.get("sense"))
            if sense is not None:
                kwargs["sense"] = sense
            throw_m = _optional_float(row, h.get("throw_m"))
            if throw_m is not None:
                kwargs["throw_m"] = throw_m
            strike_deg = _optional_float(row, h.get("strike_deg"))
            if strike_deg is not None:
                kwargs["strike_deg"] = strike_deg
            z_top = _optional_float(row, h.get("z_top"))
            if z_top is not None:
                kwargs["z_top"] = z_top
            conf = _optional_float(row, h.get("confidence"))
            if conf is not None:
                kwargs["confidence"] = conf
            evidence = _optional_str(row, h.get("evidence"))
            if evidence is not None:
                kwargs["evidence"] = evidence
            notes = _optional_str(row, h.get("notes"))
            if notes is not None:
                kwargs["notes"] = notes
            out.append(
                FaultTrace(
                    x=float(row[h["x"]]),
                    dip_deg=float(row[h["dip_deg"]]),
                    downthrown_side=row[h["downthrown_side"]].strip(),
                    **kwargs,
                )
            )
    return out
