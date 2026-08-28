# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Compatibility adapters between PCBH and legacy geology boreholes.

The legacy :class:`pycsamt.geology.borehole.Borehole` is a profile-based,
vertical log view. It cannot represent a PCBH CRS, Y coordinate, deviated
trajectory, multiple log families, structures, or most interval metadata.
Consequently, downgrade requires an explicit profile coordinate and never
mutates the source PCBH object.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Mapping, Sequence
from typing import Union

from ...geology.borehole import Borehole, Interval
from ...geology.lithology import RockDatabase
from .schema import (
    Collar,
    LogInterval,
    PCBHBorehole,
    PCBHDocument,
    Trajectory,
    VocabularyEntry,
)

__all__ = [
    "from_legacy_borehole",
    "to_legacy_borehole",
    "legacy_borehole_views",
    "legacy_conversion_losses",
]

ProfileMapping = Union[Mapping[str, float], Callable[[PCBHBorehole], float]]


def from_legacy_borehole(
    borehole: Borehole,
    *,
    collar: Collar,
    borehole_id: str | None = None,
    total_depth_md: float | None = None,
    kind: str = "unknown",
    status: str = "unknown",
) -> PCBHBorehole:
    """Promote one legacy geology borehole to a vertical PCBH borehole.

    Parameters
    ----------
    borehole : Borehole
        Existing profile-based borehole. Depth and resistivity are assumed to
        be metres and ohm metres, matching the legacy API.
    collar : Collar
        Explicit absolute XYZ collar in the target PCBH CRS. The legacy
        ``borehole.x`` is a profile coordinate and is not used as easting.
    borehole_id : str, optional
        Stable PCBH identifier. Defaults to the legacy name.
    total_depth_md : float, optional
        Total measured depth in metres. Defaults to the deepest interval.
        It is required when the legacy borehole has no intervals.
    kind, status : str, optional
        PCBH borehole classification values.

    Returns
    -------
    PCBHBorehole
        A validated vertical PCBH borehole with a ``lithology`` log.

    Raises
    ------
    TypeError
        If the source or collar has the wrong type.
    ValueError
        If depth or resistivity values cannot form a valid PCBH object.

    Notes
    -----
    The legacy profile coordinate is retained in metadata as
    ``legacy_profile_x``. Empty legacy lithology names become ``unknown``
    because PCBH intervals require a code or label.
    """
    if not isinstance(borehole, Borehole):
        raise TypeError("borehole must be a geology.Borehole")
    if not isinstance(collar, Collar):
        raise TypeError("collar must be a PCBH Collar")
    collar.validate()

    intervals = [_legacy_interval_to_pcbh(item) for item in borehole.intervals]
    deepest = max((item.to_md for item in intervals), default=0.0)
    total_depth = (
        deepest
        if total_depth_md is None
        else _finite_float(total_depth_md, "total_depth_md")
    )
    if total_depth <= 0:
        raise ValueError(
            "total_depth_md must be provided and greater than zero when "
            "the legacy borehole has no intervals"
        )
    if total_depth < deepest:
        raise ValueError(
            "total_depth_md cannot be shallower than the deepest interval"
        )

    result = PCBHBorehole(
        id=borehole_id or borehole.name,
        name=borehole.name,
        kind=kind,
        status=status,
        collar=collar,
        total_depth_md=total_depth,
        trajectory=Trajectory(method="vertical"),
        interval_logs={"lithology": intervals} if intervals else {},
        metadata={"legacy_profile_x": float(borehole.x)},
    )
    result.validate()
    return result


def to_legacy_borehole(
    borehole: PCBHBorehole,
    *,
    profile_x: float,
    vocabulary: Sequence[VocabularyEntry] = (),
    rock_db: RockDatabase | None = None,
    family: str = "lithology",
) -> Borehole:
    """Create a legacy vertical-log view of a PCBH borehole.

    Parameters
    ----------
    borehole : PCBHBorehole
        Source PCBH borehole.
    profile_x : float
        Explicit distance along the target 2-D profile, in metres. PCBH collar
        easting is never silently treated as profile distance.
    vocabulary : sequence of VocabularyEntry, optional
        Embedded dictionary for resolving interval codes to names.
    rock_db : RockDatabase, optional
        If supplied, unnamed intervals with positive resistivity are
        classified through this database. It is not used by default.
    family : str, default='lithology'
        PCBH interval-log family to expose through the single legacy log.

    Returns
    -------
    Borehole
        A new legacy object suitable for existing calibrators and plots.

    Raises
    ------
    TypeError
        If public object parameters have incompatible types.
    ValueError
        If ``profile_x`` is non-finite or vocabulary codes conflict.

    Notes
    -----
    See :func:`legacy_conversion_losses` for fields unavailable in the
    returned view. Interval bounds and resistivity remain unchanged.
    """
    if not isinstance(borehole, PCBHBorehole):
        raise TypeError("borehole must be a PCBHBorehole")
    if rock_db is not None and not isinstance(rock_db, RockDatabase):
        raise TypeError("rock_db must be a RockDatabase or None")
    borehole.validate()
    x = _finite_float(profile_x, "profile_x")
    names = _vocabulary_names(vocabulary)

    intervals = [
        Interval(
            top=float(item.from_md),
            bottom=float(item.to_md),
            lithology=_lithology_name(item, names, rock_db),
            resistivity=(
                float(item.resistivity_ohm_m)
                if item.resistivity_ohm_m is not None
                else None
            ),
        )
        for item in borehole.interval_logs.get(family, [])
    ]
    return Borehole(
        name=borehole.name,
        x=x,
        intervals=intervals,
        collar_elevation=float(borehole.collar.z),
    )


def legacy_borehole_views(
    document: PCBHDocument,
    *,
    profile_x: ProfileMapping,
    family: str = "lithology",
    rock_db: RockDatabase | None = None,
) -> list[Borehole]:
    """Build legacy views for existing calibration workflows.

    Parameters
    ----------
    document : PCBHDocument
        Valid PCBH document using metres and ohm metres.
    profile_x : mapping or callable
        Either ``{borehole_id: profile_distance_m}`` or a callable receiving a
        PCBH borehole and returning its profile distance in metres.
    family : str, default='lithology'
        Interval family exposed to the calibrator.
    rock_db : RockDatabase, optional
        Explicit fallback classifier for unnamed resistivity intervals.

    Returns
    -------
    list of Borehole
        Independent legacy views accepted by ``ModelCalibrator.fit``.

    Raises
    ------
    TypeError
        If the document or profile mapping has an incompatible type.
    ValueError
        If document units are incompatible with the legacy API.
    KeyError
        If a mapping lacks a borehole identifier.
    """
    if not isinstance(document, PCBHDocument):
        raise TypeError("document must be a PCBHDocument")
    document.validate()
    if document.units.depth != "m" or document.units.resistivity != "ohm.m":
        raise ValueError(
            "legacy borehole views require depth='m' and "
            "resistivity='ohm.m'; convert document units first"
        )
    if not isinstance(profile_x, Mapping) and not callable(profile_x):
        raise TypeError("profile_x must be a mapping or callable")

    views = []
    vocabulary = (
        document.lithologies
        if family == "lithology"
        else document.formations
        if family == "formation"
        else ()
    )
    for borehole in document.boreholes:
        x = (
            profile_x(borehole)
            if callable(profile_x)
            else profile_x[borehole.id]
        )
        views.append(
            to_legacy_borehole(
                borehole,
                profile_x=x,
                vocabulary=vocabulary,
                rock_db=rock_db,
                family=family,
            )
        )
    return views


def legacy_conversion_losses(
    borehole: PCBHBorehole,
    *,
    family: str = "lithology",
) -> tuple[str, ...]:
    """Describe information omitted from a PCBH-to-legacy view.

    Parameters
    ----------
    borehole : PCBHBorehole
        Source object to inspect.
    family : str, default='lithology'
        The one interval family retained by the legacy view.

    Returns
    -------
    tuple of str
        Stable, human-readable loss descriptions.
    """
    if not isinstance(borehole, PCBHBorehole):
        raise TypeError("borehole must be a PCBHBorehole")
    borehole.validate()
    losses = [
        "absolute collar X/Y and CRS are replaced by one profile_x",
        "borehole id, kind, status, aliases, metadata, and extensions",
        "interval codes, descriptions, provenance, confidence, and properties",
    ]
    if borehole.trajectory.method != "vertical":
        losses.append("deviated trajectory and survey orientations")
    omitted = sorted(name for name in borehole.interval_logs if name != family)
    if omitted:
        losses.append("interval log families: " + ", ".join(omitted))
    if borehole.structures:
        losses.append("structural observations")
    if any(
        value is not None
        for value in (
            borehole.collar.longitude,
            borehole.collar.latitude,
            borehole.collar.position_uncertainty_m,
        )
    ):
        losses.append("collar longitude, latitude, and position uncertainty")
    if borehole.diameter is not None:
        losses.append("borehole diameter")
    selected = borehole.interval_logs.get(family, [])
    selected_depth = max((item.to_md for item in selected), default=0.0)
    if selected_depth != borehole.total_depth_md:
        losses.append("total depth outside the selected log coverage")
    return tuple(losses)


def _legacy_interval_to_pcbh(interval: Interval) -> LogInterval:
    if not isinstance(interval, Interval):
        raise TypeError("legacy borehole intervals must be Interval objects")
    resistivity = interval.resistivity
    if resistivity is not None:
        resistivity = _finite_float(resistivity, "interval resistivity")
        if resistivity <= 0:
            raise ValueError("interval resistivity must be greater than zero")
    return LogInterval(
        from_md=_finite_float(interval.top, "interval top"),
        to_md=_finite_float(interval.bottom, "interval bottom"),
        label=interval.lithology.strip() or "unknown",
        resistivity_ohm_m=resistivity,
        data_nature="observed",
    )


def _vocabulary_names(
    vocabulary: Sequence[VocabularyEntry],
) -> dict[str, str]:
    names: dict[str, str] = {}
    for entry in vocabulary:
        if not isinstance(entry, VocabularyEntry):
            raise TypeError(
                "vocabulary entries must be VocabularyEntry objects"
            )
        entry.validate()
        if entry.code in names:
            raise ValueError(f"duplicate vocabulary code {entry.code!r}")
        names[entry.code] = entry.name
    return names


def _lithology_name(
    interval: LogInterval,
    vocabulary: Mapping[str, str],
    rock_db: RockDatabase | None,
) -> str:
    candidate = vocabulary.get(interval.code or "") or interval.label
    if candidate:
        if rock_db is None:
            return candidate
        canonical = {
            entry.name.casefold(): entry.name for entry in rock_db.entries
        }
        return canonical.get(candidate.casefold(), candidate)
    if rock_db is not None and interval.resistivity_ohm_m is not None:
        return rock_db.classify(float(interval.resistivity_ohm_m)).name
    return interval.code or "unknown"


def _finite_float(value: float, name: str) -> float:
    if isinstance(value, bool):
        raise TypeError(f"{name} must be a finite number")
    try:
        number = float(value)
    except (TypeError, ValueError) as error:
        raise TypeError(f"{name} must be a finite number") from error
    if not math.isfinite(number):
        raise ValueError(f"{name} must be a finite number")
    return number
