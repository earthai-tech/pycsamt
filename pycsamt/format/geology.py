# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""PCGL - pyCSAMT Common Geology Legend.

PCGL is a small, portable document for the resistivity-to-geology
*interpretation legend* that Map View's Interpretation overlay applies to
a 3-D block / fence / depth-slice / iso-surface view: an ordered table of
named units, each with a resistivity range, a display colour, and an
optional pattern/swatch reference, plus a little document metadata so the
legend can travel between projects and applications on its own.

PCGL is deliberately a *thin wrapper*: the table itself is a
:class:`pycsamt.geology.lithology.RockDatabase` (already the
resistivity-to-rock classification engine used across pycsamt), extended
there with ``pattern_id``/``pattern_source`` fields. This module only adds
the document envelope (id, timestamps, provenance, versioning) and
canonical JSON (de)serialization on top of it -- it introduces no parallel
class hierarchy, the same way :mod:`pycsamt.format.borehole` bridges
:class:`pycsamt.geology.Borehole` rather than reimplementing it.

Canonical encoding is UTF-8 JSON, canonical extension ``.pcgl.json``. A
plain CSV with columns ``name, rho_min, rho_max, color, description, code,
source, pattern_id, pattern_source`` (only ``name, rho_min, rho_max``
required) round-trips through :meth:`GeologyLegend.from_csv` /
:meth:`GeologyLegend.to_csv` -- the same columns
:meth:`pycsamt.geology.lithology.RockDatabase.from_csv` accepts, so an
existing rock-property CSV already loads as a legend with no pattern
assigned.

Example
-------
>>> from pycsamt.format.geology import GeologyLegend
>>> legend = GeologyLegend.from_rock_database(title="Site A interpretation")
>>> legend.db.classify(250.0).name
'Granite (weathered)'
>>> legend.issues()
[]
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Union

from ..api.property import MetadataMixin, PyCSAMTObject
from ..geology.lithology import RockDatabase, RockEntry

__all__ = [
    "PCGL_VERSION",
    "GeologyLegend",
    "GeologyLegendValidationError",
    "read_legend",
    "write_legend",
    "legend_from_csv",
    "legend_to_dict",
    "legend_from_dict",
]

PCGL_VERSION = "0.1.0"

PathLike = Union[str, Path]


class GeologyLegendValidationError(ValueError):
    """Raised when a PCGL document fails validation."""


@dataclass(repr=False)
class GeologyLegend(PyCSAMTObject, MetadataMixin):
    """A versioned, portable resistivity-to-geology legend document.

    Parameters
    ----------
    document_id : str
    created_at : str
        ISO-8601 timestamp.
    created_by : str
        Creating application, workflow step, or person
        (e.g. ``"auto-suggest"``, ``"manual"``, a user name).
    db : RockDatabase
        The legend table itself -- resistivity ranges, colours, and
        pattern references. See :class:`pycsamt.geology.lithology.
        RockDatabase`.
    pcgl_version : str
    title, description : str
    rho_unit : str
        Resistivity unit the ranges are expressed in. PCGL 0.1 only
        defines ``"ohm.m"`` (linear ohm-metres); readers MUST reject any
        other value rather than silently reinterpreting the ranges.
    """

    document_id: str
    created_at: str
    created_by: str
    db: RockDatabase
    pcgl_version: str = PCGL_VERSION
    title: str = ""
    description: str = ""
    rho_unit: str = "ohm.m"
    metadata: dict[str, Any] = field(default_factory=dict)

    # -- convenience passthrough -----------------------------------
    @property
    def entries(self) -> tuple[RockEntry, ...]:
        return self.db.entries

    def classify(self, rho_ohm_m: float, method: str = "nearest") -> RockEntry:
        return self.db.classify(rho_ohm_m, method=method)

    # -- construction -------------------------------------------------
    @classmethod
    def from_rock_database(
        cls,
        db: RockDatabase | None = None,
        *,
        document_id: str | None = None,
        created_by: str = "pycsamt",
        title: str = "",
        description: str = "",
    ) -> GeologyLegend:
        """Wrap an existing (or the built-in default) rock database as a
        PCGL document, ready to write or hand to Map View."""
        legend = cls(
            document_id=document_id or f"pcgl:{_now()}",
            created_at=_now(),
            created_by=created_by,
            db=db if db is not None else RockDatabase.default(),
            title=title,
            description=description,
        )
        legend.validate()
        return legend

    @classmethod
    def from_csv(
        cls,
        path: PathLike,
        *,
        document_id: str | None = None,
        created_by: str = "pycsamt CSV importer",
        title: str = "",
    ) -> GeologyLegend:
        """Load a legend from a CSV (see module docstring for columns)."""
        p = Path(path)
        db = RockDatabase.from_csv(p)
        return cls.from_rock_database(
            db,
            document_id=document_id or f"pcgl:{p.stem}",
            created_by=created_by,
            title=title or p.stem,
        )

    def to_csv(self, path: PathLike) -> Path:
        """Write the legend table as CSV (document metadata is dropped;
        use :func:`write_legend` to keep it)."""
        return self.db.to_csv(path)

    # -- validation -----------------------------------------------------
    def issues(self) -> list[str]:
        out: list[str] = []
        if not str(self.document_id).strip():
            out.append("document_id must be non-empty")
        try:
            parts = tuple(int(p) for p in self.pcgl_version.split("."))
        except (AttributeError, ValueError):
            parts = ()
        if len(parts) != 3 or parts[0] != 0:
            out.append(f"pcgl_version must be compatible with {PCGL_VERSION}")
        if self.rho_unit != "ohm.m":
            out.append(f"unsupported rho_unit {self.rho_unit!r} (only 'ohm.m')")
        if not isinstance(self.db, RockDatabase):
            out.append("db must be a RockDatabase")
            return out
        if len(self.db) == 0:
            out.append("legend must contain at least one entry")
        seen: set[str] = set()
        for e in self.db.entries:
            if not str(e.name).strip():
                out.append("every entry must have a non-empty name")
            if not (e.rho_min >= 0 and e.rho_max > e.rho_min):
                out.append(
                    f"entry {e.name!r}: rho_min must be >= 0 and < rho_max"
                )
            if e.pattern_id and not e.pattern_source:
                out.append(
                    f"entry {e.name!r}: pattern_id set without pattern_source"
                )
            key = e.name.strip().casefold()
            if key in seen:
                out.append(f"duplicate entry name {e.name!r}")
            seen.add(key)
        return out

    def validate(self) -> None:
        problems = self.issues()
        if problems:
            raise GeologyLegendValidationError("; ".join(problems[:6]))


# ---------------------------------------------------------------------------
# JSON I/O
# ---------------------------------------------------------------------------


def legend_to_dict(legend: GeologyLegend) -> dict[str, Any]:
    """Return a canonical JSON-safe dict for *legend*."""
    legend.validate()
    payload = legend.db.to_dict()
    return {
        "pcgl_version": legend.pcgl_version,
        "document_id": legend.document_id,
        "created_at": legend.created_at,
        "created_by": legend.created_by,
        "title": legend.title,
        "description": legend.description,
        "rho_unit": legend.rho_unit,
        "metadata": dict(legend.metadata),
        "entries": payload["entries"],
    }


def legend_from_dict(payload: dict[str, Any]) -> GeologyLegend:
    """Build and validate a :class:`GeologyLegend` from a dict."""
    if not isinstance(payload, dict):
        raise GeologyLegendValidationError("PCGL payload must be a JSON object")
    db = RockDatabase.from_dict(
        {"entries": payload.get("entries", [])},
        metadata={"origin": "pcgl"},
    )
    legend = GeologyLegend(
        document_id=str(payload.get("document_id", "pcgl:legend")),
        created_at=str(payload.get("created_at") or _now()),
        created_by=str(payload.get("created_by", "pycsamt")),
        db=db,
        pcgl_version=str(payload.get("pcgl_version", PCGL_VERSION)),
        title=str(payload.get("title", "")),
        description=str(payload.get("description", "")),
        rho_unit=str(payload.get("rho_unit", "ohm.m")),
        metadata=dict(payload.get("metadata", {}) or {}),
    )
    legend.validate()
    return legend


def read_legend(path: PathLike) -> GeologyLegend:
    """Read a ``.pcgl.json`` file."""
    raw = Path(path).read_text(encoding="utf-8")
    return legend_from_dict(json.loads(raw))


def write_legend(legend: GeologyLegend, path: PathLike) -> Path:
    """Write *legend* as canonical ``.pcgl.json``."""
    target = Path(path)
    text = json.dumps(legend_to_dict(legend), ensure_ascii=False, indent=2)
    target.write_text(text + "\n", encoding="utf-8")
    return target


def legend_from_csv(path: PathLike, **kwargs: Any) -> GeologyLegend:
    """Module-level convenience alias for :meth:`GeologyLegend.from_csv`."""
    return GeologyLegend.from_csv(path, **kwargs)


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def _now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")
