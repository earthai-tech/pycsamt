# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""PCGS - pyCSAMT Common Geological Structure.

PCGS is a small, portable document for field structural evidence --
planar measurements (strike/dip), linear measurements (trend/plunge),
and fault traces -- collected along one or several survey lines, so it
can travel between projects and applications on its own and be
recognised and loaded straight into Map View's Geology section.

Like :mod:`pycsamt.format.geology` (PCGL) and
:mod:`pycsamt.format.borehole` (PCBH), PCGS is a thin document envelope
(id, timestamps, provenance, versioning) around an existing pycsamt
class -- here, :class:`pycsamt.geology.structural.StructuralModel` --
rather than a parallel data model.

Canonical encoding is UTF-8 JSON, canonical extension ``.pcgs.json``. A
model can also be assembled from up to three CSV files (one each for
planar measurements, linear measurements, and fault traces), the same
tables :meth:`pycsamt.geology.structural.StructuralModel.from_csv`
already accepts -- see that method for the column layout, including the
optional ``line`` column for a multi-line survey.

Example
-------
>>> from pycsamt.format.structure import StructModel
>>> from pycsamt.geology.structural import FaultTrace, StructuralModel
>>> model = StructuralModel(
...     faults=[FaultTrace(x=500.0, dip_deg=70.0, downthrown_side="right")],
... )
>>> doc = StructModel.from_structural_model(model, title="Site A structure")
>>> doc.issues()
[]
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Union

from ..api.property import MetadataMixin, PyCSAMTObject
from ..geology.structural import StructuralModel

__all__ = [
    "PCGS_VERSION",
    "StructModel",
    "StructModelValidationError",
    "read_structure",
    "write_structure",
    "structure_from_csv",
    "structure_to_dict",
    "structure_from_dict",
]

PCGS_VERSION = "0.1.0"

PathLike = Union[str, Path]


class StructModelValidationError(ValueError):
    """Raised when a PCGS document fails validation."""


@dataclass(repr=False)
class StructModel(PyCSAMTObject, MetadataMixin):
    """A versioned, portable structural-geology document.

    Parameters
    ----------
    document_id : str
    created_at : str
        ISO-8601 timestamp.
    created_by : str
    model : StructuralModel
        The structural evidence itself. See
        :class:`pycsamt.geology.structural.StructuralModel`.
    pcgs_version : str
    title, description : str
    """

    document_id: str
    created_at: str
    created_by: str
    model: StructuralModel
    pcgs_version: str = PCGS_VERSION
    title: str = ""
    description: str = ""
    metadata: dict[str, Any] = field(default_factory=dict)

    @property
    def planar(self):
        return self.model.planar

    @property
    def linear(self):
        return self.model.linear

    @property
    def faults(self):
        return self.model.faults

    def __len__(self) -> int:
        return len(self.model.planar) + len(self.model.linear) + len(
            self.model.faults
        )

    @classmethod
    def from_structural_model(
        cls,
        model: StructuralModel | None = None,
        *,
        document_id: str | None = None,
        created_by: str = "pycsamt",
        title: str = "",
        description: str = "",
    ) -> StructModel:
        doc = cls(
            document_id=document_id or f"pcgs:{_now()}",
            created_at=_now(),
            created_by=created_by,
            model=model if model is not None else StructuralModel(),
            title=title,
            description=description,
        )
        doc.validate()
        return doc

    @classmethod
    def from_csv(
        cls,
        *,
        planar_path: PathLike | None = None,
        linear_path: PathLike | None = None,
        faults_path: PathLike | None = None,
        delimiter: str = ",",
        document_id: str | None = None,
        created_by: str = "pycsamt CSV importer",
        title: str = "",
    ) -> StructModel:
        """Load a document from up to three CSV files -- see
        :meth:`pycsamt.geology.structural.StructuralModel.from_csv`."""
        model = StructuralModel.from_csv(
            planar_path=planar_path,
            linear_path=linear_path,
            faults_path=faults_path,
            delimiter=delimiter,
        )
        return cls.from_structural_model(
            model,
            document_id=document_id,
            created_by=created_by,
            title=title,
        )

    def issues(self) -> list[str]:
        out: list[str] = []
        if not str(self.document_id).strip():
            out.append("document_id must be non-empty")
        try:
            parts = tuple(int(p) for p in self.pcgs_version.split("."))
        except (AttributeError, ValueError):
            parts = ()
        if len(parts) != 3 or parts[0] != 0:
            out.append(f"pcgs_version must be compatible with {PCGS_VERSION}")
        if not isinstance(self.model, StructuralModel):
            out.append("model must be a StructuralModel")
        return out

    def validate(self) -> None:
        problems = self.issues()
        if problems:
            raise StructModelValidationError("; ".join(problems[:6]))


# ---------------------------------------------------------------------------
# JSON I/O
# ---------------------------------------------------------------------------


def structure_to_dict(doc: StructModel) -> dict[str, Any]:
    """Return a canonical JSON-safe dict for *doc*."""
    doc.validate()
    payload = doc.model.to_dict()
    return {
        "pcgs_version": doc.pcgs_version,
        "document_id": doc.document_id,
        "created_at": doc.created_at,
        "created_by": doc.created_by,
        "title": doc.title,
        "description": doc.description,
        "metadata": dict(doc.metadata),
        "planar": payload["planar"],
        "linear": payload["linear"],
        "faults": payload["faults"],
    }


def structure_from_dict(payload: dict[str, Any]) -> StructModel:
    """Build and validate a :class:`StructModel` from a dict."""
    if not isinstance(payload, dict):
        raise StructModelValidationError("PCGS payload must be a JSON object")
    model = StructuralModel.from_dict(payload)
    doc = StructModel(
        document_id=str(payload.get("document_id", "pcgs:structure")),
        created_at=str(payload.get("created_at") or _now()),
        created_by=str(payload.get("created_by", "pycsamt")),
        model=model,
        pcgs_version=str(payload.get("pcgs_version", PCGS_VERSION)),
        title=str(payload.get("title", "")),
        description=str(payload.get("description", "")),
        metadata=dict(payload.get("metadata", {}) or {}),
    )
    doc.validate()
    return doc


def read_structure(path: PathLike) -> StructModel:
    """Read a ``.pcgs.json`` file."""
    raw = Path(path).read_text(encoding="utf-8")
    return structure_from_dict(json.loads(raw))


def write_structure(doc: StructModel, path: PathLike) -> Path:
    """Write *doc* as canonical ``.pcgs.json``."""
    target = Path(path)
    text = json.dumps(structure_to_dict(doc), ensure_ascii=False, indent=2)
    target.write_text(text + "\n", encoding="utf-8")
    return target


def structure_from_csv(**kwargs: Any) -> StructModel:
    """Module-level convenience alias for :meth:`StructModel.from_csv`."""
    return StructModel.from_csv(**kwargs)


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def _now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")
