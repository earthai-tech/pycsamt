# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""SEG-EDI metadata."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from typing import (
    Any,
)

__all__ = [
    "PlainMeta",
    "References",
    "Person",
    "Software",
    "Source",
    "Processing",
    "Copyright",
    "PropertiesMixin",
]


class PlainMeta:
    """Small helper for to_dict/update/clone/validate."""

    def to_dict(self) -> dict[str, Any]:
        try:
            return asdict(self)  # type: ignore[arg-type]
        except Exception:  # pragma: no cover
            return dict(self.__dict__)

    def update(self, /, **kw: Any) -> PlainMeta:
        for k, v in kw.items():
            setattr(self, k, v)
        self.validate()
        return self

    def clone(self, /, **over: Any) -> PlainMeta:
        cp = self.__class__(**self.to_dict())  # type: ignore
        for k, v in over.items():
            setattr(cp, k, v)
        cp.validate()
        return cp

    def validate(self) -> None:
        return None


# Metadata value objects (no EDI formatting, embed-friendly)
@dataclass
class References(PlainMeta):
    author: str | None = None
    title: str | None = None
    journal: str | None = None
    doi: str | None = None
    year: int | None = None
    volume: str | None = None
    pages: str | None = None

    def validate(self) -> None:
        if self.year is not None and self.year < 0:
            raise ValueError("year must be non-negative")


@dataclass
class Person(PlainMeta):
    name: str | None = None
    organization: str | None = None
    email: str | None = None
    role: str | None = None
    phone: str | None = None

    def validate(self) -> None:
        if self.email and "@" not in self.email:
            raise ValueError("email must contain '@'")


@dataclass
class Software(PlainMeta):
    name: str | None = None
    version: str | None = None
    release: str | None = None
    author: Person = field(default_factory=Person)
    url: str | None = None
    license: str | None = None
    description: str | None = None

    def validate(self) -> None:
        if self.name is not None and not str(self.name).strip():
            raise ValueError("name must be non-empty")
        if self.version is not None and not str(self.version).strip():
            raise ValueError("version must be non-empty")
        if self.url is not None:
            u = str(self.url).strip().lower()
            if not (u.startswith("http://") or u.startswith("https://")):
                raise ValueError("url should start with http(s)://")


@dataclass
class Source(PlainMeta):
    project: str | None = None
    survey: str | None = None
    sitename: str | None = None
    creationdate: str = field(
        default_factory=lambda: (
            datetime.now(timezone.utc)
            .replace(microsecond=0)
            .isoformat()
            .replace("+00:00", "Z")
        )
    )
    creatingsoftware: str = "pyCSAMT"
    author: Person = field(default_factory=Person)
    recipient: Person = field(default_factory=Person)
    archive: str | None = None
    reprocessed_by: str | None = None

    def validate(self) -> None:
        if not self.creatingsoftware:
            raise ValueError("creatingsoftware must be set")


@dataclass
class Processing(PlainMeta):
    ProcessingSoftware: Software = field(default_factory=Software)
    processedby: str | None = None
    processingtag: str | None = None
    runlist: list[str] | None = None
    remoteref: str | None = None
    remotesite: str | None = None
    signconvention: str = "exp(+i ω t)"

    @staticmethod
    def _normalize_signconv(s: str) -> str:
        if s is None:
            return "exp(+i ω t)"
        s0 = str(s).strip().lower()
        s0 = s0.replace("\\omega", "ω").replace(" w ", " ω ")
        s0 = s0.replace("w)", "ω)").replace("(w", "(ω")
        s0 = " ".join(s0.split())
        if "exp(+i" in s0:
            return "exp(+i ω t)"
        if "exp(-i" in s0:
            return "exp(-i ω t)"
        return "exp(+i ω t)"

    def validate(self) -> None:
        self.signconvention = self._normalize_signconv(self.signconvention)
        if self.runlist is not None and not isinstance(self.runlist, list):
            raise ValueError("runlist must be a list of str")
        if self.runlist is not None:
            for i, v in enumerate(self.runlist):
                if not isinstance(v, str):
                    raise ValueError(f"runlist[{i}] must be str")


@dataclass
class Copyright(PlainMeta):
    references: References = field(default_factory=References)
    conditions_of_use: str = field(
        default_factory=lambda: (
            "Data in 'data/' may be used to test pyCSAMT. "
            "Commercial redistribution needs citation and "
            "permission from original source."
        )
    )
    release_status: str | None = None
    owner: str | None = None
    contact: str | None = None
    additional_info: str | None = None

    def validate(self) -> None:
        if self.release_status is None:
            return
        allowed = {"open", "public", "proprietary"}
        if self.release_status.lower() not in allowed:
            raise ValueError(f"release_status must be in {sorted(allowed)}")


class PropertiesMixin:
    """
    Lightweight helper to manage Source/Processing/Copyright
    metadata on host classes without coupling to EDI writers.
    """

    # map flat keys to nested holders
    _ROUTE = {
        # Source
        "project": ("Source", "project"),
        "survey": ("Source", "survey"),
        "sitename": ("Source", "sitename"),
        "creationdate": ("Source", "creationdate"),
        "creatingsoftware": ("Source", "creatingsoftware"),
        "author": ("Source", "author"),
        "recipient": ("Source", "recipient"),
        "archive": ("Source", "archive"),
        "reprocessed_by": ("Source", "reprocessed_by"),
        # Processing
        "processedby": ("Processing", "processedby"),
        "processingtag": ("Processing", "processingtag"),
        "runlist": ("Processing", "runlist"),
        "remoteref": ("Processing", "remoteref"),
        "remotesite": ("Processing", "remotesite"),
        "signconvention": ("Processing", "signconvention"),
        "processingsoftware": ("Processing", "ProcessingSoftware"),
        # Copyright
        "references": ("Copyright", "references"),
        "conditions_of_use": ("Copyright", "conditions_of_use"),
        "release_status": ("Copyright", "release_status"),
        "owner": ("Copyright", "owner"),
        "contact": ("Copyright", "contact"),
        "additional_info": ("Copyright", "additional_info"),
    }

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        try:
            super().__init__(*args, **kwargs)  # type: ignore[misc]
        except Exception:
            pass
        # ensure holders exist (preserve if already attached)
        if not hasattr(self, "Source") or not isinstance(self.Source, Source):
            self.Source = Source()
        if not hasattr(self, "Processing") or not isinstance(
            self.Processing, Processing
        ):
            self.Processing = Processing()
        if not hasattr(self, "Copyright") or not isinstance(
            self.Copyright, Copyright
        ):
            self.Copyright = Copyright()

    # ----------------------------
    # ensure / attach / validate
    # ----------------------------
    def ensure_properties(self) -> PropertiesMixin:
        if not isinstance(self.Source, Source):
            self.Source = Source()
        if not isinstance(self.Processing, Processing):
            self.Processing = Processing()
        if not isinstance(self.Copyright, Copyright):
            self.Copyright = Copyright()
        return self

    def attach_properties(
        self,
        *,
        source: Source | None = None,
        processing: Processing | None = None,
        copyright: Copyright | None = None,
    ) -> PropertiesMixin:
        if source is not None:
            if not isinstance(source, Source):
                raise TypeError("source must be Source")
            self.Source = source
        if processing is not None:
            if not isinstance(processing, Processing):
                raise TypeError("processing must be Processing")
            self.Processing = processing
        if copyright is not None:
            if not isinstance(copyright, Copyright):
                raise TypeError("copyright must be Copyright")
            self.Copyright = copyright
        return self

    def validate_properties(self) -> None:
        self.Source.validate()
        self.Processing.validate()
        self.Copyright.validate()

    # ----------------------------
    # update helpers
    # ----------------------------
    def update_properties(
        self,
        *,
        source: dict[str, Any] | None = None,
        processing: dict[str, Any] | None = None,
        copyright: dict[str, Any] | None = None,
        **flat: Any,
    ) -> PropertiesMixin:
        """
        Update nested holders via dicts and/or flat keys.

        Examples
        --------
        self.update_properties(
            source={"project": "A", "survey": "S"}
        )
        self.update_properties(processedby="John")
        """
        self.ensure_properties()

        if source:
            self.Source.update(**source)
        if processing:
            self.Processing.update(**processing)
        if copyright:
            self.Copyright.update(**copyright)

        for k, v in flat.items():
            route = self._ROUTE.get(k)
            if not route:
                # ignore unknown flat keys by default
                continue
            holder_name, attr = route
            holder = getattr(self, holder_name)
            if holder_name == "Processing" and attr == "ProcessingSoftware":
                if isinstance(v, dict):
                    holder.ProcessingSoftware.update(**v)
                elif isinstance(v, Software):
                    holder.ProcessingSoftware = v
                else:
                    # treat as name string
                    holder.ProcessingSoftware.update(name=str(v))
            else:
                current = getattr(holder, attr, None)
                if isinstance(current, PlainMeta) and isinstance(v, dict):
                    current.update(**v)
                else:
                    setattr(holder, attr, v)

        self.validate_properties()
        return self

    # ----------------------------
    # accessors
    # ----------------------------
    def properties_as_dict(self) -> dict[str, Any]:
        return {
            "source": self.Source.to_dict(),
            "processing": self.Processing.to_dict(),
            "copyright": self.Copyright.to_dict(),
        }

    def copy_properties(self) -> dict[str, PlainMeta]:
        return {
            "source": self.Source.clone(),
            "processing": self.Processing.clone(),
            "copyright": self.Copyright.clone(),
        }

    def property_summary(self) -> str:
        s = self.Source
        p = self.Processing
        c = self.Copyright
        parts: list[str] = []
        if s.project or s.survey or s.sitename:
            parts.append(
                f"Source(project={s.project!r}, survey={s.survey!r}, "
                f"sitename={s.sitename!r})"
            )
        if p.processedby or p.processingtag:
            parts.append(
                f"Processing(by={p.processedby!r}, tag={p.processingtag!r})"
            )
        if c.release_status or c.owner:
            parts.append(
                f"Copyright(status={c.release_status!r}, owner={c.owner!r})"
            )
        return "; ".join(parts) if parts else "no properties"
