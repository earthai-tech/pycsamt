# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""Format-neutral electromagnetic processing metadata."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from ..api.property import PyCSAMTObject
from ._property import Software

__all__ = [
    "RemoteReferenceMeta",
    "ProcessingMeta",
    "normalize_sign_convention",
]


def normalize_sign_convention(value: str | None) -> str | None:
    """Normalize common Fourier sign-convention spellings.

    Unknown non-empty values are preserved rather than silently coerced. This
    is intentionally stricter than the historical SEG helper, which defaulted
    unrecognized values to the positive convention.
    """
    if value is None:
        return None
    raw = str(value).strip()
    if not raw:
        return None
    text = raw.lower().replace("\\omega", "ω")
    text = text.replace("omega", "ω").replace(" w ", " ω ")
    compact = "".join(text.split())
    if "exp(+i" in compact and ("ωt" in compact or "wt" in compact):
        return "exp(+i ω t)"
    if "exp(-i" in compact and ("ωt" in compact or "wt" in compact):
        return "exp(-i ω t)"
    return raw


@dataclass(repr=False)
class RemoteReferenceMeta(PyCSAMTObject):
    """Describe remote-reference processing when it was used."""

    reference_type: str | None = None
    site: str | None = None
    extra: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        for attr in ("reference_type", "site"):
            value = getattr(self, attr)
            if value is not None:
                text = str(value).strip()
                setattr(self, attr, text or None)
        self.extra = dict(self.extra or {})


@dataclass(repr=False)
class ProcessingMeta(PyCSAMTObject):
    """Processing assumptions and software associated with an EMTF product.

    This is the format-neutral canonical model. The historical
    :class:`pycsamt.seg.property.Processing` object remains unchanged for EDI
    compatibility and can be adapted later by the EDI interoperability layer.
    """

    sign_convention: str | None = None
    processed_by: str | None = None
    software: Software | None = None
    remote_reference: RemoteReferenceMeta | None = None
    processing_tag: str | None = None
    run_list: list[str] | None = None
    extra: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        self.sign_convention = normalize_sign_convention(self.sign_convention)
        for attr in ("processed_by", "processing_tag"):
            value = getattr(self, attr)
            if value is not None:
                text = str(value).strip()
                setattr(self, attr, text or None)
        if self.software is not None and not isinstance(
            self.software, Software
        ):
            raise TypeError("software must be a pycsamt.metadata.Software")
        if self.remote_reference is not None and not isinstance(
            self.remote_reference, RemoteReferenceMeta
        ):
            raise TypeError("remote_reference must be RemoteReferenceMeta")
        if self.run_list is not None:
            if not isinstance(self.run_list, list):
                raise TypeError("run_list must be a list of strings")
            normalized: list[str] = []
            for index, value in enumerate(self.run_list):
                if not isinstance(value, str):
                    raise TypeError(f"run_list[{index}] must be str")
                text = value.strip()
                if text:
                    normalized.append(text)
            self.run_list = normalized
        self.extra = dict(self.extra or {})
