# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""Reusable provenance metadata for electromagnetic data products."""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime
from typing import Any

from ..api.property import PyCSAMTObject
from ._property import Person

__all__ = ["ProvenanceMeta"]


@dataclass(repr=False)
class ProvenanceMeta(PyCSAMTObject):
    """Describe creation and stewardship of a scientific data product.

    This object is deliberately independent of EDI and XML syntax.  It maps
    naturally to EMTF provenance fields while remaining reusable by future
    formats and airborne data products.

    Parameters
    ----------
    create_time : datetime or str, optional
        Product creation time. Strings are preserved verbatim so historical
        metadata can be retained without inventing a timezone or precision.
    creating_application : str, optional
        Application that created the data product.
    creator, submitter : Person, optional
        Original creator and person responsible for submission/archiving.
    extra : dict
        Unmodelled provenance fields retained losslessly in memory.
    """

    create_time: datetime | str | None = None
    creating_application: str | None = None
    creator: Person | None = None
    submitter: Person | None = None
    extra: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        if self.creating_application is not None:
            value = str(self.creating_application).strip()
            self.creating_application = value or None
        if self.creator is not None and not isinstance(self.creator, Person):
            raise TypeError("creator must be a pycsamt.metadata.Person")
        if self.submitter is not None and not isinstance(
            self.submitter, Person
        ):
            raise TypeError("submitter must be a pycsamt.metadata.Person")
        if self.create_time is not None and not isinstance(
            self.create_time, (datetime, str)
        ):
            raise TypeError("create_time must be datetime, str, or None")
        if isinstance(self.create_time, str):
            value = self.create_time.strip()
            self.create_time = value or None
        self.extra = dict(self.extra or {})
