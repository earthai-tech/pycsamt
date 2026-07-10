# Author: LKouadio  <etanoyau@gmail.com>
# License: LGPL-3.0

"""
pycsamt.metadata

Metadata classes for citations, authorship, and software/licensing info.
"""
from __future__ import annotations

import re
from dataclasses import asdict, dataclass, field
from typing import Any


class MetadataError(Exception):
    """Base exception for metadata classes."""
    pass

@dataclass
class Reference:
    """
    Bibliographic reference information.

    Attributes
    ----------
    author : str
        Author list, e.g. "Doe, J.; Smith, A.".
    title : str
        Title of the work.
    journal : Optional[str]
        Journal or publication name.
    year : Optional[int]
        Year of publication.
    volume : Optional[str]
        Volume number or identifier.
    pages : Optional[str]
        Page range, e.g. "123--130".
    doi : Optional[str]
        Digital Object Identifier, e.g. "10.1000/xyz123".
    extra : Dict[str, Any]
        Any additional fields.
    """
    author: str
    title: str
    journal: str | None = None
    year: int | None = None
    volume: str | None = None
    pages: str | None = None
    doi: str | None = None
    extra: dict[str, Any] = field(default_factory=dict)

    DOI_PATTERN = re.compile(r"^10\.\d{4,9}/[-._;()/:A-Z0-9]+$", re.IGNORECASE)

    def __post_init__(self):
        if self.year is not None:
            if not isinstance(self.year, int) or self.year < 0:
                raise MetadataError(f"Invalid year: {self.year}")
        if self.doi:
            if not self.DOI_PATTERN.match(self.doi):
                raise MetadataError(f"Invalid DOI format: {self.doi}")

    def to_dict(self) -> dict[str, Any]:
        """
        Return all reference fields as a dict.
        """
        base = asdict(self)
        base.update(self.extra)
        base.pop('extra', None)
        return base

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> Reference:
        """
        Construct a Reference from a dict.
        """
        extra = {k: v for k, v in data.items()
                 if k not in {'author','title','journal','year','volume','pages','doi'}}
        return cls(
            author=data['author'],
            title=data['title'],
            journal=data.get('journal'),
            year=data.get('year'),
            volume=data.get('volume'),
            pages=data.get('pages'),
            doi=data.get('doi'),
            extra=extra
        )

    def to_bibtex(self, key: str) -> str:
        """
        Export this reference as a BibTeX entry.

        Parameters
        ----------
        key : str
            Citation key, e.g. "Doe2021".

        Returns
        -------
        str
            Multi-line BibTeX entry.
        """
        fields = []
        for attr in ['author','title','journal','year','volume','pages','doi']:
            val = getattr(self, attr)
            if val:
                fields.append(f"  {attr} = {{{val}}},")
        for k, v in self.extra.items():
            fields.append(f"  {k} = {{{v}}},")
        body = "\n".join(fields)
        return f"@article{{{key},\n{body}\n}}"

    def __repr__(self) -> str:
        return (
            f"<Reference title={self.title!r} "
            f"author={self.author!r} year={self.year!r}>"
        )


@dataclass
class CopyrightInfo:
    """
    Copyright and usage conditions.

    Attributes
    ----------
    release_status : str
        e.g. "open", "public", "proprietary".
    conditions_of_use : str
        Text describing usage terms.
    reference : Reference
        Citation for data source.
    extra : Dict[str, Any]
        Additional fields.
    """
    release_status: str
    conditions_of_use: str
    reference: Reference = field(default_factory=Reference)
    extra: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d['reference'] = self.reference.to_dict()
        d.update(self.extra)
        d.pop('extra', None)
        return d

    def __repr__(self) -> str:
        return (
            f"<CopyrightInfo status={self.release_status!r}>")

@dataclass
class Person:
    """
    Person contact information.

    Attributes
    ----------
    name : Optional[str]
    email : Optional[str]
    organization : Optional[str]
    organization_url : Optional[str]
    extra : Dict[str, Any]
    """
    name: str | None = None
    email: str | None = None
    organization: str | None = None
    organization_url: str | None = None
    extra: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d.update(self.extra)
        d.pop('extra', None)
        return d

    def __repr__(self) -> str:
        return f"<Person name={self.name!r} email={self.email!r}>"

@dataclass
class Software:
    """
    Software metadata.

    Attributes
    ----------
    name : str
    version : Optional[str]
    release : Optional[str]
    author : Person
    extra : Dict[str, Any]
    """
    name: str
    version: str | None = None
    release: str | None = None
    author: Person = field(default_factory=Person)
    extra: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d['author'] = self.author.to_dict()
        d.update(self.extra)
        d.pop('extra', None)
        return d

    def __repr__(self) -> str:
        return f"<Software name={self.name!r} version={self.version!r}>"

__all__ = ['Reference', 'CopyrightInfo', 'Person', 'Software']
