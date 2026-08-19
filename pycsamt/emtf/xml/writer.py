# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""Filesystem/stream writer for the EMTF XML serializer."""

from __future__ import annotations

from os import PathLike
from pathlib import Path
from typing import IO
import xml.etree.ElementTree as ET

from ...api.property import PyCSAMTObject
from ..document import EMTF
from .serializer import EMTFXMLSerializer

__all__ = ["EMTFXMLWriter", "write_emtf_xml"]


class EMTFXMLWriter(PyCSAMTObject):
    """Serialize :class:`~pycsamt.emtf.EMTF` documents as EMTF XML.

    Parameters
    ----------
    strict : bool, default=True
        Reject unsupported/ambiguous scientific content rather than silently
        dropping it.
    precision : int, default=17
        Significant digits used for floating-point response data.
    pretty : bool, default=True
        Indent the generated XML for human readability.
    """

    def __init__(
        self,
        *,
        strict: bool = True,
        precision: int = 17,
        pretty: bool = True,
    ) -> None:
        self.serializer = EMTFXMLSerializer(
            strict=strict,
            precision=precision,
        )
        self.pretty = bool(pretty)

    def to_element(self, document: EMTF) -> ET.Element:
        """Return a serialized root element."""
        root = self.serializer.to_element(document)
        if self.pretty:
            ET.indent(root, space="  ")
        return root

    def dumps(
        self,
        document: EMTF,
        *,
        xml_declaration: bool = True,
        encoding: str = "utf-8",
    ) -> str:
        """Return a Unicode EMTF XML document."""
        root = self.to_element(document)
        body = ET.tostring(root, encoding="unicode", short_empty_elements=True)
        if not xml_declaration:
            return body
        label = str(encoding).replace("_", "-").upper()
        return f'<?xml version="1.0" encoding="{label}"?>\n' + body

    def write(
        self,
        document: EMTF,
        target: str | PathLike[str] | IO[str] | IO[bytes],
        *,
        xml_declaration: bool = True,
        encoding: str = "utf-8",
    ) -> Path | IO[str] | IO[bytes]:
        """Write *document* to a filesystem path or writable stream."""
        text = self.dumps(
            document,
            xml_declaration=xml_declaration,
            encoding=encoding,
        )
        if isinstance(target, (str, PathLike)):
            path = Path(target).expanduser()
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(text, encoding=encoding)
            return path

        write = getattr(target, "write", None)
        if not callable(write):
            raise TypeError("target must be a path or writable file object")
        try:
            write(text)
        except TypeError:
            write(text.encode(encoding))
        return target


def write_emtf_xml(
    document: EMTF,
    target: str | PathLike[str] | IO[str] | IO[bytes],
    *,
    strict: bool = True,
    precision: int = 17,
    pretty: bool = True,
    xml_declaration: bool = True,
    encoding: str = "utf-8",
) -> Path | IO[str] | IO[bytes]:
    """Convenience wrapper around :class:`EMTFXMLWriter`."""
    return EMTFXMLWriter(
        strict=strict,
        precision=precision,
        pretty=pretty,
    ).write(
        document,
        target,
        xml_declaration=xml_declaration,
        encoding=encoding,
    )
