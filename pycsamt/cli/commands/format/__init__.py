# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.cli.commands.format
===========================

Click command group for the PCSF / PCSM inversion-result format — the
CLI face of :mod:`pycsamt.format`.

Importing this package registers every sub-command on the ``format``
group.  External code only needs::

    from pycsamt.cli.commands.format import fmt

Sub-command modules
-------------------
_base        Root Click group + shared detect / build / write helpers.
convert      ``pycsamt format convert``     — any source → .pcsf / .pcsm.
detect_cmd   ``pycsamt format detect``      — classify a file or folder.
info         ``pycsamt format info``        — full summary of a .pcsf / .pcsm.
validate     ``pycsamt format validate``    — structural + round-trip check.
build_pcbh   ``pycsamt format build-pcbh``  — CSV/XLSX/LAS → .pcbh.json.
build_pcgl   ``pycsamt format build-pcgl``  — CSV → .pcgl.json.
build_pcgs   ``pycsamt format build-pcgs``  — CSV(s) → .pcgs.json.
build_pcpt   ``pycsamt format build-pcpt``  — CSV/XLSX → .pcpt.json.
edi_xml      ``pycsamt format edi-to-xml`` / ``xml-to-edi`` — EDI ⇄ EMTF-XML.
"""

from . import (
    build_pcbh,  # noqa: F401
    build_pcgl,  # noqa: F401
    build_pcgs,  # noqa: F401
    build_pcpt,  # noqa: F401
    convert,  # noqa: F401
    detect_cmd,  # noqa: F401
    edi_xml,  # noqa: F401
    info,  # noqa: F401
    validate,  # noqa: F401
)
from ._base import fmt  # noqa: F401

__all__ = ["fmt"]
