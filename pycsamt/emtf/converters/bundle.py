# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""Small explicit wrappers for EMTF <-> TFBundle interoperability."""

from __future__ import annotations

from ...core.base import TFBundle
from ..document import EMTF

__all__ = ["bundle_to_emtf", "emtf_to_bundle"]


def bundle_to_emtf(bundle: TFBundle) -> EMTF:
    """Return :class:`EMTF` from the established ``TFBundle`` bridge."""
    return EMTF.from_bundle(bundle)


def emtf_to_bundle(document: EMTF) -> TFBundle:
    """Return a backward-compatible ``TFBundle`` view of an EMTF document."""
    if not isinstance(document, EMTF):
        raise TypeError("document must be a pycsamt.emtf.EMTF")
    return document.to_bundle()
