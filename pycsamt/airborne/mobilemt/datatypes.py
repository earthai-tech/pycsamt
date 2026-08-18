# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""MobileMT transfer-function datatype registration."""

from __future__ import annotations

from ...emtf.datatypes import (
    DataTypeDefinition,
    get_emtf_datatype,
    register_emtf_datatype,
)
from .constants import MOBILEMT_ADMITTANCE_CODE, MOBILEMT_ADMITTANCE_TAG

__all__ = [
    "MOBILEMT_ADMITTANCE_DEFINITION",
    "register_mobilemt_datatypes",
]


MOBILEMT_ADMITTANCE_DEFINITION = DataTypeDefinition(
    name=MOBILEMT_ADMITTANCE_CODE,
    tag=MOBILEMT_ADMITTANCE_TAG,
    data_kind="complex",
    input_kind="E",
    output_kind="H",
    units=None,
    intention="primary",
    description=(
        "MobileMT magnetic-to-electric admittance with horizontal electric "
        "inputs and three-component magnetic outputs"
    ),
)


def register_mobilemt_datatypes() -> DataTypeDefinition:
    """Register and return the MobileMT admittance definition.

    Registration is idempotent when the existing definition has the same
    semantic tag. A conflicting global ``Y``/tag registration is surfaced
    rather than silently overwritten.
    """
    existing = get_emtf_datatype(MOBILEMT_ADMITTANCE_TAG)
    if existing is not None:
        if existing.name != MOBILEMT_ADMITTANCE_CODE:
            raise ValueError(
                "mobilemt_admittance is already registered with "
                f"code {existing.name!r}"
            )
        return existing

    by_code = get_emtf_datatype(MOBILEMT_ADMITTANCE_CODE)
    if by_code is not None and by_code.tag != MOBILEMT_ADMITTANCE_TAG:
        raise ValueError(
            "EMTF datatype code 'Y' is already registered for "
            f"{by_code.tag!r}"
        )
    return register_emtf_datatype(MOBILEMT_ADMITTANCE_DEFINITION)
