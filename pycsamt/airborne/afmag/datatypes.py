# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""AFMAG-family derived EMTF datatype registration."""

from __future__ import annotations

from ...emtf.datatypes import (
    DataTypeDefinition,
    get_emtf_datatype,
    register_emtf_datatype,
)
from .constants import (
    AFMAG_AP_CODE,
    AFMAG_AP_TAG,
    AFMAG_AP_UNITS,
    AFMAG_TILT_CODE,
    AFMAG_TILT_TAG,
)

__all__ = [
    "AFMAG_TILT_DEFINITION",
    "AFMAG_AP_DEFINITION",
    "register_afmag_datatypes",
]


AFMAG_TILT_DEFINITION = DataTypeDefinition(
    name=AFMAG_TILT_CODE,
    tag=AFMAG_TILT_TAG,
    data_kind="real",
    input_kind=None,
    output_kind=None,
    units=None,
    intention="derived",
    description=(
        "Original AFMAG line-direction magnetic polarization tilt response"
    ),
)

AFMAG_AP_DEFINITION = DataTypeDefinition(
    name=AFMAG_AP_CODE,
    tag=AFMAG_AP_TAG,
    data_kind="complex",
    input_kind=None,
    output_kind=None,
    units=AFMAG_AP_UNITS,
    intention="derived",
    description=(
        "AirMt rotationally invariant complex amplification parameter"
    ),
    derived_from="interstation_transfer_functions",
)


def _register(definition: DataTypeDefinition) -> DataTypeDefinition:
    existing = get_emtf_datatype(definition.tag)
    if existing is not None:
        if existing.name != definition.name:
            raise ValueError(
                f"{definition.tag!r} is already registered with code "
                f"{existing.name!r}"
            )
        return existing

    by_code = get_emtf_datatype(definition.name)
    if by_code is not None and by_code.tag != definition.tag:
        raise ValueError(
            f"EMTF datatype code {definition.name!r} is already registered "
            f"for {by_code.tag!r}"
        )
    return register_emtf_datatype(definition)


def register_afmag_datatypes() -> tuple[DataTypeDefinition, ...]:
    """Register AFMAG derived response definitions idempotently.

    The tensor response itself reuses the existing EMTF ``TI`` interstation
    magnetic transfer-function definition and therefore is not re-registered.
    """
    return (
        _register(AFMAG_TILT_DEFINITION),
        _register(AFMAG_AP_DEFINITION),
    )
