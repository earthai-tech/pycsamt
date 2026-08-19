# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""AFMAG-family derived EMTF datatype registration.

Both derived response definitions below are registered through
:func:`~pycsamt.emtf.datatypes.ensure_emtf_datatype_registered`, the
same idempotent-registration helper used by
:mod:`pycsamt.airborne.mobilemt`, rather than each technology adapter
reimplementing its own "already registered compatibly? then no-op;
registered incompatibly? then raise" guard.
"""

from __future__ import annotations

from ...emtf.datatypes import (
    DataTypeDefinition,
    ensure_emtf_datatype_registered,
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


def register_afmag_datatypes() -> tuple[DataTypeDefinition, ...]:
    """Register AFMAG derived response definitions idempotently.

    The tensor response itself reuses the existing EMTF ``TI`` interstation
    magnetic transfer-function definition and therefore is not re-registered.

    Returns
    -------
    (DataTypeDefinition, DataTypeDefinition)
        ``(AFMAG_TILT_DEFINITION, AFMAG_AP_DEFINITION)``, each either
        newly registered or the already-registered definition sharing
        its tag; see
        :func:`~pycsamt.emtf.datatypes.ensure_emtf_datatype_registered`.

    Raises
    ------
    ValueError
        If either definition's tag or code is already registered
        under a materially different definition.
    """
    return (
        ensure_emtf_datatype_registered(AFMAG_TILT_DEFINITION),
        ensure_emtf_datatype_registered(AFMAG_AP_DEFINITION),
    )
