# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Inversion record: data, solver, regularization and misfit.

The implementation lives in the standalone ``pcsf`` package
(:mod:`pcsf.core.record`), the reference implementation of the format.
This module re-exports it under its pyCSAMT name, so that existing
imports keep working and pyCSAMT objects are the pcsf objects
themselves.

See Also
--------
pcsf.core.record : Implementation.
"""

from pcsf.core.record import (
    RECOMMENDED_FIELDS,
    REQUIRED_FIELDS,
    TERMINATIONS,
    InversionRecord,
    attach_inversion_record,
    check_inversion_record,
    get_inversion_record,
)

__all__ = [
    "InversionRecord",
    "RECOMMENDED_FIELDS",
    "REQUIRED_FIELDS",
    "TERMINATIONS",
    "attach_inversion_record",
    "check_inversion_record",
    "get_inversion_record",
]
