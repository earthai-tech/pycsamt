# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Reading and writing ``.pcsm``, the text encoding of PCSF.

The implementation lives in the standalone ``pcsf`` package
(:mod:`pcsf.core.text`), the reference implementation of the format.
This module re-exports it under its pyCSAMT name, so that existing
imports keep working and pyCSAMT objects are the pcsf objects
themselves.

See Also
--------
pcsf.core.text : Implementation.
"""

from pcsf.core.text import (
    PCSF_VERSION,
    RESISTIVITY_UNIT,
    pcsf_to_pcsm,
    pcsm_to_pcsf,
    peek_kind,
    read_pcsf_or_pcsm,
    read_pcsm,
    write_pcsm,
)

__all__ = [
    "PCSF_VERSION",
    "RESISTIVITY_UNIT",
    "pcsf_to_pcsm",
    "pcsm_to_pcsf",
    "peek_kind",
    "read_pcsf_or_pcsm",
    "read_pcsm",
    "write_pcsm",
]
