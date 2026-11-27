# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Reading and writing ``.pcsf`` (HDF5) files.

The implementation lives in the standalone ``pcsf`` package
(:mod:`pcsf.core.hdf5`), the reference implementation of the format.
This module re-exports it under its pyCSAMT name, so that existing
imports keep working and pyCSAMT objects are the pcsf objects
themselves.

See Also
--------
pcsf.core.hdf5 : Implementation.
"""

from pcsf.core.hdf5 import (
    PCSF_VERSION,
    RESISTIVITY_UNIT,
    read_pcsf,
    write_pcsf,
)

__all__ = [
    "PCSF_VERSION",
    "RESISTIVITY_UNIT",
    "read_pcsf",
    "write_pcsf",
]
