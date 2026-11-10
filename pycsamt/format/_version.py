# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""PCSF format version parsing and compatibility checks.

The implementation lives in the standalone ``pcsf`` package
(:mod:`pcsf.core.versioning`), the reference implementation of the
format. This module re-exports it under its pyCSAMT name, so that
existing imports keep working and pyCSAMT objects are the pcsf objects
themselves.

See Also
--------
pcsf.core.versioning : Implementation.
"""

from pcsf.core.versioning import (
    check_pcsf_version,
    parse_pcsf_version,
)

__all__ = [
    "check_pcsf_version",
    "parse_pcsf_version",
]
