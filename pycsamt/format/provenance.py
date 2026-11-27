# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Provenance of models predicted by learned methods.

The implementation lives in the standalone ``pcsf`` package
(:mod:`pcsf.core.provenance`), the reference implementation of the
format. This module re-exports it under its pyCSAMT name, so that
existing imports keep working and pyCSAMT objects are the pcsf objects
themselves.

See Also
--------
pcsf.core.provenance : Implementation.
"""

from pcsf.core.provenance import (
    ModelProvenance,
    compute_checkpoint_hash,
)

__all__ = [
    "ModelProvenance",
    "compute_checkpoint_hash",
]
