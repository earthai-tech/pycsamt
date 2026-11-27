# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Sampling and comparison of PCSF models in log10 resistivity.

The implementation lives in the standalone ``pcsf`` package
(:mod:`pcsf.algorithms.comparison`), the reference implementation of
the format. This module re-exports it under its pyCSAMT name, so that
existing imports keep working and pyCSAMT objects are the pcsf objects
themselves.

See Also
--------
pcsf.algorithms.comparison : Implementation.
pcsf.algorithms.sampling : Implementation.
pcsf.algorithms.uncertainty : Implementation.
"""

from pcsf.algorithms.comparison import (
    ComparisonResult,
    compare_log10,
)
from pcsf.algorithms.sampling import (
    sample_log10,
)
from pcsf.algorithms.uncertainty import (
    log_uncertainty,
)

__all__ = [
    "ComparisonResult",
    "compare_log10",
    "log_uncertainty",
    "sample_log10",
]
