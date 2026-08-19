# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""Constants used by the EMTF XML reader.

This module intentionally contains XML vocabulary only. Scientific data-type
semantics remain in :mod:`pycsamt.emtf.datatypes`.
"""

from __future__ import annotations

__all__ = [
    "EMTF_ROOT",
    "ESTIMATE_CODES",
    "PERIOD_UNIT_ALIASES",
    "FREQUENCY_UNIT_ALIASES",
]

EMTF_ROOT = "EM_TF"

ESTIMATE_CODES = {
    "VAR": "variance",
    "INVSIGCOV": "inverse_signal_covariance",
    "RESIDCOV": "residual_covariance",
}

PERIOD_UNIT_ALIASES = frozenset(
    {
        "s",
        "sec",
        "secs",
        "second",
        "seconds",
    }
)
FREQUENCY_UNIT_ALIASES = frozenset({"hz", "hertz"})
