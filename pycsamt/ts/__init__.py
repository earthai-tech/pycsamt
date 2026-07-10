# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
:mod:`pycsamt.ts` — MT field time series to EDI.

Read raw magnetotelluric field recordings (SAMTEX/LiMS
``.ts``, EMSLAB ``.asc``, EDI-embedded ``>TSERIES``, generic
ASCII), estimate band-averaged cross-power spectra, recover
the impedance tensor and tipper, and write complete SEG EDI
files.

Quick start
-----------
>>> from pycsamt.ts import read_ts, ts_to_edi
>>> rec = read_ts("data/MT/TS/kap103as.ts/kap103as.ts")
>>> out = ts_to_edi(rec, savepath="edi_from_ts", nfft=10240)

See Also
--------
pycsamt.seg.spectra.Spectra
    Cross-spectra container reused by this subpackage.
pycsamt.seg.edi.EDIFile
    EDI writer used for the final output.
"""
from __future__ import annotations

from importlib import import_module
from typing import TYPE_CHECKING

__all__ = [
    "TSData",
    "CHANNEL_ORDER",
    "read_ts",
    "read_lims",
    "read_emslab",
    "read_edi_tseries",
    "read_ascii",
    "sniff_format",
    "preprocess",
    "target_frequencies",
    "cross_spectra",
    "ts_to_spectra",
    "ts_to_z",
    "ts_to_edifile",
    "ts_to_edi",
]

_EXPORTS = {
    "TSData": ("_base", "TSData"),
    "CHANNEL_ORDER": ("_base", "CHANNEL_ORDER"),
    "read_ts": ("readers", "read_ts"),
    "read_lims": ("readers", "read_lims"),
    "read_emslab": ("readers", "read_emslab"),
    "read_edi_tseries": ("readers", "read_edi_tseries"),
    "read_ascii": ("readers", "read_ascii"),
    "sniff_format": ("readers", "sniff_format"),
    "preprocess": ("process", "preprocess"),
    "target_frequencies": ("process", "target_frequencies"),
    "cross_spectra": ("process", "cross_spectra"),
    "ts_to_spectra": ("process", "ts_to_spectra"),
    "ts_to_z": ("convert", "ts_to_z"),
    "ts_to_edifile": ("convert", "ts_to_edifile"),
    "ts_to_edi": ("convert", "ts_to_edi"),
}


def __getattr__(name: str):
    if name in _EXPORTS:
        mod, attr = _EXPORTS[name]
        m = import_module(f".{mod}", package=__name__)
        val = getattr(m, attr)
        globals()[name] = val  # cache after first access
        return val
    raise AttributeError(name)


if TYPE_CHECKING:  # pragma: no cover
    from ._base import CHANNEL_ORDER, TSData
    from .convert import ts_to_edi, ts_to_edifile, ts_to_z
    from .process import (
        cross_spectra,
        preprocess,
        target_frequencies,
        ts_to_spectra,
    )
    from .readers import (
        read_ascii,
        read_edi_tseries,
        read_emslab,
        read_lims,
        read_ts,
        sniff_format,
    )
