# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""Validation helpers for the format-neutral EMTF scientific model."""

from __future__ import annotations

from collections.abc import Iterable

import numpy as np

__all__ = [
    "normalize_channels",
    "normalize_periods",
    "normalize_tf_data",
]


def normalize_channels(channels: Iterable[str] | None) -> tuple[str, ...]:
    """Normalize channel names while preserving order."""
    if channels is None:
        return ()
    out: list[str] = []
    for value in channels:
        name = str(value).strip()
        if not name:
            raise ValueError("channel names must be non-empty")
        if name in out:
            raise ValueError(f"duplicate channel name: {name!r}")
        out.append(name)
    return tuple(out)


def normalize_periods(
    periods,
    *,
    n_periods: int | None = None,
) -> np.ndarray | None:
    """Return a validated positive 1-D period vector."""
    if periods is None:
        return None
    arr = np.asarray(periods, dtype=float)
    if arr.ndim == 0:
        arr = arr.reshape(1)
    if arr.ndim != 1:
        raise ValueError("periods must be a 1-D array")
    if arr.size and (not np.all(np.isfinite(arr)) or np.any(arr <= 0.0)):
        raise ValueError("periods must contain finite positive values")
    if n_periods is not None and arr.size != int(n_periods):
        raise ValueError(
            "period count does not match transfer-function data: "
            f"{arr.size} != {n_periods}"
        )
    return arr


def normalize_tf_data(
    data,
    *,
    n_output: int,
    n_input: int,
) -> np.ndarray:
    """Normalize a TF payload to ``(n_period, n_output, n_input)``.

    Empty channel lists represent scalar 1x1 products.  A 1-D array is
    therefore accepted only for a scalar response and is promoted to
    ``(n_period, 1, 1)``.  A 2-D array matching one output/input matrix is
    interpreted as a single period and promoted to 3-D.
    """
    arr = np.asarray(data)
    if arr.dtype.kind not in "biufc":
        raise TypeError("transfer-function data must be numeric")

    nout = max(1, int(n_output))
    nin = max(1, int(n_input))

    if arr.ndim == 0:
        if (nout, nin) != (1, 1):
            raise ValueError(
                "scalar data are only valid for a scalar 1x1 response"
            )
        arr = arr.reshape(1, 1, 1)
    elif arr.ndim == 1:
        if (nout, nin) != (1, 1):
            raise ValueError(
                "1-D data require scalar input/output channel dimensions"
            )
        arr = arr[:, None, None]
    elif arr.ndim == 2:
        if arr.shape != (nout, nin):
            raise ValueError(
                "2-D data are interpreted as one TF matrix; expected "
                f"{(nout, nin)}, got {arr.shape}"
            )
        arr = arr[None, ...]
    elif arr.ndim == 3:
        if arr.shape[1:] != (nout, nin):
            raise ValueError(
                "TF matrix shape does not match channel dimensions: "
                f"expected (*, {nout}, {nin}), got {arr.shape}"
            )
    else:
        raise ValueError(
            "transfer-function data must be scalar, 1-D, 2-D, or 3-D"
        )

    return arr
