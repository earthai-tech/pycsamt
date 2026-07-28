# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""
Anomaly peak and boundary selection utilities.
"""

import warnings
from collections.abc import Sequence
from typing import (
    Any,
    Optional,
    Union,
)

import numpy as np

__all__ = [
    "select_anomaly_peak",
    "find_position_bounds",
    "find_closest_positions",
    "find_nearest_indices",
]


def select_anomaly_peak(
    resistivity: Sequence[float],
    positions: Optional[Sequence[float]] = None,
    *,
    dipole: Optional[float] = None,
    rank: int = 1,
    user_peak: Optional[float] = None,
    anomaly_infos: Optional[dict[str, np.ndarray]] = None,
    return_bounds: bool = True,
) -> dict[str, Any]:
    """
    Select anomaly peak position and its boundaries.

    Can operate on raw resistivity profile + positions or
    on precomputed anomaly_infos dict with ranking keys.

    Parameters
    ----------
    resistivity : sequence of float
        Resistivity values along profile.
    positions : sequence of float, optional
        Spatial positions corresponding to resistivity.
        If None, will generate from dipole and length.
    dipole : float, optional
        Spacing between successive positions. Required if
        positions is None.
    rank : int, default=1
        Rank of anomaly to select (1 = strongest = min resistivity).
    user_peak : float, optional
        User-specified peak position override.
    anomaly_infos : dict, optional
        Mapping of ranking keys to arrays of resistivity
        bounds as produced by earlier boundary detection.
        Keys are like '<rank>_pk'.
    return_bounds : bool, default=True
        Whether to include pos_bounds and rho_bounds in output.

    Returns
    -------
    result : dict
        "peak_pos": float   - selected peak position
        "rho_peak": float   - resistivity at peak
        "pos_bounds": tuple - (min, max) position bounds
        "rho_bounds": tuple - (min, max) resistivity bounds
        "profile": ndarray  - full resistivity range used

    Raises
    ------
    ValueError
        If inputs insufficient or inconsistent.

    Examples
    --------
    >>> res = [100, 80, 90, 110]
    >>> pos = [0, 1, 2, 3]
    >>> select_anomaly_peak(res, pos)
    {'peak_pos': 1.0, 'rho_peak': 80.0,
     'pos_bounds': (0.0, 3.0),
     'rho_bounds': (100.0, 110.0),
     'profile': array([100.,  80.,  90., 110.])}
    """
    res = np.asarray(resistivity, dtype=float)
    n = res.size
    if n == 0:
        raise ValueError("Empty resistivity array")
    # If anomaly_infos provided, use it
    if anomaly_infos is not None:
        key = f"{rank}_pk"
        matches = [k for k in anomaly_infos if k.startswith(key)]
        if not matches:
            raise KeyError(f"No anomaly_infos key for rank {rank}")
        code = matches[0]
        profile = np.asarray(anomaly_infos[code], float)
        rho_peak = profile.min()
        # estimate peak position from index
        idx = int(np.argmin(profile))
        if positions is None:
            if dipole is None:
                raise ValueError("dipole needed without positions")
            base_pos = 0.0
            positions = base_pos + np.arange(n) * dipole
        else:
            positions = np.asarray(positions, float)
        peak_pos = positions[idx]
    else:
        # raw profile
        profile = res
        # determine positions
        if positions is None:
            if dipole is None:
                raise ValueError("Either positions or dipole must be provided")
            positions = np.arange(n) * dipole
        else:
            positions = np.asarray(positions, float)
            if positions.size != n:
                raise ValueError("positions length must match resistivity")
        # find index of ranked anomaly
        if rank != 1:
            # sort indices by ascending resistivity
            order = np.argsort(profile)
            idx = int(order[rank - 1])
        else:
            idx = int(np.argmin(profile))
        rho_peak = float(profile[idx])
        peak_pos = float(positions[idx])
    # user override
    if user_peak is not None:
        up = float(user_peak)
        mn, mx = positions.min(), positions.max()
        if not (mn <= up <= mx):
            warnings.warn(f"user_peak {up} outside [{mn},{mx}]", stacklevel=2)
        else:
            peak_pos = up
    result = {"peak_pos": peak_pos, "rho_peak": rho_peak, "profile": profile}
    if return_bounds:
        # full extents
        pos_min, pos_max = positions.min(), positions.max()
        rho_min, rho_max = profile.min(), profile.max()
        result["pos_bounds"] = (pos_min, pos_max)
        result["rho_bounds"] = (rho_min, rho_max)
    return result


def find_position_bounds(
    peak: Union[float, str],
    rho_peak: float,
    rho_range: Sequence[float],
    positions: Optional[Sequence[float]] = None,
    dipole: float = 1.0,
) -> dict[str, float]:
    """
    Compute spatial and resistivity bounds around an anomaly peak.

    Parameters
    ----------
    peak : float or str
        Peak position or 'pk<value>'.
    rho_peak : float
        Resistivity at the peak.
    rho_range : sequence of float
        Resistivity values defining the anomaly.
    positions : sequence of float, optional
        Spatial positions corresponding to `rho_range`. If omitted,
        positions are generated centered at `peak` with spacing `dipole`.
    dipole : float, default=1.0
        Spacing between consecutive measurements.

    Returns
    -------
    bounds : dict
        Dictionary with:
        - 'pos_min', 'pos_max': float spatial limits
        - 'rho_min', 'rho_max': float resistivity limits

    Raises
    ------
    ValueError
        On empty inputs or mismatched lengths.
    """

    # parse peak position
    if isinstance(peak, str):
        # extract numeric part
        num = "".join(ch for ch in peak if ch.isdigit() or ch in ".-")
        try:
            pk_pos = float(num)
        except Exception:
            raise ValueError(f"Cannot parse peak position {peak!r}")
    else:
        pk_pos = float(peak)

    # resistivity array
    rho_arr = np.asarray(rho_range, dtype=float)
    n = rho_arr.size
    if n == 0:
        raise ValueError("rho_range must be non-empty")

    # index of peak resistivity
    idxs = np.where(rho_arr == rho_peak)[0]
    idx = int(idxs[0]) if idxs.size else 0

    # compute positions array
    if positions is not None:
        pos_arr = np.asarray(positions, dtype=float)
        if pos_arr.size != n:
            raise ValueError(f"positions length {pos_arr.size} != rho_range length {n}")
    else:
        # generate centered positions
        offsets = (np.arange(n) - idx) * dipole
        pos_arr = pk_pos + offsets

    # bounds
    pos_min = float(pos_arr.min())
    pos_max = float(pos_arr.max())
    rho_min = float(rho_arr.min())
    rho_max = float(rho_arr.max())

    return {
        "pos_min": pos_min,
        "pos_max": pos_max,
        "rho_min": rho_min,
        "rho_max": rho_max,
    }


def find_closest_positions(
    reference: Sequence[float],
    targets: Sequence[float],
    *,
    return_values: bool = False,
) -> Union[list[int], tuple[list[int], list[float]]]:
    """
    For each target, find index of nearest reference element.

    Parameters
    ----------
    reference : sequence of float
        Array to search within.
    targets : sequence of float
        Values to match to closest reference.
    return_values : bool, default=False
        If True, also return matched reference values.

    Returns
    -------
    indices : list of int
    values : list of float, optional
    """
    ref = np.asarray(reference, float)
    tgt = np.asarray(targets, float)
    if ref.ndim != 1 or tgt.ndim != 1:
        raise ValueError("reference and targets must be 1D")
    idxs = []
    vals = []
    for t in tgt:
        i = int(np.argmin(np.abs(ref - t)))
        idxs.append(i)
        vals.append(ref[i])
    if return_values:
        return idxs, vals
    return idxs


def find_nearest_indices(
    reference: Sequence[float],
    targets: Sequence[float],
    *,
    side: str = "both",
    return_values: bool = False,
) -> Union[list[int], tuple[list[int], list[float]]]:
    """
    Find nearest index(es) in `reference` for each `target`.

    Parameters
    ----------
    reference : sequence of float
        1D array of reference points.
    targets : sequence of float
        Values to match in `reference`.
    side : {'both','left','right'}, default 'both'
        Matching direction:
        - 'both': closest by absolute difference.
        - 'left': nearest value <= target.
        - 'right': nearest value >= target.
    return_values : bool, default False
        If True, also return matched reference values.

    Returns
    -------
    indices : list of int
    values : list of float, optional
        Only if `return_values` is True.

    Raises
    ------
    ValueError
        If inputs are not 1D or `side` is invalid.

    Examples
    --------
    >>> ref = [0, 10, 20, 30]
    >>> find_nearest_indices(ref, [12, -5, 25])
    [1, 0, 2]
    >>> find_nearest_indices(ref, [12, -5], side="left")
    [1, 0]
    >>> find_nearest_indices(ref, [12, 35], side="right", return_values=True)
    ([2, 3], [20.0, 30.0])
    """

    ref = np.asarray(reference, dtype=float).flatten()
    tgt = np.asarray(targets, dtype=float).flatten()
    if ref.ndim != 1 or tgt.ndim != 1:
        raise ValueError("reference and targets must be 1D")
    side = side.lower()
    if side not in ("both", "left", "right"):
        raise ValueError(f"Invalid side argument: {side!r}")
    idxs: list[int] = []
    vals: list[float] = []
    for t in tgt:
        if side == "left":
            mask = ref <= t
            if np.any(mask):
                cand = ref[mask]
                diffs = t - cand
                idx0 = np.where(mask)[0][np.argmin(diffs)]
            else:
                idx0 = 0
        elif side == "right":
            mask = ref >= t
            if np.any(mask):
                cand = ref[mask]
                diffs = cand - t
                idx0 = np.where(mask)[0][np.argmin(diffs)]
            else:
                idx0 = ref.size - 1
        else:
            idx0 = int(np.argmin(np.abs(ref - t)))
        idxs.append(int(idx0))
        vals.append(float(ref[idx0]))
    return (idxs, vals) if return_values else idxs
