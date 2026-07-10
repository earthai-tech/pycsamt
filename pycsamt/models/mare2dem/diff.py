# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Compare two MARE2DEM resistivity models.

Port of ``diffMARE2DEM_Resistivity.m``.

:func:`diff_resistivity` reads two ``.resistivity`` files with
identical region layouts, applies a differencing function to their
resistivity arrays, and writes the result to a third file.

The default differencing function is the log10 difference:

.. math::

    \\Delta \\rho = \\log_{10}(\\rho_1) - \\log_{10}(\\rho_2)

A custom callable can be provided for other metrics such as
percentage difference or absolute difference.
"""

from __future__ import annotations

from pathlib import Path
from typing import Callable

import numpy as np

from .iotools.resistivity import (
    ResistivityFile,
    read_resistivity,
    write_resistivity,
)

__all__ = ["diff_resistivity"]


def diff_resistivity(
    file1: str | Path,
    file2: str | Path,
    out_file: str | Path,
    *,
    diff_fn: Callable[[np.ndarray, np.ndarray], np.ndarray] | None = None,
) -> ResistivityFile:
    """Difference two MARE2DEM resistivity files and write the result.

    Port of ``diffMARE2DEM_Resistivity.m``.

    Parameters
    ----------
    file1 : path-like
        First ``.resistivity`` file (minuend, or reference model).
    file2 : path-like
        Second ``.resistivity`` file (subtrahend, or inverted model).
    out_file : path-like
        Destination ``.resistivity`` file for the difference model.
    diff_fn : callable or None, default None
        Function ``(A, B) -> C`` applied element-wise to the
        resistivity arrays.  Both *A* and *B* have shape
        ``(n_regions, nrho)``.  The default is:

        .. code-block:: python

            lambda A, B: np.log10(A) - np.log10(B)

        Custom alternatives:

        * Absolute percentage difference:
          ``lambda A, B: np.abs((A - B) / A * 100)``
        * Linear difference:
          ``lambda A, B: A - B``

    Returns
    -------
    ResistivityFile
        The difference resistivity model (also written to *out_file*).

    Raises
    ------
    ValueError
        When the two files have different numbers of regions.

    Examples
    --------
    Default log10 difference:

    >>> from pycsamt.models.mare2dem.diff import diff_resistivity
    >>> dm = diff_resistivity(
    ...     "mare2dem_iter00.resistivity",
    ...     "mare2dem_iter20.resistivity",
    ...     "mare2dem_diff.resistivity",
    ... )
    >>> dm.num_regions
    4812

    Percentage change:

    >>> dm = diff_resistivity(
    ...     "iter00.resistivity",
    ...     "iter20.resistivity",
    ...     "pct_change.resistivity",
    ...     diff_fn=lambda A, B: np.abs((A - B) / A * 100),
    ... )
    """
    if diff_fn is None:
        diff_fn = lambda A, B: np.log10(A) - np.log10(B)  # noqa: E731

    rf1 = read_resistivity(file1)
    rf2 = read_resistivity(file2)

    if rf1.num_regions != rf2.num_regions:
        raise ValueError(
            f"Cannot diff: file1 has {rf1.num_regions} regions but "
            f"file2 has {rf2.num_regions} regions."
        )

    result = ResistivityFile(
        resistivity_file=str(out_file),
        poly_file=rf1.poly_file,
        data_file=rf1.data_file,
        settings_file=rf1.settings_file,
        version=rf1.version,
        anisotropy=rf1.anisotropy,
        target_misfit=rf1.target_misfit,
        max_iterations=rf1.max_iterations,
        bounds_transform=rf1.bounds_transform,
        global_bounds=rf1.global_bounds.copy(),
        roughness_penalty_method=rf1.roughness_penalty_method,
        yz_penalty_weights=rf1.yz_penalty_weights.copy(),
        penalty_cut_weight=rf1.penalty_cut_weight,
        roughness_with_prejudice=rf1.roughness_with_prejudice,
        beta_mgs=rf1.beta_mgs,
        debug_level=rf1.debug_level,
    )

    result.resistivity = diff_fn(rf1.resistivity, rf2.resistivity)
    result.free_parameter = rf1.free_parameter.copy()
    result.bounds = rf1.bounds.copy() if rf1.bounds is not None else np.array([])
    result.prejudice = rf1.prejudice.copy() if rf1.prejudice is not None else np.array([])

    write_resistivity(result, out_file)
    return result
