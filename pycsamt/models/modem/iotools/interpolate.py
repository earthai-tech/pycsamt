# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Grid-to-grid model and data interpolation — Python translation of MATLAB ioAscii/.

Translates:
  interpCond_3D.m  → :func:`interp_model3d`
  interpZ_3D.m     → :func:`interp_z3d`

Both functions require **SciPy** (``scipy.interpolate``).

Model interpolation
--------------------
:func:`interp_model3d` resamples a ``ModEmModel3D`` from one grid onto
another using trilinear interpolation (matching MATLAB's ``interp3`` default).
Cells outside the source domain are filled with a user-supplied background
resistivity.

Data interpolation
-------------------
:func:`interp_z3d` resamples ``ImpedanceFile`` data from existing site
locations to a new set of locations using cubic scattered-data interpolation
(matching MATLAB's ``griddata(...,'cubic')``).
"""

from __future__ import annotations

import copy

import numpy as np

from .impedance import ImpedanceFile, ZBlock

__all__ = [
    "interp_model3d",
    "interp_z3d",
]

# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------


def _cell_centers(
    widths: np.ndarray,
    origin_offset: float = 0.0,
) -> np.ndarray:
    """Return cell-centre coordinates from an array of cell widths.

    Equivalent to ``origin + cumsum([0; dw])(1:end-1) + dw/2`` in MATLAB.
    """
    nodes = origin_offset + np.concatenate([[0.0], np.cumsum(widths)])
    return nodes[:-1] + widths / 2.0


def _bg_loge(bg_rho: float | np.ndarray, nz: int) -> np.ndarray:
    """Convert background resistivity to a per-layer ln(ρ) array."""
    bg = np.asarray(bg_rho, dtype=float)
    if bg.ndim == 0:
        bg = np.full(nz, float(bg))
    elif len(bg) != nz:
        bg = np.resize(bg, nz)
    return np.log(np.where(bg > 0, bg, np.finfo(float).tiny))


# -------------------------------------------------------------------------
# Model interpolation
# -------------------------------------------------------------------------


def interp_model3d(
    source: ModEmModel3D,  # noqa: F821
    target: ModEmModel3D,  # noqa: F821
    bg_rho: float | np.ndarray = 1e3,
) -> ModEmModel3D:  # noqa: F821
    """Interpolate a 3-D resistivity model onto a new grid.

    Equivalent to ``interpCond_3D(newgrid, oldCond, bg)`` in MATLAB.
    Trilinear interpolation is used; cells outside the source domain are
    filled with *bg_rho*.

    Parameters
    ----------
    source : ModEmModel3D
        Model to interpolate *from*.
    target : ModEmModel3D
        Model whose grid specifies the destination.  Only its grid
        attributes (``x_widths``, ``y_widths``, ``z_widths``, ``n_air``)
        are used; the ``rho_loge`` array is overwritten in the returned
        copy.
    bg_rho : float or array-like
        Background resistivity in Ω·m used to fill cells outside the
        source domain.  Scalar applies to all layers; a 1-D array of
        length ``target.nz`` applies per layer.

    Returns
    -------
    ModEmModel3D
        New model on the target grid with interpolated values.

    Raises
    ------
    ImportError
        If SciPy is not installed.
    """
    try:
        from scipy.interpolate import RegularGridInterpolator
    except ImportError as exc:
        raise ImportError(
            "interp_model3d requires SciPy: pip install scipy"
        ) from exc

    # Source cell centres (include optional origin offset)
    origin_s = np.asarray(
        getattr(source, "origin", [0.0, 0.0, 0.0]), dtype=float
    )
    zc_s = _cell_centers(source.z_widths, origin_s[2])
    yc_s = _cell_centers(source.y_widths, origin_s[1])
    xc_s = _cell_centers(source.x_widths, origin_s[0])

    # Target cell centres
    origin_t = np.asarray(
        getattr(target, "origin", [0.0, 0.0, 0.0]), dtype=float
    )
    zc_t = _cell_centers(target.z_widths, origin_t[2])
    yc_t = _cell_centers(target.y_widths, origin_t[1])
    xc_t = _cell_centers(target.x_widths, origin_t[0])

    # Build interpolator: axes in (z, y, x) order matching rho_loge shape
    interp = RegularGridInterpolator(
        (zc_s, yc_s, xc_s),
        source.rho_loge,
        method="linear",
        bounds_error=False,
        fill_value=np.nan,
    )

    # Evaluate at all target cell centres
    nz_t, ny_t, nx_t = target.nz, target.ny, target.nx
    zz, yy, xx = np.meshgrid(zc_t, yc_t, xc_t, indexing="ij")  # (nz, ny, nx)
    pts = np.column_stack([zz.ravel(), yy.ravel(), xx.ravel()])
    new_rho = interp(pts).reshape(nz_t, ny_t, nx_t)

    # Fill NaNs per z-layer with background
    bg = _bg_loge(bg_rho, nz_t)
    for k in range(nz_t):
        mask = np.isnan(new_rho[k])
        if mask.any():
            new_rho[k, mask] = bg[k]

    result = copy.copy(target)
    result.rho_loge = new_rho
    result.log_type = getattr(source, "log_type", "LOGE")
    return result


# -------------------------------------------------------------------------
# Data interpolation
# -------------------------------------------------------------------------


def interp_z3d(
    imp: ImpedanceFile,
    new_site_loc: np.ndarray,
    new_site_names: list[str] | None = None,
    pct_error: float | None = None,
) -> ImpedanceFile:
    """Interpolate 3-D MT impedance data to new site locations.

    Equivalent to ``interpZ_3D(data, siteLoc, siteChar, pererr)`` in MATLAB.
    Cubic scattered interpolation is applied independently to the real and
    imaginary parts of each component and each period.

    Parameters
    ----------
    imp : ImpedanceFile
        Source impedance data.
    new_site_loc : array-like, shape (n_new, 2) or (n_new, 3)
        Target site locations.  Columns 0 and 1 are X and Y in metres.
    new_site_names : list of str, optional
        Site codes for the new locations.  Auto-generated as ``S000``,
        ``S001``, … if not given.
    pct_error : float, optional
        If provided, replace error estimates with this percentage of
        ``|Z|`` (e.g. ``pct_error=5`` → 5 % error).  Mirrors the
        *ADD_NOISE* branch of the MATLAB function.

    Returns
    -------
    ImpedanceFile
        New data file with impedances at the target locations.

    Raises
    ------
    ImportError
        If SciPy is not installed.
    """
    try:
        from scipy.interpolate import griddata
    except ImportError as exc:
        raise ImportError(
            "interp_z3d requires SciPy: pip install scipy"
        ) from exc

    new_site_loc = np.asarray(new_site_loc, dtype=float)
    if new_site_loc.ndim == 1:
        new_site_loc = new_site_loc.reshape(1, -1)

    n_new = len(new_site_loc)
    newx = new_site_loc[:, 0]
    newy = new_site_loc[:, 1]

    if new_site_names is None:
        new_site_names = [f"S{k:03d}" for k in range(n_new)]

    new_imp = copy.deepcopy(imp)

    for i, blk in enumerate(imp.blocks):
        oldx = blk.site_loc[:, 0]
        oldy = blk.site_loc[:, 1]
        old_pts = np.column_stack([oldx, oldy])
        new_pts = np.column_stack([newx, newy])

        new_Z = np.zeros((n_new, blk.n_comp), dtype=complex)
        new_Zerr = np.zeros((n_new, blk.n_comp), dtype=float)

        for j in range(blk.n_comp):
            new_Z[:, j] = griddata(
                old_pts, blk.Z[:, j].real, new_pts, method="cubic"
            ) + 1j * griddata(
                old_pts, blk.Z[:, j].imag, new_pts, method="cubic"
            )
            new_Zerr[:, j] = griddata(
                old_pts, blk.Zerr[:, j], new_pts, method="cubic"
            )

        new_imp.blocks[i] = ZBlock(
            period=blk.period,
            site_names=list(new_site_names),
            site_loc=new_site_loc.copy(),
            comp_names=blk.comp_names,
            Z=new_Z,
            Zerr=new_Zerr,
            mode=blk.mode,
        )

    if pct_error is not None:
        frac = pct_error / 100.0
        for blk in new_imp.blocks:
            blk.Zerr = np.abs(blk.Z) * frac

    return new_imp
