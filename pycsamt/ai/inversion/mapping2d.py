# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""Geometry-aware mapping between regular AI and Occam2D grids.

An AI inversion commonly predicts a regular array, whereas Occam2D
uses unequal finite-element cells grouped into inversion parameters.
This module interpolates an AI grid to physical Occam earth-cell
centres and then computes an area-weighted mean for every parameter
group defined by ``OccamModel.layers``.
"""

from __future__ import annotations

import numpy as np

from ...compat.sklearn import validate_params
from ...models.occam2d import OccamMesh, OccamModel
from .schema import GRID_MAPPING_SCHEMA

__all__ = ["map_ai_grid_to_occam"]


@validate_params(GRID_MAPPING_SCHEMA)
def map_ai_grid_to_occam(
    grid: np.ndarray,
    model: OccamModel,
    mesh: OccamMesh,
    x_coordinates: np.ndarray | None = None,
    z_coordinates: np.ndarray | None = None,
) -> np.ndarray:
    r"""Map a regular AI grid to Occam parameter order.

    AI values are bilinearly interpolated to physical Occam cell
    centres. Values within each model-parameter rectangle are then
    averaged with finite-element cell area as weight:

    .. math::

        \bar m_j
        =
        \frac{\sum_{k \in G_j} m_k A_k}
        {\sum_{k \in G_j} A_k},

    where :math:`G_j` is parameter group :math:`j` and :math:`A_k`
    is the area of mesh cell :math:`k`. Air layers are excluded, and
    the returned vector follows the layer-major order used by Occam
    startup and iteration files.

    Parameters
    ----------
    grid : array-like of float, shape (n_depth, n_horizontal)
        Finite AI values on a regular or explicitly coordinated grid.
        Values are typically log10 resistivity or predictive standard
        deviation.
    model : OccamModel
        Populated Occam model definition. Horizontal parameter codes
        give the number of mesh cells spanned, and ``n_merge`` gives
        the number of earth rows spanned vertically.
    mesh : OccamMesh
        Populated Occam finite-element mesh supplying physical cell
        widths, layer thicknesses, node positions, and air-layer count.
    x_coordinates : array-like of float, optional
        Strictly increasing horizontal coordinates of AI grid columns,
        in the same coordinate system as ``mesh.x_nodes``. The length
        must equal ``grid.shape[1]``. If omitted, uniformly spaced cell
        centres spanning the complete Occam horizontal domain are used.
    z_coordinates : array-like of float, optional
        Strictly increasing AI depth coordinates, positive downward
        from the earth surface. The length must equal ``grid.shape[0]``.
        If omitted, uniformly spaced cell centres spanning the Occam
        earth domain are used.

    Returns
    -------
    numpy.ndarray of float, shape (model.n_params,)
        Area-weighted values in Occam layer-major parameter order.

    Raises
    ------
    TypeError
        Raised when ``model`` or ``mesh`` lacks the required Occam
        geometry interface.
    ValueError
        Raised for invalid grids or coordinates, non-positive mesh
        dimensions, model groups inconsistent with the mesh, or a
        returned parameter count inconsistent with ``model.n_params``.

    Notes
    -----
    Values outside the supplied AI coordinate range use nearest-edge
    extrapolation through :func:`numpy.interp`. Publication workflows
    should normally supply coordinates spanning the complete inversion
    domain so extrapolation is unnecessary.

    The mapping averages log10 resistivity when ``grid`` is in log10
    units. It therefore preserves the parameterization optimized by
    Occam rather than arithmetic resistivity.

    When ``grid`` contains predictive standard deviation, the same
    area-weighted spatial averaging is applied. This treats uncertainty
    as a spatial field and does not assume independent AI pixels; it is
    therefore not a standard-error reduction by the number of cells.

    See Also
    --------
    pycsamt.ai.inversion.duhi2d.DUHIInverter2D
        Uses this mapper for AI means and standard deviations.
    pycsamt.models.occam2d.OccamModel
        Defines the parameter grouping traversed here.
    pycsamt.models.occam2d.OccamMesh
        Defines the physical cell areas used as weights.

    Examples
    --------
    Map a grid after reading an Occam project:

    >>> from pycsamt.ai.inversion.mapping2d import (
    ...     map_ai_grid_to_occam,
    ... )
    >>> from pycsamt.models.occam2d import OccamMesh, OccamModel
    >>> mesh = OccamMesh.read("occam_run/Occam2DMesh")
    >>> model = OccamModel.read("occam_run/Occam2DModel")
    >>> parameters = map_ai_grid_to_occam(  # doctest: +SKIP
    ...     ai_grid,
    ...     model,
    ...     mesh,
    ...     x_coordinates=ai_x,
    ...     z_coordinates=ai_z,
    ... )
    >>> parameters.shape == (model.n_params,)  # doctest: +SKIP
    True
    """
    values = np.asarray(grid, dtype=float)
    if values.ndim != 2 or values.size == 0:
        raise ValueError("grid must be a non-empty two-dimensional array")
    if not np.all(np.isfinite(values)):
        raise ValueError("grid must contain only finite values")

    _validate_geometry_interfaces(model, mesh)
    x_widths = np.asarray(mesh.x_widths, dtype=float).reshape(-1)
    z_widths = np.asarray(mesh.z_widths, dtype=float).reshape(-1)
    if (
        x_widths.size == 0
        or z_widths.size == 0
        or np.any(x_widths <= 0)
        or np.any(z_widths <= 0)
        or not np.all(np.isfinite(x_widths))
        or not np.all(np.isfinite(z_widths))
    ):
        raise ValueError("mesh widths must be finite and positive")

    n_air = int(mesh.n_airlayers)
    if n_air < 0 or n_air >= z_widths.size:
        raise ValueError("mesh must contain at least one earth layer")
    earth_widths = z_widths[n_air:]
    earth_depth = float(earth_widths.sum())
    x_extent = float(x_widths.sum())

    source_x = _source_coordinates(
        x_coordinates,
        values.shape[1],
        x_extent,
        "x_coordinates",
    )
    source_z = _source_coordinates(
        z_coordinates,
        values.shape[0],
        earth_depth,
        "z_coordinates",
    )
    target_x = np.cumsum(x_widths) - 0.5 * x_widths
    target_z = np.cumsum(earth_widths) - 0.5 * earth_widths
    cell_values = _interpolate_grid(
        values,
        source_x,
        source_z,
        target_x,
        target_z,
    )
    return _aggregate_parameter_groups(
        cell_values,
        model,
        x_widths,
        earth_widths,
    )


def _validate_geometry_interfaces(model, mesh) -> None:
    """Validate the minimum model and mesh geometry interfaces."""
    if not hasattr(model, "layers") or not hasattr(model, "n_params"):
        raise TypeError("model must provide the OccamModel interface")
    required_mesh = ("x_widths", "z_widths", "n_airlayers")
    missing = [name for name in required_mesh if not hasattr(mesh, name)]
    if missing:
        raise TypeError(
            "mesh does not provide required attributes: " + ", ".join(missing)
        )
    if not model.layers:
        raise ValueError("Occam model has no parameter layers")


def _source_coordinates(
    coordinates,
    size: int,
    extent: float,
    name: str,
) -> np.ndarray:
    """Return validated explicit or inferred source-cell centres."""
    if size < 1 or not np.isfinite(extent) or extent <= 0:
        raise ValueError("source coordinate size and extent must be positive")
    if coordinates is None:
        width = extent / size
        return (np.arange(size, dtype=float) + 0.5) * width
    result = np.asarray(coordinates, dtype=float).reshape(-1)
    if result.size != size:
        raise ValueError(f"{name} length must match its grid dimension")
    if not np.all(np.isfinite(result)) or np.any(np.diff(result) <= 0):
        raise ValueError(f"{name} must be finite and strictly increasing")
    return result


def _interpolate_grid(
    grid: np.ndarray,
    source_x: np.ndarray,
    source_z: np.ndarray,
    target_x: np.ndarray,
    target_z: np.ndarray,
) -> np.ndarray:
    """Bilinearly interpolate a regular grid to target cell centres."""
    horizontal = np.vstack(
        [np.interp(target_x, source_x, row) for row in grid]
    )
    return np.vstack(
        [
            np.interp(target_z, source_z, horizontal[:, column])
            for column in range(target_x.size)
        ]
    ).T


def _aggregate_parameter_groups(
    cell_values: np.ndarray,
    model: OccamModel,
    x_widths: np.ndarray,
    z_widths: np.ndarray,
) -> np.ndarray:
    """Area-average earth-cell values over Occam parameter groups."""
    output = []
    z_start = 0
    n_x = x_widths.size
    n_z = z_widths.size
    for layer_index, layer in enumerate(model.layers):
        try:
            n_merge = int(layer["n_merge"])
            codes = np.asarray(layer["params"], dtype=int).reshape(-1)
            n_columns = int(layer["n_cols"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(
                f"invalid Occam model layer {layer_index}"
            ) from exc
        if n_merge < 1 or n_columns != codes.size or np.any(codes < 1):
            raise ValueError(
                f"invalid parameter grouping in model layer {layer_index}"
            )
        z_end = z_start + n_merge
        if z_end > n_z:
            raise ValueError("model layers exceed the Occam earth mesh")
        if int(codes.sum()) != n_x:
            raise ValueError(
                f"model layer {layer_index} does not span the mesh width"
            )

        x_start = 0
        for span in codes:
            x_end = x_start + int(span)
            block = cell_values[z_start:z_end, x_start:x_end]
            areas = np.outer(
                z_widths[z_start:z_end],
                x_widths[x_start:x_end],
            )
            output.append(float(np.sum(block * areas) / np.sum(areas)))
            x_start = x_end
        z_start = z_end

    if z_start != n_z:
        raise ValueError("model layers do not span the Occam earth mesh")
    mapped = np.asarray(output, dtype=float)
    if mapped.size != int(model.n_params):
        raise ValueError(
            "mapped AI vector does not match Occam parameter count"
        )
    return mapped
