# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
2-D resistivity grid for MT forward modelling.

The :class:`Grid2D` dataclass is the central data object shared by the 2-D
MT forward solver and any downstream processing that needs the subsurface
resistivity structure.  It stores the model as a regular (but non-uniform)
finite-difference grid together with station positions and padding metadata.

Coordinate convention
---------------------
* **x** — horizontal distance along the profile [m], increasing to the right.
* **z** — depth below surface [m], increasing downward.
* Cell *(i, j)* has its top-left corner at ``(x_nodes[j], z_nodes[i])``
  where *i* is the row index (depth) and *j* is the column index (x).
* ``resistivity[i, j]`` is the resistivity of cell *(i, j)* [Ω·m].

Grid layout (``nx = 4``, ``nz = 3``, no padding shown)::

    x_nodes:   0    dx[0]  dx[0]+dx[1]  …
               ┌────┬────┬────┬────┐
    z_nodes[0] │(0,0)│(0,1)│(0,2)│(0,3)│  ← row 0 (shallowest)
               ├────┼────┼────┼────┤
    z_nodes[1] │(1,0)│(1,1)│(1,2)│(1,3)│
               ├────┼────┼────┼────┤
    z_nodes[2] │(2,0)│(2,1)│(2,2)│(2,3)│  ← row nz-1 (deepest)
               └────┴────┴────┴────┘

Padding
-------
The FD solver requires padding cells on both horizontal sides and at the
bottom to push the Dirichlet boundaries far enough from the model region
that edge effects are negligible.  Padding cells are constructed with
exponentially growing widths/depths.  The resistivity in padding cells is
copied from the nearest edge column (horizontal) or deepest row (vertical).

The :attr:`n_pad` attribute records how many padding cells were added on
each side so that station positions and visualisation can clip to the
*core* region (:attr:`core_slice`).

References
----------
Wannamaker, P.E., Stodt, J.A., Rijo, L. (1987). A stable finite element
solution for two-dimensional magnetotelluric modelling. *Geophys. J. Int.*,
88(1), 277-296.

de Groot-Hedlin, C., Constable, S. (1990). Occam's inversion to generate
smooth, two-dimensional models from magnetotelluric data. *Geophysics*,
55(12), 1613-1624.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Sequence, Tuple, Union

import numpy as np

__all__ = [
    "Grid2D",
    "make_padding",
]

# ─────────────────────────────────────────────────────────────────────────────
# Helper
# ─────────────────────────────────────────────────────────────────────────────

def make_padding(
    core_spacing: float,
    n_pad: int,
    factor: float = 1.3,
) -> np.ndarray:
    """Return cell widths for one padding strip.

    The first padding cell has width ``core_spacing * factor`` and each
    subsequent cell is ``factor`` times wider.

    Parameters
    ----------
    core_spacing : float
        Width/depth of the last core cell adjacent to the padding strip [m].
    n_pad : int
        Number of padding cells.
    factor : float
        Growth factor between successive padding cells.  Values in [1.2, 1.5]
        give good results for MT; 1.3 is a safe default.

    Returns
    -------
    numpy.ndarray, shape (n_pad,)
        Cell widths ordered from the core edge outward.
    """
    return core_spacing * (factor ** np.arange(1, n_pad + 1))


def _ensure_rng(seed) -> np.random.Generator:
    if isinstance(seed, np.random.Generator):
        return seed
    return np.random.default_rng(seed)


# ─────────────────────────────────────────────────────────────────────────────
# Grid2D
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class Grid2D:
    """
    Non-uniform 2-D finite-difference grid for MT forward modelling.

    Parameters
    ----------
    dx : array-like, shape (nx,)
        Cell widths in the horizontal (x) direction [m].
    dz : array-like, shape (nz,)
        Cell heights in the vertical (z, depth) direction [m].
    resistivity : array-like, shape (nz, nx)
        Cell resistivities [Ω·m], top-to-bottom, left-to-right.
    x_stations : array-like, shape (n_stations,)
        Horizontal positions of MT stations [m].  Must lie within
        ``[x_nodes[0], x_nodes[-1]]``.
    n_pad : int
        Number of padding cells on each side (left, right) and at the
        bottom that were added during construction.  Used internally to
        identify the core model region.  Set to 0 if no padding is
        present.
    name : str
        Optional label for the model.

    Attributes
    ----------
    nx : int
        Total number of cells in x (includes padding).
    nz : int
        Total number of cells in z (includes padding).
    x_nodes : ndarray, shape (nx+1,)
        x-coordinates of vertical node lines [m].
    z_nodes : ndarray, shape (nz+1,)
        z-coordinates of horizontal node lines (depth) [m].
    x_centers : ndarray, shape (nx,)
        x-coordinates of cell centres [m].
    z_centers : ndarray, shape (nz,)
        z-coordinates (depth) of cell centres [m].
    core_x_slice : slice
        Slice that selects core (non-padding) columns from any ``[:, j]``
        array.
    core_z_slice : slice
        Slice that selects core (non-padding) rows from any ``[i, :]``
        array.

    Examples
    --------
    Uniform halfspace::

        >>> g = Grid2D.halfspace(rho=100.0, nx=30, nz=20,
        ...                       x_max=5_000.0, z_max=3_000.0,
        ...                       n_stations=10, station_x_max=4_000.0)
        >>> g.resistivity.shape
        (20, 30)

    One conductive anomaly::

        >>> g = Grid2D.with_anomaly(
        ...     bg_rho=500.0,
        ...     anomaly_rho=5.0,
        ...     anomaly_bounds=(1_000.0, 3_000.0, 200.0, 800.0),
        ...     nx=40, nz=25,
        ...     x_max=6_000.0, z_max=4_000.0,
        ...     n_stations=12,
        ... )

    1-D layered model extended horizontally::

        >>> from pycsamt.forward.synthetic import LayeredModel
        >>> m = LayeredModel([100, 10, 500], [300, 800])
        >>> g = Grid2D.from_1d_layers(m, nx=40, x_max=6_000.0, n_stations=12)
    """

    dx:          np.ndarray
    dz:          np.ndarray
    resistivity: np.ndarray
    x_stations:  np.ndarray
    n_pad:       int   = 0
    name:        str   = ""

    # ─── post-init normalisation ──────────────────────────────────────────

    def __post_init__(self):
        self.dx          = np.asarray(self.dx,          dtype=float)
        self.dz          = np.asarray(self.dz,          dtype=float)
        self.resistivity = np.asarray(self.resistivity, dtype=float)
        self.x_stations  = np.asarray(self.x_stations,  dtype=float)

        if self.resistivity.shape != (self.nz, self.nx):
            raise ValueError(
                f"resistivity shape {self.resistivity.shape} does not match "
                f"(nz={self.nz}, nx={self.nx})."
            )
        if np.any(self.resistivity <= 0.0):
            raise ValueError("All resistivities must be strictly positive.")
        if np.any(self.dx <= 0.0) or np.any(self.dz <= 0.0):
            raise ValueError("All cell spacings must be strictly positive.")

        x_lo, x_hi = self.x_nodes[0], self.x_nodes[-1]
        out = (self.x_stations < x_lo) | (self.x_stations > x_hi)
        if np.any(out):
            raise ValueError(
                f"{out.sum()} station(s) lie outside the grid "
                f"[{x_lo:.1f}, {x_hi:.1f}] m."
            )

    # ─── shape properties ─────────────────────────────────────────────────

    @property
    def nx(self) -> int:
        """Number of cells in x."""
        return len(self.dx)

    @property
    def nz(self) -> int:
        """Number of cells in z (depth)."""
        return len(self.dz)

    @property
    def n_stations(self) -> int:
        """Number of MT stations."""
        return len(self.x_stations)

    @property
    def n_nodes(self) -> int:
        """Total FD node count  ``(nx+1) × (nz+1)``."""
        return (self.nx + 1) * (self.nz + 1)

    # ─── coordinate arrays ────────────────────────────────────────────────

    @property
    def x_nodes(self) -> np.ndarray:
        """x-coordinates of vertical node lines [m], shape (nx+1,)."""
        return np.concatenate([[0.0], np.cumsum(self.dx)])

    @property
    def z_nodes(self) -> np.ndarray:
        """Depth coordinates of horizontal node lines [m], shape (nz+1,)."""
        return np.concatenate([[0.0], np.cumsum(self.dz)])

    @property
    def x_centers(self) -> np.ndarray:
        """x-coordinates of cell centres [m], shape (nx,)."""
        xn = self.x_nodes
        return 0.5 * (xn[:-1] + xn[1:])

    @property
    def z_centers(self) -> np.ndarray:
        """Depth coordinates of cell centres [m], shape (nz,)."""
        zn = self.z_nodes
        return 0.5 * (zn[:-1] + zn[1:])

    @property
    def x_extent(self) -> float:
        """Total horizontal extent of the grid [m]."""
        return float(self.dx.sum())

    @property
    def z_extent(self) -> float:
        """Total depth extent of the grid [m]."""
        return float(self.dz.sum())

    # ─── padding helpers ──────────────────────────────────────────────────

    @property
    def core_x_slice(self) -> slice:
        """Column slice that selects only core (non-padding) cells."""
        p = self.n_pad
        return slice(p, self.nx - p) if p > 0 else slice(None)

    @property
    def core_z_slice(self) -> slice:
        """Row slice that selects only core (non-padding) cells."""
        p = self.n_pad
        return slice(None, self.nz - p) if p > 0 else slice(None)

    @property
    def core_resistivity(self) -> np.ndarray:
        """Resistivity array clipped to the core (non-padding) region."""
        return self.resistivity[self.core_z_slice, self.core_x_slice]

    # ─── column / row profiles for 1-D boundary conditions ───────────────

    def column_profile(
        self,
        col: int,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Return the 1-D resistivity/thickness profile of column *col*.

        Used by the FD solver to compute Dirichlet boundary conditions at
        the left and right edges from the 1-D MT response of the edge
        columns.

        Parameters
        ----------
        col : int
            Column index (0 = leftmost, nx-1 = rightmost).

        Returns
        -------
        rho : ndarray, shape (nz,)
            Layer resistivities top-to-bottom [Ω·m].
        thick : ndarray, shape (nz-1,)
            Layer thicknesses [m].  (Bottom layer is a halfspace.)
        """
        rho   = self.resistivity[:, col]
        thick = self.dz[:-1]      # nz-1 thicknesses; last cell = halfspace
        return rho, thick

    def station_cell_indices(self) -> np.ndarray:
        """Return the column index of the cell that contains each station.

        Returns
        -------
        indices : ndarray of int, shape (n_stations,)
        """
        xn = self.x_nodes
        idx = np.searchsorted(xn, self.x_stations, side="right") - 1
        return np.clip(idx, 0, self.nx - 1)

    def station_node_indices(self) -> np.ndarray:
        """Return the index of the nearest surface node for each station.

        Surface nodes are at ``z = 0`` and numbered ``0, 1, …, nx``.

        Returns
        -------
        indices : ndarray of int, shape (n_stations,)
        """
        xn = self.x_nodes
        return np.argmin(np.abs(xn[:, None] - self.x_stations[None, :]), axis=0)

    # ─── conductivity helpers ─────────────────────────────────────────────

    @property
    def conductivity(self) -> np.ndarray:
        """Cell conductivities [S/m], shape (nz, nx)."""
        return 1.0 / self.resistivity

    def harmonic_mean_x(self) -> np.ndarray:
        """Harmonic-mean conductivity at x-directed cell interfaces.

        Used by the TM-mode FD assembler.

        Returns
        -------
        ndarray, shape (nz, nx-1)
            ``sigma_hx[i, j]`` is the harmonic mean of ``sigma[i, j]``
            and ``sigma[i, j+1]``.
        """
        s = self.conductivity
        s_l, s_r = s[:, :-1], s[:, 1:]
        return 2.0 * s_l * s_r / (s_l + s_r)

    def harmonic_mean_z(self) -> np.ndarray:
        """Harmonic-mean conductivity at z-directed cell interfaces.

        Returns
        -------
        ndarray, shape (nz-1, nx)
            ``sigma_hz[i, j]`` is the harmonic mean of ``sigma[i, j]``
            and ``sigma[i+1, j]``.
        """
        s = self.conductivity
        s_u, s_d = s[:-1, :], s[1:, :]
        return 2.0 * s_u * s_d / (s_u + s_d)

    # ─── model visualisation ──────────────────────────────────────────────

    def plot(
        self,
        ax=None,
        *,
        log_scale: bool = True,
        cmap: str = "jet_r",
        clip_core: bool = True,
        show_stations: bool = True,
        figsize: Tuple[float, float] = (10, 5),
        vmin: Optional[float] = None,
        vmax: Optional[float] = None,
    ):
        """Plot the 2-D resistivity model.

        Parameters
        ----------
        ax : Axes or None
            Target axes.  Created when not given.
        log_scale : bool
            Plot log₁₀(ρ) rather than ρ.
        cmap : str
            Matplotlib colormap name.
        clip_core : bool
            Clip to the non-padding core region.
        show_stations : bool
            Mark station positions on the surface.
        figsize : (float, float)
            Figure size when *ax* is not provided.
        vmin, vmax : float or None
            Colormap limits in log₁₀(Ω·m) space when *log_scale* is True.

        Returns
        -------
        ax : Axes
        """
        import matplotlib.pyplot as plt
        from matplotlib.colors import LogNorm

        if ax is None:
            _, ax = plt.subplots(figsize=figsize)

        xs = self.core_x_slice if clip_core else slice(None)
        zs = self.core_z_slice if clip_core else slice(None)
        rho = self.resistivity[zs, xs]
        xn  = self.x_nodes[self.n_pad : self.nx + 1 - self.n_pad] if clip_core else self.x_nodes
        zn  = self.z_nodes[: self.nz + 1 - self.n_pad] if clip_core else self.z_nodes

        if log_scale:
            pc = ax.pcolormesh(
                xn, zn, np.log10(rho),
                cmap=cmap, shading="flat",
                vmin=vmin, vmax=vmax,
            )
            cbar = ax.get_figure().colorbar(pc, ax=ax, pad=0.02)
            cbar.set_label("log₁₀(ρ)  [Ω·m]", fontsize=8)
        else:
            norm = LogNorm(
                vmin=vmin or rho.min(),
                vmax=vmax or rho.max(),
            )
            pc = ax.pcolormesh(xn, zn, rho, cmap=cmap, norm=norm, shading="flat")
            cbar = ax.get_figure().colorbar(pc, ax=ax, pad=0.02)
            cbar.set_label("ρ  [Ω·m]", fontsize=8)

        if show_stations and self.n_stations > 0:
            ax.plot(
                self.x_stations, np.zeros_like(self.x_stations),
                "v", color="k", ms=5, zorder=5, label="stations",
            )

        ax.invert_yaxis()
        ax.set_xlabel("Distance  (m)")
        ax.set_ylabel("Depth  (m)")
        ax.set_title(self.name or "2-D resistivity model")
        return ax

    # ─── repr ─────────────────────────────────────────────────────────────

    def __repr__(self) -> str:
        return (
            f"Grid2D(nx={self.nx}, nz={self.nz}, "
            f"x_extent={self.x_extent:.0f} m, z_extent={self.z_extent:.0f} m, "
            f"n_stations={self.n_stations}, n_pad={self.n_pad}"
            + (f", name={self.name!r}" if self.name else "")
            + ")"
        )

    # ─────────────────────────────────────────────────────────────────────────
    # Class-method constructors
    # ─────────────────────────────────────────────────────────────────────────

    @classmethod
    def halfspace(
        cls,
        rho: float = 100.0,
        *,
        nx: int = 40,
        nz: int = 30,
        x_max: float = 10_000.0,
        z_max: float =  5_000.0,
        n_pad: int = 8,
        pad_factor: float = 1.3,
        n_stations: int = 10,
        station_x_max: Optional[float] = None,
        name: str = "halfspace",
    ) -> "Grid2D":
        """Create a uniform resistivity halfspace grid.

        Parameters
        ----------
        rho : float
            Background resistivity [Ω·m].
        nx : int
            Number of *core* cells in x (padding added on both sides).
        nz : int
            Number of *core* cells in z (padding added at the bottom).
        x_max : float
            Total horizontal extent of the *core* region [m].
        z_max : float
            Total depth of the *core* region [m].
        n_pad : int
            Number of padding cells on each side (left, right, bottom).
        pad_factor : float
            Exponential growth factor for padding cell sizes.
        n_stations : int
            Number of equally spaced stations across the core region.
        station_x_max : float or None
            Right-most station position [m].  Defaults to ``x_max``.
        name : str
            Label for the model.

        Returns
        -------
        Grid2D
        """
        dx_core = np.full(nx, x_max / nx)
        dz_core = np.full(nz, z_max / nz)

        dx_pad = make_padding(dx_core[0],  n_pad, pad_factor)
        dz_pad = make_padding(dz_core[-1], n_pad, pad_factor)

        dx_full = np.concatenate([dx_pad[::-1], dx_core, dx_pad])
        dz_full = np.concatenate([dz_core, dz_pad])

        nx_total = len(dx_full)
        nz_total = len(dz_full)
        rho_grid = np.full((nz_total, nx_total), rho)

        x_pad_offset = dx_pad[::-1].sum()
        sx_max = station_x_max if station_x_max is not None else x_max
        x_st = np.linspace(0.0, sx_max, n_stations) + x_pad_offset

        return cls(
            dx=dx_full, dz=dz_full,
            resistivity=rho_grid,
            x_stations=x_st,
            n_pad=n_pad,
            name=name,
        )

    @classmethod
    def from_1d_layers(
        cls,
        model,           # LayeredModel
        *,
        nx: int = 40,
        x_max: float = 10_000.0,
        n_pad: int = 8,
        pad_factor: float = 1.3,
        n_stations: int = 10,
        station_x_max: Optional[float] = None,
        dz_min: float = 10.0,
        name: str = "",
    ) -> "Grid2D":
        """Build a 2-D grid by extending a 1-D :class:`~pycsamt.forward.synthetic.LayeredModel`.

        Every column has the same 1-D layer stack.  This is useful as a
        background model for synthetic anomaly tests and inversion starting
        models.

        Parameters
        ----------
        model : LayeredModel
            Source 1-D model.  ``model.thickness`` defines the cell
            thicknesses; the deepest cell becomes the halfspace.
        nx : int
            Number of core cells in x.
        x_max : float
            Horizontal extent of the core region [m].
        n_pad : int
            Padding cells on each side.
        pad_factor : float
            Padding growth factor.
        n_stations : int
            Number of equally spaced surface stations.
        station_x_max : float or None
            Right-most station x [m].
        dz_min : float
            Minimum cell height [m]; avoids excessively thin cells at depth.
        name : str
            Model label.

        Returns
        -------
        Grid2D
        """
        dx_core = np.full(nx, x_max / nx)
        dx_pad  = make_padding(dx_core[0], n_pad, pad_factor)
        dx_full = np.concatenate([dx_pad[::-1], dx_core, dx_pad])

        # Build dz from model thicknesses; add a deep halfspace cell
        thick = np.asarray(model.thickness, dtype=float)
        dz_core = np.concatenate([np.maximum(thick, dz_min),
                                   [np.maximum(thick[-1] * 3.0, dz_min)]])
        nz_core = len(dz_core)
        dz_pad  = make_padding(dz_core[-1], n_pad, pad_factor)
        dz_full = np.concatenate([dz_core, dz_pad])

        nx_total = len(dx_full)
        nz_total = len(dz_full)

        # Assign resistivity: repeat 1D column across all x positions.
        # model.resistivity has n_layers values; dz_core has len(thick)+1
        # cells (last cell is the halfspace discretisation) — they match.
        rho_col = np.asarray(model.resistivity, dtype=float)
        rho_col_full = np.concatenate(
            [rho_col, np.full(n_pad, rho_col[-1])]
        )
        rho_grid = np.tile(rho_col_full[:, None], (1, nx_total))

        x_pad_offset = dx_pad[::-1].sum()
        sx_max = station_x_max if station_x_max is not None else x_max
        x_st = np.linspace(0.0, sx_max, n_stations) + x_pad_offset

        return cls(
            dx=dx_full, dz=dz_full,
            resistivity=rho_grid,
            x_stations=x_st,
            n_pad=n_pad,
            name=name or f"1D-extended ({model.n_layers} layers)",
        )

    @classmethod
    def with_anomaly(
        cls,
        bg_rho: float = 100.0,
        anomaly_rho: float = 5.0,
        anomaly_bounds: Tuple[float, float, float, float] = (2_000.0, 4_000.0, 200.0, 800.0),
        *,
        nx: int = 50,
        nz: int = 35,
        x_max: float = 8_000.0,
        z_max: float = 4_000.0,
        n_pad: int = 8,
        pad_factor: float = 1.3,
        n_stations: int = 12,
        station_x_max: Optional[float] = None,
        name: str = "",
    ) -> "Grid2D":
        """Create a background halfspace with one rectangular anomaly.

        Parameters
        ----------
        bg_rho : float
            Background resistivity [Ω·m].
        anomaly_rho : float
            Anomaly resistivity [Ω·m].
        anomaly_bounds : (x_lo, x_hi, z_lo, z_hi)
            Anomaly extent in core coordinates [m]:
            horizontal from *x_lo* to *x_hi*, depth from *z_lo* to *z_hi*.
        nx, nz : int
            Core cell counts.
        x_max, z_max : float
            Core region extents [m].
        n_pad : int
            Padding cells per side.
        pad_factor : float
            Padding growth factor.
        n_stations : int
            Number of equally spaced surface stations.
        station_x_max : float or None
            Right-most station x [m].
        name : str
            Model label.

        Returns
        -------
        Grid2D

        Examples
        --------
        Resistive body at mid-depth::

            >>> g = Grid2D.with_anomaly(
            ...     bg_rho=50.0, anomaly_rho=2000.0,
            ...     anomaly_bounds=(2000.0, 4000.0, 500.0, 1500.0),
            ... )
        """
        g = cls.halfspace(
            rho=bg_rho, nx=nx, nz=nz,
            x_max=x_max, z_max=z_max,
            n_pad=n_pad, pad_factor=pad_factor,
            n_stations=n_stations,
            station_x_max=station_x_max,
            name=name or f"bg={bg_rho} + anomaly={anomaly_rho} Ω·m",
        )

        # Identify which cells fall inside the anomaly (core coords)
        x_pad_offset = make_padding(x_max / nx, n_pad, pad_factor).sum()
        z_pad_offset = 0.0  # no top padding

        xc = g.x_centers - x_pad_offset
        zc = g.z_centers

        x_lo, x_hi, z_lo, z_hi = anomaly_bounds
        col_mask = (xc >= x_lo) & (xc <= x_hi)
        row_mask = (zc >= z_lo) & (zc <= z_hi)
        g.resistivity[np.ix_(row_mask, col_mask)] = anomaly_rho

        return g

    @classmethod
    def random(
        cls,
        *,
        nx: int = 40,
        nz: int = 30,
        x_max: float = 8_000.0,
        z_max: float = 4_000.0,
        n_pad: int = 8,
        pad_factor: float = 1.3,
        n_stations: int = 10,
        rho_min: float = 1.0,
        rho_max: float = 10_000.0,
        n_layers: int = 4,
        lateral_variation: bool = True,
        seed=None,
        name: str = "random",
    ) -> "Grid2D":
        """Generate a random layered model with optional lateral variation.

        Resistivities for ``n_layers`` horizontal layers are drawn from a
        log-uniform prior.  When *lateral_variation* is ``True``, each
        column's resistivity is perturbed by ±30 % in log₁₀ space.

        Parameters
        ----------
        nx, nz : int
            Core cell counts.
        x_max, z_max : float
            Core extents [m].
        n_pad : int
            Padding cells per side.
        rho_min, rho_max : float
            Resistivity bounds [Ω·m].
        n_layers : int
            Number of horizontal layers.
        lateral_variation : bool
            Add column-wise log-normal perturbation.
        seed : int or Generator or None
        name : str

        Returns
        -------
        Grid2D
        """
        rng = _ensure_rng(seed)

        # Build grid geometry
        g = cls.halfspace(
            rho=100.0, nx=nx, nz=nz,
            x_max=x_max, z_max=z_max,
            n_pad=n_pad, pad_factor=pad_factor,
            n_stations=n_stations,
            name=name,
        )

        nx_tot = g.nx
        nz_tot = g.nz

        # Assign 1D background
        layer_bounds = np.linspace(0, nz_tot, n_layers + 1).astype(int)
        log_lo = np.log10(rho_min)
        log_hi = np.log10(rho_max)
        rho_bg = np.ones((nz_tot, nx_tot))
        for k in range(n_layers):
            rho_k = 10.0 ** rng.uniform(log_lo, log_hi)
            rho_bg[layer_bounds[k]: layer_bounds[k + 1], :] = rho_k

        if lateral_variation:
            # Column-wise log-normal perturbation, spatially correlated in x
            log_perturb = rng.normal(0.0, 0.3, (nz_tot, nx_tot))
            # Smooth in x with a running mean over 3 cells
            log_perturb = np.apply_along_axis(
                lambda v: np.convolve(v, np.ones(3) / 3, mode="same"),
                axis=1, arr=log_perturb,
            )
            rho_bg = np.clip(rho_bg * 10.0 ** log_perturb, rho_min, rho_max)

        g.resistivity = rho_bg
        return g
