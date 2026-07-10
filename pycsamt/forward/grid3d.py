# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
3-D resistivity grid for MT forward modelling.

:class:`Grid3D` stores the subsurface resistivity as a regular
(non-uniform) finite-difference grid in x–y–z together with the
surface station positions needed by :class:`~pycsamt.forward.em3d.MT3DForward`.

Coordinate convention
---------------------
* **x** — easting [m], increasing to the right.
* **y** — northing [m], increasing into the page.
* **z** — depth [m], increasing downward.
* ``resistivity[iz, iy, ix]`` is the resistivity of the cell whose
  top-left-front corner is at ``(x_nodes[ix], y_nodes[iy], z_nodes[iz])``.

Padding
-------
The quasi-3D solver extracts 2-D XZ and YZ slices and runs
:class:`~pycsamt.forward.em2d.MT2DForward` on each.  Those 2-D solvers
need padding to push Dirichlet BCs away from the model region.
:attr:`n_pad` records the padding count added on **each** side in x and y
and at the **bottom** in z, consistent with :class:`~pycsamt.forward.grid2d.Grid2D`.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .grid2d import Grid2D, make_padding

__all__ = ["Grid3D"]


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def _ensure_rng(seed) -> np.random.Generator:
    if isinstance(seed, np.random.Generator):
        return seed
    return np.random.default_rng(seed)


# ─────────────────────────────────────────────────────────────────────────────
# Grid3D
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class Grid3D:
    """Non-uniform 3-D finite-difference grid for MT forward modelling.

    Parameters
    ----------
    dx : array-like, shape (nx,)
        Cell widths in x (easting) [m], including padding.
    dy : array-like, shape (ny,)
        Cell widths in y (northing) [m], including padding.
    dz : array-like, shape (nz,)
        Cell heights in z (depth) [m], including padding.
    resistivity : array-like, shape (nz, ny, nx)
        Cell resistivities [Ω·m], top→bottom, front→back, left→right.
    stations_xy : array-like, shape (n_stations, 2)
        Surface station (x, y) coordinates [m].  Must lie within the
        grid extent.
    n_pad : int
        Padding cells added on each side in x and y and at the bottom
        in z during construction.
    name : str
        Optional label.

    Examples
    --------
    Uniform halfspace::

        >>> g = Grid3D.halfspace(rho=100.0, nx=20, ny=20, nz=15,
        ...     x_max=8_000.0, y_max=8_000.0, z_max=4_000.0,
        ...     nx_stations=5, ny_stations=5)
        >>> g.resistivity.shape   # includes padding
        (23, 28, 28)

    3-D conductive block::

        >>> g = Grid3D.block_anomaly(
        ...     bg_rho=500.0, anomaly_rho=5.0,
        ...     bounds=(2_000.0, 6_000.0, 2_000.0, 6_000.0, 300.0, 1_500.0),
        ...     nx=25, ny=25, nz=18,
        ...     x_max=8_000.0, y_max=8_000.0, z_max=5_000.0,
        ...     nx_stations=5, ny_stations=5)
    """

    dx:          np.ndarray
    dy:          np.ndarray
    dz:          np.ndarray
    resistivity: np.ndarray
    stations_xy: np.ndarray
    n_pad:       int  = 0
    name:        str  = ""

    def __post_init__(self):
        self.dx          = np.asarray(self.dx,          dtype=float)
        self.dy          = np.asarray(self.dy,          dtype=float)
        self.dz          = np.asarray(self.dz,          dtype=float)
        self.resistivity = np.asarray(self.resistivity, dtype=float)
        self.stations_xy = np.asarray(self.stations_xy, dtype=float)

        if self.stations_xy.ndim == 1:
            self.stations_xy = self.stations_xy.reshape(-1, 2)

        if self.resistivity.shape != (self.nz, self.ny, self.nx):
            raise ValueError(
                f"resistivity shape {self.resistivity.shape} does not match "
                f"(nz={self.nz}, ny={self.ny}, nx={self.nx})."
            )
        if np.any(self.resistivity <= 0):
            raise ValueError("All resistivities must be strictly positive.")

        xlo, xhi = self.x_nodes[0], self.x_nodes[-1]
        ylo, yhi = self.y_nodes[0], self.y_nodes[-1]
        bad_x = (self.stations_xy[:, 0] < xlo) | (self.stations_xy[:, 0] > xhi)
        bad_y = (self.stations_xy[:, 1] < ylo) | (self.stations_xy[:, 1] > yhi)
        if bad_x.any() or bad_y.any():
            raise ValueError(
                f"{(bad_x | bad_y).sum()} station(s) outside the grid extent."
            )

    # ── shape ────────────────────────────────────────────────────────────────

    @property
    def nx(self) -> int: return len(self.dx)
    @property
    def ny(self) -> int: return len(self.dy)
    @property
    def nz(self) -> int: return len(self.dz)
    @property
    def n_stations(self) -> int: return len(self.stations_xy)

    # ── coordinates ──────────────────────────────────────────────────────────

    @property
    def x_nodes(self) -> np.ndarray:
        return np.concatenate([[0.0], np.cumsum(self.dx)])

    @property
    def y_nodes(self) -> np.ndarray:
        return np.concatenate([[0.0], np.cumsum(self.dy)])

    @property
    def z_nodes(self) -> np.ndarray:
        return np.concatenate([[0.0], np.cumsum(self.dz)])

    @property
    def x_centers(self) -> np.ndarray:
        xn = self.x_nodes; return 0.5 * (xn[:-1] + xn[1:])

    @property
    def y_centers(self) -> np.ndarray:
        yn = self.y_nodes; return 0.5 * (yn[:-1] + yn[1:])

    @property
    def z_centers(self) -> np.ndarray:
        zn = self.z_nodes; return 0.5 * (zn[:-1] + zn[1:])

    @property
    def x_extent(self) -> float: return float(self.dx.sum())
    @property
    def y_extent(self) -> float: return float(self.dy.sum())
    @property
    def z_extent(self) -> float: return float(self.dz.sum())

    # ── station cell lookup ───────────────────────────────────────────────────

    def _station_x_cells(self) -> np.ndarray:
        """x-cell index for every station."""
        idx = np.searchsorted(self.x_nodes, self.stations_xy[:, 0], side="right") - 1
        return np.clip(idx, 0, self.nx - 1)

    def _station_y_cells(self) -> np.ndarray:
        """y-cell index for every station."""
        idx = np.searchsorted(self.y_nodes, self.stations_xy[:, 1], side="right") - 1
        return np.clip(idx, 0, self.ny - 1)

    # ── 2-D slice extraction (used by the quasi-3D solver) ───────────────────

    def xz_slice(self, yi: int) -> tuple[Grid2D, np.ndarray]:
        """Extract an XZ (east-west) 2-D slice at y-cell *yi*.

        Returns a :class:`~pycsamt.forward.grid2d.Grid2D` whose horizontal
        axis is x (easting) and vertical axis is z (depth).  The station
        x-positions are those of all stations whose y-coordinate falls in
        cell *yi*.

        Returns
        -------
        grid2d : Grid2D
        station_indices : ndarray of int
            Global indices (into :attr:`stations_xy`) of stations in this
            y-row.  ``grid2d.x_stations[k]`` is the x-coordinate of
            ``stations_xy[station_indices[k]]``.
        """
        y_cells = self._station_y_cells()
        mask    = y_cells == yi
        indices = np.where(mask)[0]

        x_st = self.stations_xy[indices, 0] if len(indices) > 0 else np.array(
            [self.x_nodes[self.n_pad + self.nx // 2]]
        )

        rho_xz = self.resistivity[:, yi, :]   # (nz, nx)
        g2d    = Grid2D(
            dx=self.dx, dz=self.dz,
            resistivity=rho_xz,
            x_stations=x_st,
            n_pad=self.n_pad,
        )
        return g2d, indices

    def yz_slice(self, xi: int) -> tuple[Grid2D, np.ndarray]:
        """Extract a YZ (north-south) 2-D slice at x-cell *xi*.

        The ``Grid2D`` horizontal axis is mapped to y (northing); its
        ``dx`` array contains :attr:`dy` and its resistivity is
        ``resistivity[:, :, xi]``.

        Returns
        -------
        grid2d : Grid2D
        station_indices : ndarray of int
            Global indices of stations in this x-column.
            ``grid2d.x_stations[k]`` is the **y**-coordinate of the
            station (used as the horizontal axis in the yz-plane solver).
        """
        x_cells = self._station_x_cells()
        mask    = x_cells == xi
        indices = np.where(mask)[0]

        y_st = self.stations_xy[indices, 1] if len(indices) > 0 else np.array(
            [self.y_nodes[self.n_pad + self.ny // 2]]
        )

        rho_yz = self.resistivity[:, :, xi]   # (nz, ny)
        g2d    = Grid2D(
            dx=self.dy, dz=self.dz,
            resistivity=rho_yz,
            x_stations=y_st,
            n_pad=self.n_pad,
        )
        return g2d, indices

    def column_profile_3d(self, xi: int, yi: int) -> tuple[np.ndarray, np.ndarray]:
        """Return the 1-D resistivity/thickness profile at cell (xi, yi).

        Parameters
        ----------
        xi : int   Column index (x direction).
        yi : int   Row index (y direction).

        Returns
        -------
        rho : ndarray, shape (nz,)
        thick : ndarray, shape (nz-1,)
        """
        return self.resistivity[:, yi, xi], self.dz[:-1]

    # ── conductivity helpers ──────────────────────────────────────────────────

    @property
    def conductivity(self) -> np.ndarray:
        """Cell conductivities [S/m], shape (nz, ny, nx)."""
        return 1.0 / self.resistivity

    # ── core slices (clip padding for visualisation) ──────────────────────────

    @property
    def _cx(self) -> slice:
        p = self.n_pad
        return slice(p, self.nx - p) if p else slice(None)

    @property
    def _cy(self) -> slice:
        p = self.n_pad
        return slice(p, self.ny - p) if p else slice(None)

    @property
    def _cz(self) -> slice:
        p = self.n_pad
        return slice(None, self.nz - p) if p else slice(None)

    # ── visualisation ────────────────────────────────────────────────────────

    def plot(
        self,
        *,
        cmap: str = "jet_r",
        log_scale: bool = True,
        clip_core: bool = True,
        show_stations: bool = True,
        figsize: tuple[float, float] = (13, 4.5),
        vmin: float | None = None,
        vmax: float | None = None,
    ):
        """Plot orthogonal slices: XZ (mid-y), YZ (mid-x), XY (mid-z).

        Parameters
        ----------
        cmap : str
        log_scale : bool
        clip_core : bool
            Clip to non-padding region.
        show_stations : bool
        figsize : (float, float)
        vmin, vmax : float or None

        Returns
        -------
        fig : matplotlib.figure.Figure
        axes : ndarray of Axes, shape (3,)
        """
        import matplotlib.pyplot as plt

        cx, cy, cz = self._cx, self._cy, self._cz
        rho = self.resistivity

        xn = (self.x_nodes[self.n_pad: self.nx + 1 - self.n_pad]
              if clip_core and self.n_pad else self.x_nodes)
        yn = (self.y_nodes[self.n_pad: self.ny + 1 - self.n_pad]
              if clip_core and self.n_pad else self.y_nodes)
        zn = (self.z_nodes[: self.nz + 1 - self.n_pad]
              if clip_core and self.n_pad else self.z_nodes)

        # midpoint indices in core
        mid_y = self.n_pad + (self.ny - 2 * self.n_pad) // 2 if self.n_pad else self.ny // 2
        mid_x = self.n_pad + (self.nx - 2 * self.n_pad) // 2 if self.n_pad else self.nx // 2
        mid_z = (self.nz - self.n_pad) // 2 if self.n_pad else self.nz // 2

        slices = [
            (rho[cz, mid_y, cx], xn, zn, "x (m)", "z (m)",
             f"XZ slice (y = {self.y_centers[mid_y]:.0f} m)"),
            (rho[cz, cy, mid_x], yn, zn, "y (m)", "z (m)",
             f"YZ slice (x = {self.x_centers[mid_x]:.0f} m)"),
            (rho[mid_z, cy, cx], xn, yn, "x (m)", "y (m)",
             f"XY slice (z = {self.z_centers[mid_z]:.0f} m)"),
        ]

        fig, axs = plt.subplots(1, 3, figsize=figsize, constrained_layout=True)
        for ax, (data, h_nodes, v_nodes, xlb, ylb, ttl) in zip(axs, slices):
            d = np.log10(np.maximum(data, 1e-12)) if log_scale else data
            pc = ax.pcolormesh(h_nodes, v_nodes, d, cmap=cmap, shading="flat",
                               vmin=vmin, vmax=vmax)
            if ylb == "z (m)":
                ax.invert_yaxis()
            ax.set_xlabel(xlb, fontsize=8)
            ax.set_ylabel(ylb, fontsize=8)
            ax.set_title(ttl, fontsize=8, pad=4)
            cb = fig.colorbar(pc, ax=ax, pad=0.02, shrink=0.95, aspect=22)
            cb.ax.tick_params(labelsize=7)
            lbl = (r"$\log_{10}\rho$  (Ω·m)" if log_scale else r"$\rho$  (Ω·m)")
            cb.set_label(lbl, fontsize=7)

        if show_stations and self.n_stations > 0:
            # Draw stations on the top XY-slice panel (axs[2])
            axs[2].plot(
                self.stations_xy[:, 0], self.stations_xy[:, 1],
                "v", ms=5, color="k", label="stations", zorder=5,
            )

        fig.suptitle(self.name or "3-D resistivity model", fontsize=10, y=1.01)
        return fig, axs

    def __repr__(self) -> str:
        return (
            f"Grid3D(nx={self.nx}, ny={self.ny}, nz={self.nz}, "
            f"x_ext={self.x_extent:.0f} m, y_ext={self.y_extent:.0f} m, "
            f"z_ext={self.z_extent:.0f} m, "
            f"n_stations={self.n_stations}, n_pad={self.n_pad}"
            + (f", name={self.name!r}" if self.name else "") + ")"
        )

    # ─────────────────────────────────────────────────────────────────────────
    # Constructors
    # ─────────────────────────────────────────────────────────────────────────

    @classmethod
    def halfspace(
        cls,
        rho: float = 100.0,
        *,
        nx: int = 20,
        ny: int = 20,
        nz: int = 15,
        x_max: float = 8_000.0,
        y_max: float = 8_000.0,
        z_max: float = 4_000.0,
        n_pad: int = 8,
        pad_factor: float = 1.3,
        nx_stations: int = 5,
        ny_stations: int = 5,
        name: str = "halfspace",
    ) -> Grid3D:
        """Create a uniform resistivity 3-D halfspace with a regular station grid.

        Parameters
        ----------
        rho : float
            Background resistivity [Ω·m].
        nx, ny, nz : int
            Core cell counts (padding added automatically).
        x_max, y_max : float
            Core horizontal extents [m].
        z_max : float
            Core depth extent [m].
        n_pad : int
            Padding cells per side (x, y) and at the bottom (z).
        pad_factor : float
            Exponential growth factor for padding cells.
        nx_stations, ny_stations : int
            Regular station grid dimensions.
        name : str

        Returns
        -------
        Grid3D
        """
        dx_core = np.full(nx, x_max / nx)
        dy_core = np.full(ny, y_max / ny)
        dz_core = np.full(nz, z_max / nz)

        dx_pad = make_padding(dx_core[0],  n_pad, pad_factor)
        dy_pad = make_padding(dy_core[0],  n_pad, pad_factor)
        dz_pad = make_padding(dz_core[-1], n_pad, pad_factor)

        dx_full = np.concatenate([dx_pad[::-1], dx_core, dx_pad])
        dy_full = np.concatenate([dy_pad[::-1], dy_core, dy_pad])
        dz_full = np.concatenate([dz_core, dz_pad])

        nx_tot, ny_tot, nz_tot = len(dx_full), len(dy_full), len(dz_full)
        rho_grid = np.full((nz_tot, ny_tot, nx_tot), rho)

        x_pad_off = dx_pad[::-1].sum()
        y_pad_off = dy_pad[::-1].sum()

        xs = np.linspace(0.0, x_max, nx_stations) + x_pad_off
        ys = np.linspace(0.0, y_max, ny_stations) + y_pad_off
        gx, gy = np.meshgrid(xs, ys)
        stations = np.column_stack([gx.ravel(), gy.ravel()])

        return cls(
            dx=dx_full, dy=dy_full, dz=dz_full,
            resistivity=rho_grid,
            stations_xy=stations,
            n_pad=n_pad,
            name=name,
        )

    @classmethod
    def block_anomaly(
        cls,
        bg_rho: float = 100.0,
        anomaly_rho: float = 5.0,
        bounds: tuple[float, float, float, float, float, float] = (
            2_000.0, 6_000.0, 2_000.0, 6_000.0, 300.0, 1_500.0
        ),
        *,
        nx: int = 25,
        ny: int = 25,
        nz: int = 18,
        x_max: float = 8_000.0,
        y_max: float = 8_000.0,
        z_max: float = 5_000.0,
        n_pad: int = 8,
        pad_factor: float = 1.3,
        nx_stations: int = 5,
        ny_stations: int = 5,
        name: str = "",
    ) -> Grid3D:
        """Create a background halfspace with one rectangular 3-D anomaly.

        Parameters
        ----------
        bg_rho : float
            Background resistivity [Ω·m].
        anomaly_rho : float
            Anomaly resistivity [Ω·m].
        bounds : (x_lo, x_hi, y_lo, y_hi, z_lo, z_hi)
            Anomaly extents in core coordinates [m].
        nx, ny, nz : int   Core cell counts.
        x_max, y_max, z_max : float   Core extents [m].
        n_pad, pad_factor : int, float
        nx_stations, ny_stations : int
        name : str

        Returns
        -------
        Grid3D

        Examples
        --------
        3-D conductive fault zone::

            >>> g = Grid3D.block_anomaly(
            ...     bg_rho=500.0, anomaly_rho=3.0,
            ...     bounds=(2_500., 5_500., 2_500., 5_500., 200., 1_800.))
        """
        g = cls.halfspace(
            rho=bg_rho,
            nx=nx, ny=ny, nz=nz,
            x_max=x_max, y_max=y_max, z_max=z_max,
            n_pad=n_pad, pad_factor=pad_factor,
            nx_stations=nx_stations, ny_stations=ny_stations,
            name=name or f"bg={bg_rho} + block={anomaly_rho} Ω·m",
        )

        x_pad_off = make_padding(x_max / nx, n_pad, pad_factor).sum()
        y_pad_off = make_padding(y_max / ny, n_pad, pad_factor).sum()

        xc = g.x_centers - x_pad_off
        yc = g.y_centers - y_pad_off
        zc = g.z_centers

        x_lo, x_hi, y_lo, y_hi, z_lo, z_hi = bounds
        col_mask = (xc >= x_lo) & (xc <= x_hi)
        row_mask = (yc >= y_lo) & (yc <= y_hi)
        dep_mask = (zc >= z_lo) & (zc <= z_hi)

        g.resistivity[np.ix_(dep_mask, row_mask, col_mask)] = anomaly_rho
        return g

    @classmethod
    def random_layered(
        cls,
        *,
        n_layers: int = 4,
        nx: int = 20,
        ny: int = 20,
        nz: int = 15,
        x_max: float = 8_000.0,
        y_max: float = 8_000.0,
        z_max: float = 4_000.0,
        n_pad: int = 8,
        pad_factor: float = 1.3,
        nx_stations: int = 5,
        ny_stations: int = 5,
        rho_min: float = 1.0,
        rho_max: float = 10_000.0,
        lateral_variation: bool = True,
        corr_length: float = 2_000.0,
        seed=None,
        name: str = "random_3d",
    ) -> Grid3D:
        """Generate a random layered 3-D model with optional lateral variation.

        Resistivities for *n_layers* horizontal layers are drawn from a
        log-uniform prior.  When *lateral_variation* is ``True``, each
        column's resistivity is smoothly perturbed using a Gaussian random
        field with correlation length *corr_length*.

        Parameters
        ----------
        n_layers : int   Number of horizontal layers.
        nx, ny, nz, x_max, y_max, z_max, n_pad, pad_factor : as above.
        rho_min, rho_max : float   Resistivity bounds [Ω·m].
        lateral_variation : bool
        corr_length : float   GRF correlation length [m].
        seed : int or Generator or None.
        name : str

        Returns
        -------
        Grid3D
        """
        rng = _ensure_rng(seed)

        g = cls.halfspace(
            rho=100.0,
            nx=nx, ny=ny, nz=nz,
            x_max=x_max, y_max=y_max, z_max=z_max,
            n_pad=n_pad, pad_factor=pad_factor,
            nx_stations=nx_stations, ny_stations=ny_stations,
            name=name,
        )
        nx_tot, ny_tot, nz_tot = g.nx, g.ny, g.nz

        # Background layers
        layer_bounds = np.round(np.linspace(0, nz_tot, n_layers + 1)).astype(int)
        log_lo, log_hi = np.log10(rho_min), np.log10(rho_max)
        rho_3d = np.ones((nz_tot, ny_tot, nx_tot))
        for k in range(n_layers):
            rho_k = 10.0 ** rng.uniform(log_lo, log_hi)
            rho_3d[layer_bounds[k]: layer_bounds[k + 1], :, :] = rho_k

        if lateral_variation:
            # 2-D Gaussian random field in xy for each depth layer
            for iz in range(nz_tot):
                # Fast approximation: convolve white noise with a Gaussian kernel
                noise = rng.standard_normal((ny_tot, nx_tot))
                # Gaussian smoothing kernel width (cells)
                dx_mean = float(g.dx.mean())
                dy_mean = float(g.dy.mean())
                sx = max(1.0, corr_length / dx_mean)
                sy = max(1.0, corr_length / dy_mean)
                from scipy.ndimage import gaussian_filter
                smooth = gaussian_filter(noise, sigma=[sy, sx], mode="reflect")
                smooth /= max(smooth.std(), 1e-10)   # normalise to unit std
                log_perturb = 0.3 * smooth            # ±30% log perturbation
                rho_3d[iz] = np.clip(
                    rho_3d[iz] * 10.0 ** log_perturb, rho_min, rho_max
                )

        g.resistivity = rho_3d
        return g
