# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
ForwardWorker — QThread for 1-D / 2-D / 3-D MT forward modelling.

Runs the appropriate solver (MT1DForward, MT2DForward, MT3DForward) in a
background thread so the GUI stays responsive.  The caller passes a plain
dict of parameters; this worker owns all forward imports so they are never
on the main thread.

Signals
-------
progress(int)      0–100 estimated progress
finished(object)   ForwardResponse | ForwardResponse2D | ForwardResponse3D
error(str)         human-readable error on failure
"""

from __future__ import annotations

import logging
from typing import Any

import numpy as np
from PySide6.QtCore import QThread, Signal

logger = logging.getLogger(__name__)


class ForwardWorker(QThread):
    """Background thread for a single forward-modelling run."""

    progress = Signal(int)
    finished = Signal(object)
    error = Signal(str)

    def __init__(self, params: dict[str, Any], parent=None) -> None:
        super().__init__(parent)
        self._params = params

    # ── QThread entry point ───────────────────────────────────────────────────

    def run(self) -> None:
        dim = self._params.get("dim", "1D")
        try:
            if dim == "1D":
                result = self._run_1d()
            elif dim == "2D":
                result = self._run_2d()
            else:
                result = self._run_3d()
            self.finished.emit(result)
        except Exception as exc:
            logger.exception("Forward worker error")
            self.error.emit(str(exc))

    # ── 1-D solver ────────────────────────────────────────────────────────────

    def _run_1d(self):
        from pycsamt.forward import LayeredModel, MT1DForward
        from pycsamt.forward.noise import add_noise

        p = self._params
        rho = np.asarray(p["resistivity"], dtype=float)
        thick = np.asarray(p["thickness"], dtype=float)
        n_freq = int(p.get("n_freq", 30))
        f_min = max(float(p.get("f_min", 1e-3)), 1e-6)  # guard log10(0)
        f_max = max(float(p.get("f_max", 1e3)), f_min * 10)
        noise = p.get("noise", "none")

        self.progress.emit(10)
        model = LayeredModel(resistivity=rho, thickness=thick)
        freqs = np.logspace(np.log10(f_min), np.log10(f_max), n_freq)

        self.progress.emit(30)
        resp = MT1DForward(freqs).run(model)

        if noise != "none":
            self.progress.emit(70)
            resp = add_noise(resp, noise_model=noise, level=0.05)

        self.progress.emit(100)
        return resp

    # ── 2-D solver ────────────────────────────────────────────────────────────

    def _run_2d(self):
        from pycsamt.forward import Grid2D, MT2DForward

        p = self._params
        nx = int(p.get("nx", 30))
        nz = int(p.get("nz", 20))
        dx = float(p.get("dx", 500.0))
        dz_min = float(p.get("dz_min", 50.0))
        dz_max = float(p.get("dz_max", 1000.0))
        n_pad = int(p.get("n_pad", 5))
        bg_rho = float(p.get("bg_rho", 100.0))
        n_freq = int(p.get("n_freq", 15))
        f_min = max(float(p.get("f_min", 1e-3)), 1e-6)
        f_max = max(float(p.get("f_max", 1e2)), f_min * 10)
        n_sta = int(p.get("n_stations", 10))

        self.progress.emit(10)

        # Build uniform-ish core cells
        dx_core = np.full(nx, dx)
        dz_core = np.geomspace(dz_min, dz_max, nz)

        # Padding strips
        from pycsamt.forward import make_padding

        pad = make_padding(dx, n_pad)
        dx_full = np.concatenate([pad[::-1], dx_core, pad])
        pad_z = make_padding(dz_max, n_pad)
        dz_full = np.concatenate([dz_core, pad_z])

        nz_full = len(dz_full)
        nx_full = len(dx_full)

        resistivity = np.full((nz_full, nx_full), bg_rho)

        # Inject anomaly if provided
        if p.get("anomaly", False):
            # Centre x is in *survey* coordinates (0..nx*dx; the UI default
            # 7500 m is mid-survey for 30 x 500 m). It used to be applied
            # from the outer padding edge, shifting the anomaly left by the
            # whole padding width (~5.9 km with the defaults).
            ax = float(p.get("anom_x", nx * dx / 2))
            az = float(p.get("anom_z", 500.0))
            aw = float(p.get("anom_w", 2000.0))
            ah = float(p.get("anom_h", 1000.0))
            arho = float(p.get("anom_rho", 10.0))
            xn = np.concatenate([[0], np.cumsum(dx_full)])
            zn = np.concatenate([[0], np.cumsum(dz_full)])
            ix = _cells_near(xn, np.sum(pad) + ax, aw / 2)
            iz = _cells_in(zn, az, az + ah)
            if ix.size and iz.size:
                resistivity[np.ix_(iz, ix)] = arho

        # Station positions in core region (centred)
        x_start = np.sum(pad) + dx * 0.5
        x_end = np.sum(pad) + dx * (nx - 0.5)
        sta_x = np.linspace(x_start, x_end, n_sta)

        grid = Grid2D(
            dx=dx_full,
            dz=dz_full,
            resistivity=resistivity,
            x_stations=sta_x,
            n_pad=n_pad,
        )

        freqs = np.logspace(np.log10(f_min), np.log10(f_max), n_freq)

        self.progress.emit(30)
        solver = MT2DForward(freqs, grid)
        self.progress.emit(50)
        resp = solver.run()
        self.progress.emit(100)
        return resp

    # ── 3-D solver ────────────────────────────────────────────────────────────

    def _run_3d(self):
        from pycsamt.forward import Grid3D, MT3DForward

        p = self._params
        nx = int(p.get("nx", 12))
        ny = int(p.get("ny", 12))
        nz = int(p.get("nz", 10))
        dx = float(p.get("dx", 1000.0))
        dy = float(p.get("dy", 1000.0))
        n_pad = int(p.get("n_pad", 4))
        bg_rho = float(p.get("bg_rho", 100.0))
        n_freq = int(p.get("n_freq", 8))
        f_min = max(float(p.get("f_min", 1e-3)), 1e-6)
        f_max = max(float(p.get("f_max", 1e1)), f_min * 10)
        n_sx = int(p.get("n_sx", 4))
        n_sy = int(p.get("n_sy", 4))

        x_max = nx * dx
        y_max = ny * dy
        z_max = nz * 200.0  # nominal depth extent

        self.progress.emit(10)

        # Grid3D.halfspace builds the full padded grid internally
        grid = Grid3D.halfspace(
            rho=bg_rho,
            nx=nx,
            ny=ny,
            nz=nz,
            x_max=x_max,
            y_max=y_max,
            z_max=z_max,
            n_pad=n_pad,
            nx_stations=n_sx,
            ny_stations=n_sy,
        )

        self.progress.emit(20)

        # Inject anomaly if provided
        if p.get("anomaly", False):
            ax = float(p.get("anom_x", x_max / 2))
            ay = float(p.get("anom_y", y_max / 2))
            az = float(p.get("anom_z", 500.0))
            arho = float(p.get("anom_rho", 10.0))
            xn = np.concatenate([[0], np.cumsum(grid.dx)])
            yn = np.concatenate([[0], np.cumsum(grid.dy)])
            zn = np.concatenate([[0], np.cumsum(grid.dz)])
            # The UI's anomaly centre is in *survey* coordinates (0..nx*dx,
            # default = mid-survey); the node arrays start at the outer edge
            # of the padding. Without this offset the default anomaly landed
            # in the padding and -- combined with the old "cell must lie
            # entirely inside the window" test, which the wide padding cells
            # never pass -- was silently never inserted at all.
            gp = int(getattr(grid, "n_pad", n_pad))
            ax_abs = xn[gp] + ax
            ay_abs = yn[gp] + ay
            half = min(dx, dy)
            ix = _cells_near(xn, ax_abs, half)
            iy = _cells_near(yn, ay_abs, half)
            iz = _cells_in(zn, az, az + 2000.0)
            if ix.size and iy.size and iz.size:
                grid.resistivity[np.ix_(iz, iy, ix)] = arho

        freqs = np.logspace(np.log10(f_min), np.log10(f_max), n_freq)

        self.progress.emit(40)
        solver = MT3DForward(freqs, grid)
        self.progress.emit(60)
        resp = solver.run()
        self.progress.emit(100)
        return resp


def _cells_near(nodes: np.ndarray, centre: float, half: float) -> np.ndarray:
    """Indices of cells whose centre lies within ``centre ± half``.

    Falls back to the single nearest cell, so an anomaly narrower than the
    local cell size (e.g. in coarse padding) is still represented -- but
    only when *centre* lies inside the grid; outside it, nothing is returned
    rather than clamping the body onto an edge cell.
    """
    centres = 0.5 * (nodes[:-1] + nodes[1:])
    idx = np.where(np.abs(centres - centre) <= half)[0]
    if idx.size:
        return idx
    if nodes[0] <= centre <= nodes[-1]:
        return np.array([int(np.argmin(np.abs(centres - centre)))])
    return np.array([], dtype=int)


def _cells_in(nodes: np.ndarray, top: float, bottom: float) -> np.ndarray:
    """Indices of cells whose centre lies in ``[top, bottom]``.

    Falls back to the cell containing *top* when the interval is thinner
    than the local cells; empty when *top* is below the grid.
    """
    centres = 0.5 * (nodes[:-1] + nodes[1:])
    idx = np.where((centres >= top) & (centres <= bottom))[0]
    if idx.size:
        return idx
    if nodes[0] <= top < nodes[-1]:
        return np.array([int(np.searchsorted(nodes, top, side="right") - 1)])
    return np.array([], dtype=int)
