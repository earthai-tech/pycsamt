# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""MeshTools3D export — Python translation of the MATLAB ioAscii/ scripts.

Translates:
  write_meshtools3d_model.m  → :func:`write_meshtools3d`

The export creates two companion files from a ``ModEmModel3D``:

* ``<stem>.msh`` — ASCII mesh file (grid dimensions + widths)
* ``<stem>.con`` — ASCII conductivity file (σ values, one per line)

MeshTools3D mesh format (``*.msh``)
--------------------------------------
::

    nx  ny  nzEarth
    ox  oy  oz

    dx1  dx2  dx3  dx4  dx5     (5 per line)
    ...

    dy1  dy2  ...

    dz1  dz2  ...

MeshTools3D conductivity format (``*.con``)
---------------------------------------------
Values written in the order: for j=1..ny, for i=1..nx, for k=1..nzEarth
(z incremented fastest), one value per line.  Air layers are excluded.
"""

from __future__ import annotations

from pathlib import Path
from typing import Union

import numpy as np

__all__ = ["write_meshtools3d"]

PathLike = Union[str, Path]
_WIDTHS_PER_LINE = 5


def write_meshtools3d(
    model: ModEmModel3D,  # noqa: F821
    path: PathLike,
) -> tuple[Path, Path]:
    """Export a ``ModEmModel3D`` to MeshTools3D mesh and conductivity files.

    Equivalent to ``write_meshtools3d_model(fname, Cond)`` in MATLAB.
    The model is converted from LOGE resistivity to linear conductivity
    before writing.  Air layers (the first ``model.n_air`` z-layers) are
    excluded from the export.

    Parameters
    ----------
    model : ModEmModel3D
        Source model.  ``rho_loge`` values are expected as ``ln(Ω·m)``.
    path : str or Path
        Base path for output.  The ``.msh`` and ``.con`` extensions are
        appended automatically.  Any existing extension is stripped.

    Returns
    -------
    msh_path : Path
        Path to the mesh file.
    con_path : Path
        Path to the conductivity file.

    Examples
    --------
    >>> msh, con = write_meshtools3d(model, "output/my_model")
    """
    p = Path(path).with_suffix("")  # strip extension
    msh_path = p.with_suffix(".msh")
    con_path = p.with_suffix(".con")

    nx = model.nx
    ny = model.ny
    n_air = model.n_air
    nz_earth = model.nz - n_air

    origin = np.asarray(getattr(model, "origin", [0.0, 0.0, 0.0]), dtype=float)
    ox, oy, oz = (
        origin * 1000.0
    )  # km → m (matches MATLAB: 1000*Cond.grid.origin)

    dx = model.x_widths
    dy = model.y_widths
    dz = model.z_widths[n_air:]  # earth layers only

    # -- mesh file -----------------------------------------------------------
    with msh_path.open("w") as fh:
        fh.write(f"{nx} {ny} {nz_earth}\n")
        fh.write(f"{ox:G} {oy:G} {oz:G}\n\n")

        for widths in (dx, dy, dz):
            for j in range(0, len(widths), _WIDTHS_PER_LINE):
                row = widths[j : j + _WIDTHS_PER_LINE]
                fh.write(" ".join(f"{v:G}" for v in row) + "\n")
            fh.write("\n\n")

    # -- conductivity file ---------------------------------------------------
    sigma = np.exp(-model.rho_loge[n_air:])  # ln(ρ) → σ = 1/ρ
    # Write order: j (ny) outer, i (nx) middle, k (nzEarth) inner
    with con_path.open("w") as fh:
        for j in range(ny):
            for i in range(nx):
                for k in range(nz_earth):
                    fh.write(f"{sigma[k, j, i]:13.5E}\n")

    return msh_path, con_path
