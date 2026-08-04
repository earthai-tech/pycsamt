# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Subprocess wrapper for J. R. Shewchuk's Triangle mesh generator.

MARE2DEM's own source build compiles Triangle alongside the main
solver (see ``TRICOPTS`` in
:func:`pycsamt.models.mare2dem.source._generate_inc`), so a
``triangle``/``Triangle`` executable is normally a byproduct of
:meth:`~pycsamt.models.mare2dem.source.SourceManager.build`. Meshing
itself has no MPI/Fortran/LAPACK dependency, so this module never
launches the (much heavier) MARE2DEM solver -- only the standalone
Triangle refinement tool, resolved independently via
:meth:`~pycsamt.models.mare2dem.source.SourceManager.resolve_triangle_binary`.

This module writes and runs a subprocess only; it does not parse
Triangle's ``.node``/``.ele`` output itself -- see
:func:`pycsamt.models.mare2dem.iotools.poly.read_triangulation` for that.
"""

from __future__ import annotations

import subprocess
from pathlib import Path
from typing import TYPE_CHECKING, Sequence, Union

from .source import SourceManager

if TYPE_CHECKING:
    from .config import Mare2DEMConfig

PathLike = Union[str, Path]

__all__ = ["TriangleRunError", "run_triangle"]


class TriangleRunError(RuntimeError):
    """Raised when Triangle cannot be found, fails, or produces no output.

    Examples
    --------
    >>> isinstance(TriangleRunError("no triangle"), RuntimeError)
    True
    """


def run_triangle(
    poly_path: PathLike,
    *,
    executable: PathLike | Sequence[PathLike] | None = None,
    min_angle: float | None = 30.0,
    max_area: float | None = None,
    extra_args: Sequence[str] = (),
    config: "Mare2DEMConfig | None" = None,
    timeout: float | None = None,
) -> Path:
    """Refine a Triangle ``.poly`` PSLG into a quality FEM mesh.

    Parameters
    ----------
    poly_path : path-like
        Input ``.poly`` PSLG file, e.g. from
        :func:`~pycsamt.models.mare2dem.iotools.poly.write_poly`.
    executable : path-like or sequence of path-like, optional
        Explicit Triangle executable, or an ``[interpreter, script]``
        sequence (e.g. ``[sys.executable, "fake_triangle.py"]``) for a
        scripted test double standing in for a real binary. Resolved via
        :meth:`~pycsamt.models.mare2dem.source.SourceManager.resolve_triangle_binary`
        when omitted.
    min_angle : float or None, default=30.0
        Minimum triangle angle in degrees, Triangle's ``-q`` quality
        constraint. Pass ``None`` to disable it.
    max_area : float or None, optional
        Maximum triangle area, Triangle's ``-a`` constraint.
    extra_args : sequence of str, optional
        Additional flags appended verbatim (e.g. ``("-Y",)`` to forbid
        Steiner points on PSLG boundary segments).
    config : Mare2DEMConfig, optional
        Forwarded to :class:`~pycsamt.models.mare2dem.source.SourceManager`
        when resolving the executable.
    timeout : float, optional
        Subprocess timeout in seconds.

    Returns
    -------
    pathlib.Path
        Path to the refined ``<stem>.1.poly`` file -- Triangle's own
        output-naming convention when run with ``-p`` on ``<stem>.poly``.
        Its companion ``<stem>.1.node``/``<stem>.1.ele`` files sit
        alongside it; read all three together with
        :func:`~pycsamt.models.mare2dem.iotools.poly.read_triangulation`.

    Raises
    ------
    FileNotFoundError
        If ``poly_path`` does not exist.
    TriangleRunError
        If no executable is found, the process exits non-zero, or the
        expected ``.1.poly`` output is missing afterward.

    Examples
    --------
    >>> run_triangle("missing.poly")  # doctest: +SKIP
    Traceback (most recent call last):
    ...
    FileNotFoundError: poly_path not found: missing.poly
    """
    poly = Path(poly_path)
    if not poly.exists():
        raise FileNotFoundError(f"poly_path not found: {poly}")

    if executable is None:
        resolved = SourceManager(config=config).resolve_triangle_binary()
        exe_args: list[str] = [] if resolved is None else [str(resolved)]
    elif isinstance(executable, (list, tuple)):
        # A [interpreter, script] pair, e.g. for a scripted test double
        # standing in for a real Triangle binary -- a single Path cannot
        # express "run this script with this interpreter".
        exe_args = [str(part) for part in executable]
    else:
        exe_args = [str(executable)]
    if not exe_args:
        raise TriangleRunError(
            "no Triangle executable found. Build MARE2DEM's source tree "
            "(SourceManager.build(), which compiles Triangle alongside "
            "the solver) or install a 'triangle'/'Triangle' executable "
            "on PATH, then retry."
        )

    flags = "p"
    if min_angle is not None:
        flags += f"q{min_angle:g}"
    if max_area is not None:
        flags += f"a{max_area:g}"
    flags += "ze"
    args = [*exe_args, f"-{flags}", *extra_args, str(poly.name)]

    try:
        completed = subprocess.run(
            args,
            cwd=poly.parent if str(poly.parent) else ".",
            capture_output=True,
            text=True,
            timeout=timeout,
        )
    except OSError as exc:
        raise TriangleRunError(f"failed to launch Triangle: {exc}") from exc

    if completed.returncode != 0:
        raise TriangleRunError(
            f"Triangle exited {completed.returncode} for {poly}: "
            f"{completed.stdout[-2000:]}\n{completed.stderr[-2000:]}"
        )

    refined = poly.with_name(f"{poly.stem}.1.poly")
    if not refined.exists():
        raise TriangleRunError(
            f"Triangle reported success (exit 0) but {refined} was not "
            "produced -- check its stdout/stderr for warnings:\n"
            f"{completed.stdout[-2000:]}\n{completed.stderr[-2000:]}"
        )
    return refined
