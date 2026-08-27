# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Write ModEM 3-D forward-solver control files.

The module defines :class:`ModEmForwardControl`, the container
for the plain text forward-solver control file (``rFile_fwdCtrl``
in ModEM's own ``-I NLCG`` usage text) read by Mod3DMT's
``readEMsolveControl`` (``3D_MT/FWD/EMsolve3D.f90``).

This file has one purpose in this project: Mod3DMT's ``-I NLCG``
argument order is ``rFile_Model rFile_Data [rFile_invCtrl
rFile_fwdCtrl] [rFile_Cov]`` -- the covariance file is the *sixth*
positional argument and is only reached once a forward-control file
(the fifth) is also supplied and exists on disk. Without one, a
custom covariance file cannot be passed to Mod3DMT at all (it
silently falls back to its own default covariance instead, or -- if
a covariance path is forced into the wrong slot -- misparses it as
forward-solver settings and aborts with a Fortran runtime error).

:meth:`ModEmForwardControl.write` writes exactly the six required
lines, using the same numeric values compiled into Mod3DMT as its
own defaults (``IterPerDivCorDef``, ``MaxDivCorDef``,
``MaxIterDivCorDef``, ``tolEMDef``, ``tolDivCorDef`` in
``EMsolve3D.f90``), and deliberately omits the two optional
sections (a nested-boundary-condition EM-solution file, and an
explicit air-layers override) -- so that supplying this file changes
nothing about forward-solver behaviour versus omitting it entirely;
its only effect is unlocking the covariance argument slot.
"""

from __future__ import annotations

from pathlib import Path
from typing import Union

from .base import ModEmBase
from .config import ModEmConfig
from .doc import _modem_param_docs as _params

PathLike = Union[str, Path]

__all__ = ["ModEmForwardControl"]

# (attribute, label, is_integer) -- label width and format below match
# EMsolve3D.f90's readEMsolveControl: `read (unit,'(a48,i5)') string,int_val`
# or `read (unit,'(a48,g15.7)') string,float_val`.
_KEYS = [
    ("qmr_iters_per_divcor", "Number of QMR iters per divergence correction", True),
    ("max_divcor", "Maximum number of divergence correction calls", True),
    ("max_iter_divcor", "Maximum number of divergence correction iters", True),
    ("tol_em_fwd", "Misfit tolerance for EM forward solver", False),
    ("tol_em_adj", "Misfit tolerance for EM adjoint solver", False),
    ("tol_divcor", "Misfit tolerance for divergence correction", False),
]


class ModEmForwardControl(ModEmBase):
    def __init__(self, config: ModEmConfig | None = None, **kwargs):
        """Initialize forward-solver control values from a configuration."""
        super().__init__(**kwargs)
        cfg = config or ModEmConfig()
        self.config: ModEmConfig = cfg
        self.qmr_iters_per_divcor: int = cfg.qmr_iters_per_divcor
        self.max_divcor: int = cfg.max_divcor
        self.max_iter_divcor: int = cfg.max_iter_divcor
        self.tol_em_fwd: float = cfg.tol_em_fwd
        self.tol_em_adj: float = cfg.tol_em_adj
        self.tol_divcor: float = cfg.tol_divcor

    @classmethod
    def from_config(
        cls,
        config: ModEmConfig | None = None,
        **kwargs,
    ) -> ModEmForwardControl:
        """Build a forward-control object from a :class:`ModEmConfig`.

        Parameters
        ----------
        config : ModEmConfig, optional
            Configuration object supplying the QMR iteration count
            and solver tolerances. If omitted, a default
            :class:`ModEmConfig` is used -- which reproduces
            Mod3DMT's own compiled-in defaults exactly.
        **kwargs : dict
            Additional keyword arguments forwarded to
            :class:`ModEmForwardControl`, commonly ``verbose`` or
            ``logger`` inherited from :class:`ModEmBase`.

        Returns
        -------
        ModEmForwardControl

        Examples
        --------
        >>> from pycsamt.models.modem.config import ModEmConfig
        >>> from pycsamt.models.modem.forward_control import (
        ...     ModEmForwardControl,
        ... )
        >>> ctrl = ModEmForwardControl.from_config(ModEmConfig())
        >>> ctrl.qmr_iters_per_divcor
        40
        """
        return cls(config=config, **kwargs)

    def write(self, path: PathLike) -> Path:
        """Write the six required lines to a ModEM forward-control file.

        Parameters
        ----------
        path : path-like
            Destination file. Parent directories are created
            before writing. Existing files are overwritten.

        Returns
        -------
        pathlib.Path
            Path passed to the writer, converted to
            :class:`pathlib.Path`.

        Notes
        -----
        ``EMsolve3D.f90``'s ``readEMsolveControl`` parses this file
        with fixed column widths, not by splitting on ``:`` -- each
        label is read as ``a48`` (columns 1-48) and the value
        immediately after as ``i5`` (integers) or ``g15.7``
        (floats). A label field wider than 48 columns would push
        the value out of alignment, the same class of bug already
        fixed for :class:`~pycsamt.models.modem.control.ModEmControl`
        (``a36`` there, not ``a48`` -- the two file formats use
        different fixed widths; do not share the constant).

        Float values are written in ``%.6E`` scientific notation, not
        Python's default ``%g``: Fortran's ``G`` edit descriptor on
        *input* requires an explicit decimal point, or the field's own
        decimal-digit count (``.7`` here) silently re-places one --
        confirmed both by a real run (a written ``1e-07`` was read back
        as ``0.1000000E-13``) and by ModEM's own usage-text examples in
        ``UserCtrl.f90``, which write even whole numbers with a trailing
        ``.`` (``"1.0e-7"``, never ``"1e-7"``).

        Examples
        --------
        >>> from pycsamt.models.modem.config import ModEmConfig
        >>> from pycsamt.models.modem.forward_control import (
        ...     ModEmForwardControl,
        ... )
        >>> ctrl = ModEmForwardControl.from_config(ModEmConfig())
        >>> path = ctrl.write("ModEM_fwd.ctrl")
        >>> path.name
        'ModEM_fwd.ctrl'
        """
        p = Path(path)
        p.parent.mkdir(parents=True, exist_ok=True)

        _W = 48  # key field width -- matches Fortran a48 (EMsolve3D.f90)
        lines: list[str] = []
        for attr, label, is_int in _KEYS:
            val = getattr(self, attr)
            key = f"{label + ':':<{_W}}"
            if is_int:
                lines.append(f"{key}{int(val)}\n")
            else:
                lines.append(f"{key}{val:.6E}\n")

        with p.open("w") as fh:
            fh.writelines(lines)
        return p


ModEmForwardControl.__doc__ = rf"""
Represent a ModEM 3-D forward-solver control file.

``ModEmForwardControl`` stores the small key-value file
(``rFile_fwdCtrl``) that configures Mod3DMT's forward/adjoint EM
solver: QMR iterations per divergence correction, divergence
correction limits, and solver tolerances. Written with the same
values Mod3DMT already uses by default, so its only functional
effect in this project is occupying the fifth positional CLI
argument -- required by Mod3DMT's own ``-I NLCG`` argument order
before a sixth argument (the covariance file) can be supplied at
all. See :class:`~pycsamt.models.modem.runner.ModEmRunner` for how
the two are passed together.

Parameters
----------
{_params.common.config}
{_params.common.verbose}
{_params.common.logger}

Attributes
----------
config : ModEmConfig
    Configuration object used to initialize the control values.
{_params.forward.qmr_iters_per_divcor}
{_params.forward.max_divcor}
{_params.forward.max_iter_divcor}
{_params.forward.tol_em_fwd}
{_params.forward.tol_em_adj}
{_params.forward.tol_divcor}

Notes
-----
Two optional sections of the real file format are deliberately not
written: a nested-boundary-condition EM-solution filename, and an
explicit air-layers override (mirror / fixed height / read from
file). Both are read with a ``advance='no'``-then-EOF pattern in
``readEMsolveControl`` that fails gracefully when absent, and
Mod3DMT's own pre-existing air-layer initialization (confirmed via
a real run: ``Air layers setup complete according to the method :
mirror``) is left untouched by their absence.

See Also
--------
ModEmConfig
    Supplies the forward-solver values used here.
ModEmControl
    The sibling inversion-control (``.inv``) file -- a different
    fixed column width (``a36``, not ``a48``).
ModEmRunner
    Passes both this file and the covariance file to Mod3DMT.

Examples
--------
Create and write a forward-control file matching Mod3DMT's own
defaults exactly:

>>> from pycsamt.models.modem.config import ModEmConfig
>>> from pycsamt.models.modem.forward_control import ModEmForwardControl
>>> ctrl = ModEmForwardControl.from_config(ModEmConfig())
>>> path = ctrl.write("ModEM_fwd.ctrl")

References
----------
.. [ModEmForwardControl-1] Kelbert, A., Meqbel, N., Egbert, G. D., and
   Tandon, K., "ModEM: A modular system for inversion of
   electromagnetic geophysical data", Computers and
   Geosciences, 66, 40-53, 2014,
   doi:10.1016/j.cageo.2014.01.010.
"""
