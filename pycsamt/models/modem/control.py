# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Read and write ModEM inversion-control files.

The module defines :class:`ModEmControl`, the container for
the plain text ``.inv`` file used to configure ModEM nonlinear
conjugate-gradient inversions.
"""

from __future__ import annotations

from pathlib import Path

from .base import ModEmBase
from .config import ModEmConfig
from .doc import _modem_param_docs as _params

PathLike = str | Path

__all__ = ["ModEmControl"]

_KEYS = [
    ("output_stem",       "Model and data output file name"),
    ("initial_lambda",    "Initial damping factor lambda"),
    ("lambda_divisor",    "To update lambda divide by"),
    ("initial_alpha",     "Initial search step in model units"),
    ("rms_diff_tol",      "Restart when rms diff is less than"),
    ("target_rms",        "Exit search when rms is less than"),
    ("lambda_exit",       "Exit when lambda is less than"),
    ("max_iterations",    "Maximum number of iterations"),
]
_FLOAT_ATTRS = {
    "initial_lambda",
    "lambda_divisor",
    "initial_alpha",
    "rms_diff_tol",
    "target_rms",
    "lambda_exit",
}
_INT_ATTRS = {"max_iterations"}


class ModEmControl(ModEmBase):
    def __init__(self, config: ModEmConfig | None = None, **kwargs):
        """Initialize control values from a configuration."""
        super().__init__(**kwargs)
        cfg = config or ModEmConfig()
        self.config: ModEmConfig = cfg
        self.output_stem: str = cfg.output_stem
        self.initial_lambda: float = cfg.initial_lambda
        self.lambda_divisor: float = cfg.lambda_divisor
        self.initial_alpha: float = cfg.initial_alpha
        self.rms_diff_tol: float = cfg.rms_diff_tol
        self.target_rms: float = cfg.target_rms
        self.lambda_exit: float = cfg.lambda_exit
        self.max_iterations: int = cfg.max_iterations

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------
    @classmethod
    def from_config(
        cls,
        config: ModEmConfig | None = None,
        **kwargs,
    ) -> ModEmControl:
        """Build a control object from a :class:`ModEmConfig`.

        Parameters
        ----------
        config : ModEmConfig, optional
            Configuration object supplying output stem, lambda
            controls, target RMS, convergence thresholds, and
            maximum iteration count. If omitted, a default
            :class:`ModEmConfig` is used.
        **kwargs : dict
            Additional keyword arguments forwarded to
            :class:`ModEmControl`, commonly ``verbose`` or
            ``logger`` inherited from :class:`ModEmBase`.

        Returns
        -------
        ModEmControl
            Control-file container initialized from ``config``.

        Examples
        --------
        >>> from pycsamt.models.modem.config import ModEmConfig
        >>> from pycsamt.models.modem.control import ModEmControl
        >>> cfg = ModEmConfig(max_iterations=50, target_rms=1.03)
        >>> ctrl = ModEmControl.from_config(cfg)
        >>> ctrl.max_iterations
        50
        """
        return cls(config=config, **kwargs)

    # ------------------------------------------------------------------
    # I/O
    # ------------------------------------------------------------------
    @classmethod
    def read(cls, path: PathLike, **kwargs) -> ModEmControl:
        """Parse an existing ModEM inversion-control file.

        The reader scans the key-value labels recognized by
        ModEM and maps them to object attributes. Unknown lines
        are ignored so comments or executable-specific additions
        do not prevent the standard fields from loading.

        Parameters
        ----------
        path : path-like
            Path to an existing ModEM ``.inv`` control file.
            The file must contain colon-separated key-value
            rows such as ``Initial damping factor lambda`` and
            ``Maximum number of iterations``.
        **kwargs : dict
            Additional keyword arguments forwarded to
            :class:`ModEmControl`.

        Returns
        -------
        ModEmControl
            Parsed control object. Fields not found in the file
            retain their default values from ``ModEmConfig``.

        Raises
        ------
        FileNotFoundError
            If ``path`` does not exist.

        Examples
        --------
        >>> from pycsamt.models.modem.control import ModEmControl
        >>> ctrl = ModEmControl.read("ModEM.inv")
        >>> ctrl.target_rms > 0
        True
        """
        p = Path(path)
        if not p.exists():
            msg = f"ModEM control file not found: {p}"
            raise FileNotFoundError(msg)

        obj = cls(**kwargs)
        with p.open("r", errors="replace") as fh:
            for raw in fh:
                if ":" not in raw:
                    continue
                key_raw, _, val_raw = raw.partition(":")
                key = key_raw.strip().upper()
                val = val_raw.strip()
                for attr, label in _KEYS:
                    if label.upper() == key:
                        if attr in _FLOAT_ATTRS:
                            try:
                                setattr(obj, attr, float(val))
                            except ValueError:
                                pass
                        elif attr in _INT_ATTRS:
                            try:
                                setattr(obj, attr, int(float(val)))
                            except ValueError:
                                pass
                        else:
                            setattr(obj, attr, val)
                        break
        if obj.verbose:
            obj.logger.info("ModEmControl.read: loaded from %s", p)
        return obj

    def write(self, path: PathLike) -> Path:
        """Write the control object to a ModEM ``.inv`` file.

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
        Floating-point values are written with compact
        ``%.4g`` formatting. This matches the simple
        key-value style expected by ModEM while keeping files
        readable in version control.

        Examples
        --------
        >>> from pycsamt.models.modem.config import ModEmConfig
        >>> from pycsamt.models.modem.control import ModEmControl
        >>> ctrl = ModEmControl.from_config(ModEmConfig())
        >>> path = ctrl.write("control.inv")
        >>> path.name
        'control.inv'
        """

        p = Path(path)
        p.parent.mkdir(parents=True, exist_ok=True)

        _W = 44   # key field width
        lines: list[str] = []
        for attr, label in _KEYS:
            val = getattr(self, attr)
            key = f"{label:<{_W}}"
            if attr in _FLOAT_ATTRS:
                lines.append(f"{key}: {val:.4g}\n")
            elif attr in _INT_ATTRS:
                lines.append(f"{key}: {int(val)}\n")
            else:
                lines.append(f"{key}: {val}\n")

        with p.open("w") as fh:
            fh.writelines(lines)
        return p


ModEmControl.__doc__ = rf"""
Represent a ModEM inversion-control file.

``ModEmControl`` stores the small key-value ``.inv`` file used
by ModEM to control nonlinear inversion. The same container is
used for 2-D and 3-D runs. It records the output file stem,
lambda damping parameters, line-search starting step,
convergence thresholds, target misfit, and maximum number of
iterations. The object can be created directly from
:class:`ModEmConfig`, read from an existing control file, or
written as part of an :class:`InputBuilder` workflow.

The control parameters guide the regularized nonlinear solve.
A typical ModEM objective can be written as

.. math::

   \Phi(m) =
   \| W_d (F(m) - d) \|_2^2 +
   \lambda \| W_m (m - m_0) \|_2^2 ,

where :math:`m` is the current model, :math:`F(m)` is the
forward response, :math:`d` is the observed data vector,
:math:`W_d` and :math:`W_m` are weighting operators, and
:math:`\lambda` is the damping factor written to this file.
The target RMS controls when the data-fit part is considered
adequate for the assigned uncertainties.

Parameters
----------
{_params.common.config}
{_params.common.verbose}
{_params.common.logger}

Attributes
----------
config : ModEmConfig
    Configuration object used to initialize the control
    values. It is retained so callers can inspect the source
    settings used to create the control file.
{_params.control.output_stem}
{_params.control.initial_lambda}
{_params.control.lambda_divisor}
{_params.control.initial_alpha}
{_params.control.rms_diff_tol}
{_params.control.target_rms}
{_params.control.lambda_exit}
{_params.control.max_iterations}

Notes
-----
The ModEM control file is a plain text file with one
colon-separated key-value pair per line. The canonical labels
written by this class are:

* ``Model and data output file name``;
* ``Initial damping factor lambda``;
* ``To update lambda divide by``;
* ``Initial search step in model units``;
* ``Restart when rms diff is less than``;
* ``Exit search when rms is less than``;
* ``Exit when lambda is less than``;
* ``Maximum number of iterations``.

The reader is intentionally tolerant: fields that cannot be
parsed retain their current defaults, and unrecognized lines
are skipped. This allows the loader to handle control files
with comments or minor executable-specific additions.

See Also
--------
ModEmConfig
    Supplies the inversion-control values used here.
InputBuilder
    Writes a control file as part of a complete ModEM input
    set.
ModEmRunner
    Passes the written control file to the ModEM executable.
InversionResult
    Loads the control file found in a completed run directory.

Examples
--------
Create a control file from configuration values:

>>> from pycsamt.models.modem.config import ModEmConfig
>>> from pycsamt.models.modem.control import ModEmControl
>>> cfg = ModEmConfig(max_iterations=80, target_rms=1.05)
>>> ctrl = ModEmControl.from_config(cfg)
>>> ctrl.max_iterations
80

Write and read a control file:

>>> path = ctrl.write("ModEM.inv")
>>> loaded = ModEmControl.read(path)
>>> loaded.target_rms == ctrl.target_rms
True

Tune lambda controls before writing:

>>> ctrl.initial_lambda = 20.0
>>> ctrl.lambda_divisor = 50.0
>>> ctrl.write("strong_start.inv")
PosixPath('strong_start.inv')

References
----------
.. [1] Egbert, G. D., and Kelbert, A., "Computational
   recipes for electromagnetic inverse problems", Geophysical
   Journal International, 189(1), 251-267, 2012,
   doi:10.1111/j.1365-246X.2011.05347.x.
.. [2] Kelbert, A., Meqbel, N., Egbert, G. D., and Tandon,
   K., "ModEM: A modular system for inversion of
   electromagnetic geophysical data", Computers and
   Geosciences, 66, 40-53, 2014,
   doi:10.1016/j.cageo.2014.01.010.
"""
