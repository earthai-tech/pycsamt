# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""OccamStartup and OccamIter — Startup / .iter file handling.

Both file types share the ``OCCAMITER_FLEX`` format header.  The
distinction is that a Startup file has ``Iteration: 0`` and a uniform
initial model vector, while an ``.iter`` file carries a real iteration
number and the model parameters from that iteration.

OCCAMITER_FLEX format
---------------------
The file is split into two sections separated by ``Param Count: N``:

Header (keyword: value pairs, one per line)::

    Format:             OCCAMITER_FLEX
    Description:        <text>
    Model File:         <filename>
    Data File:          <filename>
    Date/Time:          <text>
    Iterations to run:  <int>
    Target Misfit:      <float>
    Roughness Type:     <int>
    Diagonal Penalties: <int>
    Stepsize Cut Count: <int>
    !Model Limits:      <text>       ← comment lines (!) are skipped
    !Model Value Steps: <text>
    Debug Level:        <int>
    Iteration:          <int>        ← 0 = Startup, N ≥ 1 = iter file
    Lagrange Value:     <float>
    Roughness Value:    <float>
    Misfit Value:       <float>
    Misfit Reached:     <int>        ← 0 or 1
    Param Count:        <int>        ← triggers switch to param block

Parameter block (space-separated floats, 4 per line)::

    log10(ρ₀)  log10(ρ₁)  log10(ρ₂)  log10(ρ₃)
    ...

Classes
-------
OccamStartup
    Reads / writes / builds the initial Startup file (Iteration: 0).
OccamIter
    Reads any ``.iter`` file (Iteration: N ≥ 1) produced by Occam.
"""

from __future__ import annotations

import datetime
from pathlib import Path
from typing import Union

import numpy as np

from .base import OccamBase
from .config import OccamConfig

PathLike = Union[str, Path]

__all__ = ["OccamStartup", "OccamIter"]

_FORMAT_TAG = "OCCAMITER_FLEX"

# -----------------------------------------------------------------------
# Shared low-level parser
# -----------------------------------------------------------------------

# Canonical key map: upper-cased header keyword → attribute name
_KEY_MAP: dict[str, str] = {
    "FORMAT": "format_str",
    "DESCRIPTION": "description",
    "MODEL FILE": "model_file",
    "DATA FILE": "data_file",
    "DATE/TIME": "datetime_str",
    "ITERATIONS TO RUN": "max_iterations",
    "TARGET MISFIT": "target_misfit",
    "ROUGHNESS TYPE": "roughness_type",
    "DIAGONAL PENALTIES": "diagonal_penalties",
    "STEPSIZE CUT COUNT": "stepsize_cut_count",
    "DEBUG LEVEL": "debug_level",
    "ITERATION": "iteration",
    "LAGRANGE VALUE": "lagrange_value",
    "ROUGHNESS VALUE": "roughness_value",
    "MISFIT VALUE": "misfit_value",
    "MISFIT REACHED": "misfit_reached",
    "PARAM COUNT": "n_params",
}

# Fields that need integer coercion
_INT_FIELDS = {
    "max_iterations",
    "roughness_type",
    "diagonal_penalties",
    "stepsize_cut_count",
    "debug_level",
    "iteration",
    "n_params",
}


def _parse_iter_flex(path: Path) -> dict:
    """Parse an OCCAMITER_FLEX file and return a plain dict.

    Parameters
    ----------
    path : Path
        Resolved path to the file.

    Returns
    -------
    dict
        Keys matching ``_KEY_MAP`` values plus
        ``"param_values"`` as an array of log10 resistivity.

    Raises
    ------
    FileNotFoundError
        If *path* does not exist.
    ValueError
        If the file lacks a valid OCCAMITER_FLEX header.
    """
    if not path.exists():
        raise FileNotFoundError(f"OCCAMITER_FLEX file not found: {path}")

    result: dict = {v: None for v in _KEY_MAP.values()}
    param_tokens: list[str] = []
    reading_params = False
    n_params = 0

    with path.open("r", errors="replace") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("!"):
                continue

            if reading_params:
                param_tokens.extend(line.split())
                continue

            # Header key: value
            if ":" not in line:
                continue
            raw_key, _, raw_val = line.partition(":")
            key_up = raw_key.strip().upper()
            val = raw_val.strip()

            attr = _KEY_MAP.get(key_up)
            if attr is None:
                continue

            if attr in _INT_FIELDS:
                try:
                    result[attr] = int(float(val))
                except ValueError:
                    result[attr] = 0
            elif attr == "misfit_reached":
                try:
                    result[attr] = bool(int(float(val)))
                except ValueError:
                    result[attr] = False
            elif attr == "format_str":
                result[attr] = val
                if val.upper() != _FORMAT_TAG:
                    raise ValueError(
                        f"Expected format '{_FORMAT_TAG}', got '{val}' in {path}"
                    )
            else:
                # float fields and string fields
                try:
                    result[attr] = float(val)
                except ValueError:
                    result[attr] = val  # store as string (description, filenames…)

            if attr == "n_params":
                n_params = result["n_params"] or 0
                reading_params = True

    if result["format_str"] is None:
        raise ValueError(f"File does not contain a valid OCCAMITER_FLEX header: {path}")

    # Parse parameter vector
    try:
        arr = np.array([float(t) for t in param_tokens[:n_params]], dtype=float)
    except ValueError:
        arr = np.array([], dtype=float)

    if arr.size < n_params and n_params > 0:
        # Pad with NaN if fewer tokens than declared (truncated file)
        pad = np.full(n_params - arr.size, np.nan)
        arr = np.concatenate([arr, pad])

    result["param_values"] = arr
    return result


# -----------------------------------------------------------------------
# OccamStartup
# -----------------------------------------------------------------------


class OccamStartup(OccamBase):
    r"""Represent an Occam2D startup control file.

    ``OccamStartup`` stores the iteration-zero
    ``OCCAMITER_FLEX`` file passed to the Occam2D executable.
    It defines run controls, file references, inversion
    options, and the initial model vector. Unlike ``.iter``
    files produced by the solver, a valid startup file has
    ``Iteration: 0``.

    The startup parameter vector is initialized as a uniform
    half-space:

    .. math::

        m_i = \log_{10}(\rho_0),
        \qquad i = 1, \ldots, N_p.

    Here :math:`\rho_0` is ``config.initial_rho`` and
    :math:`N_p` is the number of model parameters defined by
    :class:`OccamModel`.

    Parameters
    ----------
    config : OccamConfig, optional
        Configuration object providing file names, inversion
        controls, starting resistivity, target misfit,
        roughness settings, and debug level. If omitted, a
        default :class:`OccamConfig` is created.
    description : str, default "startup created by pycsamt"
        Description written to the ``Description`` header.
        Use it to record the purpose or provenance of the run.
    verbose : int or bool, default 0
        Verbosity level inherited from :class:`OccamBase`.
        Positive values enable progress messages through the
        instance logger.
    logger : logging.Logger, optional
        Logger used for progress and diagnostic messages. If
        omitted, a class-specific PyCSAMT logger is created.

    Attributes
    ----------
    format_str : str
        File format tag, always ``"OCCAMITER_FLEX"``.
    description : str
        Human-readable startup description.
    model_file : str
        Model file name referenced by the startup file.
    data_file : str
        Data file name referenced by the startup file.
    datetime_str : str
        Creation or file timestamp string.
    max_iterations : int
        Maximum number of iterations requested from Occam2D.
    target_misfit : float
        Target normalized RMS misfit.
    roughness_type : int
        Roughness penalty type written to the startup file.
    diagonal_penalties : int
        Flag controlling diagonal roughness penalties.
    stepsize_cut_count : int
        Maximum number of step-size cuts in a line search.
    debug_level : int
        Debug verbosity passed to the Fortran executable.
    iteration : int
        Always 0 for a valid Startup file.
    lagrange_value : float
        Initial Lagrange multiplier written as
        ``Lagrange Value``.
    roughness_value : float
        Initial roughness value written before the first run.
    misfit_value : float
        Initial misfit value written before the first run.
    misfit_reached : bool
        Whether the target misfit has already been reached.
        Startup files normally use ``False``.
    n_params : int
        Number of model parameters. This must match the model
        file and the length of :attr:`param_values`.
    param_values : numpy.ndarray of float, shape (n_params,)
        Initial log10-resistivity values. Values are uniform
        after :meth:`from_model` and equal to
        ``log10(config.initial_rho)``.

    Notes
    -----
    ``OccamStartup`` writes the same flexible iteration format
    that Occam later uses for ``.iter`` files. The distinction
    is semantic: startup files carry ``Iteration: 0`` and are
    input to the solver, while ``OccamIter`` files carry
    non-zero iteration numbers and are output from the solver.

    See Also
    --------
    OccamModel
        Provides the parameter count for the startup vector.
    OccamIter
        Reads iteration files produced after running Occam2D.
    OccamRunner
        Launches the executable with this startup file.

    Examples
    --------
    Build a startup file from a model definition:

    >>> from pycsamt.models.occam2d import OccamModel
    >>> from pycsamt.models.occam2d import OccamStartup
    >>> model = OccamModel.read("occam_run/Occam2DModel")
    >>> startup = OccamStartup.from_model(model)
    >>> startup.write("occam_run/Startup")

    Read an existing startup file:

    >>> from pycsamt.models.occam2d import OccamStartup
    >>> startup = OccamStartup.read("occam_run/Startup")
    >>> startup.n_params

    References
    ----------
    .. [1] Constable, S. C., Parker, R. L., and Constable,
       C. G., "Occam's inversion: A practical algorithm for
       generating smooth models from electromagnetic sounding
       data", Geophysics, 52(3), 289-300, 1987.
    .. [2] deGroot-Hedlin, C., and Constable, S.,
       "Occam's inversion to generate smooth, two-dimensional
       models from magnetotelluric data", Geophysics, 55(12),
       1613-1624, 1990.
    """

    def __init__(
        self,
        config: OccamConfig | None = None,
        description: str = "startup created by pycsamt",
        **kwargs,
    ):
        super().__init__(**kwargs)
        cfg = config or OccamConfig()
        self.config = cfg
        self.format_str = _FORMAT_TAG
        self.description = description
        self.model_file = cfg.model_file
        self.data_file = cfg.data_file
        self.datetime_str = str(datetime.datetime.now())
        self.max_iterations = cfg.max_iterations
        self.target_misfit = cfg.target_misfit
        self.roughness_type = cfg.roughness_type
        self.diagonal_penalties = cfg.diagonal_penalties
        self.stepsize_cut_count = cfg.stepsize_cut_count
        self.debug_level = cfg.debug_level
        self.iteration = 0
        self.lagrange_value = cfg.lagrange_start
        self.roughness_value = 1e10
        self.misfit_value = 1000.0
        self.misfit_reached = False
        self.n_params = 0
        self.param_values: np.ndarray = np.array([])

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------
    @classmethod
    def from_model(
        cls,
        model: OccamModel,  # noqa: F821
        config: OccamConfig | None = None,
        **kwargs,
    ) -> OccamStartup:
        r"""Build a startup object from a model definition.

        The method uses ``model.n_params`` to size the startup
        vector and fills every entry with the log10 value of
        starting half-space resistivity:

        .. math::

            m_i = \log_{10}(\rho_0),
            \qquad i = 1,\ldots,N_p.

        This produces the standard smooth-inversion initial
        model: a homogeneous half-space whose value is later
        updated by the Occam solver.

        Parameters
        ----------
        model : OccamModel
            Populated model-definition object. It must contain
            at least one parameter so the vector can be
            sized consistently with the ``Occam2DModel`` file.
        config : OccamConfig, optional
            Configuration object providing ``initial_rho``,
            model and data file names, iteration controls, and
            inversion settings. If omitted, a default
            :class:`OccamConfig` is created.
        **kwargs
            Additional keyword arguments forwarded to the
            ``OccamStartup`` constructor. Use this for
            ``description``, ``verbose``, or ``logger``.

        Returns
        -------
        OccamStartup
            Startup object with ``n_params`` and uniform
            ``param_values`` populated.

        Raises
        ------
        ValueError
            Raised when ``model.n_params`` is not positive.

        See Also
        --------
        OccamStartup.write
            Serializes the generated startup object.
        OccamModel
            Supplies the parameter count used here.

        Examples
        --------
        >>> from pycsamt.models.occam2d import OccamModel
        >>> from pycsamt.models.occam2d import OccamStartup
        >>> model = OccamModel.read("occam_run/Occam2DModel")
        >>> startup = OccamStartup.from_model(model)
        >>> startup.param_values.shape
        """
        cfg = config or OccamConfig()
        obj = cls(config=cfg, **kwargs)

        n = model.n_params
        if n <= 0:
            raise ValueError("OccamStartup.from_model: model has no parameters")

        obj.n_params = n
        obj.param_values = np.full(n, np.log10(cfg.initial_rho), dtype=float)
        obj.model_file = cfg.model_file
        obj.data_file = cfg.data_file

        if obj.verbose:
            obj.logger.info(
                "OccamStartup.from_model: %d params, rho0=%.1f Ω·m",
                n,
                cfg.initial_rho,
            )
        return obj

    # ------------------------------------------------------------------
    # I/O
    # ------------------------------------------------------------------
    @classmethod
    def read(cls, path: PathLike, **kwargs) -> OccamStartup:
        """Read an Occam2D startup file.

        The reader parses an ``OCCAMITER_FLEX`` file and then
        validates that its ``Iteration`` header is zero. Use
        :meth:`OccamIter.read` for non-zero iteration files
        produced by the solver.

        Parameters
        ----------
        path : path-like
            Path to the startup file. The value may be a
            string, :class:`pathlib.Path`, or any object
            accepted by :class:`pathlib.Path`.
        **kwargs
            Additional keyword arguments forwarded to the
            ``OccamStartup`` constructor before parsed values
            are attached.

        Returns
        -------
        OccamStartup
            Parsed startup object with header fields and
            parameter vector populated.

        Raises
        ------
        FileNotFoundError
            Raised when ``path`` does not exist.
        ValueError
            Raised when the file is not ``OCCAMITER_FLEX`` or
            has a non-zero iteration number.

        Examples
        --------
        >>> from pycsamt.models.occam2d import OccamStartup
        >>> startup = OccamStartup.read("occam_run/Startup")
        >>> startup.iteration
        """
        p = Path(path)
        parsed = _parse_iter_flex(p)

        if parsed.get("iteration", 0) != 0:
            raise ValueError(
                f"Expected Iteration: 0 (Startup) but got "
                f"Iteration: {parsed['iteration']} in {p}"
            )

        obj = cls(**kwargs)
        _apply_parsed(obj, parsed)

        if obj.verbose:
            obj.logger.info(
                "OccamStartup.read: %d params loaded from %s", obj.n_params, p
            )
        return obj

    def write(self, path: PathLike) -> Path:
        """Write this startup file in ``OCCAMITER_FLEX`` format.

        Parameters
        ----------
        path : path-like
            Destination path for the startup file. Parent
            directories are created when needed.

        Returns
        -------
        pathlib.Path
            Path to the file that was written.

        See Also
        --------
        OccamStartup.read
            Parses files written by this method.
        OccamRunner
            Passes the written startup file to the executable.
        """
        return _write_iter_flex(self, Path(path))


# -----------------------------------------------------------------------
# OccamIter
# -----------------------------------------------------------------------


class OccamIter(OccamBase):
    r"""Represent an Occam2D iteration file.

    ``OccamIter`` reads ``OCCAMITER_FLEX`` files written by
    the Occam2D executable after one or more inversion
    iterations. These files have the same structural format as
    ``Startup`` but carry ``Iteration`` values greater than
    zero and store the accepted model vector.

    The parameter vector is stored in log10-resistivity units.
    Physical resistivity is recovered as

    .. math::

        \rho_i = 10^{m_i},

    where :math:`m_i` is the stored value for parameter
    :math:`i`.

    Parameters
    ----------
    verbose : int or bool, default 0
        Verbosity level inherited from :class:`OccamBase`.
        Positive values enable progress messages through the
        instance logger.
    logger : logging.Logger, optional
        Logger used for progress and diagnostic messages. If
        omitted, a class-specific PyCSAMT logger is created.

    Attributes
    ----------
    format_str : str
        File format tag, usually ``"OCCAMITER_FLEX"``.
    description : str
        Iteration description written by Occam.
    model_file : str
        Model file referenced by this iteration.
    data_file : str
        Data file referenced by this iteration.
    datetime_str : str
        Date and time string written by Occam.
    max_iterations : int
        Iteration limit stored in the file.
    target_misfit : float
        Target normalized RMS misfit.
    roughness_type : int
        Roughness penalty type.
    diagonal_penalties : int
        Diagonal penalty flag.
    stepsize_cut_count : int
        Maximum number of line-search step-size cuts.
    debug_level : int
        Debug verbosity setting.
    iteration : int
        Iteration number. Valid ``.iter`` files use values
        greater than zero.
    lagrange_value : float
        Lagrange multiplier accepted at this iteration.
    roughness_value : float
        Model roughness value reported by Occam.
    misfit_value : float
        Normalized RMS misfit at this iteration.
    misfit_reached : bool
        ``True`` if the target misfit was achieved.
    n_params : int
        Number of model parameters.
    param_values : numpy.ndarray of float, shape (n_params,)
        Accepted log10-resistivity values for this iteration.

    See Also
    --------
    OccamStartup
        Represents the corresponding iteration-zero file.
    InversionResult
        Selects iteration files from a run directory.
    OccamResponse
        Reads the response file for the same iteration.

    Examples
    --------
    Read an iteration file and convert to resistivity:

    >>> from pycsamt.models.occam2d import OccamIter
    >>> iteration = OccamIter.read("occam_run/ITER17.iter")
    >>> rho = iteration.to_resistivity()

    Inspect log10-resistivity statistics:

    >>> iteration.log10_rho_stats

    References
    ----------
    .. [1] Constable, S. C., Parker, R. L., and Constable,
       C. G., "Occam's inversion: A practical algorithm for
       generating smooth models from electromagnetic sounding
       data", Geophysics, 52(3), 289-300, 1987.
    .. [2] deGroot-Hedlin, C., and Constable, S.,
       "Occam's inversion to generate smooth, two-dimensional
       models from magnetotelluric data", Geophysics, 55(12),
       1613-1624, 1990.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.format_str: str = _FORMAT_TAG
        self.description: str = ""
        self.model_file: str = ""
        self.data_file: str = ""
        self.datetime_str: str = ""
        self.max_iterations: int = 0
        self.target_misfit: float = 1.0
        self.roughness_type: int = 1
        self.diagonal_penalties: int = 0
        self.stepsize_cut_count: int = 8
        self.debug_level: int = 1
        self.iteration: int = 0
        self.lagrange_value: float = 0.0
        self.roughness_value: float = 0.0
        self.misfit_value: float = 0.0
        self.misfit_reached: bool = False
        self.n_params: int = 0
        self.param_values: np.ndarray = np.array([])

    # ------------------------------------------------------------------
    # I/O
    # ------------------------------------------------------------------
    @classmethod
    def read(cls, path: PathLike, **kwargs) -> OccamIter:
        """Read an Occam2D ``.iter`` file.

        The reader parses an ``OCCAMITER_FLEX`` file and
        validates that the ``Iteration`` value is non-zero.
        Startup files should be loaded with
        :meth:`OccamStartup.read`.

        Parameters
        ----------
        path : path-like
            Path to the iteration file. The value may be a
            string, :class:`pathlib.Path`, or any object
            accepted by :class:`pathlib.Path`.
        **kwargs
            Additional keyword arguments forwarded to the
            ``OccamIter`` constructor before parsed values are
            attached.

        Returns
        -------
        OccamIter
            Parsed iteration object with header fields and
            parameter vector populated.

        Raises
        ------
        FileNotFoundError
            Raised when ``path`` does not exist.
        ValueError
            Raised when the file is not ``OCCAMITER_FLEX`` or
            has ``Iteration: 0``.

        Examples
        --------
        >>> from pycsamt.models.occam2d import OccamIter
        >>> iteration = OccamIter.read("ITER17.iter")
        >>> iteration.misfit_value
        """
        p = Path(path)
        parsed = _parse_iter_flex(p)

        if parsed.get("iteration", 0) == 0:
            raise ValueError(
                f"File has Iteration: 0 — this is a Startup file. "
                f"Use OccamStartup.read() instead: {p}"
            )

        obj = cls(**kwargs)
        _apply_parsed(obj, parsed)

        if obj.verbose:
            obj.logger.info(
                "OccamIter.read: iter %d, misfit=%.4f, %d params from %s",
                obj.iteration,
                obj.misfit_value,
                obj.n_params,
                p,
            )
        return obj

    def write(self, path: PathLike) -> Path:
        """Write this iteration in ``OCCAMITER_FLEX`` format.

        Parameters
        ----------
        path : path-like
            Destination path for the iteration file. Parent
            directories are created when needed.

        Returns
        -------
        pathlib.Path
            Path to the file that was written.
        """
        return _write_iter_flex(self, Path(path))

    # ------------------------------------------------------------------
    # Derived
    # ------------------------------------------------------------------
    def to_resistivity(self) -> np.ndarray:
        r"""Return resistivity values from log10 parameters.

        Returns
        -------
        numpy.ndarray of float
            Resistivity values in ohm metres computed as
            :math:`10^m`, where :math:`m` is each element of
            :attr:`param_values`.
        """
        if not self.param_values.size:
            return np.array([])
        return 10.0**self.param_values

    @property
    def log10_rho_stats(self) -> dict:
        """Return summary statistics for log10 resistivity."""
        if not self.param_values.size:
            return {}
        v = self.param_values
        return {
            "min": float(np.nanmin(v)),
            "max": float(np.nanmax(v)),
            "mean": float(np.nanmean(v)),
            "median": float(np.nanmedian(v)),
        }


# -----------------------------------------------------------------------
# Shared writer (used by both OccamStartup.write and OccamIter.write)
# -----------------------------------------------------------------------

# Ordered write sequence: (file keyword, attribute name, value type)
_WRITE_ORDER = [
    ("Format", "format_str", "str"),
    ("Description", "description", "str"),
    ("Model File", "model_file", "str"),
    ("Data File", "data_file", "str"),
    ("Date/Time", "datetime_str", "str"),
    ("Iterations to run", "max_iterations", "int"),
    ("Target Misfit", "target_misfit", "float"),
    ("Roughness Type", "roughness_type", "int"),
    ("Diagonal Penalties", "diagonal_penalties", "int"),
    ("Stepsize Cut Count", "stepsize_cut_count", "int"),
    ("Debug Level", "debug_level", "int"),
    ("Iteration", "iteration", "int"),
    ("Lagrange Value", "lagrange_value", "float"),
    ("Roughness Value", "roughness_value", "float"),
    ("Misfit Value", "misfit_value", "float"),
    ("Misfit Reached", "misfit_reached", "int"),
    ("Param Count", "n_params", "int"),
]

_KEY_WIDTH = 20  # columns reserved for "Key:" + padding


def _write_iter_flex(obj: OccamBase, path: Path) -> Path:
    """Serialise an OCCAMITER_FLEX object to *path*."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines: list[str] = []

    for kw, attr, vtype in _WRITE_ORDER:
        val = getattr(obj, attr, None)
        if val is None:
            val = 0
        key_str = f"{kw}:"
        pad = max(1, _KEY_WIDTH - len(key_str))
        if vtype == "int":
            val_str = str(int(val))
        elif vtype == "float":
            val_str = str(float(val))
        else:
            val_str = str(val)
        lines.append(f"{key_str}{' ' * pad}{val_str}\n")

        if attr == "stepsize_cut_count":
            # Insert comment lines between Stepsize Cut Count and Debug Level
            lines.append("!Model Limits:      none\n")
            lines.append("!Model Value Steps: none\n")

    # Param block: 4 values per line
    pv = getattr(obj, "param_values", np.array([]))
    for i in range(0, len(pv), 4):
        chunk = pv[i : i + 4]
        lines.append("  " + "    ".join(f"{v:.4f}" for v in chunk) + "\n")

    with path.open("w") as fh:
        fh.writelines(lines)

    obj.path = path
    return path


# -----------------------------------------------------------------------
# Helper: apply parsed dict to any object that has the right attributes
# -----------------------------------------------------------------------


def _apply_parsed(obj: OccamBase, parsed: dict) -> None:
    """Copy all non-None values from *parsed* onto *obj*."""
    for attr, val in parsed.items():
        if val is not None and hasattr(obj, attr):
            setattr(obj, attr, val)
    # Ensure n_params is consistent with param_values length
    pv = parsed.get("param_values")
    if pv is not None and isinstance(pv, np.ndarray):
        obj.param_values = pv
        if not parsed.get("n_params"):
            obj.n_params = len(pv)
