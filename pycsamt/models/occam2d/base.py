# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Base class shared by all Occam2D model objects.

Follows the same cooperative-init / logger / verbose pattern as
``pycsamt.seg.base.Base`` so that Occam objects slot naturally into the
v2 mixin hierarchy.
"""

from __future__ import annotations

from pathlib import Path
from typing import Union

from ...log.logger import get_logger as _get_logger

PathLike = Union[str, Path]

__all__ = ["OccamBase"]


class OccamBase:
    """Provide shared behavior for Occam2D API objects.

    ``OccamBase`` is the lightweight parent used by the v2
    Occam2D file containers, builders, runners, loaders, and
    plotting helpers. It gives every subclass a consistent
    ``verbose`` flag, class-specific logger, and compact
    representation.

    The class is intentionally small. It does not prescribe
    Occam file format, numerical model, or plotting behavior.
    Subclasses keep their own domain state while inheriting
    same logging and display conventions. The initializer also
    attempts cooperative initialization through ``super``.
    This lets it participate in mixed inheritance trees.

    Parameters
    ----------
    *args : tuple
        Positional arguments forwarded to the next class in
        method-resolution order. They are used only when a
        parent class accepts cooperative initialization.
    verbose : int or bool, default 0
        Verbosity level for progress reporting. ``0`` or
        ``False`` keeps the object quiet. Positive values emit
        messages. Larger values may request diagnostics.
    logger : logging.Logger, optional
        Logger used by the object. If omitted, a PyCSAMT
        logger is created with class name. Pass an explicit
        logger to integrate Occam2D objects into
        an application-level logging configuration.
    **kwargs : dict
        Extra keyword arguments forwarded to the next class in
        the method-resolution order. Unsupported cooperative
        calls are ignored so simple subclasses remain easy to
        construct.

    Attributes
    ----------
    verbose : int
        Integer verbosity level converted from the input.
    logger : logging.Logger
        Logger associated with the concrete subclass.

    Notes
    -----
    The string representation shows up to four public values
    and omits ``logger`` and ``verbose``. This keeps file
    containers readable in notebooks and debugging sessions
    without dumping large arrays or internal state.

    See Also
    --------
    OccamData
        Data-file container inheriting this logging behavior.
    OccamRunner
        Uses the same logger interface to launch the binary.
    InversionResult
        Loads outputs and exposes concise object summaries.

    Examples
    --------
    >>> from pycsamt.models.occam2d.base import OccamBase
    >>> base = OccamBase(verbose=True)
    >>> base.verbose
    1
    >>> "OccamBase" in repr(base)
    True

    References
    ----------
    .. [1] deGroot-Hedlin, C., and Constable, S.,
       "Occam's inversion to generate smooth, two-dimensional
       models from magnetotelluric data", Geophysics, 55(12),
       1613-1624, 1990.
    """

    def __init__(self, *args, verbose: int | bool = 0, logger=None, **kwargs):
        """Initialize verbosity and logger state."""
        try:
            super().__init__(*args, **kwargs)
        except Exception:
            pass

        self.verbose: int = int(verbose)
        self.path: Path | None = None
        name = f"{self.__class__.__module__}.{self.__class__.__qualname__}"
        self.logger = logger if logger is not None else _get_logger(name)

    # ------------------------------------------------------------------
    # Repr helpers
    # ------------------------------------------------------------------
    def __repr__(self) -> str:
        """Return a compact representation of public state."""
        cls = self.__class__.__name__
        attrs = {
            k: v
            for k, v in self.__dict__.items()
            if not k.startswith("_") and k not in ("logger", "verbose")
        }
        pairs = ", ".join(f"{k}={v!r}" for k, v in list(attrs.items())[:4])
        return f"{cls}({pairs})"

    def __str__(self) -> str:
        """Return :meth:`__repr__` for user-facing display."""
        return self.__repr__()
