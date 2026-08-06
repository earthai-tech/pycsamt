# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Base class shared by MARE2DEM API objects."""

from __future__ import annotations

from pathlib import Path
from typing import Union

from ...log.logger import get_logger as _get_logger
from .doc import _mare2dem_param_docs as _params

PathLike = Union[str, Path]

__all__ = ["Mare2DEMBase"]


class Mare2DEMBase:
    def __init__(
        self,
        *args,
        verbose: int | bool = 0,
        logger=None,
        **kwargs,
    ):
        """Initialize verbosity and logger state."""
        try:
            super().__init__(*args, **kwargs)
        except Exception:
            pass

        self.verbose: int = int(verbose)
        module = self.__class__.__module__
        qualname = self.__class__.__qualname__
        name = f"{module}.{qualname}"
        self.logger = logger if logger is not None else _get_logger(name)

    def __repr__(self) -> str:
        cls = self.__class__.__name__
        attrs = {
            k: v
            for k, v in self.__dict__.items()
            if not k.startswith("_") and k not in ("logger", "verbose")
        }
        pairs = ", ".join(f"{k}={v!r}" for k, v in list(attrs.items())[:4])
        return f"{cls}({pairs})"

    def __str__(self) -> str:
        return self.__repr__()


Mare2DEMBase.__doc__ = rf"""
Provide shared behavior for MARE2DEM API objects.

``Mare2DEMBase`` is the lightweight parent used by the MARE2DEM
subpackage. It gives source managers, file containers,
configuration objects, builders, runners, result loaders, and
plotting helpers a common logger, integer verbosity level, and
concise object representation.

Parameters
----------
*args : tuple
    Positional arguments forwarded cooperatively.
{_params.common.verbose}
{_params.common.logger}
**kwargs : dict
    Extra keyword arguments forwarded cooperatively.

Attributes
----------
verbose : int
    Integer verbosity level converted from the ``verbose``
    argument.
logger : logging.Logger
    Logger associated with the concrete subclass, named from
    the subclass module and qualified class name.

Notes
-----
MARE2DEM implements 2.5-D finite-element electromagnetic
inversion for marine and land MT and CSEM data [Mare2DEMBase-1]_. The base
class does not implement any physics; it provides the
infrastructure used by all MARE2DEM wrapper classes.

See Also
--------
Mare2DEMConfig
    Configuration object for MARE2DEM run parameters.
SourceManager
    Download and compile the MARE2DEM Fortran source.
Mare2DEMRunner
    Launch the MARE2DEM binary subprocess.

Examples
--------
>>> from pycsamt.models.mare2dem.base import Mare2DEMBase
>>> base = Mare2DEMBase(verbose=True)
>>> base.verbose
1

References
----------
.. [Mare2DEMBase-1] Key, K. (2016). MARE2DEM: A 2-D inversion code for
   controlled-source electromagnetic and magnetotelluric data.
   *Geophysical Journal International*, 207(1), 571-588.
   doi:10.1093/gji/ggw290.
"""
