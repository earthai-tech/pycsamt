# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Base class shared by ModEM API objects.

The module defines the small common parent used by ModEM
configuration, file containers, builders, runners, result
loaders, and plotting helpers. It centralizes logger creation,
verbosity handling, and compact display behavior while leaving
all ModEM-specific data structures to subclasses.
"""

from __future__ import annotations

from pathlib import Path
from typing import Union

from ...log.logger import get_logger as _get_logger
from .doc import _modem_param_docs as _params

PathLike = Union[str, Path]

__all__ = ["ModEmBase"]


class ModEmBase:
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


ModEmBase.__doc__ = rf"""
Provide shared behavior for ModEM API objects.

``ModEmBase`` is the lightweight parent used by the ModEM v2
subpackage. It gives file containers, configuration objects,
builders, runners, result loaders, and plotting helpers a
common logger, integer verbosity level, and concise object
representation. The class deliberately avoids any assumption
about ModEM file format or dimensionality; subclasses keep
their own model, data, covariance, control, or plotting state.

The initializer attempts cooperative initialization with
``super``. This lets ``ModEmBase`` participate in simple
mixin hierarchies while remaining tolerant of parent classes
that do not accept the forwarded arguments. After cooperative
initialization, the object stores ``verbose`` as an integer
and creates a logger named from the concrete subclass.

Parameters
----------
*args : tuple
    Positional arguments forwarded to the next class in the
    method-resolution order. They are used only when another
    parent class supports cooperative initialization.
{_params.common.verbose}
{_params.common.logger}
**kwargs : dict
    Extra keyword arguments forwarded to the next class in the
    method-resolution order. If the cooperative call is not
    accepted, the error is ignored so direct subclasses can be
    constructed without boilerplate initialization code.

Attributes
----------
verbose : int
    Integer verbosity level converted from ``verbose``. It is
    available to subclasses for optional progress messages,
    parsing diagnostics, and run-status reporting.
logger : logging.Logger
    Logger associated with the concrete subclass. The default
    logger name includes the subclass module and qualified
    class name, making messages easy to trace in larger ModEM
    workflows.

Notes
-----
The string representation is intentionally compact. It prints
up to four public attributes and omits ``logger`` and
``verbose``. This keeps interactive summaries readable when a
subclass stores large NumPy arrays, response blocks, or model
grids.

ModEM inversions solve regularized electromagnetic inverse
problems in which a model :math:`m` is updated to reduce an
objective of the form

.. math::

   \Phi(m) =
   \| W_d (F(m) - d) \|_2^2 +
   \lambda \| W_m (m - m_0) \|_2^2 ,

where :math:`F(m)` is the forward response, :math:`d` is the
observed data vector, :math:`W_d` and :math:`W_m` are data
and model weighting operators, and :math:`\lambda` controls
the trade-off between fit and regularization. ``ModEmBase``
does not implement this objective directly, but it provides
the common object behavior used by the classes that prepare,
run, and inspect these inversions.

See Also
--------
ModEmConfig
    Configuration object holding data, model, covariance, and
    runner parameters.
InputBuilder
    Builds a complete ModEM input set from EDI-like survey
    data.
ModEmRunner
    Launches the selected ModEM executable from a working
    directory.
InversionResult
    Scans a ModEM run directory and loads logs, models, data,
    and covariance outputs.

Examples
--------
>>> from pycsamt.models.modem.base import ModEmBase
>>> base = ModEmBase(verbose=True)
>>> base.verbose
1
>>> "ModEmBase" in repr(base)
True

Subclasses inherit the same logger and representation
behavior:

>>> class TinyModEmObject(ModEmBase):
...     def __init__(self, value, **kwargs):
...         super().__init__(**kwargs)
...         self.value = value
>>> obj = TinyModEmObject(10, verbose=2)
>>> obj.verbose
2
>>> repr(obj)
'TinyModEmObject(value=10)'

References
----------
.. [1] Egbert, G. D., and Kelbert, A., "Computational
   recipes for electromagnetic inverse problems", Geophysical
   Journal International, 189(1), 251-267, 2012,
   doi:10.1111/j.1365-246X.2011.05347.x.
"""
