# Author: LKouadio <etanoyau@gmail.com>
# License: GPL-3.0

"""
Runtime guard for *optional* dependencies.

The :pyfunc:`ensure_pkg` decorator checks that a given third-party
package is importable **before** the decorated object (function,
method or class) is executed.  It can

* silently *ignore* missing packages;
* *warn* or *raise* when the requirement is unmet;
* *auto-install* (or upgrade) the package through *pip* and retry.

The check occurs **once** at decoration time: no penalty at call
time when everything is already satisfied.
"""

from __future__ import annotations

import functools
import sys
from types import ModuleType
from typing import Any, Callable, TypeVar

from ._dependency import import_optional_dependency
from .generic import ensure_package

__all__ = ["ensure_pkg"]

T = TypeVar("T", bound=Callable[..., Any])


def ensure_pkg(  # noqa: D401
    name: str,
    *,
    extra: str = "",
    dist_name: str | None = None,
    min_version: str | None = None,
    auto_install: bool = False,
    use_conda: bool = False,
    errors: str = "raise",
    verbose: int = 1,
) -> Callable[[T], T]:
    """
    Guarantee *package_name* is importable.

    Parameters
    ----------
    name : str
        Canonical ``import`` target, e.g. ``"scipy"`` or
        ``"matplotlib.pyplot"``.
    extra : str, default ``""``
        Text appended to the *ImportError* message.
    dist_name : str, default *None*
        Distribution name on *PyPI*.  When *None* we assume it is the
        same as *package_name*.
    min_version : str, default *None*
        Enforce a minimum version (PEP 440 compliant string).
    auto_install : bool, default *False*
        Attempt a **one-shot** ``pip install --upgrade`` if the import
        fails or the version is too old.
    use_conda : bool, default *False*
        Placeholder for a future *conda-forge* backend (ignored).
    errors : {"raise", "warn", "ignore"}, default ``"raise"``
        Behaviour when the requirement is unsatisfied.
    verbose : int, default ``0``
        *0* → silent — *1* → short log — *≥ 2* → very chatty.

    Notes
    -----
    The wrapped object is returned *unchanged*; there is no runtime
    overhead once the dependency is resolved.

    Examples
    --------
    >>> from pycsamt.utils.deps import ensure_pkg
    >>>
    >>> @ensure_pkg("tqdm", auto_install=True, verbose=1)
    ... def progress_bar(seq):
    ...     from tqdm import tqdm
    ...     for item in tqdm(seq):
    ...         yield item
    """
    dist_name = dist_name or name

    def _try_import() -> ModuleType | None:
        """Import helper, with optional auto-install."""
        try:
            return import_optional_dependency(
                name,
                extra=extra,
                errors=errors,
                min_version=min_version,
            )
        except ImportError as err:
            if not auto_install:
                raise
            if verbose:
                print(f"[deps] installing '{dist_name}' via pip ...",
                      file=sys.stderr)
            ok = ensure_package(
                dist_name,
                install=True,
                upgrade=True,
                verbose=max(verbose - 1, 0),
            )
            if not ok:
                raise err
            # retry once
            return import_optional_dependency(
                name,
                extra=extra,
                errors=errors,
                min_version=min_version,
            )

    # perform the check now (module import time)
    _try_import()

    def _decorator(obj: T) -> T:
        """Simple passthrough – dependency is already validated."""
        return functools.wraps(obj)(obj)

    return _decorator
