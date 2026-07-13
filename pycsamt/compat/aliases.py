# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

r"""
Lightweight utilities for API deprecation and compatibility.

This module provides helpers to create *runtime* aliases that
preserve v1.x public names while forwarding calls to the new
v2 implementations. Each alias emits a ``FutureWarning`` on
use, guiding users to migrate with clear, uniform messages.

Exports
-------
- :func:`make_compat_alias`
  Build one alias that warns and forwards.
- :func:`install_compat_aliases`
  Register many aliases from a mapping in one step.
- :func:`compat_alias`
  Decorator that injects an alias at definition time.

Design
------
- Aliases are created with ``functools.wraps`` to keep the
  signature and metadata of the new function.
- Aliases set markers for tooling:
  ``__is_compat_alias__ = True`` and ``__compat_target__``.
- Warnings use ``stacklevel=2`` so messages point to user
  call sites rather than inside the library.
- Keep central maps under ``pycsamt.compat``. Re-export at
  the package root if needed, to preserve imports.

Quickstart
----------
Single alias:

.. code-block:: python

   from pycsamt.compat.aliases import make_compat_alias
   from pycsamt.utils.em import plot_lcurve

   plot_l_curve = make_compat_alias(
       plot_lcurve,
       old_name="plot_l_curve",
       new_name="plot_lcurve",
       since="2.0.0",
       remove_in="3.0.0",
   )

Bulk install in a module:

.. code-block:: python

   from pycsamt.compat.aliases import install_compat_aliases

   _MAP = {"plot_l_curve": plot_lcurve}
   install_compat_aliases(_MAP, g=globals())

Decorator at definition:

.. code-block:: python

   from pycsamt.compat.aliases import compat_alias

   @compat_alias("plot_l_curve", remove_in="3.0.0")
   def plot_lcurve(...):
       ...

Notes
-----
Follow PEP 387 for deprecation timelines and user facing
messaging. Keep removal versions realistic and documented.

References
----------
.. [1] PEP 387 — Backwards Compatibility Policy,
       Python Software Foundation.
"""

from __future__ import annotations

import functools
import warnings
from typing import Callable

__all__ = [
    "make_compat_alias",
    "install_compat_aliases",
    "compat_alias",
]


def make_compat_alias(
    new_func: Callable,
    *,
    old_name: str,
    new_name: str | None = None,
    since: str = "2.0.0",
    remove_in: str = "2.7.0",
    extra: str | None = None,
) -> Callable:
    r"""
    Create a callable that acts as a deprecation alias for an
    existing function. The alias exposes the *old* public name,
    forwards to the *new* implementation, and emits a
    ``FutureWarning`` on each call.

    The alias preserves metadata via ``functools.wraps``. It
    also marks itself for tooling:

    * ``__is_compat_alias__ = True``
    * ``__compat_target__`` points to the wrapped function

    Parameters
    ----------
    new_func : Callable
        The canonical replacement function (the v2 target)
        that implements the behavior going forward.
    old_name : str
        The public v1.x name to preserve temporarily. The
        returned alias will be bound to this name and will
        display it in warnings.
    new_name : Optional[str], optional
        The public v2 name to show in the message. If not
        given, defaults to ``new_func.__name__``.
    since : str, default "2.0.0"
        Version in which the old name was first deprecated.
        Displayed to the user in the warning text.
    remove_in : str, default "2.17.0"
        Planned removal version for the old name. Used to
        guide user migration timelines.
    extra : Optional[str], optional
        Extra guidance appended to the message. Useful for
        noting argument changes, added options, or links to
        migration notes.

    Returns
    -------
    Callable
        A thin wrapper that warns (``FutureWarning``) and
        calls ``new_func`` with the provided arguments.

    Notes
    -----
    * The warning uses ``stacklevel=2`` so that it points to
      the user call site rather than inside the library.
    * The alias keeps ``old_name`` as its ``__name__`` for
      introspection, while signature and docs are inherited
      from ``new_func`` via ``functools.wraps``.
    * Bulk installers can detect aliases via the
      ``__is_compat_alias__`` attribute. This avoids double
      installation and duplicate warnings when combining
      local decorators with central maps.
    * Prefer placing aliases under ``pycsamt.compat`` and
      optionally re-exporting them at the package root. See
      PEP 387 for deprecation policy guidance [1]_.

    Examples
    --------
    Create a single alias and use it:

    .. code-block:: python

        def plot_lcurve(...):
            ...

        plot_l_curve = make_compat_alias(
            plot_lcurve,
            old_name="plot_l_curve",
            new_name="plot_lcurve",
            since="2.0.0",
            remove_in="2.1.0",
        )

        # Emits FutureWarning, then forwards:
        plot_l_curve(data)

    Add extra migration tips:

    .. code-block:: python

        plot_l_curve = make_compat_alias(
            plot_lcurve,
            old_name="plot_l_curve",
            extra="API gained 'style' and 'ax'.",
        )

    See Also
    --------
    compat_alias
        Decorator that injects a back-compat alias into the
        defining module without wrapping the new function.
    install_compat_aliases
        Bulk installer that registers many aliases from a
        mapping of ``old_name -> new_func``.

    References
    ----------
    .. [1] PEP 387 — Backwards Compatibility Policy,
           Python Software Foundation.

    """

    if new_name is None:
        new_name = getattr(new_func, "__name__", "<?>")

    base = (
        f"'{old_name}' is deprecated since v{since} and will be "
        f"removed in v{remove_in}. Use '{new_name}' instead."
    )

    msg = f"{base} {extra}" if extra else base

    @functools.wraps(new_func)
    def _wrapper(*args, **kwargs):
        warnings.warn(msg, FutureWarning, stacklevel=2)
        return new_func(*args, **kwargs)

    _wrapper.__name__ = old_name
    _wrapper.__is_compat_alias__ = True
    _wrapper.__compat_target__ = new_func
    return _wrapper


def install_compat_aliases(
    mapping: dict[str, Callable],
    *,
    g: dict[str, object],
    since: str = "2.0.0",
    remove_in: str = "2.1.0",
    extras: dict[str, str] | None = None,
) -> None:
    r"""
    Bulk-install deprecation aliases into a target namespace.

    For each ``old_name -> new_func`` pair in ``mapping``, this
    function creates a compatibility alias using
    :func:`make_compat_alias` and binds it into ``g`` (usually
    the caller's ``globals()``). Each alias warns with a
    ``FutureWarning`` and forwards to ``new_func``.

    Parameters
    ----------
    mapping : dict[str, Callable]
        Map from old public names (v1.x) to the new callables.
        Keys must be valid Python identifiers, as they become
        variable names in ``g``.
    g : dict[str, object]
        Target namespace to mutate. Pass ``globals()`` from the
        module where aliases should be created.
    since : str, default "2.0.0"
        Version where the old names were first deprecated.
        Shown in the warning text.
    remove_in : str, default "2.1.0"
        Planned removal version. Guides user migration plans.
    extras : Optional[dict[str, str]], optional
        Optional per-alias advice. Keys match ``old_name`` and
        values are appended to the warning (e.g., argument
        changes, new options, or a doc link).

    Returns
    -------
    None

    Notes
    -----
    - Existing symbols in ``g`` that already carry the marker
      ``__is_compat_alias__`` are skipped. This prevents
      duplicate installation and duplicate warnings when you
      mix local decorators with central maps.
    - The namespace ``g`` is mutated in place. This function
      does not return the created aliases. Extend ``__all__``
      yourself if you wish to export the old names.
    - Aliases are produced by :func:`make_compat_alias`, which
      preserves metadata via ``functools.wraps`` and warns with
      ``stacklevel=2`` so the message points at user code.
    - Place bulk maps in a dedicated ``compat`` module to keep
      deprecation logic centralized. See PEP 387 for policy
      guidance [1]_.

    Examples
    --------
    Install several aliases inside a module:

    .. code-block:: python

        from pycsamt.compat.aliases import install_compat_aliases
        from .utils.em import full_freq, tensor2d

        _ALIAS_MAP = {
            "get_full_frequency": full_freq,
            "get2dtensor": tensor2d,
        }

        _EXTRAS = {
            "get2dtensor": "Use 'tensor2d'.",
        }

        install_compat_aliases(
            _ALIAS_MAP,
            g=globals(),
            since="2.0.0",
            remove_in="3.0.0",
            extras=_EXTRAS,
        )

        __all__ += [n for n in _ALIAS_MAP if n not in __all__]

    See Also
    --------
    make_compat_alias
        Create a single alias with a uniform warning message.
    compat_alias
        Decorator that injects one alias at definition time.

    References
    ----------
    .. [1] PEP 387 — Backwards Compatibility Policy,
           Python Software Foundation.
    """

    for old_name, new_func in mapping.items():
        existing = g.get(old_name)
        if getattr(existing, "__is_compat_alias__", False):
            continue

        extra = extras.get(old_name) if extras else None
        new_name = getattr(new_func, "__name__", None)

        alias = make_compat_alias(
            new_func,
            old_name=old_name,
            new_name=new_name,
            since=since,
            remove_in=remove_in,
            extra=extra,
        )
        g[old_name] = alias


def compat_alias(
    old_name: str,
    *,
    since: str = "2.0.0",
    remove_in: str = "2.7.0",
    extra: str | None = None,
    export: bool = True,
) -> Callable:
    r"""
    Decorator factory that injects a back-compat alias into the
    defining module without wrapping the new function. The alias
    is created via :func:`make_compat_alias` and bound under
    ``old_name`` in the function's module ``globals()``. The
    decorated function remains unmodified, so direct calls do
    not warn; only calls via the alias emit a ``FutureWarning``.

    When ``export`` is ``True`` and the module defines a list
    ``__all__``, the alias name is appended to that list to
    preserve import semantics for users relying on v1 names.

    Parameters
    ----------
    old_name : str
        The v1.x public name to expose as a temporary alias.
        This identifier will be inserted into the module's
        namespace and will warn when called.
    since : str, default "2.0.0"
        Version where the old name became deprecated. Shown in
        the warning message to inform the user.
    remove_in : str, default "2.7.0"
        Planned removal version for the alias. Helps users
        plan migrations before hard removal.
    extra : Optional[str], optional
        Additional guidance appended to the warning text. Useful
        for noting argument changes or linking migration notes.
    export : bool, default True
        If ``True`` and the module defines ``__all__`` as a
        list, append ``old_name`` to ``__all__`` to re-export
        the alias.

    Returns
    -------
    Callable
        A decorator that, when applied to the new function,
        installs the alias in the defining module and returns
        the function unchanged.

    Notes
    -----
    - The decorator checks if a symbol named ``old_name`` in
      the module already has the marker
      ``__is_compat_alias__``. If so, installation is skipped.
      This prevents duplicate aliases when combined with a
      bulk installer.
    - Only the alias warns and forwards. The decorated
      function is not wrapped, avoiding double warnings when
      both names are reachable.
    - Keep alias declarations close to definitions for local
      clarity, or manage them centrally with a bulk map in a
      ``compat`` module. See PEP 387 for deprecation policy
      guidance [1]_.

    Examples
    --------
    Inject a legacy alias at definition time:

    .. code-block:: python

        from pycsamt.compat.aliases import compat_alias

        __all__ = ["plot_lcurve"]

        @compat_alias(
            "plot_l_curve",
            since="2.0.0",
            remove_in="3.0.0",
            extra="Use 'style' and 'ax' kwargs.",
            export=True,
        )
        def plot_lcurve(...):
            ...

        # Available names in the module:
        # - plot_lcurve  (no warning)
        # - plot_l_curve (warns, then forwards)

    Interaction with bulk maps:

    .. code-block:: python

        # If install_compat_aliases(...) later tries to install
        # 'plot_l_curve', it will detect the existing alias via
        # '__is_compat_alias__' and skip re-installation.

    See Also
    --------
    make_compat_alias
        Build a single alias that warns and forwards to the
        new function.
    install_compat_aliases
        Register many aliases from a mapping of
        ``old_name -> new_func``.

    References
    ----------
    .. [1] PEP 387 — Backwards Compatibility Policy,
           Python Software Foundation.
    """

    def _dec(func: Callable) -> Callable:
        g = func.__globals__
        existing = g.get(old_name)
        if getattr(existing, "__is_compat_alias__", False):
            return func

        alias = make_compat_alias(
            func,
            old_name=old_name,
            new_name=getattr(func, "__name__", None),
            since=since,
            remove_in=remove_in,
            extra=extra,
        )
        g[old_name] = alias

        if export and "__all__" in g:
            al = g["__all__"]
            if isinstance(al, list) and old_name not in al:
                al.append(old_name)

        return func

    return _dec
