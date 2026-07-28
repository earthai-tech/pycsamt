# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""
pycsamt.decorators
"""

from __future__ import annotations

import functools
import inspect
import os
import subprocess
import warnings
from typing import (
    Any,
    Callable,
    Literal,
    TypeVar,
)

import pandas as pd

from .log.logger import get_logger

logger = get_logger(__name__)


__all__ = [
    "noop",
    "Deprecated",
    "GdalDataCheck",
    "ReplaceBy",
    "isdf",
    "check_empty",
    "has_fit",
]

_T = TypeVar("_T", bound=type)
F = TypeVar("F", bound=Callable[..., Any])


def has_fit(
    error: Literal["raise", "warn", "ignore"] = "raise",
) -> Callable[[_T], _T]:
    """
    Ensure a class exposes a ``fit`` method.

    The decorator inspects the *target* class at import-time and
    applies the following logic:

    1. **If** the class already defines :pymeth:`fit`, nothing is
       changed (the decorator is a no-op).

    2. **Else if** a :pymeth:`read` method exists, a thin shim is
       injected so that::

           obj.fit(*args, **kw)

       transparently forwards the call to
       ``obj.read(*args, **kw)``.

    3. **Else** neither *fit* nor *read* are present.
       The *error* policy determines the outcome:

       ``"raise"``
           Raise :class:`AttributeError`
           (default – safest behaviour).

       ``"warn"``
           Emit :class:`RuntimeWarning` **and** install a silent
           no-op ``fit`` placeholder.

       ``"ignore"``
           Silently add the no-op alias without warnings.

    Notes
    -----
    * The decorator is *idempotent* — applying it several times to
      the same class leaves the public API intact.
    * The original :pymeth:`read` signature & docstring are preserved
      on the auto-generated *fit* alias (via
      :func:`functools.wraps`).

    Examples
    --------
    >>> @has_fit("warn")
    ... class Loader:
    ...     def read(self, src, **kw):
    ...         print(f"reading {src!r}")
    >>> loader = Loader()
    >>> loader.fit("dataset.avg")
    reading 'dataset.avg'

    Attempting to decorate a class that lacks both *read* **and**
    *fit* while ``error="raise"``::

        @has_fit()  # default is "raise"
        class Empty:
            pass

    … raises ``AttributeError: Empty defines neither 'fit' nor 'read'``.
    """

    def _decorator(cls: _T) -> _T:  # type: ignore[valid-type]
        # Prevent double-patching
        if getattr(cls, "_has_fit_alias", False):
            return cls

        # Native .fit already there → nothing to do
        if "fit" in cls.__dict__:
            return cls

        read_fn = getattr(cls, "read", None)
        if inspect.isfunction(read_fn):
            # Build a thin wrapper re-exporting .read as .fit
            @functools.wraps(read_fn)
            def _fit(self, *a: Any, **k: Any):
                return read_fn(self, *a, **k)  # type: ignore[arg-type]

            cls.fit = _fit  # type: ignore[attr-defined]
            cls._has_fit_alias = True
            return cls

        # No suitable candidate found — decide according to *error*
        msg = f"{cls.__name__} defines neither 'fit' nor 'read'"
        if error == "raise":
            raise AttributeError(msg)
        if error == "warn":
            warnings.warn(msg, RuntimeWarning, stacklevel=2)

        # Silent (or warned) no-op fallback
        def _noop(self, *a: Any, **k: Any) -> None:
            return None

        cls.fit = _noop  # type: ignore[attr-defined]
        cls._has_fit_alias = True
        return cls

    return _decorator


def noop(reason: str | None = None) -> Callable[[F], F]:
    """
    No-op decorator. Returns the original function or class unchanged.

    Parameters
    ----------
    reason : str, optional
        Placeholder message for skipping or future implementation.
    """

    def decorator(obj: F) -> F:
        return obj

    return decorator


class Deprecated:
    """
    Decorator to mark functions or classes as deprecated.
    Emits a DeprecationWarning on use and updates the docstring.

    Parameters
    ----------
    reason : str
        Explanation for deprecation and guidance for replacement.
    """

    def __init__(self, reason: str) -> None:
        if not isinstance(reason, str) or not reason.strip():
            raise ValueError("A non-empty deprecation reason must be provided.")
        self.reason = reason

    def __call__(self, obj: F) -> F:
        name = obj.__name__
        original_doc = obj.__doc__ or ""
        obj.__doc__ = f"[DEPRECATED] {self.reason}\n\n{original_doc}"

        @functools.wraps(obj)
        def wrapper(*args: Any, **kwargs: Any) -> Any:
            warnings.warn(
                f"{name} is deprecated: {self.reason}",
                category=DeprecationWarning,
                stacklevel=2,
            )
            return obj(*args, **kwargs)

        return wrapper  # type: ignore


class GdalDataCheck:
    """
    Decorator ensuring GDAL_DATA is set and valid before function execution.

    If GDAL_DATA is missing, attempts `gdal-config --datadir`. Logs status
    and optionally raises ImportError.

    Parameters
    ----------
    raise_on_missing : bool
        If True, raises ImportError when GDAL data cannot be located.
    """

    _checked: bool = False
    _available: bool = False

    def __init__(self, raise_on_missing: bool = False) -> None:
        self.raise_on_missing = raise_on_missing
        if not GdalDataCheck._checked:
            GdalDataCheck._available = self._locate_gdal_data()
            GdalDataCheck._checked = True

    def __call__(self, func: F) -> F:
        @functools.wraps(func)
        def wrapper(*args: Any, **kwargs: Any) -> Any:
            if not GdalDataCheck._available:
                msg = "GDAL data not found; geospatial operations may fail."
                logger.error(msg)
                if self.raise_on_missing:
                    raise ImportError(msg)
            return func(*args, **kwargs)

        return wrapper  # type: ignore

    @staticmethod
    def _locate_gdal_data() -> bool:
        path = os.environ.get("GDAL_DATA")
        if path and os.path.isdir(path):
            logger.info(f"GDAL_DATA found: {path}")
            return True
        try:
            proc = subprocess.run(
                ["gdal-config", "--datadir"],
                capture_output=True,
                check=True,
                text=True,
            )
            path = proc.stdout.strip()
            if os.path.isdir(path):
                os.environ["GDAL_DATA"] = path
                logger.info(f"Set GDAL_DATA to: {path}")
                return True
        except Exception:
            logger.debug("gdal-config lookup failed.")

        if os.environ.get("PYCSAMT_DOCS_BUILD") != "1":
            logger.warning(
                "GDAL_DATA is unavailable. Geospatial features will degrade."
            )
        return False


class ReplaceBy:
    """
    Decorator to alias a deprecated function or class to a new implementation.
    Emits a DeprecationWarning and redirects calls.

    Usage:
        @ReplaceBy(new_func, reason="Use new_func instead.")
        def old_func(...):
            ...
    """

    def __init__(
        self, new_obj: Callable[..., Any] | type, reason: str | None = None
    ) -> None:
        self.new_obj = new_obj
        self.reason = reason or f"Use {new_obj.__name__} instead."

    def __call__(self, old_obj: F) -> F:
        old_name = old_obj.__name__
        new_name = self.new_obj.__name__

        @functools.wraps(old_obj)
        def wrapper(*args: Any, **kwargs: Any) -> Any:
            warnings.warn(
                f"{old_name} is replaced by {new_name}: {self.reason}",
                category=DeprecationWarning,
                stacklevel=2,
            )
            return self.new_obj(*args, **kwargs)

        return wrapper  # type: ignore


def isdf(func):
    """
    Decorator that ensures the first positional argument passed to the
    decorated callable is a pandas DataFrame. If it's not, attempts to convert
    it to a DataFrame using an optional `columns` keyword argument.

    Function is designed to be flexible and efficient, suitable for
    both functions and methods.
    """

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        # Get the signature of the function
        sig = inspect.signature(func)
        params = sig.parameters
        param_list = list(params.values())

        # Check if the function has any parameters
        if not param_list:
            # No parameters to process
            return func(*args, **kwargs)

        # Determine if we're decorating a method (with 'self' or 'cls')
        is_method = False
        start_idx = 0
        if param_list[0].name in ("self", "cls"):
            is_method = True
            start_idx = 1  # Skip 'self' or 'cls'

        # Map arguments to their names
        bound_args = sig.bind_partial(*args, **kwargs)
        bound_args.apply_defaults()

        # Identify the data parameter name
        # Prefer 'data' if it's among the parameters
        data_param_name = None
        for _idx, param in enumerate(param_list[start_idx:], start=start_idx):
            if param.name == "data":
                data_param_name = "data"
                break
        else:
            # If 'data' is not a parameter, use the first positional
            # parameter after 'self'/'cls'
            if (len(param_list) > start_idx) or is_method:
                data_param_name = param_list[start_idx].name
            else:
                # No parameters left to consider
                return func(*args, **kwargs)

        # Get 'data' argument from bound arguments
        data = bound_args.arguments.get(data_param_name, None)
        columns = bound_args.arguments.get("columns", kwargs.get("columns", None))
        if isinstance(columns, str):
            columns = [columns]

        # Proceed with conversion if necessary
        if data is not None and not isinstance(data, pd.DataFrame):
            try:
                if columns and len(columns) != data.shape[1]:
                    data = pd.DataFrame(data)
                else:
                    data = pd.DataFrame(data, columns=columns)
            except Exception as e:
                raise ValueError(
                    f"Unable to convert {type(data).__name__!r} to DataFrame: {e}"
                )
            data.columns = data.columns.astype(str)
            # Update the bound arguments with the new data
            bound_args.arguments[data_param_name] = data

        # Call the original function with the updated arguments
        return func(*bound_args.args, **bound_args.kwargs)

    return wrapper


def check_empty(
    func=None,
    *,
    params=None,
    allow_none=True,
    none_as_empty=False,
    error="raise",
):
    r"""
    Validate that certain inputs are non-empty or valid
    (not :math:`\varnothing`). Can be used in two modes:

    - **No parentheses** (e.g. ``@check_empty``):
      Only the first argument (positional or keyword) of the
      decorated function is checked. This is useful for simple
      scenarios where only one primary parameter requires
      validation.

    - **With parentheses** (e.g. ``@check_empty(params=['x'])``):
      Only the parameters listed in `params` are checked, ignoring
      ``*args`` or ``**kwargs``. This mode offers precise control
      over which parameters must be verified.

    .. math::
        \text{valid\_data}(x) =
        \begin{cases}
        \emptyset & \text{if } x \text{ is invalid},\\
        x & \text{if } x \text{ is non-empty}.
        \end{cases}

    Parameters
    ----------
    func : callable, optional
        The function being decorated. If this parameter is
        provided without parentheses (e.g. ``@check_empty``), then
        the first argument to `func` is validated. If omitted
        (i.e., ``@check_empty(...)``), the returned decorator
        checks only the parameters listed in `params`.
    params : list of str, optional
        A list of parameter names that should be checked for
        emptiness. Only used when the decorator is invoked with
        parentheses (e.g.
        ``@check_empty(params=['data'], allow_none=True)``). If
        this is ``None``, no parameter is explicitly checked in
        this mode.
    allow_none : bool, default=True
        Indicates whether `None` is acceptable. If set to
        ``False``, any encounter of `None` triggers the defined
        error handling strategy in `<error>`.
    none_as_empty : bool, default=False
        If ``True``, `None` is treated as empty
        (i.e. :math:`\varnothing`). Hence, if any parameter is
        `None`, it is considered a violation of non-emptiness,
        regardless of `<allow_none>`.
    error : {'raise', 'warn', 'ignore'}, default='raise'
        Determines the behavior when an empty parameter is
        detected.
        - ``'raise'`` : Raises :class:`ValueError` (stops
                       execution).
        - ``'warn'``  : Issues a :func:`warnings.warn` (allows
                       execution to continue).
        - ``'ignore'``: Takes no action.

    Returns
    -------
    callable
        A decorator function if used with parentheses (or if
        `func` is not given). Otherwise, a wrapper function that
        checks the first argument of the decorated function.

    Notes
    -----
    - When used without parentheses, only the *first argument* of
      the decorated function is checked—no matter if it is
      positional or a keyword.
    - When used with parentheses (and `params` is specified),
      *only* those parameters named in `<params>` are checked. The
      other arguments (including ``*args`` or ``**kwargs``) are
      ignored.
    - To detect emptiness, this function attempts :func:`len`. If
      length is zero, it is considered empty. If no length
      (e.g. an integer), it is not considered empty unless it is
      `None` and `<none_as_empty>` is ``True`` or
      `<allow_none>` is ``False``.

    Examples
    --------
    Minimal usage with no parentheses:

    >>> from fusionlab.core.checks import check_empty
    >>> @check_empty
    ... def process_data(data, *args, **kwargs):
    ...     print("Data:", data)

    Using parentheses with custom parameters:

    >>> from fusionlab.core.checks import check_empty
    >>> @check_empty(params=['data'], allow_none=False)
    ... def load_data(data, path=None):
    ...     print("Loading data:", data, "from", path)

    See Also
    --------
    some_other_checker : Another hypothetical checker function for
        demonstration purposes.

    References
    ----------
    .. [1] *Doe, J.*, *Smith, A.*, "On Data Quality," Journal of
           Integrity, 2021.

    """
    # This function operates in two distinct modes:
    #
    # 1) No parentheses -> e.g. @check_empty
    #    In this scenario, we have a direct reference to the
    #    decorated function as `func`. The decorator will check the
    #    first argument in the function call.
    #
    # 2) With parentheses -> e.g. @check_empty(params=['x'])
    #    Here, `func` is None or not yet bound, and we return a
    #    decorator function that can handle multiple specified
    #    parameters listed in `params`.

    def _handle_violation(param_name):
        """
        Internal method to manage empty or None parameter
        violations based on <error> strategy.
        """
        msg = (
            f"Parameter '{param_name}' is considered empty. Please "
            "provide a valid non-empty input."
        )
        if error == "raise":
            raise ValueError(msg)
        elif error == "warn":
            warnings.warn(msg, UserWarning, stacklevel=2)
        # 'ignore' => do nothing

    def _check_val(val, param_name):
        """
        Internal checker for a single parameter <param_name>
        against emptiness or None depending on <allow_none>,
        <none_as_empty>, etc.
        """
        # If <val> is None and:
        #   - <none_as_empty> is True, or
        #   - <allow_none> is False,
        # => violation
        if val is None and (none_as_empty or not allow_none):
            _handle_violation(param_name)
        else:
            try:
                # Check length to detect empty
                if len(val) == 0:
                    _handle_violation(param_name)
            except TypeError:
                # If object has no __len__, skip.
                pass

    def decorator_with_args(func):
        """
        Decorator for the scenario: @check_empty(params=['...'])
        Only checks the parameters listed in <params>.
        """

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            # Use signature to map positional and keyword
            # arguments to parameter names. Then check only
            # those listed in <params>.
            sig = inspect.signature(func)
            bound = sig.bind(*args, **kwargs)
            bound.apply_defaults()

            # If <params> is provided, iterate over them
            if params:
                for p in params:
                    if p in bound.arguments:
                        _check_val(bound.arguments[p], p)
            return func(*args, **kwargs)

        return wrapper

    def decorator_no_args(func):
        """
        Decorator for the scenario: @check_empty
        Only checks the first argument (positional or
        first keyword).
        """

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            # If there's at least one positional arg, check it;
            # otherwise check the first keyword arg (if any).
            if args:
                first_param_name = func.__code__.co_varnames[0]
                _check_val(args[0], first_param_name)
            elif kwargs:
                first_kwarg = next(iter(kwargs))
                _check_val(kwargs[first_kwarg], first_kwarg)
            return func(*args, **kwargs)

        return wrapper

    # Distinguish between the two usage modes:
    #  - If <func> is a callable, user did @check_empty without
    #    parentheses.
    #  - Otherwise, user did @check_empty(...) with parentheses.
    if callable(func):
        return decorator_no_args(func)
    else:
        return decorator_with_args


def ensure_fit(
    error: Literal["raise", "warn", "ignore"] = "raise",
) -> Callable[[_T], _T]:
    """
    Class decorator that guarantees an alias **fit → read** when the
    target class does **not** already define *fit*.

    Behaviour
    ---------
    * If the class **already** implements ``fit()``, the decorator is a
      no-op.
    * If the class lacks ``fit()`` **but** implements ``read()``, an
      alias is injected so that calls to ``obj.fit(…)`` are transparently
      forwarded to ``obj.read(…)``.
    * If the class implements **neither** method the response depends on
      *error*:

      ``"raise"``   → :class:`AttributeError`
      ``"warn"``    → a *RuntimeWarning* is emitted and ``fit()`` becomes
                      a silent no-op.
      ``"ignore"``  → silent no-op.

    Examples
    --------
    ```python
    @ensure_fit("warn")
    class MyLoader:
        def read(self, src, **kw): ...
    ```
    """

    def _decorator(cls: _T) -> _T:  # type: ignore[valid-type]
        # nothing to do if the class already exposes `.fit`
        if hasattr(cls, "fit"):
            return cls

        # fallback to `.read`
        if hasattr(cls, "read"):

            @functools.wraps(cls.read)
            def _fit(self, *args, **kwargs):  # type: ignore[no-self-use]
                return self.read(
                    *args, **kwargs
                )  # pyright: ignore[reportGeneralTypeIssues]

            cls.fit = _fit  # type: ignore[attr-defined]
            return cls

        # no suitable method found — decide how to react
        msg = f"{cls.__name__} defines neither 'fit' nor 'read'"
        if error == "raise":
            raise AttributeError(msg)
        if error == "warn":
            warnings.warn(msg, RuntimeWarning, stacklevel=2)

        def _noop(self, *a, **k):  # type: ignore[no-self-use]
            return None

        cls.fit = _noop  # type: ignore[attr-defined]
        return cls

    return _decorator
