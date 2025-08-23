# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: GPL-3.0
from __future__ import annotations

import math
import subprocess
import importlib
import sys
import ast
import time
import re 
import warnings 
import inspect
from typing import ( 
    Any, 
    Sequence, 
    List, 
    Union, 
    Optional 
)
import numpy as np

StrLike = Union[str, np.str_]
Container = Union[StrLike, Sequence[StrLike], np.ndarray]

def make_ids(                    # noqa: D401
    items: Sequence[Any] | int,
    *,
    prefix: str | None = "",
    start: int = 0,
    zfill: int | None = None,
    skip_leading_zero: bool = False,
) -> list[str]:
    r"""
    Generate **station/site identifiers** like ``S00``, ``K1_07`` …

    Parameters
    ----------
    items : sequence | int
        *Either* a sequence whose length drives the enumeration,
        *or* an explicit count.
    prefix : str, default *""*
        String prepended to every id – leave empty to get plain
        numeric labels.
    start : int, default ``0``
        First index (use ``1`` for MATLAB-like counting).
    zfill : int, optional
        Pad the numeric part with leading zeroes so that **at
        least** *zfill* digits are printed.  If *None* the width
        is inferred from the largest index **unless**
        *skip_leading_zero* is *True*.
    skip_leading_zero : bool, default *False*
        Remove the zero padding when the computed width is one
        (gives ``"S1"`` instead of ``"S01"``).

    Returns
    -------
    list[str]
        Generated identifiers.

    Notes
    -----
    * When *items* is a sequence **its content is ignored** – only
      ``len(items)`` is used.
    * The routine honours :pycode{numpy.integer} subtype arguments
      and falls back to ``int(count)`` for exotic objects.

    Examples
    --------
    >>> make_ids(range(3), prefix="ix")
    ['ix0', 'ix1', 'ix2']

    >>> make_ids(8, prefix="L", start=1)
    ['L01', 'L02', 'L03', 'L04', 'L05', 'L06', 'L07', 'L08']

    >>> make_ids(12, prefix="line", zfill=4, start=1)
    ['line0001', 'line0002', 'line0003', ... 'line0012']
    """
    # sanitise count
    if isinstance(items, int):
        count = int(items)
    elif isinstance(items, np.integer):
        count = int(items.item())
    else:
        # any iterable
        count = len(items)
    if count < 0:
        raise ValueError("number of ids must be ≥ 0")

    # derive padding
    if zfill is None:
        width = max(1, math.floor(math.log10(count + start)) + 1)
        if skip_leading_zero and width == 1:
            width = 0
    else:
        width = max(0, int(zfill))

    # build ids
    tmpl = f"{{:0{width}d}}" if width else "{:d}"
    ids = [f"{prefix}{tmpl.format(idx)}" for idx in range(
        start, start + count)]
    return ids

def ensure_package(  # noqa: D401
    package: str,
    *,
    install: bool = True,
    upgrade: bool = True,
    extra_pip_args: Sequence[str] | None = None,
    silent: bool = False,
    capture_output: bool = True,
    timeout: float | None = 600,
    verbose: int = 0,
) -> bool:
    r"""
    *Programmatically* install / uninstall a PyPI package.

    This helper wraps :pymod:`subprocess` around
    ``python -m pip …`` so that a project can lazily pull a
    dependency the first time it is required.

    Parameters
    ----------
    package : str
        Distribution name, e.g. ``"tqdm"`` or ``"scipy>=1.11"``.
    install : bool, default *True*
        *True* → run ``pip install``  
        *False* → run ``pip uninstall -y``
    upgrade : bool, default *True*
        Add ``--upgrade`` when *install=True*.
    extra_pip_args : sequence of str, optional
        Additional CLI flags forwarded *as is* to ``pip``.
    silent : bool, default *False*
        Route both *stdout* and *stderr* to :pydata:`subprocess.DEVNULL`.
    capture_output : bool, default *True*
        Capture the child output and log it (ignored if *silent*).
    timeout : float or None, default *600*
        Bail out if the subprocess hangs longer than *timeout*
        seconds.
    verbose : int, default *0*
        * 0 → quiet  
        * 1 → minimal prints  
        * ≥2 → log complete ``pip`` command and duration.

    Returns
    -------
    bool
        *True* if the subprocess returned ``0`` (success).

    Notes
    -----
    *The function always* returns – failures never crash the parent
    process.  Inspect the returned flag and logs to decide what to do
    next.

    Examples
    --------
    >>> from pycsamt.utils.generic import ensure_package
    >>> ensure_package("tqdm", verbose=1)
    ---> installing 'tqdm' (may take a moment) ...
    True
    >>> ensure_package("tqdm", install=False, verbose=1)
    ---> uninstalling 'tqdm' ...
    True
    """
    action = "install" if install else "uninstall"
    cmd: list[str] = [sys.executable, "-m", "pip", action, package]

    if install and upgrade:
        cmd.append("--upgrade")
    if not install:
        cmd.append("-y")
    if extra_pip_args:
        cmd.extend(extra_pip_args)

    if verbose >= 2:
        print("[pip]", " ".join(cmd), file=sys.stderr)

    stdout = subprocess.DEVNULL if silent else (
        subprocess.PIPE if capture_output else None
    )
    stderr = subprocess.DEVNULL if silent else subprocess.STDOUT

    t0 = time.perf_counter()
    try:
        subprocess.check_call(cmd, stdout=stdout, stderr=stderr, timeout=timeout)
        ok = True
    except subprocess.SubprocessError as exc:
        ok = False
        if verbose:
            print(f"[pip] {action} failed → {exc}", file=sys.stderr)
    else:
        if verbose:
            took = time.perf_counter() - t0
            print(f"[pip] {action} succeeded in {took:.1f}s", file=sys.stderr)
    return ok


def strip_item(
    item_to_clean: Optional[Container],
    item: Optional[str] = None,
    multi_space: int = 12,
) -> Optional[Container]:
    """
    Strip a token repeatedly from both ends of strings.

    This utility cleans leading/trailing *repetitions* of a token
    (default: whitespace) in a flexible way. It accepts a scalar
    string, a list of strings, or a NumPy array of strings and
    returns the same container type.

    If the input (or all items within) becomes empty after
    sanitization, ``None`` is returned and a warning is issued.
    This mirrors the legacy behavior where completely blank values
    are treated as missing.

    Parameters
    ----------
    item_to_clean : {str, list of str, np.ndarray of str}, optional
        The text or collection of texts to sanitize. If ``None``,
        returns ``None``.
    item : str, optional
        The token to strip from both ends. If ``None`` or a blank
        string, standard whitespace stripping (``str.strip()``) is
        used. For multi-character tokens (e.g., ``"//"``), the
        token is removed as repeated *whole* substrings rather
        than character sets.
    multi_space : int, default=12
        Maximum repetition count to consider at each end for a
        multi-character token (upper bound in the regex quantifier).
        Must be a positive integer.

    Returns
    -------
    {str, list of str, np.ndarray of str} or None
        Sanitized output preserving the input container type, or
        ``None`` if content is effectively empty after cleaning.

    Notes
    -----
    * For ``item`` that is ``None``/blank, this function applies
      ``str.strip()`` (whitespace). For a non-blank ``item``, it
      removes repeated occurrences of that *exact* token from both
      ends using a compiled regular expression anchored at the start
      and end of the string.
    * Returning ``None`` for fully empty results keeps backward
      compatibility with legacy usage that treats blank fields as
      missing.

    Examples
    --------
    >>> strip_item("    ss_data   ")
    'ss_data'
    >>> strip_item(["  a  ", "   b"], item=None)
    ['a', 'b']
    >>> arr = np.array(["////name////", "////x"], dtype="<U16")
    >>> strip_item(arr, item="//")
    array(['name', 'x'], dtype='<U16')
    >>> strip_item(["      ", "   "]) is None
    True
    """
    # Fast exits and validations
    if item_to_clean is None:
        return None

    if not isinstance(multi_space, int) or multi_space <= 0:
        raise TypeError(
            f"'multi_space' must be a positive integer, got {multi_space!r}"
        )

    # Normalize to a list of Python strings for processing; track original type
    input_is_array = isinstance(item_to_clean, np.ndarray)
    input_is_scalar = isinstance(item_to_clean, (str, np.str_))

    if input_is_scalar:
        data_list: List[str] = [str(item_to_clean)]
        orig_dtype = None
    elif input_is_array:
        # Ensure 1D iteration over elements (even for object dtype)
        arr = np.asarray(item_to_clean)
        orig_dtype = arr.dtype
        data_list = [str(x) if x is not None else "" for x in arr.tolist()]
    elif isinstance(item_to_clean, Sequence):
        data_list = [str(x) if x is not None else "" for x in item_to_clean]
        orig_dtype = None
    else:
        raise TypeError(
            "item_to_clean must be a string, a list/sequence of strings, "
            "or a NumPy array of strings."
        )

    # Build the stripper
    token = (item or "").strip("\n\r")
    cleaned: List[str] = []

    if token == "":
        # Whitespace strip
        for s in data_list:
            cleaned.append(s.strip())
    else:
        # Strip repeated *token* (not char-set) from both ends
        # ^(?:token){1,N} | (?:token){1,N}$
        pat = re.compile(
            rf"^(?:{re.escape(token)}){{1,{multi_space}}}|"
            rf"(?:{re.escape(token)}){{1,{multi_space}}}$"
        )
        for s in data_list:
            # Remove at both ends (pattern matches both anchors)
            ns = pat.sub("", s)
            # Also trim surrounding whitespace that may remain
            cleaned.append(ns.strip())

    # Decide if the result is effectively empty
    if all(c == "" for c in cleaned):
        warnings.warn(
            "No data found for sanitization; returning None.",
            RuntimeWarning,
            stacklevel=2,
        )
        return None

    # Reconstruct original container type
    if input_is_scalar:
        return cleaned[0]
    if input_is_array:
        # Preserve dtype if possible; fallback to unicode
        try:
            return np.array(cleaned, dtype=orig_dtype)
        except Exception:
            return np.array(cleaned)
    # Sequence -> list
    return cleaned

def count_functions(
    module_name: str,
    include_class: bool = False,
    return_counts: bool = True,
    include_private: bool = False,
    include_local: bool = False,
) -> Union[int, List[str]]:
    r"""
    Count or list functions (and classes) in a module.

    Parameters
    ----------
    module_name : str
        Dotted name, e.g. ``"pkg.mod"``.
    include_class : bool, default False
        Include classes in the result.
    return_counts : bool, default True
        Return a count. If ``False``, return sorted names.
    include_private : bool, default False
        Include names that start with ``"_"``.
    include_local : bool, default False
        Include nested (local) functions.

    Returns
    -------
    int or list of str
        Count if ``return_counts=True`` else sorted names.

    Notes
    -----
    Parses the module's AST. Nested functions are excluded
    unless ``include_local=True``.

    Examples
    --------
    >>> count_functions("collections", include_class=True)
    20  # doctest: +SKIP
    """

    try:
        mod = importlib.import_module(module_name)
    except Exception as exc:
        raise ImportError(
            "Cannot import module: {}".format(module_name)
        ) from exc

    try:
        src = inspect.getsource(mod)
    except Exception as exc:
        raise ValueError(
            "Cannot read source for {}".format(module_name)
        ) from exc

    tree = ast.parse(src)

    # attach parent links (for nested detection)
    for node in ast.walk(tree):
        for child in ast.iter_child_nodes(node):
            setattr(child, "parent", node)

    def _is_nested(node: ast.AST) -> bool:
        p = getattr(node, "parent", None)
        while p is not None:
            if isinstance(p, ast.FunctionDef):
                return True
            p = getattr(p, "parent", None)
        return False

    funcs: List[str] = []
    clss: List[str] = []

    for n in ast.walk(tree):
        if isinstance(n, ast.FunctionDef):
            if not include_private and n.name.startswith("_"):
                continue
            if not include_local and _is_nested(n):
                continue
            funcs.append(n.name)
        elif include_class and isinstance(n, ast.ClassDef):
            if not include_private and n.name.startswith("_"):
                continue
            clss.append(n.name)

    names = sorted(funcs + (clss if include_class else []))
    return len(names) if return_counts else names


def get_valid_kwargs(
    callable_obj: Any,
    kwargs: dict[str, Any],
) -> dict[str, Any]:
    r"""
    Filter keyword args to those accepted by a callable.

    Given ``kwargs`` and a target callable (function, method,
    class, or callable instance), return only the keys that
    can be passed as keyword arguments. Extra keys are
    ignored with a warning. If the callable accepts
    ``**kwargs``, all keys are considered valid.

    Parameters
    ----------
    callable_obj : callable or object
        Target to inspect. If an instance is provided, its
        class or ``__call__`` is inspected.
    kwargs : dict
        Candidate keyword arguments.

    Returns
    -------
    dict
        Subset of ``kwargs`` compatible with the callable.

    Notes
    -----
    - Positional-only parameters (``/``) are not eligible as
      keyword arguments and are excluded.
    - If the signature cannot be resolved (e.g., some C
      builtins), an empty dict is returned and a warning is
      emitted.

    Examples
    --------
    >>> def f(a, b=0, *, c=1): ...
    >>> get_valid_kwargs(f, {"a": 1, "x": 9, "c": 2})
    {'a': 1, 'c': 2}
    """
    # Resolve a suitable inspection target
    target = callable_obj
    if (
        not inspect.isclass(target)
        and not inspect.isfunction(target)
        and not inspect.ismethod(target)
        and not callable(target)
    ):
        target = target.__class__

    # Try direct signature; then fall back to __call__
    sig = None
    try:
        sig = inspect.signature(target)
    except (ValueError, TypeError):
        call = getattr(target, "__call__", None)
        if call is not None:
            try:
                sig = inspect.signature(call)
            except (ValueError, TypeError):
                sig = None

    if sig is None:
        warnings.warn(
            "Unable to retrieve signature; no kwargs will be "
            "passed.",
            stacklevel=2,
        )
        return {}

    # If **kwargs present, all keys are valid
    if any(
        p.kind is inspect.Parameter.VAR_KEYWORD
        for p in sig.parameters.values()
    ):
        return dict(kwargs)

    # Eligible names: POSITIONAL_OR_KEYWORD and KEYWORD_ONLY
    eligible_kinds = {
        inspect.Parameter.POSITIONAL_OR_KEYWORD,
        inspect.Parameter.KEYWORD_ONLY,
    }
    valid_names = {
        name
        for name, p in sig.parameters.items()
        if p.kind in eligible_kinds
    }

    valid = {k: v for k, v in kwargs.items() if k in valid_names}
    invalid = [k for k in kwargs.keys() if k not in valid_names]

    if invalid:
        bad = ", ".join(repr(k) for k in invalid)
        warnings.warn(
            "Ignoring invalid keyword(s): {}".format(bad),
            stacklevel=2,
        )

    return valid


def error_policy(
    error: str | None,
    *,
    policy: str = "auto",
    base: str = "ignore",
    exception: type[Exception] = None,
    msg: str | None = None,
    valid_policies: set = None,
) -> str:
    r"""
    Manage error-handling policies like 'warn', 'raise', or 'ignore'.

    The `error_policy` function determines how to handle potential
    errors by mapping the user-provided ``error`` argument to a valid
    policy. It helps standardize responses such as warnings, raised
    exceptions, or silent ignores. The function can adapt to different
    modes, allowing for strict or flexible behavior depending on the
    ``policy`` and ``base`` settings.

    Parameters
    ----------
    error : str or None
        The user-provided error setting. Can be `'warn'`, `'raise'`,
        `'ignore'`, or `None`. If `None`, the behavior is resolved
        based on ``policy`` and ``base``.

    policy : str, default='auto'
        Determines how to interpret a `None` error setting. Valid
        options:

        - `'auto'`: Resolve `None` to the default `base` policy.
        - `'strict'`: Disallows `None` for `error`; raises an error
          if encountered.
        - `None`: Defers strictly to `base`.

    base : str, default='ignore'
        The fallback error policy when `None` is encountered and
        `policy='auto'` or `policy=None`. Must be one of `'warn'`,
        `'raise'`, or `'ignore'`.

    exception : type of Exception, default=ValueError
        The exception class to be raised if an invalid policy or
        error is encountered.

    msg : str, optional
        A custom message for the raised exception if an invalid
        `error` or `policy` is detected. If omitted, a default is
        used.

    Returns
    -------
    str
        A valid error policy: one of `'warn'`, `'raise'`, or
        `'ignore'`.

    Raises
    ------
    ValueError
        If `policy` is invalid or if `None` is not permitted by
        `policy='strict'` but is used. Also raised if `error` cannot
        be resolved to a valid policy or if `base` is invalid when
        `policy='auto'`.

    Notes
    -----
    - If `error` is already a valid policy (`'warn'`, `'raise'`,
      `'ignore'`), it is returned immediately.
    - When `error=None`, the behavior depends on the `policy` and
      `base` parameters. Setting `policy='strict'` disallows `None`
      for `error`.


    .. math::
       \\text{error\\_policy}:
       \\begin{cases}
         \\text{'warn'}, & \\text{issue a warning} \\\\
         \\text{'raise'}, & \\text{raise an exception} \\\\
         \\text{'ignore'}, & \\text{do nothing}
       \\end{cases}


    Examples
    --------
    >>> from pycsamt.utils.generic_utils import error_policy
    >>> # Basic usage:
    >>> resolved_error = error_policy('warn')
    >>> print(resolved_error)
    'warn'

    >>> # Using 'auto' policy with a default base of 'ignore'
    >>> resolved_error = error_policy(None, policy='auto',
    ...                                base='warn')
    >>> print(resolved_error)
    'warn'

    >>> # Strict policy disallows None
    >>> error_policy(None, policy='strict')
    ValueError: In strict policy, `None` is not acceptable as error.

    See Also
    --------
    gofast.utils.validator.validate_nan_policy : A function that
        validate NaN policies.
    """  # noqa: E501

    # Predefined valid policies.
    valid_policies = valid_policies or {"warn", "raise", "ignore"}

    # Default message if none is provided.
    default_msg = (
        "Invalid error policy: '{error}'. Valid options are "
        f"{valid_policies}."
    )
    if exception is None:
        exception = ValueError

    # Use custom message or default if not supplied.
    msg = msg or default_msg

    # Validate the `policy` argument.
    if policy not in {"auto", "strict", None}:
        raise ValueError(
            f"Invalid policy: '{policy}'. Valid options are "
            "'auto', 'strict', or None."
        )

    # Resolve None values for `error` according to `policy`.
    if error is None:
        if policy == "auto":
            # If policy='auto', fallback to `base` if no override is set.
            error = base or "ignore"
        elif policy == "strict":
            # If policy='strict', disallow None for `error`.
            raise ValueError(
                "In strict policy, `None` is not acceptable as an "
                "error. Please set `error` explicitly or switch "
                "policy to 'auto'."
            )
        else:
            # policy=None means strictly use `base` for resolution.
            if base not in valid_policies:
                raise ValueError(
                    f"Invalid base policy: '{base}'. Must be one of "
                    f"{valid_policies} when `error` is None and "
                    "policy is None."
                )
            error = base

    # Final check to ensure `error` is valid.
    if error not in valid_policies:
        # Raise the specified exception if the policy is invalid.
        raise exception(msg.format(error=error))

    # Return the resolved error policy.
    return error
