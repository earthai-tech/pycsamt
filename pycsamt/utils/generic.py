# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: GPL-3.0

import math
import subprocess
import sys
import time
import re 
import warnings 

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
