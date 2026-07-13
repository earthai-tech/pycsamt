# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""
Text formatting utilities for pycsamt.
"""

import re
from collections.abc import Iterable, Sequence
from typing import Any, Optional, Union

__all__ = ["smart_format", "fmt_text", "str2columns", "listing_items_format"]


def smart_format(
    items: Union[Iterable[Any], Any],
    choice: str = "and",
    sep: str = ", ",
    oxford_comma: bool = True,
    quote: bool = True,
) -> str:
    """
    Nicely format a sequence into a human-readable list string.

    Parameters
    ----------
    items : iterable or any
        Elements to format. Non-iterables or strings are returned directly.
    choice : {'and', 'or', ...}, default 'and'
        Conjunction to insert before the last element.
    sep : str, default ', '
        Separator between items.
    oxford_comma : bool, default True
        If True, include a comma before the conjunction when >2 items.
    quote : bool, default True
        If True, wrap each item in repr(); else use str().

    Returns
    -------
    str
        A string like "'a', 'b', and 'c'".

    Examples
    --------
    >>> smart_format(['a','b','c'])
    "'a', 'b', and 'c'"
    >>> smart_format(['a','b'], choice='or', oxford_comma=False)
    "'a' or 'b'"
    >>> smart_format('x')
    "'x'"
    >>> smart_format([])
    ""
    """
    # Treat strings as atomic
    if isinstance(items, str) or not isinstance(items, Iterable):
        s = repr(items) if quote else str(items)
        return s
    seq = list(items)
    if not seq:
        return ""

    # prepare representation
    def fmt(x: Any) -> str:
        return repr(x) if quote else str(x)

    formatted = [fmt(x) for x in seq]
    n = len(formatted)
    if n == 1:
        return formatted[0]
    if n == 2:
        return f"{formatted[0]} {choice} {formatted[1]}"
    # n > 2
    if oxford_comma:
        # e.g. a, b, and c
        body = sep.join(formatted[:-1])
        return f"{body}{sep}{choice} {formatted[-1]}"
    else:
        # e.g. a, b and c
        body = sep.join(formatted[:-1])
        return f"{body} {choice} {formatted[-1]}"


def fmt_text(
    features: Union[Sequence[Any], Iterable[Sequence[Any]], dict],
    headers: Optional[Sequence[str]] = None,
    inline: str = "-",
    width: int = 80,
    col_sep: str = " | ",
) -> str:
    # normalize features to list of rows
    rows = []
    if isinstance(features, dict):
        for k, v in features.items():
            if isinstance(v, Iterable) and not isinstance(v, str):
                row = [k] + list(v)
            else:
                row = [k, v]
            rows.append(row)
    elif isinstance(features, Iterable) and not isinstance(
        features, (str, bytes)
    ):
        # list of rows or single row
        first = next(iter(features), None)
        if (
            not first
            or not isinstance(first, Iterable)
            or isinstance(first, str)
        ):
            rows = [list(features)]
        else:
            rows = [list(r) for r in features]
    else:
        rows = [[features]]
    # determine headers
    ncols = max(len(r) for r in rows)
    if headers is None:
        headers = [str(i) for i in range(ncols)]
    if len(headers) < ncols:
        headers = list(headers) + [""] * (ncols - len(headers))
    # compute column widths
    cols = list(zip(*rows))
    widths = []
    for i, col in enumerate(cols):
        max_cell = max(len(str(cell)) for cell in col)
        widths.append(max(len(headers[i]), max_cell))
    # build line
    total_width = sum(widths) + len(col_sep) * (ncols - 1)
    line = inline * min(width, total_width)
    # build header
    header_cells = [f"{headers[i]:^{widths[i]}}" for i in range(ncols)]
    header = col_sep.join(header_cells)
    # build rows
    body = []
    for r in rows:
        cells = []
        for i in range(ncols):
            val = r[i] if i < len(r) else ""
            cells.append(f"{str(val):^{widths[i]}}")
        body.append(col_sep.join(cells))
    # assemble
    out = [line, header, line] + body + [line]
    return "\n".join(out)


def str2columns(
    text: str,
    regex: Optional[Union[str, re.Pattern]] = None,  # r'[#&.*@!_,;\s-]\s*'
    pattern: Optional[str] = None,
    lower: bool = False,
    unique: bool = False,
    strip_chars: Optional[str] = None,
) -> list:
    """
    Split a string into tokens using non-alphanumeric separators.

    Remove all string non-alphanumeric and some operator indicators,  and
    fetch attributes names.

    """
    txt = str(text)
    # determine splitter
    if isinstance(regex, re.Pattern):
        splitter = regex
    else:
        pat = (
            regex if isinstance(regex, str) else (pattern or r"[^0-9A-Za-z]+")
        )
        splitter = re.compile(pat)
    tokens = [tok for tok in splitter.split(txt) if tok]
    if strip_chars:
        tokens = [tok.strip(strip_chars) for tok in tokens]
    if lower:
        tokens = [tok.lower() for tok in tokens]
    if unique:
        seen = set()
        uniq = []
        for tok in tokens:
            if tok not in seen:
                seen.add(tok)
                uniq.append(tok)
        tokens = uniq
    return tokens


def listing_items_format(
    items: Any,
    /,
    begin_text: str = "",
    end_text: str = "",
    bullet: str = "-",
    enumerate_items: bool = True,
    indent: int = 3,
    inline: bool = False,
    verbose: bool = True,
) -> Optional[str]:
    """Format list by enumerate them successively with carriage return"""
    # normalize to list
    from .arrayops import is_iterable

    try:
        seq = list(items) if is_iterable(items) else [items]
    except Exception:
        seq = [items]
    # header
    out_lines = []
    if begin_text:
        out_lines.append(begin_text)
        if verbose:
            print(begin_text)
    # list items
    for idx, val in enumerate(seq, 1):
        prefix = f"{idx}. " if enumerate_items else f"{bullet} "
        line = " " * indent + prefix + str(val)
        out_lines.append(line)
        if verbose:
            print(line)
    # footer
    if end_text:
        out_lines.append(end_text)
        if verbose:
            print(end_text)
    result = "".join(out_lines)
    return None if verbose else result
