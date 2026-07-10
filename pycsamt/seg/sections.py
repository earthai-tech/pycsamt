# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
from __future__ import annotations

import os
from collections.abc import Iterator

from ..exceptions import EdIDataError
from ..log.logger import get_logger
from .base import EDIComponentBase
from .mtemap import MTEMAP
from .other import OtherSECT
from .spectra import SpectraSECT
from .time_series import TSect
from .validation import _extract_tag

logger = get_logger(__name__)

__all__ = ["SECT_REGISTRY", "iter_sections"]


# Tag -> header loader
SECT_REGISTRY: dict[str, type[EDIComponentBase]] = {
    ">=MTSECT": MTEMAP,
    ">=EMAPSECT": MTEMAP,
    ">=SPECTRASECT": SpectraSECT,
    ">=TSERIESSECT": TSect,
    ">=OTHERSECT": OtherSECT,
}


def _find_header_end(lines: list[str], start: int) -> int:
    """
    Return index of first line after header. Stop at next line
    that starts with '>' (either '>BLOCK' or '>=' new section).
    """
    j = start + 1
    n = len(lines)
    while j < n:
        s = lines[j].lstrip()
        if s.startswith(">"):
            break
        j += 1
    return j


def _make_header(
    tag: str,
    header_lines: list[str],
    body_start: int,
) -> EDIComponentBase:
    """
    Create a header instance for tag and feed it the raw lines.
    If the class exposes a .read, call it; then inject
    .start_data_lines_num for downstream IO.
    """
    cls = SECT_REGISTRY[tag]
    inst = cls()  # type: ignore[call-arg]

    # Be liberal: if header exposes .read(lines), use it.
    read = getattr(inst, "read", None)
    if callable(read):
        try:
            read(header_lines)
        except TypeError:
            # Some headers may not accept args; skip.
            read()

    # Attach where body starts for IO loaders.
    try:
        inst.start_data_lines_num = body_start
    except Exception:
        pass

    return inst


def iter_sections(
    edi_path: str,
    *,
    include: list[str] | None = None,
) -> Iterator[tuple[str, EDIComponentBase, int]]:
    """
    Yield (tag, header_obj, body_start_line) for each '>=...' block.

    Parameters
    ----------
    edi_path : str
        Path to EDI file.
    include : list of str, optional
        If given, only yield tags in this list (case-insens.).
        Tags must match keys in SECT_REGISTRY.

    Yields
    ------
    (tag, header, body_start) : tuple
        tag : str, e.g. ">=MTSECT"
        header : parsed header instance
        body_start : int, index of first line after header,
            i.e., where the first '>BLOCK' (or next section)
            begins.
    """
    if not os.path.isfile(edi_path):
        raise FileNotFoundError(f"{edi_path!r} is not a file.")

    with open(edi_path, encoding="utf-8") as f:
        lines = f.readlines()

    want = (
        {t.upper(): True for t in include} if include else None
    )

    i = 0
    n = len(lines)
    while i < n:
        ln = lines[i]
        tag = _extract_tag(ln)
        if not tag or not tag.startswith(">="):
            i += 1
            continue

        tag_u = tag.upper()
        if tag_u not in SECT_REGISTRY:
            # Unknown top-level block: skip to next tag.
            i += 1
            continue

        if want is not None and tag_u not in want:
            i += 1
            continue

        j = _find_header_end(lines, i)
        header_lines = [
            s.rstrip("\n") for s in lines[i + 1 : j]
        ]

        try:
            header = _make_header(tag_u, header_lines, j)
        except Exception as exc:
            raise EdIDataError(
                f"Failed to parse header {tag_u} at line "
                f"{i}: {exc}"
            ) from exc

        yield (tag_u, header, j)
        i = j
