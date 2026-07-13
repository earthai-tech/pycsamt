# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""SEG-EDI validators (no base-class deps)."""

from __future__ import annotations

import os
import re
from abc import ABC, abstractmethod
from collections.abc import Iterable, Sequence
from pathlib import Path
from typing import (
    Any,
)

from ..exceptions import EdIDataError, FileHandlingError

_TAG_RE = re.compile(r"^\s*>(=)?\s*([A-Za-z0-9!.]+)")
_FLOAT_RE = re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][-+]?\d+)?")
_NFREQ_RE = re.compile(r"\bNFREQ\s*=\s*(\d+)\b", re.IGNORECASE)


def _strip_norm(s: str) -> str:
    return s.strip().strip('"').strip("'")


def _to_int_or_none(v: Any) -> int | None:
    if v in (None, "", "None"):
        return None
    try:
        return int(float(v))
    except Exception:
        return None


def _to_float_or_none(v: Any) -> float | None:
    if v in (None, "", "None"):
        return None
    try:
        return float(v)
    except Exception:
        return None


def _split_comment(s: str) -> tuple[str, str | None]:
    """Return (before, after) around '//' if present."""
    i = s.find("//")
    if i < 0:
        return (s, None)
    return (s[:i], s[i + 2 :].strip())


def _norm_str(value: Any) -> str | None:
    if value is None:
        return None
    s = str(value).strip().strip('"').strip("'")
    return s if s != "" else None


def _is_tag(line: str, tag: str) -> bool:
    s = line.strip()
    if not s.startswith(">"):
        return False
    return s.upper().startswith(tag.upper())


def _extract_tag(line: str) -> str | None:
    s = line.strip()
    if not s.startswith(">"):
        return None
    # up to first space
    parts = s.split(None, 1)
    return parts[0]


def _extract_tag_in(line: str) -> str | None:
    m = _TAG_RE.match(line)
    if not m:
        return None
    eq, tag = m.groups()
    prefix = ">=" if eq else ">"
    return f"{prefix}{tag.upper()}"


def _iter_blocks(lines: Sequence[str]) -> list[str]:
    tags: list[str] = []
    for ln in lines:
        t = _extract_tag_in(ln)
        if t:
            tags.append(t)
    return tags


def _count_freq_values(
    lines: Sequence[str], start_idx: int
) -> tuple[int, int]:
    n = 0
    i = start_idx + 1
    N = len(lines)
    while i < N:
        if _extract_tag_in(lines[i]) is not None:
            break
        n += len(_FLOAT_RE.findall(lines[i]))
        i += 1
    return n, i


def _expected_nfreq_from_header(line: str) -> int | None:
    m = _NFREQ_RE.search(line)
    if m:
        return int(m.group(1))
    mm = re.search(r"//\s*(\d+)\s*$", line)
    if mm:
        return int(mm.group(1))
    return None


def _has_any_data_block(tags: Iterable[str]) -> bool:
    tagset = set(tags)
    cand = {
        ">ZXXR",
        ">ZXXI",
        ">ZYYR",
        ">ZYYI",
        ">ZXYR",
        ">ZXYI",
        ">ZYXR",
        ">ZYXI",
        ">RHOXX",
        ">RHOXY",
        ">RHOYX",
        ">RHOYY",
        ">PHSXX",
        ">PHSXY",
        ">PHSYX",
        ">PHSYY",
        ">TXR.EXP",
        ">TYR.EXP",
        ">TIPMAG",
        ">TIPPHS",
    }
    return any(t in tagset for t in cand)


class IsEdi(ABC):
    r"""
    Abstract base for SEG-EDI validation helpers.

    Subclasses implement :pyattr:`is_valid`.  The static
    method :py:meth:`_assert_edi` provides a robust, file-level
    validator that accepts either a path or an existing
    :class:`IsEdi` instance.  The check is heuristic, fast, and
    tolerant of minor formatting issues.

    The validator recognizes three EDI *families*:

    1. *Impedance-style* files that contain a ``>FREQ`` block
       and a measurement section such as ``>=MTSECT``,
       ``>=DEFINEMEAS``, ``>=EMAPSECT`` or ``>=OTHERSECT``.
    2. *Spectra-style* files that include ``>=SPECTRASECT`` and
       or at least one ``>SPECTRA`` block.
    3. *Time-series* files that include ``>=TSERIESSECT`` and
       or at least one ``>TSERIES`` block.

    In addition, the first top-level tag must be ``>HEAD`` and
    the last must be ``>END``.  On failure the method raises a
    descriptive :class:`~pycsamt.exceptions.EdIDataError`.

    Notes
    -----
    The check is *structural* rather than semantic.  It does
    not validate numerical values or cross-block consistency.
    Files are read as text using ``utf-8-sig`` with
    ``errors="replace"`` to gracefully handle odd encodings.

    Examples
    --------
    Validate an EDI file path:

    >>> from pycsamt.seg.validation import IsEdi
    >>> IsEdi._assert_edi("path/to/site.edi")
    True

    Use with an object that implements ``is_valid``:

    >>> class MyEdi(IsEdi):
    ...     @property
    ...     def is_valid(self):
    ...         return True
    ...
    >>> IsEdi._assert_edi(MyEdi())
    True

    See Also
    --------
    pycsamt.seg.mtemap.MTEMAP.from_file :
        Parse ``>=MTSECT`` / ``>=EMAPSECT`` headers.
    pycsamt.seg.spectra.SpectraSECT.from_file :
        Parse ``>=SPECTRASECT`` headers.
    pycsamt.seg.time_series.TSect.from_file :
        Parse ``>=TSERIESSECT`` headers.

    References
    ----------
    .. [1] SEG (1987). *MT/EMAP EDI Format Standard*. Society of
       Exploration Geophysicists. Available online:
       https://www.mtnet.info/docs/seg_mt_emap_1987.pdf
    """

    @property
    @abstractmethod
    def is_valid(self) -> bool:
        r"""
        Indicate whether the current instance represents a
        structurally valid EDI object.

        This property is used by
        :py:meth:`IsEdi._assert_edi` when a concrete
        :class:`IsEdi` instance is provided instead of a file
        path.

        Returns
        -------
        bool
            ``True`` when the instance is valid.  Subclasses
            decide the exact criteria.

        Notes
        -----
        Implementations should be lightweight and side-effect
        free.  Heavy validation belongs to dedicated readers
        (e.g., MTEMAP, Spectra, Time-Series parsers).
        """
        raise NotImplementedError

    @staticmethod
    def _assert_edi(
        file: str | os.PathLike | IsEdi, *, deep: bool = True
    ) -> bool:
        r"""
        Validate that ``file`` is a SEG-EDI file (or object).

        The function accepts either a path-like object or an
        :class:`IsEdi` instance.  When ``deep`` is ``False``,
        only the extension check ``.edi`` is performed.  When
        ``deep`` is ``True`` (default), the file is scanned for
        required structural tags:

        - First tag must be ``>HEAD`` and last tag ``>END``.
        - At least one of the following families must be seen:

          * *Impedance-style*: ``>FREQ`` plus a measurement
            section (``>=MTSECT``, ``>=DEFINEMEAS``,
            ``>=EMAPSECT`` or ``>=OTHERSECT``).
          * *Spectra-style*: ``>=SPECTRASECT`` and/or
            ``>SPECTRA``.
          * *Time-series*: ``>=TSERIESSECT`` and/or
            ``>TSERIES``.

        Parameters
        ----------
        file : str or os.PathLike or IsEdi
            Path to an EDI file, or an :class:`IsEdi` instance.
        deep : bool, default: True
            If ``False``, perform a light extension check.
            If ``True``, scan the file contents and tags.

        Returns
        -------
        bool
            ``True`` if the input satisfies the structural
            checks.  Otherwise an exception is raised.

        Raises
        ------
        FileHandlingError
            If ``file`` is ``None``.
        FileNotFoundError
            If the path does not exist or is not a file.
        PermissionError
            If the file cannot be read due to permissions.
        EdIDataError
            If the structure does not match any supported EDI
            family, or if ``>HEAD``/``>END`` ordering fails.

        Notes
        -----
        The validator is intentionally permissive and aims to
        support EDI variants found in the field.  It does not
        guarantee that downstream readers will succeed.

        Examples
        --------
        >>> from pycsamt.seg.validation import IsEdi
        >>> IsEdi._assert_edi("survey01.edi", deep=False)
        True
        >>> IsEdi._assert_edi("survey01.edi", deep=True)
        True
        """
        seg_url = "https://www.mtnet.info/docs/seg_mt_emap_1987.pdf"
        base_msg = f"Unrecognized SEG EDI file. See standard: {seg_url}"

        if file is None:
            raise FileHandlingError("NoneType not allowed")

        if isinstance(file, IsEdi):
            if not file.is_valid:
                raise EdIDataError("EDI object reports invalid state")
            return True

        path = Path(file)  # type: ignore[arg-type]
        if not path.exists() or not path.is_file():
            raise FileNotFoundError(f"{str(path)!r} not a file")

        if not deep:
            if path.suffix.lower().lstrip(".") == "edi":
                return True
            raise EdIDataError("SEG-EDI files usually use '.edi' extension")

        try:
            text = path.read_text(encoding="utf-8-sig", errors="replace")
        except PermissionError:
            raise PermissionError("Permission denied while reading file")

        lines = text.splitlines()

        first_tag = None
        for ln in lines:
            t = _extract_tag_in(ln)
            if t:
                first_tag = t
                break

        last_tag = None
        for ln in reversed(lines):
            t = _extract_tag_in(ln)
            if t:
                last_tag = t
                break

        if first_tag is None or last_tag is None:
            raise EdIDataError(base_msg)

        if not first_tag.startswith(">HEAD"):
            raise EdIDataError(f"{base_msg} (first block should be >HEAD)")

        if not last_tag.startswith(">END"):
            raise EdIDataError(f"{base_msg} (last block should be >END)")

        # All top-level tags (upper-cased: '>HEAD', '>=MTSECT', ...)
        tags = list(_iter_blocks(lines))

        # Impedance-style markers
        have_freq = any(t.startswith(">FREQ") for t in tags)
        have_meas_section = any(
            t
            in {
                ">=MTSECT",
                ">=DEFINEMEAS",
                ">=EMAPSECT",
                ">=OTHERSECT",
            }
            for t in tags
        ) or any(t.startswith(">=DEFINEMEAS") for t in tags)

        # Spectra-style markers
        have_spectra_sect = any(t.startswith(">=SPECTRASECT") for t in tags)
        have_spectra_blocks = any(t.startswith(">SPECTRA") for t in tags)

        # Time-series markers
        have_tseries_sect = any(t.startswith(">=TSERIESSECT") for t in tags)
        have_tseries_blocks = any(t.startswith(">TSERIES") for t in tags)

        # NEW: Other-style markers (may exist without >FREQ)
        have_other_sect = any(t.startswith(">=OTHERSECT") for t in tags)
        # You can be permissive and accept header-only OTHER
        valid_other = have_other_sect

        valid_impedance = have_freq and (
            have_meas_section or have_spectra_sect
        )
        valid_spectra = have_spectra_sect or have_spectra_blocks
        valid_tseries = have_tseries_sect or have_tseries_blocks

        if not (
            valid_impedance or valid_spectra or valid_tseries or valid_other
        ):
            raise EdIDataError(
                f"{base_msg} (need >FREQ+section OR "
                f">=SPECTRASECT/>SPECTRA OR "
                f">=TSERIESSECT/>TSERIES OR "
                f">=OTHERSECT; found: "
                f"have_freq={have_freq}, "
                f"have_meas_section={have_meas_section}, "
                f"have_spectra_sect={have_spectra_sect}, "
                f"have_spectra_blocks={have_spectra_blocks}, "
                f"have_tseries_sect={have_tseries_sect}, "
                f"have_tseries_blocks={have_tseries_blocks}, "
                f"have_other_sect={have_other_sect})"
            )

        return True
