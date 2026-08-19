# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
from __future__ import annotations

import re
from pathlib import Path
from typing import Any

from ..exceptions import EdIDataError
from ..log.logger import get_logger
from .base import EDIComponentBase
from .utils import (
    gather_measurement_key_value_with_str_parser,
)
from .validation import (
    IsEdi,
    _split_comment,
    _strip_norm,
    _to_int_or_none,
)

logger = get_logger(__name__)

__all__ = ["MTEMAP", "EMAPComponents", "EMAPMixin"]


def _kv_tokens_from_lines(lines: list[str]) -> list[str]:
    """
    Try project parser first, then robust regex fallback.
    """
    tokens = gather_measurement_key_value_with_str_parser(lines)
    if tokens:
        return tokens

    out: list[str] = []
    pat = re.compile(
        r"([A-Za-z][A-Za-z0-9_]*)\s*=\s*"
        r"("  # group value (quoted or bare)
        r'"[^"]*"|'
        r"'[^']*'|"
        r"[^,\s;]+"  # bare value until comma/space/semicolon
        r")"
    )
    for raw in lines:
        before, _ = _split_comment(raw)
        s = before.strip()
        if not s or s.startswith(">"):
            continue
        for k, v in pat.findall(s):
            vv = _strip_norm(v)
            out.append(f"{k}={vv}")
    return out


class MTEMAP(EDIComponentBase):
    r"""
    Minimal header container for ``>=MTSECT`` and
    ``>=EMAPSECT`` blocks.

    The class stores option key-values from the section
    header. It keeps measurement IDs as strings to avoid
    loss of formatting (e.g. ``"251.025"``). Count fields
    are parsed as integers where possible.

    Notes
    -----
    - If ``SECTID`` is missing or looks numeric, the value
      from ``>HEAD/DATAID`` is used as a fallback.
    - Presence of ``NDIPOLE`` or ``TYPE`` marks the block
      as EMAP when writing.
    - ``NFREQ`` in the header should match the number of
      items in the following ``>FREQ`` block [1]_.
    - The class does *not* validate that the referenced
      measurement IDs exist. That is handled elsewhere.

    Attributes
    ----------
    sectid : str or None
        Section identifier. May be replaced by ``DATAID``
        when numeric or missing.
    nfreq : int or None
        Number of frequencies expected in the section.
    maxblks : int or None
        Optional maximum count of data blocks.
    hx, hy, hz : str or None
        Magnetic measurement IDs. May be absent.
    ex, ey : str or None
        Electric measurement IDs. May be absent.
    rx, ry : str or None
        Reference measurement IDs. May be absent.
    ndipole : int or None
        EMAP dipole count. Presence implies EMAP.
    type : str or None
        Optional EMAP type tag.
    chksum : int or None
        Optional checksum reported by some writers.
    start_data_lines_num : int or None
        Index in the source file where data blocks start.
    temp_sectid : str or None
        Internal cache of ``DATAID`` used for fallback.

    See Also
    --------
    EMAPComponents
        Small container exposing only component IDs.
    EMAPMixin
        Thin facade to build an ``MTEMAP`` from a file.

    References
    ----------
    .. [1] SEG EDI (MT/EMAP) standard. See the spectra
           and MT section rules for frequency lists.
           https://www.mtnet.info/docs/seg_mt_emap_1987.pdf

    Examples
    --------
    Parse a header from an EDI file and write it back::

        from pycsamt.seg.mtemap import MTEMAP

        m = MTEMAP.from_file("site.edi")
        lines = m.write()
        print("".join(lines))
    """

    KEY_ORDER: list[str] = [
        "sectid",
        "nfreq",
        "maxblks",
        "hx",
        "hy",
        "hz",
        "ex",
        "ey",
        "rx",
        "ry",
        "ndipole",
        "type",
        "chksum",
    ]

    sectid: str | None
    nfreq: int | None
    maxblks: int | None
    hx: str | None
    hy: str | None
    hz: str | None
    ex: str | None
    ey: str | None
    rx: str | None
    ry: str | None
    ndipole: int | None
    type: str | None
    chksum: int | None

    start_data_lines_num: int | None = None
    temp_sectid: str | None = None

    def __init__(
        self,
        mt_or_emap_section_list: list[str] | None = None,
        verbose: int | bool = 0,
        logger=None,
        **kwargs: Any,
    ):
        super().__init__(verbose=verbose, logger=logger)
        self.sectid = None
        self.nfreq = None
        self.maxblks = None
        self.hx = None
        self.hy = None
        self.hz = None
        self.ex = None
        self.ey = None
        self.rx = None
        self.ry = None
        self.ndipole = None
        self.type = None
        self.chksum = None

        self._section_lines: list[str] | None = mt_or_emap_section_list

        for k, v in kwargs.items():
            setattr(self, k, v)

        if self._section_lines is not None:
            self.read(self._section_lines)

    @classmethod
    def from_file(cls, edi_path: str) -> MTEMAP:
        r"""
        Build an instance from an EDI path. The file is first
        validated, then the header lines between the section
        tag and the next block are parsed.

        Parameters
        ----------
        edi_path : str
            Path to a readable EDI file.

        Returns
        -------
        MTEMAP
            Parsed header object.

        Raises
        ------
        FileNotFoundError
            If the path does not point to a file.
        EdIDataError
            If no ``>=MTSECT`` or ``>=EMAPSECT`` is found.

        Notes
        -----
        The method also caches ``DATAID`` from ``>HEAD`` and
        uses it to replace a numeric or missing ``SECTID``.
        This helps with writers that output numeric section
        labels.

        Examples
        --------
        >>> m = MTEMAP.from_file("site.edi")
        >>> m.sectid
        'SITE_A'
        """

        p = Path(edi_path)
        IsEdi._assert_edi(p, deep=True)

        # Read once; later slices will reuse these lines.
        lines = p.read_text(
            encoding="utf-8-sig", errors="replace"
        ).splitlines()

        # Capture DATAID early for fallback SECTID.
        dataid = None
        for ln in lines:
            u = ln.upper()
            if re.match(r"^\s*DATAID\s*=", u):
                try:
                    dataid = _strip_norm(ln.split("=", 1)[1])
                except Exception:
                    pass
                break

        # Locate >=MTSECT or >=EMAPSECT
        start_idx = None
        for i, ln in enumerate(lines):
            u = ln.upper().lstrip()
            if u.startswith(">=MTSECT") or u.startswith(">=EMAPSECT"):
                start_idx = i
                break

        if start_idx is None:
            raise EdIDataError("No >=MTSECT or >=EMAPSECT found.")

        # Find end of this header block: next data/meta block.
        stop_idx = len(lines)
        for j in range(start_idx + 1, len(lines)):
            u = lines[j].upper().lstrip()
            if (
                u.startswith(">FREQ")
                or u.startswith(">ZROT")
                or u.startswith(">RHOROT")
                or u.startswith(">=")
            ):
                stop_idx = j
                break

        section_lines = [ln for ln in lines[start_idx + 1 : stop_idx]]

        tokens = _kv_tokens_from_lines(section_lines)

        # inst = cls(mt_or_emap_section_list=tokens)
        # inst.start_data_lines_num = stop_idx
        # inst.temp_sectid = dataid

        # Create instance *without* lines so __init__ doesn't
        # auto-read. Then set fallback and call read() ourselves.
        inst = cls(verbose=0, logger=None)
        inst.start_data_lines_num = stop_idx
        inst.temp_sectid = dataid
        inst.read(tokens)

        return inst

    def read(
        self,
        mt_or_emap_section_list: list[str] | None = None,
    ) -> MTEMAP:
        r"""
        Parse key-value tokens into attributes.

        Parameters
        ----------
        mt_or_emap_section_list : list of str, optional
            Items like ``["NFREQ=60", "HX=251.025"]``. If
            omitted, the list passed at construction is used.

        Returns
        -------
        MTEMAP
            The instance itself, for chaining.

        Notes
        -----
        Unknown keys are ignored. Values are normalized by
        stripping quotes and whitespace. Integer fields are
        parsed with a best-effort cast.
        """

        if mt_or_emap_section_list is not None:
            self._section_lines = mt_or_emap_section_list

        if not self._section_lines:
            raise EdIDataError("No MT/EMAP section lines to read.")

        for raw in self._section_lines:
            if "=" not in raw:
                continue
            k, v = raw.split("=", 1)
            key = _strip_norm(k).lower()
            val = _strip_norm(v)

            if key == "sectid":
                self.sectid = val
            elif key == "nfreq":
                self.nfreq = _to_int_or_none(val)
            elif key == "maxblks":
                self.maxblks = _to_int_or_none(val)
            elif key in (
                "hx",
                "hy",
                "hz",
                "ex",
                "ey",
                "rx",
                "ry",
            ):
                setattr(
                    self,
                    key,
                    val if val not in ("", "NONE") else None,
                )
            elif key == "ndipole":
                self.ndipole = _to_int_or_none(val)
            elif key == "type":
                self.type = val
            elif key == "chksum":
                self.chksum = _to_int_or_none(val)

        if (
            (self.sectid is None) or self._looks_numeric(self.sectid)
        ) and self.temp_sectid:
            self.sectid = self.temp_sectid

        return self

    def write(self, nfreq: int | None = None) -> list[str]:
        r"""
        Serialize the header to text lines.

        Parameters
        ----------
        nfreq : int, optional
            Override for ``NFREQ`` in the output.

        Returns
        -------
        list of str
            Lines including the opening section tag.

        Notes
        -----
        If either ``NDIPOLE`` or ``TYPE`` is present the
        opening tag is ``>=EMAPSECT``. Otherwise ``>=MTSECT``
        is used. For EMAP, the fields ``EX``, ``EY`` and
        ``HZ`` are omitted by convention.
        """

        is_emap = (self.ndipole is not None) or (self.type not in (None, ""))
        header = ">=EMAPSECT\n" if is_emap else ">=MTSECT\n"
        lines: list[str] = [header]

        values: dict[str, Any] = {
            "sectid": self.sectid,
            "nfreq": nfreq if nfreq is not None else self.nfreq,
            "maxblks": self.maxblks,
            "hx": self.hx,
            "hy": self.hy,
            "hz": self.hz,
            "ex": self.ex,
            "ey": self.ey,
            "rx": self.rx,
            "ry": self.ry,
            "ndipole": self.ndipole,
            "type": self.type,
            "chksum": self.chksum,
        }

        omit_for_emap = {"ex", "ey", "hz"} if is_emap else set()

        for key in self.KEY_ORDER:
            if key in omit_for_emap:
                continue
            val = values.get(key, None)
            if val in (None, "", "None"):
                continue

            if key in (
                "hx",
                "hy",
                "hz",
                "ex",
                "ey",
                "rx",
                "ry",
            ):
                fmt = "{:>12}"
            else:
                fmt = "{:>7}"

            lines.append(f"  {key.upper()}={fmt.format(str(val).upper())}\n")

        return lines

    @staticmethod
    def _looks_numeric(s: str) -> bool:
        try:
            float(s)
            return True
        except Exception:
            return False


class EMAPComponents(EDIComponentBase):
    r"""
    Lightweight view of EMAP/MT component IDs.

    The container exposes only component identifiers and
    does not include section counts or metadata. It is
    useful when you need the mapping from channel type to
    measurement ID without other header details.

    Attributes
    ----------
    hx, hy, hz : str or None
        Magnetic component measurement IDs.
    ex, ey : str or None
        Electric component measurement IDs.
    rx, ry : str or None
        Reference component measurement IDs.

    Notes
    -----
    The class does not normalize or validate that the IDs
    exist in the ``>=DEFINEMEAS`` block. It simply carries
    the strings.

    See Also
    --------
    MTEMAP
        Full section header with counts and flags.

    Examples
    --------
    Create from an existing ``MTEMAP``::

        m = MTEMAP.from_file("site.edi")
        comp = EMAPComponents.from_mtemap(m)
        comp.as_dict()
    """

    def __init__(
        self,
        *args,
        verbose: int | bool = 0,
        logger=None,
        **kws: Any,
    ):
        super().__init__(verbose=verbose, logger=logger)
        self.hx: str | None = kws.get("hx")
        self.hy: str | None = kws.get("hy")
        self.hz: str | None = kws.get("hz")
        self.ex: str | None = kws.get("ex")
        self.ey: str | None = kws.get("ey")
        self.rx: str | None = kws.get("rx")
        self.ry: str | None = kws.get("ry")

    @classmethod
    def from_mtemap(cls, m: MTEMAP) -> EMAPComponents:
        r"""
        Construct a component view from an ``MTEMAP``.

        Parameters
        ----------
        m : MTEMAP
            A parsed section header.

        Returns
        -------
        EMAPComponents
            Container populated with IDs from ``m``.
        """
        return cls(
            hx=m.hx,
            hy=m.hy,
            hz=m.hz,
            ex=m.ex,
            ey=m.ey,
            rx=m.rx,
            ry=m.ry,
        )

    def as_dict(self) -> dict[str, str | None]:
        r"""
        Return a plain mapping of component names to IDs.

        Returns
        -------
        dict
            Keys are ``"hx"``, ``"hy"``, ``"hz"``, ``"ex"``,
            ``"ey"``, ``"rx"``, ``"ry"``. Values are strings
            or ``None``.
        """

        return {
            "hx": self.hx,
            "hy": self.hy,
            "hz": self.hz,
            "ex": self.ex,
            "ey": self.ey,
            "rx": self.rx,
            "ry": self.ry,
        }


class EMAPMixin:
    r"""
    Thin facade that lets a host class load an ``MTEMAP``
    with a consistent API.

    Attach this mixin to any class that should expose a
    ``from_file`` constructor returning an ``MTEMAP``. The
    mixin does not keep state; it only forwards the call.

    See Also
    --------
    MTEMAP.from_file
        The underlying constructor that performs parsing.

    Examples
    --------
    Add the mixin and use it on a host type::

        class Host(EMAPMixin):
            pass


        hdr = Host.from_file("site.edi")
        assert hdr.nfreq is not None
    """

    @classmethod
    def from_file(cls, edi_fn: str) -> MTEMAP:
        return MTEMAP.from_file(edi_fn)
