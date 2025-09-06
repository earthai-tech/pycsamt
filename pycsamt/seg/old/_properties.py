# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""SEG-EDI metadata properties.

This module defines light-weight component classes that capture
bibliographic, copyright, provenance, and contact information
used across SEG-EDI headers and related structures.

Classes inherit from :class:`pycsamt.seg.base.EDIComponentBase`
to get consistent ``__repr__``/``__str__`` and rendering helpers.
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Optional, List


from abc import ABC, abstractmethod
from pathlib import Path
from typing import Iterable, Sequence, Tuple, Union

from ..exceptions import EdIDataError, FileHandlingError
from .base import EDIComponentBase# no need to inherit since base.py called 
# already IsEdi. This will avoid circular import 

# revise this class to best follow our new structure ... 



# -----------------------
# Internal parsing helpers
# -----------------------

_TAG_RE = re.compile(r"^\s*>(=)?\s*([A-Za-z0-9!.]+)")
_FLOAT_RE = re.compile(
    r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][-+]?\d+)?"
)  # robust float, E or D exponent
_NFREQ_RE = re.compile(r"\bNFREQ\s*=\s*(\d+)\b", re.IGNORECASE)



_EDI =[
    #Head=Infos-Frequency-Rhorot
    # Zrot and end blocks
    '>HEAD','>INFO',
    
    #Definitions Measurments Blocks
    '>=DEFINEMEAS',
    '>=MTSECT',
    '>FREQ ORDER=INC', 
    '>ZROT',
    '>RHOROT',
    '>!Comment',
    '>FREQ',
    '>HMEAS',
    '>EMEAS',
    
    #Apparents Resistivities
    #and Phase Blocks
    '>RHOXY','>PHSXY',
    '>RHOXY.VAR','>PHSXY.VAR',
    '>RHOXY.ERR','>PHSXY.ERR',
    '>RHOXY.FIT','>PHSXY.FIT', 
    '>RHOYX','>PHSYX',
    '>RHOYX.VAR','>PHSYX.VAR',
    '>RHOYX.ERR','>PHSYX.ERR',
    '>RHOYX.FIT','>PHSYX.FIT',
    '>RHOXX','>PHSXX',
    '>RHOXX.VAR','>PHSXX.VAR',
    '>RHOXX.ERR','>PHSXX.ERR',
    '>RHOXX.FIT','>PHSXX.FIT',
    '>RHOYY','>PHSYY',
    '>RHOYY.VAR','>PHSYY.VAR',
    '>RHOYY.ERR','>PHSYY.ERR',
    '>RHOYY.FIT','>PHSYY.FIT',
    '>FRHOXY','>FPHSXY',
    '>FRHOXY.VAR','>FPHSXY.VAR',
    '>FRHOXY.ERR','>FPHSXY.ERR',
    '>FRHOXY.FIT','>FPHSXY.FIT', 
    '>FRHOXX','>FPHSXX',
    '>FRHOXX.VAR','>FPHSXX.VAR',
    '>FRHOXX.ERR','>FPHSXX.ERR',
    '>FRHOXX.FIT','>FPHSXX.FIT',
    
     #Time series-Sepctra 
     #and EM/OTHERSECT
    '>TSERIESSECT', '>TSERIES', 
    '>=SPECTRASECT', '>=EMAPSECT',
    '>=OTHERSECT',
    
    #Impedance Data Blocks

    '>ZXYR','>ZXYI',
    '>ZXY.VAR', '>ZXYR.VAR',
    '>ZXYI.VAR', '>ZXY.COV', 
    '>ZYXR','>ZYXI',
    '>ZYX.VAR', '>ZYXR.VAR',
    '>ZYXI.VAR', '>ZYX.COV',
    '>ZYYR','>ZYYI',
    '>ZYY.VAR', '>ZYYR.VAR',
    '>ZYYI.VAR', '>ZYY.COV',
    '>ZXXR', '>ZXXI',
    '>ZXXR.VAR','>ZXXI.VAR',
    '>ZXX.VAR','>ZXX.COV',
    '>FZXXR','>FZXXI',
    '>FZXYR','>FZXYI', 
    
    #Continuous 1D inversion 
    '>RES1DXX', '>DEP1DXX', 
    '>RES1DXY', '>DEP1DXY',
    '>RES1DYX', '>DEP1DYX', 
    '>RES1DYY', '>DEP1DYY',
    '>FRES1DXX', '>FDEP1DXX',
    '>FRES1DXY', '>FDEP1DXY',
    
    # Coherency and 
    #Signal Data Blocks
    '>COH','>EPREDCOH',
    '>HPREDCOH','>SIGAMP',
    '>SIGNOISE', 
    
    #Tipper Data blocks 
    '>TIPMAG','>TIPPHS',
    'TIPMAG.ERR', '>TIPMAG.FIT',
    '>TIPPHS.FIT', '>TXR.EXP',
    '>TXI.EXP', '>TXVAR.EXP',
    '>TYR.EXP','>TYI.EXP',
    '>TYVAR.EXP',
    
     # Strike, Skew, and
     # Ellipticity Data Blocks 
    '>ZSTRIKE','>ZSKEW',
    '>ZELLIP', '>TSTRIKE',
    '>TSKEW','>TELLIP',
    
    #Spatial filter blocks
    '>FILWIDTH','>FILANGLE',
    '>EQUIVLEN' , '>END'
] 
#
_ZRP_COMPS =[
    #z
    [
         ['zxxr', 'zxxi', 'zxx.var'],
         ['zxyr', 'zxyi', 'zxy.var'],
         ['zyxr', 'zyxi', 'zyx.var'],
         ['zyyr', 'zyyi', 'zyy.var']
    ],
    #Rho
    [
         ['rhoxx', 'rhoxx.var','rhoxx.err', 'rhoxx.fit'],
         ['rhoxy','rhoxy.var','rhoxy.err', 'rhoxy.fit']
    ],
                                    
    [
         ['phxx','phsxx.var', 'phsxx.err', 'phsxx.fit'],
         ['phxy','phsxy.var', 'phsxy.err', 'phsxy.fit']
    ], 
    # filtered 
    [
         ['frhoxx','frhoxx.var','frhoxx.err', 'frhoxx.fit'],
         ['frhoxy','frhoxy.var', 'frhoxy.err', 'frhoxy.fit']
    ],
                                         
    [
         ['fphsxx','fphsxx.var', 'fphsxx.err', 'fphsxx.fit'],
         ['fphsxy','fphsxy.var', 'fphsxy.err', 'fphsxy.fit']
     ],
 ]
                                                    
_TIP_COMPS =[
    ['txr.exp', 'txi.exp', 'txvar.exp'],
    ['tyr.exp', 'tyi.exp', 'tyvar.exp']
    ]

__all__ = ["References", "Copyright", "Person", "Source", 
           "Software", "Processing", "IsEdi"]


@dataclass
class References(EDIComponentBase):
    """Reference (citation) metadata.

    Parameters
    ----------
    author : str, optional
        Author string (can include initials and separators).
    title : str, optional
        Title of article, report, dataset, or publication.
    journal : str, optional
        Journal or outlet name.
    doi : str, optional
        DOI string, e.g. ``"10.1234/abcd.efgh"``.
    year : int, optional
        Year of publication.
    volume : str or int, optional
        Volume number or identifier.
    pages : str, optional
        Page range, e.g. ``"234--241"`` or article id.

    Notes
    -----
    - Additional ad-hoc metadata can be added as attributes on
      the instance if needed. They will appear in the rendered
      output unless set to ``None``.
    - This class is *sectionless* by default (no ``>BLOCK``)
      because it is typically embedded inside other components.
      Its ``__str__`` prints a compact ``KEY=VALUE`` sequence.

    Examples
    --------
    >>> ref = References(
    ...     author="D. MaryE",
    ...     title="pyCSAMT: Python toolbox for CSAMT",
    ...     journal="Computers & Geosciences",
    ...     year=2021,
    ...     volume="18",
    ...     pages="234--241",
    ...     doi="10.1016/j.cageo.2021.123456",
    ... )
    >>> "journal" in str(ref)
    True
    """
    _section: Optional[str] = None  # keep as embedded/value object

    author: Optional[str] = None
    title: Optional[str] = None
    journal: Optional[str] = None
    doi: Optional[str] = None
    year: Optional[int] = None
    volume: Optional[str] = None
    pages: Optional[str] = None

    def validate(self) -> None:
        if self.year is not None and self.year < 0:
            raise ValueError("References.year must be non-negative.")


@dataclass
class Copyright(EDIComponentBase):
    """Copyright and reuse conditions for the dataset.

    Parameters
    ----------
    references : References, optional
        Citation to published work using the data. Defaults to an
        empty :class:`References` object.
    conditions_of_use : str, optional
        Short license/usage notice. If omitted, a permissive,
        non-commercial default text is provided.
    release_status : {"open", "public", "proprietary"}, optional
        High-level availability flag; not enforced but validated
        against these known values (case-insensitive).
    owner : str, optional
        Rights holder (organization or person).
    contact : str, optional
        Contact person or mailbox.
    additional_info : str, optional
        Free-text details (links, clauses, etc.).

    Notes
    -----
    The default `conditions_of_use` text clarifies that bundled
    test data can be used to exercise the software but should not
    be redistributed commercially without proper citation.

    Examples
    --------
    >>> c = Copyright(owner="University of CSAMT",
    ...               contact="Cagniard",
    ...               release_status="public")
    >>> "release_status" in repr(c)
    True
    """
    _section: Optional[str] = None  # embedded/value object

    references: References = field(default_factory=References)
    conditions_of_use: str = field(
        default_factory=lambda: (
            "Data in the bundled 'data/' directory may be used to test "
            "the pyCSAMT software. Redistribution or commercial use is "
            "not permitted unless properly cited and permitted by the "
            "original source. For additional public MT datasets, see the "
            "IRIS DMC. Data may be copied and redistributed provided the "
            "source is cited."
        )
    )
    release_status: Optional[str] = None
    owner: Optional[str] = None
    contact: Optional[str] = None
    additional_info: Optional[str] = None

    def validate(self) -> None:
        if self.release_status is None:
            return
        allowed = {"open", "public", "proprietary"}
        if self.release_status.lower() not in allowed:
            raise ValueError(
                "release_status must be one of "
                f"{sorted(allowed)} (case-insensitive)."
            )


@dataclass
class Person(EDIComponentBase):
    """Person or contact information.

    Parameters
    ----------
    name : str, optional
        Person's full name.
    organization : str, optional
        Affiliation or organization name.
    email : str, optional
        Contact e-mail address.
    role : str, optional
        Role in data lifecycle, e.g. ``"author"``, ``"submitter"``.
    phone : str, optional
        Contact phone number.

    Examples
    --------
    >>> p = Person(name="A. Cagniard", role="author")
    >>> "name=" in repr(p)
    True
    """
    _section: Optional[str] = None  # embedded/value object

    name: Optional[str] = None
    organization: Optional[str] = None
    email: Optional[str] = None
    role: Optional[str] = None
    phone: Optional[str] = None

    def validate(self) -> None:
        if self.email and "@" not in self.email:
            raise ValueError("Person.email must contain '@'.")


@dataclass
class Source(EDIComponentBase):
    """Provenance of the EDI file (how/when/by whom it was made).

    Parameters
    ----------
    project : str, optional
        Project identifier or name.
    survey : str, optional
        Survey identifier or campaign name.
    sitename : str, optional
        Site or station name.
    creationdate : str, optional
        Creation timestamp in ISO-8601 (UTC), e.g.
        ``"2025-07-10T01:37:08Z"``. If omitted, set on instantiation.
    creatingsoftware : str, optional
        Software name and version creating the file. Defaults to
        ``"pyCSAMT"``.
    author : Person, optional
        Person who created the file. Defaults to empty :class:`Person`.
    recipient : Person, optional
        Person to whom the file is submitted/archived. Defaults to
        empty :class:`Person`.
    archive : str, optional
        Archival location or system (e.g., IRIS).
    reprocessed_by : str, optional
        Free-text note for post-processing agent.

    Notes
    -----
    - ``creationdate`` is stored in UTC with ``"Z"`` suffix.
    - This class is *sectionless* to allow embedding under higher
      level headers; its ``__str__`` prints a compact one-liner.

    Examples
    --------
    >>> s = Source(project="DemoProject", sitename="E1_2")
    >>> "creatingsoftware" in str(s)
    True
    """
    _section: Optional[str] = None  # embedded/value object

    project: Optional[str] = None
    survey: Optional[str] = None
    sitename: Optional[str] = None

    creationdate: str = field(
        default_factory=lambda: datetime.now(timezone.utc)
        .replace(microsecond=0)
        .isoformat()
        .replace("+00:00", "Z")
    )
    creatingsoftware: str = "pyCSAMT"

    author: Person = field(default_factory=Person)
    recipient: Person = field(default_factory=Person)

    archive: Optional[str] = None
    reprocessed_by: Optional[str] = None

    def validate(self) -> None:
        # Minimal sanity checks
        if not self.creatingsoftware:
            raise ValueError("creatingsoftware must be a non-empty string.")


@dataclass
class Software(EDIComponentBase):
    """Software metadata.

    Parameters
    ----------
    name : str, optional
        Name of the software, e.g. ``"pyCSAMT"`` or
        ``"MTPROC"``.
    version : str, optional
        Version string. Flexible (semver preferred), e.g.
        ``"1.2.3"`` or ``"V1.0.7"``.
    release : str, optional
        Release tag or date string, e.g. ``"2025-07-10"``,
        ``"0.11.23"``, or ``"May 02 2022"``.
    author : Person, optional
        Author or maintainer contact. Defaults to an empty
        :class:`Person`.
    url : str, optional
        Project homepage or repository URL.
    license : str, optional
        License identifier, e.g. ``"LGPL-3.0"``.
    description : str, optional
        Short free-text description.

    Notes
    -----
    This class is sectionless (no ``>BLOCK``) by default and
    renders a compact one-liner via :class:`EDIComponentBase`.

    Examples
    --------
    >>> sw = Software(name="pyCSAMT", version="0.11.23",
    ...               release="2025-07-10",
    ...               license="LGPL-3.0")
    >>> "version=" in repr(sw)
    True
    """
    _section: Optional[str] = None  # embedded/value object

    name: Optional[str] = None
    version: Optional[str] = None
    release: Optional[str] = None
    author: "Person" = field(default_factory=Person)
    url: Optional[str] = None
    license: Optional[str] = None
    description: Optional[str] = None

    def validate(self) -> None:
        if self.name is not None and not str(self.name).strip():
            raise ValueError("Software.name must be a non-empty string.")
        if self.version is not None and not str(self.version).strip():
            raise ValueError("Software.version must be a non-empty string.")
        if self.url is not None:
            # minimal sanity check
            u = str(self.url).strip().lower()
            if not (u.startswith("http://") or u.startswith("https://")):
                raise ValueError(
                    "Software.url should start with http:// or https://"
                )

@dataclass
class Processing(EDIComponentBase):
    """Processing metadata for an EDI site/section.

    Parameters
    ----------
    ProcessingSoftware : Software, optional
        Software details used for processing (keeps legacy
        capitalization for compatibility). Defaults to an
        empty :class:`Software`.
    processedby : str, optional
        Name or ID of the processing operator.
    processingtag : str, optional
        Specific tag or label for this processing run, e.g.
        ``"robust-v2"`` or ``"remote-ref(25-E01)"``.
    runlist : list of str, optional
        Free-form list of run identifiers or notes.
    remoteref : str, optional
        Reference point/station used for remote reference.
    remotesite : str, optional
        Reference site name (human-readable).
    signconvention : {"exp(+i ω t)", "exp(-i ω t)"}, optional
        Complex time dependence used in impedance definition.
        Defaults to ``"exp(+i ω t)"``. Greek omega or plain
        ``"w"`` are both accepted.

    Notes
    -----
    - Normalizes sign convention variants such as
      ``exp(+i \\omega t)``, ``exp(+i w t)``, and their
      negative counterparts.
    - This class is sectionless; it is typically embedded
      inside higher-level headers (e.g. INFO).

    Examples
    --------
    >>> p = Processing(
    ...     processedby="pyCSAMT",
    ...     processingtag="robust-v2",
    ...     runlist=["run-001", "run-002"],
    ...     remoteref="RX",
    ...     remotesite="E01-360",
    ... )
    >>> "signconvention" in str(p)
    True
    """
    _section: Optional[str] = None  # embedded/value object

    ProcessingSoftware: Software = field(default_factory=Software)
    processedby: Optional[str] = None
    processingtag: Optional[str] = None
    runlist: Optional[List[str]] = None
    remoteref: Optional[str] = None
    remotesite: Optional[str] = None
    signconvention: str = "exp(+i ω t)"

    @staticmethod
    def _normalize_signconv(s: str) -> str:
        """Normalize sign convention strings to canonical forms.

        Accepts variants like:
        - ``exp(+i ω t)`` / ``exp(+i \\omega t)`` / ``exp(+i w t)``
        - ``exp(-i ω t)`` / ``exp(-i \\omega t)`` / ``exp(-i w t)``
        """
        if s is None:
            return "exp(+i ω t)"
        s0 = str(s).strip().lower()
        # normalize omega to unicode ω
        s0 = s0.replace("\\omega", "ω").replace(" w ", " ω ")
        s0 = s0.replace("w)", "ω)").replace("(w", "(ω")
        # strip duplicate spaces
        s0 = " ".join(s0.split())
        # enforce canonical plus/minus
        if "exp(+i" in s0:
            return "exp(+i ω t)"
        if "exp(-i" in s0:
            return "exp(-i ω t)"
        # fallback to default
        return "exp(+i ω t)"

    # --- lifecycle -------------------------------------------------------

    def validate(self) -> None:
        # normalize and check sign convention
        self.signconvention = self._normalize_signconv(self.signconvention)
        if self.runlist is not None and not isinstance(self.runlist, list):
            raise ValueError("Processing.runlist must be a list of strings.")
        # Ensure runlist items are strings
        if self.runlist is not None:
            for i, v in enumerate(self.runlist):
                if not isinstance(v, str):
                    raise ValueError(
                        f"Processing.runlist[{i}] must be a string, got {type(v)}"
                    )



def _extract_tag(line: str) -> Optional[str]:
    """
    Return a normalized tag (uppercased), e.g. '>HEAD', '>=DEFINEMEAS',
    or None if the line is not a tag line.
    """
    m = _TAG_RE.match(line)
    if not m:
        return None
    eq, tag = m.groups()
    prefix = ">=" if eq else ">"
    # normalize e.g. "ZXX.VAR" or "PHSXY.VAR" as-is, uppercase for matching
    return f"{prefix}{tag.upper()}"


def _iter_blocks(lines: Sequence[str]) -> List[str]:
    """Collect tag headers from the raw text (normalized)."""
    tags: List[str] = []
    for ln in lines:
        t = _extract_tag(ln)
        if t:
            tags.append(t)
    return tags


def _count_freq_values(lines: Sequence[str], start_idx: int) -> Tuple[int, int]:
    """
    From a 'FREQ' header line at 'start_idx', count how many float tokens
    follow until the next '>'-tag line. Returns (n_vals, next_tag_index).

    Robust against 'NFREQ=XX', ' // XX' on the header line; those are
    advisory, we compute the real count from the payload.
    """
    n = 0
    i = start_idx + 1
    N = len(lines)
    while i < N:
        if _extract_tag(lines[i]) is not None:
            break
        n += len(_FLOAT_RE.findall(lines[i]))
        i += 1
    return n, i


def _expected_nfreq_from_header(line: str) -> Optional[int]:
    """
    Try to read an expected NFREQ from the FREQ header, e.g.
    'FREQ  NFREQ=53 // 53' or 'FREQ //52'. Returns None if not found.
    """
    m = _NFREQ_RE.search(line)
    if m:
        return int(m.group(1))
    # Try trailing comment form like '// 52'
    mm = re.search(r"//\s*(\d+)\s*$", line)
    if mm:
        return int(mm.group(1))
    return None


def _has_any_data_block(tags: Iterable[str]) -> bool:
    """
    Check presence of at least one recognized data block:
    any of Z*, RHO*, PHS*, TIP* blocks as defined in registries.
    """
    tagset = set(tags)
    # Known canonical data starters
    candidates = {
        ">ZXXR", ">ZXXI", ">ZYYR", ">ZYYI", ">ZXYR", ">ZXYI", ">ZYXR", ">ZYXI",
        ">RHOXX", ">RHOXY", ">RHOYX", ">RHOYY",
        ">PHSXX", ">PHSXY", ">PHSYX", ">PHSYY",
        ">TXR.EXP", ">TYR.EXP", ">TIPMAG", ">TIPPHS",
    }

    # Expand registry-based candidates
    for group in _ZRP_COMPS:
        for row in group:
            for name in row:
                candidates.add(f">{name.upper()}")

    for row in _TIP_COMPS:
        for name in row:
            candidates.add(f">{name.upper()}")

    return any(t in tagset for t in candidates)



class IsEdi(ABC):
    """
    Abstract base class to assert SEG MT/EMAP EDI files or objects.

    This class provides a robust static validator :meth:`_assert_edi`
    for path-like EDI files and a polymorphic ``is_valid`` property
    (implemented by concrete EDI objects) so you can also validate
    *instances* via ``isinstance(obj, IsEdi)`` and ``obj.is_valid``.

    Notes
    -----
    * **Shallow** check (``deep=False``) only verifies the ``.edi``
      extension.
    * **Deep** check (``deep=True``) parses tags case-insensitively,
      tolerates extra whitespace/CRLF/BOM, and verifies:
        - presence of ``>HEAD`` and ``>END`` blocks
        - presence of a frequency block (``>FREQ``)
        - presence of a measurement/section header (e.g. ``>=MTSECT``
          or ``>=DEFINEMEAS``)
        - presence of at least one data block (any of ``Z*``, ``RHO*``,
          ``PHS*``, ``TIP*``)
        - optional sanity: ``NFREQ`` equals the number of numeric values
          found under the ``>FREQ`` payload.

    Examples
    --------
    Validate a file path:

    >>> IsEdi._assert_edi("E01.edi", deep=True)
    True

    Validate an object registered as a virtual subclass:

    >>> # IsEdi.register(Edi)  # somewhere in your code
    >>> isinstance(edi_obj, IsEdi) and edi_obj.is_valid
    True
    """

    @property
    @abstractmethod
    def is_valid(self) -> bool:
        """Whether this EDI object is valid."""
        raise NotImplementedError

    # ----------------------------
    # File / object level assertion
    # ----------------------------
    @staticmethod
    def _assert_edi(file: Union[str, os.PathLike, "IsEdi"],
                    deep: bool = True) -> bool:
        """
        Assert that the given path (or object) represents a valid EDI.

        Parameters
        ----------
        file : str or PathLike or IsEdi
            Path to a ``.edi`` file, or an EDI object implementing
            ``is_valid``.
        deep : bool, default=True
            When ``True``, fully parse and sanity-check the file
            contents. When ``False``, only check the file extension.

        Returns
        -------
        bool
            ``True`` if valid; otherwise raises an exception.

        Raises
        ------
        FileHandlingError
            If ``file`` is ``None`` or not a path/object we can handle.
        FileNotFoundError
            If the path does not exist.
        PermissionError
            If the file cannot be opened due to OS permissions.
        EdIDataError
            If the contents do not match an expected EDI structure.
        """
        seg_url = (
            "https://www.mtnet.info/docs/seg_mt_emap_1987.pdf"
        )
        base_msg = (
            "Unrecognized SEG EDI file. Refer to the SEG MT/EMAP "
            f"Data Interchange Standard: {seg_url}"
        )

        if file is None:
            raise FileHandlingError(
                "NoneType cannot be checked. Provide a valid path or EDI object."
            )

        # If an EDI object is passed, defer to its 'is_valid'
        if isinstance(file, IsEdi):
            if not file.is_valid:
                raise EdIDataError("EDI object reports invalid state.")
            return True

        # Path validation (accept str / PathLike)
        path = Path(file)  # type: ignore[arg-type]
        if not path.exists():
            raise FileNotFoundError(
                f"{str(path)!r} is not a file. Provide a path to an EDI file."
            )
        if not path.is_file():
            raise FileNotFoundError(
                f"{str(path)!r} is not a regular file. Provide a path to an EDI file."
            )

        suffix = path.suffix.lower().lstrip(".")
        if not deep:
            if suffix == "edi":
                return True
            raise EdIDataError(
                "SEG-EDI files normally use the '.edi' extension. "
                f"Got {path.suffix!r}. Set deep=True to validate contents."
            )

        # Deep parsing
        try:
            # Handle potential BOM with 'utf-8-sig'
            text = path.read_text(encoding="utf-8-sig", errors="replace")
        except PermissionError:
            help_url = (
                "https://stackoverflow.com/questions/36434764/"
                "permissionerror-errno-13-permission-denied"
            )
            raise PermissionError(
                "Permission denied while reading file. See guidance: "
                f"{help_url}"
            )

        # Normalize lines (keep raw to count numeric tokens)
        lines = text.splitlines()

        # Quick mandatory endpoints
        # (Be tolerant to leading BOM/whitespace and trailing whitespace.)
        # Find first and last meaningful tag
        first_tag = None
        for ln in lines:
            t = _extract_tag(ln)
            if t:
                first_tag = t
                break

        last_tag = None
        for ln in reversed(lines):
            t = _extract_tag(ln)
            if t:
                last_tag = t
                break

        if first_tag is None or last_tag is None:
            raise EdIDataError(base_msg)

        if not first_tag.startswith(">HEAD"):
            raise EdIDataError(
                f"{base_msg} (First block should be >HEAD, found {first_tag!r})"
            )
        if not last_tag.startswith(">END"):
            raise EdIDataError(
                f"{base_msg} (Last block should be >END, found {last_tag!r})"
            )

        # Collect all tag headers
        tags = _iter_blocks(lines)
        have_freq = any(t.startswith(">FREQ") for t in tags)
        have_section = any(
            t in {">=MTSECT", ">=DEFINEMEAS", ">=SPECTRASECT", ">=EMAPSECT", ">=OTHERSECT"}
            for t in tags
        )

        if not have_freq:
            raise EdIDataError(f"{base_msg} (Missing >FREQ block).")
        if not have_section:
            raise EdIDataError(
                f"{base_msg} (Missing a section header such as >=MTSECT or >=DEFINEMEAS)."
            )
        if not _has_any_data_block(tags):
            raise EdIDataError(
                f"{base_msg} (No data blocks found: expected Z*/RHO*/PHS*/TIP* entries)."
            )

        # Sanity: NFREQ vs actual FREQ payload count
        # Find the first FREQ header and validate
        for i, ln in enumerate(lines):
            tag = _extract_tag(ln)
            if tag and tag.startswith(">FREQ"):
                n_payload, _ = _count_freq_values(lines, i)
                n_expected = _expected_nfreq_from_header(ln)
                if n_expected is not None and n_expected != n_payload:
                    raise EdIDataError(
                        f"Inconsistent frequency count: header expects "
                        f"{n_expected}, but found {n_payload} values."
                    )
                break  # only check the first FREQ

        # If the file extension was not '.edi' but content is valid, still accept.
        return True
