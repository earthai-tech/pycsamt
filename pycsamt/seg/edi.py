# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Iterable, Optional, Tuple

import numpy as np
import datetime

from ..exceptions import (
    EdIDataError,
    FileHandlingError,
)
from ..core.base import CoreObject
from ..log.logger import get_logger
from ..z.resphase import ResPhase
from ..z.z import Z
from ..z.tipper import Tipper

from .validation import IsEdi
from .sections import iter_sections
from .spectra import SpectraSECT, SpectraIO, Spectra
from .time_series import TSIO, TimeSeries
from .other import OtherSECT, OtherIO
from .utils import _format_block_numbers

logger = get_logger(__name__)

__all__ = ["EDIMixin", "EDIOMixin", "EDIFile"]


class EDIMixin(CoreObject):
    r"""
    Lightweight registry and helpers used by EDI readers.

    The mixin stores parsed section objects (e.g. ``>HEAD``,
    ``>=MTSECT``, ``>=SPECTRASECT``) in a simple dictionary,
    and provides tiny utilities to manage and query them.

    It does **not** perform I/O.  Host classes remain free
    to decide *when* and *how* sections are discovered and
    populated.  This keeps the orchestration logic small and
    testable.

    Attributes
    ----------
    sections : dict[str, object]
        Case–insensitive mapping from a logical section key
        (e.g. ``"head"``, ``"mtsect"``, ``"spectra"``) to the
        parsed object that represents that section.

    Methods
    -------
    add_section(key, obj)
        Register ``obj`` under ``key``.  Keys are normalized
        to lower case.
    get_section(key)
        Retrieve a previously added section or ``None``.
    has_section(key)
        Return ``True`` if a section exists under ``key``.
    _tag2name(tag)
        Translate a raw ``>=...`` header tag to a canonical
        registry key (e.g. ``">=MTSECT" → "mtsect"``).

    Notes
    -----
    The registry is intentionally untyped so that different
    parser implementations can coexist.  For instance, a
    project may store a *header* object (``SpectraSECT``)
    and also the decoded *payload* object (``Spectra``) under
    separate keys.

    Examples
    --------
    >>> mix = EDIMixin()
    >>> mix._init_registry()
    >>> mix.add_section("HEAD", object())
    >>> mix.has_section("head")
    True
    >>> isinstance(mix.get_section("head"), object)
    True

    See Also
    --------
    EDIFile
        High level reader that uses the registry to expose
        parsed sections to callers.

    References
    ----------
    .. [1] SEG EDI MT/EMAP standard (1987).  MTNet.
    """
    def _init_registry(self) -> None:
        self.sections: Dict[str, Any] = {}

    def add_section(self, key: str, obj: Any) -> None:
        self.sections[str(key).lower()] = obj

    def get_section(self, key: str) -> Any:
        return self.sections.get(str(key).lower())

    def has_section(self, key: str) -> bool:
        return str(key).lower() in self.sections

    def _tag2name(self, tag: str) -> str:
        t = tag.upper()
        if t in {">=MTSECT", ">=EMAPSECT"}:
            return "mtsect"
        if t in {">=SPECTRA", ">=SPECTRASECT"}:
            return "spectra_sect"
        if t in {">=TSERIES", ">=TSERIESSECT"}:
            return "timeseries_sect"
        if t == ">=OTHERSECT":
            return "other"
        return t


class EDIOMixin(CoreObject):
    r"""
    Tolerant ``>BLOCK`` parser and TF (Z/Tipper) builder.

    The mixin provides two core utilities used by :class:`EDIFile`
    after headers are discovered:

    * :meth:`_scan_blocks` reads numeric blocks starting at a
      given line (e.g. ``>FREQ``, ``>ZXXR``, ``>TXR.EXP``).
    * :meth:`_build_from_comp` converts the scanned dictionary
      into :class:`pycsamt.z.z.Z` and
      :class:`pycsamt.z.tipper.Tipper` objects.

    The reader accepts both complex tensor blocks and the
    scalar families (``RHO*`` and ``PHS*``).  When complex
    blocks are missing, the impedance tensor is reconstructed
    from resistivity and phase if possible.


    Attributes
    ----------
    None

    Methods
    -------
    _scan_blocks(path, start=None, empty_val=1e32)
        Return a mapping ``key → list[float]`` by streaming
        lines until the next section or EOF.  Unknown keys are
        ignored.  Values equal to ``empty_val`` are converted
        to zeros.
    _build_from_comp(comp, z_obj, tip_obj)
        Populate ``z_obj`` and ``tip_obj`` from the components.
        Frequency order is normalized to *descending*.  Z–error
        arrays are converted from variance blocks (``.VAR``).

    Notes
    -----
    * Frequency order is unified to high→low.  This follows a
      common practice in EDI archives and simplifies plotting.
    * The tipper is optional.  If no tipper blocks are found,
      ``tip_obj`` is left untouched.
    * For ``RHO*/PHS*`` the method also carries ``*.ERR`` if
      present.  Otherwise zero errors are assumed.

    Examples
    --------
    >>> comp = {"freq": [10, 1], "zxxr": [1, 2], "zxxi": [0, 0]}
    >>> from pycsamt.z.z import Z
    >>> from pycsamt.z.tipper import Tipper
    >>> mix = EDIOMixin()
    >>> z, t = Z(), Tipper()
    >>> mix._build_from_comp(comp, z_obj=z, tip_obj=t)
    >>> z.n_freq
    2

    See Also
    --------
    EDIFile
        Uses these utilities during :meth:`read_data`.

    References
    ----------
    .. [1] SEG EDI MT/EMAP standard (1987).  MTNet.
    """
    # tolerant key sets (lower-case)
    _Z_PAIRS = ("zxx", "zxy", "zyx", "zyy")
    _RHO_PREF = ("rho", "frho")
    _PHS_PREF = ("phs", "fphs")

    _TIP_R = ("txr", "tyr")
    _TIP_I = ("txi", "tyi")
    _TIP_V = ("txvar", "tyvar")

    def _scan_blocks(
        self,
        path: Path,
        *,
        start: Optional[int] = None,
        empty_val: float = 1.0e32,
    ) -> Dict[str, list]:
        IsEdi._assert_edi(path, deep=True)

        with path.open("r", encoding="utf-8") as f:
            lines = f.readlines()

        i0 = 0 if start is None else int(start)
        if i0 <= 0:
            for i, ln in enumerate(lines):
                s = ln.lstrip()
                if s.startswith(">FREQ"):
                    i0 = i
                    break
        if i0 <= 0 or i0 >= len(lines):
            raise EdIDataError("data section not found")

        comp: Dict[str, list] = {}
        cur: Optional[str] = None

        def _is_known(key: str) -> bool:
            k = key.lower()
            if k in {"freq", "zrot", "trot"}:
                return True
            if any(
                k in (f"{c}r", f"{c}i", f"{c}.var")
                for c in self._Z_PAIRS
            ):
                return True
            if any(k.startswith(p) for p in self._RHO_PREF):
                return True
            if any(k.startswith(p) for p in self._PHS_PREF):
                return True
            base = k.split(".")[0]
            if base in {"txr", "txi", "tyr", "tyi"}:
                return True
            if base in {"txvar", "tyvar"}:
                return True
            return False

        for raw in lines[i0:]:
            s = raw.strip()
            if not s:
                continue

            if s.startswith(">!"):
                cur = None
                continue

            if s.startswith(">"):
                toks = s[1:].split()
                if not toks:
                    cur = None
                    continue
                key = toks[0].lower()
                if not _is_known(key): 
                    cur =None
                    continue
                
                base = key.split(".", 1)[0]
                # normalize tipper keys: TXR.EXP -> txr
                if base in {
                    "txr", "txi", "tyr", "tyi",
                    "txvar", "tyvar",
                }:
                    cur = base
                else:
                    cur = key

                comp.setdefault(cur, [])
                
                continue

            if cur is None:
                continue

            vals = []
            for tok in s.split():
                try:
                    v = float(tok)
                    if v == empty_val:
                        v = 0.0
                    vals.append(v)
                except Exception:
                    vals.append(0.0)
            comp[cur].extend(vals)

        return _lower_keys(comp)

    def _build_from_comp(
        self,
        comp: Dict[str, list],
        *,
        z_obj: Z,
        tip_obj: Tipper,
    ) -> None:
        if "freq" not in comp or not comp["freq"]:
            raise EdIDataError("FREQ block missing")

        f = np.asarray(comp["freq"], float)
        f, rev = _rev_if_asc(f)
        n = f.size

        # ---------------- Z via complex blocks
        z = np.zeros((n, 2, 2), complex)
        zv = np.zeros((n, 2, 2), float)

        def _arr(name: str) -> np.ndarray:
            a = np.asarray(comp.get(name, []), float)
            return _nz(a, n)

        def _fill(cc: str, i: int, j: int) -> None:
            r = _arr(f"{cc}r")
            im = _arr(f"{cc}i")
            vr = _arr(f"{cc}.var")
            z[:, i, j] = r + 1j * im
            zv[:, i, j] = np.sqrt(np.abs(vr))

        for cc, ij in zip(
            self._Z_PAIRS,
            ((0, 0), (0, 1), (1, 0), (1, 1)),
        ):
            _fill(cc, *ij)

        if rev:
            z = z[::-1]
            zv = zv[::-1]

        # ---------------- Z via ρ/φ (if Z blocks absent)
        have_z = np.any(np.abs(z) > 0)
        have_rp = any(
            any(k.startswith(p) for k in comp)
            for p in (*self._RHO_PREF, *self._PHS_PREF)
        )

        if (not have_z) and have_rp:
            rho = np.zeros((n, 2, 2), float)
            phi = np.zeros((n, 2, 2), float)
            rho_e = None
            phi_e = None

            def _get2(prefix: str, suf: str) -> np.ndarray:
                # e.g. "rho", "xy" -> try RHOXY / FRHOXY
                for p in (self._RHO_PREF if
                          prefix == "rho" else self._PHS_PREF):
                    key = f"{p}{suf}".lower()
                    if key in comp:
                        a = np.asarray(comp[key], float)
                        return _nz(a, n)
                return np.zeros(n, float)

            def _gete2(prefix: str, suf: str) -> np.ndarray:
                # ERR or .ERR forms
                for p in (self._RHO_PREF if
                          prefix == "rho" else self._PHS_PREF):
                    for tail in (".err", "err"):
                        key = f"{p}{suf}{tail}".lower()
                        if key in comp:
                            a = np.asarray(comp[key], float)
                            return _nz(a, n)
                return np.zeros(n, dtype=float)

            # XX, XY, YX, YY
            for suf, ij in zip(
                ("xx", "xy", "yx", "yy"),
                ((0, 0), (0, 1), (1, 0), (1, 1)),
            ):
                rho[:, ij[0], ij[1]] = _get2("rho", suf)
                phi[:, ij[0], ij[1]] = _get2("phs", suf)

            re = sum(
                np.any(_gete2("rho", s))
                for s in ("xx", "xy", "yx", "yy")
            )
            pe = sum(
                np.any(_gete2("phs", s))
                for s in ("xx", "xy", "yx", "yy")
            )

            if re or pe:
                rho_e = np.zeros_like(rho)
                phi_e = np.zeros_like(phi)
                for suf, ij in zip(
                    ("xx", "xy", "yx", "yy"),
                    ((0, 0), (0, 1), (1, 0), (1, 1)),
                ):
                    rho_e[:, ij[0], ij[1]] = _gete2("rho", suf)
                    phi_e[:, ij[0], ij[1]] = _gete2("phs", suf)

            if rev:
                rho = rho[::-1]
                phi = phi[::-1]
                if rho_e is not None:
                    rho_e = rho_e[::-1]
                if phi_e is not None:
                    phi_e = phi_e[::-1]

            z_obj.set_res_phase(
                rho,
                phi,
                f,
                res_err_array=rho_e,
                phase_err_array=phi_e,
            )
        else:
            z_obj._freq = f
            z_obj._z = z
            z_obj._z_err = zv
            z_obj.compute_resistivity_phase()

        zr = np.asarray(comp.get("zrot", []), float)
        zr = _nz(zr, n)
        z_obj.rotation_angle = zr[::-1] if rev else zr

        # Tipper (optional)
        have_tip = any(k.startswith(("tx", "ty")) for k in comp)
        if not have_tip:
            return

        txr = _nz(np.asarray(comp.get("txr", []), float), n)
        txi = _nz(np.asarray(comp.get("txi", []), float), n)
        tyr = _nz(np.asarray(comp.get("tyr", []), float), n)
        tyi = _nz(np.asarray(comp.get("tyi", []), float), n)

        txv = comp.get("txvar", [])
        tyv = comp.get("tyvar", [])
        txv = _nz(np.asarray(txv, float), n)
        tyv = _nz(np.asarray(tyv, float), n)

        t = np.zeros((n, 1, 2), complex)
        te = np.zeros((n, 1, 2), float)
        t[:, 0, 0] = txr + 1j * txi
        t[:, 0, 1] = tyr + 1j * tyi
        te[:, 0, 0] = np.sqrt(np.abs(txv))
        te[:, 0, 1] = np.sqrt(np.abs(tyv))

        if rev:
            t = t[::-1]
            te = te[::-1]

        tip_obj._freq = f
        tip_obj._tipper = t
        tip_obj._tipper_err = te

        tr = np.asarray(comp.get("trot", []), float)
        tr = _nz(tr, n)
        if tr.size == n:
            tip_obj.rotation_angle = tr[::-1] if rev else tr
        else:
            tip_obj.rotation_angle = z_obj.rotation_angle

        tip_obj.compute_amp_phase()
        tip_obj.compute_mag_direction()


class EDIFile(EDIMixin, EDIOMixin):
    r"""
    High–level EDI dispatcher for SEG/EMAP/CSAMT archives.

    The class discovers top–level headers (e.g. ``>=MTSECT``,
    ``>=SPECTRASECT``, ``>=TSERIESSECT``), loads the matching
    data blocks, and exposes convenient Python containers for
    impedance tensors, tippers, spectra, and time series.

    It also writes EDI files by reusing the in–memory sections,
    preserving headers when possible, and regenerating block
    payloads from the current objects.

    Parameters
    ----------
    path : str or Path, optional
        File to open.  If given, :meth:`read` is executed on
        construction.
    verbose : int, optional
        Verbosity level propagated to subcomponents.

    Attributes
    ----------
    path : Path or None
        Bound file path (if any).
    Z : pycsamt.z.z.Z
        Impedance tensor container with errors and rotations.
    Tip : pycsamt.z.tipper.Tipper
        Tipper container (optional).
    sections : dict[str, object]
        Registry populated via :class:`EDIMixin`.
    block_size : int
        Numbers per line when writing numeric payloads.
    float_fmt : str
        Float formatter used for numeric blocks.
    header_tpl : str
        Template used for logical block titles.

    Properties
    ----------
    station : str or None
        Shortcut to ``>HEAD.DATAID``.  Setter also mirrors the
        value to the MT/EMAP section id.
    processingsoftware : str or None
        Shortcut to the name of the processing software from
        ``>INFO``.

    Methods
    -------
    read(path=None)
        Load headers and sections.  Parse numeric blocks into
        :attr:`Z` and :attr:`Tip`.  Also attaches ``Spectra`` and
        ``TimeSeries`` if present.
    read_data()
        Low–level numeric parsing used by :meth:`read`.
    compose_headers()
        Serialize only the headers (no data blocks).
    write(...)
        Write a full EDI assembled from current objects and
        sections.  MT or EMAP numeric families are chosen from
        the MT/EMAP header or inferred from context.
    interpolate(new_freq, kind="slinear", ...)
        Interpolate ``Z`` on a new frequency grid.  The grid is
        rounded to two decimals for stable serialization.
    write_new_edi(edi_fn=None, Z=None, Tipper=None, ...)
        Rebuild a clean container bound to the same source, swap
        selected transfer functions, and delegate to
        :meth:`write`.

    Notes
    -----
    * Frequency order is normalized to descending on read.
      Therefore, a file written by :meth:`write` and read back
      will expose ``Z.freq`` in high→low order.
    * The interpolation routine enforces the new grid to live
      strictly inside the source span when ``bounds_error`` is
      ``True``.
    * Missing blocks are handled gracefully.  If complex
      impedances are absent, the reader tries to reconstruct
      them from ``RHO*/PHS*`` families.

    Examples
    --------
    >>> ed = EDIFile("site.edi")
    >>> ed.station
    'SITE'
    >>> ed.Z.n_freq > 0
    True
    >>> out = ed.write(savepath="outdir")
    >>> Path(out).exists()
    True
    >>> fnew = np.geomspace(ed.Z.freq.min()*1.1,
    ...                     ed.Z.freq.max()*0.9, 16)
    >>> z2 = ed.interpolate(fnew, kind="linear")
    >>> ed.write_new_edi(edi_fn="interp.edi", Z=z2)

    See Also
    --------
    EDIMixin
        Registry and convenience helpers used internally.
    EDIOMixin
        Numeric block reader and TF builder used by
        :meth:`read_data`.
    pycsamt.seg.spectra.Spectra, pycsamt.seg.time_series.TimeSeries
        High–level containers for spectra and time–series data.

    References
    ----------
    .. [1] SEG EDI MT/EMAP standard (1987).  MTNet.
    .. [2] B. Groom, R. Bailey (1989).  Decomposition of the
           magnetotelluric impedance tensor.  *Geophysics*.
    """
    def __init__(
        self,
        path: str | Path | None = None,
        *,
        verbose: int = 0,
    ) -> None:
        self.path: Optional[Path] = _as_path(path)
        self.verbose = int(verbose)

        self._init_registry()
        self._data_start: Optional[int] = None

        self.Z = Z(name=None, verbose=verbose)
        self.Res = ResPhase(name=None, verbose=verbose) 
        self.Tip = Tipper()
        
        # writer configuration
        self.block_size = 6
        self.float_fmt = "{: .6E}"      # scientific, uppercase E
        self.header_tpl = ">!****{title}****!\n"


        if self.path is not None:
            self.read()
            
    def read(
        self,
        path: str | Path | None = None,
    ) -> "EDIFile":
        if path is not None:
            self.path = _as_path(path)
        if self.path is None:
            raise FileHandlingError("path is not set")
        IsEdi._assert_edi(self.path, deep=True)
    
        for mod, key in (
            (".heads", "head"),
            (".heads", "info"),
            (".meas", "definemeas"),
        ):
            try:
                pkg = __import__(
                    f"pycsamt.seg{mod}",
                    fromlist=["*"],
                )
                try:
                    cls = getattr(pkg, key.capitalize())
                except AttributeError:
                    cls = getattr(pkg, "DefineMeas")
                obj = cls.from_file(self.path)
                self.add_section(key, obj)
            except:
                pass
    
        has_tf = False
        head = self.get_section("head")
        if head is not None:
            try:
                self.chainage = getattr(head, "chainage")
            except:
                pass

        for tag, header, j in iter_sections(str(self.path)):
            name = self._tag2name(tag)
            self.add_section(name, header)
    
            if name == "mtsect":
                has_tf = True
                self._data_start = j
    
            if isinstance(header, OtherSECT):
                try:
                    io = OtherIO.from_file(
                        str(self.path),
                        start_line=header.start_data_lines_num,
                    )
                    self.add_section("otherio", io)
                except Exception as exc:
                    logger.debug("OTHERIO skip: %s", exc)
    
            if name == "spectra_sect":
                # header from iter_sections may be minimal → enrich
                try:
                    need_ids = not getattr(header, "meas_ids", None)
                    need_map = not getattr(header, "id_to_chtype", None)
                    if need_ids or need_map:
                        rich = SpectraSECT.from_file(str(self.path))
                        header = rich
                except Exception:
                    pass
            
                self.add_section("spectra_sect", header)
      
                spec_obj = None
                try:
                    io = SpectraIO.from_file(
                        str(self.path),
                        start_line=header.start_data_lines_num,
                    )
                    self.add_section("spectra_io", io)
                
                    spec_obj = Spectra.from_io(
                        header,
                        io,
                        empty=self._sentinel(),
                        verbose=self.verbose,
                    )
                except Exception as exc:
                    logger.debug("SPECTRA load failed: %s", exc)
                
       
                # Fallback: expose a non-empty Spectra container only if needed
                if spec_obj is None:
                    spec_obj = Spectra(verbose=self.verbose)
                    sid = getattr(header, "sectid", None)
                    if sid:
                        # keep name for reporting
                        spec_obj.name = sid
                
                if getattr(spec_obj, "name", None) in (None, "", " "):
                    # prefer SPECTRA.SECTION id → >HEAD.DATAID → filename stem
                    head = self.get_section("head")
                    nm = (
                        getattr(header, "sectid", None)
                        or (getattr(head, "dataid", None) if head else None)
                        or (self.path.stem if self.path else None)
                    )
                    if nm:
                        try:
                            spec_obj.name = str(nm)
                        except Exception:
                            pass

                self.add_section("spectra", spec_obj)
                
                # assert/log once to catch this class of issues
                if getattr(spec_obj, "chan_ids", None) in (None, []):
                    logger.debug(
                        "Spectra chan_ids empty; check SPECTRASECT.")
                    
                # --- FALLBACK: build Z/Tipper from Spectra if no MT blocks
                if (not has_tf) and getattr(spec_obj, "n_freq", 0) > 0:
                    try:
                        z_from_sp, tip_from_sp = spec_obj.to_Z(estimate_error=False)
                        if z_from_sp is not None:
                            self.Z = z_from_sp
                            # ensure rho/phi are available downstream
                            self.Z.compute_resistivity_phase()
                            # we *do not* mark has_tf here; we only read raw
                            # MT blocks when MTSECT was found.
                        if tip_from_sp is not None:
                            self.Tip = tip_from_sp
                    except Exception as exc:
                        logger.debug("Z from SPECTRA failed: %s", exc)
            
            # attach TSERIES data blocks
            if name == "timeseries_sect":
                # keep the header object for metadata
                self.add_section("timeseries_sect", header)
            
                ts_obj = None
                try:
                    ts_io = TSIO.from_file(
                        str(self.path),
                        start_line=header.start_data_lines_num,
                    )
                    self.add_section("timeseries_io", ts_io)
            
                    # build high-level container from parsed IO
                    ts_obj = TimeSeries.from_io(header, ts_io)
                except Exception as exc:
                    logger.debug("TSERIES load failed: %s", exc)
            
                # fallback: still expose a TimeSeries container so
                # callers have a uniform API even with no blocks.
                if ts_obj is None:
                    try:
                        # if you have a helper:
                        #   ts_obj = TimeSeries.from_header(header)
                        # otherwise synthesize a minimal container:
                        ts_obj = TimeSeries(verbose=self.verbose)
                        sid = getattr(header, "sectid", None)
                        if sid:
                            setattr(ts_obj, "sectid", sid)
                        dt = getattr(header, "dt", None)
                        if dt is not None:
                            # keep section-level dt for time() defaults
                            setattr(ts_obj, "_sect_dt", float(dt))
                    except Exception:
                        # ultimate guardrail
                        ts_obj = TimeSeries(verbose=self.verbose)
            
                self.add_section("timeseries", ts_obj)

        # if has_tf:
        #     self.read_data()
        # only scan MT numeric blocks when we actually found them
        if self._data_start is not None:
            self.read_data()
        
        return self

    def read_data(self) -> "EDIFile":
        if self.path is None:
            raise FileHandlingError("path is not set")

        empty_val = 1.0e32
        head = self.get_section("head")
        if head is not None:
            ev = getattr(head, "empty", None)
            try:
                if ev is not None:
                    empty_val = float(ev)
            except Exception:
                pass

        comp = self._scan_blocks(
            self.path,
            start=self._data_start,
            empty_val=empty_val,
        )
        self._build_from_comp(
            comp,
            z_obj=self.Z,
            tip_obj=self.Tip,
        )
 
        # After Z is built, create and populate 
        # the ResPhase object for API consistency.
        if self.Z.freq is not None and self.Z.z is not None:
            station_name = getattr(
                head, "dataid", None
            ) if head else self.station
            
            # Instantiate and populate the ResPhase container
            self.Res = ResPhase(
                freq=self.Z.freq, name=station_name
                )
            try:
                self.Res.compute_resistivity_phase(
                    z_array=self.Z.z,
                    z_err_array=self.Z.z_err,
                    freq=self.Z.freq,
                )
            except Exception as e:
                logger.warning(
                    f"Could not compute resistivity/phase"
                    f" for {station_name}: {e}"
                )
   
        return self

    def compose_headers(self) -> str:
        out: list[str] = []
        for key in (
            "head",
            "info",
            "definemeas",
            "mtsect",
            "spectra",
            "timeseries",
            "other",
        ):
            sec = self.get_section(key)
            if sec is None:
                continue
            write = getattr(sec, "write", None)
            if callable(write):
                try:
                    out.extend(write())
                except Exception:
                    pass
        return "".join(out)
 
    def write(  # noqa: C901
        self,
        edi_fn: str | None = None,
        new_edifn: str | None = None,
        datatype: str | None = None,
        savepath: str | Path | None = None,
        add_filter_array: np.ndarray | None = None,
        synthesize_spectra: bool = False,  
        **kwargs,
    ) -> str:
        # verbosity: prefer kwarg, else instance value
        v = kwargs.pop("verbose", None)
        if v is not None:
            self.verbose = int(v)
  
        # v2 auto-detects MT vs EMAP from the section header.
        def _detect_tf_mode() -> str | None:
            m = self.get_section("mtsect")
            if m is None:
                return None
            try:
                hl = m.write()
                if hl:
                    f0 = str(hl[0]).upper()
                    if "MTSECT" in f0:
                        return "mt"
                    if "EMAPSECT" in f0:
                        return "emap"
            except Exception:
                pass
            # fallback: if tipper present assume MT, else EMAP
            tip = getattr(self.Tip, "tipper", None)
            if tip is not None and tip.size and not np.all(tip == 0.0):
                return "mt"
            return "emap"
    
        # explicit override if user passed a value
        dt_auto = _detect_tf_mode()
        dtype = (
            (datatype or dt_auto or "mt").strip().lower()
        )
    
        out_dir = (
            Path(savepath).expanduser().resolve()
            if savepath is not None
            else (Path.cwd() / "edi_out")
        )
        out_dir.mkdir(parents=True, exist_ok=True)
    
        # choose the output file name
        # - new_edifn has priority (ensure .edi)
        # - else derive from source path
        # - else synthesize from HEAD.dataid + dtype + year
        src = edi_fn or getattr(self, "path", None)
        year = datetime.datetime.utcnow().year
    
        if new_edifn:
            name = str(new_edifn)
            if not name.lower().endswith(".edi"):
                name += ".edi"
        elif isinstance(src, (str, Path)):
            name = f"new_{Path(src).name}"
        else:
            head = self.get_section("head")
            sid = getattr(head, "dataid", None) if head else None
            sid = sid or "site"
            name = f"{sid}_{dtype}.{year}.edi"
    
        out_path = out_dir / name
    
        lines: list[str] = []
    
        # 1) structural section headers (safe compose)
        blob = self.compose_headers()
        if blob:
            lines.append(blob)
    
        # helper: append string or list[str]
        def _append(b: list[str] | str | None) -> None:
            if not b:
                return
            if isinstance(b, str):
                lines.append(b)
            else:
                lines.extend(b)
    
        # helper: emit >FREQ and get ZROT/TROT & tipper flag
        def _emit_freq_and_rot() -> tuple[
            np.ndarray, np.ndarray, np.ndarray, bool
        ]:
            if self.Z.freq is None or self.Z.n_freq == 0:
                raise EdIDataError(
                    "no frequency vector for MT/EMAP"
                )
            f = np.asarray(self.Z.freq, float)
            _append(self.header_tpl.format(title="FREQUENCIES"))
            _append(self._emit_block("FREQ", f, rot=None))
    
            zrot = getattr(self.Z, "rotation_angle", None)
            if np.ndim(zrot) == 0:
                zrot = np.full(f.size, float(zrot or 0.0))
            else:
                zrot = np.asarray(zrot, float)
                if zrot.size != f.size:
                    zrot = np.full(f.size, 0.0, dtype=float)
    
            tip = getattr(self.Tip, "tipper", None)
            tip_ok = (
                tip is not None and tip.size and not np.all(tip == 0.0)
            )
            trot = getattr(self.Tip, "rotation_angle", None)
            if tip_ok:
                if np.ndim(trot) == 0:
                    trot = np.full(f.size, float(trot or 0.0))
                else:
                    trot = np.asarray(trot, float)
                    if trot.size != f.size:
                        trot = np.full(f.size, 0.0, dtype=float)
            else:
                trot = np.zeros_like(f)
    
            return f, zrot, trot, bool(tip_ok)
    
        # 2) MT/EMAP numeric blocks, only if we have MT/EMAP
        if self.has_section("mtsect") or self.Z.n_freq > 0:
            freq, zrot, trot, tip_ok = _emit_freq_and_rot()
    
            z = np.nan_to_num(self.Z.z)
            ze = (
                np.nan_to_num(self.Z.z_err)
                if (self.Z.z_err is not None)
                else np.zeros_like(z.real)
            )
    
            if dtype == "mt":
                _append(
                    self.header_tpl.format(
                        title="IMPEDANCE ROTATION ANGLES"
                    )
                )
                _append(self._emit_block("ZROT", zrot))
    
            _append(self.header_tpl.format(title="IMPEDANCES"))
            
            for tag, i, j in (
                ("ZXX", 0, 0),
                ("ZXY", 0, 1),
                ("ZYX", 1, 0),
                ("ZYY", 1, 1),
            ):
                rot_tag = "ROT=ZROT" if dtype == "mt" else "ROT=NONE"
                _append(
                    self._emit_block(
                        f"{tag}R",
                        self._sanitize_payload(z[:, i, j].real),
                        rot=rot_tag,
                    )
                )
                _append(
                    self._emit_block(
                        f"{tag}I",
                        self._sanitize_payload(z[:, i, j].imag),
                        rot=rot_tag,
                    )
                )
                _append(
                    self._emit_block(
                        f"{tag}.VAR",
                        self._sanitize_payload(
                            np.square(ze[:, i, j])
                        ),
                        rot=rot_tag,
                    )
                )
    
            if dtype == "emap" or dtype =="mt":
                _append(
                    self.header_tpl.format(
                        title="RESISTIVITIES AND PHASES"
                    )
                )
    
                def _rho_phi(
                    tag: str,
                    a: np.ndarray,
                    e: np.ndarray | None,
                ) -> None:
                    err = np.zeros_like(a) if e is None else e
                    _append(
                        self._emit_block(
                            tag,
                            self._sanitize_payload(a),
                            rot="ROT=NONE",
                        )
                    )
                    _append(
                        self._emit_block(
                            f"{tag}.ERR",
                            self._sanitize_payload(err),
                            rot="ROT=NONE",
                        )
                    )
    
                _rho_phi("RHOXX", self.Z.res_xx, self.Z.res_err_xx)
                _rho_phi("RHOXY", self.Z.res_xy, self.Z.res_err_xy)
                _rho_phi("RHOYX", self.Z.res_yx, self.Z.res_err_yx)
                _rho_phi("RHOYY", self.Z.res_yy, self.Z.res_err_yy)
    
                _rho_phi(
                    "PHSXX", self.Z.phase_xx, self.Z.phase_err_xx
                )
                _rho_phi(
                    "PHSXY", self.Z.phase_xy, self.Z.phase_err_xy
                )
                _rho_phi(
                    "PHSYX", self.Z.phase_yx, self.Z.phase_err_yx
                )
                _rho_phi(
                    "PHSYY", self.Z.phase_yy, self.Z.phase_err_yy
                )
    
                if add_filter_array is not None:
                    a = np.asarray(add_filter_array, float)
                    ok = (
                        a.ndim == 3
                        and a.shape[1:] == (2, 2)
                        and a.shape[0] == freq.size
                    )
                    if ok:
                        _append(
                            self._emit_block(
                                "FRHOXY",
                                self._sanitize_payload(a[:, 0, 1]),
                                rot="ROT=NONE",
                            )
                        )
                        _append(
                            self._emit_block(
                                "FRHOYX",
                                self._sanitize_payload(a[:, 1, 0]),
                                rot="ROT=NONE",
                            )
                        )
    
            if tip_ok and dtype == "mt":
                _append(
                    self.header_tpl.format(
                        title="TIPPER ROTATION ANGLES"
                    )
                )
                _append(self._emit_block("TROT", trot))
    
                _append(
                    self.header_tpl.format(
                        title="TIPPER PARAMETERS"
                    )
                )
                tip = np.asarray(self.Tip.tipper)
                terr = getattr(self.Tip, "_tipper_err", None)
                tvar = (
                    np.square(terr)
                    if terr is not None
                    else np.zeros_like(tip.real)
                )
                for idx, tag in enumerate(("TX", "TY")):
                    _append(
                        self._emit_block(
                            f"{tag}R.EXP",
                            self._sanitize_payload(
                                tip[:, 0, idx].real
                            ),
                            rot="ROT=TROT",
                        )
                    )
                    _append(
                        self._emit_block(
                            f"{tag}I.EXP",
                            self._sanitize_payload(
                                tip[:, 0, idx].imag
                            ),
                            rot="ROT=TROT",
                        )
                    )
                    _append(
                        self._emit_block(
                            f"{tag}VAR.EXP",
                            self._sanitize_payload(
                                tvar[:, 0, idx]
                            ),
                            rot="ROT=TROT",
                        )
                    )
    
        # 3) SPECTRA: prefer parsed header/io; else synthesize
        spec_sect = self.get_section("spectra_sect")
        spec_io = self.get_section("spectra_io")
        spec_obj = self.get_section("spectra")
    
        # synthesize when asked and no spectra loaded
        if (spec_io is None and spec_obj is None and synthesize_spectra
                and self.Z.n_freq > 0):
            try:
                sp = Spectra.from_Z(
                    self.Z,
                    include_hz=self.has_tipper,
                    tipper=self.Tip if self.has_tipper else None,
                    name=(getattr(self.get_section("head"), "dataid", None)
                          or getattr(self.Z, "name", None)),
                    verbose=self.verbose,
                )
                spec_obj = sp
                spec_sect, spec_io = sp.to_io()
                # register so compose_headers / later code can use them
                self.add_section("spectra", sp)
                if spec_sect is not None:
                    self.add_section("spectra_sect", spec_sect)
                if spec_io is not None:
                    self.add_section("spectra_io", spec_io)
            except Exception as exc:
                logger.debug("synth spectra failed: %s", exc)
        
        # existing behavior (write parsed spectra if present)
        
        if spec_sect is None and spec_io is None and spec_obj:
            try:
                spec_sect, spec_io = spec_obj.to_io()
            except Exception as exc:
                logger.debug("SPECTRA synth failed: %s", exc)
            if spec_sect is not None:
                _append(spec_sect.write())
    
        if spec_io is not None:
            try:
                _append(spec_io.write())
            except Exception as exc:
                logger.warning("Skip spectra write: %s", exc)
    
        # 4) TSERIES: prefer parsed header/io; else synthesize
        ts_sect = self.get_section("timeseries_sect")
        ts_io = self.get_section("timeseries_io")
        ts_obj = self.get_section("timeseries")
    
        if ts_sect is None and ts_io is None and ts_obj:
            try:
                ts_sect, ts_io = ts_obj.to_io()
            except Exception as exc:
                logger.debug("TSERIES synth failed: %s", exc)
            if ts_sect is not None:
                _append(ts_sect.write())
    
        if ts_io is not None:
            try:
                _append(ts_io.write())
            except Exception as exc:
                logger.warning("Skip time series write: %s", exc)
    
        # 5) OTHER: write back verbatim if captured
        other_sect = self.get_section("other")
        other_io = self.get_section("otherio")
        if other_sect is not None:
            try:
                _append(other_sect.write())
            except Exception:
                pass
        if other_io is not None:
            try:
                _append(other_io.write())
            except Exception:
                pass
    
        # 6) end marker + flush
        lines.append(">END")
        with out_path.open("w", encoding="utf-8") as fw:
            fw.writelines(lines)
    
        if getattr(self, "verbose", 0) > 0:
            logger.info("EDI written: %s", str(out_path))
    
        return str(out_path)

    def _sentinel(self) -> float:
        # pull EMPTY from >HEAD if present, else 1e32
        empty = 1.0e32
        head = self.get_section("head")
        if head is not None:
            ev = getattr(head, "empty", None)
            try:
                if ev is not None:
                    empty = float(ev)
            except Exception:
                pass
        return float(empty)
    
    
    def _sanitize_payload(
        self,
        arr: np.ndarray,
        *,
        replace_zero: bool = True,
    ) -> np.ndarray:
        # replace 0/NaN/Inf by EMPTY sentinel for data payloads
        out = np.asarray(arr, dtype=float).copy()
        if not replace_zero:
            return out
        sentinel = self._sentinel()
        bad = (~np.isfinite(out)) | (out == 0.0)
        out[bad] = sentinel
        return out


    def _emit_block(
        self,
        key: str,
        values: np.ndarray,
        *,
        rot: str | None = None,
    ) -> list[str]:
        # one header + wrapped numeric payload
        v = np.asarray(values, dtype=float).ravel()
        head = f">{key.upper()}"
        if rot:
            head += f" {rot}"
        head += f"  //{v.size}\n"
        body = _format_block_numbers(
            v,
            per_line=self.block_size,
            fmt=self.float_fmt,
            indent=2,
        )
        return [head, body, "\n"]

    
    def interpolate(
        self,
        new_freq: np.ndarray | list[float],
        *,
        kind: str = "slinear",
        bounds_error: bool = True,
        period_buffer: float | None = None,
    ) -> Z:
        # interpolate Z onto a new frequency grid. errors are
        # interpolated independently and carried over.
        try:
            from scipy.interpolate import interp1d
        except Exception as exc:  # pragma: no cover
            raise ImportError(
                "SciPy is required for interpolation."
            ) from exc
    
        nf = np.asarray(new_freq, float)
        nf = np.around(nf, 2)
    
        of = np.asarray(self.Z.freq, float)
        if of.ndim != 1 or of.size == 0:
            raise EdIDataError("Z.freq is empty or invalid.")
    
        if bounds_error:
            if nf.min() < of.min() or nf.max() > of.max():
                raise ValueError(
                    "new frequency range must lie within src."
                )
    
        n_new = int(nf.size)
        out = Z(
            z_array=np.zeros((n_new, 2, 2), complex),
            z_err_array=np.zeros((n_new, 2, 2), float),
            freq=nf,
        )
    
        for i in range(2):
            for j in range(2):
                src = np.asarray(self.Z.z[:, i, j])
                err = np.asarray(
                    self.Z.z_err[:, i, j]
                    if self.Z.z_err is not None
                    else np.zeros_like(src.real)
                )
    
                nz = np.nonzero(src)[0]
                if nz.size == 0:
                    continue
    
                f = of[nz]
                zr = src[nz].real
                zi = src[nz].imag
                ze = err[nz]
    
                keep = (nf >= f.min()) & (nf <= f.max())
                idx = np.where(keep)[0]
                if idx.size == 0:
                    continue
                fn = nf[idx]
    
                if isinstance(period_buffer, (int, float)):
                    # accept new points close (in period) to
                    # existing samples. ratio < (pb + 1).
                    pb = float(period_buffer) + 1.0
                    sel: list[int] = []
                    for k, ff in enumerate(fn):
                        dif = np.abs(np.log10(ff) - np.log10(f))
                        k0 = int(np.argmin(dif))
                        r = max(f[k0] / ff, ff / f[k0])
                        if r < pb:
                            sel.append(idx[k])
                    idx = np.array(sel, int)
                    if idx.size == 0:
                        continue
                    fn = nf[idx]
    
                order = np.argsort(f)
                f2 = f[order]
                fr = interp1d(
                    f2, zr[order], kind=kind, bounds_error=True
                )
                fi = interp1d(
                    f2, zi[order], kind=kind, bounds_error=True
                )
                fe = interp1d(
                    f2, ze[order], kind=kind, bounds_error=True
                )
    
                out.z[idx, i, j] = fr(fn) + 1j * fi(fn)
                out.z_err[idx, i, j] = fe(fn)
    
        out.compute_resistivity_phase()
        return out
    
    
    # backward-compat
    interpolate_z = interpolate
    
    def write_new_edi(  # noqa: D401
        self,
        edi_fn: str | None = None,
        Z: Z | None = None,
        Tipper: Tipper | None = None,
        *,
        Spectra: object | None = None,
        TimeSeries: object | None = None,
        sections: dict[str, object] | None = None,
        **kwargs,
    ) -> str:
        # rebuild a fresh container bound to same source file and
        # sections. swap transfer functions or aux sections if
        # provided. then delegate to write().
        if self.path is None:
            raise EdIDataError("no source EDI bound; abort.")
    
        fresh = self.__class__(self.path, verbose=self.verbose)
    
        # keep structural sections already parsed by fresh.read().
        # allow targeted overrides via `sections` mapping.
        if sections:
            for k, v in sections.items():
                fresh.add_section(k, v)
    
        # swap TFs (Z/Tipper)
        fresh.Z = Z if Z is not None else self.Z
        fresh.Tip = Tipper if Tipper is not None else self.Tip
    
        # optional high-level spectra / time-series overrides
        if Spectra is not None:
            fresh.add_section("spectra", Spectra)
            try:
                sect, io = Spectra.to_io()  # type: ignore[attr-defined]  # noqa: E501
                fresh.add_section("spectra_sect", sect)
                fresh.add_section("spectra_io", io)
            except Exception:
                # tolerate objects that don't expose .to_io()
                pass
    
        if TimeSeries is not None:
            fresh.add_section("timeseries", TimeSeries)
            try:
                tsect, tio = TimeSeries.to_io()  # type: ignore[attr-defined]  # noqa: E501
                fresh.add_section("timeseries_sect", tsect)
                fresh.add_section("timeseries_io", tio)
            except Exception:
                pass
    
        # forward to main writer; map to new_edifn param
        return fresh.write(new_edifn=edi_fn, **kwargs)
    

    @property
    def n_freq(self) -> int:
        """Number of frequencies in Z or Spectra."""
        n = int(getattr(self.Z, "n_freq", 0) or 0)
        if n:
            return n
        sp = self.get_section("spectra")
        return int(getattr(sp, "n_freq", 0) or 0)

    @property
    def station(self) -> str | None:
        """Return DATAID from >HEAD if present."""
        h = self.get_section("head")
        return getattr(h, "dataid", None) if h else None

    @station.setter
    def station(self, value: str) -> None:
        """Set DATAID and keep section IDs aligned."""
        name = str(value)
        h = self.get_section("head")
        if h is not None:
            try:
                h.dataid = name
            except Exception:
                pass

        # keep >=... section SECTID aligned when possible
        for key in (
            "mtsect",
            "spectra",
            "timeseries",
            "spectra_sect",
            "timeseries_sect",
        ):
            obj = self.get_section(key)
            if obj is not None and hasattr(obj, "sectid"):
                try:
                    setattr(obj, "sectid", name)
                except Exception:
                    pass

    @property
    def empty(self) -> float | None:
        """Return EMPTY sentinel from >HEAD if present."""
        h = self.get_section("head")
        try:
            v = getattr(h, "empty", None)
            return float(v) if v is not None else None
        except Exception:
            return None

    @empty.setter
    def empty(self, value: float) -> None:
        """Set EMPTY sentinel on >HEAD if present."""
        h = self.get_section("head")
        if h is None:
            return
        try:
            h.empty = float(value)
        except Exception:
            pass

    @property
    def dtype(self) -> str | None:
        """Infer 'mt' or 'emap' from >=MT/EMAPSECT or tipper."""
        m = self.get_section("mtsect")
        if m is not None:
            try:
                hl = m.write()
                if hl:
                    f0 = str(hl[0]).upper()
                    if "MTSECT" in f0:
                        return "mt"
                    if "EMAPSECT" in f0:
                        return "emap"
            except Exception:
                pass
        tip = getattr(self.Tip, "tipper", None)
        if tip is not None and getattr(tip, "size", 0):
            if not np.all(tip == 0.0):
                return "mt"
        return "emap" if self.has_section("mtsect") else None

    @property
    def has_tipper(self) -> bool:
        """True if non-zero tipper array present."""
        tip = getattr(self.Tip, "tipper", None)
        if tip is None:
            return False
        try:
            return bool(tip.size and not np.all(tip == 0.0))
        except Exception:
            return False
        
    @property
    def spectra_sect(self):
        return self.get_section("spectra_sect")
    
    @property
    def timeseries_sect(self):
        return self.get_section("timeseries_sect")

    @property
    def spectra(self):
        """Return high-level Spectra object if present."""
        return self.get_section("spectra")

    @property
    def spectra_io(self):
        """Return SpectraIO if present."""
        return self.get_section("spectra_io")

    @property
    def timeseries(self):
        """Return high-level TimeSeries object if present."""
        return self.get_section("timeseries")

    @property
    def timeseries_io(self):
        """Return TSIO if present."""
        return self.get_section("timeseries_io")

    @property
    def channels(self) -> list[str]:
        """Return TS channels if any, else []."""
        ts = self.get_section("timeseries")
        if ts is None:
            return []
        try:
            return list(ts.channels())
        except Exception:
            return []

    @property
    def path_str(self) -> str:
        """EDI path as string or empty if not set."""
        return str(self.path) if self.path is not None else ""

    @property
    def edi_dir(self) -> Path | None:
        """Parent directory of EDI file if path set."""
        return self.path.parent if self.path is not None else None

    @property
    def processingsoftware(self) -> str | None:
        """Return INFO.Processing.ProcessingSoftware.name."""
        info = self.get_section("info")
        if info is None:
            return None
        ps = getattr(info, "Processing", None)
        if ps is None:
            return None
        sw = getattr(ps, "ProcessingSoftware", None)
        return getattr(sw, "name", None)

    @processingsoftware.setter
    def processingsoftware(self, name: str) -> None:
        """Set INFO.Processing.ProcessingSoftware.name if any."""
        info = self.get_section("info")
        if info is None:
            return
        ps = getattr(info, "Processing", None)
        if ps is None:
            return
        sw = getattr(ps, "ProcessingSoftware", None)
        if sw is None:
            return
        try:
            sw.name = str(name)
        except Exception:
            pass
        
    def __repr__(self) -> str:  # pragma: no cover
        # concise, eval-like summary for debugging.
        p = str(self.path) if self.path else "-"
        nf = self.n_freq # int(getattr(self.Z, "n_freq", 0) or 0)
        has_mt = bool(self.has_section("mtsect"))
        has_sp = bool(self.get_section("spectra") is not None)
        has_ts = bool(self.get_section("timeseries") is not None)
        tip_ok = False
        tip = getattr(self.Tip, "tipper", None)
        if tip is not None:
            tip_ok = bool(tip.size and not np.all(tip == 0.0))
        return (
            "EDIFile("
            f"path={p!r}, n_freq={nf}, mt={has_mt}, "
            f"spectra={has_sp}, ts={has_ts}, tipper={tip_ok}"
            ")"
        )
    
    
    def __str__(self) -> str:  # pragma: no cover
        # human-friendly multi-line snapshot.
        p = str(self.path) if self.path else "-"
        head = self.get_section("head")
        sid = getattr(head, "dataid", None) if head else None
        sid = sid or "?"
        nf = self.n_freq #int(getattr(self.Z, "n_freq", 0) or 0)
    
        tags: list[str] = []
        if self.has_section("mtsect"):
            tags.append("MT/EMAP")
        if self.get_section("spectra") is not None:
            tags.append("SPECTRA")
        if self.get_section("timeseries") is not None:
            tags.append("TSERIES")
        if self.get_section("other") is not None:
            tags.append("OTHER")
    
        tip_ok = False
        tip = getattr(self.Tip, "tipper", None)
        if tip is not None:
            tip_ok = bool(tip.size and not np.all(tip == 0.0))
    
        lines = [
            "EDIFile",
            f"  path   : {p}",
            f"  station: {sid}",
            f"  n_freq : {nf}",
            f"  tipper : {'yes' if tip_ok else 'no'}",
            f"  sects  : {', '.join(tags) if tags else '-'}",
        ]
        return "\n".join(lines)
    

def _as_path(p: str | Path | None) -> Optional[Path]:
    if p is None:
        return None
    return Path(p)


def _nz(a: np.ndarray, n: int) -> np.ndarray:
    if a.size == n:
        return a
    out = np.zeros(n, dtype=a.dtype)
    out[: min(n, a.size)] = a[: min(n, a.size)]
    return out


def _lower_keys(d: Dict[str, Iterable[float]]) -> Dict[str, list]:
    out: Dict[str, list] = {}
    for k, v in d.items():
        out[k.lower()] = list(v)
    return out


def _rev_if_asc(freq: np.ndarray) -> Tuple[np.ndarray, bool]:
    if freq.size == 0:
        return freq, False
    rev = bool(freq[-1] > freq[0])
    return (freq[::-1], True) if rev else (freq, False)


