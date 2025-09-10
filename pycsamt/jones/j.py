# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

 
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union
import math
import numpy as np

from ..constants import MU_0
from ..log.logger import get_logger
from ..z.z import Z
from ..z.tipper import Tipper
from ..z.resphase import ResPhase

from .config import (
    FLOAT_FORMAT_R,
    FLOAT_FORMAT_TF,
    TENSOR_INDEX
)
from .heads import Heads
from .blocks import JBlocks, RBlock, TFBlock

logger = get_logger(__name__)

def _as_path(p: Union[str, Path]) -> Path:
    return p if isinstance(p, Path) else Path(p)

__all__= ['JMixin', "JIOMixin", "JFile"]

class JMixin:
    r"""
    Lightweight helpers shared by J-format readers.

    The mixin groups small, allocation-friendly utilities
    that many parsers need.  It keeps math, alignment and
    token-handling logic out of higher-level classes.

    Notes
    -----
    The helpers are intentionally tiny and avoid importing
    heavy libraries beyond :mod:`numpy`.  They favor pure
    functions that accept and return arrays.

    Methods
    -------
    _complex(re, im)
        Combine real and imaginary parts into a complex
        vector.  Shapes must be broadcast-compatible.

    _deg2rad(x)
        Degrees to radians conversion with float return.

    _hz_from_period(p)
        Convert period seconds to frequency in Hz.  Non-
        positive entries are mapped to ``nan``.

    _align_by_periods(p0, p1)
        Return a tuple ``(p_common, i0, i1)`` where
        ``p_common`` are periods present in both sequences,
        ordered like ``p0``; ``i0`` and ``i1`` are indices
        to select the matching rows in the original arrays.

    Examples
    --------
    >>> jm = JMixin()
    >>> jm._complex([1, -2], [0.5, 3]).dtype.kind == 'c'
    True
    >>> jm._deg2rad(180.0)
    3.141592653589793
    >>> jm._hz_from_period([1.0, 0.5])
    array([1., 2.])
    >>> pc, i0, i1 = jm._align_by_periods([1, 2, 3], [3, 1])
    >>> pc.tolist(), i0.tolist(), i1.tolist()
    ([1, 3], [0, 2], [1, 0])

    See Also
    --------
    JIOMixin : Block scanning and object building helpers.
    JFile : High-level reader that uses both mixins.

    References
    ----------
    .. [1] A. G. Jones, *Magnetotelluric data file J-format*,
           version 2.0, 1994.
    .. [2] MTNet, *J format documentation*.
    """

    _tidx: Dict[str, Tuple[int, int]] = dict(TENSOR_INDEX)

    @staticmethod
    def _deg2rad(a: np.ndarray) -> np.ndarray:
        return np.deg2rad(a, dtype=float)

    @staticmethod
    def _hz_from_period(p: np.ndarray) -> np.ndarray:
        return np.where(p > 0.0, 1.0 / p, np.nan)

    @staticmethod
    def _complex(re: np.ndarray, im: np.ndarray) -> np.ndarray:
        return re.astype(float) + 1j * im.astype(float)

    @staticmethod
    def _align_by_periods(
        a_p: np.ndarray,
        b_p: np.ndarray,
        *,
        rtol: float = 1e-8,
        atol: float = 1e-12,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        if a_p.size == 0 or b_p.size == 0:
            z = np.array([])
            zi = np.array([], dtype=int)
            zj = np.array([], dtype=int)
            return z, zi, zj

        ap = a_p.copy()
        bp = b_p.copy()
        idx_a: List[int] = []
        idx_b: List[int] = []
        for i, pa in enumerate(ap):
            j = int(np.argmin(np.abs(bp - pa)))
            if math.isclose(pa, bp[j], rel_tol=rtol, abs_tol=atol):
                idx_a.append(i)
                idx_b.append(j)
        if not idx_a:
            z = np.array([])
            zi = np.array([], dtype=int)
            zj = np.array([], dtype=int)
            return z, zi, zj

        ia = np.asarray(idx_a, dtype=int)
        ib = np.asarray(idx_b, dtype=int)
        p = ap[ia]
        order = np.argsort(p)
        return p[order], ia[order], ib[order]


class JIOMixin(JMixin):
    r"""
    Tolerant BLOCK parser and TF/R/Tipper builder.

    This mixin concentrates the I/O-oriented parts for J
    files: scanning blocks, aligning periods, normalizing
    rows, and assembling higher-level objects (``Z``,
    ``Tipper``, and ``ResPhase``).

    It aims to be robust against non-canonical header
    orders and minor format quirks.

    Notes
    -----
    
    * Data rows are normalized before assembly.  Period
      sign conventions (negative means Hz) are corrected.
    * Missing sentinels (e.g., ``-999``) are mapped to
      ``nan`` for numeric arrays, but objects are pre-
      allocated with zeros to keep shapes consistent.
    * When only rho/phi are present, impedance is rebuilt
      using :math:`|Z|=\sqrt{\mu_0 \,\omega\, \rho}` and
      :math:`\phi` for the phase.  The vacuum permeability
      :math:`\mu_0` is imported from
      :mod:`pycsamt.constants`.

    Methods
    -------
    _scan_blocks(path, *, start=None, empty_val=...)
        Parse the file into a component dictionary indexed
        by tokens like ``'ZXY'``, ``'RXX'`` or ``'TZX'``.
        Values include period, real/imag/error (TF) or rho/
        phi (+ auxiliary columns for R blocks).

    _build_from_comp(comp, *, z_obj=None, tip_obj=None)
        Assemble :class:`~pycsamt.z.z.Z`,
        :class:`~pycsamt.z.tipper.Tipper`, and
        :class:`~pycsamt.z.resphase.ResPhase` from the
        scanned components.  Returns a triple
        ``(Z|None, Tipper|None, ResPhase|None)``.

    Examples
    --------
    >>> mix = JIOMixin()
    >>> # (typically used via JFile; direct use shown here)
    >>> # comp = mix._scan_blocks(Path("data/j/site.j"))
    >>> # z, tip, rp = mix._build_from_comp(comp)

    See Also
    --------
    JMixin : Mathematical utilities used by this mixin.
    JFile : High-level reader/writer that calls these
        methods.
    pycsamt.z.z.Z : Impedance tensor container.
    pycsamt.z.tipper.Tipper : Tipper container.
    pycsamt.z.resphase.ResPhase : R–φ container.

    References
    ----------
    .. [1] A. G. Jones, *Magnetotelluric data file J-format*,
           version 2.0, 1994.
    .. [2] MTNet, *J format documentation*.
    """
    
    def _scan_blocks(
        self,
        path: Path,
        *,
        start: Optional[int] = None,  # noqa: ARG002
        empty_val: float = -999.0,  # noqa: ARG002
    ) -> Dict[str, Dict[str, Any]]:
        jb = JBlocks.from_file(path, verbose=getattr(self, "verbose", 0))
        comp: Dict[str, Dict[str, Any]] = {}

        for blk in jb.blocks:
            if isinstance(blk, RBlock):
                a = blk.to_numpy()
                comp[blk.head.dtype.kind + blk.head.dtype.comp] = {
                    "period": a["period"],
                    "rho": a["rho"],
                    "pha": a["pha"],
                    "rhomax": a["rhomax"],
                    "rhomin": a["rhomin"],
                    "phamax": a["phamax"],
                    "phamin": a["phamin"],
                    "wrho": a["wrho"],
                    "wpha": a["wpha"],
                    "rej": a["rej"].astype(bool),
                    "head": blk.head,
                    "kind": blk.head.dtype.kind,
                }
            elif isinstance(blk, TFBlock):
                a = blk.to_numpy()
                comp[blk.head.dtype.kind + blk.head.dtype.comp] = {
                    "period": a["period"],
                    "real": a["real"],
                    "imag": a["imag"],
                    "err": a["error"],
                    "w": a["weight"],
                    "rej": a["rej"].astype(bool),
                    "head": blk.head,
                    "kind": blk.head.dtype.kind,
                }
            else:
                logger.debug("Skip block: %s", type(blk).__name__)
        return comp

    def _build_from_comp(
        self,
        comp: Dict[str, Dict[str, Any]],
        *,
        z_obj: Optional[Z] = None,
        tip_obj: Optional[Tipper] = None,
    ) -> Tuple[Optional[Z], Optional[Tipper], Optional[ResPhase]]:
        """
        Assemble Z, Tipper, and Res/Phase from parsed component dicts.
        Uses zero initialization (no NaNs) for API consistency.
        """
    
        z_parts: Dict[str, Dict[str, Any]] = {
            k: d for k, d in comp.items() if k[:1] in {"Z", "Q"}
        }
        r_parts: Dict[str, Dict[str, Any]] = {
            k: d for k, d in comp.items() if k[:1] in {"R", "S"}
        }
        t_parts: Dict[str, Dict[str, Any]] = {
            k: d for k, d in comp.items() if k[:1] == "T"
        }
    
        z_out: Optional[Z] = z_obj
        rp_out: Optional[ResPhase] = None
    
        # Z from TF (preferred) 
        if z_out is None and z_parts:
            k0 = next(iter(z_parts.keys()))
            p0 = np.asarray(z_parts[k0]["period"], dtype=float)
            freq = self._hz_from_period(p0)
            n = int(p0.size)
    
            zxx = np.zeros(n, dtype=complex)
            zxy = np.zeros(n, dtype=complex)
            zyx = np.zeros(n, dtype=complex)
            zyy = np.zeros(n, dtype=complex)
    
            exx = np.zeros(n, dtype=float)
            exy = np.zeros(n, dtype=float)
            eyx = np.zeros(n, dtype=float)
            eyy = np.zeros(n, dtype=float)
    
            for key, d in z_parts.items():
                comp_code = key[-2:].upper()
                p = np.asarray(d["period"], dtype=float)
                _, ia, ib = self._align_by_periods(p0, p)
                if ia.size == 0:
                    continue
                zval = self._complex(d["real"][ib], d["imag"][ib])
                err = np.asarray(d.get("err", 0.0), float)
                err = err[ib] if err.size else np.zeros(ib.size, dtype=float)
    
                if comp_code == "XX":
                    zxx[ia] = zval
                    exx[ia] = err
                elif comp_code == "XY":
                    zxy[ia] = zval
                    exy[ia] = err
                elif comp_code == "YX":
                    zyx[ia] = zval
                    eyx[ia] = err
                elif comp_code == "YY":
                    zyy[ia] = zval
                    eyy[ia] = err
    
            z_mat = np.zeros((n, 2, 2), dtype=complex)
            z_err = np.zeros((n, 2, 2), dtype=float)
    
            z_mat[:, 0, 0] = zxx
            z_mat[:, 0, 1] = zxy
            z_mat[:, 1, 0] = zyx
            z_mat[:, 1, 1] = zyy
    
            z_err[:, 0, 0] = exx
            z_err[:, 0, 1] = exy
            z_err[:, 1, 0] = eyx
            z_err[:, 1, 1] = eyy
    
            z_out = Z(
                z_array=z_mat,
                z_err_array=z_err,
                freq=freq,
                name=z_parts[k0]["head"].station,  # type: ignore
            )
    
            # Also provide Res/Phase when Z exists.
            rp_out = ResPhase(
                freq=freq, name=z_parts[k0]["head"].station
                )  # type: ignore  # noqa: E501
            try:
                rp_out.compute_resistivity_phase(
                    z_array=z_mat,
                    z_err_array=z_err,
                    freq=freq,
                )
            except Exception:
                pass
    
        # Z from R/phi (fallback) 
        if z_out is None and r_parts:
            have_keys = [k for k in ("RXX", "RXY", "RYX", "RYY") if k in r_parts]
            if not have_keys:
                have_keys = list(r_parts.keys())
            k0 = have_keys[0]
    
            p0 = np.asarray(r_parts[k0]["period"], dtype=float)
            freq = self._hz_from_period(p0)
            n = int(p0.size)
            w = 2.0 * math.pi * freq
    
            def to_imp(d: Dict[str, Any]) -> np.ndarray:
                rho = np.asarray(d["rho"], dtype=float)
                pha = np.asarray(d["pha"], dtype=float)
                mag = np.sqrt(np.maximum(rho, 0.0) * MU_0 * w)
                ang = self._deg2rad(pha)
                return mag * (np.cos(ang) + 1j * np.sin(ang))
    
            zxx = np.zeros(n, dtype=complex)
            zxy = np.zeros(n, dtype=complex)
            zyx = np.zeros(n, dtype=complex)
            zyy = np.zeros(n, dtype=complex)
    
            for key, d in r_parts.items():
                comp_code = key[-2:].upper()
                p = np.asarray(d["period"], float)
                _, ia, ib = self._align_by_periods(p0, p)
                if ia.size == 0:
                    continue
                v = to_imp(d)[ib]
                if comp_code == "XX":
                    zxx[ia] = v
                elif comp_code == "XY":
                    zxy[ia] = v
                elif comp_code == "YX":
                    zyx[ia] = v
                elif comp_code == "YY":
                    zyy[ia] = v
    
            z_mat = np.zeros((n, 2, 2), dtype=complex)
            z_mat[:, 0, 0] = zxx
            z_mat[:, 0, 1] = zxy
            z_mat[:, 1, 0] = zyx
            z_mat[:, 1, 1] = zyy
    
            z_out = Z(z_array=z_mat, freq=freq, name=k0)
    
            rp_out = ResPhase(freq=freq, name=k0)  # type: ignore
            try:
                rp_out.compute_resistivity_phase(
                    z_array=z_mat,
                    freq=freq,
                )
            except Exception:
                pass
    
        # Tipper if present 
        tip_out: Optional[Tipper] = tip_obj
        if tip_out is None and t_parts:
            k0 = next(iter(t_parts.keys()))
            p0 = np.asarray(t_parts[k0]["period"], dtype=float)
            freq = self._hz_from_period(p0)
            n = int(p0.size)
    
            tzx = np.zeros(n, dtype=complex)
            tzy = np.zeros(n, dtype=complex)
            ezx = np.zeros(n, dtype=float)
            ezy = np.zeros(n, dtype=float)
    
            for key, d in t_parts.items():
                comp_code = key[-2:].upper()
                p = np.asarray(d["period"], float)
                _, ia, ib = self._align_by_periods(p0, p)
                if ia.size == 0:
                    continue
                v = self._complex(d["real"][ib], d["imag"][ib])
                e = np.asarray(d.get("err", 0.0), dtype=float)
                e = e[ib] if e.size else np.zeros(ib.size, dtype=float)
                if comp_code == "ZX":
                    tzx[ia] = v
                    ezx[ia] = e
                elif comp_code == "ZY":
                    tzy[ia] = v
                    ezy[ia] = e
    
            t_arr = np.column_stack((tzx, tzy))
            e_arr = np.column_stack((ezx, ezy))
            tip_out = Tipper(
                tipper_array=t_arr,
                tipper_err_array=e_arr,
                freq=freq,
            )
    
        return z_out, tip_out, rp_out


class JFile(JIOMixin):
    r"""
    High-level J dispatcher for MT/SEG archives.

    The class reads a J file, extracts headers and blocks,
    and builds analysis-ready objects for impedance (``Z``),
    resistivity/phase (``Res``), and tipper (``Tip``).  It
    also writes new J files from the in-memory state.

    Parameters
    ----------
    path : str or Path, optional
        Input J file.  If omitted, call :meth:`read` later.
    verbose : int, default=0
        Verbosity level.  Non-zero emits informational
        messages during parsing and writing.

    Attributes
    ----------
    path : Path or None
        Source path when set via ``__init__`` or
        :meth:`from_file`.
    heads : pycsamt.jones.heads.Heads or None
        Parsed banner, info and a single head triple.
    blocks : pycsamt.jones.blocks.JBlocks or None
        Parsed data blocks (R and/or TF).
    Z : pycsamt.z.z.Z or None
        Impedance tensor container, possibly rebuilt from
        rho/phi when TF are absent.
    Tip : pycsamt.z.tipper.Tipper or None
        Tipper container if ZX/ZY are present.
    Res : pycsamt.z.resphase.ResPhase or None
        Resistivity/phase view (direct or derived).
    freq : ndarray or None
        Shared frequency vector inferred from available
        objects.  Periods are available via
        :pyattr:`periods`.
    periods : ndarray or None
        Convenience view of ``1.0/freq`` when known.
    n_freq : int
        Number of frequency samples (``0`` if unknown).
    name : str or None
        Friendly site/station name.  Precedence is:
        ``Z.name`` → head.station → file stem.
    site : str or None
        Alias for the station code (if known).
    lat, lon, azimuth, az_hint, elev : float or None
        Geographic metadata proxied from headers.

    Methods
    -------
    from_file(path, *, verbose=0)
        Construct and read in one call.
    read(path=None, *, start=None)
        Parse headers and blocks, then build ``Z``,
        ``Tip`` and ``Res`` as available.
    write(j_fn=None, new_jfn=None, datatype=None,
          savepath=None, *, verbose=None, overwrite=True)
        Serialize the current state to a J file.  The
        ``datatype`` selector accepts combinations like
        ``'Z'``, ``'R'``, ``'T'``, ``'ZR'``, or ``'ALL'``.
    compose_headers()
        Render banner, head and info lines only.
    __has_read__()
        Return ``True`` once a full :meth:`read` completed.

    Examples
    --------
    >>> jf = JFile.from_file("data/j/kb0-s001.txt", verbose=0)
    >>> jf.n_freq > 0
    True
    >>> out = jf.write(new_jfn="out.j", datatype="ZR",
    ...                overwrite=True)
    >>> isinstance(out, str)
    True
    >>> jf.lat, jf.lon  # site coordinates if present
    ( ... )  # doctest: +SKIP

    Notes
    -----
    * ``write`` prefers existing uncertainties; missing
      errors are filled with zeros.  Periods are written
      from the active frequency vector.
    * When only R-blocks exist, ``Z`` is rebuilt so that
      downstream code can still compute QA metrics or
      plot tensor-based products.

    See Also
    --------
    pycsamt.jones.heads.Heads : Header and metadata view.
    pycsamt.jones.blocks.JBlocks : Low-level parsed blocks.
    pycsamt.z.z.Z : Impedance tensor class.
    pycsamt.z.tipper.Tipper : Tipper class.
    pycsamt.z.resphase.ResPhase : R–φ class.

    References
    ----------
    .. [1] A. G. Jones, *Magnetotelluric data file J-format*,
           version 2.0, 1994.
    .. [2] MTNet, *J format documentation*.
    """

    def __init__(
        self,
        path: str | Path | None = None,
        *,
        verbose: int = 0,
    ) -> None:
        self.verbose = verbose
        self.path: Optional[Path] = _as_path(path) if path else None
        self.heads: Optional[Heads] = None
        self.blocks: Optional[JBlocks] = None
        self.Z: Optional[Z] = None
        self.Tip: Optional[Tipper] = None
        self.Res: Optional[ResPhase] = None
        self._read_ok: bool = False
        
        if self.path is not None: 
            self.read(self.path)

    @classmethod
    def from_file(
        cls, path: str | Path, *, verbose: int = 0
    ) -> "JFile":
        r"""
        Construct and read a J file in one call.

        This convenience constructor mirrors ``__init__`` +
        :meth:`read`.  It resolves ``path`` to a filesystem
        location, parses headers and data blocks, and builds
        analysis-ready objects (``Z``, ``Tip``, ``Res``) when
        present or derivable.

        Parameters
        ----------
        path : str or Path
            Path to a J-format text file (Jones v2.0 style).
        verbose : int, default=0
            Verbosity flag.  When non-zero, progress/info
            messages may be emitted during parsing.

        Returns
        -------
        JFile
            Instance with :pyattr:`heads`, :pyattr:`blocks`
            and objects (:pyattr:`Z`, :pyattr:`Tip`,
            :pyattr:`Res`) populated where possible.

        Notes
        -----
        * Headers (banner + ``>KEY=VALUE`` + the first head
          triple) are parsed via :class:`~pycsamt.jones.heads.Heads`.
        * Blocks are scanned with
          :class:`~pycsamt.jones.blocks.JBlocks`, then
          assembled into ``Z``/``Tip``/``Res`` via
          :class:`JIOMixin`.

        Examples
        --------
        >>> jf = JFile.from_file("data/j/kb0-s001.txt")
        >>> jf.n_freq > 0
        True

        See Also
        --------
        JFile.read : Lower-level method if you already have
            an instance.
        JBlocks : Low-level block parser.
        Heads : Header and site metadata container.

        References
        ----------
        .. [1] A. G. Jones, *Magnetotelluric data file
               J-format*, version 2.0, 1994.
        """
        inst = cls(path, verbose=verbose)
        # inst.read(path)
        return inst

    def read(
        self,
        path: str | Path | None = None,
        *,
        start: Optional[int] = None,
    ) -> "JFile":
        r"""
        Parse a J file and build objects in memory.

        The method reads the banner, info block and the first
        head triple, then scans all following data blocks.
        Transfer functions (``Zxx``, ``Zxy``, ...) and tipper
        (``Tzx``, ``Tzy``) are assembled when present.  If only
        resistivity/phase blocks exist, a synthetic impedance
        is rebuilt from :math:`\rho` and :math:`\phi`.

        Parameters
        ----------
        path : str or Path, optional
            If given, set as the source and read from it.
            If omitted, reuse ``self.path`` set at construction.
        start : int, optional
            Line index hint to start scanning blocks.  Most
            users can leave this as ``None``.

        Returns
        -------
        JFile
            The instance itself (for chaining).

        Notes
        -----
        * Block scanning is tolerant to minor format quirks
          (blank lines, non-canonical head order where the row
          count precedes the data-type).
        * Period sign conventions are normalized (negative
          values mean input was frequency in Hz).
        * Missing sentinels (e.g., ``-999``) are mapped to
          ``nan`` in numeric arrays, while objects are pre-
          allocated with zeros to keep shapes consistent.

        Examples
        --------
        >>> jf = JFile(verbose=0)
        >>> _ = jf.read("data/j/kb0-s001.txt")
        >>> jf.Z is not None or jf.Res is not None
        True

        See Also
        --------
        JFile.from_file : Shortcut that constructs then calls
            this method.
        JBlocks : Underlying block parser.
        Z, Tipper, ResPhase : Target containers built by this
            method.

        References
        ----------
        .. [1] A. G. Jones, *Magnetotelluric data file
               J-format*, version 2.0, 1994.
        """

        if path is not None:
            self.path = _as_path(path)
        if self.path is None:
            raise ValueError("path is required")

        self.heads = Heads.from_file(
            self.path, verbose=self.verbose
        )
        self.blocks = JBlocks.from_file(
            self.path, verbose=self.verbose)

        comp = self._scan_blocks(self.path, start=start)
        z, tip, rp = self._build_from_comp(comp)
        self.Z = z
        self.Tip = tip
        self.Res = rp

        self._read_ok = True
        return self

    def write(  # noqa: C901
        self,
        j_fn: str | None = None,
        new_jfn: str | None = None,
        datatype: str | None = None,
        savepath: str | Path | None = None,
        *,
        verbose: int | None = None,
        overwrite: bool = True,
    ) -> str:
        r"""
        Serialize the current state to a J-format file.

        The writer renders a banner, info lines and one or more
        data blocks selected via ``datatype``.  When uncertainties
        are unavailable, zero-filled error columns are emitted to
        preserve column layout.  Periods are derived from the
        active frequency vector.

        Parameters
        ----------
        j_fn : str, optional
            Base filename to use.  If omitted, derive from
            :pyattr:`path` or default to ``'out.j'``.
        new_jfn : str, optional
            Replacement filename.  Takes precedence over
            ``j_fn`` when provided.
        datatype : {'Z','R','T','ZR','RT','ZT','ZRT','ALL'}, optional
            Select families to emit.  If ``None``, the writer
            auto-detects from available objects on the instance.
        savepath : str or Path, optional
            Folder where to save.  Defaults to the parent of
            :pyattr:`path` or the current directory.
        verbose : int, optional
            Override verbosity.  If ``None``, reuse
            :pyattr:`verbose`.
        overwrite : bool, default=True
            If ``False`` and the target exists, a numeric
            suffix is appended to avoid clobbering.

        Returns
        -------
        out_path : str
            The filesystem path of the written file.

        Notes
        -----
        * The station code and optional azimuth hint are taken
          from the parsed head.  Units for transfer functions
          are written as ``SI``.
        * If only ``R`` blocks exist, they can be emitted
          directly; if only ``Z`` is present, synthetic
          ``R``/``φ`` can be computed for writing when the
          selector requests it.
        * The banner defaults to ``PYCSAMT`` when no producer
          is known.  The original banner (if parsed) can be
          preserved or referenced by the caller before writing.

        Examples
        --------
        >>> jf = JFile.from_file("data/j/kb0-s001.txt")
        >>> out = jf.write(new_jfn="site_out.j",
        ...                datatype="ZR", overwrite=True)
        >>> isinstance(out, str)
        True

        See Also
        --------
        JFile.compose_headers : Render banner + headers only.
        Z, Tipper, ResPhase : Sources used by the writer.

        References
        ----------
        .. [1] A. G. Jones, *Magnetotelluric data file
               J-format*, version 2.0, 1994.
        """

        vb = self.verbose if verbose is None else int(verbose)

        def _ensure_parent(p: Path) -> None:
            p.parent.mkdir(parents=True, exist_ok=True)

        def _fmt(fmt: str, v: float) -> str:
            try:
                return fmt.format(v)
            except Exception:
                return fmt.format(float("nan"))

        def _fmt_row_r(p: float, rho: float, pha: float) -> str:
            return " ".join(
                [
                    _fmt(FLOAT_FORMAT_R, p),
                    _fmt(FLOAT_FORMAT_R, rho),
                    _fmt(FLOAT_FORMAT_R, pha),
                    _fmt(FLOAT_FORMAT_R, rho),
                    _fmt(FLOAT_FORMAT_R, rho),
                    _fmt(FLOAT_FORMAT_R, pha),
                    _fmt(FLOAT_FORMAT_R, pha),
                    _fmt(FLOAT_FORMAT_R, 1.0),
                    _fmt(FLOAT_FORMAT_R, 1.0),
                ]
            )

        def _fmt_row_tf(
            p: float,
            re: float,
            im: float,
            err: float,
            w: float = 1.0,
        ) -> str:
            return " ".join(
                [
                    _fmt(FLOAT_FORMAT_TF, p),
                    _fmt(FLOAT_FORMAT_TF, re),
                    _fmt(FLOAT_FORMAT_TF, im),
                    _fmt(FLOAT_FORMAT_TF, err),
                    _fmt(FLOAT_FORMAT_TF, w),
                ]
            )

        def _periods_from_freq(f: Optional[np.ndarray]) -> np.ndarray:
            if f is None:
                return np.array([], float)
            with np.errstate(divide="ignore", invalid="ignore"):
                p = np.where(f > 0.0, 1.0 / f, np.nan)
            return p.astype(float)

        # decide families to emit 
        have_z = (
                self.Z is not None
                and (
                    getattr(self.Z, "_z", None) is not None
                    or getattr(self.Z, "z", None) is not None
                )
            )

        have_r = self.Res is not None
        have_t = self.Tip is not None

        if datatype is None:
            flags = {"Z": have_z, "R": have_r, "T": have_t}
        else:
            key = datatype.upper().strip()
            if key in {"ALL", "ZRT"}:
                flags = {"Z": have_z, "R": have_r, "T": have_t}
            else:
                flags = {
                    "Z": "Z" in key and have_z,
                    "R": "R" in key and have_r,
                    "T": "T" in key and have_t,
                }

        # header preamble 
        lines: List[str] = []

        # banner
        if self.heads is not None and hasattr(self.heads, "banner"):
            try:
                lines.extend(
                    self.heads.banner.write(include_origin=True))
            except Exception:
                pass

        # info
        if self.heads is not None and hasattr(self.heads, "info"):
            try:
                lines.extend(self.heads.info.write())
            except Exception:
                pass

        # station id and optional azimuth hint
        station = "SITE"
        az_hint: str = ""
        if self.heads is not None and self.heads.head is not None:
            if self.heads.head.station:
                station = str(self.heads.head.station).upper()
            ah = getattr(self.heads.head, "az_hint", None)
            if isinstance(ah, (int, float)):
                az_hint = f" {ah:g}"

        #  helper: write 1 block 
        def _write_head(dtype: str, nrows: int) -> None:
            lines.append(f"{station}{az_hint}")
            lines.append(dtype)
            lines.append(str(int(nrows)))

        #  write Z 
        if flags.get("Z", False):
            z = self.Z  # already checked not None
            za = getattr(z, "_z", None)
            if za is None:
                za = getattr(z, "z", None)
            
            ze = getattr(z, "_z_err", None)
            if ze is None: 
               ze = getattr(z, "z_err", None)
        
            f = getattr(z, "freq", None)
            p = _periods_from_freq(f)

            # components map
            comps = {
                "XX": (0, 0),
                "XY": (0, 1),
                "YX": (1, 0),
                "YY": (1, 1),
            }
            for code, (i, j) in comps.items():
                if za is None:
                    continue
                try:
                    v = np.asarray(za)[:, i, j]
                except Exception:
                    continue
                
                # errors
                err = None
                if ze is None:
                    err = np.zeros_like(v.real, dtype=float)
                    
                if ze is not None:
                    try:
                        err = np.asarray(ze)[:, i, j].astype(float)
                    except Exception:
                        err = None
    
                if err is None:
                    err = np.zeros_like(v.real, dtype=float)

                # valid rows: period & value finite
                keep = np.isfinite(p) & np.isfinite(v.real) & np.isfinite(
                    v.imag
                )
                nrows = int(np.count_nonzero(keep))
                if nrows == 0:
                    continue

                _write_head(f"Z{code} SI", nrows)
                for pi, zi, ei in zip(p[keep], v[keep], err[keep]):
                    lines.append(
                        _fmt_row_tf(float(pi), float(zi.real),
                                    float(zi.imag), float(ei), 1.0)
                    )

        # write R 
        if flags.get("R", False):
            # prefer direct Res; else compute from Z
            rp = self.Res  # type: ignore[assignment]
            if rp is not None:
                f = getattr(rp, "freq", None)
                p = _periods_from_freq(f)
                # expect attributes like res_xy, phase_xy, ...
                pairs = {
                    "XX": ("res_xx", "phase_xx"),
                    "XY": ("res_xy", "phase_xy"),
                    "YX": ("res_yx", "phase_yx"),
                    "YY": ("res_yy", "phase_yy"),
                }
                for code, (rname, pname) in pairs.items():
                    rho = getattr(rp, rname, None)
                    pha = getattr(rp, pname, None)
                    if rho is None or pha is None:
                        continue
                    r = np.asarray(rho, float)
                    ph = np.asarray(pha, float)
                    keep = np.isfinite(p) & np.isfinite(r) & np.isfinite(ph)
                    nrows = int(np.count_nonzero(keep))
                    if nrows == 0:
                        continue
                    _write_head(f"R{code}", nrows)
                    for pi, ri, phii in zip(p[keep], r[keep], ph[keep]):
                        lines.append(_fmt_row_r(float(pi), float(ri),
                                                float(phii)))
            elif self.Z is not None:
                # derive from Z magnitudes
                z = self.Z
                za = getattr(z, "z")
                f = getattr(z, "freq", None)
                p = _periods_from_freq(f)
                if za is not None:
                    mu0 = 4.0e-7 * math.pi
                    w = 2.0 * math.pi * np.asarray(f, float)
                    comps = {
                        "XX": (0, 0),
                        "XY": (0, 1),
                        "YX": (1, 0),
                        "YY": (1, 1),
                    }
                    for code, (i, j) in comps.items():
                        try:
                            zc = np.asarray(za)[:, i, j]
                        except Exception:
                            continue
                        mag2 = (zc.real ** 2 + zc.imag ** 2)
                        with np.errstate(divide="ignore", invalid="ignore"):
                            rho = mag2 / (mu0 * w)
                        pha = np.degrees(np.arctan2(zc.imag, zc.real))
                        keep = (
                            np.isfinite(p)
                            & np.isfinite(rho)
                            & np.isfinite(pha)
                        )
                        nrows = int(np.count_nonzero(keep))
                        if nrows == 0:
                            continue
                        _write_head(f"R{code}", nrows)
                        for pi, ri, phii in zip(
                            p[keep], rho[keep], pha[keep]
                        ):
                            lines.append(
                                _fmt_row_r(float(pi), float(ri),
                                           float(phii))
                            )

        # write Tipper 
        if flags.get("T", False) and self.Tip is not None:
            tp = self.Tip
            ta = getattr(tp, "tipper")
            te = getattr(tp, "tipper_err", None)
            f = getattr(tp, "freq", None)
            p = _periods_from_freq(f)
            
            arr = None
            # shape: (n, 2) -> ZX, ZY
            if ta is not None:
                arr = np.asarray(ta)
                if arr.ndim == 3 and arr.shape[1:] == (1, 2):
                    arr = arr[:, 0, :]
                if arr.shape[-1] != 2:
                    arr = None
            if arr is not None:
                err = None
                if te is not None:
                    err = np.asarray(te, float)
                    if err.ndim == 3 and err.shape[1:] == (1, 2):
                        err = err[:, 0, :]
                    if err.shape[-1] != 2:
                        err = None
                if err is None:
                    err = np.zeros_like(arr.real, float)

                comps = {"ZX": 0, "ZY": 1}
                for code, k in comps.items():
                    v = arr[:, k]
                    e = err[:, k]
                    keep = (
                        np.isfinite(p) & np.isfinite(v.real)
                        & np.isfinite(v.imag)
                    )
                    nrows = int(np.count_nonzero(keep))
                    if nrows == 0:
                        continue
                    _write_head(f"T{code}", nrows)
                    for pi, zi, ei in zip(p[keep], v[keep], e[keep]):
                        lines.append(
                            _fmt_row_tf(float(pi), float(zi.real),
                                        float(zi.imag), float(ei), 1.0)
                        )

        # finalize target 
        base = None
        if new_jfn:
            base = Path(new_jfn)
        elif j_fn:
            base = Path(j_fn)
        elif self.path is not None:
            base = Path(self.path.name)
        else:
            base = Path("out.j")

        folder = None
        if savepath is not None:
            folder = Path(savepath)
        elif self.path is not None and self.path.parent.exists():
            folder = self.path.parent
        else:
            folder = Path(".")

        out_path = folder / base
        if not overwrite and out_path.exists():
            stem = out_path.stem
            suf = out_path.suffix or ".j"
            k = 1
            while True:
                cand = out_path.with_name(f"{stem}_{k}{suf}")
                if not cand.exists():
                    out_path = cand
                    break
                k += 1

        _ensure_parent(out_path)
        out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

        if vb:
            logger.info("Wrote J file: %s", out_path)

        return str(out_path)


    @property
    def freq(self) -> Optional[np.ndarray]:
        if self.Z is not None and hasattr(self.Z, "freq"):
            return self.Z.freq  # type: ignore[return-value]
        if self.Tip is not None and hasattr(self.Tip, "freq"):
            return self.Tip.freq  # type: ignore[return-value]
        if self.Res is not None and hasattr(self.Res, "freq"):
            return self.Res.freq  # type: ignore[return-value]
        return None

    @property
    def periods(self) -> Optional[np.ndarray]:
        f = self.freq
        if f is None:
            return None
        with np.errstate(divide="ignore", invalid="ignore"):
            p = np.where(f > 0.0, 1.0 / f, np.nan)
        return p

    @property
    def n_freq(self) -> int:
        f = self.freq
        return 0 if f is None else int(f.size) 

    
    @property
    def station(self) -> Optional[str]:
        if self.heads is not None and self.heads.head is not None:
            return self.heads.head.station
        if self.Z is not None and getattr(self.Z, "name", None):
            return str(self.Z.name)
        if self.path is not None:
            return self.path.stem
        return None
    
    
    @property
    def name(self) -> Optional[str]:
        # alias with a slightly different precedence
        if self.Z is not None and getattr(self.Z, "name", None):
            return str(self.Z.name)
        if self.heads is not None and self.heads.head is not None:
            return self.heads.head.station
        if self.path is not None:
            return self.path.stem
        return None
    
    
    @property
    def lat(self) -> Optional[float]:
        if self.heads is not None:
            return self.heads.latitude
        return None
    
    
    @property
    def lon(self) -> Optional[float]:
        if self.heads is not None:
            return self.heads.longitude
        return None
    
    
    @property
    def azimuth(self) -> Optional[float]:
        # prefer site AZIMUTH; fall back to header azimuth hint
        if self.heads is not None:
            az = self.heads.azimuth
            if az is not None:
                return az
            # explicit fallback to the hint
            return self.az_hint
        return None
    
    
    @property
    def az_hint(self) -> Optional[float]:
        if (
            self.heads is not None
            and self.heads.head is not None
            and hasattr(self.heads.head, "az_hint")
        ):
            return self.heads.head.az_hint  # type: ignore[attr-defined]
        return None
    
    @property
    def elev(self) -> Optional[float]:
        if self.heads is not None:
            return self.heads.elevation
        return None

    @property 
    def site (self): 
        # site / station
        # be in defensive
        site = (
            getattr(self, "station", None)
            or getattr(self, "site", None)
            or getattr(getattr(self, "heads", None), "station", None)
            or getattr(getattr(self, "head", None), "station", None)
            or getattr(getattr(self, "header", None), "dataid", None)
            or "UNKNOWN"
        )
        return site 
    
    def __has_read__(self) -> bool:
        return bool(self._read_ok)

    def compose_headers(self) -> List[str]:
        out: List[str] = []
        if self.heads is not None:
            out.extend(self.heads.write())
        return out

    def _comps(self) -> tuple[str, dict[str, bool]]:
        z = bool(getattr(self, "Z", None))
        r = bool(getattr(self, "Res", None))
        t = bool(getattr(self, "Tip", None))
        cmap = {"Z": z, "R": r, "T": t}
        cstr = ",".join([k for k, ok in cmap.items() if ok]) or "-"
        return cstr, cmap
    
    def _nfreq(self) -> int:
        f = self.freq
        return 0 if f is None else int(getattr(f, "size", len(f)))
    
    def _pname(self) -> str:
        return self.path.name if self.path else "-"
    
    def _pstr(self) -> str:
        return str(self.path) if self.path else "-"
    
    def _summary_dict(self) -> dict[str, object]:
        def _flt(x: float | None) -> str:
            return f"{x:.5f}" if isinstance(x, (int, float)) else "NA"
        
        cstr, cmap = self._comps()
        return {
            "cls": self.__class__.__name__,
            "site": self.site,
            "nf": self._nfreq(),
            "cstr": cstr,
            "cmap": cmap,
            "lat": _flt(self.lat),
            "lon": _flt(self.lon),
            "az": _flt(self.azimuth),
            "pname": self._pname(),
            "pstr": self._pstr(),
        }

    def __str__(self) -> str:  # noqa: D401
        s = self._summary_dict()
        return (
            f"{s['cls']}(site={s['site']!r}, nfreq={s['nf']}, "
            f"comps={s['cstr']}, lat={s['lat']}, lon={s['lon']}, "
            f"az={s['az']}, path={s['pname']!r})"
        )
    
    def __repr__(self) -> str:  # noqa: D401
        s = self._summary_dict()
        return (
            f"{s['cls']}(\n"
            f"  site={s['site']!r},\n"
            f"  nfreq={s['nf']},\n"
            f"  comps={s['cstr']},\n"
            f"  lat={s['lat']}, lon={s['lon']}, az={s['az']},\n"
            f"  path={s['pstr']!r},\n"
            f")"
        )

