# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0


from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import numpy as np

from ..constants import MU_0
from ..log.logger import get_logger
from ..z.resphase import ResPhase
from ..z.tipper import Tipper
from ..z.z import Z
from .blocks import JBlocks, RBlock, TFBlock
from .config import TENSOR_INDEX
from .heads import Heads

logger = get_logger(__name__)


def _as_path(p: str | Path) -> Path:
    return p if isinstance(p, Path) else Path(p)


__all__ = ["JMixin", "JIOMixin", "JFile"]


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
    >>> jm._complex([1, -2], [0.5, 3]).dtype.kind == "c"
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
    .. [JMixin-1] A. G. Jones, *Magnetotelluric data file J-format*,
           version 2.0, 1994.
    .. [JMixin-2] MTNet, *J format documentation*.
    """

    _tidx: dict[str, tuple[int, int]] = dict(TENSOR_INDEX)

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
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        if a_p.size == 0 or b_p.size == 0:
            z = np.array([])
            zi = np.array([], dtype=int)
            zj = np.array([], dtype=int)
            return z, zi, zj

        ap = a_p.copy()
        bp = b_p.copy()
        idx_a: list[int] = []
        idx_b: list[int] = []
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
    .. [JIOMixin-1] A. G. Jones, *Magnetotelluric data file J-format*,
           version 2.0, 1994.
    .. [JIOMixin-2] MTNet, *J format documentation*.
    """

    def _scan_blocks(
        self,
        j_blocks: JBlocks,
    ) -> dict[str, dict[str, Any]]:
        """
        Converts a JBlocks object into a dictionary of component data.
        """
        comp: dict[str, dict[str, Any]] = {}

        # Operate on the blocks from the passed object
        for blk in j_blocks.blocks:
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
        comp: dict[str, dict[str, Any]],
        *,
        z_obj: Z | None = None,
        tip_obj: Tipper | None = None,
    ) -> tuple[Z | None, Tipper | None, ResPhase | None]:
        """
        Assemble Z, Tipper, and Res/Phase from parsed component dicts.
        Uses zero initialization (no NaNs) for API consistency.
        """

        z_parts: dict[str, dict[str, Any]] = {
            k: d for k, d in comp.items() if k[:1] in {"Z", "Q"}
        }
        r_parts: dict[str, dict[str, Any]] = {
            k: d for k, d in comp.items() if k[:1] in {"R", "S"}
        }
        t_parts: dict[str, dict[str, Any]] = {
            k: d for k, d in comp.items() if k[:1] == "T"
        }

        z_out: Z | None = z_obj
        rp_out: ResPhase | None = None

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
                err = np.asarray(d.get("err", 0.0), dtype=float)
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
            rp_out = ResPhase(freq=freq, name=z_parts[k0]["head"].station)  # type: ignore  # noqa: E501
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
            have_keys = [
                k for k in ("RXX", "RXY", "RYX", "RYY") if k in r_parts
            ]
            if not have_keys:
                have_keys = list(r_parts.keys())
            k0 = have_keys[0]

            p0 = np.asarray(r_parts[k0]["period"], dtype=float)
            freq = self._hz_from_period(p0)
            n = int(p0.size)
            w = 2.0 * math.pi * freq

            def to_imp(d: dict[str, Any]) -> np.ndarray:
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
                p = np.asarray(d["period"], dtype=float)
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
        tip_out: Tipper | None = tip_obj
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
                name=t_parts[k0]["head"].station,
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
        :attr:`JFile.periods`.
    periods : ndarray or None
        Convenience view of ``1.0/freq`` when known.
    n_freq : int
        Number of frequency samples (``0`` if unknown).
    name : str or None
        Friendly site/station name.  Precedence is:
        ``Z.name`` -> head.station -> file stem.
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
    >>> out = jf.write(new_jfn="out.j", datatype="ZR", overwrite=True)
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
    .. [JFile-1] A. G. Jones, *Magnetotelluric data file J-format*,
           version 2.0, 1994.
    .. [JFile-2] MTNet, *J format documentation*.
    """

    def __init__(
        self,
        path: str | Path | None = None,
        *,
        verbose: int = 0,
    ) -> None:
        self.verbose = verbose
        self.path: Path | None = _as_path(path) if path else None
        self.heads: Heads | None = None
        self.blocks: JBlocks | None = None
        self.Z: Z | None = None
        self.Tip: Tipper | None = None
        self.Res: ResPhase | None = None
        self._read_ok: bool = False

        if self.path is not None:
            self.read(self.path)

    @classmethod
    def from_file(cls, path: str | Path, *, verbose: int = 0) -> JFile:
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
            Instance with :attr:`JFile.heads`, :attr:`JFile.blocks`
            and objects (:attr:`JFile.Z`, :attr:`JFile.Tip`,
            :attr:`JFile.Res`) populated where possible.

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
        .. [JFile-from-file-1] A. G. Jones, *Magnetotelluric data file
               J-format*, version 2.0, 1994.
        """
        inst = cls(path, verbose=verbose)
        # inst.read(path)
        return inst

    def read(
        self,
        path: str | Path | None = None,
        *,
        start: int | None = None,
    ) -> JFile:
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
        .. [JFile-read-1] A. G. Jones, *Magnetotelluric data file
               J-format*, version 2.0, 1994.
        """

        if path is not None:
            self.path = _as_path(path)
        if self.path is None:
            raise ValueError("path is required")

        # Parse headers and all blocks ONCE.
        self.heads = Heads.from_file(self.path, verbose=self.verbose)
        self.blocks = JBlocks.from_file(self.path, verbose=self.verbose)

        # Scan the already-parsed blocks object
        # (no re-reading from disk).
        comp = self._scan_blocks(self.blocks)

        # Build the final Z, Tip, and Res objects.
        z, tip, rp = self._build_from_comp(comp)
        self.Z = z
        self.Tip = tip
        self.Res = rp

        self._read_ok = True
        return self

    def write(
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
            :attr:`JFile.path` or default to ``'out.j'``.
        new_jfn : str, optional
            Replacement filename.  Takes precedence over
            ``j_fn`` when provided.
        datatype : {'Z','R','T','ZR','RT','ZT','ZRT','ALL'}, optional
            Select families to emit.  If ``None``, the writer
            auto-detects from available objects on the instance.
        savepath : str or Path, optional
            Folder where to save.  Defaults to the parent of
            :attr:`JFile.path` or the current directory.
        verbose : int, optional
            Override verbosity.  If ``None``, reuse
            :attr:`JFile.verbose`.
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
        >>> out = jf.write(new_jfn="site_out.j", datatype="ZR", overwrite=True)
        >>> isinstance(out, str)
        True

        See Also
        --------
        JFile.compose_headers : Render banner + headers only.
        Z, Tipper, ResPhase : Sources used by the writer.

        References
        ----------
        .. [JFile-write-1] A. G. Jones, *Magnetotelluric data file
               J-format*, version 2.0, 1994.
        """

        if not self.__has_read__():
            raise RuntimeError("Cannot write JFile; call .read() first.")

        def _ensure_parent(p: Path) -> None:
            p.parent.mkdir(parents=True, exist_ok=True)

        lines: list[str] = []
        vb = self.verbose if verbose is None else int(verbose)

        # --- 1. Formatting Helpers for Standard J-Format Output ---
        def _fmt(val: float, width: int, precision: int) -> str:
            """Formats a number with fixed width and precision."""
            if not np.isfinite(val):
                return f"{'NaN':>{width}}"
            # Space for sign pad, G for general format
            return f"{val:{width}.{precision}G}"

        def _fmt_sci(val: float, width: int, precision: int) -> str:
            """Formats a number in scientific notation."""
            if not np.isfinite(val):
                return f"{'NaN':>{width}}"
            # Format and clean up E+0 -> E+, E-0 -> E-
            s = f"{val:{width}.{precision}E}"
            return s.replace("E+0", "E+").replace("E-0", "E-")

        # --- 2. Write Main Headers (Banner & Info) ---
        if self.heads:
            lines.extend(
                self.heads.banner.write(new=True, include_origin=True)
            )
            lines.extend(self.heads.info.write())

        # --- 3. Prepare to Write Data Blocks ---
        first_block_written = False
        station_line = self.station or "UNKNOWN"
        if self.az_hint is not None:
            station_line += f"  {self.az_hint:g}"

        have = {
            "Z": self.Z is not None,
            "R": self.Res is not None,
            "T": self.Tip is not None,
        }
        key = (datatype or "ZRT").upper().strip()
        flags = {
            "Z": "Z" in key and have["Z"],
            "R": "R" in key and have["R"],
            "T": "T" in key and have["T"],
        }

        # --- 4. Write R-Blocks (Resistivity/Phase) ---
        if flags.get("R") and self.Res and self.periods is not None:
            for code, (i, j) in self._tidx.items():
                rho = self.Res.resistivity[:, i, j]
                if np.all(np.abs(rho) < 1e-12):
                    continue

                if not first_block_written:
                    lines.append(station_line)
                    first_block_written = True

                lines.append(f"R{code}")
                lines.append(str(self.n_freq))
                pha = self.Res.phase[:, i, j]

                for k in range(self.n_freq):
                    p, r_val, ph_val = (self.periods[k], rho[k], pha[k])
                    # Create nicely aligned 9-column output
                    row = (
                        f"{_fmt_sci(p, 11, 4)}"
                        f"{_fmt(r_val, 11, 4)}"
                        f"{_fmt(ph_val, 11, 2)}"
                        f"{_fmt(r_val, 11, 4)}"
                        f"{_fmt(r_val, 11, 4)}"
                        f"{_fmt(ph_val, 11, 2)}"
                        f"{_fmt(ph_val, 11, 2)}"
                        f"{_fmt(1.0, 9, 2)}"
                        f"{_fmt(1.0, 9, 2)}"
                    )
                    lines.append(row)

        # --- 5. Write T-Blocks (Tipper) ---
        if (
            flags.get("T")
            and self.Tip
            and self.Tip.tipper is not None
            and self.periods is not None
        ):
            for code, k in {"ZX": 0, "ZY": 1}.items():
                t_comp = self.Tip.tipper[:, 0, k]
                if np.all(np.abs(t_comp) < 1e-12):
                    continue
                if not first_block_written:
                    lines.append(station_line)
                    first_block_written = True

                lines.append(f"T{code}")
                lines.append(str(self.n_freq))
                err = self.Tip.tipper_err
                t_err = (
                    err[:, 0, k] if err is not None else np.zeros(self.n_freq)
                )

                for l_idx in range(self.n_freq):
                    p, t, e = (
                        self.periods[l_idx],
                        t_comp[l_idx],
                        t_err[l_idx],
                    )
                    row = (
                        f"{_fmt_sci(p, 11, 4)}"
                        f"{_fmt(t.real, 12, 4)}"
                        f"{_fmt(t.imag, 12, 4)}"
                        f"{_fmt(e, 12, 4)}"
                        f"{_fmt(1.0, 9, 2)}"
                    )
                    lines.append(row)

        # --- 6. Write Z-Blocks (Impedance) ---
        if (
            flags.get("Z")
            and self.Z
            and self.Z.z is not None
            and self.periods is not None
        ):
            for comp_code, (i, j) in self._tidx.items():
                z_comp = self.Z.z[:, i, j]
                if np.all(np.abs(z_comp) < 1e-12):
                    continue
                if not first_block_written:
                    lines.append(station_line)
                    first_block_written = True

                lines.append(f"Z{comp_code} SI")
                lines.append(str(self.n_freq))
                err = self.Z.z_err
                z_err = (
                    err[:, i, j] if err is not None else np.zeros(self.n_freq)
                )

                for k in range(self.n_freq):
                    p, z, e = (self.periods[k], z_comp[k], z_err[k])
                    row = (
                        f"{_fmt_sci(p, 11, 4)}"
                        f"{_fmt(z.real, 12, 4)}"
                        f"{_fmt(z.imag, 12, 4)}"
                        f"{_fmt(e, 12, 4)}"
                        f"{_fmt(1.0, 9, 2)}"
                    )
                    lines.append(row)

        # --- 7. Finalize and Save File ---
        path = self.path
        base = Path(new_jfn or j_fn or (path.name if path else "out.j"))
        folder = Path(
            savepath or (path.parent if path and path.parent.exists() else ".")
        )
        out_path = folder / base

        if not overwrite and out_path.exists():
            stem, suf = out_path.stem, out_path.suffix or ".j"
            c = 1
            while out_path.exists():
                out_path = out_path.with_name(f"{stem}_{c}{suf}")
                c += 1

        _ensure_parent(out_path)
        out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        if vb:
            logger.info("Wrote J file: %s", out_path)

        return str(out_path)

    @property
    def freq(self) -> np.ndarray | None:
        if self.Z is not None and hasattr(self.Z, "freq"):
            return self.Z.freq  # type: ignore[return-value]
        if self.Tip is not None and hasattr(self.Tip, "freq"):
            return self.Tip.freq  # type: ignore[return-value]
        if self.Res is not None and hasattr(self.Res, "freq"):
            return self.Res.freq  # type: ignore[return-value]
        return None

    @property
    def periods(self) -> np.ndarray | None:
        f = self.freq
        return None if f is None else self._hz_from_period(f)

    @property
    def n_freq(self) -> int:
        f = self.freq
        return 0 if f is None else int(f.size)

    @property
    def station(self) -> str | None:
        if self.heads is not None and self.heads.head is not None:
            return self.heads.head.station
        if self.Z is not None and getattr(self.Z, "name", None):
            return str(self.Z.name)
        if self.path is not None:
            return self.path.stem
        return None

    @property
    def name(self) -> str | None:
        # alias with a slightly different precedence
        if self.Z is not None and getattr(self.Z, "name", None):
            return str(self.Z.name)
        if self.heads is not None and self.heads.head is not None:
            return self.heads.head.station
        if self.path is not None:
            return self.path.stem
        return None

    @property
    def lat(self) -> float | None:
        if self.heads is not None:
            return self.heads.latitude
        return None

    @property
    def lon(self) -> float | None:
        if self.heads is not None:
            return self.heads.longitude
        return None

    @property
    def azimuth(self) -> float | None:
        # prefer site AZIMUTH; fall back to header azimuth hint
        if self.heads is not None:
            az = self.heads.azimuth
            if az is not None:
                return az
            # explicit fallback to the hint
            return self.az_hint
        return None

    @property
    def az_hint(self) -> float | None:
        if (
            self.heads is not None
            and self.heads.head is not None
            and hasattr(self.heads.head, "az_hint")
        ):
            return self.heads.head.az_hint  # type: ignore[attr-defined]
        return None

    @property
    def elev(self) -> float | None:
        if self.heads is not None:
            return self.heads.elevation
        return None

    @property
    def site(self):
        # site / station
        # be in defensive
        return self.station or "UNKNOWN"

    def __has_read__(self) -> bool:
        return bool(self._read_ok)

    def compose_headers(self) -> list[str]:
        out: list[str] = []
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
