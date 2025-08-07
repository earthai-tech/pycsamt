# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""
EDI container and orchestrator.

This module provides the v2 `Edi` class, a façade that
assembles and coordinates all EDI components (HEAD, INFO,
DEFINEMEAS, MT/EMAP, SPECTRA, TIME SERIES) and numeric data
blocks (FREQ, Z, Tipper, RHO/PHASE). It exposes a consistent
API to read, manipulate, and write SEG EDI files while
keeping per-section parsing encapsulated within their own
classes.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import datetime
import numpy as np

from ..exceptions import (
    EdIDataError,
    HeaderError,
    FileHandlingError,
)
from ..z.z import Z as EMZ
from ..z.z import Tipper
from ..log.logger import get_logger

from .base import EdiFileBase
from .properties import IsEdi, _ZRP_COMPS, _TIP_COMPS
from .heads import Head, Info
from .definemeas import DefineMeasurement
from .mtemap import MTEMAP
from .spectra import Spectra
from .time_series import TimeSeries

__all__ = ["Edi"]

logger = get_logger(__name__)


class Edi(EdiFileBase):
    """
    High-level container for SEG EDI files (v2 API).

    The `Edi` class coordinates per-section parsers and the
    reading/writing of numeric data blocks. It provides a
    uniform API:

    - :meth:`read` loads the EDI file and parses sections.
    - :meth:`read_data_blocks` parses FREQ/Z/Tipper/RHO/PHASE
      blocks and populates `Z` and `Tipper`.
    - :meth:`compose` builds the text of all structural
      sections; numeric blocks are added by write paths.
    - :meth:`write` (or `write_edifile`) writes a complete EDI.
    - :meth:`interpolate` re-samples the impedance tensor over
      a new frequency grid and returns a new `Z` object.
    - :meth:`write_new_edifile` writes a new EDI using a given
      `Z` (without mutating `self.Z`).

    Protocol
    --------
    1) **Validation**: `IsEdi._assert_edi(path, deep=True)`
       checks extension and sentinel lines (`>HEAD`, `>END`).

    2) **Section parsing**:
       - `Head.from_edi`
       - `Info.from_edi`
       - `DefineMeasurement.from_edi`
       - `MTEMAP.from_edi` (records the first data-line index)
       - Optional: `Spectra.from_edi`, `TimeSeries.from_edi`

       Each section class owns its parsing and rendering
       logic. `Edi` simply stores the resulting objects in
       an internal registry (`add_section/get_section`).

    3) **Numeric blocks**:
       - `_get_specific_comp()` scans from the data-section
         start and collects numeric payload into a dict of
         lists: `{'freq': [...], 'zxxr': [...], ...}`.
       - `_fill_data_array()` converts lists to arrays:
         `self.Z._freq`, `self.Z._z`, `self.Z._z_err`,
         `self.Z.rotation_angle`; optional `Tipper` arrays
         if present. It reverses frequency order to high→low
         and converts "empty" values (default `1.0e32`) to
         zero. Then it computes resistivity and phase.

    4) **Composition and writing**:
       - `compose()` concatenates section writers
         (`Head.write()`, `Info.write()`, etc.).
       - Write helpers append numeric block headers (with
         appropriate `ROT=ZROT/TROT/NONE`) and packed lines
         of values with controlled block width and scientific
         notation in uppercase.

    Attributes
    ----------
    path : Path | None
        EDI file path to read or write.
    sections : dict[str, object]
        Registry of parsed components; keys include
        ``head``, ``info``, ``definemeasurement``,
        ``mtemap``, optionally ``spectra``, ``timeseries``.
    Z : EMZ
        Impedance container. Populated by
        :meth:`read_data_blocks`.
    Tip : Tipper
        Tipper container. Populated if tipper blocks exist.
    block_size : int
        Number of numeric values per line in data blocks.
    data_header_template : str
        Template for data section comment headers,
        e.g. ``'>!****{0}****!\\n'``.
    float_format : str
        Numeric format for data values, e.g. ``' 15.6e'``.
    _data_start_line : int | None
        Line index where numeric data blocks begin.

    Methods (summary)
    -----------------
    read(path: str | Path | None) -> Edi
        Parse sections and then call `read_data_blocks()`.
    read_data_blocks() -> Edi
        Parse numeric datasets and populate `Z`/`Tipper`.
    compose() -> str
        Render all structural sections (no numeric data).
    write(...): str
        Write a complete EDI file to disk.
    interpolate(new_freq, kind='slinear', ...) -> EMZ
        Interpolate impedance tensor onto a new frequency
        grid and compute derived quantities.
    write_new_edi(new_edi_fn=None, new_Z=None, **kws) -> str
        Write a new EDI using a provided `Z` object without
        mutating `self.Z`.

    Examples
    --------
    Read and access metadata::

        from pycsamt.seg.edi import Edi
        edi = Edi(path="data/S00.edi").read()
        print(edi.get_section("head").dataid)
        print(edi.Z.freq.shape, edi.Z.z.shape)

    Interpolate impedance and write a new EDI::

        new_f = np.logspace(-3, 3, 24)
        z_new = edi.interpolate(new_f, kind="slinear")
        out = edi.write_new("S00_interp.edi", new_Z=z_new)

    Convert to EMAP-style write (rho/phase blocks)::

        # when composing EMAP, writer adds RHO/PHS blocks;
        # empty values are formatted as 1.0E+32 as per SEG.
        out = edi.write(new_edi_fn="S00_emap.edi", datatype="emap")

    Notes
    -----
    * Frequencies are normalized to high→low order.
    * Empty values are represented by `1.0E+32` in files and
      mapped to `0.0` in memory for convenience.
    * Rotation blocks: Z uses ``ZROT``, Tipper uses ``TROT``,
      EMAP scalar blocks use ``ROT=NONE``.
    * `Spectra` and `TimeSeries` are optional. If present,
      class methods `from_edi` and `write` will handle them.
    * The writer enforces consistent numeric formatting and
      line wrapping with the configured `block_size`.

    References
    ----------
    Wight, D. E., & Bostick Jr., F. X. (1987/1988). SEG
    MT/EMAP Data Interchange Standard. Society of Exploration
    Geophysicists.
    MTpy project: https://github.com/MTgeophysics/mtpy
    """
    def __init__(
        self,
        path: Optional[str | Path] = None,
        *,
        verbose: int = 0,
        **kwargs,
    ) -> None:
        super().__init__(path=Path(path) if path else None, **kwargs)

        self.verbose = int(verbose)

        self.block_size = 6
        self.number_fmt = " 15.6e"
        self.data_header_tpl = ">!****{title}****!\n"

        # placeholders for data sections
        self._comp_raw: dict[str, list[float]] = {}
        self._data_start_line: Optional[int] = None

        # numeric components will be built later
        self.Z = EMZ()
        self.Tip = Tipper()

        # optional sections (added when present)
        # Head/Info/DefineMeasurement/MTEMAP/Spectra/TimeSeries
        # will be registered in read()

        if self.path is not None:
            self.read()

    # Lifecycle
    def read(self, path=None ) -> "Edi":
        """
        Read self.path and populate sections. Spectra and
        TimeSeries are optional and added when found.
        """
        self.path = path or self.path 
        
        if self.path is None:
            raise FileHandlingError("No EDI path provided.")

        # Validate EDI file (extension + sentinel lines).
        IsEdi._assert_edi(str(self.path), deep=True)

        # Core sections (required)
        head = Head.from_edi(self.path)
        self.add_section("head", head)

        info = Info.from_edi(self.path)
        self.add_section("info", info)

        dmeas = DefineMeasurement.from_edi(self.path)
        self.add_section("definemeasurement", dmeas)

        mtemap = MTEMAP.from_edi(self.path)
        self.add_section("mtemap", mtemap)

        # Keep where data blocks start (if exposed by MTEMAP).
        self._data_start_line = getattr(
            mtemap, "start_data_line", None
        )

        # Optional sections
        spectra: Optional[Spectra] = None
        try:
            spectra = Spectra.from_edi(self.path)
        except (EdIDataError, HeaderError, FileHandlingError):
            spectra = None
        if spectra is not None:
            self.add_section("spectra", spectra)

        tser: Optional[TimeSeries] = None
        try:
            tser = TimeSeries.from_edi(self.path)
        except (EdIDataError, HeaderError, FileHandlingError):
            tser = None
        if tser is not None:
            self.add_section("timeseries", tser)

        # Impedance/tipper blocks and other numeric datasets are
        # parsed later in _read_data_blocks().
        self.read_data_blocks()
        
        return self
        
    def read_data_blocks(self) -> "Edi":
        """
        Parse numeric data blocks and populate Z/Tip arrays.
        """
        comp = self._get_specific_comp(
            edifile=str(self.path) if self.path else None,
            data_sect_line=self._data_start_line,
        )
        self._fill_data_array(comp)
        return self

    def compose(self) -> str:
        """
        Build the full EDI text. Will be completed as we
        refactor write paths for each section and data.
        """
        lines: list[str] = []

        head = self.get_section("head")
        if head is not None:
            lines.extend(head.write())

        info = self.get_section("info")
        if info is not None:
            lines.extend(info.write())

        dmeas = self.get_section("definemeasurement")
        if dmeas is not None:
            lines.extend(dmeas.write())

        mtemap = self.get_section("mtemap")
        if mtemap is not None:
            lines.extend(mtemap.write())

        spectra = self.get_section("spectra")
        if spectra is not None:
            lines.extend(spectra.write())

        tser = self.get_section("timeseries")
        if tser is not None:
            lines.extend(tser.write())

        # data blocks (freq/Z/tipper/rho/phase) will be
        # appended in a later step.
        return "".join(lines)
    

    # Data blocks (freq, Z, tipper, rho/phase) parsing
    def _get_specific_comp(
        self,
        edifile: Optional[str | Path] = None,
        data_sect_line: Optional[int] = None,
    ) -> dict[str, list[float]]:
        """
        Parse numeric blocks into a dict of lists:
        {'freq': [...], 'zxxr': [...], 'zxxi': [...], ...}.
        """
        if edifile is not None:
            self.path = Path(edifile)

        if self.path is None:
            raise FileHandlingError(
                "No EDI path provided."
            )

        IsEdi._assert_edi(str(self.path), deep=True)

        with self.path.open("r", encoding="utf-8") as f:
            lines = f.readlines()

        start = (
            int(data_sect_line)
            if data_sect_line is not None
            else (self._data_start_line or 0)
        )

        if start <= 0:
            # fallback: find first data block header like '>FREQ'
            for i, ln in enumerate(lines):
                if ln.lstrip().startswith(">FREQ"):
                    start = i
                    break

        if start <= 0 or start >= len(lines):
            raise EdIDataError(
                "Could not locate the data section start."
            )

        # known keys
        z_pairs = ("zxx", "zxy", "zyx", "zyy")
        rho_prefix = ("rho", "frho")
        phs_prefix = ("phs", "fphs")
        tip_keys = {
            "txr.exp", "txi.exp", "txvar.exp",
            "tyr.exp", "tyi.exp", "tyvar.exp",
        }
        rot_keys = {"zrot", "trot"}

        def _is_known_key(k: str) -> bool:
            k = k.lower()
            if k == "freq" or k in tip_keys or k in rot_keys:
                return True
            if any(k.startswith(p) for p in rho_prefix):
                return True
            if any(k.startswith(p) for p in phs_prefix):
                return True
            if any(
                k in (f"{c}r", f"{c}i", f"{c}.var")
                for c in z_pairs
            ):
                return True
            return False

        empty_val = 1.0e32
        head = self.get_section("head")
        if head is not None and getattr(head, "empty", None):
            try:
                empty_val = float(head.empty)
            except Exception:
                pass

        comp: dict[str, list[float]] = {}
        current: Optional[str] = None

        for raw in lines[start:]:
            s = raw.strip()
            if not s:
                continue

            if s.startswith(">"):
                # skip comment/data headers like >!****...****!
                if s.startswith(">!"):
                    current = None
                    continue
                # parse key token after '>'
                toks = s[1:].strip().split()
                if not toks:
                    current = None
                    continue
                key = toks[0].lower()
                if _is_known_key(key):
                    current = key
                    comp.setdefault(current, [])
                else:
                    current = None
                continue

            # accumulate numeric payload
            if current is not None:
                vals = []
                for tok in s.split():
                    try:
                        v = float(tok)
                        if v == empty_val:
                            v = 0.0
                        vals.append(v)
                    except Exception:
                        vals.append(0.0)
                comp[current].extend(vals)

        self._comp_raw = comp
        
        return comp

    def _fill_data_array(
        self,
        data_dict: Optional[dict[str, list[float]]] = None,
    ) -> None:
        """
        Build Z/Tipper arrays from parsed dict, normalize
        frequency order, and compute derived quantities.
        """
        comp = data_dict if data_dict is not None else self._comp_raw
        if not comp:
            raise EdIDataError("No parsed data to fill.")

        if "freq" not in comp or not comp["freq"]:
            raise EdIDataError("Missing 'FREQ' data block.")

        freq = np.asarray(comp["freq"], dtype=float)

        # high → low; reverse if ascending
        rev = bool(freq[-1] > freq[0])

        # Z and variance
        z = np.zeros((freq.size, 2, 2), dtype=np.complex128)
        zvar = np.zeros((freq.size, 2, 2), dtype=float)

        def _get(name: str) -> np.ndarray:
            return np.asarray(comp.get(name, []), dtype=float)

        # helper: fill one component
        def _fill(cc: str, i: int, j: int) -> None:
            r = _get(f"{cc}r")
            im = _get(f"{cc}i")
            vr = _get(f"{cc}.var")
            n = min(r.size, im.size) if im.size else r.size
            if n:
                z[:n, i, j] = r[:n] + 1j * (im[:n] if im.size else 0.0)
            if vr.size:
                zvar[:vr.size, i, j] = np.sqrt(np.abs(vr[:vr.size]))

        _fill("zxx", 0, 0)
        _fill("zxy", 0, 1)
        _fill("zyx", 1, 0)
        _fill("zyy", 1, 1)

        if rev:
            freq = freq[::-1]
            z = z[::-1]
            zvar = zvar[::-1]

        self.Z._freq = freq
        self.Z._z = z
        self.Z._z_err = zvar

        zrot = _get("zrot")
        if zrot.size == freq.size:
            self.Z.rotation_angle = zrot[::-1] if rev else zrot
        else:
            self.Z.rotation_angle = np.zeros_like(freq, dtype=float)

        self.Z.compute_resistivity_phase()

        # tipper (optional)
        has_tip = any(k.startswith("t") for k in comp.keys())
        if not has_tip:
            return

        tip = np.zeros((freq.size, 1, 2), dtype=np.complex128)
        tip_err = np.zeros((freq.size, 1, 2), dtype=float)

        txr = _get("txr.exp")
        txi = _get("txi.exp")
        tyr = _get("tyr.exp")
        tyi = _get("tyi.exp")
        txv = _get("txvar.exp")
        tyv = _get("tyvar.exp")

        n0 = min(txr.size, txi.size) if txi.size else txr.size
        n1 = min(tyr.size, tyi.size) if tyi.size else tyr.size

        if n0:
            tip[:n0, 0, 0] = txr[:n0] + 1j * (txi[:n0] if txi.size else 0.0)
        if n1:
            tip[:n1, 0, 1] = tyr[:n1] + 1j * (tyi[:n1] if tyi.size else 0.0)

        if txv.size:
            tip_err[:txv.size, 0, 0] = np.sqrt(np.abs(txv[:txv.size]))
        if tyv.size:
            tip_err[:tyv.size, 0, 1] = np.sqrt(np.abs(tyv[:tyv.size]))

        if rev:
            tip = tip[::-1]
            tip_err = tip_err[::-1]

        self.Tip._freq = freq
        self.Tip._tipper = tip
        self.Tip._tipper_err = tip_err

        trot = _get("trot")
        if trot.size == freq.size:
            self.Tip.rotation_angle = trot[::-1] if rev else trot
        elif zrot.size == freq.size:
            self.Tip.rotation_angle = (
                zrot[::-1] if rev else zrot
            )
        else:
            self.Tip.rotation_angle = np.zeros_like(freq, dtype=float)

        self.Tip.compute_amp_phase()
        self.Tip.compute_mag_direction()

    def write(
        self,
        edi_fn: str | None = None,
        new_edifilename: str | None = None,
        datatype: str | None = None,
        savepath: str | Path | None = None,
        add_filter_array: np.ndarray | None = None,
        **kwargs,
    ) -> str:
        """
        Write an EDI file assembled from the current object.
        """
        verbose = kwargs.pop("verbose", None)
        if verbose is not None:
            self.verbose = verbose
        if getattr(self, "verbose", None) is None:
            self.verbose = 0

        if edi_fn is not None:
            self.edifile = edi_fn

        if self.Head.dataid is None:
            # ensure sections are populated when writing from source file
            if self.edifile:
                self.read_edi(self.edifile)

        # output directory
        out_dir = (
            Path(savepath).expanduser().resolve()
            if savepath is not None
            else (Path.cwd() / "edi_out")
        )
        out_dir.mkdir(parents=True, exist_ok=True)

        # choose filename
        now_year = datetime.datetime.utcnow().year
        sect_lines = self.MTEMAP.write(
            nfreq=(self.Z.freq.size if self.Z.freq is not None else None)
        )
        first_line = sect_lines[0].strip()
        inferred_type = "mt" if "MTSECT" in first_line else "emap"

        if datatype is not None:
            dt = datatype.strip().lower()
            if dt not in {"mt", "emap"}:
                raise EdIDataError(
                    "datatype must be 'mt' or 'emap'."
                )
            self.typefile = dt
        else:
            self.typefile = inferred_type

        if new_edifilename:
            name = new_edifilename
            if not name.lower().endswith(".edi"):
                name += ".edi"
        else:
            base = (
                f"new_{Path(self.edifile).name}"
                if self.edifile
                else f"{self.Head.dataid or 'site'}_"
                     f"{self.typefile}.{now_year}.edi"
            )
            name = base

        out_path = out_dir / name

        # section headers
        head_lines = self.Head.write()
        info_lines = self.Info.write()
        dmeas_lines = self.DefineMeasurement.write()
        mtemap_lines = sect_lines + ["\n"]

        # frequencies
        freq_head = self.data_head_com.format(
            "FREQUENCIES".upper()
        )
        freq_lines = [freq_head] + self._write_components_blocks(
            edi_datacomp=self.Z.freq, comp_key="freq"
        )

        # Z rotation (MT)
        zrot_lines: list[str] = []
        if self.typefile == "mt":
            zrot_head = self.data_head_com.format(
                "IMPEDANCE ROTATION ANGLES".upper()
            )
            angles = self.Z.rotation_angle
            if np.ndim(angles) == 0 or (
                hasattr(angles, "size") and angles.size == 1
            ):
                angles = np.full_like(
                    self.Z.freq, float(angles), dtype=float
                )
            zrot_lines = [zrot_head] + self._write_components_blocks(
                edi_datacomp=np.asarray(angles, dtype=float),
                comp_key="zrot",
            )

        # impedance blocks
        z_head = self.data_head_com.format("IMPEDANCES".upper())
        z_data_lines = [z_head]
        z = np.nan_to_num(self.Z.z)
        zerr = np.nan_to_num(self.Z.z_err)

        for i in range(2):
            for j in range(2):
                # real
                zreal = self._write_components_blocks(
                    edi_datacomp=z[:, i, j].real,
                    comp_key=_ZRP_COMPS[0][i * 2 + j][0],
                )
                # imag
                zimag = self._write_components_blocks(
                    edi_datacomp=z[:, i, j].imag,
                    comp_key=_ZRP_COMPS[0][i * 2 + j][1],
                )
                # variance (square of error)
                zvar = self._write_components_blocks(
                    edi_datacomp=np.square(zerr[:, i, j]),
                    comp_key=_ZRP_COMPS[0][i * 2 + j][2],
                )
                z_data_lines.extend(zreal + ["\n"])
                z_data_lines.extend(zimag + ["\n"])
                z_data_lines.extend(zvar + ["\n"])

        # EMAP: resistivity/phase (+ optional filtered)
        rhophs_lines: list[str] = []
        if self.typefile == "emap":
            rhophs_head = self.data_head_com.format(
                "RESISTIVITIES AND PHASES".upper()
            )
            rhophs_lines = [rhophs_head]

            rho = np.zeros((self.Z.freq.size, 2, 4), float)
            phs = np.zeros((self.Z.freq.size, 2, 4), float)

            # xy, xx order as in _ZRP_COMPS[1] and [2]
            rho[:, 0, 0] = np.nan_to_num(self.Z.res_xx)
            rho[:, 0, 1] = np.nan_to_num(self.Z.res_err_xx) ** 2
            rho[:, 0, 2] = np.nan_to_num(self.Z.res_err_xx)
            rho[:, 1, 0] = np.nan_to_num(self.Z.res_xy)
            rho[:, 1, 1] = np.nan_to_num(self.Z.res_err_xy) ** 2
            rho[:, 1, 2] = np.nan_to_num(self.Z.res_err_xy)

            phs[:, 0, 0] = np.nan_to_num(self.Z.phase_xx)
            phs[:, 0, 1] = np.nan_to_num(self.Z.phase_err_xx) ** 2
            phs[:, 0, 2] = np.nan_to_num(self.Z.phase_err_xx)
            phs[:, 1, 0] = np.nan_to_num(self.Z.phase_xy)
            phs[:, 1, 1] = np.nan_to_num(self.Z.phase_err_xy) ** 2
            phs[:, 1, 2] = np.nan_to_num(self.Z.phase_err_xy)

            for ii in range(2):
                for jj in range(4):
                    if np.all(rho[:, ii, jj] == 0.0) and np.all(
                        phs[:, ii, jj] == 0.0
                    ):
                        continue
                    r_lines = self._write_components_blocks(
                        edi_datacomp=rho[:, ii, jj],
                        comp_key=_ZRP_COMPS[1][ii][jj],
                    )
                    p_lines = self._write_components_blocks(
                        edi_datacomp=phs[:, ii, jj],
                        comp_key=_ZRP_COMPS[2][ii][jj],
                    )
                    rhophs_lines.extend(r_lines + ["\n"])
                    rhophs_lines.extend(p_lines + ["\n"])

            # optional filtered rho
            if add_filter_array is not None:
                frho = np.zeros_like(rho)
                # example: caller supplies shape (nf, 2, 2). Map to [ii,0]
                if (
                    add_filter_array.ndim == 3
                    and add_filter_array.shape[1:] == (2, 2)
                    and add_filter_array.shape[0] == self.Z.freq.size
                ):
                    frho[:, 1, 0] = add_filter_array[:, 0, 1]
                    frho[:, 0, 1] = add_filter_array[:, 1, 0]

                    for ii in range(2):
                        for jj in range(4):
                            f_lines = self._write_components_blocks(
                                edi_datacomp=frho[:, ii, jj],
                                comp_key=_ZRP_COMPS[3][ii][jj],
                            )
                            rhophs_lines.extend(f_lines + ["\n"])

            # clean trailing extra newline from Z block when EMAP
            if z_data_lines and z_data_lines[-1] == "\n":
                z_data_lines = z_data_lines[:-1]

        # tipper (optional)
        tip_rot_lines: list[str] = []
        tip_data_lines: list[str] = []
        tip = getattr(self.Tip, "tipper", None)

        if tip is not None and tip.size and not np.all(tip == 0.0):
            trot = getattr(self.Tip, "rotation_angle", None)
            if trot is None:
                trot = np.zeros_like(self.Z.freq, dtype=float)
            if np.ndim(trot) == 0 or (
                hasattr(trot, "size") and trot.size == 1
            ):
                trot = np.full_like(self.Z.freq, float(trot), dtype=float)

            tip_rot_head = self.data_head_com.format(
                "TIPPER ROTATION ANGLES".upper()
            )
            tip_rot_lines = [tip_rot_head] + self._write_components_blocks(
                edi_datacomp=np.asarray(trot, dtype=float),
                comp_key="trot",
            )

            tip_head = self.data_head_com.format(
                "TIPPER PARAMETERS".upper()
            )
            tip_data_lines = [tip_head]

            # tx / ty real, imag, variance
            for comp_idx, keys in enumerate(_TIP_COMPS):
                # 0: real, 1: imag, 2: var
                treal = self._write_components_blocks(
                    edi_datacomp=tip[:, 0, comp_idx].real,
                    comp_key=keys[0],
                )
                timag = self._write_components_blocks(
                    edi_datacomp=tip[:, 0, comp_idx].imag,
                    comp_key=keys[1],
                )
                terr = getattr(self.Tip, "_tipper_err", None)
                tvar_arr = (
                    np.square(terr[:, 0, comp_idx])
                    if terr is not None
                    else np.zeros_like(self.Z.freq, dtype=float)
                )
                tvar = self._write_components_blocks(
                    edi_datacomp=tvar_arr,
                    comp_key=keys[2],
                )
                tip_data_lines.extend(treal)
                tip_data_lines.extend(timag)
                tip_data_lines.extend(tvar)
                tip_data_lines.append("\n")

        # assemble
        out_lines: list[str] = []
        out_lines += head_lines
        out_lines += info_lines
        out_lines += dmeas_lines
        out_lines += mtemap_lines

        if self.typefile == "mt":
            for block in [
                freq_lines,
                zrot_lines,
                z_data_lines,
                tip_rot_lines,
                tip_data_lines,
            ]:
                if block:
                    out_lines += block + ["\n"]

        if self.typefile == "emap":
            for block in [freq_lines, z_data_lines, rhophs_lines]:
                if block:
                    out_lines += block + ["\n"]

        out_lines.append(">END")

        # write to disk
        with out_path.open("w", encoding="utf-8") as fw:
            fw.writelines(out_lines)

        if self.verbose > 0:
            logger.info(
                "EDI written: %s", str(out_path)
            )

        return str(out_path)

    def _write_components_blocks(
        self,
        edi_datacomp: np.ndarray,
        comp_key: str,
        datatype: str | None = None,
    ) -> list[str]:
        """
        Build a data block for a given component.
        """
        logger.info("Write component block: %s", comp_key)

        if edi_datacomp is None:
            raise EdIDataError(
                f"No data for component '{comp_key}'."
            )

        arr = np.asarray(edi_datacomp)
        if arr.ndim != 1:
            arr = arr.ravel()

        dtype = (datatype or getattr(self, "typefile", None)
                 or "mt").lower()

        rot_tags = ("ROT=ZROT", "ROT=TROT", "ROT=NONE")
        key = comp_key.strip().lower()

        # header line
        if key == "freq":
            header = f">FREQ  //{arr.size}\n"
        elif dtype in {"mt", ">=mtsect", "mtsect"}:
            if key in {"zrot", "trot"}:
                header = f">{key.upper()}  //{arr.size}\n"
            elif key.startswith("z"):
                header = (f">{key.upper()} {rot_tags[0]}  //"
                          f"{arr.size}\n")
            elif key.startswith("t"):
                header = (f">{key.upper()} {rot_tags[1]}  //"
                          f"{arr.size}\n")
            else:
                raise EdIDataError(
                    f"Unsupported MT key '{comp_key}'."
                )
        elif dtype in {"emap", ">=emapsect", "emapsect"}:
            if key == "freq":
                header = f">FREQ  //{arr.size}\n"
            else:
                header = (f">{key.upper()} {rot_tags[2]}  //"
                          f"{arr.size}\n")
        else:
            raise EdIDataError(
                f"datatype must be 'mt' or 'emap', got "
                f"'{datatype}'."
            )

        # payload lines
        out: list[str] = [header]
        fmt = f"{{0:{self.bloc_num_format}}}".upper()

        # sentinel for "empty"
        sentinel = 1.0e32

        for i, v in enumerate(arr, start=1):
            if (key not in {"zrot", "trot", "freq"} and
                    (v == 0.0 or np.isnan(v) or np.isinf(v))):
                v = sentinel

            token = fmt.format(v)
            if i % self.block_size == 0:
                token += "\n"
            if i == arr.size:
                token += "\n"

            out.append(token)

        return out

    def interpolate_z(
        self,
        new_freq_array,
        interp_type: str = "slinear",
        bounds_error: bool = True,
        period_buffer: float | None = None,
    ) -> EMZ:
        """
        Interpolate the impedance tensor onto a new frequency
        grid and return a new Z object.
        """
        try:
            from scipy.interpolate import interp1d
        except Exception as exc:
            raise ImportError(
                "SciPy is required for interpolation. "
                "Install 'scipy' to use interpolate_z()."
            ) from exc

        if not isinstance(new_freq_array, np.ndarray):
            new_freq_array = np.array(new_freq_array, dtype=float)
        else:
            new_freq_array = new_freq_array.astype(float, copy=True)

        # keep previous behavior
        new_freq_array = np.around(new_freq_array, 2)

        if period_buffer is not None and 0.0 < period_buffer < 1.0:
            period_buffer += 1.0
            # gentle nudge to caller; avoid print
            import warnings as _warnings
            _warnings.warn(
                "period_buffer must be > 1; updated to "
                f"{period_buffer:.3g}.",
                RuntimeWarning,
            )

        old_freq = np.asarray(self.Z.freq, dtype=float)
        if old_freq.ndim != 1 or old_freq.size == 0:
            raise EdIDataError("Edi.Z.freq is empty or invalid.")

        if bounds_error:
            nfmin, nfmax = float(new_freq_array.min()), \
                           float(new_freq_array.max())
            ofmin, ofmax = float(old_freq.min()), float(old_freq.max())
            if nfmin < ofmin or nfmax > ofmax:
                raise ValueError(
                    "New frequency range must be within the old one; "
                    f"new[{nfmin:.5g}, {nfmax:.5g}] vs "
                    f"old[{ofmin:.5g}, {ofmax:.5g}]."
                )

        # build the target container
        n_new = int(new_freq_array.shape[0])
        new_Z = EMZ(
            z_array=np.zeros((n_new, 2, 2), dtype=complex),
            z_err_array=np.zeros((n_new, 2, 2), dtype=float),
            freq=new_freq_array,
        )

        # interpolate each tensor component independently on
        # the non-zero support of source values
        for i in range(2):
            for j in range(2):
                src = np.asarray(self.Z.z[:, i, j])
                src_err = np.asarray(self.Z.z_err[:, i, j])

                # non-zero complex entries (exclude pure zeros)
                nz = np.nonzero(src)[0]
                if nz.size == 0:
                    continue

                f = old_freq[nz]
                z_real = src[nz].real
                z_imag = src[nz].imag
                z_err = src_err[nz]

                # enforce interpolation only within the support
                inb = np.where(
                    (new_freq_array >= f.min()) &
                    (new_freq_array <= f.max())
                )[0]
                if inb.size == 0:
                    continue

                f_new = new_freq_array[inb]

                # optional period proximity filtering
                if isinstance(period_buffer, (int, float)):
                    keep_idx = []
                    for k, ff in enumerate(f_new):
                        # nearest original frequency (log10 metric)
                        dif = np.abs(np.log10(ff) - np.log10(f))
                        k0 = int(np.argmin(dif))
                        ratio = max(f[k0] / ff, ff / f[k0])
                        if ratio < float(period_buffer):
                            keep_idx.append(inb[k])
                    inb = np.array(keep_idx, dtype=int)
                    if inb.size == 0:
                        continue
                    f_new = new_freq_array[inb]

                # set up interpolants
                f_order = np.argsort(f)
                f_sorted = f[f_order]
                zr_sorted = z_real[f_order]
                zi_sorted = z_imag[f_order]
                ze_sorted = z_err[f_order]

                fr = interp1d(
                    f_sorted, zr_sorted, kind=interp_type,
                    bounds_error=True,
                )
                fi = interp1d(
                    f_sorted, zi_sorted, kind=interp_type,
                    bounds_error=True,
                )
                fe = interp1d(
                    f_sorted, ze_sorted, kind=interp_type,
                    bounds_error=True,
                )

                new_Z.z[inb, i, j] = fr(f_new) + 1j * fi(f_new)
                new_Z.z_err[inb, i, j] = fe(f_new)

        new_Z.compute_resistivity_phase()
        return new_Z


    def write_new_edi(
        self,
        new_edi_fn: Optional[str] = None,
        new_Z: Optional[EMZ] = None,
        new_Tipper=None,
        **kwargs,
    ) -> str:
        """
        Write a new EDI without mutating this instance.

        - If `new_Z` (and/or `new_Tipper`) is provided, those are
          used in the written file only; `self` is unchanged.
        - `new_edi_fn` is passed through to `write_edi` as the
          output filename (optional).
        """
        if not getattr(self, "edifile", None):
            raise EdIDataError(
                "No source EDI file bound to this instance; "
                "cannot rebuild a fresh writer."
            )

        # Rebuild a clean Edi from the same source so that any
        # manipulations here won't touch `self`.
        edi_obj = self.__class__(edi_filename=self.edifile)

        # Swap in new transfer functions if provided.
        edi_obj.Z = new_Z if new_Z is not None else self.Z
        if new_Tipper is not None:
            edi_obj.Tip = new_Tipper
        elif hasattr(self, "Tip"):
            edi_obj.Tip = self.Tip

        # Delegate to the main writer (v2 API).
        out_path = edi_obj.write_edi(
            new_edi_filename=new_edi_fn,
            **kwargs,
        )
        return out_path


    @property
    def station(self) -> str | None:
        return getattr(self.Head, "dataid", None)

    @station.setter
    def station(self, value: str) -> None:
        name = str(value)
        self.Head.dataid = name
        # keep section ID aligned with station name
        if hasattr(self, "MTEMAP"):
            self.MTEMAP.sectid = name

    @property
    def processingsoftware(self) -> str | None:
        ps = getattr(self.Info, "Processing", None)
        if ps is None:
            return None
        sw = getattr(ps, "ProcessingSoftware", None)
        return getattr(sw, "name", None)

    @processingsoftware.setter
    def processingsoftware(self, name: str) -> None:
        if not hasattr(self.Info, "Processing"):
            return
        ps = self.Info.Processing
        if not hasattr(ps, "ProcessingSoftware"):
            return
        ps.ProcessingSoftware.name = str(name)
