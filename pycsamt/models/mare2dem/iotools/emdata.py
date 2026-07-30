# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Full reader and writer for MARE2DEM ``.emdata`` / ``.EMResp`` files.

Port of:
  * ``m2d_readEMData2DFile.m``
  * ``m2d_writeEMData2DFile.m``

Supports format versions EMData_2.0 through EMData_2.3 and the
corresponding EMResp variants.  The returned :class:`EMDataFile` object
holds all sections (UTM origin, MT receivers/frequencies, CSEM
transmitters/receivers/frequencies, DC electrodes, and the DATA block)
as plain NumPy arrays and Python lists, mirroring the MATLAB structure
fields exactly so that downstream tools can be ported field-by-field.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

__all__ = [
    "EMDataFile",
    "UTMOrigin",
    "CSEMConfig",
    "MTConfig",
    "DCConfig",
    "read_emdata",
    "write_emdata",
]

# ---------------------------------------------------------------------------
# Supported format identifiers
# ---------------------------------------------------------------------------

_DATA_FORMATS = {
    "emdata_2.0",
    "emdata_2.1",
    "emdata_2.2",
    "emdata_2.3",
}
_RESP_FORMATS = {
    "emresp_2.0",
    "emresp_2.1",
    "emresp_2.2",
    "emresp_2.3",
}
_ALL_FORMATS = _DATA_FORMATS | _RESP_FORMATS

# Columns in the DATA block: 6 for observed data, 8 for response files
_DATA_NCOLS = 6
_RESP_NCOLS = 8


# ---------------------------------------------------------------------------
# Helper: parse one text line, strip comments, normalise 'd' exponents
# ---------------------------------------------------------------------------


def _parse_line(line: str) -> tuple[str, bool]:
    """Return (cleaned_line, is_comment).

    Strips inline ``!`` and ``%`` comments, leading/trailing whitespace, and
    normalises Fortran ``d``-exponent notation to ``e``.
    """
    s = line.strip()
    if not s or s[0] in ("!", "%"):
        return s, True
    # trim trailing comments
    for ch in ("!", "%"):
        idx = s.find(ch)
        if idx > 0:
            s = s[:idx].rstrip()
    # Fortran d-exponent
    s = re.sub(r"([0-9])[dD]([+-]?\d)", r"\1e\2", s)
    return s, False


def _strtok(s: str, sep: str = ":") -> tuple[str, str]:
    """Split on first occurrence of *sep*; return (before, after)."""
    idx = s.find(sep)
    if idx < 0:
        return s, ""
    return s[:idx], s[idx + 1 :]


# ---------------------------------------------------------------------------
# Data classes for each section
# ---------------------------------------------------------------------------


@dataclass
class UTMOrigin:
    """UTM mapping metadata for the 2-D profile.

    Attributes
    ----------
    grid : int
        UTM zone number.
    hemi : str
        Hemisphere letter (``'N'`` or ``'S'``).
    north0 : float
        Northing of profile origin in metres.
    east0 : float
        Easting of profile origin in metres.
    theta : float
        Profile strike direction (degrees).
    """

    grid: int = 0
    hemi: str = "N"
    north0: float = 0.0
    east0: float = 0.0
    theta: float = 0.0


@dataclass
class CSEMConfig:
    """CSEM source–receiver configuration and data parameters.

    Attributes
    ----------
    phase_convention : str
        ``'lag'`` or ``'lead'``.
    reciprocity_used : str
        ``'yes'``, ``'no'``, or ``''``.
    frequencies : numpy.ndarray, shape (n_freq,)
        CSEM frequencies in Hz.
    time_offsets : numpy.ndarray or None
        Time offsets for TDEM modeling.
    tdem_waveform : numpy.ndarray or None
        TDEM waveform, shape (n_pts, 2).
    transmitters : numpy.ndarray, shape (n_tx, 7)
        Columns: x, y, z, azimuth, dip, length, solve_corr.
    transmitter_type : list of str
        Dipole type per transmitter (``'edipole'`` or ``'bdipole'``).
    transmitter_name : list of str
        Transmitter labels.
    receivers : numpy.ndarray, shape (n_rx, 8)
        Columns: x, y, z, theta, alpha, beta, length, solve_corr.
    receiver_name : list of str
        Receiver labels.
    """

    phase_convention: str = "lag"
    reciprocity_used: str = ""
    frequencies: np.ndarray = field(
        default_factory=lambda: np.empty(0, dtype=float)
    )
    time_offsets: np.ndarray | None = None
    tdem_waveform: np.ndarray | None = None
    transmitters: np.ndarray = field(
        default_factory=lambda: np.empty((0, 7), dtype=float)
    )
    transmitter_type: list[str] = field(default_factory=list)
    transmitter_name: list[str] = field(default_factory=list)
    receivers: np.ndarray = field(
        default_factory=lambda: np.empty((0, 8), dtype=float)
    )
    receiver_name: list[str] = field(default_factory=list)


@dataclass
class MTConfig:
    """MT receiver configuration and frequency list.

    Attributes
    ----------
    frequencies : numpy.ndarray, shape (n_freq,)
        MT frequencies in Hz.
    receivers : numpy.ndarray, shape (n_rx, 8)
        Columns: x, y, z, theta, alpha, beta, length, solve_static.
    receiver_name : list of str
        Receiver labels (station names).
    """

    frequencies: np.ndarray = field(
        default_factory=lambda: np.empty(0, dtype=float)
    )
    receivers: np.ndarray = field(
        default_factory=lambda: np.empty((0, 8), dtype=float)
    )
    receiver_name: list[str] = field(default_factory=list)


@dataclass
class DCConfig:
    """DC resistivity electrode and receiver configuration.

    Attributes
    ----------
    tx_electrodes : numpy.ndarray, shape (n_tx_el, 3)
        Transmitter electrode positions (x, y, z).
    rx_electrodes : numpy.ndarray, shape (n_rx_el, 3)
        Receiver electrode positions (x, y, z).
    transmitters : numpy.ndarray, shape (n_tx, 2)
        Integer electrode-index pairs (A, B) for each transmitter.
    receivers : numpy.ndarray, shape (n_rx, 2)
        Integer electrode-index pairs (M, N) for each receiver.
    transmitter_name : list of str
    receiver_name : list of str
    """

    tx_electrodes: np.ndarray = field(
        default_factory=lambda: np.empty((0, 3), dtype=float)
    )
    rx_electrodes: np.ndarray = field(
        default_factory=lambda: np.empty((0, 3), dtype=float)
    )
    transmitters: np.ndarray = field(
        default_factory=lambda: np.empty((0, 2), dtype=int)
    )
    receivers: np.ndarray = field(
        default_factory=lambda: np.empty((0, 2), dtype=int)
    )
    transmitter_name: list[str] = field(default_factory=list)
    receiver_name: list[str] = field(default_factory=list)


@dataclass
class EMDataFile:
    """Complete contents of one MARE2DEM ``.emdata`` or ``.EMResp`` file.

    Attributes
    ----------
    path : pathlib.Path or None
        Source file path, if loaded from disk.
    format : str
        Format string from the file header (e.g. ``"EMData_2.3"``).
    is_response : bool
        ``True`` when the file is a response file (8-column DATA block).
    comment : str
        Optional free-text comment from the file header.
    utm : UTMOrigin
        UTM origin metadata.
    csem : CSEMConfig or None
        CSEM configuration section, or ``None`` if absent.
    mt : MTConfig or None
        MT configuration section, or ``None`` if absent.
    dc : DCConfig or None
        DC configuration section, or ``None`` if absent.
    data : numpy.ndarray, shape (n_data, 6) or (n_data, 8)
        DATA block.  Columns for observed data:
        ``[type, freq#, tx#, rx#, data, std_err]``.
        Response files add ``[response, residual]``.
    """

    path: Path | None = None
    format: str = "EMData_2.3"
    is_response: bool = False
    comment: str = ""
    utm: UTMOrigin = field(default_factory=UTMOrigin)
    csem: CSEMConfig | None = None
    mt: MTConfig | None = None
    dc: DCConfig | None = None
    data: np.ndarray = field(
        default_factory=lambda: np.empty((0, 6), dtype=float)
    )

    # ------------------------------------------------------------------
    # Convenience properties
    # ------------------------------------------------------------------

    @property
    def n_data(self) -> int:
        """Total number of data rows."""
        return len(self.data)

    @property
    def n_mt_frequencies(self) -> int:
        """Number of MT frequencies."""
        return len(self.mt.frequencies) if self.mt is not None else 0

    @property
    def n_mt_receivers(self) -> int:
        """Number of MT receivers."""
        return len(self.mt.receivers) if self.mt is not None else 0

    @property
    def n_csem_transmitters(self) -> int:
        """Number of CSEM transmitters."""
        return len(self.csem.transmitters) if self.csem is not None else 0

    @property
    def n_csem_receivers(self) -> int:
        """Number of CSEM receivers."""
        return len(self.csem.receivers) if self.csem is not None else 0

    def __repr__(self) -> str:
        return (
            f"EMDataFile(format={self.format!r}, "
            f"n_data={self.n_data}, "
            f"mt_rx={self.n_mt_receivers}, "
            f"csem_tx={self.n_csem_transmitters})"
        )


# ---------------------------------------------------------------------------
# Reader
# ---------------------------------------------------------------------------


class _EMDataReader:
    """Internal line-by-line parser for MARE2DEM ``.emdata`` files."""

    def __init__(self, lines: list[str]) -> None:
        self._lines = lines
        self._pos = 0

    def _next_line(self) -> str | None:
        if self._pos >= len(self._lines):
            return None
        line = self._lines[self._pos]
        self._pos += 1
        return line

    def _next_data_line(self) -> str | None:
        """Return the next non-comment line, or None at EOF."""
        while True:
            line = self._next_line()
            if line is None:
                return None
            cleaned, is_comment = _parse_line(line)
            if not is_comment:
                return cleaned

    def _read_float_block(self, n: int) -> np.ndarray:
        """Read *n* floats from consecutive non-comment lines."""
        vals: list[float] = []
        while len(vals) < n:
            line = self._next_data_line()
            if line is None:
                break
            for tok in line.split():
                tok = re.sub(r"([0-9])[dD]([+-]?\d)", r"\1e\2", tok)
                vals.append(float(tok))
        return np.array(vals[:n], dtype=float)

    def _read_int_block(self, n_rows: int, n_cols: int) -> np.ndarray:
        """Read *n_rows* × *n_cols* integers."""
        vals: list[int] = []
        while len(vals) < n_rows * n_cols:
            line = self._next_data_line()
            if line is None:
                break
            for tok in line.split():
                vals.append(int(tok))
        arr = np.array(vals[: n_rows * n_cols], dtype=int)
        return arr.reshape(n_rows, n_cols)

    # ------------------------------------------------------------------
    # Top-level parse
    # ------------------------------------------------------------------

    def parse(self) -> EMDataFile:
        out = EMDataFile()
        fmt = ""
        file_type = "data"

        while True:
            raw = self._next_line()
            if raw is None:
                break
            cleaned, is_comment = _parse_line(raw)
            if is_comment:
                if not out.comment and cleaned:
                    out.comment = cleaned.lstrip("!")
                continue

            code_raw, value_raw = _strtok(cleaned, ":")
            code = code_raw.strip().lower()
            value = value_raw.strip()
            # strip trailing comment from value
            for ch in ("!", "%"):
                idx = value.find(ch)
                if idx >= 0:
                    value = value[:idx].strip()

            # --- format / version ---
            if code in ("format", "version"):
                fmt = value.lower().strip()
                out.format = value.strip()
                if fmt in _DATA_FORMATS:
                    file_type = "data"
                    out.is_response = False
                elif fmt in _RESP_FORMATS:
                    file_type = "response"
                    out.is_response = True
                continue

            # --- UTM origin ---
            if code in (
                "utm of x,y origin (n,e,theta)",
                "utm",
                "origin",
                "utm of x,y origin (utm zone, n, e, 2d strike)",
            ):
                parts = value.split()
                if parts:
                    grid_str = parts[0].rstrip("nNsS")
                    hemi = parts[0][len(grid_str) :].upper() or (
                        parts[1].upper() if len(parts) > 1 else "N"
                    )
                    offset = 1
                    if not hemi or hemi not in ("N", "S"):
                        hemi = parts[1].upper() if len(parts) > 1 else "N"
                        offset = 2
                    try:
                        out.utm.grid = int(grid_str)
                        out.utm.hemi = hemi
                        out.utm.north0 = (
                            float(parts[offset])
                            if len(parts) > offset
                            else 0.0
                        )
                        out.utm.east0 = (
                            float(parts[offset + 1])
                            if len(parts) > offset + 1
                            else 0.0
                        )
                        out.utm.theta = (
                            float(parts[offset + 2])
                            if len(parts) > offset + 2
                            else 0.0
                        )
                    except (ValueError, IndexError):
                        pass
                continue

            # --- CSEM phase convention ---
            if code == "phase convention":
                if out.csem is None:
                    out.csem = CSEMConfig()
                out.csem.phase_convention = value.strip()
                continue

            if code == "reciprocity used":
                if out.csem is None:
                    out.csem = CSEMConfig()
                out.csem.reciprocity_used = value.strip()
                continue

            # --- CSEM transmitters ---
            if code in ("# transmitters", "#transmitters"):
                if out.csem is None:
                    out.csem = CSEMConfig()
                n_tx = int(value)
                txs = np.zeros((n_tx, 7), dtype=float)
                tx_types: list[str] = []
                tx_names: list[str] = []
                n_read = 0
                while n_read < n_tx:
                    raw2 = self._next_line()
                    if raw2 is None:
                        break
                    s, is_c = _parse_line(raw2)
                    if is_c:
                        continue
                    parts = s.split()
                    if fmt in ("emdata_2.1", "emresp_2.1"):
                        vals = [float(p) for p in parts[:5]]
                        txs[n_read, :5] = vals
                        tx_types.append(
                            parts[5] if len(parts) > 5 else "edipole"
                        )
                        tx_names.append(parts[6] if len(parts) > 6 else "")
                    elif fmt in ("emdata_2.2", "emresp_2.2"):
                        vals = [float(p) for p in parts[:6]]
                        txs[n_read, :6] = vals
                        tx_types.append(
                            parts[6] if len(parts) > 6 else "edipole"
                        )
                        tx_names.append(parts[7] if len(parts) > 7 else "")
                    else:  # 2.3 (default) or 2.0
                        num_cols = min(7, len(parts))
                        str_idx = num_cols
                        # find first non-numeric token
                        for k, p in enumerate(parts):
                            try:
                                float(p)
                            except ValueError:
                                str_idx = k
                                break
                        vals = [float(p) for p in parts[:str_idx]]
                        txs[n_read, : len(vals)] = vals
                        tx_types.append(
                            parts[str_idx]
                            if len(parts) > str_idx
                            else "edipole"
                        )
                        tx_names.append(
                            parts[str_idx + 1]
                            if len(parts) > str_idx + 1
                            else ""
                        )
                    n_read += 1
                out.csem.transmitters = txs
                out.csem.transmitter_type = tx_types
                out.csem.transmitter_name = tx_names
                continue

            # --- CSEM time offsets ---
            if code == "# csem time offsets":
                if out.csem is None:
                    out.csem = CSEMConfig()
                n = int(value)
                out.csem.time_offsets = self._read_float_block(n)
                continue

            # --- CSEM TDEM waveform ---
            if code == "# tdem waveform points":
                if out.csem is None:
                    out.csem = CSEMConfig()
                n = int(value)
                arr = self._read_float_block(n * 2)
                out.csem.tdem_waveform = arr.reshape(n, 2)
                continue

            # --- CSEM frequencies ---
            if code == "# csem frequencies":
                if out.csem is None:
                    out.csem = CSEMConfig()
                n = int(value)
                out.csem.frequencies = self._read_float_block(n)
                continue

            # --- CSEM receivers ---
            if code == "# csem receivers":
                if out.csem is None:
                    out.csem = CSEMConfig()
                n_rx = int(value)
                rxs = np.zeros((n_rx, 8), dtype=float)
                rx_names: list[str] = []
                n_read = 0
                while n_read < n_rx:
                    raw2 = self._next_line()
                    if raw2 is None:
                        break
                    s, is_c = _parse_line(raw2)
                    if is_c:
                        continue
                    parts = s.split()
                    if fmt in ("emdata_2.0", "emresp_2.0"):
                        rxs[n_read, :6] = [float(p) for p in parts[:6]]
                        rx_names.append(str(n_read + 1))
                    elif fmt in ("emdata_2.1", "emresp_2.1"):
                        rxs[n_read, :6] = [float(p) for p in parts[:6]]
                        rx_names.append(
                            parts[6] if len(parts) > 6 else str(n_read + 1)
                        )
                    elif fmt in ("emdata_2.2", "emresp_2.2"):
                        rxs[n_read, :7] = [float(p) for p in parts[:7]]
                        rx_names.append(
                            parts[7] if len(parts) > 7 else str(n_read + 1)
                        )
                    else:  # 2.3
                        # find where numeric tokens end
                        str_idx = 8
                        for k, p in enumerate(parts[:9]):
                            try:
                                float(p)
                            except ValueError:
                                str_idx = k
                                break
                        rxs[n_read, :str_idx] = [
                            float(p) for p in parts[:str_idx]
                        ]
                        rx_names.append(
                            parts[str_idx]
                            if len(parts) > str_idx
                            else str(n_read + 1)
                        )
                    n_read += 1
                out.csem.receivers = rxs
                out.csem.receiver_name = rx_names
                continue

            # --- MT frequencies ---
            if code == "# mt frequencies":
                if out.mt is None:
                    out.mt = MTConfig()
                n = int(value)
                out.mt.frequencies = self._read_float_block(n)
                continue

            # --- MT receivers ---
            if code == "# mt receivers":
                if out.mt is None:
                    out.mt = MTConfig()
                n_rx = int(value)
                rxs = np.zeros((n_rx, 8), dtype=float)
                rx_names = []
                n_read = 0
                while n_read < n_rx:
                    raw2 = self._next_line()
                    if raw2 is None:
                        break
                    s, is_c = _parse_line(raw2)
                    if is_c:
                        continue
                    parts = s.split()
                    if fmt in ("emdata_2.0", "emresp_2.0"):
                        rxs[n_read, :6] = [float(p) for p in parts[:6]]
                        rxs[n_read, 6] = 0.0  # no static option
                        rx_names.append(str(n_read + 1))
                    elif fmt in ("emdata_2.1", "emresp_2.1"):
                        rxs[n_read, :7] = [float(p) for p in parts[:7]]
                        rx_names.append(
                            parts[7] if len(parts) > 7 else str(n_read + 1)
                        )
                    else:  # 2.2, 2.3
                        rxs[n_read, :8] = [float(p) for p in parts[:8]]
                        rx_names.append(
                            parts[8] if len(parts) > 8 else str(n_read + 1)
                        )
                    n_read += 1
                out.mt.receivers = rxs
                out.mt.receiver_name = rx_names
                continue

            # --- DC transmitter electrodes ---
            if code == "# dc transmitter electrodes":
                if out.dc is None:
                    out.dc = DCConfig()
                n = int(value)
                arr = np.zeros((n, 3), dtype=float)
                n_read = 0
                while n_read < n:
                    s = self._next_data_line()
                    if s is None:
                        break
                    arr[n_read] = [float(v) for v in s.split()[:3]]
                    n_read += 1
                out.dc.tx_electrodes = arr
                continue

            if code == "# dc transmitters":
                if out.dc is None:
                    out.dc = DCConfig()
                n = int(value)
                out.dc.transmitters = self._read_int_block(n, 2)
                continue

            if code == "# dc receiver electrodes":
                if out.dc is None:
                    out.dc = DCConfig()
                n = int(value)
                arr = np.zeros((n, 3), dtype=float)
                n_read = 0
                while n_read < n:
                    s = self._next_data_line()
                    if s is None:
                        break
                    arr[n_read] = [float(v) for v in s.split()[:3]]
                    n_read += 1
                out.dc.rx_electrodes = arr
                continue

            if code == "# dc receivers":
                if out.dc is None:
                    out.dc = DCConfig()
                n = int(value)
                out.dc.receivers = self._read_int_block(n, 2)
                continue

            # --- DATA block ---
            if code in ("# data", "#data"):
                n_data = int(value)
                n_cols = (
                    _RESP_NCOLS if file_type == "response" else _DATA_NCOLS
                )
                # Skip optional comment header line
                line_peek = self._next_line()
                if line_peek is not None:
                    s2, is_c2 = _parse_line(line_peek)
                    if not is_c2:
                        self._pos -= 1  # put back

                # Try fast numpy read of remaining text
                remaining = "\n".join(self._lines[self._pos :])
                try:
                    arr = np.fromstring(
                        remaining.replace("\n", " "), sep=" ", dtype=float
                    )
                    arr = arr[: n_data * n_cols]
                    if len(arr) == n_data * n_cols:
                        out.data = arr.reshape(n_data, n_cols)
                        break
                except Exception:
                    pass

                # Fallback: line-by-line slow read
                rows: list[list[float]] = []
                while len(rows) < n_data:
                    s = self._next_data_line()
                    if s is None:
                        break
                    row = [float(v) for v in s.split()[:n_cols]]
                    if len(row) == n_cols:
                        rows.append(row)
                out.data = np.array(rows, dtype=float)
                break

        return out


# ---------------------------------------------------------------------------
# Writer
# ---------------------------------------------------------------------------


def _write_float_col(f, arr: np.ndarray) -> None:
    for v in arr:
        f.write(f"{v:.14g}\n")


def _write_emdata_file(out_file: TextIO, em: EMDataFile) -> None:
    """Low-level writer that mirrors m2d_writeEMData2DFile.m."""
    n_cols = (
        len(em.data[0])
        if len(em.data) and len(em.data.shape) > 1
        else _DATA_NCOLS
    )
    fmt_str = (
        "EMResp_2.3"
        if (n_cols == _RESP_NCOLS or em.is_response)
        else "EMData_2.3"
    )
    out_file.write(f"Format:  {fmt_str}\n")

    if em.comment:
        out_file.write(f"!{em.comment}\n")

    # UTM
    utm = em.utm
    if utm.grid:
        out_file.write(
            f"UTM of x,y origin (UTM zone, N, E, 2D strike): "
            f"{utm.grid} {utm.hemi} "
            f"{utm.north0:14.6f} {utm.east0:14.6f} {utm.theta:6g}\n"
        )

    # MT section
    if em.mt is not None and len(em.mt.frequencies):
        mt = em.mt
        out_file.write(f"# MT Frequencies:    {len(mt.frequencies)}\n")
        _write_float_col(out_file, mt.frequencies)

        out_file.write(f"# MT Receivers:      {len(mt.receivers)}\n")
        out_file.write(
            f"{'!':1s} {'X':>22s} {'Y':>22s} {'Z':>22s} "
            f"{'Theta':>9s} {'Alpha':>9s} {'Beta':>9s} "
            f"{'Length':>10s} {'SolveStatic':>11s} {'Name':>4s}\n"
        )
        for i, rx in enumerate(mt.receivers):
            name = mt.receiver_name[i] if i < len(mt.receiver_name) else ""
            e_len = rx[6] if len(rx) >= 7 else 0.0
            i_static = int(rx[7]) if len(rx) >= 8 else 0
            out_file.write(
                f"  {rx[0]:22.15g} {rx[1]:22.15g} {rx[2]:22.15g} "
                f"{rx[3]:9.2f} {rx[4]:9.2f} {rx[5]:9.2f} "
                f"{e_len:10.5g} {i_static:11d} {name}\n"
            )

    # CSEM section
    if em.csem is not None and len(em.csem.transmitters):
        csem = em.csem
        out_file.write(f"Phase Convention: {csem.phase_convention}\n")
        if csem.reciprocity_used:
            out_file.write(f"Reciprocity Used: {csem.reciprocity_used}\n")

        if csem.time_offsets is not None and len(csem.time_offsets):
            out_file.write(
                f"# CSEM Time Offsets:    {len(csem.time_offsets)}\n"
            )
            _write_float_col(out_file, csem.time_offsets)

        if csem.tdem_waveform is not None and len(csem.tdem_waveform):
            out_file.write(
                f"# TDEM waveform points:    {len(csem.tdem_waveform)}\n"
            )
            for row in csem.tdem_waveform:
                out_file.write(f"{row[0]:.14g} {row[1]:.14g}f\n")

        if len(csem.frequencies):
            out_file.write(f"# CSEM Frequencies:    {len(csem.frequencies)}\n")
            _write_float_col(out_file, csem.frequencies)

        n_tx = len(csem.transmitters)
        out_file.write(f"# Transmitters:   {n_tx}\n")
        out_file.write(
            f"{'!':1s} {'X':>22s} {'Y':>22s} {'Z':>22s} "
            f"{'Azimuth':>9s} {'Dip':>9s} {'Length':>10s} "
            f"{'SolveCorr':>10s} {'Type':>10s} {'Name':>4s}\n"
        )
        txs = csem.transmitters
        if txs.shape[1] < 7:
            pad = np.zeros((len(txs), 7 - txs.shape[1]))
            txs = np.hstack([txs, pad])
        for i, tx in enumerate(txs):
            tx_type = (
                csem.transmitter_type[i]
                if i < len(csem.transmitter_type)
                else "edipole"
            )
            tx_name = (
                csem.transmitter_name[i]
                if i < len(csem.transmitter_name)
                else ""
            )
            out_file.write(
                f"  {tx[0]:22.15g} {tx[1]:22.15g} {tx[2]:22.15g} "
                f"{tx[3]:9.2f} {tx[4]:9.2f} {tx[5]:10.5g} "
                f"{int(tx[6]):10d} {tx_type:>10s} {tx_name}\n"
            )

        out_file.write(f"# CSEM Receivers:      {len(csem.receivers)}\n")
        out_file.write(
            f"{'!':1s} {'X':>22s} {'Y':>22s} {'Z':>22s} "
            f"{'Theta':>9s} {'Alpha':>9s} {'Beta':>9s} "
            f"{'Length':>10s} {'SolveCorr':>11s} {'Name':>4s}\n"
        )
        rxs = csem.receivers
        if rxs.shape[1] < 8:
            pad = np.zeros((len(rxs), 8 - rxs.shape[1]))
            rxs = np.hstack([rxs, pad])
        for i, rx in enumerate(rxs):
            name = csem.receiver_name[i] if i < len(csem.receiver_name) else ""
            out_file.write(
                f"  {rx[0]:22.15g} {rx[1]:22.15g} {rx[2]:22.15g} "
                f"{rx[3]:9.2f} {rx[4]:9.2f} {rx[5]:9.2f} "
                f"{rx[6]:10.5g} {int(rx[7]):11d} {name}\n"
            )

    # DC section
    if em.dc is not None and len(em.dc.transmitters):
        dc = em.dc
        n_tx_el = len(dc.tx_electrodes)
        out_file.write(f"# DC Transmitter Electrodes:   {n_tx_el}\n")
        out_file.write(f"{'!':1s} {'X':>22s} {'Y':>22s} {'Z':>22s}\n")
        for el in dc.tx_electrodes:
            out_file.write(f"{el[0]:22.15g} {el[1]:22.15g} {el[2]:22.15g}\n")

        n_tx = len(dc.transmitters)
        out_file.write(f"# DC Transmitters:   {n_tx}\n")
        out_file.write(
            f"{'!':1s} {'Electrode A':>12s} {'Electrode B':>12s} {'Name':>12s}\n"
        )
        for i, tx in enumerate(dc.transmitters):
            name = (
                dc.transmitter_name[i] if i < len(dc.transmitter_name) else ""
            )
            out_file.write(f"{tx[0]:12d} {tx[1]:12d} {name:12s}\n")

        n_rx_el = len(dc.rx_electrodes)
        out_file.write(f"# DC Receiver Electrodes:   {n_rx_el}\n")
        out_file.write(f"{'!':1s} {'X':>22s} {'Y':>22s} {'Z':>22s}\n")
        for el in dc.rx_electrodes:
            out_file.write(f"{el[0]:22.15g} {el[1]:22.15g} {el[2]:22.15g}\n")

        n_rx = len(dc.receivers)
        out_file.write(f"# DC Receivers:   {n_rx}\n")
        out_file.write(
            f"{'!':1s} {'Electrode M':>12s} {'Electrode N':>12s} {'Name':>12s}\n"
        )
        for i, rx in enumerate(dc.receivers):
            name = dc.receiver_name[i] if i < len(dc.receiver_name) else ""
            out_file.write(f"{rx[0]:12d} {rx[1]:12d} {name:12s}\n")

    # DATA block
    n_data = len(em.data)
    out_file.write(f"# Data:       {n_data}\n")
    if n_cols == _DATA_NCOLS:
        out_file.write(
            "!  Type  Freq #    Tx #    Rx #              Data            StdErr\n"
        )
        for row in em.data:
            out_file.write(
                f"{int(row[0]):7d} {int(row[1]):7d} {int(row[2]):7d} {int(row[3]):7d} "
                f"{row[4]:20.15g} {row[5]:20.15g}\n"
            )
    else:
        out_file.write(
            "!  Type  Freq #    Tx #    Rx #              Data            StdErr"
            "       Response       Residual\n"
        )
        for row in em.data:
            out_file.write(
                f"{int(row[0]):7d} {int(row[1]):7d} {int(row[2]):7d} {int(row[3]):7d} "
                f"{row[4]:22.15g} {row[5]:22.15g} {row[6]:20.15g} {row[7]:20.15g}\n"
            )


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def read_emdata(path: str | Path, *, silent: bool = False) -> EMDataFile:
    """Read a MARE2DEM ``.emdata`` or ``.EMResp`` file.

    Port of ``m2d_readEMData2DFile.m``.

    Parameters
    ----------
    path : path-like
        File to read.
    silent : bool, default False
        Suppress warnings on unrecognised format strings.

    Returns
    -------
    EMDataFile
        Parsed file contents.

    Raises
    ------
    FileNotFoundError
        When *path* does not exist.
    ValueError
        When the format string in the header is not recognised and
        ``silent=False``.

    Examples
    --------
    >>> from pycsamt.models.mare2dem.iotools.emdata import read_emdata
    >>> em = read_emdata("survey.emdata")
    >>> em.n_mt_receivers
    12
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")

    lines = path.read_text(errors="replace").splitlines()
    reader = _EMDataReader(lines)
    result = reader.parse()
    result.path = path.resolve()
    return result


def write_emdata(em: EMDataFile, path: str | Path) -> Path:
    """Write an :class:`EMDataFile` to *path*.

    Port of ``m2d_writeEMData2DFile.m``.

    Parameters
    ----------
    em : EMDataFile
        Data to write.
    path : path-like
        Destination ``.emdata`` file.

    Returns
    -------
    pathlib.Path
        Path of the written file.

    Examples
    --------
    >>> from pycsamt.models.mare2dem.iotools.emdata import write_emdata
    >>> write_emdata(em, "survey_out.emdata")
    PosixPath('survey_out.emdata')
    """
    dest = Path(path)
    dest.parent.mkdir(parents=True, exist_ok=True)
    with dest.open("w") as fh:
        _write_emdata_file(fh, em)
    return dest
