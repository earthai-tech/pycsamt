# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
stratagem.io
============

Low-level I/O for Stratagem hardware surveys:

* :class:`StratagemRawReader` — parse raw 19-column ASCII spectral files
  (.HXX / .MDX) for quality diagnostics (SNR mask, stack counts, sensor
  metadata).  Does **not** produce EDI output; use WinGLink for that.

* :class:`EDIBatch` — load a directory of WinGLink-exported EDI files as a
  list of :class:`~pycsamt.seg.edi.EDIFile` objects with natural-sort ordering
  that respects Stratagem's three-digit station numbering (001 … 087).

Raw file format
---------------
Each component file (X*.NNN, Y*.NNN, Z*.NNN) is a whitespace-separated
ASCII table whose columns are::

    col 0   frequency bin (Hz)
    col 1   instrument constant (~2.93 for Stratagem)
    col 2   stack count (0 = measurement absent / below threshold)
    col 3-18  16 cross-spectral values (Re/Im of impedance components)

Rows where col 2 == 0 carry no useful signal and are marked False in the
:attr:`~StratagemRawReader.snr_mask_` array.
"""

from __future__ import annotations

import re
import warnings
from pathlib import Path

# Stratagem hardware files store values in scientific notation without spaces
# between adjacent numbers when the second starts with '-' (e.g.
# "2.152e+002-5.546e-001").  We use a regex to extract individual numbers.
_NUM_RE = re.compile(r"[+-]?(?:\d+\.\d+|\d+)(?:[eE][+-]?\d+)?")

# Stratagem files end each line with \r\r\n (double-CR + LF).  We read binary
# and normalise \r\n → \n then lone \r → \n before splitting, so the parser
# handles real hardware files (\r\r\n) and plain \n test fixtures equally.


import numpy as np
import pandas as pd

from ..api.property import PyCSAMTObject
from ..exceptions import FileHandlingError, NotFittedError
from ..seg.edi import EDIFile

__all__ = ["StratagemRawReader", "EDIBatch"]


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def _station_number(stem: str) -> int:
    """Extract the trailing numeric suffix from a file stem (e.g. 'X2HX.001' → 1)."""
    m = re.search(r"(\d+)$", stem)
    return int(m.group(1)) if m else 0


def _edi_sort_key(path: Path) -> tuple:
    """Natural-sort key: splits the stem into (text, int, text, …) tuples."""
    parts = re.split(r"(\d+)", path.stem.lower())
    return tuple(int(p) if p.isdigit() else p for p in parts)


def _read_19col(path: Path) -> np.ndarray:
    """Parse one Stratagem ASCII file into a float array of shape (n_rows, 19).

    Two format quirks are handled:

    * **CRLF + trailing CR**: physical lines end with ``\\r\\r\\n``.  Opening
      in text mode on Linux causes the first ``\\r`` to be treated as a
      standalone line-end, splitting every data record into two.  We read
      binary and split on ``\\r\\n`` instead.
    * **Adjacent negative numbers**: valid-measurement rows omit the space
      before a negative value (e.g. ``2.152e+002-5.546e-001``).  A regex
      extractor is used instead of ``str.split()``.

    Rows with fewer than 3 extracted numbers are skipped.  Rows with fewer
    than 19 numbers are zero-padded on the right.
    """
    rows: list[list[float]] = []
    with open(path, "rb") as fh:
        raw = fh.read()
    # Normalise \r\n → \n then lone \r → \n so both hardware files (\r\r\n)
    # and plain-\n test fixtures are handled identically.
    content = raw.replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    for line_bytes in content.split(b"\n"):
        line = line_bytes.decode("ascii", errors="replace").strip()
        if not line:
            continue
        vals = [float(m) for m in _NUM_RE.findall(line)]
        if len(vals) < 3:
            continue
        while len(vals) < 19:
            vals.append(0.0)
        rows.append(vals[:19])
    return np.array(rows, dtype=np.float64) if rows else np.zeros((0, 19))


def _parse_sensors_tbl(path: Path) -> dict[str, str]:
    """Parse SENSORS.TBL into a mapping of lower-case name → original name."""
    sensors: dict[str, str] = {}
    if not path.is_file():
        return sensors
    with open(path, errors="replace") as fh:
        for raw in fh:
            name = raw.strip()
            if name:
                sensors[name.lower()] = name
    return sensors


def _build_masks(
    files: list[Path],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (freq_grid, snr_mask, stack_counts) from a list of component files.

    *freq_grid*   : 1-D float array of length n_freqs
    *snr_mask*    : bool array (n_stations, n_freqs) — True where col-2 > 0
    *stack_counts*: int  array (n_stations, n_freqs) — raw stack count
    """
    freq_grid: np.ndarray | None = None
    snr_rows: list[np.ndarray | None] = []
    stack_rows: list[np.ndarray | None] = []

    for fp in files:
        try:
            mat = _read_19col(fp)
        except OSError:
            snr_rows.append(None)
            stack_rows.append(None)
            continue

        if mat.shape[0] == 0:
            snr_rows.append(None)
            stack_rows.append(None)
            continue

        freqs = mat[:, 0]
        stacks = mat[:, 2].astype(np.int32)

        if freq_grid is None:
            freq_grid = freqs
        snr_rows.append(stacks > 0)
        stack_rows.append(stacks)

    if freq_grid is None:
        freq_grid = np.array([], dtype=np.float64)

    n_freq = freq_grid.size
    n_sta = len(files)
    snr_mask = np.zeros((n_sta, n_freq), dtype=bool)
    stack_counts = np.zeros((n_sta, n_freq), dtype=np.int32)

    for i, (sm, sc) in enumerate(zip(snr_rows, stack_rows)):
        if sm is not None and sm.size == n_freq:
            snr_mask[i] = sm
            stack_counts[i] = sc

    return freq_grid, snr_mask, stack_counts


# ---------------------------------------------------------------------------
# StratagemRawReader
# ---------------------------------------------------------------------------


class StratagemRawReader(PyCSAMTObject):
    """Parse raw Stratagem hardware files (19-column ASCII) for QC diagnostics.

    The Stratagem AMT system records frequency-band cross-spectral data for
    each station in separate component files named ``X*.NNN``, ``Y*.NNN``, and
    ``Z*.NNN`` (where NNN is the zero-padded station number).  This reader
    extracts the frequency grid, stack counts, and SNR masks from those files.

    .. important::
       This class does **not** compute impedance tensors or produce EDI
       output.  Use WinGLink 1.0.x for raw → EDI conversion.  The masks
       produced here can later be consumed by
       :class:`~pycsamt.stratagem.qc.FrequencyFilter` to apply
       hardware-level quality information to the WinGLink EDIs.

    Parameters
    ----------
    station_dir : path-like, optional
        Directory containing the raw Stratagem component files.  May also
        be supplied to :meth:`fit`.
    component : {'X', 'Y', 'Z', 'ALL'}, default 'X'
        Which component family to read for QC masks.  ``'ALL'`` reads X, Y,
        and Z and stores per-component masks under :attr:`component_masks_`.
    verbose : int, default 0
        Verbosity level.  0 = silent; ≥1 = progress messages.

    Attributes
    ----------
    stations_ : list of str
        Station file names in acquisition order (e.g. ``['X2HX.001', …]``).
    station_numbers_ : ndarray of int, shape (n_stations,)
        Hardware station numbers extracted from the extension
        (e.g. ``[1, 2, …, 87]``).
    freqs_ : ndarray of shape (n_freqs,)
        Frequency grid in Hz, read from the X-component of the first station.
    snr_mask_ : ndarray of shape (n_stations, n_freqs), dtype bool
        True where the stack count is non-zero (measurement present).
    stack_counts_ : ndarray of shape (n_stations, n_freqs), dtype int
        Raw stack counts as recorded by the hardware.
    n_stations_ : int
    n_freqs_ : int
    sensors_ : dict
        Contents of ``SENSORS.TBL`` as ``{lower_name: original_name}``.
    component_masks_ : dict, optional
        Present only when ``component='ALL'``.  Maps ``'X'``, ``'Y'``,
        ``'Z'`` to their respective ``(snr_mask, stack_counts)`` tuples.

    Examples
    --------
    >>> rdr = StratagemRawReader("原始数据/2HX").fit()
    >>> rdr.snr_mask_.shape   # (n_stations, n_freqs)
    (87, 292)
    >>> good = rdr.snr_mask_.sum(axis=1)   # usable freqs per station
    """

    __repr_fields__ = ("station_dir", "component", "n_stations_", "n_freqs_")

    def __init__(
        self,
        station_dir: str | Path | None = None,
        *,
        component: str = "X",
        verbose: int = 0,
    ) -> None:
        self.station_dir = station_dir
        self.component = component.upper()
        self.verbose = verbose

    # ------------------------------------------------------------------
    def fit(
        self,
        station_dir: str | Path | None = None,
    ) -> StratagemRawReader:
        """Read raw Stratagem files and build QC arrays.

        Parameters
        ----------
        station_dir : path-like, optional
            Override the directory set in ``__init__``.

        Returns
        -------
        self
        """
        if station_dir is not None:
            self.station_dir = station_dir
        if self.station_dir is None:
            raise ValueError("station_dir must be provided")

        d = Path(self.station_dir).expanduser().resolve()
        if not d.is_dir():
            raise FileHandlingError(f"Directory not found: {d}")

        components = (
            ["X", "Y", "Z"] if self.component == "ALL" else [self.component]
        )

        # collect station files for the primary component (always X-first)
        # Note: Stratagem uses the extension as the station number (X2HX.001,
        # X2HX.002, …), so we sort/filter on p.name, not p.stem.
        primary = components[0]
        x_files = sorted(
            d.glob(f"{primary}*.*"),
            key=lambda p: _station_number(p.name),
        )
        # keep only files whose content looks numeric (skip calibration files)
        x_files = [f for f in x_files if _station_number(f.name) > 0]

        if not x_files:
            raise FileHandlingError(
                f"No {primary}-component station files found in {d}.\n"
                "Expected files named like 'X2HX.001', 'X2HX.002', …"
            )

        # Store file names (e.g. 'X2HX.001') so that station numbers can be
        # extracted from the extension; stems (e.g. 'X2HX') are all identical
        # across stations and therefore useless as identifiers.
        self._station_files_: list[Path] = x_files
        self.stations_: list[str] = [p.name for p in x_files]
        self.station_numbers_: np.ndarray = np.array(
            [_station_number(p.name) for p in x_files], dtype=np.int32
        )

        freq_grid, snr_mask, stack_counts = _build_masks(x_files)
        self.freqs_ = freq_grid
        self.snr_mask_ = snr_mask
        self.stack_counts_ = stack_counts
        self.n_stations_ = len(x_files)
        self.n_freqs_ = freq_grid.size

        if self.component == "ALL":
            self.component_masks_: dict[str, tuple] = {}
            for comp in components:
                comp_files = sorted(
                    d.glob(f"{comp}*.*"),
                    key=lambda p: _station_number(p.name),
                )
                comp_files = [
                    f for f in comp_files if _station_number(f.name) > 0
                ]
                if comp_files:
                    _, cm, cs = _build_masks(comp_files)
                    self.component_masks_[comp] = (cm, cs)

        self.sensors_ = _parse_sensors_tbl(d / "SENSORS.TBL")

        if self.verbose:
            print(
                f"[StratagemRawReader] {self.n_stations_} stations, "
                f"{self.n_freqs_} frequencies loaded from {d.name}"
            )

        return self

    # ------------------------------------------------------------------
    def usable_freq_counts(self) -> np.ndarray:
        """Return the number of usable frequencies per station.

        Returns
        -------
        ndarray of shape (n_stations,)
            ``snr_mask_.sum(axis=1)``
        """
        if not hasattr(self, "snr_mask_"):
            raise NotFittedError("Call fit() first.")
        return self.snr_mask_.sum(axis=1)

    def station_coverage(self) -> float:
        """Fraction of (station, frequency) cells with valid data.

        Returns
        -------
        float in [0, 1]
        """
        if not hasattr(self, "snr_mask_"):
            raise NotFittedError("Call fit() first.")
        return float(self.snr_mask_.mean())

    # ------------------------------------------------------------------
    def match_to_edis(self, edi_objects: list) -> dict[int, int]:
        """Map EDI batch indices to raw station indices by hardware number.

        Stratagem raw files are numbered ``X*.001`` … ``X*.087`` (extension
        = hardware station number).  WinGLink EDI files are named
        ``Z*002.edi`` … ``Z*087.edi`` (stem suffix = same number).
        Index-based alignment is wrong when WinGLink skipped stations or
        the two sequences start at different offsets.

        This method extracts the numeric suffix from each EDI file's stem
        and from each raw file's name, then cross-references by value.

        Parameters
        ----------
        edi_objects : list of EDIFile

        Returns
        -------
        dict[int, int]
            ``{edi_batch_index: raw_station_index}`` for every EDI that
            has a matching raw file.  EDIs with no matching raw file are
            absent from the dict.

        Examples
        --------
        >>> mapping = rdr.match_to_edis(batch.edi_objects_)
        >>> mapping[0]   # raw index for the first EDI station
        1
        """
        if not hasattr(self, "station_numbers_"):
            raise NotFittedError("Call fit() first.")

        raw_num_to_idx: dict[int, int] = {
            int(n): i for i, n in enumerate(self.station_numbers_)
        }
        result: dict[int, int] = {}
        for j, edi in enumerate(edi_objects):
            path = getattr(edi, "path", None)
            if path is not None:
                num = _station_number(Path(path).stem)
            else:
                # fall back to DATAID numeric suffix
                dataid = getattr(edi, "station", None) or ""
                num = _station_number(str(dataid))
            if num and num in raw_num_to_idx:
                result[j] = raw_num_to_idx[num]
        return result

    # ------------------------------------------------------------------
    def station_frame(self) -> pd.DataFrame:
        """Per-station coverage summary as a DataFrame.

        Returns
        -------
        pandas.DataFrame
            One row per station.  Columns: ``station``, ``station_number``,
            ``total_freqs``, ``usable_freqs``, ``coverage``,
            ``max_stacks``, ``med_stacks``.

        Examples
        --------
        >>> rdr.station_frame().sort_values("coverage").head()
        """
        if not hasattr(self, "snr_mask_"):
            raise NotFittedError("Call fit() first.")

        usable = self.snr_mask_.sum(axis=1)
        valid_stacks = np.where(
            self.stack_counts_ > 0, self.stack_counts_, np.nan
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            med_stacks = np.nanmedian(valid_stacks, axis=1)

        return pd.DataFrame(
            {
                "station": self.stations_,
                "station_number": self.station_numbers_.tolist(),
                "total_freqs": self.n_freqs_,
                "usable_freqs": usable.tolist(),
                "coverage": (usable / max(1, self.n_freqs_)).tolist(),
                "max_stacks": self.stack_counts_.max(axis=1).tolist(),
                "med_stacks": med_stacks.tolist(),
            }
        )

    # ------------------------------------------------------------------
    def freq_frame(self) -> pd.DataFrame:
        """Per-frequency coverage summary as a DataFrame.

        Returns
        -------
        pandas.DataFrame
            One row per frequency bin.  Columns: ``freq_hz``,
            ``stations_ok``, ``frac_ok``, ``med_stacks``.

        Examples
        --------
        >>> rdr.freq_frame().query("frac_ok > 0.8")
        """
        if not hasattr(self, "snr_mask_"):
            raise NotFittedError("Call fit() first.")

        stations_ok = self.snr_mask_.sum(axis=0)
        valid_stacks = np.where(
            self.stack_counts_ > 0, self.stack_counts_, np.nan
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            med_stacks = np.nanmedian(valid_stacks, axis=0)

        return pd.DataFrame(
            {
                "freq_hz": self.freqs_.tolist(),
                "stations_ok": stations_ok.tolist(),
                "frac_ok": (stations_ok / max(1, self.n_stations_)).tolist(),
                "med_stacks": med_stacks.tolist(),
            }
        )

    # ------------------------------------------------------------------
    def stack_audit(self) -> pd.DataFrame:
        """Full stack-count audit as a (stations × frequencies) DataFrame.

        Returns
        -------
        pandas.DataFrame
            Index = station file names, columns = frequency values (Hz),
            values = integer stack counts (0 = no measurement).

        Examples
        --------
        >>> audit = rdr.stack_audit()
        >>> audit.loc["X2HX.005"]          # one station across all freqs
        >>> (audit > 0).sum(axis=1).plot()  # usable freqs per station
        """
        if not hasattr(self, "stack_counts_"):
            raise NotFittedError("Call fit() first.")

        cols = np.round(self.freqs_, 6)
        return pd.DataFrame(
            self.stack_counts_,
            index=self.stations_,
            columns=cols,
        )

    # ------------------------------------------------------------------
    def plot_coverage(
        self,
        *,
        kind: str = "snr",
        cmap: str = "RdYlGn",
        figsize: tuple | None = None,
        log_freq: bool = True,
        title: str | None = None,
    ):
        """Plot hardware data coverage as a station × frequency heatmap.

        Parameters
        ----------
        kind : {'snr', 'stacks'}, default ``'snr'``
            ``'snr'`` plots the boolean presence/absence mask;
            ``'stacks'`` shows the raw stack count values.
        cmap : str, default ``'RdYlGn'``
            Matplotlib colormap.
        figsize : tuple, optional
        log_freq : bool, default True
            Use a log-frequency x-axis.
        title : str, optional

        Returns
        -------
        matplotlib.figure.Figure
        """
        if not hasattr(self, "snr_mask_"):
            raise NotFittedError("Call fit() first.")

        import matplotlib.pyplot as plt  # noqa: PLC0415 – optional dep

        data = (
            self.snr_mask_.astype(float)
            if kind == "snr"
            else self.stack_counts_
        )

        fig, ax = plt.subplots(figsize=figsize or (12, 5))

        freqs = self.freqs_
        x = np.log10(freqs) if log_freq else freqs
        extent = [x[0], x[-1], self.n_stations_ - 0.5, -0.5]

        im = ax.imshow(
            data,
            aspect="auto",
            extent=extent,
            cmap=cmap,
            origin="upper",
            vmin=0,
            vmax=(1 if kind == "snr" else int(self.stack_counts_.max())),
        )
        fig.colorbar(im, ax=ax, label="Valid" if kind == "snr" else "Stacks")

        ax.set_xlabel(
            "log₁₀ Frequency (Hz)" if log_freq else "Frequency (Hz)"
        )
        ax.set_ylabel("Station index")
        ax.set_title(
            title
            or f"Hardware coverage ({self.n_stations_} stations, {self.n_freqs_} freqs)"
        )

        fig.tight_layout()
        return fig


# ---------------------------------------------------------------------------
# EDIBatch
# ---------------------------------------------------------------------------


class EDIBatch(PyCSAMTObject):
    """Load a directory of WinGLink-exported EDI files as :class:`EDIFile` objects.

    Files are sorted with a natural-sort key so that Stratagem's three-digit
    station numbering (``001``, ``002``, …, ``087``) is preserved regardless
    of how the OS returns directory entries.

    Parameters
    ----------
    edi_dir : path-like, optional
        Directory containing ``.edi`` files.  May also be given to
        :meth:`fit`.
    pattern : str, default ``'*.edi'``
        Glob pattern used to discover EDI files.
    verbose : int, default 0

    Attributes
    ----------
    edi_paths_ : list of Path
        Sorted list of discovered EDI file paths.
    edi_objects_ : list of EDIFile
        Successfully parsed :class:`~pycsamt.seg.edi.EDIFile` objects.
    n_stations_ : int
        Number of successfully loaded EDI files.

    Examples
    --------
    >>> batch = EDIBatch("2/2EDI").fit()
    >>> len(batch)
    87
    >>> batch[0].station          # DATAID from >HEAD
    'S00'
    >>> for edi in batch:
    ...     print(edi.station)
    """

    __repr_fields__ = ("edi_dir", "pattern", "n_stations_")

    def __init__(
        self,
        edi_dir: str | Path | None = None,
        *,
        pattern: str = "*.edi",
        verbose: int = 0,
    ) -> None:
        self.edi_dir = edi_dir
        self.pattern = pattern
        self.verbose = verbose

    # ------------------------------------------------------------------
    def fit(
        self,
        edi_dir: str | Path | None = None,
    ) -> EDIBatch:
        """Discover and load EDI files from *edi_dir*.

        Parameters
        ----------
        edi_dir : path-like, optional
            Override the directory set in ``__init__``.

        Returns
        -------
        self
        """
        if edi_dir is not None:
            self.edi_dir = edi_dir
        if self.edi_dir is None:
            raise ValueError("edi_dir must be provided")

        d = Path(self.edi_dir).expanduser().resolve()
        if not d.is_dir():
            raise FileHandlingError(f"Directory not found: {d}")

        paths = sorted(d.glob(self.pattern), key=_edi_sort_key)
        if not paths:
            raise FileHandlingError(
                f"No EDI files matching '{self.pattern}' in {d}"
            )

        self.edi_paths_: list[Path] = paths
        self.edi_objects_: list[EDIFile] = []
        skipped = 0

        for p in paths:
            try:
                self.edi_objects_.append(EDIFile(p, verbose=0))
            except Exception as exc:
                skipped += 1
                if self.verbose:
                    print(f"[EDIBatch] skip {p.name}: {exc}")

        self.n_stations_ = len(self.edi_objects_)

        if self.verbose:
            msg = f"[EDIBatch] loaded {self.n_stations_} EDI files from {d.name}"
            if skipped:
                msg += f" ({skipped} skipped)"
            print(msg)

        return self

    # ------------------------------------------------------------------
    # sequence protocol
    # ------------------------------------------------------------------

    def __len__(self) -> int:
        return getattr(self, "n_stations_", 0)

    def __getitem__(self, idx: int):
        if not hasattr(self, "edi_objects_"):
            raise NotFittedError("Call fit() first.")
        return self.edi_objects_[idx]

    def __iter__(self):
        if not hasattr(self, "edi_objects_"):
            raise NotFittedError("Call fit() first.")
        return iter(self.edi_objects_)

    # ------------------------------------------------------------------
    def station_names(self) -> list[str]:
        """Return DATAID strings from each loaded EDI's ``>HEAD`` section.

        Returns
        -------
        list of str
            Missing DATAIDs are replaced with the file stem.
        """
        if not hasattr(self, "edi_objects_"):
            raise NotFittedError("Call fit() first.")
        names = []
        for edi, p in zip(self.edi_objects_, self.edi_paths_):
            name = edi.station or p.stem
            names.append(str(name))
        return names
