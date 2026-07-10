# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""OccamData — build and parse OccamDataFile.dat from EDI sources.

The data file is the primary Occam2D input.  It encodes apparent
resistivity, phase (and optionally tipper) for selected TE/TM modes
at each station and frequency.

Entry points
------------
``OccamData.from_edi(edi_collection_or_sites, ...)``
    Build a new data file object directly from EDI data.
``OccamData.read(path)``
    Parse an existing ``OccamDataFile.dat``.
``OccamData.write(path)``
    Serialise to disk in the ``OCCAM2MTDATA_1.0`` format.

File format (OCCAM2MTDATA_1.0)
-------------------------------
::

    FORMAT:           OCCAM2MTDATA_1.0
    TITLE:            <title>
    SITES:            <n>
       S00
       ...
       S<n-1>
    OFFSETS (M):
       <x0>
       ...
       <x_{n-1}>
    FREQUENCIES:      <nf>
       <f0>
       ...
       <f_{nf-1}>
    DATA BLOCKS:      <nd>
    SITE  FREQ  TYPE   DATUM    ERROR
       <site_idx>  <freq_idx>  <type_code>  <value>  <error>
       ...
"""

from __future__ import annotations

from pathlib import Path
from typing import Union

import numpy as np

from .base import OccamBase
from .config import OccamConfig
from .doc import _occam_param_docs as _params

PathLike = Union[str, Path]

__all__ = ["OccamData"]

_FORMAT_TAG = "OCCAM2MTDATA_1.0"

# Data-type codes used in the OCCAM2MTDATA format
DATA_TYPE_CODES = {
    "RhoTE":   1,
    "PhsTE":   2,
    "RhoTM":   5,
    "PhsTM":   6,
    "ReZXX":   11,
    "ImZXX":   12,
    "ReZXY":   13,
    "ImZXY":   14,
    "ReZYX":   15,
    "ImZYX":   16,
    "ReZYY":   17,
    "ImZYY":   18,
    "ReTZX":   3,
    "ImTZX":   4,
    "ReTZY":   7,
    "ImTZY":   8,
}


# -----------------------------------------------------------------------
# Helpers for from_edi
# -----------------------------------------------------------------------

def _site_attr(site, *names):
    """Return the first non-None attribute from *names* on *site*."""
    for n in names:
        v = getattr(site, n, None)
        if v is not None:
            return v
    return None


def _unique_freqs(all_freqs: list, rtol: float = 0.01) -> np.ndarray:
    """Merge near-duplicate frequencies within tolerance."""
    arr = np.sort(np.unique(np.asarray(all_freqs, dtype=float)))[::-1]
    if arr.size == 0:
        return arr
    merged = [float(arr[0])]
    for f in arr[1:]:
        if abs(f - merged[-1]) / max(abs(f), abs(merged[-1]), 1e-30) > rtol:
            merged.append(float(f))
    return np.array(merged, dtype=float)


def _normalise_source(source) -> list:
    """Return a flat list of site-like items."""
    # Sites (pycsamt.site.base.Sites): has _items → iterate directly
    if hasattr(source, "_items"):
        return list(source)
    # EDICollection: has .edic → list of EDIFile; wrap lazily
    if hasattr(source, "edic"):
        try:
            from pycsamt.site.base import Site
            return [Site(e) for e in source.edic]
        except Exception:
            return list(source.edic)
    # Assume already iterable of site-like objects
    return list(source)


def _compute_offsets(items: list):
    """Return (names, offsets_m) sorted by profile chainage.

    Try ``Profile.from_sites`` first, then simple coordinates.
    """
    # Station names
    names = []
    for it in items:
        n = getattr(it, "name", None)
        if n is None:
            n = _site_attr(it, "station", "dataid") or f"S{len(names):03d}"
        names.append(str(n))

    # Try full Profile machinery
    try:
        from pycsamt.site.profile import Profile
        prof = Profile.from_sites(items)
        offsets = [float(prof.chainages.get(n, float("nan"))) for n in names]
        if all(np.isfinite(offsets)):
            return names, offsets
    except Exception:
        pass

    # Fallback: compute Euclidean offset from lat/lon
    lats, lons = [], []
    for it in items:
        c = getattr(it, "coords", None)
        if callable(c):
            c = c()
        if c is not None and len(c) >= 2:
            lats.append(float(c[0]))
            lons.append(float(c[1]))
        else:
            lats.append(float("nan"))
            lons.append(float("nan"))

    lat0 = next((v for v in lats if np.isfinite(v)), 0.0)
    lon0 = next((v for v in lons if np.isfinite(v)), 0.0)
    _M_PER_DEG = 111_195.0
    offsets = []
    for lat, lon in zip(lats, lons):
        if np.isfinite(lat) and np.isfinite(lon):
            dx = (lon - lon0) * _M_PER_DEG * np.cos(np.radians(lat0))
            dy = (lat - lat0) * _M_PER_DEG
            offsets.append(float(np.sqrt(dx * dx + dy * dy)) * np.sign(dx + dy))
        else:
            offsets.append(float(len(offsets)) * 1000.0)   # 1 km spacing

    return names, offsets


# -----------------------------------------------------------------------
# Low-level parser
# -----------------------------------------------------------------------

def _parse_data(path: Path) -> dict:
    """Parse an OCCAM2MTDATA_1.0 file.

    Returns
    -------
    dict with keys:
        format_str, title, sites, offsets, frequencies, and
        data_rows. Each row is [site_idx, freq_idx, type_code,
        datum, error].

    Raises
    ------
    FileNotFoundError
    ValueError
        If the Format tag is absent or wrong.
    """
    if not path.exists():
        raise FileNotFoundError(f"OccamDataFile not found: {path}")

    with path.open("r", errors="replace") as fh:
        lines = [line.rstrip("\n") for line in fh]

    result: dict = {
        "format_str":  None,
        "title":       "",
        "sites":       [],
        "offsets":     [],
        "frequencies": [],
        "data_rows":   [],
    }

    i = 0
    N = len(lines)
    n_sites = 0
    n_freq  = 0

    while i < N:
        raw     = lines[i]
        stripped = raw.strip()
        i       += 1

        if not stripped:
            continue

        if ":" not in stripped:
            continue

        raw_key, _, raw_val = stripped.partition(":")
        key_up = raw_key.strip().upper()
        val    = raw_val.strip()

        if key_up == "FORMAT":
            result["format_str"] = val
            if val.upper() != _FORMAT_TAG:
                raise ValueError(
                    f"Expected format '{_FORMAT_TAG}', got '{val}' in {path}"
                )

        elif key_up == "TITLE":
            result["title"] = val

        elif key_up == "SITES":
            n_sites = int(val)
            sites: list[str] = []
            while len(sites) < n_sites and i < N:
                s = lines[i].strip()
                i += 1
                if s:
                    sites.append(s)
            result["sites"] = sites

        elif key_up.startswith("OFFSETS"):
            offsets: list[float] = []
            while len(offsets) < n_sites and i < N:
                s = lines[i].strip()
                i += 1
                if not s:
                    continue
                try:
                    offsets.append(float(s))
                except ValueError:
                    i -= 1   # not a float → put line back
                    break
            result["offsets"] = offsets

        elif key_up == "FREQUENCIES":
            n_freq = int(val)
            freqs: list[float] = []
            while len(freqs) < n_freq and i < N:
                s = lines[i].strip()
                i += 1
                if not s:
                    continue
                try:
                    freqs.append(float(s))
                except ValueError:
                    i -= 1
                    break
            result["frequencies"] = freqs

        elif key_up.startswith("DATA BLOCKS"):
            # Skip the column-header line ("SITE  FREQ  TYPE   DATUM    ERROR")
            while i < N:
                s = lines[i].strip()
                i += 1
                if s and s.upper().startswith("SITE"):
                    break
            # Read data rows until end of file
            data_rows: list[list] = []
            while i < N:
                s = lines[i].strip()
                i += 1
                if not s:
                    continue
                tokens = s.split()
                if len(tokens) == 5:
                    try:
                        data_rows.append([
                            int(tokens[0]), int(tokens[1]), int(tokens[2]),
                            float(tokens[3]), float(tokens[4]),
                        ])
                    except ValueError:
                        pass
            result["data_rows"] = data_rows

    if result["format_str"] is None:
        raise ValueError(
            f"File does not contain a valid OCCAM2MTDATA_1.0 header: {path}"
        )

    return result


# -----------------------------------------------------------------------
# OccamData
# -----------------------------------------------------------------------

class OccamData(OccamBase):
    def __init__(
        self,
        title: str = "pycsamt Occam2D data file",
        config: OccamConfig | None = None,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.format_str:  str          = _FORMAT_TAG
        self.title:       str          = title
        self.config:      OccamConfig  = config or OccamConfig()
        self.sites:       list[str]    = []
        self.offsets:     np.ndarray   = np.array([])
        self.frequencies: np.ndarray   = np.array([])
        self.data_blocks: np.ndarray   = np.empty((0, 5))

    # ------------------------------------------------------------------
    # Construction from EDI
    # ------------------------------------------------------------------
    @classmethod
    def from_edi(
        cls,
        source,
        modes: list[str] | None = None,
        config: OccamConfig | None = None,
        title: str = "pycsamt Occam2D data file",
        **kwargs,
    ) -> OccamData:
        cfg   = config or OccamConfig()
        modes = modes or cfg.modes

        # ---- normalise source to a list of duck-typed site items ---------
        items = _normalise_source(source)
        if not items:
            raise ValueError("OccamData.from_edi: no sites in source")

        # ---- profile sort → chainages (offsets in metres) ----------------
        names, offsets = _compute_offsets(items)

        # sort by chainage
        order = np.argsort(offsets)
        names   = [names[i]   for i in order]
        offsets = [offsets[i] for i in order]
        items   = [items[i]   for i in order]

        # ---- global frequency list ---------------------------------------
        all_freqs: list[float] = []
        for it in items:
            f = _site_attr(it, "freq")
            if f is not None:
                all_freqs.extend(np.asarray(f, dtype=float).ravel().tolist())

        if not all_freqs:
            raise ValueError("OccamData.from_edi: no frequency data found in source")

        freqs = _unique_freqs(all_freqs)
        if cfg.freq_min is not None:
            freqs = freqs[freqs >= cfg.freq_min]
        if cfg.freq_max is not None:
            freqs = freqs[freqs <= cfg.freq_max]
        freqs = np.sort(freqs)[::-1]   # descending (Occam convention)

        if freqs.size == 0:
            raise ValueError(
                "OccamData.from_edi: no frequencies in the specified band "
                f"[{cfg.freq_min}, {cfg.freq_max}]"
            )

        # ---- mode → (row, col, rho_code, phs_code) ----------------------
        _MODE_MAP = {
            "TE": (0, 1, 1, 2),
            "TM": (1, 0, 5, 6),
        }
        mode_spec = {m: _MODE_MAP[m] for m in modes if m in _MODE_MAP}
        if not mode_spec:
            raise ValueError(f"OccamData.from_edi: no valid modes in {modes}")

        # ---- build data rows ---------------------------------------------
        _ln10 = float(np.log(10.0))
        data_rows: list[list] = []

        for si, (_name, site) in enumerate(zip(names, items)):
            s_freq    = _site_attr(site, "freq")
            s_rho     = _site_attr(site, "rho")
            s_phs     = _site_attr(site, "phase")
            s_rho_err = _site_attr(site, "rho_err", "z_err")
            s_phs_err = _site_attr(site, "phase_err")

            if s_freq is None or s_rho is None or s_phs is None:
                continue

            s_freq = np.asarray(s_freq, dtype=float)
            s_rho  = np.asarray(s_rho,  dtype=float)
            s_phs  = np.asarray(s_phs,  dtype=float)

            rho_err_arr = np.asarray(s_rho_err, dtype=float) if s_rho_err is not None else None
            phs_err_arr = np.asarray(s_phs_err, dtype=float) if s_phs_err is not None else None

            for fi, f in enumerate(freqs):
                # Match global frequency within 1 % tolerance
                fdiff = np.abs(s_freq - f) / np.maximum(np.abs(s_freq), np.abs(f))
                mi    = int(np.argmin(fdiff))
                if fdiff[mi] > 0.01:
                    continue

                for _mode_name, (ri, ci, rho_code, phs_code) in mode_spec.items():
                    rho_val = float(s_rho[mi, ri, ci]) if s_rho.ndim == 3 else float(s_rho[mi])
                    phs_val = float(s_phs[mi, ri, ci]) if s_phs.ndim == 3 else float(s_phs[mi])

                    if not (np.isfinite(rho_val) and rho_val > 0):
                        continue
                    if not np.isfinite(phs_val):
                        continue

                    # Normalise TM phase to 1st quadrant: phase_yx + 180°
                    if rho_code == 5:
                        phs_val = phs_val + 180.0

                    log_rho = float(np.log10(rho_val))

                    # Rho error (stored as delta_log10_rho)
                    if rho_err_arr is not None and rho_err_arr.ndim == 3:
                        re = float(rho_err_arr[mi, ri, ci])
                        rel_err = (re / rho_val) if (np.isfinite(re) and re > 0) else cfg.error_floor_rho
                    else:
                        rel_err = cfg.error_floor_rho
                    rho_err = max(rel_err, cfg.error_floor_rho) / _ln10

                    # Phase error (stored in degrees)
                    if phs_err_arr is not None and phs_err_arr.ndim == 3:
                        pe = float(phs_err_arr[mi, ri, ci])
                        phs_err = max(abs(pe), cfg.error_floor_phase) if np.isfinite(pe) else cfg.error_floor_phase
                    else:
                        phs_err = cfg.error_floor_phase

                    data_rows.append([si + 1, fi + 1, rho_code, log_rho, rho_err])
                    data_rows.append([si + 1, fi + 1, phs_code, phs_val, phs_err])

        obj = cls(title=title, config=cfg, **kwargs)
        obj.sites       = names
        obj.offsets     = np.array(offsets, dtype=float)
        obj.frequencies = freqs
        if data_rows:
            obj.data_blocks = np.array(data_rows, dtype=float)

        if obj.verbose:
            obj.logger.info(
                "OccamData.from_edi: %d sites, %d freqs, %d data blocks",
                obj.n_sites, obj.n_frequencies, obj.n_data,
            )
        return obj

    # ------------------------------------------------------------------
    # I/O
    # ------------------------------------------------------------------
    @classmethod
    def read(cls, path: PathLike, **kwargs) -> OccamData:
        p      = Path(path)
        parsed = _parse_data(p)
        obj    = cls(**kwargs)

        obj.format_str  = parsed["format_str"]
        obj.title       = parsed["title"]
        obj.sites       = parsed["sites"]
        obj.offsets     = np.array(parsed["offsets"], dtype=float)
        obj.frequencies = np.array(parsed["frequencies"], dtype=float)

        rows = parsed["data_rows"]
        if rows:
            obj.data_blocks = np.array(rows, dtype=float)
        else:
            obj.data_blocks = np.empty((0, 5))

        if obj.verbose:
            obj.logger.info(
                "OccamData.read: %d sites, %d freqs, %d data blocks from %s",
                obj.n_sites, obj.n_frequencies, obj.n_data, p,
            )
        return obj

    def write(self, path: PathLike) -> Path:
        p = Path(path)
        p.parent.mkdir(parents=True, exist_ok=True)

        lines: list[str] = [
            f"FORMAT:           {self.format_str}\n",
            f"TITLE:            {self.title}\n",
            f"SITES:            {self.n_sites}\n",
        ]
        for s in self.sites:
            lines.append(f"   {s}\n")

        lines.append("OFFSETS (M):\n")
        for off in self.offsets:
            lines.append(f"   {off:.1f}\n")

        lines.append(f"FREQUENCIES:      {self.n_frequencies}\n")
        for f in self.frequencies:
            lines.append(f"   {f:.7e}\n")

        lines.append(f"DATA BLOCKS:      {self.n_data}\n")
        lines.append("SITE  FREQ  TYPE   DATUM    ERROR\n")
        for row in self.data_blocks:
            s_i = int(row[0])
            f_i = int(row[1])
            t_c = int(row[2])
            dat = row[3]
            err = row[4]
            lines.append(f"{s_i:5d} {f_i:5d} {t_c:5d} {dat:12.4f} {err:12.4f}\n")

        with p.open("w") as fh:
            fh.writelines(lines)
        self.path = p
        return p

    # ------------------------------------------------------------------
    # Convenience
    # ------------------------------------------------------------------
    @property
    def n_sites(self) -> int:
        return len(self.sites)

    @property
    def n_frequencies(self) -> int:
        return len(self.frequencies)

    @property
    def n_data(self) -> int:
        return len(self.data_blocks)

    @property
    def type_codes(self) -> np.ndarray:
        """Unique data-type codes present in this dataset."""
        if not self.data_blocks.size:
            return np.array([], dtype=int)
        return np.unique(self.data_blocks[:, 2].astype(int))


OccamData.__doc__ = rf"""
Represent an Occam2D magnetotelluric data file.

``OccamData`` stores the station list, profile offsets,
global frequency table, data-type codes, datum values, and
uncertainty values written to ``OccamDataFile.dat``. The
object is both a container for parsed files and the product
of EDI conversion by :meth:`from_edi`.

The Occam2D data file uses one row for each datum. Apparent
resistivity is stored in logarithmic form, while phase is
stored in degrees:

.. math::

    d_\rho = \log_{{10}}(\rho_a),
    \qquad
    \sigma_d = \frac{{\sigma_\rho}}{{\ln(10)}} .

For PyCSAMT EDI arrays, TE mode is taken from
:math:`Z_{{xy}}` and TM mode is taken from :math:`Z_{{yx}}`.
TM phase is shifted by :math:`180^\circ` so passive-MT
:math:`Z_{{yx}}` phases are written in the first quadrant.

Parameters
----------
{_params.data.title}
{_params.common.config}
{_params.common.verbose}
{_params.common.logger}

Attributes
----------
format_str : str
    Occam format tag. The current writer uses
    ``"OCCAM2MTDATA_1.0"``.
title : str
    Free-text title written to the data-file header.
config : OccamConfig
    Configuration used for default modes, frequency limits,
    and error floors during EDI conversion.
sites : list of str
    Station names ordered along the profile. The same order is
    used by mesh construction and all one-based site indices.
offsets : numpy.ndarray of float, shape (n_sites,)
    Station chainages in metres. Values are sorted from low to
    high during :meth:`from_edi`.
frequencies : numpy.ndarray of float, shape (n_frequencies,)
    Global frequency table in hertz, sorted from high to low
    as expected by Occam2D.
data_blocks : numpy.ndarray of float, shape (n_data, 5)
    Data rows with columns ``site_index``, ``freq_index``,
    ``type_code``, ``datum``, and ``error``. Indices are
    one-based because they are written directly to Occam
    files.

Notes
-----
Occam2D type codes distinguish both data kind and component.
The common MT rows are ``1`` for ``RhoTE``, ``2`` for
``PhsTE``, ``5`` for ``RhoTM``, and ``6`` for ``PhsTM``.
Additional impedance and tipper codes are exposed through
``DATA_TYPE_CODES`` for readers and future writers.

See Also
--------
OccamConfig
    Supplies default modes, frequency bounds, and error
    floors.
OccamMesh.from_data
    Builds mesh geometry from station offsets in
    ``OccamData``.
OccamResponse
    Reads modeled responses and residuals for the same rows.
InputBuilder
    Coordinates writing data, mesh, model, and startup files.

Examples
--------
Build a data file from EDI sites and write it to disk:

>>> from pycsamt.models.occam2d import OccamData
>>> from pycsamt.site import Sites
>>> sites = Sites.from_dir("edi")
>>> data = OccamData.from_edi(sites, modes=["TE", "TM"])
>>> data.write("occam_run/OccamDataFile.dat")

Create a synthetic container for tests or scripted workflows:

>>> import numpy as np
>>> from pycsamt.models.occam2d import OccamData
>>> data = OccamData(title="synthetic profile")
>>> data.sites = ["S00", "S01"]
>>> data.offsets = np.array([0.0, 1000.0])
>>> data.frequencies = np.array([100.0, 10.0])
>>> data.data_blocks = np.array([[1, 1, 1, 2.0, 0.05]])

Read an existing Occam data file:

>>> from pycsamt.models.occam2d import OccamData
>>> data = OccamData.read("occam_run/OccamDataFile.dat")
>>> data.n_sites, data.n_frequencies, data.n_data

References
----------
.. [1] deGroot-Hedlin, C., and Constable, S.,
   "Occam's inversion to generate smooth, two-dimensional
   models from magnetotelluric data", Geophysics, 55(12),
   1613-1624, 1990.
.. [2] Constable, S. C., Parker, R. L., and Constable,
   C. G., "Occam's inversion: A practical algorithm for
   generating smooth models from electromagnetic sounding
   data", Geophysics, 52(3), 289-300, 1987.
"""

OccamData.from_edi.__func__.__doc__ = rf"""
Build an Occam data object from EDI-derived stations.

This constructor normalizes the accepted input source to a
list of site-like objects, estimates station chainages,
merges all available frequencies into a common descending
frequency table, applies frequency limits, and writes TE/TM
apparent-resistivity and phase rows using Occam type codes.

The EDI arrays are interpreted with the convention

.. math::

    \mathrm{{TE}} = Z_{{xy}},
    \qquad
    \mathrm{{TM}} = Z_{{yx}} .

For each accepted apparent-resistivity value :math:`\rho_a`,
the stored datum is :math:`\log_{{10}}(\rho_a)`. The
relative resistivity floor is converted to log10 uncertainty
by :math:`\sigma_d=\sigma_\rho/\ln(10)`. Phase rows use
degree errors and the TM phase is shifted by
:math:`180^\circ`.

Parameters
----------
{_params.data.source}
{_params.data.modes}
{_params.common.config}
{_params.data.title}
**kwargs
    Additional keyword arguments passed to the ``OccamData``
    constructor. This is commonly used for ``verbose`` or
    ``logger`` when progress messages are desired.

Returns
-------
OccamData
    Populated data object ready to be written as an
    ``OCCAM2MTDATA_1.0`` file.

Raises
------
ValueError
    Raised when the source has no sites, no frequency arrays,
    no frequencies remain after filtering, or no requested
    mode is supported.

See Also
--------
OccamConfig
    Provides default modes, frequency bounds, and error
    floors.
OccamData.write
    Serializes the returned object to ``OccamDataFile.dat``.
OccamMesh.from_data
    Uses the returned station offsets to build mesh geometry.

Examples
--------
Build TE and TM rows from a site collection:

>>> from pycsamt.models.occam2d import OccamData
>>> from pycsamt.site import Sites
>>> sites = Sites.from_dir("edi")
>>> data = OccamData.from_edi(sites, modes=["TE", "TM"])

Restrict the frequency range through ``OccamConfig``:

>>> from pycsamt.models.occam2d import OccamConfig
>>> from pycsamt.models.occam2d import OccamData
>>> cfg = OccamConfig(freq_min=0.1, freq_max=1000.0)
>>> data = OccamData.from_edi(sites, config=cfg)

Use only TM data with stronger phase floor:

>>> cfg = OccamConfig(modes=["TM"], error_floor_phase=1.0)
>>> data = OccamData.from_edi(sites, config=cfg)

References
----------
.. [1] deGroot-Hedlin, C., and Constable, S.,
   "Occam's inversion to generate smooth, two-dimensional
   models from magnetotelluric data", Geophysics, 55(12),
   1613-1624, 1990.
"""

OccamData.read.__func__.__doc__ = rf"""
Read an existing ``OCCAM2MTDATA_1.0`` file.

The parser reads the format tag, title, station names,
offsets, frequency table, and numeric data block. Site and
frequency indices are kept in the one-based form used by
Occam2D so a read-write round trip preserves the original
file structure.

Parameters
----------
{_params.common.path}
**kwargs
    Additional keyword arguments forwarded to the
    ``OccamData`` constructor before parsed values are
    attached. Use this for ``config``, ``verbose``, or
    ``logger``.

Returns
-------
OccamData
    Parsed data-file container with arrays populated from
    ``path``.

Raises
------
FileNotFoundError
    Raised when ``path`` does not exist.
ValueError
    Raised when the format tag is missing or not
    ``"OCCAM2MTDATA_1.0"``.

Examples
--------
>>> from pycsamt.models.occam2d import OccamData
>>> data = OccamData.read("occam_run/OccamDataFile.dat")
>>> data.type_codes
"""

OccamData.write.__doc__ = rf"""
Write this object as an ``OCCAM2MTDATA_1.0`` file.

The writer serializes the current title, station list,
offsets, frequency table, and data rows using the Occam2D
text layout.
Parent directories are created before writing. The method does
not modify the object, so it can be used repeatedly for
round-trip checks or alternative run directories.

Parameters
----------
{_params.common.path}

Returns
-------
pathlib.Path
    Path to the file that was written.

See Also
--------
OccamData.read
    Parses a file written by this method.
InputBuilder.build
    Calls this method as part of complete input generation.

Examples
--------
>>> from pycsamt.models.occam2d import OccamData
>>> data = OccamData.read("source/OccamDataFile.dat")
>>> written = data.write("copy/OccamDataFile.dat")
"""
