# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from ..log.logger import get_logger
from .edi import EDIFile

logger = get_logger(__name__)

__all__ = ["XAMixin", "EDIAcc", "build_dataset"]


def build_dataset(
    edis: Iterable[EDIFile],
    *,
    drop_empty: bool = True,
) -> xr.Dataset:
    r"""
    Build a multi-site xarray Dataset from an iterable of EDIFile.

    This function iterates through a collection of parsed
    :class:`~pycsamt.seg.edi.EDIFile` objects, converts each one
    into a single-site :class:`xarray.Dataset`, and then
    concatenates them into a unified, multi-site dataset.

    Parameters
    ----------
    edis : Iterable of :class:`~pycsamt.seg.edi.EDIFile`
        An iterable (e.g., a list or a
        :class:`~pycsamt.seg.collection.EDICollection`) of parsed EDI
        objects.
    drop_empty : bool, default=True
        If ``True``, any :class:`EDIFile` object that contains no
        frequency data will be skipped and excluded from the final
        dataset.

    Returns
    -------
    xr.Dataset
        A single dataset containing data from all valid EDI files.
        The dataset is indexed by a ``site`` dimension, and
        site-specific metadata (latitude, longitude, etc.) are
        stored as non-dimensional coordinates aligned with this
        dimension.

    Notes
    -----
    The resulting dataset is structured with dimensions for sites,
    frequencies, and tensor components. This structure is ideal for
    vectorized computations and advanced plotting across multiple
    sites.

    This function correctly handles site-specific metadata by
    assigning it to coordinates, preventing data loss during
    concatenation, which is a common pitfall when storing metadata
    in global attributes.

    See Also
    --------
    EDICollection.to_xarray : A convenient wrapper around this function.
    EDIFile : The per-item reader that provides the source data.
    EDIAcc : An accessor for interacting with the created dataset.

    Examples
    --------
    >>> from pycsamt.seg import EDICollection, build_dataset
    >>> # Create a collection of EDI files
    >>> edi_collection = EDICollection.from_sources("data/edis/")
    >>> # Build the xarray dataset
    >>> ds = build_dataset(edi_collection)
    >>> print(ds)
    <xarray.Dataset>
    Dimensions:      (site: 2, freq: 60, ...)
    Coordinates:
      * site         (site) object 'S01' 'S02'
      * freq         (freq) float64 320.0 286.9 ...
        ...
        lat          (site) float64 26.05 26.05
        lon          (site) float64 -10.33 -10.33
    Data variables:
        z            (site, freq, output_ch, input_ch) complex128 ...
        z_err        (site, freq, output_ch, input_ch) float64 ...
        ...
    """

    datasets = []
    metadata_records = []
    for ed in edis:
        try:
            ds = _ds_from_edi(ed)
            is_empty_tf = ds.sizes.get("freq", 0) == 0
            has_sp_or_ts = any(k in ds.data_vars for k in ("spec_vals", "ts"))
            if drop_empty and is_empty_tf and not has_sp_or_ts:
                continue
            datasets.append(ds)
            metadata_records.append(_meta_from_edi(ed))
        except Exception as exc:
            logger.debug("Skip %s: %s", ed, exc)

    if not datasets:
        return xr.Dataset(
            coords={"site": [], "freq": [], "output_ch": [], "input_ch": []}
        )

    full_ds = xr.concat(datasets, dim="site", join="outer")

    if metadata_records:
        meta_df = pd.DataFrame(metadata_records).set_index("site")
        for col in meta_df.columns:
            full_ds = full_ds.assign_coords({col: ("site", meta_df[col])})

    return full_ds


class XAMixin:
    r"""
    A mixin that adds convenient xarray exports to collection classes.

    This mixin provides a bridge between collection-like objects
    (such as :class:`~pycsamt.seg.collection.EDICollection`) and
    the powerful, multi-dimensional data structures offered by
    the `xarray` library. Any class that is iterable over
    :class:`~pycsamt.seg.edi.EDIFile` instances can inherit from this
    mixin to gain methods for data conversion and metadata
    extraction.

    Methods
    -------
    to_xarray(drop_empty=True)
        Converts the entire collection into a single, comprehensive
        :class:`xarray.Dataset`. This method leverages
        :func:`build_dataset` to handle the conversion and
        concatenation of multiple EDI files.
    meta_table()
        Extracts only the site-level metadata (e.g., coordinates,
        filenames, data quality flags) from the collection and
        returns it as a clean, tabular :class:`xarray.Dataset`,
        omitting the bulky transfer function data.

    Notes
    -----
    - This mixin is designed to be lightweight and does not impose any
      specific storage or indexing strategy on the host class; it only
      requires that the host class implements the `__iter__` method
      to yield :class:`~pycsamt.seg.edi.EDIFile` objects.
    - The site identifiers used in the resulting datasets follow the
      same robust inference rules as :func:`build_dataset`.

    See Also
    --------
    build_dataset : The core function that performs the conversion.
    EDICollection : A primary user of this mixin.
    EDIAcc : The accessor for interacting with the created dataset.

    Examples
    --------
    To use this mixin, simply inherit from it in your collection class.

    >>> from pycsamt.seg.edi import EDIFile
    >>> class MyEDICollection(XAMixin):
    ...     def __init__(self, items):
    ...         self._items = list(items)
    ...
    ...     def __iter__(self):
    ...         return iter(self._items)
    >>> # Assume "site1.edi" and "site2.edi" exist
    >>> edi_files = [
    ...     EDIFile("data/edis/S01.edi"),
    ...     EDIFile("data/edis/S02.edi"),
    ... ]
    >>> collection = MyEDICollection(edi_files)
    >>>
    >>> # Convert the entire collection to an xarray Dataset
    >>> ds = collection.to_xarray()
    >>> print(ds.site.values)
    ['S01' 'S02']
    >>>
    >>> # Get a summary table of just the metadata
    >>> metadata_ds = collection.meta_table()
    >>> print(metadata_ds[["lat", "lon"]])
    <xarray.Dataset>
    Dimensions:  (site: 2)
    Coordinates:
      * site     (site) object 'S01' 'S02'
    Data variables:
        lat      (site) float64 26.05 26.05
        lon      (site) float64 -10.33 -10.33
    """

    def to_xarray(
        self,
        *,
        drop_empty: bool = True,
    ) -> xr.Dataset:
        # host must be iterable over EDIFile
        return build_dataset(self, drop_empty=drop_empty)

    def meta_table(self) -> xr.Dataset:
        """Extracts site-level metadata into a new Dataset."""
        rows: list[dict[str, object]] = []
        for ed in self:
            rows.append(_meta_from_edi(ed))
        if not rows:
            return xr.Dataset(coords={"site": []})

        meta_df = pd.DataFrame(rows).set_index("site")
        return xr.Dataset.from_dataframe(meta_df)


@xr.register_dataset_accessor("edi")
class EDIAcc:
    r"""
    An xarray accessor for convenient interaction with EDI datasets.

    This accessor is registered under the ``.edi`` namespace and
    provides domain-specific methods and properties for datasets
    created by :func:`build_dataset`. It simplifies common data
    selection and visualization tasks that are specific to MT/EM
    (magnetotelluric/electromagnetic) data.

    Properties
    ----------
    stations : list[str]
        A list of all unique station or site names present in the
        dataset's ``site`` coordinate.

    Methods
    -------
    get(site)
        Selects and returns a new :class:`xarray.Dataset` containing
        data for only a single site, specified by its name. The
        selection is case-insensitive.
    band(fmin=None, fmax=None)
        Filters the dataset to a specific frequency range. Returns a
        new dataset containing only the data within the inclusive
        frequency bounds.
    plot_apparent_resistivity(site, **kwargs)
        Generates a standard plot of apparent resistivity and phase
        curves for the off-diagonal tensor components (XY and YX) of
        a specified site.
    attrs()
        Returns a dictionary of the dataset's global attributes.

    See Also
    --------
    build_dataset : The function used to create datasets compatible
                    with this accessor.

    Examples
    --------
    >>> from pycsamt.seg import EDICollection, build_dataset
    >>> edi_collection = EDICollection.from_sources("data/edis/")
    >>> ds = build_dataset(edi_collection)
    >>>
    >>> # Get a list of all station names
    >>> print(ds.edi.stations)
    ['S01', 'S02', 'S03', ...]
    >>>
    >>> # Select data for a single station (case-insensitive)
    >>> site_data = ds.edi.get("s01")
    >>>
    >>> # Filter the data to a specific frequency band (e.g., 1 to 100 Hz)
    >>> filtered_ds = ds.edi.band(fmin=1.0, fmax=100.0)
    >>>
    >>> # Create a standard plot for a site
    >>> fig, axes = ds.edi.plot_apparent_resistivity(site="S01")
    >>> # fig.show() # Uncomment to display plot
    """

    def __init__(self, ds: xr.Dataset) -> None:
        self._ds = ds

    @property
    def stations(self) -> list[str]:
        s = self._ds.coords.get("site", None)
        return [] if s is None else [str(v) for v in s.data]

    def get(self, site: str) -> xr.Dataset:
        """Selects data for a single site (case-insensitive)."""
        # Find the correctly cased site name from the coordinates
        site_coord = self._ds.coords["site"].values
        try:
            match = next(s for s in site_coord if s.upper() == site.upper())
            return self._ds.sel(site=match)
        except StopIteration:
            raise KeyError(f"Site '{site}' not found in dataset.")

    def plot_apparent_resistivity(
        self,
        site: str,
        components: list[str] = None,
        phase_mod: int | None = None,
        figsize: tuple[int, int] = (8, 8),
        show_grid: bool = True,
        grid_props: dict | None = None,
        savefig: str | None = None,
        **plot_kwargs,
    ):
        r"""
        Generates a standard plot of apparent resistivity and phase.

        This method provides a flexible interface for visualizing MT
        (magnetotelluric) data, allowing customization of components,
        phase wrapping, and plot aesthetics.

        Parameters
        ----------
        site : str
            The site identifier to plot.
        components : list of str, default=["xy", "yx"]
            A list of tensor components to plot (e.g., "xy", "yx", "xx").
            The selection is case-insensitive.
        phase_mod : int, optional
            If provided, wraps the phase to a specific quadrant. For
            example, ``phase_mod=90`` will display phases in the [0, 90]
            degree range, useful for visualizing data in a single quadrant.
        figsize : tuple[int, int], default=(8, 6)
            The figure size for the plot.
        show_grid : bool, default=True
            Whether to display a grid on both subplots.
        grid_props : dict, optional
            Additional properties to customize the grid lines (e.g.,
            ``{'color': 'grey', 'linestyle': '--', 'linewidth': 0.5}``).
        savefig : str, optional
            If a path is provided, the plot will be saved to that file.
        **plot_kwargs :
            Additional keyword arguments passed directly to xarray's
            ``.plot.line()`` method for customizing the lines.

        Returns
        -------
        fig : matplotlib.figure.Figure
            The matplotlib Figure object.
        axes : np.ndarray of matplotlib.axes.Axes
            An array containing the two subplot Axes objects.

        Examples
        --------
        >>> # Basic plot of off-diagonal components
        >>> fig, axes = ds.edi.plot_apparent_resistivity(site="S01")
        >>> # fig.show()

        >>> # Plot all components and save the figure
        >>> fig, axes = ds.edi.plot_apparent_resistivity(
        ...     site="S01",
        ...     components=["xy", "yx", "xx", "yy"],
        ...     savefig="S01_all_components.png",
        ... )

        >>> # Plot with phase wrapped to the first quadrant and custom styling
        >>> fig, axes = ds.edi.plot_apparent_resistivity(
        ...     site="S01",
        ...     phase_mod=90,
        ...     grid_props={"color": "red", "linestyle": ":"},
        ...     marker="o",  # passed to plot.line
        ... )
        """
        import matplotlib.pyplot as plt
        import matplotlib.ticker as mticker

        if components is None:
            components = ["xy", "yx"]
        ds_site = self.get(site)
        fig, axes = plt.subplots(2, 1, sharex=True, figsize=figsize)

        comp_map = {
            "xy": ("Hx", "Hy"),
            "yx": ("Hy", "Hx"),
            "xx": ("Hx", "Hx"),
            "yy": ("Hy", "Hy"),
        }

        default_grid_props = {
            "which": "both",
            "linestyle": "--",
            "linewidth": 0.5,
        }
        if grid_props:
            default_grid_props.update(grid_props)

        for comp in components:
            comp_lower = comp.lower()
            if comp_lower not in comp_map:
                logger.warning(f"Component '{comp}' is not valid. Skipping.")
                continue

            output_ch, input_ch = comp_map[comp_lower]

            # --- Resistivity Plot (Log-Log) ---
            rho_data = ds_site["rho"].sel(output_ch=output_ch, input_ch=input_ch)
            rho_data.plot.line(
                ax=axes[0],
                xscale="log",
                yscale="log",
                label=f"$\\rho_{{{comp_lower}}}$",
                **plot_kwargs,
            )

            # --- Phase Plot (Semi-Log) ---
            phi_data = ds_site["phi"].sel(output_ch=output_ch, input_ch=input_ch)
            if phase_mod is not None and isinstance(phase_mod, int):
                phi_data = phi_data % phase_mod

            phi_data.plot.line(
                ax=axes[1],
                xscale="log",
                label=f"$\\phi_{{{comp_lower}}}$",
                **plot_kwargs,
            )

        # --- Aesthetics and Formatting ---
        axes[0].set_ylabel("Apparent Resistivity (Ω·m)")
        axes[0].set_xlabel("")  # Remove x-label from top plot

        axes[1].set_ylabel("Phase (degrees)")
        axes[1].set_xlabel("Frequency (Hz)")

        # Grid formatting
        if show_grid:
            axes[0].grid(**default_grid_props)
            axes[1].grid(**default_grid_props)

        # Use scientific notation for log axes
        axes[0].xaxis.set_major_formatter(mticker.LogFormatterSciNotation())
        axes[0].yaxis.set_major_formatter(mticker.LogFormatterSciNotation())

        axes[0].legend()
        axes[1].legend()
        fig.suptitle(f"Site: {site}", fontsize=14)
        plt.tight_layout(rect=[0, 0, 1, 0.96])  # Adjust for suptitle

        if savefig:
            plt.savefig(savefig, dpi=300)

        return fig, axes

    def attrs(self) -> dict[str, object]:
        """Returns the global attributes of the Dataset."""
        return dict(self._ds.attrs)

    def band(
        self,
        fmin: float | None = None,
        fmax: float | None = None,
    ) -> xr.Dataset:
        ds = self._ds
        if "freq" not in ds.coords:
            return ds
        cond = np.ones(ds.sizes["freq"], bool)
        fv = ds["freq"].data
        if fmin is not None:
            cond &= fv >= float(fmin)
        if fmax is not None:
            cond &= fv <= float(fmax)
        return ds.isel(freq=np.where(cond)[0])

    # --- spectra helpers ---
    def has_spectra(self) -> bool:
        return "spec_vals" in self._ds

    def spectra(self) -> xr.Dataset:
        keep = [k for k in self._ds.data_vars if k.startswith("spec_")]
        return self._ds[keep] if keep else self._ds

    # --- time-series helpers ---
    def has_timeseries(self) -> bool:
        return "ts" in self._ds

    def timeseries(self) -> xr.Dataset:
        keep = [k for k in ("ts", "time", "dt", "npts")]
        keep = [k for k in keep if k in self._ds]
        return self._ds[keep] if keep else self._ds


def _site_id_from_edi(ed: EDIFile) -> str:
    """Infer site ID with fallbacks."""
    return ed.station or "unknown_site"


def _site_location(
    ed: EDIFile,
) -> tuple[float | None, float | None, float | None]:
    """Site (lat, lon, elev), falling back to DEFINEMEAS REFLAT/REFLONG/REFELEV
    when the HEAD section carries no LAT/LONG/ELEV (common for BIRRP/JONES
    processed EDI, e.g. kap03lmt)."""
    head = ed.get_section("head")
    dm = ed.get_section("definemeas")

    lat = getattr(head, "lat", None) if head else None
    if lat is None:
        lat = getattr(dm, "reflat", None) if dm else None

    lon = getattr(head, "long", None) if head else None
    if lon is None:
        lon = getattr(dm, "reflong", None) if dm else None

    elev = getattr(head, "elev", None) if head else None
    if elev is None:
        elev = getattr(dm, "refelev", None) if dm else None

    return lat, lon, elev


def _meta_from_edi(ed: EDIFile) -> dict[str, object]:
    """Extract metadata from an EDIFile object for a single site."""
    p = ed.path
    software = ed.processingsoftware
    lat, lon, elev = _site_location(ed)

    return {
        "site": _site_id_from_edi(ed),
        "path": str(p) if isinstance(p, Path) else None,
        "filename": p.name if isinstance(p, Path) else None,
        "dataid": ed.station,
        "lat": lat,
        "lon": lon,
        "elev": elev,
        "has_tip": ed.has_tipper,
        "nfreq": ed.n_freq,
        "has_spec": ed.get_section("spectra") is not None,
        "has_ts": ed.get_section("timeseries") is not None,
        "software": software,
    }


def _get_tensor_or_zeros(
    obj: object | None, attr: str, n_freq: int, dtype: np.dtype
) -> np.ndarray:
    """Safely get a (n_freq, 2, 2) tensor array or return zeros."""
    val = getattr(obj, attr, None) if obj else None
    if val is None:
        return np.zeros((n_freq, 2, 2), dtype=dtype)

    arr = np.asarray(val)
    if arr.ndim == 3 and arr.shape == (n_freq, 2, 2):
        return arr.astype(dtype, copy=False)

    return np.zeros((n_freq, 2, 2), dtype=dtype)


def _get_tipper_or_zeros(
    obj: object | None, attr: str, n_freq: int, dtype: np.dtype
) -> np.ndarray:
    """Safely get a (n_freq, 2) tipper array or return zeros."""
    val = getattr(obj, attr, None) if obj else None
    if val is None:
        return np.zeros((n_freq, 2), dtype=dtype)

    arr = np.asarray(val)
    # Handle the standard (n_freq, 1, 2) shape
    if arr.ndim == 3 and arr.shape == (n_freq, 1, 2):
        return arr[:, 0, :].astype(dtype, copy=False)
    # Handle cases where it might already be 2D
    if arr.ndim == 2 and arr.shape == (n_freq, 2):
        return arr.astype(dtype, copy=False)

    return np.zeros((n_freq, 2), dtype=dtype)


def _site_id(ed: EDIFile) -> str:
    # preference: station -> HEAD.dataid
    # -> Spectra.name -> path.stem -> "site"
    head = ed.get_section("head")
    spec = ed.get_section("spectra")
    for cand in (
        getattr(ed, "station", None),
        getattr(head, "dataid", None) if head else None,
        getattr(spec, "name", None) if spec else None,
    ):
        if cand:
            return str(cand)
    p = getattr(ed, "path", None)
    return p.stem if isinstance(p, Path) else "site"


def _meta(ed: EDIFile) -> dict[str, object]:
    head = ed.get_section("head")
    info = ed.get_section("info")
    spec = ed.get_section("spectra")

    p = getattr(ed, "path", None)

    dataid = getattr(head, "dataid", None) if head else None
    if dataid is None:
        dataid = getattr(spec, "name", None) or getattr(ed, "station", None)
        if dataid is None:
            p = getattr(ed, "path", None)
            dataid = p.stem if isinstance(p, Path) else None

    # --- robust software extraction ---
    software = None
    proc = getattr(info, "Processing", None)
    if proc is not None and hasattr(proc, "ProcessingSoftware"):
        sw = proc.ProcessingSoftware
        # accept either an object with .name or a plain string
        software = getattr(sw, "name", sw)

    lat, lon, elev = _site_location(ed)

    return {
        "path": str(p) if isinstance(p, Path) else None,
        "filename": p.name if isinstance(p, Path) else None,
        "dataid": dataid,
        "lat": lat,
        "lon": lon,
        "elev": elev,
        "has_tip": bool(
            getattr(ed.Tip, "tipper", None) is not None
            and getattr(ed.Tip.tipper, "size", 0) > 0
        ),
        "nfreq": int(getattr(ed.Z, "n_freq", 0) or 0),
        "has_spec": ed.get_section("spectra") is not None,
        "has_ts": ed.get_section("timeseries") is not None,
        "software": software,
    }


def _pad2d(
    seqs: list[np.ndarray],
    *,
    pad: float = np.nan,
) -> np.ndarray:
    n = len(seqs)
    m = max((s.size for s in seqs), default=0)
    out = np.full((n, m), pad, float)
    for i, a in enumerate(seqs):
        if a.size:
            out[i, : a.size] = a
    return out


def _spec_pack(ed: EDIFile) -> dict[str, xr.DataArray]:
    spec = ed.get_section("spectra")
    if spec is None:
        return {}

    # Try direct attributes first
    freq = getattr(spec, "freq", None)
    vals = getattr(spec, "values", None)
    bw = getattr(spec, "bw", None)
    avgt = getattr(spec, "avgt", None)
    rs = getattr(spec, "rotspec", None)

    if freq is None or vals is None:
        try:
            sect, io = spec.to_io()  # may not exist on all impls
            blks = getattr(io, "blocks", [])
            # prefer block.options['freq'] if present
            freq = [
                (
                    getattr(b, "freq", None)
                    if getattr(b, "freq", None) is not None
                    else b.options.get("freq", None)
                )
                for b in blks
            ]
            vals = [np.asarray(getattr(b, "values", []), float) for b in blks]
            bw = [getattr(b, "bw", None) for b in blks]
            avgt = [getattr(b, "avgt", None) for b in blks]
            rs = [getattr(b, "rotspec", None) for b in blks]
        except Exception:
            return {}  # don’t break the whole dataset

    # Sanitize freq
    try:
        freq_clean = [float(x) for x in freq if x is not None]
    except Exception:
        return {}  # invalid freq list

    f = np.asarray(freq_clean, float)
    if f.size == 0:
        return {}

    # Values: ragged → padded
    try:
        vlist = [np.asarray(v, float) for v in vals]
    except Exception:
        return {}

    spec_vals = _pad2d(vlist, pad=np.nan)
    spec_len = np.array([v.size for v in vlist], int)
    sidx = np.arange(spec_vals.shape[1], dtype=int)

    out: dict[str, xr.DataArray] = {
        "spec_vals": xr.DataArray(
            spec_vals,
            dims=("freq", "sidx"),
            coords={"freq": f, "sidx": sidx},
        ),
        "spec_len": xr.DataArray(spec_len, dims=("freq",), coords={"freq": f}),
    }

    # Optional per-freq metadata (only if lengths match)
    def _ok_1d(a):
        try:
            a = np.asarray(a, float)
            return a.size == f.size, a
        except Exception:
            return False, None

    for name, arr in (
        ("spec_bw", bw),
        ("spec_avgt", avgt),
        ("spec_rotspec", rs),
    ):
        ok, arrv = _ok_1d(arr)
        if ok:
            out[name] = xr.DataArray(arrv, dims=("freq",), coords={"freq": f})

    return out


def _ts_pack(ed: EDIFile) -> dict[str, xr.DataArray]:
    ts = ed.get_section("timeseries")
    if ts is None:
        return {}

    # channels, dt per channel, sequences
    try:
        ch = list(ts.channels())
    except Exception:
        return {}

    if not ch:
        return {}

    series = [np.asarray(ts.get(c), float) for c in ch]
    npts = [a.size for a in series]
    mat = _pad2d(series, pad=np.nan)  # (ch_len→columns)
    # we want (sample, ch)
    mat = mat.T

    dt = np.array([float(ts.dt_map.get(c, np.nan)) for c in ch])
    nmax = mat.shape[0]
    # time grid per column with padding as nan
    T = np.full((nmax, len(ch)), np.nan, float)
    for j, (a, d) in enumerate(zip(series, dt)):
        if np.isfinite(d) and a.size:
            T[: a.size, j] = np.arange(a.size, float) * d

    out: dict[str, xr.DataArray] = {}
    out["ts"] = xr.DataArray(
        mat,
        dims=("sample", "ch"),
        coords={"sample": np.arange(nmax), "ch": ch},
    )
    out["time"] = xr.DataArray(
        T,
        dims=("sample", "ch"),
        coords={"sample": np.arange(nmax), "ch": ch},
    )
    out["dt"] = xr.DataArray(dt, dims=("ch",), coords={"ch": ch})
    out["npts"] = xr.DataArray(np.array(npts, int), dims=("ch",), coords={"ch": ch})
    return out


def _ds_from_edi(ed: EDIFile) -> xr.Dataset:
    """
    Creates a single-site xarray Dataset from one EDIFile object.

    This function is refactored for clarity and completeness, leveraging
    helper functions to handle data extraction and optional data blocks
    like spectra and time-series.
    """
    sid = _site_id_from_edi(ed)

    # --- Determine primary frequency array ---
    # Prefer Z.freq, but fall back to Spectra.freq if Z is empty.
    f = np.asarray(getattr(ed.Z, "freq", []), dtype=float)
    if f.size == 0:
        spec = ed.get_section("spectra")
        if spec:
            f = np.asarray(getattr(spec, "freq", []), dtype=float)
    n_freq = f.size

    # --- Extract Transfer Function Data using helpers ---
    z = _get_tensor_or_zeros(ed.Z, "z", n_freq, dtype=complex)
    z_err = _get_tensor_or_zeros(ed.Z, "z_err", n_freq, dtype=float)
    rho = _get_tensor_or_zeros(ed.Z, "resistivity", n_freq, dtype=float)
    phi = _get_tensor_or_zeros(ed.Z, "phase", n_freq, dtype=float)
    rho_err = _get_tensor_or_zeros(ed.Z, "resistivity_err", n_freq, dtype=float)
    phi_err = _get_tensor_or_zeros(ed.Z, "phase_err", n_freq, dtype=float)

    zrot_val = getattr(ed.Z, "rotation_angle", np.zeros(n_freq))
    zrot = (
        np.asarray(zrot_val)
        if zrot_val.size == n_freq
        else np.zeros(n_freq, dtype=np.float64)
    )

    tip = _get_tipper_or_zeros(ed.Tip, "tipper", n_freq, np.complex128)
    tip_err = _get_tipper_or_zeros(ed.Tip, "tipper_err", n_freq, np.float64)

    # --- Create the base Dataset ---
    ds = xr.Dataset(
        data_vars={
            "z": (("freq", "output_ch", "input_ch"), z),
            "z_err": (("freq", "output_ch", "input_ch"), z_err),
            "rho": (("freq", "output_ch", "input_ch"), rho),
            "phi": (("freq", "output_ch", "input_ch"), phi),
            "rho_err": (("freq", "output_ch", "input_ch"), rho_err),
            "phi_err": (("freq", "output_ch", "input_ch"), phi_err),
            "zrot": (("freq",), zrot),
            "tip": (("freq", "tcomp"), tip),
            "tip_err": (("freq", "tcomp"), tip_err),
        },
        coords={
            "freq": f,
            "output_ch": ["Hx", "Hy"],
            "input_ch": ["Hx", "Hy"],
            "tcomp": ["Tx", "Ty"],
        },
    ).expand_dims(site=[sid])

    # --- Add Spectra and Time-Series data ---
    try:
        sp_data = _spec_pack(ed)
        if sp_data:
            ds = ds.merge(xr.Dataset(sp_data))
    except NameError:
        # _spec_pack not defined, skipping.
        pass

    try:
        ts_data = _ts_pack(ed)
        if ts_data:
            ds = ds.merge(xr.Dataset(ts_data))
    except NameError:
        # _ts_pack not defined, skipping.
        pass

    return ds


def _ds_from_edi_v1(ed: EDIFile) -> xr.Dataset:
    sid = _site_id(ed)

    # prefer Z freq; if empty, fall back to Spectra freq
    f = np.asarray(getattr(ed.Z, "freq", []), float)
    if f.ndim != 1 or f.size == 0:
        sp = ed.get_section("spectra")
        sf = (
            np.asarray(getattr(sp, "freq", []), float)
            if sp is not None
            else np.asarray([], float)
        )
        f = sf if (sf.ndim == 1 and sf.size > 0) else np.asarray([], float)

    def _block_or_zeros(a, n: int, dtype) -> np.ndarray:
        if a is None:
            return np.zeros((n, 2, 2), dtype)
        arr = np.asarray(a)
        # reject scalar/object junk or wrong shape
        if arr.ndim != 3 or arr.shape != (n, 2, 2):
            try:
                arr = arr.reshape(n, 2, 2)
            except Exception:
                return np.zeros((n, 2, 2), dtype)
        return arr.astype(dtype, copy=False)

    # Z as complex (freq, i, j)
    z = _block_or_zeros(getattr(ed.Z, "z", None), f.size, complex)
    ze = _block_or_zeros(getattr(ed.Z, "z_err", None), f.size, float)

    zrot = np.asarray(getattr(ed.Z, "rotation_angle", np.zeros(f.size))).astype(float)
    if zrot.size != f.size:
        zrot = np.zeros(f.size, float)

    # tipper (freq, 1, 2) -> (freq, comp)
    tip = getattr(ed.Tip, "tipper", None)
    if tip is None or getattr(tip, "size", 0) == 0:
        tip = np.zeros((f.size, 1, 2), complex)
    tip = np.asarray(tip)
    tip_da = tip[:, 0, :] if tip.ndim == 3 else np.zeros((f.size, 2), complex)

    terr = getattr(ed.Tip, "_tipper_err", None)
    if terr is None:
        terr = np.zeros_like(tip, float)
    terr = np.asarray(terr)
    terr_da = terr[:, 0, :] if terr.ndim == 3 else np.zeros((f.size, 2), float)

    # --- helper: coerce to 1D length-n safely
    def _as1d(v, n: int) -> np.ndarray:
        if v is None:
            return np.zeros(n, float)
        a = np.asarray(v, float)
        if a.ndim == 0:
            return np.full(n, float(a))
        a = a.ravel()
        out = np.zeros(n, float)
        k = min(n, a.size)
        out[:k] = a[:k]
        return out

    # --- helper: try to read Z property without crashing
    def _get1d(name: str, n: int) -> np.ndarray:
        try:
            v = getattr(ed.Z, name)  # may raise ResistivityError
        except Exception:
            v = None
        return _as1d(v, n)

    # Optional scalar forms (ρ, φ), safe when Z is absent
    rho = np.stack(
        [
            _get1d("res_xx", f.size),
            _get1d("res_xy", f.size),
            _get1d("res_yx", f.size),
            _get1d("res_yy", f.size),
        ],
        axis=1,
    ).reshape(f.size, 2, 2)

    phi = np.stack(
        [
            _get1d("phase_xx", f.size),
            _get1d("phase_xy", f.size),
            _get1d("phase_yx", f.size),
            _get1d("phase_yy", f.size),
        ],
        axis=1,
    ).reshape(f.size, 2, 2)

    # base TF dataset
    ds = xr.Dataset(
        data_vars={
            "z": (("freq", "i", "j"), z),
            "z_err": (("freq", "i", "j"), ze),
            "zrot": (("freq",), zrot),
            "tip": (("freq", "tcomp"), tip_da),
            "tip_err": (("freq", "tcomp"), terr_da),
            "rho": (("freq", "i", "j"), rho),
            "phi": (("freq", "i", "j"), phi),
        },
        coords={
            "freq": f,
            "i": [0, 1],
            "j": [0, 1],
            "tcomp": ["tx", "ty"],
        },
        attrs=_meta(ed),
    ).expand_dims(site=[sid])

    # spectra/time-series if present
    sp = _spec_pack(ed)
    if sp:
        for k, v in sp.items():
            ds[k] = v

    ts = _ts_pack(ed)
    if ts:
        for k, v in ts.items():
            ds[k] = v

    return ds
