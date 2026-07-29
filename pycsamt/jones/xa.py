# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from ..log.logger import get_logger
from ..utils.handlers import columns_manager
from .j import JFile

logger = get_logger(__name__)

__all__ = ["XAJMixin", "JFileAcc", "build_jdataset"]


def build_jdataset(
    jfiles: Iterable[JFile],
    *,
    drop_empty: bool = True,
) -> xr.Dataset:
    r"""
    Build a multi-site xarray Dataset from JFile objects.

    This function iterates through a collection of parsed
    :class:`~pycsamt.jones.j.JFile` objects, converts each one
    into a single-site :class:`xarray.Dataset`, and then
    concatenates them into a unified, multi-site dataset.

    Parameters
    ----------
    jfiles : Iterable of :class:`~pycsamt.jones.j.JFile`
        An iterable (e.g., a list or a
        :class:`~pycsamt.jones.collection.JCollection`)
        of parsed J-file objects.
    drop_empty : bool, default=True
        If ``True``, any :class:`JFile` object that contains no
        frequency data will be skipped and excluded from the
        final dataset.

    Returns
    -------
    xr.Dataset
        A single dataset containing data from all valid J-files.
        The dataset is indexed by a ``site`` dimension, and
        site-specific metadata (latitude, longitude, etc.) are
        stored as non-dimensional coordinates aligned with this
        dimension.

    Notes
    -----
    The resulting dataset is structured with dimensions for sites,
    frequencies, and tensor components. This structure is ideal
    for vectorized computations and advanced plotting across
    multiple sites.

    Data variables include impedance (`z`), tipper (`tip`),
    apparent resistivity (`rho`), phase (`phi`), their associated
    errors, and crucial data quality flags (`z_rej`, `rho_rej`).

    See Also
    --------
    JCollection.to_xarray : A convenient wrapper around this function.
    JFile : The per-item reader that provides the source data.
    JFileAcc : An accessor for interacting with the created dataset.

    Examples
    --------
    >>> from pycsamt.jones import JCollection, build_jdataset
    >>> # Create a collection of J-files from a directory
    >>> j_collection = JCollection.from_sources("data/jc/*")
    >>> # Build the xarray dataset
    >>> ds = build_jdataset(j_collection)
    >>> print(ds)
    <xarray.Dataset>
    Dimensions:    (site: 4, freq: 17, ...)
    Coordinates:
      * site       (site) object 'NIA000' 'NIA001' ...
        ...
    Data variables:
        z          (site, freq, output_ch, input_ch) complex128 ...
        z_err      (site, freq, output_ch, input_ch) float64 ...
        ...
    """
    datasets = []
    metadata_records = []

    jfiles = columns_manager(jfiles, empty_as_none=False)
    for jf in jfiles:
        try:
            ds = _ds_from_jfile(jf)
            if drop_empty and ds.sizes.get("freq", 0) == 0:
                continue
            datasets.append(ds)
            # Collect metadata for each site
            metadata_records.append(_meta_from_jfile(jf))
        except Exception as exc:
            logger.debug("Skip %s: %s", jf, exc)

    if not datasets:
        # Return an empty dataset with expected coordinates
        return xr.Dataset(
            coords={"site": [], "freq": [], "output_ch": [], "input_ch": []}
        )

    # Concatenate the data variables
    full_ds = xr.concat(datasets, dim="site", join="outer")

    # Handle Metadata
    # Create a DataFrame from metadata and assign as non-dimensional coords
    meta_df = pd.DataFrame(metadata_records).set_index("site")
    for col in meta_df.columns:
        # Use .assign_coords to add metadata to the 'site' dimension
        full_ds = full_ds.assign_coords({col: ("site", meta_df[col])})

    return full_ds


class XAJMixin:
    r"""
    A mixin to add xarray export capabilities to J-file collections.

    This class is designed to be mixed into a collection class that
    is iterable and yields :class:`~pycsamt.jones.j.JFile` objects,
    such as :class:`~pycsamt.jones.collection.JCollection`. It
    provides high-level methods to convert the entire collection
    into an :class:`xarray.Dataset`.

    Methods
    -------
    to_xarray()
        Converts the collection into a full :class:`xarray.Dataset`.
    meta_table()
        Extracts only the site-level metadata into a compact Dataset.

    See Also
    --------
    JCollection : An example of a class that can use this mixin.
    build_jdataset : The underlying function that performs the conversion.

    Examples
    --------
    >>> from pycsamt.jones import JCollection
    >>> from pycsamt.jones.xa import XAJMixin
    >>>
    >>> # Create a new class that inherits from JCollection and the mixin
    >>> class JCollectionWithXA(JCollection, XAJMixin):
    ...     pass
    >>>
    >>> # Use the new class to load data
    >>> j_collection_xa = JCollectionWithXA.from_sources("data/jc/")
    >>>
    >>> # Now you can directly call the .to_xarray() method
    >>> ds = j_collection_xa.to_xarray()
    >>> print(ds.dims)
    Frozen({'site': 4, 'freq': 17, 'output_ch': 2, 'input_ch': 2, 'tcomp': 2})
    """

    def to_xarray(
        self,
        *,
        drop_empty: bool = True,
    ) -> xr.Dataset:
        """
        Converts the entire collection into an xarray Dataset.

        This method is a convenient wrapper around the
        :func:`build_jdataset` function.
        """

        # Host must be iterable over JFile
        return build_jdataset(self, drop_empty=drop_empty)

    def meta_table(self) -> xr.Dataset:
        """
        Extracts site-level metadata into a new Dataset.

        This method provides a quick way to get a summary of all sites
        in the collection without loading the bulky tensor data.
        """
        rows: list[dict[str, object]] = []
        for jf in self:
            m = _meta_from_jfile(jf)
            m["site"] = _site_id_from_jfile(jf)
            rows.append(m)
        if not rows:
            return xr.Dataset(coords={"site": []})

        cols = {}
        for k in rows[0].keys():
            cols[k] = ("site", [r.get(k) for r in rows])

        ds = xr.Dataset(coords={"site": cols.pop("site")[1]})
        for k, v in cols.items():
            ds[k] = v
        return ds


@xr.register_dataset_accessor("jfile")
class JFileAcc:
    r"""
    An xarray accessor for convenient interaction
    with J-format datasets.

    This accessor is registered under the ``.jfile`` namespace and
    provides a set of domain-specific methods and properties for
    :class:`xarray.Dataset` objects that were created by
    :func:`build_jdataset`. It simplifies common tasks like
    selecting sites, filtering by frequency, and accessing metadata.

    Properties
    ----------
    stations : list[str]
        A list of all station/site names in the dataset.

    Methods
    -------
    get(site)
        Selects a single site by name and returns a new Dataset.
    components()
        Lists the available impedance tensor component names.
    band(fmin, fmax)
        Filters the dataset to a specific frequency band.

    Notes
    -----
    The accessor pattern is a powerful feature of xarray that
    allows external libraries to add their own functionality to
    xarray objects without modifying the xarray codebase itself.

    See Also
    --------
    build_jdataset :
        The function that creates compatible datasets.
    xarray.register_dataset_accessor :
        The decorator used to create accessors.

    Examples
    --------
    >>> from pycsamt.jones import JCollection, build_jdataset
    >>> j_collection = JCollection.from_sources("data/jc/")
    >>> ds = build_jdataset(j_collection)
    >>>
    >>> # Use the .jfile accessor
    >>> print(ds.jfile.stations)
    ['NIA000', 'NIA001', 'NIA002', 'NIA003']
    >>>
    >>> # Get data for a single site
    >>> site_nia001 = ds.jfile.get("NIA001")
    >>> print(site_nia001.dims)
    Frozen({'freq': 17, 'output_ch': 2, 'input_ch': 2, 'tcomp': 2})
    """

    def __init__(self, ds: xr.Dataset) -> None:
        self._ds = ds

    @property
    def stations(self) -> list[str]:
        s = self._ds.coords.get("site", None)
        return [] if s is None else [str(v) for v in s.data]

    def get(self, site: str) -> xr.Dataset:
        """
        Selects data for a single site.

        This is a convenience wrapper around ``.sel(site=...)`` which
        also cleans up the resulting dataset by removing the now-singleton
        'site' dimension.
        """
        return self._ds.sel(site=str(site))

    def components(self) -> list[str]:
        """Returns a list of available impedance component names."""
        comps = []
        if "output_ch" in self._ds.coords and "input_ch" in self._ds.coords:
            for och in self._ds["output_ch"].values:
                for ich in self._ds["input_ch"].values:
                    comps.append(f"z{och.lower()}{ich.lower()}")
        return comps or ["zxx", "zxy", "zyx", "zyy"]

    def band(
        self,
        fmin: float | None = None,
        fmax: float | None = None,
    ) -> xr.Dataset:
        """
        Filters the dataset to a specified frequency band.

        Parameters
        ----------
        fmin : float, optional
            The minimum frequency to include (inclusive).
        fmax : float, optional
            The maximum frequency to include (inclusive).

        Returns
        -------
        xr.Dataset
            A new dataset containing only the data within the specified
            frequency range.
        """
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
            degree range.
        figsize : tuple[int, int], default=(8, 8)
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

        if show_grid:
            axes[0].grid(**default_grid_props)
            axes[1].grid(**default_grid_props)

        axes[0].xaxis.set_major_formatter(mticker.LogFormatterSciNotation())
        axes[0].yaxis.set_major_formatter(mticker.LogFormatterSciNotation())

        axes[0].legend()
        axes[1].legend()
        fig.suptitle(f"Site: {site}", fontsize=14)
        plt.tight_layout(rect=[0, 0, 1, 0.96])

        if savefig:
            plt.savefig(savefig, dpi=300)

        return fig, axes

    def attrs(self) -> dict[str, object]:
        """
        Returns the global attributes of the Dataset.

        .. deprecated:: 2.1.0
            Metadata is now primarily stored as non-dimensional
            coordinates. Use direct coordinate access for site-specific
            metadata (e.g., `ds.coords['lat']`). This method is
            maintained for backward compatibility.
        """
        return dict(self._ds.attrs)


def _site_id_from_jfile(jf: JFile) -> str:
    """Infer site ID with fallbacks."""
    return jf.site or "unknown_site"


def _meta_from_jfile(jf: JFile) -> dict[str, object]:
    """Extract metadata from a JFile object."""
    p = jf.path
    software = None
    if jf.heads and jf.heads.banner:
        software = jf.heads.banner.software

    return {
        "site": _site_id_from_jfile(jf),  # Added for indexing
        "path": str(p) if isinstance(p, Path) else None,
        "filename": p.name if isinstance(p, Path) else None,
        "dataid": jf.site,
        "lat": jf.lat,
        "lon": jf.lon,
        "elev": jf.elev,
        "azimuth": jf.azimuth,
        "has_tip": jf.Tip is not None,
        "nfreq": jf.n_freq,
        "has_z": jf.Z is not None,
        "has_r": jf.Res is not None,
        "software": software,
    }


def _ds_from_jfile(jf: JFile) -> xr.Dataset:
    """
    Creates a single-site xarray Dataset from one JFile object.

    This version is refactored for clarity, robustness, and includes
    critical data quality (rejection) flags.
    """
    sid = _site_id_from_jfile(jf)
    f = np.asarray(getattr(jf, "freq", []), float)
    n_freq = f.size

    # --- Transfer Function Data (Z) ---
    z = _get_tensor_or_zeros(jf.Z, "z", n_freq, np.complex128)
    z_err = _get_tensor_or_zeros(jf.Z, "z_err", n_freq, np.float64)
    z_rej = _get_rejection_flags(jf, "Z", n_freq)

    zrot_val = getattr(jf.Z, "rotation_angle", np.zeros(n_freq))
    zrot = (
        np.asarray(zrot_val).astype(np.float64)
        if zrot_val.size == n_freq
        else np.zeros(n_freq, dtype=np.float64)
    )

    # --- Resistivity and Phase Data (R/S) ---
    rho = _get_tensor_or_zeros(jf.Res, "resistivity", n_freq, np.float64)
    phi = _get_tensor_or_zeros(jf.Res, "phase", n_freq, np.float64)
    rho_err = _get_tensor_or_zeros(jf.Res, "resistivity_err", n_freq, np.float64)
    phi_err = _get_tensor_or_zeros(jf.Res, "phase_err", n_freq, np.float64)
    rho_rej = _get_rejection_flags(jf, "R", n_freq)

    # --- Tipper Data (T) ---
    tip_val = getattr(jf.Tip, "tipper", None)
    tip = (
        np.asarray(tip_val)
        if tip_val is not None
        else np.zeros((n_freq, 1, 2), dtype=np.complex128)
    )
    tip_da = (
        tip[:, 0, :]
        if tip.ndim == 3 and tip.shape[1] == 1
        else np.zeros((n_freq, 2), dtype=np.complex128)
    )

    tip_err_val = getattr(jf.Tip, "tipper_err", None)
    tip_err = (
        np.asarray(tip_err_val)
        if tip_err_val is not None
        else np.zeros_like(tip_da, dtype=np.float64)
    )
    tip_err_da = (
        tip_err[:, 0, :]
        if tip_err.ndim == 3 and tip_err.shape[1] == 1
        else np.zeros((n_freq, 2), dtype=np.float64)
    )

    ds = xr.Dataset(
        data_vars={
            "z": (("freq", "output_ch", "input_ch"), z),
            "z_err": (("freq", "output_ch", "input_ch"), z_err),
            "z_rej": (("freq", "output_ch", "input_ch"), z_rej),
            "zrot": (("freq",), zrot),
            "rho": (("freq", "output_ch", "input_ch"), rho),
            "rho_err": (("freq", "output_ch", "input_ch"), rho_err),
            "rho_rej": (("freq", "output_ch", "input_ch"), rho_rej),
            "phi": (("freq", "output_ch", "input_ch"), phi),
            "phi_err": (("freq", "output_ch", "input_ch"), phi_err),
            "tip": (("freq", "tcomp"), tip_da),
            "tip_err": (("freq", "tcomp"), tip_err_da),
        },
        coords={
            "freq": f,
            # Using more descriptive coordinate
            # names improves clarity
            "output_ch": ["Hx", "Hy"],
            "input_ch": ["Hx", "Hy"],
            "tcomp": ["Tx", "Ty"],
        },
    ).expand_dims(site=[sid])

    return ds


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

    # Fallback for safety
    return np.zeros((n_freq, 2, 2), dtype=dtype)


def _get_rejection_flags(jf: JFile, kind_prefix: str, n_freq: int) -> np.ndarray:
    """
    Extracts rejection flags from the original JBlocks for a given kind.

    This is necessary because the high-level Z/Res objects do not
    preserve the per-component rejection flags.
    """
    rej_tensor = np.zeros((n_freq, 2, 2), dtype=bool)
    if jf.blocks is None or jf.freq is None:
        return rej_tensor

    # Map component suffix to tensor index, e.g., "XY" -> (0, 1)
    comp_map = {"XX": (0, 0), "XY": (0, 1), "YX": (1, 0), "YY": (1, 1)}

    for block in jf.blocks.select(kind=kind_prefix):
        comp = block.comp
        if comp not in comp_map:
            continue

        i, j = comp_map[comp]
        block_data = block.to_numpy()

        # Align periods between the block and the final JFile object
        p_common, idx_jfile, idx_block = JFile._align_by_periods(
            jf.freq, block_data["period"]
        )

        if p_common.size > 0:
            rej_tensor[idx_jfile, i, j] = block_data["rej"][idx_block]

    return rej_tensor
