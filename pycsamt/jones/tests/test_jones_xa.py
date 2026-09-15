# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

import numpy as np
import pytest
import xarray as xr

from pycsamt.jones.collection import JCollection
from pycsamt.jones.j import JFile
from pycsamt.jones.xa import (
    XAJMixin,
    build_jdataset,
)

# Mark all tests in this file as requiring the 'xarray' package
pytestmark = pytest.mark.requires_xarray


@pytest.fixture(scope="module")
def single_jfile(j_single_file):
    """Parsed JFile object from a single file fixture."""
    return JFile.from_file(j_single_file)


@pytest.fixture(scope="module")
def jcollection(jc_files):
    """JCollection object from multiple file fixtures."""
    return JCollection.from_sources(jc_files)


def test_build_dataset_single_file(single_jfile):
    """Test building a dataset from a single JFile."""
    ds = build_jdataset([single_jfile])

    assert isinstance(ds, xr.Dataset)
    assert "site" in ds.dims
    assert ds.sizes["site"] == 1
    assert ds.coords["site"].values[0] == "KB0001"

    # Test data variables and shapes with new coordinate names
    assert "z" in ds.data_vars
    assert ds["z"].shape == (1, single_jfile.n_freq, 2, 2)
    assert ds["z"].dims == ("site", "freq", "output_ch", "input_ch")
    assert np.allclose(ds["z"].isel(site=0).values, single_jfile.Z.z)

    # Test for rejection flags
    assert "z_rej" in ds.data_vars
    assert ds["z_rej"].dtype == bool
    assert "rho_rej" in ds.data_vars
    assert ds["rho_rej"].dtype == bool


def test_build_dataset_from_collection(jcollection):
    """Test building a dataset from a JCollection."""
    ds = build_jdataset(jcollection)

    assert isinstance(ds, xr.Dataset)
    assert ds.sizes["site"] == len(jcollection)
    assert "NIA000" in ds.coords["site"].values

    # Test for metadata as non-dimensional coordinates
    assert "lat" in ds.coords
    assert "lon" in ds.coords
    assert "azimuth" in ds.coords
    assert ds["lat"].dims == ("site",)

    # Check a specific metadata value
    s01_lat = jcollection.get("NIA001", "lat")
    assert np.isclose(ds["lat"].sel(site="NIA001").item(), s01_lat)


def test_build_dataset_empty():
    """Test empty dataset creation from empty input."""
    ds = build_jdataset([])
    assert isinstance(ds, xr.Dataset)
    assert ds.sizes["site"] == 0
    assert "output_ch" in ds.coords
    assert "input_ch" in ds.coords


class JCollWithMixin(JCollection, XAJMixin):
    """Test class combining JCollection with the XAMixin."""

    pass


@pytest.fixture(scope="module")
def collection_with_mixin(jc_files):
    """Instance of the test collection with the XA mixin."""
    return JCollWithMixin.from_sources(jc_files)


def test_xamixin_to_xarray(collection_with_mixin):
    """Test the to_xarray method from the mixin."""
    ds = collection_with_mixin.to_xarray()
    assert isinstance(ds, xr.Dataset)
    assert ds.sizes["site"] == len(collection_with_mixin)
    assert "NIA001" in ds.coords["site"].values


def test_xamixin_meta_table(collection_with_mixin):
    """Test the meta_table method from the mixin."""
    meta_ds = collection_with_mixin.meta_table()
    assert isinstance(meta_ds, xr.Dataset)
    # meta_table now returns a Dataset with Data variables
    assert "z" not in meta_ds.data_vars
    assert "rho" not in meta_ds.data_vars
    assert "lat" in meta_ds.data_vars  # Stored as a data variable
    assert meta_ds.sizes["site"] == len(collection_with_mixin)


@pytest.fixture(scope="module")
def main_dataset(jcollection):
    """A dataset built from the collection for accessor tests."""
    return build_jdataset(jcollection)


def test_jfileacc_stations(main_dataset):
    """Test the .jfile.stations accessor property."""
    stations = main_dataset.jfile.stations
    assert isinstance(stations, list)
    assert len(stations) == main_dataset.sizes["site"]
    assert "NIA001" in stations


def test_jfileacc_get_site(main_dataset):
    """Test the .jfile.get() accessor method."""
    site_ds = main_dataset.jfile.get("NIA000")
    assert isinstance(site_ds, xr.Dataset)
    # After sel, the dimension 'site' is dropped.
    assert "site" not in site_ds.dims
    assert site_ds.coords["dataid"] == "NIA000"


def test_jfileacc_band(main_dataset):
    """Test the .jfile.band() frequency slicing method."""
    fmin, fmax = 1, 10
    band_ds = main_dataset.jfile.band(fmin=fmin, fmax=fmax)
    assert isinstance(band_ds, xr.Dataset)
    assert np.all(band_ds.coords["freq"].values >= fmin)
    assert np.all(band_ds.coords["freq"].values <= fmax)


def test_jfileacc_attrs_compatibility(main_dataset):
    """Test the .jfile.attrs accessor for backward compatibility."""
    # Note: main metadata is now in coords. attrs() might be deprecated
    # in the future or its behavior clarified. For now, test it works.
    attrs = main_dataset.jfile.attrs()
    assert isinstance(attrs, dict)


def test_jfileacc_components(main_dataset):
    """Test the component naming from the accessor."""
    comps = main_dataset.jfile.components()
    assert isinstance(comps, list)
    assert "zhxhy" in comps
    assert "zhyhx" in comps


# ─────────────────────────────────────────────────────────────────────────
# Additional coverage: build_jdataset drop_empty / exception skip
# ─────────────────────────────────────────────────────────────────────────


class _EmptyFreqJF:
    """Stand-in JFile with zero frequencies -> dropped by drop_empty."""

    site = "EMPTY"
    path = None
    heads = None
    lat = None
    lon = None
    elev = None
    azimuth = None
    Tip = None
    Z = None
    Res = None
    blocks = None
    n_freq = 0
    freq = np.array([])


class _BoomJF:
    """Stand-in JFile whose conversion always raises -> skipped + logged."""

    site = "BOOM"

    def __getattr__(self, name):
        raise RuntimeError("boom")


def test_build_jdataset_drops_empty_freq_when_flag_set():
    ds = build_jdataset([_EmptyFreqJF()], drop_empty=True)
    assert ds.sizes["site"] == 0


def test_build_jdataset_skips_items_that_raise():
    ds = build_jdataset([_BoomJF()])
    assert ds.sizes["site"] == 0


def test_xamixin_meta_table_empty_collection():
    class _EmptyColl(list, XAJMixin):
        pass

    ds = _EmptyColl().meta_table()
    assert isinstance(ds, xr.Dataset)
    assert ds.sizes.get("site", 0) == 0


def test_jfileacc_band_without_freq_coord_returns_unchanged():
    from pycsamt.jones.xa import JFileAcc

    ds = xr.Dataset(coords={"site": ["S1"]})
    acc = JFileAcc(ds)
    out = acc.band(fmin=1, fmax=10)
    assert out is ds


def test_jfileacc_plot_apparent_resistivity_smoke(main_dataset):
    import matplotlib

    matplotlib.use("Agg")
    site = main_dataset.jfile.stations[0]
    fig, axes = main_dataset.jfile.plot_apparent_resistivity(
        site, components=["xy", "bogus"]
    )
    assert fig is not None
    assert len(axes) == 2


# ─────────────────────────────────────────────────────────────────────────
# _get_tensor_or_zeros / _get_rejection_flags direct coverage
# ─────────────────────────────────────────────────────────────────────────


def test_get_tensor_or_zeros_falls_back_on_shape_mismatch():
    from pycsamt.jones.xa import _get_tensor_or_zeros

    class _Obj:
        z = np.zeros((3, 3))  # wrong shape for n_freq=2

    out = _get_tensor_or_zeros(_Obj(), "z", 2, np.complex128)
    assert out.shape == (2, 2, 2)
    assert np.all(out == 0)


def test_get_rejection_flags_returns_zeros_when_blocks_or_freq_missing():
    from pycsamt.jones.xa import _get_rejection_flags

    class _NoBlocksJF:
        blocks = None
        freq = np.array([1.0])

    out = _get_rejection_flags(_NoBlocksJF(), "Z", 1)
    assert out.shape == (1, 2, 2)
    assert not out.any()

    class _NoFreqJF:
        blocks = object()
        freq = None

    out2 = _get_rejection_flags(_NoFreqJF(), "Z", 1)
    assert not out2.any()


def test_get_rejection_flags_skips_unknown_component(single_jfile):
    from pycsamt.jones.xa import _get_rejection_flags

    class _Block:
        comp = "ZZ"  # not in comp_map -> skipped

    class _Blocks:
        def select(self, kind):
            return [_Block()]

    class _JFLike:
        blocks = _Blocks()
        freq = single_jfile.freq

    out = _get_rejection_flags(_JFLike(), "Z", single_jfile.n_freq)
    assert not out.any()
