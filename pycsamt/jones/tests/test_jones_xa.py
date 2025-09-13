# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

import pytest
import xarray as xr
import numpy as np

from pycsamt.jones.j import JFile
from pycsamt.jones.collection import JCollection
from pycsamt.jones.xa import (
    build_jdataset,
    XAJMixin,
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
    assert np.allclose(
        ds["z"].isel(site=0).values, single_jfile.Z.z
    )

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
    assert 'NIA000' in ds.coords["site"].values
    
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
    assert 'NIA001' in ds.coords["site"].values

def test_xamixin_meta_table(collection_with_mixin):
    """Test the meta_table method from the mixin."""
    meta_ds = collection_with_mixin.meta_table()
    assert isinstance(meta_ds, xr.Dataset)
    # meta_table now returns a Dataset with Data variables
    assert "z" not in meta_ds.data_vars
    assert "rho" not in meta_ds.data_vars
    assert "lat" in meta_ds.data_vars # Stored as a data variable
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
    assert 'NIA001' in stations

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
