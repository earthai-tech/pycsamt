# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later

from __future__ import annotations

from collections.abc import Iterable, Iterator
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pytest

from pycsamt.seg.collection import EDICollection
from pycsamt.seg.edi import EDIFile
from pycsamt.seg.xa import XAMixin, _site_location, build_dataset

# Mark all tests in this file as requiring the 'xarray' package
pytestmark = pytest.mark.requires_xarray


class _DummyColl(XAMixin):
    """Minimal iterable collection to test XAMixin."""

    def __init__(self, items: Iterable[EDIFile]) -> None:
        self._items: list[EDIFile] = list(items)

    def __iter__(self) -> Iterator[EDIFile]:
        return iter(self._items)


@pytest.mark.usefixtures("edi_imp_file")
def test_build_dataset_single_imp(edi_imp_file: Path) -> None:
    """Test building a dataset from a single EDIFile with impedance."""
    ed = EDIFile(edi_imp_file)
    ds = build_dataset([ed], drop_empty=False)

    # Core coordinates and dimensions should be present
    assert "site" in ds.coords
    assert "freq" in ds.coords
    # Verify the new, descriptive coordinate names
    assert set(("output_ch", "input_ch")).issubset(ds.coords)

    # Core variables present with expected shapes
    assert "z" in ds and "z_err" in ds and "zrot" in ds
    nfreq = ds.sizes["freq"]
    nsite = ds.sizes.get("site", 1)

    assert ds["z"].shape == (nsite, nfreq, 2, 2)
    assert ds["z_err"].shape == (nsite, nfreq, 2, 2)
    assert ds["tip"].shape == (nsite, nfreq, 2)
    assert np.iscomplexobj(ds["z"].data)

    # Verify metadata is stored as a non-dimensional coordinate
    assert "lat" in ds.coords
    assert "lon" in ds.coords
    assert ds["lat"].dims == ("site",)
    assert ds["lon"].dims == ("site",)

    # Check a specific metadata value (falls back to DEFINEMEAS REFLAT
    # when HEAD carries no LAT, as with this BIRRP/JONES-processed file)
    expected_lat, _, _ = _site_location(ed)
    assert np.isclose(ds["lat"].item(), expected_lat)


@pytest.mark.usefixtures("edi_collection")
def test_edi_accessor_basic(edi_collection: EDICollection) -> None:
    """Test the basic functionality of the .edi accessor."""
    ds = build_dataset(edi_collection)

    # .stations property
    stations = ds.edi.stations
    print("stations===", stations)
    assert isinstance(stations, list)
    assert len(stations) == len(edi_collection)
    assert "S01" in stations

    # .get() method (case-insensitive)
    sub = ds.edi.get("s01")
    assert "site" in sub.coords
    assert "site" not in sub.sizes  # Dimension is dropped
    assert str(sub.site.values) == "S01"

    # .band() frequency selection
    f = ds["freq"].values
    if f.size >= 3:
        fmin, fmax = f[1], f[-2]
        slim = ds.edi.band(fmin=fmin, fmax=fmax)
        assert 0 < slim.sizes["freq"] < ds.sizes["freq"]
    else:
        slim = ds.edi.band()
        assert slim.sizes["freq"] == ds.sizes["freq"]


@pytest.mark.usefixtures("edi_imp_file")
def test_edi_accessor_plot(edi_imp_file: Path) -> None:
    """Test that the plotting method on the accessor is callable."""
    ed = EDIFile(edi_imp_file)
    ds = build_dataset([ed])
    site_name = ds.edi.stations[0]

    try:
        fig, axes = ds.edi.plot_apparent_resistivity(site=site_name)
        assert isinstance(fig, plt.Figure)
        assert isinstance(axes, np.ndarray)
    finally:
        plt.close("all")  # Ensure plots are closed after test


@pytest.mark.usefixtures("edi_spe_file")
def test_build_dataset_with_spectra(edi_spe_file: Path) -> None:
    """Test building a dataset from an EDI file containing spectra."""
    ed = EDIFile(edi_spe_file)
    ds = build_dataset([ed], drop_empty=False)

    assert "site" in ds.coords
    assert ds.sizes["site"] >= 1

    # Check for spectra variables if the accessor finds them
    if ds.edi.has_spectra():
        assert "spec_vals" in ds
        assert "spec_len" in ds
        assert "freq" in ds["spec_vals"].coords

        sp = ds.edi.spectra()
        assert "spec_vals" in sp and "spec_len" in sp
    else:
        pytest.skip("No spectra variables present in this EDI sample.")


@pytest.mark.usefixtures("edi_imp_file", "edi_spe_file")
def test_xamixin_collection_meta(
    edi_imp_file: Path,
    edi_spe_file: Path,
) -> None:
    """Test the XAMixin methods on a dummy collection."""
    ed1 = EDIFile(edi_imp_file)
    ed2 = EDIFile(edi_spe_file)
    coll = _DummyColl([ed1, ed2])

    ds = coll.to_xarray(drop_empty=False)
    assert "site" in ds.coords
    assert ds.sizes["site"] == len(ds.edi.stations)
    assert ds.sizes["site"] == 2

    meta = coll.meta_table()
    assert "site" in meta.coords
    assert meta.sizes["site"] == ds.sizes["site"]
    for k in ("filename", "nfreq", "has_tip", "has_spec", "has_ts"):
        assert k in meta
