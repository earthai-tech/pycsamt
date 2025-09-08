# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Iterator, List

import numpy as np
import pytest
import xarray as xr

from pycsamt.seg.edi import EDIFile
from pycsamt.seg.xa import build_dataset, XAMixin


class _DummyColl(XAMixin):
    """Minimal iterable collection to test XAMixin."""

    def __init__(self, items: Iterable[EDIFile]) -> None:
        self._items: List[EDIFile] = list(items)

    def __iter__(self) -> Iterator[EDIFile]:
        return iter(self._items)


@pytest.mark.usefixtures("edi_imp_file")
def test_build_dataset_single_imp(edi_imp_file: Path) -> None:
    ed = EDIFile(edi_imp_file)
    ds = build_dataset([ed], drop_empty=False)

    # core coords / dims present
    assert "site" in ds.coords
    assert "freq" in ds.coords
    assert set(("i", "j")).issubset(ds.coords)

    # core variables present with expected shapes
    assert "z" in ds and "z_err" in ds and "zrot" in ds
    nfreq = ds.sizes["freq"]
    nsite = ds.sizes.get("site", 1)

    # variables include 'site' as the leading dim
    assert ds["z"].shape == (nsite, nfreq, 2, 2)
    assert ds["z_err"].shape == (nsite, nfreq, 2, 2)
    assert ds["tip"].shape == (nsite, nfreq, 2)
    assert np.iscomplexobj(ds["z"].data)

    # per-site convenience: select one site, then check per-site shapes
    ds1 = ds.isel(site=0)
    assert ds1["tip"].shape == (nfreq, 2)
    assert ds1["tip_err"].shape == (nfreq, 2)


@pytest.mark.usefixtures("edi_imp_file")
def test_edi_accessor_basic(edi_imp_file: Path) -> None:
    ed = EDIFile(edi_imp_file)
    ds = build_dataset([ed], drop_empty=False)

    # stations list and get()
    stations = ds.edi.stations
    assert isinstance(stations, list) and len(stations) == 1

    # get() drops the 'site' dimension (scalar coord remains)
    sub = ds.edi.get(stations[0])
    assert "site" in sub.coords
    assert "site" not in sub.sizes

    # components + z_as_comp mapping (select one site to get (freq, c))
    zc = ds.isel(site=0).edi.z_as_comp()
    comps = ["zxx", "zxy", "zyx", "zyy"]
    assert zc.dims == ("c", "freq")
    assert list(zc.coords["c"].data) == comps
    assert zc.shape[1] == ds.sizes["freq"]

    # band selection (use inner range if possible)
    f = ds["freq"].values
    if f.size >= 3:
        fmin = float(np.nanmin(f) + 1e-9)
        fmax = float(np.nanmax(f) - 1e-9)
        slim = ds.edi.band(fmin=fmin, fmax=fmax)
        assert 0 < slim.sizes["freq"] <= ds.sizes["freq"]
    else:
        slim = ds.edi.band()
        assert slim.sizes["freq"] == ds.sizes["freq"]

    # spectra/timeseries helpers behave (presence-agnostic)
    assert isinstance(ds.edi.has_spectra(), bool)
    assert isinstance(ds.edi.spectra(), xr.Dataset)
    assert isinstance(ds.edi.has_timeseries(), bool)
    assert isinstance(ds.edi.timeseries(), xr.Dataset)


@pytest.mark.usefixtures("edi_spe_file")
def test_build_dataset_with_spectra(edi_spe_file: Path) -> None:
    ed = EDIFile(edi_spe_file)
    ds = build_dataset([ed], drop_empty=False)

    # Always at least 1 site (even if freq is empty)
    assert "site" in ds.coords
    assert ds.sizes["site"] >= 1

    # Only check spectra variables if present in this sample
    if ds.edi.has_spectra():
        assert "spec_vals" in ds
        assert "spec_len" in ds
        assert "freq" in ds["spec_vals"].coords

        sp = ds.edi.spectra()
        for k in ("spec_vals", "spec_len"):
            assert k in sp
    else:
        pytest.skip("No spectra variables present in this EDI sample.")


@pytest.mark.usefixtures("edi_imp_file", "edi_spe_file")
def test_xamixin_collection_meta(
    edi_imp_file: Path,
    edi_spe_file: Path,
) -> None:
    ed1 = EDIFile(edi_imp_file)
    ed2 = EDIFile(edi_spe_file)
    coll = _DummyColl([ed1, ed2])

    ds = coll.to_xarray(drop_empty=False)
    assert "site" in ds.coords
    # Number of sites equals number of stations exposed by accessor
    assert ds.sizes["site"] == len(ds.edi.stations)
    assert ds.sizes["site"] >= 1

    meta = coll.meta_table()
    # meta is a dataset with site coord and known columns
    assert "site" in meta.coords
    assert meta.sizes["site"] == ds.sizes["site"]
    for k in ("filename", "nfreq", "has_tip", "has_spec", "has_ts"):
        assert k in meta
