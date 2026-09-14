# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later

from __future__ import annotations

from collections.abc import Iterable, Iterator
from pathlib import Path
from types import SimpleNamespace

import matplotlib.pyplot as plt
import numpy as np
import pytest
import xarray as xr

from pycsamt.seg.collection import EDICollection
from pycsamt.seg.edi import EDIFile
from pycsamt.seg.xa import (
    XAMixin,
    _get_tensor_or_zeros,
    _get_tipper_or_zeros,
    _meta_from_edi,
    _pad2d,
    _site_location,
    _spec_pack,
    _ts_pack,
    build_dataset,
)

# Mark all tests in this file as requiring the 'xarray' package
pytestmark = pytest.mark.requires_xarray


class _DummyColl(XAMixin):
    """Minimal iterable collection to test XAMixin."""

    def __init__(self, items: Iterable[EDIFile]) -> None:
        self._items: list[EDIFile] = list(items)

    def __iter__(self) -> Iterator[EDIFile]:
        return iter(self._items)


class _FakeEDI:
    """Minimal duck-typed EDI-like object for direct helper unit tests."""

    def __init__(self, sections: dict[str, object] | None = None):
        self._sections = sections or {}

    def get_section(self, name: str):
        return self._sections.get(name)


# ─────────────────────────────────────────────────────────────────────────
# Pure helper functions, exercised directly
# ─────────────────────────────────────────────────────────────────────────


def test_get_tensor_or_zeros_returns_zeros_when_missing():
    out = _get_tensor_or_zeros(None, "z", 3, dtype=complex)
    assert out.shape == (3, 2, 2)
    assert np.all(out == 0)


def test_get_tensor_or_zeros_returns_zeros_for_wrong_shape():
    obj = SimpleNamespace(z=np.ones((4, 2, 2)))
    out = _get_tensor_or_zeros(obj, "z", 3, dtype=float)
    assert out.shape == (3, 2, 2)
    assert np.all(out == 0)


def test_get_tensor_or_zeros_returns_actual_values():
    obj = SimpleNamespace(z=np.ones((3, 2, 2)) * 5.0)
    out = _get_tensor_or_zeros(obj, "z", 3, dtype=float)
    assert np.all(out == 5.0)


def test_get_tipper_or_zeros_none_returns_zeros():
    out = _get_tipper_or_zeros(None, "tipper", 2, dtype=complex)
    assert out.shape == (2, 2)
    assert np.all(out == 0)


def test_get_tipper_or_zeros_standard_3d_shape():
    obj = SimpleNamespace(tipper=np.ones((2, 1, 2)) * 3.0)
    out = _get_tipper_or_zeros(obj, "tipper", 2, dtype=float)
    assert out.shape == (2, 2)
    assert np.all(out == 3.0)


def test_get_tipper_or_zeros_already_2d_shape():
    obj = SimpleNamespace(tipper=np.ones((2, 2)) * 7.0)
    out = _get_tipper_or_zeros(obj, "tipper", 2, dtype=float)
    assert out.shape == (2, 2)
    assert np.all(out == 7.0)


def test_get_tipper_or_zeros_unrecognized_shape_falls_back_to_zeros():
    obj = SimpleNamespace(tipper=np.ones((5, 5)))
    out = _get_tipper_or_zeros(obj, "tipper", 2, dtype=float)
    assert out.shape == (2, 2)
    assert np.all(out == 0)


def test_pad2d_pads_ragged_sequences():
    out = _pad2d([np.array([1.0, 2.0]), np.array([3.0])])
    assert out.shape == (2, 2)
    assert out[0].tolist() == [1.0, 2.0]
    assert out[1, 0] == 3.0
    assert np.isnan(out[1, 1])


def test_pad2d_empty_input():
    out = _pad2d([])
    assert out.shape == (0, 0)


def test_site_location_falls_back_to_definemeas():
    head = SimpleNamespace(lat=None, long=None, elev=None)
    dm = SimpleNamespace(reflat=1.0, reflong=2.0, refelev=3.0)
    ed = _FakeEDI({"head": head, "definemeas": dm})
    assert _site_location(ed) == (1.0, 2.0, 3.0)


def test_site_location_prefers_head_values():
    head = SimpleNamespace(lat=10.0, long=20.0, elev=30.0)
    ed = _FakeEDI({"head": head})
    assert _site_location(ed) == (10.0, 20.0, 30.0)


def test_site_location_no_sections_returns_all_none():
    ed = _FakeEDI({})
    assert _site_location(ed) == (None, None, None)


def test_meta_from_edi_uses_helpers():
    ed = SimpleNamespace(
        path=Path("/tmp/S01.edi"),
        processingsoftware="pyCSAMT",
        station="S01",
        has_tipper=True,
        n_freq=5,
        get_section=lambda name: None,
    )
    meta = _meta_from_edi(ed)
    assert meta["site"] == "S01"
    assert meta["filename"] == "S01.edi"
    assert meta["software"] == "pyCSAMT"
    assert meta["has_spec"] is False
    assert meta["has_ts"] is False


# ─────────────────────────────────────────────────────────────────────────
# _spec_pack / _ts_pack edge cases
# ─────────────────────────────────────────────────────────────────────────


def test_spec_pack_no_spectra_section_returns_empty():
    ed = _FakeEDI({})
    assert _spec_pack(ed) == {}


def test_spec_pack_direct_attributes():
    spec = SimpleNamespace(
        freq=[1.0, 2.0],
        values=[[1.0, 2.0], [3.0]],
        bw=[0.1, 0.2],
        avgt=[1.0, 1.0],
        rotspec=[0.0, 0.0],
    )
    ed = _FakeEDI({"spectra": spec})
    out = _spec_pack(ed)
    assert "spec_vals" in out
    assert "spec_bw" in out
    assert out["spec_vals"].shape == (2, 2)


def test_spec_pack_invalid_freq_list_returns_empty():
    spec = SimpleNamespace(freq=["not-a-number"], values=[[1.0]])
    ed = _FakeEDI({"spectra": spec})
    assert _spec_pack(ed) == {}


def test_spec_pack_empty_freq_returns_empty():
    spec = SimpleNamespace(freq=[], values=[])
    ed = _FakeEDI({"spectra": spec})
    assert _spec_pack(ed) == {}


def test_spec_pack_invalid_values_returns_empty():
    spec = SimpleNamespace(freq=[1.0], values=[object()])
    ed = _FakeEDI({"spectra": spec})
    assert _spec_pack(ed) == {}


def test_spec_pack_to_io_fallback_raises_returns_empty():
    class BadSpec:
        freq = None
        values = None

        def to_io(self):
            raise RuntimeError("boom")

    ed = _FakeEDI({"spectra": BadSpec()})
    assert _spec_pack(ed) == {}


def test_spec_pack_mismatched_optional_metadata_is_skipped():
    spec = SimpleNamespace(
        freq=[1.0, 2.0],
        values=[[1.0], [2.0]],
        bw=[0.1],  # wrong length -> skipped
        avgt=None,
        rotspec=None,
    )
    ed = _FakeEDI({"spectra": spec})
    out = _spec_pack(ed)
    assert "spec_bw" not in out


def test_spec_pack_optional_metadata_unconvertible_is_skipped():
    spec = SimpleNamespace(
        freq=[1.0, 2.0],
        values=[[1.0], [2.0]],
        bw=["not", "numeric"],  # same length but not float-convertible
        avgt=None,
        rotspec=None,
    )
    ed = _FakeEDI({"spectra": spec})
    out = _spec_pack(ed)
    assert "spec_bw" not in out


def test_ts_pack_no_timeseries_section_returns_empty():
    ed = _FakeEDI({})
    assert _ts_pack(ed) == {}


def test_ts_pack_channels_raises_returns_empty():
    class BadTS:
        def channels(self):
            raise RuntimeError("boom")

    ed = _FakeEDI({"timeseries": BadTS()})
    assert _ts_pack(ed) == {}


def test_ts_pack_empty_channels_returns_empty():
    ts = SimpleNamespace(channels=lambda: [])
    ed = _FakeEDI({"timeseries": ts})
    assert _ts_pack(ed) == {}


def test_ts_pack_builds_arrays():
    ts = SimpleNamespace(
        channels=lambda: ["HX", "HY"],
        get=lambda c: {"HX": [1.0, 2.0, 3.0], "HY": [4.0, 5.0]}[c],
        dt_map={"HX": 0.1, "HY": 0.2},
    )
    ed = _FakeEDI({"timeseries": ts})
    out = _ts_pack(ed)
    assert set(out) == {"ts", "time", "dt", "npts"}
    assert out["ts"].shape == (3, 2)
    assert out["npts"].values.tolist() == [3, 2]


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


# ─────────────────────────────────────────────────────────────────────────
# build_dataset: empty / drop_empty / error-skip branches
# ─────────────────────────────────────────────────────────────────────────


def test_build_dataset_no_items_returns_empty_dataset():
    ds = build_dataset([])
    assert ds.sizes.get("site", 0) == 0
    assert "freq" in ds.coords


def test_build_dataset_drops_empty_items_by_default():
    empty_ed = SimpleNamespace(
        station="EMPTY1",
        path=None,
        processingsoftware=None,
        Z=SimpleNamespace(freq=np.array([]), n_freq=0),
        Tip=SimpleNamespace(tipper=None),
        has_tipper=False,
        n_freq=0,
        get_section=lambda name: None,
    )
    ds = build_dataset([empty_ed])  # drop_empty=True (default)
    assert ds.sizes.get("site", 0) == 0


def test_build_dataset_skips_items_that_raise():
    class Broken:
        def __getattr__(self, name):
            raise RuntimeError("simulated failure")

    ds = build_dataset([Broken()])
    assert ds.sizes.get("site", 0) == 0


def test_build_dataset_falls_back_to_spectra_freq_when_z_empty() -> None:
    spec = SimpleNamespace(
        freq=np.array([1.0, 2.0, 3.0]),
        values=[[1.0], [2.0], [3.0]],
        bw=None,
        avgt=None,
        rotspec=None,
    )
    ed = SimpleNamespace(
        station="SPEC1",
        path=None,
        processingsoftware=None,
        Z=SimpleNamespace(freq=np.array([]), n_freq=0),
        Tip=SimpleNamespace(tipper=None),
        has_tipper=False,
        n_freq=0,
        get_section=lambda name: spec if name == "spectra" else None,
    )
    ds = build_dataset([ed], drop_empty=False)
    assert ds.sizes.get("site", 0) == 1
    assert ds.sizes["freq"] == 3


# ─────────────────────────────────────────────────────────────────────────
# XAMixin.meta_table on an empty collection
# ─────────────────────────────────────────────────────────────────────────


def test_meta_table_empty_collection_returns_empty_dataset():
    coll = _DummyColl([])
    meta = coll.meta_table()
    assert meta.sizes.get("site", 0) == 0


# ─────────────────────────────────────────────────────────────────────────
# EDIAcc: get() KeyError, band() branches, attrs/timeseries helpers,
# plot_apparent_resistivity edge cases
# ─────────────────────────────────────────────────────────────────────────


@pytest.mark.usefixtures("edi_imp_file")
def test_edi_accessor_get_unknown_site_raises_keyerror(
    edi_imp_file: Path,
) -> None:
    ed = EDIFile(edi_imp_file)
    ds = build_dataset([ed], drop_empty=False)
    with pytest.raises(KeyError):
        ds.edi.get("no-such-site")


@pytest.mark.usefixtures("edi_imp_file")
def test_edi_accessor_band_filters_by_range(edi_imp_file: Path) -> None:
    ed = EDIFile(edi_imp_file)
    ds = build_dataset([ed], drop_empty=False)
    f = ds["freq"].values
    assert f.size >= 2
    fmin, fmax = float(f.min()), float(f[len(f) // 2])
    limited = ds.edi.band(fmin=fmin, fmax=fmax)
    assert limited.sizes["freq"] <= ds.sizes["freq"]
    assert limited.sizes["freq"] > 0


@pytest.mark.usefixtures("edi_imp_file")
def test_edi_accessor_band_no_freq_coord_returns_dataset_unchanged(
    edi_imp_file: Path,
) -> None:
    ed = EDIFile(edi_imp_file)
    ds = build_dataset([ed], drop_empty=False)
    meta = xr.Dataset(coords={"site": ds["site"].values})
    result = meta.edi.band(fmin=1.0)
    assert result is meta


@pytest.mark.usefixtures("edi_imp_file")
def test_edi_accessor_attrs_and_timeseries_helpers(edi_imp_file: Path) -> None:
    ed = EDIFile(edi_imp_file)
    ds = build_dataset([ed], drop_empty=False)
    ds.attrs["note"] = "hello"
    assert ds.edi.attrs() == {"note": "hello"}
    assert ds.edi.has_timeseries() is False
    # no ts/time/dt/npts vars present -> returns dataset unchanged
    assert ds.edi.timeseries() is ds


@pytest.mark.usefixtures("edi_imp_file")
def test_plot_apparent_resistivity_invalid_component_is_skipped(
    edi_imp_file: Path,
) -> None:
    ed = EDIFile(edi_imp_file)
    ds = build_dataset([ed], drop_empty=False)
    site_name = ds.edi.stations[0]
    try:
        fig, axes = ds.edi.plot_apparent_resistivity(
            site=site_name, components=["xy", "bogus"],
        )
        assert isinstance(fig, plt.Figure)
    finally:
        plt.close("all")


@pytest.mark.usefixtures("edi_imp_file")
def test_plot_apparent_resistivity_custom_grid_and_phase_mod_and_savefig(
    edi_imp_file: Path, tmp_path: Path,
) -> None:
    ed = EDIFile(edi_imp_file)
    ds = build_dataset([ed], drop_empty=False)
    site_name = ds.edi.stations[0]
    out_file = tmp_path / "plot.png"
    try:
        fig, axes = ds.edi.plot_apparent_resistivity(
            site=site_name,
            phase_mod=90,
            grid_props={"color": "red"},
            savefig=str(out_file),
        )
        assert isinstance(fig, plt.Figure)
        assert out_file.exists()
    finally:
        plt.close("all")
