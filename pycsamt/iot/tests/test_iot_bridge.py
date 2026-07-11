"""Tests for :mod:`pycsamt.iot.bridge` — the bridge between IoT field
telemetry and the pyCSAMT data model (Z / EDI / Site) in both directions.
"""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.iot import (
    DeploymentConfig,
    FieldSession,
    StationConfig,
    build_edifile,
    deployment_from_edis,
    edi_survey_table,
    field_session_from_edis,
    field_session_to_edifiles,
    impedance_to_z,
    read_edi_survey,
    z_to_edi,
)
from pycsamt.z.z import Z

FREQ = np.logspace(4, 0, 10)


def _scalar_zxy(freq: np.ndarray = FREQ) -> np.ndarray:
    """A simple, well-behaved off-diagonal impedance for tests."""
    return (1 + 1j) * np.sqrt(freq)


# ---------------------------------------------------------------------------
# impedance_to_z
# ---------------------------------------------------------------------------
def test_impedance_to_z_scalar_offdiag():
    zxy = _scalar_zxy()
    z = impedance_to_z(zxy, FREQ, station="S01", method="amt")
    assert isinstance(z, Z)
    assert z.z.shape == (FREQ.size, 2, 2)
    assert np.allclose(z.z[:, 0, 1], zxy)
    assert np.allclose(z.z[:, 1, 0], -zxy)          # antisymmetric placement
    assert np.allclose(z.z[:, 0, 0], 0)
    assert z.name == "S01"
    assert z.resistivity.shape == (FREQ.size, 2, 2)


@pytest.mark.parametrize(
    "component,idx",
    [("xy", (0, 1)), ("yx", (1, 0))],
)
def test_impedance_to_z_scalar_component_placement(component, idx):
    zxy = _scalar_zxy()
    z = impedance_to_z(zxy, FREQ, component=component)
    assert np.allclose(z.z[:, idx[0], idx[1]], zxy)
    other = (1, 0) if idx == (0, 1) else (0, 1)
    assert np.allclose(z.z[:, other[0], other[1]], 0)


def test_impedance_to_z_full_tensor():
    zxy = _scalar_zxy()
    tensor = np.zeros((FREQ.size, 2, 2), dtype=complex)
    tensor[:, 0, 1] = zxy
    tensor[:, 1, 0] = -zxy
    z = impedance_to_z(tensor, FREQ)
    assert np.allclose(z.z, tensor)


def test_impedance_to_z_window_aggregation_sets_error():
    rng = np.random.default_rng(0)
    zxy = _scalar_zxy()
    windows = zxy[None, :] * (1 + 0.02 * rng.standard_normal((8, FREQ.size)))
    z = impedance_to_z(windows, FREQ)
    # mean of windows is close to the noise-free curve
    assert np.allclose(z.z[:, 0, 1], zxy, rtol=0.05)
    # spread across windows becomes a non-zero absolute error
    assert z.z_err is not None
    assert np.all(z.z_err[:, 0, 1] > 0)


def test_impedance_to_z_window_tensor_4d():
    zxy = _scalar_zxy()
    base = np.zeros((FREQ.size, 2, 2), dtype=complex)
    base[:, 0, 1] = zxy
    base[:, 1, 0] = -zxy
    windows = np.stack([base, base, base], axis=0)
    z = impedance_to_z(windows, FREQ)
    assert np.allclose(z.z, base)
    assert z.z_err is not None                       # error grid allocated


def test_impedance_to_z_single_tensor_one_freq():
    tensor = np.array([[0 + 0j, 1 + 1j], [-1 - 1j, 0 + 0j]])
    z = impedance_to_z(tensor, [1.0])
    assert z.z.shape == (1, 2, 2)
    assert np.allclose(z.z[0], tensor)


def test_impedance_to_z_explicit_error_vector():
    zxy = _scalar_zxy()
    err = np.full(FREQ.size, 0.5)
    z = impedance_to_z(zxy, FREQ, impedance_err=err)
    assert np.allclose(z.z_err[:, 0, 1], 0.5)
    assert np.allclose(z.z_err[:, 1, 0], 0.5)


def test_impedance_to_z_median_aggregate():
    zxy = _scalar_zxy()
    windows = np.stack([zxy, zxy * 2, zxy * 3], axis=0)
    z = impedance_to_z(windows, FREQ, aggregate="median")
    assert np.allclose(z.z[:, 0, 1], zxy * 2)        # median of {1,2,3}x


def test_impedance_to_z_rejects_bad_freq():
    with pytest.raises(ValueError):
        impedance_to_z(_scalar_zxy(), np.zeros(FREQ.size))
    with pytest.raises(ValueError):
        impedance_to_z(_scalar_zxy(), [])


def test_impedance_to_z_rejects_shape_mismatch():
    with pytest.raises(ValueError):
        impedance_to_z(np.ones(5), FREQ)             # 1-D length mismatch
    with pytest.raises(ValueError):
        impedance_to_z(np.ones((FREQ.size, 3)), FREQ)  # not windows/tensor


def test_impedance_to_z_invalid_component():
    with pytest.raises(ValueError):
        impedance_to_z(_scalar_zxy(), FREQ, component="det")


# ---------------------------------------------------------------------------
# forward: Z -> EDI
# ---------------------------------------------------------------------------
def test_build_edifile_structure():
    z = impedance_to_z(_scalar_zxy(), FREQ, station="S01")
    ed = build_edifile(z, station="S01", lat=6.5, lon=3.4, elevation=120.0)
    assert ed.station == "S01"
    assert ed.n_freq == FREQ.size
    head = ed.get_section("head")
    assert head.lat == pytest.approx(6.5)
    assert head.lon == pytest.approx(3.4)


def test_build_edifile_requires_z():
    with pytest.raises(TypeError):
        build_edifile(np.ones((FREQ.size, 2, 2)), station="S01")


def test_z_to_edi_roundtrip(tmp_path):
    from pycsamt.seg.edi import EDIFile

    z = impedance_to_z(_scalar_zxy(), FREQ, station="S01")
    path = z_to_edi(
        z, station="S01", lat=6.5, lon=3.4, elevation=120.0,
        savepath=str(tmp_path), method="amt",
    )
    assert path.endswith(".edi")
    back = EDIFile(path)
    assert back.station == "S01"
    assert back.n_freq == FREQ.size
    assert np.allclose(
        np.sort_complex(back.Z.z[:, 0, 1]),
        np.sort_complex(z.z[:, 0, 1]),
    )


def test_field_session_to_edifiles_uses_geometry(tmp_path):
    from pycsamt.seg.edi import EDIFile

    z1 = impedance_to_z(_scalar_zxy(), FREQ, station="S01")
    z2 = impedance_to_z(_scalar_zxy() * 2, FREQ, station="S02")
    session = FieldSession("SURV1", method="amt")
    session.add_station(
        StationConfig("S01", lat=6.5, lon=3.4, elevation=120.0)
    )
    session.add_station(
        StationConfig("S02", lat=6.6, lon=3.5, elevation=130.0)
    )
    paths = field_session_to_edifiles(
        session, {"S01": z1, "S02": z2}, write=True, savepath=str(tmp_path)
    )
    assert set(paths) == {"S01", "S02"}
    back = EDIFile(paths["S02"])
    assert back.get_section("head").lat == pytest.approx(6.6)
    assert back.get_section("head").lon == pytest.approx(3.5)


def test_field_session_to_edifiles_returns_objects():
    session = FieldSession("SURV1")
    session.add_station(StationConfig("S01", lat=6.5, lon=3.4))
    # value given as an (impedance_array, freq) pair rather than a Z
    out = field_session_to_edifiles(session, {"S01": (_scalar_zxy(), FREQ)})
    assert out["S01"].station == "S01"          # in-memory EDIFile, not path


# ---------------------------------------------------------------------------
# reverse: EDI survey -> acquisition planning
# ---------------------------------------------------------------------------
@pytest.fixture()
def survey_dir(tmp_path):
    """Write a small two-station EDI survey and return its directory."""
    for sid, lat, lon in (("A01", 6.5, 3.4), ("A02", 6.6, 3.5)):
        z = impedance_to_z(_scalar_zxy(), FREQ, station=sid)
        z_to_edi(z, station=sid, lat=lat, lon=lon, elevation=100.0,
                 savepath=str(tmp_path))
    return str(tmp_path)


def test_read_edi_survey(survey_dir):
    records = read_edi_survey(survey_dir)
    assert len(records) == 2
    stations = {r["station"] for r in records}
    assert stations == {"A01", "A02"}
    rec = next(r for r in records if r["station"] == "A01")
    assert rec["lat"] == pytest.approx(6.5)
    assert rec["n_freq"] == FREQ.size
    assert rec["f_max_hz"] == pytest.approx(FREQ.max())
    assert rec["channels"]


def test_edi_survey_table(survey_dir):
    tbl = edi_survey_table(survey_dir)
    assert tbl.shape[0] == 2
    for col in ("station", "lat", "lon", "f_min_hz", "f_max_hz", "n_channels"):
        assert col in tbl.columns


def test_field_session_from_edis(survey_dir):
    session = field_session_from_edis(
        survey_dir, survey_id="REOCCUPY", method="amt", operator="crew"
    )
    assert isinstance(session, FieldSession)
    assert session.n_stations == 2
    assert session.n_devices == 2
    site = next(s for s in session.to_sites() if s.station_id == "A01")
    assert site.coords == (pytest.approx(6.5), pytest.approx(3.4),
                           pytest.approx(100.0))
    # a sample-rate hint is derived from the highest recovered frequency
    device = session.devices["A01-node"]
    assert device.sample_rate_hz is not None
    assert device.sample_rate_hz > FREQ.max()


def test_deployment_from_edis(survey_dir):
    dep = deployment_from_edis(
        survey_dir, survey_id="WILLY", capabilities=["telemetry"],
    )
    assert isinstance(dep, DeploymentConfig)
    assert dep.n_devices == 2
    assert dep.has_capability("telemetry")


def test_reverse_accepts_single_file(survey_dir):
    import os

    one = sorted(
        os.path.join(survey_dir, f)
        for f in os.listdir(survey_dir) if f.endswith(".edi")
    )[0]
    records = read_edi_survey(one)
    assert len(records) == 1


# ---------------------------------------------------------------------------
# z_to_site (optional geospatial stack)
# ---------------------------------------------------------------------------
def test_z_to_site_optional():
    pytest.importorskip("pyproj")
    from pycsamt.iot import z_to_site

    z = impedance_to_z(_scalar_zxy(), FREQ, station="S01")
    site = z_to_site(z, station="S01", lat=6.5, lon=3.4)
    assert site is not None
