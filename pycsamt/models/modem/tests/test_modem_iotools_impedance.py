"""Tests for pycsamt.models.modem.iotools.impedance."""

import numpy as np
import pytest

from pycsamt.models.modem.iotools.impedance import (
    ImpedanceFile,
    ZBlock,
    convert_z2d,
    convert_z3d,
    read_z2d_old,
    read_z3d_old,
    write_z2d_list,
    write_z2d_old,
    write_z3d_list,
    write_z3d_old,
)

# ---------------------------------------------------------------------------
# Synthetic data builders
# ---------------------------------------------------------------------------


def _z3d_imp(n_sites=3, n_comp=4, periods=(10.0, 100.0)):
    """Build an ImpedanceFile with simple 3D data."""
    rng = np.random.default_rng(7)
    comps = ["ZXX", "ZXY", "ZYX", "ZYY"][:n_comp]
    blocks = []
    for T in periods:
        Z = rng.standard_normal((n_sites, n_comp)) + 1j * rng.standard_normal(
            (n_sites, n_comp)
        )
        Zerr = np.abs(rng.standard_normal((n_sites, n_comp))) * 0.1
        loc = rng.uniform(-5000, 5000, (n_sites, 3))
        loc[:, 2] = 0.0  # surface
        blocks.append(
            ZBlock(
                period=T,
                site_names=[f"S{k:02d}" for k in range(n_sites)],
                site_loc=loc,
                comp_names=comps,
                Z=Z,
                Zerr=Zerr,
            )
        )
    return ImpedanceFile(
        description="Synthetic test data",
        units="[V/m]/[T]",
        sign=-1,
        blocks=blocks,
        origin=(100.0, 200.0, 0.0),
    )


def _z2d_imp(n_sites=4, periods_te=(10.0,), periods_tm=(10.0,)):
    """Build an ImpedanceFile with 2D TE+TM data."""
    rng = np.random.default_rng(13)
    blocks = []
    for T in periods_te:
        Z = rng.standard_normal((n_sites, 1)) + 1j * rng.standard_normal(
            (n_sites, 1)
        )
        Zerr = np.abs(rng.standard_normal((n_sites, 1))) * 0.05
        ys = np.linspace(-10000, 10000, n_sites)
        loc = np.column_stack([np.zeros(n_sites), ys, np.zeros(n_sites)])
        blocks.append(
            ZBlock(
                period=T,
                site_names=[f"S{k:02d}" for k in range(n_sites)],
                site_loc=loc,
                comp_names=["ZTE"],
                Z=Z,
                Zerr=Zerr,
                mode="TE",
            )
        )
    for T in periods_tm:
        Z = rng.standard_normal((n_sites, 1)) + 1j * rng.standard_normal(
            (n_sites, 1)
        )
        Zerr = np.abs(rng.standard_normal((n_sites, 1))) * 0.05
        ys = np.linspace(-10000, 10000, n_sites)
        loc = np.column_stack([np.zeros(n_sites), ys, np.zeros(n_sites)])
        blocks.append(
            ZBlock(
                period=T,
                site_names=[f"S{k:02d}" for k in range(n_sites)],
                site_loc=loc,
                comp_names=["ZTM"],
                Z=Z,
                Zerr=Zerr,
                mode="TM",
            )
        )
    return ImpedanceFile(
        description="Synthetic 2D test",
        units="[V/m]/[T]",
        sign=-1,
        blocks=blocks,
    )


# ===========================================================================
# ZBlock / ImpedanceFile dataclass
# ===========================================================================


def test_zblock_n_sites():
    blk = _z3d_imp().blocks[0]
    assert blk.n_sites == 3


def test_zblock_n_comp():
    blk = _z3d_imp().blocks[0]
    assert blk.n_comp == 4


def test_impedancefile_periods():
    imp = _z3d_imp(periods=(1.0, 10.0, 100.0))
    np.testing.assert_array_equal(imp.periods, [1.0, 10.0, 100.0])


def test_impedancefile_n_periods():
    assert _z3d_imp(periods=(5.0, 50.0)).n_periods == 2


# ===========================================================================
# Old 3D format roundtrip
# ===========================================================================


class TestZ3dOld:
    def test_write_creates_file(self, tmp_path):
        imp = _z3d_imp()
        out = write_z3d_old(imp, tmp_path / "z3d.imp")
        assert out.exists()

    def test_header_lines(self, tmp_path):
        imp = _z3d_imp()
        out = write_z3d_old(imp, tmp_path / "z3d.imp")
        text = out.read_text()
        assert "Description:" in text
        assert "Units:" in text
        assert "Sign convention:" in text

    def test_roundtrip_n_periods(self, tmp_path):
        imp = _z3d_imp(periods=(1.0, 10.0, 100.0))
        out = write_z3d_old(imp, tmp_path / "z3d.imp")
        back = read_z3d_old(out)
        assert back.n_periods == 3

    def test_roundtrip_n_sites(self, tmp_path):
        imp = _z3d_imp(n_sites=5)
        write_z3d_old(imp, tmp_path / "z3d.imp")
        back = read_z3d_old(tmp_path / "z3d.imp")
        assert back.blocks[0].n_sites == 5

    def test_roundtrip_period_values(self, tmp_path):
        imp = _z3d_imp(periods=(7.3, 73.0))
        write_z3d_old(imp, tmp_path / "z3d.imp")
        back = read_z3d_old(tmp_path / "z3d.imp")
        np.testing.assert_allclose(back.periods, [7.3, 73.0], rtol=1e-5)

    def test_roundtrip_z_real(self, tmp_path):
        imp = _z3d_imp(n_sites=3, n_comp=4, periods=(10.0,))
        write_z3d_old(imp, tmp_path / "z3d.imp")
        back = read_z3d_old(tmp_path / "z3d.imp")
        np.testing.assert_allclose(
            back.blocks[0].Z.real, imp.blocks[0].Z.real, rtol=1e-5
        )

    def test_roundtrip_z_imag(self, tmp_path):
        imp = _z3d_imp(n_sites=3, n_comp=4, periods=(10.0,))
        write_z3d_old(imp, tmp_path / "z3d.imp")
        back = read_z3d_old(tmp_path / "z3d.imp")
        np.testing.assert_allclose(
            back.blocks[0].Z.imag, imp.blocks[0].Z.imag, rtol=1e-5
        )

    def test_roundtrip_zerr(self, tmp_path):
        imp = _z3d_imp(n_sites=3, n_comp=4, periods=(10.0,))
        write_z3d_old(imp, tmp_path / "z3d.imp")
        back = read_z3d_old(tmp_path / "z3d.imp")
        np.testing.assert_allclose(
            back.blocks[0].Zerr, imp.blocks[0].Zerr, rtol=1e-5
        )

    def test_roundtrip_site_loc(self, tmp_path):
        imp = _z3d_imp(n_sites=3, periods=(10.0,))
        write_z3d_old(imp, tmp_path / "z3d.imp")
        back = read_z3d_old(tmp_path / "z3d.imp")
        np.testing.assert_allclose(
            back.blocks[0].site_loc, imp.blocks[0].site_loc, atol=0.01
        )

    def test_roundtrip_site_names(self, tmp_path):
        imp = _z3d_imp(n_sites=3, periods=(10.0,))
        write_z3d_old(imp, tmp_path / "z3d.imp")
        back = read_z3d_old(tmp_path / "z3d.imp")
        assert back.blocks[0].site_names == imp.blocks[0].site_names

    def test_roundtrip_comp_names(self, tmp_path):
        imp = _z3d_imp(n_comp=4, periods=(10.0,))
        write_z3d_old(imp, tmp_path / "z3d.imp")
        back = read_z3d_old(tmp_path / "z3d.imp")
        assert back.blocks[0].comp_names == ["ZXX", "ZXY", "ZYX", "ZYY"]

    def test_roundtrip_description(self, tmp_path):
        imp = _z3d_imp()
        write_z3d_old(imp, tmp_path / "z3d.imp")
        back = read_z3d_old(tmp_path / "z3d.imp")
        assert back.description == imp.description

    def test_roundtrip_sign(self, tmp_path):
        imp = _z3d_imp()
        imp.sign = +1
        write_z3d_old(imp, tmp_path / "z3d.imp")
        back = read_z3d_old(tmp_path / "z3d.imp")
        assert back.sign == +1

    def test_off_diagonal_comp_names(self, tmp_path):
        """2-component (off-diagonal) roundtrip."""
        imp = _z3d_imp(n_comp=2, periods=(10.0,))
        write_z3d_old(imp, tmp_path / "z3d.imp")
        back = read_z3d_old(tmp_path / "z3d.imp")
        assert back.blocks[0].n_comp == 2

    def test_missing_file_raises(self):
        with pytest.raises(FileNotFoundError):
            read_z3d_old("/nonexistent/z3d.imp")


# ===========================================================================
# Old 2D format roundtrip
# ===========================================================================


class TestZ2dOld:
    def test_write_creates_file(self, tmp_path):
        imp = _z2d_imp()
        assert write_z2d_old(imp, tmp_path / "z2d.imp").exists()

    def test_roundtrip_n_periods(self, tmp_path):
        imp = _z2d_imp(periods_te=(1.0, 10.0), periods_tm=(1.0, 10.0))
        write_z2d_old(imp, tmp_path / "z2d.imp")
        back = read_z2d_old(tmp_path / "z2d.imp")
        assert back.n_periods == 4

    def test_roundtrip_mode(self, tmp_path):
        imp = _z2d_imp()
        write_z2d_old(imp, tmp_path / "z2d.imp")
        back = read_z2d_old(tmp_path / "z2d.imp")
        modes = [b.mode for b in back.blocks]
        assert "TE" in modes and "TM" in modes

    def test_roundtrip_z_values(self, tmp_path):
        imp = _z2d_imp(n_sites=4, periods_te=(5.0,), periods_tm=(5.0,))
        write_z2d_old(imp, tmp_path / "z2d.imp")
        back = read_z2d_old(tmp_path / "z2d.imp")
        np.testing.assert_allclose(
            back.blocks[0].Z.real, imp.blocks[0].Z.real, rtol=1e-5
        )

    def test_roundtrip_zerr(self, tmp_path):
        imp = _z2d_imp(n_sites=4, periods_te=(5.0,), periods_tm=(5.0,))
        write_z2d_old(imp, tmp_path / "z2d.imp")
        back = read_z2d_old(tmp_path / "z2d.imp")
        np.testing.assert_allclose(
            back.blocks[0].Zerr, imp.blocks[0].Zerr, rtol=1e-5
        )

    def test_missing_file_raises(self):
        with pytest.raises(FileNotFoundError):
            read_z2d_old("/nonexistent/z2d.imp")


# ===========================================================================
# Current list format writers
# ===========================================================================


class TestListFormat:
    def test_z3d_list_creates_file(self, tmp_path):
        imp = _z3d_imp()
        assert write_z3d_list(imp, tmp_path / "out.dat").exists()

    def test_z3d_list_header_hash(self, tmp_path):
        imp = _z3d_imp()
        out = write_z3d_list(imp, tmp_path / "out.dat")
        lines = out.read_text().splitlines()
        assert lines[0].startswith("#")

    def test_z3d_list_data_type(self, tmp_path):
        imp = _z3d_imp(n_comp=4)
        out = write_z3d_list(imp, tmp_path / "out.dat")
        assert "Full_Impedance" in out.read_text()

    def test_z3d_list_off_diagonal(self, tmp_path):
        imp = _z3d_imp(n_comp=2)
        out = write_z3d_list(imp, tmp_path / "out.dat")
        assert "Off_Diagonal_Impedance" in out.read_text()

    def test_z3d_list_data_lines_count(self, tmp_path):
        imp = _z3d_imp(n_sites=3, n_comp=4, periods=(10.0,))
        out = write_z3d_list(imp, tmp_path / "out.dat")
        data_lines = [
            ln
            for ln in out.read_text().splitlines()
            if ln and not ln.startswith("#") and not ln.startswith(">")
        ]
        # 3 sites × 4 components = 12 data lines
        assert len(data_lines) == 12

    def test_z2d_list_creates_file(self, tmp_path):
        imp = _z2d_imp()
        assert write_z2d_list(imp, tmp_path / "out.dat").exists()

    def test_z2d_list_te_tm_blocks(self, tmp_path):
        imp = _z2d_imp()
        out = write_z2d_list(imp, tmp_path / "out.dat")
        text = out.read_text()
        assert "TE_Impedance" in text
        assert "TM_Impedance" in text


# ===========================================================================
# Convert functions
# ===========================================================================


class TestConvert:
    def test_convert_z3d(self, tmp_path):
        imp = _z3d_imp()
        old = write_z3d_old(imp, tmp_path / "z3d.imp")
        new = convert_z3d(old, tmp_path / "z3d.dat")
        assert new.exists()
        assert "Full_Impedance" in new.read_text()

    def test_convert_z2d(self, tmp_path):
        imp = _z2d_imp()
        old = write_z2d_old(imp, tmp_path / "z2d.imp")
        new = convert_z2d(old, tmp_path / "z2d.dat")
        assert new.exists()
        assert "TE_Impedance" in new.read_text()

    def test_convert_z3d_preserves_periods(self, tmp_path):
        imp = _z3d_imp(periods=(7.0, 70.0, 700.0))
        old = write_z3d_old(imp, tmp_path / "z3d.imp")
        new = convert_z3d(old, tmp_path / "z3d.dat")
        text = new.read_text()
        for T in [7.0, 70.0, 700.0]:
            # period appears as %12.6E in each data line
            assert any(f"{T:.6E}" in ln for ln in text.splitlines()), (
                f"Period {T} not found in output"
            )
