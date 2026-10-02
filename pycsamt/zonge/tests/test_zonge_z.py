# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later
"""
Pytest suite for the Z (impedance) component.
"""

import numpy as np
import pandas as pd
import pytest

from pycsamt.constants import MU_0, PI
from pycsamt.exceptions import AvgDataError
from pycsamt.zonge.z import Z


@pytest.fixture
def sample_avg_data() -> pd.DataFrame:
    """
    Provides a sample tidy DataFrame mimicking a parsed AVG file.
    """
    data = {
        "station": [100, 100, 100, 100, 200, 200],
        "freq": [1024, 1024, 512, 512, 1024, 1024],
        "comp": ["ExHy", "EyHx", "ExHy", "EyHx", "ExHy", "EyHx"],
        "rho": [50.0, 150.0, 45.0, 140.0, 60.0, 160.0],
        "phase": [1000, 1200, 950, 1150, 1050, 1250],  # mrad
        "pc_rho": [5.0, 6.0, 4.5, 5.5, 7.0, 8.0],  # %
        "sphz": [20, 25, 18, 22, 28, 30],  # mrad
    }
    return pd.DataFrame(data)


class TestZ:
    """Test suite for the Z impedance component."""

    def test_read_success(self, sample_avg_data):
        """
        Verify successful reading and frame setup.
        """
        z_comp = Z()
        z_comp.read(sample_avg_data)
        assert not z_comp.frame.empty
        expected_cols = {
            "station",
            "freq",
            "comp",
            "rho",
            "phase",
            "pc_rho",
            "s_phz",
        }
        assert expected_cols.issubset(z_comp.frame.columns)
        assert z_comp.frame.shape[0] == 6

    def test_read_missing_required_columns(self):
        """
        Test that AvgDataError is raised for missing columns.
        """
        df_bad = pd.DataFrame({"freq": [100], "comp": ["ExHy"]})
        z_comp = Z()
        with pytest.raises(AvgDataError, match="missing required columns"):
            z_comp.read(df_bad)

    def test_z_properties(self, sample_avg_data):
        """
        Test calculation of z, z_real, and z_imag properties.
        """
        z_comp = Z()
        z_comp.read(sample_avg_data)

        # Manual calculation for the first row
        rho = 50.0
        phase_mrad = 1000
        freq = 1024
        omega = 2 * PI * freq
        phase_rad = phase_mrad * 1e-3
        z_mag_exp = np.sqrt(rho * omega * MU_0)
        z_exp = z_mag_exp * np.exp(1j * phase_rad)

        z_calc = z_comp.z.iloc[0]
        z_real_calc = z_comp.z_real.iloc[0]
        z_imag_calc = z_comp.z_imag.iloc[0]

        assert np.isclose(z_calc, z_exp)
        assert np.isclose(z_real_calc, z_exp.real)
        assert np.isclose(z_imag_calc, z_exp.imag)

    def test_z_err_property(self, sample_avg_data):
        """
        Test calculation of propagated error z_err.
        """
        z_comp = Z()
        z_comp.read(sample_avg_data)

        # Manual calculation for the first row
        rho = 50.0
        pc_rho = 5.0
        sphz_mrad = 20
        freq = 1024
        omega = 2 * PI * freq
        z_mag_exp = np.sqrt(rho * omega * MU_0)
        rel_err_rho = pc_rho / 100.0
        dphi_rad = sphz_mrad * 1e-3

        term_rho_sq = (0.5 * rel_err_rho) ** 2
        term_phi_sq = dphi_rad**2
        z_err_exp = z_mag_exp * np.sqrt(term_rho_sq + term_phi_sq)

        z_err_calc = z_comp.z_err.iloc[0]
        assert np.isclose(z_err_calc, z_err_exp)

    def test_component_properties(self, sample_avg_data):
        """
        Test component-specific properties like z_xy and z_xy_err.
        """
        z_comp = Z()
        z_comp.read(sample_avg_data)

        # First, verify that the 'comp' column was correctly
        # normalized to uppercase tokens during the read process.
        # This is the likely root cause of the original failure.
        expected_comps = ["EXHY", "EYHX", "EXHY", "EYHX", "EXHY", "EYHX"]
        assert (
            z_comp.frame["comp"].tolist() == expected_comps
        ), "Component normalization in Z.read() failed."

        # Now, test the properties which depend on this
        # Z_xy corresponds to EXHY
        z_xy = z_comp.z_xy
        z_xy_err = z_comp.z_xy_err

        assert len(z_xy) == 3
        assert len(z_xy_err) == 3

        # Verify first value corresponds to the first ExHy row
        assert np.isclose(z_xy.iloc[0], z_comp.z.iloc[0])
        assert np.isclose(z_xy_err.iloc[0], z_comp.z_err.iloc[0])

        # Z_yx corresponds to EyHx
        z_yx = z_comp.z_yx
        assert len(z_yx) == 3
        assert np.isclose(z_yx.iloc[0], z_comp.z.iloc[1])

    def test_to_tensor_single_station(self, sample_avg_data):
        """
        Test to_tensor for a single station.
        """
        z_comp = Z()
        z_comp.read(sample_avg_data)

        T, freqs, stations = z_comp.to_tensor(station=100)

        assert T.shape == (2, 2, 2)  # (n_freq, 2, 2)
        assert stations.size == 0
        assert np.allclose(freqs, [512, 1024])

        # Check values for freq=1024 (index 1)
        z_xy_val = z_comp.z.iloc[0]
        z_yx_val = z_comp.z.iloc[1]

        assert np.isclose(T[1, 0, 1], z_xy_val)  # ExHy -> (0, 1)
        assert np.isclose(T[1, 1, 0], z_yx_val)  # EyHx -> (1, 0)
        assert np.isnan(T[1, 0, 0])  # Zxx is NaN

    def test_to_tensor_multi_station(self, sample_avg_data):
        """
        Test to_tensor for multiple stations.
        """
        z_comp = Z()
        z_comp.read(sample_avg_data)

        T, freqs, stations = z_comp.to_tensor()

        assert T.shape == (2, 2, 2, 2)  # (n_st, n_freq, 2, 2)
        assert np.allclose(stations, [100, 200])
        assert np.allclose(freqs, [512, 1024])

        # Check station 200 (index 1), freq 1024 (index 1)
        z_xy_st200 = z_comp.z.iloc[4]
        assert np.isclose(T[1, 1, 0, 1], z_xy_st200)

    def test_to_xarray(self, sample_avg_data):
        """
        Test conversion to an xarray.DataArray.
        """
        z_comp = Z()
        z_comp.read(sample_avg_data)

        da = z_comp.to_xarray()

        assert da.name == "z"
        assert da.dims == ("station", "freq", "e", "h")
        assert da.shape == (2, 2, 2, 2)
        assert np.allclose(da.coords["station"], [100, 200])

        # Check a value
        val = da.sel(station=100, freq=1024, e="Ex", h="Hy").item()
        assert np.isclose(val, z_comp.z.iloc[0])


class TestZCoverage:
    """Additional coverage for edge cases and less-used branches."""

    def test_unread_instance_empty_properties(self):
        z_comp = Z()
        assert z_comp.z.empty
        assert z_comp.z_real.empty
        assert z_comp.z_imag.empty
        assert z_comp.z_err.empty
        assert z_comp.z_xy.empty
        assert z_comp.z_xy_err.empty

    def test_read_verbose_logs_missing_optional_columns(self, sample_avg_data):
        df = sample_avg_data.drop(columns=["pc_rho"])
        z_comp = Z(verbose=True)
        z_comp.read(df)
        assert "pc_rho" in z_comp.frame.columns
        assert z_comp.frame["pc_rho"].isna().all()

    def test_read_fills_missing_station_and_comp(self):
        df = pd.DataFrame(
            {"rho": [10.0], "phase": [500.0], "freq": [1024.0]}
        )
        z_comp = Z()
        z_comp.read(df)
        assert "station" in z_comp.frame.columns
        assert z_comp.frame["comp"].iloc[0] == "EXHY"

    def test_z_err_without_any_error_columns_returns_all_nan(self):
        z_comp = Z()
        z_comp._frame = pd.DataFrame(
            {
                "station": [1.0],
                "freq": [1024.0],
                "comp": ["EXHY"],
                "rho": [10.0],
            }
        )
        err = z_comp.z_err
        assert len(err) == 1
        assert err.isna().all()

    def test_z_err_with_only_rho_error(self):
        z_comp = Z()
        z_comp._frame = pd.DataFrame(
            {
                "station": [1.0],
                "freq": [1024.0],
                "comp": ["EXHY"],
                "rho": [10.0],
                "pc_rho": [5.0],
            }
        )
        err = z_comp.z_err
        assert not err.empty
        assert np.isfinite(err.iloc[0])

    def test_z_err_with_only_phase_error(self):
        z_comp = Z()
        z_comp._frame = pd.DataFrame(
            {
                "station": [1.0],
                "freq": [1024.0],
                "comp": ["EXHY"],
                "rho": [10.0],
                "s_phz": [20.0],
            }
        )
        err = z_comp.z_err
        assert not err.empty
        assert np.isfinite(err.iloc[0])

    def test_all_component_err_properties(self, sample_avg_data):
        z_comp = Z()
        z_comp.read(sample_avg_data)
        # exercise every component accessor at least once
        assert z_comp.z_xx.empty
        assert z_comp.z_yy.empty
        assert z_comp.z_xx_err.empty
        assert z_comp.z_yy_err.empty
        assert len(z_comp.z_yx_err) == 3

    def test_to_tensor_var_z_real_z_imag_z_err(self, sample_avg_data):
        z_comp = Z()
        z_comp.read(sample_avg_data)

        T_real, freqs, _ = z_comp.to_tensor(var="z_real", station=100)
        assert T_real.shape == (2, 2, 2)

        T_imag, _, _ = z_comp.to_tensor(var="z_imag", station=100)
        assert T_imag.shape == (2, 2, 2)

        T_err, _, _ = z_comp.to_tensor(var="z_err", station=100)
        assert T_err.shape == (2, 2, 2)

    def test_to_tensor_fallback_to_base_for_other_vars(self, sample_avg_data):
        z_comp = Z()
        z_comp.read(sample_avg_data)
        T, freqs, stations = z_comp.to_tensor(var="rho", station=100)
        assert T.shape == (2, 2, 2)
        assert np.isclose(T[1, 0, 1], 50.0)  # ExHy rho at freq=1024

    def test_to_xarray_merges_meta_and_explicit_attrs(self, sample_avg_data):
        z_comp = Z()
        z_comp.read(sample_avg_data, meta={"survey": "K2"})
        da = z_comp.to_xarray(station=100, attrs={"note": "test"})
        assert da.attrs["survey"] == "K2"
        assert da.attrs["note"] == "test"
        assert da.dims == ("freq", "e", "h")

    def test_write_empty_and_nonempty(self, sample_avg_data):
        z_empty = Z()
        out_empty = z_empty.write()
        assert "$Z (Impedance) Block" in out_empty[0]

        z_comp = Z()
        z_comp.read(sample_avg_data)
        out = z_comp.write()
        assert isinstance(out, list)
        assert any("$Z (Impedance) Block" in line for line in out)

    def test_str_and_repr_empty_and_nonempty(self, sample_avg_data):
        z_empty = Z()
        assert str(z_empty) == "Z(status=empty)"
        assert repr(z_empty) == "Z(status=empty)"

        z_comp = Z()
        z_comp.read(sample_avg_data)
        s = str(z_comp)
        assert s.startswith("Z(rows=")
        assert "stations=2" in s
        assert repr(z_comp) == s


if __name__ == "__main__":  # pragma: no-cover
    pytest.main([__file__])
