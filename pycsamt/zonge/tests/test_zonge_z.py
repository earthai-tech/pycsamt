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
        assert z_comp.frame["comp"].tolist() == expected_comps, (
            "Component normalization in Z.read() failed."
        )

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


if __name__ == "__main__":  # pragma: no-cover
    pytest.main([__file__])
