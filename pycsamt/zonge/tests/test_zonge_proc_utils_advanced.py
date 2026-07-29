# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later
"""
Pytest suite for the advanced analytical functions in proc_utils.
"""

import numpy as np
import pandas as pd
import pytest

from pycsamt.zonge.proc_utils import (
    get_reference_frequency,
    get_skew,
    get_strike,
)


# --- Test Data Fixture --------------------------------------------
@pytest.fixture(scope="module")
def dimensional_data_fixture() -> pd.DataFrame:
    """
    Provides a DataFrame with a full impedance tensor and QC data
    suitable for testing strike, skew, and ref. frequency.
    """
    data = {
        # FIX: Provide a full tensor for station 100 at two freqs
        "station": [100, 100, 100, 100, 100, 100, 100, 100],
        "freq": [1024, 1024, 1024, 1024, 512, 512, 512, 512],
        "comp": ["ExHx", "ExHy", "EyHx", "EyHy"] * 2,
        "rho": [10, 100, 120, 12, 15, 110, 130, 18],
        "phase": [45, 45, 45, 45, 45, 45, 45, 45],
        "pc_rho": [5, 5, 5, 5, 25, 25, 25, 25],
        "z": [
            # Zxx, Zxy, Zyx, Zyy
            1 + 1j,
            10 + 2j,
            8 + 3j,
            1.5 + 1j,  # Tensor for 1024 Hz
            2 + 2j,
            12 + 3j,
            9 + 4j,
            2.5 + 2j,  # Tensor for 512 Hz
        ],
    }
    return pd.DataFrame(data)


# --- Tests for Advanced Processing Utilities ----------------------


class TestAdvancedProcUtils:
    def test_get_reference_frequency_auto(self, dimensional_data_fixture):
        """Test auto-detection of the reference frequency."""
        df = dimensional_data_fixture
        # Station 100 is clean (pc_rho < 20), highest freq is 1024
        # Station 200 is noisy (pc_rho > 20)
        ref_freq = get_reference_frequency(df, mode="auto")
        assert ref_freq == 1024.0

    def test_get_reference_frequency_manual(self, dimensional_data_fixture):
        """Test manual override of the reference frequency."""
        df = dimensional_data_fixture
        ref_freq = get_reference_frequency(df, mode=512.0)
        assert ref_freq == 512.0

    def test_get_skew(self, dimensional_data_fixture):
        """Test the calculation of Swift's skew."""
        df = dimensional_data_fixture
        skew_df = get_skew(df)
        assert isinstance(skew_df, pd.DataFrame)
        assert "skew" in skew_df.columns
        # FIX: Correct the expected length
        assert len(skew_df) == 2  # 1 station x 2 freqs

        # Manual calculation for station 100, freq 1024
        zxx, zxy = 1 + 1j, 10 + 2j
        zyx, zyy = 8 + 3j, 1.5 + 1j
        expected_skew = np.abs(zxx + zyy) / np.abs(zxy - zyx)

        calculated_skew = skew_df[(skew_df.station == 100) & (skew_df.freq == 1024)][
            "skew"
        ].iloc[0]

        assert pd.notna(calculated_skew)
        assert np.isclose(calculated_skew, expected_skew)

    def test_get_strike(self, dimensional_data_fixture):
        """Test the calculation of geoelectric strike."""
        df = dimensional_data_fixture
        strike_df = get_strike(df)
        assert isinstance(strike_df, pd.DataFrame)
        assert "strike_angle" in strike_df.columns
        # FIX: Correct the expected length
        assert len(strike_df) == 2

        # Manual calculation for station 100, freq 1024
        zxx, zxy = 1 + 1j, 10 + 2j
        zyx, zyy = 8 + 3j, 1.5 + 1j
        expected_strike_rad = 0.5 * np.arctan2((zxy + zyx).real, (zxx - zyy).real)
        expected_strike_deg = np.rad2deg(expected_strike_rad)

        calculated_strike = strike_df[
            (strike_df.station == 100) & (strike_df.freq == 1024)
        ]["strike_angle"].iloc[0]

        assert pd.notna(calculated_strike)
        assert np.isclose(calculated_strike, expected_strike_deg)


if __name__ == "__main__":  # pragma: no-cover
    pytest.main([__file__])
