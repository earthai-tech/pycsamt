# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later
"""
Pytest suite for the processing utility functions in proc_utils.
"""

import numpy as np
import pandas as pd
import pytest

from pycsamt.constants import MU_0, PI
from pycsamt.zonge.proc_utils import (
    ama,
    flma,
    interpolate_to_log_space,
    smooth_rho_from_phase,
    tma,
)


# --- Test Data Fixture --------------------------------------------
@pytest.fixture(scope="module")
def processing_data_fixture() -> pd.DataFrame:
    """
    Provides a DataFrame suitable for testing all processing
    utility functions, with 5 data points per station.
    """
    stations = np.repeat([100, 200, 300], 5)
    freqs = np.tile([1024, 512, 256, 128, 64], 3)

    data = {
        "station": stations,
        "freq": freqs,
        "comp": ["ExHy"] * 15,
        "rho": np.linspace(50, 200, 15),
        "phase": np.linspace(1000, 1400, 15),
    }
    df = pd.DataFrame(data)
    # Add complex impedance for FLMA and AMA tests
    omega = 2 * PI * df["freq"]
    z_mag = np.sqrt(df["rho"] * omega * MU_0)
    phase_rad = df["phase"] * 1e-3
    df["z"] = z_mag * np.exp(1j * phase_rad)
    return df


# --- Tests for Processing Utilities -------------------------------


class TestProcessingUtils:
    def test_tma_series_input(self, processing_data_fixture):
        """Test TMA filter with a pandas Series."""
        rho_series = processing_data_fixture["rho"]
        smoothed = tma(rho_series, window_size=3)
        assert isinstance(smoothed, pd.Series)
        assert len(smoothed) == len(rho_series)
        # Check that smoothing occurred
        assert pd.notna(smoothed.iloc[1])

    def test_tma_dataframe_input(self, processing_data_fixture):
        """Test TMA filter with a DataFrame and column name."""
        result_df = tma("rho", data=processing_data_fixture, window_size=3)
        assert isinstance(result_df, pd.DataFrame)
        assert "rho_tma" in result_df.columns
        assert pd.notna(result_df["rho_tma"].iloc[1])

    def test_flma_series_input(self, processing_data_fixture):
        """Test FLMA filter with pandas Series."""
        df = processing_data_fixture
        smoothed = flma(df["z"], stations=df["station"], dipole_length=100)
        assert isinstance(smoothed, pd.Series)
        assert len(smoothed) == len(df)
        assert np.iscomplexobj(smoothed)

    def test_ama_array_input(self, processing_data_fixture):
        """Test AMA filter with NumPy arrays."""
        df = processing_data_fixture
        smoothed = ama(
            df["z"].to_numpy(),
            stations=df["station"].to_numpy(),
            dipole_length=100,
            frequency=1024,
        )
        assert isinstance(smoothed, np.ndarray)
        assert len(smoothed) == len(df)
        assert np.iscomplexobj(smoothed)

    def test_interpolate_to_log_space(self, processing_data_fixture):
        """Test interpolation to a regular frequency grid."""
        df = processing_data_fixture
        interp_df = interpolate_to_log_space(
            df, num_points=10, interp_kind="slinear"
        )
        assert isinstance(interp_df, pd.DataFrame)
        # 3 stations * 10 frequencies = 30 rows
        assert len(interp_df) == 30
        assert interp_df["freq"].nunique() == 10
        log_diffs = np.diff(np.log10(interp_df["freq"].unique()))
        assert np.allclose(log_diffs, log_diffs[0])

    def test_smooth_rho_from_phase(self, processing_data_fixture):
        """Test resistivity smoothing from phase data."""
        df = processing_data_fixture
        smoothed_df = smooth_rho_from_phase(df)
        assert isinstance(smoothed_df, pd.DataFrame)
        assert "rho_smoothed" in smoothed_df.columns
        # With 5 points per station, the result should be valid
        assert pd.notna(smoothed_df["rho_smoothed"]).all()
        assert not np.allclose(
            smoothed_df["rho"], smoothed_df["rho_smoothed"]
        )


if __name__ == "__main__":  # pragma: no-cover
    pytest.main([__file__])
