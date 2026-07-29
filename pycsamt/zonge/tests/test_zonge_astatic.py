# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later
"""
Pytest suite for the ASTATIC processing class.
"""

import numpy as np
import pandas as pd
import pytest

from pycsamt.exceptions import NotReadError
from pycsamt.zonge.avg import AMTAVG
from pycsamt.zonge.processing import ASTATIC

# --- Test Data Fixture --------------------------------------------


@pytest.fixture(scope="module")
def loaded_amtavg_instance() -> AMTAVG:
    """
    Provides a fully loaded AMTAVG instance suitable for processing.
    """
    data = {
        "station": [100, 100, 200, 200, 300, 300, 400, 400, 500, 500],
        "freq": [1024, 512, 1024, 512, 1024, 512, 1024, 512, 1024, 512],
        "comp": ["ExHy"] * 10,
        "ARes.mag": [50, 100, 25, 80, 200, 150, 40, 70, 60, 120],
        "Z.phz": [1000, 1100, 1050, 1150, 950, 1050, 1020, 1120, 980, 1080],
        "E.mag": [1] * 10,
        "B.mag": [1] * 10,
        "E.%err": [1] * 10,
        "B.%err": [1] * 10,
        "E.phz": [1050, 1150, 1100, 1200, 1000, 1100, 1070, 1170, 1030, 1130],
        "ARes.%err": [1] * 10,
        "E.perr": [1] * 10,
        "B.perr": [1] * 10,
        "Z.perr": [1] * 10,
    }
    df = pd.DataFrame(data)
    return AMTAVG(verbose=True).read(df)


# --- Tests for ASTATIC Class --------------------------------------


class TestASTATIC:
    def test_init_with_loaded_avg(self, loaded_amtavg_instance):
        """Test initialization directly with a loaded AVG object."""
        processor = ASTATIC().read(loaded_amtavg_instance)
        assert processor.avg is loaded_amtavg_instance
        assert processor.avg.__has_read__()

    def test_init_empty_then_read(self, loaded_amtavg_instance):
        """Test initializing empty and then using the read method."""
        processor = ASTATIC(verbose=True)
        assert processor.avg is None
        processor.read(loaded_amtavg_instance)
        assert processor.avg is loaded_amtavg_instance

    def test_init_with_unloaded_avg_raises_error(self):
        """Test that ASTATIC raises an error if given an unloaded AVG."""
        unloaded_avg = ASTATIC()
        with pytest.raises(NotReadError):
            # let try to apply correction statistic shift
            unloaded_avg.correct_static_shift(reference_freq=10.24, dipole_length=20)

    def test_read_with_filepath(self, modern_data_file):
        """Test the read method with a direct file path."""
        processor = ASTATIC(verbose=True)
        processor.read(modern_data_file)
        assert isinstance(processor.avg, AMTAVG)
        assert processor.avg.__has_read__()
        assert processor.avg._source_path.name == "K2.AVG"

    def test_correct_static_shift_tma(self, loaded_amtavg_instance):
        """Test the TMA static shift correction method."""
        avg_obj = loaded_amtavg_instance
        original_rho = avg_obj.df["rho"].copy()

        processor = ASTATIC().read(avg_obj)
        stats_df = processor.correct_static_shift(
            reference_freq=1024, filter_method="tma", window_size=3
        )

        assert isinstance(stats_df, pd.DataFrame)
        assert "shift_factor" in stats_df.columns
        # Check that the main DataFrame's rho values have changed
        assert not avg_obj.df["rho"].equals(original_rho)
        # Check that the static corrected column was created
        assert "rho_sc" in avg_obj.df.columns
        # The corrected value should equal the new rho
        assert np.allclose(avg_obj.df["rho"], avg_obj.df["rho_sc"])

    def test_correct_capacitive_coupling(self, loaded_amtavg_instance):
        """Test the capacitive coupling correction method."""
        avg_obj = loaded_amtavg_instance
        original_emag = avg_obj.df["emag"].copy()

        processor = ASTATIC().read(avg_obj)
        # Use simple scalar values for the test
        processor.correct_capacitive_coupling(
            contact_resistance=5000.0, setup_length=100.0
        )

        assert not avg_obj.df["emag"].equals(original_emag)
        # Check if Z was recomputed (rho will change as a result)
        assert "rho" in avg_obj.df.columns


if __name__ == "__main__":  # pragma: no-cover
    pytest.main([__file__])
