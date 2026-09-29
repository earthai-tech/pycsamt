# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later
"""
Pytest suite for the ASTATIC processing class.
"""

import numpy as np
import pandas as pd
import pytest

from pycsamt.exceptions import NotReadError, ProcessingError
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


# --- Coverage: constructor, _resolve_param branches, flma/ama, and  ---
# --- update_components=False / short-station skip branches.        ---


def _basic_data():
    return {
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


class TestASTATICCoverage:
    def test_constructor_reads_avg_data_directly(self):
        avg = AMTAVG(verbose=False).read(pd.DataFrame(_basic_data()))
        processor = ASTATIC(avg_data=avg, verbose=True)
        assert processor.avg is avg

    def test_capacitive_coupling_string_and_series_params(self):
        data = _basic_data()
        data["Contact"] = [5000.0] * 10
        data["SetupLen"] = [100.0] * 10
        avg = AMTAVG(verbose=False).read(pd.DataFrame(data))
        processor = ASTATIC().read(avg)

        out = processor.correct_capacitive_coupling(
            contact_resistance="contact",
            setup_length="setuplen",
            update_components=False,
        )
        assert list(out.columns) == ["emag", "ephz"]

        # Series param path
        zc_series = pd.Series([5000.0] * 10)
        out2 = processor.correct_capacitive_coupling(
            contact_resistance=zc_series,
            setup_length=100.0,
            update_components=False,
        )
        assert not out2.empty

    def test_capacitive_coupling_missing_string_column_raises(self):
        avg = AMTAVG(verbose=False).read(pd.DataFrame(_basic_data()))
        processor = ASTATIC().read(avg)
        with pytest.raises(ProcessingError, match="not found"):
            processor.correct_capacitive_coupling(
                contact_resistance="does_not_exist",
                setup_length=100.0,
            )

    def test_capacitive_coupling_update_components_false_leaves_avg_unchanged(self):
        avg = AMTAVG(verbose=False).read(pd.DataFrame(_basic_data()))
        original_emag = avg.df["emag"].copy()
        processor = ASTATIC().read(avg)
        processor.correct_capacitive_coupling(
            contact_resistance=5000.0,
            setup_length=100.0,
            update_components=False,
        )
        # main object untouched since update_components=False
        assert avg.df["emag"].equals(original_emag)

    def test_static_shift_flma_filter(self, loaded_amtavg_instance):
        processor = ASTATIC().read(loaded_amtavg_instance)
        # dipole_length=150 keeps every station's Hanning window at >=4
        # points (spacing is 100), avoiding scipy's hann(2) == [0, 0]
        # zero-weight edge case in `flma`.
        stats_df = processor.correct_static_shift(
            reference_freq=1024,
            filter_method="flma",
            dipole_length=150.0,
            update_components=False,
        )
        assert "shift_factor" in stats_df.columns
        assert not stats_df.empty

    def test_static_shift_ama_filter(self):
        # A dedicated high-resistivity, closely-spaced fixture keeps the
        # adaptive skin-depth window covering all stations on every
        # iteration, avoiding scipy's hann(2) == [0, 0] zero-weight case
        # in `ama` that a low-resistivity / widely-spaced profile hits.
        data = {
            "station": [0, 0, 10, 10, 20, 20],
            "freq": [1024, 512, 1024, 512, 1024, 512],
            "comp": ["ExHy"] * 6,
            "ARes.mag": [1.0e6] * 6,
            "Z.phz": [900] * 6,
            "E.mag": [1] * 6,
            "B.mag": [1] * 6,
            "E.%err": [1] * 6,
            "B.%err": [1] * 6,
            "E.phz": [900] * 6,
            "ARes.%err": [1] * 6,
            "E.perr": [1] * 6,
            "B.perr": [1] * 6,
            "Z.perr": [1] * 6,
        }
        avg = AMTAVG(verbose=False).read(pd.DataFrame(data))
        processor = ASTATIC().read(avg)
        stats_df = processor.correct_static_shift(
            reference_freq=1024,
            filter_method="ama",
            dipole_length=1.0,
            update_components=False,
        )
        assert "shift_factor" in stats_df.columns
        assert not stats_df.empty
        assert np.allclose(stats_df["shift_factor"], 1.0, atol=1e-6)

    def test_static_shift_unknown_filter_raises(self, loaded_amtavg_instance):
        processor = ASTATIC().read(loaded_amtavg_instance)
        with pytest.raises(ValueError, match="Unknown filter method"):
            processor.correct_static_shift(
                reference_freq=1024, filter_method="bogus"
            )

    def test_static_shift_skips_station_with_single_frequency(self):
        data = _basic_data()
        # station 600 has only a single row -> must be skipped (len<2)
        extra = pd.DataFrame(
            {
                "station": [600],
                "freq": [1024],
                "comp": ["ExHy"],
                "ARes.mag": [90],
                "Z.phz": [1000],
                "E.mag": [1],
                "B.mag": [1],
                "E.%err": [1],
                "B.%err": [1],
                "E.phz": [1050],
                "ARes.%err": [1],
                "E.perr": [1],
                "B.perr": [1],
                "Z.perr": [1],
            }
        )
        df = pd.concat([pd.DataFrame(data), extra], ignore_index=True)
        avg = AMTAVG(verbose=False).read(df)
        processor = ASTATIC().read(avg)
        stats_df = processor.correct_static_shift(
            reference_freq=1024, filter_method="tma", update_components=False
        )
        # station 600 has a single frequency sample and must be skipped
        # rather than crash the final DataFrame assembly (see fixed bug).
        assert 600 not in stats_df["station"].tolist()
        assert set(stats_df["station"]) == {100, 200, 300, 400, 500}


if __name__ == "__main__":  # pragma: no-cover
    pytest.main([__file__])
