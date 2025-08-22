# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later
"""
Pytest suite for the Tipper component and calculation method.
"""
import pytest
import pandas as pd
import numpy as np

from pycsamt.zonge.avg import AMTAVG
from pycsamt.zonge.tipper import Tipper
from pycsamt.exceptions import AvgDataError

# --- Test Data Fixture --------------------------------------------

@pytest.fixture(scope="module")
def tipper_data_fixture():
    """
    Provides an AMTAVG instance and Hz data suitable for testing
    the Tipper calculation.
    """
    # Create a DataFrame with Hx and Hy components at two freqs
    # Use modern, dot-notation names to mimic a real file
    data = {
        "station": [101, 101, 101, 101],
        "freq":    [1024, 1024, 512, 512],
        "comp":    ["ExHx", "ExHy", "ExHx", "ExHy"],
        "B.mag":   [1.0, 2.0, 1.5, 2.5],
        "B.phz":   [0, 90000, 0, 90000], # 90 deg in mrad -> imag
        "ARes.mag": [100]*4, # Dummy data for other components
        "Z.phz":   [45]*4,
    }
    df = pd.DataFrame(data)

    # Create an AMTAVG instance from this data
    avg = AMTAVG(verbose=True)
    avg.read(df)

    # Create a corresponding complex Hz Series.
    # We'll define known Tipper values (Tx, Ty) to generate Hz.
    # Let Tx = 0.1 + 0.2j and Ty = 0.3 + 0.4j
    tx, ty = (0.1 + 0.2j), (0.3 + 0.4j)

    # Reconstruct Hx and Hy from the STANDARDIZED dataframe
    # h_complex = avg.info.df["hmag"] * np.exp(
    #     1j * avg.info.df["hphz"] * 1e-3
    # )
    # hx = h_complex.where(avg.info.df["comp"] == "ExHx", 0)
    # hy = h_complex.where(avg.info.df["comp"] == "ExHy", 0)
    
    # # Group Hx and Hy by station and freq
    # h_vectors = pd.DataFrame({'hx': hx, 'hy': hy}).groupby(
    #     [df['station'], df['freq']]).sum()

    # # Calculate Hz = Tx*Hx + Ty*Hy
    # hz_values = (
    #     tx * h_vectors['hx'] + ty * h_vectors['hy']
    # ).values
    # Create an Hz series that aligns with the original DataFrame
    # hz_series = pd.Series(
    #     np.repeat(hz_values, 2), index=df.index, dtype='complex128'
    #     )
    
    # Reconstruct Hx and Hy from the STANDARDIZED dataframe
    h_complex = avg.info.df["hmag"] * np.exp(
        1j * avg.info.df["hphz"] * 1e-3
    )
    df_h = avg.info.df.copy()
    df_h["hx"] = h_complex.where(df_h["comp"] == "ExHx", 0)
    df_h["hy"] = h_complex.where(df_h["comp"] == "ExHy", 0)
 
    # Calculate Hz row-by-row, not from the sum ---
    hz_series = tx * df_h["hx"] + ty * df_h["hy"]
    hz_series = hz_series.astype('complex128')
   
    return avg, hz_series

# --- Tests for Tipper Class ---------------------------------------

class TestTipper:
    def test_tipper_initialization(self):
        """Test that the Tipper component initializes correctly."""
        tipper = Tipper()
        assert tipper._frame.empty

    def test_read_from_dataframe(self):
        """Test populating the Tipper from a valid DataFrame."""
        df = pd.DataFrame({
            "station": [101], "freq": [1024],
            "tx": [0.1+0.2j], "ty": [0.3+0.4j]
        })
        tipper = Tipper()
        tipper.read(df)
        assert not tipper.frame.empty
        assert "tx" in tipper.frame.columns
        assert np.isclose(tipper.tx.iloc[0], 0.1 + 0.2j)

    def test_read_missing_columns(self):
        """Test that read handles missing tx/ty by creating NaNs."""
        df = pd.DataFrame({"station": [101], "freq": [1024]})
        tipper = Tipper(verbose=True)
        tipper.read(df)
        assert "tx" in tipper.frame.columns
        assert "ty" in tipper.frame.columns
        assert pd.isna(tipper.frame["tx"].iloc[0])

    def test_to_xarray(self):
        """Test exporting Tipper data to xarray."""
        df = pd.DataFrame({
            "station": [101, 101], "freq": [1024, 512],
            "tx": [0.1+0.2j, 0.11+0.21j],
            "ty": [0.3+0.4j, 0.31+0.41j]
        })
        tipper = Tipper()
        tipper.read(df)
        ds = tipper.to_xarray()
        assert 'xarray.core.dataset.Dataset' in str(type(ds))
        assert "tx" in ds.data_vars
        assert "ty" in ds.data_vars
        assert ds.dims["station"] == 1
        assert ds.dims["freq"] == 2

# --- Tests for AMTAVG.calculate_tipper ----------------------------

class TestCalculateTipper:
    def test_calculation_and_update(self, tipper_data_fixture):
        """Test successful Tipper calculation and component update."""
        avg, hz_data = tipper_data_fixture
        tipper_df = avg.calculate_tipper(
            hz_data, update_components=True
        )

        assert isinstance(tipper_df, pd.DataFrame)
        assert "tx" in tipper_df.columns
        assert "ty" in tipper_df.columns
        assert len(tipper_df) == 2 # One row per frequency

        # Check numerical correctness
        # Expected values are Tx = 0.1 + 0.2j, Ty = 0.3 + 0.4j
        assert np.isclose(tipper_df["tx"].iloc[0], 0.1 + 0.2j)
        assert np.isclose(tipper_df["ty"].iloc[0], 0.3 + 0.4j)

        # Check that the component was created and populated
        assert hasattr(avg.info, "tipper")
        assert isinstance(avg.info.tipper, Tipper)
        assert not avg.info.tipper.frame.empty

    def test_calculation_no_update(self, tipper_data_fixture):
        """Test calculation without updating components."""
        avg, hz_data = tipper_data_fixture
        avg.calculate_tipper(hz_data, update_components=False)

        # The tipper component should not exist on the info object
        assert not hasattr(avg.info, "tx")

    def test_missing_hmag_raises_error(self, tipper_data_fixture):
        """Test that an error is raised if hmag is missing."""
        avg, hz_data = tipper_data_fixture
        # Create a copy and drop the required column
        df_no_hmag = avg.info.df.drop(columns=["hmag"])
        avg.info.read(df_no_hmag, avg.info.meta)

        with pytest.raises(AvgDataError, match="Horizontal magnetic"):
            avg.calculate_tipper(hz_data)

if __name__=='__main__': # pragma: no-cover 
   pytest.main( [__file__])