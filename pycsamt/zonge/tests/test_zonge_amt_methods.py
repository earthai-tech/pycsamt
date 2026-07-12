# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later
"""
Pytest suite for the analytical methods of the AMTAVG class.
"""

import numpy as np
import pandas as pd
import pytest

from pycsamt.constants import MU_0, PI
from pycsamt.zonge.avg import AMTAVG


# --- Test Data Fixture --------------------------------------------
@pytest.fixture(scope="module")
def analytical_avg_data() -> pd.DataFrame:
    """
    Provides a DataFrame with data suitable for testing
    statistics, phase unwrapping, and tensor rotation.
    """
    # Use modern, dot-notation column names to mimic a real file
    data = {
        "station": [100, 100, 100, 100, 200, 200, 200, 200],
        "freq":    [1024, 1024,  512,  512, 1024,  512,  256,  128],
        "comp":    ["ExHy", "ExHy", "ExHy", "ExHy", "ExHy", "EyHx",
                    "ExHy", "EyHx"],
        "ARes.mag": [50.0, 52.0, 45.0, 47.0, 60.0, 160.0, 55.0, 155.0],
        "Z.phz":   [1000, 1000, 3200, 3200, -500, 1200, 1800, 4000],
        "E.mag":   [20, 22, 18, 20, 30, 80, 28, 78],
        "B.mag":   [0.1, 0.11, 0.09, 0.1, 0.15, 0.2, 0.14, 0.19],

        "E.phz":   [1100, 1100, 3300, 3300, -400, 1300, 1900, 4100],
        "B.phz":   [100, 100, 100, 100, 100, 100, 100, 100],
        # Add  QC columns
        "E.%err":  [3.0, 3.2, 2.8, 2.9, 4.0, 4.2, 4.1, 4.3],
        "B.%err":  [4.0, 4.2, 3.8, 3.9, 5.0, 5.2, 5.1, 5.3],

        "ARes.%err": [5.0, 6.0, 4.5, 5.5, 7.0, 8.0, 7.5, 8.5],
        "Z.perr":  [20, 25, 18, 22, 28, 30, 29, 31],
        # Dummy Zxx/Zyy data for rotation test
        "zxx_real": [1, 1, 1, 1, 2, 2, 2, 2],
        "zyy_real": [1, 1, 1, 1, 2, 2, 2, 2],
    }
    df = pd.DataFrame(data)
    # Create full Z tensor for rotation
    # Use resistivity and phase to create a consistent z_complex
    omega = 2 * PI * df['freq']
    z_mag = np.sqrt(df['ARes.mag'] * omega * MU_0)
    phase_rad = df['Z.phz'] * 1e-3
    z_complex = z_mag * np.exp(1j * phase_rad)

    df['zxx'] = z_complex * df['zxx_real']
    df['zyy'] = z_complex * df['zyy_real']
    df['zxy'] = z_complex * 5
    df['zyx'] = z_complex * 10

    return df

@pytest.fixture
def amtavg_for_methods(analytical_avg_data):
    """
    Provides a fresh AMTAVG instance for each test function.
    """
    avg = AMTAVG(verbose=True)
    avg.read(analytical_avg_data)
    # Manually add Z components for rotation test
    # The Z component calculates its own z_complex, so we need to
    # add the off-diagonal elements for the rotation test.
    df = avg.z._frame
    df['zxx'] = analytical_avg_data['zxx']
    df['zyy'] = analytical_avg_data['zyy']
    df['zxy'] = analytical_avg_data['zxy']
    df['zyx'] = analytical_avg_data['zyx']

    return avg

# --- Tests for AMTAVG Methods -------------------------------------

class TestAMTAVGMethods:

    def test_calculate_statistics(self, amtavg_for_methods):
        """Test the calculation of statistical QC metrics."""
        stats_df = amtavg_for_methods.calculate_statistics(
            update_components=True
        )
        assert not stats_df.empty
        assert "rho_mean" in stats_df.columns
        assert "rho_std" in stats_df.columns
        assert "rho_cvar" in stats_df.columns

        st100_f1024 = stats_df[
            (stats_df.station == 100) & (stats_df.freq == 1024)
        ]
        assert np.isclose(st100_f1024["rho_mean"].iloc[0], 51.0)
        assert np.isclose(
            st100_f1024["rho_cvar"].iloc[0], 100 * (2**0.5) / 51.0
        )

        main_df = amtavg_for_methods.info.df
        assert "pc_rho" in main_df.columns
        assert pd.notna(main_df["pc_rho"].iloc[0])

    def test_unwrap_phase(self, amtavg_for_methods):
        """Test the phase unwrapping functionality."""
        original_phase = amtavg_for_methods.df["phase"].copy()
        unwrapped_df = amtavg_for_methods.unwrap_phase(
            update_components=True
        )

        assert not unwrapped_df["phase"].equals(original_phase)

        st100_phase = amtavg_for_methods.df[
            amtavg_for_methods.df.station == 100
        ].sort_values("freq")["phase"]
        assert np.all(np.abs(np.diff(st100_phase)) < 1000 * PI)

    def test_unwrap_phase_to_degree(self, amtavg_for_methods):
        """Test phase unwrapping with conversion to degrees."""
        amtavg_for_methods.unwrap_phase(
            todeg=True, update_components=True
        )
        assert amtavg_for_methods.info.meta["Unit.Phase"] == "deg"
        assert -360 < amtavg_for_methods.df["phase"].iloc[0] < 360

    def test_rotate_tensor(self):
        """90-degree rotation swaps the off-diagonal components.

        With Zxx = Zyy = 0 (scalar-CSAMT data), a 90-degree rotation
        gives Zxy' = -Zyx and Zyx' = -Zxy, so the apparent resistivity
        of each sounding's XY row must take the value of its YX row
        and vice versa.
        """
        data = {
            "station":   [100, 100, 100, 100, 200, 200, 200, 200],
            "freq":      [1024, 1024, 512, 512, 1024, 1024, 512, 512],
            "comp":      ["ExHy", "EyHx"] * 4,
            "ARes.mag":  [50.0, 150.0, 45.0, 145.0,
                          60.0, 160.0, 55.0, 155.0],
            "Z.phz":     [1000, 1200, 3200, 3400,
                          -500, -300, 1800, 2000],
            "E.mag":     [20, 80, 18, 78, 30, 82, 28, 84],
            "B.mag":     [0.1, 0.2, 0.09, 0.19, 0.15, 0.21, 0.14, 0.2],
            "E.phz":     [1100, 1300, 3300, 3500, -400, -200, 1900, 2100],
            "B.phz":     [100] * 8,
            "E.%err":    [3.0, 4.2, 2.8, 4.3, 4.0, 4.4, 4.1, 4.5],
            "B.%err":    [4.0, 5.2, 3.8, 5.3, 5.0, 5.4, 5.1, 5.5],
            "ARes.%err": [5.0, 8.0, 4.5, 8.5, 7.0, 9.0, 7.5, 9.5],
            "Z.perr":    [20, 30, 18, 31, 28, 32, 29, 33],
        }
        df = pd.DataFrame(data)
        avg = AMTAVG(verbose=0)
        avg.read(df)

        before = avg.df[["station", "freq", "comp", "rho"]].copy()
        rotated = avg.rotate(angle=90, update_components=True)

        # the returned frame holds one row per (station, freq) pair
        assert rotated.shape[0] == 4
        assert {"z_xy_rot", "z_yx_rot"} <= set(rotated.columns)

        after = avg.df[["station", "freq", "comp", "rho"]].copy()
        # comp labels keep their file casing; normalise for lookups
        before["comp"] = before["comp"].str.upper()
        after["comp"] = after["comp"].str.upper()
        key = ["station", "freq"]
        for (st, fr), grp in after.groupby(key):
            orig = before[
                (before["station"] == st) & (before["freq"] == fr)
            ].set_index("comp")["rho"]
            new = grp.set_index("comp")["rho"]
            # |Z| is sign-insensitive: rho swaps between XY and YX
            assert np.isclose(new["EXHY"], orig["EYHX"])
            assert np.isclose(new["EYHX"], orig["EXHY"])


if __name__=='__main__': # pragma: no-cover
   pytest.main( [__file__])
