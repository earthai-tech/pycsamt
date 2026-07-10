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

    def test_rotate_tensor(self, amtavg_for_methods):
        """Test the 90-degree rotation of the impedance tensor."""
        z_xy_orig = amtavg_for_methods.z.z_xy.copy()
        z_yx_orig = amtavg_for_methods.z.z_yx.copy()

        amtavg_for_methods.rotate(angle=90, update_components=True)

        z_yx_rotated = amtavg_for_methods.z.z_yx
        assert np.allclose(z_yx_rotated, z_xy_orig)

        z_xy_rotated = amtavg_for_methods.z.z_xy
        assert np.allclose(z_xy_rotated, -z_yx_orig)

        assert "rho" in amtavg_for_methods.df.columns
        new_rho_val = amtavg_for_methods.df["rho"].iloc[0]
        expected_rho = (
            np.abs(z_xy_rotated.iloc[0])**2
            / (2 * PI * 1024 * MU_0)
        )
        assert np.isclose(new_rho_val, expected_rho)


if __name__=='__main__': # pragma: no-cover
   pytest.main( [__file__])
