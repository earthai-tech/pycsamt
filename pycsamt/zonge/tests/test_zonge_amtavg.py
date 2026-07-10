# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later
"""
Pytest suite for the AMTAVG class.
"""
import numpy as np
import pandas as pd
import pytest

from pycsamt.zonge.avg import AMTAVG


@pytest.fixture(scope="module")
def full_modern_avg_content():
    """
    Provides a comprehensive modern (kind-2) AVG file content
    with all necessary columns for AMTAVG testing.
    """
    return (
        '\\ AMTAVG 7.76: "FULL.fld", Dated 2024-01-01\n'
        "$Survey.Type=CSAMT\n"
        "$Unit.E=nV/Am\n$Unit.B=pT/A\n$Unit.Phase=mrad\n"
        "\n"
        "$Rx.Stn=101\n$Rx.Cmp=ExHy\n"
        "Z.mwgt,Freq,Tx.Amp,ARes.mag,Z.phz,E.mag,B.mag,"
        "E.%err,B.%err,ARes.%err,E.perr,B.perr,Z.perr\n"
        "1,1024,10,1200,450,220,0.09,3,4,5,22,23,16\n"
        "1,512,10,1400,650,260,0.13,2,5,6,20,45,18\n"
        "\n"
        "$Rx.Stn=101\n$Rx.Cmp=EyHx\n"
        "Z.mwgt,Freq,Tx.Amp,ARes.mag,Z.phz,E.mag,B.mag,"
        "E.%err,B.%err,ARes.%err,E.perr,B.perr,Z.perr\n"
        "1,1024,10,800,300,150,0.08,2,3,4,15,20,12\n"
        "1,512,10,900,400,180,0.11,3,4,5,18,30,14\n"
    )

@pytest.fixture(scope="module")
def amtavg_instance(tmp_path_factory, full_modern_avg_content):
    """
    Creates a temporary modern AVG file and returns an initialized
    AMTAVG instance from it. This fixture is module-scoped for
    efficiency as the loading process is the same for all tests.
    """
    path = (
        tmp_path_factory.mktemp("data") / "full_modern.avg"
    )
    path.write_text(full_modern_avg_content)
    return AMTAVG.from_file(path, verbose=True)

# --- Tests for AMTAVG ---------------------------------------------
class TestAMTAVG:
    def test_inheritance(self):
        """Test that AMTAVG inherits from AVG."""
        from pycsamt.zonge.avg import AVG
        assert issubclass(AMTAVG, AVG)

    def test_z_tensor_properties(self, amtavg_instance):
        """Test impedance tensor component properties."""
        assert isinstance(amtavg_instance.z_xy, pd.Series)
        assert len(amtavg_instance.z_xy) == 2
        assert len(amtavg_instance.z_yx) == 2
        # Zxx and Zyy should be empty as they are not in the file
        assert amtavg_instance.z_xx.empty
        assert amtavg_instance.z_yy.empty
        # Check error properties
        assert isinstance(amtavg_instance.z_xy_err, pd.Series)
        assert len(amtavg_instance.z_xy_err) == 2

    def test_resistivity_tensor_properties(self, amtavg_instance):
        """Test resistivity tensor component properties."""
        assert isinstance(amtavg_instance.res_xy, pd.Series)
        assert len(amtavg_instance.res_xy) == 2
        assert amtavg_instance.res_xy.iloc[0] == 1200
        assert amtavg_instance.res_yx.iloc[0] == 800
        assert amtavg_instance.res_xx.empty

    def test_phase_tensor_properties(self, amtavg_instance):
        """Test phase tensor component properties."""
        assert isinstance(amtavg_instance.phase_xy, pd.Series)
        assert len(amtavg_instance.phase_xy) == 2
        assert amtavg_instance.phase_xy.iloc[0] == 450
        assert amtavg_instance.phase_yx.iloc[0] == 300
        assert amtavg_instance.phase_yy.empty

    def test_compute_resistivity_phase(self, amtavg_instance):
        """Test the calculation of rho and phi from Z."""
        rho, phi, rho_err, phi_err = (
            amtavg_instance.compute_resistivity_phase()
        )
        assert isinstance(rho, pd.Series)
        assert len(rho) == 4
        assert np.isclose(rho.iloc[0], 1200, rtol=1e-2)
        assert isinstance(phi_err, pd.Series)

    def test_set_resistivity_phase(self, amtavg_instance):
        """Test updating rho/phi and recalculating Z."""
        # Create new data with the same index as the original df
        df = amtavg_instance.df
        new_rho = pd.Series(np.full(len(df), 999.0), index=df.index)
        new_phi = pd.Series(np.full(len(df), 100.0), index=df.index)

        # Set new values
        amtavg_instance.set_resistivity_phase(new_rho, new_phi)

        # Verify that the internal components have been updated
        assert amtavg_instance.resistivity.frame["rho"].iloc[0] == 999.0
        assert amtavg_instance.phase.frame["phase"].iloc[0] == 100.0
        # Check if Z was recomputed (it will be different now)
        assert np.isclose(
            np.abs(amtavg_instance.z.z.iloc[0]), 2.84, rtol=1e-2
        )

    def test_get_tensor_by_station(self, amtavg_instance):
        """Test fetching a single station's tensor as xarray."""
        station_id = 101
        station_tensor = amtavg_instance.get_tensor_by_station(
            station_id, var="rho"
        )
        assert "xarray.core.dataarray.DataArray" in str(type(station_tensor))
        assert station_tensor.shape == (2, 2, 2)  # (freq, e, h)
        assert "station" not in station_tensor.dims
        assert station_tensor.coords["freq"].size == 2

    def test_get_tensor_by_station_not_found(self, amtavg_instance):
        """Test error handling for a non-existent station."""
        with pytest.raises(KeyError):
            amtavg_instance.get_tensor_by_station(999)

    def test_repr_and_str(self, amtavg_instance):
        """Test the string representations of the AMTAVG object."""
        assert "AMTAVG" in str(amtavg_instance)
        assert "source='full_modern.avg'" in str(amtavg_instance)
        assert "AMTAVG.from_file(" in repr(amtavg_instance)
        assert "full_modern.avg" in repr(amtavg_instance)

if __name__=='__main__': # pragma: no-cover
   pytest.main( [__file__])
