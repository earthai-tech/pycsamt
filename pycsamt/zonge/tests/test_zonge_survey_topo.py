# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later
"""
Pytest suite for the Topography and Station classes.
"""

import numpy as np
import pandas as pd
import pytest

from pycsamt.zonge.survey import Topography

# --- Test Data Fixtures -------------------------------------------


@pytest.fixture
def sample_stn_content():
    """Provides content for a sample .stn file."""
    return (
        "! Sample .stn file\n"
        "station,easting,northing,elevation\n"
        "100,500100,4000100,200\n"
        "200,500200,4000105,205\n"
        "300,500300,4000110,210\n"
    )


@pytest.fixture
def sample_stn_file(tmp_path, sample_stn_content):
    """Creates a temporary .stn file and returns its path."""
    f = tmp_path / "sample.stn"
    f.write_text(sample_stn_content)
    return f


# --- Tests for Topography Class -----------------------------------
class TestTopography:
    def test_initialization(self):
        """Test that Topography initializes correctly."""
        topo = Topography()
        assert topo._frame.empty

    def test_read_from_stn_file(self, sample_stn_file):
        """Test reading data from a .stn file."""
        topo = Topography().read(sample_stn_file)
        assert not topo._frame.empty
        assert len(topo._frame) == 3
        assert list(topo._frame.columns) == [
            "station",
            "easting",
            "northing",
            "elevation",
        ]
        assert topo.stations[0] == 100

    def test_generate_utm(self):
        """Test generating a synthetic survey line from UTM coords."""
        start_coord = (500000, 4000000)
        n_stations = 10
        step = 50
        azimuth = 90  # East
        topo = Topography.generate(
            start_coord=start_coord,
            n_stations=n_stations,
            step=step,
            azimuth=azimuth,
        )
        assert isinstance(topo, Topography)
        assert len(topo._frame) == n_stations
        assert topo.easting[0] == start_coord[0]
        assert topo.northing[0] == start_coord[1]
        # Check that the line extends east
        assert np.all(np.diff(topo.easting) > 0)
        assert np.allclose(np.diff(topo.northing), 0)

    def test_generate_ll(self):
        """Test generating a line from Lat/Lon and converting."""
        start_coord = (40.0, -110.0)  # Lat, Lon
        topo = Topography.generate(
            start_coord=start_coord,
            n_stations=5,
            step=100,
            azimuth=0,  # North
            coord_type="ll",
        )
        assert "easting" in topo._frame.columns
        # Check that the UTM coordinates are reasonable
        assert topo.easting[0] > 100000
        assert topo.northing[0] > 1000000

    def test_get_step(self, sample_stn_file):
        """Test the calculation of inter-station distances."""
        topo = Topography().read(sample_stn_file)
        steps = topo.get_step()
        assert isinstance(steps, pd.Series)
        assert len(steps) == 3
        # The first step is always 0
        assert steps.iloc[0] == 0
        # Check the distance between the first two points
        expected_dist = np.hypot(100, 5)
        assert np.isclose(steps.iloc[1], expected_dist)

    def test_correct_coords(self):
        """Test the coordinate regularization method."""
        # Create a slightly irregular line
        data = {
            "station": [1, 2, 3, 4],
            "easting": [100, 200, 305, 400],
            "northing": [100, 102, 98, 101],
            "elevation": [10, 11, 12, 11],
        }
        df = pd.DataFrame(data)
        topo = Topography(data=df)

        # Correct the coordinates in place
        topo.correct_coords(step=100)

        # Check that the new coordinates form a perfect line
        # i.e., the slope between consecutive points is constant
        dx = np.diff(topo.easting)
        dy = np.diff(topo.northing)
        slopes = dy / dx
        assert np.allclose(slopes, slopes[0])

    def test_to_grid(self, stn_file_k1):
        """Test interpolation onto a regular grid."""
        topo = Topography().read(stn_file_k1)
        grid_x, grid_y, grid_z = topo.to_grid(resolution=10)
        assert grid_x.shape == (10, 10)
        assert grid_z.shape == (10, 10)
        # Check that the grid has interpolated values
        assert not np.all(np.isnan(grid_z))


if __name__ == "__main__":  # pragma: no-cover
    pytest.main([__file__])
