# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later
"""
Configuration for the pytest suite.

This file defines shared fixtures for locating and accessing test
data files, ensuring that tests can be run consistently from any
directory.
"""
import pytest
from pathlib import Path

def get_project_root() -> Path:
    """Find the project root by looking for the 'pycsamt' dir."""
    # Starts from the current file and goes up the hierarchy
    current_path = Path(__file__).resolve()
    for parent in current_path.parents:
        if (parent / 'pycsamt').is_dir():
            return parent
    raise FileNotFoundError("Could not find project root.")

@pytest.fixture(scope="session")
def project_root() -> Path:
    """Fixture to provide the project's root directory."""
    return get_project_root()

@pytest.fixture(scope="session")
def data_path(project_root) -> Path:
    """Fixture to provide the path to the test data directory."""
    return project_root / "data" / "avg"

@pytest.fixture(scope="session")
def legacy_data_file(data_path):
    """
    Fixture for the path to the legacy K1.AVG data file.
    Skips the test if the file is not found.
    """
    file_path = data_path / "K1.AVG"
    if not file_path.exists():
        pytest.skip(f"Legacy data file not found at: {file_path}")
    return file_path

@pytest.fixture(scope="session")
def modern_data_file(data_path):
    """
    Fixture for the path to the modern K2.AVG data file.
    Skips the test if the file is not found.
    """
    file_path = data_path / "K2.AVG"
    if not file_path.exists():
        pytest.skip(f"Modern data file not found at: {file_path}")
    return file_path
