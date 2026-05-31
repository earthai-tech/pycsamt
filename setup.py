"""
setup.py — legacy stub (pyCSAMT v2)

All package metadata and build configuration has moved to pyproject.toml (PEP 621).
This file is kept only so that tools that require a setup.py entry point can find
it; it raises an error immediately to prevent accidental use.

    pip install .          # uses pyproject.toml via pip >= 21.3
    pip install -e ".[dev]"
    python -m build        # uses pyproject.toml via build
"""
import sys

sys.exit(
    "\n"
    "ERROR: setup.py is no longer used in pyCSAMT v2.\n"
    "Please install via pip:\n\n"
    "    pip install .\n"
    "    pip install -e '.[dev]'\n\n"
    "All metadata is now in pyproject.toml.\n"
)
