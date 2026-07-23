# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Shared fixtures for pycsamt.app.agent_master tests."""

from __future__ import annotations

import pytest


@pytest.fixture(scope="session")
def agent_app():
    """Single fully-wired agent_master Dash app (all callbacks registered)."""
    from pycsamt.app.agent_master.app import create_app

    return create_app()
