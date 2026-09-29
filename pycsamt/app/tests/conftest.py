# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
This directory now only holds tests for code that lives directly in
``pycsamt/app/`` itself -- the shared, backend-neutral overlay modules
(``_borehole.py``, ``_geology.py``, ``_patterns.py``, ``_points.py``,
``_structure.py``) that Map View's Geology/Borehole rails build on. None
of them touch Qt, Dash, or torch/scipy together, so no special fixture is
needed here.

Everything that used to live in this directory alongside these five files
-- the PySide6 desktop suite, the Map View (Dash) suite, and the web app
(Dash) suite -- has moved to sit next to the code it tests:
``pycsamt/app/desktop/tests``, ``pycsamt/app/mapview/tests``, and
``pycsamt/app/web/tests`` respectively, each with its own conftest.py.
"""

from __future__ import annotations
