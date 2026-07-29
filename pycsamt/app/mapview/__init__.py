# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""pyCSAMT Map View platform — a focused, interactive map viewer.

The platform is a thin Dash UI over :class:`pycsamt.map.MapView`; every
view (station map, profile, pseudosection, 3-D) is rendered by the
library, so code and GUI stay in sync.

Entry points
------------
>>> from pycsamt.app.mapview import launch
>>> launch()  # doctest: +SKIP

or seed it with data already loaded in code::

    from pycsamt.map import MapView

    MapView.from_folder("data/survey").launch()
"""

from __future__ import annotations

from .app import create_app, launch

__all__ = ["create_app", "launch"]
