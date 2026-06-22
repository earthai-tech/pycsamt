# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""pycsamt pipe — processing pipeline CLI commands.

Importing this package registers all sub-commands onto the ``pipe``
Click group via side-effect imports.
"""

from ._base    import pipe   # noqa: F401  (re-exported for cli/_base.py)
from .         import run    # noqa: F401  registers pipe run
from .         import steps  # noqa: F401  registers pipe steps
from .         import presets# noqa: F401  registers pipe presets
from .         import init   # noqa: F401  registers pipe init
from .         import show   # noqa: F401  registers pipe show
