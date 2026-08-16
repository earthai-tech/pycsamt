# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""pycsamt pipe — processing pipeline CLI commands.

Importing this package registers all sub-commands onto the ``pipe``
Click group via side-effect imports.
"""

from . import (
    history,  # noqa: F401  registers pipe history
    init,  # noqa: F401  registers pipe init
    plugins,  # noqa: F401  registers pipe plugins
    presets,  # noqa: F401  registers pipe presets
    run,  # noqa: F401  registers pipe run
    show,  # noqa: F401  registers pipe show
    steps,  # noqa: F401  registers pipe steps
)
from ._base import (
    pipe,  # noqa: F401  (re-exported for cli/_base.py)
)
