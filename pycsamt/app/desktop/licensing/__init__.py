# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
licensing — pycsamt-desktop's license/trial gate.

``LicenseManager``/``TrialState``/``LicenseStatus`` (``interfaces.py``) are
the stable shapes every caller depends on. ``OfflineLicenseManager``
(``manager.py``) is the real Phase 11 implementation: a locally persisted,
Ed25519-signed license key (``signing.py``) combined with a local 30-day
trial tracker (``trial.py``), bound to a lightweight machine fingerprint
(``fingerprint.py``). Everything is offline -- no network call is ever made
to check a license or trial. ``NullLicenseManager`` remains available for
tests and any code that deliberately wants an unconditional trial.

Importing this package is deliberately side-effect-free and Qt-independent
(Qt/branding imports inside ``manager.py``/``trial.py`` are deferred into
function bodies) so it can be imported from the splash screen's bootstrap
path before ``QApplication`` exists.
"""

from .interfaces import LicenseManager, LicenseStatus, TrialState
from .manager import OfflineLicenseManager, get_default_manager
from .null_manager import NullLicenseManager

__all__ = [
    "LicenseManager",
    "LicenseStatus",
    "TrialState",
    "NullLicenseManager",
    "OfflineLicenseManager",
    "get_default_manager",
]
