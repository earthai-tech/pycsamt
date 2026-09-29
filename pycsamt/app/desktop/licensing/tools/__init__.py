# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Vendor-side tooling for issuing pycsamt-desktop license keys.

Nothing under here is imported by the running desktop application --
``issue_license.py`` needs the private half of the signing keypair, which
is deliberately kept outside this repository and never ships in the app
bundle. This subpackage exists purely so earthai-tech has a
version-controlled, tested implementation of the *format* it signs
against, without the secret itself ever being committed.
"""
