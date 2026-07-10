# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Entry point: python -m pycsamt.app.web"""

from __future__ import annotations

import sys


def main(argv: list[str] | None = None) -> None:
    """Launch the pyCSAMT web application from the command line."""
    from pycsamt.app.web.app import main as _main

    _main(argv=sys.argv[1:] if argv is None else argv)

if __name__ == "__main__":
    main()
