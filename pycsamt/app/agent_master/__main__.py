# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""python -m pycsamt.app.agent_master"""

from __future__ import annotations

import argparse
import sys


def _parse() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="pycsamt-agent-master",
        description=(
            "Launch the pyCSAMT Agent Master GUI."
        ),
    )
    p.add_argument(
        "--host",
        default="127.0.0.1",
        help="Bind address (default 127.0.0.1)",
    )
    p.add_argument(
        "--port",
        type=int,
        default=8765,
        help="HTTP port (default 8765)",
    )
    p.add_argument(
        "--debug",
        action="store_true",
        help="Enable Dash debug / hot-reload",
    )
    p.add_argument(
        "--no-browser",
        dest="no_browser",
        action="store_true",
        help="Do not open browser automatically",
    )
    return p.parse_args()


def main() -> int:
    args = _parse()
    from .app import launch
    launch(
        host=args.host,
        port=args.port,
        debug=args.debug,
        open_browser=not args.no_browser,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
