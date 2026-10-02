# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""python -m pycsamt.app.mapview"""

from __future__ import annotations

import argparse
import sys


def _parse() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="pycsamt-mapview",
        description="Launch the pyCSAMT Map View platform.",
    )
    p.add_argument(
        "--host",
        default="127.0.0.1",
        help="Bind address (default 127.0.0.1)",
    )
    p.add_argument(
        "--port",
        type=int,
        default=8770,
        help="HTTP port (default 8770)",
    )
    p.add_argument(
        "--data",
        default=None,
        help="Optional EDI folder to preload (one or more lines).",
    )
    p.add_argument(
        "--pcsf",
        default=None,
        help="Optional .pcsf/.pcsm inversion model to preload.",
    )
    p.add_argument(
        "--state",
        default=None,
        help="Optional scene JSON to open with --pcsf (written by the "
        "desktop's 'Open in Map View').",
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
        help="Do not open a browser automatically",
    )
    return p.parse_args()


def main() -> int:
    args = _parse()
    from .app import launch

    view = None
    state = None
    if args.state:
        import json
        from pathlib import Path

        try:
            state = json.loads(Path(args.state).read_text(encoding="utf-8"))
        except Exception as exc:  # open the model anyway
            print(f"pycsamt-mapview: ignoring --state: {exc}", file=sys.stderr)
    if args.pcsf:
        from pycsamt.map import MapView

        view = MapView.from_pcsf(args.pcsf, fetch_elevation=False)
        elev = (state or {}).get("elevations")
        if elev:  # the topography the desktop scene used
            view = view.with_elevations(elev)
    elif args.data:
        from pycsamt.map import MapView

        view = MapView.from_folder(args.data)

    launch(
        host=args.host,
        port=args.port,
        debug=args.debug,
        open_browser=not args.no_browser,
        view=view,
        state=state,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
