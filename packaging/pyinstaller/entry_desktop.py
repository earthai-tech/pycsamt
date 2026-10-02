# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""PyInstaller entry script for pycsamt-desktop.

A separate, tiny script (rather than pointing PyInstaller at
``pycsamt/app/desktop/__main__.py`` directly) so the frozen build's entry
point is decoupled from the package layout and stays a stable target for
``pycsamt_desktop.spec`` -- same rationale as the converter app's
``entry_converter.py``.

Dual-purpose: the desktop app also launches Agent Master (a separate Dash
server) as a subprocess. When running from source that subprocess is
``sys.executable -m pycsamt.app.agent_master``; a frozen build has no real
Python interpreter behind ``sys.executable`` for ``-m`` to work against, so
``agent_master_bridge.launch_agent_master()`` instead re-execs *this same*
frozen binary with the ``--agent-master-server`` sentinel flag. This entry
script checks for that flag first and dispatches to Agent Master's own
entry point instead of the desktop UI's -- see
``pycsamt/app/desktop/agent_master_bridge.py``'s module docstring for the
full picture.
"""

from __future__ import annotations

import os

# Belt-and-suspenders: the same guard also lives in
# pycsamt/app/desktop/__main__.py and pycsamt/app/agent_master/__main__.py
# (both reachable below); set here too since this is the frozen build's
# actual process entry point, before either of those modules is imported.
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

import sys

# The windowed build (``console=False``) starts with ``sys.stdout`` and
# ``sys.stderr`` set to None. Library code that writes progress there --
# tqdm bars, ``stream.write(...)`` -- then fails (e.g. Occam1D refused to
# build its inputs: "stream must provide a callable write method"). Give
# both a real sink before anything is imported.
for _name in ("stdout", "stderr"):
    if getattr(sys, _name) is None:
        setattr(sys, _name, open(os.devnull, "w", encoding="utf-8"))


def main() -> None:
    from pycsamt.app.desktop.agent_master_bridge import (
        AGENT_MASTER_SERVER_FLAG,
    )

    if AGENT_MASTER_SERVER_FLAG in sys.argv[1:]:
        # Strip the sentinel so agent_master's own argparse-based CLI
        # (pycsamt/app/agent_master/__main__.py) sees only its own flags.
        sys.argv = [sys.argv[0]] + [
            a for a in sys.argv[1:] if a != AGENT_MASTER_SERVER_FLAG
        ]
        from pycsamt.app.agent_master.__main__ import (
            main as agent_master_main,
        )

        sys.exit(agent_master_main())

    from pycsamt.app.desktop.mapview_bridge import MAPVIEW_SERVER_FLAG

    if MAPVIEW_SERVER_FLAG in sys.argv[1:]:
        # "Open in Map View" from the 3-D window (see mapview_bridge).
        sys.argv = [sys.argv[0]] + [
            a for a in sys.argv[1:] if a != MAPVIEW_SERVER_FLAG
        ]
        from pycsamt.app.mapview.__main__ import main as mapview_main

        sys.exit(mapview_main())

    from pycsamt.app.desktop.__main__ import main as desktop_main

    desktop_main()


if __name__ == "__main__":
    main()
