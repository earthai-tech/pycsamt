# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Open the Map View web application from the desktop app.

Map View (``pycsamt-mapview``, Dash) is the full 3-D/map platform:
fence, block, depth-slice and iso views, basemaps, geology, boreholes,
structure.  The desktop's 3-D window offers **Open in Map View**, which
starts it in its own process, seeded with the same ``.pcsf``/``.pcsm``
model, and opens the browser.  Like Agent Master
(:mod:`~pycsamt.app.desktop.agent_master_bridge`), a frozen desktop
build re-executes itself with a sentinel flag, because ``sys.executable``
is not a Python interpreter there.

The launch can also carry the desktop scene (``state``, see
:func:`~pycsamt.app.desktop.controllers.pcsf_scene.mapview_state`): it is
written to a JSON file passed as ``--state``, and Map View opens straight
on the 3-D view with the same mode, colours, depth, topography, overlays,
spin and camera -- nothing to rebuild.

A running Map View cannot be re-seeded with another file, so each launch
uses a free port (from 8770 upwards).
"""

from __future__ import annotations

import json
import os
import socket
import tempfile
import subprocess
import sys
import threading
import time
import webbrowser
from dataclasses import dataclass

MAPVIEW_SERVER_FLAG = "--mapview-server"
DEFAULT_HOST = "127.0.0.1"
DEFAULT_PORT = 8770

__all__ = ["MAPVIEW_SERVER_FLAG", "MapViewLaunch", "free_port",
           "launch_mapview", "mapview_command", "write_state"]


@dataclass
class MapViewLaunch:
    url: str
    process: subprocess.Popen | None = None


def free_port(host: str = DEFAULT_HOST, start: int = DEFAULT_PORT,
              tries: int = 50) -> int:
    """First port from *start* nobody is listening on."""
    for port in range(start, start + tries):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            if s.connect_ex((host, port)) != 0:
                return port
    raise OSError(f"no free port in {start}-{start + tries - 1}")


def mapview_command(pcsf: str | None, host: str, port: int,
                    state_file: str | None = None) -> list[str]:
    head = ([sys.executable, MAPVIEW_SERVER_FLAG]
            if getattr(sys, "frozen", False)
            else [sys.executable, "-m", "pycsamt.app.mapview"])
    cmd = head + ["--host", host, "--port", str(int(port)), "--no-browser"]
    if pcsf:
        cmd += ["--pcsf", str(pcsf)]
    if state_file:
        cmd += ["--state", str(state_file)]
    return cmd


def write_state(state: dict) -> str:
    """Write a Map View seed state to a temporary JSON file."""
    fd, path = tempfile.mkstemp(prefix="pycsamt_mapview_", suffix=".json")
    with os.fdopen(fd, "w", encoding="utf-8") as fh:
        json.dump(state, fh)
    return path


def _open_when_ready(url: str, host: str, port: int,
                     timeout: float = 60.0) -> None:
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            if s.connect_ex((host, port)) == 0:
                webbrowser.open(url)
                return
        time.sleep(0.5)


def launch_mapview(pcsf: str | None = None, *, host: str = DEFAULT_HOST,
                   port: int | None = None, open_browser: bool = True,
                   state: dict | None = None) -> MapViewLaunch:
    """Start Map View (optionally seeded with *pcsf* and the scene
    *state*) and open it."""
    port = free_port(host) if port is None else int(port)
    url = f"http://{host}:{port}"
    env = dict(os.environ)
    env.setdefault("PYTHONDONTWRITEBYTECODE", "1")
    kwargs: dict = {"stdout": subprocess.DEVNULL,
                    "stderr": subprocess.DEVNULL,
                    "stdin": subprocess.DEVNULL, "env": env}
    if os.name == "nt":
        kwargs["creationflags"] = (subprocess.CREATE_NO_WINDOW
                                   | subprocess.CREATE_NEW_PROCESS_GROUP)
    else:
        kwargs["start_new_session"] = True
    state_file = write_state(state) if state else None
    proc = subprocess.Popen(mapview_command(pcsf, host, port, state_file),
                            **kwargs)
    if open_browser:
        threading.Thread(target=_open_when_ready, args=(url, host, port),
                         daemon=True).start()
    return MapViewLaunch(url=url, process=proc)
