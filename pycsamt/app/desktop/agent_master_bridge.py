# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Bridge from the desktop app to the Agent Master web interface.

Agent Master always runs as a separate subprocess (a Dash server), never
in-process. When running from source that subprocess is
``sys.executable -m pycsamt.app.agent_master``; a PyInstaller-frozen
``pycsamt-desktop`` build has no real Python interpreter behind
``sys.executable`` for ``-m`` to work against, so the frozen branch instead
re-execs the same frozen binary with ``AGENT_MASTER_SERVER_FLAG`` --
``packaging/pyinstaller/entry_desktop.py`` recognizes that flag and
dispatches to Agent Master's own entry point instead of the desktop app's.
"""

from __future__ import annotations

import os
import subprocess
import sys
import threading
import time
import webbrowser
from dataclasses import dataclass
from urllib.error import URLError
from urllib.request import urlopen

DEFAULT_HOST = "127.0.0.1"
DEFAULT_PORT = 8765
_PROCESS: subprocess.Popen | None = None
_LOCK = threading.Lock()

#: Sentinel CLI flag recognized by ``packaging/pyinstaller/entry_desktop.py``
#: when re-exec'ing a frozen ``pycsamt-desktop`` binary as an Agent Master
#: server subprocess (see ``launch_agent_master``'s frozen-build branch).
AGENT_MASTER_SERVER_FLAG = "--agent-master-server"


@dataclass(frozen=True)
class AgentMasterLaunch:
    """Result of a desktop-to-Agent-Master launch request."""

    url: str
    started: bool
    process: subprocess.Popen | None = None


def agent_master_url(
    host: str = DEFAULT_HOST,
    port: int = DEFAULT_PORT,
) -> str:
    """Return the local Agent Master URL."""
    return f"http://{host}:{int(port)}"


def is_agent_master_running(
    host: str = DEFAULT_HOST,
    port: int = DEFAULT_PORT,
    *,
    timeout: float = 0.45,
) -> bool:
    """Return whether Agent Master appears to answer HTTP requests."""
    try:
        with urlopen(agent_master_url(host, port), timeout=timeout) as resp:
            return 200 <= int(resp.status) < 500
    except (OSError, URLError, TimeoutError, ValueError):
        return False


def launch_agent_master(
    host: str = DEFAULT_HOST,
    port: int = DEFAULT_PORT,
    *,
    open_browser: bool = True,
    handoff: str = "",
) -> AgentMasterLaunch:
    """Start Agent Master if needed and open it in the default browser.

    *handoff* is a token from
    :func:`pycsamt.app.agent_master._handoff.write_handoff`: the page opens
    on ``/?handoff=<token>`` and starts on the desktop's survey (works
    whether the server is already running or not).
    """
    global _PROCESS
    base = agent_master_url(host, port)
    url = f"{base}/?handoff={handoff}" if handoff else base
    with _LOCK:
        if is_agent_master_running(host, port):
            if open_browser:
                webbrowser.open(url)
            return AgentMasterLaunch(url=url, started=False)

        if _PROCESS is not None and _PROCESS.poll() is None:
            if open_browser:
                threading.Thread(
                    target=_open_when_ready,
                    args=(host, port, url),
                    daemon=True,
                ).start()
            return AgentMasterLaunch(
                url=url,
                started=False,
                process=_PROCESS,
            )

        if getattr(sys, "frozen", False):
            # A PyInstaller-frozen ``pycsamt-desktop`` binary has no real
            # Python interpreter behind ``sys.executable`` -- it *is* the
            # frozen app, so ``-m pycsamt.app.agent_master`` cannot work
            # the way it does when running from source. Re-exec the same
            # frozen binary with a sentinel flag instead; the packaged
            # entry point (packaging/pyinstaller/entry_desktop.py) checks
            # for it before launching the desktop app itself and
            # dispatches to Agent Master's own server entry point.
            cmd = [
                sys.executable,
                AGENT_MASTER_SERVER_FLAG,
                "--host",
                host,
                "--port",
                str(int(port)),
                "--no-browser",
            ]
        else:
            cmd = [
                sys.executable,
                "-m",
                "pycsamt.app.agent_master",
                "--host",
                host,
                "--port",
                str(int(port)),
                "--no-browser",
            ]
        env = dict(os.environ)
        env.setdefault("PYTHONDONTWRITEBYTECODE", "1")
        kwargs: dict[str, object] = {
            "stdout": subprocess.DEVNULL,
            "stderr": subprocess.DEVNULL,
            "stdin": subprocess.DEVNULL,
            "env": env,
        }
        if os.name == "nt":
            kwargs["creationflags"] = (
                subprocess.CREATE_NO_WINDOW
                | subprocess.CREATE_NEW_PROCESS_GROUP
            )
        else:
            kwargs["start_new_session"] = True

        proc = subprocess.Popen(cmd, **kwargs)
        _PROCESS = proc

    if open_browser:
        threading.Thread(
            target=_open_when_ready,
            args=(host, port, url),
            daemon=True,
        ).start()
    return AgentMasterLaunch(url=url, started=True, process=proc)



def _open_when_ready(
    host: str,
    port: int,
    url: str,
    *,
    attempts: int = 40,
    interval: float = 0.25,
) -> None:
    """Open *url* after the local Dash server starts responding."""
    for _ in range(attempts):
        if is_agent_master_running(host, port):
            webbrowser.open(url)
            return
        time.sleep(interval)
    webbrowser.open(url)
