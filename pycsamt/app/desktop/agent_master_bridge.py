# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Bridge from the desktop app to the Agent Master web interface."""

from __future__ import annotations

from dataclasses import dataclass
import os
import subprocess
import sys
import threading
import time
from urllib.error import URLError
from urllib.request import urlopen
import webbrowser

DEFAULT_HOST = "127.0.0.1"
DEFAULT_PORT = 8765
_PROCESS: subprocess.Popen | None = None
_LOCK = threading.Lock()


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
) -> AgentMasterLaunch:
    """Start Agent Master if needed and open it in the default browser."""
    global _PROCESS
    url = agent_master_url(host, port)
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
