# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Inversion Studio library: recent runs and settings presets (Qt-free).

* Recent runs -- ``~/.pycsamt/inversion_runs.json``; one entry per run
  folder (the latest run in a folder replaces the earlier entry), listed
  newest first and only while the folder still exists.
* Presets -- a few built-in starting points per engine plus the user's own
  (``~/.pycsamt/inversion_presets.json``).  A preset only holds the values
  it changes; the rest keep the engine defaults.
"""

from __future__ import annotations

import json
import time
from dataclasses import dataclass
from pathlib import Path

__all__ = [
    "BUILTIN_PRESETS",
    "Preset",
    "delete_preset",
    "presets",
    "recent_runs",
    "record_run",
    "runs_path",
    "save_preset",
    "user_presets_path",
]

_MAX_RUNS = 50


def runs_path() -> Path:
    return Path.home() / ".pycsamt" / "inversion_runs.json"


def user_presets_path() -> Path:
    return Path.home() / ".pycsamt" / "inversion_presets.json"


def _read(path: Path, default):
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return default


def _write(path: Path, data) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2), encoding="utf-8")


# ── Recent runs ───────────────────────────────────────────────────────────────


def record_run(engine: str, workdir, *, stations: int = 0,
               final_rms: float | None = None, status: str = "done",
               label: str = "", path: Path | None = None) -> dict:
    """Add (or refresh) the run in *workdir* at the top of the list."""
    p = path or runs_path()
    folder = str(Path(workdir).resolve())
    runs = [r for r in _read(p, []) if r.get("workdir") != folder]
    entry = {
        "engine": engine, "workdir": folder, "stations": int(stations),
        "final_rms": None if final_rms is None else float(final_rms),
        "status": status, "label": label or Path(folder).name,
        "when": time.strftime("%Y-%m-%d %H:%M"),
    }
    runs.insert(0, entry)
    _write(p, runs[:_MAX_RUNS])
    return entry


def recent_runs(limit: int = 30, *, path: Path | None = None) -> list[dict]:
    """Newest first, skipping folders that no longer exist."""
    runs = [r for r in _read(path or runs_path(), [])
            if isinstance(r, dict) and Path(r.get("workdir", "")).is_dir()]
    return runs[:limit]


# ── Presets ───────────────────────────────────────────────────────────────────


@dataclass(frozen=True)
class Preset:
    name: str
    values: dict
    builtin: bool = True
    description: str = ""


BUILTIN_PRESETS: dict[str, list[Preset]] = {
    "occam1d": [
        Preset("Quick look", {"n_layers": 25, "max_iterations": 15},
               description="Few layers and iterations — a first check."),
        Preset("Standard", {"n_layers": 40, "max_iterations": 30},
               description="The Occam1D defaults."),
        Preset("Deep crust", {"n_layers": 60, "depth_max": 50000.0,
                              "first_thickness": 20.0,
                              "max_iterations": 40},
               description="Long-period MT down to 50 km."),
    ],
    "occam2d": [
        Preset("Fast preview", {"n_layers": 20, "max_iterations": 10,
                                "target_misfit": 2.0},
               description="Coarse mesh, 10 iterations."),
        Preset("Standard", {"n_layers": 30, "max_iterations": 100},
               description="The Occam2D defaults."),
        Preset("Fine", {"n_layers": 45, "cell_size_vertical_top": 5.0,
                        "depth_scale": 1.12, "max_iterations": 150},
               description="Thinner layers for shallow targets."),
    ],
    "modem2d": [
        Preset("Test run", {"max_iterations": 5, "nz_2d": 30},
               description="Five iterations to check the setup."),
        Preset("Standard", {"max_iterations": 100},
               description="The ModEM defaults."),
    ],
    "modem3d": [
        Preset("Coarse test", {"nx": 12, "ny": 12, "nz": 20,
                               "max_iterations": 5},
               description="Small grid, five iterations."),
        Preset("Production", {"nx": 30, "ny": 30, "nz": 40,
                              "max_iterations": 150, "use_mpi": True},
               description="Finer grid with MPI."),
    ],
    "mare2dem": [
        Preset("Quick test", {"max_iterations": 3, "n_procs": 4},
               description="Three iterations to check the setup."),
        Preset("Standard", {"max_iterations": 30},
               description="Typical MT inversion."),
    ],
}


def presets(engine: str, *, path: Path | None = None) -> list[Preset]:
    """Built-in presets, then the user's own, for *engine*."""
    user = _read(path or user_presets_path(), {}).get(engine, {})
    out = list(BUILTIN_PRESETS.get(engine, []))
    for name, values in sorted(user.items()):
        out.append(Preset(name, dict(values), builtin=False))
    return out


def save_preset(engine: str, name: str, values: dict, *,
                path: Path | None = None) -> None:
    p = path or user_presets_path()
    data = _read(p, {})
    data.setdefault(engine, {})[name] = values
    _write(p, data)


def delete_preset(engine: str, name: str, *, path: Path | None = None) -> None:
    p = path or user_presets_path()
    data = _read(p, {})
    if data.get(engine, {}).pop(name, None) is not None:
        _write(p, data)
