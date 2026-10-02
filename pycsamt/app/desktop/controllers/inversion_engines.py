# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Inversion engines for the desktop Inversion Studio (Qt-free).

One :class:`Engine` per classical inversion code, all with the same shape,
so the window never special-cases an engine:

========  =========================================================
fields    the settings the user can edit (the window builds forms)
build     write the solver's input files into a work folder and return
          a :class:`BuildInfo` (the mesh preview is drawn from it)
task      a callable that runs the solver, reporting console lines,
          iterations and stages through a :class:`RunReporter`
history   RMS per iteration parsed from the work folder's log (polled
          while an external solver runs)
detect /  recognise and open a finished run folder (from this app, the
load      library, or the solver itself) for the results view
views     the result plots an engine offers; ``render`` draws one or
          raises :class:`PlotUnavailable` with the reason
========  =========================================================

Engines: Occam1D (pure Python, per-station), Occam2D, ModEM 2-D, ModEM 3-D
and MARE2DEM.
"""

from __future__ import annotations

import math
import re
import time
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from pycsamt.app.desktop.controllers.correction_views import PlotUnavailable

__all__ = [
    "ENGINES",
    "BuildInfo",
    "Engine",
    "Field",
    "LoadedRun",
    "PlotUnavailable",
    "RunReporter",
    "detect_engine",
    "engine",
    "site_list",
]


# ── Settings fields ───────────────────────────────────────────────────────────


@dataclass(frozen=True)
class Field:
    """One editable setting; the window turns it into a widget."""

    key: str
    label: str
    kind: str  # "int" | "float" | "bool" | "choice" | "text"
    default: Any
    lo: float | None = None
    hi: float | None = None
    step: float | None = None
    decimals: int = 2
    choices: tuple = ()  # (value, label) pairs for "choice"
    help: str = ""
    section: str = "Settings"  # "Data" | "Mesh" | "Settings" | "Run"
    advanced: bool = False
    unit: str = ""
    auto_zero: bool = False  # 0 means "automatic"


def defaults(fields: list[Field]) -> dict:
    return {f.key: f.default for f in fields}


# ── Run plumbing ──────────────────────────────────────────────────────────────


class RunReporter:
    """What a running task reports back (the worker turns it into signals).

    The base class only records; the desktop worker overrides the hooks.
    """

    def __init__(self) -> None:
        self.lines: list[str] = []
        self.history: list[tuple[str, int, float]] = []
        self.stage_text = ""
        self._cancelled = False

    def log(self, line: str) -> None:
        self.lines.append(line)

    def stage(self, text: str) -> None:
        self.stage_text = text

    def iteration(self, n: int, rms: float, label: str = "") -> None:
        """One accepted iteration (``label`` = station for Occam1D)."""
        self.history.append((label, int(n), float(rms)))

    def cancel(self) -> None:
        self._cancelled = True

    def cancelled(self) -> bool:
        return self._cancelled


@dataclass
class BuildInfo:
    """Result of writing an engine's input files."""

    engine: str
    workdir: Path
    files: dict[str, Path] = field(default_factory=dict)
    summary: list[tuple[str, str]] = field(default_factory=list)
    mesh: dict = field(default_factory=dict)
    values: dict = field(default_factory=dict)
    stations: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)


@dataclass
class LoadedRun:
    """A finished (or partial) run opened for the results view."""

    engine: str
    path: Path
    result: Any = None
    rms: list[tuple[int, float]] = field(default_factory=list)
    iterations: list[int] = field(default_factory=list)
    iteration: int | None = None
    stations: list[str] = field(default_factory=list)
    station: str | None = None
    summary: list[tuple[str, str]] = field(default_factory=list)
    extra: dict = field(default_factory=dict)

    @property
    def final_rms(self) -> float | None:
        return self.rms[-1][1] if self.rms else None


# ── Site helpers ──────────────────────────────────────────────────────────────


def site_list(sites) -> list:
    """Site objects of a Sites collection (iteration yields ``Site``)."""
    if sites is None:
        return []
    try:
        return list(sites)
    except TypeError:
        return []


def site_names(sites) -> list[str]:
    return [str(getattr(s, "name", "")) for s in site_list(sites)]


def select_sites(sites, names: list[str] | None):
    if sites is None or not names:
        return sites
    return sites.select(names=list(names))


def _figure_of(obj, fig):
    """Figure produced by a library plot call (``obj`` may be a Figure,
    an Axes, a tuple, or ``None`` when it drew into *fig*)."""
    from matplotlib.figure import Figure

    if isinstance(obj, Figure):
        return obj
    if isinstance(obj, tuple):
        for item in obj:
            if isinstance(item, Figure):
                return item
            if hasattr(item, "figure") and isinstance(item.figure, Figure):
                return item.figure
    if hasattr(obj, "figure") and isinstance(getattr(obj, "figure"), Figure):
        return obj.figure
    return fig


def _constrained(fig) -> None:
    """Titles/labels never overlap (mesh previews stack two panels)."""
    try:
        fig.set_layout_engine("constrained")
    except Exception:  # matplotlib < 3.6
        fig.set_constrained_layout(True)


def _need(cond, title, reason, guidance=""):
    if not cond:
        raise PlotUnavailable(title, reason, guidance)


def _grid_axes(ax, x_nodes, z_nodes, *, stations=None, xlim=None,
               zmax=None, color="#4c6a92", title=""):
    """Draw mesh lines (x_nodes, z_nodes positive down) on *ax*."""
    x = np.asarray(x_nodes, float)
    z = np.asarray(z_nodes, float)
    lo, hi = (xlim if xlim is not None else (x.min(), x.max()))
    zb = zmax if zmax is not None else z.max()
    xs = x[(x >= lo) & (x <= hi)]
    zs = z[z <= zb]
    ax.vlines(xs, zs.min() if zs.size else 0, zb, color=color, lw=0.35)
    ax.hlines(zs, lo, hi, color=color, lw=0.35)
    if stations is not None and len(stations):
        ax.plot(stations, np.full(len(stations), zs.min() if zs.size else 0),
                "v", color="#c92a2a", ms=7, mec="k", mew=0.5,
                clip_on=False, zorder=5, label="stations")
    ax.set_xlim(lo, hi)
    ax.set_ylim(zb, min(zs.min(), 0) if zs.size else 0)
    ax.set_xlabel("Distance (m)")
    ax.set_ylabel("Depth (m)")
    if title:
        ax.set_title(title, fontsize=10)


def _core_centre(nodes) -> float:
    """Centre of the finest (core) cells of a padded grid axis."""
    x = np.asarray(nodes, float)
    w = np.diff(x)
    if w.size == 0:
        return 0.0
    core = np.flatnonzero(w <= w.min() * 1.001)
    return 0.5 * (x[core[0]] + x[core[-1] + 1])


def _align(stations, nodes):
    """Station positions shifted so their span is centred on the core
    cells -- where every builder places them, whatever the file frame."""
    if stations is None or not np.size(stations):
        return stations
    s = np.asarray(stations, float)
    return s - 0.5 * (s.min() + s.max()) + _core_centre(nodes)


def _core_depth(z_nodes, stations, min_layers=10) -> float:
    """Depth shown in a "core region" panel: ~2x the station span, never
    fewer than *min_layers* layers."""
    z = np.asarray(z_nodes, float)
    k = min(min_layers, z.size - 1)
    span = float(np.ptp(stations)) if stations is not None and         np.size(stations) > 1 else 0.0
    return float(min(z[-1], max(z[k], 2.0 * span)))


def _pad_window(stations, x_nodes, frac=0.25):
    s = np.asarray(stations, float)
    if s.size == 0:
        return None
    span = max(s.max() - s.min(), 1.0)
    x = np.asarray(x_nodes, float)
    return (max(x.min(), s.min() - frac * span),
            min(x.max(), s.max() + frac * span))


# ── Engine base ───────────────────────────────────────────────────────────────


class Engine:
    """Common interface; see the module docstring."""

    key = ""
    label = ""
    dim = "2D"  # "1D" | "2D" | "3D"
    binary_key: str | None = None  # Solver Builder key; None = built in
    description = ""

    # --- settings -------------------------------------------------------
    def fields(self) -> list[Field]:
        return []

    def defaults(self) -> dict:
        return defaults(self.fields())

    # --- build ----------------------------------------------------------
    def check_sites(self, sites) -> str | None:
        n = len(site_list(sites))
        if n == 0:
            return "Load EDI data first (File ▸ Open)."
        return None

    def build(self, sites, workdir, values: dict) -> BuildInfo:
        raise NotImplementedError

    def plot_mesh(self, fig, build: BuildInfo) -> None:
        raise PlotUnavailable("No mesh preview", f"{self.label} has none.")

    # --- run ------------------------------------------------------------
    def make_task(self, build: BuildInfo, values: dict,
                  binary: str | None) -> Callable[[RunReporter], Any]:
        raise NotImplementedError

    def history(self, workdir) -> list[tuple[int, float]]:
        return []

    # --- results --------------------------------------------------------
    def detect(self, path) -> bool:
        return False

    def load(self, path, iteration=None, station=None) -> LoadedRun:
        raise NotImplementedError

    def views(self, run: LoadedRun) -> list[tuple[str, str]]:
        return []

    def render(self, key: str, fig, run: LoadedRun):
        raise PlotUnavailable("Unknown view", key)


def _rms_view(fig, run: LoadedRun, target: float | None = None):
    _need(run.rms, "No convergence history",
          "This run folder has no readable iteration log.",
          "The solver writes it while it iterates; an interrupted run may "
          "not have one.")
    ax = fig.add_subplot(111)
    it = [a for a, _ in run.rms]
    rms = [b for _, b in run.rms]
    ax.plot(it, rms, "o-", color="#1864ab", lw=1.6, ms=4)
    if target:
        ax.axhline(target, color="#2a7f3f", ls="--", lw=1, label="target")
        ax.legend(loc="upper right", fontsize=8)
    ax.set_xlabel("Iteration")
    ax.set_ylabel("RMS misfit")
    ax.set_title(f"Convergence — final RMS {rms[-1]:.3f}", fontsize=10)
    ax.grid(True, alpha=0.3)
    return fig


# ══════════════════════════════════════════════════════════════════════════════
# Occam1D — pure Python, one independent inversion per station
# ══════════════════════════════════════════════════════════════════════════════

_O1D_MANIFEST = "occam1d_manifest.json"
_O1D_RESTART = "native-restart.json"


class Occam1DEngine(Engine):
    key = "occam1d"
    label = "Occam1D"
    dim = "1D"
    binary_key = None
    description = ("Smooth layered-earth inversion of each station "
                   "(Constable et al., 1987). Built in — no binary needed.")

    def fields(self):
        return [
            Field("mode", "Response", "choice", "auto",
                  choices=(("auto", "Auto"), ("determinant", "Determinant"),
                           ("xy", "XY (TE)"), ("yx", "YX (TM)")),
                  section="Data",
                  help="Impedance used for each 1-D sounding. Auto uses "
                       "the determinant, or XY for scalar (CSAMT) data "
                       "that has no YX component."),
            Field("error_floor_rho", "Error floor ρ", "float", 0.05,
                  0.001, 1.0, 0.01, 3, section="Data",
                  help="Relative apparent-resistivity error floor."),
            Field("error_floor_phase", "Error floor φ", "float", 2.0,
                  0.01, 30.0, 0.5, 2, section="Data", unit="°"),
            Field("n_layers", "Layers", "int", 40, 5, 200, section="Mesh"),
            Field("first_thickness", "First layer", "float", 5.0, 0.1,
                  1e4, 1.0, 1, section="Mesh", unit="m"),
            Field("depth_max", "Max depth", "float", 5000.0, 10.0, 1e6,
                  100.0, 0, section="Mesh", unit="m"),
            Field("starting_resistivity", "Starting ρ", "float", 100.0,
                  0.1, 1e6, 10.0, 1, unit="Ω·m"),
            Field("max_iterations", "Max iterations", "int", 30, 1, 500),
            Field("target_misfit", "Target RMS", "float", 1.0, 0.1, 50.0,
                  0.1, 2),
            Field("roughness_type", "Roughness", "choice", 1,
                  choices=((1, "First derivative"),
                           (2, "Second derivative"))),
            Field("lagrange_start", "Start log₁₀ μ", "float", 5.0, -5.0,
                  15.0, 0.5, 2, advanced=True),
            Field("stepsize_cut_count", "Step cuts", "int", 8, 1, 50,
                  advanced=True),
        ]

    def check_sites(self, sites):
        return super().check_sites(sites)

    def _config(self, values, sites=None):
        from pycsamt.models.occam1d import Occam1DConfig

        keys = ("mode", "error_floor_rho", "error_floor_phase", "n_layers",
                "first_thickness", "depth_max", "starting_resistivity",
                "max_iterations", "target_misfit", "roughness_type",
                "lagrange_start", "stepsize_cut_count")
        kw = {k: values[k] for k in keys if k in values}
        if kw.get("mode", "auto") == "auto":
            kw["mode"] = _auto_1d_mode(sites)
        if values.get("freq_min"):
            kw["freq_min"] = float(values["freq_min"])
        if values.get("freq_max"):
            kw["freq_max"] = float(values["freq_max"])
        return Occam1DConfig(**kw)

    def build(self, sites, workdir, values):
        from pycsamt.models.occam1d import Occam1DBatch

        wd = Path(workdir)
        cfg = self._config(values, sites)
        batch = Occam1DBatch(sites, wd, config=cfg, verbose=0)
        batch.build_all(continue_on_error=True)
        if not batch.builders:
            first = next(iter(batch.failures.values()), "no usable data")
            raise ValueError(
                f"No station could be prepared ({first}). Check the "
                "frequency window and the Response choice — scalar CSAMT "
                "data has only the XY component.")
        dirs = {b.data.station: b.workdir for b in batch.builders}
        info = BuildInfo(self.key, wd, values=dict(values),
                         stations=list(dirs))
        info.files = {name: d / cfg.data_file for name, d in dirs.items()}
        info.warnings = [f"{k}: {v}" for k, v in batch.failures.items()]
        depth = np.asarray(batch.builders[0].model.depth, float) \
            if batch.builders else np.array([])
        freqs = []
        for b in batch.builders:
            freqs.extend(np.asarray(getattr(b.data, "frequency",
                                            getattr(b.data, "frequencies",
                                                    [])), float).tolist())
        info.mesh = {"depth": depth, "dirs": dirs,
                     "freqs": np.asarray(freqs, float),
                     "rho0": float(values.get("starting_resistivity", 100))}
        info.summary = [
            ("Stations", f"{len(dirs)} built"
             + (f", {len(batch.failures)} failed" if batch.failures else "")),
            ("Layers", str(depth.size)),
            ("Depth range", f"{depth[1] if depth.size > 1 else 0:.0f} – "
                            f"{depth[-1] if depth.size else 0:,.0f} m"),
            ("Folder", str(wd)),
        ]
        return info

    def plot_mesh(self, fig, build):
        _constrained(fig)
        depth = np.asarray(build.mesh.get("depth", []), float)
        _need(depth.size > 1, "No layers", "The model has no layers.")
        ax = fig.add_subplot(121)
        tops = depth[1:]
        ax.hlines(tops, 0, 1, color="#4c6a92", lw=0.6)
        ax.set_yscale("log")
        ax.set_ylim(tops.max() * 1.1, max(tops.min() * 0.9, 1e-2))
        ax.set_xticks([])
        ax.set_ylabel("Depth (m)")
        ax.set_title(f"{depth.size} layers", fontsize=10)
        f = build.mesh.get("freqs")
        rho0 = build.mesh.get("rho0", 100.0)
        if f is not None and np.size(f):
            f = np.asarray(f, float)
            f = f[f > 0]
            if f.size:
                d_lo = 503.0 * math.sqrt(rho0 / f.max())
                d_hi = 503.0 * math.sqrt(rho0 / f.min())
                ax.axhspan(d_lo, d_hi, color="#3b5bdb", alpha=0.10,
                           label=f"skin depth @ {rho0:g} Ω·m")
                ax.legend(loc="lower right", fontsize=7)
        ax2 = fig.add_subplot(122)
        thick = np.diff(depth)
        ax2.plot(thick, depth[1:], "o-", ms=2.5, color="#1864ab")
        ax2.set_xscale("log")
        ax2.set_yscale("log")
        ax2.invert_yaxis()
        ax2.set_xlabel("Thickness (m)")
        ax2.set_title("Layer thickness", fontsize=10)
        ax2.grid(True, which="both", alpha=0.25)
        fig.suptitle(f"Occam1D layer mesh — {len(build.stations)} stations",
                     fontsize=11)

    def make_task(self, build, values, binary):
        dirs = dict(build.mesh["dirs"])

        def task(rep: RunReporter):
            done, failed = {}, {}
            for i, (name, wd) in enumerate(dirs.items()):
                if rep.cancelled():
                    break
                rep.stage(f"Station {i + 1}/{len(dirs)} · {name}")
                rep.log(f"── {name} ──")
                try:
                    done[name] = _o1d_invert_station(wd, rep, name)
                except InterruptedError:
                    rep.log(f"{name}: cancelled")
                    break
                except Exception as exc:  # keep going, like invert_all
                    failed[name] = f"{type(exc).__name__}: {exc}"
                    rep.log(f"{name}: FAILED — {failed[name]}")
            if rep.cancelled():
                raise InterruptedError("Run cancelled by the user.")
            if not done:
                raise RuntimeError("No station converged to a result: "
                                   + "; ".join(failed.values()))
            return {"workdir": build.workdir, "done": done,
                    "failed": failed}

        return task

    # --- results --------------------------------------------------------
    @staticmethod
    def _station_dirs(path: Path) -> list[Path]:
        p = Path(path)
        if (p / _O1D_MANIFEST).is_file():
            return [p]
        return sorted(d for d in p.iterdir()
                      if d.is_dir() and (d / _O1D_MANIFEST).is_file()) \
            if p.is_dir() else []

    def detect(self, path):
        return bool(self._station_dirs(Path(path)))

    def load(self, path, iteration=None, station=None):
        dirs = {d.name: d for d in self._station_dirs(Path(path))}
        _need(dirs, "Not an Occam1D run", f"No station folders in {path}.")
        stations = list(dirs)
        name = station if station in dirs else stations[0]
        st = _o1d_load_station(dirs[name])
        run = LoadedRun(self.key, Path(path), result=st,
                        stations=stations, station=name)
        its = st["result"].iterations if st["result"] else ()
        run.rms = [(it.number, float(it.rms)) for it in its]
        run.iterations = [it.number for it in its]
        run.iteration = (iteration if iteration in run.iterations
                         else (run.iterations[-1] if run.iterations
                               else None))
        run.extra = {"dirs": dirs, "target": st["target"]}
        run.summary = [
            ("Stations", str(len(stations))),
            ("Station", name),
            ("Iterations", str(max(len(run.iterations) - 1, 0))),
            ("Final RMS", f"{run.final_rms:.3f}" if run.final_rms else "–"),
            ("Status", st["status"]),
        ]
        return run

    def views(self, run):
        return [("model", "Model"), ("response", "Data fit"),
                ("convergence", "Convergence"), ("section", "Stitched 2-D"),
                ("summary", "Summary")]

    def render(self, key, fig, run):
        st = run.result
        if key == "section":
            return _o1d_section(fig, run)
        _need(st and st["result"] is not None, "No inversion history",
              f"Station {run.station} has input files but no "
              f"{_O1D_RESTART} — it has not been inverted yet.",
              "Run it from Build & Run.")
        res = _o1d_slice(st["result"], run.iteration)
        inv = st["inversion"]
        from pycsamt.api.occam1d import resolve_occam1d_style

        style = resolve_occam1d_style(None)
        makers = {
            "model": inv._plot_native_model,
            "response": inv._plot_native_response,
            "convergence": inv._plot_native_convergence,
            "summary": inv._plot_native_summary,
        }
        _need(key in makers, "Unknown view", key)
        return makers[key](res, style)


def _auto_1d_mode(sites) -> str:
    """'determinant' when every site has both off-diagonal impedances,
    'xy' for scalar data (e.g. CSAMT with Zxy only)."""
    for s in site_list(sites):
        try:
            z = np.asarray(s.z)
            if z.ndim == 3 and not np.isfinite(z[:, 1, 0]).any():
                return "xy"
        except Exception:
            continue
    return "determinant"


def _o1d_objects(wd: Path):
    import json

    from pycsamt.models.occam1d import (
        Occam1DConfig,
        Occam1DData,
        Occam1DInversion,
        Occam1DModel,
        Occam1DStartup,
    )

    manifest = json.loads((wd / _O1D_MANIFEST).read_text(encoding="utf8"))
    cfg = Occam1DConfig(**manifest["config"])
    data = Occam1DData.read(wd / cfg.data_file)
    model = Occam1DModel.read(wd / cfg.model_file)
    startup = Occam1DStartup.read(wd / cfg.startup_file)
    inv = Occam1DInversion(data, model, config=cfg, startup=startup,
                           verbose=0)
    return cfg, data, inv


def _o1d_invert_station(wd: Path, rep: RunReporter, name: str) -> dict:
    """One station, like the library's ``_invert_one_station`` but with a
    live per-iteration callback and cooperative cancel."""
    cfg, data, inv = _o1d_objects(Path(wd))

    def on_iter(it):
        rep.iteration(it.number, float(it.rms), name)
        rep.log(f"{name}  iter {it.number:3d}  RMS {it.rms:8.4f}  "
                f"roughness {it.roughness:10.4g}  μ {it.multiplier:.3g}")

    result = inv.run(callback=on_iter, cancel=rep.cancelled)
    inv.restart(result).write(Path(wd) / _O1D_RESTART)
    inv.export_result(Path(wd) / "model-text", result)
    rep.log(f"{name}: {result.convergence.value} — {result.message} "
            f"(RMS {result.final.rms:.3f})")
    return {"rms": float(result.final.rms),
            "converged": bool(result.converged),
            "iterations": int(result.n_iterations)}


def _o1d_load_station(wd: Path) -> dict:
    from pycsamt.models.occam1d import Occam1DInversionResult, Occam1DRestart

    cfg, data, inv = _o1d_objects(wd)
    out = {"inversion": inv, "data": data, "result": None,
           "target": float(cfg.target_misfit), "status": "not inverted",
           "depth": np.asarray(inv.model.depth, float)}
    rp = wd / _O1D_RESTART
    if rp.is_file():
        rs = Occam1DRestart.read(rp)
        res = Occam1DInversionResult(
            iterations=rs.iterations,
            convergence=rs.previous_convergence,
            target_rms=rs.target_rms,
            message=rs.previous_message or "",
            rejected_candidates=rs.rejected_candidates,
            failed_iterations=rs.failed_iterations,
        )
        out["result"] = res
        out["status"] = (rs.previous_convergence.value
                         if rs.previous_convergence else "finished")
    return out


def _o1d_slice(result, iteration):
    """The result truncated at *iteration* (history is contiguous)."""
    if iteration is None or iteration >= result.iterations[-1].number:
        return result
    from dataclasses import replace

    keep = tuple(it for it in result.iterations if it.number <= iteration)
    return replace(result, iterations=keep)


def _o1d_section(fig, run: LoadedRun):
    """Final models of every station side by side (station order)."""
    dirs = run.extra.get("dirs", {})
    cols, names, depth = [], [], None
    for name, wd in dirs.items():
        try:
            st = _o1d_load_station(wd)
        except Exception:
            continue
        if st["result"] is None:
            continue
        depth = st["depth"] if depth is None else depth
        if st["depth"].shape != depth.shape:
            continue
        cols.append(np.log10(st["result"].final.resistivity))
        names.append(name)
    _need(len(cols) >= 2, "Not enough inverted stations",
          f"{len(cols)} station(s) have results; a stitched section needs "
          "at least two.")
    img = np.column_stack(cols)
    ax = fig.add_subplot(111)
    edges_z = np.append(depth, depth[-1] * 1.25)
    edges_z[0] = max(depth[1] * 0.5, 1e-1) if depth.size > 1 else 1.0
    x = np.arange(len(names) + 1) - 0.5
    mesh = ax.pcolormesh(x, edges_z, img, cmap="jet_r", shading="flat",
                         vmin=np.nanpercentile(img, 2),
                         vmax=np.nanpercentile(img, 98))
    ax.set_yscale("log")
    ax.invert_yaxis()
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names, rotation=60, fontsize=7)
    ax.set_ylabel("Depth (m)")
    ax.set_title("Occam1D models stitched along the profile", fontsize=10)
    cb = fig.colorbar(mesh, ax=ax, pad=0.02)
    cb.set_label("log₁₀ ρ (Ω·m)")
    return fig


# ══════════════════════════════════════════════════════════════════════════════
# Occam2D
# ══════════════════════════════════════════════════════════════════════════════


class Occam2DEngine(Engine):
    key = "occam2d"
    label = "Occam2D"
    dim = "2D"
    binary_key = "occam2d"
    description = ("Smooth 2-D TE/TM inversion (deGroot-Hedlin & "
                   "Constable, 1990) on a regular finite-element mesh.")

    def fields(self):
        return [
            Field("modes", "Modes", "choice", "TE+TM",
                  choices=(("TE+TM", "TE + TM"), ("TE", "TE only"),
                           ("TM", "TM only")), section="Data"),
            Field("error_floor_rho", "Error floor ρ", "float", 0.05, 0.001,
                  1.0, 0.01, 3, section="Data"),
            Field("error_floor_phase", "Error floor φ", "float", 0.5, 0.01,
                  30.0, 0.1, 2, section="Data", unit="°"),
            Field("n_layers", "Layers", "int", 30, 5, 200, section="Mesh"),
            Field("max_depth", "Max depth", "float", 1500.0, 10.0, 1e6,
                  100.0, 0, section="Mesh", unit="m"),
            Field("cell_size_horizontal", "Cell width", "float", 100.0, 1.0,
                  1e5, 10.0, 0, section="Mesh", unit="m"),
            Field("cell_size_vertical_top", "First layer", "float", 10.0,
                  0.1, 1e4, 1.0, 1, section="Mesh", unit="m"),
            Field("depth_scale", "Depth growth", "float", 1.2, 1.01, 2.0,
                  0.01, 2, section="Mesh"),
            Field("n_airlayers", "Air layers", "int", 5, 1, 30,
                  section="Mesh", advanced=True),
            Field("n_padding_x", "Padding cells", "int", 7, 1, 30,
                  section="Mesh", advanced=True),
            Field("initial_rho", "Starting ρ", "float", 100.0, 0.1, 1e6,
                  10.0, 1, unit="Ω·m"),
            Field("max_iterations", "Max iterations", "int", 100, 1, 1000),
            Field("target_misfit", "Target RMS", "float", 1.0, 0.1, 50.0,
                  0.1, 2),
            Field("roughness_type", "Roughness", "choice", 1,
                  choices=((1, "Standard"), (2, "Deep-weighted"))),
            Field("lagrange_start", "Start log₁₀ μ", "float", 5.0, -5.0,
                  15.0, 0.5, 2, advanced=True),
            Field("diagonal_penalties", "Diagonal penalties", "bool",
                  False, advanced=True),
            Field("stepsize_cut_count", "Step cuts", "int", 8, 1, 50,
                  advanced=True),
        ]

    def check_sites(self, sites):
        n = len(site_list(sites))
        if n < 2:
            return "A 2-D inversion needs at least two stations on a line."
        return None

    def _config(self, values):
        from pycsamt.models.occam2d import OccamConfig

        modes = {"TE+TM": ["TE", "TM"], "TE": ["TE"], "TM": ["TM"]}[
            values.get("modes", "TE+TM")]
        keys = ("error_floor_rho", "error_floor_phase", "n_layers",
                "max_depth", "n_airlayers", "cell_size_horizontal",
                "cell_size_vertical_top", "depth_scale", "n_padding_x",
                "max_iterations", "target_misfit", "roughness_type",
                "stepsize_cut_count", "initial_rho", "lagrange_start")
        kw = {k: values[k] for k in keys if k in values}
        kw["diagonal_penalties"] = int(bool(values.get("diagonal_penalties")))
        if values.get("freq_min"):
            kw["freq_min"] = float(values["freq_min"])
        if values.get("freq_max"):
            kw["freq_max"] = float(values["freq_max"])
        return OccamConfig(modes=modes, **kw)

    def build(self, sites, workdir, values):
        from pycsamt.models.occam2d import InputBuilder
        from pycsamt.models.occam2d.data import OccamData
        from pycsamt.models.occam2d.mesh import OccamMesh

        wd = Path(workdir)
        cfg = self._config(values)
        InputBuilder(source=sites, workdir=wd, config=cfg).build()
        mesh = OccamMesh.read(wd / cfg.mesh_file)
        info = BuildInfo(self.key, wd, values=dict(values),
                         stations=site_names(sites))
        info.files = {k: wd / getattr(cfg, f"{k}_file")
                      for k in ("data", "mesh", "model", "startup")}
        st_x = None
        try:
            data = OccamData.read(wd / cfg.data_file)
            st_x = np.asarray(data.offsets, float)
        except Exception:
            pass
        n_air = int(getattr(mesh, "n_airlayers", 0) or 0)
        info.mesh = {"x_nodes": np.asarray(mesh.x_nodes, float),
                     "z_nodes": np.asarray(mesh.z_nodes, float),
                     "stations": st_x, "n_air": n_air}
        info.summary = [
            ("Stations", str(len(info.stations))),
            ("Mesh", f"{mesh.n_xcells} × {mesh.n_zcells} cells"),
            ("Depth", f"{mesh.z_nodes[-1]:,.0f} m"),
            ("Folder", str(wd)),
        ]
        return info

    def plot_mesh(self, fig, build):
        _constrained(fig)
        m = build.mesh
        x = m["x_nodes"] - m["x_nodes"][0]
        z = m["z_nodes"]
        st = _align(m.get("stations"), x)
        ax1 = fig.add_subplot(211)
        _grid_axes(ax1, x, z, stations=st,
                   xlim=_pad_window(st, x) if st is not None and np.size(st)
                   else None,
                   zmax=_core_depth(z, st), title="Core region (stations ▼)")
        ax2 = fig.add_subplot(212)
        _grid_axes(ax2, x, z, stations=st, title="Full mesh with padding")
        fig.suptitle(f"Occam2D mesh — {len(x) - 1} × {len(z) - 1} cells",
                     fontsize=11)

    def make_task(self, build, values, binary):
        wd = build.workdir
        max_it = int(values.get("max_iterations", 100))
        target = float(values.get("target_misfit", 1.0))

        def task(rep: RunReporter):
            from pycsamt.models.occam2d import InversionResult, OccamRunner

            runner = OccamRunner(workdir=wd, binary_path=binary)
            rep.stage("Running Occam2D")
            code = runner.run(max_iter=max_it, target_misfit=target,
                              auto_compile=False, on_output=rep.log,
                              cancel=rep.cancelled)
            if code not in (0, None):
                raise RuntimeError(f"Occam2D exited with code {code}.")
            return InversionResult(workdir=str(wd))

        return task

    def history(self, workdir):
        from pycsamt.models.occam2d.log import OccamLog

        logs = sorted(Path(workdir).glob("*.logfile"))
        if not logs:
            return []
        try:
            log = OccamLog.read(logs[0])
        except Exception:
            return []
        return [(int(i), float(r)) for i, r in zip(log.iterations, log.rms)
                if np.isfinite(r)]

    def detect(self, path):
        return _detect_backend(path) == "occam2d"

    def load(self, path, iteration=None, station=None):
        from pycsamt.models.occam2d import InversionResult

        res = InversionResult(workdir=str(path), iteration=iteration)
        run = LoadedRun(self.key, Path(path), result=res)
        run.rms = self.history(path)
        run.iterations = sorted({int(''.join(c for c in p.stem if
                                             c.isdigit()) or 0)
                                 for p in res.iter_files})
        run.iteration = iteration or (run.iterations[-1]
                                      if run.iterations else None)
        run.summary = [("Iterations", str(res.n_iterations)),
                       ("Final RMS", _fmt(res.final_rms)),
                       ("Model", f"{np.shape(res.rho_2d)}")]
        return run

    def views(self, run):
        return [("model", "Model"), ("section", "Section"),
                ("response", "Data fit"), ("pseudo", "Pseudo-section"),
                ("convergence", "Convergence"), ("mesh", "Mesh")]

    def render(self, key, fig, run):
        res = run.result
        if key == "section":
            raise PlotUnavailable("section")  # window's SectionPanel
        if key == "convergence":
            return _rms_view(fig, run, 1.0)
        if key == "mesh":
            _need(res.mesh is not None, "No mesh file",
                  "The run folder has no Occam2D mesh file.")
            ax = fig.add_subplot(111)
            x = np.asarray(res.mesh.x_nodes, float)
            _grid_axes(ax, x - x[0], res.mesh.z_nodes,
                       title="Occam2D mesh")
            return fig
        if key == "pseudo":
            # Scalar (CSAMT) data has TE only; PlotPseudo defaults to TM.
            errors = []
            for mode in ("TM", "TE"):
                fig.clear()
                try:
                    return _figure_of(res.plot_pseudo(
                        ax=fig.add_subplot(111), mode=mode), fig)
                except RuntimeError as exc:
                    errors.append(str(exc))
            raise PlotUnavailable("No pseudo-section", "; ".join(errors))
        ax = fig.add_subplot(111)
        fn = {"model": res.plot_model, "response": res.plot_response}[key]
        return _figure_of(fn(ax=ax), fig)


# ══════════════════════════════════════════════════════════════════════════════
# ModEM (2-D and 3-D)
# ══════════════════════════════════════════════════════════════════════════════


_NAN_RMS = re.compile(r"rms=\s*nan", re.I)


class ModEMEngine(Engine):
    binary_key = "modem2d"

    def __init__(self, mode: str) -> None:
        self.mode = mode
        self.key = f"modem{mode}"
        self.label = f"ModEM {mode.upper()}"
        self.dim = mode.upper()
        self.binary_key = f"modem{mode}"
        self.description = (
            "Non-linear conjugate-gradient inversion (Egbert & Kelbert, "
            f"2012) — {'Mod2DMT' if mode == '2d' else 'Mod3DMT'}.")

    def fields(self):
        common_data = [
            Field("error_floor_z", "Error floor |Z|", "float", 0.05, 0.001,
                  1.0, 0.01, 3, section="Data",
                  help="Relative impedance error floor."),
        ]
        run = [
            Field("initial_rho", "Starting ρ", "float", 100.0, 0.1, 1e6,
                  10.0, 1, unit="Ω·m"),
            Field("max_iterations", "Max iterations", "int", 100, 1, 1000),
            Field("target_rms", "Target RMS", "float", 1.05, 0.1, 50.0,
                  0.05, 2),
            Field("initial_lambda", "Initial λ", "float", 10.0, 1e-4, 1e6,
                  1.0, 3),
            Field("lambda_divisor", "λ divisor", "float", 100.0, 1.1, 1e4,
                  1.0, 1, advanced=True),
            Field("rms_diff_tol", "RMS tolerance", "float", 5e-4, 1e-7,
                  0.1, 1e-4, 5, advanced=True),
            Field("use_mpi", "Use MPI", "bool", False, section="Run"),
            Field("n_procs", "MPI processes", "int", 4, 1, 256,
                  section="Run"),
        ]
        if self.mode == "2d":
            # Mod2DMT only knows TE_/TM_Impedance blocks ("Unknown data
            # type: Full_Impedance" otherwise -- the library default, meant
            # for 3-D).  pycsamt writes one mode per data file.
            common_data.append(Field(
                "component_type", "Mode", "choice", "TE_Impedance",
                choices=(("TE_Impedance", "TE (Zxy)"),
                         ("TM_Impedance", "TM (Zyx)")), section="Data",
                help="Impedance mode inverted by Mod2DMT."))
            mesh = [
                Field("nx_2d", "Core cells (x)", "int", 100, 5, 2000,
                      section="Mesh"),
                Field("nz_2d", "Layers (z)", "int", 50, 5, 500,
                      section="Mesh"),
                Field("cell_size_h_2d", "Cell width", "float", 100.0, 1.0,
                      1e5, 10.0, 0, section="Mesh", unit="m"),
                Field("cell_size_v_top_2d", "First layer", "float", 10.0,
                      0.1, 1e4, 1.0, 1, section="Mesh", unit="m"),
                Field("depth_scale_2d", "Depth growth", "float", 1.2, 1.01,
                      2.0, 0.01, 2, section="Mesh"),
                Field("n_airlayers_2d", "Air layers", "int", 5, 1, 30,
                      section="Mesh", advanced=True),
                Field("n_padding_x_2d", "Padding cells", "int", 7, 1, 30,
                      section="Mesh", advanced=True),
            ]
            smooth = []
        else:
            common_data.append(Field(
                "component_type", "Components", "choice", "Full_Impedance",
                choices=(("Full_Impedance", "Full impedance"),
                         ("Off_Diagonal_Impedance", "Off-diagonal Z")),
                section="Data"))
            mesh = [
                Field("nx", "Core cells (x)", "int", 20, 4, 400,
                      section="Mesh"),
                Field("ny", "Core cells (y)", "int", 20, 4, 400,
                      section="Mesh"),
                Field("nz", "Layers (z)", "int", 30, 5, 300,
                      section="Mesh"),
                Field("cell_size_h", "Cell width", "float", 500.0, 1.0,
                      1e5, 50.0, 0, section="Mesh", unit="m"),
                Field("cell_size_v_top", "First layer", "float", 10.0, 0.1,
                      1e4, 1.0, 1, section="Mesh", unit="m"),
                Field("depth_scale", "Depth growth", "float", 1.2, 1.01,
                      2.0, 0.01, 2, section="Mesh"),
                Field("n_airlayers", "Air layers", "int", 5, 1, 30,
                      section="Mesh", advanced=True),
                Field("n_padding_xy", "Padding cells", "int", 7, 1, 30,
                      section="Mesh", advanced=True),
            ]
            smooth = [
                Field("smooth_x", "Smoothing x", "float", 0.1, 0.0, 1.0,
                      0.05, 2),
                Field("smooth_y", "Smoothing y", "float", 0.1, 0.0, 1.0,
                      0.05, 2),
                Field("smooth_z", "Smoothing z", "float", 0.1, 0.0, 1.0,
                      0.05, 2),
                Field("n_smooth_iter", "Smoothing passes", "int", 2, 1, 10,
                      advanced=True),
            ]
        return common_data + mesh + run[:4] + smooth + run[4:]

    def check_sites(self, sites):
        n = len(site_list(sites))
        if self.mode == "2d" and n < 2:
            return "A 2-D inversion needs at least two stations on a line."
        if self.mode == "3d" and n < 3:
            return "A 3-D inversion needs at least three stations."
        return None

    def _config(self, values, binary=None):
        from dataclasses import fields as dc_fields

        from pycsamt.models.modem import ModEmConfig

        known = {f.name for f in dc_fields(ModEmConfig)}
        kw = {k: v for k, v in values.items() if k in known}
        kw["mode"] = self.mode
        if values.get("freq_min"):
            kw["freq_min"] = float(values["freq_min"])
        if values.get("freq_max"):
            kw["freq_max"] = float(values["freq_max"])
        cfg = ModEmConfig(**kw)
        if binary:
            if self.mode == "3d":
                cfg.binary_3d = binary
            else:
                cfg.binary_2d = binary
        return cfg

    def build(self, sites, workdir, values):
        from pycsamt.models.modem import (
            InputBuilder,
            ModEmData,
            ModEmModel2D,
            ModEmModel3D,
        )

        wd = Path(workdir)
        cfg = self._config(values)
        files = InputBuilder(config=cfg).build(sites, wd)
        info = BuildInfo(self.key, wd, files=dict(files),
                         values=dict(values), stations=site_names(sites))
        reader = ModEmModel3D if self.mode == "3d" else ModEmModel2D
        model = reader.read(files["model"])
        coords = {}
        try:
            coords = ModEmData.read(files["data"]).site_coords
        except Exception:
            pass
        xyz = np.array([v for v in coords.values()], float) \
            if coords else np.empty((0, 3))
        info.mesh = {"x": np.asarray(model.x_nodes, float),
                     "z": np.asarray(model.z_nodes, float),
                     "stations": xyz}
        if self.mode == "3d":
            info.mesh["y"] = np.asarray(model.y_nodes, float)
            shape = f"{model.nx} × {model.ny} × {model.nz}"
        else:
            shape = f"{model.nx} × {model.nz}"
        info.summary = [("Stations", str(len(info.stations))),
                        ("Grid", f"{shape} cells"),
                        ("Depth", f"{model.z_nodes[-1]:,.0f} m"),
                        ("Folder", str(wd))]
        return info

    def plot_mesh(self, fig, build):
        _constrained(fig)
        m = build.mesh
        st = m.get("stations")
        if self.mode == "2d":
            x, z = m["x"], m["z"]
            sx = st[:, 1] if st is not None and len(st) else None
            if sx is not None and np.ptp(sx) == 0:
                sx = st[:, 0]
            sx = _align(sx, x)
            ax1 = fig.add_subplot(211)
            _grid_axes(ax1, x, z, stations=sx,
                       xlim=_pad_window(sx, x) if sx is not None else None,
                       zmax=_core_depth(z, sx),
                       title="Core region (stations ▼)")
            ax2 = fig.add_subplot(212)
            _grid_axes(ax2, x, z, stations=sx, title="Full grid")
            fig.suptitle(f"ModEM 2-D grid — {len(x) - 1} × {len(z) - 1}",
                         fontsize=11)
            return
        x, y, z = m["x"], m["y"], m["z"]
        if len(st):
            st = np.column_stack([_align(st[:, 0], x), _align(st[:, 1], y),
                                  st[:, 2]])
        ax1 = fig.add_subplot(121)
        win_x = _pad_window(st[:, 0], x, 0.6) if len(st) else None
        win_y = _pad_window(st[:, 1], y, 0.6) if len(st) else None
        lo_y, hi_y = win_y or (y.min(), y.max())
        lo_x, hi_x = win_x or (x.min(), x.max())
        ax1.vlines(y[(y >= lo_y) & (y <= hi_y)], lo_x, hi_x,
                   color="#4c6a92", lw=0.35)
        ax1.hlines(x[(x >= lo_x) & (x <= hi_x)], lo_y, hi_y,
                   color="#4c6a92", lw=0.35)
        if len(st):
            ax1.plot(st[:, 1], st[:, 0], "v", color="#c92a2a", ms=6,
                     mec="k", mew=0.5)
        ax1.set_xlim(lo_y, hi_y)
        ax1.set_ylim(lo_x, hi_x)
        ax1.set_aspect("equal")
        ax1.set_xlabel("East (m)")
        ax1.set_ylabel("North (m)")
        ax1.set_title("Plan view (core)", fontsize=10)
        ax2 = fig.add_subplot(122)
        _grid_axes(ax2, y, z, stations=st[:, 1] if len(st) else None,
                   xlim=(lo_y, hi_y),
                   zmax=_core_depth(z, st[:, 1] if len(st) else None),
                   title="E–W section")
        fig.suptitle(f"ModEM 3-D grid — {len(x) - 1} × {len(y) - 1} × "
                     f"{len(z) - 1}", fontsize=11)

    def make_task(self, build, values, binary):
        cfg = self._config(values, binary)
        files = build.files
        mode = self.mode

        def task(rep: RunReporter):
            from pycsamt.models.modem import ModEmRunner

            from pycsamt.models._process import ProcessCancelled

            rep.stage(f"Running {'Mod3DMT' if mode == '3d' else 'Mod2DMT'}")
            nan = {"run": 0}

            def on_line(line):
                rep.log(line)
                # ModEM never stops by itself once its misfit is NaN (the
                # line search retries forever): watch for it.
                if _NAN_RMS.search(line):
                    nan["run"] += 1
                elif "rms=" in line:
                    nan["run"] = 0

            def stop():
                return rep.cancelled() or nan["run"] >= 30

            try:
                res = ModEmRunner(build.workdir, config=cfg).run(
                    files["model"], files["data"], files["control"],
                    files.get("covariance"),
                    fwd_control=files.get("fwd_control"), mode=mode,
                    on_output=on_line, cancel=stop)
            except ProcessCancelled:
                if rep.cancelled():
                    raise
                raise RuntimeError(
                    "ModEM stopped: the misfit is NaN (the forward solution "
                    "failed). Check the data errors and the mesh — e.g. "
                    "co-located stations or missing components.") from None
            if res is None or not getattr(res, "n_iter", 0):
                # A solver that stops on bad input often still exits 0.
                tail = " / ".join(x.strip() for x in rep.lines[-3:]
                                  if x.strip())
                raise RuntimeError(
                    "ModEM finished without completing an iteration."
                    + (f" Last output: {tail}" if tail else ""))
            return res

        return task

    def history(self, workdir):
        from pycsamt.models.modem import ModEmLog

        logs = sorted(Path(workdir).glob("*.log"))
        for lp in logs:
            try:
                log = ModEmLog.read(lp)
            except Exception:
                continue
            its = getattr(log, "iterations", None)
            rms = np.asarray(getattr(log, "rms", []), float)
            if rms.size:
                its = (np.arange(rms.size) if its is None or
                       len(its) != rms.size else its)
                return [(int(i), float(r)) for i, r in zip(its, rms)
                        if np.isfinite(r)]
        return []

    def detect(self, path):
        if _detect_backend(path) != "modem":
            return False
        mode = _modem_mode(path)
        return mode in (None, self.mode)

    def load(self, path, iteration=None, station=None):
        from pycsamt.models.modem import InversionResult

        res = InversionResult(path)
        run = LoadedRun(self.key, Path(path), result=res)
        run.rms = [(int(i), float(r)) for i, r in
                   zip(res.iteration_numbers, res.rms_history)]
        run.iterations = [int(i) for i in res.iteration_numbers]
        run.iteration = run.iterations[-1] if run.iterations else None
        run.summary = [("Mode", str(res.mode).upper()),
                       ("Iterations", str(res.n_iter)),
                       ("Final RMS", _fmt(res.final_rms))]
        return run

    def views(self, run):
        if self.mode == "2d":  # one profile: no plan-view misfit map
            return [("model", "Model"), ("response", "Data fit"),
                    ("pseudo", "Pseudo-section"),
                    ("convergence", "Convergence")]
        return [("model", "Model"), ("section", "Section"),
                ("response", "Data fit"), ("pseudo", "Pseudo-section"),
                ("convergence", "Convergence"), ("misfitmap", "Misfit map")]

    def render(self, key, fig, run):
        from pycsamt.models.modem import plot as mp

        res = run.result
        if key == "convergence":
            return _rms_view(fig, run, getattr(res.config, "target_rms",
                                               None))
        if key == "model":
            _need(res.model_final is not None, "No final model",
                  "No ModEM model file was found in the run folder.")
            return (mp.PlotDepthMap(result=res).plot() if self.mode == "3d"
                    else mp.PlotModel2D(result=res).plot())
        if key == "section":
            return mp.PlotSection(result=res).plot()
        _need(res.data_obs is not None, "No observed data",
              "The run folder has no ModEM data file.")
        if key in ("response", "misfitmap"):
            _need(res.data_pred is not None, "No predicted data",
                  "ModEM writes the predicted data (<stem>_NLCG_NNN.dat) "
                  "after each iteration; none is in this folder yet.")
        if key == "response":
            if self.mode == "2d":
                # PlotDataFit/PlotResponse expect impedance-tensor
                # components; Mod2DMT data are TE/TM blocks.
                return _modem_2d_fit(fig, res)
            return mp.PlotDataFit(result=res, max_stations=4).plot()
        if key == "misfitmap":
            return mp.PlotMisfitMap(result=res).plot()
        if key == "pseudo":
            return mp.PlotPseudo(result=res).plot()
        raise PlotUnavailable("Unknown view", key)


def _modem_rows(data):
    """{site: (period, Z, err)} arrays from ModEM data blocks."""
    out: dict[str, list] = {}
    if data is None:
        return {}
    names = list(getattr(data, "site_names", []))
    for blk in getattr(data, "blocks", []) or []:
        for row in blk.get("rows", []):
            period, si, _x, _y, _z, _c, re_, im, err = row[:9]
            name = names[si] if isinstance(si, int) and si < len(names)                 else str(si)
            out.setdefault(name, []).append((period, complex(re_, im), err))
    return {k: np.array(v, dtype=object) for k, v in out.items()}


def _modem_2d_fit(fig, res):
    """Apparent resistivity and phase, observed vs predicted (2-D ModEM)."""
    obs = _modem_rows(res.data_obs)
    pred = _modem_rows(getattr(res, "data_pred", None))
    _need(obs, "No observed data", "The ModEM data file has no rows.")
    names = list(obs)[:4]
    _constrained(fig)
    for j, name in enumerate(names):
        o = obs[name]
        T = o[:, 0].astype(float)
        Z = o[:, 1].astype(complex)
        E = o[:, 2].astype(float)
        rho = 0.2 * T * np.abs(Z) ** 2
        phi = np.mod(np.degrees(np.angle(Z)), 180.0)
        ax1 = fig.add_subplot(2, len(names), j + 1)
        ax2 = fig.add_subplot(2, len(names), len(names) + j + 1, sharex=ax1)
        rerr = 2.0 * rho * E / np.maximum(np.abs(Z), 1e-30)
        ax1.errorbar(T, rho, yerr=rerr, fmt="o", ms=3, color="#c92a2a",
                     label="observed")
        ax2.errorbar(T, phi, yerr=np.degrees(E / np.maximum(np.abs(Z),
                                                            1e-30)),
                     fmt="o", ms=3, color="#c92a2a")
        if name in pred:
            p = pred[name]
            order = np.argsort(p[:, 0].astype(float))
            Tp = p[order, 0].astype(float)
            Zp = p[order, 1].astype(complex)
            ax1.plot(Tp, 0.2 * Tp * np.abs(Zp) ** 2, "-", color="#1864ab",
                     label="model")
            ax2.plot(Tp, np.mod(np.degrees(np.angle(Zp)), 180.0), "-",
                     color="#1864ab")
        ax1.set_xscale("log")
        ax1.set_yscale("log")
        ax1.set_title(name, fontsize=9)
        ax2.set_xlabel("Period (s)")
        if j == 0:
            ax1.set_ylabel("ρa (Ω·m)")
            ax2.set_ylabel("φ (°)")
            ax1.legend(fontsize=7)
        for ax in (ax1, ax2):
            ax.grid(True, which="both", alpha=0.25)
    return fig


def _modem_mode(path) -> str | None:
    """'2d' / '3d' from the model file header in a ModEM folder."""
    for p in sorted(Path(path).glob("*.rho")) + sorted(Path(path).glob(
            "*.ws")):
        try:
            head = p.read_text(errors="replace").splitlines()[:3]
        except Exception:
            continue
        for line in head[1:]:
            parts = line.split()
            if len(parts) >= 3 and parts[0].isdigit() and parts[1].isdigit():
                return "3d" if (len(parts) >= 4 and parts[2].isdigit()) \
                    else "2d"
    return None


# ══════════════════════════════════════════════════════════════════════════════
# MARE2DEM
# ══════════════════════════════════════════════════════════════════════════════


class Mare2DEMEngine(Engine):
    key = "mare2dem"
    label = "MARE2DEM"
    dim = "2D"
    binary_key = "mare2dem"
    description = ("Adaptive finite-element 2-D inversion (Key, 2016) — "
                   "Linux/WSL, MPI.")

    def fields(self):
        return [
            Field("output_modes", "Modes", "choice", "all",
                  choices=(("all", "TE + TM"), ("te", "TE only"),
                           ("tm", "TM only")), section="Data"),
            Field("error_floor", "Error floor |Z|", "float", 0.05, 0.001,
                  1.0, 0.01, 3, section="Data"),
            Field("cell_y", "Cell width", "float", 0.0, 0.0, 1e5, 10.0, 0,
                  section="Mesh", unit="m", auto_zero=True,
                  help="0 = half the median station spacing."),
            Field("z_first", "First layer", "float", 0.0, 0.0, 1e4, 5.0, 1,
                  section="Mesh", unit="m", auto_zero=True,
                  help="0 = a quarter of the shallowest skin depth."),
            Field("depth_max", "Max depth", "float", 0.0, 0.0, 1e7, 500.0,
                  0, section="Mesh", unit="m", auto_zero=True,
                  help="0 = 1.5 × the deepest skin depth."),
            Field("growth", "Depth growth", "float", 1.15, 1.0, 2.0, 0.01,
                  2, section="Mesh"),
            Field("initial_rho", "Starting ρ", "float", 100.0, 0.1, 1e6,
                  10.0, 1, unit="Ω·m"),
            Field("max_iterations", "Max iterations", "int", 30, 1, 500),
            Field("target_rms", "Target RMS", "float", 1.0, 0.1, 50.0, 0.1,
                  2),
            Field("use_mpi", "Use MPI", "bool", True, section="Run"),
            Field("n_procs", "MPI processes", "int", 4, 2, 256,
                  section="Run",
                  help="MARE2DEM uses one manager and n−1 workers."),
        ]

    def check_sites(self, sites):
        n = len(site_list(sites))
        if n < 2:
            return "A 2-D inversion needs at least two stations on a line."
        return None

    def _config(self, values, binary=None):
        from pycsamt.models.mare2dem import Mare2DEMConfig

        return Mare2DEMConfig(
            max_iterations=int(values.get("max_iterations", 30)),
            target_rms=float(values.get("target_rms", 1.0)),
            initial_rho=float(values.get("initial_rho", 100.0)),
            use_mpi=bool(values.get("use_mpi", True)),
            n_procs=int(values.get("n_procs", 4)),
            binary=binary or "MARE2DEM",
        )

    def build(self, sites, workdir, values):
        from pycsamt.models.mare2dem import build_mt_inputs

        def opt(k):
            v = float(values.get(k) or 0)
            return v if v > 0 else None

        wd = Path(workdir)
        out = build_mt_inputs(
            sites, wd, self._config(values),
            error_floor=float(values.get("error_floor", 0.05)),
            output_modes=str(values.get("output_modes", "all")),
            cell_y=opt("cell_y"), z_first=opt("z_first"),
            depth_max=opt("depth_max"),
            growth=float(values.get("growth", 1.15)),
        )
        info = BuildInfo(self.key, wd, files=dict(out["files"]),
                         values=dict(values), stations=site_names(sites))
        info.mesh = {k: out[k] for k in ("y_edges", "z_edges", "rx_y",
                                         "rx_z", "padding")}
        info.mesh["run_stem"] = out["run_stem"]
        ye, ze = out["y_edges"], out["z_edges"]
        info.summary = [
            ("Receivers", str(len(out["rx_y"]))),
            ("Core grid", f"{len(ye) - 1} × {len(ze) - 1} "
                          f"({out['n_parameters']} free cells)"),
            ("Depth", f"{ze[-1]:,.0f} m"),
            ("Padding", f"{out['padding'] / 1000:,.0f} km"),
            ("Folder", str(wd)),
        ]
        rx = np.round(np.asarray(out["rx_y"], float), 3)
        if np.unique(rx).size < rx.size:
            info.warnings.append(
                "Some stations share the same position along the profile; "
                "MARE2DEM can fail on coincident receivers.")
        return info

    def plot_mesh(self, fig, build):
        _constrained(fig)
        m = build.mesh
        ye, ze = m["y_edges"], m["z_edges"]
        ax1 = fig.add_subplot(211)
        _grid_axes(ax1, ye, ze, stations=m["rx_y"],
                   zmax=min(ze[-1], max(ze[-1] * 0.3, ze[min(8, len(ze) - 1)])),
                   title="Core grid, near surface (receivers ▼)")
        ax2 = fig.add_subplot(212)
        _grid_axes(ax2, ye, ze, stations=m["rx_y"], title="Whole core grid")
        fig.suptitle(
            f"MARE2DEM starting mesh — {len(ye) - 1} × {len(ze) - 1} "
            "regions (refined adaptively by the solver)", fontsize=11)

    def make_task(self, build, values, binary):
        cfg = self._config(values, binary)
        stem = build.mesh["run_stem"]

        def task(rep: RunReporter):
            from pycsamt.models.mare2dem import Mare2DEMRunner

            rep.stage("Running MARE2DEM")
            res = Mare2DEMRunner(build.workdir, config=cfg).run(
                stem, on_output=rep.log, cancel=rep.cancelled)
            if not res or not res.n_iterations:
                tail = " / ".join(rep.lines[-3:]) if rep.lines else ""
                raise RuntimeError(
                    "MARE2DEM finished without completing an iteration."
                    + (f" Last output: {tail}" if tail else ""))
            return res

        return task

    def history(self, workdir):
        from pycsamt.models.mare2dem.log import Mare2DEMLog
        from pycsamt.models.mare2dem.validation import is_log_file

        for p in sorted(Path(workdir).iterdir()):
            try:
                if p.is_file() and is_log_file(p):
                    return _m2d_rms(Mare2DEMLog(p))
            except Exception:
                continue
        return []

    def detect(self, path):
        return _detect_backend(path) == "mare2dem"

    def load(self, path, iteration=None, station=None):
        from pycsamt.models.mare2dem import InversionResult

        res = InversionResult(path)
        run = LoadedRun(self.key, Path(path), result=res)
        run.rms = _m2d_rms(res.log) if res.log is not None else []
        run.iterations = [i for i, _ in run.rms]
        run.iteration = run.iterations[-1] if run.iterations else None
        run.summary = [("Iterations", str(res.n_iterations)),
                       ("Final RMS", _fmt(res.final_rms)),
                       ("Converged", "yes" if res.converged else "no")]
        return run

    def views(self, run):
        return [("model", "Model"), ("response", "Data fit"),
                ("convergence", "Convergence")]

    def render(self, key, fig, run):
        from pycsamt.models.mare2dem import plot as m2p

        res = run.result
        if key == "convergence":
            return _rms_view(fig, run, 1.0)
        if key == "model":
            _need(res.model is not None, "No model",
                  "No .resistivity file was found in the run folder.")
            ax = fig.add_subplot(111)
            return _figure_of(m2p.PlotModel(res).plot(ax=ax), fig)
        if key == "response":
            _need(res.response is not None, "No response",
                  "MARE2DEM writes the model response (.resp) after each "
                  "iteration; none is in this folder yet.")
            ax = fig.add_subplot(111)
            return _figure_of(m2p.PlotResponse(res).plot(ax=ax), fig)
        raise PlotUnavailable("Unknown view", key)


def _m2d_rms(log) -> list[tuple[int, float]]:
    """(iteration, RMS) pairs of a Mare2DEMLog (``rms_history()`` is a
    method there, not a property)."""
    hist = log.rms_history
    hist = hist() if callable(hist) else hist
    out = []
    for i, r in enumerate(hist or [], 1):
        try:
            v = float(r)
        except (TypeError, ValueError):
            continue
        if math.isfinite(v):
            out.append((i, v))
    return out


# ── Registry & detection ──────────────────────────────────────────────────────


def _fmt(v) -> str:
    try:
        v = float(v)
    except (TypeError, ValueError):
        return "–"
    return f"{v:.3f}" if math.isfinite(v) else "–"


def _detect_backend(path) -> str | None:
    try:
        from pycsamt.format import detect_source

        kind = detect_source(Path(path))
    except Exception:
        return None
    return kind.backend if kind.category == "solver" else None


ENGINES: dict[str, Engine] = {
    e.key: e for e in (Occam1DEngine(), Occam2DEngine(), ModEMEngine("2d"),
                       ModEMEngine("3d"), Mare2DEMEngine())
}


def engine(key: str) -> Engine:
    return ENGINES[key]


def detect_engine(path) -> str | None:
    """Engine key for a run folder, or ``None`` if not recognised."""
    p = Path(path)
    if not p.is_dir():
        return None
    for key in ("occam1d", "occam2d", "mare2dem", "modem2d", "modem3d"):
        try:
            if ENGINES[key].detect(p):
                return key
        except Exception:
            continue
    return None


def sites_to_1d_features(sites, freqs) -> tuple[np.ndarray | None,
                                                   list[str]]:
    """Real-data features for the AI 1-D inverter.

    Returns ``X`` of shape ``(n_sites, 2 * n_freq)`` laid out like the
    training set of :func:`pycsamt.forward.batch.generate_dataset`
    (``[log10 rho_a ..., phase ...]``), plus the station names.  Each
    sounding is the Berdichevsky average of the off-diagonal responses
    (rho = sqrt(rho_xy * rho_yx), mean of the first-quadrant phases), or XY
    alone for scalar data; it is interpolated in log-frequency onto *freqs*
    (values outside the measured band take the nearest measured value).

    The desktop used to call ``site.interpolate_rho_a``, which does not
    exist, inside a bare ``except`` -- every AI 1-D "prediction" silently
    ran on five synthetic training samples instead of the loaded stations.
    """
    f_out = np.log10(np.asarray(freqs, float))
    rows, names = [], []
    for s in site_list(sites):
        try:
            f = np.asarray(s.freq, float)
            rho = np.asarray(s.rho, float)
            ph = np.asarray(s.phase, float)
        except Exception:
            continue
        if f.ndim != 1 or rho.ndim != 3 or f.size < 2:
            continue
        rxy, ryx = rho[:, 0, 1], rho[:, 1, 0]
        pxy = np.mod(ph[:, 0, 1], 180.0)
        pyx = np.mod(ph[:, 1, 0], 180.0)
        both = np.isfinite(ryx) & (ryx > 0) & np.isfinite(pyx)
        r = np.where(both, np.sqrt(np.abs(rxy * ryx)), rxy)
        p = np.where(both, 0.5 * (pxy + pyx), pxy)
        ok = np.isfinite(r) & (r > 0) & np.isfinite(p) & (f > 0)
        if ok.sum() < 2:
            continue
        lf = np.log10(f[ok])
        order = np.argsort(lf)
        lr = np.interp(f_out, lf[order], np.log10(r[ok][order]))
        pp = np.interp(f_out, lf[order], p[ok][order])
        rows.append(np.concatenate([lr, pp]))
        names.append(str(getattr(s, "name", len(names))))
    if not rows:
        return None, []
    return np.asarray(rows, float), names


def elapsed_eta(history: list, started: float, max_iter: int,
                now: float | None = None) -> tuple[float, float | None]:
    """Elapsed seconds and an ETA from the mean iteration time."""
    now = time.monotonic() if now is None else now
    elapsed = now - started
    n = len(history)
    if n < 2 or max_iter <= 0:
        return elapsed, None
    per = elapsed / n
    return elapsed, max(per * (max_iter - n), 0.0)
