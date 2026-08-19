# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""Validated post-inversion access for native Occam1D runs.

The result layer discovers native files by their contents, selects one
iteration deterministically, verifies scientific relationships, and exposes
stable text and image products. Parsers own syntax; this module owns
run-level consistency and provenance.
"""

from __future__ import annotations

import csv
import json
import math
import re
from dataclasses import asdict, dataclass
from enum import Enum
from pathlib import Path
from typing import Any

import numpy as np

from ...compat.sklearn import validate_params
from .base import Occam1DBase
from .data import Occam1DData
from .log import Occam1DLog
from .model import Occam1DModel
from .response import Occam1DResponse
from .schema import ITERATION_SELECTOR_SCHEMA
from .startup import Occam1DStartup
from .validation import Occam1DFileType, detect_file_type

__all__ = [
    "Occam1DResult",
    "Occam1DResultFiles",
    "Occam1DResultState",
]

_ITERATION_RE = re.compile(r"(?:iter|resp)[^0-9]*([0-9]+)", re.I)
_MANIFEST_NAME = "occam1d_manifest.json"


class Occam1DResultState(str, Enum):
    """Lifecycle of a run-level result aggregate."""

    DISCOVERING = "discovering"
    LOADING = "loading"
    READY = "ready"
    FAILED = "failed"


@dataclass(frozen=True)
class Occam1DResultFiles:
    """Immutable paths participating in one selected result.

    Parameters
    ----------
    workdir : pathlib.Path
        Absolute run directory.
    data, model : pathlib.Path
        Required native sounding and layer-model files.
    startup : pathlib.Path or None
        Initial iteration-zero controls, when present.
    iteration : pathlib.Path or None
        Selected solver iteration. It is ``None`` for startup-only runs.
    response, log, manifest : pathlib.Path or None
        Optional modeled response, convergence log, and build manifest.
    """

    workdir: Path
    data: Path
    model: Path
    startup: Path | None
    iteration: Path | None
    response: Path | None
    log: Path | None
    manifest: Path | None

    def to_dict(self) -> dict[str, str | None]:
        """Return JSON-friendly absolute paths."""
        return {
            name: str(value) if value is not None else None
            for name, value in asdict(self).items()
        }


def _filename_number(path: Path) -> int | None:
    """Return an iteration suffix when a filename clearly carries one."""
    match = _ITERATION_RE.search(path.stem)
    return int(match.group(1)) if match else None


def _safe_prefix(value: Any) -> str:
    """Return a filesystem-safe, non-empty export prefix."""
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value).strip())
    text = text.strip("._")
    if not text:
        raise ValueError("prefix must contain a filesystem-safe character.")
    return text


class Occam1DResult(Occam1DBase):
    r"""Load and validate one completed or startup-only Occam1D run.

    Parameters
    ----------
    workdir : path-like
        Directory containing native Occam1D files. Files are identified from
        format signatures, not solely from their names.
    iteration : {"best", "latest"} or int, default="best"
        Selection policy. ``best`` uses the lowest finite log RMS whose
        iteration file exists; otherwise it falls back to the latest parsed
        iteration. An integer must exist exactly once.
    strict : bool, default=True
        Reject ambiguous required files, duplicate iteration numbers, and
        malformed recognized files. With ``False``, deterministic choices
        are made and recoverable issues are recorded in :attr:`warnings`.
    verbose, logger, metadata, stream
        Shared :class:`Occam1DBase` options.

    Attributes
    ----------
    data : Occam1DData
        Observed sounding data. Resistivity is in ohm metres and phase in
        degrees.
    model : Occam1DModel
        Layer geometry, with depth and thickness in metres.
    iteration_data : Occam1DStartup
        Selected log10-resistivity parameter vector.
    resistivity : ndarray of shape (n_layers,)
        Selected physical layer resistivities in ohm metres.
    files : Occam1DResultFiles
        Immutable provenance for every selected file.

    Raises
    ------
    FileNotFoundError
        If the directory or a required data/model/startup file is absent.
    ValueError
        If selection is ambiguous or cross-file scientific invariants fail.

    Notes
    -----
    A solver response is optional because an interrupted run may have a valid
    iteration but no matching response. Missing fit metrics remain ``nan``;
    they are never converted to physical zero.
    """

    @validate_params(ITERATION_SELECTOR_SCHEMA)
    def __init__(
        self,
        workdir,
        iteration="best",
        *,
        strict=True,
        **kwargs,
    ):
        super().__init__(**kwargs)
        if not isinstance(strict, bool):
            raise TypeError("strict must be a bool.")
        self.workdir = self._normalize_path(workdir)
        self.strict = strict
        self.result_state = Occam1DResultState.DISCOVERING
        self._requested_iteration = iteration
        self._inventory: dict[str, tuple[Path, ...]] = {}
        try:
            self._load()
        except Exception as error:
            self.result_state = Occam1DResultState.FAILED
            self.mark_invalid(str(error))
            raise

    def _load(self) -> None:
        """Discover, parse, select, and cross-validate run components."""
        if not self.workdir.exists():
            raise FileNotFoundError(
                f"Occam1D run directory not found: {self.workdir}"
            )
        if not self.workdir.is_dir():
            raise NotADirectoryError(
                f"Occam1D workdir is not a directory: {self.workdir}"
            )
        self._inventory = self._discover()
        self.result_state = Occam1DResultState.LOADING

        data_path = self._one_required(Occam1DFileType.DATA)
        model_path = self._one_required(Occam1DFileType.MODEL)
        startup_path = self._one_optional(Occam1DFileType.STARTUP)
        log_path = self._one_optional(Occam1DFileType.LOG)
        self.data = Occam1DData.read(data_path)
        self.model = Occam1DModel.read(model_path)
        self.log = Occam1DLog.read(log_path) if log_path else None

        iterations = self._parse_iterations()
        self.iter_files = [path for _, path, _ in iterations]
        self.available_iterations = tuple(item[0] for item in iterations)
        self.iteration = self._select_iteration(iterations)
        selected = next(
            (item for item in iterations if item[0] == self.iteration), None
        )
        if selected is None:
            if startup_path is None:
                raise FileNotFoundError(
                    "No native Occam1D iteration or Startup file was found."
                )
            self.iteration_file = None
            self.iteration_data = Occam1DStartup.read(startup_path)
        else:
            _, self.iteration_file, self.iteration_data = selected

        self._validate_parameter_binding()
        self.resistivity = self.iteration_data.physical_resistivity
        self.response_file = self._select_response(self.iteration)
        self.response = (
            Occam1DResponse.read(self.response_file)
            if self.response_file is not None
            else None
        )
        if self.response is not None:
            self.response.validate_against(self.data)

        manifest = self.workdir / _MANIFEST_NAME
        self.files = Occam1DResultFiles(
            workdir=self.workdir,
            data=data_path,
            model=model_path,
            startup=startup_path,
            iteration=self.iteration_file,
            response=self.response_file,
            log=log_path,
            manifest=manifest if manifest.is_file() else None,
        )
        self.path = self.workdir
        self.result_state = Occam1DResultState.READY
        self.logger.info(
            "Loaded Occam1D result %s at iteration %d.",
            self.workdir,
            self.iteration,
        )
        self._emit(
            f"Loaded Occam1D iteration {self.iteration} from "
            f"{self.workdir}"
        )

    def _discover(self) -> dict[str, tuple[Path, ...]]:
        """Classify direct child files by native content signature."""
        found: dict[str, list[Path]] = {}
        for path in sorted(self.workdir.iterdir(), key=lambda item: item.name):
            if not path.is_file() or path.name == _MANIFEST_NAME:
                continue
            try:
                kind = detect_file_type(path)
            except (OSError, UnicodeError) as error:
                if self.strict:
                    raise ValueError(
                        f"Cannot inspect Occam1D candidate {path}: {error}"
                    ) from error
                self.add_warning(f"Skipped unreadable result file: {path}")
                continue
            if kind != Occam1DFileType.UNKNOWN:
                found.setdefault(kind, []).append(path.resolve())
        return {kind: tuple(paths) for kind, paths in found.items()}

    def _one_required(self, kind: str) -> Path:
        """Return one unambiguous required path of ``kind``."""
        path = self._one_optional(kind)
        if path is None:
            raise FileNotFoundError(
                f"Run directory has no recognized Occam1D {kind} file: "
                f"{self.workdir}"
            )
        return path

    def _one_optional(self, kind: str) -> Path | None:
        """Return one optional path, enforcing ambiguity policy."""
        paths = self._inventory.get(kind, ())
        if not paths:
            return None
        if len(paths) > 1:
            names = ", ".join(path.name for path in paths)
            message = f"Multiple Occam1D {kind} files found: {names}."
            if self.strict:
                raise ValueError(message)
            self.add_warning(message + f" Using {paths[0].name}.")
        return paths[0]

    def _parse_iterations(
        self,
    ) -> list[tuple[int, Path, Occam1DStartup]]:
        """Parse iteration files and enforce unique internal identifiers."""
        parsed = []
        numbers: dict[int, Path] = {}
        for path in self._inventory.get(Occam1DFileType.ITER, ()):
            try:
                item = Occam1DStartup.read(path)
            except (OSError, TypeError, ValueError) as error:
                if self.strict:
                    raise ValueError(
                        f"Cannot parse Occam1D iteration {path}: {error}"
                    ) from error
                self.add_warning(f"Skipped invalid iteration file: {path}")
                continue
            previous = numbers.get(item.iteration)
            if previous is not None:
                message = (
                    f"Iteration {item.iteration} occurs in both "
                    f"{previous.name} and {path.name}."
                )
                if self.strict:
                    raise ValueError(message)
                self.add_warning(message + f" Using {previous.name}.")
                continue
            suffix = _filename_number(path)
            if suffix is not None and suffix != item.iteration:
                self.add_warning(
                    f"Filename {path.name} suggests iteration {suffix}, but "
                    f"its native header declares {item.iteration}."
                )
            numbers[item.iteration] = path
            parsed.append((item.iteration, path, item))
        return sorted(parsed, key=lambda value: (value[0], value[1].name))

    def _select_iteration(self, iterations) -> int:
        """Resolve the public iteration selection policy."""
        available = tuple(item[0] for item in iterations)
        requested = self._requested_iteration
        if isinstance(requested, int) and not isinstance(requested, bool):
            if requested not in available:
                shown = ", ".join(map(str, available)) or "none"
                raise ValueError(
                    f"Iteration {requested} is unavailable; available "
                    f"iterations: {shown}."
                )
            return requested
        if requested == "best" and self.log is not None:
            candidates = [
                (record.rms, record.iteration)
                for record in self.log.records
                if record.iteration in available and math.isfinite(record.rms)
            ]
            if candidates:
                return min(candidates, key=lambda value: (value[0], value[1]))[
                    1
                ]
        return max(available) if available else 0

    def _validate_parameter_binding(self) -> None:
        """Verify selected parameters and referenced filenames."""
        item = self.iteration_data
        if item.parameters.size != self.model.n_layers:
            raise ValueError(
                "Iteration parameter count does not match the layer model: "
                f"{item.parameters.size} != {self.model.n_layers}."
            )
        if not np.all(np.isfinite(item.parameters)):
            raise ValueError("Selected iteration contains non-finite values.")
        for label, reference, actual in (
            ("data", item.data_file, self.data.path),
            ("model", item.model_file, self.model.path),
        ):
            if Path(reference).name.lower() != actual.name.lower():
                self.add_warning(
                    f"Iteration references {label} file {reference!r}, but "
                    f"the recognized file is {actual.name!r}."
                )

    def _select_response(self, iteration: int) -> Path | None:
        """Select the uniquely named response for ``iteration``."""
        paths = self._inventory.get(Occam1DFileType.RESPONSE, ())
        matches = [
            path for path in paths if _filename_number(path) == iteration
        ]
        if len(matches) > 1:
            names = ", ".join(path.name for path in matches)
            message = (
                f"Multiple responses match iteration {iteration}: {names}."
            )
            if self.strict:
                raise ValueError(message)
            self.add_warning(message + f" Using {matches[0].name}.")
        if matches:
            return matches[0]
        if paths:
            self.add_warning(
                f"No response file matches selected iteration {iteration}."
            )
        return None

    @property
    def is_ready(self) -> bool:
        """Whether all required run-level invariants were validated."""
        return self.result_state is Occam1DResultState.READY

    @property
    def n_iterations(self) -> int:
        """Number of unique, parsed solver iteration files."""
        return len(self.available_iterations)

    @property
    def final_rms(self) -> float:
        """Normalized RMS for the selected iteration, or ``nan``."""
        if self.response is not None and math.isfinite(self.response.rms):
            return self.response.rms
        if self.log is not None:
            try:
                return float(self.log.get_iteration(self.iteration).rms)
            except KeyError:
                pass
        return float("nan")

    @property
    def converged(self) -> bool:
        """Whether selected RMS reaches its iteration target."""
        return bool(
            math.isfinite(self.final_rms)
            and self.final_rms <= self.iteration_data.target_misfit
        )

    def export_model(self, path) -> Path:
        """Export physical layers to CSV with explicit SI units."""
        target = self._prepare_output_file(path)
        with target.open("w", newline="", encoding="utf8") as stream:
            writer = csv.writer(stream)
            writer.writerow(
                ["layer", "top_depth_m", "thickness_m", "resistivity_ohm_m"]
            )
            writer.writerows(
                (index, depth, thickness, rho)
                for index, (depth, thickness, rho) in enumerate(
                    zip(
                        self.model.depth,
                        self.model.thickness,
                        self.resistivity,
                    ),
                    start=1,
                )
            )
        return target

    def export_response(self, path) -> Path:
        """Export physical response values and weighted residuals as CSV."""
        if self.response is None:
            raise RuntimeError(
                "No response file is available for the selected iteration."
            )
        target = self._prepare_output_file(path)
        observed, modeled = self.response.physical_values()
        errors = self.response.physical_errors()
        with target.open("w", newline="", encoding="utf8") as stream:
            writer = csv.writer(stream)
            writer.writerow(
                [
                    "frequency_index",
                    "type_code",
                    "observed_physical",
                    "modeled_physical",
                    "error_physical",
                    "weighted_residual",
                ]
            )
            writer.writerows(
                zip(
                    self.response.frequency_index,
                    self.response.type_code,
                    observed,
                    modeled,
                    errors,
                    self.response.residuals,
                )
            )
        return target

    def export_iterations(self, path) -> Path:
        """Export normalized convergence history as UTF-8 CSV."""
        if self.log is None:
            raise RuntimeError("No convergence log is available.")
        return self.log.to_csv(path)

    def export_summary(self, path) -> Path:
        """Write :meth:`summary` to a UTF-8 text file."""
        target = self._prepare_output_file(path)
        target.write_text(self.summary(), encoding="utf8")
        return target

    def export_metadata(self, path) -> Path:
        """Export JSON-safe result diagnostics and provenance."""
        target = self._prepare_output_file(path)
        target.write_text(
            json.dumps(self.diagnostics(), indent=2, sort_keys=True) + "\n",
            encoding="utf8",
        )
        return target

    def export_text(self, directory, prefix=None) -> dict[str, Path]:
        """Write the stable ``model-text`` product family."""
        output = self._normalize_path(directory)
        output.mkdir(parents=True, exist_ok=True)
        name = _safe_prefix(prefix or self.data.station)
        paths = {
            "model": self.export_model(output / f"{name}_model.csv"),
            "summary": self.export_summary(output / f"{name}_summary.txt"),
        }
        if self.response is not None:
            paths["response"] = self.export_response(
                output / f"{name}_response.csv"
            )
        if self.log is not None:
            paths["iterations"] = self.export_iterations(
                output / f"{name}_iterations.csv"
            )
        return paths

    def plot_model(self, **kwargs):
        """Return the selected physical resistivity figure."""
        from .plot import PlotModel

        depth_max = kwargs.pop("depth_max", None)
        show_station = kwargs.pop("show_station", True)
        return PlotModel(self, **kwargs).plot(
            depth_max=depth_max,
            show_station=show_station,
        )

    def plot_response(self, **kwargs):
        """Return observed and modeled sounding-response panels."""
        from .plot import PlotResponse

        return PlotResponse(self, **kwargs).plot()

    def plot_convergence(self, **kwargs):
        """Return normalized RMS and roughness convergence history."""
        from .plot import PlotConvergence

        target = kwargs.pop("target", None)
        show_roughness = kwargs.pop("show_roughness", True)
        return PlotConvergence(self, **kwargs).plot(
            target=target,
            show_roughness=show_roughness,
        )

    plot_misfit = plot_convergence

    def plot_summary(self, **kwargs):
        """Return the combined model, response, and convergence figure."""
        from .plot import PlotSummary

        depth_max = kwargs.pop("depth_max", None)
        target = kwargs.pop("target", None)
        return PlotSummary(self, **kwargs).plot(
            depth_max=depth_max,
            target=target,
        )

    def save_main_images(
        self,
        directory,
        prefix=None,
        dpi=180,
        *,
        fmt="png",
        style=None,
        style_overrides=None,
    ) -> dict[str, Path]:
        """Save four style-aware ``model-image`` products.

        Parameters ``style`` and ``style_overrides`` follow
        :mod:`pycsamt.api.occam1d`; the latter accepts nested keys such as
        ``model__color`` and ``response_legend``.
        """
        import matplotlib.pyplot as plt

        if not isinstance(dpi, int) or isinstance(dpi, bool) or dpi < 1:
            raise ValueError("dpi must be a positive integer.")
        if not isinstance(fmt, str) or not fmt.strip(". "):
            raise ValueError("fmt must be a non-empty format name.")
        output = self._normalize_path(directory)
        output.mkdir(parents=True, exist_ok=True)
        name = _safe_prefix(prefix or self.data.station)
        extension = fmt.lower().strip(". ")
        factories = {
            "model": self.plot_model,
            "response": self.plot_response,
            "convergence": self.plot_convergence,
            "summary": self.plot_summary,
        }
        paths = {}
        for kind, factory in factories.items():
            figure = factory(
                dpi=dpi,
                style=style,
                style_overrides=style_overrides,
            )
            path = output / f"{name}_occam1d_{kind}.{extension}"
            try:
                figure.savefig(
                    path,
                    dpi=dpi,
                    bbox_inches="tight",
                    format=extension,
                )
            finally:
                plt.close(figure)
            paths[kind] = path
        return paths

    def summary(self) -> str:
        """Return a human-readable selected-inversion report."""
        rms = (
            f"{self.final_rms:.6g}"
            if math.isfinite(self.final_rms)
            else "N/A"
        )
        return (
            "Occam1D inversion result\n"
            f"  workdir   : {self.workdir}\n"
            f"  station   : {self.data.station}\n"
            f"  mode      : {self.data.mode}\n"
            f"  iteration : {self.iteration}\n"
            f"  available : {list(self.available_iterations)}\n"
            f"  RMS       : {rms}\n"
            f"  converged : {self.converged}\n"
            f"  layers    : {self.model.n_layers}\n"
            f"  depth max : {self.model.depth[-1]:g} m\n"
            f"  response  : {self.response is not None}\n"
        )

    def diagnostics(self) -> dict[str, Any]:
        """Return JSON-safe science, selection, and file provenance."""
        values = super().diagnostics()
        values.update(
            {
                "result_state": self.result_state.value,
                "workdir": str(self.workdir),
                "strict": self.strict,
                "requested_iteration": self._requested_iteration,
                "selected_iteration": self.iteration,
                "available_iterations": list(self.available_iterations),
                "n_iterations": self.n_iterations,
                "station": self.data.station,
                "mode": self.data.mode,
                "n_layers": self.model.n_layers,
                "final_rms": self.final_rms,
                "target_misfit": self.iteration_data.target_misfit,
                "converged": self.converged,
                "has_response": self.response is not None,
                "has_log": self.log is not None,
                "files": self.files.to_dict(),
            }
        )
        return values
