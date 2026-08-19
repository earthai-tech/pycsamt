# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Transactional construction of single and batch Occam1D input sets."""

from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import tempfile
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Any

from ...api.view.progress import get_progress_bar
from .base import Occam1DBase
from .config import Occam1DConfig
from .data import Occam1DData
from .model import Occam1DModel
from .processing import extract_sounding, normalise_sites
from .startup import Occam1DStartup

__all__ = [
    "Occam1DBatch",
    "Occam1DBuildState",
    "Occam1DInputBuilder",
    "Occam1DInputPaths",
]

_MANIFEST_NAME = "occam1d_manifest.json"


class Occam1DBuildState(str, Enum):
    """Lifecycle states for input preparation and filesystem commit."""

    IDLE = "idle"
    PREPARING = "preparing"
    PREPARED = "prepared"
    COMMITTING = "committing"
    READY = "ready"
    FAILED = "failed"


@dataclass(frozen=True)
class Occam1DInputPaths:
    """Resolved files belonging to one complete Occam1D input set."""

    workdir: Path
    data: Path
    model: Path
    startup: Path
    manifest: Path

    def to_dict(self) -> dict[str, str]:
        """Return JSON-friendly absolute path strings."""
        return {
            "workdir": str(self.workdir),
            "data": str(self.data),
            "model": str(self.model),
            "startup": str(self.startup),
            "manifest": str(self.manifest),
        }


def _safe_station_name(value) -> str:
    """Return a portable, non-traversing station directory name."""
    name = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("._")
    return name or "site"


def _station_name(site, fallback="site") -> str:
    """Return a stable station label across supported site dialects."""
    for name in ("name", "station", "dataid"):
        value = getattr(site, name, None)
        if value is not None and str(value).strip():
            return str(value).strip()
    return str(fallback)


def _thread_limit_context(n_threads):
    """Return a context manager capping BLAS threads inside joblib workers.

    Uses ``joblib.parallel_config(inner_max_num_threads=...)`` when
    available; degrades to a no-op ``contextlib.nullcontext`` otherwise (an
    older joblib still runs correctly, just with the oversubscription risk
    this exists to avoid).
    """
    import contextlib

    try:
        from joblib import parallel_config
    except ImportError:
        return contextlib.nullcontext()
    try:
        return parallel_config(
            backend="loky", inner_max_num_threads=n_threads
        )
    except TypeError:
        return contextlib.nullcontext()


def _invert_one_station(
    workdir,
    *,
    export_text=True,
    text_dirname="model-text",
    export_images=False,
    image_dirname="model-image",
    dpi=180,
    inversion_kwargs=None,
):
    """Run one already-built station's inversion to completion.

    Module-level and disk-driven by design: it re-reads data, model,
    startup, and config from the files :meth:`Occam1DBatch.build_all`
    already committed to ``workdir`` instead of accepting live objects, so
    it stays picklable for a joblib process pool (Occam1D objects carry a
    logger, which does not pickle) and behaves identically whether it runs
    in-process or in a worker.

    Returns
    -------
    tuple
        ``(station, summary, error)``. ``error`` is ``None`` on success and
        a human-readable message otherwise; ``summary`` is ``None`` on
        failure.
    """
    from .config import Occam1DConfig
    from .data import Occam1DData
    from .inversion import Occam1DInversion
    from .model import Occam1DModel
    from .startup import Occam1DStartup

    workdir = Path(workdir)
    station = workdir.name
    try:
        manifest = json.loads(
            (workdir / _MANIFEST_NAME).read_text(encoding="utf8")
        )
        config = Occam1DConfig(**manifest["config"])
        data = Occam1DData.read(workdir / config.data_file)
        model = Occam1DModel.read(workdir / config.model_file)
        startup = Occam1DStartup.read(workdir / config.startup_file)
        station = data.station
        inversion = Occam1DInversion(
            data,
            model,
            config=config,
            startup=startup,
            verbose=0,
            **(inversion_kwargs or {}),
        )
        result = inversion.run()
        restart = inversion.restart(result)
        restart.write(workdir / "native-restart.json")
        if export_text:
            text_dir = workdir / text_dirname
            inversion.export_result(text_dir, result)
        if export_images:
            import matplotlib

            matplotlib.use("Agg")
            image_dir = workdir / image_dirname
            inversion.save_main_images(image_dir, result, dpi=dpi)
    except Exception as error:  # noqa: BLE001 - reported, not swallowed
        return station, None, f"{type(error).__name__}: {error}"
    summary = {
        "station": station,
        "workdir": str(workdir),
        "status": result.convergence.value,
        "converged": result.converged,
        "iterations": result.n_iterations,
        "initial_rms": result.initial.rms,
        "final_rms": result.final.rms,
        "target_rms": result.target_rms,
    }
    return station, summary, None


def _sha256(path: Path) -> str:
    """Return the hexadecimal SHA-256 digest of one file."""
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(65536), b""):
            digest.update(block)
    return digest.hexdigest()


class Occam1DInputBuilder(Occam1DBase):
    """Prepare and transactionally write one Occam1D inversion input set.

    Parameters
    ----------
    source : object or collection
        Site, EDI-like object, path, or site collection accepted by
        :func:`normalise_sites`.
    workdir : path-like, default="."
        Destination for data, model, startup, and manifest files.
    config : Occam1DConfig, optional
        Base configuration. A defensive validated copy is retained; the
        caller's object is never mutated.
    station : str or int, optional
        Station name or non-negative source index. Required when ``source``
        contains more than one site.
    verbose, logger, path, metadata, stream
        Shared options inherited from :class:`Occam1DBase`.

    Attributes
    ----------
    build_state : Occam1DBuildState
        Preparation/commit state, independent of base path lifecycle.
    data, model, startup
        Prepared scientific objects, or ``None`` before successful prepare.
    paths : Occam1DInputPaths
        Expected resolved output paths.
    manifest : dict
        Checksummed reproducibility metadata after commit.

    Notes
    -----
    :meth:`prepare` performs no output writes. :meth:`commit` writes into a
    temporary staging directory, reparses all native files, then replaces the
    destination files. Existing files are backed up during replacement and
    restored if a later replacement fails. Unrelated workdir files are never
    removed.

    :meth:`build` is the convenient prepare-and-commit operation retained for
    normal use.

    Examples
    --------
    >>> builder = Occam1DInputBuilder(site, "occam1d-inversion")  # noqa: F821
    >>> builder.prepare(mode="xy")  # doctest: +SKIP
    >>> builder.commit()  # doctest: +SKIP
    >>> builder.is_ready  # doctest: +SKIP
    True
    """

    def __init__(
        self,
        source,
        workdir=".",
        config=None,
        *,
        station=None,
        **kwargs,
    ):
        super().__init__(**kwargs)
        if config is not None and not isinstance(config, Occam1DConfig):
            raise TypeError("config must be an Occam1DConfig or None.")
        if isinstance(station, bool) or not isinstance(
            station, (str, int, type(None))
        ):
            raise TypeError("station must be a name, index, or None.")
        if isinstance(station, int) and station < 0:
            raise ValueError("station index must be non-negative.")
        self.source = source
        self.workdir = self._normalize_path(workdir)
        self.config = (config or Occam1DConfig()).copy()
        self.station = station
        self.data: Occam1DData | None = None
        self.model: Occam1DModel | None = None
        self.startup: Occam1DStartup | None = None
        self.selected_site = None
        self.build_state = Occam1DBuildState.IDLE
        self.manifest: dict[str, Any] = {}

    @property
    def paths(self) -> Occam1DInputPaths:
        """Resolved expected files for the effective configuration."""
        return Occam1DInputPaths(
            workdir=self.workdir,
            data=self.workdir / self.config.data_file,
            model=self.workdir / self.config.model_file,
            startup=self.workdir / self.config.startup_file,
            manifest=self.workdir / _MANIFEST_NAME,
        )

    def _select_site(self):
        """Resolve one unambiguous site from the configured source."""
        sites = normalise_sites(self.source)
        if not sites:
            raise ValueError("source did not produce any usable sites.")
        if self.station is None:
            if len(sites) != 1:
                raise ValueError(
                    "source contains multiple sites; select station by name "
                    "or index."
                )
            return sites[0]
        if isinstance(self.station, int):
            if self.station >= len(sites):
                raise IndexError("station index is outside the source.")
            return sites[self.station]
        wanted = self.station.casefold()
        matches = [
            site
            for site in sites
            if _station_name(site, fallback="").casefold() == wanted
        ]
        if not matches:
            raise KeyError(f"Unknown station: {self.station!r}")
        if len(matches) > 1:
            raise ValueError(
                f"Station name {self.station!r} is not unique in source."
            )
        return matches[0]

    def prepare(self, **overrides):
        """Build validated scientific objects without writing output files.

        Configuration overrides are validated on a copy. The builder adopts
        that effective copy only after all three objects are successfully
        constructed.
        """
        self.build_state = Occam1DBuildState.PREPARING
        try:
            unknown = set(overrides).difference(
                self.config.__dataclass_fields__
            )
            if unknown:
                names = ", ".join(sorted(unknown))
                raise TypeError(
                    f"Unknown Occam1D configuration fields: {names}."
                )
            effective = self.config.updated(**overrides)
            if _MANIFEST_NAME in effective.native_filenames:
                raise ValueError(
                    "Native input filenames cannot use reserved manifest "
                    f"name {_MANIFEST_NAME!r}."
                )
            site = self._select_site()
            data = extract_sounding(site, effective)
            model = Occam1DModel.build(
                effective.n_layers,
                effective.first_thickness,
                effective.depth_max,
                resistivity=effective.starting_resistivity,
            )
            startup = Occam1DStartup.from_model(model, effective)
            startup.apply_to_model(model)
        except Exception:
            self.build_state = Occam1DBuildState.FAILED
            self.logger.exception("Could not prepare Occam1D input objects.")
            raise
        self.config = effective
        self.selected_site = site
        self.data = data
        self.model = model
        self.startup = startup
        self.manifest = {}
        self.build_state = Occam1DBuildState.PREPARED
        self.logger.debug(
            "Prepared Occam1D inputs for %s: %d data, %d layers.",
            data.station,
            data.n_data,
            model.n_layers,
        )
        return self

    def commit(self, *, overwrite=True):
        """Validate staged files and commit them to ``workdir``.

        Parameters
        ----------
        overwrite : bool, default=True
            Replace existing managed input/manifest files. ``False`` fails
            before staging if any managed destination exists.
        """
        if not isinstance(overwrite, bool):
            raise TypeError("overwrite must be a boolean.")
        if self.build_state is not Occam1DBuildState.PREPARED:
            raise RuntimeError("Call prepare() successfully before commit().")
        if self.workdir.exists() and not self.workdir.is_dir():
            raise NotADirectoryError(
                f"Occam1D workdir is not a directory: {self.workdir}"
            )
        final = self.paths
        managed = (final.data, final.model, final.startup, final.manifest)
        if not overwrite:
            existing = [path for path in managed if path.exists()]
            if existing:
                names = ", ".join(path.name for path in existing)
                raise FileExistsError(f"Managed output files exist: {names}.")

        self.build_state = Occam1DBuildState.COMMITTING
        self.workdir.parent.mkdir(parents=True, exist_ok=True)
        try:
            with tempfile.TemporaryDirectory(
                prefix=".occam1d-stage-", dir=self.workdir.parent
            ) as directory:
                stage = Path(directory)
                staged = self._write_stage(stage)
                self._verify_stage(staged)
                manifest = self._build_manifest(staged)
                staged.manifest.write_text(
                    json.dumps(manifest, indent=2, sort_keys=True) + "\n",
                    encoding="utf8",
                )
                self.workdir.mkdir(parents=True, exist_ok=True)
                self._commit_files(staged, final, overwrite=overwrite)
        except Exception:
            self.build_state = Occam1DBuildState.FAILED
            self.logger.exception(
                "Could not commit Occam1D inputs to %s.", self.workdir
            )
            raise
        self.data._bind_path(final.data)
        self.model._bind_path(final.model)
        self.startup._bind_path(final.startup)
        self.manifest = manifest
        self._bind_path(self.workdir)
        self.build_state = Occam1DBuildState.READY
        self.logger.info(
            "Built Occam1D inputs for %s in %s.",
            self.data.station,
            self.workdir,
        )
        self._emit(
            f"Occam1D inputs for {self.data.station} written to "
            f"{self.workdir}"
        )
        return self

    def _write_stage(self, stage: Path) -> Occam1DInputPaths:
        """Serialize prepared objects inside one staging directory."""
        paths = Occam1DInputPaths(
            workdir=stage,
            data=stage / self.config.data_file,
            model=stage / self.config.model_file,
            startup=stage / self.config.startup_file,
            manifest=stage / _MANIFEST_NAME,
        )
        self.data.write(paths.data)
        self.model.write(paths.model)
        self.startup.write(paths.startup)
        return paths

    @staticmethod
    def _verify_stage(paths: Occam1DInputPaths) -> None:
        """Reparse staged files and cross-check model/startup compatibility."""
        data = Occam1DData.read(paths.data)
        model = Occam1DModel.read(paths.model)
        startup = Occam1DStartup.read(paths.startup)
        startup.apply_to_model(model)
        if startup.data_file != paths.data.name:
            raise ValueError(
                "Staged startup data reference does not match staged data."
            )
        if startup.model_file != paths.model.name:
            raise ValueError(
                "Staged startup model reference does not match staged model."
            )
        if data.n_data < 1:
            raise ValueError("Staged data contain no inversion observations.")

    def _build_manifest(self, paths: Occam1DInputPaths) -> dict[str, Any]:
        """Build deterministic configuration and checksum provenance."""
        return {
            "schema": "pycsamt.occam1d.input-manifest/v1",
            "station": self.data.station,
            "mode": self.data.mode,
            "config": self.config.to_dict(),
            "data": {
                "n_frequencies": self.data.n_frequencies,
                "n_observations": self.data.n_data,
                "frequency_bounds_hz": list(self.data.frequency_bounds),
            },
            "model": {
                "n_layers": self.model.n_layers,
                "depth_max_m": self.model.depth_max,
            },
            "files": {
                self.config.data_file: _sha256(paths.data),
                self.config.model_file: _sha256(paths.model),
                self.config.startup_file: _sha256(paths.startup),
            },
        }

    @staticmethod
    def _commit_files(staged, final, *, overwrite):
        """Replace managed files and restore prior versions on failure."""
        pairs = (
            (staged.data, final.data),
            (staged.model, final.model),
            (staged.startup, final.startup),
            (staged.manifest, final.manifest),
        )
        backup = staged.workdir / ".backup"
        backup.mkdir()
        previous = {}
        for _, destination in pairs:
            if destination.exists():
                if destination.is_dir():
                    raise IsADirectoryError(
                        f"Managed output path is a directory: {destination}"
                    )
                if not overwrite:
                    raise FileExistsError(
                        f"Managed output file exists: {destination}"
                    )
                saved = backup / destination.name
                shutil.copy2(destination, saved)
                previous[destination] = saved
        committed = []
        try:
            for source, destination in pairs:
                os.replace(source, destination)
                committed.append(destination)
        except Exception:
            for destination in reversed(committed):
                saved = previous.get(destination)
                if saved is not None and saved.exists():
                    os.replace(saved, destination)
                elif destination.exists() and destination.is_file():
                    destination.unlink()
            raise

    def build(self, *, overwrite=True, **overrides):
        """Prepare and commit the complete Occam1D input set."""
        return self.prepare(**overrides).commit(overwrite=overwrite)

    @property
    def is_ready(self) -> bool:
        """Whether objects and committed managed files are complete."""
        if self.build_state is not Occam1DBuildState.READY:
            return False
        objects = (self.data, self.model, self.startup)
        files = (
            self.paths.data,
            self.paths.model,
            self.paths.startup,
            self.paths.manifest,
        )
        return all(item is not None for item in objects) and all(
            path.is_file() for path in files
        )

    def summary(self) -> str:
        """Return a human-readable input-build summary."""
        if self.data is None or self.model is None:
            return (
                f"Occam1DInputBuilder(state={self.build_state.value!r}, "
                f"workdir={str(self.workdir)!r})"
            )
        return (
            "Occam1D input summary\n"
            f"  state      : {self.build_state.value}\n"
            f"  station    : {self.data.station}\n"
            f"  mode       : {self.data.mode}\n"
            f"  frequencies: {self.data.n_frequencies}\n"
            f"  data       : {self.data.n_data}\n"
            f"  layers     : {self.model.n_layers}\n"
            f"  depth max  : {self.model.depth_max:g} m\n"
            f"  workdir    : {self.workdir}\n"
        )

    def diagnostics(self) -> dict[str, Any]:
        """Extend lifecycle diagnostics with build and manifest state."""
        values = super().diagnostics()
        values.update(
            {
                "workdir": str(self.workdir),
                "build_state": self.build_state.value,
                "station_selector": self.station,
                "selected_station": (
                    self.data.station if self.data is not None else None
                ),
                "is_ready": self.is_ready,
                "paths": self.paths.to_dict(),
                "config": self.config.to_dict(),
                "manifest": dict(self.manifest),
            }
        )
        return values


class Occam1DBatch(Occam1DBase):
    """Build isolated Occam1D directories for a station collection.

    Each child receives its own configuration copy. Sanitized station names
    are planned before any write; collisions raise instead of merging two
    stations into one directory.
    """

    def __init__(
        self,
        source,
        workdir="occam1d-inversion",
        config=None,
        **kwargs,
    ):
        super().__init__(**kwargs)
        if config is not None and not isinstance(config, Occam1DConfig):
            raise TypeError("config must be an Occam1DConfig or None.")
        self.source = source
        self.workdir = self._normalize_path(workdir)
        self.config = (config or Occam1DConfig()).copy()
        self.builders: list[Occam1DInputBuilder] = []
        self.failures: dict[str, str] = {}
        self.build_state = Occam1DBuildState.IDLE

    def build_all(
        self,
        *,
        overwrite=True,
        continue_on_error=False,
        **overrides,
    ):
        """Build every source station and return the populated batch.

        Parameters
        ----------
        overwrite : bool, default=True
            Child managed-file replacement policy.
        continue_on_error : bool, default=False
            Continue remaining stations after a child failure and record its
            exception in :attr:`failures`. The batch is not ready when any
            child fails.
        **overrides
            Validated configuration overrides applied independently per child.
        """
        if not isinstance(overwrite, bool):
            raise TypeError("overwrite must be a boolean.")
        if not isinstance(continue_on_error, bool):
            raise TypeError("continue_on_error must be a boolean.")
        sites = normalise_sites(self.source)
        if not sites:
            raise ValueError("source did not produce any usable sites.")
        unknown = set(overrides).difference(self.config.__dataclass_fields__)
        if unknown:
            names = ", ".join(sorted(unknown))
            raise TypeError(
                f"Unknown Occam1D configuration fields: {names}."
            )
        effective = self.config.updated(**overrides)
        plan = []
        used = {}
        for index, site in enumerate(sites, start=1):
            name = _station_name(site, fallback=index)
            directory = _safe_station_name(name)
            if directory.casefold() in used:
                prior = used[directory.casefold()]
                raise ValueError(
                    "Station directory collision after sanitization: "
                    f"{prior!r} and {name!r} both map to {directory!r}."
                )
            used[directory.casefold()] = name
            plan.append((name, directory, site))

        self.builders = []
        self.failures = {}
        self.build_state = Occam1DBuildState.PREPARING
        total = len(plan)
        for index, (name, directory, site) in enumerate(plan, start=1):
            builder = Occam1DInputBuilder(
                site,
                self.workdir / directory,
                effective,
                verbose=0,
                logger=self.logger,
            )
            try:
                builder.build(overwrite=overwrite)
            except Exception as error:
                self.failures[name] = f"{type(error).__name__}: {error}"
                self.logger.exception(
                    "Occam1D batch build failed for station %s.", name
                )
                if not continue_on_error:
                    self.build_state = Occam1DBuildState.FAILED
                    raise
            else:
                self.builders.append(builder)
                self._emit(
                    f"Occam1D build {index}/{total}: {name}", level=1
                )
        self.config = effective
        self.build_state = (
            Occam1DBuildState.FAILED
            if self.failures
            else Occam1DBuildState.READY
        )
        if not self.failures:
            self._bind_path(self.workdir)
        self.logger.info(
            "Built %d/%d Occam1D station directories.",
            len(self.builders),
            total,
        )
        return self

    def invert_all(
        self,
        *,
        n_jobs=1,
        export_text=True,
        text_dirname="model-text",
        export_images=False,
        image_dirname="model-image",
        dpi=180,
        continue_on_error=True,
        **inversion_kwargs,
    ):
        """Run the nonlinear inversion for every built station.

        Each station is fully independent -- no lateral coupling exists in
        Occam1D -- so this is the natural place to spend extra cores.
        ``n_jobs != 1`` runs stations through a :mod:`joblib` process pool
        (soft dependency; falls back to a sequential loop with a warning
        when joblib is not installed). Processes, not threads, because a
        station inversion is CPU-bound Python-level work (candidate search,
        bookkeeping), not I/O -- the GIL would otherwise serialize it.

        Every station re-reads its own data/model/startup/config from the
        files :meth:`build_all` already wrote, so live Occam1D objects
        (which hold a non-picklable logger) never cross a process boundary.

        Parameters
        ----------
        n_jobs : int, default=1
            Worker processes. ``1`` runs sequentially in this process.
            ``-1`` uses all CPU cores.
        export_text, export_images : bool
            Write :meth:`~.inversion.Occam1DInversion.export_result` /
            :meth:`~.inversion.Occam1DInversion.save_main_images` products
            into each station's ``workdir``. Images are off by default:
            rendering is typically slower than the inversion itself and
            better run once, after inspecting results.
        text_dirname, image_dirname : str
            Directory names created inside each station's ``workdir``.
        dpi : int, default=180
            Resolution used when ``export_images=True``.
        continue_on_error : bool, default=True
            Continue remaining stations after one fails and record its
            error in the returned mapping's ``"failures"`` key.
        **inversion_kwargs
            Forwarded to every :class:`~.inversion.Occam1DInversion`
            constructor (e.g. ``solver_policy``, ``rms_tolerance``). Must be
            picklable when ``n_jobs != 1``.

        Returns
        -------
        dict
            ``{"results": {station: summary}, "failures": {station: message}}``.
        """
        if not self.builders:
            raise RuntimeError("The batch has no completed station builders.")
        if not isinstance(continue_on_error, bool):
            raise TypeError("continue_on_error must be a boolean.")
        workdirs = [str(builder.workdir) for builder in self.builders]
        jobs = dict(
            export_text=export_text,
            text_dirname=text_dirname,
            export_images=export_images,
            image_dirname=image_dirname,
            dpi=dpi,
            inversion_kwargs=inversion_kwargs or None,
        )

        results: dict[str, Any] = {}
        failures: dict[str, str] = {}
        total = len(workdirs)
        with get_progress_bar(
            total=total,
            desc="Occam1D batch inversion",
            unit="station",
            verbose=self.verbose,
        ) as bar:
            outcomes = None
            if n_jobs != 1:
                try:
                    import inspect

                    from joblib import Parallel, delayed
                except ImportError:
                    self.add_warning(
                        "joblib is not installed; running the batch "
                        "inversion sequentially instead of n_jobs="
                        f"{n_jobs}.",
                        log=True,
                    )
                else:
                    tasks = [
                        delayed(_invert_one_station)(workdir, **jobs)
                        for workdir in workdirs
                    ]
                    # Each worker process would otherwise let OpenBLAS/MKL
                    # spawn its own multi-threaded pool for the per-iteration
                    # linear solves, oversubscribing the machine once
                    # n_jobs processes each try to use every core. One BLAS
                    # thread per worker process is the standard fix.
                    thread_limit = _thread_limit_context(1)
                    supports_generator = "return_as" in inspect.signature(
                        Parallel.__init__
                    ).parameters
                    with thread_limit:
                        if supports_generator:
                            outcomes = []
                            for outcome in Parallel(
                                n_jobs=n_jobs,
                                prefer="processes",
                                return_as="generator_unordered",
                            )(tasks):
                                outcomes.append(outcome)
                                bar.update(1)
                        else:
                            # No incremental updates: this joblib is too
                            # old to stream results as they complete.
                            outcomes = Parallel(
                                n_jobs=n_jobs, prefer="processes"
                            )(tasks)
                            bar.update(total)
            if outcomes is None:
                outcomes = []
                for workdir in workdirs:
                    outcomes.append(_invert_one_station(workdir, **jobs))
                    bar.update(1)

        for station, summary, error in outcomes:
            if error is not None:
                failures[station] = error
                self.logger.error(
                    "Occam1D inversion failed for station %s: %s",
                    station,
                    error,
                )
                if not continue_on_error:
                    raise RuntimeError(
                        f"Occam1D inversion failed for station {station}: "
                        f"{error}"
                    )
            else:
                results[station] = summary
        self.logger.info(
            "Inverted %d/%d Occam1D stations (%d failures).",
            len(results),
            total,
            len(failures),
        )
        return {"results": results, "failures": failures}

    @property
    def is_ready(self) -> bool:
        """Whether all planned stations committed complete input sets."""
        return (
            self.build_state is Occam1DBuildState.READY
            and bool(self.builders)
            and not self.failures
            and all(builder.is_ready for builder in self.builders)
        )

    def export_images(self, directory="model-image", *, iteration="best"):
        """Save main images for each completed inversion result.

        A relative ``directory`` is resolved beneath the batch workdir.
        Input building alone does not create responses or convergence logs;
        call this only after the station inversions have completed.
        """
        from .results import Occam1DResult

        if not self.builders:
            raise RuntimeError("The batch has no completed station builders.")
        output = Path(directory)
        if not output.is_absolute():
            output = self.workdir / output
        paths = {}
        for builder in self.builders:
            result = Occam1DResult(builder.workdir, iteration=iteration)
            paths[result.data.station] = result.save_main_images(output)
        return paths

    def summary(self) -> str:
        """Return a concise batch-build summary."""
        return (
            f"Occam1DBatch(state={self.build_state.value!r}, "
            f"completed={len(self.builders)}, failures={len(self.failures)}, "
            f"workdir={str(self.workdir)!r})"
        )

    def diagnostics(self) -> dict[str, Any]:
        """Extend lifecycle diagnostics with batch outcomes."""
        values = super().diagnostics()
        values.update(
            {
                "workdir": str(self.workdir),
                "build_state": self.build_state.value,
                "is_ready": self.is_ready,
                "completed": len(self.builders),
                "stations": [
                    builder.data.station for builder in self.builders
                ],
                "failures": dict(self.failures),
                "config": self.config.to_dict(),
            }
        )
        return values
