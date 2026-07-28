# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Occam2D adapter backend.

This module connects :mod:`pycsamt.inversion` to the Occam2D workflow helpers
under :mod:`pycsamt.models.occam2d`. It prepares Occam2D input directories,
validates runner files, optionally executes the external binary when
``run_external=True``, and loads finished Occam2D results into the common
:class:`pycsamt.inversion.results.InversionResult` container.
"""

from __future__ import annotations

import shlex
from pathlib import Path
from typing import Any

from ..base import BaseInversionBackend
from ..results import InversionResult

__all__ = ["Occam2DBackend"]


class Occam2DBackend(BaseInversionBackend):
    name = "occam2d"
    supports = (
        ("mt", "2d"),
        ("amt", "2d"),
        ("csamt", "2d"),
        ("emap", "2d"),
    )

    def run(self, data: Any | None = None) -> InversionResult:
        """Prepare, validate, optionally run, and load Occam2D outputs.

        Parameters
        ----------
        data : object, mapping, path-like, or sequence, optional
            Input source used by ``pycsamt.models.occam2d.InputBuilder``. When
            omitted, ``self.config.data`` is used. If no source is supplied,
            configured file names are resolved inside ``workdir`` and the
            backend validates whether the runner can be launched.

        Returns
        -------
        InversionResult
            Backend-neutral result. The ``files`` mapping contains data, mesh,
            model, and startup file paths when known. ``metadata`` stores the
            command, startup file, binary name, external execution flag, exit
            code, and whether a runner was executed.

        Raises
        ------
        ImportError
            If :mod:`pycsamt.models.occam2d` is unavailable.
        NotImplementedError
            If the configured method/dimension pair is not supported.

        Notes
        -----
        External execution is disabled unless ``self.config.run_external`` is
        ``True``. With the default value, this method prepares or validates the
        directory and reports the command that would be used when all required
        files exist.

        Status values follow the adapter lifecycle:

        ``"prepared"``
            Files were prepared or configured, but the runner is not ready.
        ``"ready"``
            Required files exist and an Occam2D command was assembled.
        ``"executed"``
            The external runner returned exit code ``0``.
        ``"loaded"``
            A finished Occam2D result with ``rho_2d`` was loaded.
        ``"needs_review"``
            Preparation, runner execution, or result loading needs user review.

        Examples
        --------
        >>> from pycsamt.inversion.backends.occam2d import Occam2DBackend
        >>> from pycsamt.inversion.config import InversionConfig
        >>> cfg = InversionConfig(
        ...     method="mt",
        ...     dimension="2d",
        ...     backend="occam2d",
        ...     workdir="occam_run",
        ...     backend_options={
        ...         "files": {
        ...             "data": "OccamData.dat",
        ...             "mesh": "OccamMesh",
        ...             "model": "OccamModel",
        ...             "startup": "Startup",
        ...         }
        ...     },
        ... )
        >>> result = Occam2DBackend(cfg).run()  # doctest: +SKIP
        """
        self.check_supported()
        cfg = self.config
        source = cfg.data if data is None else data
        workdir = Path(cfg.workdir)
        workdir.mkdir(parents=True, exist_ok=True)
        warnings: list[str] = []
        files: dict[str, str] = {}
        native = None
        status = "prepared"
        command: str | None = None
        executed = False
        exit_code: int | None = None

        try:
            from ...models.occam2d import (
                InputBuilder,
                OccamConfig,
                OccamRunner,
            )
            from ...models.occam2d import (
                InversionResult as OccamResult,
            )
        except ImportError as exc:
            raise ImportError(
                "Occam2D backend requires pycsamt.models.occam2d."
            ) from exc

        occam_cfg = _occam_config(OccamConfig, cfg)
        if source is not None:
            try:
                builder = InputBuilder(source, workdir=workdir, config=occam_cfg)
                build_options = _build_options(cfg.backend_options)
                build_options.setdefault("error_floor_rho", cfg.error_floor)
                build_options.setdefault("error_floor_phase", cfg.phase_error)
                builder.build(**build_options)
                native = builder
                files = _builder_files(builder)
            except Exception as exc:
                warnings.append(f"Occam2D preparation failed: {exc}")
                status = "needs_review"
        else:
            files = _configured_files(workdir, occam_cfg, cfg.backend_options)

        runner_options = dict(cfg.backend_options.get("runner", {}))
        missing = _missing_runner_files(files)
        runner = OccamRunner(
            workdir,
            binary_path=runner_options.get(
                "binary_path", cfg.backend_options.get("binary_path")
            ),
            startup_file=runner_options.get("startup_file", occam_cfg.startup_file),
        )
        native = native or runner
        if missing:
            warnings.append(
                "Occam2D runner not ready; missing files: " + ", ".join(missing)
            )
        else:
            binary_name = str(
                runner_options.get(
                    "binary_path",
                    cfg.backend_options.get("binary_path", occam_cfg.binary_name),
                )
            )
            command = _command_string(binary_name, runner.startup_file)
            status = "ready" if status == "prepared" else status

            if cfg.run_external:
                try:
                    exit_code = runner.run(
                        max_iter=runner_options.get("max_iter", None),
                        target_misfit=runner_options.get("target_misfit", None),
                        auto_compile=bool(runner_options.get("auto_compile", False)),
                    )
                    native = runner
                    executed = True
                    if exit_code == 0:
                        status = "executed"
                    else:
                        status = "needs_review"
                        warnings.append(
                            f"Occam2D exited with code {exit_code}; see runner logs."
                        )
                except Exception as exc:
                    warnings.append(f"Occam2D runner failed: {exc}")
                    status = "needs_review"

        rms = float("nan")
        try:
            loaded = OccamResult(workdir, config=occam_cfg)
            if getattr(loaded, "rho_2d", None) is not None:
                native = loaded
                status = "loaded"
                rms = float(getattr(loaded, "final_rms", float("nan")))
        except Exception as exc:
            warnings.append(f"Occam2D result loading skipped: {exc}")

        return InversionResult(
            method=cfg.method,
            dimension=cfg.dimension,
            backend=self.name,
            status=status,
            data=source,
            rms=rms,
            workdir=str(workdir),
            files=files,
            native=native,
            warnings=warnings,
            metadata={
                **cfg.metadata,
                "engine": "occam2d",
                "run_external": cfg.run_external,
                "executed": executed,
                "exit_code": exit_code,
                "command": command,
                "binary": occam_cfg.binary_name,
                "startup_file": occam_cfg.startup_file,
            },
        )


def _occam_config(config_cls: Any, cfg: Any) -> Any:
    raw = cfg.backend_options.get("config")
    if raw is None:
        values = {
            key: value
            for key, value in cfg.backend_options.items()
            if key in getattr(config_cls, "__dataclass_fields__", {})
        }
        return config_cls(**values)
    if isinstance(raw, config_cls):
        return raw
    if isinstance(raw, dict):
        return config_cls(**raw)
    return raw


def _build_options(options: dict[str, Any]) -> dict[str, Any]:
    excluded = {"config", "runner", "files", "binary_path"}
    config_fields = {
        "data_file",
        "mesh_file",
        "model_file",
        "startup_file",
        "binary_name",
    }
    return {
        key: value
        for key, value in options.items()
        if key not in excluded and key not in config_fields
    }


def _builder_files(builder: Any) -> dict[str, str]:
    files: dict[str, str] = {}
    defaults = {
        "data": builder.config.data_file,
        "mesh": builder.config.mesh_file,
        "model": builder.config.model_file,
        "startup": builder.config.startup_file,
    }
    for key, obj in {
        "data": builder.data,
        "mesh": builder.mesh,
        "model": builder.model,
        "startup": builder.startup,
    }.items():
        path = getattr(obj, "path", None)
        if path is None:
            path = builder.workdir / defaults[key]
        files[key] = str(path)
    return files


def _configured_files(
    workdir: Path, occam_cfg: Any, options: dict[str, Any]
) -> dict[str, str]:
    raw = dict(options.get("files", {}))
    defaults = {
        "data": occam_cfg.data_file,
        "mesh": occam_cfg.mesh_file,
        "model": occam_cfg.model_file,
        "startup": occam_cfg.startup_file,
    }
    for role, name in defaults.items():
        raw.setdefault(role, name)
    return {
        role: str(_resolve_workdir_path(path, workdir)) for role, path in raw.items()
    }


def _missing_runner_files(files: dict[str, str]) -> list[str]:
    missing: list[str] = []
    for key in ("data", "mesh", "model", "startup"):
        path = files.get(key)
        if path is None:
            missing.append(key)
        elif not Path(path).exists():
            missing.append(f"{key}={path}")
    return missing


def _command_string(binary: str, startup_file: str) -> str:
    return " ".join(shlex.quote(part) for part in [binary, startup_file])


def _resolve_workdir_path(path: Any, workdir: Path) -> Path:
    p = Path(path)
    return p if p.is_absolute() else workdir / p


Occam2DBackend.__doc__ = r"""
Adapt PyCSAMT inversion workflows to Occam2D.

``Occam2DBackend`` is the external-run adapter for smooth 2-D EM inversion
with Occam2D-style workflows. It uses
:mod:`pycsamt.models.occam2d` to prepare input files, assemble runner commands,
optionally launch the external executable, and load completed model sections.
The actual inversion physics are provided by Occam2D, while PyCSAMT keeps a
consistent configuration and result interface.

The adapter supports MT, AMT, CSAMT, and EMAP in 2-D mode. Its main purpose is
to make external Occam2D projects reproducible from
:class:`pycsamt.inversion.config.InversionConfig` and to return a common
:class:`pycsamt.inversion.results.InversionResult`.

Parameters
----------
config : pycsamt.inversion.config.InversionConfig
    Inversion configuration. Important fields are ``method``, ``dimension``,
    ``workdir``, ``data``, ``backend_options``, ``run_external``,
    ``error_floor``, ``phase_error``, and ``metadata``.

Attributes
----------
name : str
    Backend registry name, always ``"occam2d"``.
supports : tuple of tuple
    Supported ``(method, dimension)`` pairs.

Notes
-----
``backend_options`` can contain:

``config``
    An ``OccamConfig`` instance or dictionary used to configure file names and
    binary settings.
``files``
    Mapping from roles ``"data"``, ``"mesh"``, ``"model"``, and
    ``"startup"`` to existing or expected paths. Relative paths are resolved
    inside ``workdir``.
``runner``
    Runner-only options such as ``binary_path``, ``startup_file``,
    ``max_iter``, ``target_misfit``, and ``auto_compile``.
``binary_path``
    Convenience runner executable path when not supplied inside
    ``runner``.

When a data source is supplied, ``InputBuilder`` receives the source and the
configured error floors. Without a data source, the adapter treats
``backend_options["files"]`` as an existing Occam2D project and only validates
whether it is ready to run or load.

Examples
--------
Validate an existing Occam2D directory without executing the binary::

    >>> from pycsamt.inversion.backends.occam2d import Occam2DBackend
    >>> from pycsamt.inversion.config import InversionConfig
    >>> cfg = InversionConfig(
    ...     method="mt",
    ...     dimension="2d",
    ...     backend="occam2d",
    ...     workdir="occam_run",
    ...     backend_options={
    ...         "files": {
    ...             "data": "OccamData.dat",
    ...             "mesh": "OccamMesh",
    ...             "model": "OccamModel",
    ...             "startup": "Startup",
    ...         }
    ...     },
    ... )
    >>> result = Occam2DBackend(cfg).run()  # doctest: +SKIP
    >>> result.status in {"ready", "loaded", "needs_review"}  # doctest: +SKIP
    True

Run through the high-level workflow and allow external execution::

    >>> from pycsamt.inversion.workflow import run_inversion
    >>> result = run_inversion(
    ...     method="csamt",
    ...     dimension="2d",
    ...     backend="occam2d",
    ...     workdir="occam_csamt",
    ...     run_external=True,
    ...     backend_options={
    ...         "runner": {
    ...             "binary_path": "Occam2D",
    ...             "startup_file": "Startup",
    ...             "target_misfit": 1.0,
    ...         },
    ...         "files": {
    ...             "data": "OccamData.dat",
    ...             "mesh": "OccamMesh",
    ...             "model": "OccamModel",
    ...             "startup": "Startup",
    ...         },
    ...     },
    ... )  # doctest: +SKIP

See Also
--------
pycsamt.inversion.workflow.InversionWorkflow
    High-level entry point that instantiates this backend.
pycsamt.inversion.backends.modem.ModEMBackend
    Similar external-run adapter for ModEM.
pycsamt.inversion.results.InversionResult
    Backend-neutral result returned by :meth:`run`.

References
----------
.. [1] deGroot-Hedlin, C. and Constable, S. (1990). Occam's inversion to
       generate smooth, two-dimensional models from magnetotelluric data.
       *Geophysics*, 55(12), 1613-1624.
.. [2] Constable, S. C., Parker, R. L. and Constable, C. G. (1987). Occam's
       inversion: A practical algorithm for generating smooth models from
       electromagnetic sounding data. *Geophysics*, 52(3), 289-300.
"""
