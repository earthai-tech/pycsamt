# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""ModEM adapter backend.

This module connects :mod:`pycsamt.inversion` to the ModEM workflow helpers
under :mod:`pycsamt.models.modem`. It does not implement the ModEM
electromagnetic solver itself. Instead, it prepares or locates ModEM input
files, validates the external-run lifecycle, optionally calls the runner when
``run_external=True``, and loads finished results back into a common
:class:`pycsamt.inversion.results.InversionResult`.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from ..base import BaseInversionBackend
from ..results import InversionResult

__all__ = ["ModEMBackend"]


class ModEMBackend(BaseInversionBackend):
    name = "modem"
    supports = (
        ("mt", "2d"),
        ("mt", "3d"),
        ("amt", "2d"),
        ("amt", "3d"),
    )

    def run(self, data: Any | None = None) -> InversionResult:
        """Prepare, validate, optionally run, and load ModEM outputs.

        Parameters
        ----------
        data : object, mapping, path-like, or sequence, optional
            Input source used by ``pycsamt.models.modem.InputBuilder`` when it
            is available. When omitted, ``self.config.data`` is used. If no
            source is supplied, configured file names are resolved inside the
            work directory and the backend only validates runner readiness.

        Returns
        -------
        InversionResult
            Backend-neutral result. The ``files`` mapping contains prepared or
            configured ModEM file paths, ``metadata["command"]`` records the
            command that would be executed, and ``native`` stores the builder,
            runner, loaded ModEM result, or another backend-native object.

        Raises
        ------
        ImportError
            If :mod:`pycsamt.models.modem` is unavailable.
        NotImplementedError
            If the configured method/dimension pair is not supported by this
            adapter.

        Notes
        -----
        External execution is disabled unless ``self.config.run_external`` is
        ``True``. With the default value, the method prepares or validates the
        ModEM directory and reports the command that would be executed when the
        runner files are complete.

        The returned status describes the lifecycle stage:

        ``"prepared"``
            Input files were built or configured, but runner readiness was not
            fully established.
        ``"ready"``
            Required files exist and a command string could be assembled.
        ``"executed"``
            The external runner completed without loading a result object.
        ``"loaded"``
            A ModEM result object was found and attached to ``native``.
        ``"needs_review"``
            Preparation, runner execution, or validation raised a recoverable
            issue. Details are stored in ``warnings``.

        Examples
        --------
        >>> from pycsamt.inversion.backends.modem import ModEMBackend
        >>> from pycsamt.inversion.config import InversionConfig
        >>> cfg = InversionConfig(
        ...     method="mt",
        ...     dimension="3d",
        ...     backend="modem",
        ...     workdir="modem_run",
        ...     backend_options={
        ...         "files": {
        ...             "data": "ModEM_Data.dat",
        ...             "model": "ModEM_Model.rho",
        ...             "control": "ModEM_Control.inv",
        ...             "covariance": "covariance.cov",
        ...         }
        ...     },
        ... )
        >>> result = ModEMBackend(cfg).run()  # doctest: +SKIP
        """
        self.check_supported()
        cfg = self.config
        source = cfg.data if data is None else data
        workdir = Path(cfg.workdir)
        workdir.mkdir(parents=True, exist_ok=True)
        warnings: list[str] = []
        native = None
        file_map: dict[str, str] = {}
        status = "prepared"
        command: str | None = None
        executed = False

        try:
            from ...models import modem
        except ImportError as exc:
            raise ImportError(
                "ModEM backend requires pycsamt.models.modem."
            ) from exc

        modem_cfg = _modem_config(modem, cfg)
        builder_cls = getattr(modem, "InputBuilder", None)
        if builder_cls is not None and source is not None:
            try:
                builder = builder_cls(config=modem_cfg)
                build_options = _build_options(cfg.backend_options)
                files = builder.build(
                    source, workdir=workdir, **build_options
                )
                native = builder
                file_map = {key: str(value) for key, value in files.items()}
                status = "prepared"
            except Exception as exc:
                warnings.append(f"ModEM preparation failed: {exc}")
                status = "needs_review"
        else:
            file_map = _configured_files(
                workdir, modem_cfg, cfg.backend_options
            )
            warnings.append(
                "ModEM InputBuilder is not available; created workdir only."
            )

        runner_cls = getattr(modem, "ModEmRunner", None)
        runner = None
        runner_options = dict(cfg.backend_options.get("runner", {}))
        if runner_cls is not None:
            runner = runner_cls(workdir=workdir, config=modem_cfg)
            native = native or runner
            missing = _missing_runner_files(
                file_map, required=_required_runner_keys(modem_cfg)
            )
            if missing:
                warnings.append(
                    "ModEM runner not ready; missing files: "
                    + ", ".join(missing)
                )
            else:
                command = runner.command(
                    _relative_to_workdir(file_map["model"], workdir),
                    _relative_to_workdir(file_map["data"], workdir),
                    _relative_to_workdir(file_map["control"], workdir),
                    covariance=(
                        _relative_to_workdir(file_map["covariance"], workdir)
                        if "covariance" in file_map
                        else None
                    ),
                    mode=modem_cfg.mode,
                    use_mpi=runner_options.get("use_mpi", None),
                    n_procs=runner_options.get("n_procs", None),
                )
                status = "ready" if status == "prepared" else status
                if cfg.run_external:
                    try:
                        loaded = runner.run(
                            _relative_to_workdir(file_map["model"], workdir),
                            _relative_to_workdir(file_map["data"], workdir),
                            _relative_to_workdir(
                                file_map["control"], workdir
                            ),
                            covariance=(
                                _relative_to_workdir(
                                    file_map["covariance"], workdir
                                )
                                if "covariance" in file_map
                                else None
                            ),
                            mode=modem_cfg.mode,
                            use_mpi=runner_options.get("use_mpi", None),
                            n_procs=runner_options.get("n_procs", None),
                            extra_args=runner_options.get("extra_args", None),
                            timeout=runner_options.get("timeout", None),
                            load_result=bool(
                                runner_options.get("load_result", True)
                            ),
                        )
                        native = loaded or runner
                        executed = True
                        status = (
                            "loaded" if loaded is not None else "executed"
                        )
                    except Exception as exc:
                        warnings.append(f"ModEM runner failed: {exc}")
                        status = "needs_review"
        else:
            warnings.append(
                "ModEM ModEmRunner is not available; execution disabled."
            )

        result_cls = getattr(modem, "InversionResult", None)
        rms = float("nan")
        if result_cls is not None and status != "needs_review":
            try:
                loaded = result_cls(workdir, config=modem_cfg)
                if _has_loaded_result(loaded):
                    native = loaded
                    status = "loaded"
                    rms = float(
                        getattr(
                            loaded,
                            "final_rms",
                            getattr(loaded, "rms", float("nan")),
                        )
                    )
            except Exception as exc:
                warnings.append(f"ModEM result loading skipped: {exc}")

        return InversionResult(
            method=cfg.method,
            dimension=cfg.dimension,
            backend=self.name,
            status=status,
            data=source,
            rms=rms,
            workdir=str(workdir),
            files=file_map,
            native=native,
            warnings=warnings,
            metadata={
                **cfg.metadata,
                "engine": "modem",
                "mode": modem_cfg.mode,
                "run_external": cfg.run_external,
                "executed": executed,
                "command": command,
                "binary": modem_cfg.binary_name,
                "use_mpi": modem_cfg.use_mpi,
                "n_procs": modem_cfg.n_procs,
            },
        )


def _modem_config(modem: Any, cfg: Any) -> Any:
    cfg_cls = modem.ModEmConfig
    raw = cfg.backend_options.get("config")
    if raw is None:
        values = {
            key: value
            for key, value in cfg.backend_options.items()
            if key in getattr(cfg_cls, "__dataclass_fields__", {})
        }
        values.setdefault("mode", cfg.dimension)
        return cfg_cls(**values)
    if isinstance(raw, cfg_cls):
        raw.mode = cfg.dimension
        return raw
    if isinstance(raw, dict):
        values = dict(raw)
        values.setdefault("mode", cfg.dimension)
        return cfg_cls(**values)
    if hasattr(raw, "mode"):
        raw.mode = cfg.dimension
    return raw


def _build_options(options: dict[str, Any]) -> dict[str, Any]:
    excluded = {"config", "runner", "files"}
    modem_fields = {
        "mode",
        "binary_2d",
        "binary_3d",
        "use_mpi",
        "n_procs",
        "mpi_command",
    }
    return {
        key: value
        for key, value in options.items()
        if key not in excluded and key not in modem_fields
    }


def _configured_files(
    workdir: Path, modem_cfg: Any, options: dict[str, Any]
) -> dict[str, str]:
    raw = dict(options.get("files", {}))
    aliases = {
        "data_file": "data",
        "model_file": "model",
        "control_file": "control",
        "covariance_file": "covariance",
    }
    for option_key, role in aliases.items():
        if option_key in options and role not in raw:
            raw[role] = options[option_key]
    defaults = {
        "data": modem_cfg.data_file,
        "model": modem_cfg.model_file,
        "control": modem_cfg.control_file,
    }
    if modem_cfg.is_3d:
        defaults["covariance"] = modem_cfg.covariance_file
    for role, name in defaults.items():
        raw.setdefault(role, name)
    return {
        role: str(_resolve_workdir_path(path, workdir))
        for role, path in raw.items()
    }


def _required_runner_keys(modem_cfg: Any) -> tuple[str, ...]:
    if modem_cfg.is_3d:
        return ("model", "data", "control", "covariance")
    return ("model", "data", "control")


def _missing_runner_files(
    files: dict[str, str], *, required: tuple[str, ...]
) -> list[str]:
    missing: list[str] = []
    for key in required:
        path = files.get(key)
        if path is None:
            missing.append(key)
        elif not Path(path).exists():
            missing.append(f"{key}={path}")
    return missing


def _relative_to_workdir(path: str, workdir: Path) -> str:
    p = Path(path)
    try:
        return str(p.resolve().relative_to(workdir.resolve()))
    except Exception:
        return str(p)


def _resolve_workdir_path(path: Any, workdir: Path) -> Path:
    p = Path(path)
    return p if p.is_absolute() else workdir / p


def _has_loaded_result(result: Any) -> bool:
    return bool(
        getattr(result, "models", None)
        or getattr(result, "log", None)
        or getattr(result, "data", None)
        or getattr(result, "rho_2d", None) is not None
    )


ModEMBackend.__doc__ = r"""
Adapt PyCSAMT inversion workflows to ModEM.

``ModEMBackend`` is the external-run adapter for 2-D and 3-D MT/AMT
inversions. It shares the same
:class:`pycsamt.inversion.workflow.InversionWorkflow` interface as the
runnable in-process backends, but the numerical inversion is performed by the
ModEM ecosystem through :mod:`pycsamt.models.modem`.

The backend handles four lifecycle steps:

1. normalize ModEM configuration from ``backend_options``;
2. prepare input files with ``InputBuilder`` when a data source is supplied;
3. validate required runner files and assemble the external command; and
4. optionally run ModEM and load finished result files.

Parameters
----------
config : pycsamt.inversion.config.InversionConfig
    Inversion configuration. Important fields are ``method``, ``dimension``,
    ``workdir``, ``data``, ``backend_options``, ``run_external``, and
    ``metadata``. The adapter supports MT and AMT in 2-D or 3-D mode.

Attributes
----------
name : str
    Backend registry name, always ``"modem"``.
supports : tuple of tuple
    Supported ``(method, dimension)`` pairs.

Notes
-----
``backend_options`` can contain:

``config``
    A ``ModEmConfig`` instance or dictionary used to configure file names,
    binary names, MPI settings, and mode-specific options.
``files``
    Mapping from roles such as ``"data"``, ``"model"``, ``"control"``, and
    ``"covariance"`` to existing or expected file paths. Relative paths are
    resolved inside ``workdir``.
``runner``
    Runner-only options including ``use_mpi``, ``n_procs``, ``extra_args``,
    ``timeout``, and ``load_result``.
``data_file``, ``model_file``, ``control_file``, ``covariance_file``
    Convenience aliases for file names when a full ``files`` mapping is not
    supplied.

By default the backend prepares and validates files but does not execute the
external binary. Set ``run_external=True`` in
:class:`pycsamt.inversion.config.InversionConfig` to allow execution.

Examples
--------
Validate an existing 3-D ModEM run directory without execution::

    >>> from pycsamt.inversion.backends.modem import ModEMBackend
    >>> from pycsamt.inversion.config import InversionConfig
    >>> cfg = InversionConfig(
    ...     method="mt",
    ...     dimension="3d",
    ...     backend="modem",
    ...     workdir="modem_run",
    ...     backend_options={
    ...         "files": {
    ...             "data": "ModEM_Data.dat",
    ...             "model": "ModEM_Model.rho",
    ...             "control": "ModEM_Control.inv",
    ...             "covariance": "covariance.cov",
    ...         }
    ...     },
    ... )
    >>> result = ModEMBackend(cfg).run()  # doctest: +SKIP
    >>> result.status in {"ready", "loaded", "needs_review"}  # doctest: +SKIP
    True

Use the high-level inversion workflow and request execution::

    >>> from pycsamt.inversion.workflow import run_inversion
    >>> result = run_inversion(
    ...     method="amt",
    ...     dimension="2d",
    ...     backend="modem",
    ...     workdir="modem_2d",
    ...     run_external=True,
    ...     backend_options={
    ...         "runner": {"use_mpi": False, "timeout": 600},
    ...         "files": {
    ...             "data": "ModEM_Data.dat",
    ...             "model": "ModEM_Model.rho",
    ...             "control": "ModEM_Control.inv",
    ...         },
    ...     },
    ... )  # doctest: +SKIP

See Also
--------
pycsamt.inversion.workflow.InversionWorkflow
    High-level entry point that instantiates this backend.
pycsamt.inversion.backends.occam2d.Occam2DBackend
    Similar external-run adapter for Occam2D.
pycsamt.inversion.results.InversionResult
    Backend-neutral result returned by :meth:`run`.

References
----------
.. [1] Egbert, G. D. and Kelbert, A. (2012). Computational recipes for
       electromagnetic inverse problems. *Geophysical Journal International*,
       189(1), 251-267.
.. [2] Kelbert, A., Meqbel, N., Egbert, G. D. and Tandon, K. (2014). ModEM:
       A modular system for inversion of electromagnetic geophysical data.
       *Computers & Geosciences*, 66, 40-53.
"""
