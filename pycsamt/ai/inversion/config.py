# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Configuration for pyCSAMT 1-D AI-based EM inversion.

The module exposes :class:`InversionConfig`, a dataclass that collects every
tuneable parameter for :class:`~pycsamt.ai.inversion.inv1d.EMInverter1D` —
architecture, training loop, regularisation, checkpointing, and output.

The recommended workflow mirrors the pattern used by ``ModEmConfig`` and
``OccamConfig``:

1. Call :meth:`InversionConfig.write_template` to generate a fully annotated
   source-of-truth file (Python, JSON, or YAML).
2. Edit the file to reflect the desired architecture and training budget.
3. Load the edited file with :meth:`InversionConfig.from_file`.
4. Optionally call :meth:`InversionConfig.validate` to catch range errors.
5. Call :meth:`InversionConfig.to_inverter` to get a ready-to-fit
   :class:`~pycsamt.ai.inversion.inv1d.EMInverter1D`.
6. Call ``inv.fit(dataset, **cfg.to_fit_kwargs())`` to train.

Quick start
-----------
Generate a default template, edit it, train::

    from pycsamt.ai.inversion.config import InversionConfig
    from pycsamt.forward.batch import ForwardDataset

    # 1 — write annotated source-of-truth file
    InversionConfig.write_template("my_inversion.yml")

    # 2 — edit my_inversion.yml …

    # 3 — load and train
    cfg = InversionConfig.from_file("my_inversion.yml")
    cfg.validate()

    ds = ForwardDataset.load("mt1d_train.npz")
    inv = cfg.to_inverter()
    inv.fit(ds, **cfg.to_fit_kwargs())
    inv.save(cfg.checkpoint_path())

Snapshot a fitted inverter for reproducibility::

    cfg = InversionConfig.from_inverter(inv)
    cfg.write_template("snapshot.yml")
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field, fields
from pathlib import Path
from typing import Any, Dict, Optional, Union

from ...models.config_io import (
    ConfigParameter,
    read_config_file,
    write_config_template,
)

__all__ = ["InversionConfig"]


# ── parameter schema ────────────────────────────────────────────────────────

_INVERSION_CONFIG_SCHEMA: list[ConfigParameter] = [
    # ── Architecture ─────────────────────────────────────────────────────────
    ConfigParameter(
        "arch",
        "Neural network architecture for 1-D inversion.  "
        "Accepted values: 'resnet' — residual network (Liu 2021 style, "
        "best accuracy for MT datasets with > 5 000 samples); "
        "'cnn1d' — 1-D convolutional network (Puzyrev 2019 style, "
        "faster to train, good for smaller datasets); "
        "'fcn' — fully-connected network (Moghadas 2020 style, "
        "handles variable input length).",
        "Architecture",
    ),
    ConfigParameter(
        "n_layers",
        "Number of earth layers the network inverts for, including the "
        "halfspace.  The output vector has length 2*n_layers - 1 "
        "(n_layers log-resistivities + n_layers-1 thicknesses).  "
        "Must match the n_layers used during dataset generation.",
        "Architecture",
    ),
    ConfigParameter(
        "solver",
        "EM method the inverter targets.  Accepted values: 'mt1d', "
        "'csamt1d', 'tem1d'.  Determines which default frequency grid "
        "is used when coercing Z objects or ForwardResponse inputs to "
        "numpy feature arrays for prediction.",
        "Architecture",
    ),
    ConfigParameter(
        "device",
        "Compute device for network training and inference.  "
        "Accepted values: 'cpu', 'cuda', 'mps', or null/None.  "
        "When null the backend auto-detects the best available device "
        "(CUDA > MPS > CPU).  Set explicitly to 'cpu' to force CPU-only "
        "execution on a machine with a GPU.",
        "Architecture",
    ),
    ConfigParameter(
        "include_phase",
        "Include impedance phase in the MT/CSAMT input feature vector.  "
        "When True the feature vector length is 2*n_freqs (log10(rho_a) "
        "concatenated with phase in degrees).  When False only log10(rho_a) "
        "is used (feature vector length n_freqs).  Must match the value "
        "used during dataset generation.",
        "Architecture",
    ),
    ConfigParameter(
        "log_thickness",
        "Apply log10 to layer thicknesses in the training targets.  "
        "Strongly recommended when thicknesses span more than two orders of "
        "magnitude.  Must match the setting used when the dataset was "
        "generated if loading a pre-built dataset.",
        "Architecture",
    ),
    ConfigParameter(
        "augment_noise",
        "Standard deviation of Gaussian noise added to each training batch "
        "on-the-fly during every epoch.  Acts as a form of input "
        "regularisation and improves generalisation to field data noise "
        "levels that differ from the training noise.  Set to 0.0 to disable.",
        "Architecture",
    ),
    # ── Training ─────────────────────────────────────────────────────────────
    ConfigParameter(
        "epochs",
        "Maximum number of training epochs.  Early stopping (controlled by "
        "'patience') usually terminates training before this limit is "
        "reached.  A value of 200–500 is typical for datasets with "
        "10 000–50 000 samples.",
        "Training",
    ),
    ConfigParameter(
        "batch_size",
        "Mini-batch size for stochastic gradient descent.  Larger batches "
        "reduce variance in the gradient estimate and often converge in fewer "
        "epochs, but require more GPU memory.  256 is a good default for "
        "most GPU configurations.",
        "Training",
    ),
    ConfigParameter(
        "lr",
        "Initial learning rate for the Adam optimiser.  The "
        "ReduceLROnPlateau scheduler halves this value automatically when "
        "validation loss stagnates.  1e-3 is a robust starting point for "
        "most architectures and dataset sizes.",
        "Training",
    ),
    ConfigParameter(
        "weight_decay",
        "L2 regularisation coefficient (Adam weight decay).  Controls "
        "over-fitting for large networks.  1e-5 is a conservative default "
        "that rarely hurts performance.  Increase to 1e-4 or 1e-3 if "
        "training loss converges but validation loss diverges.",
        "Training",
    ),
    ConfigParameter(
        "patience",
        "Early-stopping patience: number of epochs without a validation-loss "
        "improvement larger than 'min_delta' before training is halted and "
        "the best checkpoint is restored.  20–50 epochs is standard; "
        "increase for large datasets or low learning rates.",
        "Training",
    ),
    ConfigParameter(
        "min_delta",
        "Minimum absolute decrease in validation loss to count as an "
        "improvement for early stopping.  Setting this too small allows "
        "the patience counter to reset on noise fluctuations.  "
        "1e-5 is appropriate for masked MSE loss on log-scaled targets.",
        "Training",
    ),
    ConfigParameter(
        "val_frac",
        "Fraction of the training dataset held out for validation.  "
        "Must be in (0, 1).  0.1 (ten percent) is standard.  With very "
        "large datasets (> 100 000 samples) 0.05 is sufficient.",
        "Training",
    ),
    ConfigParameter(
        "grad_clip",
        "Gradient-norm clipping threshold.  Gradients with a norm exceeding "
        "this value are rescaled before the optimiser step, preventing "
        "exploding gradients in deep networks.  1.0 is a safe default.  "
        "Set to null/None to disable clipping.",
        "Training",
    ),
    ConfigParameter(
        "seed",
        "Random seed for the train / validation split.  Does not affect "
        "dataset generation (controlled by ForwardConfig.seed).  "
        "Set to null/None for a non-reproducible split.",
        "Training",
    ),
    # ── Checkpointing ────────────────────────────────────────────────────────
    ConfigParameter(
        "checkpoint_dir",
        "Directory where model checkpoints (.npz files) are saved.  "
        "Created automatically if it does not exist.  Set to null/None to "
        "skip automatic checkpoint saving (the fitted inverter is still "
        "returned in memory).",
        "Checkpointing",
    ),
    ConfigParameter(
        "checkpoint_name",
        "Base file name (without extension) for checkpoint files.  "
        "The full path is checkpoint_dir/checkpoint_name.npz.  "
        "For example 'mt1d_resnet_5l' produces "
        "'checkpoints/mt1d_resnet_5l.npz'.",
        "Checkpointing",
    ),
    ConfigParameter(
        "save_best",
        "Automatically save the best-validation-loss checkpoint to "
        "checkpoint_path() after training completes.  When False the "
        "caller is responsible for calling inv.save().",
        "Checkpointing",
    ),
    # ── Output ───────────────────────────────────────────────────────────────
    ConfigParameter(
        "verbose",
        "Print per-epoch training summaries and a final report after "
        "training completes.  Set to False for batch scripts.",
        "Output",
    ),
]


# ── dataclass ────────────────────────────────────────────────────────────────

@dataclass
class InversionConfig:
    """Collect settings that define a 1-D AI-based EM inversion run.

    ``InversionConfig`` is the configuration object for
    :class:`~pycsamt.ai.inversion.inv1d.EMInverter1D`.  It covers four
    concern areas: network architecture, training hyperparameters,
    regularisation, and checkpoint management.

    The recommended workflow:

    1. Generate a template with :meth:`write_template`.
    2. Edit the values in the generated file.
    3. Load the edited file with :meth:`from_file`.
    4. Optionally call :meth:`validate` to catch range errors.
    5. Call :meth:`to_inverter` to instantiate a ready-to-fit
       :class:`~pycsamt.ai.inversion.inv1d.EMInverter1D`.
    6. Pass :meth:`to_fit_kwargs` to ``inv.fit(dataset, **cfg.to_fit_kwargs())``.

    Parameters
    ----------
    arch : {'resnet', 'cnn1d', 'fcn'}
        Network architecture.
    n_layers : int
        Number of earth layers (including halfspace).
    solver : {'mt1d', 'csamt1d', 'tem1d'}
        Forward solver this inverter targets.
    device : str or None
        Compute device; ``None`` auto-detects (CUDA > MPS > CPU).
    include_phase : bool
        Include impedance phase in the input feature vector.
    log_thickness : bool
        Apply log10 to thickness targets during training.
    augment_noise : float
        On-the-fly per-epoch noise augmentation level.
    epochs : int
        Maximum training epochs.
    batch_size : int
        Mini-batch size.
    lr : float
        Initial Adam learning rate.
    weight_decay : float
        Adam L2 regularisation coefficient.
    patience : int
        Early-stopping patience (epochs without improvement).
    min_delta : float
        Minimum validation-loss decrease to count as an improvement.
    val_frac : float
        Fraction of data used for validation.
    grad_clip : float or None
        Gradient-norm clipping threshold; ``None`` disables clipping.
    seed : int or None
        Random seed for train/val split.
    checkpoint_dir : str or None
        Directory for checkpoint files; ``None`` disables auto-saving.
    checkpoint_name : str
        Base file name for checkpoints (without extension).
    save_best : bool
        Auto-save the best checkpoint after training.
    verbose : bool
        Print training progress.

    Examples
    --------
    Default configuration (ResNet, 5 layers, MT1D)::

        >>> cfg = InversionConfig()
        >>> cfg.arch
        'resnet'

    Deep ResNet for a crystalline-crust survey::

        >>> cfg = InversionConfig(
        ...     arch="resnet",
        ...     n_layers=6,
        ...     solver="mt1d",
        ...     epochs=300,
        ...     lr=5e-4,
        ...     seed=0,
        ... )

    Round-trip template::

        >>> path = InversionConfig.write_template("inv_config.yml")
        >>> cfg = InversionConfig.from_file(path)
        >>> cfg.solver
        'mt1d'

    Snapshot a fitted inverter::

        >>> cfg = InversionConfig.from_inverter(inv)      # doctest: +SKIP
        >>> cfg.write_template("run_snapshot.py")         # doctest: +SKIP
    """

    # ── Architecture ─────────────────────────────────────────────────────────
    arch:           str            = "resnet"
    n_layers:       int            = 5
    solver:         str            = "mt1d"
    device:         Optional[str]  = None
    include_phase:  bool           = True
    log_thickness:  bool           = True
    augment_noise:  float          = 0.02

    # ── Training ─────────────────────────────────────────────────────────────
    epochs:        int             = 100
    batch_size:    int             = 256
    lr:            float           = 1e-3
    weight_decay:  float           = 1e-5
    patience:      int             = 20
    min_delta:     float           = 1e-5
    val_frac:      float           = 0.1
    grad_clip:     Optional[float] = 1.0
    seed:          Optional[int]   = None

    # ── Checkpointing ────────────────────────────────────────────────────────
    checkpoint_dir:  Optional[str] = "checkpoints"
    checkpoint_name: str           = "em_inverter"
    save_best:       bool          = True

    # ── Output ───────────────────────────────────────────────────────────────
    verbose: bool = True

    # ─────────────────────────────────────────────────────────────────────────
    # Validation
    # ─────────────────────────────────────────────────────────────────────────

    def validate(self) -> None:
        """Check parameter ranges and raise :class:`ValueError` on errors.

        Raises
        ------
        ValueError
            Descriptive message pointing to the offending parameter.
        """
        _VALID_ARCH = {"resnet", "cnn1d", "fcn"}
        if self.arch not in _VALID_ARCH:
            raise ValueError(
                f"arch must be one of {_VALID_ARCH!r}, got {self.arch!r}."
            )

        _VALID_SOLVERS = {"mt1d", "csamt1d", "tem1d"}
        if self.solver not in _VALID_SOLVERS:
            raise ValueError(
                f"solver must be one of {_VALID_SOLVERS!r}, "
                f"got {self.solver!r}."
            )

        if self.n_layers < 2:
            raise ValueError("n_layers must be at least 2 (1 layer + halfspace).")

        if self.device is not None:
            _VALID_DEVICES = {"cpu", "cuda", "mps"}
            if self.device not in _VALID_DEVICES and not self.device.startswith("cuda:"):
                raise ValueError(
                    f"device must be one of {_VALID_DEVICES!r} or 'cuda:N', "
                    f"got {self.device!r}."
                )

        if self.epochs < 1:
            raise ValueError("epochs must be at least 1.")
        if self.batch_size < 1:
            raise ValueError("batch_size must be at least 1.")
        if self.lr <= 0.0:
            raise ValueError("lr must be strictly positive.")
        if self.weight_decay < 0.0:
            raise ValueError("weight_decay must be non-negative.")
        if self.patience < 1:
            raise ValueError("patience must be at least 1.")
        if self.min_delta < 0.0:
            raise ValueError("min_delta must be non-negative.")
        if not 0.0 < self.val_frac < 1.0:
            raise ValueError("val_frac must be in (0, 1).")
        if self.grad_clip is not None and self.grad_clip <= 0.0:
            raise ValueError("grad_clip must be strictly positive when set.")
        if self.augment_noise < 0.0:
            raise ValueError("augment_noise must be non-negative.")

    # ─────────────────────────────────────────────────────────────────────────
    # Assemblers
    # ─────────────────────────────────────────────────────────────────────────

    def to_inverter(self) -> "EMInverter1D":
        """Instantiate a :class:`~pycsamt.ai.inversion.inv1d.EMInverter1D`.

        Returns an untrained inverter configured according to the
        architecture and feature settings stored in this config.  Call
        ``inv.fit(dataset, **cfg.to_fit_kwargs())`` to train it.

        Returns
        -------
        EMInverter1D

        Examples
        --------
        >>> cfg = InversionConfig(arch="cnn1d", n_layers=4, epochs=50)
        >>> inv = cfg.to_inverter()
        >>> type(inv).__name__
        'EMInverter1D'
        """
        from .inv1d import EMInverter1D
        return EMInverter1D(
            arch=self.arch,
            n_layers=self.n_layers,
            solver=self.solver,
            device=self.device,
            log_thickness=self.log_thickness,
            include_phase=self.include_phase,
            augment_noise=self.augment_noise,
        )

    def to_fit_kwargs(self) -> Dict[str, Any]:
        """Assemble keyword arguments for :meth:`EMInverter1D.fit`.

        The returned dict is ready to be unpacked directly::

            inv = cfg.to_inverter()
            inv.fit(dataset, **cfg.to_fit_kwargs())

        Returns
        -------
        dict
            Keys: ``epochs``, ``batch_size``, ``lr``, ``patience``,
            ``val_frac``, ``grad_clip``, ``seed``, ``verbose``.

        Notes
        -----
        ``weight_decay`` and ``min_delta`` are ``EMTrainer`` parameters
        not currently exposed through ``EMInverter1D.fit``.  They are
        stored in ``InversionConfig`` for documentation and round-trip
        reproducibility but are not included in the returned dict.
        """
        return dict(
            epochs=self.epochs,
            batch_size=self.batch_size,
            lr=self.lr,
            patience=self.patience,
            val_frac=self.val_frac,
            grad_clip=self.grad_clip,
            seed=self.seed,
            verbose=self.verbose,
        )

    def checkpoint_path(self) -> Optional[Path]:
        """Return the full checkpoint file path, or ``None`` if disabled.

        Returns
        -------
        pathlib.Path or None
        """
        if self.checkpoint_dir is None:
            return None
        p = Path(self.checkpoint_dir).expanduser()
        p.mkdir(parents=True, exist_ok=True)
        return p / f"{self.checkpoint_name}.npz"

    # ─────────────────────────────────────────────────────────────────────────
    # Snapshot / from_inverter
    # ─────────────────────────────────────────────────────────────────────────

    @classmethod
    def from_inverter(cls, inv: "EMInverter1D") -> "InversionConfig":
        """Snapshot a fitted (or unfitted) inverter's architecture settings.

        Creates an ``InversionConfig`` whose architecture and feature
        fields match those of *inv*.  Training hyperparameters are reset
        to their defaults because the inverter does not record them after
        training.

        Use this to generate a reproducible record of a training run::

            cfg = InversionConfig.from_inverter(inv)
            cfg.write_template("run_snapshot.py")

        Parameters
        ----------
        inv : EMInverter1D
            Source inverter (fitted or unfitted).

        Returns
        -------
        InversionConfig
        """
        return cls(
            arch=inv.arch,
            n_layers=inv.n_layers,
            solver=inv.solver,
            device=inv.device,
            log_thickness=getattr(inv, "log_thickness", True),
            include_phase=getattr(inv, "include_phase", True),
            augment_noise=getattr(inv, "augment_noise", 0.02),
        )

    # ─────────────────────────────────────────────────────────────────────────
    # Config file I/O
    # ─────────────────────────────────────────────────────────────────────────

    def to_template(
        self,
        path: Union[str, Path] = "inversion_config.py",
        *,
        fmt: Optional[str] = None,
    ) -> Path:
        """Write this configuration to an annotated source-of-truth file.

        Parameters
        ----------
        path : path-like, default "inversion_config.py"
            Destination file.  The suffix selects the output format
            (``.py``, ``.json``, ``.yml``).
        fmt : {"py", "json", "yml", "yaml"}, optional
            Explicit format override.

        Returns
        -------
        pathlib.Path
        """
        return write_config_template(
            path,
            self,
            _INVERSION_CONFIG_SCHEMA,
            fmt=fmt,
            title="PyCSAMT AI inversion configuration",
        )

    @classmethod
    def write_template(
        cls,
        path: Union[str, Path] = "inversion_config.py",
        *,
        fmt: Optional[str] = None,
    ) -> Path:
        """Generate a documented source-of-truth configuration file.

        Creates a file with default parameter values and an inline comment
        for every parameter.  Edit the file, then load with
        :meth:`from_file`.

        Parameters
        ----------
        path : path-like, default "inversion_config.py"
            Destination file.
        fmt : {"py", "json", "yml", "yaml"}, optional
            Explicit format override.

        Returns
        -------
        pathlib.Path

        Examples
        --------
        >>> from pycsamt.ai.inversion.config import InversionConfig
        >>> path = InversionConfig.write_template("my_inv.yml")
        >>> path.suffix
        '.yml'
        """
        return cls().to_template(path, fmt=fmt)

    @classmethod
    def from_file(
        cls,
        path: Union[str, Path],
        *,
        strict: bool = True,
    ) -> "InversionConfig":
        """Load a configuration from a source-of-truth file.

        Parameters
        ----------
        path : path-like
            Python, JSON, YML, or YAML file generated by
            :meth:`write_template` or following the same structure.
        strict : bool, default True
            If ``True``, unknown keys raise :class:`ValueError`.
            If ``False``, unknown keys are silently ignored.

        Returns
        -------
        InversionConfig

        Examples
        --------
        >>> InversionConfig.write_template("inv_config.json")
        PosixPath('inv_config.json')
        >>> cfg = InversionConfig.from_file("inv_config.json")
        >>> cfg.arch
        'resnet'
        """
        values = read_config_file(path, cls, strict=strict)
        return cls(**values)

    #: Alias — matches the convention used by ModEmConfig and OccamConfig.
    read = from_file

    # ─────────────────────────────────────────────────────────────────────────
    # repr / summary
    # ─────────────────────────────────────────────────────────────────────────

    def summary(self) -> str:
        """Return a human-readable multi-line summary of the configuration."""
        phase_s = "yes" if self.include_phase else "no"
        log_th_s = "yes" if self.log_thickness else "no"
        clip_s = str(self.grad_clip) if self.grad_clip is not None else "off"
        ckpt = str(self.checkpoint_path()) if self.checkpoint_path() else "disabled"
        lines = [
            "InversionConfig",
            "  ── Architecture ──",
            f"  {'arch':<22s} = {self.arch!r}",
            f"  {'n_layers':<22s} = {self.n_layers}",
            f"  {'solver':<22s} = {self.solver!r}",
            f"  {'device':<22s} = {self.device!r}  (None → auto)",
            f"  {'include_phase':<22s} = {phase_s}",
            f"  {'log_thickness':<22s} = {log_th_s}",
            f"  {'augment_noise':<22s} = {self.augment_noise}",
            "  ── Training ──",
            f"  {'epochs':<22s} = {self.epochs}",
            f"  {'batch_size':<22s} = {self.batch_size}",
            f"  {'lr':<22s} = {self.lr}",
            f"  {'weight_decay':<22s} = {self.weight_decay}",
            f"  {'patience':<22s} = {self.patience}  (min_delta={self.min_delta})",
            f"  {'val_frac':<22s} = {self.val_frac}",
            f"  {'grad_clip':<22s} = {clip_s}",
            f"  {'seed':<22s} = {self.seed!r}",
            "  ── Checkpointing ──",
            f"  {'checkpoint':<22s} = {ckpt}",
            f"  {'save_best':<22s} = {self.save_best}",
        ]
        return "\n".join(lines)

    def __repr__(self) -> str:
        return self.summary()
