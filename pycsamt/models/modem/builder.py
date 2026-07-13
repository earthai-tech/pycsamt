# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Build complete ModEM input sets from survey data.

This module provides :class:`InputBuilder`, the high-level
factory that converts EDI-like survey objects or an existing
:class:`~pycsamt.models.modem.data.ModEmData` instance into a
self-contained ModEM working directory.
"""

from __future__ import annotations

from pathlib import Path
from typing import Union

from .base import ModEmBase
from .config import ModEmConfig
from .control import ModEmControl
from .covariance import ModEmCovariance
from .data import ModEmData
from .doc import _modem_param_docs as _params
from .model2d import ModEmModel2D
from .model3d import ModEmModel3D

PathLike = Union[str, Path]

__all__ = ["InputBuilder"]


class InputBuilder(ModEmBase):
    def __init__(
        self,
        config: ModEmConfig | None = None,
        **kwargs,
    ):
        """Initialize an input builder with a run configuration."""
        super().__init__(**kwargs)
        self.config: ModEmConfig = config or ModEmConfig()
        self.data: ModEmData | None = None
        self.model: ModEmModel2D | ModEmModel3D | None = None
        self.covariance: ModEmCovariance | None = None
        self.control: ModEmControl | None = None

    # ------------------------------------------------------------------
    # Main entry point
    # ------------------------------------------------------------------

    def build(
        self,
        source,
        workdir: PathLike = ".",
        *,
        data_filename: str = "data.dat",
        model_filename: str | None = None,
        cov_filename: str = "covariance.cov",
        ctrl_filename: str = "control.inv",
    ) -> dict[str, Path]:
        """Write the complete ModEM input set to ``workdir``.

        The method is the normal entry point for preparing a
        new ModEM inversion from EDI-like survey data. It first
        converts ``source`` to :class:`ModEmData`, then builds
        a uniform half-space starting model with the geometry
        defined by ``self.config``. For 3-D runs it also derives
        a covariance file from the model grid. Finally, it
        writes the inversion-control file.

        The generated file set is:

        * observed data, usually ``data.dat``;
        * starting model, ``m0.ws`` in 3-D or ``m0.rho`` in 2-D;
        * covariance file for 3-D runs only;
        * inversion-control file, usually ``control.inv``.

        Parameters
        ----------
        source : iterable of site-like objects
            Survey source accepted by :meth:`ModEmData.from_edi`.
            Each item should provide a station name, coordinates,
            frequency samples, a complex impedance tensor, and
            impedance errors. The builder consumes the iterable
            immediately, so generators are supported but cannot be
            reused after the call.
        workdir : path-like, default "."
            Output directory that will contain the ModEM input
            files. The directory is created if it does not exist.
            Relative paths are resolved from the current Python
            working directory.
        data_filename : str, default "data.dat"
            Name of the observed-data file written in ``workdir``.
            Use this to keep several component selections or
            frequency bands in the same parent folder.
        model_filename : str, optional
            Name of the starting-model file. If omitted,
            ``"m0.ws"`` is used for 3-D runs and ``"m0.rho"``
            is used for 2-D runs.
        cov_filename : str, default "covariance.cov"
            Name of the covariance file written for 3-D runs.
            The value is ignored for 2-D runs because the current
            2-D builder path does not write a covariance file.
        ctrl_filename : str, default "control.inv"
            Name of the inversion-control file written in
            ``workdir``.

        Returns
        -------
        dict
            Mapping from output role to resolved
            :class:`pathlib.Path`. The mapping always contains
            ``"data"``, ``"model"``, and ``"control"``. It also
            contains ``"covariance"`` for 3-D runs.

        Raises
        ------
        ValueError
            Raised by :meth:`ModEmData.from_edi` when ``source``
            is empty or cannot provide the required station and
            impedance information.

        Examples
        --------
        Build a standard 3-D input set from EDI-like stations:

        >>> from pycsamt.models.modem.builder import InputBuilder
        >>> from pycsamt.models.modem.config import ModEmConfig
        >>> cfg = ModEmConfig(mode="3d", initial_rho=100.0)
        >>> builder = InputBuilder(config=cfg)
        >>> files = builder.build(sites, workdir="modem_3d")
        >>> sorted(files)
        ['control', 'covariance', 'data', 'model']

        Build a 2-D TE inversion input set with custom names:

        >>> cfg = ModEmConfig(
        ...     mode="2d",
        ...     component_type="TE_Impedance",
        ... )
        >>> builder = InputBuilder(config=cfg)
        >>> files = builder.build(
        ...     sites,
        ...     workdir="modem_2d_te",
        ...     data_filename="te_data.dat",
        ...     model_filename="te_start.rho",
        ...     ctrl_filename="te_control.inv",
        ... )
        >>> "covariance" in files
        False
        """
        cfg = self.config
        wd = Path(workdir)
        wd.mkdir(parents=True, exist_ok=True)

        is_3d = cfg.is_3d

        # -- default model filename based on mode ----------------------
        if model_filename is None:
            model_filename = "m0.ws" if is_3d else "m0.rho"

        # 1. Data ------------------------------------------------------
        self.data = ModEmData.from_edi(list(source), config=cfg)
        data_path = self.data.write(wd / data_filename)

        # 2. Starting model --------------------------------------------
        if is_3d:
            self.model = ModEmModel3D.halfspace(self.data, config=cfg)
        else:
            self.model = ModEmModel2D.halfspace(self.data, config=cfg)
        model_path = self.model.write(wd / model_filename)

        # 3. Covariance (3-D only) -------------------------------------
        cov_path: Path | None = None
        if is_3d:
            self.covariance = ModEmCovariance.from_model(
                self.model,
                config=cfg,
            )
            cov_path = self.covariance.write(wd / cov_filename)

        # 4. Control ---------------------------------------------------
        self.control = ModEmControl.from_config(cfg)
        ctrl_path = self.control.write(wd / ctrl_filename)

        files: dict[str, Path] = {
            "data": data_path,
            "model": model_path,
            "control": ctrl_path,
        }
        if cov_path is not None:
            files["covariance"] = cov_path

        if self.verbose:
            self.logger.info(
                "InputBuilder.build: wrote %d files to %s",
                len(files),
                wd,
            )

        return files

    # ------------------------------------------------------------------
    # Convenience: build from a ModEmData that is already assembled
    # ------------------------------------------------------------------

    def build_from_data(
        self,
        data: ModEmData,
        workdir: PathLike = ".",
        *,
        model_filename: str | None = None,
        cov_filename: str = "covariance.cov",
        ctrl_filename: str = "control.inv",
    ) -> dict[str, Path]:
        """Write run files from an existing :class:`ModEmData`.

        This method is useful when the observed-data object has
        already been assembled, filtered, edited, or written by a
        caller. It uses the supplied data to build the starting
        model, derives a covariance file for 3-D workflows, and
        writes the inversion-control file. It does not write the
        data file itself, so callers should write ``data`` to the
        run directory when the ModEM executable will need it.

        Parameters
        ----------
        data : ModEmData
            Populated data object used to derive station extents,
            periods, and model geometry. The object is retained as
            ``self.data`` after the method returns.
        workdir : path-like, default "."
            Output directory for the generated starting model,
            optional covariance file, and control file. The
            directory is created if it does not exist.
        model_filename : str, optional
            Name of the starting-model file. If omitted,
            ``"m0.ws"`` is used for 3-D runs and ``"m0.rho"``
            is used for 2-D runs.
        cov_filename : str, default "covariance.cov"
            Name of the covariance file written for 3-D runs.
            Ignored for 2-D workflows.
        ctrl_filename : str, default "control.inv"
            Name of the inversion-control file written in
            ``workdir``.

        Returns
        -------
        dict
            Mapping from output role to resolved
            :class:`pathlib.Path`. The mapping contains
            ``"model"`` and ``"control"`` for all modes, and
            also ``"covariance"`` for 3-D runs.

        Examples
        --------
        Reuse a data object that was prepared elsewhere:

        >>> from pycsamt.models.modem.builder import InputBuilder
        >>> from pycsamt.models.modem.data import ModEmData
        >>> data = ModEmData.from_edi(sites, config=cfg)
        >>> builder = InputBuilder(config=cfg)
        >>> files = builder.build_from_data(data, workdir="run")
        >>> sorted(files)
        ['control', 'covariance', 'model']
        """
        cfg = self.config
        wd = Path(workdir)
        wd.mkdir(parents=True, exist_ok=True)

        is_3d = cfg.is_3d
        if model_filename is None:
            model_filename = "m0.ws" if is_3d else "m0.rho"

        self.data = data
        if is_3d:
            self.model = ModEmModel3D.halfspace(data, config=cfg)
        else:
            self.model = ModEmModel2D.halfspace(data, config=cfg)
        model_path = self.model.write(wd / model_filename)

        cov_path: Path | None = None
        if is_3d:
            self.covariance = ModEmCovariance.from_model(
                self.model,
                config=cfg,
            )
            cov_path = self.covariance.write(wd / cov_filename)

        self.control = ModEmControl.from_config(cfg)
        ctrl_path = self.control.write(wd / ctrl_filename)

        files: dict[str, Path] = {
            "model": model_path,
            "control": ctrl_path,
        }
        if cov_path is not None:
            files["covariance"] = cov_path
        return files


InputBuilder.__doc__ = rf"""
Build and write a complete ModEM input set.

``InputBuilder`` is the main preparation object for the ModEM
v2 workflow. It turns a survey source into the files required
by the ModEM executable: an observed-data file, a half-space
starting model, a covariance file for 3-D runs, and an
inversion-control file. The builder does not launch ModEM.
Instead, it creates a consistent working directory that can be
passed to :class:`~pycsamt.models.modem.runner.ModEmRunner`.

The build sequence is deterministic:

1. Convert the survey source to :class:`ModEmData`.
2. Build a 2-D or 3-D half-space starting model.
3. Derive a 3-D covariance file when ``config.mode == "3d"``.
4. Write the ModEM inversion-control file.

For a uniform starting resistivity :math:`\rho_0`, model
writers store logarithmic resistivity values, commonly

.. math::

   m_0 = \ln(\rho_0).

These files are then used by ModEM to solve a regularized
nonlinear inverse problem for the subsurface resistivity
distribution [1]_, [2]_.

Parameters
----------
{_params.common.config}
{_params.common.verbose}
{_params.common.logger}

Attributes
----------
config : ModEmConfig
    Mutable run configuration used by all build steps. It
    selects dimensionality, component type, error floors,
    model geometry, covariance smoothing, control values, and
    default executable/file names.
data : ModEmData or None
    Data object generated by :meth:`build` or supplied through
    :meth:`build_from_data`. It stores station coordinates,
    periods, component rows, complex values, and errors.
model : ModEmModel2D or ModEmModel3D or None
    Starting model generated from ``data`` and ``config``.
    The concrete class depends on ``config.mode``.
covariance : ModEmCovariance or None
    Covariance object generated only for 3-D runs. It stores
    smoothing controls and active-cell masks derived from the
    starting model.
control : ModEmControl or None
    Inversion-control object generated from ``config``.

Build Parameters
----------------
{_params.data.source}
{_params.common.workdir}
{_params.builder.data_filename}
{_params.builder.model_filename}
{_params.builder.cov_filename}
{_params.builder.ctrl_filename}

Notes
-----
``InputBuilder`` writes a half-space model by design. This is
the standard starting point for many ModEM inversions and gives
the solver a stable initial model. Users who need a custom
starting model can write it separately and pass that file to
the runner while still using this builder for data, covariance,
and control-file generation.

The 2-D builder path writes data, model, and control files.
The 3-D path additionally writes a covariance file because
ModEM 3-D inversions use explicit smoothing and mask
definitions.

See Also
--------
ModEmData.from_edi
    Convert EDI-like station objects into ModEM data rows.
ModEmModel2D.halfspace
    Build a 2-D half-space model from station geometry.
ModEmModel3D.halfspace
    Build a 3-D half-space model from station geometry.
ModEmCovariance.from_model
    Derive 3-D smoothing and active-cell masks from a model.
ModEmControl.from_config
    Build inversion controls from :class:`ModEmConfig`.
ModEmRunner
    Execute ModEM with the generated input files.

Examples
--------
Build the default 3-D ModEM files:

>>> from pycsamt.models.modem.builder import InputBuilder
>>> from pycsamt.models.modem.config import ModEmConfig
>>> cfg = ModEmConfig(mode="3d", initial_rho=100.0)
>>> builder = InputBuilder(config=cfg)
>>> files = builder.build(sites, workdir="modem_run")
>>> sorted(files)
['control', 'covariance', 'data', 'model']

Build a 2-D TE input set:

>>> cfg = ModEmConfig(
...     mode="2d",
...     component_type="TE_Impedance",
...     nz_2d=60,
... )
>>> builder = InputBuilder(config=cfg)
>>> files = builder.build(sites, workdir="modem_te")
>>> "covariance" in files
False

Reuse an already assembled data object:

>>> from pycsamt.models.modem.data import ModEmData
>>> data = ModEmData.from_edi(sites, config=cfg)
>>> data.write("modem_te/data.dat")
>>> files = builder.build_from_data(data, workdir="modem_te")
>>> sorted(files)
['control', 'model']

Use custom file names for a test run:

>>> cfg = ModEmConfig(mode="3d", output_stem="trial")
>>> builder = InputBuilder(config=cfg)
>>> files = builder.build(
...     sites,
...     workdir="trial_run",
...     data_filename="obs_trial.dat",
...     model_filename="start_trial.ws",
...     cov_filename="smooth_trial.cov",
...     ctrl_filename="trial.inv",
... )
>>> files["control"].name
'trial.inv'

References
----------
.. [1] Egbert, G. D., and Kelbert, A., "Computational
   recipes for electromagnetic inverse problems", Geophysical
   Journal International, 189(1), 251-267, 2012,
   doi:10.1111/j.1365-246X.2011.05347.x.
.. [2] Kelbert, A., Meqbel, N., Egbert, G. D., and Tandon,
   K., "ModEM: A modular system for inversion of
   electromagnetic geophysical data", Computers and
   Geosciences, 66, 40-53, 2014,
   doi:10.1016/j.cageo.2014.01.010.
"""
