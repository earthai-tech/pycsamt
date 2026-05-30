# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""InputBuilder — assemble a complete ModEM input set from EDI data.

Given a collection of EDI sites (or objects with compatible attributes)
and a configuration, :class:`InputBuilder` writes all input files required
to launch a ModEM inversion into a target directory:

* ``data.dat``      — observed data
* ``m0.ws`` / ``m0.rho`` — half-space starting model (3D / 2D)
* ``covariance.cov`` — covariance file  (3D only)
* ``control.inv``   — inversion control parameters

Quick start
-----------
>>> from pycsamt.models.modem.builder import InputBuilder
>>> from pycsamt.models.modem.config  import ModEmConfig
>>> cfg     = ModEmConfig(mode="3d", initial_rho=100.0)
>>> builder = InputBuilder(config=cfg)
>>> files   = builder.build(edi_list, workdir="./modem_run")
>>> print(files)          # dict of written paths
"""

from __future__ import annotations

from pathlib import Path
from typing  import Dict, List, Optional, Union

from .base       import ModEmBase
from .config     import ModEmConfig
from .data       import ModEmData
from .model2d    import ModEmModel2D
from .model3d    import ModEmModel3D
from .covariance import ModEmCovariance
from .control    import ModEmControl

PathLike = Union[str, Path]

__all__ = ["InputBuilder"]


class InputBuilder(ModEmBase):
    """Build and write a complete ModEM input set.

    Parameters
    ----------
    config : ModEmConfig, optional
        Run configuration.  Defaults to 3-D with ``Full_Impedance``.

    Attributes
    ----------
    config : ModEmConfig
    data : ModEmData or None
        Populated after :meth:`build`.
    model : ModEmModel2D or ModEmModel3D or None
    covariance : ModEmCovariance or None  (3-D only)
    control : ModEmControl or None
    """

    def __init__(
        self,
        config: Optional[ModEmConfig] = None,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.config:     ModEmConfig            = config or ModEmConfig()
        self.data:       Optional[ModEmData]    = None
        self.model:      Optional[Union[ModEmModel2D, ModEmModel3D]] = None
        self.covariance: Optional[ModEmCovariance] = None
        self.control:    Optional[ModEmControl] = None

    # ------------------------------------------------------------------
    # Main entry point
    # ------------------------------------------------------------------

    def build(
        self,
        source,
        workdir: PathLike = ".",
        *,
        data_filename:  str = "data.dat",
        model_filename: Optional[str] = None,
        cov_filename:   str = "covariance.cov",
        ctrl_filename:  str = "control.inv",
    ) -> Dict[str, Path]:
        """Write all ModEM input files to *workdir*.

        Parameters
        ----------
        source : iterable of site-like objects
            Each item must expose the same interface accepted by
            :meth:`ModEmData.from_edi`:
            ``name``, ``coords``, ``freq``, ``z``, ``z_err``.
        workdir : path-like
            Output directory (created if absent).
        data_filename, model_filename, cov_filename, ctrl_filename : str
            Override default output file names.

        Returns
        -------
        dict
            Mapping ``{'data', 'model', 'covariance', 'control'}`` → resolved
            :class:`~pathlib.Path` (``'covariance'`` absent for 2-D runs).
        """
        cfg = self.config
        wd  = Path(workdir)
        wd.mkdir(parents=True, exist_ok=True)

        is_3d = cfg.is_3d

        # -- default model filename based on mode ----------------------
        if model_filename is None:
            model_filename = "m0.ws" if is_3d else "m0.rho"

        # 1. Data ---------------------------------------------------------
        self.data = ModEmData.from_edi(list(source), config=cfg)
        data_path = self.data.write(wd / data_filename)

        # 2. Starting model -----------------------------------------------
        if is_3d:
            self.model = ModEmModel3D.halfspace(self.data, config=cfg)
        else:
            self.model = ModEmModel2D.halfspace(self.data, config=cfg)
        model_path = self.model.write(wd / model_filename)

        # 3. Covariance (3-D only) ----------------------------------------
        cov_path: Optional[Path] = None
        if is_3d:
            self.covariance = ModEmCovariance.from_model(self.model, config=cfg)
            cov_path = self.covariance.write(wd / cov_filename)

        # 4. Control -------------------------------------------------------
        self.control = ModEmControl.from_config(cfg)
        ctrl_path    = self.control.write(wd / ctrl_filename)

        files: Dict[str, Path] = {
            "data":    data_path,
            "model":   model_path,
            "control": ctrl_path,
        }
        if cov_path is not None:
            files["covariance"] = cov_path

        if self.verbose:
            self.logger.info(
                "InputBuilder.build: wrote %d files to %s",
                len(files), wd,
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
        model_filename: Optional[str] = None,
        cov_filename:   str = "covariance.cov",
        ctrl_filename:  str = "control.inv",
    ) -> Dict[str, Path]:
        """Write model, covariance, and control files from an existing
        :class:`ModEmData` object (data file must already be written).

        Returns
        -------
        dict
            ``{'model', 'control'}`` and optionally ``'covariance'``.
        """
        cfg = self.config
        wd  = Path(workdir)
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

        cov_path: Optional[Path] = None
        if is_3d:
            self.covariance = ModEmCovariance.from_model(self.model, config=cfg)
            cov_path = self.covariance.write(wd / cov_filename)

        self.control  = ModEmControl.from_config(cfg)
        ctrl_path     = self.control.write(wd / ctrl_filename)

        files: Dict[str, Path] = {"model": model_path, "control": ctrl_path}
        if cov_path is not None:
            files["covariance"] = cov_path
        return files
