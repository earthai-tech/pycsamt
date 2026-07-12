# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Build complete MARE2DEM input sets from survey data.

:class:`InputBuilder` is the high-level factory that:

1. Writes the ``.emdata`` data file via
   :func:`~pycsamt.models.mare2dem.survey.make_data_file`.
2. Writes a starting ``.resistivity`` model file.
3. Writes the ``.settings`` parallel-decomposition file.

The full FEM mesh generation (Triangle ``.poly`` file) is not yet
implemented — it will be added when ``Mamba2D.m`` is fully ported.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Union

from .base import Mare2DEMBase
from .config import Mare2DEMConfig
from .doc import _mare2dem_param_docs as _params
from .iotools.emdata import EMDataFile, read_emdata
from .iotools.resistivity import (
    ResistivityFile,
    write_resistivity,
)
from .iotools.settings import SettingsFile, write_settings
from .survey import (
    CSEMSurveyConfig,
    MTSurveyConfig,
    make_data_file,
)

PathLike = Union[str, Path]

__all__ = ["InputBuilder"]


class InputBuilder(Mare2DEMBase):
    """Prepare a MARE2DEM working directory from survey parameters.

    Parameters
    ----------
    config : Mare2DEMConfig, optional
        Configuration controlling initial resistivity, iteration
        limits, target misfit, and output file names.
    **kwargs :
        Forwarded to :class:`Mare2DEMBase`.
    """

    def __init__(
        self,
        config: Mare2DEMConfig | None = None,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.config: Mare2DEMConfig = config or Mare2DEMConfig()
        self._em: EMDataFile | None = None

    # ------------------------------------------------------------------
    # Settings file
    # ------------------------------------------------------------------

    def write_settings(
        self,
        workdir: PathLike = ".",
        *,
        filename: str | None = None,
        **sf_kwargs,
    ) -> Path:
        """Write the MARE2DEM ``.settings`` file.

        Parameters
        ----------
        workdir : path-like, default "."
            Target directory.
        filename : str, optional
            Override the settings filename.
        **sf_kwargs :
            Extra keyword arguments forwarded to :class:`SettingsFile`.

        Returns
        -------
        pathlib.Path
            Path of the written file.
        """
        cfg = self.config
        dest = Path(workdir)
        dest.mkdir(parents=True, exist_ok=True)
        fname = filename or cfg.settings_file
        sf = SettingsFile(**sf_kwargs)
        path = write_settings(sf, dest / fname)
        if self.verbose:
            self.logger.info("InputBuilder: wrote settings to %s", path)
        return path

    # ------------------------------------------------------------------
    # Resistivity model file
    # ------------------------------------------------------------------

    def write_resistivity(
        self,
        workdir: PathLike = ".",
        *,
        filename: str | None = None,
        poly_file: str | None = None,
    ) -> Path:
        """Write a homogeneous half-space ``.resistivity`` file.

        Parameters
        ----------
        workdir : path-like, default "."
            Target directory.
        filename : str, optional
            Override the resistivity filename.
        poly_file : str, optional
            Poly file reference written into the model header.

        Returns
        -------
        pathlib.Path
            Path of the written file.
        """
        cfg = self.config
        dest = Path(workdir)
        dest.mkdir(parents=True, exist_ok=True)
        fname = filename or cfg.resistivity_file
        log10_rho = math.log10(max(cfg.initial_rho, 1e-10))
        import numpy as np
        rf = ResistivityFile(
            resistivity_file=str(dest / fname),
            poly_file=poly_file or fname.replace(".resistivity", ".poly"),
            data_file=cfg.data_file,
            settings_file=cfg.settings_file,
            target_misfit=cfg.target_rms,
            max_iterations=cfg.max_iterations,
        )
        rf.resistivity = np.array([[log10_rho]])
        rf.free_parameter = np.array([[1]])
        rf.bounds = np.array([[-2.0, 5.0]])
        rf.prejudice = np.zeros((1, 2))
        path = write_resistivity(rf, dest / fname)
        if self.verbose:
            self.logger.info("InputBuilder: wrote resistivity to %s", path)
        return path

    # ------------------------------------------------------------------
    # Main entry point
    # ------------------------------------------------------------------

    def build(
        self,
        source,
        workdir: PathLike = ".",
        *,
        topo=0.0,
        mt: MTSurveyConfig | None = None,
        csem: CSEMSurveyConfig | None = None,
        data_filename: str | None = None,
        model_filename: str | None = None,
        settings_filename: str | None = None,
    ) -> dict[str, Path]:
        """Write a MARE2DEM input set to *workdir*.

        Parameters
        ----------
        source : path-like, EMDataFile, or None
            Existing ``.emdata`` file (copied to workdir) **or**
            ``None`` (generate from *mt*/*csem* config objects).
        workdir : path-like, default "."
            Target directory.
        topo : float or array-like
            Topography for receiver/transmitter placement (used only
            when *source* is ``None``).
        mt : MTSurveyConfig, optional
            MT survey config (used when *source* is ``None``).
        csem : CSEMSurveyConfig, optional
            CSEM survey config (used when *source* is ``None``).
        data_filename : str, optional
            Override data file name.
        model_filename : str, optional
            Override resistivity model file name.
        settings_filename : str, optional
            Override settings file name.

        Returns
        -------
        dict[str, pathlib.Path]
            Keys: ``"data"``, ``"model"``, ``"settings"``.
        """
        import shutil

        cfg = self.config
        dest = Path(workdir)
        dest.mkdir(parents=True, exist_ok=True)
        result: dict[str, Path] = {}

        _data_fname = data_filename or cfg.data_file
        _model_fname = model_filename or cfg.resistivity_file
        _settings_fname = settings_filename or cfg.settings_file

        # ---- data file ----
        if isinstance(source, (str, Path)) and Path(source).exists():
            data_path = dest / _data_fname
            shutil.copy2(str(source), str(data_path))
            self._em = read_emdata(data_path)
        elif isinstance(source, EMDataFile):
            from .iotools.emdata import write_emdata
            data_path = write_emdata(source, dest / _data_fname)
            self._em = source
        elif mt is not None or csem is not None:
            data_path = dest / _data_fname
            self._em = make_data_file(
                data_path, topo, mt=mt, csem=csem
            )
        else:
            data_path = dest / _data_fname
            self.logger.warning(
                "InputBuilder.build: no data source provided — data file not written."
            )
        result["data"] = dest / _data_fname

        # ---- resistivity model ----
        model_path = self.write_resistivity(dest, filename=_model_fname)
        result["model"] = model_path

        # ---- settings ----
        settings_path = self.write_settings(dest, filename=_settings_fname)
        result["settings"] = settings_path

        if self.verbose:
            self.logger.info(
                "InputBuilder.build: wrote %d files to %s",
                len(result),
                dest,
            )
        return result


InputBuilder.__doc__ = rf"""
Prepare a MARE2DEM working directory from survey parameters.

``InputBuilder`` produces the three required input files for a MARE2DEM
inversion run:

* the ``.emdata`` observed-data file;
* the starting ``.resistivity`` model;
* the ``.settings`` parallel-decomposition control file.

Parameters
----------
config : Mare2DEMConfig, optional
    Configuration for inversion parameters and file names.
{_params.common.verbose}
{_params.common.logger}

Examples
--------
Build from an existing ``.emdata`` file:

>>> from pycsamt.models.mare2dem import Mare2DEMConfig, InputBuilder
>>> cfg = Mare2DEMConfig(initial_rho=1.0, max_iterations=100)
>>> builder = InputBuilder(config=cfg)
>>> files = builder.build("survey.emdata", workdir="./run")

Build from MTSurveyConfig:

>>> import numpy as np
>>> from pycsamt.models.mare2dem.survey import MTSurveyConfig
>>> mt = MTSurveyConfig(
...     frequencies=np.logspace(-3, 3, 10),
...     rx_y=np.linspace(-5000, 5000, 20),
...     rx_type="marine", lTE=True, lTM=True,
... )
>>> files = builder.build(None, workdir="./run", mt=mt)

See Also
--------
Mare2DEMRunner
    Launch MARE2DEM on the written input files.
make_data_file
    Low-level data file generator.
"""
