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

__all__ = ["InputBuilder", "build_mt_inputs", "mt_core_grid"]


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
        import numpy as np

        rf = ResistivityFile(
            resistivity_file=str(dest / fname),
            poly_file=poly_file or fname.replace(".resistivity", ".poly"),
            data_file=cfg.data_file,
            settings_file=cfg.settings_file,
            target_misfit=cfg.target_rms,
            max_iterations=cfg.max_iterations,
        )
        # The Rho column is linear ohm-m (it used to receive log10 rho, so
        # 100 ohm-m became 2 ohm-m); 0 0 bounds defer to the global ones.
        rf.resistivity = np.array([[max(cfg.initial_rho, 1e-10)]])
        rf.free_parameter = np.array([[1]])
        rf.bounds = np.zeros((1, 2))
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
            self._em = make_data_file(data_path, topo, mt=mt, csem=csem)
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


# ---------------------------------------------------------------------------
# EDI → complete MT inversion input set (data + mesh + model + settings)
# ---------------------------------------------------------------------------


def mt_core_grid(
    rx_y,
    frequencies,
    *,
    initial_rho: float = 100.0,
    cell_y: float | None = None,
    n_extra_y: int = 2,
    z_first: float | None = None,
    depth_max: float | None = None,
    growth: float = 1.15,
):
    """Return the core-grid cell edges under a line of MT receivers.

    Parameters
    ----------
    rx_y : array-like
        Along-profile receiver positions in metres.
    frequencies : array-like
        Data frequencies in Hz; they set the depth range through skin
        depths ``δ = 503·sqrt(ρ / f)`` of the starting resistivity.
    initial_rho : float, default 100
        Starting resistivity in ohm metres (skin-depth estimate only).
    cell_y : float, optional
        Lateral cell width. Default: half the median receiver spacing.
    n_extra_y : int, default 2
        Core cells added beyond the first and last receiver.
    z_first : float, optional
        First cell thickness. Default: a quarter of the shallowest skin
        depth, at least 5 m.
    depth_max : float, optional
        Bottom of the core grid. Default: 1.5 × the deepest skin depth.
    growth : float, default 1.15
        Geometric thickness growth factor with depth.

    Returns
    -------
    y_edges, z_edges : numpy.ndarray
        Cell edges in metres (``z`` positive down, starting at 0).
    """
    import numpy as np

    y = np.unique(np.asarray(rx_y, dtype=float))
    f = np.asarray(frequencies, dtype=float)
    f = f[np.isfinite(f) & (f > 0)]
    if y.size == 0 or f.size == 0:
        raise ValueError("need at least one receiver and one frequency")
    rho = max(float(initial_rho), 1e-3)
    if cell_y is None:
        spacing = float(np.median(np.diff(y))) if y.size > 1 else 500.0
        cell_y = max(spacing / 2.0, 10.0)
    y0 = y.min() - n_extra_y * cell_y
    y1 = y.max() + n_extra_y * cell_y
    n_y = max(int(math.ceil((y1 - y0) / cell_y)), 1)
    y_edges = y0 + cell_y * np.arange(n_y + 1)

    delta_min = 503.0 * math.sqrt(rho / f.max())
    delta_max = 503.0 * math.sqrt(rho / f.min())
    dz = max(z_first if z_first is not None else delta_min / 4.0, 5.0)
    bottom = depth_max if depth_max is not None else 1.5 * delta_max
    z_edges = [0.0]
    while z_edges[-1] < bottom:
        z_edges.append(z_edges[-1] + dz)
        dz *= max(float(growth), 1.0)
    return y_edges, np.asarray(z_edges)


def build_mt_inputs(
    sites,
    workdir: PathLike,
    config: Mare2DEMConfig | None = None,
    *,
    stem: str = "mare2dem",
    error_floor: float = 0.05,
    output_modes: str = "all",
    cell_y: float | None = None,
    z_first: float | None = None,
    depth_max: float | None = None,
    growth: float = 1.15,
    padding: float | None = None,
    **edi_kwargs,
) -> dict:
    """Write a complete, runnable MARE2DEM MT inversion from EDI sites.

    :meth:`InputBuilder.write_resistivity` writes a one-region stub with
    no ``.poly`` mesh, which MARE2DEM cannot run.  This builds the whole
    set: the ``.emdata`` from the sites
    (:func:`~pycsamt.models.mare2dem.edi.make_mt_data_from_edi`), then a
    regular core grid under the receivers (:func:`mt_core_grid`) turned
    into ``.poly`` + ``.resistivity`` + ``.settings`` by
    :func:`~pycsamt.models.mare2dem.grid_to_m2d.grid_to_mare2dem` (fixed
    air, fixed padding, one free parameter per core cell).

    Parameters
    ----------
    sites : Sites or path-like
        Anything accepted by ``make_mt_data_from_edi``.
    workdir : path-like
        Output directory (created).
    config : Mare2DEMConfig, optional
        Supplies ``initial_rho``, ``target_rms`` and ``max_iterations``.
    stem : str, default "mare2dem"
        Base name; the model is ``<stem>.0.resistivity``.
    error_floor : float, default 0.05
        Relative impedance error floor (TE and TM).
    output_modes : str, default "all"
        Forwarded to ``make_mt_data_from_edi``.
    cell_y, z_first, depth_max, growth
        Core-grid controls, see :func:`mt_core_grid`.
    padding : float, optional
        Lateral and vertical padding around the core grid in metres.
        Default: twice the larger core dimension, at least 50 km, so the
        model boundaries stay far from the receivers at every period.
    **edi_kwargs
        Extra ``make_mt_data_from_edi`` options.

    Returns
    -------
    dict
        ``files`` ({data, poly, model, settings}), ``run_stem`` (argument
        for :meth:`Mare2DEMRunner.run`), ``y_edges``, ``z_edges``,
        ``rx_y``, ``rx_z``, ``frequencies`` and ``n_parameters``.
    """
    import numpy as np

    from .edi import make_mt_data_from_edi
    from .grid_to_m2d import grid_to_mare2dem

    cfg = config or Mare2DEMConfig()
    dest = Path(workdir)
    dest.mkdir(parents=True, exist_ok=True)
    data_path = dest / f"{stem}.emdata"
    make_mt_data_from_edi(
        sites, data_path, output_modes=output_modes,
        error_floor_te=error_floor, error_floor_tm=error_floor,
        **edi_kwargs,
    )
    em = read_emdata(data_path)
    if em.mt is None or not len(em.mt.receivers):
        raise ValueError(f"{data_path.name} contains no MT receivers")
    rx = np.asarray(em.mt.receivers, dtype=float)
    freqs = np.asarray(em.mt.frequencies, dtype=float)
    rho0 = max(float(cfg.initial_rho), 1e-3)
    y_edges, z_edges = mt_core_grid(
        rx[:, 1], freqs, initial_rho=rho0, cell_y=cell_y,
        z_first=z_first, depth_max=depth_max, growth=growth,
    )
    yc = 0.5 * (y_edges[:-1] + y_edges[1:])
    zc = 0.5 * (z_edges[:-1] + z_edges[1:])
    Y, Z = np.meshgrid(yc, zc)
    if padding is None:
        core = max(y_edges[-1] - y_edges[0], z_edges[-1])
        padding = max(50000.0, 2.0 * core)
    files = grid_to_mare2dem(
        Y, Z, np.full(Y.shape, rho0),
        padding_y=padding, padding_z=padding, out_dir=dest,
        model_name=stem, data_file=data_path.name,
        target_misfit=float(cfg.target_rms),
        max_iterations=int(cfg.max_iterations),
    )
    return {
        "files": {"data": data_path, "poly": files["poly"],
                  "model": files["resistivity"],
                  "settings": files["settings"]},
        "run_stem": files["resistivity"].name[: -len(".resistivity")],
        "y_edges": y_edges,
        "z_edges": z_edges,
        "rx_y": rx[:, 1],
        "rx_z": rx[:, 2],
        "frequencies": freqs,
        "n_parameters": int(Y.size),
        "padding": float(padding),
    }
