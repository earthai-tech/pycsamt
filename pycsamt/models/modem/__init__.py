# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""pycsamt.models.modem — Python interface to the ModEM MT inversion code.

ModEM (Modular EM inversion system) supports 2-D and 3-D magnetotelluric
inversion using a nonlinear conjugate-gradient (NLCG) algorithm.

References
----------
Egbert, G.D. & Kelbert, A. (2012). Computational recipes for electromagnetic
inverse problems. Geophysical Journal International, 189(1), 251-267.

Kelbert, A., Meqbel, N., Egbert, G.D., Tandon, K. (2014). ModEM: A modular
system for inversion of electromagnetic geophysical data.
Computers & Geosciences, 66, 40-53.

Quick start — 3D
-----------------
>>> from pycsamt.models.modem import ModEmConfig, InputBuilder
>>> cfg = ModEmConfig(mode='3d', initial_rho=100.0)
>>> builder = InputBuilder(config=cfg)
>>> builder.build(edi_source, workdir='./modem_run')
>>> from pycsamt.models.modem import ModEmRunner
>>> ModEmRunner(workdir='./modem_run', config=cfg).run()
>>> from pycsamt.models.modem import InversionResult
>>> result = InversionResult(workdir='./modem_run')
>>> result.plot_model().savefig('model.png')
"""

from .builder import InputBuilder
from .config import ModEmConfig
from .control import ModEmControl
from .covariance import ModEmCovariance
from .data import ModEmData
from .iotools import (
    ImpedanceFile,
    ZBlock,
    convert_z2d,
    convert_z3d,
    imp_units_factor,
    interp_model3d,
    interp_z3d,
    linear_to_loge,
    log10_to_loge,
    loge_to_linear,
    loge_to_log10,
    read_mackie2d,
    read_mackie3d,
    read_z2d_old,
    read_z3d_old,
    skin_depth,
    write_mackie2d,
    write_mackie3d,
    write_meshtools3d,
    write_z2d_list,
    write_z2d_old,
    write_z3d_list,
    write_z3d_old,
)
from .log import ModEmLog
from .model2d import ModEmModel2D
from .model3d import ModEmModel3D
from .plot import (
    PlotMisfit,
    PlotModel2D,
    PlotModel3D,
    PlotPseudo,
    PlotResponse,
)
from .results import InversionResult
from .runner import ModEmRunner
from .validation import (
    ModEmFileType,
    detect_file_type,
    is_control_file,
    is_covariance_file,
    is_data_file,
    is_log_file,
    is_model_2d_file,
    is_model_3d_file,
    is_model_file,
)

__all__ = [
    # config
    "ModEmConfig",
    # validation
    "ModEmFileType",
    "detect_file_type",
    "is_data_file",
    "is_model_file",
    "is_model_2d_file",
    "is_model_3d_file",
    "is_covariance_file",
    "is_control_file",
    "is_log_file",
    # I/O objects
    "ModEmData",
    "ModEmModel2D",
    "ModEmModel3D",
    "ModEmCovariance",
    "ModEmControl",
    "ModEmLog",
    "InversionResult",
    "ModEmRunner",
    # plots
    "PlotMisfit",
    "PlotModel2D",
    "PlotModel3D",
    "PlotResponse",
    "PlotPseudo",
    # builder
    "InputBuilder",
    # iotools — interpolation
    "interp_model3d",
    "interp_z3d",
    # iotools — export
    "write_meshtools3d",
    # iotools — utils
    "skin_depth",
    "imp_units_factor",
    "loge_to_log10",
    "log10_to_loge",
    "loge_to_linear",
    "linear_to_loge",
    # iotools — Mackie format
    "read_mackie2d",
    "write_mackie2d",
    "read_mackie3d",
    "write_mackie3d",
    # iotools — impedance old format
    "ZBlock",
    "ImpedanceFile",
    "read_z3d_old",
    "write_z3d_old",
    "read_z2d_old",
    "write_z2d_old",
    "write_z3d_list",
    "write_z2d_list",
    "convert_z3d",
    "convert_z2d",
]
