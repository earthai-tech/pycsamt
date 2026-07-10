# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""ModEM iotools — Python translations of the MATLAB ioAscii/ scripts."""

from .export import write_meshtools3d
from .impedance import (
    ImpedanceFile,
    ZBlock,
    convert_z2d,
    convert_z3d,
    read_z2d_old,
    read_z3d_old,
    write_z2d_list,
    write_z2d_old,
    write_z3d_list,
    write_z3d_old,
)
from .interpolate import interp_model3d, interp_z3d
from .mackie import (
    read_mackie2d,
    read_mackie3d,
    write_mackie2d,
    write_mackie3d,
)
from .utils import (
    imp_units_factor,
    linear_to_loge,
    log10_to_loge,
    loge_to_linear,
    loge_to_log10,
    skin_depth,
)

__all__ = [
    # interpolation
    "interp_model3d",
    "interp_z3d",
    # export
    "write_meshtools3d",
    # utils
    "skin_depth",
    "imp_units_factor",
    "loge_to_log10",
    "log10_to_loge",
    "loge_to_linear",
    "linear_to_loge",
    # Mackie format
    "read_mackie2d",
    "write_mackie2d",
    "read_mackie3d",
    "write_mackie3d",
    # Impedance old format
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
