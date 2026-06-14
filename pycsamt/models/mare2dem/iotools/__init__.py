# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""MARE2DEM file I/O tools.

Python ports of the MATLAB ``a_util/io/`` and related utilities.

Submodules
----------
codes
    Data-type code lookup table (``m2d_getDataCodeLookupTable``).
emdata
    ``.emdata`` / ``.EMResp`` reader and writer
    (``m2d_readEMData2DFile`` / ``m2d_writeEMData2DFile``).
resistivity
    ``.resistivity`` reader and writer
    (``m2d_readResistivity`` / ``m2d_writeResistivity``).
poly
    Triangle ``.poly`` PSLG reader and writer
    (``m2d_readPoly`` / ``m2d_writePoly``).
settings
    ``.settings`` parallel-decomposition writer
    (``m2d_writeSettingsFile``).
group_rms
    Group-level RMS log reader (``m2d_read_group_rms_log``).
data_group
    ``.emdata_group`` data-group file reader and writer
    (``m2d_readDataGroupFile`` / ``m2d_writeDataGroupFile``).
most_recent
    Most-recently-modified file finder (``m2d_getMostRecent``).
"""

from .codes import (
    DATA_CODES,
    code_label,
    code_component,
    code_representation,
    is_mt_code,
    is_csem_code,
)
from .emdata import (
    EMDataFile,
    UTMOrigin,
    CSEMConfig,
    MTConfig,
    DCConfig,
    read_emdata,
    write_emdata,
)
from .resistivity import (
    ResistivityFile,
    read_resistivity,
    write_resistivity,
)
from .poly import (
    PolyFile,
    read_poly,
    write_poly,
)
from .settings import (
    SettingsFile,
    write_settings,
)
from .group_rms import (
    GroupRMSLog,
    read_group_rms_log,
)
from .data_group import (
    DataGroupFile,
    read_data_group,
    write_data_group,
)
from .most_recent import get_most_recent

__all__ = [
    # codes
    "DATA_CODES",
    "code_label",
    "code_component",
    "code_representation",
    "is_mt_code",
    "is_csem_code",
    # emdata
    "EMDataFile",
    "UTMOrigin",
    "CSEMConfig",
    "MTConfig",
    "DCConfig",
    "read_emdata",
    "write_emdata",
    # resistivity
    "ResistivityFile",
    "read_resistivity",
    "write_resistivity",
    # poly
    "PolyFile",
    "read_poly",
    "write_poly",
    # settings
    "SettingsFile",
    "write_settings",
    # group RMS log
    "GroupRMSLog",
    "read_group_rms_log",
    # data group
    "DataGroupFile",
    "read_data_group",
    "write_data_group",
    # most recent file
    "get_most_recent",
]
