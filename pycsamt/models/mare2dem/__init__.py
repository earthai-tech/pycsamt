# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""pycsamt.models.mare2dem — Python interface to the MARE2DEM EM inversion code.

MARE2DEM (Modeling with Adaptively Refined Elements for 2.5-D
Electromagnetic inversion) implements finite-element MT and CSEM
forward modelling and Occam-style regularized inversion. It is an
MPI-parallel Fortran code developed by Kerry Key at LDEO, Columbia
University.

Because the compiled binary is ~49 MB, the source tree is **not**
bundled in the repository. Use :class:`SourceManager` to download and
build it before running inversions.

Ported MATLAB scripts (complete)
---------------------------------
**I/O** (``a_util/io/`` + ``a_util/data/``):

* ``.emdata`` / ``.EMResp`` reader + writer
* ``.resistivity`` reader + writer (all anisotropy modes)
* Triangle ``.poly`` PSLG reader + writer
* ``.settings`` parallel-decomposition writer
* ``.emdata_group`` data-group file reader + writer
* Group-RMS CSV log reader
* Most-recently-modified file finder

**Data management**:

* Data-type code lookup table
* High-level ``.emdata`` builder from MT / CSEM survey configs
* ZMM impedance file reader + MT data-file builder
* Synthetic noise addition + ``make_synthetic_data`` wrapper
* Multi-file merge utility

**Geometry** (``a_util/geom/``, ``a_util/mapping/``):

* Topography interpolation + slope angles
* Line-segment intersections (bounding-box pre-filtered)
* Douglas-Peucker polyline simplification
* PSLG polygon simplification (collinear node removal)
* Area-weighted triangle region centroids
* Survey-profile line orientation
* Station-to-profile projection
* UTM ↔ Lon/Lat conversion (``pyproj`` with pure-Python WGS-84 fallback)
* Survey area-of-interest estimator
* Triangle FEM region flood-fill assignment

**Model construction**:

* 2-D resistivity grid → MARE2DEM ``.poly`` + ``.resistivity``
* Topography import + profile projection

**Model comparison**:

* Log10 (or custom) difference of two ``.resistivity`` files

**Plotting**:

* RMS convergence curve
* Survey map (Rx/Tx positions in UTM)
* Receiver geometry QC (6-panel)
* Transmitter geometry QC
* ``.poly`` PSLG mesh plot
* Resistivity section stub (pending full mesh parser)

Quick start
-----------
Step 1 — Download and compile (once per machine)::

    from pycsamt.models.mare2dem import SourceManager

    sm = SourceManager(verbose=1)
    sm.download()  # git-clone from Bitbucket
    sm.build()  # requires Intel oneAPI + MKL

Step 2 — Create a data file from MT survey parameters::

    import numpy as np
    from pycsamt.models.mare2dem import MTSurveyConfig, make_data_file

    mt = MTSurveyConfig(
        frequencies=np.logspace(-3, 3, 20),
        rx_y=np.linspace(-5000, 5000, 20),
        rx_type="marine",
        lTE=True,
        lTM=True,
    )
    em = make_data_file("survey.emdata", topo=-1000.0, mt=mt)

Step 3 — Prepare a full input set and run::

    from pycsamt.models.mare2dem import (
        Mare2DEMConfig,
        InputBuilder,
        Mare2DEMRunner,
    )

    cfg = Mare2DEMConfig(initial_rho=1.0, n_procs=8)
    InputBuilder(cfg).build("survey.emdata", workdir="./run")
    result = Mare2DEMRunner("./run", cfg).run("mare2dem")
    result.print_summary()

References
----------
Key, K. (2016). MARE2DEM: A 2-D inversion code for controlled-source
electromagnetic and magnetotelluric data. *Geophysical Journal
International*, 207(1), 571–588. doi:10.1093/gji/ggw290.
"""

# ---- configuration ----
from .builder import InputBuilder
from .config import Mare2DEMConfig

# ---- compatibility wrappers (original stub API) ----
from .data import EMData

# ---- model diff ----
from .diff import diff_resistivity

# ---- EDI → MT data builder ----
from .edi import make_mt_data_from_edi, stations_from_edi

# ---- geometry utilities ----
from .geom import (
    ELLIPSOIDS,
    do_rects_overlap,
    dp_simplify,
    estimate_area_of_interest,
    get_centroids,
    get_intersections,
    get_line_orientation,
    get_triangle_regions,
    lonlat_to_utm,
    parse_topo,
    project_onto_line,
    simplify_poly,
    topo_depth,
    topo_slope,
    triangle_areas,
    triangle_centroids,
    utm_to_lonlat,
)

# ---- model construction from grid ----
from .grid_to_m2d import grid_to_mare2dem

# ---- topography import ----
from .import_topo import TopoConfig, TopoProfile, import_topo

# ---- I/O tools (full parsers) ----
from .iotools import (
    # data codes
    DATA_CODES,
    CSEMConfig,
    # data group
    DataGroupFile,
    DCConfig,
    # emdata
    EMDataFile,
    # group RMS log
    GroupRMSLog,
    MTConfig,
    # poly
    PolyFile,
    # resistivity
    ResistivityFile,
    # settings
    SettingsFile,
    UTMOrigin,
    code_component,
    code_label,
    code_representation,
    # most recent
    get_most_recent,
    is_csem_code,
    is_mt_code,
    read_data_group,
    read_emdata,
    read_group_rms_log,
    read_poly,
    read_resistivity,
    write_data_group,
    write_emdata,
    write_poly,
    write_resistivity,
    write_settings,
)

# ---- log parsers ----
from .log import (
    GroupRMSLog,
    IterationRecord,
    Mare2DEMLog,
    read_group_rms_log,
)

# ---- file merge ----
from .merge import merge_data_files, merge_emdata
from .mesh import PolyMesh, ResistivityModel

# ---- synthetic noise ----
from .noise import (
    NoiseConfig,
    add_synthetic_noise,
    make_synthetic_data,
)

# ---- plots ----
from .plot import (
    PlotConvergence,
    PlotModel,
    PlotResponse,
    PlotRxParams,
    PlotSurveyLayout,
    PlotTxParams,
    plot_poly,
)

# ---- run lifecycle ----
from .results import InversionResult
from .runner import Mare2DEMRunner

# ---- source management ----
from .source import SourceManager

# ---- survey data builder ----
from .survey import (
    CSEMSurveyConfig,
    MTSurveyConfig,
    make_data_file,
)

# ---- file-type detection ----
from .validation import (
    Mare2DEMFileType,
    detect_file_type,
    is_emdata_file,
    is_log_file,
    is_resistivity_file,
    is_response_file,
    is_settings_file,
)

# ---- ZMM reader + MT data builder ----
from .zmm import (
    ZMMStation,
    make_mt_data_from_stations,
    make_mt_data_from_zmm,
    read_zmm,
)

__all__ = [
    # config
    "Mare2DEMConfig",
    # source management
    "SourceManager",
    # file-type detection
    "Mare2DEMFileType",
    "detect_file_type",
    "is_emdata_file",
    "is_resistivity_file",
    "is_settings_file",
    "is_log_file",
    "is_response_file",
    # data codes
    "DATA_CODES",
    "code_label",
    "code_component",
    "code_representation",
    "is_mt_code",
    "is_csem_code",
    # emdata I/O
    "EMDataFile",
    "UTMOrigin",
    "CSEMConfig",
    "MTConfig",
    "DCConfig",
    "read_emdata",
    "write_emdata",
    # resistivity I/O
    "ResistivityFile",
    "read_resistivity",
    "write_resistivity",
    # poly I/O
    "PolyFile",
    "read_poly",
    "write_poly",
    # settings I/O
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
    # geometry
    "parse_topo",
    "topo_depth",
    "topo_slope",
    "get_intersections",
    "do_rects_overlap",
    "dp_simplify",
    "simplify_poly",
    "get_centroids",
    "triangle_centroids",
    "triangle_areas",
    "lonlat_to_utm",
    "utm_to_lonlat",
    "ELLIPSOIDS",
    "estimate_area_of_interest",
    "get_triangle_regions",
    "get_line_orientation",
    "project_onto_line",
    # survey builder
    "MTSurveyConfig",
    "CSEMSurveyConfig",
    "make_data_file",
    # ZMM
    "ZMMStation",
    "read_zmm",
    "make_mt_data_from_zmm",
    "make_mt_data_from_stations",
    # EDI → emdata
    "stations_from_edi",
    "make_mt_data_from_edi",
    # noise
    "NoiseConfig",
    "add_synthetic_noise",
    "make_synthetic_data",
    # merge
    "merge_data_files",
    "merge_emdata",
    # grid to MARE2DEM
    "grid_to_mare2dem",
    # topo import
    "TopoConfig",
    "TopoProfile",
    "import_topo",
    # model diff
    "diff_resistivity",
    # compatibility wrappers
    "EMData",
    "ResistivityModel",
    "PolyMesh",
    # logs
    "Mare2DEMLog",
    "IterationRecord",
    # run lifecycle
    "InversionResult",
    "Mare2DEMRunner",
    "InputBuilder",
    # plots
    "PlotConvergence",
    "PlotSurveyLayout",
    "PlotRxParams",
    "PlotTxParams",
    "plot_poly",
    "PlotModel",
    "PlotResponse",
]
