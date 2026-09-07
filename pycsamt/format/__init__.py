"""pyCSAMT Common Subsurface Format (PCSF) — backend-neutral inversion results.

See ``PYCSAMT-PCSF-INVERSION-FORMAT-PLAN.md`` at the repository root
for the full design and phase-by-phase status.

Public second-level namespaces (importable as ``pycsamt.format.<name>``
right after ``import pycsamt.format``, not only via an explicit
``from pycsamt.format.<name> import ...``):

- :mod:`pycsamt.format.adapters` — per-backend converters
  (``adapters.occam2d_to_pcsf``, ``adapters.modem3d_to_pcsf``,
  ``adapters.mare2dem_to_pcsf``), plus a solver-agnostic
  ``adapters.grid2d_to_pcsf``/``adapters.grid3d_to_pcsf``/
  ``adapters.mesh_to_pcsf`` for any AI/DL inversion result (see
  :mod:`pycsamt.format.adapters.generic` and
  :mod:`pycsamt.format.provenance`).
- :mod:`pycsamt.format.schema` — the dataclasses re-exported at this
  top level (``Grid2DGeometry``, ``PCSFModel``, ...).
- :mod:`pycsamt.format.io` — ``read_pcsf``/``write_pcsf``, also
  re-exported here.
- :mod:`pycsamt.format.text` — ``read_pcsm``/``write_pcsm``/
  ``pcsf_to_pcsm``/``pcsm_to_pcsf``: PCSM, the hand-editable ASCII
  sibling of a ``.pcsf`` file, also re-exported here.
- :mod:`pycsamt.format.multiline`, :mod:`pycsamt.format.topography`,
  :mod:`pycsamt.format.pointcloud`, :mod:`pycsamt.format.regrid` —
  likewise re-exported here.
- :mod:`pycsamt.format.pointset` — PCPT, a tiny ``.pcpt.json`` format
  for targets / points of interest (``read_points``, ``points_from_csv``,
  ``points_from_xlsx``), re-exported here.
- :mod:`pycsamt.format.geology` — PCGL, a tiny ``.pcgl.json`` format for
  a resistivity-to-geology interpretation legend (``GeologyLegend``,
  ``read_legend``, ``write_legend``, ``legend_from_csv``), re-exported
  here.
- :mod:`pycsamt.format.structure` — PCGS, a tiny ``.pcgs.json`` format
  for field structural evidence -- planar/linear measurements and fault
  traces (``StructModel``, ``read_structure``, ``write_structure``,
  ``structure_from_csv``), re-exported here.
"""

from . import adapters, borehole
from .borehole import (
    align_pcbh_to_pcsf,
    borehole_from_las,
    boreholes_from_csv,
    boreholes_from_csv_directory,
    build_render_model,
    desurvey,
    embed_pcbh,
    extract_pcbh,
    read_pcbh,
    reference_pcbh,
    write_csv_directory,
    write_geojson,
    write_gltf,
    write_las_subset,
    write_pcbh,
    write_vtp,
)
from .detect import SourceKind, describe_source, detect_source
from .geology import (
    PCGL_VERSION,
    GeologyLegend,
    GeologyLegendValidationError,
    legend_from_csv,
    legend_from_dict,
    legend_to_dict,
    read_legend,
    write_legend,
)
from .io import read_pcsf, write_pcsf
from .structure import (
    PCGS_VERSION,
    StructModel,
    StructModelValidationError,
    read_structure,
    structure_from_csv,
    structure_from_dict,
    structure_to_dict,
    write_structure,
)
from .multiline import (
    build_multiline_pcsf,
    line_offsets_from_stations,
    multiline_pcsf_to_profiles,
    stack_lines_to_common_grid,
)
from .pointcloud import PointCloud, pcsf_to_point_cloud
from .pointset import (
    PCPT_VERSION,
    Point,
    PointSet,
    PointSetValidationError,
    point_set_from_dict,
    point_set_to_dict,
    points_from_csv,
    points_from_xlsx,
    read_points,
    write_points,
)
from .provenance import ModelProvenance, compute_checkpoint_hash
from .regrid import mesh_to_grid2d
from .schema import (
    DERIVATION_METHODS,
    GEOMETRY_KINDS,
    PCSF_VERSION,
    RESISTIVITY_ENCODINGS,
    RESISTIVITY_UNIT,
    TOPOGRAPHY_KINDS,
    DerivedVolume,
    Grid2DGeometry,
    Grid3DGeometry,
    LineEntry,
    MultilineGeometry,
    PCSFModel,
    StationTable,
    TopographyPerStation,
    TopographyRaster,
    UnstructuredMeshGeometry,
)
from .text import (
    pcsf_to_pcsm,
    pcsm_to_pcsf,
    peek_kind,
    read_pcsf_or_pcsm,
    read_pcsm,
    write_pcsm,
)
from .topo_source import (
    TopoAttribution,
    TopoTable,
    attribute_topo,
    read_topo_file,
    resolve_topo,
    topo_from_sites,
)
from .topography import (
    topography_from_elevation_file,
    topography_from_grid,
    topography_from_map_data,
    topography_raster_to_grid,
    topography_to_elev_map,
)

__all__ = [
    "adapters",
    "borehole",
    "read_pcbh",
    "write_pcbh",
    "desurvey",
    "boreholes_from_csv",
    "boreholes_from_csv_directory",
    "write_csv_directory",
    "borehole_from_las",
    "write_las_subset",
    "build_render_model",
    "embed_pcbh",
    "reference_pcbh",
    "extract_pcbh",
    "align_pcbh_to_pcsf",
    "write_geojson",
    "write_vtp",
    "write_gltf",
    "SourceKind",
    "detect_source",
    "describe_source",
    "PCSF_VERSION",
    "RESISTIVITY_UNIT",
    "GEOMETRY_KINDS",
    "RESISTIVITY_ENCODINGS",
    "TOPOGRAPHY_KINDS",
    "DERIVATION_METHODS",
    "Grid2DGeometry",
    "Grid3DGeometry",
    "UnstructuredMeshGeometry",
    "LineEntry",
    "DerivedVolume",
    "MultilineGeometry",
    "StationTable",
    "TopographyPerStation",
    "TopographyRaster",
    "PCSFModel",
    "read_pcsf",
    "write_pcsf",
    "read_pcsm",
    "write_pcsm",
    "pcsf_to_pcsm",
    "pcsm_to_pcsf",
    "read_pcsf_or_pcsm",
    "peek_kind",
    "build_multiline_pcsf",
    "multiline_pcsf_to_profiles",
    "line_offsets_from_stations",
    "stack_lines_to_common_grid",
    "topography_from_map_data",
    "topography_from_elevation_file",
    "topography_to_elev_map",
    "topography_from_grid",
    "topography_raster_to_grid",
    "TopoTable",
    "TopoAttribution",
    "read_topo_file",
    "topo_from_sites",
    "attribute_topo",
    "resolve_topo",
    "PointCloud",
    "pcsf_to_point_cloud",
    "PCPT_VERSION",
    "Point",
    "PointSet",
    "PointSetValidationError",
    "read_points",
    "write_points",
    "points_from_csv",
    "points_from_xlsx",
    "point_set_to_dict",
    "point_set_from_dict",
    "PCGL_VERSION",
    "GeologyLegend",
    "GeologyLegendValidationError",
    "read_legend",
    "write_legend",
    "legend_from_csv",
    "legend_to_dict",
    "legend_from_dict",
    "PCGS_VERSION",
    "StructModel",
    "StructModelValidationError",
    "read_structure",
    "write_structure",
    "structure_from_csv",
    "structure_to_dict",
    "structure_from_dict",
    "ModelProvenance",
    "compute_checkpoint_hash",
    "mesh_to_grid2d",
]
