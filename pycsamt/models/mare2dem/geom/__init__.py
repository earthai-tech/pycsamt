# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""MARE2DEM geometry utilities.

Python ports of the MATLAB ``a_util/geom/`` and
``a_util/mapping/`` utilities.

Submodules
----------
topo
    Topography interpolation and slope computation
    (``m2d_parseTopo``).
intersections
    Line-segment intersection tests (``m2d_getIntersections``).
simplify
    Douglas-Peucker polyline simplification (``m2d_dpsimplify``).
centroids
    Area-weighted triangle region centroids (``m2d_getCentroids``).
utm
    UTM ↔ Longitude/Latitude conversion (``LonLat2UTM``,
    ``UTM2LonLat``).
area_of_interest
    Estimate display / mesh y-z bounds from survey geometry
    (``m2d_estimateAreaOfInterest``).
triangle_regions
    Assign FEM triangles to PSLG-bounded regions
    (``m2d_getTriangleRegions``).
line_orientation
    Survey profile orientation from UTM stations
    (``m2d_getLineOrientation``).
simplify_poly
    PSLG polygon simplification by removing collinear interior nodes
    (``m2d_simplify_poly``).
"""

from .topo import parse_topo, topo_depth, topo_slope
from .intersections import get_intersections, do_rects_overlap
from .simplify import dp_simplify
from .centroids import get_centroids, triangle_centroids, triangle_areas
from .utm import lonlat_to_utm, utm_to_lonlat, ELLIPSOIDS
from .area_of_interest import estimate_area_of_interest
from .triangle_regions import get_triangle_regions
from .line_orientation import get_line_orientation, project_onto_line
from .simplify_poly import simplify_poly

__all__ = [
    # topo
    "parse_topo",
    "topo_depth",
    "topo_slope",
    # intersections
    "get_intersections",
    "do_rects_overlap",
    # polyline simplify (Douglas-Peucker)
    "dp_simplify",
    # centroids
    "get_centroids",
    "triangle_centroids",
    "triangle_areas",
    # utm
    "lonlat_to_utm",
    "utm_to_lonlat",
    "ELLIPSOIDS",
    # area of interest
    "estimate_area_of_interest",
    # triangle region assignment
    "get_triangle_regions",
    # line orientation
    "get_line_orientation",
    "project_onto_line",
    # PSLG polygon simplify
    "simplify_poly",
]
