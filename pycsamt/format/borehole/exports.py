# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""GeoJSON, VTP, and glTF exports for PCBH visualization subsets."""

from __future__ import annotations

import base64
import json
import math
import struct
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from ...api.property import PyCSAMTObject
from .render import build_render_model
from .schema import PCBHDocument

__all__ = [
    "ExportLoss",
    "PCBHExportReport",
    "write_geojson",
    "write_vtp",
    "write_gltf",
]


@dataclass(frozen=True, repr=False)
class ExportLoss(PyCSAMTObject):
    """One field or semantic category omitted from a view export."""

    code: str
    message: str


@dataclass(frozen=True, repr=False)
class PCBHExportReport(PyCSAMTObject):
    """Summary and explicit losses for a non-canonical export."""

    format: str
    boreholes: int
    features: int
    source_crs: str
    target_crs: str | None
    losses: tuple[ExportLoss, ...]


@dataclass(frozen=True)
class _TubeMesh:
    positions: np.ndarray
    normals: np.ndarray
    colors: np.ndarray
    measured_depth: np.ndarray
    material_ids: np.ndarray
    triangles: np.ndarray
    materials: tuple[tuple[str, str], ...]


def write_geojson(
    document: PCBHDocument,
    path: str | Path,
    *,
    target_crs: str = "EPSG:4326",
    include_z: bool = True,
) -> tuple[Path, PCBHExportReport]:
    """Write WGS84 collar Points and trajectory LineStrings as GeoJSON."""
    document.validate()
    transformer = _transformer(document.crs.horizontal, target_crs)
    render = build_render_model(document)
    features = []
    for hole in render.boreholes:
        collar = hole.collar
        cx, cy = _xy(transformer, collar.position[0], collar.position[1])
        collar_coords = [cx, cy]
        if include_z:
            collar_coords.append(collar.position[2])
        common = dict(collar.metadata)
        common["pcbh_document_id"] = document.document_id
        features.append(
            {
                "type": "Feature",
                "id": f"{hole.borehole_id}:collar",
                "geometry": {"type": "Point", "coordinates": collar_coords},
                "properties": {**common, "feature_type": "collar"},
            }
        )
        coordinates = []
        for point in hole.centerline.points:
            x, y = _xy(transformer, point.x, point.y)
            coordinate = [x, y]
            if include_z:
                coordinate.append(point.z)
            coordinates.append(coordinate)
        features.append(
            {
                "type": "Feature",
                "id": f"{hole.borehole_id}:trajectory",
                "geometry": {
                    "type": "LineString",
                    "coordinates": coordinates,
                },
                "properties": {
                    **common,
                    "feature_type": "trajectory",
                    "z_reference": document.crs.vertical,
                    "z_unit": document.crs.coordinate_unit,
                },
            }
        )
    payload = {"type": "FeatureCollection", "features": features}
    output = _prepare_path(path, ".geojson")
    output.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    losses = _common_losses(document, family=None)
    if include_z:
        losses.append(
            ExportLoss(
                "geojson_z_semantics",
                "GeoJSON Z retains PCBH elevation/reference; RFC 7946 only "
                "standardizes horizontal WGS84 coordinates.",
            )
        )
    return output, PCBHExportReport(
        "geojson",
        len(document.boreholes),
        len(features),
        document.crs.horizontal,
        target_crs,
        tuple(losses),
    )


def write_vtp(
    document: PCBHDocument,
    path: str | Path,
    *,
    family: str = "lithology",
    sides: int = 8,
) -> tuple[Path, PCBHExportReport]:
    """Write interval tubes as ASCII VTK XML PolyData (``.vtp``)."""
    mesh = _tube_mesh(document, family=family, sides=sides)
    root = ET.Element(
        "VTKFile",
        type="PolyData",
        version="1.0",
        byte_order="LittleEndian",
    )
    poly = ET.SubElement(root, "PolyData")
    piece = ET.SubElement(
        poly,
        "Piece",
        NumberOfPoints=str(len(mesh.positions)),
        NumberOfPolys=str(len(mesh.triangles)),
    )
    points = ET.SubElement(piece, "Points")
    _xml_array(points, "Float64", None, mesh.positions, components=3)
    point_data = ET.SubElement(piece, "PointData", Scalars="material_id")
    _xml_array(point_data, "Float64", "measured_depth", mesh.measured_depth)
    _xml_array(point_data, "Int32", "material_id", mesh.material_ids)
    _xml_array(point_data, "UInt8", "RGB", mesh.colors, components=3)
    polys = ET.SubElement(piece, "Polys")
    _xml_array(polys, "Int64", "connectivity", mesh.triangles.reshape(-1))
    offsets = np.arange(1, len(mesh.triangles) + 1, dtype=np.int64) * 3
    _xml_array(polys, "Int64", "offsets", offsets)
    field = ET.SubElement(piece, "FieldData")
    material_json = json.dumps(dict(mesh.materials), separators=(",", ":"))
    data = ET.SubElement(
        field,
        "DataArray",
        type="String",
        Name="pcbh_materials_json",
        NumberOfTuples="1",
        format="ascii",
    )
    data.text = material_json
    output = _prepare_path(path, ".vtp")
    ET.ElementTree(root).write(output, encoding="utf-8", xml_declaration=True)
    losses = _common_losses(document, family=family)
    return output, PCBHExportReport(
        "vtp",
        len(document.boreholes),
        len(mesh.triangles),
        document.crs.horizontal,
        document.crs.horizontal,
        tuple(losses),
    )


def write_gltf(
    document: PCBHDocument,
    path: str | Path,
    *,
    family: str = "lithology",
    sides: int = 8,
) -> tuple[Path, PCBHExportReport]:
    """Write browser-ready glTF 2.0 (``.gltf``) or binary GLB tubes."""
    output = Path(path)
    if output.suffix.lower() not in {".gltf", ".glb"}:
        output = output.with_suffix(".gltf")
    output.parent.mkdir(parents=True, exist_ok=True)
    mesh = _tube_mesh(document, family=family, sides=sides)
    gltf, binary = _gltf_payload(mesh, document)
    if output.suffix.lower() == ".glb":
        output.write_bytes(_glb_bytes(gltf, binary))
        fmt = "glb"
    else:
        gltf["buffers"][0]["uri"] = (
            "data:application/octet-stream;base64,"
            + base64.b64encode(binary).decode("ascii")
        )
        output.write_text(
            json.dumps(gltf, ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8",
        )
        fmt = "gltf"
    return output, PCBHExportReport(
        fmt,
        len(document.boreholes),
        len(mesh.triangles),
        document.crs.horizontal,
        document.crs.horizontal,
        tuple(_common_losses(document, family=family)),
    )


def _tube_mesh(
    document: PCBHDocument,
    *,
    family: str,
    sides: int,
) -> _TubeMesh:
    if isinstance(sides, bool) or not isinstance(sides, int) or sides < 3:
        raise ValueError("tube sides must be an integer >= 3")
    model = build_render_model(document, family=family)
    positions = []
    normals = []
    colors = []
    depths = []
    material_ids = []
    triangles = []
    material_map: dict[str, int] = {}
    for hole in model.boreholes:
        segments = hole.interval_segments or (hole.centerline,)
        for segment in segments:
            color = segment.color
            material_id = material_map.setdefault(color, len(material_map))
            radius = getattr(segment, "display_radius", hole.display_radius)
            points = segment.points
            start = len(positions)
            rgb = _hex_rgb(color)
            for index, point in enumerate(points):
                tangent = _tangent(points, index)
                u, v = _normal_basis(tangent)
                for side in range(sides):
                    angle = 2.0 * math.pi * side / sides
                    normal = math.cos(angle) * u + math.sin(angle) * v
                    positions.append(
                        np.array([point.x, point.y, point.z]) + radius * normal
                    )
                    normals.append(normal)
                    colors.append(rgb)
                    depths.append(point.md)
                    material_ids.append(material_id)
            for ring in range(len(points) - 1):
                a = start + ring * sides
                b = a + sides
                for side in range(sides):
                    nxt = (side + 1) % sides
                    triangles.extend(
                        [(a + side, b + side, b + nxt),
                         (a + side, b + nxt, a + nxt)]
                    )
    if not positions:
        raise ValueError(f"no {family!r} intervals or trajectories to export")
    materials = tuple(
        (str(index), color)
        for color, index in sorted(
            material_map.items(), key=lambda item: item[1]
        )
    )
    return _TubeMesh(
        np.asarray(positions, dtype=np.float64),
        np.asarray(normals, dtype=np.float64),
        np.asarray(colors, dtype=np.uint8),
        np.asarray(depths, dtype=np.float64),
        np.asarray(material_ids, dtype=np.int32),
        np.asarray(triangles, dtype=np.uint32),
        materials,
    )


def _gltf_payload(mesh: _TubeMesh, document: PCBHDocument):
    arrays = [
        mesh.positions.astype("<f4"),
        mesh.normals.astype("<f4"),
        mesh.colors.astype(np.float32) / 255.0,
        mesh.triangles.reshape(-1).astype("<u4"),
    ]
    binary = bytearray()
    views = []
    accessors = []
    targets = [34962, 34962, 34962, 34963]
    types = ["VEC3", "VEC3", "VEC3", "SCALAR"]
    components = [5126, 5126, 5126, 5125]
    for array, target, accessor_type, component in zip(
        arrays, targets, types, components
    ):
        while len(binary) % 4:
            binary.append(0)
        offset = len(binary)
        raw = array.tobytes(order="C")
        binary.extend(raw)
        views.append(
            {"buffer": 0, "byteOffset": offset, "byteLength": len(raw),
             "target": target}
        )
        accessor: dict[str, Any] = {
            "bufferView": len(views) - 1,
            "componentType": component,
            "count": int(array.shape[0]),
            "type": accessor_type,
        }
        if len(accessors) == 0:
            accessor["min"] = array.min(axis=0).tolist()
            accessor["max"] = array.max(axis=0).tolist()
        accessors.append(accessor)
    gltf = {
        "asset": {"version": "2.0", "generator": "pyCSAMT PCBH"},
        "scene": 0,
        "scenes": [{"nodes": [0]}],
        "nodes": [{"mesh": 0, "name": document.document_id}],
        "meshes": [{"primitives": [{
            "attributes": {"POSITION": 0, "NORMAL": 1, "COLOR_0": 2},
            "indices": 3,
            "material": 0,
        }]}],
        "materials": [{
            "name": "PCBH interval colors",
            "pbrMetallicRoughness": {
                "baseColorFactor": [1.0, 1.0, 1.0, 1.0],
                "metallicFactor": 0.0,
                "roughnessFactor": 0.8,
            },
            "doubleSided": True,
        }],
        "buffers": [{"byteLength": len(binary)}],
        "bufferViews": views,
        "accessors": accessors,
        "extras": {
            "pcbh_crs": document.crs.horizontal,
            "pcbh_vertical_crs": document.crs.vertical,
            "pcbh_materials": dict(mesh.materials),
        },
    }
    return gltf, bytes(binary)


def _glb_bytes(gltf: dict[str, Any], binary: bytes) -> bytes:
    json_chunk = json.dumps(gltf, separators=(",", ":")).encode("utf-8")
    json_chunk += b" " * ((-len(json_chunk)) % 4)
    binary += b"\x00" * ((-len(binary)) % 4)
    total = 12 + 8 + len(json_chunk) + 8 + len(binary)
    return (
        struct.pack("<4sII", b"glTF", 2, total)
        + struct.pack("<I4s", len(json_chunk), b"JSON")
        + json_chunk
        + struct.pack("<I4s", len(binary), b"BIN\x00")
        + binary
    )


def _transformer(source: str, target: str):
    try:
        from pyproj import CRS, Transformer
    except ImportError as error:
        if source != target:
            raise ImportError(
                "pyproj is required for GeoJSON reprojection"
            ) from error
        return None
    source_crs = CRS.from_user_input(source)
    target_crs = CRS.from_user_input(target)
    if target_crs.to_epsg() != 4326:
        raise ValueError("RFC 7946 GeoJSON target CRS must be EPSG:4326")
    if source_crs == target_crs:
        return None
    return Transformer.from_crs(source_crs, target_crs, always_xy=True)


def _xy(transformer, x, y):
    if transformer is None:
        return float(x), float(y)
    result = transformer.transform(x, y)
    return float(result[0]), float(result[1])


def _tangent(points, index):
    left = points[max(0, index - 1)]
    right = points[min(len(points) - 1, index + 1)]
    vector = np.array([right.x - left.x, right.y - left.y, right.z - left.z])
    norm = np.linalg.norm(vector)
    return vector / norm if norm else np.array([0.0, 0.0, -1.0])


def _normal_basis(tangent):
    helper = np.array([0.0, 0.0, 1.0])
    if abs(float(np.dot(tangent, helper))) > 0.9:
        helper = np.array([1.0, 0.0, 0.0])
    u = np.cross(tangent, helper)
    u /= np.linalg.norm(u)
    return u, np.cross(tangent, u)


def _hex_rgb(color):
    value = color.lstrip("#")
    return tuple(int(value[index:index + 2], 16) for index in (0, 2, 4))


def _xml_array(parent, vtk_type, name, values, *, components=None):
    attributes = {"type": vtk_type, "format": "ascii"}
    if name is not None:
        attributes["Name"] = name
    if components is not None:
        attributes["NumberOfComponents"] = str(components)
    data = ET.SubElement(parent, "DataArray", **attributes)
    array = np.asarray(values).reshape(-1)
    data.text = " ".join(str(value) for value in array)


def _common_losses(document, *, family):
    losses = [
        ExportLoss(
            "noncanonical_view",
            "Export is a visualization subset; retain PCBH JSON as authority.",
        )
    ]
    if family is not None:
        other = sorted(
            {
                name
                for hole in document.boreholes
                for name in hole.interval_logs
                if name != family
            }
        )
        if other:
            losses.append(
                ExportLoss(
                    "interval_families_omitted",
                    "Omitted interval families: " + ", ".join(other),
                )
            )
    if any(hole.structures for hole in document.boreholes):
        losses.append(
            ExportLoss(
                "structures_omitted",
                "Structural observations are not encoded in tube geometry.",
            )
        )
    return losses


def _prepare_path(path, suffix):
    output = Path(path)
    if output.suffix.lower() != suffix:
        output = output.with_suffix(suffix)
    output.parent.mkdir(parents=True, exist_ok=True)
    return output
