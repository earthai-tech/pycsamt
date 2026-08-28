"""Repeatable PCBH parser and rendering micro-benchmark."""

from __future__ import annotations

import copy
import json
from dataclasses import asdict, dataclass
from importlib import resources
from time import perf_counter

from .jsonio import pcbh_from_dict, pcbh_to_dict
from .render import build_render_model


@dataclass(frozen=True)
class PCBHBenchmarkResult:
    """Measured timings and output size for a synthetic project."""

    boreholes: int
    json_bytes: int
    render_points: int
    parse_seconds: float
    render_seconds: float

    def to_dict(self) -> dict[str, int | float]:
        """Return a JSON-serializable result."""
        return asdict(self)


def benchmark_pcbh(boreholes: int = 1000) -> PCBHBenchmarkResult:
    """Benchmark canonical parsing and rendering of vertical boreholes."""
    if not isinstance(boreholes, int) or boreholes < 1:
        raise ValueError("boreholes must be a positive integer")
    example = resources.files(__package__).joinpath(
        "examples/minimal-vertical.pcbh.json"
    )
    root = json.loads(example.read_text(encoding="utf-8"))
    prototype = root["boreholes"][0]
    root["boreholes"] = []
    for index in range(boreholes):
        hole = copy.deepcopy(prototype)
        hole["id"] = hole["name"] = f"BH-{index + 1:06d}"
        hole["collar"]["x"] += float(index % 100) * 25.0
        hole["collar"]["y"] += float(index // 100) * 25.0
        root["boreholes"].append(hole)

    started = perf_counter()
    document = pcbh_from_dict(root)
    payload = json.dumps(pcbh_to_dict(document), separators=(",", ":"))
    parse_seconds = perf_counter() - started
    started = perf_counter()
    rendered = build_render_model(document)
    render_seconds = perf_counter() - started
    points = sum(len(hole.centerline.points) for hole in rendered.boreholes)
    return PCBHBenchmarkResult(
        boreholes=boreholes,
        json_bytes=len(payload.encode("utf-8")),
        render_points=points,
        parse_seconds=parse_seconds,
        render_seconds=render_seconds,
    )
