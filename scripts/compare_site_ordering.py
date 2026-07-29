"""Compare station ordering strategies for an EDI directory."""

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def _dms(value: str) -> float:
    parts = value.strip().strip('"').split(":")
    deg, minute, second = map(float, parts)
    sign = -1.0 if deg < 0 else 1.0
    return sign * (abs(deg) + minute / 60.0 + second / 3600.0)


def _natural_key(value: str) -> tuple:
    return tuple(
        int(p) if p.isdigit() else p.lower() for p in re.split(r"(\d+)", value)
    )


def _read_head(path: Path) -> dict[str, object]:
    text = path.read_text(encoding="utf-8", errors="replace")
    head = text.split(">INFO", 1)[0]
    fields = {}
    for key in ("DATAID", "LAT", "LONG", "ELEV"):
        match = re.search(rf"(?mi)^\s*{key}\s*=\s*([^\r\n]+)", head)
        fields[key] = match.group(1).strip().strip('"') if match else ""
    return {
        "file": path.name,
        "station": fields["DATAID"] or path.stem,
        "lat": _dms(str(fields["LAT"])),
        "lon": _dms(str(fields["LONG"])),
        "elev": float(fields["ELEV"]),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    rows = [_read_head(p) for p in args.directory.glob("*.edi")]
    lat = np.array([r["lat"] for r in rows], dtype=float)
    lon = np.array([r["lon"] for r in rows], dtype=float)
    lat0, lon0 = float(lat.mean()), float(lon.mean())
    x = (lon - lon0) * 111_000.0 * math.cos(math.radians(lat0))
    y = (lat - lat0) * 111_000.0
    xy = np.column_stack((x, y))
    _, _, vh = np.linalg.svd(xy - xy.mean(axis=0), full_matrices=False)
    axis = vh[0]
    # Fix PCA's arbitrary sign so the smallest natural filename starts the line.
    ref = min(
        range(len(rows)), key=lambda i: _natural_key(str(rows[i]["file"]))
    )
    chainage = xy @ axis
    if chainage[ref] > np.median(chainage):
        axis *= -1
        chainage *= -1
    chainage -= chainage.min()
    azimuth = math.degrees(math.atan2(axis[0], axis[1])) % 360.0
    for i, row in enumerate(rows):
        row["x_m"] = x[i]
        row["y_m"] = y[i]
        row["chainage_m"] = chainage[i]

    strategies = {
        "Input / filename lexical": sorted(
            range(len(rows)), key=lambda i: str(rows[i]["file"])
        ),
        "Station name lexical": sorted(
            range(len(rows)), key=lambda i: str(rows[i]["station"])
        ),
        "Filename natural numeric": sorted(
            range(len(rows)), key=lambda i: _natural_key(str(rows[i]["file"]))
        ),
        "Longitude": sorted(range(len(rows)), key=lambda i: lon[i]),
        "Latitude": sorted(range(len(rows)), key=lambda i: lat[i]),
        "Profile chainage (PCA)": sorted(
            range(len(rows)), key=lambda i: chainage[i]
        ),
    }

    def metrics(order: list[int]) -> tuple[float, int]:
        pts = xy[order]
        length = float(np.linalg.norm(np.diff(pts, axis=0), axis=1).sum())
        dc = np.diff(chainage[order])
        reversals = int(np.sum(dc < -1e-9))
        return length, reversals

    with (args.output / "ordering_summary.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "strategy",
                "path_length_m",
                "chainage_reversals",
                "station_sequence",
            ]
        )
        for name, order in strategies.items():
            length, reversals = metrics(order)
            writer.writerow(
                [
                    name,
                    f"{length:.2f}",
                    reversals,
                    " | ".join(str(rows[i]["station"]) for i in order),
                ]
            )

    with (args.output / "stations_by_chainage.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(sorted(rows, key=lambda r: float(r["chainage_m"])))

    fig, axes = plt.subplots(2, 3, figsize=(17, 11), constrained_layout=True)
    for ax, (name, order) in zip(axes.flat, strategies.items()):
        pts = xy[order]
        ax.plot(pts[:, 0], pts[:, 1], "-o", ms=4, lw=1.2)
        for step, idx in enumerate(order, 1):
            ax.annotate(
                str(step),
                xy[idx],
                xytext=(3, 3),
                textcoords="offset points",
                fontsize=7,
            )
        length, reversals = metrics(order)
        ax.set_title(f"{name}\npath={length:.0f} m; reversals={reversals}")
        ax.set_xlabel("East offset (m)")
        ax.set_ylabel("North offset (m)")
        ax.axis("equal")
        ax.grid(alpha=0.25)
    fig.suptitle(
        f"L22PLT station ordering comparison — PCA azimuth {azimuth:.1f}°",
        fontsize=15,
    )
    fig.savefig(args.output / "site_ordering_comparison.png", dpi=180)
    print(f"stations={len(rows)} azimuth={azimuth:.3f} output={args.output}")


if __name__ == "__main__":
    main()
