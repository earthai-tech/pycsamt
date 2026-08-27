"""Create browsable PCSM examples from the bundled PCSF conversion outputs.

Run from the repository root::

    python examples/pcsm_conversion/run_demo.py

The Occam2D and MARE2DEM text projections are retained in full. Concise,
clearly marked, non-parseable ModEM and MARE2DEM excerpts are also written
for documentation. The real ModEM volume contains 590,400 cells and expands
to about 19.7 MB as plain text; users can retain the complete or
gzip-compressed projection with ``--full-modem``.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from pycsamt.format import pcsf_to_pcsm, read_pcsm

ROOT = Path(__file__).resolve().parents[2]
PCSF_OUTPUT = ROOT / "examples" / "pcsf_conversion_demo" / "output"
AI_PCSF_OUTPUT = ROOT / "examples" / "ai_inversion_pcsf" / "output"
OUTPUT = Path(__file__).resolve().parent / "output"


def _block_excerpt(lines: list[str], keyword: str, keep: int = 2) -> list[str]:
    """Return a marked head/tail excerpt of one END_KEYWORD block."""
    start = next(
        i
        for i, line in enumerate(lines)
        if line == keyword or line.startswith(f"{keyword}  #")
    )
    end_keyword = f"END_{keyword}"
    end = next(
        i for i in range(start + 1, len(lines)) if lines[i] == end_keyword
    )
    body = lines[start + 1 : end]
    return [
        lines[start],
        *body[:keep],
        f"# ... {max(0, len(body) - 2 * keep)} data rows omitted ...",
        *body[-keep:],
        lines[end],
    ]


def _raw_line_excerpt(
    lines: list[str], keyword: str, *, max_chars: int = 340
) -> list[str]:
    """Return one raw (non-block) ``KEYWORD value`` line, wrapped with a
    truncation marker if it exceeds *max_chars* -- for a single-line
    field such as ``METADATA_JSON`` that has no ``END_KEYWORD``."""
    line = next(line for line in lines if line.startswith(f"{keyword} "))
    if len(line) <= max_chars:
        return [line]
    return [line[:max_chars] + " ... (truncated)"]


def _write_excerpt(
    full_path: Path,
    excerpt_path: Path,
    keywords: tuple[str, ...],
    regeneration_hint: str,
    *,
    raw_keywords: tuple[str, ...] = (),
    keep: int = 2,
) -> None:
    """Write real PCSM headers plus shortened representative blocks.

    *raw_keywords* names single-line (non-``END_KEYWORD``) fields, such
    as ``METADATA_JSON``, appended verbatim (truncated only if very
    long) after the block sections. *keep* is lowered to 1 for a small
    demo block where the default 2-head/2-tail excerpt would otherwise
    print the same row twice.
    """
    lines = full_path.read_text(encoding="utf8").splitlines()
    first_block = next(
        i
        for i, line in enumerate(lines)
        if line == keywords[0] or line.startswith(f"{keywords[0]}  #")
    )
    header = lines[:first_block]
    sections = [header]
    for keyword in keywords:
        sections.append(_block_excerpt(lines, keyword, keep=keep))
    for keyword in raw_keywords:
        sections.append(_raw_line_excerpt(lines, keyword))
    text = [
        "# DOCUMENTATION EXCERPT ONLY -- ellipses make this file",
        "# non-parseable. Generate the complete model with:",
        f"# {regeneration_hint}",
        "",
    ]
    for section in sections:
        text.extend(section)
        text.append("")
    excerpt_path.write_text("\n".join(text), encoding="utf8")


def main(*, full_modem: bool = False) -> None:
    """Generate and verify the three real-backend PCSM demonstrations."""
    OUTPUT.mkdir(parents=True, exist_ok=True)
    pairs = (
        ("occam2d_no_topo.pcsf", "occam2d.pcsm"),
        ("mare2dem.pcsf", "mare2dem.pcsm"),
    )
    for source_name, target_name in pairs:
        target = pcsf_to_pcsm(PCSF_OUTPUT / source_name, OUTPUT / target_name)
        restored = read_pcsm(target)
        size = target.stat().st_size
        print(f"{target.relative_to(ROOT)}: {restored.kind}, {size:,} bytes")

    mare_full = OUTPUT / "mare2dem.pcsm"
    mare_excerpt = OUTPUT / "mare2dem_excerpt.pcsm"
    _write_excerpt(
        mare_full,
        mare_excerpt,
        (
            "NODES",
            "CONNECTIVITY",
            "REGION_IDS",
            "RESISTIVITY",
            "RESISTIVITY_BY_REGION",
        ),
        "python examples/pcsm_conversion/run_demo.py",
    )
    print(
        f"{mare_excerpt.relative_to(ROOT)}: documented excerpt, "
        f"{mare_excerpt.stat().st_size:,} bytes"
    )

    # occam2d_with_bln_topo.pcsf carries real (well, illustrative-but-
    # populated) stations/lon/lat -- converting it to .pcsm and back
    # confirms that data (added via occam2d_to_pcsf's topo= parameter,
    # see pcsf_conversion_demo) survives the PCSF<->PCSM round trip
    # too, not just a raw .pcsf write.
    topo_source = PCSF_OUTPUT / "occam2d_with_bln_topo.pcsf"
    if topo_source.exists():
        topo_target = pcsf_to_pcsm(
            topo_source, OUTPUT / "occam2d_with_topo.pcsm"
        )
        restored_topo = read_pcsm(topo_target)
        has_lonlat = restored_topo.stations.lon is not None
        size = topo_target.stat().st_size
        print(
            f"{topo_target.relative_to(ROOT)}: {restored_topo.kind}, "
            f"{size:,} bytes, "
            f"stations/lon+lat present={has_lonlat} "
            f"({len(restored_topo.stations.name)} stations)"
        )

    modem_full = OUTPUT / "modem3d.pcsm"
    pcsf_to_pcsm(PCSF_OUTPUT / "modem3d_no_topo.pcsf", modem_full)
    _write_excerpt(
        modem_full,
        OUTPUT / "modem3d_excerpt.pcsm",
        (
            "X_COORDS",
            "Y_COORDS",
            "Z_COORDS",
            "X_NODES",
            "Y_NODES",
            "Z_NODES",
            "RESISTIVITY",
            "RESISTIVITY_NATIVE",
            "STATIONS",
        ),
        "python examples/pcsm_conversion/run_demo.py --full-modem",
    )
    if full_modem:
        restored = read_pcsm(modem_full)
        size = modem_full.stat().st_size
        name = modem_full.relative_to(ROOT)
        print(f"{name}: {restored.kind}, {size:,} bytes")
        pcsf_to_pcsm(
            PCSF_OUTPUT / "modem3d_no_topo.pcsf", OUTPUT / "modem3d.pcsm.gz"
        )
    else:
        modem_full.unlink()
    excerpt = OUTPUT / "modem3d_excerpt.pcsm"
    size = excerpt.stat().st_size
    print(f"{excerpt.relative_to(ROOT)}: documented excerpt, {size:,} bytes")

    # Standalone AI/DL adapter (pycsamt.format.adapters.generic) -- run
    # examples/ai_inversion_pcsf/run_demo.py first to produce these two
    # source .pcsf files.
    unet_source = AI_PCSF_OUTPUT / "toy_unet_grid2d.pcsf"
    if unet_source.exists():
        unet_full = OUTPUT / "toy_unet_grid2d.pcsm"
        pcsf_to_pcsm(unet_source, unet_full)
        unet_excerpt = OUTPUT / "toy_unet_grid2d_excerpt.pcsm"
        _write_excerpt(
            unet_full,
            unet_excerpt,
            ("X_COORDS", "Z_COORDS", "RESISTIVITY", "RESISTIVITY_NATIVE", "STATIONS"),
            "python examples/ai_inversion_pcsf/run_demo.py",
            raw_keywords=("METADATA_JSON",),
            keep=1,
        )
        unet_full.unlink()
        print(
            f"{unet_excerpt.relative_to(ROOT)}: documented excerpt, "
            f"{unet_excerpt.stat().st_size:,} bytes"
        )

    gcn_source = AI_PCSF_OUTPUT / "toy_gcn_mesh.pcsf"
    if gcn_source.exists():
        gcn_full = OUTPUT / "toy_gcn_mesh.pcsm"
        pcsf_to_pcsm(gcn_source, gcn_full)
        gcn_excerpt = OUTPUT / "toy_gcn_mesh_excerpt.pcsm"
        _write_excerpt(
            gcn_full,
            gcn_excerpt,
            ("NODES", "CONNECTIVITY", "REGION_IDS", "RESISTIVITY", "RESISTIVITY_BY_NODE"),
            "python examples/ai_inversion_pcsf/run_demo.py",
            raw_keywords=("METADATA_JSON",),
        )
        gcn_full.unlink()
        print(
            f"{gcn_excerpt.relative_to(ROOT)}: documented excerpt, "
            f"{gcn_excerpt.stat().st_size:,} bytes"
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--full-modem",
        action="store_true",
        help="retain the full plain and gzip-compressed ModEM PCSM files",
    )
    args = parser.parse_args()
    main(full_modem=args.full_modem)
