"""Round-trip real EDI stations through EMTF-XML and back to EDI.

Exercises both public entry points introduced for XML support:

* the ``Site``/``Sites`` wrapper API (``Sites.write_xml``, ``to_sites``
  discovering a directory of ``*.xml`` files, ``Sites.write``), and
* the underlying ``pycsamt.emtf.document.EMTF`` document directly
  (``EMTF.from_edi``, ``.write_xml``, ``EMTF.from_xml``,
  ``pycsamt.emtf.converters.edi.write_edi``),

and checks that both agree with each other and with the original EDI
source numerically (periods, impedance, tipper, variance).

Run from the repository root::

    python examples/emtf_xml_roundtrip_demo/run_demo.py
"""

from __future__ import annotations

import argparse
import json
import warnings
from pathlib import Path

import numpy as np

from pycsamt.emtf.converters.edi import edi_to_emtf, write_edi
from pycsamt.emtf.document import EMTF
from pycsamt.seg.edi import EDIFile
from pycsamt.site.base import Site, Sites, to_sites

ROOT = Path(__file__).resolve().parents[2]
EDI_DIRECTORY = ROOT / "data" / "gv_data" / "gv_final_edi"
DEFAULT_OUTPUT = Path(__file__).resolve().parent / "output"
DEFAULT_STATIONS = ["gv100", "gv101", "gv102"]


def _arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Round-trip real EDI stations through EMTF-XML and back, "
            "via both the Site/Sites API and the raw EMTF document API."
        )
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help="Demo output root (default: examples/emtf_xml_roundtrip_demo/output).",
    )
    parser.add_argument(
        "--station",
        action="append",
        default=None,
        help=(
            "EDI stem (from data/gv_data/gv_final_edi/) to include; "
            f"repeat for several. Default: {', '.join(DEFAULT_STATIONS)}."
        ),
    )
    return parser.parse_args()


def _sources(selected: list[str] | None) -> list[Path]:
    stems = selected or DEFAULT_STATIONS
    paths = [EDI_DIRECTORY / f"{stem}.edi" for stem in stems]
    missing = [p for p in paths if not p.exists()]
    if missing:
        raise FileNotFoundError(f"Missing EDI source(s): {missing}")
    return paths


# ---------------------------------------------------------------------------
# Numeric comparison
# ---------------------------------------------------------------------------


def _by_period(tf: EMTF) -> np.ndarray:
    return np.argsort(tf.periods)


def _max_abs(a: np.ndarray | None, b: np.ndarray | None) -> float | None:
    if a is None or b is None:
        return None
    return float(np.nanmax(np.abs(np.asarray(a) - np.asarray(b))))


def _compare(label: str, tf_ref: EMTF, tf_other: EMTF) -> dict:
    order_ref = _by_period(tf_ref)
    order_other = _by_period(tf_other)
    report = {
        "label": label,
        "n_periods_ref": int(tf_ref.periods.size),
        "n_periods_other": int(tf_other.periods.size),
        "max_abs_diff_periods": _max_abs(
            tf_ref.periods[order_ref], tf_other.periods[order_other]
        ),
        "max_abs_diff_z": _max_abs(
            tf_ref.z[order_ref] if tf_ref.z is not None else None,
            tf_other.z[order_other] if tf_other.z is not None else None,
        ),
    }
    if tf_ref.tipper is not None and tf_other.tipper is not None:
        report["max_abs_diff_tipper"] = _max_abs(
            tf_ref.tipper[order_ref], tf_other.tipper[order_other]
        )
    try:
        v_ref = tf_ref.get_transfer_function("impedance").estimates[
            "variance"
        ].data
        v_other = tf_other.get_transfer_function("impedance").estimates[
            "variance"
        ].data
        report["max_abs_diff_variance"] = _max_abs(
            v_ref[order_ref], v_other[order_other]
        )
    except Exception:
        report["max_abs_diff_variance"] = None
    return report


# ---------------------------------------------------------------------------
# Path A: via Site / Sites
# ---------------------------------------------------------------------------


def run_via_sites(
    sources: list[Path], output: Path
) -> tuple[Sites, Sites, list[str], list[Path]]:
    """EDI -> XML -> EDI entirely through the Site/Sites wrapper API."""

    print("\n=== Path A: Site/Sites API ===")
    sites = Sites([EDIFile(p) for p in sources])
    print(f"Loaded {len(sites)} station(s) via Sites([EDIFile(...), ...])")

    xml_dir = output / "xml"
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        xml_paths = sites.write_xml(xml_dir)
    warn_msgs = [str(w.message) for w in caught]
    print(f"Sites.write_xml -> {len(xml_paths)} file(s) in {xml_dir}")
    for p in xml_paths:
        print(f"  wrote {p.relative_to(ROOT)}")

    # Discover the *.xml directory back into a fresh Sites via the new
    # XML-aware to_sites() coercion -- exactly what ensure_sites() now
    # does automatically for any XML (or mixed EDI+XML) directory.
    sites_from_xml = to_sites(xml_dir)
    print(
        f"to_sites({xml_dir.name}/) rediscovered "
        f"{len(sites_from_xml)} XML-native site(s), "
        f"backend={sorted({s.backend for s in sites_from_xml})}"
    )

    edi_roundtrip_dir = output / "edi_roundtrip_via_site"
    edi_paths = sites_from_xml.write(edi_roundtrip_dir, exist_ok=True)
    print(f"Sites.write -> {len(edi_paths)} file(s) in {edi_roundtrip_dir}")
    for p in edi_paths:
        print(f"  wrote {p.relative_to(ROOT)}")

    return sites, sites_from_xml, warn_msgs, edi_paths


# ---------------------------------------------------------------------------
# Path B: via the raw EMTF document API
# ---------------------------------------------------------------------------


def run_via_emtf(
    sources: list[Path], output: Path
) -> tuple[list[EMTF], list[EMTF], list[str], list[Path]]:
    """EDI -> XML -> EDI directly through pycsamt.emtf.document.EMTF."""

    print("\n=== Path B: raw EMTF document API ===")
    xml_dir = output / "xml_via_emtf"
    xml_dir.mkdir(parents=True, exist_ok=True)
    edi_roundtrip_dir = output / "edi_roundtrip_via_emtf"
    edi_roundtrip_dir.mkdir(parents=True, exist_ok=True)

    docs_original: list[EMTF] = []
    docs_roundtrip: list[EMTF] = []
    warn_msgs: list[str] = []
    edi_paths: list[Path] = []
    for src in sources:
        edi = EDIFile(src)
        tf = edi_to_emtf(edi)
        docs_original.append(tf)

        xml_path = xml_dir / f"{src.stem}.xml"
        tf.write_xml(xml_path)
        print(f"EMTF.write_xml -> {xml_path.relative_to(ROOT)}")

        tf2 = EMTF.from_xml(xml_path, strict=False)
        docs_roundtrip.append(tf2)

        edi_path = edi_roundtrip_dir / f"{src.stem}.edi"
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            write_edi(tf2, edi_path, on_loss="warn")
        warn_msgs.extend(str(w.message) for w in caught)
        print(f"write_edi(EMTF, ...) -> {edi_path.relative_to(ROOT)}")
        edi_paths.append(edi_path)

    return docs_original, docs_roundtrip, warn_msgs, edi_paths


def _reread_from_disk(path: Path) -> EMTF:
    """Parse a *written* EDI file back off disk (not the in-memory
    object that wrote it) -- this is the check that actually catches a
    writer silently emitting a placeholder instead of real content."""

    return edi_to_emtf(EDIFile(path))


def main() -> None:
    args = _arguments()
    output = args.output.expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    sources = _sources(args.station)
    print(f"Source EDI files ({len(sources)}):")
    for p in sources:
        print(f"  {p.relative_to(ROOT)}")

    sites_original, sites_roundtrip, warn_a, edi_paths_a = run_via_sites(
        sources, output
    )
    docs_original_b, docs_roundtrip_b, warn_b, edi_paths_b = run_via_emtf(
        sources, output
    )

    print("\n=== Numeric verification ===")
    per_station = []
    for src, site_orig, site_rt, tf_orig_b, tf_rt_b, edi_path_a, edi_path_b in zip(
        sources,
        sites_original,
        sites_roundtrip,
        docs_original_b,
        docs_roundtrip_b,
        edi_paths_a,
        edi_paths_b,
    ):
        tf_orig_a = site_orig.tf
        tf_rt_a = site_rt.tf
        # Re-parse the files actually written to disk -- comparing only
        # in-memory objects would miss a writer silently emitting a
        # placeholder instead of real EDI content (a real bug this demo
        # caught in Sites.write(), since fixed).
        tf_from_disk_a = _reread_from_disk(edi_path_a)
        tf_from_disk_b = _reread_from_disk(edi_path_b)
        record = {
            "station": src.stem,
            "source": str(src.relative_to(ROOT)),
            "via_site_roundtrip_vs_original": _compare(
                "Site: EDI -> XML -> EDI vs. original", tf_orig_a, tf_rt_a
            ),
            "via_emtf_roundtrip_vs_original": _compare(
                "EMTF: EDI -> XML -> EDI vs. original", tf_orig_b, tf_rt_b
            ),
            "site_path_vs_emtf_path_agreement": _compare(
                "Site path vs. raw EMTF path (roundtripped)", tf_rt_a, tf_rt_b
            ),
            "on_disk_edi_via_site_vs_original": _compare(
                "on-disk written EDI (Site path) vs. original",
                site_orig.tf,
                tf_from_disk_a,
            ),
            "on_disk_edi_via_emtf_vs_original": _compare(
                "on-disk written EDI (EMTF path) vs. original",
                tf_orig_b,
                tf_from_disk_b,
            ),
        }
        per_station.append(record)
        print(f"\n{src.stem}:")
        for key in (
            "via_site_roundtrip_vs_original",
            "via_emtf_roundtrip_vs_original",
            "site_path_vs_emtf_path_agreement",
            "on_disk_edi_via_site_vs_original",
            "on_disk_edi_via_emtf_vs_original",
        ):
            c = record[key]
            print(
                f"  {c['label']}: "
                f"max|dZ|={c['max_abs_diff_z']!r}  "
                f"max|dT|={c.get('max_abs_diff_tipper')!r}  "
                f"max|dVar|={c.get('max_abs_diff_variance')!r}"
            )

    summary = {
        "schema": "pycsamt.emtf_xml_roundtrip.demo/v1",
        "source_directory": str(EDI_DIRECTORY.relative_to(ROOT)),
        "stations": per_station,
        "data_loss_warnings_via_site": warn_a,
        "data_loss_warnings_via_emtf": warn_b,
    }
    summary_path = output / "demo-summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf8"
    )
    print(f"\nDemo completed: {summary_path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
