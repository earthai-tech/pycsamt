"""Execute and verify the MARE2DEM preparation tutorial's pycon blocks.

Mirrors ``generate_tutorial_modem.py``: parses every
``.. code-block:: pycon`` region out of
``docs/source/tutorials/prepare_mare2dem_inversion.rst``, in document order,
and executes the examples in one shared namespace -- exactly like a single
pasted REPL session. Any ``print(...)`` output is compared byte-exact
against what the page shows. Because the page's own plotting code runs for
real, the same call produces the actual PNGs referenced by the page's
``.. figure::`` directives; this script copies them from the run's
``figures/`` folder into
``docs/source/images/tutorials/prepare_mare2dem_inversion/``.
"""

from __future__ import annotations

import doctest
import io
import re
import shutil
import sys
from contextlib import redirect_stdout
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

ROOT = Path(__file__).resolve().parents[2]
RST_PATH = ROOT / "docs" / "source" / "tutorials" / "prepare_mare2dem_inversion.rst"
IMAGE_DIR = (
    ROOT / "docs" / "source" / "images" / "tutorials" / "prepare_mare2dem_inversion"
)

_BLOCK_RE = re.compile(
    r"^\.\. code-block:: pycon\n(?:^ {3}:linenos:\n)?\n(?P<body>(?:^ {3}.*\n|\n)+)",
    re.MULTILINE,
)

_FIGURE_COPIES = {
    "receiver_profile.png": "mare2dem_receiver_profile.png",
    "poly_mesh.png": "mare2dem_poly_mesh.png",
}


def _extract_blocks(text: str) -> list[str]:
    blocks = []
    for m in _BLOCK_RE.finditer(text):
        body = m.group("body")
        lines = [ln[3:] if ln.startswith("   ") else ln for ln in body.splitlines()]
        blocks.append("\n".join(lines))
    return blocks


def main() -> int:
    sys.path.insert(0, str(ROOT))
    text = RST_PATH.read_text(encoding="utf-8")
    blocks = _extract_blocks(text)
    print(f"Found {len(blocks)} pycon blocks")

    parser = doctest.DocTestParser()
    namespace: dict = {"__name__": "__main__"}
    failures = 0
    total_examples = 0

    import os

    os.chdir(ROOT)

    for block_idx, block in enumerate(blocks):
        examples = parser.get_examples(block)
        for example in examples:
            total_examples += 1
            buf = io.StringIO()
            try:
                with redirect_stdout(buf):
                    exec(compile(example.source, "<tutorial>", "single"), namespace)
            except Exception as exc:  # noqa: BLE001
                failures += 1
                print(f"[block {block_idx}] EXCEPTION for: {example.source.strip()!r}")
                print(f"    -> {type(exc).__name__}: {exc}")
                continue
            got = buf.getvalue()
            if example.want and got != example.want:
                failures += 1
                print(f"[block {block_idx}] MISMATCH for: {example.source.strip()!r}")
                print(f"    expected: {example.want!r}")
                print(f"    got     : {got!r}")

    print(f"Executed {total_examples} examples, {failures} mismatches/exceptions")

    figures_dir = namespace["figure_dir"]
    IMAGE_DIR.mkdir(parents=True, exist_ok=True)
    for src_name, dst_name in _FIGURE_COPIES.items():
        src = Path(figures_dir) / src_name
        if not src.exists():
            print(f"MISSING FIGURE: {src}")
            failures += 1
            continue
        shutil.copyfile(src, IMAGE_DIR / dst_name)
        print(f"copied {src} -> {IMAGE_DIR / dst_name}")

    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
