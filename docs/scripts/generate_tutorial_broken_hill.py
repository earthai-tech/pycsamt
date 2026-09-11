"""Execute and verify the Broken Hill MT tutorial's pycon blocks.

Mirrors ``generate_tutorial_modem.py``: parses every
``.. code-block:: pycon`` region out of
``docs/source/tutorials/model_broken_hill_mt_3d.rst``, in document order,
and executes the examples in one shared namespace -- exactly like a
single pasted REPL session. Any ``print(...)`` output is compared
byte-exact against what the page shows. Because the page's own plotting
code runs for real, the same calls produce the actual PNGs referenced by
the page's ``.. figure::`` directives; this script copies them from the
run's ``figures/`` folder into
``docs/source/images/tutorials/model_broken_hill_mt_3d/``.
"""

from __future__ import annotations

import doctest
import io
import os
import re
import shutil
import sys
import warnings
from contextlib import redirect_stdout
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[2]
RST_PATH = (
    ROOT / "docs" / "source" / "tutorials" / "model_broken_hill_mt_3d.rst"
)
IMAGE_DIR = (
    ROOT
    / "docs"
    / "source"
    / "images"
    / "tutorials"
    / "model_broken_hill_mt_3d"
)

_BLOCK_RE = re.compile(
    r"^\.\. code-block:: pycon\n(?:^ {3}:linenos:\n)?\n"
    r"(?P<body>(?:^ {3}.*\n|\n)+)",
    re.MULTILINE,
)

_FIGURE_COPIES = {
    "phase_tensor_grid.png": "phase_tensor_grid.png",
    "misfit_map.png": "misfit_map.png",
    "depth_slices.png": "depth_slices.png",
    "conductance.png": "conductance.png",
    "section.png": "section.png",
}


def _extract_blocks(text: str) -> list[str]:
    blocks = []
    for m in _BLOCK_RE.finditer(text):
        body = m.group("body")
        lines = [
            ln[3:] if ln.startswith("   ") else ln
            for ln in body.splitlines()
        ]
        blocks.append("\n".join(lines))
    return blocks


def main() -> int:
    sys.path.insert(0, str(ROOT))
    os.chdir(ROOT)
    text = RST_PATH.read_text(encoding="utf-8")
    blocks = _extract_blocks(text)
    print(f"Found {len(blocks)} pycon blocks")

    parser = doctest.DocTestParser()
    namespace: dict = {"__name__": "__main__"}
    failures = 0
    total_examples = 0

    for block_idx, block in enumerate(blocks):
        for example in parser.get_examples(block):
            total_examples += 1
            buf = io.StringIO()
            try:
                with redirect_stdout(buf):
                    exec(  # noqa: S102 - executing verified docs code
                        compile(example.source, "<tutorial>", "single"),
                        namespace,
                    )
            except Exception as exc:  # noqa: BLE001
                failures += 1
                print(
                    f"[block {block_idx}] EXCEPTION for: "
                    f"{example.source.strip()!r}"
                )
                print(f"    -> {type(exc).__name__}: {exc}")
                continue
            got = buf.getvalue()
            if example.want and got != example.want:
                failures += 1
                print(
                    f"[block {block_idx}] MISMATCH for: "
                    f"{example.source.strip()!r}"
                )
                print(f"    expected: {example.want!r}")
                print(f"    got     : {got!r}")

    print(
        f"Executed {total_examples} examples, "
        f"{failures} mismatches/exceptions"
    )

    figures_dir = Path(namespace["figure_dir"])
    IMAGE_DIR.mkdir(parents=True, exist_ok=True)
    for src_name, dst_name in _FIGURE_COPIES.items():
        src = figures_dir / src_name
        if not src.exists():
            print(f"MISSING FIGURE: {src}")
            failures += 1
            continue
        shutil.copyfile(src, IMAGE_DIR / dst_name)
        print(f"copied {src.name} -> {IMAGE_DIR.name}/{dst_name}")

    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
