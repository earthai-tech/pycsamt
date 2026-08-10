"""Validate that every public API catalogue entry has a Sphinx target.

Run this after an HTML documentation build::

    python docs/scripts/validate_api_catalog.py

The validator reads the package list from ``docs/source/api/index.rst``, uses
the same discovery code as the catalogue directive, and checks the generated
``objects.inv``.  A nonzero exit status prevents plain, unlinked object names
from reaching the published long table.
"""

from __future__ import annotations

import argparse
import importlib.util
import os
import sys
from pathlib import Path
from typing import Any

from sphinx.util.inventory import InventoryFile


REPO_ROOT = Path(__file__).resolve().parents[2]
os.environ.setdefault("PYCSAMT_DOCS_BUILD", "1")
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
API_INDEX = REPO_ROOT / "docs" / "source" / "api" / "index.rst"
CATALOG_EXTENSION = (
    REPO_ROOT / "docs" / "source" / "_ext" / "public_api_catalog.py"
)
DEFAULT_INVENTORY = REPO_ROOT / "docs" / "build" / "html" / "objects.inv"


def _load_catalog_module() -> Any:
    spec = importlib.util.spec_from_file_location(
        "pycsamt_public_api_catalog", CATALOG_EXTENSION
    )
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load catalogue extension: {CATALOG_EXTENSION}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def catalogue_packages(path: Path = API_INDEX) -> list[str]:
    """Return package names inside the ``public-api-index`` directive."""
    packages: list[str] = []
    in_directive = False
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(".. public-api-index::"):
            in_directive = True
            continue
        if not in_directive:
            continue
        if not line.strip():
            continue
        if line.startswith("   "):
            value = line.strip()
            if value and not value.startswith("#"):
                packages.append(value)
            continue
        break
    return packages


def catalogue_objects(packages: list[str]) -> list[str]:
    """Discover the exact class and function names emitted by the table."""
    catalog = _load_catalog_module()
    names: set[str] = set()
    for package in packages:
        members = catalog._members(package)
        names.update(members["Classes"])
        names.update(members["Functions"])
    return sorted(names, key=str.casefold)


def inventory_objects(path: Path) -> set[str]:
    """Load indexed Python classes and functions from a Sphinx inventory."""
    with path.open("rb") as stream:
        inventory = InventoryFile.load(stream, "", lambda _base, uri: uri)
    return set(inventory.get("py:class", {})) | set(
        inventory.get("py:function", {})
    )


def unresolved_objects(inventory: Path = DEFAULT_INVENTORY) -> list[str]:
    """Return catalogue objects absent from *inventory*."""
    expected = catalogue_objects(catalogue_packages())
    indexed = inventory_objects(inventory)
    return [name for name in expected if name not in indexed]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--inventory",
        type=Path,
        default=DEFAULT_INVENTORY,
        help="path to the built Sphinx objects.inv",
    )
    args = parser.parse_args(argv)
    inventory = args.inventory.resolve()
    if not inventory.is_file():
        parser.error(
            f"inventory does not exist: {inventory}; build the HTML docs first"
        )

    missing = unresolved_objects(inventory)
    if missing:
        print(
            f"ERROR: {len(missing)} public API catalogue object(s) have no "
            "Sphinx target:",
            file=sys.stderr,
        )
        for name in missing:
            print(f"  - {name}", file=sys.stderr)
        return 1

    count = len(catalogue_objects(catalogue_packages()))
    print(f"OK: all {count} public API catalogue objects resolve")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
