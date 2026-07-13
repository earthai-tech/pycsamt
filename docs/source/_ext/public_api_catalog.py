"""Sphinx directive for grouped public API autosummary catalogues."""

from __future__ import annotations

import importlib
import inspect
import pkgutil
from collections.abc import Iterable

from sphinx.util.docutils import SphinxDirective


def _public_names(module: object) -> Iterable[str]:
    """Return the declared public API, falling back to visible attributes."""
    declared = getattr(module, "__all__", None)
    if declared:
        return sorted(
            {str(name) for name in declared if not str(name).startswith("_")}
        )
    return sorted(name for name in vars(module) if not name.startswith("_"))


def _members(module_name: str) -> dict[str, list[str]]:
    """Classify public modules, classes, and functions for *module_name*."""
    module = importlib.import_module(module_name)
    result: dict[str, list[str]] = {
        "Modules": [],
        "Classes": [],
        "Functions": [],
    }

    package_path = getattr(module, "__path__", None)
    if package_path is not None:
        for info in pkgutil.iter_modules(package_path):
            if not info.name.startswith("_") and info.name not in {
                "tests",
                "test",
            }:
                result["Modules"].append(f"{module_name}.{info.name}")

    for name in _public_names(module):
        try:
            value = getattr(module, name)
        except (AttributeError, ImportError):
            continue
        qualified = f"{module_name}.{name}"
        if inspect.isclass(value):
            owner = getattr(value, "__module__", "")
            owner_is_public = owner.startswith("pycsamt.") and not any(
                part.startswith("_") for part in owner.split(".")[1:]
            )
            if owner_is_public:
                qualified = f"{owner}.{value.__qualname__}"
            result["Classes"].append(qualified)
        elif inspect.isroutine(value):
            owner = getattr(value, "__module__", "")
            owner_is_public = owner.startswith("pycsamt.") and not any(
                part.startswith("_") for part in owner.split(".")[1:]
            )
            if owner_is_public:
                qualified = f"{owner}.{value.__qualname__}"
            result["Functions"].append(qualified)

    for values in result.values():
        values[:] = sorted(set(values), key=str.casefold)
    return result


class PublicAPICatalogue(SphinxDirective):
    """Insert grouped autosummary tables for one public package."""

    required_arguments = 1
    has_content = False

    def run(self):
        module_name = self.arguments[0].strip()
        lines: list[str] = []
        for heading, names in _members(module_name).items():
            if not names:
                continue
            lines.extend(
                [
                    heading,
                    "~" * len(heading),
                    "",
                    ".. autosummary::",
                    "   :nosignatures:",
                    "",
                    *(f"   {name}" for name in names),
                    "",
                ]
            )
        self.state_machine.insert_input(lines, self.get_source_info()[0])
        return []


class PublicAPIIndex(SphinxDirective):
    """Insert one flat autosummary table for several public packages."""

    has_content = True

    def run(self):
        names: set[str] = set()
        for raw_name in self.content:
            module_name = raw_name.strip()
            if not module_name or module_name.startswith("#"):
                continue
            members = _members(module_name)
            names.update(members["Classes"])
            names.update(members["Functions"])

        lines = [
            ".. autosummary::",
            "   :nosignatures:",
            "",
            *(f"   {name}" for name in sorted(names, key=str.casefold)),
            "",
        ]
        self.state_machine.insert_input(lines, self.get_source_info()[0])
        return []


def setup(app):
    app.add_directive("public-api-catalogue", PublicAPICatalogue)
    app.add_directive("public-api-index", PublicAPIIndex)
    return {"parallel_read_safe": True, "parallel_write_safe": True}
