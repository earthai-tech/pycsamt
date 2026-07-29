from __future__ import annotations

import importlib
import inspect
from pathlib import Path


def test_all_public_plot_functions_accept_axes_parameter():
    root = Path(__file__).resolve().parents[1]
    missing = []

    for path in sorted(root.glob("*.py")):
        if path.name.startswith("_") or path.name == "__init__.py":
            continue

        module_name = f"pycsamt.emtools.{path.stem}"
        module = importlib.import_module(module_name)

        for name, obj in vars(module).items():
            if not (
                name.startswith("plot")
                and inspect.isfunction(obj)
                and obj.__module__ == module_name
            ):
                continue

            params = inspect.signature(obj).parameters
            if not any(
                param == "ax" or param.startswith("ax_") or param in {"axes", "axs"}
                for param in params
            ):
                missing.append(f"{module_name}.{name}")

    assert missing == []
