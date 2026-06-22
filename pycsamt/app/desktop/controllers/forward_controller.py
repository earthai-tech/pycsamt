# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
ForwardController — model library persistence and geological preset helpers.

The library is stored as JSON at ``~/.pycsamt/forward_models.json`` and
contains a flat list of serialised 1-D LayeredModel records (resistivity +
thickness + name + dim).  2-D / 3-D parameter sets are stored as raw dicts.
"""
from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any

import numpy as np

logger = logging.getLogger(__name__)

_LIBRARY_PATH = Path.home() / ".pycsamt" / "forward_models.json"

# All 13 geology priors available in pycsamt.forward
GEOLOGY_PRESET_NAMES: list[str] = [
    "sedimentary",
    "crystalline",
    "geothermal",
    "marine",
    "permafrost",
    "basement",
    "coastal",
    "evaporite",
    "hydrothermal",
    "laterite",
    "mineralized",
    "porphyry",
    "volcanic",
]


class ForwardController:
    """Manages the forward model library and geological preset generation."""

    def __init__(self) -> None:
        self._library: list[dict] = []
        self.load_library()

    # ── Library persistence ───────────────────────────────────────────────────

    def load_library(self) -> None:
        if _LIBRARY_PATH.exists():
            try:
                data = json.loads(_LIBRARY_PATH.read_text(encoding="utf-8"))
                self._library = data.get("models", [])
            except Exception:
                logger.warning("Could not load forward model library; starting empty.")
                self._library = []

    def save_library(self) -> None:
        try:
            _LIBRARY_PATH.parent.mkdir(parents=True, exist_ok=True)
            _LIBRARY_PATH.write_text(
                json.dumps({"models": self._library}, indent=2),
                encoding="utf-8",
            )
        except Exception:
            logger.warning("Could not save forward model library.")

    # ── Library CRUD ──────────────────────────────────────────────────────────

    @property
    def model_names(self) -> list[str]:
        return [m["name"] for m in self._library]

    def save_model(self, name: str, record: dict[str, Any]) -> None:
        """Add or overwrite a model record by name."""
        record = dict(record)
        record["name"] = name
        existing = next(
            (i for i, m in enumerate(self._library) if m["name"] == name), None
        )
        if existing is not None:
            self._library[existing] = record
        else:
            self._library.append(record)
        self.save_library()

    def delete_model(self, name: str) -> bool:
        before = len(self._library)
        self._library = [m for m in self._library if m["name"] != name]
        if len(self._library) < before:
            self.save_library()
            return True
        return False

    def rename_model(self, old_name: str, new_name: str) -> bool:
        for m in self._library:
            if m["name"] == old_name:
                m["name"] = new_name
                self.save_library()
                return True
        return False

    def get_model_record(self, name: str) -> dict | None:
        return next((m for m in self._library if m["name"] == name), None)

    # ── Serialisation helpers ─────────────────────────────────────────────────

    @staticmethod
    def model_to_record(
        name: str,
        dim: str,
        resistivity: np.ndarray,
        thickness: np.ndarray,
    ) -> dict:
        """Serialise a 1-D LayeredModel to a plain dict."""
        return {
            "name": name,
            "dim": dim,
            "resistivity": resistivity.tolist(),
            "thickness": thickness.tolist(),
        }

    @staticmethod
    def record_to_arrays(record: dict) -> tuple[np.ndarray, np.ndarray]:
        """Return (resistivity, thickness) arrays from a library record."""
        return (
            np.asarray(record["resistivity"], dtype=float),
            np.asarray(record["thickness"],   dtype=float),
        )

    # ── Geological presets ────────────────────────────────────────────────────

    @staticmethod
    def build_preset_1d(name: str, seed: int | None = None) -> dict:
        """
        Generate a random 1-D model from the named geological prior.

        Returns a dict with keys ``resistivity``, ``thickness``.
        """
        from pycsamt.forward import LayeredModel
        model = LayeredModel.from_geology(name, seed=seed)
        return {
            "resistivity": model.resistivity.tolist(),
            "thickness":   model.thickness.tolist(),
        }
