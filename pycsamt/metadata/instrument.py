# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.metadata.instrument
============================

Structured acquisition-system descriptor for MT / AMT / CSAMT / TEM surveys.

:class:`InstrumentMeta` captures the recording system and sensor configuration
in a form that can be round-tripped to/from EDI ``>HEAD`` fields and JSON/YAML.
Pre-built presets for common field systems (Phoenix, Metronix, LEMI, Zonge)
are available via :func:`known_system`.

Quick start
-----------
::

    from pycsamt.metadata.instrument import InstrumentMeta, SensorSpec

    # Manually construct
    inst = InstrumentMeta(
        system="Phoenix V8",
        serial="V8-20473",
        magnetic_sensor=SensorSpec(
            sensor_type="induction_coil",
            model="MTC-150H",
            frequency_range=(0.00001, 10000.0),
        ),
        electric_sensor=SensorSpec(sensor_type="electrode", model="Pb-PbCl2"),
    )
    print(inst.summary())

    # Build from a preset
    inst2 = InstrumentMeta.from_preset("metronix_adu07")

    # Round-trip JSON
    import json

    s = inst.to_json()
    inst3 = InstrumentMeta.from_json(s)

    # Push fields to an EDI Head dict
    head_dict = (
        inst.to_head_fields()
    )  # {"acqby": "Phoenix V8 / V8-20473", ...}

API Reference
-------------
.. autosummary::
   :toctree:

   SensorSpec
   InstrumentMeta
   known_system
   list_presets
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from typing import Any, Literal

try:
    import yaml as _yaml

    _HAS_YAML = True
except ModuleNotFoundError:
    _HAS_YAML = False

__all__ = [
    "SensorSpec",
    "InstrumentMeta",
    "KNOWN_SYSTEMS",
    "known_system",
    "list_presets",
]

# ---------------------------------------------------------------------------
# SensorType alias
# ---------------------------------------------------------------------------
SensorType = Literal["induction_coil", "fluxgate", "electrode", "other"]


# ---------------------------------------------------------------------------
# SensorSpec
# ---------------------------------------------------------------------------
@dataclass
class SensorSpec:
    """Specification for a single sensor (magnetic or electric).

    Parameters
    ----------
    sensor_type : {"induction_coil", "fluxgate", "electrode", "other"}
        Physical transducer principle.
    model : str
        Manufacturer model string, e.g. ``"MTC-150H"``.
    frequency_range : tuple[float, float] or None
        (f_min, f_max) in Hz that the sensor is rated for.
    notes : str
        Free-text annotation (calibration file name, mounting details, etc.).

    Examples
    --------
    >>> s = SensorSpec("induction_coil", "MTC-150H", (1e-5, 1e4))
    >>> s.covers(0.1)
    True
    >>> s.covers(1e6)
    False
    """

    sensor_type: SensorType = "induction_coil"
    model: str = ""
    frequency_range: tuple[float, float] | None = None
    notes: str = ""

    # ------------------------------------------------------------------
    def covers(self, frequency_hz: float) -> bool:
        """Return True if *frequency_hz* lies within the rated band."""
        if self.frequency_range is None:
            return True  # assume universal if unspecified
        f_min, f_max = self.frequency_range
        return f_min <= frequency_hz <= f_max

    # ------------------------------------------------------------------
    def __str__(self) -> str:
        fr = (
            f"{self.frequency_range[0]:.3g}–{self.frequency_range[1]:.3g} Hz"
            if self.frequency_range
            else "full band"
        )
        label = self.model if self.model else self.sensor_type
        return f"<SensorSpec {label!r} [{fr}]>"

    __repr__ = __str__

    # ------------------------------------------------------------------
    # (de)serialisation helpers used by InstrumentMeta
    # ------------------------------------------------------------------
    def _to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        if d["frequency_range"] is not None:
            d["frequency_range"] = list(d["frequency_range"])
        return d

    @classmethod
    def _from_dict(cls, d: dict[str, Any]) -> SensorSpec:
        fr = d.get("frequency_range")
        return cls(
            sensor_type=d.get("sensor_type", "induction_coil"),
            model=d.get("model", ""),
            frequency_range=tuple(fr) if fr is not None else None,
            notes=d.get("notes", ""),
        )


# ---------------------------------------------------------------------------
# InstrumentMeta
# ---------------------------------------------------------------------------
@dataclass
class InstrumentMeta:
    """Structured descriptor for an EM acquisition system.

    Parameters
    ----------
    system : str
        Human-readable system name, e.g. ``"Phoenix V8"`` or
        ``"Metronix ADU-07e"``.
    serial : str or None
        Logger serial number or site-specific tag.
    magnetic_sensor : SensorSpec or None
        Magnetic-field transducer (H field) specification.
    electric_sensor : SensorSpec or None
        Electric-field transducer (E field) specification.
    software_version : str
        Processing / acquisition software version string.
    notes : str
        Free-text annotation (calibration date, operator remarks, …).

    Notes
    -----
    The class deliberately stays thin — only fields that appear (or could
    appear) in a standard EDI ``>HEAD`` block or in a compact YAML survey
    manifest are stored here.  Calibration tables belong in separate files
    referenced by :attr:`notes`.

    Examples
    --------
    >>> from pycsamt.metadata.instrument import InstrumentMeta
    >>> inst = InstrumentMeta(system="Phoenix V8", serial="V8-00123")
    >>> inst.to_head_fields()
    {'acqby': 'Phoenix V8 / V8-00123', 'progvers': 'pyCSAMT'}
    """

    system: str = ""
    serial: str | None = None
    magnetic_sensor: SensorSpec | None = None
    electric_sensor: SensorSpec | None = None
    software_version: str = ""
    notes: str = ""

    # ------------------------------------------------------------------
    # Convenience properties
    # ------------------------------------------------------------------
    @property
    def label(self) -> str:
        """Short identifier: ``"<system> / <serial>"`` or just system."""
        if self.serial:
            return f"{self.system} / {self.serial}"
        return self.system

    # ------------------------------------------------------------------
    # Summaries
    # ------------------------------------------------------------------
    def summary(self) -> str:
        """Return a multi-line human-readable summary."""
        lines = [
            f"System       : {self.system or '—'}",
            f"Serial       : {self.serial or '—'}",
        ]
        if self.software_version:
            lines.append(f"Software     : {self.software_version}")
        if self.magnetic_sensor:
            lines.append(f"Mag. sensor  : {self.magnetic_sensor}")
        if self.electric_sensor:
            lines.append(f"Elec. sensor : {self.electric_sensor}")
        if self.notes:
            lines.append(f"Notes        : {self.notes}")
        return "\n".join(lines)

    def __repr__(self) -> str:
        return (
            f"<InstrumentMeta system={self.system!r} serial={self.serial!r}>"
        )

    # ------------------------------------------------------------------
    # EDI HEAD integration
    # ------------------------------------------------------------------
    def to_head_fields(self) -> dict[str, str]:
        """Return a dict of EDI ``>HEAD`` key-value pairs.

        The mapping is::

            acqby    ← system / serial  (or just system)
            progvers ← software_version (or the package default)

        Returns
        -------
        dict
            Ready to pass to :meth:`~pycsamt.seg.heads.Head.update`.
        """
        from pycsamt import (
            __version__ as _ver,  # lazy import to avoid cycles
        )

        fields: dict[str, str] = {}
        if self.system:
            fields["acqby"] = self.label
        fields["progvers"] = self.software_version or f"pyCSAMT {_ver}"
        return fields

    @classmethod
    def from_head(cls, head: Any) -> InstrumentMeta:
        """Construct from a :class:`~pycsamt.seg.heads.Head` object.

        Only ``acqby`` and ``progvers`` are extracted; sensor details
        cannot be inferred from the EDI header alone.

        Parameters
        ----------
        head :
            A :class:`pycsamt.seg.heads.Head` instance (or any object with
            ``acqby`` and ``progvers`` string attributes).

        Returns
        -------
        InstrumentMeta
        """
        acqby = getattr(head, "acqby", None) or ""
        progvers = getattr(head, "progvers", None) or ""

        # Split "System / Serial" heuristic
        parts = [p.strip() for p in acqby.split("/", 1)]
        system = parts[0]
        serial = parts[1] if len(parts) == 2 else None

        return cls(
            system=system,
            serial=serial,
            software_version=progvers,
        )

    # ------------------------------------------------------------------
    # Presets
    # ------------------------------------------------------------------
    @classmethod
    def from_preset(cls, name: str) -> InstrumentMeta:
        """Build from a named preset in :data:`KNOWN_SYSTEMS`.

        Parameters
        ----------
        name : str
            Case-insensitive preset key, e.g. ``"phoenix_v8"``.

        Returns
        -------
        InstrumentMeta

        Raises
        ------
        KeyError
            If *name* is not in :data:`KNOWN_SYSTEMS`.

        Examples
        --------
        >>> inst = InstrumentMeta.from_preset("lemi_424")
        >>> inst.system
        'LEMI-424'
        """
        key = name.lower().replace("-", "_").replace(" ", "_")
        if key not in KNOWN_SYSTEMS:
            available = ", ".join(sorted(KNOWN_SYSTEMS))
            raise KeyError(f"Unknown preset {name!r}. Available: {available}")
        d = KNOWN_SYSTEMS[key]
        return cls._from_dict(d)

    # ------------------------------------------------------------------
    # Serialisation
    # ------------------------------------------------------------------
    def to_dict(self) -> dict[str, Any]:
        """Serialise to a plain dict (JSON-safe)."""
        d: dict[str, Any] = {
            "system": self.system,
            "serial": self.serial,
            "software_version": self.software_version,
            "notes": self.notes,
        }
        if self.magnetic_sensor is not None:
            d["magnetic_sensor"] = self.magnetic_sensor._to_dict()
        if self.electric_sensor is not None:
            d["electric_sensor"] = self.electric_sensor._to_dict()
        return d

    @classmethod
    def _from_dict(cls, d: dict[str, Any]) -> InstrumentMeta:
        mag = d.get("magnetic_sensor")
        elec = d.get("electric_sensor")
        return cls(
            system=d.get("system", ""),
            serial=d.get("serial"),
            magnetic_sensor=SensorSpec._from_dict(mag) if mag else None,
            electric_sensor=SensorSpec._from_dict(elec) if elec else None,
            software_version=d.get("software_version", ""),
            notes=d.get("notes", ""),
        )

    def to_json(self, indent: int = 2) -> str:
        """Serialise to a JSON string."""
        return json.dumps(self.to_dict(), indent=indent)

    @classmethod
    def from_json(cls, s: str) -> InstrumentMeta:
        """Deserialise from a JSON string."""
        return cls._from_dict(json.loads(s))

    def to_yaml(self) -> str:
        """Serialise to a YAML string (requires *PyYAML*)."""
        if not _HAS_YAML:
            raise ImportError(
                "PyYAML is required for YAML serialisation: pip install pyyaml"
            )
        return _yaml.dump(self.to_dict(), default_flow_style=False)

    @classmethod
    def from_yaml(cls, s: str) -> InstrumentMeta:
        """Deserialise from a YAML string (requires *PyYAML*)."""
        if not _HAS_YAML:
            raise ImportError(
                "PyYAML is required for YAML deserialisation: pip install pyyaml"
            )
        return cls._from_dict(_yaml.safe_load(s))

    def save(self, path: str, fmt: str = "json") -> None:
        """Write to *path* as JSON or YAML.

        Parameters
        ----------
        path : str
            Destination file path.
        fmt : {"json", "yaml"}
        """
        from pathlib import Path

        p = Path(path)
        if fmt == "json":
            p.write_text(self.to_json(), encoding="utf-8")
        elif fmt in ("yaml", "yml"):
            p.write_text(self.to_yaml(), encoding="utf-8")
        else:
            raise ValueError(
                f"Unknown format {fmt!r}; choose 'json' or 'yaml'."
            )

    @classmethod
    def load(cls, path: str) -> InstrumentMeta:
        """Load from a JSON or YAML file (format detected by extension)."""
        from pathlib import Path

        p = Path(path)
        text = p.read_text(encoding="utf-8")
        if p.suffix.lower() in (".yaml", ".yml"):
            return cls.from_yaml(text)
        return cls.from_json(text)


# ---------------------------------------------------------------------------
# Pre-built system presets
# ---------------------------------------------------------------------------
#: Registry of known acquisition systems.
#: Each value is a dict compatible with :meth:`InstrumentMeta._from_dict`.
KNOWN_SYSTEMS: dict[str, dict[str, Any]] = {
    # ── Phoenix Geophysics ──────────────────────────────────────────────
    "phoenix_v8": {
        "system": "Phoenix V8",
        "serial": None,
        "software_version": "SSMT2000",
        "magnetic_sensor": {
            "sensor_type": "induction_coil",
            "model": "MTC-150H",
            "frequency_range": [1e-4, 2e3],
            "notes": "Broadband MT coil; low-noise below 1 Hz",
        },
        "electric_sensor": {
            "sensor_type": "electrode",
            "model": "Pb-PbCl2 porous pot",
            "frequency_range": None,
            "notes": "Non-polarising; standard 50–100 m dipole",
        },
        "notes": "Standard Phoenix V8 MT/AMT system",
    },
    "phoenix_mtx": {
        "system": "Phoenix MTX",
        "serial": None,
        "software_version": "SSMT2000",
        "magnetic_sensor": {
            "sensor_type": "induction_coil",
            "model": "MTC-155",
            "frequency_range": [1e-5, 1e4],
            "notes": "Ultra-broadband Phoenix coil",
        },
        "electric_sensor": {
            "sensor_type": "electrode",
            "model": "Pb-PbCl2 porous pot",
            "frequency_range": None,
            "notes": "",
        },
        "notes": "Phoenix MTX broadband system",
    },
    # ── Metronix ────────────────────────────────────────────────────────
    "metronix_adu07": {
        "system": "Metronix ADU-07e",
        "serial": None,
        "software_version": "ADU-07e firmware",
        "magnetic_sensor": {
            "sensor_type": "induction_coil",
            "model": "MFS-06e",
            "frequency_range": [1e-4, 4096.0],
            "notes": "Metronix broadband coil; pairs with ADU-07e",
        },
        "electric_sensor": {
            "sensor_type": "electrode",
            "model": "Pb-PbCl2 porous pot",
            "frequency_range": None,
            "notes": "",
        },
        "notes": "Metronix ADU-07e broadband MT system",
    },
    "metronix_adu08": {
        "system": "Metronix ADU-08e",
        "serial": None,
        "software_version": "ADU-08e firmware",
        "magnetic_sensor": {
            "sensor_type": "induction_coil",
            "model": "MFS-07e",
            "frequency_range": [1e-5, 4096.0],
            "notes": "High-resolution Metronix coil for quiet sites",
        },
        "electric_sensor": {
            "sensor_type": "electrode",
            "model": "Pb-PbCl2 porous pot",
            "frequency_range": None,
            "notes": "",
        },
        "notes": "Metronix ADU-08e high-resolution system",
    },
    # ── LEMI ────────────────────────────────────────────────────────────
    "lemi_424": {
        "system": "LEMI-424",
        "serial": None,
        "software_version": "",
        "magnetic_sensor": {
            "sensor_type": "induction_coil",
            "model": "LEMI-120",
            "frequency_range": [3e-4, 1e3],
            "notes": "Low-frequency induction coil for long-period MT",
        },
        "electric_sensor": {
            "sensor_type": "electrode",
            "model": "Pb-PbCl2",
            "frequency_range": None,
            "notes": "",
        },
        "notes": "LEMI-424 low-frequency broadband system",
    },
    # ── Zonge ───────────────────────────────────────────────────────────
    "zonge_gdp32": {
        "system": "Zonge GDP-32",
        "serial": None,
        "software_version": "",
        "magnetic_sensor": {
            "sensor_type": "induction_coil",
            "model": "ANT/6",
            "frequency_range": [1e-3, 1e4],
            "notes": "Zonge broadband coil, common in CSAMT surveys",
        },
        "electric_sensor": {
            "sensor_type": "electrode",
            "model": "Pb-PbCl2",
            "frequency_range": None,
            "notes": "",
        },
        "notes": "Zonge GDP-32 CSAMT/MT system",
    },
    # ── Geometrics / OYO ────────────────────────────────────────────────
    "geometrics_stratagem": {
        "system": "Geometrics Stratagem EH4",
        "serial": None,
        "software_version": "",
        "magnetic_sensor": {
            "sensor_type": "induction_coil",
            "model": "EH4 coil",
            "frequency_range": [10.0, 1e5],
            "notes": "AMT band; controlled source / natural field",
        },
        "electric_sensor": {
            "sensor_type": "electrode",
            "model": "Pb-PbCl2",
            "frequency_range": None,
            "notes": "",
        },
        "notes": "EH4 tensor AMT system (10 Hz – 100 kHz)",
    },
    # ── Generic fluxgate (TEM / airborne) ───────────────────────────────
    "generic_fluxgate": {
        "system": "Generic fluxgate magnetometer",
        "serial": None,
        "software_version": "",
        "magnetic_sensor": {
            "sensor_type": "fluxgate",
            "model": "",
            "frequency_range": [0.0, 10.0],
            "notes": "Low-frequency or DC fluxgate; used in TEM & airborne",
        },
        "electric_sensor": None,
        "notes": "Generic fluxgate configuration",
    },
}


# ---------------------------------------------------------------------------
# Module-level helpers
# ---------------------------------------------------------------------------
def list_presets() -> list:
    """Return sorted list of available preset keys.

    Examples
    --------
    >>> from pycsamt.metadata.instrument import list_presets
    >>> list_presets()
    ['geometrics_stratagem', 'generic_fluxgate', 'lemi_424', ...]
    """
    return sorted(KNOWN_SYSTEMS.keys())


def known_system(name: str) -> InstrumentMeta:
    """Shortcut for :meth:`InstrumentMeta.from_preset`.

    Parameters
    ----------
    name : str
        Preset key (case-insensitive, spaces/hyphens treated as underscores).

    Returns
    -------
    InstrumentMeta

    Examples
    --------
    >>> from pycsamt.metadata.instrument import known_system
    >>> inst = known_system("phoenix_v8")
    >>> inst.system
    'Phoenix V8'
    """
    return InstrumentMeta.from_preset(name)
