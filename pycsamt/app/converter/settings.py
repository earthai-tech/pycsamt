# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.app.converter.settings
================================

Persisted defaults for the converter app, backed by ``QSettings``
(``earthai-tech`` / ``pycsamt-converter`` -- an INI file under the OS
user-config directory on Windows/Linux, or the native registry/plist on
Windows/macOS depending on Qt's default format).

Every value has a plain-Python default, so this module -- and every
caller of :func:`load_settings` -- works even before PySide6 is
imported anywhere else, and tests can call :func:`load_settings` /
:func:`save_settings` without a running ``QApplication`` as long as a
``QCoreApplication`` name/org has been set (Qt requires that for
``QSettings`` to resolve a file path).
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, fields


@dataclass
class ConverterSettings:
    """User-configurable defaults applied to new conversion jobs."""

    # Georeferencing defaults (Inversion -> PCSF/PCSM, topo attachment).
    default_epsg: int | None = None
    default_utm_zone: str = ""
    # ModEM adapter defaults.
    station_z_convention: str = "auto"
    air_threshold_ohm_m: float = 1e8
    # EDI <-> EMTF-XML defaults.
    on_loss: str = "warn"
    prefer_spectra: bool = True
    xml_strict: bool = True
    # Output behaviour.
    default_pcsf_format: str = "pcsf"
    log10_view: bool = False
    open_output_folder_after_convert: bool = False
    overwrite_without_asking: bool = False
    # Appearance.
    theme: str = "system"
    # Remembered folders (one per page, keyed by page name).
    last_dirs: dict[str, str] | None = None

    def __post_init__(self) -> None:
        if self.last_dirs is None:
            self.last_dirs = {}


_ORG = "earthai-tech"
_APP = "pycsamt-converter"


def _qsettings():
    from PySide6.QtCore import QSettings  # noqa: PLC0415

    return QSettings(_ORG, _APP)


def load_settings() -> ConverterSettings:
    """Read persisted settings, falling back to defaults for anything unset."""
    qs = _qsettings()
    defaults = ConverterSettings()
    values: dict = {}
    for f in fields(ConverterSettings):
        default = getattr(defaults, f.name)
        raw = qs.value(f.name, default)
        values[f.name] = _coerce(raw, default)
    last_dirs = qs.value("last_dirs", {})
    values["last_dirs"] = dict(last_dirs) if last_dirs else {}
    return ConverterSettings(**values)


def save_settings(settings: ConverterSettings) -> None:
    """Persist *settings* via ``QSettings``."""
    qs = _qsettings()
    for key, value in asdict(settings).items():
        qs.setValue(key, value)
    qs.sync()


def _coerce(raw, default):
    if raw is None:
        return default
    if isinstance(default, bool):
        if isinstance(raw, bool):
            return raw
        return str(raw).strip().lower() in {"1", "true", "yes", "on"}
    if isinstance(default, int) and not isinstance(default, bool):
        try:
            return int(raw)
        except (TypeError, ValueError):
            return default
    if isinstance(default, float):
        try:
            return float(raw)
        except (TypeError, ValueError):
            return default
    return raw
