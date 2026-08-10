"""Package-wide contour-line configuration for pyCSAMT plots."""

from __future__ import annotations

import copy
from collections.abc import Generator, Mapping
from contextlib import contextmanager
from dataclasses import dataclass
from typing import Any

_PRESETS = {"default", "subtle", "review", "publication", "off"}


@dataclass
class ContourStyle:
    """Rendering and optional labeling settings for contour overlays."""

    enabled: bool = True
    levels: int | tuple[float, ...] = 7
    colors: Any = "#202020"
    linewidths: float = 0.8
    linestyles: str = "solid"
    alpha: float = 0.8
    labels: bool = False
    label_fmt: str = "%g"
    label_fontsize: float = 7.0
    label_inline: bool = True

    def contour_kwargs(self) -> dict[str, Any]:
        """Return keyword arguments accepted by Matplotlib ``contour``."""
        return {
            "levels": self.levels,
            "colors": self.colors,
            "linewidths": self.linewidths,
            "linestyles": self.linestyles,
            "alpha": self.alpha,
        }

    def label_kwargs(self) -> dict[str, Any]:
        """Return keyword arguments accepted by Matplotlib ``clabel``."""
        return {
            "fmt": self.label_fmt,
            "fontsize": self.label_fontsize,
            "inline": self.label_inline,
        }


class PyCSAMTContour:
    """Package-wide contour presets and live default style."""

    def __init__(self) -> None:
        self._init_defaults()

    def _init_defaults(self) -> None:
        self.subtle = ContourStyle(linewidths=0.35, alpha=0.45)
        self.review = ContourStyle()
        self.publication = ContourStyle(
            colors="#111111",
            linewidths=0.65,
            alpha=0.95,
            labels=True,
            label_fontsize=6.5,
        )
        self.off = ContourStyle(enabled=False)
        self.default = copy.deepcopy(self.review)

    def style_for(self, preset: str = "default") -> ContourStyle:
        """Return the contour style associated with *preset*."""
        name = str(preset).lower()
        if name not in _PRESETS:
            raise ValueError(
                f"contour preset must be one of {sorted(_PRESETS)}."
            )
        return getattr(self, name)

    def use_preset(self, preset: str) -> None:
        """Copy a named preset into the live default style."""
        if str(preset).lower() == "default":
            return
        self.default = copy.deepcopy(self.style_for(preset))

    def configure(self, **kw: Any) -> None:
        """Configure styles with ``preset__attribute`` paths.

        A key without a preset prefix configures the live ``default`` style,
        so ``linewidths=1.0`` is equivalent to
        ``default__linewidths=1.0``.
        """
        for path, value in kw.items():
            parts = path.split("__")
            if len(parts) == 1:
                obj = self.default
                attr = parts[0]
            else:
                obj = self
                for part in parts[:-1]:
                    obj = getattr(obj, part)
                attr = parts[-1]
            if not hasattr(obj, attr):
                raise AttributeError(f"unknown contour setting: {path!r}")
            setattr(obj, attr, value)

    def resolve(
        self,
        enabled: bool | None = None,
        overrides: Mapping[str, Any] | None = None,
    ) -> tuple[bool, dict[str, Any], dict[str, Any]]:
        """Resolve enabled state, contour kwargs, and label kwargs."""
        style = self.default
        is_enabled = style.enabled if enabled is None else bool(enabled)
        contour_kwargs = style.contour_kwargs()
        contour_kwargs.update(dict(overrides or {}))
        label_kwargs = style.label_kwargs() if style.labels else {}
        return is_enabled, contour_kwargs, label_kwargs

    @contextmanager
    def context(
        self,
        preset: str | None = None,
        **kw: Any,
    ) -> Generator[PyCSAMTContour, None, None]:
        """Temporarily apply a preset or overrides, then restore state."""
        snapshot = self._snapshot()
        try:
            if preset:
                self.use_preset(preset)
            if kw:
                self.configure(**kw)
            yield self
        finally:
            self._restore(snapshot)

    def reset(self) -> None:
        """Restore every contour preset and the live default."""
        self._init_defaults()

    def summary(self) -> str:
        """Return a compact description of the live contour style."""
        style = self.default
        return "\n".join(
            [
                "PyCSAMTContour",
                f"  enabled     = {style.enabled}",
                f"  levels      = {style.levels!r}",
                f"  colors      = {style.colors!r}",
                f"  linewidths  = {style.linewidths}",
                f"  linestyles  = {style.linestyles!r}",
                f"  alpha       = {style.alpha}",
                f"  labels      = {style.labels}",
            ]
        )

    def __repr__(self) -> str:  # noqa: D105
        return self.summary()

    def _snapshot(self) -> dict[str, ContourStyle]:
        return {
            name: copy.deepcopy(getattr(self, name))
            for name in _PRESETS
        }

    def _restore(self, snapshot: dict[str, ContourStyle]) -> None:
        for name, style in snapshot.items():
            setattr(self, name, style)


PYCSAMT_CONTOUR = PyCSAMTContour()


def configure_contour(**kw: Any) -> None:
    """Configure the live :data:`PYCSAMT_CONTOUR` style."""
    PYCSAMT_CONTOUR.configure(**kw)


def use_contour(preset: str) -> None:
    """Copy a named contour preset into the live default style."""
    PYCSAMT_CONTOUR.use_preset(preset)


def reset_contour() -> None:
    """Restore package contour defaults."""
    PYCSAMT_CONTOUR.reset()


__all__ = [
    "ContourStyle",
    "PYCSAMT_CONTOUR",
    "PyCSAMTContour",
    "configure_contour",
    "reset_contour",
    "use_contour",
]
