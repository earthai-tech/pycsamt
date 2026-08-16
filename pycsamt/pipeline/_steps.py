"""Step and StepResult — user-facing wrappers around StepSpec entries.

A :class:`Step` binds a registry entry (identified by code or name) with a set
of user-supplied parameter overrides.  It is the atomic unit the user places
into a :class:`~pycsamt.pipeline.Pipeline`.

A :class:`StepResult` is the immutable record produced after a step executes.
The pipeline collects one per step and bundles them into a
:class:`~pycsamt.pipeline.PipelineResult`.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from ._registry import StepSpec, lookup_step

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _to_figure(obj: Any) -> Any:
    """Return a :class:`matplotlib.figure.Figure` from *obj*.

    emtools QC functions return either a ``Figure`` or an ``Axes``.  This
    normalises both so the pipeline always receives a ``Figure`` it can pass
    to ``savefig`` and ``plt.close``.
    """
    try:
        import matplotlib.figure as _mf

        if isinstance(obj, _mf.Figure):
            return obj
        get_fig = getattr(obj, "get_figure", None)
        if callable(get_fig):
            f = get_fig()
            if isinstance(f, _mf.Figure):
                return f
    except Exception:
        pass
    return None


# ---------------------------------------------------------------------------
# Step
# ---------------------------------------------------------------------------


class Step:
    """A configured pipeline step.

    Parameters
    ----------
    code_or_name:
        Registry code (e.g. ``"NR001"``) **or** snake-case name
        (e.g. ``"notch_powerline"``).  Both forms are accepted.
    **params:
        Keyword arguments forwarded to the underlying transform function.
        Values here are merged on top of the registry defaults.

    Examples
    --------
    >>> from pycsamt.emtools.pipe import Step
    >>> Step("NR001", mains_hz=50)
    Step [NR001] Power-line Harmonic Notch  (mains_hz=50, n_harm=30, tol_hz=0.08)
    >>> Step("notch_powerline")
    Step [NR001] Power-line Harmonic Notch  (mains_hz=50, n_harm=30, tol_hz=0.08)
    """

    def __init__(self, code_or_name: str, **params: Any) -> None:
        self.spec: StepSpec = lookup_step(code_or_name)
        self.params: dict = {**self.spec.defaults, **params}

    # ------------------------------------------------------------------
    # Core execution
    # ------------------------------------------------------------------

    def transform(self, sites: Any) -> Any:
        """Apply the step's transform function to *sites*.

        If the step has ``returns_sites=False`` (diagnostic-only), the input
        *sites* is returned unchanged.
        """
        if not self.spec.returns_sites:
            # Diagnostic step — run the function for side-effects but
            # pass sites through unchanged.
            try:
                self.spec.get_fn()(sites, **self.params)
            except Exception as exc:
                warnings.warn(
                    f"Diagnostic step {self.spec.code!r} raised "
                    f"{type(exc).__name__}: {exc}",
                    stacklevel=2,
                )
            return sites
        return self.spec.get_fn()(sites, **self.params)

    def generate_qc_plots(self, sites: Any) -> list:
        """Call QC plot functions on *sites* and return ``(name, Figure)`` pairs.

        Return values are normalised to :class:`matplotlib.figure.Figure` so
        the pipeline can call ``savefig`` regardless of whether the underlying
        emtools function returned a ``Figure`` or an ``Axes`` object.
        Failures are silently skipped so that a broken QC plot never kills
        a successful processing run.
        """
        figs = []
        for fn_name, fn in self.spec.get_qc_fns():
            try:
                result = fn(sites)
                if result is not None:
                    fig = _to_figure(result)
                    if fig is not None:
                        figs.append((fn_name, fig))
            except Exception:
                pass
        return figs

    # ------------------------------------------------------------------
    # Display
    # ------------------------------------------------------------------

    def __repr__(self) -> str:
        params_str = ", ".join(f"{k}={v!r}" for k, v in self.params.items())
        code = self.spec.code
        label = self.spec.label
        if params_str:
            return f"Step [{code}] {label}  ({params_str})"
        return f"Step [{code}] {label}"

    def to_dict(self) -> dict:
        """Serialise step to a plain dict (used by to_yaml / to_json)."""
        return {
            "code": self.spec.code,
            "name": self.spec.name,
            "params": dict(self.params),
        }


# ---------------------------------------------------------------------------
# StepResult
# ---------------------------------------------------------------------------


@dataclass
class StepResult:
    """Immutable record produced after one pipeline step runs.

    Attributes
    ----------
    step_idx:
        1-based position in the pipeline.
    step_name:
        User-supplied label for this step in the pipeline.
    step_code:
        Registry code (e.g. ``"NR001"``).
    step_label:
        Human-readable label from the registry.
    params:
        Parameters that were passed to the transform function.
    elapsed_sec:
        Wall-clock time in seconds for this step (transform + QC plots).
    plots:
        Paths of figures saved to disk during this step.
    n_sites_in:
        Number of sites fed into this step.
    n_sites_out:
        Number of sites coming out of this step.
    error:
        If the step raised an exception and execution continued
        (``on_step_error != "raise"``), the exception is stored here.
        ``None`` means the step completed without error.
    cached:
        ``True`` when this step's result was replayed from the pipeline's
        step cache (``Pipeline.run(..., cache=...)``) instead of actually
        being recomputed.  Only ever ``True`` when caching was enabled for
        this run.
    """

    step_idx: int
    step_name: str
    step_code: str
    step_label: str
    params: dict
    elapsed_sec: float
    plots: list[Path] = field(default_factory=list)
    n_sites_in: int = 0
    n_sites_out: int = 0
    error: Exception | None = None
    cached: bool = False

    @property
    def ok(self) -> bool:
        """``True`` when the step completed without error."""
        return self.error is None

    def summary_line(self) -> str:
        """Return a one-line text summary of this step result."""
        status = "OK" if self.ok else "ERR"
        cached_str = "  (cached)" if self.cached else ""
        return (
            f"  [{self.step_idx:2d}] {self.step_name:<20s} "
            f"[{self.step_code}]  {status}  "
            f"{self.elapsed_sec:.2f}s  "
            f"sites {self.n_sites_in}→{self.n_sites_out}  "
            f"plots={len(self.plots)}{cached_str}"
        )

    def __repr__(self) -> str:
        return (
            f"StepResult(idx={self.step_idx}, code={self.step_code!r}, "
            f"ok={self.ok}, elapsed={self.elapsed_sec:.2f}s)"
        )


__all__ = ["Step", "StepResult"]
