"""pyCSAMT processing pipeline.

A declarative, automated MT processing engine that chains emtools operations
from raw Sites to publication-ready output.

This package is the **single entry point** for all pipeline-related symbols.
The processing engine lives in the private sub-modules of this package
(``_registry``, ``_steps``, ``_pipeline``, …).  Runtime configuration is
provided by :mod:`pycsamt.api.pipe` and re-exported here so users only need
one import location.

Typical usage
-------------
Build and run from code::

    from pycsamt.pipeline import Pipeline, Step

    pipe = Pipeline(
        [
            ("notch", Step("NR001", mains_hz=50)),
            ("band", Step("FREQ001")),
            ("align", Step("FREQ004")),
            ("correct_ss", Step("SS001")),
            ("rotate", Step("TZ001")),
        ]
    )
    print(pipe)
    result = pipe.run(sites, outdir="willy_results/")

    # shortest path — Pipeline and Step are also on pycsamt directly
    from pycsamt import Pipeline, Step

Load from a config file::

    pipe = Pipeline.from_yaml("config/workflow.yaml")
    pipe = Pipeline.from_preset("publication_ready")

Discover available steps and presets::

    from pycsamt.pipeline import list_steps, preset_catalogue

    list_steps()  # all built-in StepSpec objects
    list_steps("noise_removal")  # by category
    print(preset_catalogue())  # named presets

Configure pipeline output globally::

    from pycsamt.pipeline import configure_pipe

    configure_pipe(plot_dpi=300, plot_fmt="pdf", output_root="results/")

    # or temporarily with a context manager
    from pycsamt.pipeline import PYCSAMT_PIPE

    with PYCSAMT_PIPE.context(show_progress=False):
        result = pipe.run(sites)

Cache step outputs so an interrupted run resumes instead of recomputing::

    result = pipe.run(sites, outdir="results/", cache=True)

Observe a run live, or log it for later comparison::

    result = pipe.run(sites, on_step=lambda sr: print(sr.summary_line()))
    result = pipe.run(sites, history=True)  # -> load_history() later

A finished ``PipelineResult`` renders inline in Jupyter (``_repr_html_``),
and ``configure_pipe(progress_style="rich")`` (or CLI ``--live``) renders a
live-updating status table in the terminal while a run is in progress.
"""

from __future__ import annotations

# ── Pipeline API configuration (relative — sibling package api/pipe) ─────────
from ..api.pipe import (
    PYCSAMT_PIPE,
    PipelineAPIConfig,
    configure_pipe,
    reset_pipe,
)

# ── Pipeline engine (relative — lives in this package) ───────────────────────
from ._base import PipelineBase
from ._config import load_json, load_py, load_yaml
from ._pipeline import Pipeline, PipelineResult
from ._presets import (
    PRESETS,
    Preset,
    get_preset,
    get_preset_for_method,
    list_presets,
    preset_catalogue,
)
from ._cache import StepCache, fingerprint_sites
from ._history import load_history
from ._plugins import ENTRY_POINT_GROUP, PluginLoadResult, discover_plugins
from .ai_steps import register_ai_steps
from ._registry import (
    STEP_REGISTRY,
    StepSpec,
    categories,
    list_steps,
    lookup_step,
    register_step,
    step_codes,
    step_names,
    unregister_step,
)
from ._steps import Step, StepResult
from .plot import (
    plot_pipeline_dashboard,
    plot_pipeline_status,
    plot_pipeline_timing,
    plot_site_count_flow,
)
from .stratagem import (
    STRATAGEM_PRESETS,
    StratagemPipeline,
    StratagemPreset,
    get_stratagem_preset,
    list_stratagem_presets,
    run_stratagem_preset,
    stratagem_preset_catalogue,
)

__all__ = [
    # ── engine core
    "PipelineBase",
    "Pipeline",
    "PipelineResult",
    "Step",
    "StepResult",
    # ── registry
    "StepSpec",
    "STEP_REGISTRY",
    "lookup_step",
    "list_steps",
    "step_codes",
    "step_names",
    "categories",
    "register_step",
    "unregister_step",
    # ── plugins
    "discover_plugins",
    "PluginLoadResult",
    "ENTRY_POINT_GROUP",
    # ── opt-in AI steps
    "register_ai_steps",
    # ── step cache
    "StepCache",
    "fingerprint_sites",
    # ── run history
    "load_history",
    # ── presets
    "Preset",
    "PRESETS",
    "get_preset",
    "get_preset_for_method",
    "list_presets",
    "preset_catalogue",
    # ── config loaders
    "load_yaml",
    "load_json",
    "load_py",
    # ── plotting
    "plot_pipeline_dashboard",
    "plot_pipeline_status",
    "plot_pipeline_timing",
    "plot_site_count_flow",
    # ── stratagem pipeline
    "StratagemPreset",
    "STRATAGEM_PRESETS",
    "StratagemPipeline",
    "get_stratagem_preset",
    "list_stratagem_presets",
    "run_stratagem_preset",
    "stratagem_preset_catalogue",
    # ── runtime configuration
    "PYCSAMT_PIPE",
    "PipelineAPIConfig",
    "configure_pipe",
    "reset_pipe",
]
