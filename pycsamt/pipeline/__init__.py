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

    pipe = Pipeline([
        ("notch",      Step("NR001", mains_hz=50)),
        ("band",       Step("FREQ001")),
        ("align",      Step("FREQ004")),
        ("correct_ss", Step("SS001")),
        ("rotate",     Step("TZ001")),
    ])
    print(pipe)
    result = pipe.run(sites, outdir="willy_results/")

    # shortest path — Pipeline and Step are also on pycsamt directly
    from pycsamt import Pipeline, Step

Load from a config file::

    pipe = Pipeline.from_yaml("config/workflow.yaml")
    pipe = Pipeline.from_preset("publication_ready")

Discover available steps and presets::

    from pycsamt.pipeline import list_steps, preset_catalogue
    list_steps()                  # all 33 StepSpec objects
    list_steps("noise_removal")   # by category
    print(preset_catalogue())     # named presets

Configure pipeline output globally::

    from pycsamt.pipeline import configure_pipe
    configure_pipe(plot_dpi=300, plot_fmt="pdf", output_root="results/")

    # or temporarily with a context manager
    from pycsamt.pipeline import PYCSAMT_PIPE
    with PYCSAMT_PIPE.context(show_progress=False):
        result = pipe.run(sites)
"""

from __future__ import annotations

# ── Pipeline engine (relative — lives in this package) ───────────────────────
from ._base import PipelineBase
from ._pipeline import Pipeline, PipelineResult
from ._steps import Step, StepResult
from ._registry import (
    StepSpec,
    STEP_REGISTRY,
    lookup_step,
    list_steps,
    step_codes,
    step_names,
    categories,
)
from ._presets import (
    Preset,
    PRESETS,
    get_preset,
    list_presets,
    preset_catalogue,
)
from ._config import load_yaml, load_json, load_py
from .stratagem import (
    StratagemPreset,
    STRATAGEM_PRESETS,
    StratagemPipeline,
    get_stratagem_preset,
    list_stratagem_presets,
    run_stratagem_preset,
    stratagem_preset_catalogue,
)

# ── Pipeline API configuration (relative — sibling package api/pipe) ─────────
from ..api.pipe import (
    PYCSAMT_PIPE,
    PipelineAPIConfig,
    configure_pipe,
    reset_pipe,
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
    # ── presets
    "Preset",
    "PRESETS",
    "get_preset",
    "list_presets",
    "preset_catalogue",
    # ── config loaders
    "load_yaml",
    "load_json",
    "load_py",
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
