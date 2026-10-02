"""Before/after QC plots in the pipeline (``call_qc_fn`` / ``generate_qc_plots``).

Regression: every comparison QC plot registered for the static-shift
(SS001-SS004) and frequency-edit (FREQ003/008/009) steps needs the step's
input *and* output, but was called as ``fn(sites)``.  The ``TypeError`` was
swallowed by ``generate_qc_plots``, so 13 of the 81 registered QC figures
were silently never produced — in Pipeline.run, the CLI and the desktop.
"""

from __future__ import annotations

import importlib
import inspect
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest

from pycsamt.pipeline import Pipeline, Step, list_steps
from pycsamt.pipeline._steps import _required_positional, call_qc_fn

_WILLY = (Path(__file__).parents[3] / "data" / "AMT" / "WILLY_DATA"
          / "L18PLT")


@pytest.fixture(scope="module")
def sites():
    if not (_WILLY.exists() and any(_WILLY.glob("*.edi"))):
        pytest.skip("WILLY L18PLT data not available")
    from pycsamt.emtools import ensure_sites

    return ensure_sites(str(_WILLY))


@pytest.fixture(autouse=True)
def _close():
    yield
    plt.close("all")


def test_every_registered_qc_function_has_a_supported_shape():
    """A new QC plot with an unsupported signature must fail loudly here,
    not silently at run time."""
    for spec in list_steps():
        for mod, name in spec.qc_defs or []:
            fn = getattr(importlib.import_module(mod), name)
            req = _required_positional(fn)
            assert len(req) <= 2, (spec.code, name, req)
            if len(req) == 2:
                names = [r.lower() for r in req]
                ok = (names[0].startswith("logrho")
                      or "before" in names[0] and "after" in names[1])
                assert ok, (spec.code, name, req)


@pytest.mark.parametrize(
    "code,expected",
    [("SS001", 2), ("SS002", 2), ("SS004", 3), ("FREQ003", 2)],
)
def test_comparison_plots_render_with_before(sites, code, expected):
    step = Step(code)
    out = step.transform(sites)
    assert step.generate_qc_plots(out) == [] or code.startswith("FREQ00")
    figs = step.generate_qc_plots(out, before=sites)
    assert len(figs) == expected


def test_call_qc_fn_skips_comparison_without_before(sites):
    from pycsamt.emtools.ss import plot_ss_delta_psection

    assert call_qc_fn(plot_ss_delta_psection, sites) is None


def test_ss_logrho_arrays_shapes(sites):
    from pycsamt.emtools import correct_ss_ama, ss_logrho_arrays

    after = correct_ss_ama(sites, inplace=False, verbose=0)
    b, a, freqs, labels = ss_logrho_arrays(sites, after)
    assert b.shape == a.shape == (len(labels), freqs.size)
    assert np.isfinite(b).any() and len(labels) == len(sites)


def test_pipeline_run_saves_static_shift_qc_figures(sites, tmp_path):
    pipe = Pipeline([("ss", Step("SS001"))], name="qc_check")
    result = pipe.run(sites, outdir=tmp_path, save_edis=False,
                      save_report=False)
    names = {p.stem for p in result.plots}
    assert any("plot_ss_delta_psection" in n for n in names)
    assert any("plot_ss_summary" in n for n in names)


def test_signature_introspection_ignores_keyword_defaults():
    def f(a, b=1, *, c=2):
        return a

    assert _required_positional(f) == ["a"]
    assert inspect.isfunction(call_qc_fn)


def test_run_still_accepts_old_style_qc_steps():
    """Custom/plugin steps with the pre-``before`` signature keep working."""
    from pycsamt.pipeline._pipeline import _qc_plots

    class OldStyle:
        def generate_qc_plots(self, sites):
            return [("old", sites)]

    class NewStyle:
        def generate_qc_plots(self, sites, before=None):
            return [("new", (before, sites))]

    assert _qc_plots(OldStyle(), "after", "before") == [("old", "after")]
    assert _qc_plots(NewStyle(), "after", "before") == [
        ("new", ("before", "after"))]
