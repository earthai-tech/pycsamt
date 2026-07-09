r"""
Build and run a pipeline
========================

This is the core workflow: assemble an ordered list of
:class:`~pycsamt.pipeline.Step` s, run them on a survey with
:meth:`Pipeline.run`, and read the :class:`~pycsamt.pipeline.PipelineResult`
— per-step timings, the run summary, and the on-disk output package (cleaned
EDIs, QC plots, and a reproducible config).
"""

# %%
# Assemble the pipeline
# ---------------------
# Each entry is a ``(label, Step(code, **params))`` pair; the steps run in
# order. This chain de-duplicates and band-limits frequencies, aligns the
# frequency axis, notches the 50 Hz power line, then snaps a QC summary.

from _pipe_data import demo_sites, scratch_dir, quiet_logs
from pycsamt.pipeline import Pipeline, Step, configure_pipe, reset_pipe

sites = demo_sites(n=8)
print(f"loaded {len(sites)} stations from WILLY_DATA L22")

pipe = Pipeline(
    [
        ("drop_dup",    Step("FREQ002")),
        ("select_band", Step("FREQ001", band_hz=(0.01, 10_000))),
        ("align",       Step("FREQ004")),
        ("notch",       Step("NR001", mains_hz=50)),
        ("qc_snap",     Step("QC001")),
    ],
    name="willy_l22_demo",
)
print(pipe)

# %%
# Run it
# ------
# :meth:`Pipeline.run` executes every step, writes outputs to ``outdir``,
# and returns the result. ``configure_pipe`` just turns off the progress bar
# and sets a light plot DPI for the docs build.

configure_pipe(show_progress=False, plot_dpi=72, plot_fmt="png")
outdir = scratch_dir()

with quiet_logs():
    result = pipe.run(sites, outdir=outdir,
                      save_plots=True, save_edis=True, save_report=True)

print(result.summary())

# %%
# Read the result
# ---------------
# :class:`~pycsamt.pipeline.PipelineResult` exposes the run at every level:
# overall status, and a :class:`~pycsamt.pipeline.StepResult` per step with
# its timing and in/out station counts — the audit trail of what happened.

print(f"ok = {result.ok}   errors = {result.n_errors}   "
      f"elapsed = {result.elapsed_sec:.2f} s\n")
print(f"{'step':<12}{'code':<9}{'ok':<5}{'seconds':>8}{'sites':>10}")
for sr in result.step_results:
    print(f"{sr.step_name:<12}{sr.step_code:<9}{str(sr.ok):<5}"
          f"{sr.elapsed_sec:>8.3f}{f'{sr.n_sites_in}->{sr.n_sites_out}':>10}")

# %%
# The output package
# ------------------
# The run leaves a self-contained folder: the cleaned EDIs, the QC figures
# each plotting step produced, and ``pipeline.yaml`` — the recipe needed to
# reproduce the run (see :doc:`plot_4_config_reproducibility`).

from pathlib import Path

files = sorted(p.relative_to(outdir).as_posix() for p in outdir.rglob("*") if p.is_file())
print(f"{len(files)} files written to the output directory:")
for f in files:
    print("  ", f)

# %%
# What a step produced
# --------------------
# Since ``save_plots=True``, the QC step wrote real diagnostic figures. Here
# is one of them, loaded straight from the output package — this is the
# pipeline's own output, not a plot made for this page.

import matplotlib.pyplot as plt
import matplotlib.image as mpimg

pngs = sorted(outdir.rglob("*.png"))
pick = next((p for p in pngs if "snr" in p.name.lower()), pngs[0])
img = mpimg.imread(pick)
fig, ax = plt.subplots(figsize=(9, 5.0), constrained_layout=True)
ax.imshow(img)
ax.set_axis_off()
ax.set_title(f"pipeline output: {pick.name}", fontsize=10)

reset_pipe()

# %%
# **Takeaway.** A pipeline turns a list of step codes into a fully
# documented run — cleaned data, QC figures, timings, and a reproducible
# config — in one call. The next examples package common chains as
# :doc:`presets <plot_3_presets>` and reproduce runs from
# :doc:`config files <plot_4_config_reproducibility>`.
