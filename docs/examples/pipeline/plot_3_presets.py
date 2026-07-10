r"""
Reproducible presets
=====================

Typing out the same step chain every time is error-prone, so
:mod:`pycsamt.pipeline` ships **presets** — named, curated workflows you can
run in one line. This example browses the catalogue and runs one end to end.
"""

# %%
# The preset catalogue
# --------------------
# :func:`~pycsamt.pipeline.preset_catalogue` prints every preset with its
# purpose and step list — the menu of ready-made workflows.

from pycsamt.pipeline import (
    Pipeline,
    configure_pipe,
    get_preset,
    list_presets,
    preset_catalogue,
)

print(preset_catalogue())

# %%
# Presets as objects
# ------------------
# :func:`~pycsamt.pipeline.list_presets` returns
# :class:`~pycsamt.pipeline.Preset` objects; each carries a ``name``,
# ``description`` and its ordered ``steps``.

presets = list_presets()
print(f"{len(presets)} presets:\n")
for p in presets:
    print(f"  {p.name:<18} {len(p.steps)} steps  — {p.description.splitlines()[0]}")

# %%
# Inspect one preset
# ------------------
# :func:`~pycsamt.pipeline.get_preset` fetches one by name. Its ``steps`` are
# exactly the ``(label, Step)`` pairs you would otherwise write by hand.

preset = get_preset("basic_qc")
print(f"preset {preset.name!r}: {preset.description}\n")
for label, step in preset.steps:
    print(f"  {label:<14} {step}")

# %%
# Run a preset
# ------------
# A preset drops straight into :class:`~pycsamt.pipeline.Pipeline` — its
# steps *are* a pipeline definition — so running it is one call on the data.

from _pipe_data import demo_sites, quiet_logs, scratch_dir

sites = demo_sites(n=8)
configure_pipe(show_progress=False, plot_dpi=72)
pipe = Pipeline(preset.steps, name=preset.name)

with quiet_logs():
    result = pipe.run(sites, outdir=scratch_dir(),
                      save_plots=False, save_edis=True, save_report=True)
print(result.summary())

# %%
# Comparing the presets
# ---------------------
# A quick chart of how much processing each preset applies — from the light
# ``basic_qc`` to the full end-to-end workflows.

import matplotlib.pyplot as plt

names = [p.name for p in presets]
sizes = [len(p.steps) for p in presets]
order = sorted(range(len(presets)), key=lambda i: sizes[i])
fig, ax = plt.subplots(figsize=(8, 4.2), constrained_layout=True)
ax.barh([names[i] for i in order], [sizes[i] for i in order], color="#fbb040")
ax.set_xlabel("number of steps")
ax.set_title("pycsamt.pipeline presets, by workflow length")
ax.margins(x=0.08)

# %%
# **Takeaway.** Presets make a standard workflow a one-liner and a shared
# vocabulary across a project. To pin an *exact* run — preset plus your own
# tweaks — serialise it to a config file: see
# :doc:`plot_4_config_reproducibility`.
