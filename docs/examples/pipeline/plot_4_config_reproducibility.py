r"""
Config-driven pipelines and reproducibility
===========================================

The strongest reproducibility guarantee is a pipeline you can serialise to a
file and rebuild byte-for-byte later. :class:`~pycsamt.pipeline.Pipeline`
round-trips to YAML/JSON/Python, and every ``run`` also drops a
``pipeline.yaml`` beside its outputs — so any result can be traced back to
the exact recipe that made it.
"""

# %%
# A pipeline serialises to YAML
# -----------------------------
# :meth:`Pipeline.to_yaml_string` writes the whole definition — step order,
# codes, and every parameter — as human-readable YAML. This is the artefact
# you commit to version control.

from pycsamt.pipeline import Pipeline, Step

pipe = Pipeline(
    [
        ("select_band", Step("FREQ001", band_hz=(0.01, 10_000))),
        ("notch",       Step("NR001", mains_hz=50)),
        ("static_shift", Step("SS001")),
        ("qc_snap",     Step("QC001")),
    ],
    name="documented_workflow",
)
yaml_text = pipe.to_yaml_string()
print(yaml_text)

# %%
# ... and rebuilds exactly
# ------------------------
# :meth:`Pipeline.from_yaml` reconstructs the pipeline from that file. The
# reloaded object has the identical step sequence and parameters — the basis
# for "run the same processing on next year's survey".

from _pipe_data import scratch_dir

cfg_path = scratch_dir() / "workflow.yaml"
pipe.to_yaml(cfg_path)
reloaded = Pipeline.from_yaml(cfg_path)

original = [(lbl, s.spec.code) for lbl, s in pipe]
restored = [(lbl, s.spec.code) for lbl, s in reloaded]
print("round-trip identical:", original == restored)
for (lbl, code) in restored:
    print(f"  {lbl:<14} {code}")

# %%
# Config files can also compose a preset
# --------------------------------------
# A YAML config may name a ``preset`` and then *append* extra steps, so you
# can start from a curated workflow and add project-specific tweaks — all
# captured in one file. The config format is::
#
#     name: my_survey
#     preset: basic_qc        # seed steps from a preset (optional)
#     output_dir: results/
#     steps:                  # appended after the preset's steps
#       - static_shift: {code: SS001}
#       - qc_snap:      {code: QC001}
#
# :meth:`Pipeline.from_yaml`, :meth:`~pycsamt.pipeline.Pipeline.from_json`,
# and :meth:`~pycsamt.pipeline.Pipeline.from_py` all read this shape.

# %%
# The workflow, drawn from the config
# -----------------------------------
# A tiny diagram of the reloaded pipeline makes the recipe legible at a
# glance — the sequence of steps the config encodes.

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

fig, ax = plt.subplots(figsize=(10, 2.4), constrained_layout=True)
ax.set_axis_off(); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
n = len(restored)
xs = [(i + 0.5) / n for i in range(n)]
for i, (lbl, code) in enumerate(restored):
    ax.add_patch(FancyBboxPatch((xs[i] - 0.4 / n, 0.35), 0.8 / n, 0.34,
                 boxstyle="round,pad=0.01,rounding_size=0.02",
                 fc="#eef3fb", ec="#3e65b0", lw=1.4))
    ax.text(xs[i], 0.57, code, ha="center", va="center", fontsize=9, weight="bold")
    ax.text(xs[i], 0.45, lbl, ha="center", va="center", fontsize=7.5, color="#374151")
    if i < n - 1:
        ax.add_patch(FancyArrowPatch((xs[i] + 0.42 / n, 0.52),
                     (xs[i + 1] - 0.42 / n, 0.52), arrowstyle="-|>",
                     mutation_scale=13, color="#64748b", lw=1.3))
ax.set_title(f"pipeline '{reloaded.name}' — rebuilt from workflow.yaml", fontsize=10)

# %%
# **Takeaway.** Because the pipeline *is* the config, a result and its recipe
# never drift apart: commit the ``pipeline.yaml``, and the run is
# reproducible anywhere. This is what makes preset-plus-tweaks workflows
# auditable across a project and over time.
