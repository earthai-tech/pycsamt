r"""
Stratagem instrument presets
============================

The Zonge Stratagem workflow has its own end-to-end presets in
:mod:`pycsamt.pipeline.stratagem` — full pipelines that go from raw
instrument files (with survey coordinates and reprojection) all the way to
processed, renamed EDIs. This example browses those presets; running one
needs raw Stratagem data plus a coordinate file, so here we focus on the
recipes themselves.
"""

# %%
# The Stratagem catalogue
# -----------------------
# :func:`~pycsamt.pipeline.stratagem_preset_catalogue` prints the available
# instrument workflows and what each does.

from pycsamt.pipeline import (
    list_stratagem_presets, stratagem_preset_catalogue,
)

print(stratagem_preset_catalogue())

# %%
# Presets as objects
# ------------------
# :func:`~pycsamt.pipeline.list_stratagem_presets` returns
# :class:`~pycsamt.pipeline.StratagemPreset` objects. Each bundles
# survey-wide defaults (projection, coordinate handling) with an ordered
# step list — a complete acquisition-to-EDI recipe.

presets = list_stratagem_presets()
print(f"{len(presets)} Stratagem presets:\n")
for p in presets:
    print(f"  {p.name:<18} {len(p.steps)} steps — {p.description.splitlines()[0]}")

# %%
# Inside a preset
# ---------------
# The ``publication_ready`` preset is the fullest — its steps and
# survey defaults show exactly what a finished Stratagem run applies.

pub = next(p for p in presets if "pub" in p.name.lower())
print(f"preset {pub.name!r}: {pub.description}\n")
print("survey defaults:")
for k, v in pub.survey_defaults.items():
    print(f"  {k:<16} {v}")
print("\nsteps:")
for name, params in pub.steps:
    print(f"  {name:<20} {params if params else ''}")

# %%
# How you would run it
# --------------------
# With raw data in hand it is a single call —
# :func:`~pycsamt.pipeline.run_stratagem_preset` — pointing at the raw
# directory and coordinate file::
#
#     from pycsamt.pipeline import run_stratagem_preset
#     result = run_stratagem_preset(
#         "publication_ready",
#         raw_dir="field/stratagem_raw/",
#         outdir="processed/",
#         epsg=32649, utm_zone="49N",
#     )
#
# For step-by-step control, build a
# :class:`~pycsamt.pipeline.StratagemPipeline` directly and call ``run``.

# %%
# Preset lengths at a glance
# --------------------------

import matplotlib.pyplot as plt

names = [p.name for p in presets]
sizes = [len(p.steps) for p in presets]
order = sorted(range(len(presets)), key=lambda i: sizes[i])
fig, ax = plt.subplots(figsize=(7.5, 3.4), constrained_layout=True)
ax.barh([names[i] for i in order], [sizes[i] for i in order], color="#f15a29")
ax.set_xlabel("number of steps")
ax.set_title("Stratagem presets, by workflow length")
ax.margins(x=0.1)

# %%
# **Takeaway.** The Stratagem presets package the whole instrument workflow —
# import, reprojection, processing, export — as reproducible recipes, the
# same philosophy as the generic :doc:`presets <plot_3_presets>` applied to a
# specific acquisition system.
