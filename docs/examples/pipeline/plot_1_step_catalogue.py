r"""
The processing-step catalogue
=============================

A pipeline is built from **steps**, each identified by a short code (like
``NR001`` for the power-line notch). Before assembling a workflow it helps
to see what is available: :mod:`pycsamt.pipeline` ships a registry of 47
steps across eight categories, and this example browses it.
"""

# %%
# Categories and step counts
# --------------------------
# :func:`~pycsamt.pipeline.categories` lists the processing families;
# :func:`~pycsamt.pipeline.list_steps` returns their
# :class:`~pycsamt.pipeline.StepSpec` entries. Together they answer "what
# can this pipeline do?".

from pycsamt.pipeline import categories, list_steps, step_codes, lookup_step

cats = categories()
print(f"{len(step_codes())} steps across {len(cats)} categories\n")
counts = {}
for cat in cats:
    specs = list_steps(cat)
    counts[cat] = len(specs)
    print(f"  {cat:<15} {len(specs):>2} steps")

# %%
# Browsing one category
# ---------------------
# Passing a category name lists just its steps. Here are the noise-removal
# steps — the largest family — with their code, name, and one-line label.

print("noise_removal steps:\n")
for spec in list_steps("noise_removal"):
    print(f"  [{spec.code}]  {spec.name:<24}  {spec.label}")

# %%
# Inspecting one step
# -------------------
# :func:`~pycsamt.pipeline.lookup_step` returns the full spec for a code,
# including its **default parameters** — the values you would override when
# adding it to a pipeline (e.g. ``Step("NR001", mains_hz=60)`` for 60 Hz
# mains).

spec = lookup_step("NR001")
print(f"code      : {spec.code}")
print(f"name      : {spec.name}")
print(f"label     : {spec.label}")
print(f"category  : {spec.category}")
print(f"defaults  : {spec.defaults}")
print(f"returns_sites: {spec.returns_sites}")

# %%
# A catalogue at a glance
# -----------------------
# One small chart makes the shape of the toolbox obvious — where the
# processing effort is concentrated (noise removal and frequency handling
# dominate).

import matplotlib.pyplot as plt

order = sorted(counts, key=counts.get)
fig, ax = plt.subplots(figsize=(8, 4.2), constrained_layout=True)
ax.barh(order, [counts[c] for c in order], color="#3e65b0")
for i, c in enumerate(order):
    ax.text(counts[c] + 0.1, i, str(counts[c]), va="center", fontsize=9)
ax.set_xlabel("number of steps")
ax.set_title(f"pycsamt.pipeline step registry — {len(step_codes())} steps")
ax.margins(x=0.08)

# %%
# **Next.** With the catalogue in hand, :doc:`plot_2_build_and_run` assembles
# a few of these steps into a pipeline and runs it on a real survey line.
