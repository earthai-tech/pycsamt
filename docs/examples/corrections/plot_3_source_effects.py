r"""
Source effects (CSAMT)
======================

A controlled-source survey adds a distortion that natural-source MT never
sees: near the transmitter, at low frequency, the field is **near-field** —
the plane-wave assumption breaks down and apparent resistivity rises with a
tell-tale 45° slope. :mod:`pycsamt.emtools` detects this source overprint and
corrects it. This wave diagnoses the overprint first, then removes it.
"""

# %%
# Detect the overprint
# --------------------
# :func:`~pycsamt.emtools.detect_source_overprint` flags, per station, the
# low-frequency band where the response departs from plane-wave behaviour —
# the diagnostic you check before deciding to correct.

import numpy as np

from _corr_data import demo_line, curves, plot_before_after
from pycsamt.emtools import (
    detect_source_overprint, source_overprint_table, correct_near_field,
)

S = demo_line("L18PLT")
raw = curves(S, "rho")

flags = detect_source_overprint(S, recursive=False)
print("detect_source_overprint ->", type(flags).__name__)
print(flags.head(10).to_string(index=False))

# %%
# A survey-wide overprint table
# -----------------------------
# :func:`~pycsamt.emtools.source_overprint_table` summarises the diagnostic
# across the whole line — how many stations are affected and over what period
# band — so you can judge whether a near-field correction is warranted.

table = source_overprint_table(S, recursive=False)
print(table.head(12).to_string(index=False))

# %%
# Correct the near field
# ----------------------
# :func:`~pycsamt.emtools.correct_near_field` removes the near-field rise
# given the transmitter-receiver offset (``source_offset``, metres). It
# rescales the affected low-frequency band back toward the plane-wave
# response, leaving the far-field high frequencies untouched.

S_nf = correct_near_field(S, source_offset=5000.0, recursive=False)
nf = curves(S_nf, "rho")

stations = list(raw)
pick = [stations[2], stations[len(stations) // 2], stations[-3]]
plot_before_after(raw, nf, pick, quantity="rho",
                  labels=("raw", "near-field corrected"),
                  title="Near-field / source-overprint correction "
                        "(source_offset = 5 km)")

# %%
# Where the correction acted
# --------------------------
# Imaging the change station-by-period confirms the correction is confined to
# the low-frequency (long-period) near-field band and leaves the far field
# alone — the signature of a physically-sensible source correction.

names = [s for s in raw if np.array_equal(raw[s][0], nf[s][0])]
periods = raw[names[0]][0]
delta = np.column_stack([np.log10(nf[s][1]) - np.log10(raw[s][1]) for s in names])
vmax = float(np.nanpercentile(np.abs(delta), 98)) or 0.1
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(10, 4.6), constrained_layout=True)
im = ax.pcolormesh(np.arange(len(names)), np.log10(periods), delta,
                   cmap="PuOr_r", vmin=-vmax, vmax=vmax, shading="auto")
ax.set_xticks(range(len(names)))
ax.set_xticklabels(names, rotation=90, fontsize=6)
ax.set_ylabel(r"$\log_{10}$ period (s)")
ax.set_title(r"$\Delta\log_{10}\rho_a$ from the near-field correction")
fig.colorbar(im, ax=ax, label=r"$\Delta\log_{10}\rho_a$")

# %%
# **Takeaway.** Always *detect* the source overprint before correcting, set
# ``source_offset`` from the real survey geometry, and confirm the change is
# confined to the low-frequency band. With the source distortion removed, the
# tensor is ready to be rotated onto strike —
# :doc:`tensor rotation <plot_4_tensor_rotation>`.

# sphinx_gallery_thumbnail_number = 1
