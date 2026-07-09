r"""
Processing a tipper survey (spectra to induction arrows)
=======================================================

The corrections so far worked on the impedance tensor; this one follows the
**vertical-field tipper** through a complete run. MT data begins life as raw
cross-spectra, from which both the impedance *and* the tipper are derived —
so we open at the spectra level, then process a real tipper-bearing line
(**KAP03**, a SAMTEX long-period profile whose station coordinates live in
the hidden ``REFLAT``/``REFLONG`` fields), watching the induction response
clean up and writing tipper-preserving EDIs at the end.
"""

# %%
# Where tipper comes from: the raw spectra
# ----------------------------------------
# :meth:`Spectra.from_file <pycsamt.seg.spectra.Spectra.from_file>` reads the
# ``>=SPECTRASECT`` cross-power block. :func:`~pycsamt.emtools.plot_psd`
# shows the channel power spectra — the signal levels every transfer function
# is estimated from.

import os
from pathlib import Path

import numpy as np

from pycsamt.seg.spectra import Spectra
from pycsamt.emtools import (
    spectra_summary, plot_psd, plot_tipper_from_spectra,
)

ROOT = Path(os.environ.get("PYCSAMT_DOCS_REPO_ROOT", "."))
sp = Spectra.from_file(str(ROOT / "data" / "MT" / "SPECTRA" / "spectra01.edi"))
print("spectra summary:")
print(spectra_summary(sp))
_ax = plot_psd(sp)

# %%
# The tipper, recovered from the spectra
# --------------------------------------
# :func:`~pycsamt.emtools.plot_tipper_from_spectra` solves the vertical-field
# transfer function ``(Tx, Ty)`` straight from the cross-spectra — the tipper
# at its origin, before it is written into an EDI. This is the quantity the
# rest of the workflow conditions.

_tzx = plot_tipper_from_spectra(sp)

# %%
# Load the tipper line
# --------------------
# KAP03 ships the derived impedance **and** tipper. Its 26 stations span
# 25 s - 17,067 s, and unlike the AMT lines it carries a genuine vertical
# field, so every induction arrow below is real tipper.

from pycsamt.emtools._core import ensure_sites, _iter_items
from pycsamt.emtools import (
    drop_duplicates, select_band, drop_low_confidence_frequencies,
    notch_powerline, smooth_logfreq,
    plot_induction_arrows, plot_induction_section, plot_induction_rose,
)
from pycsamt.site.export import write_sites

K = ensure_sites(str(ROOT / "data" / "MT" / "kap03lmt_edis"))


def n_freq(sites):
    return np.array([len(np.asarray(e.freq)) for e in _iter_items(sites)])


print(f"KAP03: {len(n_freq(K))} stations, {n_freq(K).max()} frequencies "
      f"(period 25-17067 s), real tipper channel")

# %%
# Raw induction arrows
# --------------------
# :func:`~pycsamt.emtools.plot_induction_arrows` draws the real induction
# vectors at several periods. Under the Parkinson convention they point
# *away* from conductors, so a coherent arrow pattern is real structure and a
# noisy one flags tipper that needs cleaning.

ax = plot_induction_arrows(K)

# %%
# Raw tipper section
# ------------------
# :func:`~pycsamt.emtools.plot_induction_section` images tipper magnitude
# across station and period — effectively a *spectral* view of the vertical
# field. The raw section is speckled where the tipper estimate is noisy.

ax = plot_induction_section(K, component="real")

# %%
# Process the survey (tipper carried through)
# -------------------------------------------
# The correction functions accept ``also="both"``, so they clean the
# **tipper alongside the impedance**. We trim the frequency axis, then notch
# and smooth — each step returning a new ``Sites`` whose tipper is processed,
# not dropped.

s1 = drop_duplicates(K, recursive=False)
s2 = select_band(s1, fmin=1e-4, fmax=5e-2, recursive=False)
s3 = drop_low_confidence_frequencies(s2, threshold=0.5, recursive=False)
s4 = notch_powerline(s3, also="both", recursive=False)
final = smooth_logfreq(s4, win=5, also="both", recursive=False)
print(f"frequencies: {n_freq(K).max()} raw -> {n_freq(final).max()} processed "
      f"(band-trimmed and confidence-pruned)")
print("chain: raw -> drop_dup -> select_band -> drop_low_conf -> "
      "notch(both) -> smooth(both)")

# %%
# Processed tipper section
# ------------------------
# The same section on the processed data: band-limited to the trustworthy
# periods and smoothed, so the coherent induction signal stands clear of the
# noise the raw section carried.

ax = plot_induction_section(final, component="real")

# %%
# Processed induction rose
# ------------------------
# :func:`~pycsamt.emtools.plot_induction_rose` folds the arrow azimuths into
# a rose — after cleaning, a tighter preferred direction confirms a coherent
# regional conductor rather than scatter.

ax = plot_induction_rose(final, component="real")

# %%
# Write tipper-preserving EDIs
# ----------------------------
# :func:`~pycsamt.site.export.write_sites` serialises the processed line back
# to EDIs with the cleaned tipper intact — the deliverable for a
# tipper-inclusive inversion.

import tempfile

outdir = Path(tempfile.mkdtemp(prefix="kap03_tipper_"))
paths = write_sites(final, outdir, exist_ok=True)
reloaded = ensure_sites(str(outdir))
print(f"wrote {len(paths)} tipper-preserving EDIs; re-loaded "
      f"{len(n_freq(reloaded))} with {n_freq(reloaded).max()} frequencies")

# %%
# **Takeaway.** Tipper rides through the same ``Sites -> Sites`` pipeline as
# the impedance: process it with ``also="both"``, watch the induction section
# and rose clean up, and write EDIs that keep the vertical field. Combined
# with the impedance waves (static shift, rotation), this completes a
# tipper-inclusive processing run — spectra in, sanitised induction arrows
# out.

# sphinx_gallery_thumbnail_number = 3
