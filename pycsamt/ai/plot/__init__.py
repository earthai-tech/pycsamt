# -*- coding: utf-8 -*-
"""
pycsamt.ai.plot
===============

Unified, publication-quality visualisation for AI-based EM results.

All functions use the :class:`EMStyle` context manager to enforce a
consistent look across all plots:

* **Resistivity sections** — ``RdYlBu_r`` diverging colormap with
  log-normalisation; colorbar on the right labelled in Ω·m or
  log₁₀(Ω·m).
* **Uncertainty** — ``YlOrRd`` sequential colormap for ensemble spread.
* **Typography** — axes 11 pt, title 13 pt, ticks 9 pt.
* **Figure sizes** — single-column 3.5 in, double-column 7.0 in.
* **Background** — white with light-grey grid (``0.92``).

Planned modules (Phase 2 / Phase 4)
-------------------------------------
_style.py  (Phase 2)
    :class:`EMStyle` — ``matplotlib.rc_context``-based context manager
    that applies the shared style sheet.  All other plot modules use it.

compare.py  (Phase 2)
    :func:`plot_compare` — side-by-side true vs predicted 1-D
    resistivity profiles for N sites.  Each panel shows the true model
    (dashed) and the predicted model (solid) with an inset showing the
    data misfit.

convergence.py  (Phase 2)
    :func:`plot_convergence` — training / validation loss curves with
    shaded 1-σ band.  Marks the early-stopping epoch and the minimum
    validation loss.

section.py  (Phase 4)
    :func:`plot_section` — 2-D resistivity section (depth vs distance)
    with optional ensemble spread rendered as an alpha channel.

diagnostics.py  (Phase 4)
    :func:`plot_diagnostics` — calibration curve, residual distribution,
    per-frequency RMSE heat-map, and IoU matrix for binary
    anomaly/salt delineation tasks.
"""
__all__: list = []
