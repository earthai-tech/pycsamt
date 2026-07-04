.. _emtools_advanced:

pycsamt.emtools.advanced — Novel MT Visualizations
=====================================================

:mod:`pycsamt.emtools.advanced` is, in its own module docstring's
words, a set of "novel visualisations unique to pycsamt v2 — none of
these plots exist in MTPy, MARE2DEM, ModEM or other standard MT
packages." Sixteen functions span single-station tensor diagnostics
(Mohr circles, Argand trajectories, Bode consistency, polar apparent
resistivity, the phase-tensor "period clock"), survey-wide
dimensionality and distortion views (a ternary diagram, a
six-axis distortion radar), full pseudosection summaries (sensitivity
depth, apparent anisotropy, dimensionality depth, rotation invariants,
survey fingerprint, MT composite, SNR), and a strike-stability ribbon
and a geographic coherence network.

Functions
---------

- :func:`~pycsamt.emtools.advanced.plot_impedance_mohr_circles`
- :func:`~pycsamt.emtools.advanced.plot_zt_argand`
- :func:`~pycsamt.emtools.advanced.plot_rho_phase_bode`
- :func:`~pycsamt.emtools.advanced.plot_apparent_resistivity_polar`
- :func:`~pycsamt.emtools.advanced.plot_pt_period_clock`
- :func:`~pycsamt.emtools.advanced.plot_dimensionality_ternary`
- :func:`~pycsamt.emtools.advanced.plot_distortion_radar`
- :func:`~pycsamt.emtools.advanced.plot_sensitivity_depth_section`
- :func:`~pycsamt.emtools.advanced.plot_apparent_anisotropy_section`
- :func:`~pycsamt.emtools.advanced.plot_dimensionality_depth_profile`
- :func:`~pycsamt.emtools.advanced.plot_z_invariants_section`
- :func:`~pycsamt.emtools.advanced.plot_survey_fingerprint`
- :func:`~pycsamt.emtools.advanced.plot_mt_composite_section`
- :func:`~pycsamt.emtools.advanced.plot_snr_section`
- :func:`~pycsamt.emtools.advanced.plot_strike_stability_bands`
- :func:`~pycsamt.emtools.advanced.plot_tf_coherence_network`

Full signatures and parameter documentation live in the
:doc:`API reference <../api/emtools>`.

Worked Example
--------------

Built on **L18PLT** (``data/AMT/WILLY_DATA/``), with **KAP03**
(``data/MT/kap03lmt_edis/``) brought in once to demonstrate a real
bug. Three real bugs turned up while building this example and are
fixed along the way:
:func:`~pycsamt.emtools.advanced.plot_apparent_anisotropy_section`'s
``show_pt_arrows`` option was documented as overlaying phase-tensor
strike arrows but never implemented — the parameter existed and did
nothing;
:func:`~pycsamt.emtools.advanced.plot_z_invariants_section`'s fourth
panel, documented as a distinct anisotropy proxy, computed the exact
same formula as its first panel, an unintentional numerical duplicate;
and :func:`~pycsamt.emtools.advanced.plot_tf_coherence_network` crashed
on any survey whose EDI files carry no per-station coordinates in the
standard header field — the identical class of bug already found and
fixed in :func:`pycsamt.emtools.tensor.plot_phase_tensor_map` (see
:doc:`tensor`), here surfacing in a second, independent function.

The example moves from single-station tensor diagnostics — Mohr
circles, an Argand-space trajectory, a Bode rho/phase consistency
check, a polar apparent-resistivity diagram, and the phase-tensor
period clock — through survey-wide dimensionality and distortion views
(a ternary diagram whose 98.3% mean 3-D membership is the *fourth*
independent confirmation across this gallery of the same strong 3-D
distortion already seen via skew in :doc:`qc`, static shift in
:doc:`ss`, and the phase-tensor dimensionality grid in :doc:`tensor`;
and a six-proxy distortion radar), a Bostick sensitivity-depth section,
the apparent-anisotropy section (exposing the first bug above), a
dimensionality depth profile, the rotation-invariants section (exposing
the second bug above), a compact survey fingerprint, a five-row MT
composite section, an SNR pseudosection, multi-method strike-stability
bands (confirming L18PLT's known lack of a tipper channel), and
finishes with the geographic coherence network (exposing the third bug
above).

.. include:: auto_examples/plot_advanced.rst
