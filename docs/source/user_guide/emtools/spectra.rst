.. _emtools_spectra:

Cross-Spectra Analysis
======================

``pycsamt.emtools.spectra`` works below the processed impedance tensor.
Most ``emtools`` workflows start from an EDI ``Z`` tensor. This module
starts from the raw cross-power spectral matrix stored in a
``pycsamt.seg.spectra.Spectra`` object: one complex Hermitian
``(n_channel, n_channel)`` matrix for each frequency in an EDI
``>=SPECTRASECT`` block.

Use this module when you need to inspect spectra before they become
impedance and tipper estimates:

* channel power spectral density (PSD);
* inter-channel squared coherence;
* coherence-derived SNR;
* frequency-band selection from coherence;
* the full cross-power matrix at one frequency;
* impedance and tipper recovered from spectra;
* PSD and coherence pseudo-sections across multiple spectra stations.

Full function signatures and parameter defaults are maintained in the
:doc:`API reference <../../api/emtools>`. The analysis and plotting
functions are exported through the public ``pycsamt.emtools`` namespace.
The ``Spectra`` container itself lives in ``pycsamt.seg.spectra``.

When To Use This Page
---------------------

Use ``spectra`` tools only when your EDI file contains a
``>=SPECTRASECT`` block. Many EDI files contain only final impedance and
tipper estimates. Those files are valid for most pyCSAMT workflows, but
they do not contain the raw cross-power matrices required here.

The bundled spectra examples are under ``data/MT/SPECTRA/``:

.. code-block:: text

   data/MT/SPECTRA/spectra01.edi
   data/MT/SPECTRA/spectra02.edi

They are de-identified field spectra files. Their identifying metadata
was sanitized, while the numeric spectra blocks were kept for
documentation and examples.

Workflow Map
------------

.. list-table::
   :header-rows: 1
   :widths: 26 36 38

   * - Task
     - Use this
     - Output
   * - Compute coherence cube
     - ``coherence_matrix``
     - Array with shape ``(n_freq, n_chan, n_chan)``.
   * - Build PSD table
     - ``psd_table``
     - Tidy table with station, frequency, channel, and PSD.
   * - Build coherence table
     - ``coherence_table``
     - Tidy table with channel pairs and squared coherence.
   * - Estimate SNR from coherence
     - ``snr_table``
     - Coherence table plus ``snr`` and ``snr_db``.
   * - Restrict frequency band
     - ``band_select``
     - New ``Spectra`` object sliced to ``[f_min, f_max]``.
   * - Create coherence mask
     - ``mask_low_coherence``
     - Boolean mask of frequencies that pass coherence criteria.
   * - Summarize spectra
     - ``spectra_summary``
     - Per-frequency table with PSD and mean coherence.
   * - Plot single-station spectra
     - ``plot_psd``, ``plot_coherence``, ``plot_spectra_matrix``
     - PSD curves, pairwise coherence, and matrix image.
   * - Recover transfer functions
     - ``plot_z_from_spectra``, ``plot_tipper_from_spectra``
     - Apparent resistivity/phase and induction tipper from spectra.
   * - Plot station sections
     - ``plot_psd_section``, ``plot_coherence_section``
     - Station-period pseudo-sections from multiple spectra objects.

Loading Spectra
---------------

Load spectra with ``Spectra.from_file``. This is the one required
three-level import on this page because ``Spectra`` is a segmentation
container, not an ``emtools`` function.

.. code-block:: python
   :linenos:

   from pathlib import Path

   from pycsamt.seg.spectra import Spectra

   spectra_dir = Path("data/MT/SPECTRA")
   sp1 = Spectra.from_file(spectra_dir / "spectra01.edi")
   sp2 = Spectra.from_file(spectra_dir / "spectra02.edi")

   print(sp1.name, sp1.n_freq, sp1.n_chan)
   print(sp1.freq.min(), sp1.freq.max())
   print(sp1.id_to_chtype)

A ``Spectra`` object exposes the core arrays used by this module:

.. code-block:: python
   :linenos:

   print(sp1.freq.shape)       # frequency vector
   print(sp1.S.shape)          # cross-power matrix: (n_freq, n_chan, n_chan)
   print(sp1.bw.shape)         # bandwidth metadata
   print(sp1.avgt.shape)       # averaging-time metadata
   print(sp1.chan_ids)         # channel identifiers

The matrix ``sp1.S[k]`` is the full channel-by-channel cross-power
matrix at one frequency.

Coherence Matrix
----------------

``coherence_matrix`` computes squared coherence for every channel pair:

.. math::

   \gamma^2_{ij}
   =
   \frac{|S_{ij}|^2}{S_{ii} S_{jj}}

The result is real-valued and bounded between 0 and 1.

.. code-block:: python
   :linenos:

   import numpy as np

   from pycsamt.emtools import coherence_matrix

   coh = coherence_matrix(sp1)

   print(coh.shape)
   print(np.nanmin(coh), np.nanmax(coh))
   print(np.diagonal(coh, axis1=1, axis2=2)[0])

The diagonal is 1 because each channel is perfectly coherent with
itself. The off-diagonal values are the useful ones for judging whether
two channels carry a stable relationship at a given frequency.

PSD Table
---------

``psd_table`` extracts the diagonal of the cross-power matrix as a tidy
table. It accepts one ``Spectra`` object, a list of spectra, or a
dictionary of station names to spectra.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import psd_table

   psd = psd_table(sp1)
   print(psd.head())

   psd_norm = psd_table({"spectra01": sp1, "spectra02": sp2}, normalize=True)
   print(psd_norm.groupby(["station", "channel"])["psd"].max().head())

Expected columns:

.. code-block:: text

   station, freq, period, channel, psd

Set ``normalize=True`` when you want to compare channel shapes rather
than absolute units. Electric and magnetic channels often have very
different physical units and raw PSD magnitudes.

Coherence Table
---------------

``coherence_table`` converts the coherence cube into a tidy table. By
default it includes all upper-triangle channel pairs. Pass ``pairs`` to
focus on physically meaningful pairs.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import coherence_table

   # Example channel-index pairs used by the bundled spectra example:
   # EX-HY and EY-HX.
   mt_pairs = [(3, 1), (4, 0)]

   coherence = coherence_table(sp1, pairs=mt_pairs)
   print(coherence.head())
   print(coherence.groupby("pair")["coherence"].describe())

Expected columns:

.. code-block:: text

   station, freq, period, ch_i, ch_j, pair, coherence

Use channel indices after checking ``sp.id_to_chtype``. Spectra files can
carry duplicate reference channels or project-specific channel ordering.

Coherence-Derived SNR
---------------------

``snr_table`` estimates SNR from squared coherence:

.. math::

   \mathrm{SNR} = \frac{\gamma^2}{1 - \gamma^2}

.. math::

   \mathrm{SNR}_{dB} = 10\log_{10}(\mathrm{SNR})

.. code-block:: python
   :linenos:

   from pycsamt.emtools import snr_table

   mt_pairs = [(3, 1), (4, 0)]

   snr = snr_table(sp1, pairs=mt_pairs)
   print(snr[["station", "freq", "pair", "coherence", "snr", "snr_db"]].head())
   print(snr.groupby("pair")["snr_db"].mean())

This is not the same table as the impedance-error SNR used in
``remove_noise``. In this spectra workflow, SNR is derived from channel
coherence before transfer functions are estimated.

Band Selection
--------------

``band_select`` returns a new ``Spectra`` object restricted to a
frequency interval. It slices all spectra arrays and metadata together.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import band_select

   high_band = band_select(sp1, f_min=100.0, f_max=10400.0)

   print(sp1.n_freq, high_band.n_freq)
   print(high_band.freq.min(), high_band.freq.max())

Use band selection when a frequency range is known to be more reliable
or when you need a common band across stations.

Coherence Masks
---------------

``mask_low_coherence`` returns a boolean mask over frequencies. It does
not modify the ``Spectra`` object. A ``True`` value means the frequency
passes the coherence criterion.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import mask_low_coherence

   mt_pairs = [(3, 1), (4, 0)]

   pass_any = mask_low_coherence(
       sp1,
       pairs=mt_pairs,
       threshold=0.5,
       require_all=False,
   )

   pass_all = mask_low_coherence(
       sp1,
       pairs=mt_pairs,
       threshold=0.5,
       require_all=True,
   )

   print(f"any pair passes: {pass_any.sum()} / {pass_any.size}")
   print(f"all pairs pass: {pass_all.sum()} / {pass_all.size}")

Use ``require_all=True`` when every requested channel pair must be
coherent before a frequency is accepted. Use ``False`` for a looser
screen where at least one pair is enough.

Spectra Summary
---------------

``spectra_summary`` produces one row per frequency. It combines
frequency metadata, channel PSD values, and mean off-diagonal coherence.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import spectra_summary

   summary = spectra_summary(sp1)
   print(summary.head())
   print(summary[["freq", "period", "mean_coherence"]].head())

Use this table for quick reporting and for finding frequency ranges
where average coherence collapses.

PSD Plot
--------

``plot_psd`` draws the auto-spectrum for selected channels.

.. code-block:: python
   :linenos:

   import matplotlib.pyplot as plt

   from pycsamt.emtools import plot_psd

   fig, ax = plt.subplots(figsize=(9, 5))
   plot_psd(
       sp1,
       channels=None,
       log_psd=True,
       title=f"{sp1.name} PSD",
       ax=ax,
   )

Pass a channel list when you want only a subset:

.. code-block:: python
   :linenos:

   plot_psd(sp1, channels=[0, 1, 3, 4], log_psd=True)

The x-axis follows the global pyCSAMT plotting control, which is usually
period-oriented for MT-style figures.

Coherence Plot
--------------

``plot_coherence`` creates one axis per channel pair and draws the
threshold line.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import plot_coherence

   mt_pairs = [(3, 1), (4, 0)]

   axes = plot_coherence(
       sp1,
       pairs=mt_pairs,
       threshold=0.5,
       show_threshold=True,
       title=f"{sp1.name} MT-pair coherence",
   )

Use this plot before deciding on a band cut. A single mean coherence
number can hide whether failures are isolated, broad-band, or confined
to one channel pair.

Full Spectral Matrix
--------------------

``plot_spectra_matrix`` visualizes the complete cross-power matrix at
one frequency. The diagonal cells are auto-spectra. Off-diagonal cells
are cross-spectra.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import plot_spectra_matrix

   fig = plot_spectra_matrix(
       sp1,
       freq_idx=0,
       quantity="abs",
       log_scale=True,
       title="Cross-power matrix",
   )

Available ``quantity`` values are ``"abs"``, ``"real"``, ``"imag"``, and
``"phase"``. The log absolute view is usually the best first look because
electric and magnetic channel powers can span many orders of magnitude.

Recovering Impedance From Spectra
---------------------------------

``plot_z_from_spectra`` calls ``Spectra.to_Z`` internally and plots
apparent resistivity and phase from the spectra-derived impedance.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import plot_z_from_spectra

   fig = plot_z_from_spectra(
       sp1,
       e_labels=("EX", "EY"),
       h_labels=("HX", "HY"),
       ridge=None,
       estimate_error=False,
       show_error=True,
   )

For programmatic access, call ``to_Z`` on the ``Spectra`` object:

.. code-block:: python
   :linenos:

   z_obj, tipper_obj = sp1.to_Z(
       e_labels=("EX", "EY"),
       h_labels=("HX", "HY"),
       estimate_error=False,
   )

   print(z_obj.z.shape)
   print(z_obj.resistivity[:, 0, 1].min(), z_obj.resistivity[:, 0, 1].max())

Use ``ridge`` when the magnetic cross-power block is poorly conditioned.
Keep ``estimate_error=False`` if degrees-of-freedom metadata is
incomplete and you do not need uncertainty envelopes.

Recovering Tipper From Spectra
------------------------------

``plot_tipper_from_spectra`` recovers induction tipper components from
the same spectra object when an ``HZ`` channel is available.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import plot_tipper_from_spectra

   axes = plot_tipper_from_spectra(
       sp1,
       h_labels=("HX", "HY"),
       ridge=None,
       estimate_error=False,
       show_error=True,
   )

If no vertical magnetic channel is available, the function returns axes
with a no-tipper message rather than failing.

Multiple-Station Sections
-------------------------

``plot_psd_section`` and ``plot_coherence_section`` accept a list or
dictionary of ``Spectra`` objects. They build a common log-frequency grid
over the overlapping frequency range and interpolate each station onto
that grid.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import plot_coherence_section, plot_psd_section

   spectra_sites = {
       "spectra01": sp1,
       "spectra02": sp2,
   }

   ax_psd = plot_psd_section(
       spectra_sites,
       channel=3,
       log_psd=True,
       title="EX PSD section",
   )

   ax_coh = plot_coherence_section(
       spectra_sites,
       pair=(3, 1),
       threshold=0.5,
       show_threshold=True,
       title="EX-HY coherence section",
   )

The section functions use only the shared frequency overlap. They do not
extrapolate outside a station's spectra band.

Practical QC Recipe
-------------------

A compact spectra QC sequence is:

.. code-block:: python
   :linenos:

   from pycsamt.emtools import (
       band_select,
       coherence_table,
       mask_low_coherence,
       spectra_summary,
   )

   mt_pairs = [(3, 1), (4, 0)]

   full_mask = mask_low_coherence(
       sp1,
       pairs=mt_pairs,
       threshold=0.5,
       require_all=True,
   )
   print(f"full band pass: {full_mask.sum()} / {full_mask.size}")

   clean = band_select(sp1, f_min=100.0, f_max=10400.0)
   clean_mask = mask_low_coherence(
       clean,
       pairs=mt_pairs,
       threshold=0.5,
       require_all=True,
   )
   print(f"selected band pass: {clean_mask.sum()} / {clean_mask.size}")

   coh = coherence_table(clean, pairs=mt_pairs)
   summary = spectra_summary(clean)

   print(coh.groupby("pair")["coherence"].describe())
   print(summary[["freq", "mean_coherence"]].head())

This sequence turns a visual coherence problem into a reproducible band
selection: define pairs, apply threshold, slice the spectra, and verify
the selected band.

Pitfalls
--------

Do not use these functions on impedance-only EDIs. They require a
``Spectra`` object with a real cross-power matrix.

Do not assume channel indices are universal. Always inspect
``sp.id_to_chtype`` before hard-coding pairs such as ``(3, 1)``.

Do not confuse spectra coherence SNR with impedance-error SNR. The
spectra ``snr_table`` is derived from squared coherence. The
noise-removal SNR table is based on impedance errors.

Do not compare raw PSD values across electric and magnetic channels as
if they shared units. Use log scaling or normalization when comparing
shape.

Worked Example
--------------

The example below loads the bundled de-identified spectra EDI files,
plots PSD and coherence, builds PSD/coherence/SNR tables, selects a
cleaner band, displays the full spectral matrix, recovers impedance and
tipper from spectra, and builds multi-station PSD/coherence sections.

.. literalinclude:: ../../../examples/emtools/plot_spectra.py
   :language: python
   :linenos:
