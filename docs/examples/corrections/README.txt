.. _corrections_gallery:

Data corrections
----------------

Correcting the impedance tensor — before inversion or interpretation — is
one of pyCSAMT's defining strengths. The package ships a catalogue of
correction methods across six families; this section works through them as
**processing waves**, each an example that takes the raw survey line to a
cleaner one with publication-quality before/after figures.

The waves, applied to the bundled **WILLY_DATA L18PLT** line:

* **Static shift** — remove the frequency-independent ρ\ :sub:`a` offset with
  the AMA spatial average and the Hanning-EMAP filter;
* **Noise removal** — power-line notch, log-frequency smoothing, and robust
  ρ/φ trend smoothing;
* **Confidence editing** — keep, recover, mask, or drop station-frequency
  rows using auditable confidence thresholds;
* **Confidence-gated EMAP** — blend spatial filtering only where the
  confidence score says the data need help;
* **Source effects** — detect and correct CSAMT near-field / source
  overprint;
* **Tensor rotation** — rotate onto geoelectric strike and antisymmetrise
  for 2-D inversion;
* **Galvanic distortion** — estimate and audit a Groom--Bailey-style real
  distortion matrix before deciding whether to apply it;
* **A publication workflow** — chain the corrections into one final,
  inversion-ready dataset;
* **A step-by-step walkthrough** — take L18PLT and L22PLT from raw to
  *sanitised EDIs*, dropping bad frequencies and conditioning the survivors,
  with every step returning a new, cleaner dataset.
* **A pre-inversion case study** — start from collected EDIs, make explicit
  processing decisions, correct the line, export corrected EDIs, and reload
  them as a final handoff check.

Every wave uses the real ``pycsamt.emtools`` correction functions — the same
ones the desktop and web apps expose on their *Correction* page — so the
scripts double as a recipe for scripting the corrections yourself. Further
waves (coordinate-geometry conditioning, Stratagem instrument chains) will
be added over time. See the :doc:`emtools API reference </api/emtools>`.
