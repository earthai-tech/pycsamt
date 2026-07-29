r"""
From IoT edge QC to an MT sounding: a complete acquisition pipeline
===================================================================

This is an end-to-end study that connects the two halves of pyCSAMT — the
IoT acquisition layer and the geophysical processing engine — on a single
**real** recording (``data/MT/TS/kap103as.ts``). The workflow is:

1. **acquire** the raw long-period MT time series;
2. run IoT **edge QC** (coverage, dropout, resolvable band) and reach an
   accept/warn decision;
3. **process** the accepted recording with :func:`pycsamt.ts.ts_to_z` to
   recover the impedance tensor, apparent resistivity, and phase;
4. read the classic **MT sounding** back, with the QC-resolvable band
   drawn on it so you can see which periods the data actually support; and
5. emit a **provenance** trail that ties the raw-file hash to the QC
   decision, the processing steps, and the written EDI; and
6. **sign** the manifest and use the data-model **bridge** to seed a
   re-occupation session straight from the EDI just written.

It is the "IoT acquisition feeds directly into the processing pipeline"
claim, made concrete on real data.
"""

# sphinx_gallery_thumbnail_number = 1

# %%
# 1. Acquire and QC the raw recording
# -----------------------------------
# The IoT edge layer inspects the raw series before any processing: how
# complete is it, and over what band is it resolvable?

from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from _ts_sample import load_ts_sample, real_ts_path

from pycsamt.iot import (
    EdgeProcessingConfig,
    EdgeProcessor,
    ProvenanceRecord,
    build_acquisition_manifest,
    detect_sensor_dropout,
    edi_survey_table,
    estimate_frequency_coverage,
    export_reproducibility_bundle,
    field_session_from_edis,
    verify_manifest,
)
from pycsamt.ts import ts_to_edi, ts_to_z

C_XY, C_YX = "#0072B2", "#D55E00"  # Okabe-Ito: two components, two hues
BAND_FILL = "#009E73"


def style_axis(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, which="both", color="#000000", alpha=0.07, lw=0.7)
    ax.set_axisbelow(True)


ts_path = real_ts_path()
record, ts_is_real = load_ts_sample()
if not ts_is_real:
    print("NOTE: bundled TS recording absent — using a synthetic sample.")
sample_rate = 1.0 / float(record.dt)
channels = [c.lower() for c in record.channels()]

block = np.column_stack(
    [np.asarray(record.get(c.upper()), float) for c in channels]
)
edge = EdgeProcessor(
    EdgeProcessingConfig(finite_threshold=0.9, warn_finite_threshold=0.999)
).process(block, channel_names=channels)
coverage_band = estimate_frequency_coverage(record.get("EX"), sample_rate)
dropout = detect_sensor_dropout(record.get("EX"))

print(
    f"station={record.station}  fs={sample_rate:g} Hz  n={record.n_samples:,}"
)
print(
    f"edge decision: {edge.decision.value}  "
    f"(coverage={edge.metrics['finite_coverage']:.4f}, "
    f"nan={dropout['nan_fraction']:.2%})"
)
print(
    f"resolvable band: {coverage_band.f_low_hz * 1e3:.2f}"
    f"-{coverage_band.f_high_hz * 1e3:.2f} mHz"
)

# %%
# 2. Process the accepted recording to impedance
# ----------------------------------------------
# Only because edge QC did not *reject* the recording do we spend compute
# on it. ``ts_to_z`` estimates band-averaged cross-spectra and recovers the
# impedance tensor, apparent resistivity, and phase.

if edge.decision.value == "reject":  # pragma: no cover - QC gate
    raise SystemExit("edge QC rejected the recording; skipping processing.")

z, tipper, spectra = ts_to_z(
    record, nfft=8192, per_decade=6, estimate_error=True
)
freq = np.asarray(z.freq, dtype=float)
period = 1.0 / freq
rho = {"xy": z.resistivity[:, 0, 1], "yx": z.resistivity[:, 1, 0]}
# yx phase is in the third quadrant; bring it to the first (standard MT
# convention, +180 deg) so both curves are comparable.
phase = {"xy": z.phase[:, 0, 1], "yx": z.phase[:, 1, 0] + 180.0}
print(
    f"processed {z.n_freq} frequencies  "
    f"({period.min():.1f}-{period.max():.0f} s period)"
)

# %%
# 3. The MT sounding, annotated with the QC-resolvable band
# ---------------------------------------------------------
# Apparent resistivity (log-log) and phase (semi-log) versus period, for
# the two off-diagonal components. The shaded span marks the band the edge
# QC deemed resolvable — a direct link from acquisition QC to the periods
# an interpreter should trust.

p_lo, p_hi = 1.0 / coverage_band.f_high_hz, 1.0 / coverage_band.f_low_hz

fig, (ax_rho, ax_phi) = plt.subplots(
    2,
    1,
    figsize=(8.5, 7.2),
    sharex=True,
    constrained_layout=True,
    gridspec_kw={"height_ratios": [2, 1]},
)
for comp, color in (("xy", C_XY), ("yx", C_YX)):
    good = np.isfinite(period) & np.isfinite(rho[comp]) & (rho[comp] > 0)
    ax_rho.loglog(
        period[good],
        rho[comp][good],
        "o-",
        ms=5,
        lw=1.6,
        color=color,
        label=rf"$\rho_{{{comp}}}$",
    )
    good_p = np.isfinite(period) & np.isfinite(phase[comp])
    ax_phi.semilogx(
        period[good_p],
        phase[comp][good_p],
        "o-",
        ms=5,
        lw=1.6,
        color=color,
        label=rf"$\phi_{{{comp}}}$",
    )

for ax in (ax_rho, ax_phi):
    ax.axvspan(p_lo, p_hi, color=BAND_FILL, alpha=0.10, lw=0)
    style_axis(ax)
ax_rho.set(
    ylabel=r"apparent resistivity ($\Omega\!\cdot\!$m)",
    title=f"MT sounding from IoT-QC'd time series - {record.station}",
)
ax_rho.legend(frameon=False, ncol=2, loc="best")
ax_phi.set(xlabel="period (s)", ylabel="phase (deg)", ylim=(0, 90))
ax_phi.set_yticks([0, 30, 45, 60, 90])
ax_phi.legend(frameon=False, ncol=2, loc="best")
ax_rho.annotate(
    "QC-resolvable band",
    xy=((p_lo * p_hi) ** 0.5, ax_rho.get_ylim()[1]),
    xytext=(0, -14),
    textcoords="offset points",
    ha="center",
    fontsize=9,
    color="#2f7d5b",
)

# %%
# 4. Write the EDI and a provenance bundle
# ----------------------------------------
# The processed result is written to a SEG EDI file, and a provenance
# manifest records the raw-file hash, the QC decision, the accepted band,
# and the exact processing steps — a reproducible chain from field bytes to
# interpreted curve.

out_dir = Path(tempfile.gettempdir()) / "kap103_acquisition"
out_dir.mkdir(exist_ok=True)
edi_path = ts_to_edi(
    record,
    out="kap103_from_ts.edi",
    savepath=str(out_dir),
    nfft=8192,
    per_decade=6,
)
print(
    f"wrote EDI: {Path(edi_path).name} ({os.path.getsize(edi_path):,} bytes)"
)

provenance = ProvenanceRecord(
    station_id=record.station,
    sample_rate_hz=sample_rate,
    lat=record.lat,
    lon=record.lon,
    accepted_band_hz=(coverage_band.f_low_hz, coverage_band.f_high_hz),
    processing_steps=[
        f"iot_edge_qc -> {edge.decision.value}",
        "ts_to_z(nfft=8192, per_decade=6)",
        "resistivity_phase",
        "ts_to_edi",
    ],
    field_notes="Real SAMTEX/LiMS long-period MT series.",
)
if ts_is_real:
    provenance.add_raw_file(str(ts_path))
manifest = build_acquisition_manifest(
    record.station,
    records=[provenance],
    method="mt",
)
bundle = export_reproducibility_bundle(manifest, str(out_dir / "provenance"))
print(f"provenance steps: {provenance.processing_steps}")
print(f"manifest content hash: {manifest.as_dict()['content_hash'][:16]}...")
print(
    f"bundle: {Path(bundle['manifest']).name} "
    f"+ {len(bundle['audits'])} audit(s)"
)

# %%
# 5. Sign the manifest and seed a re-occupation from the EDI
# ----------------------------------------------------------
# An HMAC signature makes the audit trail tamper-evident: a reviewer
# holding the shared key can confirm the manifest was not altered. Then the
# data-model bridge reads the EDI we just wrote back into an IoT
# ``FieldSession``, so the station's recorded geometry is ready to be
# re-occupied in a follow-up campaign — closing the loop from acquisition
# to data model and back.

key = "survey-key-2026"  # in practice, a per-survey secret
signed_path = manifest.write_signed(
    str(out_dir / "provenance" / "manifest.signed.json"), key
)
signed = json.loads(Path(signed_path).read_text())
print(
    f"signed manifest verifies: {verify_manifest(signed, key)}  "
    f"(wrong key: {verify_manifest(signed, 'nope')})"
)

reoccupy = field_session_from_edis(edi_path, survey_id="REOCCUPY", method="mt")
print(
    f"re-occupation session: {reoccupy.n_stations} station(s), "
    f"{reoccupy.n_devices} node(s)"
)
print(edi_survey_table(edi_path).to_string(index=False))

plt.show()
