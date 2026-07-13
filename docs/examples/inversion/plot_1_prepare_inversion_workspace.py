r"""
Prepare an inversion workspace from corrected EDIs
=================================================

After correction, the next practical task is not usually "run inversion"
immediately.  A careful interpreter first creates an **inversion workspace**:
a folder that contains the corrected EDIs, station metadata, selected
frequency band, data/error tables, notes, and a machine-readable processing
policy.

This example shows that preparation step in a solver-neutral way.  The output
can later be adapted to ModEM, Occam2D, MARE2DEM, a custom 1-D inversion, or
any internal workflow.

The example uses the bundled WILLY ``L18PLT`` line so it can run during the
documentation build.  In a real project, replace ``corrected_edi_dir`` with
the folder produced by your correction workflow, for example the exported
folder from the pre-inversion correction case study.
"""

# %%
# 1. Imports and project paths
# ----------------------------
# Keep imports at the top so users can copy this page into a processing
# notebook or a field-project script.

import csv
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

# sphinx-gallery executes examples without __file__ (the gallery
# runner sets the working directory to this example's folder).
try:
    EXAMPLE_DIR = Path(__file__).resolve().parent
except NameError:
    EXAMPLE_DIR = Path.cwd()

import matplotlib.pyplot as plt
import numpy as np


def repo_root():
    root = os.environ.get("PYCSAMT_DOCS_REPO_ROOT")
    return Path(root) if root else EXAMPLE_DIR.parents[2]


ROOT = repo_root()
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pycsamt.emtools import drop_duplicates, ensure_sites, select_band
from pycsamt.emtools._core import _get_z_block, _iter_items, _name
from pycsamt.site import SitesReport, write_sites

# Replace this path with your corrected EDI folder.
corrected_edi_dir = ROOT / "data" / "AMT" / "WILLY_DATA" / "L18PLT"

# The workspace is deliberately inside the example folder for this gallery.
# In a real project this would be something like:
# ``PROJECT_ROOT / "inversion" / "L18_2026_07_prepared"``.
workspace = EXAMPLE_DIR / "workspaces" / "l18_prepared_workspace"

# %%
# 2. Decide the inversion contract before writing files
# -----------------------------------------------------
# A workspace is valuable because it records decisions.  These choices should
# be made explicitly, before exporting any inversion files:
#
# * What is the input data source?
# * Which frequency band is trusted?
# * Which tensor components will be inverted?
# * What error floor will be used when the EDI error is missing or too small?
# * Is topography included now, or only later during model/mesh preparation?

policy = {
    "project": "WILLY L18 demonstration",
    "line": "L18PLT",
    "source_type": "corrected_edi_folder",
    "source_path": str(corrected_edi_dir),
    "frequency_band_hz": {
        "fmin": 1e-3,
        "fmax": 1e3,
    },
    "components": ["xy", "yx"],
    "error_floor": {
        "relative_z": 0.05,
        "absolute_z": 1e-12,
        "comment": "Use max(EDI error, relative_z * |Z|, absolute_z).",
    },
    "topography": {
        "include_station_elevation": True,
        "mesh_topography": "defer_to_next_example",
    },
    "target_backend": "solver_neutral",
    "created_utc": datetime.now(timezone.utc).isoformat(),
}

print("Workspace policy:")
print(json.dumps(policy, indent=2))

# %%
# 3. Load and audit the corrected EDIs
# ------------------------------------
# The loader should find one impedance-bearing station per EDI.  If this
# audit fails, do not build inversion files yet; go back to correction/export.

sites = ensure_sites(corrected_edi_dir, recursive=False, verbose=0)
summary = SitesReport(sites).to_dataframe(api=False)

print(f"Stations loaded: {len(summary)}")
print(
    "Frequency rows per station: "
    f"{summary['nfreq'].min()}-{summary['nfreq'].max()}"
)
print("Station summary:")
print(
    summary[["station", "nfreq", "lat", "lon", "elev"]]
    .head(8)
    .to_string(index=False)
)

if summary.empty:
    raise RuntimeError(
        "No stations were loaded from the corrected EDI folder."
    )
if summary["nfreq"].min() <= 0:
    raise RuntimeError("At least one station has no valid frequency rows.")

# %%
# 4. Normalize the frequency band
# -------------------------------
# Frequency preparation for inversion is usually more conservative than for a
# quick plot.  We remove duplicates and select the intended band.  This step
# does not decide error floors or mesh design; it only makes the input line
# consistent enough for the workspace.

prepared = drop_duplicates(sites, recursive=False)
prepared = select_band(
    prepared,
    fmin=policy["frequency_band_hz"]["fmin"],
    fmax=policy["frequency_band_hz"]["fmax"],
    recursive=False,
)
prepared_summary = SitesReport(prepared).to_dataframe(api=False)

print(
    "Prepared frequency rows per station: "
    f"{prepared_summary['nfreq'].min()}-{prepared_summary['nfreq'].max()}"
)

# %%
# 5. Create the workspace layout
# ------------------------------
# A predictable layout makes handoff easier.  The external solver may require
# different filenames later, but these folders keep the preparation artefacts
# organized.

folders = {
    "root": workspace,
    "edi": workspace / "01_corrected_edis",
    "tables": workspace / "02_tables",
    "model": workspace / "03_model_placeholder",
    "run": workspace / "04_run_files",
    "figures": workspace / "05_figures",
    "notes": workspace / "06_notes",
}

for path in folders.values():
    path.mkdir(parents=True, exist_ok=True)

print("Workspace folders:")
for key, path in folders.items():
    print(f"  {key:>7}: {path}")

# %%
# 6. Export corrected EDIs into the workspace
# -------------------------------------------
# Even when the input folder already contains corrected EDIs, copying/exporting
# them into the inversion workspace freezes the exact dataset used for the
# inversion.  That matters when several correction experiments exist.

edi_manifest_path = folders["tables"] / "edi_manifest.csv"
edi_paths = write_sites(
    prepared,
    folders["edi"],
    template="{index:03d}_{station}.edi",
    exist_ok=True,
    manifest_csv=edi_manifest_path,
)

print(f"Exported EDI files: {len(edi_paths)}")
print("EDI manifest:", edi_manifest_path)

# %%
# 7. Build station and frequency manifests
# ----------------------------------------
# The station manifest is for humans and scripts: it records coordinates,
# elevation, and row counts.  The frequency manifest tells us which frequencies
# are common, sparse, or absent across the line.


def station_rows(sites_obj):
    rows = []
    for i, ed in enumerate(_iter_items(sites_obj)):
        _Z, z, fr = _get_z_block(ed)
        fr = np.asarray(fr, dtype=float) if fr is not None else np.array([])
        rows.append(
            {
                "station": _name(ed, i),
                "n_frequency": int(fr.size),
                "fmin_hz": float(np.nanmin(fr)) if fr.size else np.nan,
                "fmax_hz": float(np.nanmax(fr)) if fr.size else np.nan,
                "period_min_s": float(1.0 / np.nanmax(fr))
                if fr.size
                else np.nan,
                "period_max_s": float(1.0 / np.nanmin(fr))
                if fr.size
                else np.nan,
                "has_finite_z": bool(
                    z is not None and np.isfinite(np.asarray(z)).all()
                ),
            }
        )
    return rows


def frequency_rows(sites_obj):
    station_freqs = []
    for ed in _iter_items(sites_obj):
        _Z, _z, fr = _get_z_block(ed)
        if fr is not None:
            station_freqs.append(np.asarray(fr, dtype=float))
    if not station_freqs:
        return []
    all_freq = np.sort(np.unique(np.concatenate(station_freqs)))
    rows = []
    for freq in all_freq:
        present = sum(
            np.isclose(fr, freq, rtol=1e-6, atol=1e-12).any()
            for fr in station_freqs
        )
        rows.append(
            {
                "frequency_hz": float(freq),
                "period_s": float(1.0 / freq),
                "n_station_present": int(present),
                "coverage_fraction": float(present / len(station_freqs)),
            }
        )
    return rows


def write_csv(path, rows):
    rows = list(rows)
    if not rows:
        return
    with Path(path).open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


station_manifest = folders["tables"] / "station_manifest.csv"
frequency_manifest = folders["tables"] / "frequency_manifest.csv"

write_csv(station_manifest, station_rows(prepared))
write_csv(frequency_manifest, frequency_rows(prepared))

print("Station manifest:", station_manifest)
print("Frequency manifest:", frequency_manifest)

# %%
# 8. Write a solver-neutral impedance data table
# ----------------------------------------------
# This table is not a final ModEM/Occam file.  It is a clean intermediate
# contract.  A later example can translate it to a backend-specific format.
# Each row carries station, frequency, component, complex impedance, and the
# error value that should be used unless a backend requires a different
# convention.


def impedance_rows(sites_obj, components, relative_floor, absolute_floor):
    ij = {"xx": (0, 0), "xy": (0, 1), "yx": (1, 0), "yy": (1, 1)}
    rows = []
    for i, ed in enumerate(_iter_items(sites_obj)):
        station = _name(ed, i)
        Z, z, fr = _get_z_block(ed)
        if Z is None or z is None or fr is None:
            continue
        z = np.asarray(z, dtype=np.complex128)
        fr = np.asarray(fr, dtype=float)
        z_err = getattr(Z, "z_err", None)
        z_err = np.asarray(z_err, dtype=float) if z_err is not None else None
        for row_index, freq in enumerate(fr):
            if not np.isfinite(freq) or freq <= 0:
                continue
            for component in components:
                a, b = ij[component]
                value = z[row_index, a, b]
                if not np.isfinite(value):
                    continue
                edi_error = (
                    float(abs(z_err[row_index, a, b]))
                    if z_err is not None
                    and z_err.shape[:3] == z.shape
                    and np.isfinite(z_err[row_index, a, b])
                    else np.nan
                )
                floor = max(
                    float(relative_floor) * abs(value), float(absolute_floor)
                )
                error = max(
                    edi_error if np.isfinite(edi_error) else 0.0, floor
                )
                rows.append(
                    {
                        "station": station,
                        "frequency_hz": float(freq),
                        "period_s": float(1.0 / freq),
                        "component": component,
                        "z_real": float(np.real(value)),
                        "z_imag": float(np.imag(value)),
                        "z_abs": float(abs(value)),
                        "edi_error": edi_error,
                        "error_floor": float(floor),
                        "error_used": float(error),
                    }
                )
    return rows


data_table = folders["tables"] / "impedance_data_table.csv"
rows = impedance_rows(
    prepared,
    policy["components"],
    policy["error_floor"]["relative_z"],
    policy["error_floor"]["absolute_z"],
)
write_csv(data_table, rows)

print(f"Impedance rows written: {len(rows)}")
print("Data table:", data_table)

# %%
# 9. Write README, policy, and run notes
# --------------------------------------
# These files make the folder self-explanatory months later.  A future solver
# adapter can read ``inversion_policy.json`` directly.

policy_path = workspace / "inversion_policy.json"
policy_path.write_text(json.dumps(policy, indent=2), encoding="utf-8")

workspace_readme = workspace / "README.md"
workspace_readme.write_text(
    "\n".join(
        [
            "# pyCSAMT inversion workspace",
            "",
            "This folder was prepared by `plot_1_prepare_inversion_workspace.py`.",
            "",
            "## Contents",
            "",
            "- `01_corrected_edis/`: corrected EDI files frozen for inversion",
            "- `02_tables/`: station, frequency, EDI, and impedance manifests",
            "- `03_model_placeholder/`: starting model files will be added later",
            "- `04_run_files/`: backend-specific control files will be added later",
            "- `05_figures/`: QC figures produced during preparation",
            "- `06_notes/`: human-readable processing notes",
            "",
            "The workspace is solver-neutral.  Convert `impedance_data_table.csv`",
            "to ModEM, Occam2D, MARE2DEM, or another backend in a later step.",
            "",
        ]
    ),
    encoding="utf-8",
)

run_notes = folders["notes"] / "handoff_notes.md"
run_notes.write_text(
    "\n".join(
        [
            "# Inversion handoff notes",
            "",
            "Before running inversion:",
            "",
            "1. Confirm the frequency band is appropriate for the target depth.",
            "2. Confirm strike/dimensionality assumptions from diagnostics.",
            "3. Review station spacing before designing mesh cells.",
            "4. Choose backend-specific error floors and data modes.",
            "5. Record the external solver version and command line.",
            "",
            "This example has not created a solver-specific control file yet.",
            "",
        ]
    ),
    encoding="utf-8",
)

print("Policy JSON:", policy_path)
print("Workspace README:", workspace_readme)
print("Run notes:", run_notes)

# %%
# 10. Plot workspace diagnostics
# ------------------------------
# The first figure checks whether every station has enough rows after band
# selection.  The second figure shows frequency coverage across the line.

station_table = station_rows(prepared)
frequency_table = frequency_rows(prepared)

labels = [row["station"] for row in station_table]
nfreq = np.array([row["n_frequency"] for row in station_table], dtype=float)

fig, ax = plt.subplots(figsize=(10.5, 4.0))
ax.bar(np.arange(len(labels)), nfreq, color="#2563eb", alpha=0.85)
ax.set_xticks(np.arange(len(labels)))
ax.set_xticklabels(labels, rotation=90, fontsize=7)
ax.set_ylabel("Frequency rows")
ax.set_title("Prepared inversion workspace: station row counts")
ax.grid(axis="y", alpha=0.25)
fig.tight_layout()
fig.savefig(folders["figures"] / "station_frequency_counts.png", dpi=160)

freq = np.array([row["frequency_hz"] for row in frequency_table], dtype=float)
coverage = np.array(
    [row["coverage_fraction"] for row in frequency_table], dtype=float
)

fig, ax = plt.subplots(figsize=(8.5, 4.2))
ax.semilogx(freq, coverage, "o-", color="#16a34a")
ax.set_xlabel("Frequency (Hz)")
ax.set_ylabel("Station coverage fraction")
ax.set_ylim(0.0, 1.05)
ax.set_title("Prepared inversion workspace: frequency coverage")
ax.grid(True, which="both", alpha=0.25)
fig.tight_layout()
fig.savefig(folders["figures"] / "frequency_coverage.png", dpi=160)

# %%
# 11. Reload the workspace EDIs
# -----------------------------
# The final sanity check is simple but essential: reload the exported EDIs.
# If this fails, the inversion handoff is not ready.

reloaded = ensure_sites(folders["edi"], recursive=False, verbose=0)
reload_summary = SitesReport(reloaded).to_dataframe(api=False)

print(f"Reloaded stations: {len(reload_summary)}")
print(
    "Reloaded frequency rows per station: "
    f"{reload_summary['nfreq'].min()}-{reload_summary['nfreq'].max()}"
)

if len(reload_summary) != len(summary):
    raise RuntimeError(
        "Reloaded station count differs from the input station count."
    )

# %%
# 12. What comes next?
# --------------------
# This workspace is now ready for backend-specific preparation.  The next
# gallery examples can build on it:
#
# * design a starting model and depth grid;
# * convert ``impedance_data_table.csv`` to a ModEM/Occam data file;
# * validate error floors and dimensionality assumptions;
# * plot inversion convergence and final sections.

# sphinx_gallery_thumbnail_number = 1
