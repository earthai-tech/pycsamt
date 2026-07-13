r"""
Validate inversion inputs before running
========================================

Before launching an inversion, stop and validate the workspace.  This is the
last inexpensive moment to catch mistakes such as:

* missing tables or stale files;
* non-finite impedance values;
* zero or unrealistically small errors;
* incomplete station-period coverage;
* impossible phase values;
* model grids that are not monotonic;
* starting-model arrays whose shape does not match the grid.

This example validates the solver-neutral workspace created by the previous
inversion examples.  The checks are intentionally explicit and readable so a
user can copy them into a project-specific pre-run gate.
"""

# %%
# 1. Imports and workspace paths
# ------------------------------

import csv
import json
import os
import sys
from pathlib import Path

# sphinx-gallery executes examples without __file__ (the gallery
# runner sets the working directory to this example's folder).
try:
    EXAMPLE_DIR = Path(__file__).resolve().parent
except NameError:
    EXAMPLE_DIR = Path.cwd()

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def repo_root():
    root = os.environ.get("PYCSAMT_DOCS_REPO_ROOT")
    return Path(root) if root else EXAMPLE_DIR.parents[2]


ROOT = repo_root()
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pycsamt.emtools import ensure_sites
from pycsamt.site import SitesReport

workspace = EXAMPLE_DIR / "workspaces" / "l18_prepared_workspace"
edi_dir = workspace / "01_corrected_edis"
table_dir = workspace / "02_tables"
model_dir = workspace / "03_model_placeholder"
figure_dir = workspace / "05_figures"

figure_dir.mkdir(parents=True, exist_ok=True)

# %%
# 2. Validation helpers
# ---------------------
# Each check returns a small dictionary with a status:
#
# * ``PASS`` — safe to proceed;
# * ``WARN`` — not necessarily wrong, but review before launching;
# * ``FAIL`` — do not run inversion until fixed.


def check(name, status, message, **details):
    return {
        "check": name,
        "status": status,
        "message": message,
        **details,
    }


def require_columns(frame, required):
    missing = [column for column in required if column not in frame.columns]
    return missing


def write_csv(path, rows):
    rows = list(rows)
    if not rows:
        return
    keys = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with Path(path).open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


results = []

# %%
# 3. File-existence checks
# ------------------------
# If required files are missing, the workspace is incomplete.  The user should
# rerun the previous preparation examples or rebuild the workspace manually.

required_files = {
    "policy": table_dir / "inversion_data_error_policy.json",
    "complex_table": table_dir / "inversion_impedance_complex_table.csv",
    "rho_phase_table": table_dir / "inversion_rho_phase_table.csv",
    "station_manifest": table_dir / "station_manifest.csv",
    "frequency_manifest": table_dir / "frequency_manifest.csv",
    "x_edges": model_dir / "x_edges_m.csv",
    "z_edges": model_dir / "z_edges_m.csv",
    "starting_model": model_dir / "starting_log10_resistivity.csv",
    "model_policy": model_dir / "starting_model_policy.json",
}

for label, path in required_files.items():
    results.append(
        check(
            f"file:{label}",
            "PASS" if path.exists() else "FAIL",
            f"{path.name} {'exists' if path.exists() else 'is missing'}",
            path=str(path),
        )
    )

if any(row["status"] == "FAIL" for row in results):
    print("One or more required files are missing.  See validation report.")

# %%
# 4. Load tables
# --------------

policy = json.loads(required_files["policy"].read_text(encoding="utf-8"))
complex_df = pd.read_csv(required_files["complex_table"])
rho_phase_df = pd.read_csv(required_files["rho_phase_table"])
station_df = pd.read_csv(required_files["station_manifest"])
frequency_df = pd.read_csv(required_files["frequency_manifest"])

print(f"Complex impedance rows: {len(complex_df)}")
print(f"Rho/phase rows: {len(rho_phase_df)}")
print(f"Stations in manifest: {len(station_df)}")
print(f"Frequencies in manifest: {len(frequency_df)}")

# %%
# 5. Schema and finite-value checks
# ---------------------------------
# These are hard failures for most inversion engines.  A table with missing
# columns or non-finite data should not be converted to a backend format.

complex_required = [
    "station",
    "frequency_hz",
    "period_s",
    "component",
    "z_real",
    "z_imag",
    "error_used",
    "error_source",
]
rho_phase_required = [
    "station",
    "frequency_hz",
    "period_s",
    "component",
    "log10_rho_a",
    "log10_rho_error",
    "phase_deg",
    "phase_error_deg",
]

missing_complex = require_columns(complex_df, complex_required)
missing_rho_phase = require_columns(rho_phase_df, rho_phase_required)

results.append(
    check(
        "schema:complex_impedance",
        "PASS" if not missing_complex else "FAIL",
        "complex impedance table has required columns",
        missing=";".join(missing_complex),
    )
)
results.append(
    check(
        "schema:rho_phase",
        "PASS" if not missing_rho_phase else "FAIL",
        "rho/phase table has required columns",
        missing=";".join(missing_rho_phase),
    )
)

for table_name, frame, numeric_columns in [
    (
        "complex_impedance",
        complex_df,
        ["frequency_hz", "period_s", "z_real", "z_imag", "error_used"],
    ),
    (
        "rho_phase",
        rho_phase_df,
        [
            "frequency_hz",
            "period_s",
            "log10_rho_a",
            "log10_rho_error",
            "phase_deg",
            "phase_error_deg",
        ],
    ),
]:
    finite = np.isfinite(frame[numeric_columns].to_numpy(dtype=float)).all()
    results.append(
        check(
            f"finite:{table_name}",
            "PASS" if finite else "FAIL",
            f"{table_name} numeric fields are finite",
        )
    )

# %%
# 6. Error-floor validation
# -------------------------
# Errors must be positive, finite, and not so small that the solver overfits.
# A high floor-controlled fraction is a warning: it may be correct, but it
# means the EDI errors were often missing or below the chosen floor.

positive_complex_error = (complex_df["error_used"] > 0).all()
positive_logrho_error = (rho_phase_df["log10_rho_error"] > 0).all()
positive_phase_error = (rho_phase_df["phase_error_deg"] > 0).all()

results.extend(
    [
        check(
            "error:complex_positive",
            "PASS" if positive_complex_error else "FAIL",
            "complex impedance errors are positive",
        ),
        check(
            "error:logrho_positive",
            "PASS" if positive_logrho_error else "FAIL",
            "log-rho errors are positive",
        ),
        check(
            "error:phase_positive",
            "PASS" if positive_phase_error else "FAIL",
            "phase errors are positive",
        ),
    ]
)

floor_fraction = float((complex_df["error_source"] != "edi_error").mean())
floor_status = (
    "PASS"
    if floor_fraction < 0.35
    else "WARN"
    if floor_fraction < 0.75
    else "FAIL"
)
results.append(
    check(
        "error:floor_fraction",
        floor_status,
        "fraction of impedance rows controlled by an error floor",
        value=floor_fraction,
    )
)

phase_min = float(policy["derived_error"]["phase_min_deg"])
logrho_min = float(policy["derived_error"]["logrho_min"])
phase_at_floor = float(
    (rho_phase_df["phase_error_deg"] <= phase_min + 1e-12).mean()
)
logrho_at_floor = float(
    (rho_phase_df["log10_rho_error"] <= logrho_min + 1e-12).mean()
)

results.append(
    check(
        "error:derived_floor_usage",
        "WARN" if max(phase_at_floor, logrho_at_floor) > 0.8 else "PASS",
        "derived rho/phase floor usage is not excessive",
        phase_at_floor=phase_at_floor,
        logrho_at_floor=logrho_at_floor,
    )
)

# %%
# 7. Phase and component validation
# ---------------------------------
# Phase values outside the expected wrapping interval are a hard failure.
# Missing components are usually a hard failure for a planned TE/TM inversion.

lo, hi = policy["phase_wrap_degrees"]
phase_ok = (
    (rho_phase_df["phase_deg"] >= lo) & (rho_phase_df["phase_deg"] <= hi)
).all()
results.append(
    check(
        "phase:range",
        "PASS" if phase_ok else "FAIL",
        f"phase values lie within [{lo}, {hi}] degrees",
    )
)

expected_components = set(policy["components"])
actual_components = set(complex_df["component"].astype(str))
missing_components = sorted(expected_components - actual_components)
results.append(
    check(
        "components:selected",
        "PASS" if not missing_components else "FAIL",
        "all requested components are present",
        missing=";".join(missing_components),
    )
)

# %%
# 8. Coverage validation
# ----------------------
# Count how many selected components exist at each station-period pair.  For
# this policy we expect two components: ``xy`` and ``yx``.

coverage = (
    complex_df.groupby(["station", "period_s"])["component"]
    .nunique()
    .reset_index(name="n_component")
)
expected_n_component = len(expected_components)
coverage_fraction = float(
    (coverage["n_component"] == expected_n_component).mean()
)
coverage_status = (
    "PASS"
    if coverage_fraction > 0.95
    else "WARN"
    if coverage_fraction > 0.80
    else "FAIL"
)

results.append(
    check(
        "coverage:station_period_components",
        coverage_status,
        "station-period rows contain the requested component set",
        coverage_fraction=coverage_fraction,
    )
)

n_station = complex_df["station"].nunique()
n_period = complex_df["period_s"].nunique()
results.extend(
    [
        check(
            "coverage:station_count",
            "PASS" if n_station >= 5 else "WARN",
            "enough stations are present for a line inversion",
            n_station=int(n_station),
        ),
        check(
            "coverage:period_count",
            "PASS" if n_period >= 8 else "WARN",
            "enough periods/frequencies are present",
            n_period=int(n_period),
        ),
    ]
)

# %%
# 9. EDI reload validation
# ------------------------
# The frozen EDI folder should still load.  This catches failed exports,
# accidental deletion, or writer problems before a solver sees the data.

try:
    reloaded = ensure_sites(edi_dir, recursive=False, verbose=0)
    reload_summary = SitesReport(reloaded).to_dataframe(api=False)
    reload_ok = len(reload_summary) == n_station
    results.append(
        check(
            "edi:reload",
            "PASS" if reload_ok else "FAIL",
            "workspace EDI folder reloads with expected station count",
            reloaded_station_count=int(len(reload_summary)),
            table_station_count=int(n_station),
        )
    )
except Exception as exc:
    results.append(
        check(
            "edi:reload",
            "FAIL",
            "workspace EDI folder could not be reloaded",
            error=str(exc),
        )
    )

# %%
# 10. Model-grid validation
# -------------------------
# The starting model must match the grid.  For cell-centred model arrays:
#
# * number of columns = len(x_edges) - 1;
# * number of rows = len(z_edges) - 1;
# * x and z edges must be strictly increasing.

x_edges = pd.read_csv(required_files["x_edges"])["x_edge_m"].to_numpy(
    dtype=float
)
z_edges = pd.read_csv(required_files["z_edges"])["z_edge_m"].to_numpy(
    dtype=float
)
model = np.loadtxt(required_files["starting_model"], delimiter=",")
model = np.atleast_2d(model)

x_ok = np.all(np.diff(x_edges) > 0)
z_ok = np.all(np.diff(z_edges) > 0)
shape_ok = model.shape == (len(z_edges) - 1, len(x_edges) - 1)
finite_model = np.isfinite(model).all()

results.extend(
    [
        check(
            "grid:x_monotonic",
            "PASS" if x_ok else "FAIL",
            "horizontal grid edges are strictly increasing",
        ),
        check(
            "grid:z_monotonic",
            "PASS" if z_ok else "FAIL",
            "vertical grid edges are strictly increasing",
        ),
        check(
            "model:shape",
            "PASS" if shape_ok else "FAIL",
            "starting model shape matches grid cell count",
            model_shape=str(model.shape),
            expected_shape=str((len(z_edges) - 1, len(x_edges) - 1)),
        ),
        check(
            "model:finite",
            "PASS" if finite_model else "FAIL",
            "starting model values are finite",
        ),
    ]
)

bottom_depth = float(z_edges[-1])
if "frequency_hz" in frequency_df and "coverage_fraction" in frequency_df:
    freq_coverage_ok = (frequency_df["coverage_fraction"] > 0.0).all()
else:
    freq_coverage_ok = False
results.append(
    check(
        "frequency_manifest:coverage",
        "PASS" if freq_coverage_ok else "WARN",
        "frequency manifest has positive station coverage for every row",
    )
)

results.append(
    check(
        "model:bottom_depth",
        "PASS" if bottom_depth > 500.0 else "WARN",
        "model bottom is deeper than a minimal demonstration depth",
        bottom_depth_m=bottom_depth,
    )
)

# %%
# 11. Visual validation dashboard
# -------------------------------
# A compact dashboard makes the validation report easier to review.

status_order = {"PASS": 0, "WARN": 1, "FAIL": 2}
status_counts = {
    status: sum(row["status"] == status for row in results)
    for status in ["PASS", "WARN", "FAIL"]
}

stations = list(dict.fromkeys(complex_df["station"].astype(str)))
periods = np.sort(complex_df["period_s"].unique())
matrix = np.full((len(periods), len(stations)), np.nan)
for _, row in coverage.iterrows():
    i = int(np.where(periods == row["period_s"])[0][0])
    j = stations.index(str(row["station"]))
    matrix[i, j] = row["n_component"]

fig, axs = plt.subplots(2, 2, figsize=(13.0, 8.0))

axs[0, 0].bar(
    status_counts.keys(),
    status_counts.values(),
    color=["#16a34a", "#f59e0b", "#dc2626"],
)
axs[0, 0].set_title("Validation statuses")
axs[0, 0].set_ylabel("Number of checks")
axs[0, 0].grid(axis="y", alpha=0.25)

axs[0, 1].hist(
    complex_df["relative_error_used"], bins=30, color="#2563eb", alpha=0.82
)
axs[0, 1].set_title("Relative impedance errors")
axs[0, 1].set_xlabel("error / |Z|")
axs[0, 1].grid(axis="y", alpha=0.25)

im = axs[1, 0].imshow(
    matrix,
    aspect="auto",
    origin="lower",
    interpolation="nearest",
    cmap="viridis",
    vmin=0,
    vmax=expected_n_component,
    extent=(
        -0.5,
        len(stations) - 0.5,
        np.log10(periods.min()),
        np.log10(periods.max()),
    ),
)
axs[1, 0].set_xticks(np.arange(len(stations)))
axs[1, 0].set_xticklabels(stations, rotation=90, fontsize=6)
axs[1, 0].set_ylabel(r"$\log_{10}T$ (s)")
axs[1, 0].set_title("Component coverage")
fig.colorbar(im, ax=axs[1, 0], pad=0.02, label="components")

axs[1, 1].plot(
    0.5 * (z_edges[:-1] + z_edges[1:]) / 1000.0,
    np.nanmedian(model, axis=1),
    "o-",
)
axs[1, 1].set_xlabel("Depth (km)")
axs[1, 1].set_ylabel(r"median $\log_{10}\rho$")
axs[1, 1].set_title("Starting model depth trend")
axs[1, 1].grid(alpha=0.25)

fig.tight_layout()
fig.savefig(figure_dir / "inversion_input_validation_dashboard.png", dpi=160)

# %%
# 12. Write validation reports
# ----------------------------

report_json = table_dir / "inversion_input_validation_report.json"
report_csv = table_dir / "inversion_input_validation_report.csv"

report = {
    "workspace": str(workspace),
    "status_counts": status_counts,
    "checks": results,
}

report_json.write_text(json.dumps(report, indent=2), encoding="utf-8")
write_csv(report_csv, results)

print("Validation status counts:")
print(json.dumps(status_counts, indent=2))
print("Validation report JSON:", report_json)
print("Validation report CSV:", report_csv)

if status_counts["FAIL"]:
    print("Validation failed.  Fix FAIL rows before running inversion.")
elif status_counts["WARN"]:
    print(
        "Validation passed with warnings.  Review WARN rows before running inversion."
    )
else:
    print(
        "Validation passed.  Inputs are ready for backend-specific conversion."
    )

# %%
# 13. How to use this in a real project
# -------------------------------------
# Use this page as a pre-run gate.  For a production inversion, make FAIL rows
# stop the run automatically.  WARN rows should be reviewed by a human and
# recorded in the project notes.  Only after this validation passes should you
# create backend-specific files and launch ModEM, Occam2D, MARE2DEM, or a
# custom solver.

# sphinx_gallery_thumbnail_number = 1
