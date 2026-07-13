r"""
Write inversion data and error tables
=====================================

An inversion data file is more than a copy of the impedance tensor.  It is a
contract between processed data and the inversion engine:

* which stations and frequencies are used;
* which components are inverted;
* what value is supplied for each datum;
* what error controls the weight of that datum.

This example writes robust, solver-neutral inversion tables from corrected
EDIs.  It also writes backend-friendly long tables that can be translated to
ModEM, Occam2D, MARE2DEM, or another solver in a later step.

The central lesson is error handling.  If the EDI error is missing or too
small, we apply explicit floors so the inversion does not overfit individual
points.  Every row records which error source was used.
"""

# %%
# 1. Imports and paths
# --------------------

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


def repo_root():
    root = os.environ.get("PYCSAMT_DOCS_REPO_ROOT")
    return Path(root) if root else EXAMPLE_DIR.parents[2]


ROOT = repo_root()
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pycsamt.emtools import drop_duplicates, ensure_sites, select_band
from pycsamt.emtools._core import _get_z_block, _iter_items, _name

edi_dir = ROOT / "data" / "AMT" / "WILLY_DATA" / "L18PLT"
workspace = EXAMPLE_DIR / "workspaces" / "l18_prepared_workspace"
table_dir = workspace / "02_tables"
figure_dir = workspace / "05_figures"

for path in (table_dir, figure_dir):
    path.mkdir(parents=True, exist_ok=True)

# %%
# 2. Data policy
# --------------
# Define the inversion table policy before reading values.  The most important
# choices are the component list and the error floors.
#
# ``relative_z`` protects against unrealistically tiny impedance errors.
# ``absolute_z`` protects against zero-amplitude values.
# ``phase_min_deg`` prevents phase rows from receiving impossible precision.
# ``logrho_min`` is a lower bound for log10 apparent-resistivity errors.

policy = {
    "frequency_band_hz": {"fmin": 1e-3, "fmax": 1e3},
    "components": ["xy", "yx"],
    "impedance_error": {
        "relative_z": 0.05,
        "absolute_z": 1e-12,
    },
    "derived_error": {
        "logrho_min": 0.05,
        "phase_min_deg": 2.0,
    },
    "phase_wrap_degrees": [-180.0, 180.0],
    "units": {
        "frequency": "Hz",
        "period": "s",
        "impedance": "EDI/native impedance units",
        "apparent_resistivity": "ohm.m",
        "phase": "degree",
    },
}

print("Data/error table policy:")
print(json.dumps(policy, indent=2))

# %%
# 3. Load and normalize the line
# ------------------------------
# In production this input should be the corrected EDI folder prepared by the
# correction workflow.  Here we use the bundled line and apply the same
# duplicate/band normalization used in the workspace example.

sites = ensure_sites(edi_dir, recursive=False, verbose=0)
sites = drop_duplicates(sites, recursive=False)
sites = select_band(
    sites,
    fmin=policy["frequency_band_hz"]["fmin"],
    fmax=policy["frequency_band_hz"]["fmax"],
    recursive=False,
)

# %%
# 4. Helper functions
# -------------------
# Keep the conversion rules explicit.  Different inversion codes use different
# conventions, but these intermediate tables are transparent and easy to
# translate.

IJ = {"xx": (0, 0), "xy": (0, 1), "yx": (1, 0), "yy": (1, 1)}


def apparent_resistivity(z_value, frequency_hz):
    return 0.2 * abs(z_value) ** 2 / max(float(frequency_hz), 1e-24)


def phase_degrees(z_value):
    phase = float(np.degrees(np.angle(z_value)))
    lo, hi = policy["phase_wrap_degrees"]
    width = hi - lo
    return float(((phase - lo) % width) + lo)


def impedance_error_value(z_value, edi_error):
    rel = policy["impedance_error"]["relative_z"] * abs(z_value)
    abs_floor = policy["impedance_error"]["absolute_z"]
    floor = max(rel, abs_floor)
    edi = float(abs(edi_error)) if np.isfinite(edi_error) else np.nan
    if np.isfinite(edi) and edi >= floor:
        return edi, floor, "edi_error"
    if np.isfinite(edi):
        return floor, floor, "relative_floor"
    return floor, floor, "missing_edi_error"


def derived_errors(z_value, z_error):
    ratio = abs(z_error) / max(abs(z_value), 1e-24)
    logrho_err = max(
        2.0 * ratio / np.log(10.0), policy["derived_error"]["logrho_min"]
    )
    phase_err = max(
        np.degrees(ratio), policy["derived_error"]["phase_min_deg"]
    )
    return float(logrho_err), float(phase_err)


def write_csv(path, rows):
    rows = list(rows)
    if not rows:
        return
    with Path(path).open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


# %%
# 5. Build impedance, apparent-resistivity, and phase rows
# --------------------------------------------------------
# The impedance table is the primary table.  The rho/phase table is derived
# from it and is useful for 2-D MT codes that expect log apparent resistivity
# and phase rather than complex impedance.

impedance_rows = []
rho_phase_rows = []

for station_index, ed in enumerate(_iter_items(sites)):
    station = _name(ed, station_index)
    Z, z, fr = _get_z_block(ed)
    if Z is None or z is None or fr is None:
        continue

    z = np.asarray(z, dtype=np.complex128)
    fr = np.asarray(fr, dtype=float)
    z_err = getattr(Z, "z_err", None)
    z_err = np.asarray(z_err, dtype=float) if z_err is not None else None

    for frequency_index, frequency_hz in enumerate(fr):
        if not np.isfinite(frequency_hz) or frequency_hz <= 0:
            continue
        for component in policy["components"]:
            i, j = IJ[component]
            value = z[frequency_index, i, j]
            if not np.isfinite(value):
                continue
            edi_error = (
                float(z_err[frequency_index, i, j])
                if z_err is not None
                and z_err.shape == z.shape
                and np.isfinite(z_err[frequency_index, i, j])
                else np.nan
            )
            error_used, error_floor, source = impedance_error_value(
                value, edi_error
            )
            logrho_error, phase_error = derived_errors(value, error_used)
            rho = apparent_resistivity(value, frequency_hz)
            phase = phase_degrees(value)

            common = {
                "station": station,
                "station_index": int(station_index),
                "frequency_index": int(frequency_index),
                "frequency_hz": float(frequency_hz),
                "period_s": float(1.0 / frequency_hz),
                "component": component,
            }
            impedance_rows.append(
                {
                    **common,
                    "z_real": float(np.real(value)),
                    "z_imag": float(np.imag(value)),
                    "z_abs": float(abs(value)),
                    "edi_error": edi_error,
                    "error_floor": float(error_floor),
                    "error_used": float(error_used),
                    "error_source": source,
                    "relative_error_used": float(
                        error_used / max(abs(value), 1e-24)
                    ),
                }
            )
            rho_phase_rows.append(
                {
                    **common,
                    "rho_a_ohmm": float(rho),
                    "log10_rho_a": float(np.log10(max(rho, 1e-24))),
                    "log10_rho_error": float(logrho_error),
                    "phase_deg": float(phase),
                    "phase_error_deg": float(phase_error),
                    "derived_from": "complex_impedance",
                }
            )

print(f"Impedance rows: {len(impedance_rows)}")
print(f"Rho/phase rows: {len(rho_phase_rows)}")

# %%
# 6. Write solver-neutral tables
# ------------------------------
# These files preserve the most information.  Backend-specific writers should
# read these rather than going back to the EDIs.

impedance_table = table_dir / "inversion_impedance_complex_table.csv"
rho_phase_table = table_dir / "inversion_rho_phase_table.csv"
policy_file = table_dir / "inversion_data_error_policy.json"

write_csv(impedance_table, impedance_rows)
write_csv(rho_phase_table, rho_phase_rows)
policy_file.write_text(json.dumps(policy, indent=2), encoding="utf-8")

print("Complex impedance table:", impedance_table)
print("Rho/phase table:", rho_phase_table)
print("Policy file:", policy_file)

# %%
# 7. Write backend-friendly long tables
# -------------------------------------
# These are not strict official ModEM/Occam files yet.  They are deliberately
# simple long tables with names that match the concepts those solvers need.
#
# * The ModEM-style table keeps real and imaginary impedance values.
# * The Occam-style table stores TE/TM log-rho and phase rows.

mode_map = {
    "xy": ("TE", "ZXY"),
    "yx": ("TM", "ZYX"),
}

modem_rows = []
for row in impedance_rows:
    mode, modem_component = mode_map[row["component"]]
    modem_rows.append(
        {
            "station": row["station"],
            "period_s": row["period_s"],
            "mode": mode,
            "component": modem_component,
            "real": row["z_real"],
            "imag": row["z_imag"],
            "error": row["error_used"],
            "error_source": row["error_source"],
        }
    )

occam_rows = []
for row in rho_phase_rows:
    mode, _component = mode_map[row["component"]]
    occam_rows.append(
        {
            "station": row["station"],
            "period_s": row["period_s"],
            "mode": mode,
            "data_type": "log10_rho",
            "value": row["log10_rho_a"],
            "error": row["log10_rho_error"],
        }
    )
    occam_rows.append(
        {
            "station": row["station"],
            "period_s": row["period_s"],
            "mode": mode,
            "data_type": "phase_deg",
            "value": row["phase_deg"],
            "error": row["phase_error_deg"],
        }
    )

modem_table = table_dir / "backend_modem_style_impedance_long.csv"
occam_table = table_dir / "backend_occam_style_rho_phase_long.csv"
write_csv(modem_table, modem_rows)
write_csv(occam_table, occam_rows)

print("ModEM-style long table:", modem_table)
print("Occam-style long table:", occam_table)

# %%
# 8. Error audit summaries
# ------------------------
# Count which error rule controlled each row.  A healthy table usually uses a
# mix: real EDI errors where trustworthy, floors where errors are missing or
# unrealistically tiny.  If almost everything is floor-controlled, revisit the
# field errors or choose a more conservative floor.

sources = {}
for row in impedance_rows:
    sources[row["error_source"]] = sources.get(row["error_source"], 0) + 1

print("Error source counts:")
for key, value in sources.items():
    print(f"  {key}: {value}")

relative_errors = np.array(
    [row["relative_error_used"] for row in impedance_rows], dtype=float
)
logrho_errors = np.array(
    [row["log10_rho_error"] for row in rho_phase_rows], dtype=float
)
phase_errors = np.array(
    [row["phase_error_deg"] for row in rho_phase_rows], dtype=float
)

fig, axs = plt.subplots(1, 3, figsize=(13.0, 4.0))

axs[0].hist(relative_errors, bins=30, color="#2563eb", alpha=0.82)
axs[0].axvline(
    policy["impedance_error"]["relative_z"], color="#dc2626", lw=1.6
)
axs[0].set_xlabel("|Z| relative error used")
axs[0].set_ylabel("Rows")
axs[0].set_title("Complex impedance errors")

axs[1].hist(logrho_errors, bins=30, color="#16a34a", alpha=0.82)
axs[1].axvline(policy["derived_error"]["logrho_min"], color="#dc2626", lw=1.6)
axs[1].set_xlabel("log10 rho error")
axs[1].set_title("Derived rho errors")

axs[2].hist(phase_errors, bins=30, color="#7c3aed", alpha=0.82)
axs[2].axvline(
    policy["derived_error"]["phase_min_deg"], color="#dc2626", lw=1.6
)
axs[2].set_xlabel("phase error (deg)")
axs[2].set_title("Derived phase errors")

for ax in axs:
    ax.grid(axis="y", alpha=0.25)

fig.tight_layout()
fig.savefig(figure_dir / "inversion_error_floor_audit.png", dpi=160)

# %%
# 9. Data coverage audit
# ----------------------
# A data table can have good row count but poor distribution.  This heat map
# shows how many selected components exist for every station-period location.

stations = list(dict.fromkeys(row["station"] for row in impedance_rows))
periods = np.sort(np.unique([row["period_s"] for row in impedance_rows]))
coverage = np.zeros((len(periods), len(stations)), dtype=float)

station_index = {station: i for i, station in enumerate(stations)}
period_index = {period: i for i, period in enumerate(periods)}

for row in impedance_rows:
    coverage[
        period_index[row["period_s"]], station_index[row["station"]]
    ] += 1

fig, ax = plt.subplots(figsize=(11.0, 5.0))
im = ax.imshow(
    coverage,
    aspect="auto",
    origin="lower",
    interpolation="nearest",
    cmap="viridis",
    extent=(
        -0.5,
        len(stations) - 0.5,
        np.log10(periods.min()),
        np.log10(periods.max()),
    ),
)
ax.set_xticks(np.arange(len(stations)))
ax.set_xticklabels(stations, rotation=90, fontsize=7)
ax.set_ylabel(r"$\log_{10}T$ (s)")
ax.set_title("Inversion data coverage by station and period")
cbar = fig.colorbar(im, ax=ax, pad=0.02)
cbar.set_label("number of selected components")
fig.tight_layout()
fig.savefig(figure_dir / "inversion_data_coverage.png", dpi=160)

# %%
# 10. Practical checks before inversion
# -------------------------------------
# The table is now ready for a backend writer, but these checks are worth
# reading before you run any solver.

n_floor = sum(row["error_source"] != "edi_error" for row in impedance_rows)
floor_fraction = n_floor / max(len(impedance_rows), 1)
phase_bad = np.sum(
    np.abs([row["phase_deg"] for row in rho_phase_rows]) > 180.0
)

print(f"Floor-controlled impedance rows: {floor_fraction:.1%}")
print(f"Phase rows outside +/-180 deg: {phase_bad}")
print(f"Unique stations: {len(stations)}")
print(f"Unique periods: {len(periods)}")

if phase_bad:
    raise RuntimeError(
        "Phase wrapping failed; phase values exceed +/-180 degrees."
    )
if len(stations) < 3:
    raise RuntimeError(
        "Too few stations for a meaningful 2-D inversion table."
    )

# %%
# 11. What comes next?
# --------------------
# A later example can convert these long tables into exact backend syntax:
#
# * ModEM data file with station coordinates and impedance units;
# * Occam2D data file using TE/TM log-rho and phase rows;
# * MARE2DEM-style data blocks;
# * custom 1-D station-by-station inversion input.
#
# The important part is already done: each datum now has an explicit value,
# uncertainty, component, frequency, and audit trail.

# sphinx_gallery_thumbnail_number = 1
