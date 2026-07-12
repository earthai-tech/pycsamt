r"""
Create ModEM/Occam-style run folders
====================================

At this point the inversion workspace contains corrected EDIs, data/error
tables, a starting model grid, and a validation report.  The next practical
step is to assemble **run folders** for one or more inversion backends.

This example creates two backend-style folders:

* ``modem_2d`` — a ModEM-style folder with impedance data, model arrays,
  a real ``ModEM.inv`` control file, and launch notes;
* ``occam2d`` — an Occam2D-style folder with rho/phase data, mesh/model
  arrays, startup/config notes, and launch notes.

The example does **not** launch ModEM or Occam.  Documentation builds should
not require external Fortran executables.  Instead, it prepares an auditable
folder that a user can inspect, adapt, version-control, and run manually.
"""

# %%
# 1. Imports and workspace paths
# ------------------------------

import csv
import json
import os
import shutil
import sys
from dataclasses import asdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def repo_root():
    root = os.environ.get("PYCSAMT_DOCS_REPO_ROOT")
    return Path(root) if root else Path(__file__).resolve().parents[3]


ROOT = repo_root()
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pycsamt.models.modem import ModEmConfig, ModEmControl
from pycsamt.models.occam2d import OccamConfig


workspace = Path(__file__).resolve().parent / "workspaces" / "l18_prepared_workspace"
table_dir = workspace / "02_tables"
model_dir = workspace / "03_model_placeholder"
run_root = workspace / "04_run_files"
figure_dir = workspace / "05_figures"

modem_dir = run_root / "modem_2d"
occam_dir = run_root / "occam2d"

for path in (modem_dir, occam_dir, figure_dir):
    path.mkdir(parents=True, exist_ok=True)

# %%
# 2. Load validated workspace artefacts
# -------------------------------------
# A run folder should be built from validated artefacts, not directly from raw
# field data.  If the validation report contains FAIL rows, stop here.

validation_report = table_dir / "inversion_input_validation_report.json"
if not validation_report.exists():
    raise RuntimeError(
        "Validation report missing. Run plot_4_validate_inversion_inputs.py first."
    )

report = json.loads(validation_report.read_text(encoding="utf-8"))
status_counts = report["status_counts"]
print("Validation status counts:")
print(json.dumps(status_counts, indent=2))

if status_counts.get("FAIL", 0):
    raise RuntimeError("Workspace validation has FAIL rows; do not create run folders.")

complex_table = table_dir / "inversion_impedance_complex_table.csv"
rho_phase_table = table_dir / "inversion_rho_phase_table.csv"
data_policy_file = table_dir / "inversion_data_error_policy.json"
model_policy_file = model_dir / "starting_model_policy.json"

complex_df = pd.read_csv(complex_table)
rho_phase_df = pd.read_csv(rho_phase_table)
data_policy = json.loads(data_policy_file.read_text(encoding="utf-8"))
model_policy = json.loads(model_policy_file.read_text(encoding="utf-8"))

print(f"Complex impedance rows: {len(complex_df)}")
print(f"Rho/phase rows: {len(rho_phase_df)}")

# %%
# 3. Common helpers
# -----------------


def copy_required(src, dst):
    src = Path(src)
    dst = Path(dst)
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    return dst


def write_text(path, text):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


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


def json_default(obj):
    if isinstance(obj, Path):
        return str(obj)
    return str(obj)


# %%
# 4. Build a ModEM-style 2-D run folder
# -------------------------------------
# ``ModEmControl`` writes a real ModEM ``.inv`` key-value control file.  The
# data and model files here are still lightweight CSV/JSON artefacts; a later
# backend writer can convert them to exact ModEM ASCII syntax.

modem_config = ModEmConfig(
    mode="2d",
    component_type="Off_Diagonal_Impedance",
    error_floor_z=data_policy["impedance_error"]["relative_z"],
    error_floor_z_floor=data_policy["impedance_error"]["absolute_z"],
    freq_min=data_policy["frequency_band_hz"]["fmin"],
    freq_max=data_policy["frequency_band_hz"]["fmax"],
    nx_2d=model_policy["horizontal"]["n_cells"],
    nz_2d=model_policy["vertical"]["n_layers"],
    cell_size_h_2d=model_policy["horizontal"]["inner_dx_m"],
    cell_size_v_top_2d=model_policy["vertical"]["first_dz_m"],
    depth_scale_2d=model_policy["vertical"]["growth"],
    n_padding_x_2d=model_policy["horizontal"]["padding_cells_each_side"],
    initial_rho=model_policy["resistivity"]["background_ohmm"],
    target_rms=1.05,
    max_iterations=80,
    output_stem="L18_modem2d",
)

modem_control = ModEmControl.from_config(modem_config)
modem_control_path = modem_control.write(modem_dir / modem_config.control_file)

modem_data = modem_dir / modem_config.data_file
modem_model = modem_dir / modem_config.model_file
modem_x_edges = modem_dir / "x_edges_m.csv"
modem_z_edges = modem_dir / "z_edges_m.csv"

copy_required(table_dir / "backend_modem_style_impedance_long.csv", modem_data)
copy_required(model_dir / "starting_resistivity_ohmm.csv", modem_model)
copy_required(model_dir / "x_edges_m.csv", modem_x_edges)
copy_required(model_dir / "z_edges_m.csv", modem_z_edges)
copy_required(validation_report, modem_dir / "validation_report.json")

modem_manifest = {
    "backend": "ModEM-style",
    "dimensionality": "2D",
    "status": "prepared_not_launched",
    "control_file": modem_config.control_file,
    "data_file": modem_config.data_file,
    "model_file": modem_config.model_file,
    "validation_report": "validation_report.json",
    "config": asdict(modem_config),
    "notes": [
        "Data file is a pyCSAMT long-table staging file, not final ModEM ASCII syntax.",
        "Use a backend writer or project adapter to convert before running ModEM.",
        "The ModEM.inv file is a real key-value control file created by ModEmControl.",
    ],
}

write_text(
    modem_dir / "run_manifest.json",
    json.dumps(modem_manifest, indent=2, default=json_default),
)

modem_command = (
    f"{modem_config.binary_name} {modem_config.control_file} "
    f"{modem_config.model_file} {modem_config.data_file}"
)

write_text(
    modem_dir / "RUN_NOTES.md",
    "\n".join(
        [
            "# ModEM-style run folder",
            "",
            "This folder is prepared but not launched.",
            "",
            "## Candidate command",
            "",
            "```bash",
            modem_command,
            "```",
            "",
            "Before running, convert the staging CSV data/model files to the exact",
            "ModEM ASCII format required by your executable and confirm impedance",
            "units/sign convention.",
            "",
        ]
    ),
)

print("ModEM-style run folder:", modem_dir)
print("ModEM control file:", modem_control_path)

# %%
# 5. Build an Occam2D-style run folder
# ------------------------------------
# Occam2D commonly works with apparent resistivity and phase rows.  We use the
# rho/phase long table and write a documented configuration snapshot.

occam_config = OccamConfig(
    modes=["TE", "TM"],
    error_floor_rho=data_policy["derived_error"]["logrho_min"],
    error_floor_phase=data_policy["derived_error"]["phase_min_deg"],
    freq_min=data_policy["frequency_band_hz"]["fmin"],
    freq_max=data_policy["frequency_band_hz"]["fmax"],
    n_layers=model_policy["vertical"]["n_layers"],
    cell_size_horizontal=model_policy["horizontal"]["inner_dx_m"],
    cell_size_vertical_top=model_policy["vertical"]["first_dz_m"],
    depth_scale=model_policy["vertical"]["growth"],
    n_padding_x=model_policy["horizontal"]["padding_cells_each_side"],
    initial_rho=model_policy["resistivity"]["background_ohmm"],
    max_iterations=80,
    target_misfit=1.0,
)

copy_required(table_dir / "backend_occam_style_rho_phase_long.csv", occam_dir / occam_config.data_file)
copy_required(model_dir / "x_edges_m.csv", occam_dir / occam_config.mesh_file)
copy_required(model_dir / "starting_log10_resistivity.csv", occam_dir / occam_config.model_file)
copy_required(model_dir / "z_edges_m.csv", occam_dir / "z_edges_m.csv")
copy_required(validation_report, occam_dir / "validation_report.json")

occam_startup_text = "\n".join(
    [
        "# Occam2D-style startup staging file",
        "# This is a readable pyCSAMT staging file, not final Fortran syntax.",
        f"data_file: {occam_config.data_file}",
        f"mesh_file: {occam_config.mesh_file}",
        f"model_file: {occam_config.model_file}",
        f"target_misfit: {occam_config.target_misfit}",
        f"max_iterations: {occam_config.max_iterations}",
        f"initial_rho: {occam_config.initial_rho}",
        f"lagrange_start: {occam_config.lagrange_start}",
        f"roughness_type: {occam_config.roughness_type}",
    ]
)
write_text(occam_dir / occam_config.startup_file, occam_startup_text)

occam_manifest = {
    "backend": "Occam2D-style",
    "status": "prepared_not_launched",
    "data_file": occam_config.data_file,
    "mesh_file": occam_config.mesh_file,
    "model_file": occam_config.model_file,
    "startup_file": occam_config.startup_file,
    "validation_report": "validation_report.json",
    "config": asdict(occam_config),
    "notes": [
        "Data file is a pyCSAMT rho/phase long-table staging file.",
        "Mesh/model/startup files are readable staging files.",
        "Use Occam2D InputBuilder or a project adapter for exact Occam2D syntax.",
    ],
}

write_text(
    occam_dir / "run_manifest.json",
    json.dumps(occam_manifest, indent=2, default=json_default),
)

occam_command = f"{occam_config.binary_name} {occam_config.startup_file}"
write_text(
    occam_dir / "RUN_NOTES.md",
    "\n".join(
        [
            "# Occam2D-style run folder",
            "",
            "This folder is prepared but not launched.",
            "",
            "## Candidate command",
            "",
            "```bash",
            occam_command,
            "```",
            "",
            "Before running, convert the staging files to exact Occam2D input",
            "syntax or rebuild them with `pycsamt.models.occam2d.InputBuilder`.",
            "",
        ]
    ),
)

print("Occam2D-style run folder:", occam_dir)

# %%
# 6. Run-folder audit
# -------------------
# Check that each folder contains the files a user expects before handing it
# to a solver-specific conversion or launch step.

expected = {
    "modem_2d": [
        modem_config.control_file,
        modem_config.data_file,
        modem_config.model_file,
        "x_edges_m.csv",
        "z_edges_m.csv",
        "validation_report.json",
        "run_manifest.json",
        "RUN_NOTES.md",
    ],
    "occam2d": [
        occam_config.data_file,
        occam_config.mesh_file,
        occam_config.model_file,
        occam_config.startup_file,
        "z_edges_m.csv",
        "validation_report.json",
        "run_manifest.json",
        "RUN_NOTES.md",
    ],
}

audit_rows = []
for folder_name, files in expected.items():
    folder = run_root / folder_name
    for filename in files:
        path = folder / filename
        audit_rows.append(
            {
                "folder": folder_name,
                "file": filename,
                "exists": path.exists(),
                "size_bytes": path.stat().st_size if path.exists() else 0,
            }
        )

audit_path = run_root / "run_folder_audit.csv"
write_csv(audit_path, audit_rows)

audit_df = pd.DataFrame(audit_rows)
print("Run-folder audit:")
print(audit_df.to_string(index=False))

if not audit_df["exists"].all():
    raise RuntimeError("One or more expected run-folder files are missing.")

# %%
# 7. Visual summary
# -----------------

fig, axs = plt.subplots(1, 2, figsize=(11.0, 4.2))

row_counts = [
    len(pd.read_csv(modem_dir / modem_config.data_file)),
    len(pd.read_csv(occam_dir / occam_config.data_file)),
]
axs[0].bar(["ModEM-style", "Occam2D-style"], row_counts, color=["#2563eb", "#16a34a"])
axs[0].set_ylabel("Data rows")
axs[0].set_title("Rows staged for each backend")
axs[0].grid(axis="y", alpha=0.25)

status = audit_df.groupby("folder")["exists"].mean()
axs[1].bar(status.index, status.values, color="#7c3aed")
axs[1].set_ylim(0, 1.05)
axs[1].set_ylabel("Completeness fraction")
axs[1].set_title("Run-folder completeness")
axs[1].grid(axis="y", alpha=0.25)

fig.tight_layout()
fig.savefig(figure_dir / "run_folder_creation_summary.png", dpi=160)

# %%
# 8. What is safe to do next?
# ---------------------------
# These folders are now organized, documented, and validated.  The next step
# is backend-specific conversion:
#
# * for ModEM, convert the impedance long table and model grid to exact ModEM
#   ASCII data/model/covariance syntax, then run the command in ``RUN_NOTES``;
# * for Occam2D, either convert the staged rho/phase table to ``OccamDataFile``
#   syntax or rebuild from EDIs using :class:`pycsamt.models.occam2d.InputBuilder`;
# * after a solver run, parse the log and plot RMS convergence before
#   interpreting the final resistivity section.

# sphinx_gallery_thumbnail_number = 1
