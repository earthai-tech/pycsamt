"""Publish an AI/DL inversion result as a citable, reproducible PCSF file.

Every other PCSF example (`examples/pcsf_conversion_demo`) converts a
*specific solver's* in-memory result object (Occam2D's ``InversionResult``,
ModEM's ``ModEmModel3D``, MARE2DEM's ``TriMesh``). This demo needs none of
that: it builds two synthetic "AI predictions" from plain ``numpy`` arrays
via :mod:`pycsamt.format.adapters.generic`, exactly the way a UNet, a GCN,
a ResNet, or any third-party AI/DL inversion tool -- with zero dependency
on pycsamt's own :mod:`pycsamt.ai` subpackage -- would do it:

* A toy **UNet**-style prediction: a regular 2-D resistivity section
  (log10-encoded, the common training-stability choice for DL models),
  converted with :func:`pycsamt.format.adapters.generic.grid2d_to_pcsf`.
* A toy **GCN**-style prediction: per-*node* resistivity on a small
  triangular mesh -- the natural output shape of a graph model -- converted
  with :func:`pycsamt.format.adapters.generic.mesh_to_pcsf`, which projects
  it onto PCSF's canonical per-triangle ``resistivity`` (documented
  arithmetic mean of each triangle's 3 vertex values).

Both predictions are pure synthetic data (seeded, reproducible), not a
real trained model -- this demo is about the *file format's* reproducibility
story, not a scientific inversion result. Each carries a
:class:`pycsamt.format.provenance.ModelProvenance` block recording a toy
"checkpoint" file's real SHA-256 hash, so a reviewer downloading the
(illustrative) checkpoint alongside a shared ``.pcsf`` could verify they
have the exact weights the result claims -- the actual point of a
community-shareable AI-inversion format.

Run from the repository root::

    python examples/ai_inversion_pcsf/run_demo.py
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from pycsamt.format import read_pcsf, read_pcsm, write_pcsf, write_pcsm
from pycsamt.format.adapters.generic import grid2d_to_pcsf, mesh_to_pcsf
from pycsamt.format.provenance import ModelProvenance, compute_checkpoint_hash

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).resolve().parent / "output"


def _print_pcsf_tree(path: Path) -> None:
    import h5py

    print(f"\n--- {path.relative_to(ROOT)} ({path.stat().st_size:,} bytes) ---")

    def _visit(name, obj):
        if isinstance(obj, h5py.Dataset):
            print(f"  {name}  dataset  shape={obj.shape}  dtype={obj.dtype}")
        else:
            print(f"  {name}/")

    with h5py.File(path, "r") as fh:
        print("  [root attrs]", dict(fh.attrs))
        fh.visititems(_visit)


def _write_toy_checkpoint(name: str, seed: int) -> Path:
    """A small deterministic binary file standing in for real trained
    weights -- just enough to give :func:`compute_checkpoint_hash` a real
    file to hash, so ``checkpoint_sha256`` below is a genuinely
    independently-verifiable value, not a placeholder string."""
    path = OUTPUT / name
    rng = np.random.default_rng(seed)
    path.write_bytes(rng.bytes(4096))
    return path


# ---------------------------------------------------------------------
# Toy UNet -> grid2d
# ---------------------------------------------------------------------


def build_toy_unet_grid2d() -> dict:
    print("\n=== Toy UNet prediction -> grid2d ===")
    rng = np.random.default_rng(42)

    n_x, n_z = 24, 14
    x = np.linspace(0.0, 460.0, n_x)
    z = np.linspace(5.0, 280.0, n_z)

    # A plausible layered-earth log10(rho) section a UNet trained on
    # inversion sections might predict: conductive near-surface layer,
    # resistive basement, plus small per-cell prediction noise.
    depth_trend = 1.3 + 1.1 * (1.0 - np.exp(-z / 90.0))
    log10_rho = depth_trend[:, None] * np.ones((n_z, n_x))
    log10_rho += rng.normal(scale=0.05, size=(n_z, n_x))

    station_names = [f"AI{i:02d}" for i in range(n_x)]
    station_elevations = {
        name: float(1200.0 + 3.0 * np.sin(i / 3.0))
        for i, name in enumerate(station_names)
    }
    station_lonlat = {
        name: (11.0 + 0.001 * i, 46.0 + 0.0005 * i)
        for i, name in enumerate(station_names)
    }

    checkpoint_path = _write_toy_checkpoint("toy_unet_checkpoint.bin", seed=1)
    provenance = ModelProvenance(
        architecture="UNet",
        framework="pytorch",
        framework_version="2.3.0",
        checkpoint=str(checkpoint_path.relative_to(ROOT)),
        checkpoint_sha256=compute_checkpoint_hash(checkpoint_path),
        training_data="synthetic_layered_earth_v1 (illustrative, not a real survey)",
        hyperparameters={"lr": 1e-3, "epochs": 50, "batch_size": 16},
        random_seed=42,
        authors=["pycsamt example"],
        notes="Synthetic demo prediction -- not a real trained model.",
    )
    history = {
        "epoch": list(range(6)),
        "train_loss": [0.42, 0.31, 0.22, 0.17, 0.14, 0.12],
        "val_loss": [0.46, 0.35, 0.27, 0.22, 0.19, 0.18],
    }

    model = grid2d_to_pcsf(
        log10_rho, x=x, z=z, encoding="log10",
        station_names=station_names, station_x=list(x),
        station_elevations=station_elevations, station_lonlat=station_lonlat,
        history=history, provenance=provenance,
        source_backend="unet",
        created_by="ai_inversion_pcsf demo",
        description="Toy UNet log10(rho) prediction -- synthetic layered earth",
    )
    model.validate()

    pcsf_path = write_pcsf(model, OUTPUT / "toy_unet_grid2d.pcsf")
    pcsm_path = write_pcsm(model, OUTPUT / "toy_unet_grid2d.pcsm")
    print(f"wrote {pcsf_path.relative_to(ROOT)} and {pcsm_path.relative_to(ROOT)}")

    restored_pcsf = read_pcsf(pcsf_path)
    restored_pcsm = read_pcsm(pcsm_path)
    np.testing.assert_allclose(restored_pcsf.resistivity, model.resistivity)
    np.testing.assert_allclose(restored_pcsm.resistivity, model.resistivity)
    assert (
        restored_pcsf.metadata["model_provenance"]["checkpoint_sha256"]
        == provenance.checkpoint_sha256
    )
    print(
        "round-trip check: resistivity bit-exact through .pcsf and .pcsm; "
        "model_provenance survives both encodings"
    )
    print(f"checkpoint_sha256 (independently re-computable): {provenance.checkpoint_sha256}")

    _print_pcsf_tree(pcsf_path)
    return _summarize(pcsf_path, model)


# ---------------------------------------------------------------------
# Toy GCN -> mesh_unstructured
# ---------------------------------------------------------------------


def _structured_triangle_mesh(
    n_x: int, n_z: int, dx: float, dz: float
) -> tuple[np.ndarray, np.ndarray]:
    """A plain structured grid, split into 2 triangles per cell.

    Deterministic and dependency-free (no ``triangle`` package needed) --
    good enough to demonstrate the mesh/GCN path without pretending it is
    a real unstructured survey mesh.
    """
    xs = np.arange(n_x) * dx
    zs = np.arange(n_z) * dz
    xv, zv = np.meshgrid(xs, zs)
    nodes = np.column_stack([xv.ravel(), zv.ravel()])

    def node_index(i: int, j: int) -> int:
        return j * n_x + i

    triangles = []
    for j in range(n_z - 1):
        for i in range(n_x - 1):
            a = node_index(i, j)
            b = node_index(i + 1, j)
            c = node_index(i, j + 1)
            d = node_index(i + 1, j + 1)
            triangles.append([a, b, d])
            triangles.append([a, d, c])
    return nodes, np.asarray(triangles, dtype=np.int64)


def build_toy_gcn_mesh() -> dict:
    print("\n=== Toy GCN prediction -> mesh_unstructured ===")
    rng = np.random.default_rng(7)

    nodes, connectivity = _structured_triangle_mesh(n_x=9, n_z=7, dx=20.0, dz=15.0)
    depth = nodes[:, 1]
    # A GCN's natural output shape: one value per graph node/mesh vertex.
    node_log10_rho = 1.5 + 0.01 * depth + rng.normal(scale=0.03, size=depth.shape)

    checkpoint_path = _write_toy_checkpoint("toy_gcn_checkpoint.bin", seed=2)
    provenance = ModelProvenance(
        architecture="GCN",
        framework="pytorch-geometric",
        framework_version="2.5.2",
        checkpoint=str(checkpoint_path.relative_to(ROOT)),
        checkpoint_sha256=compute_checkpoint_hash(checkpoint_path),
        training_data="synthetic_mesh_earth_v1 (illustrative, not a real survey)",
        hyperparameters={"lr": 5e-4, "layers": 4, "hidden_dim": 64},
        random_seed=7,
        authors=["pycsamt example"],
        notes="Synthetic demo prediction -- not a real trained model.",
    )

    model = mesh_to_pcsf(
        nodes, connectivity, resistivity_by_node=node_log10_rho,
        encoding="log10", plane="xz", provenance=provenance,
        source_backend="gcn",
        created_by="ai_inversion_pcsf demo",
        description="Toy GCN per-node log10(rho) prediction on a small mesh",
    )
    model.validate()

    pcsf_path = write_pcsf(model, OUTPUT / "toy_gcn_mesh.pcsf")
    pcsm_path = write_pcsm(model, OUTPUT / "toy_gcn_mesh.pcsm")
    print(f"wrote {pcsf_path.relative_to(ROOT)} and {pcsm_path.relative_to(ROOT)}")

    restored_pcsf = read_pcsf(pcsf_path)
    restored_pcsm = read_pcsm(pcsm_path)
    np.testing.assert_allclose(restored_pcsf.resistivity_by_node, model.resistivity_by_node)
    np.testing.assert_allclose(restored_pcsm.resistivity_by_node, model.resistivity_by_node)
    np.testing.assert_allclose(
        restored_pcsf.resistivity,
        model.resistivity_by_node[connectivity].mean(axis=1),
    )
    print(
        "round-trip check: resistivity_by_node bit-exact through .pcsf and "
        ".pcsm; per-triangle resistivity matches the documented "
        "arithmetic-mean-of-3-vertices projection"
    )

    _print_pcsf_tree(pcsf_path)
    return _summarize(pcsf_path, model)


# ---------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------


def _summarize(path: Path, model) -> dict:
    return {
        "file": str(path.relative_to(ROOT)),
        "size_bytes": path.stat().st_size,
        "kind": model.kind,
        "source_backend": model.source_backend,
        "resistivity_shape": list(model.resistivity.shape),
        "has_resistivity_by_node": model.resistivity_by_node is not None,
        "n_stations": len(model.stations.name) if model.stations is not None else 0,
        "model_provenance": model.metadata.get("model_provenance"),
    }


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)

    records = [build_toy_unet_grid2d(), build_toy_gcn_mesh()]

    summary = {
        "schema": "pycsamt.ai_inversion_pcsf.demo/v1",
        "plan": "PYCSAMT-PCSF-INVERSION-FORMAT-PLAN.md (Phase 9)",
        "note": (
            "Both files are built with pycsamt.format.adapters.generic "
            "from plain numpy arrays -- no pycsamt.ai or solver result "
            "object involved -- the same path a third-party AI/DL "
            "inversion tool would use."
        ),
        "files": records,
    }
    summary_path = OUTPUT / "demo-summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf8"
    )
    print(f"\nDemo completed: {summary_path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
