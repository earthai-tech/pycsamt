# AI/DL inversion -> PCSF: a shareable, reproducible result format

Turns a bare AI/DL inversion prediction — a UNet, a GCN, a ResNet, or any
third party's own model — into a citable, reproducible `.pcsf`/`.pcsm`
file, the same backend-neutral pyCSAMT Common Subsurface Format
`examples/pcsf_conversion_demo` produces from real Occam2D/ModEM/MARE2DEM
runs. See `PYCSAMT-PCSF-INVERSION-FORMAT-PLAN.md` (Phase 9) and
`pycsamt/format/SPEC.md` (section 6) for the full design.

Unlike every other adapter in `pycsamt.format.adapters`, the functions used
here (`pycsamt.format.adapters.generic`) need **no solver result object at
all** — only plain `numpy` arrays. That is the point: a third-party AI/DL
tool with zero dependency on pycsamt's own `pycsamt.ai` subpackage can
still publish a result others can reload, inspect, and independently
verify.

Run from the repository root:

```bash
python examples/ai_inversion_pcsf/run_demo.py
```

## What this demonstrates

| Model type | Function used | PCSF geometry kind | Encoding |
|---|---|---|---|
| Toy **UNet** (regular 2-D section) | `grid2d_to_pcsf` | `grid2d` | `log10` -> canonical linear |
| Toy **GCN** (per-node mesh prediction) | `mesh_to_pcsf` | `mesh_unstructured` | `log10` -> canonical linear |

Both predictions are synthetic, seeded data — this demo is about the
**file format's** reproducibility story, not a real trained model or a
real survey.

### The GCN case: per-node resistivity

A graph model naturally predicts one value per mesh vertex, not per cell.
`mesh_to_pcsf`'s `resistivity_by_node` parameter accepts that directly and
projects it onto PCSF's canonical per-triangle `resistivity` as the
documented arithmetic mean of each triangle's 3 vertex values — both the
raw per-node prediction and the derived per-cell array are kept in the
file (`model/resistivity_by_node` and `model/resistivity`), so nothing is
lost in the projection.

### Model provenance: what makes this actually reproducible

Reloading the array back is only half of "reproducible" — the other half
is knowing *what produced it*. Each file carries a
`pycsamt.format.provenance.ModelProvenance` block (architecture, framework
+ version, hyperparameters, random seed, and a **checkpoint SHA-256**),
stored at `metadata['model_provenance']`:

```python
from pycsamt.format import read_pcsf

model = read_pcsf("examples/ai_inversion_pcsf/output/toy_unet_grid2d.pcsf")
model.metadata["model_provenance"]
# {'architecture': 'UNet', 'framework': 'pytorch', 'framework_version': '2.3.0',
#  'checkpoint': 'examples/ai_inversion_pcsf/output/toy_unet_checkpoint.bin',
#  'checkpoint_sha256': '...', 'hyperparameters': {'lr': 0.001, ...}, ...}
```

This demo's "checkpoint" is a small deterministic toy binary file (not a
real trained model), but the hash is real: anyone who downloads a
published checkpoint alongside a shared `.pcsf` file can re-hash it with
`pycsamt.format.provenance.compute_checkpoint_hash` and confirm it matches
`checkpoint_sha256` — a concrete, checkable claim rather than just a
filename.

## What gets written

```
output/
  toy_unet_grid2d.pcsf         grid2d, log10-encoded UNet-style prediction
  toy_unet_grid2d.pcsm         same model, hand-editable ASCII sibling
  toy_gcn_mesh.pcsf            mesh_unstructured, per-node GCN-style prediction
  toy_gcn_mesh.pcsm            same model, hand-editable ASCII sibling
  toy_unet_checkpoint.bin      toy "trained weights" (hashed, not real weights)
  toy_gcn_checkpoint.bin       toy "trained weights" (hashed, not real weights)
  demo-summary.json            machine-readable file-by-file summary
```

## Reading a file back — no pycsamt.ai import needed

```python
from pycsamt.format import read_pcsf

model = read_pcsf("examples/ai_inversion_pcsf/output/toy_gcn_mesh.pcsf")
model.kind                       # "mesh_unstructured"
model.resistivity.shape          # per-triangle, linear ohm.m
model.resistivity_by_node.shape  # the model's own native per-node output
model.metadata["model_provenance"]["architecture"]  # "GCN"
```

Or with `h5py` directly, no pycsamt import at all:

```python
import h5py
with h5py.File("examples/ai_inversion_pcsf/output/toy_gcn_mesh.pcsf") as f:
    f.visititems(lambda name, obj: print(name, getattr(obj, "shape", "")))
```

## Extending this to a real model

Swap the synthetic array-building in `run_demo.py` for your own model's
real output and a real checkpoint path:

```python
from pycsamt.format.adapters.generic import grid2d_to_pcsf
from pycsamt.format.provenance import ModelProvenance, compute_checkpoint_hash
from pycsamt.format import write_pcsf

provenance = ModelProvenance(
    architecture="ResNet18",
    framework="tensorflow", framework_version="2.16.1",
    checkpoint="resnet18_mt_v3.h5",
    checkpoint_sha256=compute_checkpoint_hash("resnet18_mt_v3.h5"),
    random_seed=123,
)
model = grid2d_to_pcsf(
    my_predicted_log10_rho, x=my_x, z=my_z, encoding="log10",
    provenance=provenance, source_backend="resnet",
)
write_pcsf(model, "my_result.pcsf")
```

Anything with a numpy array output fits one of `grid2d_to_pcsf`,
`grid3d_to_pcsf`, or `mesh_to_pcsf` — see
`pycsamt/format/adapters/generic.py` for the full parameter set
(uncertainty/sensitivity arrays, station elevation/lon-lat via the same
smart `topo=` resolver every other adapter uses, and `origin`/rotation for
real-world placement).
