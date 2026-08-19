# Native Occam1D EDI demo

This example builds and inverts the three EDI soundings under `data/3edis`
using the native pyCSAMT Occam1D forward model, analytic Jacobian, bounded
regularized solver, restart checkpoints, and progress API.

Run from the repository root:

```bash
python examples/occam1_demo/run_demo.py
```

Run only one sounding:

```bash
python examples/occam1_demo/run_demo.py --station new_E1_1
```

Each station receives:

- `occam1d-inversion/`: native data, model, startup, manifest, and restart;
- `model-text/`: final model, response, convergence, audit, and metadata;
- `model-image/`: model, response, convergence, and summary PNG figures.

`output/demo-summary.json` is the machine-readable index of every product.

## Customizing the figures

Occam1D figures use the shared API in `pycsamt.api.occam1d`. The default
model/layer legend is intentionally hidden. Every major artist can be changed
or disabled without altering the inversion:

```python
from pycsamt.api import PYCSAMT_OCCAM1D

custom = PYCSAMT_OCCAM1D.default.copy(
    observed__marker="x",
    observed__color="black",
    predicted__color="purple",
    phase_predicted__linestyle="--",
    model__color="darkgreen",
    model__linewidth=3.0,
    iteration__marker="s",
    iteration__color="navy",
    target__color="crimson",
    target__linestyle=":",
    model_legend=False,
    response_legend=True,
)

inversion.save_main_images("model-image", style=custom)
```

Use `target__visible=False`, `roughness__visible=False`, or
`observed__visible=False` to hide individual layers. Named presets are
`"pycsamt"`, `"publication"`, and `"minimal"`.
