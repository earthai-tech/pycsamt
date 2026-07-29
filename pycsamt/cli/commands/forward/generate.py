# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt forward generate — batch-generate synthetic (model, response) pairs.

Wraps :func:`~pycsamt.forward.batch.generate_dataset` with a convenient
CLI interface for producing training datasets without writing any Python.

Usage
-----
::

    # quick run — 1 000 MT samples, Gaussian noise, saved as .npz
    pycsamt forward generate --n-samples 1000 --output mt1d_train.npz

    # large dataset, field-realistic noise, parallel
    pycsamt forward generate \\
        --solver mt1d \\
        --n-samples 50000 \\
        --n-layers 3:8 \\
        --freqs 0.001:10000 \\
        --n-freqs 30 \\
        --noise-type field \\
        --noise-level 0.03 \\
        --seed 42 \\
        --jobs 4 \\
        --output mt1d_field_noise.npz

    # TEM dataset
    pycsamt forward generate \\
        --solver tem1d --n-samples 2000 \\
        --noise-type multiplicative \\
        --output tem1d_train.npz
"""

from __future__ import annotations

import sys

import click

from ....api.cli.config import configure_cli
from ....api.cli.options import (
    n_jobs_option,
    no_color_option,
    verbose_option,
)
from ....api.cli.params import FreqRange
from ._base import forward

# ---------------------------------------------------------------------------
# generate command
# ---------------------------------------------------------------------------


@forward.command("generate")
@click.option(
    "--solver",
    type=click.Choice(["mt1d", "tem1d", "csamt1d"], case_sensitive=False),
    default="mt1d",
    show_default=True,
    help="1-D EM solver used to compute responses.",
)
@click.option(
    "--n-samples",
    "-n",
    type=click.IntRange(min=1),
    default=1_000,
    show_default=True,
    metavar="INT",
    help="Number of (model, response) pairs to generate.",
)
@click.option(
    "--n-layers",
    default="3:7",
    show_default=True,
    metavar="INT|MIN:MAX",
    help=(
        "Layers per sample — a fixed integer (e.g. 5) or "
        "a colon-separated uniform range (e.g. 3:7)."
    ),
)
@click.option(
    "--freqs",
    "freq_range",
    type=FreqRange(),
    default="0.001:10000",
    show_default=True,
    metavar="FMIN:FMAX",
    help="Frequency range in Hz for the forward model.",
)
@click.option(
    "--n-freqs",
    type=click.IntRange(min=2),
    default=30,
    show_default=True,
    metavar="INT",
    help="Number of log-spaced frequencies between FMIN and FMAX.",
)
@click.option(
    "--noise-type",
    type=click.Choice(
        ["none", "gaussian", "field", "multiplicative"],
        case_sensitive=False,
    ),
    default="gaussian",
    show_default=True,
    help="Noise model applied to synthetic responses.",
)
@click.option(
    "--noise-level",
    type=float,
    default=0.05,
    show_default=True,
    metavar="FLOAT",
    help="Relative noise amplitude (standard deviation as a fraction).",
)
@click.option(
    "--seed",
    type=int,
    default=None,
    metavar="INT",
    help="Global random seed for fully reproducible datasets.",
)
@click.option(
    "--output",
    "-o",
    default="forward_dataset.npz",
    show_default=True,
    metavar="FILE",
    help=(
        "Output file path.  .npz saves a compressed NumPy archive; "
        ".csv saves the feature matrix only."
    ),
)
@n_jobs_option
@verbose_option
@no_color_option
@click.pass_context
def generate(
    ctx: click.Context,
    solver: str,
    n_samples: int,
    n_layers: str,
    freq_range: tuple[float, float],
    n_freqs: int,
    noise_type: str,
    noise_level: float,
    seed: int | None,
    output: str,
    n_jobs: int,
    verbose: int,
    no_color: bool,
) -> None:
    """Batch-generate synthetic (model, response) training datasets.

    Randomly samples :class:`~pycsamt.forward.synthetic.LayeredModel`
    objects and runs the chosen 1-D EM solver on each, optionally adding
    noise.  The resulting feature matrix (X) and target matrix (y) are
    saved to a compressed NumPy archive or CSV file.

    \b
    Dataset layout (saved in the .npz):
      X   — (n_samples, 2*n_freqs)  log10(ρ_a) ‖ phase for MT/CSAMT;
                                     log10(|dBz/dt|) for TEM
      y   — (n_samples, 2*max_layers-1)  log10(ρ) ‖ thickness
      meta — structured array: n_layers, noise_level per sample

    \b
    Examples:
      pycsamt forward generate --n-samples 5000 --output mt1d_train.npz

      pycsamt forward generate \\
          --solver mt1d --n-samples 50000 --n-layers 3:8 \\
          --noise-type field --noise-level 0.03 --seed 42 --jobs 4 \\
          --output mt1d_field.npz

      pycsamt forward generate \\
          --solver tem1d --n-samples 2000 \\
          --noise-type multiplicative \\
          --output tem1d_train.npz
    """
    configure_cli(log__level=verbose, log__color=not no_color)

    import numpy as np  # noqa: PLC0415

    from ....forward import generate_dataset  # noqa: PLC0415

    # ── parse --n-layers ─────────────────────────────────────────────────────
    try:
        if ":" in str(n_layers):
            lo_s, hi_s = str(n_layers).split(":", 1)
            n_layers_arg: int | tuple[int, int] = (int(lo_s), int(hi_s))
        else:
            n_layers_arg = int(n_layers)
    except (ValueError, TypeError):
        raise click.BadParameter(
            "Expected an integer (e.g. 5) or a range (e.g. 3:7).",
            param_hint="--n-layers",
        )

    freqs = np.logspace(
        np.log10(freq_range[0]),
        np.log10(freq_range[1]),
        n_freqs,
    )

    noise_type_arg = None if noise_type.lower() == "none" else noise_type.lower()
    save_path = output if output.endswith(".npz") else None

    if verbose >= 1:
        nl_str = (
            f"{n_layers_arg[0]}–{n_layers_arg[1]}"
            if isinstance(n_layers_arg, tuple)
            else str(n_layers_arg)
        )
        click.echo(
            f"Generating {n_samples} samples  "
            f"solver={solver}  n_layers={nl_str}  "
            f"freqs={freq_range[0]:.3g}–{freq_range[1]:.3g} Hz  "
            f"noise={noise_type}({noise_level:.2f})  jobs={n_jobs}",
            err=True,
        )

    try:
        ds = generate_dataset(
            solver=solver,
            n_samples=n_samples,
            freqs=freqs,
            n_layers=n_layers_arg,
            noise_level=noise_level if noise_type_arg else 0.0,
            noise_type=noise_type_arg or "gaussian",
            seed=seed,
            n_jobs=n_jobs,
            output=save_path,
        )
    except Exception as exc:  # noqa: BLE001
        click.echo(f"Error: {exc}", err=True)
        sys.exit(1)

    # ── save CSV when .npz was not requested ─────────────────────────────────
    if not output.endswith(".npz"):
        import pandas as pd  # noqa: PLC0415

        df = pd.DataFrame(
            ds.X,
            columns=[f"feature_{i}" for i in range(ds.X.shape[1])],
        )
        df.to_csv(output, index=False)

    click.echo(f"Dataset saved → {output}  (X: {ds.X.shape}, y: {ds.y.shape})")
