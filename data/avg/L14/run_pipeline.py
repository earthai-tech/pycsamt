"""Run the ``full_processing`` pipeline preset on the L14 EDIs.

Loads every ``.edi`` file under ``L14/edi/`` (produced by
``avg_to_edi.py``), runs :class:`~pycsamt.pipeline.Pipeline`'s
``full_processing`` preset (notch -> drop-duplicates -> band-select ->
grid-align -> skew-mask -> strike-rotate -> static-shift -> QC snapshot),
and writes the corrected EDIs, QC plots, and an HTML/text report to
``L14/output/``.

Usage (any cwd)::

    python data/avg/L14/run_pipeline.py
"""

from __future__ import annotations

from pathlib import Path

from pycsamt.emtools import ensure_sites
from pycsamt.pipeline import Pipeline, configure_pipe

HERE = Path(__file__).resolve().parent
EDI_DIR = HERE / "edi"
OUTPUT_DIR = HERE / "output"


def main() -> None:
    paths = sorted(EDI_DIR.glob("*.edi"))
    if not paths:
        raise SystemExit(
            f"No .edi files found in {EDI_DIR} -- run avg_to_edi.py first."
        )
    sites = ensure_sites([str(p) for p in paths])
    print(f"loaded {len(paths)} stations from {EDI_DIR}")

    pipe = Pipeline.from_preset(
        "full_processing", pipeline_name="L14_full_processing"
    )
    print(pipe)

    configure_pipe(show_progress=False)
    result = pipe.run(
        sites,
        outdir=OUTPUT_DIR,
        save_plots=True,
        save_edis=True,
        save_report=True,
    )

    print(result.summary())
    print(f"{'step':<15}{'code':<9}{'ok':<6}{'seconds':>8}{'sites':>10}")
    for sr in result.step_results:
        print(
            f"{sr.step_name:<15}{sr.step_code:<9}{str(sr.ok):<6}"
            f"{sr.elapsed_sec:>8.2f}{f'{sr.n_sites_in}->{sr.n_sites_out}':>10}"
        )


if __name__ == "__main__":
    main()
