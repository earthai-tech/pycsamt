"""Named pipeline presets for common MT processing workflows.

Each preset is a :class:`Preset` object that pairs a description with an
ordered list of :class:`~pycsamt.emtools.pipe._steps.Step` instances.  Users
can load a preset via :meth:`~pycsamt.emtools.pipe.Pipeline.from_preset` and
then add, remove, or override steps as needed.

Available presets
-----------------
``"basic_qc"``
    Minimal denoising + frequency cleanup.  Good for a first-pass inspection.
``"noise_reduction"``
    Stacked noise-removal chain for high-EMI environments.
``"full_processing"``
    Standard end-to-end workflow: noise → frequency → skew gate → SS → strike.
``"tensor_analysis"``
    Tensor-only cleanup (rotation, antisymmetry, sigma-clip, balance).
``"dimensionality_filter"``
    Classify, mask, and project to 2-D.
``"publication_ready"``
    Full chain suitable for Computers & Geosciences quality output.
``"stratagem_mt"``
    emtools-level preset for Stratagem AMT data already loaded as Sites.
    Static-shift → AMT band select → powerline notch → Hampel → QC.
    Use :mod:`pycsamt.pipeline.stratagem` when you also need coordinate
    injection and renaming from raw EDI + CSV inputs.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from ._steps import Step

# ---------------------------------------------------------------------------
# Preset dataclass
# ---------------------------------------------------------------------------

@dataclass
class Preset:
    """A named, ordered collection of pipeline steps.

    Attributes
    ----------
    name:
        Identifier used in :meth:`~pycsamt.emtools.pipe.Pipeline.from_preset`.
    description:
        One-line description shown when the user prints the preset catalogue.
    steps:
        Ordered list of ``(label, Step)`` tuples — same format as a
        :class:`~pycsamt.emtools.pipe.Pipeline` step list.
    """

    name: str
    description: str
    steps: list[tuple[str, Step]] = field(default_factory=list)

    def __repr__(self) -> str:
        codes = " → ".join(s.spec.code for _, s in self.steps)
        return f"Preset({self.name!r}  [{codes}])"


# ---------------------------------------------------------------------------
# Preset definitions
# ---------------------------------------------------------------------------

PRESETS: dict[str, Preset] = {

    "basic_qc": Preset(
        name="basic_qc",
        description=(
            "Minimal denoising + frequency cleanup.  "
            "Good for quick-look inspection."
        ),
        steps=[
            ("notch",            Step("NR001")),
            ("drop_duplicates",  Step("FREQ002")),
            ("select_band",      Step("FREQ001")),
            ("align_grid",       Step("FREQ004")),
            ("qc_snapshot",      Step("QC001")),
        ],
    ),

    "noise_reduction": Preset(
        name="noise_reduction",
        description=(
            "Stacked noise-removal chain for high-EMI environments."
        ),
        steps=[
            ("notch",         Step("NR001")),
            ("hampel",        Step("NR004")),
            ("spatial_med",   Step("NR005")),
            ("shrink_trend",  Step("NR003")),
            ("mask_incoher",  Step("NR010")),
            ("qc_snapshot",   Step("QC001")),
        ],
    ),

    "full_processing": Preset(
        name="full_processing",
        description=(
            "Standard end-to-end workflow: "
            "noise → frequency → skew gate → static-shift → strike rotation."
        ),
        steps=[
            ("notch",           Step("NR001")),
            ("drop_dup",        Step("FREQ002")),
            ("select_band",     Step("FREQ001")),
            ("align_grid",      Step("FREQ004")),
            ("mask_skew",       Step("SK001")),
            ("rotate_strike",   Step("TZ001")),
            ("correct_ss",      Step("SS001")),
            ("qc_snapshot",     Step("QC001")),
        ],
    ),

    "tensor_analysis": Preset(
        name="tensor_analysis",
        description=(
            "Tensor-only cleanup: strike rotation, antisymmetry, "
            "sigma-clip, and off-diagonal balance."
        ),
        steps=[
            ("rotate_strike",  Step("TZ001")),
            ("antisymm",       Step("TZ002")),
            ("sigma_clip",     Step("TZ003")),
            ("balance",        Step("TZ004")),
            ("qc_snapshot",    Step("QC001")),
        ],
    ),

    "dimensionality_filter": Preset(
        name="dimensionality_filter",
        description=(
            "Classify 1-D / 2-D / 3-D regions, mask unwanted class, "
            "and project to 2-D."
        ),
        steps=[
            ("classify_dim",  Step("DIM001")),
            ("mask_dim",      Step("DIM002")),
            ("project_2d",    Step("DIM003")),
            ("qc_snapshot",   Step("QC001")),
        ],
    ),

    "publication_ready": Preset(
        name="publication_ready",
        description=(
            "Full chain for publication-quality output "
            "(C&G paper standard): noise removal, frequency editing, "
            "static-shift correction, strike rotation, and skew gating."
        ),
        steps=[
            ("notch",          Step("NR001")),
            ("drop_dup",       Step("FREQ002")),
            ("select_band",    Step("FREQ001")),
            ("align_grid",     Step("FREQ004")),
            ("correct_ss",     Step("SS001")),
            ("rotate_strike",  Step("TZ001")),
            ("antisymm",       Step("TZ002")),
            ("mask_skew",      Step("SK001")),
            ("qc_snapshot",    Step("QC001")),
        ],
    ),

    # ── Stratagem-specific preset ─────────────────────────────────────────
    "stratagem_mt": Preset(
        name="stratagem_mt",
        description=(
            "emtools-level preset for Stratagem AMT data already loaded "
            "as a Sites object.  Applies AMA static-shift correction, "
            "selects the standard AMT band (10 Hz – 100 kHz), removes "
            "powerline harmonics and outliers, and masks incoherent "
            "frequency bins.  For the full raw EDI + GPS workflow "
            "(coordinate injection, hardware SNR masking, renaming) use "
            "pycsamt.pipeline.stratagem.StratagemPipeline instead."
        ),
        steps=[
            # static shift first — requires complete frequency range
            ("correct_ss",    Step("SS001")),
            # trim to AMT acquisition band
            ("select_band",   Step("FREQ001", band_hz=(10.0, 1e5))),
            # drop duplicate frequencies from multi-band acquisition
            ("drop_dup",      Step("FREQ002")),
            # powerline harmonics (50 Hz grid, 30 harmonics covers AMT range)
            ("notch",         Step("NR001", mains_hz=50.0, n_harm=30,
                                   tol_hz=0.08)),
            # Hampel outlier filter on impedance tensor
            ("hampel",        Step("NR004")),
            # mask frequencies that are incoherent across the profile
            ("mask_incoher",  Step("NR010")),
            ("qc_snapshot",   Step("QC001")),
        ],
    ),
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def list_presets() -> list[Preset]:
    """Return all available :class:`Preset` objects."""
    return list(PRESETS.values())


def get_preset(name: str) -> Preset:
    """Return the :class:`Preset` for *name*.

    Raises
    ------
    KeyError
        When *name* is not registered.
    """
    if name not in PRESETS:
        available = sorted(PRESETS)
        raise KeyError(
            f"Unknown preset {name!r}.  "
            f"Available presets: {available}"
        )
    return PRESETS[name]


def preset_catalogue() -> str:
    """Return a formatted catalogue of all presets."""
    lines = ["Available pipeline presets", "─" * 60]
    for p in PRESETS.values():
        codes = " → ".join(s.spec.code for _, s in p.steps)
        lines.append(f"  {p.name:<22s}  {p.description}")
        lines.append(f"  {'':22s}  {codes}")
        lines.append("")
    return "\n".join(lines)


__all__ = [
    "Preset",
    "PRESETS",
    "get_preset",
    "list_presets",
    "preset_catalogue",
]
