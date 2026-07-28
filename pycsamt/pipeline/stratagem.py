# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pipeline.stratagem
==================

Stratagem AMT pipeline — extends the standard :class:`~._pipeline.Pipeline`
with hardware-survey pre-processing and post-processing hooks.

Two levels of API
-----------------

**Level 1 — :class:`StratagemPipeline`** (``Pipeline`` subclass):

    Follows exactly the same ``run(sites, *, outdir=...)`` interface as
    :class:`~._pipeline.Pipeline`.  The *sites* argument accepts every input
    type that :func:`~pycsamt.emtools._core.ensure_sites` handles:

    * :class:`~pycsamt.site.base.Sites` object (used as-is)
    * ``str`` or :class:`pathlib.Path` pointing to a **single EDI file**
    * ``str`` or :class:`pathlib.Path` pointing to an **EDI directory**
    * A single :class:`~pycsamt.seg.edi.EDIFile` object
    * A list of :class:`~pycsamt.seg.edi.EDIFile` objects

    Optional pre-processing before the emtools steps:

    * **Coordinate injection** — when *coord_file* is set, GPS coordinates
      from a CSV / XLS / XLSX file are injected into the EDI ``>HEAD``
      sections.
    * **Hardware SNR masking** — when *raw_dir* is set, zero-stack frequency
      rows from the Stratagem hardware files are masked.

    Optional post-processing after the emtools steps:

    * **Rename** — when *rename_basename* is set, the ``processed/`` EDI
      files are renamed with the given basename (``T2.000.edi``, …).

**Level 2 — :func:`run_stratagem_preset` / :class:`StratagemPreset`**:

    Convenience wrappers around :class:`~pycsamt.stratagem.survey.StratagemSurvey`
    that use the full Stratagem workflow (coord injection + SS + noise + rename).
    Useful for quick one-liner scripts that replace the legacy
    ``stratagem_edi_process_script.py``.

Usage examples
--------------
One-liner (Level 2, backwards-compatible)::

    from pycsamt.pipeline.stratagem import run_stratagem_preset

    sv = run_stratagem_preset(
        "basic",
        edi_dir="2/2EDI",
        coord_file="2.csv",
        outdir="2/processed",
        epsg=32649,
        rename_basename="T2.",
    )

Pipeline API with a directory of EDIs (Level 1)::

    from pycsamt.pipeline.stratagem import StratagemPipeline

    pipe = StratagemPipeline.from_preset(
        "stratagem_mt",
        coord_file="2.csv",
        epsg=32649,
        rename_basename="T2.",
    )
    result = pipe.run("2/2EDI", outdir="2/processed")

Pipeline API with a single EDI file::

    result = pipe.run("2/2EDI/Z2HX002.edi", outdir="tmp/")

Pipeline API with an already-loaded Sites::

    from pycsamt.emtools._core import ensure_sites

    sites = ensure_sites("2/2EDI")
    result = pipe.run(sites, outdir="2/processed")

Custom step list::

    from pycsamt.pipeline import Step

    pipe = StratagemPipeline(
        [
            ("correct_ss", Step("SS001")),
            ("select_band", Step("FREQ001", band_hz=(10.0, 1e5))),
            ("notch", Step("NR001")),
            ("hampel", Step("NR004")),
            ("qc", Step("QC001")),
        ],
        coord_file="2.csv",
        raw_dir="原始数据/2HX",
        epsg=32649,
        rename_basename="T2.",
    )
    result = pipe.run("2/2EDI", outdir="2/processed")
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from ._pipeline import Pipeline, PipelineResult
from ._steps import Step

__all__ = [
    "StratagemPreset",
    "STRATAGEM_PRESETS",
    "StratagemPipeline",
    "get_stratagem_preset",
    "list_stratagem_presets",
    "run_stratagem_preset",
    "stratagem_preset_catalogue",
]


# ---------------------------------------------------------------------------
# StratagemPreset  (Level 2 — StratagemSurvey-based presets)
# ---------------------------------------------------------------------------


@dataclass
class StratagemPreset:
    """Named configuration bundle for :class:`~pycsamt.stratagem.survey.StratagemSurvey`.

    Used by :func:`run_stratagem_preset`.  Each preset describes an ordered
    sequence of :class:`~pycsamt.stratagem.survey.StratagemSurvey` method
    calls with default keyword arguments.

    Attributes
    ----------
    name : str
    description : str
    survey_defaults : dict
        Default keyword arguments for
        :class:`~pycsamt.stratagem.survey.StratagemSurvey.__init__`.
    steps : list of (str, dict)
        Ordered ``(method_name, kwargs)`` pairs executed on the survey.
    """

    name: str
    description: str
    survey_defaults: dict = field(default_factory=dict)
    steps: list[tuple[str, dict]] = field(default_factory=list)

    def __repr__(self) -> str:
        step_names = " → ".join(n for n, _ in self.steps)
        return f"StratagemPreset({self.name!r}  [{step_names}])"


STRATAGEM_PRESETS: dict[str, StratagemPreset] = {
    "basic": StratagemPreset(
        name="basic",
        description=(
            "Direct equivalent of the legacy stratagem_edi_process_script.py: "
            "inject GPS coords → AMA static-shift → frequency trim → "
            "noise removal → export → rename."
        ),
        survey_defaults={"epsg": 32649, "utm_zone": "49N"},
        steps=[
            (
                "remove_static_shift",
                {"half_window": 3, "weights": "tri"},
            ),
            (
                "drop_frequencies",
                {
                    "fmin": 10.0,
                    "snr_thresh": 2.5,
                    "min_frac": 0.4,
                    "use_hardware_mask": True,
                },
            ),
            (
                "remove_noises",
                {
                    "mains_hz": 50.0,
                    "n_harm": 30,
                    "hampel_win": 3,
                    "smooth": False,
                },
            ),
            ("export", {}),
            ("rename", {}),
        ],
    ),
    "full_processing": StratagemPreset(
        name="full_processing",
        description=(
            "QC report → AMA static-shift → hardware-aware frequency filter "
            "→ noise removal with smoothing → export → rename."
        ),
        survey_defaults={"epsg": 32649, "utm_zone": "49N"},
        steps=[
            (
                "run_qc",
                {
                    "include_skew": True,
                    "min_frac_ok": 0.6,
                    "min_snr_med": 2.0,
                    "max_skew_med": 6.0,
                },
            ),
            (
                "remove_static_shift",
                {
                    "half_window": 3,
                    "weights": "tri",
                    "max_skew": 6.0,
                },
            ),
            (
                "drop_frequencies",
                {
                    "fmin": 10.0,
                    "fmax": 1e5,
                    "snr_thresh": 2.5,
                    "min_frac": 0.4,
                    "use_hardware_mask": True,
                },
            ),
            (
                "remove_noises",
                {
                    "mains_hz": 50.0,
                    "n_harm": 30,
                    "hampel_win": 3,
                    "smooth": True,
                    "smooth_win": 3,
                },
            ),
            ("export", {}),
            ("rename", {}),
        ],
    ),
    "publication_ready": StratagemPreset(
        name="publication_ready",
        description=(
            "C&G standard: full QC, hardware SNR masking, tight AMT band "
            "(15 Hz – 100 kHz), aggressive smoothing → export → rename."
        ),
        survey_defaults={"epsg": 32649, "utm_zone": "49N"},
        steps=[
            (
                "run_qc",
                {
                    "include_skew": True,
                    "min_frac_ok": 0.5,
                    "min_snr_med": 3.0,
                    "max_skew_med": 5.0,
                },
            ),
            (
                "remove_static_shift",
                {
                    "half_window": 5,
                    "weights": "tri",
                    "max_skew": 5.0,
                },
            ),
            (
                "drop_frequencies",
                {
                    "fmin": 15.0,
                    "fmax": 1e5,
                    "snr_thresh": 3.0,
                    "min_frac": 0.5,
                    "use_hardware_mask": True,
                },
            ),
            (
                "remove_noises",
                {
                    "mains_hz": 50.0,
                    "n_harm": 30,
                    "tol_hz": 0.05,
                    "hampel_win": 3,
                    "hampel_nsig": 3.0,
                    "smooth": True,
                    "smooth_win": 3,
                },
            ),
            ("export", {}),
            ("rename", {}),
        ],
    ),
}


def get_stratagem_preset(name: str) -> StratagemPreset:
    """Return the :class:`StratagemPreset` for *name*."""
    if name not in STRATAGEM_PRESETS:
        raise KeyError(
            f"Unknown Stratagem preset {name!r}.  "
            f"Available: {sorted(STRATAGEM_PRESETS)}"
        )
    return STRATAGEM_PRESETS[name]


def list_stratagem_presets() -> list[StratagemPreset]:
    """Return all :class:`StratagemPreset` objects."""
    return list(STRATAGEM_PRESETS.values())


def stratagem_preset_catalogue() -> str:
    """Return a formatted catalogue of all Stratagem presets."""
    lines = ["Stratagem pipeline presets", "─" * 64]
    for p in STRATAGEM_PRESETS.values():
        steps_str = " → ".join(n for n, _ in p.steps)
        lines.append(f"  {p.name:<22s}  {p.description}")
        lines.append(f"  {'':22s}  {steps_str}")
        lines.append("")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def _coerce_to_sites(sites: Any) -> Any:
    """Normalise any Stratagem-compatible input to a :class:`~pycsamt.site.base.Sites`.

    Accepted inputs:

    * :class:`~pycsamt.site.base.Sites`                 — returned unchanged
    * ``str`` / :class:`pathlib.Path` (dir)             — loads all EDIs
    * ``str`` / :class:`pathlib.Path` (single .edi)     — loads one EDI
    * :class:`~pycsamt.seg.edi.EDIFile`                 — wrapped in Sites
    * list of :class:`~pycsamt.seg.edi.EDIFile`         — wrapped in Sites
    """
    from ..emtools._core import ensure_sites  # noqa: PLC0415
    from ..seg.edi import EDIFile  # noqa: PLC0415

    # Wrap a bare EDIFile so ensure_sites doesn't confuse it with a path-like
    if isinstance(sites, EDIFile):
        sites = [sites]

    return ensure_sites(sites)


def _inject_coordinates(
    sites: Any,
    coord_file: str | Path,
    epsg: int,
    utm_zone: str,
    order: str,
) -> Any:
    """Apply coordinate injection and return the updated Sites.

    When the number of EDIs does not match the number of rows in
    *coord_file* (e.g. the caller passed a single EDI but the CSV has
    83 rows), injection is skipped with a warning and the original
    Sites is returned unchanged.
    """
    import warnings  # noqa: PLC0415

    from ..emtools._core import (  # noqa: PLC0415
        _iter_items,
        ensure_sites,
    )
    from ..exceptions import ValidationError  # noqa: PLC0415
    from ..stratagem.gis_correct import (
        CoordinateInjector,  # noqa: PLC0415
    )

    edis = [getattr(ed, "edi", ed) for ed in _iter_items(sites)]
    try:
        inj = CoordinateInjector(
            epsg=epsg,
            utm_zone=utm_zone,
            order=order,
        ).fit(edis, coord_file)
        return ensure_sites(inj.edi_objects_)
    except ValidationError as exc:
        warnings.warn(
            f"[StratagemPipeline] coordinate injection skipped: {exc}  "
            "Proceeding without injecting coordinates.",
            UserWarning,
            stacklevel=3,
        )
        return sites


def _apply_hardware_mask(sites: Any, raw_dir: str | Path) -> Any:
    """Apply hardware SNR mask and return the updated Sites."""
    from ..emtools._core import (  # noqa: PLC0415
        _iter_items,
        ensure_sites,
    )
    from ..stratagem.io import (
        StratagemRawReader,  # noqa: PLC0415
    )
    from ..stratagem.qc import (
        FrequencyFilter,  # noqa: PLC0415
    )

    rdr = StratagemRawReader(raw_dir).fit()
    edis = [getattr(ed, "edi", ed) for ed in _iter_items(sites)]
    ff = FrequencyFilter(use_hardware_mask=True, snr_thresh=0.0, min_frac=0.0)
    ff.fit(edis, raw_reader=rdr)
    return ensure_sites(ff.edi_objects_)


def _rename_processed(
    processed_dir: Path,
    rename_dir: Path,
    basename: str,
    overwrite: bool = False,
) -> list[Path]:
    """Rename EDI files in *processed_dir* to *rename_dir* with *basename*."""
    from ..stratagem.rename import EDIRenamer  # noqa: PLC0415

    if not processed_dir.is_dir():
        return []
    rename_dir.mkdir(parents=True, exist_ok=True)
    rn = EDIRenamer(basename=basename, overwrite=overwrite)
    rn.fit(processed_dir, rename_dir)
    return rn.dst_paths()


# ---------------------------------------------------------------------------
# StratagemPipeline  (Level 1 — proper Pipeline subclass)
# ---------------------------------------------------------------------------


class StratagemPipeline(Pipeline):
    """A :class:`~._pipeline.Pipeline` extended with Stratagem pre/post-processing.

    Follows exactly the same ``run(sites, *, outdir=...)`` interface as the
    standard pipeline.  The *sites* argument is normalised via
    :func:`~pycsamt.emtools._core.ensure_sites`, so it accepts:

    * :class:`~pycsamt.site.base.Sites`
    * ``str`` / :class:`pathlib.Path` pointing to a **directory** of EDIs
    * ``str`` / :class:`pathlib.Path` pointing to a **single EDI** file
    * A single :class:`~pycsamt.seg.edi.EDIFile` object
    * A list of :class:`~pycsamt.seg.edi.EDIFile` objects

    Optional pre-processing (applied before emtools steps)
    -------------------------------------------------------
    * **Coordinate injection** — when *coord_file* is given, GPS coordinates
      from a CSV / XLS / XLSX table are written into each EDI ``>HEAD``
      section via :class:`~pycsamt.stratagem.gis_correct.CoordinateInjector`.
    * **Hardware SNR mask** — when *raw_dir* is given, zero-stack frequency
      rows in the Stratagem raw files are masked in the impedance tensor via
      :class:`~pycsamt.stratagem.qc.FrequencyFilter`.

    Optional post-processing (applied after emtools steps)
    -------------------------------------------------------
    * **Rename** — when *rename_basename* is given, the output EDI files in
      ``<outdir>/processed/`` are copied to ``<outdir>/renamed/`` with
      a standardised naming convention using
      :class:`~pycsamt.stratagem.rename.EDIRenamer`.

    Parameters
    ----------
    steps : list of (str, Step) or list of Step
        emtools processing steps (same format as :class:`~._pipeline.Pipeline`).
    coord_file : path-like, optional
        GPS coordinate table (CSV / XLS / XLSX).
    raw_dir : path-like, optional
        Directory of raw Stratagem hardware files for hardware SNR masking.
    epsg : int, default 32649
        EPSG code of the projected CRS in *coord_file*.
    utm_zone : str, default ``'49N'``
    order : str, default ``'auto'``
        Station ordering for
        :class:`~pycsamt.stratagem.gis_correct.StationLocator`.
    rename_basename : str, optional
        When given, output EDIs are renamed ``{basename}000.edi``, …
    rename_dir : path-like, optional
        Destination for renamed files.  Defaults to ``<outdir>/renamed/``.
    name : str, default ``'stratagem_pipeline'``

    Examples
    --------
    Load a directory and inject coordinates::

        from pycsamt.pipeline.stratagem import StratagemPipeline

        pipe = StratagemPipeline.from_preset(
            "stratagem_mt",
            coord_file="2.csv",
            epsg=32649,
            rename_basename="T2.",
        )
        result = pipe.run("2/2EDI", outdir="2/processed")
        print(result.summary())

    Single EDI file::

        result = pipe.run("2/2EDI/Z2HX002.edi", outdir="tmp/single/")

    Already-loaded Sites::

        from pycsamt.emtools._core import ensure_sites

        S = ensure_sites("2/2EDI")
        result = pipe.run(S, outdir="2/processed")

    Custom steps with hardware mask::

        from pycsamt.pipeline import Step

        pipe = StratagemPipeline(
            [
                ("ss", Step("SS001")),
                ("band", Step("FREQ001", band_hz=(10, 1e5))),
            ],
            coord_file="2.csv",
            raw_dir="原始数据/2HX",
            epsg=32649,
        )
        result = pipe.run("2/2EDI", outdir="out/")
    """

    def __init__(
        self,
        steps: list[tuple[str, Step]] | list[Step],
        *,
        coord_file: str | Path | None = None,
        raw_dir: str | Path | None = None,
        epsg: int = 32649,
        utm_zone: str = "49N",
        order: str = "auto",
        rename_basename: str | None = None,
        rename_dir: str | Path | None = None,
        name: str = "stratagem_pipeline",
    ) -> None:
        super().__init__(steps, name=name)
        self.coord_file = coord_file
        self.raw_dir = raw_dir
        self.epsg = epsg
        self.utm_zone = utm_zone
        self.order = order
        self.rename_basename = rename_basename
        self.rename_dir = rename_dir

    # ------------------------------------------------------------------
    def run(
        self,
        sites: Any,
        *,
        outdir: Any = None,
        save_plots: bool = True,
        save_edis: bool = True,
        save_report: bool = True,
        api: Any = None,
        rename_basename: str | None = None,
        rename_dir: str | Path | None = None,
        overwrite: bool = False,
    ) -> PipelineResult:
        """Run the pipeline on *sites*.

        Parameters
        ----------
        sites : Sites | str | Path | EDIFile | list[EDIFile]
            Input data.  Any form accepted by
            :func:`~pycsamt.emtools._core.ensure_sites`.
        outdir : path-like, optional
            Root output directory (same as :meth:`Pipeline.run`).
        save_plots : bool, default True
        save_edis : bool, default True
        save_report : bool, default True
        api : PipelineAPIConfig, optional
        rename_basename : str, optional
            Override the *rename_basename* set in ``__init__``.
        rename_dir : path-like, optional
            Override the *rename_dir* set in ``__init__``.
        overwrite : bool, default False
            Overwrite existing renamed files.

        Returns
        -------
        :class:`~._pipeline.PipelineResult`
        """
        # ── 1. normalise input ────────────────────────────────────────
        S = _coerce_to_sites(sites)

        # ── 2. optional coord injection ───────────────────────────────
        if self.coord_file is not None:
            S = _inject_coordinates(
                S,
                self.coord_file,
                self.epsg,
                self.utm_zone,
                self.order,
            )

        # ── 3. optional hardware SNR mask ─────────────────────────────
        if self.raw_dir is not None:
            S = _apply_hardware_mask(S, self.raw_dir)

        # ── 4. standard emtools pipeline ──────────────────────────────
        result = super().run(
            S,
            outdir=outdir,
            save_plots=save_plots,
            save_edis=save_edis,
            save_report=save_report,
            api=api,
        )

        # ── 5. optional rename ────────────────────────────────────────
        _basename = rename_basename or self.rename_basename
        if _basename and result.outdir is not None:
            from ..api.pipe import (
                PYCSAMT_PIPE,  # noqa: PLC0415
            )

            cfg = api or PYCSAMT_PIPE
            processed_sub = getattr(cfg, "processed_subdir", "processed")
            processed_dir = result.outdir / processed_sub

            _rdir = rename_dir or self.rename_dir
            dst_dir = (
                Path(_rdir).expanduser().resolve()
                if _rdir is not None
                else result.outdir / "renamed"
            )
            renamed = _rename_processed(
                processed_dir, dst_dir, _basename, overwrite=overwrite
            )
            result.processed_paths.extend(renamed)

        return result

    # ------------------------------------------------------------------
    @classmethod
    def from_preset(
        cls,
        name: str = "stratagem_mt",
        *,
        coord_file: str | Path | None = None,
        raw_dir: str | Path | None = None,
        epsg: int = 32649,
        utm_zone: str = "49N",
        order: str = "auto",
        rename_basename: str | None = None,
        rename_dir: str | Path | None = None,
        pipeline_name: str | None = None,
    ) -> StratagemPipeline:
        """Build a :class:`StratagemPipeline` from a named emtools preset.

        Parameters
        ----------
        name : str, default ``'stratagem_mt'``
            Any preset in :data:`~._presets.PRESETS`
            (e.g. ``'stratagem_mt'``, ``'full_processing'``,
            ``'publication_ready'``).
        coord_file, raw_dir, epsg, utm_zone, order
            Forwarded to the constructor.
        rename_basename, rename_dir
            Forwarded to the constructor.
        pipeline_name : str, optional
            Override the pipeline label.

        Returns
        -------
        StratagemPipeline

        Examples
        --------
        >>> pipe = StratagemPipeline.from_preset(
        ...     "stratagem_mt",
        ...     coord_file="2.csv",
        ...     epsg=32649,
        ...     rename_basename="T2.",
        ... )
        >>> result = pipe.run("2/2EDI", outdir="2/processed")
        """
        from ._presets import get_preset  # noqa: PLC0415

        preset = get_preset(name)
        return cls(
            list(preset.steps),
            coord_file=coord_file,
            raw_dir=raw_dir,
            epsg=epsg,
            utm_zone=utm_zone,
            order=order,
            rename_basename=rename_basename,
            rename_dir=rename_dir,
            name=pipeline_name or preset.name,
        )


# ---------------------------------------------------------------------------
# run_stratagem_preset  (Level 2 — StratagemSurvey-based convenience fn)
# ---------------------------------------------------------------------------


def run_stratagem_preset(
    preset: str,
    edi_dir: str | Path,
    coord_file: str | Path,
    outdir: str | Path,
    *,
    raw_dir: str | Path | None = None,
    epsg: int = 32649,
    utm_zone: str = "49N",
    rename_basename: str = "S",
    rename_dir: str | Path | None = None,
    step_overrides: dict[str, dict] | None = None,
    overwrite: bool = False,
    verbose: int = 0,
) -> Any:
    """Execute a named Stratagem preset in one call (convenience wrapper).

    Uses :class:`~pycsamt.stratagem.survey.StratagemSurvey` internally for
    the full EDI-directory + GPS-CSV workflow.  For the pipeline-style API
    (``run(sites, outdir=...)``) use :class:`StratagemPipeline` directly.

    Parameters
    ----------
    preset : {'basic', 'full_processing', 'publication_ready'}
    edi_dir : path-like
    coord_file : path-like
    outdir : path-like
    raw_dir : path-like, optional
    epsg : int, default 32649
    utm_zone : str, default ``'49N'``
    rename_basename : str, default ``'S'``
    rename_dir : path-like, optional
    step_overrides : dict, optional
        Per-step parameter overrides (see :class:`StratagemPipeline`).
    overwrite : bool, default False
    verbose : int, default 0

    Returns
    -------
    :class:`~pycsamt.stratagem.survey.StratagemSurvey`
        The completed survey object.

    Examples
    --------
    Replicate the legacy script::

        sv = run_stratagem_preset(
            "basic",
            edi_dir="2/2EDI",
            coord_file="2.csv",
            outdir="2/processed",
            epsg=32649,
            rename_basename="T2.",
        )

    With hardware files::

        sv = run_stratagem_preset(
            "full_processing",
            edi_dir="2/2EDI",
            coord_file="2.csv",
            outdir="2/processed",
            raw_dir="原始数据/2HX",
            epsg=32649,
            rename_basename="T2.",
        )
    """
    from ..stratagem.survey import (
        StratagemSurvey,  # noqa: PLC0415
    )

    preset_cfg = get_stratagem_preset(preset)
    defaults = preset_cfg.survey_defaults.copy()
    defaults.update(epsg=epsg, utm_zone=utm_zone)

    sv = StratagemSurvey(
        edi_dir=edi_dir,
        coord_file=coord_file,
        raw_dir=raw_dir,
        verbose=verbose,
        epsg=defaults.get("epsg", epsg),
        utm_zone=defaults.get("utm_zone", utm_zone),
    ).fit()

    out_base = Path(outdir).expanduser().resolve()
    export_dir = out_base / "corrected"
    ren_dir = (
        Path(rename_dir).expanduser().resolve()
        if rename_dir is not None
        else out_base / "renamed"
    )

    for method_name, preset_kw in preset_cfg.steps:
        kw = {**preset_kw, **(step_overrides or {}).get(method_name, {})}
        if method_name == "export":
            sv.export(export_dir, overwrite=overwrite)
        elif method_name == "rename":
            sv.rename(
                basename=rename_basename,
                dst_path=ren_dir,
                overwrite=overwrite,
            )
        else:
            method = getattr(sv, method_name, None)
            if method is not None:
                method(**kw)
            elif verbose:
                print(f"[run_stratagem_preset] unknown step {method_name!r}, skip")

    return sv
