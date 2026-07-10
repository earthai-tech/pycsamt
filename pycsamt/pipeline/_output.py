"""Output-directory management for the pyCSAMT pipeline.

:class:`OutputDir` creates and manages the canonical directory tree produced
by a pipeline run::

    <root>/
    ├── processed/          ← final EDI files (config: processed_subdir)
    ├── plots/              ← QC figures, one sub-folder per step
    │   ├── 01_notch/
    │   │   ├── harmonic_waterfall.png
    │   │   └── snr_gain_profile.png
    │   └── 02_select_band/
    │       └── band_microstrips.png
    ├── pipeline.yaml       ← reproduced config for reproducibility
    ├── report.html
    └── summary.txt

The directory is created lazily when :meth:`OutputDir.setup` is called.
Nothing is written until that point so that a pipeline that only runs in
memory (``outdir=None``) produces no filesystem side-effects.
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    import matplotlib.figure


class OutputDir:
    """Manages the output directory tree for a single pipeline run.

    Parameters
    ----------
    root:
        Root path.  If ``None``, falls back to
        :attr:`~pycsamt.api.pipe.PipelineAPIConfig.output_root`.
    api:
        A :class:`~pycsamt.api.pipe.PipelineAPIConfig` instance.  When
        ``None``, the global :data:`~pycsamt.api.pipe.PYCSAMT_PIPE` singleton
        is used.
    """

    def __init__(
        self,
        root: str | Path | None,
        *,
        api: Any = None,
    ) -> None:
        self._api = api
        cfg = self._cfg()
        self.root = Path(root) if root is not None else Path(cfg.output_root)
        self.processed_dir = self.root / cfg.processed_subdir
        self.plots_dir = self.root / cfg.plots_subdir
        self._ready = False

    # ------------------------------------------------------------------
    # Setup
    # ------------------------------------------------------------------

    def setup(self) -> None:
        """Create the directory tree.  Safe to call multiple times."""
        self.root.mkdir(parents=True, exist_ok=True)
        self.processed_dir.mkdir(parents=True, exist_ok=True)
        self.plots_dir.mkdir(parents=True, exist_ok=True)
        self._ready = True

    # ------------------------------------------------------------------
    # Per-step plot directory
    # ------------------------------------------------------------------

    def step_plot_dir(self, idx: int, step_name: str) -> Path:
        """Return (and create) the plot sub-directory for step *idx*."""
        folder = self.plots_dir / f"{idx:02d}_{step_name}"
        folder.mkdir(parents=True, exist_ok=True)
        return folder

    # ------------------------------------------------------------------
    # Saving figures
    # ------------------------------------------------------------------

    def save_figure(
        self,
        fig: matplotlib.figure.Figure,
        fn_name: str,
        step_idx: int,
        step_name: str,
        *,
        api: Any = None,
    ) -> Path | None:
        """Save *fig* into the step's plot sub-directory.

        Parameters
        ----------
        fig:
            Matplotlib figure to save.
        fn_name:
            Name used as the filename stem (derived from the QC function name).
        step_idx:
            1-based step position.
        step_name:
            Step label (used for the sub-folder name).
        api:
            Config override.

        Returns
        -------
        Path or None
            The path where the figure was saved, or ``None`` on failure.
        """
        cfg = api or self._cfg()
        plot_dir = self.step_plot_dir(step_idx, step_name)
        path = plot_dir / f"{fn_name}.{cfg.plot_fmt}"
        try:
            fig.savefig(str(path), dpi=cfg.plot_dpi, bbox_inches="tight")
            return path
        except Exception as exc:
            warnings.warn(
                f"Failed to save figure {fn_name!r}: {exc}",
                stacklevel=2,
            )
            return None

    # ------------------------------------------------------------------
    # Writing processed EDIs
    # ------------------------------------------------------------------

    def write_edis(self, sites: Any) -> list[Path]:
        """Write processed EDI files to the ``processed/`` sub-directory.

        Uses :func:`pycsamt.site.export.write_sites`.

        Returns
        -------
        list[Path]
            Paths of written EDI files.
        """
        try:
            from pycsamt.site.export import write_sites
            return write_sites(
                sites,
                self.processed_dir,
                exist_ok=True,
            )
        except Exception as exc:
            warnings.warn(
                f"EDI export failed: {exc}.  "
                "Processed EDI files were not written.",
                stacklevel=2,
            )
            return []

    # ------------------------------------------------------------------
    # Config snapshot
    # ------------------------------------------------------------------

    def save_pipeline_config(self, yaml_str: str) -> Path | None:
        """Write *yaml_str* as ``pipeline.yaml`` in the root directory."""
        path = self.root / "pipeline.yaml"
        try:
            path.write_text(yaml_str, encoding="utf-8")
            return path
        except Exception as exc:
            warnings.warn(f"Could not save pipeline.yaml: {exc}", stacklevel=2)
            return None

    def save_text(self, text: str, filename: str) -> Path | None:
        """Write *text* to *filename* in the root directory."""
        path = self.root / filename
        try:
            path.write_text(text, encoding="utf-8")
            return path
        except Exception as exc:
            warnings.warn(f"Could not save {filename}: {exc}", stacklevel=2)
            return None

    # ------------------------------------------------------------------
    # Internal
    # ------------------------------------------------------------------

    def _cfg(self) -> Any:
        if self._api is not None:
            return self._api
        from pycsamt.api.pipe import PYCSAMT_PIPE
        return PYCSAMT_PIPE

    def __repr__(self) -> str:
        status = "ready" if self._ready else "pending"
        return f"OutputDir({self.root}  [{status}])"


__all__ = ["OutputDir"]
