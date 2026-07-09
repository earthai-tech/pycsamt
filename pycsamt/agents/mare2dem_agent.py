"""
pycsamt.agents.mare2dem_agent
==============================

:class:`Mare2DEMAgent` — Orchestrate a complete MARE2DEM 2.5-D EM
inversion workflow.

MARE2DEM (Modeling with Adaptively Refined Elements for 2.5-D
Electromagnetic inversion) is a 2.5-D FEM MT and CSEM inversion code
(Key 2016).  Unlike Occam2D and ModEM, the MARE2DEM binary is NOT
bundled with pycsamt and must be downloaded and compiled from source
once per machine before any inversion run.

This agent handles the complete lifecycle:

1. **Source management** (optional) — detect whether the MARE2DEM
   binary already exists; if not, offer to download and compile using
   :class:`~pycsamt.models.mare2dem.SourceManager`.

2. **Input preparation** — write the ``.emdata`` data file, starting
   ``.resistivity`` model, and ``.settings`` parallel-decomposition
   file using :class:`~pycsamt.models.mare2dem.InputBuilder`.
   Accepts either an existing ``.emdata`` path *or* an inline survey
   specification via :class:`~pycsamt.models.mare2dem.MTSurveyConfig`
   / :class:`~pycsamt.models.mare2dem.CSEMSurveyConfig`.

3. **Inversion run** (optional) — launch the MARE2DEM binary through
   :class:`~pycsamt.models.mare2dem.Mare2DEMRunner` when ``mode`` is
   ``"run"`` and the binary is available.

4. **Result reporting** — scan the run directory with
   :class:`~pycsamt.models.mare2dem.InversionResult` and return the
   log, final model, and RMS convergence history.  When an LLM API key
   is provided, a brief interpretation of the inversion outcome is
   appended.

References
----------
Key, K. (2016). MARE2DEM: A 2-D inversion code for controlled-source
electromagnetic and magnetotelluric data.  *Geophysical Journal
International*, 207(1), 571–588.  doi:10.1093/gji/ggw290.
"""

from __future__ import annotations

import os
import time
from pathlib import Path
from typing import Any

from ._base import AgentResult, BaseAgent

_SYSTEM_PROMPT = """\
You are an expert in 2.5-D EM inversion using MARE2DEM.
Given a summary of the survey geometry, input files, and inversion
outcome, write 4–5 sentences that:
1. Assess the receiver and transmitter coverage relative to the target depth.
2. Comment on the starting resistivity model and Occam regularisation settings.
3. Evaluate the convergence: is the final RMS close to target? Any stagnation?
4. Recommend parallel decomposition settings (tx/rx per group) for the cluster size.
5. Flag any data issues (low-amplitude offsets, asymmetric MT responses, topo artefacts).
Reply in plain English. Be specific about MARE2DEM parameters.
"""


class Mare2DEMAgent(BaseAgent):
    """Orchestrate a complete MARE2DEM 2.5-D EM inversion workflow.

    Parameters
    ----------
    api_key : str or None
        LLM API key (Anthropic / OpenAI / Gemini).
    model : str or None
        LLM model identifier; defaults to the provider's recommended model.
    llm_provider : str
        ``"claude"`` (default), ``"openai"``,
        ``"gemini"``, or ``"deepseek"``.
    n_procs : int
        Number of MPI processes for MARE2DEM (default 4).
    use_mpi : bool
        Whether to prefix the binary with ``mpirun -np n_procs``
        (default ``True``).  Set to ``False`` only for single-process
        debug builds.
    initial_rho : float
        Starting half-space resistivity in Ω·m (default 1.0).
    target_rms : float
        Target normalised RMS misfit (default 1.0).
    max_iterations : int
        Maximum Occam inversion iterations (default 150).

    Input keys
    ----------
    ``sites`` / ``path`` : Sites or str, optional
        EDI source (directory, files, or a loaded ``Sites``).  The agent
        converts the impedances to a MARE2DEM ``.emdata`` file via
        :func:`~pycsamt.models.mare2dem.edi.make_mt_data_from_edi`
        (TE = Zxy, TM = Zyx).  Used when ``emdata``/``mt``/``csem`` are
        not supplied — this is the pathway the workflow orchestrator
        uses.
    ``error_floor`` : float, optional
        Relative error floor applied to the TE/TM apparent
        resistivities built from EDI data (default 0.05 = 5 %).
    ``output_modes`` : str, optional
        Data types written from EDI data: ``"all"`` (default), ``"TE"``,
        ``"TM"``, ``"TE+tipper"``, or ``"all impedance"``.
    ``emdata`` : str or path-like, optional
        Path to an existing ``.emdata`` data file.  When supplied the
        agent copies it to the output directory and uses it directly.
    ``resistivity`` : str or path-like, optional
        Path to an existing starting ``.resistivity`` file.  When
        omitted a homogeneous half-space model is written.
    ``mt`` : dict, optional
        Keyword arguments for
        :class:`~pycsamt.models.mare2dem.MTSurveyConfig` when building
        the data file from survey parameters rather than a pre-existing
        file.
    ``csem`` : dict, optional
        Keyword arguments for
        :class:`~pycsamt.models.mare2dem.CSEMSurveyConfig`.
    ``topo`` : float or array-like, optional
        Topography for receiver/transmitter placement. Flat scalar
        (e.g. ``-1000.0`` m for a flat seafloor) or an ``(n, 2)``
        array of ``[y, z]`` pairs.  Default ``0.0`` (flat surface).
    ``output_dir`` : str
        Directory that will receive all MARE2DEM input files and run
        output.  Created if absent (default ``"pycsamt_mare2dem"``).
    ``mode`` : {"prepare", "run", "report"}
        Workflow mode:

        * ``"prepare"`` — write input files only (default).
        * ``"run"`` — write files then launch MARE2DEM.
        * ``"report"`` — scan an existing run directory and return
          results without (re-)running.

    ``source_dir`` : str or path-like, optional
        Explicit path to the MARE2DEM source tree.  Passed to
        :class:`~pycsamt.models.mare2dem.SourceManager`; auto-resolved
        when omitted.
    ``download_source`` : bool
        When ``True`` and the binary is missing, :meth:`execute`
        automatically runs ``SourceManager.download()`` and
        ``SourceManager.build()``.  Defaults to ``False`` (raises a
        warning instead).
    ``n_procs`` : int, optional
        Per-call override for :attr:`n_procs`.
    ``max_iterations`` : int, optional
        Per-call override.
    ``target_rms`` : float, optional
        Per-call override.
    ``initial_rho`` : float, optional
        Per-call override.

    Output data keys
    ----------------
    ``data_path``          : Path or None
    ``resistivity_path``   : Path or None
    ``settings_path``      : Path or None
    ``binary_found``       : bool
    ``source_downloaded``  : bool
    ``n_mt_receivers``     : int
    ``n_csem_transmitters``: int
    ``n_data``             : int
    ``final_rms``          : float or None
    ``n_iterations``       : int
    ``converged``          : bool
    ``result``             : InversionResult or None
    ``output_dir``         : str

    Examples
    --------
    Prepare input files from an existing data file:

    >>> agent  = Mare2DEMAgent(n_procs=8)
    >>> result = agent.execute({
    ...     "emdata":     "/data/survey.emdata",
    ...     "output_dir": "/run/mare2dem",
    ... })
    >>> print(result["data_path"])
    /run/mare2dem/mare2dem.emdata

    Build from MT survey parameters (flat seafloor at −1000 m):

    >>> import numpy as np
    >>> result = agent.execute({
    ...     "mt": {
    ...         "frequencies": list(np.logspace(-3, 3, 20)),
    ...         "rx_y":        list(np.linspace(-5000, 5000, 20)),
    ...         "rx_type":     "marine",
    ...         "lTE":         True,
    ...         "lTM":         True,
    ...     },
    ...     "topo":       -1000.0,
    ...     "output_dir": "/run/mare2dem_mt",
    ... })

    Run a full inversion (requires compiled MARE2DEM binary):

    >>> result = agent.execute({
    ...     "emdata":     "survey.emdata",
    ...     "output_dir": "./run",
    ...     "mode":       "run",
    ...     "n_procs":    16,
    ... })
    >>> print(result["final_rms"])
    0.98

    Report results from a completed run directory:

    >>> result = agent.execute({
    ...     "output_dir": "./run",
    ...     "mode":       "report",
    ... })
    >>> result["converged"]
    True
    """

    SYSTEM_PROMPT = _SYSTEM_PROMPT

    def __init__(
        self,
        *,
        api_key: str | None = None,
        model: str | None = None,
        llm_provider: str = "claude",
        n_procs: int = 4,
        use_mpi: bool = True,
        initial_rho: float = 1.0,
        target_rms: float = 1.0,
        max_iterations: int = 150,
    ) -> None:
        super().__init__(
            "Mare2DEMAgent",
            api_key=api_key,
            model=model,
            llm_provider=llm_provider,
            section_preset="inversion",
        )
        self.n_procs        = n_procs
        self.use_mpi        = use_mpi
        self.initial_rho    = initial_rho
        self.target_rms     = target_rms
        self.max_iterations = max_iterations

    # ------------------------------------------------------------------
    # main entry point
    # ------------------------------------------------------------------

    def execute(self, input_data: dict[str, Any]) -> AgentResult:
        self._last_cost = 0.0
        t0 = time.time()
        warnings: list[str] = []

        # ── resolve parameters ─────────────────────────────────────────
        mode           = str(input_data.get("mode", "prepare")).lower()
        output_dir     = input_data.get("output_dir", "pycsamt_mare2dem")
        source_dir     = input_data.get("source_dir")
        auto_download  = bool(input_data.get("download_source", False))
        n_procs        = int(input_data.get("n_procs",        self.n_procs))
        max_iter       = int(input_data.get("max_iterations", self.max_iterations))
        target_rms     = float(input_data.get("target_rms",   self.target_rms))
        initial_rho    = float(input_data.get("initial_rho",  self.initial_rho))
        topo           = input_data.get("topo", 0.0)
        emdata_src     = input_data.get("emdata")
        resist_src     = input_data.get("resistivity")
        mt_kwargs      = input_data.get("mt")
        csem_kwargs    = input_data.get("csem")
        sites_src      = input_data.get("sites")
        if sites_src is None:
            sites_src = input_data.get("path")

        os.makedirs(output_dir, exist_ok=True)

        if mode not in {"prepare", "run", "report"}:
            return AgentResult.failed(
                f"Unknown mode {mode!r}. Choose 'prepare', 'run', or 'report'.",
                elapsed=time.time() - t0,
            )

        # ── import MARE2DEM model layer ────────────────────────────────
        try:
            from ..models.mare2dem import (
                Mare2DEMConfig,
                SourceManager,
                InputBuilder,
                Mare2DEMRunner,
                InversionResult,
                MTSurveyConfig,
                CSEMSurveyConfig,
                make_data_file,
                read_emdata,
            )
        except ImportError as exc:
            return AgentResult.failed(
                f"pycsamt.models.mare2dem not available: {exc}",
                elapsed=time.time() - t0,
            )

        # ── build config ───────────────────────────────────────────────
        cfg = Mare2DEMConfig(
            n_procs=n_procs,
            use_mpi=self.use_mpi,
            initial_rho=initial_rho,
            target_rms=target_rms,
            max_iterations=max_iter,
            source_dir=source_dir,
        )

        # ── source management ──────────────────────────────────────────
        source_downloaded = False
        binary_found      = False

        sm = SourceManager(config=cfg)
        binary_found = sm.is_built()
        if not binary_found:
            if auto_download:
                try:
                    self._log.info("Mare2DEMAgent: downloading MARE2DEM source...")
                    sm.download()
                    self._log.info("Mare2DEMAgent: building MARE2DEM binary...")
                    sm.build()
                    binary_found      = sm.is_built()
                    source_downloaded = binary_found
                    if binary_found:
                        self._log.info("Mare2DEMAgent: MARE2DEM binary ready.")
                    else:
                        warnings.append(
                            "Source downloaded but build may have failed. "
                            "Check compiler output."
                        )
                except Exception as exc:
                    warnings.append(
                        f"Auto-build failed: {exc}. "
                        "Install Intel oneAPI + MKL then run "
                        "SourceManager().download(); SourceManager().build()."
                    )
            else:
                warnings.append(
                    "MARE2DEM binary not found. "
                    "Set download_source=True to auto-build, or run: "
                    "SourceManager().download(); SourceManager().build()"
                )

        if mode == "report":
            return self._report(
                output_dir, cfg, binary_found, source_downloaded,
                warnings, t0,
            )

        # ── prepare input files ────────────────────────────────────────
        data_path: Path | None     = None
        resist_path: Path | None   = None
        settings_path: Path | None = None
        n_mt_rx = n_csem_tx = n_data = 0

        try:
            builder = InputBuilder(config=cfg, verbose=0)

            # ---- data file ----
            if emdata_src is not None:
                # copy / link existing file
                files = builder.build(
                    emdata_src, output_dir,
                    data_filename=cfg.data_file,
                    model_filename=cfg.resistivity_file,
                    settings_filename=cfg.settings_file,
                )
            elif mt_kwargs is not None or csem_kwargs is not None:
                # build from survey parameters
                mt_cfg = MTSurveyConfig(**mt_kwargs) if mt_kwargs else None
                csem_cfg = CSEMSurveyConfig(**csem_kwargs) if csem_kwargs else None
                files = builder.build(
                    None, output_dir,
                    topo=topo,
                    mt=mt_cfg,
                    csem=csem_cfg,
                    data_filename=cfg.data_file,
                    model_filename=cfg.resistivity_file,
                    settings_filename=cfg.settings_file,
                )
            elif sites_src is not None:
                # EDI pathway (workflow orchestrator): impedances → .emdata
                from ..models.mare2dem.edi import make_mt_data_from_edi
                err_floor = float(input_data.get("error_floor", 0.05))
                out_modes = str(input_data.get("output_modes", "all"))
                cr_value = input_data.get(
                    "confidence_weighting",
                    input_data.get("cr_weighting", False),
                )
                confidence_weighting = (
                    cr_value.strip().lower() in {"1", "true", "yes", "on"}
                    if isinstance(cr_value, str)
                    else bool(cr_value)
                )
                data_file = Path(output_dir) / cfg.data_file
                make_mt_data_from_edi(
                    sites_src, data_file,
                    output_modes=out_modes,
                    error_floor_te=err_floor,
                    error_floor_tm=err_floor,
                    confidence_weighting=confidence_weighting,
                    confidence_method=str(
                        input_data.get("confidence_method", "composite")
                    ),
                    confidence_weights=input_data.get("confidence_weights"),
                    confidence_min=float(
                        input_data.get("confidence_min", 0.05)
                    ),
                    confidence_power=float(
                        input_data.get("confidence_power", 1.0)
                    ),
                    topo=topo,
                )
                files = {
                    "data": data_file,
                    "model": builder.write_resistivity(
                        output_dir, filename=cfg.resistivity_file
                    ),
                    "settings": builder.write_settings(
                        output_dir, filename=cfg.settings_file
                    ),
                }
            else:
                return AgentResult.failed(
                    "No data source supplied. Provide 'sites'/'path' (EDI), "
                    "'emdata' (path), or 'mt'/'csem' (survey dicts).",
                    elapsed=time.time() - t0,
                )

            data_path     = files.get("data")
            resist_path   = files.get("model")
            settings_path = files.get("settings")

            # override resistivity if caller supplied one
            if resist_src is not None:
                import shutil
                dest = Path(output_dir) / cfg.resistivity_file
                shutil.copy2(str(resist_src), str(dest))
                resist_path = dest
                self._log.info("Mare2DEMAgent: used provided resistivity: %s", dest)

            # read back statistics
            if data_path and data_path.exists():
                try:
                    em = read_emdata(data_path)
                    n_mt_rx   = em.n_mt_receivers
                    n_csem_tx = em.n_csem_transmitters
                    n_data    = em.n_data
                except Exception:
                    pass

            self._log.info(
                "Mare2DEMAgent: files written to %s "
                "(mt_rx=%d, csem_tx=%d, n_data=%d)",
                output_dir, n_mt_rx, n_csem_tx, n_data,
            )

        except Exception as exc:
            return AgentResult.failed(
                f"InputBuilder failed: {exc}",
                hint="Check that the emdata file or mt/csem survey dicts are valid.",
                elapsed=time.time() - t0,
            )

        # ── run MARE2DEM ───────────────────────────────────────────────
        result_obj = None
        final_rms  = None
        n_iter     = 0
        converged  = False

        if mode == "run":
            if not binary_found:
                warnings.append(
                    "mode='run' requested but MARE2DEM binary not found. "
                    "Inversion skipped. Use download_source=True or build manually."
                )
            else:
                try:
                    runner = Mare2DEMRunner(output_dir, config=cfg, verbose=1)
                    stem = Path(cfg.resistivity_file).stem
                    self._log.info(
                        "Mare2DEMAgent: launching MARE2DEM (%d MPI processes)...",
                        n_procs,
                    )
                    result_obj = runner.run(stem, load_result=True)
                    if result_obj is not None:
                        final_rms = result_obj.final_rms
                        n_iter    = result_obj.n_iterations
                        converged = result_obj.converged
                        self._log.info(
                            "Mare2DEMAgent: inversion complete — "
                            "iterations=%d, final_rms=%.4f, converged=%s",
                            n_iter, final_rms or 0, converged,
                        )
                except Exception as exc:
                    warnings.append(
                        f"MARE2DEM run failed: {exc}. "
                        "Check the log file in the output directory."
                    )

        # ── LLM interpretation ─────────────────────────────────────────
        interp: str | None = None
        if self.api_key:
            prompt = (
                f"MARE2DEM 2.5-D EM inversion summary:\n"
                f"  Mode: {mode}, output dir: {output_dir}\n"
                f"  Survey: {n_mt_rx} MT receivers, "
                f"{n_csem_tx} CSEM transmitters, {n_data} data points\n"
                f"  Starting resistivity: {initial_rho:.4g} Ω·m\n"
                f"  Target RMS: {target_rms:.3g}, Max iterations: {max_iter}\n"
                f"  MPI processes: {n_procs}\n"
                f"  Binary found: {binary_found}\n"
                f"  Final RMS: {final_rms}, Iterations: {n_iter}, "
                f"Converged: {converged}\n"
                f"  Warnings: {warnings[:3] if warnings else 'none'}\n\n"
                "Provide specific advice on the MARE2DEM inversion setup and results."
            )
            interp = self.query_llm(prompt, max_tokens=250)

        # ── assemble result ────────────────────────────────────────────
        files_ok = sum(
            1 for p in [data_path, resist_path, settings_path]
            if p and p.exists()
        )
        elapsed = time.time() - t0

        status = "success"
        if not data_path or not data_path.exists():
            status = "needs_review"
        if mode == "run" and not converged and final_rms is None:
            status = "needs_review"

        summary_parts = [f"MARE2DEM {mode}:"]
        if n_mt_rx:
            summary_parts.append(f"{n_mt_rx} MT receivers")
        if n_csem_tx:
            summary_parts.append(f"{n_csem_tx} CSEM transmitters")
        if n_data:
            summary_parts.append(f"{n_data} data points")
        summary_parts.append(f"{files_ok}/3 files written")
        if mode == "run" and final_rms is not None:
            summary_parts.append(
                f"final RMS={final_rms:.3f} (converged={converged})"
            )
        summary = ". ".join(summary_parts) + f". Output: {output_dir}."

        return AgentResult(
            status=status,
            summary=summary,
            data={
                "data_path":           data_path,
                "resistivity_path":    resist_path,
                "settings_path":       settings_path,
                "binary_found":        binary_found,
                "source_downloaded":   source_downloaded,
                "n_mt_receivers":      n_mt_rx,
                "n_csem_transmitters": n_csem_tx,
                "n_data":              n_data,
                "final_rms":           final_rms,
                "n_iterations":        n_iter,
                "converged":           converged,
                "result":              result_obj,
                "output_dir":          output_dir,
            },
            warnings=warnings,
            llm_interpretation=interp,
            elapsed_seconds=elapsed,
            cost_estimate_usd=self._last_cost,
        )

    # ------------------------------------------------------------------
    # report mode helper
    # ------------------------------------------------------------------

    def _report(
        self,
        output_dir: str,
        cfg: Any,
        binary_found: bool,
        source_downloaded: bool,
        warnings: list[str],
        t0: float,
    ) -> AgentResult:
        """Scan *output_dir* for MARE2DEM results and return a report."""
        from ..models.mare2dem import InversionResult

        try:
            result_obj = InversionResult(output_dir, config=cfg)
        except Exception as exc:
            return AgentResult.failed(
                f"Could not scan run directory: {exc}",
                elapsed=time.time() - t0,
            )

        final_rms = result_obj.final_rms
        n_iter    = result_obj.n_iterations
        converged = result_obj.converged

        n_mt_rx   = result_obj.data.n_data if result_obj.data else 0
        n_csem_tx = (
            result_obj.data.n_csem_transmitters
            if result_obj.data else 0
        )
        n_data    = result_obj.data.n_data if result_obj.data else 0

        self._log.info(
            "Mare2DEMAgent (report): iterations=%d, final_rms=%s, converged=%s",
            n_iter, final_rms, converged,
        )

        # LLM interpretation
        interp: str | None = None
        if self.api_key:
            prompt = (
                f"MARE2DEM inversion result report:\n"
                f"  Run directory: {output_dir}\n"
                f"  Completed iterations: {n_iter}\n"
                f"  Final normalised RMS: {final_rms}\n"
                f"  Converged: {converged}\n"
                f"  Model loaded: {result_obj.model is not None}\n"
                f"  Response loaded: {result_obj.response is not None}\n"
                f"  Warnings: {warnings[:3] if warnings else 'none'}\n\n"
                "Interpret the inversion quality and recommend next steps."
            )
            interp = self.query_llm(prompt, max_tokens=200)

        elapsed = time.time() - t0
        status = "success" if result_obj.log else "needs_review"
        summary = (
            f"MARE2DEM report: {n_iter} iterations, "
            f"final RMS={final_rms}, converged={converged}. "
            f"Dir: {output_dir}."
        )

        return AgentResult(
            status=status,
            summary=summary,
            data={
                "data_path":           None,
                "resistivity_path":    None,
                "settings_path":       None,
                "binary_found":        binary_found,
                "source_downloaded":   source_downloaded,
                "n_mt_receivers":      n_mt_rx,
                "n_csem_transmitters": n_csem_tx,
                "n_data":              n_data,
                "final_rms":           final_rms,
                "n_iterations":        n_iter,
                "converged":           converged,
                "result":              result_obj,
                "output_dir":          output_dir,
            },
            warnings=warnings,
            llm_interpretation=interp,
            elapsed_seconds=elapsed,
            cost_estimate_usd=self._last_cost,
        )


__all__ = ["Mare2DEMAgent"]
