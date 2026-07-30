"""
python -m pycsamt.agents — CLI entry point.

Usage
-----
Preview a workflow without running it::

    python -m pycsamt.agents preview "Load EDIs from /data/L22PLT, \\
        QC them, period range 1e-4 to 1 s"

Show the agent catalogue::

    python -m pycsamt.agents list

Show pricing tables::

    python -m pycsamt.agents pricing

Browse the pre-trained model zoo::

    python -m pycsamt.agents zoo

Download a specific pre-trained checkpoint::

    python -m pycsamt.agents zoo mt1d-resnet-5layer-v1

Launch the Gradio web interface::

    python -m pycsamt.agents web [--share] [--port=7860]
"""

from __future__ import annotations

import builtins
import sys


def print(*args, **kwargs):  # noqa: A001 — module-local shadow
    """Console-safe print: legacy Windows consoles use cp1252, which
    cannot encode the box-drawing/emoji characters in this CLI's
    output. Degrade lossily instead of crashing."""
    try:
        builtins.print(*args, **kwargs)
    except UnicodeEncodeError:
        enc = getattr(sys.stdout, "encoding", None) or "ascii"
        cleaned = [
            str(a).encode(enc, errors="replace").decode(enc) for a in args
        ]
        builtins.print(*cleaned, **kwargs)


def _cmd_preview(request: str) -> None:
    from . import (
        AgentCoordinator,
        ContextInputAgent,
        MTLoaderAgent,
    )

    ctx = ContextInputAgent()
    loader = MTLoaderAgent()

    coord = AgentCoordinator("preview_workflow", verbose=True)
    coord.add_step(
        "parse",
        ctx,
        description="Parse natural-language request into workflow config",
    )
    coord.add_step(
        "load",
        loader,
        input_fn=lambda r: {
            "path": (r["parse"].get("config") or {}).get("data_path", "")
        },
        description="Load MT data and run per-station QC scan",
    )

    # show config extracted from request
    ctx_result = ctx.execute({"request": request})
    print("\n── Extracted config ──────────────────────────────────────")
    import json

    print(json.dumps(ctx_result.get("config") or {}, indent=2, default=str))
    if ctx_result.warnings:
        print("\n── Warnings ──────────────────────────────────────────────")
        for w in ctx_result.warnings:
            print(f"  ⚠  {w}")

    print("\n── Workflow dry-run ──────────────────────────────────────")
    coord.execute({"request": request}, dry_run=True)


def _cmd_list() -> None:
    from . import __all__ as exports

    agents = [n for n in exports if n.endswith("Agent")]
    descriptions = {
        "BaseAgent": "Abstract base class; inherit to build custom agents",
        "AgentResult": "Standardised return dataclass (status/data/cost/…)",
        "AgentCoordinator": "Chain agents into workflows with checkpointing",
        "ContextInputAgent": "NL request → structured workflow config",
        "MTLoaderAgent": "EDI/AVG/J → Sites + per-station QC table",
        "DataQCAgent": "SNR, dead-band flags, outlier detection",
        "StaticShiftAgent": "Static-shift detection & correction",
        "PhaseAnalysisAgent": "PT, strike, dimensionality, Mohr, Argand",
        "ForwardModelAgent": "1-D / 2-D / 3-D MT forward modelling",
        "InversionPrepAgent": "Write Occam2D / ModEM data files",
        "InversionEvaluationAgent": "RMS, residual PT, misfit assessment",
        "InterpretationAgent": "Resistivity → lithology interpretation",
        "ReportAgent": "Generate markdown / HTML / PDF report",
        "CodeGenerationAgent": "Emit reproducible Python scripts",
        "WorkflowOrchestratorAgent": "NL → classify + chain agents automatically",
        "DenoisingAgent": "RPCA / Hampel / EMAP / AI denoising",
        "AIInversionAgent": "1-D AI inversion (ResNet / CNN / FCN)",
        "Occam2DAgent": "Write Occam2D data + mesh + startup files",
        "ModEmAgent": "Write ModEM3D impedance data file",
        "AnomalyDetectionAgent": "Unsupervised CAE anomaly flagging",
        "Inv2DAgent": "U-Net 2-D profile AI inversion",
        "Inv3DAgent": "GCN 3-D spatial AI inversion",
        "EnsembleAgent": "Ensemble 1-D inversion + uncertainty bands",
        "JointInversionAgent": "DRCNN multi-modal joint inversion",
        "ModelZooAgent": "Browse / download / run pre-trained models",
        "TensorRotationAgent": "Rotate impedance tensors + write corrected EDIs",
        "EDIExportAgent": "Export processed Sites to EDI files on disk",
        "TipperAnalysisAgent": "Induction arrows, tipper amplitude/phase maps",
        "SensitivityAgent": "Bostick DOI, vertical resolution, sensitivity kernels",
        "FrequencyDecimationAgent": "Intelligent log-spaced period selection for inversion",
        "InversionComparisonAgent": "Side-by-side section comparison with Pearson ρ + RMSE",
        "ResistivityMapAgent": "Horizontal depth-slice maps from 1-D inversions",
        "BatchSurveyAgent": "Parallel pipeline over multiple MT profiles",
        "InversionBackendAgent": "Drive pycsamt.inversion physics backends",
    }
    ordered = ["BaseAgent", "AgentResult", "AgentCoordinator"] + [
        a
        for a in agents
        if a not in {"BaseAgent", "AgentResult", "AgentCoordinator"}
    ]
    print(
        "\n── pycsamt.agents — full catalogue ───────────────────────────────"
    )
    for name in ordered:
        desc = descriptions.get(name, "")
        if desc:
            print(f"  {name:<32s}  {desc}")
    print()


def _cmd_pricing() -> None:
    from ._pricing import PROVIDER_RATES, format_cost

    print("\n── LLM pricing (USD / 1 M tokens) ───────────────────────")
    for provider, models in PROVIDER_RATES.items():
        print(f"\n  {provider.upper()}")
        for model, rates in models.items():
            in_k = format_cost(rates["input"] / 1000)
            out_k = format_cost(rates["output"] / 1000)
            print(f"    {model:<42s}  in={in_k}/1K  out={out_k}/1K")
    print()


def _cmd_zoo(args: list[str]) -> None:
    """Browse or download from the model zoo.

    Usage
    -----
    List all models::

        python -m pycsamt.agents zoo

    Download a checkpoint::

        python -m pycsamt.agents zoo mt1d-resnet-5layer-v1
        python -m pycsamt.agents zoo mt1d-resnet-5layer-v1 --force
    """
    from pycsamt.ai._zoo import (
        download_checkpoint,
        list_pretrained,
    )

    if not args or args[0].startswith("-"):
        # list mode
        models = list_pretrained()
        print(
            f"\n── pycsamt model zoo — {len(models)} pre-trained models ──────────────"
        )
        print(f"  {'Name':<38s}  Description")
        print(f"  {'-' * 37}  {'-' * 55}")
        for name, desc in models.items():
            print(f"  {name:<38s}  {desc[:55]}")
        print(
            "\nTo download: python -m pycsamt.agents zoo <model_name>\n"
            "To predict:  use ModelZooAgent(action='predict') in Python.\n"
        )
        return

    model_name = args[0]
    force = "--force" in args
    print(f"\nDownloading '{model_name}' …")
    try:
        path = download_checkpoint(model_name, force=force, verbose=True)
        print(f"Checkpoint saved: {path}")
    except KeyError as exc:
        print(f"Error: {exc}")
        sys.exit(1)
    except RuntimeError as exc:
        print(f"Download failed: {exc}")
        sys.exit(1)


def main(argv: list[str] | None = None) -> None:
    # explicit empty argv must mean "no arguments" (help), not fall
    # through to sys.argv
    argv = sys.argv[1:] if argv is None else argv

    if not argv or argv[0] in ("-h", "--help", "help"):
        print(__doc__)
        return

    cmd = argv[0].lower()

    if cmd == "preview":
        if len(argv) < 2:
            print("Usage: python -m pycsamt.agents preview '<request>'")
            sys.exit(1)
        _cmd_preview(" ".join(argv[1:]))

    elif cmd == "list":
        _cmd_list()

    elif cmd == "pricing":
        _cmd_pricing()

    elif cmd == "zoo":
        _cmd_zoo(argv[1:])

    elif cmd == "web":
        from .web import launch

        share = "--share" in argv
        port = 7860
        for a in argv[1:]:
            if a.startswith("--port="):
                port = int(a.split("=", 1)[1])
        print(f"Launching pycsamt web app on http://localhost:{port} …")
        launch(share=share, server_port=port)

    else:
        print(
            f"Unknown command: {cmd!r}.  "
            f"Use 'preview', 'list', 'pricing', 'zoo', or 'web'."
        )
        sys.exit(1)


if __name__ == "__main__":
    main()
