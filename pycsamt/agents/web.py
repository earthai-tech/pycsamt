"""
pycsamt.agents.web
===================

Gradio web interface for the pycsamt agent system.

Launch
------
From Python::

    from pycsamt.agents.web import launch

    launch()  # opens http://localhost:7860

From the CLI::

    python -m pycsamt.agents web

Tabs
----
**Chat (Orchestrator)**
    Type a free-text request. The :class:`WorkflowOrchestratorAgent`
    classifies the workflow, builds the chain, and runs it.

**Load & QC**
    Browse to an EDI directory.  :class:`MTLoaderAgent` + :class:`DataQCAgent`
    display the per-station quality table and confidence section figure.

**Static Shift**
    Choose correction method.  :class:`StaticShiftAgent` shows before/after
    comparison.

**Phase Analysis**
    :class:`PhaseAnalysisAgent` runs the full tensor suite and displays
    figures inline.

**Forward Modelling**
    Enter a layered model (resistivities and thicknesses).
    :class:`ForwardModelAgent` plots the synthetic response.

**AI Inversion (1-D)**
    Choose architecture and epochs.  :class:`AIInversionAgent` trains and
    predicts, displaying the model section.

**2-D AI Inversion**
    :class:`Inv2DAgent` (U-Net) inverts the full profile at once, capturing
    lateral continuity.

**Ensemble Inversion**
    :class:`EnsembleAgent` — multiple independent models produce mean + σ
    uncertainty bands with conformal calibration.

**Joint Inversion**
    :class:`JointInversionAgent` (DRCNN) fuses MT with a secondary modality
    (TEM, CSAMT, or synthesised low-frequency sub-band).

**Model Zoo**
    :class:`ModelZooAgent` — list and load pre-trained checkpoints from the
    earthai-tech model zoo; run zero-shot prediction on any EDI dataset.

**Report**
    Consolidate all results into a downloadable markdown / HTML report.

Requires ``gradio >= 4.0``.  Install with::

    pip install gradio
"""

from __future__ import annotations

import os
import tempfile
from typing import Any

# ── lazy gradio import ────────────────────────────────────────────────────────


def _require_gradio():
    try:
        import gradio as gr

        return gr
    except ImportError:
        raise ImportError(
            "Gradio is required for the web interface. "
            "Install with: pip install gradio"
        ) from None


# ── shared state ──────────────────────────────────────────────────────────────

_SESSION: dict[str, Any] = {}


def _get_agent(name: str, api_key: str = "", llm_provider: str = "claude"):
    """Instantiate or retrieve a cached agent."""
    key = f"{name}_{llm_provider}"
    if key not in _SESSION or (_SESSION.get(f"{key}_key") != api_key):
        from . import (
            AIInversionAgent,
            DataQCAgent,
            EnsembleAgent,
            ForwardModelAgent,
            Inv2DAgent,
            Inv3DAgent,
            JointInversionAgent,
            ModelZooAgent,
            MTLoaderAgent,
            PhaseAnalysisAgent,
            ReportAgent,
            StaticShiftAgent,
        )
        from .orchestrator import WorkflowOrchestratorAgent

        kw = dict(
            api_key=api_key or None,
            llm_provider=llm_provider,
        )
        registry = {
            "orchestrator": WorkflowOrchestratorAgent(**kw),
            "loader": MTLoaderAgent(**kw),
            "qc": DataQCAgent(**kw),
            "static_shift": StaticShiftAgent(**kw),
            "phase_analysis": PhaseAnalysisAgent(**kw),
            "forward": ForwardModelAgent(**kw),
            "ai_inversion": AIInversionAgent(**kw),
            "inv2d": Inv2DAgent(**kw),
            "inv3d": Inv3DAgent(**kw),
            "ensemble": EnsembleAgent(**kw),
            "joint": JointInversionAgent(**kw),
            "zoo": ModelZooAgent(**kw),
            "report": ReportAgent(**kw),
        }
        _SESSION.update(
            {f"{n}_{llm_provider}": v for n, v in registry.items()}
        )
        _SESSION[f"{key}_key"] = api_key
    return _SESSION.get(key)


# ── tab handlers ──────────────────────────────────────────────────────────────


def _chat_run(
    request: str,
    data_path: str,
    output_dir: str,
    api_key: str,
    llm_provider: str,
):
    """Chat tab: route request through WorkflowOrchestratorAgent."""
    orch = _get_agent("orchestrator", api_key, llm_provider)
    r = orch.execute(
        {
            "request": request,
            "data_path": data_path,
            "output_dir": output_dir or tempfile.mkdtemp(prefix="pycsamt_"),
            "dry_run": False,
        }
    )
    log = (
        f"**Status:** {r.status}\n\n"
        f"**Workflow:** {r.get('workflow_type', '?')}\n\n"
        f"**Summary:** {r.summary}\n\n"
    )
    if r.warnings:
        log += (
            "**Warnings:**\n"
            + "\n".join(f"- {w}" for w in r.warnings[:5])
            + "\n\n"
        )
    if r.llm_interpretation:
        log += f"**LLM Interpretation:**\n{r.llm_interpretation}\n"
    return log


def _load_qc_run(data_path: str, api_key: str, llm_provider: str):
    """Load & QC tab."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if not data_path or not os.path.exists(data_path):
        return "⚠ Path does not exist.", None, None

    out_dir = tempfile.mkdtemp(prefix="pycsamt_qc_")

    loader = _get_agent("loader", api_key, llm_provider)
    qc_ag = _get_agent("qc", api_key, llm_provider)

    lr = loader.execute({"path": data_path, "output_dir": out_dir})
    if lr.status == "failed":
        return f"⚠ Loader failed: {lr.error}", None, None

    qr = qc_ag.execute({"sites": lr["sites"], "output_dir": out_dir})

    # summary text
    stats = lr.get("summary_stats") or {}
    txt = (
        f"**{lr['n_stations']} stations loaded.** "
        f"Mean QC score: {stats.get('mean_qc_score', '?'):.0f}/100. "
        f"Flagged: {qr['n_flagged']}.\n\n"
    )
    if qr.get("flagged_stations"):
        txt += f"Flagged stations: {', '.join(qr['flagged_stations'])}\n"
    if qr.warnings:
        txt += "\n".join(f"- {w}" for w in qr.warnings[:4])

    # figures
    figs = qr.get("figures", {})
    fig1 = fig2 = None
    if "confidence_section" in figs:
        p1 = os.path.join(out_dir, "conf_section.png")
        figs["confidence_section"].savefig(p1, dpi=100, bbox_inches="tight")
        plt.close(figs["confidence_section"])
        fig1 = p1
    if "confidence_profile" in figs:
        p2 = os.path.join(out_dir, "conf_profile.png")
        figs["confidence_profile"].savefig(p2, dpi=100, bbox_inches="tight")
        plt.close(figs["confidence_profile"])
        fig2 = p2

    _SESSION["_last_sites"] = lr["sites"]
    _SESSION["_last_load"] = lr
    _SESSION["_last_qc"] = qr
    return txt, fig1, fig2


def _ss_run(method: str, api_key: str, llm_provider: str):
    """Static-shift tab."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sites = _SESSION.get("_last_sites")
    if sites is None:
        return "⚠ Load data first (Load & QC tab).", None

    out_dir = tempfile.mkdtemp(prefix="pycsamt_ss_")
    ag = _get_agent("static_shift", api_key, llm_provider)
    r = ag.execute({"sites": sites, "method": method, "output_dir": out_dir})

    _SESSION["_last_ss"] = r

    txt = f"**{r.summary}**\n\n"
    ds = r.get("delta_stats") or {}
    if ds:
        txt += (
            f"Mean correction: {ds.get('mean', 0):.3f} log₁₀.  "
            f"Max: {ds.get('max', 0):.3f}.  "
            f"Stations shifted: {ds.get('n_shifted', 0)}.\n"
        )
    if r.warnings:
        txt += "\n".join(f"- {w}" for w in r.warnings[:4])

    fig_path = None
    figs = r.get("figures", {})
    if "ss_comparison" in figs:
        p = os.path.join(out_dir, "ss_comparison.png")
        figs["ss_comparison"].savefig(p, dpi=100, bbox_inches="tight")
        plt.close(figs["ss_comparison"])
        fig_path = p

    return txt, fig_path


def _pt_run(api_key: str, llm_provider: str):
    """Phase analysis tab."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ss_r = _SESSION.get("_last_ss")
    sites = (ss_r.get("corrected_sites") if ss_r else None) or _SESSION.get(
        "_last_sites"
    )
    if sites is None:
        return "⚠ Load data first.", None, None, None

    out_dir = tempfile.mkdtemp(prefix="pycsamt_pt_")
    ag = _get_agent("phase_analysis", api_key, llm_provider)
    r = ag.execute({"sites": sites, "output_dir": out_dir})

    _SESSION["_last_pt"] = r

    txt = (
        f"**{r.summary}**\n\n"
        f"Strike: {r.get('strike_consensus', float('nan')):.1f}° "
        f"± {r.get('strike_iqr', float('nan')):.1f}°\n"
        f"1-D/2-D/3-D: {r.get('n_1d', 0)}/{r.get('n_2d', 0)}/{r.get('n_3d', 0)}"
    )
    if r.llm_interpretation:
        txt += f"\n\n**LLM:** {r.llm_interpretation}"

    figs = r.get("figures", {})
    saved = []
    for name in ["pt_psection", "strike_analysis", "survey_fingerprint"]:
        if name in figs:
            p = os.path.join(out_dir, f"{name}.png")
            figs[name].savefig(p, dpi=100, bbox_inches="tight")
            plt.close(figs[name])
            saved.append(p)
    while len(saved) < 3:
        saved.append(None)
    return txt, saved[0], saved[1], saved[2]


def _fwd_run(rhos_str: str, ths_str: str, api_key: str, llm_provider: str):
    """Forward modelling tab."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = tempfile.mkdtemp(prefix="pycsamt_fwd_")
    try:
        rhos = [float(x.strip()) for x in rhos_str.split(",")]
        ths = [float(x.strip()) for x in ths_str.split(",")]
    except Exception:
        return (
            "⚠ Invalid resistivity or thickness format. Use comma-separated numbers.",
            None,
        )

    ag = _get_agent("forward", api_key, llm_provider)
    r = ag.execute(
        {
            "model": {"resistivity": rhos, "thickness": ths},
            "output_dir": out_dir,
        }
    )

    txt = f"**{r.summary}**"
    if r.warnings:
        txt += "\n" + "\n".join(f"- {w}" for w in r.warnings[:3])

    fig_path = None
    figs = r.get("figures", {})
    for name in ["response_and_model", "response", "model"]:
        if name in figs:
            p = os.path.join(out_dir, f"{name}.png")
            figs[name].savefig(p, dpi=100, bbox_inches="tight")
            plt.close(figs[name])
            fig_path = p
            break

    return txt, fig_path


def _ai_inv_run(arch: str, epochs: int, api_key: str, llm_provider: str):
    """AI inversion tab."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sites = (_SESSION.get("_last_ss", {}) or {}).get(
        "corrected_sites"
    ) or _SESSION.get("_last_sites")
    if sites is None:
        return "⚠ Load data first.", None

    out_dir = tempfile.mkdtemp(prefix="pycsamt_ai_")
    ag = AIInversionAgent_lazy(
        arch=arch,
        epochs=int(epochs),
        api_key=api_key or None,
        llm_provider=llm_provider,
    )
    r = ag.execute({"sites": sites, "output_dir": out_dir})

    txt = f"**{r.summary}**"
    if r.llm_interpretation:
        txt += f"\n\n**LLM:** {r.llm_interpretation}"
    if r.warnings:
        txt += "\n" + "\n".join(f"- {w}" for w in r.warnings[:3])

    fig_path = None
    figs = r.get("figures", {})
    if "ai_section" in figs:
        p = os.path.join(out_dir, "ai_section.png")
        figs["ai_section"].savefig(p, dpi=100, bbox_inches="tight")
        plt.close(figs["ai_section"])
        fig_path = p

    return txt, fig_path


def AIInversionAgent_lazy(**kwargs):
    from .ai_inversion import AIInversionAgent

    return AIInversionAgent(**kwargs)


def _inv2d_run(epochs: int, api_key: str, llm_provider: str):
    """2-D AI inversion tab (U-Net)."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sites = (_SESSION.get("_last_ss", {}) or {}).get(
        "corrected_sites"
    ) or _SESSION.get("_last_sites")
    if sites is None:
        return "⚠ Load data first (Load & QC tab).", None

    out_dir = tempfile.mkdtemp(prefix="pycsamt_inv2d_")
    ag = _get_agent("inv2d", api_key, llm_provider)
    # update epochs without rebuilding entire registry
    ag.epochs = int(epochs)
    r = ag.execute({"sites": sites, "output_dir": out_dir})

    _SESSION["_last_inv2d"] = r
    txt = f"**{r.summary}**"
    if r.llm_interpretation:
        txt += f"\n\n**LLM:** {r.llm_interpretation}"
    if r.warnings:
        txt += "\n" + "\n".join(f"- {w}" for w in r.warnings[:3])

    fig_path = None
    figs = r.get("figures", {})
    if "inv2d_section" in figs:
        p = os.path.join(out_dir, "inv2d_section.png")
        figs["inv2d_section"].savefig(p, dpi=100, bbox_inches="tight")
        plt.close(figs["inv2d_section"])
        fig_path = p

    return txt, fig_path


def _ensemble_run(
    n_estimators: int, epochs: int, api_key: str, llm_provider: str
):
    """Ensemble inversion tab."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sites = (_SESSION.get("_last_ss", {}) or {}).get(
        "corrected_sites"
    ) or _SESSION.get("_last_sites")
    if sites is None:
        return "⚠ Load data first.", None, None

    out_dir = tempfile.mkdtemp(prefix="pycsamt_ens_")
    ag = _get_agent("ensemble", api_key, llm_provider)
    ag.n_estimators = int(n_estimators)
    ag.epochs = int(epochs)
    r = ag.execute({"sites": sites, "output_dir": out_dir})

    _SESSION["_last_ensemble"] = r
    rms = r.get("rms_global", float("nan"))
    cov = r.get("coverage", float("nan"))
    txt = (
        f"**{r.summary}**\n\nRMS: {rms:.3f}   90% coverage: {cov:.1%}"
        if not (rms != rms)
        else f"**{r.summary}**"
    )
    if r.llm_interpretation:
        txt += f"\n\n**LLM:** {r.llm_interpretation}"

    f1 = f2 = None
    figs = r.get("figures", {})
    for key, pname in [
        ("uncertainty_section", "ens_section"),
        ("uncertainty_profile", "ens_profile"),
    ]:
        if key in figs:
            p = os.path.join(out_dir, f"{pname}.png")
            figs[key].savefig(p, dpi=100, bbox_inches="tight")
            plt.close(figs[key])
            if f1 is None:
                f1 = p
            else:
                f2 = p

    return txt, f1, f2


def _joint_run(modality: str, epochs: int, api_key: str, llm_provider: str):
    """Joint inversion tab."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sites = (_SESSION.get("_last_ss", {}) or {}).get(
        "corrected_sites"
    ) or _SESSION.get("_last_sites")
    if sites is None:
        return "⚠ Load data first.", None

    out_dir = tempfile.mkdtemp(prefix="pycsamt_joint_")
    ag = _get_agent("joint", api_key, llm_provider)
    ag.modalities = ["mt", modality]
    ag.epochs = int(epochs)
    r = ag.execute({"sites": sites, "output_dir": out_dir})

    _SESSION["_last_joint"] = r
    txt = f"**{r.summary}**"
    if r.llm_interpretation:
        txt += f"\n\n**LLM:** {r.llm_interpretation}"

    fig_path = None
    figs = r.get("figures", {})
    if "joint_section" in figs:
        p = os.path.join(out_dir, "joint_section.png")
        figs["joint_section"].savefig(p, dpi=100, bbox_inches="tight")
        plt.close(figs["joint_section"])
        fig_path = p

    return txt, fig_path


def _inv3d_run(
    epochs: int, radius_km: float, n_mc: int, api_key: str, llm_provider: str
):
    """3-D GCN inversion tab."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sites = (_SESSION.get("_last_ss", {}) or {}).get(
        "corrected_sites"
    ) or _SESSION.get("_last_sites")
    if sites is None:
        return "⚠ Load data first (Load & QC tab).", None, None

    out_dir = tempfile.mkdtemp(prefix="pycsamt_inv3d_")
    ag = _get_agent("inv3d", api_key, llm_provider)
    ag.epochs = int(epochs)
    ag.radius = float(radius_km) * 1000.0
    ag.n_mc = int(n_mc)
    r = ag.execute({"sites": sites, "output_dir": out_dir})

    _SESSION["_last_inv3d"] = r
    rms_str = f"{r.get('rms_global', float('nan')):.3f}"
    n_edges = r.get("n_edges", "?")
    txt = f"**{r.summary}**\n\nGraph edges: {n_edges}   RMS: {rms_str}"
    if r.llm_interpretation:
        txt += f"\n\n**LLM:** {r.llm_interpretation}"
    if r.warnings:
        txt += "\n" + "\n".join(f"- {w}" for w in r.warnings[:4])

    f1 = f2 = None
    figs = r.get("figures", {})
    for key, pname in [
        ("depth_slices", "inv3d_slices"),
        ("resistivity_section", "inv3d_section"),
        ("uncertainty_map", "inv3d_uncertainty"),
    ]:
        if key in figs:
            p = os.path.join(out_dir, f"{pname}.png")
            figs[key].savefig(p, dpi=100, bbox_inches="tight")
            plt.close(figs[key])
            if f1 is None:
                f1 = p
            elif f2 is None:
                f2 = p

    return txt, f1, f2


def _zoo_list_run(api_key: str, llm_provider: str):
    """Model zoo — list available models."""
    ag = _get_agent("zoo", api_key, llm_provider)
    r = ag.execute({"action": "list"})
    rows = r.get("details", [])
    table = "| Model | Arch | Layers | Solver | Description |\n|---|---|---|---|---|\n"
    for row in rows:
        table += (
            f"| `{row['name']}` | {row['arch']} | {row['n_layers']} "
            f"| {row['solver']} | {row['description'][:55]} |\n"
        )
    return table


def _zoo_predict_run(model_name: str, api_key: str, llm_provider: str):
    """Model zoo — predict with a named pre-trained model."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if not model_name.strip():
        return "⚠ Enter a model name (from the table above).", None

    sites = (_SESSION.get("_last_ss", {}) or {}).get(
        "corrected_sites"
    ) or _SESSION.get("_last_sites")
    if sites is None:
        return "⚠ Load data first (Load & QC tab).", None

    out_dir = tempfile.mkdtemp(prefix="pycsamt_zoo_")
    ag = _get_agent("zoo", api_key, llm_provider)
    r = ag.execute(
        {
            "action": "predict",
            "model_name": model_name.strip(),
            "sites": sites,
            "output_dir": out_dir,
        }
    )

    txt = f"**{r.summary}**"
    if r.llm_interpretation:
        txt += f"\n\n**LLM:** {r.llm_interpretation}"
    if r.warnings:
        txt += "\n" + "\n".join(f"- {w}" for w in r.warnings[:4])

    fig_path = None
    figs = r.get("figures", {})
    for key in ["ai_section", "convergence"]:
        if key in figs:
            p = os.path.join(out_dir, f"{key}.png")
            figs[key].savefig(p, dpi=100, bbox_inches="tight")
            plt.close(figs[key])
            fig_path = p
            break

    return txt, fig_path


def _report_run(title: str, api_key: str, llm_provider: str):
    """Report tab."""
    out_dir = tempfile.mkdtemp(prefix="pycsamt_report_")

    results = {}
    for key, sess_key in [
        ("load", "_last_load"),
        ("qc", "_last_qc"),
        ("static_shift", "_last_ss"),
        ("phase_analysis", "_last_pt"),
    ]:
        v = _SESSION.get(sess_key)
        if v is not None:
            results[key] = v

    if not results:
        return "⚠ Run at least Load & QC first.", None

    ag = _get_agent("report", api_key, llm_provider)
    r = ag.execute(
        {
            "results": results,
            "output_dir": out_dir,
            "title": title or "pycsamt MT Survey Report",
        }
    )

    return f"**{r.summary}**", r.get("report_path_md")


# ── main build function ───────────────────────────────────────────────────────


def build_app():
    """Build and return the Gradio Blocks app."""
    gr = _require_gradio()

    _LOGO = "# 🌐 pycsamt — AI-assisted MT/AMT Processing"
    _LLM_HINT = "Leave blank to run without LLM (regex-only mode)."

    with gr.Blocks(
        title="pycsamt Agent System", theme=gr.themes.Soft()
    ) as app:
        gr.Markdown(_LOGO)

        # ── shared LLM config (sidebar-style row) ─────────────────────────────
        with gr.Row():
            api_key_box = gr.Textbox(
                label="LLM API key",
                type="password",
                placeholder=_LLM_HINT,
                scale=3,
            )
            llm_prov_box = gr.Dropdown(
                [
                    "claude",
                    "openai",
                    "gemini",
                    "deepseek",
                    "minimax",
                ],
                value="claude",
                label="Provider",
                scale=1,
            )

        with gr.Tabs():
            # ── Tab 1: Chat / Orchestrator ─────────────────────────────────────
            with gr.TabItem("💬 Chat (Orchestrator)"):
                gr.Markdown(
                    "Describe your MT workflow in plain English. "
                    "The orchestrator selects and chains the right agents."
                )
                req_box = gr.Textbox(
                    label="Request",
                    lines=3,
                    placeholder="Load EDIs from /data/L22PLT, "
                    "run QC and phase tensor analysis…",
                )
                dp_box = gr.Textbox(label="Data path (optional override)")
                out_box = gr.Textbox(
                    label="Output directory", value="pycsamt_output"
                )
                chat_btn = gr.Button("Run workflow", variant="primary")
                chat_out = gr.Markdown()

                chat_btn.click(
                    fn=_chat_run,
                    inputs=[
                        req_box,
                        dp_box,
                        out_box,
                        api_key_box,
                        llm_prov_box,
                    ],
                    outputs=chat_out,
                )

            # ── Tab 2: Load & QC ───────────────────────────────────────────────
            with gr.TabItem("📂 Load & QC"):
                gr.Markdown("Load EDI / AVG / J files and run frequency QC.")
                lqc_path = gr.Textbox(label="EDI directory or file path")
                lqc_btn = gr.Button("Load + QC", variant="primary")
                lqc_txt = gr.Markdown()
                lqc_f1 = gr.Image(label="Confidence section", type="filepath")
                lqc_f2 = gr.Image(label="Confidence profile", type="filepath")

                lqc_btn.click(
                    fn=_load_qc_run,
                    inputs=[lqc_path, api_key_box, llm_prov_box],
                    outputs=[lqc_txt, lqc_f1, lqc_f2],
                )

            # ── Tab 3: Static Shift ────────────────────────────────────────────
            with gr.TabItem("⚡ Static Shift"):
                gr.Markdown(
                    "Detect and correct galvanic static shift.  "
                    "Load data first in the **Load & QC** tab."
                )
                ss_meth = gr.Dropdown(
                    ["ama", "loess", "refmedian", "bilateral"],
                    value="ama",
                    label="Correction method",
                )
                ss_btn = gr.Button("Correct", variant="primary")
                ss_txt = gr.Markdown()
                ss_fig = gr.Image(
                    label="Before / after comparison", type="filepath"
                )

                ss_btn.click(
                    fn=_ss_run,
                    inputs=[ss_meth, api_key_box, llm_prov_box],
                    outputs=[ss_txt, ss_fig],
                )

            # ── Tab 4: Phase Analysis ──────────────────────────────────────────
            with gr.TabItem("🔬 Phase Analysis"):
                gr.Markdown(
                    "Phase tensor, strike, and dimensionality analysis.  "
                    "Uses corrected sites from Static Shift if available."
                )
                pt_btn = gr.Button("Run analysis", variant="primary")
                pt_txt = gr.Markdown()
                pt_f1 = gr.Image(label="PT pseudosection", type="filepath")
                pt_f2 = gr.Image(label="Strike analysis", type="filepath")
                pt_f3 = gr.Image(label="Survey fingerprint", type="filepath")

                pt_btn.click(
                    fn=_pt_run,
                    inputs=[api_key_box, llm_prov_box],
                    outputs=[pt_txt, pt_f1, pt_f2, pt_f3],
                )

            # ── Tab 5: Forward Modelling ───────────────────────────────────────
            with gr.TabItem("📡 Forward Model"):
                gr.Markdown(
                    "Run a 1-D MT forward model. "
                    "Enter resistivities and thicknesses as comma-separated numbers."
                )
                fwd_rhos = gr.Textbox(
                    label="Resistivities (Ω·m)",
                    value="200, 20, 5000, 150",
                )
                fwd_ths = gr.Textbox(
                    label="Thicknesses (m)",
                    value="300, 800, 3000",
                )
                fwd_btn = gr.Button("Run forward", variant="primary")
                fwd_txt = gr.Markdown()
                fwd_fig = gr.Image(label="Response + model", type="filepath")

                fwd_btn.click(
                    fn=_fwd_run,
                    inputs=[fwd_rhos, fwd_ths, api_key_box, llm_prov_box],
                    outputs=[fwd_txt, fwd_fig],
                )

            # ── Tab 6: AI Inversion ────────────────────────────────────────────
            with gr.TabItem("🤖 AI Inversion"):
                gr.Markdown(
                    "Train a 1-D MT neural inverter on synthetic data "
                    "then predict for each loaded station."
                )
                ai_arch = gr.Dropdown(
                    ["resnet", "cnn1d", "fcn"],
                    value="resnet",
                    label="Architecture",
                )
                ai_epochs = gr.Slider(
                    10, 200, value=30, step=10, label="Training epochs"
                )
                ai_btn = gr.Button("Train + Predict", variant="primary")
                ai_txt = gr.Markdown()
                ai_fig = gr.Image(
                    label="Predicted model section", type="filepath"
                )

                ai_btn.click(
                    fn=_ai_inv_run,
                    inputs=[ai_arch, ai_epochs, api_key_box, llm_prov_box],
                    outputs=[ai_txt, ai_fig],
                )

            # ── Tab 7: 2-D AI Inversion ───────────────────────────────────────
            with gr.TabItem("🗺 2-D AI Inversion"):
                gr.Markdown(
                    "U-Net profile inversion — captures lateral continuity.  "
                    "Uses corrected sites from Static Shift if available."
                )
                inv2d_epochs = gr.Slider(
                    10, 100, value=30, step=5, label="Training epochs"
                )
                inv2d_btn = gr.Button("Run U-Net inversion", variant="primary")
                inv2d_txt = gr.Markdown()
                inv2d_fig = gr.Image(
                    label="2-D resistivity section", type="filepath"
                )

                inv2d_btn.click(
                    fn=_inv2d_run,
                    inputs=[inv2d_epochs, api_key_box, llm_prov_box],
                    outputs=[inv2d_txt, inv2d_fig],
                )

            # ── Tab 8: 3-D GCN Inversion ──────────────────────────────────────
            with gr.TabItem("🌐 3-D GCN Inversion"):
                gr.Markdown(
                    "Graph-convolutional 3-D spatial inversion.  "
                    "Station proximity is encoded as a graph; GCN message-passing "
                    "enforces spatial coherence across the survey."
                )
                inv3d_epochs = gr.Slider(
                    10, 100, value=30, step=5, label="Training epochs"
                )
                inv3d_radius = gr.Slider(
                    1.0,
                    50.0,
                    value=5.0,
                    step=1.0,
                    label="Adjacency radius (km)",
                )
                inv3d_nmc = gr.Slider(
                    0,
                    50,
                    value=20,
                    step=5,
                    label="MC dropout samples (0 = skip)",
                )
                inv3d_btn = gr.Button("Run GCN inversion", variant="primary")
                inv3d_txt = gr.Markdown()
                inv3d_f1 = gr.Image(label="Depth slices", type="filepath")
                inv3d_f2 = gr.Image(
                    label="Resistivity section / uncertainty", type="filepath"
                )

                inv3d_btn.click(
                    fn=_inv3d_run,
                    inputs=[
                        inv3d_epochs,
                        inv3d_radius,
                        inv3d_nmc,
                        api_key_box,
                        llm_prov_box,
                    ],
                    outputs=[inv3d_txt, inv3d_f1, inv3d_f2],
                )

            # ── Tab 9: Ensemble Inversion ──────────────────────────────────────
            with gr.TabItem("📊 Ensemble Inversion"):
                gr.Markdown(
                    "Train *N* independent 1-D inverters and report mean ± σ "
                    "uncertainty bands with conformal calibration."
                )
                ens_n = gr.Slider(
                    2, 10, value=3, step=1, label="Number of estimators"
                )
                ens_ep = gr.Slider(
                    10, 100, value=20, step=5, label="Epochs per estimator"
                )
                ens_btn = gr.Button("Train ensemble", variant="primary")
                ens_txt = gr.Markdown()
                ens_f1 = gr.Image(
                    label="Mean + uncertainty section", type="filepath"
                )
                ens_f2 = gr.Image(
                    label="Uncertainty profile (first station)",
                    type="filepath",
                )

                ens_btn.click(
                    fn=_ensemble_run,
                    inputs=[ens_n, ens_ep, api_key_box, llm_prov_box],
                    outputs=[ens_txt, ens_f1, ens_f2],
                )

            # ── Tab 10: Joint Inversion ────────────────────────────────────────
            with gr.TabItem("🔗 Joint Inversion"):
                gr.Markdown(
                    "DRCNN multi-modal joint inversion.  The MT dataset is fused "
                    "with a secondary modality (or a synthesised low-frequency proxy)."
                )
                joint_mod = gr.Dropdown(
                    ["tem", "csamt", "gravity", "seismic"],
                    value="tem",
                    label="Secondary modality",
                )
                joint_ep = gr.Slider(
                    10, 100, value=20, step=5, label="Training epochs"
                )
                joint_btn = gr.Button("Run joint inversion", variant="primary")
                joint_txt = gr.Markdown()
                joint_fig = gr.Image(
                    label="Joint resistivity section", type="filepath"
                )

                joint_btn.click(
                    fn=_joint_run,
                    inputs=[joint_mod, joint_ep, api_key_box, llm_prov_box],
                    outputs=[joint_txt, joint_fig],
                )

            # ── Tab 11: Model Zoo ──────────────────────────────────────────────
            with gr.TabItem("🏛 Model Zoo"):
                gr.Markdown(
                    "Browse and deploy pre-trained EM inverters from the "
                    "**earthai-tech** model zoo.  "
                    "Weights are cached in `~/.pycsamt/model_zoo/`."
                )
                zoo_list_btn = gr.Button(
                    "Refresh model list", variant="secondary"
                )
                zoo_table = gr.Markdown()

                gr.Markdown("---")
                gr.Markdown(
                    "**Zero-shot prediction** — load a named model and "
                    "predict on the currently loaded sites."
                )
                zoo_name_box = gr.Textbox(
                    label="Model name",
                    placeholder="mt1d-resnet-5layer-v1",
                )
                zoo_pred_btn = gr.Button(
                    "Download + Predict", variant="primary"
                )
                zoo_txt = gr.Markdown()
                zoo_fig = gr.Image(
                    label="Predicted model section", type="filepath"
                )

                zoo_list_btn.click(
                    fn=_zoo_list_run,
                    inputs=[api_key_box, llm_prov_box],
                    outputs=zoo_table,
                )
                zoo_pred_btn.click(
                    fn=_zoo_predict_run,
                    inputs=[zoo_name_box, api_key_box, llm_prov_box],
                    outputs=[zoo_txt, zoo_fig],
                )

            # ── Tab 12: Report ─────────────────────────────────────────────────
            with gr.TabItem("📄 Report"):
                gr.Markdown(
                    "Assemble all results into a markdown / HTML report."
                )
                rep_title = gr.Textbox(
                    label="Report title", value="MT Survey Report"
                )
                rep_btn = gr.Button("Generate report", variant="primary")
                rep_txt = gr.Markdown()
                rep_file = gr.File(label="Download report (.md)")

                rep_btn.click(
                    fn=_report_run,
                    inputs=[rep_title, api_key_box, llm_prov_box],
                    outputs=[rep_txt, rep_file],
                )

    return app


# ── public entry points ───────────────────────────────────────────────────────


def launch(
    *,
    share: bool = False,
    server_port: int = 7860,
    server_name: str = "0.0.0.0",
    **kwargs: Any,
) -> None:
    """Launch the Gradio web app.

    Parameters
    ----------
    share : bool
        Create a public Gradio share link.
    server_port : int
    server_name : str
    **kwargs
        Passed to :meth:`gradio.Blocks.launch`.
    """
    app = build_app()
    app.launch(
        share=share,
        server_port=server_port,
        server_name=server_name,
        **kwargs,
    )


__all__ = ["build_app", "launch"]
