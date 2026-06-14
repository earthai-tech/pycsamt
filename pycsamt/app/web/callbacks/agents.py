# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""callbacks/agents.py — AI agent execution panel."""

from __future__ import annotations

import traceback

from dash import Input, Output, State, no_update

from pycsamt.app.desktop.agent_registry import default_params, get_entry
from pycsamt.app.web.cache import cache_get, has_diskcache
from pycsamt.app.web.layout import IDs
from pycsamt.app.web.utils import empty_src, fig_to_src


def register_agents(app) -> None:
    _bg = has_diskcache()
    _cb_kwargs: dict = dict(prevent_initial_call=True)
    if _bg:
        _cb_kwargs["background"] = True
        _cb_kwargs["running"] = [
            (Output(IDs.BTN_RUN_AGENT, "disabled"), True, False),
            (Output(IDs.AGENT_SPINNER, "children"), " ",  ""),
        ]

    @app.callback(
        Output(IDs.AGENT_LOG,     "children"),
        Output(IDs.AGENT_RESULT,  "src"),
        Output(IDs.AGENT_SUMMARY, "children"),
        Output(IDs.TOAST_ERROR,   "is_open",  allow_duplicate=True),
        Output(IDs.TOAST_BODY,    "children", allow_duplicate=True),
        Input(IDs.BTN_RUN_AGENT,  "n_clicks"),
        State(IDs.AGENT_SELECT,   "value"),
        State(IDs.SESSION_ID,     "data"),
        **_cb_kwargs,
    )
    def run_agent(_n_clicks, agent_name, session_id):
        if not agent_name:
            return "No agent selected.", no_update, no_update, False, ""

        if not session_id:
            msg = "Session not initialised — refresh the page."
            return f"⚠ {msg}", no_update, no_update, True, msg

        sites = cache_get(session_id)
        if sites is None:
            msg = "Load survey data before running an agent."
            return f"⚠ {msg}", no_update, no_update, True, msg

        entry      = get_entry(agent_name)
        params     = default_params(agent_name)
        log_lines: list[str] = [f"Running: {agent_name}…"]
        result_src  = empty_src(dark=True)
        summary_txt = ""

        try:
            import inspect
            import matplotlib
            import matplotlib.pyplot as plt
            import pycsamt.emtools as et

            matplotlib.use("Agg")

            if entry["type"] == "processing":
                fn_name  = entry["fn_name"]
                fn       = getattr(et, fn_name)
                accepted = set(inspect.signature(fn).parameters.keys())
                kwargs   = {k: v for k, v in params.items() if k in accepted}
                log_lines.append(f"et.{fn_name}(sites, **params)…")
                result = fn(sites, **kwargs, verbose=0)

                if hasattr(result, "savefig"):
                    result_src = fig_to_src(result)
                elif hasattr(result, "get_figure"):
                    result_src = fig_to_src(result.get_figure())
                else:
                    rp = entry.get("result_plot")
                    if rp and hasattr(et, rp):
                        fig_or_ax = getattr(et, rp)(sites, verbose=0)
                        if hasattr(fig_or_ax, "savefig"):
                            result_src = fig_to_src(fig_or_ax)
                        elif hasattr(fig_or_ax, "get_figure"):
                            result_src = fig_to_src(fig_or_ax.get_figure())

                log_lines.append("Completed.")
                summary_txt = f"✓ {agent_name} completed."

            else:
                import pycsamt.agents as ag
                cls            = getattr(ag, entry["class_name"])
                agent_instance = cls(**{
                    k: v for k, v in params.items()
                    if k in cls.__init__.__code__.co_varnames
                })
                log_lines.append("Executing LLM agent…")
                result = agent_instance.execute({"sites": sites})
                log_lines.append(f"Status: {result.status}")
                if result.warnings:
                    log_lines.extend(f"⚠ {w}" for w in result.warnings)
                summary_txt = result.summary or ""
                if result.llm_interpretation:
                    summary_txt += "\n\n" + result.llm_interpretation[:400]

        except Exception as exc:
            log_lines.append(f"✗ Error: {exc}")
            log_lines.append(traceback.format_exc()[-400:])
            summary_txt = f"Agent failed: {exc}"
            return "\n".join(log_lines), result_src, summary_txt, True, str(exc)

        return "\n".join(log_lines), result_src, summary_txt, False, ""
