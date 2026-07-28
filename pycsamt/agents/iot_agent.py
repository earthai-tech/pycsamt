# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""
pycsamt.agents.iot_agent
========================

:class:`IoTFieldAgent` — AI-assisted monitoring for IoT-enabled AMT / MT /
CSAMT / CSEM field acquisition.

Where :class:`~pycsamt.agents.qc.DataQCAgent` reasons about *processed* EDI
data, this agent reasons about the **live field campaign**: the telemetry
stream coming off edge nodes while the survey is still running.  It is a thin,
defensive orchestration layer over :mod:`pycsamt.iot` that turns whatever
acquisition evidence the caller has — a :class:`~pycsamt.iot.session.FieldSession`,
a bag of telemetry packets, or an existing EDI survey to re-occupy — into:

* a single :class:`~pycsamt.iot.monitoring.MonitoringStatus` (OK / WARN /
  CRITICAL) with the concrete issues that tripped it,
* per-topic telemetry, per-station, and pipeline hand-off tables,
* optional clock-sync and power-budget summaries,
* an optional reproducible :class:`~pycsamt.iot.provenance.AcquisitionManifest`
  (written and, on request, signed),
* dashboard / edge-QC / power / sync figures rendered with the shared pycsamt
  plotting configuration, and
* a plain-English field-ops interpretation from the configured LLM.

Every sub-step is wrapped so that a failure in one diagnostic degrades to a
warning rather than sinking the whole assessment — the agent always returns a
usable :class:`~pycsamt.agents._base.AgentResult`.

Wraps :mod:`pycsamt.iot`:

* :func:`~pycsamt.iot.assess_telemetry` / :class:`~pycsamt.iot.TelemetryMonitor`
* :func:`~pycsamt.iot.telemetry_summary`, :func:`~pycsamt.iot.packet_table`,
  :func:`~pycsamt.iot.monitoring_status_table`
* :func:`~pycsamt.iot.batch_assess_sync`, :func:`~pycsamt.iot.sync_status_table`
* :func:`~pycsamt.iot.estimate_deployment_energy`
* :func:`~pycsamt.iot.build_acquisition_manifest`
* :func:`~pycsamt.iot.field_session_from_edis`,
  :func:`~pycsamt.iot.deployment_report`
* :func:`~pycsamt.iot.plot_field_dashboard`,
  :func:`~pycsamt.iot.plot_edge_qc_summary`,
  :func:`~pycsamt.iot.plot_power_budget`,
  :func:`~pycsamt.iot.plot_sync_quality`
"""

from __future__ import annotations

import time
from typing import Any

from ._base import AgentResult, BaseAgent

_SYSTEM_PROMPT = """\
You are an expert IoT field-operations analyst for pycsamt v2, supervising a
live AMT/MT/CSAMT/CSEM acquisition through its edge telemetry stream.
Given a survey monitoring summary, write 3-5 sentences that:
1. State the overall health of the acquisition (healthy / degraded / critical).
2. Call out the specific stations, channels, or devices that need field
   attention, and the dominant failure mode (packet loss, high latency,
   low battery, clock drift, edge rejections, sensor dropout).
3. Explain the most likely field cause (power budget, GPS/clock, comms link,
   contact resistance, powerline or near-field interference).
4. Recommend the single most important field action to take right now.
Reply in plain English. No bullet points, no markdown headings.
"""


class IoTFieldAgent(BaseAgent):
    r"""Assess and monitor an IoT-enabled EM field acquisition.

    Parameters
    ----------
    api_key, model, llm_provider : str, optional
        LLM configuration.  When ``api_key`` is ``None`` the agent runs fully
        offline and ``llm_interpretation`` is ``None``.
    method : str, optional
        EM method hint (``"amt"``, ``"mt"``, ``"csamt"``, ``"csem"``,
        ``"tdem"``, ``"tem"``).  Used to seed sessions built from packets and
        to tag the manifest.  When omitted it is inferred from the telemetry.
    write_manifest : bool, default False
        When ``True`` an :class:`~pycsamt.iot.provenance.AcquisitionManifest`
        is built for every run (and written when ``manifest_path`` /
        ``output_dir`` is available).

    Input keys
    ----------
    One acquisition source is required (tried in this order):

    ``session`` : FieldSession or mapping
        A live session, or a mapping from
        :meth:`~pycsamt.iot.FieldSession.to_dict`.
    ``packets`` : iterable of TelemetryPacket or mapping
        Raw telemetry to fold into a fresh session.
    ``path`` / ``edis`` / ``sites`` : EDI survey
        An existing survey to seed a re-occupation session via
        :func:`~pycsamt.iot.field_session_from_edis` (no live packets).
    ``deployment`` : DeploymentConfig, optional
        Declared device capabilities; tabulated via
        :func:`~pycsamt.iot.deployment_report`.

    Other optional keys:

    ``survey_id`` : str
        Identifier for sessions built from packets/EDIs (default
        ``"iot_survey"``).
    ``now`` : float
        Reference epoch seconds for live latency / gap calculations.
    ``output_dir`` : str
        Directory for saved figures and (when requested) the manifest.
    ``manifest`` : bool
        Force building the acquisition manifest for this run.
    ``manifest_path`` : str
        Explicit path to write the manifest JSON.
    ``sign_key`` : str or bytes
        When given with a manifest, the written manifest is HMAC-signed.
    ``sync_references`` : mapping
        Per-device reference clocks for :func:`~pycsamt.iot.batch_assess_sync`.
    ``energy_configs`` : iterable of EnergyConfig
        Device energy budgets for :func:`~pycsamt.iot.estimate_deployment_energy`.
    ``figures`` : bool, default True
        Whether to render dashboard/edge/power/sync figures.
    ``api`` : bool
        Passed through to pycsamt table builders (API-object vs raw frame).

    Output data keys
    ----------------
    ``session``           the resolved :class:`~pycsamt.iot.FieldSession`
    ``status``            :class:`~pycsamt.iot.MonitoringStatus`
    ``status_table``      one-row monitoring-status table
    ``telemetry_summary`` packet counts by device and topic
    ``packet_table``      full telemetry packet table
    ``station_table``     registered stations
    ``pipeline_input``    acquisition hand-off dict for the processing flow
    ``deployment_table``  device-capability table (when ``deployment`` given)
    ``sync_table``        clock-sync table (when ``sync_references`` given)
    ``power_table``       energy-budget table (when ``energy_configs`` given)
    ``manifest``          :class:`~pycsamt.iot.AcquisitionManifest` (optional)
    ``manifest_path``     written manifest path (optional)
    ``signature``         manifest HMAC signature dict (when ``sign_key`` given)
    ``level``             monitoring level string (``ok`` / ``warn`` / ``critical``)
    ``issues``            list of issue strings
    ``n_packets`` / ``n_stations`` / ``n_devices`` : int
    ``figures``           dict of matplotlib Figure objects
    ``figure_paths``      dict of saved figure paths (when ``output_dir`` set)

    Examples
    --------
    >>> agent = IoTFieldAgent()
    >>> res = agent.execute({"packets": packets, "output_dir": "/out/iot"})
    >>> res["level"]
    'warn'
    >>> res["status"].issues
    ['battery_min_v below 11.5 V on 2 device(s)']
    >>> res["figures"]["dashboard"]
    <Figure ...>
    """

    SYSTEM_PROMPT = _SYSTEM_PROMPT

    def __init__(
        self,
        *,
        api_key: str | None = None,
        model: str | None = None,
        llm_provider: str = "claude",
        method: str | None = None,
        write_manifest: bool = False,
    ) -> None:
        super().__init__(
            "IoTFieldAgent",
            api_key=api_key,
            model=model,
            llm_provider=llm_provider,
            section_preset="pseudosection",
        )
        self.method = method
        self.write_manifest = write_manifest

    # ── main entry point ──────────────────────────────────────────────────────
    def execute(self, input_data: dict[str, Any]) -> AgentResult:
        self._last_cost = 0.0
        t0 = time.time()
        warnings: list[str] = []

        api = input_data.get("api")
        now = input_data.get("now")

        # ── resolve a FieldSession from whatever the caller supplied ──────────
        session, src = self._resolve_session(input_data, warnings)
        if session is None:
            return AgentResult.failed(
                "No acquisition source: provide 'session', 'packets', or an "
                "EDI 'path'/'edis'/'sites'.",
                hint="e.g. execute({'packets': [...]}) or "
                "execute({'path': '/data/EDIs'}).",
                elapsed=time.time() - t0,
            )

        n_packets = session.n_packets
        n_stations = session.n_stations
        n_devices = session.n_devices

        data: dict[str, Any] = {
            "session": session,
            "source": src,
            "n_packets": n_packets,
            "n_stations": n_stations,
            "n_devices": n_devices,
        }

        # ── monitoring status ─────────────────────────────────────────────────
        status = None
        try:
            status = session.assess(now=now)
            data["status"] = status
            data["level"] = (
                status.level.value
                if hasattr(status.level, "value")
                else str(status.level)
            )
            data["issues"] = list(status.issues)
        except Exception as exc:  # noqa: BLE001
            warnings.append(f"assess: {exc}")
            data["level"] = "unknown"
            data["issues"] = []

        # ── telemetry / station / pipeline tables ─────────────────────────────
        from ..iot import (
            monitoring_status_table,
            packet_table,
            telemetry_summary,
        )

        if status is not None:
            self._try(
                warnings,
                "monitoring_status_table",
                lambda: data.__setitem__(
                    "status_table", monitoring_status_table(status, api=api)
                ),
            )
        if n_packets:
            self._try(
                warnings,
                "telemetry_summary",
                lambda: data.__setitem__(
                    "telemetry_summary",
                    telemetry_summary(session.packets, api=api),
                ),
            )
            self._try(
                warnings,
                "packet_table",
                lambda: data.__setitem__(
                    "packet_table", packet_table(session.packets, api=api)
                ),
            )
        self._try(
            warnings,
            "station_table",
            lambda: data.__setitem__("station_table", session.station_table(api=api)),
        )
        self._try(
            warnings,
            "to_pipeline_input",
            lambda: data.__setitem__("pipeline_input", session.to_pipeline_input()),
        )

        # ── optional deployment capability table ──────────────────────────────
        deployment = input_data.get("deployment")
        if deployment is not None:
            from ..iot import deployment_report

            self._try(
                warnings,
                "deployment_report",
                lambda: data.__setitem__(
                    "deployment_table", deployment_report(deployment, api=api)
                ),
            )

        # ── optional clock-sync summary ───────────────────────────────────────
        sync_refs = input_data.get("sync_references")
        if sync_refs:
            from ..iot import batch_assess_sync

            self._try(
                warnings,
                "batch_assess_sync",
                lambda: data.__setitem__(
                    "sync_table",
                    batch_assess_sync(sync_refs, config=session.sync_config, api=api),
                ),
            )

        # ── optional power-budget summary ─────────────────────────────────────
        energy_cfgs = input_data.get("energy_configs")
        if energy_cfgs:
            from ..iot import estimate_deployment_energy

            self._try(
                warnings,
                "estimate_deployment_energy",
                lambda: data.__setitem__(
                    "power_table",
                    estimate_deployment_energy(energy_cfgs, api=api),
                ),
            )

        # ── optional acquisition manifest (provenance) ────────────────────────
        self._maybe_manifest(session, input_data, data, warnings)

        # ── figures ───────────────────────────────────────────────────────────
        if input_data.get("figures", True):
            self._render_figures(
                session,
                input_data.get("output_dir"),
                data,
                warnings,
                now=now,
            )

        # ── LLM field-ops interpretation ──────────────────────────────────────
        interp = self._interpret(status, data, warnings)

        # ── assemble result ───────────────────────────────────────────────────
        level = data.get("level", "unknown")
        if level == "ok":
            result_status = "success"
        elif level == "unknown":
            result_status = "needs_review"
        else:  # warn / critical
            result_status = "needs_review"

        n_fig = len(data.get("figures", {}))
        n_issue = len(data.get("issues", []))
        summary = (
            f"IoT field assessment: level={level}, "
            f"{n_packets} packet(s) across {n_stations} station(s)/"
            f"{n_devices} device(s), {n_issue} issue(s), "
            f"{n_fig} figure(s)."
        )

        return AgentResult(
            status=result_status,
            summary=summary,
            data=data,
            warnings=warnings,
            llm_interpretation=interp,
            elapsed_seconds=time.time() - t0,
            cost_estimate_usd=self._last_cost,
        )

    # ── source resolution ─────────────────────────────────────────────────────
    def _resolve_session(
        self,
        input_data: dict[str, Any],
        warnings: list[str],
    ) -> tuple[Any, str]:
        """Return ``(FieldSession | None, source_label)`` from *input_data*."""
        from ..iot import FieldSession

        survey_id = str(input_data.get("survey_id") or "iot_survey")

        # 1. an explicit session (object or serialised mapping)
        session = input_data.get("session")
        if session is not None:
            if isinstance(session, FieldSession):
                return session, "session"
            if isinstance(session, dict):
                try:
                    return FieldSession.from_dict(session), "session_dict"
                except Exception as exc:  # noqa: BLE001
                    warnings.append(f"FieldSession.from_dict: {exc}")

        # 2. a bag of telemetry packets → fresh session
        packets = input_data.get("packets")
        if packets:
            try:
                sess = FieldSession(survey_id, method=self.method)
                sess.add_packets(packets)
                return sess, "packets"
            except Exception as exc:  # noqa: BLE001
                warnings.append(f"build session from packets: {exc}")

        # 3. an EDI survey to re-occupy → seeded session (no live packets)
        edis = (
            input_data.get("edis") or input_data.get("path") or input_data.get("sites")
        )
        if edis is not None:
            from ..iot import field_session_from_edis

            try:
                sess = field_session_from_edis(
                    edis,
                    survey_id=survey_id,
                    operator=input_data.get("operator"),
                    method=self.method,
                )
                warnings.append(
                    "Seeded a re-occupation session from EDIs; no live "
                    "telemetry, so monitoring reflects planned stations only."
                )
                return sess, "edis"
            except Exception as exc:  # noqa: BLE001
                warnings.append(f"field_session_from_edis: {exc}")

        return None, "none"

    # ── manifest ──────────────────────────────────────────────────────────────
    def _maybe_manifest(
        self,
        session: Any,
        input_data: dict[str, Any],
        data: dict[str, Any],
        warnings: list[str],
    ) -> None:
        want = (
            self.write_manifest
            or bool(input_data.get("manifest"))
            or bool(input_data.get("manifest_path"))
            or bool(input_data.get("sign_key"))
        )
        if not want:
            return

        manifest = None
        try:
            manifest = session.to_manifest()
            data["manifest"] = manifest
        except Exception as exc:  # noqa: BLE001
            warnings.append(f"to_manifest: {exc}")
            return

        # resolve a write path (explicit, else inside output_dir)
        path = input_data.get("manifest_path")
        if not path and input_data.get("output_dir"):
            import os

            path = os.path.join(
                input_data["output_dir"],
                f"{session.survey_id}_manifest.json",
            )
        if not path:
            return

        sign_key = input_data.get("sign_key")
        try:
            if sign_key:
                data["signature"] = manifest.sign(sign_key)
                data["manifest_path"] = manifest.write_signed(path)
            else:
                data["manifest_path"] = manifest.write(path)
        except Exception as exc:  # noqa: BLE001
            warnings.append(f"write manifest: {exc}")

    # ── figures ───────────────────────────────────────────────────────────────
    def _render_figures(
        self,
        session: Any,
        output_dir: str | None,
        data: dict[str, Any],
        warnings: list[str],
        *,
        now: float | None,
    ) -> None:
        from ..iot import (
            plot_edge_qc_summary,
            plot_field_dashboard,
            plot_power_budget,
            plot_sync_quality,
        )

        figures: dict[str, Any] = {}
        fig_paths: dict[str, str] = {}

        specs = (
            ("dashboard", lambda: plot_field_dashboard(session, now=now)),
            ("edge_qc", lambda: plot_edge_qc_summary(session)),
            ("power_budget", lambda: plot_power_budget(session)),
            ("sync_quality", lambda: plot_sync_quality(session)),
        )
        for name, fn in specs:
            try:
                fig = self._as_fig(fn())
                if fig is None:
                    continue
                figures[name] = fig
                p = self._save_figure(
                    fig, output_dir, f"iot_{name}", warnings_list=warnings
                )
                if p:
                    fig_paths[name] = p
            except Exception as exc:  # noqa: BLE001
                warnings.append(f"plot {name}: {exc}")

        data["figures"] = figures
        data["figure_paths"] = fig_paths

    # ── LLM interpretation ─────────────────────────────────────────────────────
    def _interpret(
        self,
        status: Any,
        data: dict[str, Any],
        warnings: list[str],
    ) -> str | None:
        if not self.api_key or status is None:
            return None
        prompt = (
            "IoT field acquisition monitoring summary:\n"
            f"  Level: {data.get('level')}\n"
            f"  Packets: {data.get('n_packets')} across "
            f"{data.get('n_stations')} station(s) / "
            f"{data.get('n_devices')} device(s)\n"
            f"  Packet success rate: {getattr(status, 'packet_success_rate', '?')}\n"
            f"  Edge acceptance rate: {getattr(status, 'edge_acceptance_rate', '?')}\n"
            f"  Mean latency (s): {getattr(status, 'mean_latency_s', '?')}\n"
            f"  Max gap (s): {getattr(status, 'max_gap_s', '?')}\n"
            f"  Min battery (V): {getattr(status, 'battery_min_v', '?')}\n"
            f"  Max clock offset (ms): {getattr(status, 'clock_offset_max_ms', '?')}\n"
            f"  Methods: {getattr(status, 'methods', [])}\n"
            f"  Issues: {list(getattr(status, 'issues', []))[:6]}\n\n"
            "Assess the acquisition health and recommend the single most "
            "important field action."
        )
        try:
            return self.query_llm(prompt, max_tokens=260)
        except Exception as exc:  # noqa: BLE001
            warnings.append(f"llm interpretation: {exc}")
            return None

    # ── small helpers ──────────────────────────────────────────────────────────
    @staticmethod
    def _try(warnings: list[str], label: str, fn: Any) -> None:
        """Run *fn*, appending ``"<label>: <err>"`` to *warnings* on failure."""
        try:
            fn()
        except Exception as exc:  # noqa: BLE001
            warnings.append(f"{label}: {exc}")

    @staticmethod
    def _as_fig(obj: Any) -> Any:
        """Normalise a plot return (Figure or Axes) to a Figure."""
        if obj is None:
            return None
        if hasattr(obj, "get_figure"):
            fig = obj.get_figure()
            if fig is not None:
                return fig
        return obj


__all__ = ["IoTFieldAgent"]
