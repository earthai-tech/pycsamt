"""Visualisation helpers for IoT-enabled field acquisition.

The plotting layer turns telemetry, station metadata, edge QC, power, and
clock information into compact operational figures. These figures are not
geophysical inversions; they are acquisition dashboards that help explain
what happened before EDI/impedance processing.
"""

from __future__ import annotations

import math
from collections import defaultdict
from collections.abc import Iterable, Mapping
from typing import (
    Any,
)

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from .core import PacketKind, TelemetryPacket
from .edge import EdgeProcessingResult
from .power import (
    EnergyConfig,
    EnergyEstimate,
    estimate_energy_budget,
)
from .session import FieldSession
from .sync import SyncStatus

__all__ = [
    "plot_field_dashboard",
    "plot_edge_qc_summary",
    "plot_power_budget",
    "plot_sync_quality",
]


_LEVEL_COLORS = {
    "ok": "#2ca25f",
    "sustaining": "#2ca25f",
    "excellent": "#2ca25f",
    "good": "#74c476",
    "warning": "#fdae6b",
    "fair": "#fdae6b",
    "critical": "#de2d26",
    "poor": "#de2d26",
    "reject": "#de2d26",
    "no_data": "#9e9e9e",
    "unknown": "#9e9e9e",
}


def plot_field_dashboard(
    session: FieldSession | Mapping[str, Any],
    *,
    now: float | None = None,
    figsize: tuple[float, float] = (13.0, 8.0),
    station_axis: str = "auto",
    title: str | None = None,
    output_path: str | None = None,
    close: bool = False,
) -> Any:
    """Plot a compact IoT acquisition dashboard.

    Parameters
    ----------
    session : FieldSession or mapping
        Field session, or a mapping produced by
        :meth:`pycsamt.iot.FieldSession.to_dict`.
    now : float, optional
        Reference timestamp used for live latency/status calculations.
    figsize : tuple
        Matplotlib figure size in inches.
    station_axis : {"auto", "profile", "map"}
        Station layout. ``"profile"`` uses profile chainage/index;
        ``"map"`` uses longitude/latitude when all stations have
        coordinates; ``"auto"`` chooses map only when coordinates exist.
    title : str, optional
        Figure title. Defaults to the survey id.
    output_path : str, optional
        If given, save the figure to this path.
    close : bool
        Close the figure before returning it. Useful for batch report
        generation after saving.

    Returns
    -------
    matplotlib.figure.Figure
        The dashboard figure. The computed data are also attached as
        ``fig.pycsamt_iot_dashboard`` for reproducible report workflows.
    """
    sess = _as_session(session)
    dashboard = _dashboard_data(sess, now=now)

    fig, axes = plt.subplots(
        2,
        2,
        figsize=figsize,
        constrained_layout=True,
    )
    fig.suptitle(title or f"IoT field dashboard: {sess.survey_id}", fontsize=14)

    _plot_station_panel(axes[0, 0], dashboard, station_axis=station_axis)
    _plot_acceptance_panel(axes[0, 1], dashboard)
    _plot_operations_panel(axes[1, 0], dashboard)
    _plot_timeline_panel(axes[1, 1], dashboard)

    handles = [
        Line2D([0], [0], marker="o", linestyle="", color=color, label=label)
        for label, color in [
            ("ok", _LEVEL_COLORS["ok"]),
            ("warning", _LEVEL_COLORS["warning"]),
            ("critical/reject", _LEVEL_COLORS["critical"]),
            ("unknown", _LEVEL_COLORS["unknown"]),
        ]
    ]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False)
    fig.pycsamt_iot_dashboard = dashboard  # type: ignore[attr-defined]
    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
    if close:
        plt.close(fig)
    return fig


def plot_edge_qc_summary(
    edge: EdgeProcessingResult
    | TelemetryPacket
    | FieldSession
    | Mapping[str, Any]
    | Iterable[EdgeProcessingResult | TelemetryPacket | Mapping[str, Any]],
    *,
    figsize: tuple[float, float] = (12.0, 7.5),
    title: str = "Edge QC summary",
    output_path: str | None = None,
    close: bool = False,
) -> Any:
    """Plot edge quality-control decisions and channel metrics.

    Parameters
    ----------
    edge : EdgeProcessingResult, TelemetryPacket, FieldSession, mapping, or iterable
        Edge-processing result(s), QC telemetry packet(s), a field session,
        or serialised mappings. Session inputs are filtered to QC packets.
    figsize : tuple
        Matplotlib figure size in inches.
    title : str
        Figure title.
    output_path : str, optional
        If given, save the figure to this path.
    close : bool
        Close the figure before returning it.

    Returns
    -------
    matplotlib.figure.Figure
        The QC summary figure. The normalised rows are attached as
        ``fig.pycsamt_iot_edge_qc``.
    """
    rows = _edge_qc_rows(edge)
    fig, axes = plt.subplots(
        2,
        2,
        figsize=figsize,
        constrained_layout=True,
    )
    fig.suptitle(title, fontsize=14)
    _plot_qc_decisions(axes[0, 0], rows)
    _plot_qc_coverage(axes[0, 1], rows)
    _plot_qc_spikes(axes[1, 0], rows)
    _plot_qc_reasons(axes[1, 1], rows)
    fig.pycsamt_iot_edge_qc = rows  # type: ignore[attr-defined]
    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
    if close:
        plt.close(fig)
    return fig


def plot_power_budget(
    power: EnergyConfig
    | EnergyEstimate
    | TelemetryPacket
    | FieldSession
    | Mapping[str, Any]
    | Iterable[EnergyConfig | EnergyEstimate | TelemetryPacket | Mapping[str, Any]],
    *,
    figsize: tuple[float, float] = (12.0, 7.5),
    title: str = "IoT power budget",
    output_path: str | None = None,
    close: bool = False,
) -> Any:
    """Plot IoT energy budget, runtime, and power-state summaries.

    Parameters
    ----------
    power : EnergyConfig, EnergyEstimate, TelemetryPacket, FieldSession, mapping, or iterable
        Power budget input(s). ``EnergyConfig`` objects are estimated before
        plotting. Session inputs are filtered to ``PacketKind.POWER`` packets.
    figsize : tuple
        Matplotlib figure size in inches.
    title : str
        Figure title.
    output_path : str, optional
        If given, save the figure to this path.
    close : bool
        Close the figure before returning it.

    Returns
    -------
    matplotlib.figure.Figure
        The power-budget figure. Normalised rows are attached as
        ``fig.pycsamt_iot_power_budget``.
    """
    rows = _power_rows(power)
    fig, axes = plt.subplots(
        2,
        2,
        figsize=figsize,
        constrained_layout=True,
    )
    fig.suptitle(title, fontsize=14)
    _plot_power_load_harvest(axes[0, 0], rows)
    _plot_power_runtime(axes[0, 1], rows)
    _plot_power_breakdown(axes[1, 0], rows)
    _plot_power_states(axes[1, 1], rows)
    fig.pycsamt_iot_power_budget = rows  # type: ignore[attr-defined]
    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
    if close:
        plt.close(fig)
    return fig


def plot_sync_quality(
    sync: SyncStatus
    | TelemetryPacket
    | FieldSession
    | Mapping[str, Any]
    | Iterable[SyncStatus | TelemetryPacket | Mapping[str, Any]],
    *,
    figsize: tuple[float, float] = (12.0, 7.5),
    title: str = "Clock synchronisation quality",
    tolerance_ms: float | None = 1.0,
    max_drift_ppm: float | None = None,
    max_jitter_ms: float | None = None,
    output_path: str | None = None,
    close: bool = False,
) -> Any:
    """Plot clock offset, drift, jitter, GPS lock, and quality grades.

    Parameters
    ----------
    sync : SyncStatus, TelemetryPacket, FieldSession, mapping, or iterable
        Synchronisation status input(s). Session inputs are filtered to
        ``PacketKind.SYNC`` packets.
    figsize : tuple
        Matplotlib figure size in inches.
    title : str
        Figure title.
    tolerance_ms, max_drift_ppm, max_jitter_ms : float, optional
        Optional visual threshold lines.
    output_path : str, optional
        If given, save the figure to this path.
    close : bool
        Close the figure before returning it.

    Returns
    -------
    matplotlib.figure.Figure
        The synchronisation figure. Normalised rows are attached as
        ``fig.pycsamt_iot_sync_quality``.
    """
    rows = _sync_rows(sync)
    fig, axes = plt.subplots(
        2,
        2,
        figsize=figsize,
        constrained_layout=True,
    )
    fig.suptitle(title, fontsize=14)
    _plot_sync_offset(axes[0, 0], rows, tolerance_ms=tolerance_ms)
    _plot_sync_drift_jitter(
        axes[0, 1],
        rows,
        max_drift_ppm=max_drift_ppm,
        max_jitter_ms=max_jitter_ms,
    )
    _plot_sync_quality_counts(axes[1, 0], rows)
    _plot_sync_reference_points(axes[1, 1], rows)
    fig.pycsamt_iot_sync_quality = rows  # type: ignore[attr-defined]
    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
    if close:
        plt.close(fig)
    return fig


def _as_session(session: FieldSession | Mapping[str, Any]) -> FieldSession:
    if isinstance(session, FieldSession):
        return session
    if isinstance(session, Mapping):
        return FieldSession.from_dict(session)
    raise TypeError("session must be a FieldSession or session mapping.")


def _sync_rows(sync: Any) -> list[dict[str, Any]]:
    items = _sync_items(sync)
    rows: list[dict[str, Any]] = []
    for index, item in enumerate(items):
        rows.append(_sync_row(item, index=index))
    return rows


def _sync_items(sync: Any) -> list[Any]:
    if isinstance(sync, FieldSession):
        return [packet for packet in sync.packets if packet.kind is PacketKind.SYNC]
    if isinstance(sync, (SyncStatus, TelemetryPacket)):
        return [sync]
    if isinstance(sync, Mapping):
        if "packets" in sync and "survey_id" in sync:
            return _sync_items(FieldSession.from_dict(sync))
        if "payload" in sync:
            return [TelemetryPacket(**dict(sync))]
        return [SyncStatus(**dict(sync))]
    if isinstance(sync, Iterable) and not isinstance(sync, (str, bytes)):
        out: list[Any] = []
        for item in sync:
            out.extend(_sync_items(item))
        return out
    raise TypeError(
        "sync must be a SyncStatus, sync packet, FieldSession, mapping, "
        "or iterable of those objects."
    )


def _sync_row(item: Any, *, index: int) -> dict[str, Any]:
    if isinstance(item, SyncStatus):
        row = item.as_dict()
    elif isinstance(item, TelemetryPacket):
        payload = dict(item.payload or {})
        row = dict(payload)
        row.setdefault("device_id", item.device_id)
    else:
        row = dict(item)
    row.setdefault("device_id", f"device-{index + 1}")
    row.setdefault("reference", "gps")
    row.setdefault("quality", "unknown")
    row["offset_ms"] = _float_or_nan(row.get("offset_ms", row.get("clock_offset_ms")))
    row["drift_ppm"] = _float_or_nan(row.get("drift_ppm"))
    row["jitter_ms"] = _float_or_nan(row.get("jitter_ms"))
    row["n_reference_points"] = int(_finite_or_zero(row.get("n_reference_points")))
    row["gps_lock"] = _optional_bool(row.get("gps_lock"))
    row["within_tolerance"] = _optional_bool(row.get("within_tolerance"))
    row["quality"] = str(row.get("quality") or "unknown").lower()
    row["reference"] = str(row.get("reference") or "gps").lower()
    return row


def _plot_sync_offset(
    ax: Any,
    rows: list[Mapping[str, Any]],
    *,
    tolerance_ms: float | None,
) -> None:
    ax.set_title("Clock offset")
    if not rows:
        _empty_panel(ax, "No sync data")
        return
    labels = _sync_labels(rows)
    values = [_float_or_nan(row.get("offset_ms")) for row in rows]
    if not any(math.isfinite(v) for v in values):
        _empty_panel(ax, "No offset data")
        return
    x = np.arange(len(rows), dtype=float)
    colors = [_sync_quality_color(str(row.get("quality"))) for row in rows]
    plot_values = [v if math.isfinite(v) else 0.0 for v in values]
    ax.bar(x, plot_values, color=colors, edgecolor="black", linewidth=0.4)
    if tolerance_ms is not None:
        tol = abs(float(tolerance_ms))
        ax.axhline(tol, color="#de2d26", lw=1.0, ls="--")
        ax.axhline(-tol, color="#de2d26", lw=1.0, ls="--")
    ax.axhline(0.0, color="#333333", lw=0.8)
    ax.set_ylabel("Offset (ms)")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.grid(True, axis="y", alpha=0.25)


def _plot_sync_drift_jitter(
    ax: Any,
    rows: list[Mapping[str, Any]],
    *,
    max_drift_ppm: float | None,
    max_jitter_ms: float | None,
) -> None:
    ax.set_title("Drift and jitter")
    if not rows:
        _empty_panel(ax, "No sync data")
        return
    labels = _sync_labels(rows)
    x = np.arange(len(rows), dtype=float)
    width = 0.38
    drift = [_finite_or_zero(row.get("drift_ppm")) for row in rows]
    jitter = [_finite_or_zero(row.get("jitter_ms")) for row in rows]
    if not any(v != 0 for v in drift + jitter):
        _empty_panel(ax, "No drift/jitter data")
        return
    ax.bar(x - width / 2, drift, width, label="drift ppm", color="#756bb1")
    ax.bar(x + width / 2, jitter, width, label="jitter ms", color="#6baed6")
    if max_drift_ppm is not None:
        ax.axhline(float(max_drift_ppm), color="#756bb1", lw=1.0, ls="--")
        ax.axhline(-float(max_drift_ppm), color="#756bb1", lw=1.0, ls="--")
    if max_jitter_ms is not None:
        ax.axhline(float(max_jitter_ms), color="#6baed6", lw=1.0, ls=":")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_ylabel("ppm / ms")
    ax.legend(frameon=False)
    ax.grid(True, axis="y", alpha=0.25)


def _plot_sync_quality_counts(ax: Any, rows: list[Mapping[str, Any]]) -> None:
    ax.set_title("Quality grades")
    if not rows:
        _empty_panel(ax, "No quality data")
        return
    grades = [str(row.get("quality", "unknown")) for row in rows]
    counts = {grade: grades.count(grade) for grade in sorted(set(grades))}
    labels = list(counts.keys())
    colors = [_sync_quality_color(label) for label in labels]
    ax.bar(labels, list(counts.values()), color=colors, edgecolor="black")
    ax.set_ylabel("Device count")
    ax.tick_params(axis="x", rotation=25)
    ax.grid(True, axis="y", alpha=0.25)


def _plot_sync_reference_points(
    ax: Any,
    rows: list[Mapping[str, Any]],
) -> None:
    ax.set_title("Reference support and GPS lock")
    if not rows:
        _empty_panel(ax, "No reference data")
        return
    labels = _sync_labels(rows)
    points = [_finite_or_zero(row.get("n_reference_points")) for row in rows]
    gps = [row.get("gps_lock") for row in rows]
    colors = [
        _LEVEL_COLORS["ok"]
        if value is True
        else _LEVEL_COLORS["critical"]
        if value is False
        else _LEVEL_COLORS["unknown"]
        for value in gps
    ]
    x = np.arange(len(rows), dtype=float)
    ax.bar(x, points, color=colors, edgecolor="black", linewidth=0.4)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_ylabel("Reference points")
    ax.grid(True, axis="y", alpha=0.25)
    handles = [
        Line2D(
            [0],
            [0],
            marker="s",
            linestyle="",
            color=_LEVEL_COLORS["ok"],
            label="GPS lock",
        ),
        Line2D(
            [0],
            [0],
            marker="s",
            linestyle="",
            color=_LEVEL_COLORS["critical"],
            label="GPS lost",
        ),
        Line2D(
            [0],
            [0],
            marker="s",
            linestyle="",
            color=_LEVEL_COLORS["unknown"],
            label="unknown",
        ),
    ]
    ax.legend(handles=handles, frameon=False, fontsize=8, loc="best")


def _sync_labels(rows: list[Mapping[str, Any]]) -> list[str]:
    return [
        str(row.get("device_id") or f"device-{idx + 1}") for idx, row in enumerate(rows)
    ]


def _sync_quality_color(quality: str) -> str:
    quality = str(quality or "unknown").lower()
    return _LEVEL_COLORS.get(quality, "#9e9e9e")


def _optional_bool(value: Any) -> bool | None:
    if value is None or value == "":
        return None
    if isinstance(value, str):
        text = value.strip().lower()
        if text in {"1", "true", "yes", "y", "ok", "locked"}:
            return True
        if text in {"0", "false", "no", "n", "lost", "none"}:
            return False
    return bool(value)


def _power_rows(power: Any) -> list[dict[str, Any]]:
    items = _power_items(power)
    rows: list[dict[str, Any]] = []
    for index, item in enumerate(items):
        rows.append(_power_row(item, index=index))
    return rows


def _power_items(power: Any) -> list[Any]:
    if isinstance(power, FieldSession):
        return [packet for packet in power.packets if packet.kind is PacketKind.POWER]
    if isinstance(power, (EnergyConfig, EnergyEstimate, TelemetryPacket)):
        return [power]
    if isinstance(power, Mapping):
        if "packets" in power and "survey_id" in power:
            return _power_items(FieldSession.from_dict(power))
        if "payload" in power:
            return [TelemetryPacket(**dict(power))]
        if "battery_wh" in power and "active_power_w" in power:
            return [EnergyConfig(**dict(power))]
        return [EnergyEstimate(**dict(power))]
    if isinstance(power, Iterable) and not isinstance(power, (str, bytes)):
        out: list[Any] = []
        for item in power:
            out.extend(_power_items(item))
        return out
    raise TypeError(
        "power must be an EnergyConfig, EnergyEstimate, power packet, "
        "FieldSession, mapping, or iterable of those objects."
    )


def _power_row(item: Any, *, index: int) -> dict[str, Any]:
    device_id: str | None = None
    if isinstance(item, EnergyConfig):
        device_id = item.device_id
        estimate = estimate_energy_budget(item)
        row = estimate.as_dict()
    elif isinstance(item, EnergyEstimate):
        row = item.as_dict()
    elif isinstance(item, TelemetryPacket):
        device_id = item.device_id
        row = dict(item.payload or {})
    else:
        row = dict(item)
    row = dict(row)
    row.setdefault("device_id", device_id or f"device-{index + 1}")
    row.setdefault("state", "unknown")
    row.setdefault("issues", "")
    for key in (
        "average_power_w",
        "net_wh_per_day",
        "load_wh_per_day",
        "harvest_wh_per_day",
        "usable_battery_wh",
        "telemetry_wh_per_day",
        "edge_wh_per_day",
        "auxiliary_wh_per_day",
        "reserve_wh",
        "energy_margin_wh_per_day",
    ):
        row[key] = _float_or_nan(row.get(key))
    for key in ("runtime_hours", "runtime_days", "autonomy_days_no_harvest"):
        row[key] = _float_allow_inf(row.get(key))
    row["state"] = str(row.get("state") or "unknown").lower()
    row["issues"] = _split_reasons(row.get("issues"))
    return row


def _plot_power_load_harvest(ax: Any, rows: list[Mapping[str, Any]]) -> None:
    ax.set_title("Daily load and harvest")
    if not rows:
        _empty_panel(ax, "No power data")
        return
    labels = _power_labels(rows)
    x = np.arange(len(rows), dtype=float)
    width = 0.38
    load = [_finite_or_zero(row.get("load_wh_per_day")) for row in rows]
    harvest = [_finite_or_zero(row.get("harvest_wh_per_day")) for row in rows]
    ax.bar(x - width / 2, load, width, label="load", color="#756bb1")
    ax.bar(x + width / 2, harvest, width, label="harvest", color="#74c476")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_ylabel("Wh/day")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(frameon=False)


def _plot_power_runtime(ax: Any, rows: list[Mapping[str, Any]]) -> None:
    ax.set_title("Runtime and no-harvest autonomy")
    if not rows:
        _empty_panel(ax, "No runtime data")
        return
    labels = _power_labels(rows)
    x = np.arange(len(rows), dtype=float)
    width = 0.38
    runtime = [_plot_runtime_value(row.get("runtime_days")) for row in rows]
    autonomy = [
        _plot_runtime_value(row.get("autonomy_days_no_harvest")) for row in rows
    ]
    colors = [_power_state_color(str(row.get("state"))) for row in rows]
    ax.bar(x - width / 2, runtime, width, label="runtime", color=colors)
    ax.bar(
        x + width / 2,
        autonomy,
        width,
        label="no-harvest autonomy",
        color="#9ecae1",
    )
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_ylabel("Days")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(frameon=False)
    for idx, row in enumerate(rows):
        if math.isinf(_float_or_nan(row.get("runtime_days"))):
            ax.text(idx - width / 2, runtime[idx], "inf", ha="center", fontsize=8)


def _plot_power_breakdown(ax: Any, rows: list[Mapping[str, Any]]) -> None:
    ax.set_title("Daily load breakdown")
    if not rows:
        _empty_panel(ax, "No load data")
        return
    labels = _power_labels(rows)
    x = np.arange(len(rows), dtype=float)
    telemetry = np.asarray(
        [_finite_or_zero(row.get("telemetry_wh_per_day")) for row in rows]
    )
    edge = np.asarray([_finite_or_zero(row.get("edge_wh_per_day")) for row in rows])
    aux = np.asarray([_finite_or_zero(row.get("auxiliary_wh_per_day")) for row in rows])
    load = np.asarray([_finite_or_zero(row.get("load_wh_per_day")) for row in rows])
    base = np.maximum(load - telemetry - edge - aux, 0.0)
    bottom = np.zeros(len(rows))
    for values, label, color in [
        (base, "base", "#bcbddc"),
        (telemetry, "telemetry", "#6baed6"),
        (edge, "edge", "#fd8d3c"),
        (aux, "auxiliary", "#969696"),
    ]:
        ax.bar(x, values, bottom=bottom, label=label, color=color)
        bottom = bottom + values
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_ylabel("Wh/day")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(frameon=False, fontsize=8)


def _plot_power_states(ax: Any, rows: list[Mapping[str, Any]]) -> None:
    ax.set_title("Power states and issues")
    if not rows:
        _empty_panel(ax, "No states")
        return
    states = [str(row.get("state", "unknown")) for row in rows]
    counts = {state: states.count(state) for state in sorted(set(states))}
    labels = list(counts.keys())
    colors = [_power_state_color(label) for label in labels]
    ax.bar(labels, list(counts.values()), color=colors, edgecolor="black")
    ax.set_ylabel("Device count")
    ax.tick_params(axis="x", rotation=25)
    issue_counts: dict[str, int] = defaultdict(int)
    for row in rows:
        for issue in row.get("issues") or []:
            issue_counts[str(issue)] += 1
    if issue_counts:
        text = "\n".join(
            f"{key}: {value}"
            for key, value in sorted(
                issue_counts.items(), key=lambda kv: (-kv[1], kv[0])
            )[:5]
        )
        ax.text(
            0.98,
            0.95,
            text,
            transform=ax.transAxes,
            va="top",
            ha="right",
            fontsize=8,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#cccccc"),
        )
    ax.grid(True, axis="y", alpha=0.25)


def _power_labels(rows: list[Mapping[str, Any]]) -> list[str]:
    return [
        str(row.get("device_id") or f"device-{i + 1}") for i, row in enumerate(rows)
    ]


def _finite_or_zero(value: Any) -> float:
    out = _float_or_nan(value)
    return out if math.isfinite(out) else 0.0


def _plot_runtime_value(value: Any) -> float:
    out = _float_allow_inf(value)
    if math.isinf(out):
        return 1.1 * max(1.0, out if math.isfinite(out) else 1.0)
    return out if math.isfinite(out) else 0.0


def _power_state_color(state: str) -> str:
    return _LEVEL_COLORS.get(str(state or "unknown").lower(), "#9e9e9e")


def _float_allow_inf(value: Any) -> float:
    try:
        return float(value)
    except Exception:
        return float("nan")


def _edge_qc_rows(edge: Any) -> list[dict[str, Any]]:
    items = _edge_items(edge)
    rows: list[dict[str, Any]] = []
    for index, item in enumerate(items):
        if isinstance(item, EdgeProcessingResult):
            rows.extend(_rows_from_edge_result(item, result_index=index))
        else:
            rows.extend(_rows_from_qc_packet(item, result_index=index))
    return rows


def _edge_items(edge: Any) -> list[Any]:
    if isinstance(edge, FieldSession):
        return [packet for packet in edge.packets if packet.kind is PacketKind.QC]
    if isinstance(edge, EdgeProcessingResult):
        return [edge]
    if isinstance(edge, TelemetryPacket):
        return [edge]
    if isinstance(edge, Mapping):
        if "packets" in edge and "survey_id" in edge:
            return _edge_items(FieldSession.from_dict(edge))
        if "payload" in edge:
            return [TelemetryPacket(**dict(edge))]
        if "metrics" in edge:
            return [EdgeProcessingResult(**dict(edge))]
        raise TypeError("edge mapping must describe a session, packet, or result.")
    if isinstance(edge, Iterable) and not isinstance(edge, (str, bytes)):
        out: list[Any] = []
        for item in edge:
            out.extend(_edge_items(item))
        return out
    raise TypeError(
        "edge must be an EdgeProcessingResult, QC packet, FieldSession, "
        "mapping, or iterable of those objects."
    )


def _rows_from_edge_result(
    result: EdgeProcessingResult,
    *,
    result_index: int,
) -> list[dict[str, Any]]:
    metrics = dict(result.metrics or {})
    base = dict(
        result_index=result_index,
        station=None,
        decision=result.decision.value,
        accepted=bool(result.accepted),
        finite_coverage=_float_or_nan(metrics.get("finite_coverage")),
        spike_fraction=_float_or_nan(metrics.get("spike_fraction_max")),
        rms=_float_or_nan(metrics.get("rms")),
        reasons=list(result.reasons or []),
        warnings=_split_reasons(metrics.get("warnings")),
    )
    rows: list[dict[str, Any]] = []
    for channel in result.channels:
        row = dict(base)
        row.update(
            channel=channel.channel,
            finite_coverage=channel.finite_coverage,
            spike_fraction=channel.spike_fraction,
            rms=channel.rms,
            channel_accepted=channel.accepted,
            channel_reasons=list(channel.reasons or []),
        )
        rows.append(row)
    if not rows:
        row = dict(base)
        row.update(
            channel="window",
            channel_accepted=bool(result.accepted),
            channel_reasons=list(result.reasons or []),
        )
        rows.append(row)
    return rows


def _rows_from_qc_packet(
    packet: TelemetryPacket | Mapping[str, Any],
    *,
    result_index: int,
) -> list[dict[str, Any]]:
    pkt = (
        packet
        if isinstance(packet, TelemetryPacket)
        else TelemetryPacket(**dict(packet))
    )
    payload = dict(pkt.payload or {})
    metrics = dict(payload.get("metrics") or {})
    decision = str(payload.get("decision", "") or "").lower()
    if not decision:
        accepted = _accepted_from_payload(payload)
        decision = (
            "accept"
            if accepted is True
            else "reject"
            if accepted is False
            else "unknown"
        )
    accepted = _accepted_from_payload(payload)
    if accepted is None:
        accepted = decision in {"accept", "ok", "pass", "warning"}
    station = _payload_first(payload, "station", "site", "station_id")
    reasons = _split_reasons(payload.get("reasons"))
    warnings = _split_reasons(metrics.get("warnings") or payload.get("warnings"))
    base = dict(
        result_index=result_index,
        station=station,
        decision=decision,
        accepted=bool(accepted),
        finite_coverage=_float_or_nan(
            metrics.get("finite_coverage", payload.get("finite_coverage"))
        ),
        spike_fraction=_float_or_nan(
            metrics.get("spike_fraction_max", payload.get("spike_fraction"))
        ),
        rms=_float_or_nan(metrics.get("rms", payload.get("rms"))),
        reasons=reasons,
        warnings=warnings,
    )
    channels = payload.get("channels")
    rows: list[dict[str, Any]] = []
    if isinstance(channels, list) and channels and isinstance(channels[0], Mapping):
        for channel in channels:
            row = dict(base)
            row.update(
                channel=str(channel.get("channel", "channel")).lower(),
                finite_coverage=_float_or_nan(
                    channel.get("finite_coverage", base["finite_coverage"])
                ),
                spike_fraction=_float_or_nan(
                    channel.get("spike_fraction", base["spike_fraction"])
                ),
                rms=_float_or_nan(channel.get("rms", base["rms"])),
                channel_accepted=bool(channel.get("accepted", accepted)),
                channel_reasons=_split_reasons(channel.get("reasons")),
            )
            rows.append(row)
    else:
        channel_names = (
            channels
            if isinstance(channels, list)
            else ([channels] if isinstance(channels, str) else ["window"])
        )
        for channel in channel_names:
            row = dict(base)
            row.update(
                channel=str(channel).lower(),
                channel_accepted=bool(accepted),
                channel_reasons=list(reasons),
            )
            rows.append(row)
    return rows


def _split_reasons(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, str):
        return [part for part in value.split(";") if part]
    if isinstance(value, Iterable) and not isinstance(value, (bytes, str)):
        return [str(part) for part in value if str(part)]
    return [str(value)]


def _plot_qc_decisions(ax: Any, rows: list[Mapping[str, Any]]) -> None:
    ax.set_title("QC decisions")
    if not rows:
        _empty_panel(ax, "No QC rows")
        return
    decisions = [str(row.get("decision", "unknown")) for row in rows]
    counts = {
        decision: decisions.count(decision) for decision in sorted(set(decisions))
    }
    labels = list(counts.keys())
    colors = [_decision_color(label) for label in labels]
    ax.bar(
        labels,
        list(counts.values()),
        color=colors,
        edgecolor="black",
        linewidth=0.5,
    )
    ax.set_ylabel("Channel/window count")
    ax.tick_params(axis="x", rotation=25)
    for idx, value in enumerate(counts.values()):
        ax.text(idx, value + 0.05, str(value), ha="center", fontsize=8)


def _plot_qc_coverage(ax: Any, rows: list[Mapping[str, Any]]) -> None:
    ax.set_title("Finite coverage by channel")
    _plot_qc_metric(
        ax,
        rows,
        key="finite_coverage",
        ylabel="Finite coverage",
        ylim=(0, 1.05),
        thresholds=(0.85, 0.95),
    )


def _plot_qc_spikes(ax: Any, rows: list[Mapping[str, Any]]) -> None:
    ax.set_title("Spike fraction by channel")
    _plot_qc_metric(
        ax,
        rows,
        key="spike_fraction",
        ylabel="Spike fraction",
        ylim=(0, None),
        thresholds=(0.05,),
        high_bad=True,
    )


def _plot_qc_metric(
    ax: Any,
    rows: list[Mapping[str, Any]],
    *,
    key: str,
    ylabel: str,
    ylim: tuple[float, float | None],
    thresholds: tuple[float, ...] = (),
    high_bad: bool = False,
) -> None:
    values = [_float_or_nan(row.get(key)) for row in rows]
    if not any(math.isfinite(v) for v in values):
        _empty_panel(ax, f"No {ylabel.lower()} data")
        return
    labels = _qc_labels(rows)
    plot_values = [v if math.isfinite(v) else 0.0 for v in values]
    colors = [_decision_color(str(row.get("decision", "unknown"))) for row in rows]
    x = np.arange(len(rows), dtype=float)
    ax.bar(x, plot_values, color=colors, edgecolor="black", linewidth=0.4)
    for threshold in thresholds:
        color = "#de2d26" if high_bad else "#2ca25f"
        ax.axhline(threshold, color=color, lw=1.0, ls="--")
    ax.set_ylabel(ylabel)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right")
    if ylim[1] is None:
        top = max(max(plot_values) * 1.25, max(thresholds or (0.05,)) * 1.2, 0.1)
        ax.set_ylim(ylim[0], top)
    else:
        ax.set_ylim(*ylim)
    ax.grid(True, axis="y", alpha=0.25)


def _plot_qc_reasons(ax: Any, rows: list[Mapping[str, Any]]) -> None:
    ax.set_title("Reasons and warnings")
    reason_counts: dict[str, int] = defaultdict(int)
    for row in rows:
        for reason in list(row.get("reasons") or []) + list(
            row.get("channel_reasons") or []
        ):
            reason_counts[str(reason)] += 1
        for warning in row.get("warnings") or []:
            reason_counts[f"warn:{warning}"] += 1
    if not reason_counts:
        _empty_panel(ax, "No rejection reasons")
        return
    labels = sorted(reason_counts, key=lambda k: (-reason_counts[k], k))[:10]
    values = [reason_counts[label] for label in labels]
    y = np.arange(len(labels), dtype=float)
    ax.barh(y, values, color="#9ecae1", edgecolor="black", linewidth=0.4)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlabel("Count")
    ax.grid(True, axis="x", alpha=0.25)


def _qc_labels(rows: list[Mapping[str, Any]]) -> list[str]:
    labels = []
    for row in rows:
        station = row.get("station")
        channel = row.get("channel", "window")
        idx = row.get("result_index")
        prefix = f"{station}:" if station else f"{idx}:"
        labels.append(f"{prefix}{channel}")
    return labels


def _decision_color(decision: str) -> str:
    decision = str(decision or "unknown").lower()
    if decision in {"accept", "ok", "pass"}:
        return _LEVEL_COLORS["ok"]
    if decision in {"warning", "repeat"}:
        return _LEVEL_COLORS["warning"]
    if decision in {"reject", "critical", "fail"}:
        return _LEVEL_COLORS["critical"]
    return _LEVEL_COLORS["unknown"]


def _dashboard_data(
    session: FieldSession,
    *,
    now: float | None,
) -> dict[str, Any]:
    pipeline = session.to_pipeline_input()
    stations = [dict(row) for row in pipeline.get("stations", [])]
    packets = [_packet_row(packet) for packet in session.packets]
    status = session.assess(now=now)
    latest = _latest_station_metrics(packets, session)
    for row in stations:
        sid = row.get("station_id")
        row.update(latest.get(sid, {}))
        row["health_level"] = _station_level(row)
    return dict(
        survey_id=session.survey_id,
        method=session.method or pipeline.get("method"),
        n_devices=session.n_devices,
        n_stations=session.n_stations,
        n_packets=session.n_packets,
        stations=stations,
        packets=packets,
        monitoring=status.as_dict(),
        issues=list(status.issues),
    )


def _packet_row(packet: TelemetryPacket) -> dict[str, Any]:
    row = packet.as_dict()
    payload = dict(row.get("payload") or {})
    row["kind"] = (
        row["kind"].value
        if isinstance(row.get("kind"), PacketKind)
        else str(row.get("kind"))
    )
    row["station"] = _payload_first(
        payload,
        "station",
        "site",
        "station_id",
        "station_name",
    )
    row["accepted"] = _accepted_from_payload(payload)
    row["battery_v"] = _float_or_nan(
        _payload_first(payload, "battery_v", "battery_voltage_v")
    )
    row["clock_offset_ms"] = _float_or_nan(
        _payload_first(payload, "clock_offset_ms", "offset_ms")
    )
    row["runtime_days"] = _float_or_nan(payload.get("runtime_days"))
    row["power_state"] = str(payload.get("state", "") or "").lower()
    row["sync_quality"] = str(payload.get("quality", "") or "").lower()
    return row


def _payload_first(payload: Mapping[str, Any], *keys: str) -> Any:
    for key in keys:
        if key in payload and payload[key] is not None:
            return payload[key]
    return None


def _accepted_from_payload(payload: Mapping[str, Any]) -> bool | None:
    value = _payload_first(payload, "accepted", "edge_accepted", "qc_accepted")
    if value is not None:
        return _as_bool(value)
    decision = _payload_first(payload, "decision", "edge_decision")
    if decision is None:
        return None
    return str(decision).strip().lower() in {
        "accept",
        "ok",
        "pass",
        "warning",
    }


def _as_bool(value: Any) -> bool:
    if isinstance(value, str):
        return value.strip().lower() in {"1", "true", "yes", "y", "ok"}
    return bool(value)


def _float_or_nan(value: Any) -> float:
    try:
        out = float(value)
    except Exception:
        return float("nan")
    return out if math.isfinite(out) else float("nan")


def _latest_station_metrics(
    packets: Iterable[Mapping[str, Any]],
    session: FieldSession,
) -> dict[str, dict[str, Any]]:
    latest: dict[str, dict[str, Any]] = defaultdict(dict)
    by_device = {
        device_id: device.station
        for device_id, device in session.devices.items()
        if device.station
    }
    ordered = sorted(packets, key=lambda r: float(r.get("timestamp", 0.0)))
    for row in ordered:
        sid = row.get("station") or by_device.get(str(row.get("device_id")))
        if sid is None:
            continue
        bucket = latest[str(sid)]
        for key in (
            "battery_v",
            "clock_offset_ms",
            "runtime_days",
            "power_state",
            "sync_quality",
        ):
            value = row.get(key)
            if isinstance(value, float) and not math.isfinite(value):
                continue
            if value not in (None, ""):
                bucket[key] = value
    return dict(latest)


def _station_level(row: Mapping[str, Any]) -> str:
    acceptance = row.get("acceptance_rate")
    try:
        rate = float(acceptance)
    except Exception:
        rate = float("nan")
    power = str(row.get("power_state", "") or "").lower()
    sync = str(row.get("sync_quality", "") or "").lower()
    if power in {"critical"} or sync in {"poor"}:
        return "critical"
    if math.isfinite(rate):
        if rate < 0.85:
            return "critical"
        if rate < 0.95:
            return "warning"
    if power in {"warning"} or sync in {"fair"}:
        return "warning"
    if math.isfinite(rate) or power or sync:
        return "ok"
    return "unknown"


def _plot_station_panel(ax: Any, data: Mapping[str, Any], *, station_axis: str) -> None:
    stations = list(data.get("stations", []))
    ax.set_title("Station health")
    if not stations:
        _empty_panel(ax, "No stations")
        return
    use_map = _use_map_axis(stations, station_axis)
    xs, ys, labels, colors, sizes = [], [], [], [], []
    profiles = {row.get("profile") for row in stations if row.get("profile")}
    profile_index = {p: i for i, p in enumerate(sorted(profiles))}
    for idx, row in enumerate(stations):
        if use_map:
            lat, lon = _lat_lon(row.get("coords"))
            x = lon
            y = lat
        else:
            x = _float_or_nan(row.get("position_m"))
            if not math.isfinite(x):
                x = float(idx)
            profile = row.get("profile")
            y = float(profile_index.get(profile, 0))
        labels.append(str(row.get("station_id")))
        colors.append(_LEVEL_COLORS.get(row.get("health_level"), "#9e9e9e"))
        rate = _float_or_nan(row.get("acceptance_rate"))
        sizes.append(80.0 + 220.0 * (rate if math.isfinite(rate) else 0.25))
        xs.append(x)
        ys.append(y)
    ax.scatter(xs, ys, s=sizes, c=colors, edgecolor="black", linewidth=0.8)
    for x, y, label in zip(xs, ys, labels):
        ax.text(x, y, f" {label}", va="center", fontsize=8)
    ax.set_xlabel("Longitude" if use_map else "Profile position / station index")
    ax.set_ylabel("Latitude" if use_map else "Profile")
    if not use_map and profile_index:
        ax.set_yticks(list(profile_index.values()))
        ax.set_yticklabels(list(profile_index.keys()))
    ax.grid(True, alpha=0.25)


def _use_map_axis(stations: list[Mapping[str, Any]], station_axis: str) -> bool:
    mode = station_axis.lower()
    if mode not in {"auto", "profile", "map"}:
        raise ValueError("station_axis must be 'auto', 'profile', or 'map'.")
    if mode == "profile":
        return False
    has_coords = []
    for row in stations:
        coords = row.get("coords")
        lat, lon = _lat_lon(coords)
        has_coords.append(math.isfinite(lat) and math.isfinite(lon))
    return all(has_coords) if mode == "auto" else True


def _lat_lon(coords: Any) -> tuple[float, float]:
    try:
        seq = list(coords)
        return _float_or_nan(seq[0]), _float_or_nan(seq[1])
    except Exception:
        return float("nan"), float("nan")


def _plot_acceptance_panel(ax: Any, data: Mapping[str, Any]) -> None:
    stations = list(data.get("stations", []))
    ax.set_title("Edge QC acceptance")
    if not stations:
        _empty_panel(ax, "No QC data")
        return
    labels = [str(row.get("station_id")) for row in stations]
    values = [_float_or_nan(row.get("acceptance_rate")) for row in stations]
    plot_values = [v if math.isfinite(v) else 0.0 for v in values]
    colors = [_LEVEL_COLORS.get(row.get("health_level"), "#9e9e9e") for row in stations]
    ax.bar(labels, plot_values, color=colors, edgecolor="black", linewidth=0.5)
    ax.axhline(0.85, color="#de2d26", lw=1.0, ls="--", label="0.85")
    ax.axhline(0.95, color="#2ca25f", lw=1.0, ls=":", label="0.95")
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Acceptance rate")
    ax.tick_params(axis="x", rotation=45)
    ax.legend(frameon=False, loc="lower right")
    for idx, value in enumerate(values):
        text = "n/a" if not math.isfinite(value) else f"{value:.0%}"
        ax.text(idx, plot_values[idx] + 0.03, text, ha="center", fontsize=8)


def _plot_operations_panel(ax: Any, data: Mapping[str, Any]) -> None:
    stations = list(data.get("stations", []))
    ax.set_title("Power and synchronisation")
    if not stations:
        _empty_panel(ax, "No operations data")
        return
    labels = [str(row.get("station_id")) for row in stations]
    battery = [_float_or_nan(row.get("battery_v")) for row in stations]
    runtime = [_float_or_nan(row.get("runtime_days")) for row in stations]
    x = np.arange(len(labels), dtype=float)
    if any(math.isfinite(v) for v in battery):
        ax.plot(x, battery, marker="o", color="#3182bd", label="battery V")
        ax.set_ylabel("Battery voltage")
    elif any(math.isfinite(v) for v in runtime):
        ax.plot(x, runtime, marker="o", color="#756bb1", label="runtime d")
        ax.set_ylabel("Runtime days")
    else:
        _empty_panel(ax, "No battery/runtime packets")
        return
    for idx, row in enumerate(stations):
        sync = str(row.get("sync_quality", "") or "")
        power = str(row.get("power_state", "") or "")
        notes = " / ".join(v for v in (power, sync) if v)
        if notes:
            ax.text(
                idx,
                ax.get_ylim()[0],
                notes,
                rotation=90,
                va="bottom",
                ha="center",
                fontsize=7,
            )
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45)
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(frameon=False, loc="best")


def _plot_timeline_panel(ax: Any, data: Mapping[str, Any]) -> None:
    packets = list(data.get("packets", []))
    ax.set_title("Telemetry timeline")
    if not packets:
        _empty_panel(ax, "No packets")
        return
    kinds = sorted({str(row.get("kind")) for row in packets})
    kind_index = {kind: i for i, kind in enumerate(kinds)}
    times = np.asarray([_float_or_nan(row.get("timestamp")) for row in packets])
    finite_times = times[np.isfinite(times)]
    if finite_times.size and finite_times.max() > finite_times.min():
        xvals = (times - finite_times.min()) / 60.0
        xlabel = "Minutes since first packet"
    else:
        xvals = np.arange(len(packets), dtype=float)
        xlabel = "Packet index"
    yvals = [kind_index[str(row.get("kind"))] for row in packets]
    colors = [
        _LEVEL_COLORS["critical"]
        if row.get("accepted") is False
        else _LEVEL_COLORS["ok"]
        if row.get("accepted") is True
        else "#6baed6"
        for row in packets
    ]
    ax.scatter(xvals, yvals, c=colors, s=45, edgecolor="black", linewidth=0.4)
    ax.set_yticks(list(kind_index.values()))
    ax.set_yticklabels(kinds)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Packet kind")
    ax.grid(True, alpha=0.25)
    status = data.get("monitoring", {})
    level = status.get("level", "unknown")
    issues = data.get("issues", [])
    issue_text = ", ".join(issues[:3]) if issues else "no issues"
    ax.text(
        0.01,
        0.98,
        f"status: {level}\n{issue_text}",
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=8,
        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#cccccc"),
    )


def _empty_panel(ax: Any, text: str) -> None:
    ax.text(0.5, 0.5, text, ha="center", va="center", transform=ax.transAxes)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_alpha(0.2)
