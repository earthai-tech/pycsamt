from __future__ import annotations

import math
from collections.abc import Iterable
from dataclasses import dataclass, field
from enum import Enum
from typing import Any

import pandas as pd

from ..api.property import MetadataMixin, PyCSAMTObject
from ..api.view import maybe_wrap_frame
from .core import DeviceConfig, PacketKind, TelemetryPacket


class PowerState(str, Enum):
    """Energy status for an IoT field device."""

    SUSTAINING = "sustaining"
    OK = "ok"
    WARNING = "warning"
    CRITICAL = "critical"


def _as_positive(value: Any, name: str) -> float:
    out = float(value)
    if out <= 0:
        raise ValueError(f"{name} must be positive.")
    return out


def _as_nonnegative(value: Any, name: str) -> float:
    out = float(value)
    if out < 0:
        raise ValueError(f"{name} must be >= 0.")
    return out


def _as_probability(value: Any, name: str) -> float:
    out = float(value)
    if not 0.0 <= out <= 1.0:
        raise ValueError(f"{name} must be between 0 and 1.")
    return out


def _as_optional_positive(value: Any, name: str) -> float | None:
    if value is None:
        return None
    return _as_positive(value, name)


@dataclass
class EnergyConfig(PyCSAMTObject, MetadataMixin):
    """Power-budget inputs for a field IoT device."""

    battery_wh: float
    active_power_w: float
    sleep_power_w: float = 0.05
    duty_cycle: float = 1.0
    solar_wh_per_day: float = 0.0
    reserve_fraction: float = 0.0
    regulator_efficiency: float = 1.0
    charge_efficiency: float = 1.0
    telemetry_power_w: float = 0.0
    telemetry_seconds_per_day: float = 0.0
    edge_power_w: float = 0.0
    edge_duty_cycle: float = 0.0
    auxiliary_wh_per_day: float = 0.0
    min_runtime_days: float | None = None
    device_id: str | None = None
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        """Validate and normalise power-budget inputs."""
        self.battery_wh = _as_positive(self.battery_wh, "battery_wh")
        self.active_power_w = _as_nonnegative(
            self.active_power_w,
            "active_power_w",
        )
        self.sleep_power_w = _as_nonnegative(
            self.sleep_power_w,
            "sleep_power_w",
        )
        self.duty_cycle = _as_probability(self.duty_cycle, "duty_cycle")
        self.solar_wh_per_day = _as_nonnegative(
            self.solar_wh_per_day,
            "solar_wh_per_day",
        )
        self.reserve_fraction = _as_probability(
            self.reserve_fraction,
            "reserve_fraction",
        )
        self.regulator_efficiency = _as_probability(
            self.regulator_efficiency,
            "regulator_efficiency",
        )
        self.charge_efficiency = _as_probability(
            self.charge_efficiency,
            "charge_efficiency",
        )
        if self.regulator_efficiency <= 0:
            raise ValueError("regulator_efficiency must be > 0.")
        if self.charge_efficiency <= 0:
            raise ValueError("charge_efficiency must be > 0.")
        self.telemetry_power_w = _as_nonnegative(
            self.telemetry_power_w,
            "telemetry_power_w",
        )
        self.telemetry_seconds_per_day = _as_nonnegative(
            self.telemetry_seconds_per_day,
            "telemetry_seconds_per_day",
        )
        self.edge_power_w = _as_nonnegative(self.edge_power_w, "edge_power_w")
        self.edge_duty_cycle = _as_probability(
            self.edge_duty_cycle,
            "edge_duty_cycle",
        )
        self.auxiliary_wh_per_day = _as_nonnegative(
            self.auxiliary_wh_per_day,
            "auxiliary_wh_per_day",
        )
        self.min_runtime_days = _as_optional_positive(
            self.min_runtime_days,
            "min_runtime_days",
        )
        if self.device_id is not None:
            self.device_id = str(self.device_id).strip() or None
        if not isinstance(self.metadata, dict):
            self.metadata = dict(self.metadata or {})

    @property
    def usable_battery_wh(self) -> float:
        """Battery energy available after reserve is held back."""
        return self.battery_wh * (1.0 - self.reserve_fraction)

    @property
    def base_average_power_w(self) -> float:
        """Average active/sleep power before auxiliary loads."""
        return (
            self.duty_cycle * self.active_power_w
            + (1.0 - self.duty_cycle) * self.sleep_power_w
        )

    @property
    def telemetry_wh_per_day(self) -> float:
        """Daily energy used by radio transmission windows."""
        return self.telemetry_power_w * self.telemetry_seconds_per_day / 3600.0

    @property
    def edge_wh_per_day(self) -> float:
        """Daily energy used by edge processing overhead."""
        return 24.0 * self.edge_power_w * self.edge_duty_cycle

    @property
    def harvest_wh_per_day(self) -> float:
        """Daily usable harvested energy after charge efficiency."""
        return self.solar_wh_per_day * self.charge_efficiency


@dataclass
class DevicePowerProfile(PyCSAMTObject):
    """Named power profile for a field node or recorder."""

    name: str
    active_power_w: float
    sleep_power_w: float = 0.05
    telemetry_power_w: float = 0.0
    edge_power_w: float = 0.0
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        """Validate and normalise profile fields."""
        self.name = str(self.name).strip()
        if not self.name:
            raise ValueError("name cannot be empty.")
        self.active_power_w = _as_nonnegative(
            self.active_power_w,
            "active_power_w",
        )
        self.sleep_power_w = _as_nonnegative(
            self.sleep_power_w,
            "sleep_power_w",
        )
        self.telemetry_power_w = _as_nonnegative(
            self.telemetry_power_w,
            "telemetry_power_w",
        )
        self.edge_power_w = _as_nonnegative(self.edge_power_w, "edge_power_w")
        if not isinstance(self.metadata, dict):
            self.metadata = dict(self.metadata or {})

    def apply(
        self,
        *,
        battery_wh: float,
        duty_cycle: float = 1.0,
        solar_wh_per_day: float = 0.0,
        **kwargs: Any,
    ) -> EnergyConfig:
        """Build an :class:`EnergyConfig` from this profile."""
        return EnergyConfig(
            battery_wh=battery_wh,
            active_power_w=self.active_power_w,
            sleep_power_w=self.sleep_power_w,
            duty_cycle=duty_cycle,
            solar_wh_per_day=solar_wh_per_day,
            telemetry_power_w=self.telemetry_power_w,
            edge_power_w=self.edge_power_w,
            **kwargs,
        )


@dataclass
class EnergyEstimate(PyCSAMTObject):
    """Estimated runtime and power draw."""

    average_power_w: float
    runtime_hours: float
    runtime_days: float
    net_wh_per_day: float
    load_wh_per_day: float = 0.0
    harvest_wh_per_day: float = 0.0
    usable_battery_wh: float = 0.0
    telemetry_wh_per_day: float = 0.0
    edge_wh_per_day: float = 0.0
    auxiliary_wh_per_day: float = 0.0
    reserve_wh: float = 0.0
    energy_margin_wh_per_day: float = 0.0
    autonomy_days_no_harvest: float = 0.0
    state: PowerState | str = PowerState.OK
    issues: list[str] = field(default_factory=list)

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        """Validate and normalise estimate fields."""
        self.average_power_w = float(self.average_power_w)
        self.runtime_hours = float(self.runtime_hours)
        self.runtime_days = float(self.runtime_days)
        self.net_wh_per_day = float(self.net_wh_per_day)
        self.load_wh_per_day = float(self.load_wh_per_day)
        self.harvest_wh_per_day = float(self.harvest_wh_per_day)
        self.usable_battery_wh = float(self.usable_battery_wh)
        self.telemetry_wh_per_day = float(self.telemetry_wh_per_day)
        self.edge_wh_per_day = float(self.edge_wh_per_day)
        self.auxiliary_wh_per_day = float(self.auxiliary_wh_per_day)
        self.reserve_wh = float(self.reserve_wh)
        self.energy_margin_wh_per_day = float(self.energy_margin_wh_per_day)
        self.autonomy_days_no_harvest = float(self.autonomy_days_no_harvest)
        self.state = (
            self.state
            if isinstance(self.state, PowerState)
            else (PowerState(str(self.state)))
        )
        self.issues = [str(issue) for issue in list(self.issues or [])]

    @property
    def sustaining(self) -> bool:
        """Return whether harvested energy covers daily load."""
        return self.state == PowerState.SUSTAINING

    def as_dict(self) -> dict[str, Any]:
        """Return a serialisable estimate dictionary."""
        return dict(
            average_power_w=self.average_power_w,
            runtime_hours=self.runtime_hours,
            runtime_days=self.runtime_days,
            net_wh_per_day=self.net_wh_per_day,
            load_wh_per_day=self.load_wh_per_day,
            harvest_wh_per_day=self.harvest_wh_per_day,
            usable_battery_wh=self.usable_battery_wh,
            telemetry_wh_per_day=self.telemetry_wh_per_day,
            edge_wh_per_day=self.edge_wh_per_day,
            auxiliary_wh_per_day=self.auxiliary_wh_per_day,
            reserve_wh=self.reserve_wh,
            energy_margin_wh_per_day=self.energy_margin_wh_per_day,
            autonomy_days_no_harvest=self.autonomy_days_no_harvest,
            state=self.state.value
            if isinstance(self.state, PowerState)
            else str(self.state),
            issues=";".join(self.issues),
        )

    def to_packet(
        self,
        device: DeviceConfig,
        *,
        timestamp: float,
        survey_id: str | None = None,
        qos: int = 0,
        retained: bool = False,
    ) -> TelemetryPacket:
        """Encode this estimate as a power telemetry packet."""
        return TelemetryPacket.from_device(
            device,
            timestamp=timestamp,
            payload=self.as_dict(),
            kind=PacketKind.POWER,
            survey_id=survey_id,
            qos=qos,
            retained=retained,
        )


def estimate_energy_budget(config: EnergyConfig) -> EnergyEstimate:
    """Estimate runtime from battery, duty cycle, and optional solar input."""
    config.validate()
    avg_base = config.base_average_power_w
    load_daily = (
        24.0 * avg_base / config.regulator_efficiency
        + config.telemetry_wh_per_day
        + config.edge_wh_per_day
        + config.auxiliary_wh_per_day
    )
    avg = load_daily / 24.0
    avg = max(avg, 1e-12)
    net_daily = load_daily - config.harvest_wh_per_day
    if net_daily <= 0:
        runtime_hours = float("inf")
    else:
        runtime_hours = 24.0 * config.usable_battery_wh / net_daily
    no_harvest_days = (
        config.usable_battery_wh / load_daily if load_daily > 0 else float("inf")
    )
    issues = _estimate_issues(config, runtime_hours, net_daily)
    state = _estimate_state(config, runtime_hours, net_daily, issues)
    return EnergyEstimate(
        average_power_w=avg,
        runtime_hours=runtime_hours,
        runtime_days=runtime_hours / 24.0,
        net_wh_per_day=net_daily,
        load_wh_per_day=load_daily,
        harvest_wh_per_day=config.harvest_wh_per_day,
        usable_battery_wh=config.usable_battery_wh,
        telemetry_wh_per_day=config.telemetry_wh_per_day,
        edge_wh_per_day=config.edge_wh_per_day,
        auxiliary_wh_per_day=config.auxiliary_wh_per_day,
        reserve_wh=config.battery_wh - config.usable_battery_wh,
        energy_margin_wh_per_day=-net_daily,
        autonomy_days_no_harvest=no_harvest_days,
        state=state,
        issues=issues,
    )


def _estimate_issues(
    config: EnergyConfig,
    runtime_hours: float,
    net_wh_per_day: float,
) -> list[str]:
    issues: list[str] = []
    if net_wh_per_day > 0:
        issues.append("daily_energy_deficit")
    if (
        config.min_runtime_days is not None
        and not math.isinf(runtime_hours)
        and runtime_hours / 24.0 < config.min_runtime_days
    ):
        issues.append("runtime_below_minimum")
    if config.usable_battery_wh <= 0:
        issues.append("no_usable_battery_after_reserve")
    return issues


def _estimate_state(
    config: EnergyConfig,
    runtime_hours: float,
    net_wh_per_day: float,
    issues: list[str],
) -> PowerState:
    if net_wh_per_day <= 0:
        return PowerState.SUSTAINING
    if "no_usable_battery_after_reserve" in issues:
        return PowerState.CRITICAL
    if config.min_runtime_days is None:
        return PowerState.OK
    runtime_days = runtime_hours / 24.0
    if runtime_days >= config.min_runtime_days:
        return PowerState.OK
    if runtime_days >= 0.5 * config.min_runtime_days:
        return PowerState.WARNING
    return PowerState.CRITICAL


def power_summary_table(
    estimates: EnergyEstimate | Iterable[EnergyEstimate],
    *,
    device_ids: Iterable[str] | None = None,
    api: bool | None = None,
) -> Any:
    """Return energy estimates as a pyCSAMT table."""
    items = [estimates] if isinstance(estimates, EnergyEstimate) else (list(estimates))
    ids = list(device_ids or [])
    rows = []
    for idx, estimate in enumerate(items):
        estimate.validate()
        row = estimate.as_dict()
        if idx < len(ids):
            row["device_id"] = ids[idx]
        rows.append(row)
    df = pd.DataFrame.from_records(rows)
    return maybe_wrap_frame(
        df,
        api=api,
        name="iot_power_summary",
        kind="iot.power",
        source=items,
        description="IoT energy budget and runtime estimates.",
    )


def estimate_deployment_energy(
    configs: Iterable[EnergyConfig],
    *,
    api: bool | None = None,
) -> Any:
    """Estimate and tabulate energy budgets for multiple devices."""
    rows = []
    for config in configs:
        estimate = estimate_energy_budget(config)
        row = estimate.as_dict()
        row["device_id"] = config.device_id
        rows.append(row)
    df = pd.DataFrame.from_records(rows)
    return maybe_wrap_frame(
        df,
        api=api,
        name="iot_deployment_energy",
        kind="iot.power.deployment",
        source=configs,
        description="Deployment-scale IoT energy budget estimates.",
    )


__all__ = [
    "PowerState",
    "DevicePowerProfile",
    "EnergyConfig",
    "EnergyEstimate",
    "estimate_deployment_energy",
    "estimate_energy_budget",
    "power_summary_table",
]
