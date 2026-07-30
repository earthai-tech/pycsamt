"""Per-method acquisition profiles for method-aware IoT QC.

:class:`~pycsamt.iot.monitoring.EMMethod` names the survey method, but on
its own it does not change any threshold: an operator still has to set the
expected frequency band, the required channels, and so on by hand. This
module carries the canonical, method-specific knowledge that turns
``EMMethod`` into concrete QC defaults.

A :class:`MethodProfile` records, per method, the typical acquisition band,
the channels a valid record must carry, a nominal sample rate, and whether
the method uses a controlled source or is sensitive to powerline noise.
:func:`method_profile` looks one up, :func:`target_bands_for_method` breaks
the band into decade sub-bands for coverage checks, and
:meth:`~pycsamt.iot.monitoring.MonitoringConfig.for_method` seeds a monitor
from a profile.

The band edges are representative defaults, not hard standards; override
them per survey when acquisition parameters differ.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

from .monitoring import EMMethod, _normalise_enum

__all__ = [
    "MethodProfile",
    "METHOD_PROFILES",
    "method_profile",
    "target_bands_for_method",
]


@dataclass(frozen=True)
class MethodProfile:
    """Canonical acquisition characteristics for one EM method."""

    method: EMMethod
    frequency_band_hz: tuple[float, float] | None
    required_channels: tuple[str, ...]
    default_sample_rate_hz: float | None = None
    controlled_source: bool = False
    powerline_sensitive: bool = False
    description: str = ""

    @property
    def is_natural_source(self) -> bool:
        """True for passive methods (no transmitter)."""
        return not self.controlled_source

    def as_dict(self) -> dict[str, Any]:
        return dict(
            method=self.method.value,
            frequency_band_hz=(
                list(self.frequency_band_hz)
                if self.frequency_band_hz
                else None
            ),
            required_channels=list(self.required_channels),
            default_sample_rate_hz=self.default_sample_rate_hz,
            controlled_source=self.controlled_source,
            powerline_sensitive=self.powerline_sensitive,
            description=self.description,
        )


# Representative acquisition profiles. Frequency edges follow common
# practice (e.g. Zonge/Phoenix instrument bands); tune per survey.
METHOD_PROFILES: dict[EMMethod, MethodProfile] = {
    EMMethod.AMT: MethodProfile(
        method=EMMethod.AMT,
        frequency_band_hz=(1.0, 10_000.0),
        required_channels=("ex", "ey", "hx", "hy"),
        default_sample_rate_hz=24_000.0,
        controlled_source=False,
        powerline_sensitive=True,
        description="Audio-frequency natural-source magnetotellurics.",
    ),
    EMMethod.MT: MethodProfile(
        method=EMMethod.MT,
        frequency_band_hz=(1.0e-4, 1_000.0),
        required_channels=("ex", "ey", "hx", "hy", "hz"),
        default_sample_rate_hz=1_000.0,
        controlled_source=False,
        powerline_sensitive=True,
        description="Broadband/long-period natural-source magnetotellurics.",
    ),
    EMMethod.CSAMT: MethodProfile(
        method=EMMethod.CSAMT,
        frequency_band_hz=(0.125, 8_192.0),
        required_channels=("ex", "ey", "hx", "hy"),
        default_sample_rate_hz=32_768.0,
        controlled_source=True,
        powerline_sensitive=True,
        description="Controlled-source audio-frequency magnetotellurics.",
    ),
    EMMethod.CSEM: MethodProfile(
        method=EMMethod.CSEM,
        frequency_band_hz=(0.01, 100.0),
        required_channels=("ex", "ey", "hx", "hy"),
        default_sample_rate_hz=1_000.0,
        controlled_source=True,
        powerline_sensitive=True,
        description=(
            "Controlled-source electromagnetics: a dipole source and a "
            "receiver array recording amplitude/phase versus offset."
        ),
    ),
    EMMethod.TDEM: MethodProfile(
        method=EMMethod.TDEM,
        frequency_band_hz=None,  # time-domain: gated transients
        required_channels=("hz",),
        default_sample_rate_hz=100_000.0,
        controlled_source=True,
        powerline_sensitive=False,
        description="Time-domain (transient) electromagnetics.",
    ),
    EMMethod.TEM: MethodProfile(
        method=EMMethod.TEM,
        frequency_band_hz=None,
        required_channels=("hz",),
        default_sample_rate_hz=100_000.0,
        controlled_source=True,
        powerline_sensitive=False,
        description="Transient electromagnetics (alias of TDEM).",
    ),
    EMMethod.UNKNOWN: MethodProfile(
        method=EMMethod.UNKNOWN,
        frequency_band_hz=None,
        required_channels=(),
        default_sample_rate_hz=None,
        controlled_source=False,
        powerline_sensitive=False,
        description="Unspecified method; no method-derived defaults.",
    ),
}


def method_profile(method: EMMethod | str) -> MethodProfile:
    """Return the :class:`MethodProfile` for *method*.

    Parameters
    ----------
    method : EMMethod or str
        Method enum member or a name/value such as ``"csamt"``.

    Returns
    -------
    MethodProfile
        The registered profile, or the ``UNKNOWN`` profile for a method
        with no registered defaults.
    """
    key = _normalise_enum(method, EMMethod, "method")
    return METHOD_PROFILES.get(key, METHOD_PROFILES[EMMethod.UNKNOWN])


def target_bands_for_method(
    method: EMMethod | str,
) -> list[tuple[float, float]]:
    """Return per-decade sub-bands spanning a method's frequency band.

    These feed coverage checks such as
    :func:`~pycsamt.iot.edge_amt.estimate_frequency_coverage`, which report
    the fraction of target bands a recording actually resolves. A method
    without a defined band (e.g. time-domain) returns an empty list.
    """
    profile = method_profile(method)
    band = profile.frequency_band_hz
    if band is None:
        return []
    lo, hi = band
    lo_exp = math.floor(math.log10(lo))
    hi_exp = math.ceil(math.log10(hi))
    edges = [10.0**e for e in range(lo_exp, hi_exp + 1)]
    # Clamp the outer edges to the actual band limits.
    edges[0] = lo
    edges[-1] = hi
    return [
        (edges[i], edges[i + 1])
        for i in range(len(edges) - 1)
        if edges[i] < edges[i + 1]
    ]
