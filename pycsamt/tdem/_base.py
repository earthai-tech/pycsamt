# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""Core data model for a single TEM sounding station."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from pycsamt.api.property import PyCSAMTObject

__all__ = ["TEMSounding"]


@dataclass(repr=False)
class TEMSounding(PyCSAMTObject):
    r"""
    Container for a single time-domain EM (TEM) sounding.

    Stores the measured transient decay data together with all
    geometric and acquisition parameters needed to compute
    apparent resistivity and to synthesise a frequency-domain
    EDI file.

    Parameters
    ----------
    time_gates : array-like of float
        Centre times of the off-time measurement windows in
        seconds.  Shape ``(n_gates,)``.
    data : array-like of float
        Measured decay values at each time gate.  Units depend on
        ``data_type`` (see below).  Shape ``(n_gates,)``.
    current : float
        Transmitter current in Amperes.
    tx_area : float
        Transmitter loop area in m².  For a square loop of side
        *L* pass ``L**2``; for a circular loop of radius *r* pass
        ``np.pi * r**2``.
    data_type : str, optional
        Describes ``data`` units.

        ``"dBdt"``
            :math:`\partial B_z / \partial t` normalised by the
            receiver coil parameters
            (:math:`V / (n_{rx} A_{rx})`), in T s⁻¹.
        ``"dHdt"``
            :math:`\partial H_z / \partial t` in A m⁻¹ s⁻¹.
        ``"voltage"``
            Raw induced voltage in V.  Must also supply ``rx_area``
            and ``rx_turns`` so the converter can normalise.
        ``"normalized_voltage"``
            Voltage normalised by transmitter moment
            (V / (A·m²)).

        Default is ``"dBdt"``.
    tx_turns : int, optional
        Number of transmitter loop turns.  Default 1.
    rx_area : float, optional
        Receiver coil effective area in m².  Required when
        ``data_type="voltage"``.  Default 1.0.
    rx_turns : int, optional
        Number of receiver coil turns.  Required when
        ``data_type="voltage"``.  Default 1.
    offset : float, optional
        Horizontal transmitter–receiver separation in metres.
        ``0.0`` means coincident / central-loop configuration
        (the most common field geometry).  Default 0.0.
    loop_shape : str, optional
        Transmitter loop geometry: ``"square"``, ``"circle"``, or
        ``"rectangular"``.  Affects the geometric correction
        factor for *non*-coincident configurations.  Default
        ``"square"``.
    loop_dims : tuple of float, optional
        Loop dimensions:

        * ``"square"``   → ``(side_length,)``
        * ``"circle"``   → ``(radius,)``
        * ``"rectangular"`` → ``(length, width)``

    station_name : str, optional
        Site identifier used in the output EDI ``>HEAD`` block.
    x : float, optional
        Easting / longitude of the station.  Default 0.0.
    y : float, optional
        Northing / latitude of the station.  Default 0.0.
    elevation : float, optional
        Station elevation in metres.  Default 0.0.
    error : array-like of float, optional
        Absolute uncertainty on each ``data`` value, same units
        as ``data``.  When provided, it propagates to the EDI
        impedance variance arrays.
    waveform : object, optional
        Transmitter waveform descriptor (a
        :class:`~pycsamt.tdem.waveform.SquareWaveform` or
        :class:`~pycsamt.tdem.waveform.RampWaveform` instance).
        Required only for the rigorous Fourier-domain transform.

    Attributes
    ----------
    moment : float
        Transmitter magnetic moment
        :math:`M = I \cdot n_{tx} \cdot A_{tx}` in A·m².

    Notes
    -----
    For a **coincident-loop** system (offset = 0) the late-time
    apparent resistivity at gate *k* is:

    .. math::

        \rho_a(t_k) = \frac{\mu_0}{\pi}
            \left(\frac{\mu_0^2 M}
                  {20\,|\partial B_z/\partial t(t_k)|\,t_k^{5/2}}
            \right)^{2/3}

    (Ward & Hohmann 1988, eq. 4.96; Nabighian & Macnae 1991).

    The equivalent pseudo-frequency is mapped as:

    .. math::

        f_\text{eq}(t) = \frac{1}{2\pi t}

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.tdem import TEMSounding
    >>> t = np.logspace(-5, -2, 30)        # 10 µs … 10 ms
    >>> M = 8.0 * 100.**2                  # magnetic moment
    >>> import math
    >>> MU0 = 4 * math.pi * 1e-7
    >>> dBdt = M * MU0**2.5 / (10 * math.sqrt(math.pi) * 100.**1.5 * t**2.5)
    >>> s = TEMSounding(t, dBdt, current=8.0, tx_area=100.**2)
    >>> s.moment
    640000.0
    """

    time_gates: np.ndarray
    data: np.ndarray
    current: float
    tx_area: float

    data_type: str = "dBdt"
    tx_turns: int = 1
    rx_area: float = 1.0
    rx_turns: int = 1
    offset: float = 0.0
    loop_shape: str = "square"
    loop_dims: tuple[float, ...] = field(default_factory=lambda: (100.0,))

    station_name: str = ""
    x: float = 0.0
    y: float = 0.0
    elevation: float = 0.0

    error: np.ndarray | None = None
    waveform: object | None = None
    rx_position: tuple[float, float] | None = None
    """2-D receiver offset ``(rx, ry)`` in metres from the loop centre.
    When set, it overrides the scalar ``offset`` for the Biot-Savart
    in-loop correction.  ``None`` (default) falls back to ``(offset, 0)``."""
    rx_offset: tuple[float, float] | None = None
    """Alias for ``rx_position`` kept for reader and transform callers."""
    verbose: int = 0
    logger: object | None = None

    def __post_init__(self) -> None:
        self.time_gates = np.asarray(self.time_gates, dtype=float)
        self.data = np.asarray(self.data, dtype=float)
        if self.error is not None:
            self.error = np.asarray(self.error, dtype=float)
        if self.rx_position is None and self.rx_offset is not None:
            self.rx_position = self.rx_offset
        if self.rx_position is not None and self.offset <= 0.0:
            rx, ry = self.rx_position
            self.offset = float(np.hypot(rx, ry))

        if self.time_gates.ndim != 1:
            raise ValueError("time_gates must be 1-D")
        if self.data.shape != self.time_gates.shape:
            raise ValueError(
                f"data shape {self.data.shape} must match "
                f"time_gates shape {self.time_gates.shape}"
            )
        if (
            self.error is not None
            and self.error.shape != self.time_gates.shape
        ):
            raise ValueError("error shape must match time_gates shape")

        valid_dtypes = {"dBdt", "dHdt", "voltage", "normalized_voltage"}
        if self.data_type not in valid_dtypes:
            raise ValueError(
                f"data_type must be one of {valid_dtypes}, "
                f"got '{self.data_type}'"
            )

    @property
    def moment(self) -> float:
        r"""Transmitter magnetic moment :math:`M = I n_{tx} A_{tx}` (A·m²)."""
        return float(self.current) * int(self.tx_turns) * float(self.tx_area)

    @property
    def n_gates(self) -> int:
        """Number of time gates."""
        return int(self.time_gates.size)

    def dBdt(self) -> np.ndarray:
        r"""
        Return :math:`\partial B_z / \partial t` in T s⁻¹.

        Converts from ``data_type`` to normalised dB/dt.
        Raises ``ValueError`` for ``"normalized_voltage"`` (needs
        further physical context).
        """
        if self.data_type == "dBdt":
            return self.data.copy()

        MU0 = 4.0 * np.pi * 1e-7
        if self.data_type == "dHdt":
            return self.data * MU0
        if self.data_type == "voltage":
            nA = float(self.rx_turns) * float(self.rx_area)
            if nA == 0.0:
                raise ValueError(
                    "rx_area and rx_turns must be non-zero for "
                    "voltage → dBdt conversion"
                )
            return self.data / nA
        raise ValueError(
            "Cannot convert 'normalized_voltage' to dBdt automatically; "
            "please supply rx_area and rx_turns and use data_type='voltage'."
        )

    @classmethod
    def from_arrays(
        cls,
        time_gates,
        data,
        *,
        current: float,
        loop_side: float | None = None,
        loop_radius: float | None = None,
        tx_area: float | None = None,
        **kwargs,
    ) -> TEMSounding:
        """
        Convenience constructor — supply either ``loop_side``,
        ``loop_radius``, or ``tx_area`` directly.
        """
        if tx_area is None:
            if loop_side is not None:
                tx_area = loop_side ** 2
                kwargs.setdefault("loop_shape", "square")
                kwargs.setdefault("loop_dims", (float(loop_side),))
            elif loop_radius is not None:
                tx_area = np.pi * loop_radius ** 2
                kwargs.setdefault("loop_shape", "circle")
                kwargs.setdefault("loop_dims", (float(loop_radius),))
            else:
                raise ValueError(
                    "Supply one of loop_side, loop_radius, or tx_area"
                )
        return cls(
            time_gates=np.asarray(time_gates, float),
            data=np.asarray(data, float),
            current=float(current),
            tx_area=float(tx_area),
            **kwargs,
        )

    def summary(self) -> dict:
        """Return a compact dict of key acquisition parameters."""
        return {
            "station":     self.station_name or "(unnamed)",
            "n_gates":     self.n_gates,
            "t_min_ms":    float(self.time_gates.min()) * 1e3,
            "t_max_ms":    float(self.time_gates.max()) * 1e3,
            "moment":      self.moment,
            "data_type":   self.data_type,
            "loop_shape":  self.loop_shape,
            "offset_m":    self.offset,
            "has_error":   self.error is not None,
        }

    def __repr__(self) -> str:
        s = self.station_name or "(unnamed)"
        return (
            f"TEMSounding(station={s!r}, n_gates={self.n_gates}, "
            f"t=[{self.time_gates.min():.2e}, {self.time_gates.max():.2e}] s, "
            f"moment={self.moment:.3g} A·m²)"
        )

    def __str__(self) -> str:
        """Return the compact sounding representation."""
        return self.__repr__()
