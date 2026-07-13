# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""High-level survey data-file builder for MARE2DEM.

Port of ``m2d_makeDataFile.m``.

:func:`make_data_file` creates a MARE2DEM ``.emdata`` file for forward
modeling or inversion from:

* survey topography (flat scalar or ``[y, z]`` array);
* MT receiver positions and frequency list;
* CSEM receiver positions, transmitter positions, and frequency list.

The function handles the logic for placing receivers/transmitters
relative to the topographic surface (``'land'``, ``'marine'``,
``'amphibious'``), computing slope angles, and building the DATA block
that maps each datum to transmitter, receiver, and frequency indices.

Receiver placement rules (matching the MATLAB code exactly):

* ``'land'``       — receiver is 0.1 m *below* the topographic surface;
                     electric dipole tilted along the slope; magnetic
                     dipole horizontal (beta = 0).
* ``'marine'``     — receiver is 0.1 m *above* the topographic surface;
                     dipole aligned slope-parallel.
* ``'amphibious'`` — ``'land'`` where z > 0, ``'marine'`` where z ≤ 0.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from .geom.topo import parse_topo
from .iotools.emdata import (
    CSEMConfig,
    EMDataFile,
    MTConfig,
    UTMOrigin,
    write_emdata,
)

__all__ = [
    "MTSurveyConfig",
    "CSEMSurveyConfig",
    "make_data_file",
]


# ---------------------------------------------------------------------------
# Survey configuration dataclasses
# ---------------------------------------------------------------------------


@dataclass
class MTSurveyConfig:
    """MT receiver and data-selection configuration.

    Parameters
    ----------
    frequencies : array-like
        MT frequencies in Hz.
    rx_y : array-like
        Receiver y positions along the 2-D profile.
    rx_type : str
        Receiver placement mode: ``'land'``, ``'marine'``,
        or ``'amphibious'``.
    rx_z : array-like or None
        Override receiver depth positions. When supplied, slope-
        based tilt angles are not computed.
    rx_z_offset : float or array-like or None
        Offset from topography (< 0 above, > 0 below). Overrides
        the ``rx_type`` depth rule.
    rx_beta : float or array-like or None
        Override receiver beta (y-tilt) angles in degrees.
    rx_name : list of str or None
        Station names.
    lTE : bool
        Include log10 apparent resistivity and phase (TE mode).
    lTM : bool
        Include log10 apparent resistivity and phase (TM mode).
    lZDet : bool
        Include impedance-determinant apparent resistivity and phase.
    lTipper : bool
        Include TE tipper amplitude and phase.
    lTipperRealImag : bool
        Include TE tipper real and imaginary parts.
    lMTFields : bool
        Include TE and TM field vectors (Ex, Ey, Ez, Hx, Hy, Hz).
    """

    frequencies: np.ndarray = field(
        default_factory=lambda: np.empty(0, dtype=float)
    )
    rx_y: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=float))
    rx_type: str = "marine"
    rx_z: np.ndarray | None = None
    rx_z_offset: float | np.ndarray | None = None
    rx_beta: float | np.ndarray | None = None
    rx_name: list[str] | None = None

    lTE: bool = True
    lTM: bool = True
    lZDet: bool = False
    lTipper: bool = False
    lTipperRealImag: bool = False
    lMTFields: bool = False


@dataclass
class CSEMSurveyConfig:
    """CSEM source–receiver and data-selection configuration.

    Parameters
    ----------
    frequencies : array-like
        CSEM frequencies in Hz.
    tx_y : array-like
        Transmitter y positions.
    rx_y : array-like or None
        Receiver y positions. Mutually exclusive with *rx_r*.
    rx_r : array-like or None
        Towed receiver offsets from transmitter positions.
    tx_z : array-like or None
        Override transmitter depth positions.
    tx_z_offset : float or array-like or None
        Transmitter offset relative to topography.
    rx_z : array-like or None
        Override receiver depth positions.
    rx_z_offset : float or array-like or None
        Receiver offset relative to topography.
    rx_type : str
        Receiver placement mode.
    tx_type : str
        Dipole type: ``'edipole'`` or ``'bdipole'``.
    tx_azimuth : float or array-like or None
        Transmitter azimuth (degrees from x towards y).
    tx_dip : float or array-like or None
        Transmitter dip (degrees positive down).
    tx_length : float or array-like or None
        Electric dipole length in metres.
    rx_length : float or array-like or None
        Receiver dipole length in metres.
    rx_beta : float or array-like or None
        Receiver beta (y-tilt) angles in degrees.
    tx_name : list of str or None
    rx_name : list of str or None
    phase_convention : str
        ``'lag'`` (default) or ``'lead'``.
    min_range : float or None
        Minimum Tx–Rx range to include in data file.
    max_range : float or None
        Maximum Tx–Rx range.
    lEx : bool
        Include Ex log10 amplitude and phase.
    lEy : bool
    lEz : bool
    lBx : bool
    lBy : bool
    lBz : bool
    """

    frequencies: np.ndarray = field(
        default_factory=lambda: np.empty(0, dtype=float)
    )
    tx_y: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=float))
    rx_y: np.ndarray | None = None
    rx_r: np.ndarray | None = None
    tx_z: np.ndarray | None = None
    tx_z_offset: float | np.ndarray | None = None
    rx_z: np.ndarray | None = None
    rx_z_offset: float | np.ndarray | None = None
    rx_type: str = "marine"
    tx_type: str = "edipole"
    tx_azimuth: float | np.ndarray | None = None
    tx_dip: float | np.ndarray | None = None
    tx_length: float | np.ndarray | None = None
    rx_length: float | np.ndarray | None = None
    rx_beta: float | np.ndarray | None = None
    tx_name: list[str] | None = None
    rx_name: list[str] | None = None
    phase_convention: str = "lag"
    min_range: float | None = None
    max_range: float | None = None

    lEx: bool = True
    lEy: bool = False
    lEz: bool = False
    lBx: bool = False
    lBy: bool = False
    lBz: bool = False


# ---------------------------------------------------------------------------
# Internal geometry helpers
# ---------------------------------------------------------------------------


def _broadcast(val: Any, n: int) -> np.ndarray:
    """Broadcast a scalar or 1-D array to length *n*."""
    arr = np.asarray(val, dtype=float).ravel()
    if len(arr) == 1:
        return np.full(n, arr[0])
    return arr


def _place_receivers(
    topo: Any,
    y_rx: np.ndarray,
    rx_type: str,
    *,
    rx_z: np.ndarray | None = None,
    rx_z_offset: Any = None,
    rx_beta: Any = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return (x, z, theta=0, alpha=0, beta) for receiver positions.

    Returns
    -------
    x_rx, z_rx, theta_rx (zeros), alpha_rx (zeros), beta_rx
    """
    n = len(y_rx)
    z_topo, slope, _ = parse_topo(topo, y_rx)

    if rx_z is not None:
        z_rx_arr = _broadcast(rx_z, n)
        beta_arr = np.zeros(n)
    else:
        if rx_z_offset is not None:
            dz = _broadcast(rx_z_offset, n)
            z_rx_arr = z_topo + dz
        else:
            dz_map = {"land": 0.1, "marine": -0.1}
            if rx_type in dz_map:
                z_rx_arr = z_topo + dz_map[rx_type]
            else:  # amphibious
                dz_arr = np.where(z_topo > 0, -0.1, 0.1)
                z_rx_arr = z_topo + dz_arr
        beta_arr = slope.copy()

    if rx_beta is not None:
        beta_arr = _broadcast(rx_beta, n)

    x_rx = np.zeros(n)
    theta_rx = np.zeros(n)
    alpha_rx = np.zeros(n)
    return x_rx, z_rx_arr, theta_rx, alpha_rx, beta_arr


def _place_transmitters(
    topo: Any,
    y_tx: np.ndarray,
    *,
    tx_z: np.ndarray | None = None,
    tx_z_offset: Any = None,
    tx_azimuth: Any = None,
    tx_dip: Any = None,
    tx_length: Any = None,
    tx_type: str = "edipole",
) -> tuple[np.ndarray, list[str]]:
    """Return transmitters array (n, 6) and type list."""
    n = len(y_tx)
    z_topo, slope, _ = parse_topo(topo, y_tx)

    if tx_z is not None:
        z_tx = _broadcast(tx_z, n)
    elif tx_z_offset is not None:
        z_tx = z_topo + _broadcast(tx_z_offset, n)
    else:
        z_tx = z_topo

    x_tx = np.zeros(n)

    if tx_azimuth is not None:
        azimuth = _broadcast(tx_azimuth, n)
    else:
        azimuth = np.full(n, 90.0)  # inline default

    if tx_dip is not None:
        dip = _broadcast(tx_dip, n)
    else:
        dip = np.zeros(n)

    if tx_length is not None:
        length = _broadcast(tx_length, n)
    else:
        length = np.zeros(n)

    txs = np.column_stack([x_tx, y_tx, z_tx, azimuth, dip, length])
    types = [tx_type] * n
    return txs, types


# ---------------------------------------------------------------------------
# MT data builder
# ---------------------------------------------------------------------------


def _build_mt_config(
    topo: Any,
    mt: MTSurveyConfig,
) -> tuple[MTConfig, np.ndarray]:
    """Build MTConfig and DATA rows for MT forward modeling."""
    y_rx = np.asarray(mt.rx_y, dtype=float).ravel()
    n_rx = len(y_rx)
    x_rx, z_rx, theta, alpha, beta = _place_receivers(
        topo,
        y_rx,
        mt.rx_type,
        rx_z=mt.rx_z,
        rx_z_offset=mt.rx_z_offset,
        rx_beta=mt.rx_beta,
    )
    freqs = np.asarray(mt.frequencies, dtype=float)

    # Land tilted Rx: add extra H-field receivers with beta=0
    irx_mag = np.arange(n_rx)
    extra_rx: list[np.ndarray] = []
    extra_names: list[str] = []
    if mt.rx_type in ("land", "amphibious"):
        tilted = np.where(beta != 0)[0]
        for idx in tilted:
            hrx = np.array(
                [x_rx[idx], y_rx[idx], z_rx[idx], 0.0, 0.0, 0.0, 0.0, 0.0]
            )
            extra_rx.append(hrx)
            base_name = (
                mt.rx_name[idx]
                if mt.rx_name and idx < len(mt.rx_name)
                else f"MT{idx + 1}"
            )
            extra_names.append(f"{base_name}_B")
            irx_mag[idx] = n_rx + len(extra_rx) - 1

    rx_arr = np.column_stack(
        [x_rx, y_rx, z_rx, theta, alpha, beta, np.zeros(n_rx), np.zeros(n_rx)]
    )
    if extra_rx:
        rx_arr = np.vstack([rx_arr, np.array(extra_rx)])

    rx_names: list[str]
    if mt.rx_name:
        rx_names = list(mt.rx_name)
    else:
        rx_names = [f"MT{i + 1}" for i in range(n_rx)]
    rx_names.extend(extra_names)

    mt_cfg = MTConfig(
        frequencies=freqs,
        receivers=rx_arr,
        receiver_name=rx_names,
    )

    # Build DATA rows
    nf = len(freqs)
    rows: list[list[float]] = []
    data_codes_map = {
        "lTE": [(123, True), (104, False)],  # True = use irx_mag for H
        "lTM": [(125, True), (106, False)],
        "lZDet": [(129, True), (110, False)],
        "lTipperRealImag": [(133, True), (134, True)],
        "lTipper": [(135, True), (136, True)],
        "lMTFields": [
            (151, False),
            (152, False),
            (153, False),
            (154, False),
            (155, False),
            (156, False),
            (161, False),
            (162, False),
            (163, False),
            (164, False),
            (165, False),
            (166, False),
        ],
    }

    for ifreq in range(1, nf + 1):
        for irx in range(1, n_rx + 1):
            itx = int(irx_mag[irx - 1]) + 1
            for flag, entries in data_codes_map.items():
                if not getattr(mt, flag, False):
                    continue
                for code, use_mag in entries:
                    h_idx = itx if use_mag else irx
                    e_idx = irx
                    if code in (123, 104, 113, 114, 133, 134, 135, 136):
                        rows.append([code, ifreq, h_idx, e_idx, 0.0, 0.0])
                    elif code in (125, 106, 115, 116):
                        rows.append([code, ifreq, h_idx, e_idx, 0.0, 0.0])
                    else:
                        rows.append([code, ifreq, e_idx, e_idx, 0.0, 0.0])

    data = (
        np.array(rows, dtype=float) if rows else np.empty((0, 6), dtype=float)
    )
    return mt_cfg, data


# ---------------------------------------------------------------------------
# CSEM data builder
# ---------------------------------------------------------------------------


def _build_csem_config(
    topo: Any,
    csem: CSEMSurveyConfig,
) -> tuple[CSEMConfig, np.ndarray]:
    """Build CSEMConfig and DATA rows for CSEM forward modeling."""
    y_tx = np.asarray(csem.tx_y, dtype=float).ravel()
    nf = len(np.asarray(csem.frequencies, dtype=float))
    freqs = np.asarray(csem.frequencies, dtype=float)

    # Transmitters
    txs, tx_types = _place_transmitters(
        topo,
        y_tx,
        tx_z=csem.tx_z,
        tx_z_offset=csem.tx_z_offset,
        tx_azimuth=csem.tx_azimuth,
        tx_dip=csem.tx_dip,
        tx_length=csem.tx_length,
        tx_type=csem.tx_type,
    )
    tx_names = csem.tx_name or [f"Tx{i + 1}" for i in range(len(y_tx))]
    # pad to 7 columns
    if txs.shape[1] < 7:
        txs = np.hstack([txs, np.zeros((len(txs), 7 - txs.shape[1]))])

    # Receivers (towed or fixed)
    towed = csem.rx_r is not None and len(csem.rx_r) > 0
    if towed:
        r_rx = np.asarray(csem.rx_r, dtype=float).ravel()
        n_rx_per_tx = len(r_rx)
        tow_sign = np.sign(np.mean(np.diff(y_tx))) if len(y_tx) > 1 else 1
        y_rx = np.repeat(y_tx, n_rx_per_tx) - tow_sign * np.tile(
            r_rx, len(y_tx)
        )
    else:
        y_rx = np.asarray(csem.rx_y, dtype=float).ravel()

    n_rx = len(y_rx)
    x_rx, z_rx, theta_rx, alpha_rx, beta_rx = _place_receivers(
        topo,
        y_rx,
        csem.rx_type,
        rx_z=csem.rx_z,
        rx_z_offset=csem.rx_z_offset,
        rx_beta=csem.rx_beta,
    )

    # Extra H-field receivers for tilted land stations
    irx_mag = np.arange(n_rx)
    extra_rx: list[np.ndarray] = []
    has_b = csem.lBx or csem.lBy or csem.lBz
    if has_b and csem.rx_type in ("land", "amphibious"):
        tilted = np.where(beta_rx != 0)[0]
        for idx in tilted:
            hrx = np.array(
                [x_rx[idx], y_rx[idx], z_rx[idx], 0.0, 0.0, 0.0, 0.0, 0.0]
            )
            extra_rx.append(hrx)
            irx_mag[idx] = n_rx + len(extra_rx) - 1
    elif not has_b:
        # drop duplicate H receivers
        pass

    rx_arr = np.column_stack(
        [
            x_rx,
            y_rx,
            z_rx,
            theta_rx,
            alpha_rx,
            beta_rx,
            np.zeros(n_rx),
            np.zeros(n_rx),
        ]
    )
    if extra_rx:
        rx_arr = np.vstack([rx_arr, np.array(extra_rx)])

    # optional dipole length column
    if csem.rx_length is not None:
        rl = _broadcast(csem.rx_length, n_rx)
        rx_arr[:n_rx, 6] = rl

    rx_names = csem.rx_name or [f"Rx{i + 1}" for i in range(n_rx)]

    n_tx = len(txs)

    # Range limits per frequency
    min_range = np.zeros(nf)
    max_range = np.full(nf, 1e300)
    if csem.min_range is not None:
        min_range[:] = csem.min_range
    if csem.max_range is not None:
        max_range[:] = csem.max_range

    # Build DATA block
    csem_data_flags = {
        "lEx": [(27, 22)],
        "lEy": [(28, 24)],
        "lEz": [(29, 26)],
        "lBx": [(37, 32)],
        "lBy": [(38, 34)],
        "lBz": [(39, 36)],
    }
    rows: list[list[float]] = []

    for ifreq in range(1, nf + 1):
        for itx in range(1, n_tx + 1):
            for _irx_local in range(n_rx):
                irx = _irx_local + 1  # 1-based
                if towed:
                    n_rx_per_tx_local = n_rx // n_tx
                    irx = (
                        (_irx_local % n_rx_per_tx_local)
                        + (itx - 1) * n_rx_per_tx_local
                        + 1
                    )

                # range check
                tx_pos = txs[itx - 1, :3]
                rx_pos = rx_arr[irx - 1, :3]
                rng = np.linalg.norm(tx_pos - rx_pos)
                if rng < min_range[ifreq - 1] or rng > max_range[ifreq - 1]:
                    continue

                for flag, code_pairs in csem_data_flags.items():
                    if not getattr(csem, flag, False):
                        continue
                    for amp_code, phase_code in code_pairs:
                        if flag in ("lBx", "lBy", "lBz"):
                            irx_b = int(irx_mag[irx - 1]) + 1
                        else:
                            irx_b = irx
                        rows.append([amp_code, ifreq, itx, irx_b, 0.0, 0.0])
                        rows.append([phase_code, ifreq, itx, irx_b, 0.0, 0.0])

    data = (
        np.array(rows, dtype=float) if rows else np.empty((0, 6), dtype=float)
    )

    csem_cfg = CSEMConfig(
        phase_convention=csem.phase_convention,
        frequencies=freqs,
        transmitters=txs,
        transmitter_type=tx_types,
        transmitter_name=tx_names,
        receivers=rx_arr,
        receiver_name=rx_names,
    )
    return csem_cfg, data


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------


def make_data_file(
    out_file: str | Path,
    topo: Any,
    *,
    mt: MTSurveyConfig | None = None,
    csem: CSEMSurveyConfig | None = None,
    file_type: str = "forward",
    comment: str = "Created by pycsamt make_data_file",
    utm: UTMOrigin | None = None,
) -> EMDataFile:
    """Create a MARE2DEM ``.emdata`` file from survey parameters.

    Port of ``m2d_makeDataFile.m``.

    Parameters
    ----------
    out_file : path-like
        Output ``.emdata`` filename.
    topo : float or array-like (n, 2)
        Topography: scalar depth (flat) or ``[y, z]`` table.
    mt : MTSurveyConfig, optional
        MT receiver/data configuration.
    csem : CSEMSurveyConfig, optional
        CSEM transmitter/receiver/data configuration.
    file_type : str, default "forward"
        ``"forward"`` or ``"inversion"`` (inversion data not yet
        implemented).
    comment : str
        Optional comment written on the second line of the file.
    utm : UTMOrigin, optional
        UTM mapping metadata.

    Returns
    -------
    EMDataFile
        The constructed data object (also written to *out_file*).

    Examples
    --------
    MT forward dataset at 10 frequencies, 20 land stations:

    >>> import numpy as np
    >>> from pycsamt.models.mare2dem.survey import MTSurveyConfig, make_data_file
    >>> freqs = np.logspace(-3, 3, 10)
    >>> rx_y  = np.linspace(-10000, 10000, 20)
    >>> mt_cfg = MTSurveyConfig(frequencies=freqs, rx_y=rx_y, rx_type="land", lTE=True, lTM=True)
    >>> em = make_data_file("survey.emdata", topo=0.0, mt=mt_cfg)
    >>> em.n_mt_receivers
    20
    """
    em = EMDataFile(comment=comment)
    em.utm = utm or UTMOrigin()

    all_data: list[np.ndarray] = []

    if mt is not None:
        mt_cfg, mt_data = _build_mt_config(topo, mt)
        em.mt = mt_cfg
        if len(mt_data):
            all_data.append(mt_data)

    if csem is not None:
        csem_cfg, csem_data = _build_csem_config(topo, csem)
        em.csem = csem_cfg
        if len(csem_data):
            all_data.append(csem_data)

    em.data = (
        np.vstack(all_data) if all_data else np.empty((0, 6), dtype=float)
    )

    write_emdata(em, out_file)
    return em
