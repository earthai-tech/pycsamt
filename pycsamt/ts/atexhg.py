# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Reader for ``<ATEXHG>`` binary MT/AMT time series (``.atts``).

The ``<ATEXHG>`` container is written by some Chinese AMT/MT
acquisition systems and read by the ``prMT`` processing package
(it was met in a 2026 survey whose EDI files carry
``PROGVERS="prMT version 1016"``). No public specification is
known; the layout below was **inferred** from 123 field files and
is checked on every read (magic bytes, sample count against file
size, channel records). Fields whose meaning is unknown are kept
verbatim in :attr:`AtexhgHeader.raw`.

Layout (little-endian)
----------------------
==========  ======  ==================================================
offset      type    content
==========  ======  ==================================================
0x00        8 B     magic ``b"<ATEXHG>"``
0x08        u16     number of recorded channels ``n``
0x0A        u16     bits per sample (32 = IEEE float32)
0x10        f64     sampling rate (Hz)
0x20        u32     unknown (599 in all files seen)
0x24        u32     unknown (999750 in all files seen)
0x28        u32     total samples, all channels (``n * n_samples``)
0x58        54 B    channel record *i* at ``0x58 + 54 i``:
                    ``+2`` u16 kind (1 electric, 2 magnetic),
                    ``+4`` u16 component index (1 x, 2 y, 3 z),
                    ``+6`` f32 dipole length (m, electric only),
                    ``+10`` 6 B name, ``+30`` i32 code (varies
                    between stations; meaning unknown, kept as
                    ``gain_code``), ``+38`` f32 azimuth (deg)
end - 4N    f32     samples, channel-interleaved, ``N`` = total
==========  ======  ==================================================

Samples are the recorded voltages (mV). Magnetic channels are the
raw induction-coil output and must be calibrated (see
:func:`pycsamt.ts.calibrate.apply_sensor_response`) before they
are in nT; electric channels are divided by the dipole length only
when converting to field units downstream.
"""

from __future__ import annotations

import struct
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from ..exceptions import FileHandlingError
from ._base import TSData

__all__ = ["AtexhgHeader", "is_atexhg", "read_atexhg_header",
           "read_atexhg", "gain_from_code"]

MAGIC = b"<ATEXHG>"
_REC0 = 0x58
_RECLEN = 54
_KIND = {1: "E", 2: "H"}
_COMP = {1: "X", 2: "Y", 3: "Z"}


@dataclass
class AtexhgHeader:
    """Parsed ``<ATEXHG>`` header.

    Attributes
    ----------
    n_channels : int
        Number of recorded channels.
    bits : int
        Bits per sample (only 32, float32, is supported).
    sampling_rate : float
        Sampling rate in Hz.
    n_samples : int
        Samples per channel.
    data_offset : int
        Byte offset of the first sample.
    channels : list of dict
        One dict per channel with keys ``id`` (``"EX"``, ``"HY"``
        ...), ``kind``, ``index``, ``name``, ``dipole_m``,
        ``azimuth_deg`` and ``gain_code``.
    raw : dict
        Header words of unknown meaning, kept verbatim.
    """

    n_channels: int
    bits: int
    sampling_rate: float
    n_samples: int
    data_offset: int
    channels: list[dict] = field(default_factory=list)
    raw: dict = field(default_factory=dict)


def is_atexhg(path: str | Path) -> bool:
    """Return ``True`` when *path* starts with the ATEXHG magic."""
    try:
        with Path(path).open("rb") as f:
            return f.read(len(MAGIC)) == MAGIC
    except OSError:
        return False


def read_atexhg_header(path: str | Path) -> AtexhgHeader:
    """Parse and validate the header of an ``.atts`` file.

    Parameters
    ----------
    path : str or Path
        ``<ATEXHG>`` file.

    Returns
    -------
    AtexhgHeader

    Raises
    ------
    FileHandlingError
        If the magic bytes are wrong, the sample format is not
        float32, a channel record is malformed, or the sample count
        in the header does not fit the file size.
    """
    p = Path(path)
    size = p.stat().st_size
    with p.open("rb") as f:
        head = f.read(4096)
    if head[: len(MAGIC)] != MAGIC:
        raise FileHandlingError(f"{p}: not an <ATEXHG> file.")
    n_ch, bits = struct.unpack_from("<HH", head, 0x08)
    if bits != 32:
        raise FileHandlingError(
            f"{p}: {bits}-bit samples are not supported (float32 "
            "only).")
    if n_ch < 1 or _REC0 + _RECLEN * n_ch > len(head):
        raise FileHandlingError(f"{p}: invalid channel count {n_ch}.")
    fs = struct.unpack_from("<d", head, 0x10)[0]
    if not np.isfinite(fs) or fs <= 0:
        raise FileHandlingError(f"{p}: invalid sampling rate {fs}.")
    w20, w24, total = struct.unpack_from("<III", head, 0x20)
    if total == 0 or total % n_ch:
        raise FileHandlingError(
            f"{p}: total sample count {total} is not a multiple of "
            f"{n_ch} channels.")
    offset = size - 4 * total
    if offset < _REC0 + _RECLEN * n_ch:
        raise FileHandlingError(
            f"{p}: header announces {total} samples, more than the "
            f"file of {size} bytes can hold.")
    channels, seen = [], set()
    for i in range(n_ch):
        o = _REC0 + _RECLEN * i
        if head[o:o + 2] != b"\xff\xff":
            raise FileHandlingError(
                f"{p}: channel record {i + 1} lacks its marker.")
        kind, idx = struct.unpack_from("<HH", head, o + 2)
        if kind not in _KIND or idx not in _COMP:
            raise FileHandlingError(
                f"{p}: channel {i + 1} has unknown kind/index "
                f"({kind}, {idx}).")
        cid = _KIND[kind] + _COMP[idx]
        if cid in seen:
            raise FileHandlingError(f"{p}: duplicate channel {cid}.")
        seen.add(cid)
        channels.append({
            "id": cid,
            "kind": _KIND[kind],
            "index": idx,
            "name": head[o + 10:o + 16].split(b"\x00")[0]
            .decode("latin-1").strip(),
            "dipole_m": float(struct.unpack_from("<f", head,
                                                 o + 6)[0]),
            "gain_code": int(struct.unpack_from("<i", head,
                                                o + 30)[0]),
            "azimuth_deg": float(struct.unpack_from("<f", head,
                                                    o + 38)[0]),
        })
    return AtexhgHeader(
        n_channels=n_ch, bits=bits, sampling_rate=float(fs),
        n_samples=total // n_ch, data_offset=offset,
        channels=channels, raw={"word_0x20": w20, "word_0x24": w24})


def gain_from_code(code: int) -> float:
    """Channel gain implied by an ATEXHG per-channel code.

    Positive codes are gains (``5`` -> x5), negative codes
    attenuations (``-2`` -> x1/2). **Inferred**, not specified: on
    123 field files the station-to-station scale between impedances
    estimated from these samples and the vendor-processed (prMT)
    impedances followed ``G_H / G_E`` with this rule to within
    ~0.05 decade (codes 5/25 on E, -2/1/5/25 on H).

    Raises
    ------
    ValueError
        For code 0, which has no interpretation.

    Examples
    --------
    >>> from pycsamt.ts.atexhg import gain_from_code
    >>> gain_from_code(25), gain_from_code(-2)
    (25.0, 0.5)
    """
    c = int(code)
    if c == 0:
        raise ValueError("gain code 0 has no interpretation.")
    return float(c) if c > 0 else 1.0 / abs(c)


def read_atexhg(
    path: str | Path,
    *,
    station: str | None = None,
    apply_gain: bool = False,
    verbose: int = 0,
) -> TSData:
    """Read an ``<ATEXHG>`` (``.atts``) time-series file.

    Parameters
    ----------
    path : str or Path
        Input file.
    station : str, optional
        Station name; defaults to the file stem.
    apply_gain : bool, default False
        Divide each channel by :func:`gain_from_code` of its code,
        giving voltages at the sensor rather than after the field
        gain. Needed to compare stations recorded with different
        gain settings (otherwise impedances jump by factors such as
        5 or 50 between stations).
    verbose : int, default 0
        Verbosity forwarded to :class:`TSData`.

    Returns
    -------
    TSData
        Channels ``EX``, ``EY``, ``HX``, ``HY`` (and ``EZ``/``HZ``
        if recorded) as float64 arrays of recorded voltage in mV;
        ``dt = 1 / sampling_rate``; ``azim``, ``dipole`` (electric
        channels), ``units`` (``"mV"``) and, in ``meta``,
        ``gain_code`` per channel, ``gain_applied``,
        ``calibrated=False`` and the parsed :class:`AtexhgHeader`.

    Raises
    ------
    FileHandlingError
        See :func:`read_atexhg_header`.
    ValueError
        If ``apply_gain`` is set and a channel code is 0.

    Notes
    -----
    The magnetic channels are **not** in nT: they are the coil
    output voltage. In the field files examined, the vendor's
    processing scaled them by a frequency-independent constant only
    (consistent with flat-response sensors in 0.1-860 Hz); with a
    calibration table, apply it with
    :func:`pycsamt.ts.calibrate.apply_sensor_response` before
    estimating impedances.

    Examples
    --------
    >>> from pycsamt.ts import read_ts
    >>> ts = read_ts("L1-450.atts")             # doctest: +SKIP
    >>> ts.channels(), 1 / ts.dt                # doctest: +SKIP
    (['EX', 'EY', 'HX', 'HY'], 4000.0)
    """
    p = Path(path)
    hdr = read_atexhg_header(p)
    raw = np.fromfile(p, dtype="<f4", offset=hdr.data_offset)
    raw = raw.reshape(hdr.n_samples, hdr.n_channels)
    data = {c["id"]: raw[:, j].astype(float)
            for j, c in enumerate(hdr.channels)}
    if apply_gain:
        for c in hdr.channels:
            data[c["id"]] /= gain_from_code(c["gain_code"])
    return TSData(
        data=data,
        dt=1.0 / hdr.sampling_rate,
        station=station or p.stem,
        verbose=verbose,
        azim={c["id"]: c["azimuth_deg"] for c in hdr.channels},
        dipole={c["id"]: c["dipole_m"] for c in hdr.channels
                if c["kind"] == "E"},
        units={c["id"]: "mV" for c in hdr.channels},
        source=str(p),
        format="atexhg",
        gain_code={c["id"]: c["gain_code"] for c in hdr.channels},
        gain_applied=bool(apply_gain),
        calibrated=False,
        header=hdr,
    )
