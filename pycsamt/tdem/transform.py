# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""
TEM time-domain → frequency-domain transforms.

Two transform classes are provided:

:class:`LateTimeTransform`
    Fast, approximate method based on the late-time apparent
    resistivity formula and a pseudo-frequency mapping.
    Suitable for routine 1-D and constrained 2-D interpretation.

:class:`FourierTransform`
    Rigorous conversion via Hankel / sine digital-filter
    coefficients (Meju 1996, Christensen 1990).
    **Not yet implemented** — stub only.

:class:`TEMtoEDI`
    High-level dispatcher that wraps either transform and
    produces an :class:`~pycsamt.seg.collection.EDICollection`.
"""

from __future__ import annotations

import warnings
from typing import Dict, List, Optional, Sequence, Union

import numpy as np

from ._base import TEMSounding

__all__ = ["LateTimeTransform", "FourierTransform", "TEMtoEDI"]

MU0 = 4.0 * np.pi * 1e-7          # H/m (magnetic permeability of free space)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _rho_a_late_time(
    dBdt: np.ndarray,
    t: np.ndarray,
    moment: float,
) -> np.ndarray:
    r"""
    Late-time apparent resistivity for a coincident / central-loop system.

    Derived from the late-time asymptotic response of a horizontal
    circular loop source over a uniform half-space:

    .. math::

        \frac{\partial B_z}{\partial t}
            = -\frac{M \mu_0^{5/2}\,\sigma^{3/2}}
                    {10\,\sqrt{\pi}\,t^{5/2}}

    Solving for apparent resistivity :math:`\rho_a = 1/\sigma`:

    .. math::

        \rho_a(t) =
            \left(\frac{M\,\mu_0^{5/2}}
                  {10\,\sqrt{\pi}\,|\partial B_z/\partial t|\,t^{5/2}}
            \right)^{2/3}

    Parameters
    ----------
    dBdt : np.ndarray
        :math:`\partial B_z/\partial t` in T s⁻¹.  Shape ``(n,)``.
    t : np.ndarray
        Time gates in seconds.  Shape ``(n,)``.
    moment : float
        Transmitter magnetic moment :math:`M = I\,n\,A` in A m².

    Returns
    -------
    np.ndarray
        Apparent resistivity in Ω m.  Shape ``(n,)``.

    References
    ----------
    .. [1] Ward & Hohmann (1988), "Electromagnetic Theory for Geophysical
       Applications", Chapter 4, eq. 4.96.
    .. [2] Nabighian & Macnae (1991), "Time Domain Electromagnetic
       Prospecting Methods", *Electromagnetic Methods in Applied
       Geophysics*, Vol. 2, pp. 427–520.
    .. [3] Christiansen et al. (2009), "An open-source code for TEM
       forward modeling", *Computers & Geosciences*, 35, 1475–1480.
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        arg = (moment * MU0 ** 2.5) / (
            10.0 * np.sqrt(np.pi) * np.abs(dBdt) * t ** 2.5
        )
        rho = arg ** (2.0 / 3.0)
    rho = np.where(np.isfinite(rho) & (rho > 0.0), rho, np.nan)
    return rho


def _pseudo_freq(t: np.ndarray, convention: str = "skin_depth") -> np.ndarray:
    r"""
    Map time gates to equivalent pseudo-frequencies.

    ``"skin_depth"`` (default):
        :math:`f = 1 / (2\pi t)` — skin-depth equivalence.

    ``"diffusion"`` :
        :math:`f = 1 / (4\pi t)` — diffusion-depth equivalence.

    Parameters
    ----------
    t : np.ndarray
        Time gates in seconds.
    convention : str
        ``"skin_depth"`` or ``"diffusion"``.
    """
    if convention == "skin_depth":
        return 1.0 / (2.0 * np.pi * t)
    if convention == "diffusion":
        return 1.0 / (4.0 * np.pi * t)
    raise ValueError(f"Unknown convention '{convention}'. Use 'skin_depth' or 'diffusion'.")


def _phase_from_rho(
    rho_a: np.ndarray,
    freq: np.ndarray,
    mode: str = "homogeneous",
) -> np.ndarray:
    r"""
    Estimate MT phase from the apparent resistivity curve.

    ``"homogeneous"``
        Assumes a uniform half-space: :math:`\phi_{xy} = 45°`
        everywhere.

    ``"weidelt"``
        Uses the Weidelt–Boehl dispersion relation (Weidelt 1972,
        Boehl et al. 1983):

        .. math::

            \phi(\omega) \approx \frac{\pi}{4}
                + \frac{1}{2}
                  \frac{\partial \ln\rho_a}{\partial \ln\omega}

        Requires at least 3 valid (non-NaN) data points.

    Parameters
    ----------
    rho_a : np.ndarray
        Apparent resistivity in Ω m.  Shape ``(n,)``.
    freq : np.ndarray
        Frequencies in Hz.  Shape ``(n,)``.
    mode : str
        ``"homogeneous"`` or ``"weidelt"``.

    Returns
    -------
    np.ndarray
        Phase of Z_xy in degrees (positive, 0–90 range for 1-D earth).
    """
    if mode == "homogeneous":
        return np.full_like(rho_a, 45.0)

    if mode == "weidelt":
        log_rho = np.log(rho_a)
        log_omega = np.log(2.0 * np.pi * freq)

        # numerical derivative using central differences
        d_log_rho = np.gradient(log_rho, log_omega)

        # Weidelt relation: φ = 45° + 90°/2 × d(ln ρ)/d(ln ω)
        phase_rad = np.pi / 4.0 + 0.5 * d_log_rho
        phase_deg = np.degrees(phase_rad)

        # clip to physically reasonable range [0°, 90°]
        phase_deg = np.clip(phase_deg, 0.0, 90.0)

        # propagate NaNs
        mask = ~np.isfinite(rho_a)
        phase_deg[mask] = np.nan

        return phase_deg

    raise ValueError(f"Unknown phase mode '{mode}'. Use 'homogeneous' or 'weidelt'.")


def _build_z_array(
    rho_a: np.ndarray,
    phase_xy_deg: np.ndarray,
    freq: np.ndarray,
) -> np.ndarray:
    r"""
    Construct complex impedance tensor array from apparent
    resistivity and phase.

    Assumes 1-D earth:
    - :math:`Z_{xx} = Z_{yy} = 0`
    - :math:`Z_{yx} = -Z_{xy}` (anti-symmetric)

    The relationship between apparent resistivity, angular
    frequency and impedance magnitude is:

    .. math::

        |Z_{xy}(\omega)| = \sqrt{\rho_a(\omega)\,\omega\,\mu_0}

    Parameters
    ----------
    rho_a : np.ndarray
        Apparent resistivity in Ω m.  Shape ``(n,)``.
    phase_xy_deg : np.ndarray
        Phase of :math:`Z_{xy}` in degrees.  Shape ``(n,)``.
    freq : np.ndarray
        Frequencies in Hz.  Shape ``(n,)``.

    Returns
    -------
    np.ndarray
        Complex Z tensor of shape ``(n, 2, 2)``.
        Order: ``[Zxx, Zxy; Zyx, Zyy]``.
    """
    omega = 2.0 * np.pi * freq
    Z_mag = np.sqrt(rho_a * omega * MU0)          # |Z_xy| in mV/km/nT (SI: Ω)

    phi_rad = np.radians(phase_xy_deg)
    Z_xy = Z_mag * (np.cos(phi_rad) + 1j * np.sin(phi_rad))

    n = len(freq)
    Z_arr = np.zeros((n, 2, 2), dtype=complex)
    Z_arr[:, 0, 1] = Z_xy           # Zxy
    Z_arr[:, 1, 0] = -Z_xy          # Zyx = -Zxy  (1-D assumption)

    return Z_arr


def _z_error_from_rho_error(
    rho_a: np.ndarray,
    rho_err: np.ndarray,
    freq: np.ndarray,
) -> np.ndarray:
    r"""
    Propagate apparent-resistivity uncertainty to impedance uncertainty.

    From :math:`|Z| = \sqrt{\rho_a \omega \mu_0}`:

    .. math::

        \delta|Z| = \frac{\delta\rho_a}{2\rho_a} \cdot |Z|

    Returns
    -------
    np.ndarray
        Z error tensor, shape ``(n, 2, 2)``, dtype float.
        Only the off-diagonal elements are non-zero.
    """
    omega = 2.0 * np.pi * freq
    Z_mag = np.sqrt(rho_a * omega * MU0)
    with np.errstate(divide="ignore", invalid="ignore"):
        dZ = np.where(rho_a > 0.0, rho_err / (2.0 * rho_a) * Z_mag, np.nan)

    n = len(freq)
    Z_err = np.zeros((n, 2, 2))
    Z_err[:, 0, 1] = dZ
    Z_err[:, 1, 0] = dZ

    return Z_err


# ---------------------------------------------------------------------------
# LateTimeTransform
# ---------------------------------------------------------------------------

class LateTimeTransform:
    r"""
    Convert TEM soundings to frequency-domain apparent impedance
    using the late-time apparent-resistivity approximation.

    This is the standard industry approach (Meju 1996; Spies &
    Frischknecht 1991).  It is valid when measurement times are
    significantly later than the transmitter ramp-off time.

    Parameters
    ----------
    freq_convention : str, optional
        How to map each time gate to an equivalent frequency:

        * ``"skin_depth"`` (default): :math:`f = 1 / (2\pi t)`
        * ``"diffusion"`` : :math:`f = 1 / (4\pi t)`

    phase_mode : str, optional
        How to estimate the MT phase from the apparent resistivity:

        * ``"homogeneous"`` (default): :math:`\phi_{xy} = 45°`
        * ``"weidelt"`` : dispersion relation
          (Weidelt 1972, Boehl et al. 1983)

    drop_nan : bool, optional
        If ``True`` (default), remove gates where the apparent
        resistivity could not be computed (zero / negative data,
        NaN, or zero time).

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.tdem import TEMSounding
    >>> from pycsamt.tdem.transform import LateTimeTransform
    >>> t = np.logspace(-5, -2, 30)
    >>> M = 8.0 * 100.0**2          # I=8 A, square loop 100 m
    >>> dBdt = (4e-7/np.pi) * (4e-7**2 * M / (20 * 10.0 ** (3./2.) * t**(5./2.))) ** (2./3.)
    >>> snd = TEMSounding(t, dBdt, current=8.0, tx_area=100.**2)
    >>> tr = LateTimeTransform()
    >>> result = tr.transform(snd)
    >>> result["freq"].shape
    (30,)
    """

    def __init__(
        self,
        freq_convention: str = "skin_depth",
        phase_mode: str = "homogeneous",
        drop_nan: bool = True,
    ) -> None:
        self.freq_convention = freq_convention
        self.phase_mode = phase_mode
        self.drop_nan = drop_nan

    def transform(self, sounding: TEMSounding) -> Dict:
        r"""
        Transform one :class:`~pycsamt.tdem.TEMSounding` to
        frequency-domain data arrays.

        Parameters
        ----------
        sounding : TEMSounding
            Input TEM sounding.

        Returns
        -------
        dict with keys:

        ``freq``
            Pseudo-frequencies in Hz.  Shape ``(n,)``.
        ``Z``
            Complex impedance tensor.  Shape ``(n, 2, 2)``.
        ``Z_err``
            Impedance uncertainty tensor (all NaN if no error
            in input).  Shape ``(n, 2, 2)``.
        ``rho_a``
            Apparent resistivity in Ω m.  Shape ``(n,)``.
        ``phase_xy``
            Phase of Z_xy in degrees.  Shape ``(n,)``.
        ``station_name``
            Station identifier string.
        ``x``, ``y``, ``elevation``
            Station coordinates.
        """
        dBdt = sounding.dBdt()
        t = sounding.time_gates
        M = sounding.moment

        rho_a = _rho_a_late_time(dBdt, t, M)
        freq = _pseudo_freq(t, self.freq_convention)

        if self.drop_nan:
            mask = np.isfinite(rho_a) & np.isfinite(freq) & (freq > 0.0)
            rho_a = rho_a[mask]
            freq = freq[mask]
            t = t[mask]
            dBdt = dBdt[mask]
            err_raw = sounding.error[mask] if sounding.error is not None else None
        else:
            err_raw = sounding.error

        # sort by ascending frequency (= descending time)
        order = np.argsort(freq)
        freq = freq[order]
        rho_a = rho_a[order]
        if err_raw is not None:
            err_raw = err_raw[order]

        # propagate ρ_a error → Z_err
        if err_raw is not None:
            # dBdt error → rho_a error via the 2/3-power law:
            # δρ_a / ρ_a = (2/3) × δ(dBdt) / |dBdt|
            rho_err = (2.0 / 3.0) * rho_a * np.abs(err_raw / (dBdt[order] + 1e-300))
        else:
            rho_err = None

        phase_xy = _phase_from_rho(rho_a, freq, mode=self.phase_mode)
        Z_arr = _build_z_array(rho_a, phase_xy, freq)

        if rho_err is not None:
            Z_err = _z_error_from_rho_error(rho_a, rho_err, freq)
        else:
            Z_err = np.full_like(Z_arr, np.nan, dtype=float)

        return {
            "freq":         freq,
            "Z":            Z_arr,
            "Z_err":        Z_err,
            "rho_a":        rho_a,
            "phase_xy":     phase_xy,
            "station_name": sounding.station_name,
            "x":            sounding.x,
            "y":            sounding.y,
            "elevation":    sounding.elevation,
        }

    def transform_many(
        self,
        soundings: Sequence[TEMSounding],
    ) -> List[Dict]:
        """Transform a list of soundings.  Returns a list of result dicts."""
        return [self.transform(s) for s in soundings]


# ---------------------------------------------------------------------------
# FourierTransform  (stub — not yet implemented)
# ---------------------------------------------------------------------------

class FourierTransform:
    r"""
    Rigorous TDEM → MT impedance conversion via Hankel / sine
    digital-filter transform (Meju 1996; Christensen 1990).

    This method uses the exact mathematical relationship between the
    TDEM step-off response and the MT impedance for a 1-D layered
    earth:

    .. math::

        Z(\omega) = i\omega\mu_0\,K(\omega)

    where :math:`K(\omega)` is recovered from the measured
    :math:`\partial B_z/\partial t(t)` via the inverse Fourier /
    digital-filter convolution.  Requires accurate knowledge of the
    transmitter waveform (``sounding.waveform``).

    .. note::
        **Not yet implemented.**  The digital-filter coefficient set
        (Anderson 1982, Christensen 1990) and waveform deconvolution
        will be added in a future release.  Use
        :class:`LateTimeTransform` for now.

    Parameters
    ----------
    n_filters : int, optional
        Number of digital-filter coefficients.  Default 61.
    waveform_correction : bool, optional
        If ``True``, apply waveform deconvolution before the
        transform.  Default ``True``.

    Raises
    ------
    NotImplementedError
        Always, until the implementation is complete.
    """

    def __init__(
        self,
        n_filters: int = 61,
        waveform_correction: bool = True,
    ) -> None:
        self.n_filters = n_filters
        self.waveform_correction = waveform_correction

    def transform(self, sounding: TEMSounding) -> Dict:
        raise NotImplementedError(
            "FourierTransform is not yet implemented. "
            "Use LateTimeTransform for the approximate method."
        )

    def transform_many(self, soundings: Sequence[TEMSounding]) -> List[Dict]:
        raise NotImplementedError(
            "FourierTransform is not yet implemented. "
            "Use LateTimeTransform for the approximate method."
        )


# ---------------------------------------------------------------------------
# TEMtoEDI — high-level dispatcher
# ---------------------------------------------------------------------------

class TEMtoEDI:
    r"""
    Convert one or more TEM soundings to a
    :class:`~pycsamt.seg.collection.EDICollection`.

    This is the main user-facing class for the TEM→frequency-domain
    pipeline.  It wraps either :class:`LateTimeTransform` (default)
    or :class:`FourierTransform` and writes one synthetic EDI file per
    sounding into ``out_dir``.

    Parameters
    ----------
    method : str, optional
        Transform method:

        * ``"late_time"`` (default) — :class:`LateTimeTransform`
        * ``"fourier"`` — :class:`FourierTransform` (not yet available)

    freq_convention : str, optional
        Passed to :class:`LateTimeTransform`.
        ``"skin_depth"`` (default) or ``"diffusion"``.
    phase_mode : str, optional
        Passed to :class:`LateTimeTransform`.
        ``"homogeneous"`` (default) or ``"weidelt"``.
    out_dir : str or Path, optional
        Directory where synthetic EDI files are saved.  Created if
        absent.  Default ``"edi_out/tem"``.
    verbose : int, optional
        Verbosity level.  Default 0.

    Examples
    --------
    Single sounding:

    >>> from pycsamt.tdem import TEMSounding, TEMtoEDI
    >>> import numpy as np
    >>> t = np.logspace(-5, -2, 25)
    >>> dBdt = 5e-5 * t ** (-5./2.)
    >>> snd = TEMSounding(t, dBdt, current=8.0, tx_area=1e4,
    ...                   station_name="S01", x=1000.0, y=500.0)
    >>> conv = TEMtoEDI(method="late_time", out_dir="/tmp/edi_tem")
    >>> coll = conv.transform(snd)

    Profile of soundings:

    >>> soundings = [snd1, snd2, snd3]
    >>> coll = conv.transform_many(soundings)
    """

    def __init__(
        self,
        method: str = "late_time",
        freq_convention: str = "skin_depth",
        phase_mode: str = "homogeneous",
        out_dir: Union[str, "Path"] = "edi_out/tem",
        verbose: int = 0,
    ) -> None:
        self.method = method.lower()
        self.out_dir = out_dir
        self.verbose = verbose

        if self.method == "late_time":
            self._transformer = LateTimeTransform(
                freq_convention=freq_convention,
                phase_mode=phase_mode,
            )
        elif self.method == "fourier":
            self._transformer = FourierTransform()
        else:
            raise ValueError(
                f"Unknown method '{method}'. Use 'late_time' or 'fourier'."
            )

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    def transform(self, sounding: TEMSounding):
        r"""
        Convert a single TEM sounding and return an
        :class:`~pycsamt.seg.collection.EDICollection` with one site.

        Parameters
        ----------
        sounding : TEMSounding

        Returns
        -------
        EDICollection
        """
        result = self._transformer.transform(sounding)
        edi = self._result_to_edifile(result)
        return self._make_collection([edi])

    def transform_many(self, soundings: Sequence[TEMSounding]):
        r"""
        Convert a list of TEM soundings and return an
        :class:`~pycsamt.seg.collection.EDICollection`.

        Parameters
        ----------
        soundings : sequence of TEMSounding

        Returns
        -------
        EDICollection
        """
        results = self._transformer.transform_many(soundings)
        edis = [self._result_to_edifile(r) for r in results]
        return self._make_collection(edis)

    def save(
        self,
        sounding_or_soundings,
        path=None,
    ):
        r"""
        Transform and write EDI files to ``out_dir`` (or ``path``).

        Returns the list of written file paths.
        """
        from pathlib import Path as _Path

        out = _Path(self.out_dir if path is None else path)
        out.mkdir(parents=True, exist_ok=True)

        if isinstance(sounding_or_soundings, TEMSounding):
            soundings = [sounding_or_soundings]
        else:
            soundings = list(sounding_or_soundings)

        written = []
        for snd in soundings:
            result = self._transformer.transform(snd)
            edi = self._result_to_edifile(result)
            name = (result["station_name"] or "tem_site") + ".edi"
            fp = str(out / name)
            edi.write(new_edifn=name, savepath=str(out))
            written.append(fp)
            if self.verbose:
                print(f"  Wrote {fp}  ({result['freq'].size} freq)")

        return written

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _result_to_edifile(self, result: Dict):
        """
        Build an :class:`~pycsamt.seg.edi.EDIFile` from a transform
        result dict, populating Z, freq, and HEAD/INFO metadata.
        """
        from ..seg.edi import EDIFile
        from ..z.z import Z

        edi = EDIFile()

        # populate Z object
        edi.Z.freq = result["freq"]
        edi.Z.z = result["Z"]
        edi.Z.z_err = result.get("Z_err")

        # populate HEAD section via the head section object
        self._set_head(edi, result)

        return edi

    @staticmethod
    def _set_head(edi, result: Dict) -> None:
        """Populate the EDI HEAD / INFO sections with TEM provenance."""
        try:
            from ..seg.heads import Head
            head = Head()
            head.dataid = result["station_name"] or "TEM_SITE"
            head.lon = float(result.get("x", 0.0))
            head.lat = float(result.get("y", 0.0))
            head.elev = float(result.get("elevation", 0.0))
            head.acqby = "pyCSAMT.tdem"
            head.fileby = "pyCSAMT.tdem"
            edi.add_section("head", head)
        except Exception:
            pass  # HEAD is optional; transform result is still usable

        try:
            from ..seg.heads import Info
            info = Info()
            info.info_list = [
                "Data origin: TEM (time-domain EM) transformed to MT equivalent",
                "Transform: pycsamt.tdem.transform.TEMtoEDI",
                "Assumption: 1-D earth, coincident/central loop",
                "Z_xx = Z_yy = 0 (1-D assumption); Z_yx = -Z_xy",
            ]
            edi.add_section("info", info)
        except Exception:
            pass

    @staticmethod
    def _make_collection(edis: list):
        """Wrap a list of EDIFile objects in an EDICollection."""
        from ..seg.collection import EDICollection
        return EDICollection(items=edis)
