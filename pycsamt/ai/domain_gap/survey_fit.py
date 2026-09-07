# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Fit plausible :class:`~.simulator.CorruptionConfig` ranges from real data.

This module is the bridge between a real field survey (AMT, CSAMT, MT,
or otherwise) and the numpy-only :mod:`pycsamt.ai.domain_gap.simulator`.
It is deliberately split in two:

:func:`survey_data_from_sites`
    Converts EDI/``Sites``/``APISurvey`` input into a canonical
    :class:`~pycsamt.ai.data.contracts.SurveyData`, using only the
    quantities the M1 data contract defines. This path stays numpy-only
    once the bridge itself has run.

:func:`fit_corruption_config`
    Derives plausible noise/error-floor/dropout ranges purely from a
    :class:`~pycsamt.ai.data.contracts.SurveyData`'s own declared errors and
    coverage. No EDI or pandas dependency is required at this stage.

:func:`fit_distortion_priors_from_sites`
    Calls the heavier, pandas-based
    :mod:`pycsamt.emtools.gb`/:mod:`pycsamt.emtools.ss` diagnostics directly
    on real sites to estimate empirical static-shift and galvanic-distortion
    spreads. This is the one part of M3 that genuinely depends on a real
    survey's own QC diagnostics rather than literature defaults.

The frequency grid required by :class:`SurveyData` is *not* silently
interpolated across stations with different sampling: mismatched
frequency grids raise a clear error, per the plan's non-negotiable
principles.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

import numpy as np

from ...emtools._core import ensure_sites
from ...emtools.gb import groom_bailey_table
from ...emtools.ss import estimate_ss_ama
from ..data.contracts import SurveyData
from .simulator import CorruptionConfig

__all__ = [
    "survey_data_from_sites",
    "fit_corruption_config",
    "fit_distortion_priors_from_sites",
    "StructuredMissingness",
    "fit_structured_missingness",
]

_COMPONENT_ORDER = ("xx", "xy", "yx", "yy")
_COMP_IJ = {"xx": (0, 0), "xy": (0, 1), "yx": (1, 0), "yy": (1, 1)}
_DEG_THRESHOLD_M = 1e-5  # ~1 m of latitude/longitude spread

_MV_KM_NT_TO_SI = 4.0e-4 * np.pi
"""Convert EDI-native impedance ([mV/km]/[nT]) to SI (V/A = Ohm).

``pycsamt.site.Site.z`` returns the impedance tensor exactly as
stored in the EDI ``Z`` block, in the SEG/EMTF field convention of
millivolts per kilometre per nanotesla -- not SI. This is the same
factor ``Site.rho``'s own ``rho = 0.2 * |Z|**2 / f`` formula encodes
implicitly (``0.2 == _MV_KM_NT_TO_SI ** 2 / (2 * pi * mu0)`` after
substituting into the SI apparent-resistivity formula
``rho = |Z_SI|**2 / (omega * mu0)``), and the same factor
``pycsamt.forward.maxwell.modem3d`` already documents using when it
converts ModEM's own ``[mV/km]/[nT]`` output back to SI.
"""


def _project_local_metres(
    lats: np.ndarray, lons: np.ndarray, *, station_spacing: float
) -> np.ndarray:
    """Project (lat, lon) degrees to local metres; fall back to a grid.

    Mirrors the equirectangular approximation used by
    :func:`pycsamt.ai.inversion._sites_bridge.sites_to_coords_3d`: adequate
    for the short baselines of a single survey line, not a general-purpose
    geodetic projection.
    """
    n = lats.size
    ok = np.isfinite(lats) & np.isfinite(lons)
    lat_range = float(lats[ok].max() - lats[ok].min()) if ok.any() else 0.0
    lon_range = float(lons[ok].max() - lons[ok].min()) if ok.any() else 0.0
    if ok.sum() >= 2 and (
        lat_range > _DEG_THRESHOLD_M or lon_range > _DEG_THRESHOLD_M
    ):
        lat_ref = float(np.nanmean(lats[ok]))
        lon_ref = float(np.nanmean(lons[ok]))
        cos_lat = float(np.cos(np.radians(lat_ref)))
        x = np.where(ok, (lons - lon_ref) * cos_lat * 111320.0, 0.0)
        y = np.where(ok, (lats - lat_ref) * 111320.0, 0.0)
        return np.column_stack([x, y])
    side = int(np.ceil(np.sqrt(n)))
    xs = np.array([(i % side) * station_spacing for i in range(n)])
    ys = np.array([(i // side) * station_spacing for i in range(n)])
    return np.column_stack([xs, ys])


def survey_data_from_sites(
    sites: Any,
    *,
    crs: str | None = None,
    freq_rtol: float = 1e-6,
    station_spacing: float = 500.0,
    recursive: bool = True,
    on_dup: str = "replace",
    verbose: int = 0,
    metadata: dict[str, Any] | None = None,
) -> SurveyData:
    """Bridge EDI/``Sites``/``APISurvey`` input to canonical ``SurveyData``.

    Parameters
    ----------
    sites : Any
        Anything accepted by :func:`pycsamt.emtools._core.ensure_sites`: a
        filesystem path/glob/directory, ``EDIFile``/``EDICollection``,
        ``Site``/``Sites``, ``APISurvey``, or an iterable of these.
    crs : str, optional
        Coordinate reference system identifier to record. Station positions
        are always projected with a local equirectangular approximation
        (see :func:`pycsamt.ai.inversion._sites_bridge.sites_to_coords_3d`);
        pass a CRS string only if it genuinely describes that projection.
    freq_rtol : float, default=1e-6
        Relative tolerance used when checking that every station shares the
        same frequency grid.
    station_spacing : float, default=500.0
        Forwarded to the coordinate bridge as a uniform-grid fallback
        spacing, used only when no station reports finite coordinates.
    recursive, on_dup, verbose
        Forwarded to ``ensure_sites``.
    metadata : dict, optional
        Extra provenance recorded on the returned survey.

    Returns
    -------
    SurveyData
        Canonical survey with full ``xx, xy, yx, yy`` components, with
        impedance and its declared error converted from ``Site``'s
        EDI-native ``[mV/km]/[nT]`` convention to SI (V/A), matching
        :class:`~pycsamt.ai.data.contracts.SurveyData`'s default
        :class:`~pycsamt.ai.data.contracts.ImpedanceConvention`.

    Raises
    ------
    ValueError
        If no station has usable impedance data, or stations do not share
        a common frequency grid within ``freq_rtol``.

    Notes
    -----
    This function performs no frequency interpolation: a survey whose
    stations were sampled on different frequency grids must be resolved by
    an explicit, survey-matched frequency selector (an M1 concern) before
    reaching this bridge.

    Examples
    --------
    >>> survey = survey_data_from_sites(
    ...     "data/AMT/WILLY_DATA/L18PLT", recursive=False, verbose=0
    ... )  # doctest: +SKIP
    >>> survey.components  # doctest: +SKIP
    ('xx', 'xy', 'yx', 'yy')
    """
    collection = ensure_sites(
        sites, recursive=recursive, on_dup=on_dup, verbose=verbose
    )
    stations = list(collection)
    if not stations:
        raise ValueError("no stations were found in the given sites input.")

    reference_freq: np.ndarray | None = None
    reference_name = None
    kept: list[Any] = []
    for site in stations:
        freq = getattr(site, "freq", None)
        z = getattr(site, "z", None)
        if freq is None or z is None:
            continue
        freq = np.asarray(freq, dtype=float)
        if reference_freq is None:
            reference_freq = freq
            reference_name = site.name
        elif freq.shape != reference_freq.shape or not np.allclose(
            freq, reference_freq, rtol=freq_rtol, atol=0.0
        ):
            raise ValueError(
                f"station {getattr(site, 'name', '?')!r} does not share the "
                f"frequency grid of {reference_name!r}; a survey-matched "
                "frequency selector must resolve this before bridging to "
                "SurveyData."
            )
        kept.append(site)

    if reference_freq is None:
        raise ValueError("no station exposed both freq and z arrays.")

    n_station = len(kept)
    n_frequency = reference_freq.size
    impedance = np.full((n_station, n_frequency, 4), np.nan, dtype=complex)
    error = np.full((n_station, n_frequency, 4), np.nan, dtype=float)
    names: list[str] = []

    for row, site in enumerate(kept):
        names.append(str(site.name))
        z = np.asarray(site.z) * _MV_KM_NT_TO_SI
        z_err = getattr(site, "z_err", None)
        z_err = (
            None
            if z_err is None
            else np.asarray(z_err, dtype=float) * _MV_KM_NT_TO_SI
        )
        if z.ndim == 3 and z.shape[1:] == (2, 2):
            for col, name in enumerate(_COMPONENT_ORDER):
                i, j = _COMP_IJ[name]
                impedance[row, :, col] = z[:, i, j]
                if z_err is not None:
                    error[row, :, col] = z_err[:, i, j]
        elif z.ndim == 2 and z.shape[1] == 4:
            impedance[row] = z
            if z_err is not None:
                error[row] = z_err
        else:
            raise ValueError(
                f"station {site.name!r} has an unsupported z shape {z.shape}."
            )

    lats = np.full(n_station, np.nan)
    lons = np.full(n_station, np.nan)
    elevation = np.full(n_station, np.nan)
    for row, site in enumerate(kept):
        coords = getattr(site, "coords", None)
        if coords is None:
            continue
        try:
            if len(coords) >= 2:
                lats[row] = float(coords[0])
                lons[row] = float(coords[1])
            if len(coords) >= 3:
                elevation[row] = float(coords[2])
        except (TypeError, ValueError):
            pass
    xy = _project_local_metres(lats, lons, station_spacing=station_spacing)
    coordinates = np.column_stack([xy, elevation])

    has_error = np.any(np.isfinite(error))
    return SurveyData(
        impedance=impedance,
        frequencies_hz=reference_freq,
        station_names=tuple(names),
        components=_COMPONENT_ORDER,
        coordinates_m=coordinates,
        impedance_error=error if has_error else None,
        crs=crs,
        metadata=metadata or {},
    )


def fit_corruption_config(
    survey: SurveyData,
    *,
    severity_scale: float = 1.0,
) -> CorruptionConfig:
    """Derive plausible noise/dropout ranges from a real survey's QC.

    Only quantities already present in the canonical
    :class:`~pycsamt.ai.data.contracts.SurveyData` contract are used: the
    ``impedance_error``-to-``|Z|`` ratio for heteroscedastic noise and error
    floor, and :meth:`~pycsamt.ai.data.contracts.SurveyData.coverage` for
    dropout rates.

    Parameters
    ----------
    survey : SurveyData
        Real (or realistically corrupted) survey to profile.
    severity_scale : float, default=1.0
        Multiplier applied to every fitted range/rate, letting a caller
        derive a milder or harsher preset from the same empirical fit.

    Returns
    -------
    CorruptionConfig
        Configuration whose noise range spans the interquartile range of
        the observed relative error, whose error floor is the fifth
        percentile of that ratio, and whose dropout rates equal the
        observed missing fractions. Distortion and outlier parameters are
        left at zero; see :func:`fit_distortion_priors_from_sites` for
        those.

    Raises
    ------
    ValueError
        If ``survey`` has no declared ``impedance_error`` to profile.

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.ai.data.contracts import SurveyData
    >>> z = np.full((4, 6, 2), 100 + 50j)
    >>> err = np.full((4, 6, 2), 3.0)
    >>> survey = SurveyData(
    ...     z,
    ...     np.linspace(1000, 1, 6),
    ...     ["A", "B", "C", "D"],
    ...     ["xy", "yx"],
    ...     np.zeros((4, 2)),
    ...     impedance_error=err,
    ... )
    >>> config = fit_corruption_config(survey)
    >>> config.noise_level_range[0] >= 0.0
    True
    """
    if severity_scale <= 0.0 or not np.isfinite(severity_scale):
        raise ValueError("severity_scale must be finite and positive.")
    if survey.impedance_error is None:
        raise ValueError(
            "survey has no impedance_error to profile; supply a survey with "
            "declared errors, or build a CorruptionConfig from literature "
            "defaults instead."
        )
    valid = survey.valid
    ratio = np.abs(survey.impedance_error[valid]) / np.maximum(
        np.abs(survey.impedance[valid]), 1e-24
    )
    ratio = ratio[np.isfinite(ratio)]
    if ratio.size == 0:
        raise ValueError(
            "no valid observations with finite error ratio were found."
        )

    lo = float(np.percentile(ratio, 25)) * severity_scale
    hi = float(np.percentile(ratio, 75)) * severity_scale
    floor = float(np.percentile(ratio, 5)) * severity_scale
    lo, hi = sorted((lo, hi))

    coverage = survey.coverage()
    station_dropout = float(
        np.clip((1.0 - coverage.by_station).mean() * severity_scale, 0.0, 1.0)
    )
    frequency_dropout = float(
        np.clip(
            (1.0 - coverage.by_frequency).mean() * severity_scale, 0.0, 1.0
        )
    )
    random_dropout = float(
        np.clip((1.0 - coverage.overall) * severity_scale, 0.0, 1.0)
    )

    return CorruptionConfig(
        noise_level_range=(lo, hi),
        error_floor_fraction=floor,
        station_dropout_rate=station_dropout,
        frequency_dropout_rate=frequency_dropout,
        random_dropout_rate=random_dropout,
    )


def fit_distortion_priors_from_sites(
    sites: Any,
    *,
    recursive: bool = True,
    on_dup: str = "replace",
    verbose: int = 0,
    **kwargs: Any,
) -> dict[str, float]:
    """Estimate empirical static-shift and distortion spreads from real EDI.

    This is the one M3 entry point that genuinely depends on the heavier,
    pandas-based EM diagnostics in :mod:`pycsamt.emtools.gb` and
    :mod:`pycsamt.emtools.ss`, run directly on real sites (e.g. a WILLY
    line) rather than on the numpy-only :class:`SurveyData` contract.

    Parameters
    ----------
    sites : Any
        Anything accepted by :func:`pycsamt.emtools._core.ensure_sites`.
    recursive, on_dup, verbose
        Forwarded to the underlying diagnostics.
    **kwargs
        Forwarded to :func:`pycsamt.emtools.gb.groom_bailey_table`.

    Returns
    -------
    dict
        ``static_shift_log10_sigma``, ``distortion_gain_log10_sigma``,
        ``distortion_twist_deg_sigma``, ``distortion_shear_sigma``, and
        ``distortion_anisotropy_sigma``, each the population standard
        deviation of the corresponding per-station fitted parameter
        across stations with a successful fit. A parameter is ``0.0`` when
        fewer than two stations produced a usable fit.

    Examples
    --------
    >>> priors = fit_distortion_priors_from_sites(
    ...     "data/AMT/WILLY_DATA/L18PLT", recursive=False, verbose=0
    ... )  # doctest: +SKIP
    >>> sorted(priors)  # doctest: +SKIP
    ['distortion_anisotropy_sigma', 'distortion_gain_log10_sigma', \
'distortion_shear_sigma', 'distortion_twist_deg_sigma', \
'static_shift_log10_sigma']
    """

    zero = {
        "static_shift_log10_sigma": 0.0,
        "distortion_gain_log10_sigma": 0.0,
        "distortion_twist_deg_sigma": 0.0,
        "distortion_shear_sigma": 0.0,
        "distortion_anisotropy_sigma": 0.0,
    }

    gb_table = groom_bailey_table(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        verbose=verbose,
        api=False,
        **kwargs,
    )
    ok = gb_table[gb_table["status"] == "ok"] if len(gb_table) else gb_table
    if len(ok) >= 2:
        zero["distortion_gain_log10_sigma"] = float(
            np.std(np.log10(ok["gain"].to_numpy()))
        )
        zero["distortion_twist_deg_sigma"] = float(
            np.std(ok["twist_deg"].to_numpy())
        )
        zero["distortion_shear_sigma"] = float(np.std(ok["shear"].to_numpy()))
        zero["distortion_anisotropy_sigma"] = float(
            np.std(ok["anisotropy"].to_numpy())
        )

    ss_table = estimate_ss_ama(
        sites, recursive=recursive, on_dup=on_dup, verbose=verbose, api=False
    )
    if len(ss_table) >= 2:
        # delta_log10_rho is a shift on rho_a; Z scales as sqrt(rho_a).
        zero["static_shift_log10_sigma"] = 0.5 * float(
            np.std(ss_table["delta_log10_rho"].to_numpy())
        )

    return zero


def _sum_zero_contrast(n: int) -> np.ndarray:
    """Return the ``(n, n-1)`` sum-to-zero contrast matrix ``C``.

    ``effect = C @ theta`` satisfies ``sum(effect) == 0`` for every
    ``theta`` in ``R**(n-1)``, the standard "effect coding" trick for
    fitting a categorical main effect without an arbitrary reference
    level.
    """
    if n < 1:
        raise ValueError("n must be positive.")
    contrast = np.zeros((n, max(n - 1, 0)))
    if n > 1:
        contrast[: n - 1, :] = np.eye(n - 1)
        contrast[n - 1, :] = -1.0
    return contrast


@dataclass(frozen=True)
class StructuredMissingness:
    """Fitted station+frequency+component logistic missingness model.

    ``logit(pi_sjc) = alpha + a_s + b_j + q_c``, with the station,
    frequency, and component effects ``a``, ``b``, ``q`` each constrained
    to sum to zero for identifiability (an effect relative to its own
    group's mean, not to an arbitrary reference level). This replaces the
    three independent scalar dropout rates of
    :func:`fit_corruption_config`/:func:`~.simulator.apply_dropout` -- which
    are mathematically forced to be numerically identical whenever each is
    the mean of the same overall coverage deficit taken along a different
    axis -- with an explicit model of *where* missingness concentrates.
    Fit by :func:`fit_structured_missingness`.

    Parameters
    ----------
    alpha : float
        Intercept: logit of the fitted overall mean dropout probability.
    station_effect, frequency_effect, component_effect : tuple of float
        Per-station, per-frequency, and per-component effects, each
        summing to zero by construction.
    station_names, frequencies_hz, components : tuple
        Identities the effects are indexed by, kept for provenance. The
        fitted probability field is tied to this specific grid; it is not
        meant to be reused against a differently shaped survey.
    l2 : float, default=1.0
        L2 penalty weight used during fitting (not applied to ``alpha``).

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.ai.data.contracts import SurveyData
    >>> z = np.ones((3, 4, 1), dtype=complex)
    >>> z[0, :2] = np.nan
    >>> survey = SurveyData(
    ...     z, [4.0, 3.0, 2.0, 1.0], ["A", "B", "C"], ["xy"],
    ...     np.zeros((3, 2)),
    ... )
    >>> model = fit_structured_missingness(survey)
    >>> model.dropout_probability().shape
    (3, 4, 1)
    """

    alpha: float
    station_effect: tuple[float, ...]
    frequency_effect: tuple[float, ...]
    component_effect: tuple[float, ...]
    station_names: tuple[str, ...]
    frequencies_hz: tuple[float, ...]
    components: tuple[str, ...]
    l2: float = 1.0

    def dropout_probability(self) -> np.ndarray:
        """Return the fitted per-cell Bernoulli dropout probability.

        Returns
        -------
        ndarray, shape (n_station, n_frequency, n_component)
            ``sigmoid(alpha + a_s + b_j + q_c)``, ready to pass to
            :func:`~.simulator.apply_structured_dropout`.

        Examples
        --------
        >>> model = StructuredMissingness(
        ...     alpha=0.0, station_effect=(0.0,), frequency_effect=(0.0,),
        ...     component_effect=(0.0,), station_names=("A",),
        ...     frequencies_hz=(1.0,), components=("xy",),
        ... )
        >>> model.dropout_probability()
        array([[[0.5]]])
        """
        a = np.asarray(self.station_effect, dtype=float)[:, None, None]
        b = np.asarray(self.frequency_effect, dtype=float)[None, :, None]
        q = np.asarray(self.component_effect, dtype=float)[None, None, :]
        z = self.alpha + a + b + q
        return 1.0 / (1.0 + np.exp(-z))

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-serializable representation.

        Returns
        -------
        dict
            All fields as plain Python types, with a schema discriminator.

        Examples
        --------
        >>> model = StructuredMissingness(
        ...     alpha=0.0, station_effect=(0.0,), frequency_effect=(0.0,),
        ...     component_effect=(0.0,), station_names=("A",),
        ...     frequencies_hz=(1.0,), components=("xy",),
        ... )
        >>> model.to_dict()["schema_version"]
        1
        """
        return {
            "schema_version": 1,
            "alpha": self.alpha,
            "station_effect": list(self.station_effect),
            "frequency_effect": list(self.frequency_effect),
            "component_effect": list(self.component_effect),
            "station_names": list(self.station_names),
            "frequencies_hz": list(self.frequencies_hz),
            "components": list(self.components),
            "l2": self.l2,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> StructuredMissingness:
        """Restore a serialized model.

        Parameters
        ----------
        data : mapping
            State previously returned by :meth:`to_dict`.

        Returns
        -------
        StructuredMissingness
            Reconstructed, validated model.

        Examples
        --------
        >>> model = StructuredMissingness(
        ...     alpha=0.0, station_effect=(0.0,), frequency_effect=(0.0,),
        ...     component_effect=(0.0,), station_names=("A",),
        ...     frequencies_hz=(1.0,), components=("xy",),
        ... )
        >>> StructuredMissingness.from_dict(model.to_dict()) == model
        True
        """
        if data.get("schema_version") != 1:
            raise ValueError(
                "unsupported StructuredMissingness schema version."
            )
        return cls(
            alpha=float(data["alpha"]),
            station_effect=tuple(float(v) for v in data["station_effect"]),
            frequency_effect=tuple(
                float(v) for v in data["frequency_effect"]
            ),
            component_effect=tuple(
                float(v) for v in data["component_effect"]
            ),
            station_names=tuple(data["station_names"]),
            frequencies_hz=tuple(float(v) for v in data["frequencies_hz"]),
            components=tuple(data["components"]),
            l2=float(data["l2"]),
        )


def fit_structured_missingness(
    survey: SurveyData,
    *,
    l2: float = 1.0,
    max_iter: int = 500,
) -> StructuredMissingness:
    r"""Fit a station+frequency+component logistic missingness model.

    Replaces the three scalar dropout rates of :func:`fit_corruption_config`
    with an explicit model of where a survey's missingness is
    concentrated. With ``Y_sjc = 1 - survey.valid[s, j, c]`` the observed
    missingness indicator and ``pi_sjc = sigmoid(alpha + a_s + b_j + q_c)``,
    the fit solves the penalized Bernoulli maximum-likelihood problem

    .. math::

        \hat\vartheta = \arg\min_\vartheta \Bigl[
        -\sum_{sjc}\bigl\{Y_{sjc}\log\pi_{sjc}
        +(1-Y_{sjc})\log(1-\pi_{sjc})\bigr\}
        +\lambda\bigl(\lVert a\rVert_2^2+\lVert b\rVert_2^2
        +\lVert q\rVert_2^2\bigr)\Bigr],

    with ``a``, ``b``, ``q`` each parameterized through a sum-to-zero
    contrast (:func:`_sum_zero_contrast`) so the fit is identified without
    an arbitrary reference level, optimized by L-BFGS-B using the
    closed-form gradient of the penalized negative log-likelihood
    (``d(nll)/dz_sjc = pi_sjc - Y_sjc``, the standard logistic-regression
    result).

    Parameters
    ----------
    survey : SurveyData
        Survey whose ``valid`` mask defines the real missingness pattern
        to fit. Should carry a genuine, disclosed missing-data mask (e.g.
        a coverage-aware, snap-or-missing survey), not one already
        interpolated/filled -- interpolating over every gap removes the
        very signal this fit depends on.
    l2 : float, default=1.0
        Ridge penalty weight :math:`\lambda` on the station, frequency,
        and component effects (the intercept is not penalized).
    max_iter : int, default=500
        Maximum L-BFGS-B iterations.

    Returns
    -------
    StructuredMissingness
        Fitted model. Its
        :meth:`~StructuredMissingness.dropout_probability` gives the full
        ``(n_station, n_frequency, n_component)`` probability field for
        :func:`~.simulator.apply_structured_dropout`.

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.ai.data.contracts import SurveyData
    >>> z = np.ones((3, 4, 1), dtype=complex)
    >>> z[0, :2] = np.nan
    >>> survey = SurveyData(
    ...     z, [4.0, 3.0, 2.0, 1.0], ["A", "B", "C"], ["xy"],
    ...     np.zeros((3, 2)),
    ... )
    >>> model = fit_structured_missingness(survey)
    >>> abs(sum(model.station_effect)) < 1e-6
    True
    >>> abs(sum(model.frequency_effect)) < 1e-6
    True
    """
    from scipy.optimize import minimize

    if l2 < 0.0 or not np.isfinite(l2):
        raise ValueError("l2 must be finite and non-negative.")

    y = (~survey.valid).astype(float)
    n_station, n_frequency, n_component = survey.shape

    c_a = _sum_zero_contrast(n_station)
    c_b = _sum_zero_contrast(n_frequency)
    c_q = _sum_zero_contrast(n_component)
    n_a, n_b, n_q = c_a.shape[1], c_b.shape[1], c_q.shape[1]

    def unpack(theta: np.ndarray) -> tuple[float, np.ndarray, np.ndarray, np.ndarray]:
        alpha = theta[0]
        theta_a = theta[1 : 1 + n_a]
        theta_b = theta[1 + n_a : 1 + n_a + n_b]
        theta_q = theta[1 + n_a + n_b :]
        return alpha, theta_a, theta_b, theta_q

    def neg_log_posterior(theta: np.ndarray) -> tuple[float, np.ndarray]:
        alpha, theta_a, theta_b, theta_q = unpack(theta)
        a = c_a @ theta_a
        b = c_b @ theta_b
        q = c_q @ theta_q
        z = alpha + a[:, None, None] + b[None, :, None] + q[None, None, :]
        pi = 1.0 / (1.0 + np.exp(-z))
        eps = 1e-12
        nll = -np.sum(
            y * np.log(pi + eps) + (1.0 - y) * np.log(1.0 - pi + eps)
        )
        penalty = l2 * (np.sum(a**2) + np.sum(b**2) + np.sum(q**2))
        loss = nll + penalty

        g = pi - y
        d_alpha = np.sum(g)
        d_theta_a = c_a.T @ g.sum(axis=(1, 2)) + 2.0 * l2 * (c_a.T @ a)
        d_theta_b = c_b.T @ g.sum(axis=(0, 2)) + 2.0 * l2 * (c_b.T @ b)
        d_theta_q = c_q.T @ g.sum(axis=(0, 1)) + 2.0 * l2 * (c_q.T @ q)
        grad = np.concatenate(([d_alpha], d_theta_a, d_theta_b, d_theta_q))
        return float(loss), grad

    theta0 = np.zeros(1 + n_a + n_b + n_q)
    result = minimize(
        neg_log_posterior,
        theta0,
        jac=True,
        method="L-BFGS-B",
        options={"maxiter": max_iter},
    )
    alpha, theta_a, theta_b, theta_q = unpack(result.x)
    a = c_a @ theta_a
    b = c_b @ theta_b
    q = c_q @ theta_q

    return StructuredMissingness(
        alpha=float(alpha),
        station_effect=tuple(float(v) for v in a),
        frequency_effect=tuple(float(v) for v in b),
        component_effect=tuple(float(v) for v in q),
        station_names=tuple(survey.station_names),
        frequencies_hz=tuple(float(v) for v in survey.frequencies_hz),
        components=tuple(survey.components),
        l2=float(l2),
    )
