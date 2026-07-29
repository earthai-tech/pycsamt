# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
stratagem.process
=================

Data-processing classes for Stratagem AMT surveys.

:class:`StaticShiftCorrector`
    Estimates and removes the static-shift effect from apparent-resistivity
    curves using the Adaptive Moving-Average (AMA) spatial filter
    (Kouadio et al., 2024; Torres-Verdín & Bostick, 1992).  Delegates to
    :func:`~pycsamt.emtools.ss.estimate_ss_ama` and
    :func:`~pycsamt.emtools.ss.apply_ss_factors`.

:class:`NoiseRemover`
    Multi-stage noise-removal pipeline: powerline notch filter →
    Hampel outlier filter → optional log-frequency smoothing.
    Delegates to :mod:`pycsamt.emtools.remove_noise`.

Both classes follow the ``fit() → out()`` pattern: ``fit()`` applies
corrections **in-place** on the Z arrays of the supplied EDIFile objects,
stores the result in :attr:`edi_objects_`, and returns ``self`` for
chaining.  ``out()`` either returns the processed objects or writes them
to a directory.

.. note::
   Because emtools processing functions modify impedance tensors in-place,
   passing the same EDIFile list through multiple processors sequentially
   is safe.  If you need to preserve the originals, pass ``copy=True`` to
   :meth:`fit`.
"""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path

import pandas as pd

from ..api.property import MetadataMixin, PyCSAMTObject
from ..emtools.remove_noise import (
    hampel_filter_freq,
    notch_powerline,
    smooth_logfreq,
)
from ..emtools.ss import apply_ss_factors, estimate_ss_ama
from ..exceptions import NotFittedError

__all__ = ["NoiseRemover", "StaticShiftCorrector"]


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def _write_edis(
    edi_objects: list,
    savepath: Path,
    overwrite: bool,
    verbose: int,
    tag: str,
) -> list[Path]:
    """Write a list of EDIFile objects to *savepath*."""
    savepath.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    for edi in edi_objects:
        fname = (
            edi.path.name
            if getattr(edi, "path", None) is not None
            else f"{edi.station or 'station'}.edi"
        )
        out_path = savepath / fname
        if out_path.exists() and not overwrite:
            written.append(out_path)
            continue
        try:
            edi.write(new_edifn=fname, savepath=str(savepath))
            written.append(out_path)
        except Exception as exc:
            if verbose:
                print(f"[{tag}] write failed {fname}: {exc}")
    if verbose:
        print(f"[{tag}] wrote {len(written)} files → {savepath}")
    return written


# ---------------------------------------------------------------------------
# StaticShiftCorrector
# ---------------------------------------------------------------------------


class StaticShiftCorrector(PyCSAMTObject, MetadataMixin):
    """Estimate and remove static-shift from Stratagem AMT impedance data.

    Implements the AMA (Adaptive Moving-Average) spatial filter to estimate
    per-station static-shift factors and correct the impedance tensor
    amplitudes accordingly.

    The correction is applied **in-place** on ``EDIFile.Z.z``.  Use
    ``copy=True`` in :meth:`fit` to preserve originals.

    Parameters
    ----------
    sort_by : {'lon', 'lat', 'name'}, default ``'lon'``
        Spatial ordering of stations for the AMA spatial average.
        Use ``'lon'`` for E-W profiles, ``'lat'`` for N-S profiles.
    half_window : int, default 3
        Number of neighbour stations on each side used in the AMA
        spatial average.
    weights : {'tri', 'gauss', 'uniform'}, default ``'tri'``
        Distance-weighting scheme for AMA neighbours.
    pband : tuple of (float, float), optional
        Period range ``(T_min, T_max)`` in seconds used when estimating
        the shift factor.  Useful for restricting the estimation to a
        band free of near-surface distortion.
    max_skew : float or None, default 6.0
        Phase-tensor skew threshold: stations with median |β| above this
        are excluded from the spatial average (strong 3-D distortion).
        Set to ``None`` to disable.
    verbose : int, default 0

    Attributes
    ----------
    factors_ : pandas.DataFrame
        Per-station shift factors with columns ``station``,
        ``delta_log10_rho``, ``fac_rho``, ``fac_z``, ``n_used``.
    edi_objects_ : list of EDIFile
        Corrected EDI objects (in-place modified unless ``copy=True``).

    Examples
    --------
    >>> from pycsamt.stratagem.process import StaticShiftCorrector
    >>> sc = StaticShiftCorrector(sort_by="lon", half_window=3).fit(edis)
    >>> sc.factors_.head()
    >>> paths = sc.out("2/2EDISS")
    """

    __repr_fields__ = (
        "sort_by",
        "half_window",
        "weights",
        "max_skew",
        "n_stations_",
    )

    def __init__(
        self,
        *,
        sort_by: str = "lon",
        half_window: int = 3,
        weights: str = "tri",
        pband: tuple | None = None,
        max_skew: float | None = 6.0,
        verbose: int = 0,
    ) -> None:
        self.sort_by = sort_by
        self.half_window = half_window
        self.weights = weights
        self.pband = pband
        self.max_skew = max_skew
        self.verbose = verbose

    # ------------------------------------------------------------------
    def fit(
        self,
        edi_objects: list,
        *,
        copy: bool = False,
    ) -> StaticShiftCorrector:
        """Estimate and apply static-shift corrections.

        Parameters
        ----------
        edi_objects : list of EDIFile
        copy : bool, default False
            Deep-copy Z data before correcting.

        Returns
        -------
        self
        """
        if copy:
            edi_objects = [deepcopy(e) for e in edi_objects]

        self.n_stations_ = len(edi_objects)

        # estimate shift factors (api=False forces plain DataFrame).
        # estimate_ss_ama calls build_phase_tensor_table unconditionally;
        # if Z data contains NaN (e.g. after FrequencyFilter), that can raise
        # LinAlgError.  We catch it and fall back to unit factors with a warning.
        try:
            self.factors_ = estimate_ss_ama(
                edi_objects,
                sort_by=self.sort_by,
                half_window=self.half_window,
                weights=self.weights,
                pband=self.pband,
                max_skew=self.max_skew,
                api=False,
                verbose=self.verbose,
            )
        except Exception as exc:
            print(
                f"[StaticShiftCorrector] WARNING — AMA estimation failed "
                f"({type(exc).__name__}: {exc}).\n"
                f"  Tip: run StaticShiftCorrector BEFORE FrequencyFilter so "
                f"that Z data is complete during spatial averaging.\n"
                f"  No static-shift correction applied."
            )
            stations = [edi.station or f"S{i:03d}" for i, edi in enumerate(edi_objects)]
            self.factors_ = pd.DataFrame(
                {
                    "station": stations,
                    "delta_log10_rho": 0.0,
                    "fac_rho": 1.0,
                    "fac_z": 1.0,
                    "n_used": 0,
                }
            )
        else:
            # apply correction in-place on Z only when estimation succeeded
            apply_ss_factors(
                edi_objects,
                self.factors_,
                key="fac_z",
                inplace=True,
                verbose=self.verbose,
            )

        self.edi_objects_ = edi_objects

        if self.verbose:
            med_fac = float(self.factors_["fac_z"].median())
            print(
                f"[StaticShiftCorrector] corrected {self.n_stations_} stations, "
                f"median fac_z={med_fac:.3f}"
            )

        return self

    # ------------------------------------------------------------------
    def out(
        self,
        savepath: str | Path | None = None,
        *,
        overwrite: bool = False,
    ) -> list | list[Path]:
        """Return corrected EDI objects or write them to *savepath*.

        Parameters
        ----------
        savepath : path-like, optional
            When ``None`` returns the list of :class:`~pycsamt.seg.edi.EDIFile`
            objects; otherwise writes files and returns a list of
            :class:`pathlib.Path`.
        overwrite : bool, default False

        Returns
        -------
        list of EDIFile or list of Path
        """
        if not hasattr(self, "edi_objects_"):
            raise NotFittedError("Call fit() first.")
        if savepath is None:
            return self.edi_objects_
        return _write_edis(
            self.edi_objects_,
            Path(savepath).expanduser().resolve(),
            overwrite,
            self.verbose,
            "StaticShiftCorrector",
        )


# ---------------------------------------------------------------------------
# NoiseRemover
# ---------------------------------------------------------------------------


class NoiseRemover(PyCSAMTObject):
    """Multi-stage noise-removal pipeline for Stratagem AMT data.

    Applies three sequential filters to the impedance tensor:

    1. **Powerline notch** — masks (interpolates) the mains frequency and
       its harmonics.  Controlled by ``mains_hz`` and ``n_harm``.
    2. **Hampel outlier filter** — identifies and replaces frequency-domain
       spike outliers using a median-absolute-deviation test.
    3. **Log-frequency smoothing** (optional) — applies a triangular or
       Gaussian kernel along the log-frequency axis.

    All corrections are applied **in-place** on ``EDIFile.Z.z``.

    Parameters
    ----------
    mains_hz : float, default 50.0
        Mains frequency (Hz).  Use 60.0 for North American data.
    n_harm : int, default 30
        Number of powerline harmonics to notch.
    tol_hz : float, default 0.08
        Frequency tolerance (Hz) around each harmonic for the notch
        filter.
    notch_mode : {'interp', 'zero', 'nan'}, default ``'interp'``
        How to handle notched bins: ``'interp'`` interpolates across them
        (recommended), ``'nan'`` flags them as missing.
    hampel_win : int, default 3
        Half-window size for the Hampel outlier filter (in frequency
        bins).
    hampel_nsig : float, default 3.0
        Outlier threshold in units of median absolute deviation.
    smooth : bool, default False
        Enable log-frequency smoothing (stage 3).
    smooth_win : int, default 3
        Smoothing half-window.  Values above 4 may trigger a known
        shape issue in :func:`~pycsamt.emtools.remove_noise.smooth_logfreq`
        for short frequency vectors; keep ≤ 3 unless you have verified
        your data.
    verbose : int, default 0

    Attributes
    ----------
    edi_objects_ : list of EDIFile
        Denoised EDI objects (in-place modified unless ``copy=True``).

    Examples
    --------
    >>> from pycsamt.stratagem.process import NoiseRemover
    >>> nr = NoiseRemover(mains_hz=50.0, smooth=True, smooth_win=3)
    >>> nr.fit(edis)
    >>> paths = nr.out("2/2EDID")
    """

    __repr_fields__ = (
        "mains_hz",
        "n_harm",
        "hampel_win",
        "smooth",
        "n_stations_",
    )

    def __init__(
        self,
        *,
        mains_hz: float = 50.0,
        n_harm: int = 30,
        tol_hz: float = 0.08,
        notch_mode: str = "interp",
        hampel_win: int = 3,
        hampel_nsig: float = 3.0,
        smooth: bool = False,
        smooth_win: int = 3,
        verbose: int = 0,
    ) -> None:
        self.mains_hz = mains_hz
        self.n_harm = n_harm
        self.tol_hz = tol_hz
        self.notch_mode = notch_mode
        self.hampel_win = hampel_win
        self.hampel_nsig = hampel_nsig
        self.smooth = smooth
        self.smooth_win = smooth_win
        self.verbose = verbose

    # ------------------------------------------------------------------
    def fit(
        self,
        edi_objects: list,
        *,
        copy: bool = False,
    ) -> NoiseRemover:
        """Apply the noise-removal pipeline.

        Parameters
        ----------
        edi_objects : list of EDIFile
        copy : bool, default False
            Deep-copy Z data before processing.

        Returns
        -------
        self
        """
        if copy:
            edi_objects = [deepcopy(e) for e in edi_objects]

        self.n_stations_ = len(edi_objects)

        # ── stage 1: powerline notch ──────────────────────────────────
        notch_powerline(
            edi_objects,
            mains_hz=self.mains_hz,
            n_harm=self.n_harm,
            tol_hz=self.tol_hz,
            mode=self.notch_mode,
            inplace=True,
            verbose=self.verbose,
        )

        # ── stage 2: Hampel outlier filter ───────────────────────────
        hampel_filter_freq(
            edi_objects,
            win=self.hampel_win,
            nsig=self.hampel_nsig,
            inplace=True,
            verbose=self.verbose,
        )

        # ── stage 3: log-frequency smoothing (optional) ───────────────
        if self.smooth:
            try:
                smooth_logfreq(
                    edi_objects,
                    win=self.smooth_win,
                    inplace=True,
                    verbose=self.verbose,
                )
            except (ValueError, IndexError) as exc:
                if self.verbose:
                    print(
                        f"[NoiseRemover] smooth_logfreq skipped (win="
                        f"{self.smooth_win}): {exc}"
                    )

        self.edi_objects_ = edi_objects

        if self.verbose:
            print(
                f"[NoiseRemover] processed {self.n_stations_} stations "
                f"(notch + hampel" + (" + smooth)" if self.smooth else ")")
            )

        return self

    # ------------------------------------------------------------------
    def out(
        self,
        savepath: str | Path | None = None,
        *,
        overwrite: bool = False,
    ) -> list | list[Path]:
        """Return denoised EDI objects or write them to *savepath*.

        Parameters
        ----------
        savepath : path-like, optional
        overwrite : bool, default False

        Returns
        -------
        list of EDIFile or list of Path
        """
        if not hasattr(self, "edi_objects_"):
            raise NotFittedError("Call fit() first.")
        if savepath is None:
            return self.edi_objects_
        return _write_edis(
            self.edi_objects_,
            Path(savepath).expanduser().resolve(),
            overwrite,
            self.verbose,
            "NoiseRemover",
        )
