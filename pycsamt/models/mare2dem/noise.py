# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Synthetic Gaussian noise addition for MARE2DEM response data.

Port of:
  * ``m2d_addSyntheticNoise.m``
  * ``m2d_makeSyntheticData.m``

Workflow for synthetic inversion studies:

1. Run a forward model with MARE2DEM to obtain predicted responses in
   an ``.EMResp`` file (8-column DATA block).
2. Call :func:`make_synthetic_data` to add noise and write a clean
   ``.emdata`` data file suitable for inversion.
3. Invert the synthetic data to test resolution.

Noise-level specification
-------------------------
Pass a :class:`NoiseConfig` object describing relative and absolute
noise levels.  The rules match the original MATLAB code exactly:

* **log10 apparent resistivity** (codes 123, 125, 129, 27–29, 37–39):
  ``std = 0.4343 × rel_noise``  (propagated from linear to log10 scale).
* **Apparent resistivity** (codes 103, 105, 109, 21–25, 31–35):
  ``std = rel_noise × |data|``.
* **Phase** (codes 104, 106, 110, 22–26, 32–36):
  ``std = rel_noise × 180 / π``.
* **Tipper Real/Imag** (codes 133, 134):
  ``std = rel_noise × amplitude``; NaN for amplitude < abs_noise_floor.
"""

from __future__ import annotations

import copy
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .iotools.emdata import (
    EMDataFile,
    read_emdata,
    write_emdata,
)

__all__ = ["NoiseConfig", "add_synthetic_noise", "make_synthetic_data"]

# ---------------------------------------------------------------------------
# Noise-level configuration
# ---------------------------------------------------------------------------


@dataclass
class NoiseConfig:
    """Per-data-type noise specifications.

    Attributes
    ----------
    mt_rel_noise : float
        Relative noise for MT apparent resistivity (e.g. 0.05 = 5 %).
        Phase noise is ``mt_rel_noise / 2``.
    mt_rel_noise_tipper : float
        Relative noise for tipper real/imaginary components.
    mt_abs_noise_tipper : float
        Absolute noise floor for tipper.  Tipper data with amplitude
        below this value are dropped (NaN).
    csem_rel_noise_e : float
        Relative noise for CSEM electric-field responses.
    csem_rel_noise_b : float
        Relative noise for CSEM magnetic-field responses.
    max_sigma_factor : float, default 2.0
        Clip noise at this many standard deviations to avoid
        unrealistically large random draws (matches MATLAB ``abs(randn) > 2``
        clipping).
    """

    mt_rel_noise: float = 0.05
    mt_rel_noise_tipper: float = 0.01
    mt_abs_noise_tipper: float = 0.01
    csem_rel_noise_e: float = 0.05
    csem_rel_noise_b: float = 0.05
    max_sigma_factor: float = 2.0


# ---------------------------------------------------------------------------
# Low-level noise primitives (matching MATLAB sub-functions exactly)
# ---------------------------------------------------------------------------


def _apply_noise(
    data: np.ndarray,
    std: np.ndarray,
    *,
    rng: np.random.Generator,
    clip: float = 2.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Add clipped Gaussian noise and return (noisy_data, std_out)."""
    noise = rng.standard_normal(data.shape)
    too_large = np.abs(noise) > clip
    noise[too_large] = np.sign(noise[too_large]) * clip
    return data + std * noise, std


def _logamp_noise(
    log_amp: np.ndarray,
    rel_error: float,
    rng: np.random.Generator,
    clip: float,
) -> tuple[np.ndarray, np.ndarray]:
    std = np.full_like(log_amp, 0.4343 * rel_error)
    return _apply_noise(log_amp, std, rng=rng, clip=clip)


def _amp_noise(
    amp: np.ndarray,
    rel_error: float,
    rng: np.random.Generator,
    clip: float,
) -> tuple[np.ndarray, np.ndarray]:
    std = rel_error * np.abs(amp)
    return _apply_noise(amp, std, rng=rng, clip=clip)


def _phase_noise(
    phase: np.ndarray,
    rel_error: float,
    rng: np.random.Generator,
    clip: float,
) -> tuple[np.ndarray, np.ndarray]:
    std = np.full_like(phase, rel_error * 180.0 / np.pi)
    return _apply_noise(phase, std, rng=rng, clip=clip)


def _tipper_real_imag_noise(
    real: np.ndarray,
    imag: np.ndarray,
    rel_error: float,
    abs_error: float,
    rng: np.random.Generator,
    clip: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    amp = np.sqrt(real**2 + imag**2)
    std = rel_error * amp
    std[amp < abs_error] = np.nan
    r_noisy, r_std = _apply_noise(real, std, rng=rng, clip=clip)
    i_noisy, i_std = _apply_noise(imag, std, rng=rng, clip=clip)
    return r_noisy, r_std, i_noisy, i_std


# ---------------------------------------------------------------------------
# Main API
# ---------------------------------------------------------------------------


def add_synthetic_noise(
    em: EMDataFile,
    noise: NoiseConfig,
    *,
    seed: int | None = None,
) -> EMDataFile:
    """Add synthetic Gaussian noise to MARE2DEM response data.

    Port of ``m2d_addSyntheticNoise.m``.

    The function reads predicted responses from column 7 of the 8-column
    DATA block (``is_response=True``), adds noise according to *noise*,
    and stores the result in columns 5 and 6 (data, std_err).  The
    output DATA block has 6 columns (observed-data format).

    Parameters
    ----------
    em : EMDataFile
        Response file with 8-column DATA block
        (``Format: EMResp_2.3``).
    noise : NoiseConfig
        Noise specification.
    seed : int or None, default None
        Random-number generator seed for reproducibility.

    Returns
    -------
    EMDataFile
        New ``EMDataFile`` with 6-column DATA block containing the
        noisy synthetic data and standard errors.

    Raises
    ------
    ValueError
        When ``em.data`` does not have 8 columns (not a response file).

    Examples
    --------
    >>> from pycsamt.models.mare2dem import read_emdata
    >>> from pycsamt.models.mare2dem.noise import (
    ...     NoiseConfig,
    ...     add_synthetic_noise,
    ... )
    >>> resp = read_emdata("forward.EMResp")
    >>> nc = NoiseConfig(mt_rel_noise=0.05)
    >>> noisy = add_synthetic_noise(resp, nc, seed=42)
    >>> noisy.n_data
    200
    """
    if em.data is None or em.data.shape[1] != 8:
        raise ValueError(
            "add_synthetic_noise requires an 8-column response file "
            "(EMResp format).  The data column is column index 6."
        )

    rng = np.random.default_rng(seed)
    clip = noise.max_sigma_factor

    data = em.data.copy()
    # Initialise data and std columns to NaN so that any code not handled
    # below is automatically dropped by the valid-mask at the end.
    data[:, 4] = np.nan
    data[:, 5] = np.nan
    codes = data[:, 0].astype(int)
    resp = em.data[:, 6].copy()  # predicted response column (from original)

    # ---- MT ----
    # log10 ApRes (TE, TM, det)
    mask = np.isin(codes, [123, 125, 129])
    if mask.any():
        d, s = _logamp_noise(resp[mask], noise.mt_rel_noise, rng, clip)
        data[mask, 4] = d
        data[mask, 5] = s

    # linear ApRes
    mask = np.isin(codes, [103, 105, 109])
    if mask.any():
        d, s = _amp_noise(resp[mask], noise.mt_rel_noise, rng, clip)
        data[mask, 4] = d
        data[mask, 5] = s

    # Phase (MT)
    mask = np.isin(codes, [104, 106, 110])
    if mask.any():
        d, s = _phase_noise(resp[mask], noise.mt_rel_noise / 2, rng, clip)
        data[mask, 4] = d
        data[mask, 5] = s

    # Tipper Real/Imag
    mask_re = codes == 133
    mask_im = codes == 134
    if mask_re.any() and mask_im.any():
        r, rs, i, is_ = _tipper_real_imag_noise(
            resp[mask_re],
            resp[mask_im],
            noise.mt_rel_noise_tipper,
            noise.mt_abs_noise_tipper,
            rng,
            clip,
        )
        data[mask_re, 4] = r
        data[mask_re, 5] = rs
        data[mask_im, 4] = i
        data[mask_im, 5] = is_

    # ---- CSEM ----
    # log10 amplitude E
    mask = np.isin(codes, [27, 28, 29])
    if mask.any():
        d, s = _logamp_noise(resp[mask], noise.csem_rel_noise_e, rng, clip)
        data[mask, 4] = d
        data[mask, 5] = s

    # log10 amplitude B
    mask = np.isin(codes, [37, 38, 39])
    if mask.any():
        d, s = _logamp_noise(resp[mask], noise.csem_rel_noise_b, rng, clip)
        data[mask, 4] = d
        data[mask, 5] = s

    # linear amplitude E
    mask = np.isin(codes, [21, 23, 25])
    if mask.any():
        d, s = _amp_noise(resp[mask], noise.csem_rel_noise_e, rng, clip)
        data[mask, 4] = d
        data[mask, 5] = s

    # linear amplitude B
    mask = np.isin(codes, [31, 33, 35])
    if mask.any():
        d, s = _amp_noise(resp[mask], noise.csem_rel_noise_b, rng, clip)
        data[mask, 4] = d
        data[mask, 5] = s

    # phase E
    mask = np.isin(codes, [22, 24, 26])
    if mask.any():
        d, s = _phase_noise(resp[mask], noise.csem_rel_noise_e, rng, clip)
        data[mask, 4] = d
        data[mask, 5] = s

    # phase B
    mask = np.isin(codes, [32, 34, 36])
    if mask.any():
        d, s = _phase_noise(resp[mask], noise.csem_rel_noise_b, rng, clip)
        data[mask, 4] = d
        data[mask, 5] = s

    # drop NaN rows (below tipper noise floor)
    valid = ~np.isnan(data[:, 4])
    data = data[valid, :6]  # keep only first 6 columns → observed-data format

    result = copy.copy(em)
    result.data = data
    result.is_response = False
    result.format = "EMData_2.3"
    return result


def make_synthetic_data(
    in_file: str | Path,
    out_file: str | Path,
    noise: NoiseConfig,
    *,
    seed: int | None = None,
) -> EMDataFile:
    """Read a MARE2DEM response file, add noise, and write a data file.

    Port of ``m2d_makeSyntheticData.m``.

    Parameters
    ----------
    in_file : path-like
        MARE2DEM response file (``.EMResp``) with 8-column DATA block.
    out_file : path-like
        Output ``.emdata`` file with 6-column noisy DATA block.
    noise : NoiseConfig
        Noise specification.
    seed : int or None, default None
        RNG seed for reproducibility.

    Returns
    -------
    EMDataFile
        The noisy data object (also written to *out_file*).

    Examples
    --------
    >>> from pycsamt.models.mare2dem.noise import (
    ...     NoiseConfig,
    ...     make_synthetic_data,
    ... )
    >>> nc = NoiseConfig(mt_rel_noise=0.05, csem_rel_noise_e=0.05)
    >>> em = make_synthetic_data(
    ...     "forward.EMResp", "synthetic.emdata", nc, seed=0
    ... )
    """
    resp = read_emdata(in_file)
    noisy = add_synthetic_noise(resp, noise, seed=seed)
    write_emdata(noisy, out_file)
    return noisy
