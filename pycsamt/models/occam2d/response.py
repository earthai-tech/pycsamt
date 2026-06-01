# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""OccamResponse — read Occam2D response (.resp) files.

A response file is produced by the Fortran binary after each iteration.
It contains seven whitespace-delimited columns per data point — no header:

::

    site  freq  type  err_floor  observed  modeled  residual
    1.0   1.0   5.0   0.0        4.8243    4.9036   -0.12094E-01
    ...

Columns
-------
0  site_index   (1-based float → int)
1  freq_index   (1-based float → int)
2  type_code    (int; same codes as OccamData)
3  error_floor  (float; often 0.0)
4  observed     (float; observed datum)
5  modeled      (float; forward-model prediction)
6  residual     (float; weighted residual)

Entry point
-----------
``OccamResponse.read(path)``
    Parse a ``.resp`` file and populate observed / predicted arrays.
"""

from __future__ import annotations

from pathlib import Path
from typing import Union

import numpy as np

from .base import OccamBase

PathLike = Union[str, Path]

__all__ = ["OccamResponse"]

_N_COLS = 7   # expected number of columns per data row


# -----------------------------------------------------------------------
# Low-level parser
# -----------------------------------------------------------------------

def _parse_response(path: Path) -> np.ndarray:
    """Return a (n_data, 7) float array from a ``.resp`` file.

    Raises
    ------
    FileNotFoundError
    ValueError
        If no valid data rows can be parsed.
    """
    if not path.exists():
        raise FileNotFoundError(f"Occam response file not found: {path}")

    rows: list[list[float]] = []
    with path.open("r", errors="replace") as fh:
        for raw in fh:
            tokens = raw.strip().split()
            if len(tokens) == _N_COLS:
                try:
                    rows.append([float(t) for t in tokens])
                except ValueError:
                    pass

    if not rows:
        raise ValueError(f"No valid data rows found in response file: {path}")

    return np.array(rows, dtype=float)


# -----------------------------------------------------------------------
# OccamResponse
# -----------------------------------------------------------------------

class OccamResponse(OccamBase):
    r"""Represent an Occam2D response file.

    ``OccamResponse`` stores the forward response written by
    the Occam2D executable for one inversion iteration. The
    response table has one row per datum and seven columns:
    site index, frequency index, type code, error-floor value,
    observed datum, modeled datum, and weighted residual.

    The response residual is already weighted by the data
    uncertainty used by Occam. The global RMS misfit is
    therefore computed directly from the residual column:

    .. math::

        \mathrm{RMS}
        =
        \sqrt{
            \frac{1}{N}
            \sum_{i=1}^{N} r_i^2
        },

    where :math:`r_i` is the weighted residual for datum
    :math:`i` and :math:`N` is the number of response rows.
    Values near one are often consistent with data errors that
    are neither under-estimated nor over-estimated [1]_.

    Parameters
    ----------
    verbose : int or bool, default 0
        Verbosity level inherited from :class:`OccamBase`.
        Positive values enable progress messages through the
        instance logger.
    logger : logging.Logger, optional
        Logger used for progress and diagnostic messages. If
        omitted, a class-specific PyCSAMT logger is created.

    Attributes
    ----------
    data : numpy.ndarray of float, shape (n_data, 7)
        Full raw table from the ``.resp`` file. Columns are
        ``site_index``, ``freq_index``, ``type_code``,
        ``error_floor``, ``observed``, ``modeled``, and
        ``residual``.
    observed : numpy.ndarray of float, shape (n_data,)
        Observed data values from column 4. Values follow the
        datum convention of ``OccamData``: log10 apparent
        resistivity for rho rows and degrees for phase rows.
    modeled : numpy.ndarray of float, shape (n_data,)
        Forward-model predictions from column 5, ordered in
        the same row order as ``observed``.
    residuals : numpy.ndarray of float, shape (n_data,)
        Weighted residuals from column 6. These values are the
        residuals used to compute :attr:`rms`.
    rms : float
        Root-mean-square weighted residual for all response
        rows. Empty objects use ``0.0``.

    Notes
    -----
    Response files do not include a header. The parser accepts
    any line with seven numeric columns and skips non-numeric
    lines. Site and frequency indices are one-based to match
    the Occam data file.

    See Also
    --------
    OccamData
        Defines the observed data rows and type codes.
    InversionResult
        Loads the response for a selected iteration.
    Plot2D.response
        Visualizes observed and modeled response curves.

    Examples
    --------
    Read a response file and inspect its global RMS:

    >>> from pycsamt.models.occam2d import OccamResponse
    >>> response = OccamResponse.read("occam_run/RESP17.resp")
    >>> response.rms

    Summarize misfit by station and frequency index:

    >>> response.misfit_per_site()
    >>> response.misfit_per_frequency()

    Select phase rows by Occam type code:

    >>> phase_tm = response.data[response.data[:, 2] == 6]
    >>> phase_tm.shape[0]

    References
    ----------
    .. [1] deGroot-Hedlin, C., and Constable, S.,
       "Occam's inversion to generate smooth, two-dimensional
       models from magnetotelluric data", Geophysics, 55(12),
       1613-1624, 1990.
    .. [2] Constable, S. C., Parker, R. L., and Constable,
       C. G., "Occam's inversion: A practical algorithm for
       generating smooth models from electromagnetic sounding
       data", Geophysics, 52(3), 289-300, 1987.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.data:      np.ndarray = np.empty((0, _N_COLS))
        self.observed:  np.ndarray = np.array([])
        self.modeled:   np.ndarray = np.array([])
        self.residuals: np.ndarray = np.array([])
        self.rms:       float      = 0.0

    # ------------------------------------------------------------------
    # I/O
    # ------------------------------------------------------------------
    @classmethod
    def read(
        cls,
        path: PathLike,
        data_fn: PathLike | None = None,
        **kwargs,
    ) -> OccamResponse:
        """Read an Occam2D response file.

        The reader parses seven-column numeric rows from a
        response file produced by the Occam2D executable. The
        first three columns are stored in the raw table as
        floats because the file itself is numeric text, but
        convenience properties expose site, frequency, and
        type codes as integers.

        Parameters
        ----------
        path : path-like
            Path to the response file. The value may be a
            string, :class:`pathlib.Path`, or any object
            accepted by :class:`pathlib.Path`.
        data_fn : path-like, optional
            Optional data-file path reserved for consistency
            checks between observed data rows and response
            rows. It is currently accepted for API stability
            but is not used by the parser.
        **kwargs
            Additional keyword arguments forwarded to the
            ``OccamResponse`` constructor. Use this for
            ``verbose`` or ``logger``.

        Returns
        -------
        OccamResponse
            Parsed response container with raw data, observed
            values, modeled values, residuals, and global RMS
            populated.

        Raises
        ------
        FileNotFoundError
            Raised when ``path`` does not exist.
        ValueError
            Raised when no valid seven-column numeric response
            rows can be parsed.

        See Also
        --------
        OccamResponse.misfit_per_site
            Computes station-index RMS values from residuals.
        OccamResponse.misfit_per_frequency
            Computes frequency-index RMS values.

        Examples
        --------
        >>> from pycsamt.models.occam2d import OccamResponse
        >>> response = OccamResponse.read("RESP17.resp")
        >>> response.n_data
        >>> response.type_codes
        """
        p   = Path(path)
        arr = _parse_response(p)
        obj = cls(**kwargs)

        obj.data      = arr
        obj.observed  = arr[:, 4]
        obj.modeled   = arr[:, 5]
        obj.residuals = arr[:, 6]
        obj.rms       = float(np.sqrt(np.mean(obj.residuals ** 2)))

        if obj.verbose:
            obj.logger.info(
                "OccamResponse.read: %d data points, rms=%.4f from %s",
                obj.n_data, obj.rms, p,
            )
        return obj

    # ------------------------------------------------------------------
    # Convenience
    # ------------------------------------------------------------------
    @property
    def n_data(self) -> int:
        return len(self.data)

    @property
    def site_indices(self) -> np.ndarray:
        """1-based site indices (int)."""
        return self.data[:, 0].astype(int) if self.data.size else np.array([], dtype=int)

    @property
    def freq_indices(self) -> np.ndarray:
        """1-based frequency indices (int)."""
        return self.data[:, 1].astype(int) if self.data.size else np.array([], dtype=int)

    @property
    def type_codes(self) -> np.ndarray:
        """Unique data-type codes present in this response."""
        if not self.data.size:
            return np.array([], dtype=int)
        return np.unique(self.data[:, 2].astype(int))

    # ------------------------------------------------------------------
    # Derived per-site / per-frequency misfit
    # ------------------------------------------------------------------
    def misfit_per_site(self) -> dict[int, float]:
        r"""Return RMS misfit for each site index.

        The returned values are computed from the weighted
        residual column:

        .. math::

            \mathrm{RMS}_s =
            \sqrt{
                \frac{1}{N_s}
                \sum_{i \in s} r_i^2
            }.

        Returns
        -------
        dict of int to float
            Mapping from one-based site index to RMS weighted
            residual. Empty response objects return an empty
            dictionary.

        Examples
        --------
        >>> from pycsamt.models.occam2d import OccamResponse
        >>> response = OccamResponse.read("RESP17.resp")
        >>> per_site = response.misfit_per_site()
        >>> per_site[1]
        """
        if not self.data.size:
            return {}
        result: dict[int, float] = {}
        for s in np.unique(self.site_indices):
            mask = self.site_indices == s
            result[int(s)] = float(np.sqrt(np.mean(self.residuals[mask] ** 2)))
        return result

    def misfit_per_frequency(self) -> dict[int, float]:
        r"""Return RMS misfit for each frequency index.

        The returned values group residuals by Occam's
        one-based frequency index:

        .. math::

            \mathrm{RMS}_f =
            \sqrt{
                \frac{1}{N_f}
                \sum_{i \in f} r_i^2
            }.

        Returns
        -------
        dict of int to float
            Mapping from one-based frequency index to RMS
            weighted residual. Empty response objects return
            an empty dictionary.

        Examples
        --------
        >>> from pycsamt.models.occam2d import OccamResponse
        >>> response = OccamResponse.read("RESP17.resp")
        >>> per_freq = response.misfit_per_frequency()
        >>> max(per_freq.values())
        """
        if not self.data.size:
            return {}
        result: dict[int, float] = {}
        for f in np.unique(self.freq_indices):
            mask = self.freq_indices == f
            result[int(f)] = float(np.sqrt(np.mean(self.residuals[mask] ** 2)))
        return result
