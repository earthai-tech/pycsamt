# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""Occam2D model-prejudice file support.

The sparse ``OCCAM2MTPREJ_2.0`` format assigns a preferred
log10-resistivity and penalty weight to selected Occam model
parameters. Each record contains three values::

    parameter_index  encoded_prejudice  prejudice_weight

Parameter indices are one-based. The public
:class:`OccamPrejudice` API exposes physical target values rather
than solver-encoded values. This distinction is required because
the bundled Occam2D solver adds ``prewt**2`` to the penalty Hessian
but adds ``prewt*premod`` to the right-hand side. The writer uses

.. math::

    p_j = w_j m_j^{\mathrm{target}},

where :math:`p_j` is the native prejudice value and :math:`w_j` is
the native weight. This keeps the quadratic penalty centred on the
requested target as the weight changes.

Entry points
------------
``OccamPrejudice.from_dense(target_values, weights)``
    Build a sparse prejudice object from dense model vectors.
``OccamPrejudice.read(path)``
    Read and decode an ``OCCAM2MTPREJ_2.0`` file.
``OccamPrejudice.write(path)``
    Validate and write a solver-native prejudice file.
"""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path
from typing import Union

import numpy as np

from ...compat.sklearn import validate_params
from .base import OccamBase
from .config import OccamConfig
from .schema import (
    OCCAM_PREJUDICE_DENSE_SCHEMA,
    OCCAM_PREJUDICE_INIT_SCHEMA,
    OCCAM_PREJUDICE_PARAMETER_COUNT_SCHEMA,
    OCCAM_PREJUDICE_READ_SCHEMA,
    OCCAM_PREJUDICE_WRITE_SCHEMA,
)

PathLike = Union[str, Path]

__all__ = ["OccamPrejudice"]

_FORMAT_TAG = "OCCAM2MTPREJ_2.0"
_KEY_WIDTH = 18


def _as_array(values: Iterable, dtype) -> np.ndarray:
    """Convert a possibly lazy iterable to a flat NumPy array."""
    if isinstance(values, np.ndarray):
        return np.asarray(values, dtype=dtype).reshape(-1)
    return np.asarray(list(values), dtype=dtype).reshape(-1)


class OccamPrejudice(OccamBase):
    r"""Represent a sparse Occam2D model-prejudice file.

    ``OccamPrejudice`` stores selected Occam parameter indices,
    their preferred log10-resistivity values, and non-negative
    penalty weights. The object follows the same container and I/O
    conventions as :class:`OccamData`, :class:`OccamModel`, and
    :class:`OccamStartup`.

    For a target :math:`m_j^{\mathrm{target}}` and weight :math:`w_j`,
    the intended local penalty is proportional to

    .. math::

        w_j^2
        \left(m_j-m_j^{\mathrm{target}}\right)^2.

    The bundled solver assembles its prejudice terms using
    ``prewt**2`` in the Hessian and ``prewt*premod`` in the
    right-hand side. Consequently, :meth:`write` encodes the native
    prejudice field as ``target_values * weights``. :meth:`read`
    reverses that encoding so users always work with physical target
    values.

    Parameters
    ----------
    parameter_indices : iterable of int, optional
        One-based Occam model-parameter indices. Values must be
        positive and unique. If omitted, an empty prejudice object
        is created.
    target_values : iterable of float, optional
        Preferred model values in log10 resistivity, ordered like
        ``parameter_indices``. Values must be finite. If omitted,
        an empty array is used.
    weights : iterable of float, optional
        Non-negative prejudice weights, ordered like
        ``parameter_indices``. Values must be finite. A zero weight
        leaves the associated target inactive.
    config : OccamConfig, optional
        Configuration object associated with the Occam2D project.
        The prejudice format has no configuration fields of its own,
        but retaining the project configuration matches the other
        Occam2D file containers. If omitted, a default configuration
        is created.
    verbose : int or bool, default 0
        Verbosity level inherited from :class:`OccamBase`. Positive
        values enable progress messages through the instance logger.
    logger : logging.Logger, optional
        Logger used for progress and diagnostic messages. If omitted,
        a class-specific PyCSAMT logger is created.

    Attributes
    ----------
    format_str : str
        File-format identifier. The supported value is
        ``"OCCAM2MTPREJ_2.0"``.
    parameter_indices : numpy.ndarray of int, shape (n_prejudiced,)
        One-based indices of prejudiced Occam model parameters.
    target_values : numpy.ndarray of float, shape (n_prejudiced,)
        Decoded physical targets in log10 resistivity.
    weights : numpy.ndarray of float, shape (n_prejudiced,)
        Non-negative native penalty weights.
    config : OccamConfig
        Occam2D project configuration associated with the object.
    path : pathlib.Path or None
        Most recent path read or written. Inherited from
        :class:`OccamBase`.

    Notes
    -----
    The public targets are not the raw second column stored in the
    native file when a weight differs from one. Directly editing that
    column without applying the encoding can shift the effective
    penalty centre.

    Zero-weight records are accepted for controlled experiments and
    round-trip fidelity. :meth:`from_dense` omits them by default to
    keep production prejudice files sparse.

    See Also
    --------
    OccamModel
        References a prejudice file through ``prejudice_file``.
    OccamStartup
        Stores the model vector to which prejudice penalties apply.
    InputBuilder
        Creates the remaining Occam2D input files.

    Examples
    --------
    Create and write a sparse prejudice definition:

    >>> from pycsamt.models.occam2d import OccamPrejudice
    >>> prejudice = OccamPrejudice(
    ...     parameter_indices=[2, 5],
    ...     target_values=[1.5, 2.3],
    ...     weights=[2.0, 0.5],
    ... )
    >>> prejudice.write("occam_run/DUHIPrejudice")

    Build sparse records from dense model vectors:

    >>> prejudice = OccamPrejudice.from_dense(
    ...     target_values=[1.5, 2.0, 2.5],
    ...     weights=[4.0, 0.0, 1.0],
    ... )
    >>> prejudice.parameter_indices.tolist()
    [1, 3]

    Read a file and inspect its decoded target values:

    >>> restored = OccamPrejudice.read("occam_run/DUHIPrejudice")
    >>> restored.n_prejudiced
    2

    References
    ----------
    .. [OccamPrejudice-1] deGroot-Hedlin, C., and Constable, S.,
       "Occam's inversion to generate smooth, two-dimensional
       models from magnetotelluric data", Geophysics, 55(12),
       1613-1624, 1990.
    .. [OccamPrejudice-2] Constable, S. C., Parker, R. L., and
       Constable, C. G., "Occam's inversion: A practical algorithm
       for generating smooth models from electromagnetic sounding
       data", Geophysics, 52(3), 289-300, 1987.
    """

    @validate_params(OCCAM_PREJUDICE_INIT_SCHEMA)
    def __init__(
        self,
        parameter_indices: Iterable[int] | None = None,
        target_values: Iterable[float] | None = None,
        weights: Iterable[float] | None = None,
        config: OccamConfig | None = None,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.format_str: str = _FORMAT_TAG
        self.config: OccamConfig = config or OccamConfig()
        self.parameter_indices: np.ndarray = _as_array(
            [] if parameter_indices is None else parameter_indices,
            int,
        )
        self.target_values: np.ndarray = _as_array(
            [] if target_values is None else target_values,
            float,
        )
        self.weights: np.ndarray = _as_array(
            [] if weights is None else weights,
            float,
        )
        self.validate()

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------
    @classmethod
    @validate_params(OCCAM_PREJUDICE_DENSE_SCHEMA)
    def from_dense(
        cls,
        target_values: Iterable[float],
        weights: Iterable[float],
        *,
        include_zero_weight: bool = False,
        config: OccamConfig | None = None,
        **kwargs,
    ) -> OccamPrejudice:
        """Build sparse prejudice records from dense model vectors.

        Dense input vectors are interpreted in Occam parameter order.
        Their zero-based array positions are converted to one-based
        solver indices. By default, entries with zero weight are
        omitted from the sparse result.

        Parameters
        ----------
        target_values : iterable of float
            Preferred log10-resistivity value for every Occam model
            parameter.
        weights : iterable of float
            Non-negative prejudice weight for every Occam model
            parameter. The vector must have the same length as
            ``target_values``.
        include_zero_weight : bool, default False
            If ``True``, preserve inactive zero-weight records.
            Otherwise, omit them from the sparse object.
        config : OccamConfig, optional
            Occam2D project configuration attached to the returned
            object. If omitted, a default configuration is created.
        **kwargs
            Additional keyword arguments forwarded to the
            ``OccamPrejudice`` constructor. Use this for ``verbose``
            or ``logger``.

        Returns
        -------
        OccamPrejudice
            Sparse prejudice object in one-based Occam parameter
            order.

        Raises
        ------
        ValueError
            Raised when the dense vectors have different lengths or
            contain invalid target or weight values.

        See Also
        --------
        OccamPrejudice.write
            Encodes and writes the sparse result.
        OccamPrejudice.validate_parameter_count
            Checks the indices against a model parameter count.

        Examples
        --------
        >>> from pycsamt.models.occam2d import OccamPrejudice
        >>> prejudice = OccamPrejudice.from_dense(
        ...     [1.5, 2.0, 2.5],
        ...     [4.0, 0.0, 1.0],
        ... )
        >>> prejudice.parameter_indices.tolist()
        [1, 3]
        """
        targets = _as_array(target_values, float)
        penalties = _as_array(weights, float)
        if targets.size != penalties.size:
            raise ValueError(
                "target_values and weights must have equal length"
            )

        keep = np.ones(targets.size, dtype=bool)
        if not include_zero_weight:
            keep = penalties > 0

        return cls(
            parameter_indices=np.arange(1, targets.size + 1, dtype=int)[keep],
            target_values=targets[keep],
            weights=penalties[keep],
            config=config,
            **kwargs,
        )

    # ------------------------------------------------------------------
    # Validation
    # ------------------------------------------------------------------
    def validate(self) -> OccamPrejudice:
        """Validate the current sparse prejudice records.

        Validation checks vector lengths, one-based unique indices,
        finite targets, and finite non-negative weights. The method
        does not require at least one record; empty prejudice objects
        are valid containers.

        Returns
        -------
        OccamPrejudice
            The validated object. Returning ``self`` supports fluent
            preparation workflows.

        Raises
        ------
        ValueError
            Raised when the arrays have different lengths, an index
            is non-positive or duplicated, a target is non-finite, or
            a weight is non-finite or negative.

        See Also
        --------
        OccamPrejudice.validate_parameter_count
            Performs model-size validation after record validation.
        OccamPrejudice.write
            Calls this method before serialization.

        Examples
        --------
        >>> prejudice = OccamPrejudice([1], [2.0], [1.0])
        >>> prejudice.validate() is prejudice
        True
        """
        sizes = {
            self.parameter_indices.size,
            self.target_values.size,
            self.weights.size,
        }
        if len(sizes) != 1:
            raise ValueError(
                "parameter_indices, target_values, and weights "
                "must have equal length"
            )
        if np.any(self.parameter_indices < 1):
            raise ValueError(
                "Occam parameter indices must be one-based and positive"
            )
        if np.unique(self.parameter_indices).size != self.n_prejudiced:
            raise ValueError("parameter_indices must be unique")
        if not np.all(np.isfinite(self.target_values)):
            raise ValueError("target_values must be finite")
        invalid_weights = not np.all(np.isfinite(self.weights)) or np.any(
            self.weights < 0
        )
        if invalid_weights:
            raise ValueError("weights must be finite and non-negative")
        if self.format_str != _FORMAT_TAG:
            raise ValueError(
                f"unsupported prejudice format: {self.format_str!r}"
            )
        return self

    @validate_params(OCCAM_PREJUDICE_PARAMETER_COUNT_SCHEMA)
    def validate_parameter_count(self, n_params: int) -> OccamPrejudice:
        """Validate prejudice indices against an Occam model size.

        Parameters
        ----------
        n_params : int
            Total number of parameters declared by the Occam model
            and startup files. It must be positive.

        Returns
        -------
        OccamPrejudice
            The validated object.

        Raises
        ------
        TypeError
            Raised when ``n_params`` is not an integer.
        ValueError
            Raised when ``n_params`` is not positive or a prejudice
            index exceeds it.

        See Also
        --------
        OccamModel.n_params
            Supplies the expected model parameter count.
        OccamPrejudice.validate
            Checks record-level constraints first.

        Examples
        --------
        >>> prejudice = OccamPrejudice([1, 3], [2.0, 2.5], [1, 1])
        >>> prejudice.validate_parameter_count(3) is prejudice
        True
        """
        self.validate()
        if not isinstance(n_params, (int, np.integer)):
            raise TypeError("n_params must be an integer")
        if n_params < 1:
            raise ValueError("n_params must be positive")
        if self.n_prejudiced:
            largest = int(self.parameter_indices.max())
            if largest > int(n_params):
                raise ValueError(
                    "prejudice references a parameter beyond n_params"
                )
        return self

    # ------------------------------------------------------------------
    # I/O
    # ------------------------------------------------------------------
    @classmethod
    @validate_params(OCCAM_PREJUDICE_READ_SCHEMA)
    def read(
        cls,
        path: PathLike,
        config: OccamConfig | None = None,
        **kwargs,
    ) -> OccamPrejudice:
        """Read and decode an ``OCCAM2MTPREJ_2.0`` file.

        The native second column is divided by its positive weight to
        recover the public physical target. For a zero-weight record,
        the raw value is retained because the record has no numerical
        influence and cannot be uniquely decoded.

        Parameters
        ----------
        path : path-like
            Path to the Occam prejudice file. The value may be a
            string, :class:`pathlib.Path`, or any object accepted by
            :class:`pathlib.Path`.
        config : OccamConfig, optional
            Occam2D project configuration attached to the returned
            object. If omitted, a default configuration is created.
        **kwargs
            Additional keyword arguments forwarded to the
            ``OccamPrejudice`` constructor. Use this for ``verbose``
            or ``logger``.

        Returns
        -------
        OccamPrejudice
            Parsed container with decoded target values and ``path``
            set to the source file.

        Raises
        ------
        FileNotFoundError
            Raised when ``path`` does not exist.
        ValueError
            Raised when the header is invalid, the declared record
            count does not match the file, a record does not have
            three columns, or parsed values fail validation.

        See Also
        --------
        OccamPrejudice.write
            Applies the inverse encoding during serialization.
        OccamPrejudice.from_dense
            Creates sparse records from dense vectors.

        Examples
        --------
        >>> from pycsamt.models.occam2d import OccamPrejudice
        >>> prejudice = OccamPrejudice.read("DUHIPrejudice")
        >>> prejudice.path.name
        'DUHIPrejudice'
        """
        source = Path(path)
        if not source.exists():
            raise FileNotFoundError(
                f"Occam prejudice file not found: {source}"
            )

        lines = source.read_text(
            encoding="ascii", errors="replace"
        ).splitlines()
        if len(lines) < 2:
            raise ValueError(f"Incomplete Occam prejudice file: {source}")

        format_value = cls._header_value(lines[0], "Format")
        if format_value != _FORMAT_TAG:
            raise ValueError(
                f"Expected format '{_FORMAT_TAG}', got "
                f"'{format_value}' in {source}"
            )

        count_value = cls._header_value(lines[1], "Param Count")
        try:
            count = int(count_value)
        except ValueError as exc:
            raise ValueError(
                f"Invalid prejudice parameter count in {source}"
            ) from exc
        if count < 0:
            raise ValueError(f"Negative prejudice parameter count in {source}")

        rows = [line.split() for line in lines[2:] if line.strip()]
        if len(rows) != count or any(len(row) != 3 for row in rows):
            raise ValueError(
                "prejudice record count or column layout does not "
                f"match the header in {source}"
            )

        try:
            indices = np.asarray([int(row[0]) for row in rows])
            encoded = np.asarray([float(row[1]) for row in rows])
            weights = np.asarray([float(row[2]) for row in rows])
        except ValueError as exc:
            raise ValueError(
                f"Invalid numeric prejudice record in {source}"
            ) from exc

        targets = encoded.copy()
        positive = weights > 0
        targets[positive] = encoded[positive] / weights[positive]
        obj = cls(
            parameter_indices=indices,
            target_values=targets,
            weights=weights,
            config=config,
            **kwargs,
        )
        obj.path = source

        if obj.verbose:
            obj.logger.info(
                "OccamPrejudice.read: %d records from %s",
                obj.n_prejudiced,
                source,
            )
        return obj

    @validate_params(OCCAM_PREJUDICE_WRITE_SCHEMA)
    def write(self, path: PathLike) -> Path:
        """Write this object in ``OCCAM2MTPREJ_2.0`` format.

        The writer validates the current records, converts each public
        target to the solver-native value ``target * weight``, creates
        parent directories, and stores the destination on :attr:`path`.
        A zero-weight record retains its public target in the native
        column for round-trip fidelity; the solver ignores that record.

        Parameters
        ----------
        path : path-like
            Destination path for the prejudice file. The value may be
            a string, :class:`pathlib.Path`, or any object accepted by
            :class:`pathlib.Path`.

        Returns
        -------
        pathlib.Path
            Path to the file that was written.

        Raises
        ------
        ValueError
            Raised when the current records fail validation.

        See Also
        --------
        OccamPrejudice.read
            Reads and decodes files written by this method.
        OccamModel.prejudice_file
            References the written file from the Occam model.

        Examples
        --------
        >>> prejudice = OccamPrejudice([1], [1.5], [2.0])
        >>> written = prejudice.write("run/DUHIPrejudice")
        >>> written.name
        'DUHIPrejudice'
        """
        self.validate()
        output = Path(path)
        output.parent.mkdir(parents=True, exist_ok=True)
        lines = [
            self._header_line("Format", self.format_str),
            self._header_line("Param Count", self.n_prejudiced),
        ]
        for index, encoded, weight in zip(
            self.parameter_indices,
            self.native_values,
            self.weights,
        ):
            lines.append(f"{int(index):d} {encoded:.12g} {weight:.12g}\n")

        output.write_text("".join(lines), encoding="ascii")
        self.path = output
        if self.verbose:
            self.logger.info(
                "OccamPrejudice.write: %d records to %s",
                self.n_prejudiced,
                output,
            )
        return output

    # ------------------------------------------------------------------
    # Convenience
    # ------------------------------------------------------------------
    @property
    def n_prejudiced(self) -> int:
        """Return the number of sparse prejudice records.

        Returns
        -------
        int
            Number of entries in ``parameter_indices``, equivalent to
            the native ``Param Count`` header value.

        Examples
        --------
        >>> OccamPrejudice([1, 3], [2.0, 2.5], [1, 1]).n_prejudiced
        2
        """
        return int(self.parameter_indices.size)

    @property
    def native_values(self) -> np.ndarray:
        """Return solver-encoded prejudice values without writing.

        Positive-weight targets are multiplied by their weights. A
        zero-weight entry retains its target value for round-trip
        fidelity. The returned array is newly allocated and modifying
        it does not change this object.

        Returns
        -------
        numpy.ndarray of float, shape (n_prejudiced,)
            Values written in the second native file column.

        See Also
        --------
        OccamPrejudice.target_values
            Decoded physical target values.
        OccamPrejudice.write
            Serializes these encoded values.

        Examples
        --------
        >>> prejudice = OccamPrejudice([1], [1.5], [2.0])
        >>> prejudice.native_values.tolist()
        [3.0]
        """
        encoded = self.target_values.copy()
        positive = self.weights > 0
        encoded[positive] *= self.weights[positive]
        return encoded

    # ------------------------------------------------------------------
    # Formatting helpers
    # ------------------------------------------------------------------
    @staticmethod
    def _header_line(keyword: str, value) -> str:
        """Return one fixed-width Occam prejudice header line."""
        return f"{f'{keyword}:':<{_KEY_WIDTH}}{value}\n"

    @staticmethod
    def _header_value(line: str, keyword: str) -> str:
        """Return and validate one Occam prejudice header value."""
        if ":" not in line:
            raise ValueError(f"Missing '{keyword}:' prejudice header")
        raw_key, raw_value = line.split(":", 1)
        if raw_key.strip().upper() != keyword.upper():
            raise ValueError(
                f"Expected '{keyword}:' prejudice header, got "
                f"'{raw_key.strip()}:'"
            )
        return raw_value.strip()
