# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""Format-neutral transfer-function reading and writing entry points."""

from __future__ import annotations

from typing import Any

from .formats import (
    TransferFunctionFormatError,
    detect_tf_format,
    get_tf_format,
    get_tf_format_for_target,
)

__all__ = ["read_transfer_function", "write_transfer_function"]


def _buffer_nonseekable(source: Any) -> Any:
    """Buffer a non-seekable stream so detection cannot consume it."""
    read = getattr(source, "read", None)
    if not callable(read):
        return source
    seekable = getattr(source, "seekable", None)
    if callable(seekable):
        try:
            if seekable():
                return source
        except Exception:
            pass
    seek = getattr(source, "seek", None)
    tell = getattr(source, "tell", None)
    if callable(seek) and callable(tell):
        try:
            pos = tell()
            seek(pos)
            return source
        except Exception:
            pass
    return read()


def read_transfer_function(
    source: Any,
    *,
    format: str | None = None,
    **kwargs: Any,
) -> Any:
    """Read a supported electromagnetic transfer-function source.

    Parameters
    ----------
    source : path-like, str, bytes, or readable file object
        Input transfer-function source. Historical EDI currently requires a
        filesystem path because :class:`pycsamt.seg.EDIFile` owns that parser.
        EMTF XML additionally accepts inline XML, bytes, and file-like input.
    format : str, optional
        Explicit registered format name or alias. When omitted, content-based
        detection is used.
    **kwargs
        Passed to the selected backend reader.

    Returns
    -------
    object
        The established backend scientific object for the detected format.
        EDI returns ``EDIFile`` and EMTF XML returns ``EMTF``. Explicit
        conversion between them is available through ``EMTF.from_edi()`` and
        ``EMTF.to_edi()``.
    """
    prepared = _buffer_nonseekable(source)
    if format is None:
        canonical = detect_tf_format(prepared)
        spec = get_tf_format(canonical)
    else:
        spec = get_tf_format(format)
    return spec.reader(prepared, **kwargs)


def write_transfer_function(
    obj: Any,
    target: Any,
    *,
    format: str | None = None,
    **kwargs: Any,
) -> Any:
    """Write a transfer-function object through a registered serializer.

    Parameters
    ----------
    obj : object
        Scientific/backend object accepted by the selected writer. EMTF XML
        accepts :class:`pycsamt.emtf.EMTF`; EDI accepts ``EMTF`` or the
        established :class:`pycsamt.seg.EDIFile`.
    target : path-like or writable stream
        Destination. A path ending in ``.xml`` can select EMTF XML
        automatically. Streams require ``format=`` because no extension is
        available before serialization.
    format : str, optional
        Explicit registered format name or alias.
    **kwargs
        Passed to the selected backend writer.

    Notes
    -----
    EMTF -> EDI may be information-losing because standard EDI cannot carry
    every EMTF covariance, metadata, or arbitrary transfer-function type.
    The EDI writer exposes an ``on_loss`` policy through ``**kwargs``.
    """
    spec = (
        get_tf_format(format)
        if format
        else get_tf_format_for_target(target)
    )
    if spec.writer is None:
        raise TransferFunctionFormatError(
            f"transfer-function format {spec.name!r} is not writable yet"
        )
    return spec.writer(obj, target, **kwargs)
