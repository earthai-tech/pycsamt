"""Unified readers, writers, and configuration-file I/O."""

from .config import Config
from .formats import (
    TransferFunctionFormat,
    TransferFunctionFormatError,
    detect_tf_format,
    get_tf_format,
    get_tf_format_for_target,
    list_tf_formats,
    register_tf_format,
)
from .parsers import read_any, write_any
from .transfer import read_transfer_function, write_transfer_function

read_tf = read_transfer_function
write_tf = write_transfer_function

__all__ = [
    "Config",
    "read_any",
    "write_any",
    "TransferFunctionFormat",
    "TransferFunctionFormatError",
    "register_tf_format",
    "get_tf_format",
    "get_tf_format_for_target",
    "list_tf_formats",
    "detect_tf_format",
    "read_transfer_function",
    "write_transfer_function",
    "read_tf",
    "write_tf",
]
