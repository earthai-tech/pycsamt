# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""CLI configuration API for pyCSAMT.

Mirrors the ``configure_* / reset_*`` pattern used throughout
``pycsamt.api``.  Import from here to configure CLI defaults
programmatically or from a script.

Example
-------
::

    from pycsamt.api.cli import configure_cli, reset_cli, PYCSAMT_CLI

    # set session-wide defaults
    configure_cli(log__level=1, output__format="json")

    # inspect current settings
    print(PYCSAMT_CLI)

    # restore defaults
    reset_cli()
"""

from .config import (
    BuildConfig,
    LogConfig,
    OutputConfig,
    PYCSAMT_CLI,
    PyCSAMTCLI,
    configure_cli,
    reset_cli,
)
from .options import (
    common_options,
    format_option,
    n_jobs_option,
    no_cache_option,
    no_color_option,
    output_dir_option,
    overwrite_option,
    verbose_option,
)
from .params import (
    EDIDir,
    EDIPath,
    FreqRange,
    StationList,
)

__all__ = [
    # config
    "BuildConfig",
    "LogConfig",
    "OutputConfig",
    "PYCSAMT_CLI",
    "PyCSAMTCLI",
    "configure_cli",
    "reset_cli",
    # options
    "common_options",
    "format_option",
    "n_jobs_option",
    "no_cache_option",
    "no_color_option",
    "output_dir_option",
    "overwrite_option",
    "verbose_option",
    # params
    "EDIDir",
    "EDIPath",
    "FreqRange",
    "StationList",
]
