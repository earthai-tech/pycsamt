# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""
pycsamt/log/_config.py

Initialize and configure logging for pycsamt
using the YAML config file `p.configlog.yml`.
"""
import os

from .logger import configure_logging


def init_logging(config_path: str = None, use_default: bool = False) -> None:
    """
    Initialize pycsamt logging configuration.

    Parameters
    ----------
    config_path : str, optional
        Path to the YAML logging configuration file. If not provided, defaults
        to `p.configlog.yml` in this directory.
    use_default : bool, default False
        If True, ignore the YAML file and configure logging using basicConfig.
    """
    # Determine config file location
    if not config_path or not config_path.strip():
        base_dir = os.path.dirname(os.path.abspath(__file__))
        config_path = os.path.join(base_dir, 'p.configlog.yml')

    # Delegate to logger module
    configure_logging(config_path=config_path, use_default=use_default)


# Auto-configure on import
init_logging()
