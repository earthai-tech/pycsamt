# -*- coding: utf-8 -*-
# Author: LKouadio  <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt/log/logger.py

Centralized logging configuration for pycsamt.
Supports loading from a YAML config file (`p.configlog.yml`) or falling back to
basic default logging configuration.
"""
import os
import logging
import logging.config
import yaml


def configure_logging(
    config_path: str = None,
    use_default: bool = False
) -> None:
    """
    Configure logging for the pycsamt package.

    Parameters
    ----------
    config_path : str, optional
        Path to a YAML logging configuration file. If not provided, the function
        looks for `p.configlog.yml` in the same directory as this module.
    use_default : bool, default False
        If True, ignore any config file and use Python's basicConfig with INFO level.
    """
    if use_default:
        # Basic console-only configuration
        logging.basicConfig(
            level=logging.INFO,
            format="%(asctime)s %(levelname)-8s [%(name)s] %(message)s",
            datefmt="%Y-%m-%dT %H:%M:%S %z",
        )
        logging.getLogger(__name__).info("Logging configured with default/basicConfig")
        return

    # Determine YAML path
    if config_path is None or not config_path.strip():
        module_dir = os.path.dirname(os.path.abspath(__file__))
        config_path = os.path.join(module_dir, 'p.configlog.yml')

    # Attempt to load YAML config
    if os.path.exists(config_path):
        try:
            with open(config_path, 'rt', encoding='utf-8') as f:
                cfg = yaml.safe_load(f)

            # Ensure any file-handler directories exist
            handlers = cfg.get('handlers', {})
            for handler in handlers.values():
                filename = handler.get('filename')
                if filename:
                    log_dir = os.path.dirname(filename)
                    if log_dir and not os.path.exists(log_dir):
                        os.makedirs(log_dir, exist_ok=True)

            # Apply configuration
            logging.config.dictConfig(cfg)
            logging.getLogger(__name__).info(
                f"Loaded logging configuration from {config_path}" )
        except Exception as e:
            # Fallback to basicConfig on error
            logging.basicConfig(
                level=logging.INFO,
                format="%(asctime)s %(levelname)-8s [%(name)s] %(message)s",
                datefmt="%Y-%m-%dT %H:%M:%S %z",
            )
            logging.getLogger(__name__).exception(
                f"Failed to load YAML config; using basicConfig: {e}" )
    else:
        # Config file not found
        logging.basicConfig(
            level=logging.INFO,
            format="%(asctime)s %(levelname)-8s [%(name)s] %(message)s",
            datefmt="%Y-%m-%dT %H:%M:%S %z",
        )
        logging.getLogger(__name__).warning(
            f"Logging config file not found at {config_path}; using basicConfig" )


def get_logger(
    name: str = None
) -> logging.Logger:
    """
    Retrieve a named logger after ensuring logging is configured.

    Parameters
    ----------
    name : str, optional
        Name of the logger; `__name__` by default.

    Returns
    -------
    logging.Logger
    """
    logger_name = name or __name__
    return logging.getLogger(logger_name)


# Auto-configure on import (optional)
try:
    configure_logging()
except Exception:
    # If anything goes wrong, fallback to basic
    logging.basicConfig(level=logging.INFO)
    logging.getLogger(__name__).exception(
        "Unexpected error during logging setup; using basicConfig." )
