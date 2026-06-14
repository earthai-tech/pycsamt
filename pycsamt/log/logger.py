# -*- coding: utf-8 -*-
# Author: LKouadio  <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Centralized logging configuration for pycsamt.
Supports loading from a YAML config file (`p.configlog.yml`) or falling back to
basic default logging configuration.
"""
import os
import logging
import logging.config
import yaml
import warnings 


def get_data_home(data_home: str = None) -> str:
    """
    Return the pycsamt_data directory, creating it if needed.
    Priority:
      1) explicit argument `data_home`
      2) env var PYCSAMT_DATA
      3) ~/pycsamt_data
    """
    if data_home is None:
        data_home = os.environ.get(
            "PYCSAMT_DATA", os.path.join("~", "pycsamt_data"))
    data_home = os.path.expanduser(data_home)
    try:
        os.makedirs(data_home, exist_ok=True)
    except OSError as e:
        warnings.warn(
            f"Could not create pycsamt data home {data_home}: {e}")
    return data_home

def configure_logging(
    config_path: str = None,
    use_default: bool = False
) -> None:
    """
    Configure logging, but if any handler has a *relative* filename,
    remap it into <data_home>/logs/<that-filename>.
    """
    if os.environ.get("PYCSAMT_DOCS_BUILD") == "1":
        logging.basicConfig(
            level=logging.WARNING,
            format="%(levelname)-8s [%(name)s] %(message)s",
        )
        return

    if use_default:
        logging.basicConfig(
            level=logging.INFO,
            format="%(asctime)s %(levelname)-8s [%(name)s] %(message)s",
            datefmt="%Y-%m-%dT %H:%M:%S %z",
        )
        logging.getLogger(__name__).info(
            "Logging configured with default/basicConfig")
        return

    # 1) find YAML...
    if not config_path:
        module_dir = os.path.dirname(os.path.abspath(__file__))
        config_path = os.path.join(module_dir, "p.configlog.yml")

    # 2) load it (or fallback)
    if os.path.exists(config_path):
        try:
            with open(config_path, "rt", encoding="utf-8") as f:
                cfg = yaml.safe_load(f)
        except Exception as e:
            logging.basicConfig(level=logging.INFO)
            logging.getLogger(__name__).exception(
                f"Failed to parse YAML at {config_path}, using basicConfig: {e}"
            )
            return
    else:
        logging.basicConfig(level=logging.INFO)
        logging.getLogger(__name__).warning(
            f"Logging config not found at {config_path}; using basicConfig"
        )
        return

    # 3) compute our data_home/logs
    data_home = get_data_home()
    logs_home = os.path.join(data_home, "logs")
    os.makedirs(logs_home, exist_ok=True)

    # 4) rewrite any relative filenames in handlers → point into logs_home
    handlers = cfg.get("handlers", {})
    for name, h in handlers.items():
        fn = h.get("filename")
        if fn:
            # if it’s already absolute, leave it; otherwise we join:
            if not os.path.isabs(fn):
                fn = os.path.join(logs_home, fn)
                h["filename"] = fn
            # ensure its directory exists
            os.makedirs(os.path.dirname(fn), exist_ok=True)

    # 5) apply the dictConfig
    logging.config.dictConfig(cfg)
    logging.getLogger(__name__).info(
        f"Loaded logging configuration from {config_path}")

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
