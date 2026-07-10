"""
pycsamt/config.py

Utilities for managing the pycsamt data directory.
"""
import os
import shutil
import warnings

import yaml


def get_data_home(data_home: str = None) -> str:
    """
    Return the pycsamt_data directory, creating it if needed.

    Parameters
    ----------
    data_home : str, optional
        Explicit path to the desired data directory. If ``None``, the
        function checks the ``PYCSAMT_DATA`` environment variable,
        then falls back to ``~/pycsamt_data``. Tilde ('~') is expanded
        to the user's home directory.

    Returns
    -------
    data_home : str
        The absolute path to the pycsamt data directory.

    Notes
    -----
    - Directory is created if it does not exist.
    - Permission errors during creation will emit a warning.

    Examples
    --------
    >>> from pycsamt.config import get_data_home
    >>> # Use default location
    >>> path = get_data_home()
    >>> # Override with explicit path
    >>> path2 = get_data_home(
    ...     "D:/my_pycsamt_data"
    ... )
    """
    if data_home is None:
        data_home = os.environ.get(
            "PYCSAMT_DATA",
            os.path.join("~", "pycsamt_data")
        )
    data_home = os.path.expanduser(data_home)
    try:
        os.makedirs(data_home, exist_ok=True)
    except OSError as e:
        warnings.warn(
            f"Could not create pycsamt data home {data_home}: {e}", stacklevel=2
        )
    return data_home


def remove_data(data_home: str = None) -> None:
    """
    Remove the pycsamt data directory and its contents.

    Parameters
    ----------
    data_home : str, optional
        Path to the data directory to remove. If ``None``, the path
        is determined by :func:`get_data_home`.

    Returns
    -------
    None

    Notes
    -----
    - Deletes all contents under the data directory.
    - Use with caution: removal is permanent.

    Examples
    --------
    >>> from pycsamt.config import remove_data
    >>> # Remove default data directory
    >>> remove_data()
    """
    data_dir = get_data_home(data_home)
    if os.path.exists(data_dir):
        shutil.rmtree(data_dir)
    else:
        warnings.warn(
            f"pycsamt data directory not found: {data_dir}", stacklevel=2
        )


def get_config_home(config_home: str = None) -> str:
    """
    Return the pycsamt configuration directory, creating it.

    Parameters
    ----------
    config_home : str, optional
        Explicit path to the desired config directory. If
        ``None``, checks the ``PYCSAMT_CONFIG`` environment
        variable, then falls back to ``<data_home>/config``.
        Tilde ('~') is expanded.

    Returns
    -------
    config_home : str
        The absolute path to the config directory.

    Notes
    -----
    - Directory is created if it does not exist.
    - Use this to store YAML/JSON settings files.

    Examples
    --------
    >>> from pycsamt.config import get_config_home
    >>> cfg_dir = get_config_home()
    >>> # Custom config path
    >>> cfg2 = get_config_home("/etc/pycsamt/config")
    """
    if config_home is None:
        base = get_data_home()
        config_home = os.environ.get(
            "PYCSAMT_CONFIG",
            os.path.join(base, "config")
        )
    config_home = os.path.expanduser(config_home)
    try:
        os.makedirs(config_home, exist_ok=True)
    except OSError as e:
        warnings.warn(
            f"Could not create config home {config_home}: {e}", stacklevel=2
        )
    return config_home


def load_config(
    filename: str,
    config_dir: str = None
) -> dict:
    """
    Load a YAML configuration file from the config directory.

    Parameters
    ----------
    filename : str
        Name of the YAML file to load (e.g., "settings.yml").
    config_dir : str, optional
        Directory to look in. If ``None``, uses
        :func:`get_config_home`.

    Returns
    -------
    config : dict
        Parsed configuration mapping.

    Raises
    ------
    FileNotFoundError
        If the file does not exist in the config directory.
    yaml.YAMLError
        If the file cannot be parsed as valid YAML.

    Examples
    --------
    >>> from pycsamt.config import load_config
    >>> cfg = load_config("defaults.yml")
    >>> print(cfg.get("threshold"))
    """
    directory = config_dir or get_config_home()
    path = os.path.join(directory, filename)
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"Config file not found: {path}"
        )
    with open(path, encoding="utf-8") as f:
        return yaml.safe_load(f)
