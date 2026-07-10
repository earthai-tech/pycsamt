# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""
pycsamt/log/_config.py

Initialize and configure logging for pycsamt
using the YAML config file `p.configlog.yml`.
"""
import os

from .logger import configure_logging


def init_logging(
    config_path: str = None,
    use_default: bool = False
) -> None:
    """
    Initialize pycsamt logging configuration.

    Parameters
    ----------
    config_path : str, optional
        Explicit path to the YAML logging configuration file.
        If ``None``, the code will look at the
        ``PYCSAMT_LOG_CONFIG`` environment variable; if that
        is unset, falls back to the bundled
        ``p.configlog.yml`` in this package directory.
    use_default : bool, default False
        If True, ignore the YAML file and configure
        logging using basicConfig.
    """
    # Determine which config file to load
    if not config_path or not config_path.strip():
        # allow override via environment
        env_path = os.environ.get("PYCSAMT_LOG_CONFIG")
        if env_path and env_path.strip():
            config_path = env_path
        else:
            # bundled default
            pkg_dir = os.path.dirname(os.path.abspath(__file__))
            config_path = os.path.join(pkg_dir, 'p.configlog.yml')

    configure_logging(
        config_path=config_path,
        use_default=use_default
    )


# Auto-configure on import
init_logging()
