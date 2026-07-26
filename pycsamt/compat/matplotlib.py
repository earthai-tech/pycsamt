# File: pycsamt/compat/matplotlib.py
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later

"""
A compatibility module to handle API changes across different
versions of Matplotlib.

Notably, ``matplotlib.cm.get_cmap`` was deprecated in Matplotlib 3.7 and
removed in 3.11. :func:`get_cmap` here is a version-agnostic replacement
that also accepts the old ``lut`` argument (number of discrete colours),
so ``matplotlib.cm.get_cmap(name, N)`` call sites can migrate unchanged.
"""

import warnings

import matplotlib
from packaging.version import parse

# Get the installed Matplotlib version
_MPL_VERSION = parse(matplotlib.__version__)

__all__ = ["get_cmap", "is_valid_cmap"]


def is_valid_cmap(cmap, allow_none=False, **kw):
    r"""Check if a colormap identifier is valid.

    This function purely validates whether a given identifier can be
    resolved to a Matplotlib colormap. It does not retrieve the
    colormap object itself.

    Parameters
    ----------
    cmap : any
        The colormap identifier to validate, typically a string.
    allow_none : bool, optional
        If True, `None` is considered a valid input and the
        function will return True. If False (default), `None` is
        considered invalid.

    Returns
    -------
    bool
        True if the `cmap` is a valid, retrievable colormap name
        or if `cmap` is None and `allow_none` is True. Otherwise,
        returns False.
    """
    is_valid = False
    if cmap is None:
        return allow_none

    if not isinstance(cmap, str):
        return False

    try:
        # The most reliable way to check for existence is to try
        # getting it, using the modern API first.
        if _MPL_VERSION >= parse("3.6"):
            is_valid = matplotlib.colormaps.get(cmap)
        else:
            is_valid = matplotlib.cm.get_cmap(cmap)

        if is_valid is None:
            return False

        return True
    except (ValueError, KeyError):
        # ValueError is for old MPL, KeyError for new MPL.
        return False


def get_cmap(
    name,
    default="viridis",
    allow_none=False,
    error=None,
    failsafe="continuous",
    lut=None,
    **kw,
):
    r"""Robustly retrieve a Matplotlib colormap with fallbacks.

    This function ensures a valid colormap object is always returned,
    preventing runtime errors from invalid names. It uses a
    cascading fallback system.

    Parameters
    ----------
    name : str or None
        The desired colormap name.
    default : str, optional
        The fallback colormap if `name` is invalid.
        Defaults to 'viridis'.
    allow_none : bool, optional
        If True, a `name` of `None` will return `None` without
        any warnings or errors. Defaults to False.
    failsafe : {'continuous', 'discrete'}, optional
        Specifies the type of ultimate fallback colormap to use if
        both `name` and `default` are invalid.
        - 'continuous': Use 'viridis' (default).
        - 'discrete': Use 'tab10'.
    lut : int, optional
        Number of discrete colours to resample the colormap to — the
        drop-in for the removed ``matplotlib.cm.get_cmap(name, lut)``
        second argument. ``None`` (default) returns the full colormap.

    Returns
    -------
    matplotlib.colors.Colormap or None
        A valid colormap instance, or `None` if `allow_none` is
        True and the input `name` is `None`.
    """
    result = None
    # For API consistency, acknowledge the old 'error' parameter
    # but inform the user that it's no longer used.
    if error is not None:
        warnings.warn(
            "The 'error' parameter is deprecated for get_cmap and is ignored. "
            "This function now always returns a valid colormap by using "
            "fallbacks.",
            FutureWarning,
            stacklevel=2,
        )

        # but does nothing

    # Private helper to prevent repeating the retrieval code
    def _retrieve(cmap_name):
        """Retrieve the colormap object using the correct API, resampled
        to ``lut`` discrete colours when requested."""
        if _MPL_VERSION >= parse("3.6"):
            cmap = matplotlib.colormaps.get(cmap_name)
            if cmap is not None and lut is not None:
                resample = getattr(cmap, "resampled", None)
                if resample is None:  # Matplotlib 3.5 registry objects
                    resample = cmap._resample
                cmap = resample(int(lut))
            return cmap
        # Legacy (<3.6): cm.get_cmap accepts the lut directly.
        if lut is not None:
            return matplotlib.cm.get_cmap(cmap_name, int(lut))
        return matplotlib.cm.get_cmap(cmap_name)

    # 1. Handle explicit `None` input first.
    if name is None:
        if allow_none:
            return None
        # If None is not allowed, treat it as an invalid name
        # and proceed to the fallback logic below.
    # 2. Try to validate and retrieve the primary name.
    elif is_valid_cmap(name):
        result = _retrieve(name)

    if result is not None:
        return result
    # 3. If we are here, 'name' was invalid. Warn and fall back to default.
    warnings.warn(
        f"Colormap '{name}' not found. Falling back to default '{default}'.",
        UserWarning,
        stacklevel=2,
    )
    if is_valid_cmap(default):
        result = _retrieve(default)

    if result is not None:
        return result

    # apply failure safe here
    # 4. If the default is also invalid, warn and use the ultimate failsafe.
    # 4. If default is also invalid, determine and use the ultimate failsafe.
    if failsafe == "discrete":
        failsafe_cmap = "tab10"
    else:
        failsafe_cmap = "viridis"
        if failsafe != "continuous":
            warnings.warn(
                f"Invalid `failsafe` value '{failsafe}'. Defaulting to "
                f"'continuous' type ('{failsafe_cmap}').",
                UserWarning,
                stacklevel=2,
            )

    warnings.warn(
        f"Default colormap '{default}' also not found. "
        f"Falling back to failsafe '{failsafe_cmap}'.",
        UserWarning,
        stacklevel=2,
    )
    return _retrieve("viridis")


def _install_legacy_cm_get_cmap():
    """Restore the removed ``matplotlib.cm.get_cmap`` compatibility hook."""

    if not hasattr(matplotlib.cm, "get_cmap"):
        # Optional GUI dependencies may still resolve the historical public
        # name. Our adapter already dispatches through the modern registry.
        matplotlib.cm.get_cmap = get_cmap


_install_legacy_cm_get_cmap()
