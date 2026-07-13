# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later

"""
This module provides foundational base classes and configuration
objects for the Zonge subpackage.
"""

from __future__ import annotations

import copy
import inspect
from typing import Any

from ..log.logger import get_logger

__all__ = ["Zonge"]


class Zonge:
    """
    A foundational base class for Zonge data objects.

    This class provides common attributes and methods, such as
    verbosity control, a standardized logger, and parameter
    management, that can be inherited by other classes in the
    Zonge subpackage.

    Parameters
    ----------
    verbose : bool, default False
        Controls the level of detail in logging output. If ``True``,
        informational and debugging messages will be displayed.
    **kws : dict
        Additional keyword arguments for future compatibility.
    """

    def __init__(self, verbose: bool = False, **kws):
        self.verbose = verbose
        self._logger = get_logger(self.__class__.__name__)

    def get_params(self, deep: bool = True) -> dict[str, Any]:
        """
        Get parameters for this object.

        This method retrieves the parameters that were set in the
        object's ``__init__`` method.

        Parameters
        ----------
        deep : bool, default True
            If True, will recursively return the parameters for
            contained sub-objects.

        Returns
        -------
        dict
            A dictionary of parameter names mapped to their values.
        """
        out = dict()
        for key in self._get_param_names():
            value = getattr(self, key)
            if deep and hasattr(value, "get_params"):
                deep_items = value.get_params().items()
                out.update((key + "__" + k, v) for k, v in deep_items)
            out[key] = value
        return out

    def set_params(self, **params) -> Zonge:
        """
        Set the parameters of this object.

        This method allows for updating the object's parameters
        after it has been instantiated.

        Parameters
        ----------
        **params : dict
            A dictionary of parameter names and their new values.

        Returns
        -------
        self : object
            The instance with updated parameters.
        """
        if not params:
            return self
        valid_params = self.get_params(deep=True)
        for key, value in params.items():
            if key not in valid_params:
                raise ValueError(
                    f"Invalid parameter {key!r} for estimator "
                    f"{self}. Check the list of available "
                    f"parameters with `get_params().keys()`."
                )
            setattr(self, key, value)
        return self

    @classmethod
    def _get_param_names(cls) -> list[str]:
        """Get parameter names for the estimator."""
        init_signature = inspect.signature(cls.__init__)
        parameters = [
            p.name
            for p in init_signature.parameters.values()
            if p.name != "self" and p.kind != p.VAR_KEYWORD
        ]
        return sorted(parameters)

    def __copy__(self):
        """Create a shallow copy of the object."""
        cls = self.__class__
        result = cls.__new__(cls)
        result.__dict__.update(self.__dict__)
        return result

    def __deepcopy__(self, memo):
        """Create a deep copy of the object."""
        cls = self.__class__
        result = cls.__new__(cls)
        memo[id(self)] = result
        for k, v in self.__dict__.items():
            setattr(result, k, copy.deepcopy(v, memo))

        return result

    def __str__(self) -> str:
        """Provide a concise, human-readable representation."""
        attrs = {
            k: v for k, v in self.__dict__.items() if not k.startswith("_")
        }
        attrs_str = ", ".join(f"{k}={v!r}" for k, v in attrs.items())
        return f"<{self.__class__.__name__}({attrs_str})>"

    def __repr__(self) -> str:
        """Provide an unambiguous developer representation."""
        return f"{self.__class__.__name__}(verbose={self.verbose})"
