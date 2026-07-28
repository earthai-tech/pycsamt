# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Model containers used by :mod:`pycsamt.inversion`."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from ..api.property import MetadataMixin, PyCSAMTObject
from .doc import _inversion_param_docs

__all__ = ["StartingModel", "ReferenceModel"]


@dataclass
class StartingModel(PyCSAMTObject, MetadataMixin):
    resistivities: Any
    thicknesses: Any
    name: str = ""
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.resistivities = np.asarray(self.resistivities, dtype=float)
        self.thicknesses = np.asarray(self.thicknesses, dtype=float)
        self.validate()

    @classmethod
    def default(cls, n_layers: int = 4) -> StartingModel:
        """Return a conservative layered starting model.

        Parameters
        ----------
        n_layers : int, default 4
            Number of layers, including the bottom halfspace.

        Returns
        -------
        StartingModel
            Default model with 100 ohm-m layers and geometrically increasing
            thicknesses between 100 m and 2000 m.

        Raises
        ------
        ValueError
            If ``n_layers < 2``.

        Examples
        --------
        >>> from pycsamt.inversion.model import StartingModel
        >>> start = StartingModel.default(n_layers=3)
        >>> start.resistivities.tolist()
        [100.0, 100.0, 100.0]
        >>> start.thicknesses.size
        2
        """
        n_layers = int(n_layers)
        if n_layers < 2:
            raise ValueError("n_layers must be >= 2.")
        resistivities = np.full(n_layers, 100.0, dtype=float)
        thicknesses = np.geomspace(100.0, 2000.0, n_layers - 1)
        return cls(resistivities, thicknesses, name="default")

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> StartingModel:
        """Build from singular or plural mapping keys.

        Parameters
        ----------
        data : mapping
            Model mapping. Accepted keys are ``resistivity`` or
            ``resistivities`` for layer resistivity, ``thickness`` or
            ``thicknesses`` for layer thickness, plus optional ``name`` and
            ``metadata``.

        Returns
        -------
        StartingModel
            Validated layered-earth model.

        Examples
        --------
        >>> from pycsamt.inversion.model import StartingModel
        >>> start = StartingModel.from_dict(
        ...     {
        ...         "resistivity": [100.0, 500.0],
        ...         "thickness": [250.0],
        ...         "name": "two_layer",
        ...     }
        ... )
        >>> start.name
        'two_layer'
        """
        return cls(
            data.get("resistivity", data.get("resistivities")),
            data.get("thickness", data.get("thicknesses")),
            name=str(data.get("name", "")),
            metadata=dict(data.get("metadata", {})),
        )

    @classmethod
    def coerce(cls, value: Any, *, n_layers: int = 4) -> StartingModel:
        """Return *value* as a :class:`StartingModel`.

        Parameters
        ----------
        value : StartingModel, mapping, object, or None
            Input model. Existing ``StartingModel`` instances are returned
            unchanged. Mappings are read by :meth:`from_dict`. Generic objects
            may expose ``resistivity``/``resistivities`` and
            ``thickness``/``thicknesses`` attributes.
        n_layers : int, default 4
            Number of layers used only when *value* is ``None`` and a default
            model must be created.

        Returns
        -------
        StartingModel
            Validated layered-earth model.

        Examples
        --------
        >>> from pycsamt.inversion.model import StartingModel
        >>> StartingModel.coerce(None, n_layers=2).n_layers
        2
        >>> StartingModel.coerce(
        ...     {
        ...         "resistivities": [80.0, 250.0],
        ...         "thicknesses": [500.0],
        ...     }
        ... ).depths.tolist()
        [0.0, 500.0]
        """
        if value is None:
            return cls.default(n_layers=n_layers)
        if isinstance(value, cls):
            return value
        if isinstance(value, dict):
            return cls.from_dict(value)
        resistivity = getattr(value, "resistivity", None)
        if resistivity is None:
            resistivity = getattr(value, "resistivities", None)
        thickness = getattr(value, "thickness", None)
        if thickness is None:
            thickness = getattr(value, "thicknesses", None)
        return cls(resistivity, thickness, name=getattr(value, "name", ""))

    @property
    def n_layers(self) -> int:
        return int(self.resistivities.size)

    @property
    def depths(self) -> np.ndarray:
        """Top-of-layer depths in metres.

        Examples
        --------
        >>> from pycsamt.inversion.model import StartingModel
        >>> StartingModel([100.0, 300.0, 800.0], [50.0, 150.0]).depths.tolist()
        [0.0, 50.0, 200.0]
        """
        return np.r_[0.0, np.cumsum(self.thicknesses)]

    def to_layered_model(self):
        """Return the existing :class:`pycsamt.forward.LayeredModel`.

        This adapter lets inversion starting/recovered models reuse the forward
        modelling API without duplicating layer containers.

        Returns
        -------
        pycsamt.forward.LayeredModel
            Forward-model-compatible layered earth.

        Examples
        --------
        >>> from pycsamt.inversion.model import StartingModel
        >>> layered = StartingModel([100.0, 300.0], [500.0]).to_layered_model()
        >>> layered.n_layers
        2
        """
        from ..forward import LayeredModel

        return LayeredModel(
            resistivity=self.resistivities.copy(),
            thickness=self.thicknesses.copy(),
            name=self.name,
        )

    def validate(self) -> None:
        """Validate layer shape and positivity.

        Raises
        ------
        ValueError
            If resistivities or thicknesses are not 1-D, if fewer than two
            resistivity layers are supplied, if thickness count does not equal
            ``len(resistivities) - 1``, or if any value is non-positive.

        Examples
        --------
        >>> from pycsamt.inversion.model import StartingModel
        >>> StartingModel([100.0, 300.0], [500.0]).validate() is None
        True
        """
        if self.resistivities.ndim != 1 or self.thicknesses.ndim != 1:
            raise ValueError("resistivities and thicknesses must be 1-D.")
        if self.resistivities.size < 2:
            raise ValueError("at least two layers are required.")
        if self.thicknesses.size != self.resistivities.size - 1:
            raise ValueError("len(thicknesses) must be len(resistivities)-1.")
        if np.any(self.resistivities <= 0):
            raise ValueError("resistivities must be strictly positive.")
        if np.any(self.thicknesses <= 0):
            raise ValueError("thicknesses must be strictly positive.")


ReferenceModel = StartingModel

StartingModel.__doc__ = f"""
Starting or recovered layered-earth model.

Resistivities are linear ohm-m values. The final layer is a halfspace, so
``len(thicknesses)`` must equal ``len(resistivities) - 1``.

``ReferenceModel`` is currently an alias of ``StartingModel``. Use it when the
same layer container is being supplied as a regularization reference rather than
as an optimizer starting point.

Parameters
----------
{_inversion_param_docs.model.resistivities}
{_inversion_param_docs.model.thicknesses}
{_inversion_param_docs.model.name}
{_inversion_param_docs.model.metadata}

Notes
-----
The model is deliberately small and backend-neutral. Built-in solvers optimize
``log10(resistivity)`` and ``log10(thickness)`` internally, but this public
container stores physical linear values in ohm metres and metres.

{_inversion_param_docs.model.examples}

{_inversion_param_docs.model.references}
"""
