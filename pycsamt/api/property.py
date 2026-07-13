"""Common object and metadata helpers for pyCSAMT APIs."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from copy import copy
from dataclasses import asdict, fields, is_dataclass
from typing import Any, ClassVar

__all__ = [
    "MetadataMixin",
    "PyCSAMTObject",
]


class PyCSAMTObject:
    """Root object for lightweight pyCSAMT API behavior.

    ``PyCSAMTObject`` provides a small, dependency-light base for
    readable object display and shallow introspection. It is not an
    electromagnetic numerical class and does not assume EDI, MT,
    TEM, ModEM, Occam, or AI semantics.

    Subclasses may override ``__repr__`` or ``__str__`` freely.
    For smaller customization, define ``__repr_fields__`` or
    ``__repr_exclude__`` on the subclass.

    Attributes
    ----------
    __repr_fields__ : iterable of str or None
        Explicit field order for ``__repr__``. If ``None``,
        dataclass fields, ``__dict__`` keys, and ``__slots__``
        names are discovered automatically.
    __repr_exclude__ : set of str
        Attribute names excluded from automatic displays.
    __repr_max_fields__ : int
        Maximum number of fields shown in ``__repr__``.
    __repr_max_string__ : int
        Maximum string length before shortening.

    Examples
    --------
    >>> class Obj(PyCSAMTObject):
    ...     def __init__(self):
    ...         self.name = "S001"
    ...         self.verbose = 2
    ...         self.values = [1, 2, 3]
    >>> repr(Obj())
    "Obj(name='S001', values=list(len=3, sample=[1, 2, 3]))"
    """

    __repr_fields__: ClassVar[Iterable[str] | None] = None
    __repr_exclude__: ClassVar[set[str]] = {"logger", "verbose"}
    __repr_max_fields__: ClassVar[int] = 8
    __repr_max_string__: ClassVar[int] = 48

    def __repr__(self) -> str:
        """Return a compact developer-facing object display."""
        fields_ = list(self._repr_fields())
        max_fields = int(self.__repr_max_fields__)
        shown = fields_[:max_fields]
        items = []
        for name in shown:
            try:
                value = getattr(self, name)
            except Exception:
                continue
            items.append(f"{name}={self._repr_value(value)}")
        if len(fields_) > max_fields:
            items.append("...")
        return f"{self._display_name()}({', '.join(items)})"

    def __str__(self) -> str:
        """Return a user-facing compact summary."""
        text = self.summary()
        return text if isinstance(text, str) else repr(text)

    def summary(self, *, max_fields: int | None = None) -> str:
        """Return a short one-line object summary."""
        old = self.__repr_max_fields__
        if max_fields is None:
            return repr(self)
        try:
            self.__repr_max_fields__ = int(max_fields)
            return repr(self)
        finally:
            self.__repr_max_fields__ = old

    def to_dict(
        self,
        *,
        public_only: bool = True,
        max_depth: int = 1,
    ) -> dict[str, Any]:
        """Return a shallow dictionary representation."""
        return self._to_dict(self, public_only=public_only, depth=max_depth)

    def update(self, /, **kwargs: Any) -> PyCSAMTObject:
        """Set attributes from keyword arguments and validate."""
        for name, value in kwargs.items():
            setattr(self, name, value)
        self.validate()
        return self

    def clone(self, /, **overrides: Any) -> PyCSAMTObject:
        """Return a shallow copy with optional attribute overrides."""
        obj = copy(self)
        for name, value in overrides.items():
            setattr(obj, name, value)
        obj.validate()
        return obj

    def validate(self) -> None:
        """Validate object state.

        Subclasses can override this hook. The base implementation
        intentionally accepts all states.
        """
        return None

    def _display_name(self) -> str:
        """Return the class name used by display helpers."""
        return self.__class__.__name__

    def _repr_fields(self) -> Iterable[str]:
        """Return field names to include in automatic displays."""
        explicit = self.__repr_fields__
        if explicit is not None:
            return self._filter_repr_fields(explicit)

        names: list[str] = []
        if is_dataclass(self):
            names.extend(field.name for field in fields(self))
        if hasattr(self, "__dict__"):
            names.extend(self.__dict__.keys())
        names.extend(self._slot_names())
        return self._filter_repr_fields(dict.fromkeys(names))

    def _filter_repr_fields(
        self,
        names: Iterable[str],
    ) -> Iterable[str]:
        """Filter private and excluded repr field names."""
        exclude = set(self.__repr_exclude__)
        return (
            name
            for name in names
            if isinstance(name, str)
            and not name.startswith("_")
            and name not in exclude
        )

    @classmethod
    def _slot_names(cls) -> list[str]:
        """Return slot names declared across the class MRO."""
        names: list[str] = []
        for klass in reversed(cls.__mro__):
            slots = getattr(klass, "__slots__", ())
            if isinstance(slots, str):
                slots = (slots,)
            names.extend(
                slot
                for slot in slots
                if slot not in {"__dict__", "__weakref__"}
            )
        return names

    def _repr_value(self, value: Any) -> str:
        """Summarize one value for ``__repr__``."""
        return self._short_value(
            value,
            max_string=int(self.__repr_max_string__),
        )

    @classmethod
    def _short_value(cls, value: Any, *, max_string: int = 48) -> str:
        """Return a safe compact representation of ``value``."""
        if value is None or isinstance(value, (bool, int, float, complex)):
            return repr(value)
        if isinstance(value, str):
            if len(value) > max_string:
                return repr(value[: max_string - 3] + "...")
            return repr(value)
        try:
            import numpy as np

            if isinstance(value, np.ndarray):
                return (
                    "ndarray(shape="
                    f"{tuple(value.shape)}, dtype={value.dtype})"
                )
        except Exception:
            pass
        if isinstance(value, Mapping):
            keys = list(value.keys())
            sample = ", ".join(repr(key) for key in keys[:3])
            if len(keys) > 3:
                sample += ", ..."
            return f"dict(len={len(keys)}, keys=[{sample}])"
        if isinstance(value, (bytes, bytearray)):
            return f"{type(value).__name__}(len={len(value)})"
        if isinstance(value, (list, tuple, set, frozenset)):
            seq = list(value)
            sample = ", ".join(
                cls._short_value(item, max_string=max_string)
                for item in seq[:3]
            )
            if len(seq) > 3:
                sample += ", ..."
            return (
                f"{type(value).__name__}(len={len(seq)}, sample=[{sample}])"
            )
        shape = getattr(value, "shape", None)
        dtype = getattr(value, "dtype", None)
        if shape is not None:
            return f"array_like(shape={shape}, dtype={dtype})"
        return repr(value)

    @classmethod
    def _to_dict(
        cls,
        obj: Any,
        *,
        public_only: bool,
        depth: int,
    ) -> Any:
        """Convert common object shapes to dictionaries."""
        if depth < 0:
            return {"_": "max_depth"}
        if obj is None or isinstance(obj, (bool, int, float, complex, str)):
            return obj
        try:
            import numpy as np

            if isinstance(obj, np.ndarray):
                return {
                    "type": "ndarray",
                    "shape": tuple(obj.shape),
                    "dtype": str(obj.dtype),
                }
        except Exception:
            pass
        if is_dataclass(obj):
            exclude = getattr(obj, "__repr_exclude__", set())
            try:
                return {
                    key: value
                    for key, value in asdict(obj).items()
                    if key not in exclude
                }
            except Exception:
                out = {}
                for field in fields(obj):
                    if field.name in exclude:
                        continue
                    value = getattr(obj, field.name, None)
                    out[field.name] = cls._to_dict(
                        value,
                        public_only=public_only,
                        depth=depth - 1,
                    )
                return out
        if isinstance(obj, Mapping):
            return {
                str(key): cls._to_dict(
                    value,
                    public_only=public_only,
                    depth=depth - 1,
                )
                for key, value in list(obj.items())[:32]
            }
        if isinstance(obj, (list, tuple, set, frozenset)):
            return [
                cls._to_dict(
                    value,
                    public_only=public_only,
                    depth=depth - 1,
                )
                for value in list(obj)[:32]
            ]
        names: list[str] = []
        if hasattr(obj, "__dict__"):
            names.extend(obj.__dict__.keys())
        names.extend(
            obj.__class__._slot_names()
            if isinstance(obj, PyCSAMTObject)
            else []
        )
        out: dict[str, Any] = {}
        exclude = getattr(obj, "__repr_exclude__", set())
        for name in dict.fromkeys(names):
            if public_only and str(name).startswith("_"):
                continue
            if name in exclude:
                continue
            try:
                value = getattr(obj, name)
            except Exception:
                continue
            out[str(name)] = cls._to_dict(
                value,
                public_only=public_only,
                depth=depth - 1,
            )
        return out


class MetadataMixin:
    """Small mixin for attaching free-form object metadata.

    The mixin is intentionally separate from :class:`PyCSAMTObject`
    so the universal root class remains lightweight.
    """

    metadata: dict[str, Any]

    def ensure_metadata(self) -> dict[str, Any]:
        """Return the metadata mapping, creating it when needed."""
        current = getattr(self, "metadata", None)
        if not isinstance(current, dict):
            current = {}
            self.metadata = current
        return current

    def update_metadata(self, /, **metadata: Any) -> MetadataMixin:
        """Update metadata in place and return ``self``."""
        self.ensure_metadata().update(metadata)
        return self

    def metadata_dict(self) -> dict[str, Any]:
        """Return a shallow copy of attached metadata."""
        return dict(self.ensure_metadata())
