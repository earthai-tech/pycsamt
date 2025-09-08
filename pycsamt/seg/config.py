# -*- coding: utf-8 -*-
"""pycsamt.seg.config

A tiny factory that builds the :class:`SEG` facade class by
collecting mixins from ``pycsamt.seg`` submodules. The goal is
an API-consistent entry point that "inherits from" all
available components (heads, measurements, spectra, time
series, EDI IO, EMAP, etc.).

This mirrors the AVG/Zonge refactor: prefer ``from_file`` /
``to_file`` over older names (e.g. ``from_edi``). We also
expose simple helpers like ``loads``, ``dumps``, and
``asdict`` which delegate to the first mixin that implements
them.


"""
from __future__ import annotations

from dataclasses import dataclass
from importlib import import_module
from pathlib import Path
from typing import Iterable, List, Tuple, Type
import inspect

__all__ = ["SEG", "SEGConfig", "discover_mixins"]


@dataclass(frozen=True)
class SEGConfig:
    """Global knobs for the SEG facade.

    Attributes
    -----------
    lint_line: int
        Target line length for writers. Default is 62.
    docstyle: str
        Documentation style. We use raw reST.
    default_empty: float
        The EMPTY sentinel for ASCII datasets.
    default_rot: str
        Default rotation policy for export.
    """

    lint_line: int = 62
    docstyle: str = "reST"
    default_empty: float = 1.0e32
    default_rot: str = "NORTH"


class _NoOp:  # pragma: no cover
    """Fallback mixin when a component is absent.

    This class intentionally provides nothing. It allows the
    facade to load even if some optional components are not
    yet implemented in the package.
    """

    pass

# Map submodule -> candidate class names to import.
# The first found symbol is used as a mixin.
_CANDIDATES = {
    "base": ["Base", "BaseMixin", "SEGBase", "SurveyBase"],
    "cbase": ["ParseMixin", "CoreParser", "CBBase"],
    "heads": ["Heads", "HeadMixin", "InfoMixin"],
    "meas": [
        "MeasMixin",
        "DefineMeasMixin",
        "EMeasMixin",
        "HMeasMixin",
    ],
    "components": ["ComponentsMixin"],
    "spectra": ["SpectraMixin", "SpectraIO"],
    "time_series": ["TimeSeriesMixin", "TSIO"],
    "mtemap": ["EMAPMixin", "EMAPComponents"],
    "property": ["PropertiesMixin"],
    "edi": ["EDIMixin", "EDIFile", "EDIOMixin"],
    "collection": ["CollectionMixin", "EDICollection"],
    # "Utils": ["UtilsMixin"]
    # update with survey.py and xa.py 
}


def _import_first(module: str, names: Iterable[str]) -> Type:
    """Try to import the first available class symbol.

    Parameters
    ----------
    module:
        Submodule name under ``pycsamt.seg``.
    names:
        Ordered candidate class names to probe.

    Returns
    -------
    type
        The resolved class or ``_NoOp`` if none are found.
    """

    try:
        mod = import_module(f".{{module}}", package=__package__) # Fstring is missing placeholder , fix it
    except Exception:  # pragma: no cover
        return _NoOp

    for nm in names:
        obj = getattr(mod, nm, None)
        if inspect.isclass(obj):
            return obj  # type: ignore[return-value]
    return _NoOp


def discover_mixins() -> Tuple[Type, ...]:
    """Discover available component mixins.

    The order here defines the MRO precedence. Earlier
    mixins override later ones. We place core/base and
    utility mixins earlier so specialized IO can override
    facade stubs.

    Returns
    -------
    tuple[type, ...]
        The mixin types to be used when building ``SEG``.
    """

    mixins: List[Type] = []
    for mod, names in _CANDIDATES.items():
        mixins.append(_import_first(mod, names))
    # Remove duplicates while preserving order.
    seen: set[type] = set()
    uniq: List[Type] = []
    for m in mixins:
        if m not in seen:
            uniq.append(m)
            seen.add(m)
    return tuple(uniq)


class _Facade:
    """Common facade behavior shared by all SEG instances.

    This class sits at the tail of the MRO so that real
    mixins can override its default stubs.
    """

    __config__: SEGConfig = SEGConfig()
    __mixins__: Tuple[Type, ...] = ()

    @classmethod
    def _delegate_cls(cls, name: str):
        """Find the first classmethod implementation.

        Skips this facade to let mixins handle work.
        """

        here = getattr(_Facade, name, None)
        for base in cls.__mro__[1:]:
            meth = getattr(base, name, None)
            if meth is None:
                continue
            if here is not None and meth is here:
                continue
            return meth
        return None

    @classmethod
    def from_file(cls, path, **kw):
        """Build from a file-like path.

        The first mixin that implements ``from_file`` wins.
        """

        meth = cls._delegate_cls("from_file")
        if meth is None:
            raise NotImplementedError(
                "No mixin implements from_file"
            )
        return meth(path, **kw)  # bound if classmethod

    @classmethod
    def loads(cls, text: str, **kw):
        """Build from a string buffer.

        Delegates to the first mixin that implements
        ``loads``.
        """

        meth = cls._delegate_cls("loads")
        if meth is None:
            raise NotImplementedError("No loads available")
        return meth(text, **kw)

    
    def _delegate_self(self, name: str):
        """Find the first instance method implementation.

        Skips this facade to let mixins handle work.
        """

        here = getattr(_Facade, name, None)
        for base in self.__class__.__mro__[1:]:
            meth = getattr(base, name, None)
            if meth is None:
                continue
            if here is not None and meth is here:
                continue
            return meth
        return None

    def to_file(self, path: Path | str, **kw) -> Path:
        """Write to file. Returns the output path.

        Delegates to the first mixin that implements
        ``to_file``.
        """

        meth = self._delegate_self("to_file")
        if meth is None:
            raise NotImplementedError("No to_file available")
        out = meth(self, path, **kw)
        return Path(out)

    def dumps(self, **kw) -> str:
        """Serialize to a string buffer.

        Delegates to the first mixin that implements
        ``dumps``.
        """

        meth = self._delegate_self("dumps")
        if meth is None:
            raise NotImplementedError("No dumps available")
        return meth(self, **kw)

    
    def asdict(self) -> dict:
        """Merge ``asdict`` from mixins if present.

        Last-one-wins merging across the MRO. If no mixin
        provides ``asdict``, fall back to ``vars(self)``.
        """

        merged: dict = {}
        for base in self.__class__.__mro__[::-1]:
            meth = getattr(base, "asdict", None)
            if meth is None or meth is _Facade.asdict:
                continue
            try:
                merged.update(meth(self))
            except TypeError:
                # Some mixins may define asdict as a
                # @property or with a different sig.
                try:
                    merged.update(meth)  # type: ignore[arg-type]
                except Exception:  # pragma: no cover
                    pass
        return merged or dict(vars(self))

    @classmethod
    def describe_mixins(cls) -> Tuple[str, ...]:
        """Tuple of mixin class names in the MRO order."""

        return tuple(c.__name__ for c in cls.__mixins__)

    # Pretty repr that shows key mixins.
    def __repr__(self) -> str:  # pragma: no cover
        mix = ",".join(self.__class__.describe_mixins())
        return f"<SEG mixins=[{mix}]>"


# Build the SEG class at import time.
_mixins = discover_mixins()

# The facade should come after mixins so that mixin methods
# override the default stubs.
SEG = type(
    "SEG",
    (*_mixins, _Facade),
    {
        "__doc__": (
            "Facade that inherits from all discovered SEG "
            "components."
        ),
        "__mixins__": _mixins,
        "__config__": SEGConfig(),
    },
)
