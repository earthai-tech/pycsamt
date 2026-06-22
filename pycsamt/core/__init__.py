# -*- coding: utf-8 -*-

from ._registry import (
    RegistryError,
    Record,
    Manifest,
    ManifestStore,
    FileManifestStore,
    Registry,
    guess_kind,
)
from .base import CoreObject, MTBase, TFBundle, to_edi
from .registry import (
    Packer,
    register_packer,
    get_packer,
    list_packers,
    pack_to_file,
    unpack_from_file,
    RegistryAPI,
)
from .mixins import (
    bundle_from_edi,
    BundleMixin,
    BundleContainerMixin,
)
from .config import (
    StationNamePolicy,
    CoreConfig,
    get_config,
    configure,
    reset_config,
    config_context,
    to_dict,
    register_adapter,
    get_adapter,
    list_adapters,
)

# Transformers were moved from pycsamt.core to the top-level
# pycsamt.transformers package.  They are re-exported here lazily for
# backwards compatibility: pycsamt.transformers.jedi imports back into
# pycsamt.core (base/config), so an eager import at core init time would
# create a partially-initialised circular import.  PEP 562 __getattr__
# defers the import until the symbol is actually accessed.
_LAZY_TRANSFORMERS = {"AVGtoEDI", "JtoEDI", "TransformerMixin"}


def __getattr__(name):
    if name in _LAZY_TRANSFORMERS:
        from .. import transformers as _t
        return getattr(_t, name)
    raise AttributeError(
        f"module {__name__!r} has no attribute {name!r}"
    )


__all__ = [
    # Transformers / mixins
    "TransformerMixin",
    "AVGtoEDI",
    "JtoEDI",

    # Base exports
    "CoreObject",
    "MTBase",
    "TFBundle",
    "to_edi",

    # Registry (low-level)
    "RegistryError",
    "Record",
    "Manifest",
    "ManifestStore",
    "FileManifestStore",
    "Registry",
    "guess_kind",

    # Registry (high-level) + packers
    "Packer",
    "register_packer",
    "get_packer",
    "list_packers",
    "pack_to_file",
    "unpack_from_file",
    "RegistryAPI",

    # Bundle helpers / mixins
    "bundle_from_edi",
    "BundleMixin",
    "BundleContainerMixin",

    # Config & adapter registry
    "StationNamePolicy",
    "CoreConfig",
    "get_config",
    "configure",
    "reset_config",
    "config_context",
    "to_dict",
    "register_adapter",
    "get_adapter",
    "list_adapters",
]
