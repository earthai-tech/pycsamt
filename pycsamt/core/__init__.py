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
from .transformers import AVGtoEDI, JtoEDI
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
