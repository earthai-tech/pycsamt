# -*- coding: utf-8 -*-

# Core base utilities / classes
from .base import CoreObject, MTBase, TFBundle, to_edi

# Registry primitives (low-level record/manifest/store + helpers)
from ._registry import (
    RegistryError,
    Record,
    Manifest,
    ManifestStore,
    FileManifestStore,
    Registry,
    guess_kind,
)

# High-level Registry API and packing helpers
from .registry import (
    Packer,
    register_packer,
    get_packer,
    list_packers,
    pack_to_file,
    unpack_from_file,
    RegistryAPI,
)

# Bundle helpers / mixins
from .mixins import (
    bundle_from_edi,
    BundleMixin,
    BundleContainerMixin,
)

# Core configuration & adapter registry
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
    # Submodules
    # "base",
    # "transformers",

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
