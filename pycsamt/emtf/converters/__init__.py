"""Format adapters for the scientific EMTF object model."""

from .bundle import bundle_to_emtf, emtf_to_bundle
from .edi import (
    DataLossWarning,
    EMTFEDIConversionError,
    edi_to_emtf,
    emtf_to_edi,
    write_edi,
)
from .spectra import (
    SpectraChannelMap,
    SpectraCovarianceWarning,
    SpectraRecoveryError,
    SpectraRecoveryResult,
    recover_spectra_transfer_functions,
    resolve_spectra_channels,
    spectra_to_emtf,
)

__all__ = [
    "DataLossWarning",
    "EMTFEDIConversionError",
    "edi_to_emtf",
    "emtf_to_edi",
    "write_edi",
    "bundle_to_emtf",
    "emtf_to_bundle",
    "SpectraCovarianceWarning",
    "SpectraRecoveryError",
    "SpectraChannelMap",
    "SpectraRecoveryResult",
    "resolve_spectra_channels",
    "recover_spectra_transfer_functions",
    "spectra_to_emtf",
]
