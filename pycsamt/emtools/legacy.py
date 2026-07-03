from __future__ import annotations

from ..compat.aliases import install_compat_aliases

from ..utils.em import (
    check_em_kind,
    extract_z_list,
    parse_tensor,
    compute_qc,
    full_freq,
    tensor2d,
    align_tensor,
    export_edis,
    plot_confidence,
    plot_strike,
    plot_tensors,
    plot_station_tensors,
    wrap_phase,
    plot_lcurve,
)

__all__ = [
    "check_em_kind",
    "extract_z_list",
    "parse_tensor",
    "compute_qc",
    "full_freq",
    "tensor2d",
    "align_tensor",
    "export_edis",
    "plot_confidence",
    "plot_strike",
    "plot_tensors",
    "plot_station_tensors",
    "wrap_phase",
    "plot_lcurve",
]

_ALIAS_MAP = {
    "get_full_frequency": full_freq,
    "qc": compute_qc,
    "plot_tensors2": plot_station_tensors,
    "plot_l_curve": plot_lcurve,
    "get2dtensor": tensor2d
}

_EXTRAS = {
    "plot_tensors2": "Use 'plot_station_tensors'.",
    "get_full_frequency": "Use 'full_freq'.",
    "plot_l_curve": "Use 'plot_lcurve'.",
    "get2dtensor": "Use 'tensor2d'."
}

install_compat_aliases(
    _ALIAS_MAP,
    g=globals(),
    since="2.0.0",
    remove_in="2.17.0",
    extras=_EXTRAS,
)

__all__ += [n for n in _ALIAS_MAP if n not in __all__]

del _ALIAS_MAP
del _EXTRAS
