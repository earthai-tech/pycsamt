"""Public API helpers shared across pyCSAMT packages."""

from .control import (
    PYCSAMT_CONTROL,
    FrequencyAxisControl,
    PhaseViewControl,
    PyCSAMTControl,
    RhoViewControl,
    configure_control,
    reset_control,
    wrap_phase,
)
from .property import MetadataMixin, PyCSAMTObject
from .style import (
    PYCSAMT_STYLE,
    CorrectionStyle,
    MTComponentStyle,
    MultilineStyle,
    PhaseTensorEllipseStyle,
    PyCSAMTStyle,
    RawDataStyle,
    RoseStyle,
    configure_style,
    reset_style,
    use_style,
)
from .plot import (
    PLOT_CONFIG,
    PlotConfig,
    load_plot_config,
    reset_plot_config,
    save_fig,
    set_dpi,
    set_fmt,
    set_savedir,
    write_default_config,
)

__all__ = [
    "MetadataMixin",
    "PyCSAMTObject",
    # control API
    "FrequencyAxisControl",
    "PYCSAMT_CONTROL",
    "PhaseViewControl",
    "PyCSAMTControl",
    "RhoViewControl",
    "configure_control",
    "reset_control",
    "wrap_phase",
    # style API
    "PyCSAMTStyle",
    "MultilineStyle",
    "MTComponentStyle",
    "CorrectionStyle",
    "PhaseTensorEllipseStyle",
    "RawDataStyle",
    "RoseStyle",
    "PYCSAMT_STYLE",
    "use_style",
    "reset_style",
    "configure_style",
    # plot-export API
    "PlotConfig",
    "PLOT_CONFIG",
    "save_fig",
    "set_fmt",
    "set_dpi",
    "set_savedir",
    "reset_plot_config",
    "load_plot_config",
    "write_default_config",
]
