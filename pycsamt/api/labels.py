"""Shared scientific labels used by pycsamt plots."""

LOG10_PERIOD_LABEL = r"$\log_{10}(T)$ (s)"
PERIOD_LABEL = "Period (s)"
FREQUENCY_LABEL = "Frequency (Hz)"
STATION_LABEL = "Station"


def period_axis_label(kind: str = "logperiod") -> str:
    """Return the package-standard label for a period/frequency axis."""
    key = str(kind).strip().lower().replace("_", "")
    if key in {"logperiod", "log10period", "log10t", "logt"}:
        return LOG10_PERIOD_LABEL
    if key in {"period", "t"}:
        return PERIOD_LABEL
    if key in {"frequency", "freq", "f"}:
        return FREQUENCY_LABEL
    return str(kind)
