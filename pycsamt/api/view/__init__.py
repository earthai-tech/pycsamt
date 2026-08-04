"""User-facing result views for pyCSAMT public APIs."""

from .config import (
    PYCSAMT_API_VIEW,
    APIViewConfig,
    configure_api_view,
    reset_api_view,
)
from .frame import (
    APIFrame,
    FrameProfile,
    api_frame,
    maybe_wrap_frame,
    wrap_frame,
)
from .io import read_edi, read_edis, read_sites
from .progress import (
    PYCSAMT_PROGRESS,
    ProgressAPIConfig,
    ProgressConfig,
    configure_progress,
    get_progress_bar,
    iter_progress,
    progress_enabled,
    reset_progress,
)
from .result import APIResult, wrap_result
from .survey import APISurvey
from .tables import (
    geology_dataframe,
    geology_table,
    quality_dataframe,
    quality_table,
    sites_summary,
)

__all__ = [
    "APIFrame",
    "APIResult",
    "APISurvey",
    "APIViewConfig",
    "FrameProfile",
    "PYCSAMT_API_VIEW",
    "PYCSAMT_PROGRESS",
    "ProgressAPIConfig",
    "ProgressConfig",
    "api_frame",
    "configure_api_view",
    "configure_progress",
    "get_progress_bar",
    "iter_progress",
    "maybe_wrap_frame",
    "progress_enabled",
    "read_edi",
    "read_edis",
    "read_sites",
    "reset_api_view",
    "reset_progress",
    "wrap_frame",
    "wrap_result",
    "geology_dataframe",
    "geology_table",
    "quality_dataframe",
    "quality_table",
    "sites_summary",
]
