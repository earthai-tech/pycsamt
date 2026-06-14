"""User-facing result views for pyCSAMT public APIs."""

from .config import (
    APIViewConfig,
    PYCSAMT_API_VIEW,
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
from .progress import ProgressConfig, iter_progress, progress_enabled
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
    "ProgressConfig",
    "api_frame",
    "configure_api_view",
    "iter_progress",
    "maybe_wrap_frame",
    "progress_enabled",
    "read_edi",
    "read_edis",
    "read_sites",
    "reset_api_view",
    "wrap_frame",
    "wrap_result",
    "geology_dataframe",
    "geology_table",
    "quality_dataframe",
    "quality_table",
    "sites_summary",
]
