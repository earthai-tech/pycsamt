"""
pycsamt.ai.inversion
====================

High-level AI-based EM inversion workflows.
"""
from ._sites_bridge import (
    SiteObs1D,
    SiteObs2D,
    obs_to_features_1d,
    sites_to_coords_3d,
    sites_to_features_1d,
    sites_to_obs_1d,
    sites_to_obs_2d,
    sites_to_panel_2d,
)
from .calibration import (
    ConformalPredictor,
    PosteriorCalibrator,
)
from .config import InversionConfig
from .ensemble import EnsembleInverter
from .hybrid1d import HybridInverter1D
from .hybrid2d import HybridInverter2D
from .hybrid3d import HybridInverter3D
from .inv1d import EMInverter1D
from .inv2d import EMInverter2D
from .inv3d import GCNInverter3D
from .joint import JointInverter
from .pinn1d import PINNInverter1D
from .pinn2d import PINNInverter2D
from .pinn3d import PINNInverter3D
from .run_config import RunConfig

__all__ = [
    # Option 1 — supervised AI
    "EMInverter1D",
    "EMInverter2D",
    "GCNInverter3D",
    "JointInverter",
    "EnsembleInverter",
    "ConformalPredictor",
    "PosteriorCalibrator",
    "InversionConfig",
    "RunConfig",
    # Option 2 — physics-informed AI
    "PINNInverter1D",
    "PINNInverter2D",
    "PINNInverter3D",
    # Option 3 — hybrid AI + physics
    "HybridInverter1D",
    "HybridInverter2D",
    "HybridInverter3D",
    # Sites bridge utilities
    "SiteObs1D",
    "SiteObs2D",
    "sites_to_obs_1d",
    "sites_to_obs_2d",
    "sites_to_features_1d",
    "sites_to_panel_2d",
    "sites_to_coords_3d",
    "obs_to_features_1d",
]
