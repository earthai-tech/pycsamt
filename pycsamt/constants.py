# -*- coding: utf-8 -*-
#       Author: LKouadio <etanoyau@gmail.com>
#       License: LGPL-3.0
"""
pycsamt.constants

Fundamental mathematical & physical constants for EM geophysics
(AMT | MT | CSAMT).

Examples
--------
>>> from pycsamt.constants import MU_0, RHO_FACTOR, DEG2RAD
>>> Z_mod  = 1.20              # |Z| [V A⁻¹]
>>> freq   = 1_000.0           # 1 kHz
>>> rho_a  = (Z_mod**2) * RHO_FACTOR / freq     # apparent ρ [Ω·m]
"""

from __future__ import annotations

import math
from typing import Final


__all__ = [
    # maths
    "PI", "TAU", "DEG2RAD", "RAD2DEG", "LN10", "MRAD",
    "_RAD_THR", "_DEG_SCALE",
    # numerical helpers
    "EPS_TOL",
    # fundamental EM
    "MU_0", "EPSILON_0", "C_0", "ETA_0",
    # geophysical helpers
    "RHO_FACTOR",
]
# ---------------------------------------------------------------------
# Pure mathematics
# ---------------------------------------------------------------------

PI:        Final[float] = math.pi
TAU:       Final[float] = math.tau                    # 2 π
DEG2RAD:   Final[float] = PI / 180.0                  # °
RAD2DEG:   Final[float] = 180.0 / PI                  # rad
LN10:      Final[float] = math.log(10.0)
MRAD:      Final[float] = 1.0e3                       # 1 mrad = 10⁻³ rad

# ------------------------------------------------------------------
# Small-angle / auto-unit heuristics
# ------------------------------------------------------------------

# 5° expressed in *radians* – used in a few auto-unit check helpers
_RAD_THR:   Final[float] = math.radians(5.0)
# factor to convert **radians → degrees** via multiplication
_DEG_SCALE: Final[float] = RAD2DEG

# ------------------------------------------------------------------
# Numeric tolerances
# ------------------------------------------------------------------

# Generic *ε* for “close-to-zero” comparisons (unit-less)
EPS_TOL:    Final[float] = 1.0e-9

# ---------------------------------------------------------------------
# Electromagnetic (vacuum) – CODATA 2018
# ---------------------------------------------------------------------

MU_0:       Final[float] = 4.0e-7 * PI  # H m⁻¹ – permeability # ≈ 1.2566370614 × 10-6
EPSILON_0:  Final[float] = 8.854_187_817e-12          # F m⁻¹ – permittivity
C_0:        Final[float] = 1.0 / math.sqrt(MU_0 * EPSILON_0)  # m s⁻¹
ETA_0:      Final[float] = math.sqrt(MU_0 / EPSILON_0)        # Ω

# ---------------------------------------------------------------------
# Handy factor for apparent-resistivity computations
# ---------------------------------------------------------------------
#   ρₐ = |Z|² / (μ₀ 2π f)   →   |Z|² * RHO_FACTOR / f
RHO_FACTOR: Final[float] = 1.0 / (MU_0 * 2.0 * PI)

