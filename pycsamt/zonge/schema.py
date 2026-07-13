# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later
"""
pycsamt.zonge.schema

Defines the data schema and column name mappings for Zonge
AVG files. This module serves as the single source of truth for
translating between legacy (kind-1), modern (kind-2), and the
internal canonical column names used throughout the package.
"""

from __future__ import annotations

from collections import defaultdict
from typing import Literal

__all__ = [
    "_CANONICAL_MAP",
    "_CANON_TO_MODERN",
    "_CANON_TO_LEGACY",
    "_FLEXIBLE_LOOKUP",
    "_CSAVGW_ORDERED",
    "ALL_ALIASES",
    "QC_ALIASES",
    "UNIT_MAP",
    "get_unit",
    "get_aliases",
]


_CANONICAL_MAP: dict[str, str] = {
    # Coords $ Identifiers ---
    "Station": "station",
    "Stn": "station",
    "Freq": "freq",
    "Freq.": "freq",
    "Comp": "comp",
    "skp": "skp",
    # tipper,
    "tx": "tx",
    "ty": "ty",
    # --- Measurements
    "Amps": "amps",
    "Tx.Amp": "amps",
    "Emag": "emag",
    "E.mag": "emag",
    "Ephz": "ephz",
    "E.phz": "ephz",
    "Hmag": "hmag",
    "B.mag": "hmag",
    "H.mag": "hmag",
    "Hphz": "hphz",
    "B.phz": "hphz",
    "Resistivity": "rho",
    "ARes.mag": "rho",
    "Phase": "phase",
    "Z.phz": "phase",
    "Z.mag": "zmag",
    "SRes": "rho_sc",
    "TMARES": "rho_sc",
    "TMARES/SRES": "rho_sc",
    # Quality Control (QC)
    "%Emag": "pc_emag",
    "E.%err": "pc_emag",
    "sEphz": "s_ephz",
    "E.perr": "s_ephz",
    "%Hmag": "pc_hmag",
    "B.%err": "pc_hmag",
    "H.%err": "pc_hmag",
    "sHphz": "s_hphz",
    "B.perr": "s_hphz",
    "H.perr": "s_hphz",
    "%Rho": "pc_rho",
    "ARes.%err": "pc_rho",
    "sPhz": "s_phz",
    "Z.perr": "s_phz",
    # Weights & Other Modern Fields
    "Z.mwgt": "z_mwgt",
    "Z.pwgt": "z_pwgt",
    "E.wgt": "e_wgt",
    "H.wgt": "h_wgt",
    "B.wgt": "h_wgt",
    "Choer": "coh",
    "Gdp.Blk": "gdp_blk",
    "Gdp.Chn": "gdp_chn",
    "Gdp.Time": "gdp_time",
    "|Z|": "zabs",
    "Z.%err": "z.%err",
}

_CANON_TO_MODERN: dict[str, str] = {
    # survey logistics
    "station": "Station",
    "freq": "Freq",
    "comp": "Comp",
    # transmitter / field data
    "amps": "Tx.Amp",
    "emag": "E.mag",
    "ephz": "E.phz",
    "hmag": "B.mag",
    "hphz": "B.phz",
    "zmag": "Z.mag",
    "phase": "Z.phz",
    "rho": "ARes.mag",
    # quality metrics
    "rho_sc": "SRes",
    "pc_emag": "E.%err",
    "s_ephz": "E.perr",
    "pc_hmag": "B.%err",
    "s_hphz": "B.perr",
    "pc_rho": "ARes.%err",
    "s_phz": "Z.perr",
    "z_mwgt": "Z.mwgt",
    "z_pwgt": "Z.pwgt",
    "e_wgt": "E.wgt",
    "h_wgt": "H.wgt",
    "coh": "Choer",
    "gdp_blk": "Gdp.Blk",
    "gdp_chn": "Gdp.Chn",
    "gdp_time": "Gdp.Time",
    "zabs": "|Z|",
    # alternatives_config,
    "e.%err": "E.%err",
    "e.perr": "E.perr",
    "h.%err": "B.%err",
    "h.perr": "B.perr",
    "rho.%err": "ARes.%err",
    "phase.%err": "Z.perr",
    "z.%err": "Z.%err",
    "z.perr": "Z.perr",
    "z.mwgt": "Z.mwgt",
    "z.pwgt": "Z.pwgt",
    "e.wgt": "E.wgt",
    "b.wgt": "B.wgt",
    # e.g Tipper,
    "tx": "Tx",
    "ty": "Ty",
}

_CANON_TO_LEGACY: dict[str, str] = {
    "station": "Station",
    "freq": "Freq",
    "comp": "Comp",
    "amps": "Amps",
    "emag": "Emag",
    "ephz": "Ephz",
    "hmag": "Hmag",
    "hphz": "Hphz",
    "rho": "Resistivity",
    "phase": "Phase",
    "pc_emag": "%Emag",
    "s_ephz": "sEphz",
    "pc_hmag": "%Hmag",
    "s_hphz": "sHphz",
    "pc_rho": "%Rho",
    "s_phz": "sPhz",
    "skp": "skp",
}

_CSAVGW_ORDERED = [
    "Z.mwgt",
    "Z.pwgt",
    "Freq",
    "Tx.Amp",
    "E.mag",
    "E.phz",
    "B.mag",
    "B.phz",
    "Z.mag",
    "Z.phz",
    "ARes.mag",
    "SRes",
    "E.wgt",
    "H.wgt",
    "E.%err",
    "E.perr",
    "B.%err",
    "B.perr",
    "Z.%err",
    "Z.perr",
    "ARes.%err",
    # Expected CSAVGW order;
    # extras will be appended after these.
    "Choer",
    "Gdp.Blk",
    "Gdp.Chn",
    "Gdp.Time",
    "|Z|",
    "use",
]
# 3. Dynamically Built Alias Lookups
_canon_to_aliases: dict[str, list[str]] = defaultdict(list)
for alias, canon in _CANONICAL_MAP.items():
    _canon_to_aliases[canon].append(alias)

ALL_ALIASES: dict[str, tuple[str, ...]] = {
    canon: tuple(sorted(aliases))
    for canon, aliases in _canon_to_aliases.items()
}

QC_ALIASES: dict[str, tuple[str, ...]] = {
    "pc_emag": tuple(
        sorted([k for k, v in _CANONICAL_MAP.items() if v == "pc_emag"])
    ),
    "s_ephz": tuple(
        sorted([k for k, v in _CANONICAL_MAP.items() if v == "s_ephz"])
    ),
    "pc_hmag": tuple(
        sorted([k for k, v in _CANONICAL_MAP.items() if v == "pc_hmag"])
    ),
    "s_hphz": tuple(
        sorted([k for k, v in _CANONICAL_MAP.items() if v == "s_hphz"])
    ),
    "pc_rho": tuple(
        sorted([k for k, v in _CANONICAL_MAP.items() if v == "pc_rho"])
    ),
    "s_phz": tuple(
        sorted([k for k, v in _CANONICAL_MAP.items() if v == "s_phz"])
    ),
}

_FLEX_CANON_NAMES = {
    "pc_emag",
    "s_ephz",
    "pc_hmag",
    "s_hphz",
    "pc_rho",
    "s_phz",
    "z_mwgt",
    "z_pwgt",
    "e_wgt",
    "h_wgt",
}


def _create_flexible_lookup() -> dict[str, str]:
    """
    Creates a lookup dict that maps normalized variations of
    aliases ONLY for QC and weight columns to their canonical name.
    """
    lookup = {}
    for raw_alias, canon_value in _CANONICAL_MAP.items():
        if canon_value in _FLEX_CANON_NAMES:
            # Normalize the raw alias for broader matching
            norm_key = raw_alias.lower()
            variations = {
                norm_key.replace(".", ""),
                norm_key.replace("_", ""),
                norm_key.replace("%", ""),
            }
            for var in variations:
                lookup[var] = canon_value

    # Remove any keys that are also core data names
    # to prevent collisions (e.g., 'rho' should not map to 'pc_rho').
    core_data_keys = {
        "emag",
        "hmag",
        "rho",
        "phase",
        "ephz",
        "hphz",
        "zmag",
        "rho_sc",
    }
    for key in core_data_keys:
        lookup.pop(key, None)

    return lookup


# Create the flexible map once at module level
_FLEXIBLE_LOOKUP = _create_flexible_lookup()

UNIT_MAP: dict[str, tuple[str, str]] = {
    # Canonical Name -> (Simple Unit, Formatted Label for Plots)
    "rho": ("ohm.m", r"App. Res. ($\Omega \cdot m$)"),
    "rho_sc": ("ohm.m", r"Corrected Res. ($\Omega \cdot m$)"),
    "phase": ("mrad", "Phase (mrad)"),
    "phase_deg": ("deg", r"Phase ($\degree$)"),
    "emag": ("nV/Am", "E-Magnitude (nV/Am)"),
    "hmag": ("pT/A", "H-Magnitude (pT/A)"),
    "ephz": ("mrad", "E-Phase (mrad)"),
    "ephz_deg": ("deg", r"E-Phase ($\degree$)"),
    "hphz": ("mrad", "H-Phase (mrad)"),
    "hphz_deg": ("deg", r"H-Phase ($\degree$)"),
    "zmag": ("km/s", "Impedance Mag. (km/s)"),
    "zabs": ("ohm", r"$|Z|$ ($\Omega$)"),
    "pc_rho": ("%", r"$\rho_a$ Error (%)"),
    "pc_emag": ("%", r"$|E|$ Error (%)"),
    "pc_hmag": ("%", r"$|H|$ Error (%)"),
    "s_phz": ("mrad", r"$\sigma_{\phi}$ (mrad)"),
    "s_phz_deg": ("deg", r"$\sigma_{\phi}$ ($\degree$)"),
    "s_ephz": ("mrad", r"$\sigma_{E\phi}$ (mrad)"),
    "s_ephz_deg": ("deg", r"$\sigma_{E\phi}$ ($\degree$)"),
    "s_hphz": ("mrad", r"$\sigma_{H\phi}$ (mrad)"),
    "s_hphz_deg": ("deg", r"$\sigma_{H\phi}$ ($\degree$)"),
    "skew": ("", "Skew"),
    "strike_angle": ("degrees", "Strike Angle (degrees)"),
    "tx": ("", "Tipper Tx"),
    "ty": ("", "Tipper Ty"),
    "freq": ("Hz", "Frequency (Hz)"),
    "station": ("m", "Station Distance (m)"),
}


def get_unit(canonical_name: str, formatted: bool = True) -> str | None:
    """
    Fetch the unit string for a canonical variable name.
    """
    unit_tuple = UNIT_MAP.get(canonical_name)
    if unit_tuple:
        return unit_tuple[1] if formatted else unit_tuple[0]
    return None


def get_aliases(
    canonical_name: str,
    *,
    kind: Literal["qc", "all"] | None = "all",
    custom_aliases: dict[str, tuple[str, ...]] | None = None,
    normalize: bool = True,
) -> tuple[str, ...]:
    r"""Fetch all known aliases for a canonical column name.

    Parameters
    ----------
    canonical_name : str
        The internal, standardized name for the column.
    kind : {'qc', 'all'}, default 'all'
        - 'qc': Search only within quality control aliases.
        - 'all': Search within all known aliases.
    custom_aliases : dict, optional
        A dictionary to temporarily add or override aliases.
        Keys are canonical names, values are tuples of aliases.
    normalize : bool, default True
        If ``True``, the `canonical_name` and keys in the
        `custom_aliases` dictionary are treated as
        case-insensitive.

    Returns
    -------
    tuple[str, ...]
        A tuple of all raw (legacy and modern) names that map
        to the given canonical name.

    Examples
    --------
    >>> get_aliases('pc_emag', kind='qc')
    ('%Emag', 'E.%err')
    >>> get_aliases('rho')
    ('ARes.mag', 'Resistivity')
    >>> get_aliases(
    ...     'rho',
    ...     custom_aliases={'rho': ('resistivity_ohm_m',)}
    ... )
    ('ARes.mag', 'Resistivity', 'resistivity_ohm_m')
    """
    active_qc = QC_ALIASES.copy()
    active_all = ALL_ALIASES.copy()

    def _normalize_key(key):
        return key.lower() if normalize else key

    if custom_aliases:
        # Normalize keys of custom dict if requested
        norm_custom = (
            {_normalize_key(k): v for k, v in custom_aliases.items()}
            if normalize
            else custom_aliases
        )

        for canon, aliases in norm_custom.items():
            # Normalize canonical name for lookup
            lookup_canon = _normalize_key(canon)
            # Update all aliases map
            existing = list(active_all.get(lookup_canon, ()))
            existing.extend(aliases)
            active_all[lookup_canon] = tuple(sorted(set(existing)))

            # Update QC aliases map if applicable
            if kind == "qc" and lookup_canon in active_qc:
                existing_qc = list(active_qc.get(lookup_canon, ()))
                existing_qc.extend(aliases)
                active_qc[lookup_canon] = tuple(sorted(set(existing_qc)))

    # Select the correct map to use
    target_map = active_qc if kind == "qc" else active_all
    # Normalize the final lookup key
    lookup_name = _normalize_key(canonical_name)

    return target_map.get(lookup_name, ())
