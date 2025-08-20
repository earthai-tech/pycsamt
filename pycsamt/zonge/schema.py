# -*- coding: utf-8 -*-
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
from typing import Dict, List, Literal, Optional, Tuple

__all__ = [
    "_CANONICAL_MAP", 
    "_CANON_TO_MODERN", 
    "_CANON_TO_LEGACY",
    "ALL_ALIASES", 
    "QC_ALIASES",
    "get_aliases"
]


_CANONICAL_MAP: Dict[str, str] = {
    # Coords $ Identifiers ---
    'Station': 'station', 
    'Stn': 'station',
    'Freq': 'freq', 
    'Freq.': 'freq',
    'Comp': 'comp',
    'skp': 'skp',
    
    # --- Measurements 
    'Amps': 'amps', 
    'Tx.Amp': 'amps',
    'Emag': 'emag',
    'E.mag': 'emag',
    'Ephz': 'ephz',
    'E.phz': 'ephz',
    'Hmag': 'hmag', 
    'B.mag': 'hmag', 
    'H.mag': 'hmag',
    'Hphz': 'hphz', 
    'B.phz': 'hphz',
    'Resistivity': 'rho',
    'ARes.mag': 'rho',
    'Phase': 'phase', 
    'Z.phz': 'phase',
    'Z.mag': 'zmag',
    'SRes': 'rho_sc', 
    'TMARES': 'rho_sc',
    'TMARES/SRES': 'rho_sc',
    
    # Quality Control (QC) 
    '%Emag': 'pc_emag',
    'E.%err': 'pc_emag',
    'sEphz': 's_ephz',
    'E.perr': 's_ephz',
    '%Hmag': 'pc_hmag',
    'B.%err': 'pc_hmag', 
    'H.%err': 'pc_hmag',
    'sHphz': 's_hphz', 
    'B.perr': 's_hphz', 
    'H.perr': 's_hphz',
    '%Rho': 'pc_rho', 
    'ARes.%err': 'pc_rho',
    'sPhz': 's_phz', 
    'Z.perr': 's_phz',
    
    # Weights & Other Modern Fields 
    'Z.mwgt': 'z_mwgt',
    'Z.pwgt': 'z_pwgt',
    'E.wgt': 'e_wgt',
    'H.wgt': 'h_wgt',
    'B.wgt': 'h_wgt',
    'Choer': 'coh',
    'Gdp.Blk': 'gdp_blk',
    'Gdp.Chn': 'gdp_chn',
    'Gdp.Time': 'gdp_time',
    '|Z|': 'zabs',
    'Z.%err': 'z.%err',
    
}


_CANON_TO_MODERN: Dict[str, str] = {
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
    "e.%err":   "E.%err",
    "e.perr":   "E.perr",
    "h.%err":   "B.%err",
    "h.perr":   "B.perr",
    "rho.%err": "ARes.%err",
    "phase.%err": "Z.perr",
    "z.%err":   "Z.%err",
    "z.perr":   "Z.perr",
    "phase.%err": "Z.perr",
    "z.%err":   "Z.%err",
    "z.mwgt": "Z.mwgt",
    "z.pwgt": "Z.pwgt",
    "e.wgt":  "E.wgt",
    "b.wgt":  "B.wgt",
}

_CANON_TO_LEGACY: Dict[str, str] = {
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

# 3. Dynamically Built Alias Lookups
_canon_to_aliases: Dict[str, List[str]] = defaultdict(list)
for alias, canon in _CANONICAL_MAP.items():
    _canon_to_aliases[canon].append(alias)

ALL_ALIASES: Dict[str, Tuple[str, ...]] = {
    canon: tuple(sorted(aliases))
    for canon, aliases in _canon_to_aliases.items()
}

QC_ALIASES: Dict[str, Tuple[str, ...]] = {
    "pc_emag": tuple(
        sorted([k for k, v in _CANONICAL_MAP.items() if v == 'pc_emag'])
    ),
    "s_ephz": tuple(
        sorted([k for k, v in _CANONICAL_MAP.items() if v == 's_ephz'])
    ),
    "pc_hmag": tuple(
        sorted([k for k, v in _CANONICAL_MAP.items() if v == 'pc_hmag'])
    ),
    "s_hphz": tuple(
        sorted([k for k, v in _CANONICAL_MAP.items() if v == 's_hphz'])
    ),
    "pc_rho": tuple(
        sorted([k for k, v in _CANONICAL_MAP.items() if v == 'pc_rho'])
    ),
    "s_phz": tuple(
        sorted([k for k, v in _CANONICAL_MAP.items() if v == 's_phz'])
    ),
}


def get_aliases(
    canonical_name: str,
    *,
    kind: Optional[Literal["qc", "all"]] = "all",
    custom_aliases: Optional[Dict[str, Tuple[str, ...]]] = None,
    normalize: bool = True,
) -> Tuple[str, ...]:
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
        norm_custom = {
            _normalize_key(k): v for k, v in custom_aliases.items()
        } if normalize else custom_aliases

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
                active_qc[lookup_canon] = tuple(
                    sorted(set(existing_qc))
                )

    # Select the correct map to use
    target_map = active_qc if kind == "qc" else active_all
    # Normalize the final lookup key
    lookup_name = _normalize_key(canonical_name)

    return target_map.get(lookup_name, ())