# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""
Unit conversion utilities.
"""

import re
from typing import Optional, Union

__all__ = [
    'convert_value',
    'convert_temperature',
    'convert_time',
    'convert'
]

# Metric prefixes and their factors
_PREFIXES = {
    'y': 1e-24, 'z': 1e-21, 'a': 1e-18, 'f': 1e-15,
    'p': 1e-12, 'n': 1e-9,  'u': 1e-6, 'µ': 1e-6,
    'm': 1e-3, 'c': 1e-2,  'd': 1e-1, '': 1.0,
    'da':1e1, 'h': 1e2,  'k': 1e3,  'M': 1e6,
    'G': 1e9, 'T': 1e12, 'P': 1e15, 'E': 1e18,
    'Z': 1e21, 'Y': 1e24
}

def convert_value(
    value: Union[str, float, int],
    /,
    target_unit: str = 'm'
) -> float:
    """
    Convert a numeric value with unit suffix to specified unit.

    Parameters
    ----------
    value : str, int, or float
        Input value, e.g.  "20mm", "5kg", or numeric.
    target_unit : str, default='m'
        Desired output unit, e.g. 'm', 'km', 'g', 'kg'.

    Returns
    -------
    converted : float
        Numeric value expressed in `target_unit`.

    Raises
    ------
    ValueError
        If parsing fails or unsupported unit detected.

    Notes
    -----
    - Supports standard metric prefixes (y, z, a, f, p, n, µ, u,
      m, c, d, da, h, k, M, G, T, P, E, Z, Y).
    - Base unit (e.g. 'm' for length, 'g' for mass) must match
      between input and `target_unit`.

    Examples
    --------
    >>> from pycsamt.utils.conversion import convert_value
    >>> convert_value('20mm')
    0.02
    >>> convert_value('3.5km', target_unit='m')
    3500.0
    >>> convert_value('1.2kg', target_unit='g')
    1200.0
    >>> convert_value(5, target_unit='m')
    5.0
    >>> convert_value('100', target_unit='cm')
    10000.0
    >>> convert_value('10m', target_unit='g')  # mismatched base
    Traceback (most recent call last):
      ...
    ValueError: Incompatible units 'm' and 'g'
    """
    # Coerce to string for parsing
    s = str(value).strip()
    # Regex: capture number and optional unit suffix
    m = re.fullmatch(
        r"([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)"
        r"\s*([a-zA-Zµ]+)?", s
    )
    if not m:
        raise ValueError(f"Cannot parse value {value!r}")
    num_str, unit_in = m.groups()
    num = float(num_str)
    # If no suffix, assume already in target_unit
    if not unit_in:
        return num
    # Split input unit into prefix and base
    for pre in sorted(_PREFIXES, key=len, reverse=True):
        if unit_in.startswith(pre) and len(pre) < len(unit_in):
            base_in = unit_in[len(pre):]
            factor_in = _PREFIXES[pre]
            break
    else:
        # no prefix: empty prefix, entire unit is base
        base_in = unit_in
        factor_in = 1.0
    # Split target_unit likewise
    tu = target_unit.strip()
    for pre_t in sorted(_PREFIXES, key=len, reverse=True):
        if tu.startswith(pre_t) and len(pre_t) < len(tu):
            base_t = tu[len(pre_t):]
            factor_t = _PREFIXES[pre_t]
            break
    else:
        base_t = tu
        factor_t = 1.0
    # Check base units match
    if base_in != base_t:
        raise ValueError(
            f"Incompatible units '{unit_in}' and '{tu}'"
        )
    # Convert: numeric * (in_factor / out_factor)
    return num * factor_in / factor_t


def convert_temperature(
    value: Union[str, float, int],
    /,
    unit_from: str = 'C',
    unit_to: str = 'K'
) -> float:
    """
    Convert temperature between Celsius, Fahrenheit, Kelvin.

    Parameters
    ----------
    value : str or number
        e.g. '32F', '100 C', 273.15.
    unit_from : {'C','F','K'}, default 'C'
        Source unit.
    unit_to : {'C','F','K'}, default 'K'
        Target unit.

    Returns
    -------
    float
        Converted temperature.

    Raises
    ------
    ValueError
        On parsing failure or unsupported units.

    Examples
    --------
    >>> convert_temperature('0C','C','K')
    273.15
    >>> convert_temperature(100,'C','F')
    212.0
    """
    s = str(value).strip()
    m = re.fullmatch(
        r"([+-]?\d+(?:\.\d+)?)(?:\s*([cCfFkK]))?", s
    )
    if not m:
        raise ValueError(f"Cannot parse temperature {value!r}")
    num, suf = m.groups()
    t = float(num)
    # determine units
    u_from = (suf or unit_from)[0].upper()
    u_to = unit_to[0].upper()
    # to Celsius
    if u_from == 'C':
        c = t
    elif u_from == 'K':
        c = t - 273.15
    elif u_from == 'F':
        c = (t - 32) * 5/9
    else:
        raise ValueError(f"Unsupported unit_from '{unit_from}'")
    # from Celsius to target
    if u_to == 'C':
        return c
    if u_to == 'K':
        return c + 273.15
    if u_to == 'F':
        return c * 9/5 + 32
    raise ValueError(f"Unsupported unit_to '{unit_to}'")


def convert_time(
    value: Union[str, float, int],
    /,
    unit_from: str = 's',
    unit_to: str = 'h'
) -> float:
    """
    Convert time between seconds, minutes, hours, days.

    Parameters
    ----------
    value : str or number
        e.g. '60s', '2min', 3600.
    unit_from : {'s','min','h','d'}, default 's'
        Source unit.
    unit_to : {'s','min','h','d'}, default 'h'
        Target unit.

    Returns
    -------
    float
        Converted time.

    Raises
    ------
    ValueError
        On parsing failure or unsupported units.

    Examples
    --------
    >>> convert_time('3600s','s','h')
    1.0
    >>> convert_time(2,'h','min')
    120.0
    """
    s = str(value).strip()
    m = re.fullmatch(
        r"([+-]?\d+(?:\.\d+)?)(?:\s*(s|min|h|d))?", s,
        flags=re.IGNORECASE
    )
    if not m:
        raise ValueError(f"Cannot parse time {value!r}")
    num, suf = m.groups()
    t = float(num)
    # normalize units
    u_from = (suf or unit_from).lower()
    u_to = unit_to.lower()
    # to seconds
    if u_from == 's':
        sec = t
    elif u_from == 'min':
        sec = t * 60
    elif u_from == 'h':
        sec = t * 3600
    elif u_from == 'd':
        sec = t * 86400
    else:
        raise ValueError(f"Unsupported unit_from '{unit_from}'")
    # seconds to target
    if u_to == 's':
        return sec
    if u_to == 'min':
        return sec / 60
    if u_to == 'h':
        return sec / 3600
    if u_to == 'd':
        return sec / 86400
    raise ValueError(f"Unsupported unit_to '{unit_to}'")


def convert(
    value: Union[str, float, int],
    /,
    unit_from: str,
    unit_to: str,
    *,
    category: Optional[str] = None,
    round_result: Optional[int] = None,
    **kwargs
) -> float:
    """
    Universal conversion dispatch for length, mass, time, temp.

    Parameters
    ----------
    value : str or number
        Input may include unit suffix.
    unit_from : str
        Source unit, e.g. 'm','km','C','F','s','h'.
    unit_to : str
        Target unit.
    category : {'length','mass','time','temperature'}, optional
        Force conversion category.
    round_result : int, optional
        Decimals to round the output.
    **kwargs
        Passed to the specific converter.

    Returns
    -------
    float
        Converted value.

    Raises
    ------
    ValueError
        If units incompatible or unknown.

    Examples
    --------
    >>> convert(1,'km','m')
    1000.0
    >>> convert('20mm','mm','m')
    0.02
    >>> convert(100,'C','F',category='temperature')
    212.0
    >>> convert('3600s','s','h')
    1.0
    """
    uf = unit_from.strip()
    ut = unit_to.strip()
    cat = category.lower() if category else None
    # temperature
    if cat == 'temperature' or (
        uf[:1].upper() in ('C','F','K') and
        ut[:1].upper() in ('C','F','K')
    ):
        res = convert_temperature(value, uf, ut)
    # time
    elif cat == 'time' or (
        uf.lower() in ('s','min','h','d') and
        ut.lower() in ('s','min','h','d')
    ):
        res = convert_time(value, uf, ut)
    # metric
    else:
        sval = str(value).strip()
        if re.fullmatch(r"[+-]?\d+(?:\.\d+)?", sval):
            sval = f"{sval}{uf}"
        res = convert_value(sval, ut)
    if isinstance(round_result, int):
        res = round(res, round_result)
    return res
