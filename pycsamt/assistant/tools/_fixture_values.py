# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Dependency-free checks copied into the isolated fixture container."""
import math


def validate_values(value, rule):
    """Raise on invalid JSON factors or a finite array of the declared shape."""
    if rule['kind'] == 'positive_factors':
        if not isinstance(value, dict) or not value:
            raise ValueError('expected a nonempty factor mapping')
        values = list(value.values())
    else:
        def flatten(v, dims):
            if not dims:
                return [v]
            if not isinstance(v, list) or len(v) != dims[0]:
                raise ValueError('array dimensions differ from contract')
            return [x for row in v for x in flatten(row, dims[1:])]
        values = flatten(value, rule['shape'])
    if not all(type(x) in (int, float) and math.isfinite(x) and
               (rule['kind'] != 'positive_factors' or x > 0) for x in values):
        raise ValueError('expected finite numeric values (positive for factors)')
