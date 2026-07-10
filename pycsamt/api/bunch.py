# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later
"""
The `structure` module includes data structures such as
`Bunch` and `FlexDict` that provide flexible and efficient ways
to organize, store, and manipulate structured data.
"""

from .util import format_dict_result

__all__ = ['Bunch', 'FlexDict']


class Bunch(dict):
    """
    A container object that extends dictionaries by enabling
    attribute-like access to its items.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def __getattr__(self, key):
        try:
            return self[key]
        except KeyError:
            message = (
                f"'Bunch' object has no attribute '{key}'"
            )
            raise AttributeError(message)

    def __setattr__(self, key, value):
        self[key] = value

    def __setstate__(self, state):
        pass

    def __dir__(self):
        return super().__dir__() + list(self.keys())

    def __repr__(self):
        keys = ', '.join(list(self.keys()))
        return f"<Bunch object with keys: {keys}>"

    def __str__(self):
        dict_o = {
            k: v for k, v in self.items() if '__' not in str(k)
        }
        return format_dict_result(
            dict_o, dict_name="Bunch", include_message=True
        )


class FlexDict(dict):
    """
    A dictionary subclass that provides flexible attribute-style
    access to its items.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.__dict__ = self

    def __getattr__(self, item):
        try:
            return self[item]
        except KeyError:
            message = (
                f"'FlexDict' object has no attribute '{item}'"
            )
            raise AttributeError(message)

    def __setattr__(self, key, value):
        special_symbols = [
            '**', '%%', '&&', '||', '$$'
        ]
        for symbol in special_symbols:
            if symbol in key:
                key = key.split(symbol)[0]
                break
        self[key] = value

    def __setstate__(self, state):
        self.update(state)
        self.__dict__ = self

    def __dir__(self):
        return list(self.keys())

    def __repr__(self):
        keys = ', '.join(self.keys())
        return f"<FlexDict with keys: {keys}>"
