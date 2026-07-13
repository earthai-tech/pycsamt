"""Reusable docstring fragments for PyCSAMT APIs."""

from __future__ import annotations

import re

__all__ = [
    "DocstringComponents",
]


class DocstringComponents:
    """Store docstring fragments with dot-style access."""

    regexp = re.compile(r"\n((\n|.)+)\n\s*", re.MULTILINE)

    def __init__(
        self,
        comp_dict: dict[str, str],
        strip_whitespace: bool = True,
    ):
        if strip_whitespace:
            entries = {}
            for key, value in comp_dict.items():
                match = re.match(self.regexp, value)
                entries[key] = value if match is None else match.group(1)
        else:
            entries = comp_dict.copy()

        self.entries = entries

    def __getattr__(self, attr: str) -> str:
        if attr in self.entries:
            return self.entries[attr]
        raise AttributeError(attr)

    @classmethod
    def from_nested_components(
        cls, **kwargs: DocstringComponents
    ) -> DocstringComponents:
        return cls(kwargs, strip_whitespace=False)
