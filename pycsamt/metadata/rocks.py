__all__ = [
    "GEO_ROCK_RESISTIVITY",
    "ROCK_HATCH_PATTERNS",
    "RockProperties",
]

# Global resistivity ranges for common geological units
GEO_ROCK_RESISTIVITY: dict[str, list[float]] = {
    "hard rock": [1e99, 1e6],
    "igneous rock": [1e6, 1e3],
    "duricrust": [5.1e3, 5.1e2],
    "gravel/sand": [1e4, 7.943],
    "conglomerate": [1e4, 8.913e1],
    "dolomite/limestone": [1e5, 1e3],
    "permafrost": [1e5, 4.169e2],
    "metamorphic rock": [5.1e2, 1e1],
    "tills": [8.1e2, 8.512e1],
    "standstone conglomerate": [1e4, 8.318e1],
    "lignite/coal": [7.762e2, 1e1],
    "shale": [5.012e1, 3.20e1],
    "clay": [1e2, 5.012e1],
    "saprolite": [6.310e2, 3.020e1],
    "sedimentary rock": [1e4, 1e0],
    "fresh water": [3.1e2, 1e0],
    "salt water": [1e0, 1.41e0],
    "massive sulphide": [1e0, 1e-2],
    "sea water": [1.231e-1, 1e-1],
    "ore minerals": [1e0, 1e-4],
    "graphite": [3.1623e-2, 3.162e-3],
}

# Hatching patterns and colors for geological units
ROCK_HATCH_PATTERNS: dict[str, tuple[str, tuple[float, float, float]]] = {
    "hard rock": (".+++++.", (0.25, 0.5, 0.5)),
    "igneous rock": (".o.o.", (1.0, 1.0, 1.0)),
    "duricrust": ("+.+", (1.0, 0.2, 0.36)),
    "gravel": ("oO", (0.75, 0.86, 0.12)),
    "sand": ("....", (0.23, 0.36, 0.45)),
    "conglomerate": (".O.", (0.55, 0.0, 0.36)),
    "dolomite": (".-.", (0.0, 0.75, 0.23)),
    "limestone": ("//.", (0.52, 0.23, 0.125)),
    "permafrost": ("o.", (0.2, 0.26, 0.75)),
    "metamorphic rock": ("*o.", (0.2, 0.2, 0.3)),
    "tills": ("-.", (0.7, 0.6, 0.9)),
    "standstone ": ("..", (0.5, 0.6, 0.9)),
    "lignite coal": ("+/.", (0.5, 0.5, 0.4)),
    "coal": ("*.", (0.8, 0.9, 0.0)),
    "shale": ("=", (0.0, 0.0, 0.7)),
    "clay": ("=.", (0.9, 0.8, 0.8)),
    "saprolite": ("*/", (0.3, 1.2, 0.4)),
    "sedimentary rock": ("...", (0.25, 0.0, 0.25)),
    "fresh water": (".-.", (0.0, 1.0, 0.2)),
    "salt water": ("o.-", (0.2, 1.0, 0.2)),
    "massive sulphide": (".+O", (1.0, 0.5, 0.5)),
    "sea water": (".--", (0.0, 1.0, 0.0)),
    "ore minerals": ("--|", (0.8, 0.2, 0.2)),
    "graphite": (".++.", (0.2, 0.7, 0.7)),
}


class RockProperties:
    """
    Access and query geological rock metadata:
      * resistivity ranges
      * plotting hatch patterns

    Examples
    --------
    >>> rp = RockProperties()
    >>> rp.get_resistivity("shale")
    [50.12, 32.0]
    >>> pattern, color = rp.get_pattern("shale")
    "=", (0.0, 0.0, 0.7)
    """

    def __init__(self):
        self._ranges = GEO_ROCK_RESISTIVITY
        self._patterns = ROCK_HATCH_PATTERNS

    @property
    def resistivity_ranges(self) -> dict[str, list[float]]:
        """
        All resistivity ranges keyed by rock name.
        """
        return self._ranges.copy()

    @property
    def hatch_patterns(
        self,
    ) -> dict[str, tuple[str, tuple[float, float, float]]]:
        """
        All hatch patterns keyed by rock name.
        """
        return self._patterns.copy()

    def get_resistivity(self, rock: str) -> list[float]:
        """
        Return [max, min] resistivity for a given rock.

        Raises KeyError if rock not found.
        """
        return self._ranges[rock]

    def get_pattern(self, rock: str) -> tuple[str, tuple[float, float, float]]:
        """
        Return (hatch, color) tuple for a given rock.

        Raises KeyError if rock not found.
        """
        return self._patterns[rock]

    def find_matching_rocks(self, keyword: str) -> list[str]:
        """
        Return list of rock names containing keyword (case‑insensitive).
        """
        kw = keyword.lower()
        return [r for r in self._ranges if kw in r.lower()]
