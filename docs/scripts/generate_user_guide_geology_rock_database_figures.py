"""Generate the rock-database figure for the geology rock_database guide."""

from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from pycsamt.geology import RockDatabase, RockEntry  # noqa: E402

IMAGES = ROOT / "docs/source/images/user_guide/geology"

# Same three-entry regional table built directly on the page.
_REGIONAL_ENTRIES = [
    RockEntry(name="Laterite", rho_min=80, rho_max=600, color="#B5651D", source="Site report 2024"),
    RockEntry(name="Saprolite", rho_min=20, rho_max=300, color="#C9A66B", source="Site report 2024"),
    RockEntry(name="Fresh basement", rho_min=3000, rho_max=200000, color="#4A4A4A", source="Site report 2024"),
]

_TEST_RHOS = [45.0, 150.0, 5000.0]


def make_rock_database_source_comparison() -> None:
    """Write the "same resistivity, two tables" comparison figure.

    Classifies three resistivity values against both the built-in
    ``RockDatabase.default()`` and the small regional table built
    directly on the page, showing that swapping the active table
    changes the classified name even though the query value never
    changes.
    """
    db_default = RockDatabase.default()
    db_regional = RockDatabase(_REGIONAL_ENTRIES)

    fig, ax = plt.subplots(figsize=(9, 4.2))

    y_regional, y_default = 0.0, 1.0

    for e in _REGIONAL_ENTRIES:
        ax.fill_between(
            [e.rho_min, e.rho_max], y_regional - 0.18, y_regional + 0.18,
            color=e.color, edgecolor="0.25", linewidth=0.6, alpha=0.9,
        )
        ax.text(
            e.rho_mid, y_regional, e.name, ha="center", va="center", fontsize=8,
            color="white" if e.name != "Laterite" else "black", fontweight="bold",
        )

    for rho in _TEST_RHOS:
        d = db_default.classify(rho)
        ax.plot(
            [rho, rho], [y_regional + 0.22, y_default - 0.06],
            color="0.4", linestyle=":", linewidth=1.0, zorder=1,
        )
        ax.scatter([rho], [y_default], color=d.color, edgecolor="0.2", s=140, zorder=3)
        ax.text(rho, y_default + 0.16, f"{rho:.0f} Ohm.m", ha="center", fontsize=8, color="0.25")
        ax.text(rho, y_default - 0.20, d.name, ha="center", fontsize=8, color="0.15")

    ax.set_xscale("log")
    ax.set_xlim(10, 3e5)
    ax.set_ylim(-0.4, 1.5)
    ax.set_yticks([y_regional, y_default])
    ax.set_yticklabels(
        ["Regional\n(3-entry file:// table)", "Built-in default()\n(49 entries)"],
        fontsize=9,
    )
    ax.set_xlabel(r"Resistivity ($\Omega\,\mathrm{m}$, log scale)")
    ax.set_title("The same resistivity value, classified against two different tables")
    ax.grid(axis="x", which="both", alpha=0.25)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)

    fig.tight_layout()
    fig.savefig(
        IMAGES / "rock_database_source_comparison.png", dpi=200, bbox_inches="tight"
    )
    plt.close(fig)


if __name__ == "__main__":
    IMAGES.mkdir(parents=True, exist_ok=True)
    make_rock_database_source_comparison()
