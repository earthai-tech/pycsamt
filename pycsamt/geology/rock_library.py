# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Built-in rock/fluid resistivity table for :mod:`pycsamt.geology.lithology`.

This module holds *only data*: :data:`BUILTIN_ROCKS`, a flat list of
``(name, rho_min, rho_max, color, description, code, source)`` tuples
consumed by :meth:`pycsamt.geology.lithology.RockDatabase.default`. It carries
no classification logic, so it can grow -- new rows, new subsections, a
locally curated regional variant kept alongside it -- without touching
:mod:`pycsamt.geology.lithology` at all. Contributors extending the table for
a specific commodity, climate, or region should add rows here, grouped under
a labelled subsection with a literature citation in ``source``, rather than
editing the classification engine.

Resistivity ranges are ohm metres, linear scale, drawn from widely cited
compilations rather than from any single fetchable dataset -- no public,
machine-readable API currently maps rock name directly to a resistivity
range (see :mod:`pycsamt.geology.rock_providers` for how a project-specific
or organisation-internal source can be plugged in instead). Primary
references:

* Palacky, G.J. (1988), "Resistivity characteristics of geologic targets",
  in Nabighian, M.N. (ed.), *Electromagnetic Methods in Applied Geophysics*,
  Vol. 1, Society of Exploration Geophysicists.
* Telford, W.M., Geldart, L.P., Sheriff, R.E. (1990), *Applied Geophysics*,
  2nd ed., Cambridge University Press.
* Keller, G.V., Frischknecht, F.C. (1966), *Electrical Methods in
  Geophysical Prospecting*, Pergamon Press.
* Slichter, L.B., Telkes, M. (1942), "Electrical properties of rocks and
  minerals", in *Handbook of Physical Constants*, Geological Society of
  America.

Ranges are reproduced as round, order-of-magnitude-faithful brackets
consistent with these sources rather than as falsely precise per-sample
values -- individual field measurements routinely span the full width of a
bracket and beyond, which is the same reason
:doc:`/user_guide/interpretation/lithology` treats classification as a
hypothesis rather than a measurement.
"""

from __future__ import annotations

__all__ = ["BUILTIN_ROCKS"]

_PALACKY = "Palacky (1988)"
_TELFORD = "Telford, Geldart & Sheriff (1990)"
_KELLER = "Keller & Frischknecht (1966)"
_SLICHTER = "Slichter & Telkes (1942)"

BUILTIN_ROCKS: list[tuple[str, float, float, str, str, int, str]] = [
    # (name, rho_min, rho_max, color, description, code, source)
    # ------------------------------------------------------------------
    # Highly conductive
    # ------------------------------------------------------------------
    (
        "Sulfide ore body",
        0.001,
        0.1,
        "#2C3E50",
        "Massive sulfides, pyrite, chalcopyrite",
        1,
        f"{_TELFORD}; {_SLICHTER}",
    ),
    (
        "Graphite / coal",
        0.001,
        1.0,
        "#17202A",
        "Carbonaceous shale, graphitic schist",
        2,
        f"{_PALACKY}; {_TELFORD}",
    ),
    (
        "Saline water / brine",
        0.05,
        1.0,
        "#1ABC9C",
        "Saturated halite or formation brine",
        3,
        f"{_PALACKY}; {_KELLER}",
    ),
    # ------------------------------------------------------------------
    # Conductive
    # ------------------------------------------------------------------
    (
        "Clay",
        1,
        20,
        "#A9780C",
        "Smectite, illite, kaolinite clays",
        4,
        f"{_PALACKY}; {_TELFORD}",
    ),
    (
        "Alluvium (wet)",
        1,
        50,
        "#F4D03F",
        "Water-saturated unconsolidated sediment",
        5,
        _PALACKY,
    ),
    (
        "Aquifer",
        5,
        200,
        "#27AE60",
        "Fractured / porous water-bearing layer",
        6,
        _PALACKY,
    ),
    # ------------------------------------------------------------------
    # Moderate
    # ------------------------------------------------------------------
    (
        "Fractured zone",
        10,
        500,
        "#28B463",
        "Open fractures, fault zones",
        7,
        _PALACKY,
    ),
    (
        "Sand (wet)",
        20,
        200,
        "#F5CBA7",
        "Saturated sand / gravel",
        8,
        _TELFORD,
    ),
    (
        "Shale",
        5,
        100,
        "#7F8C8D",
        "Compacted argillaceous sediment",
        9,
        f"{_PALACKY}; {_TELFORD}",
    ),
    (
        "Granite (weathered)",
        50,
        2000,
        "#5DADE2",
        "MWG — moderate-to-strong weathering",
        10,
        _PALACKY,
    ),
    (
        "Basalt (weathered)",
        10,
        1000,
        "#85C1E9",
        "Tropical weathering profile on basalt",
        11,
        _PALACKY,
    ),
    # ------------------------------------------------------------------
    # Resistive
    # ------------------------------------------------------------------
    (
        "Sand (dry)",
        200,
        3000,
        "#F9E4B7",
        "Dry aeolian or vadose-zone sand",
        12,
        _TELFORD,
    ),
    (
        "Sandstone",
        50,
        5000,
        "#CA6F1E",
        "Consolidated siliciclastic",
        13,
        f"{_PALACKY}; {_TELFORD}",
    ),
    (
        "Limestone",
        500,
        10_000,
        "#D5D8DC",
        "Carbonate platform, reef limestone",
        14,
        f"{_PALACKY}; {_TELFORD}",
    ),
    (
        "Dolomite",
        500,
        20_000,
        "#BFC9CA",
        "Dolostones, evaporite-associated",
        15,
        _TELFORD,
    ),
    (
        "Schist",
        200,
        5000,
        "#8E44AD",
        "Pelitic metamorphic",
        16,
        _TELFORD,
    ),
    (
        "Marble",
        500,
        10_000,
        "#D7BDE2",
        "Calcareous metamorphic",
        17,
        _TELFORD,
    ),
    (
        "Gneiss",
        500,
        50_000,
        "#5B2C6F",
        "High-grade metamorphic basement",
        18,
        _TELFORD,
    ),
    # ------------------------------------------------------------------
    # Highly resistive
    # ------------------------------------------------------------------
    (
        "Granite (fresh)",
        5000,
        1_000_000,
        "#1A5276",
        "Unweathered plutonic rock",
        19,
        f"{_PALACKY}; {_TELFORD}",
    ),
    (
        "Basalt (fresh)",
        1000,
        100_000,
        "#154360",
        "Unweathered volcanic",
        20,
        _PALACKY,
    ),
    (
        "Gabbro",
        1000,
        100_000,
        "#0B2641",
        "Mafic intrusive",
        21,
        _TELFORD,
    ),
    (
        "Quartzite",
        1000,
        100_000,
        "#E8DAEF",
        "High-grade silicate metamorphic",
        22,
        _TELFORD,
    ),
    (
        "Igneous (basement)",
        3000,
        1_000_000,
        "#1B2631",
        "Undifferentiated crystalline basement",
        23,
        _PALACKY,
    ),
    (
        "Evaporite (dry)",
        1000,
        1_000_000,
        "#F0F3F4",
        "Dry halite, anhydrite",
        24,
        f"{_PALACKY}; {_TELFORD}",
    ),
    (
        "Air / void",
        1e6,
        1e12,
        "#FFFFFF",
        "Air-filled cavities, dry caves",
        25,
        _KELLER,
    ),
    # ------------------------------------------------------------------
    # Water and pore fluids (distinct from the generic brine entry above)
    # ------------------------------------------------------------------
    (
        "Freshwater (groundwater)",
        10,
        100,
        "#5DADE2",
        "Typical potable meteoric groundwater",
        26,
        f"{_KELLER}; {_PALACKY}",
    ),
    (
        "Brackish water",
        1,
        10,
        "#48C9B0",
        "Estuarine or partially mixed pore water",
        27,
        _KELLER,
    ),
    (
        "Seawater",
        0.15,
        0.5,
        "#117864",
        "Normal-salinity marine water, ~35 g/L",
        28,
        _KELLER,
    ),
    # ------------------------------------------------------------------
    # Named ore minerals (each individually spans several decades and
    # overlaps every other entry in this subsection; see the discussion
    # of the generic "Sulfide ore body" entry in the lithology guide)
    # ------------------------------------------------------------------
    (
        "Pyrite",
        0.01,
        100,
        "#212F3C",
        "FeS2; disseminated to massive",
        29,
        f"{_TELFORD}; {_SLICHTER}",
    ),
    (
        "Magnetite",
        0.001,
        1000,
        "#1C2833",
        "Fe3O4; disseminated to massive",
        30,
        f"{_TELFORD}; {_SLICHTER}",
    ),
    (
        "Galena",
        0.001,
        300,
        "#34495E",
        "PbS; vein to massive",
        31,
        f"{_TELFORD}; {_SLICHTER}",
    ),
    (
        "Hematite",
        1,
        100_000,
        "#943126",
        "Fe2O3; oxide-facies iron formation",
        32,
        f"{_TELFORD}; {_SLICHTER}",
    ),
    # ------------------------------------------------------------------
    # Cryosphere
    # ------------------------------------------------------------------
    (
        "Permafrost",
        1000,
        100_000,
        "#AED6F1",
        "Perennially frozen ground; ice-bonded sediment",
        33,
        _KELLER,
    ),
    (
        "Glacier ice",
        10_000,
        1e8,
        "#EBF5FB",
        "Temperate to polar glacier ice",
        34,
        _KELLER,
    ),
    # ------------------------------------------------------------------
    # Additional sediments and sedimentary rocks
    # ------------------------------------------------------------------
    (
        "Silt (wet)",
        10,
        200,
        "#D7CCC8",
        "Saturated fine-grained clastic sediment",
        35,
        _PALACKY,
    ),
    (
        "Gravel (dry)",
        1000,
        10_000,
        "#E5C9A6",
        "Unsaturated coarse clastic sediment",
        36,
        _TELFORD,
    ),
    (
        "Glacial till",
        100,
        3000,
        "#9E8465",
        "Unsorted glacial diamict",
        37,
        _PALACKY,
    ),
    (
        "Mudstone / claystone",
        5,
        100,
        "#6E6560",
        "Fine-grained, poorly fissile argillite",
        38,
        _TELFORD,
    ),
    (
        "Siltstone",
        20,
        500,
        "#8D7B68",
        "Lithified silt-grade clastic",
        39,
        _TELFORD,
    ),
    (
        "Conglomerate",
        1000,
        10_000,
        "#B9770E",
        "Coarse, poorly sorted clastic",
        40,
        _TELFORD,
    ),
    (
        "Coal",
        10,
        1000,
        "#1C2833",
        "Sub-bituminous to anthracite; excludes graphitized seams",
        41,
        _TELFORD,
    ),
    # ------------------------------------------------------------------
    # Additional igneous rocks
    # ------------------------------------------------------------------
    (
        "Andesite",
        100,
        5000,
        "#7FB3D5",
        "Intermediate volcanic",
        42,
        _TELFORD,
    ),
    (
        "Rhyolite",
        100,
        10_000,
        "#AED6F1",
        "Felsic volcanic",
        43,
        _TELFORD,
    ),
    (
        "Diorite",
        1000,
        10_000,
        "#2874A6",
        "Intermediate intrusive",
        44,
        _TELFORD,
    ),
    (
        "Peridotite / dunite",
        1000,
        100_000,
        "#0E4B39",
        "Ultramafic intrusive",
        45,
        _TELFORD,
    ),
    (
        "Serpentinite",
        100,
        10_000,
        "#186A3B",
        "Hydrothermally altered ultramafic rock",
        46,
        _PALACKY,
    ),
    # ------------------------------------------------------------------
    # Additional metamorphic rocks
    # ------------------------------------------------------------------
    (
        "Slate",
        500,
        5000,
        "#5D6D7E",
        "Low-grade pelitic metamorphic, well-cleaved",
        47,
        _TELFORD,
    ),
    (
        "Phyllite",
        300,
        3000,
        "#707B7C",
        "Low-to-medium-grade pelitic metamorphic",
        48,
        _TELFORD,
    ),
    (
        "Amphibolite",
        1000,
        10_000,
        "#4A235A",
        "Mafic-derived medium-to-high-grade metamorphic",
        49,
        _TELFORD,
    ),
]
