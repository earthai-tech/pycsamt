# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""Registry of EMTF primary and derived electromagnetic data types.

The registry mirrors the semantic role of the ``DATATYPES`` directory in
EMTF FCU v4.1 while remaining independent of XML serialization.  Entries are
Python value objects so future formats (including passive airborne EM) can
reuse the same scientific vocabulary without importing the XML reader.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from ..api.property import PyCSAMTObject

__all__ = [
    "DataTypeDefinition",
    "register_emtf_datatype",
    "ensure_emtf_datatype_registered",
    "get_emtf_datatype",
    "list_emtf_datatypes",
]

_VALID_DATA_KINDS = frozenset({"real", "complex"})
_VALID_INTENTIONS = frozenset({"primary", "derived"})


@dataclass(frozen=True, repr=False)
class DataTypeDefinition(PyCSAMTObject):
    """Describe one electromagnetic transfer-function data type.

    Parameters
    ----------
    name : str
        Short EMTF/FCU code, for example ``"Z"`` or ``"T"``.
    tag : str
        Stable semantic key, for example ``"impedance"``.
    data_kind : {"real", "complex"}
        Numerical representation of the data.
    input_kind, output_kind : str or None
        Broad channel families used by EMTF (typically ``"H"`` or ``"E"``).
        Scalar derived products may omit these values.
    units : str or None
        Native units recorded by the format definition.
    intention : {"primary", "derived"}
        Whether the quantity is measured/estimated as a primary TF or derived
        from another response.
    description : str
        Human-readable description.
    derived_from : str or None
        Semantic tag of the parent quantity for derived products.
    aliases : tuple of str
        Alternate historical codes accepted for lookup.  These are useful
        because the 2020 EMTF paper and FCU v4.1 differ for a few names
        (for example Q/P versus ZI/TI).
    see_also : tuple of str
        Related semantic tags.
    """

    name: str
    tag: str
    data_kind: str
    input_kind: str | None = None
    output_kind: str | None = None
    units: str | None = None
    intention: str = "primary"
    description: str = ""
    derived_from: str | None = None
    aliases: tuple[str, ...] = field(default_factory=tuple)
    see_also: tuple[str, ...] = field(default_factory=tuple)

    def __post_init__(self) -> None:
        name = str(self.name).strip().upper()
        tag = str(self.tag).strip().lower()
        data_kind = str(self.data_kind).strip().lower()
        intention = str(self.intention).strip().lower()
        if not name:
            raise ValueError("EMTF data-type name must be non-empty")
        if not tag:
            raise ValueError("EMTF data-type tag must be non-empty")
        if data_kind not in _VALID_DATA_KINDS:
            raise ValueError(
                "data_kind must be one of " f"{sorted(_VALID_DATA_KINDS)}"
            )
        if intention not in _VALID_INTENTIONS:
            raise ValueError(
                "intention must be one of " f"{sorted(_VALID_INTENTIONS)}"
            )

        object.__setattr__(self, "name", name)
        object.__setattr__(self, "tag", tag)
        object.__setattr__(self, "data_kind", data_kind)
        object.__setattr__(self, "intention", intention)
        object.__setattr__(
            self,
            "aliases",
            tuple(
                str(v).strip().upper()
                for v in self.aliases
                if str(v).strip()
            ),
        )
        object.__setattr__(
            self,
            "see_also",
            tuple(
                str(v).strip().lower()
                for v in self.see_also
                if str(v).strip()
            ),
        )

    @property
    def is_primary(self) -> bool:
        """Return ``True`` for a primary transfer-function type."""
        return self.intention == "primary"

    @property
    def is_derived(self) -> bool:
        """Return ``True`` for a derived data type."""
        return self.intention == "derived"


_BY_TAG: dict[str, DataTypeDefinition] = {}
_BY_NAME: dict[str, DataTypeDefinition] = {}


def register_emtf_datatype(
    definition: DataTypeDefinition,
    *,
    overwrite: bool = False,
) -> DataTypeDefinition:
    """Register and return an EMTF data-type definition.

    The semantic ``tag`` is the unique registry identity.  Short names can be
    shared by compatible variants (notably full and off-diagonal impedance),
    therefore the first short-name registration remains the default lookup
    unless ``overwrite=True`` is explicitly requested.
    """
    if not isinstance(definition, DataTypeDefinition):
        raise TypeError("definition must be a DataTypeDefinition")

    tag_key = definition.tag.lower()
    if tag_key in _BY_TAG and not overwrite:
        raise ValueError(f"EMTF data type already registered: {tag_key}")
    _BY_TAG[tag_key] = definition

    names = (definition.name, *definition.aliases)
    for name in names:
        key = name.upper()
        if overwrite or key not in _BY_NAME:
            _BY_NAME[key] = definition
    return definition


def ensure_emtf_datatype_registered(
    definition: DataTypeDefinition,
) -> DataTypeDefinition:
    """Idempotently register *definition*, tolerating a compatible repeat.

    :func:`register_emtf_datatype` raises whenever its ``tag`` is
    already registered, unless the caller passes ``overwrite=True`` --
    a reasonable default for one-shot registration, but awkward for a
    technology adapter's module-level registration function, which may
    legitimately run more than once in a process (a package re-import
    under test isolation, or more than one caller defensively calling
    the same ``register_<technology>_datatypes()``). This helper
    treats a definition already registered under the same
    ``tag``/``name`` pair as a no-op and returns the existing (first-
    registered) instance unchanged, while a genuine collision -- the
    same ``tag`` or ``name`` claimed by a materially different
    definition -- still raises rather than being silently papered
    over. :mod:`pycsamt.airborne.mobilemt` and
    :mod:`pycsamt.airborne.afmag` both register their derived
    datatypes through this helper for exactly that reason.

    Parameters
    ----------
    definition : DataTypeDefinition
        Definition to register.

    Returns
    -------
    DataTypeDefinition
        *definition* itself if newly registered, or the
        already-registered definition sharing its ``tag``.

    Raises
    ------
    ValueError
        If *definition*'s ``tag`` is already registered under a
        different ``name``, or its ``name`` is already registered
        (as a primary code or alias) under a different ``tag``.
    """
    existing = get_emtf_datatype(definition.tag)
    if existing is not None:
        if existing.name != definition.name:
            raise ValueError(
                f"{definition.tag!r} is already registered with code "
                f"{existing.name!r}"
            )
        return existing

    by_code = get_emtf_datatype(definition.name)
    if by_code is not None and by_code.tag != definition.tag:
        raise ValueError(
            f"EMTF datatype code {definition.name!r} is already "
            f"registered for {by_code.tag!r}"
        )
    return register_emtf_datatype(definition)


def get_emtf_datatype(key: str) -> DataTypeDefinition | None:
    """Return a registered definition by semantic tag, code, or alias."""
    if not isinstance(key, str) or not key.strip():
        return None
    raw = key.strip()
    return _BY_TAG.get(raw.lower()) or _BY_NAME.get(raw.upper())


def list_emtf_datatypes(
    *,
    intention: str | None = None,
) -> dict[str, DataTypeDefinition]:
    """Return a copy of registered data types keyed by semantic tag."""
    if intention is None:
        return dict(_BY_TAG)
    wanted = str(intention).strip().lower()
    if wanted not in _VALID_INTENTIONS:
        raise ValueError(
            "intention must be one of " f"{sorted(_VALID_INTENTIONS)}"
        )
    return {
        tag: definition
        for tag, definition in _BY_TAG.items()
        if definition.intention == wanted
    }


def _register_defaults() -> None:
    # Primary definitions follow EMTF FCU v4.1. Q/P are retained as aliases
    # for the interstation names used in the 2020 publication table.
    primary = (
        DataTypeDefinition(
            name="Z",
            tag="impedance",
            data_kind="complex",
            input_kind="H",
            output_kind="E",
            units="[mV/km]/[nT]",
            description="MT impedance",
        ),
        DataTypeDefinition(
            name="T",
            tag="tipper",
            data_kind="complex",
            input_kind="H",
            output_kind="H",
            units="[]",
            description="Vertical field transfer functions (tipper)",
        ),
        DataTypeDefinition(
            name="ZI",
            tag="interstation_impedance",
            data_kind="complex",
            input_kind="H",
            output_kind="E",
            units="[mV/km]/[nT]",
            description="MT interstation impedance",
            aliases=("Q",),
        ),
        DataTypeDefinition(
            name="TI",
            tag="interstation_transfer_functions",
            data_kind="complex",
            input_kind="H",
            output_kind="H",
            units="[]",
            description="Interstation magnetic transfer functions",
            aliases=("P",),
        ),
        DataTypeDefinition(
            name="Z",
            tag="off_diagonal_impedance",
            data_kind="complex",
            input_kind="H",
            output_kind="E",
            units="[mV/km]/[nT]",
            description="MT impedance - off-diagonal components only",
            see_also=("impedance",),
        ),
    )

    derived = (
        DataTypeDefinition(
            name="PT",
            tag="phase_tensor",
            data_kind="real",
            input_kind="H",
            output_kind="E",
            units="[]",
            intention="derived",
            description="MT phase tensor",
            derived_from="impedance",
        ),
        DataTypeDefinition(
            name="RHO",
            tag="apparent_resistivity",
            data_kind="real",
            input_kind="H",
            output_kind="E",
            units="[Ohm m]",
            intention="derived",
            description="Apparent resistivity computed from MT impedance",
            derived_from="impedance",
            see_also=("impedance_phase",),
        ),
        DataTypeDefinition(
            name="PHS",
            tag="impedance_phase",
            data_kind="real",
            input_kind="H",
            output_kind="E",
            units="degrees",
            intention="derived",
            description="Impedance phase in degrees",
            derived_from="impedance",
            see_also=("apparent_resistivity",),
        ),
        DataTypeDefinition(
            name="EFF",
            tag="effective_impedance",
            data_kind="complex",
            units="[mV/km]/[nT]",
            intention="derived",
            description="Square root of the impedance determinant",
            derived_from="impedance",
            aliases=("ZEFF",),
            see_also=("impedance_determinant",),
        ),
        DataTypeDefinition(
            name="DET",
            tag="impedance_determinant",
            data_kind="complex",
            units="([mV/km]/[nT])^2",
            intention="derived",
            description="Rotation-invariant determinant of MT impedance",
            derived_from="impedance",
            aliases=("ZDET",),
            see_also=("effective_impedance",),
        ),
        DataTypeDefinition(
            name="ZSTRIKE",
            tag="impedance_strike",
            data_kind="real",
            units="degrees",
            intention="derived",
            description="MT impedance strike",
            derived_from="impedance",
            see_also=("impedance_skew", "impedance_ellipticity"),
        ),
        DataTypeDefinition(
            name="ZSKEW",
            tag="impedance_skew",
            data_kind="real",
            units="[]",
            intention="derived",
            description="MT impedance skew",
            derived_from="impedance",
            see_also=("impedance_strike", "impedance_ellipticity"),
        ),
        DataTypeDefinition(
            name="ZELLIP",
            tag="impedance_ellipticity",
            data_kind="real",
            units="[]",
            intention="derived",
            description="MT impedance ellipticity",
            derived_from="impedance",
            see_also=("impedance_skew", "impedance_strike"),
        ),
        DataTypeDefinition(
            name="TIPMAG",
            tag="tipper_magnitude",
            data_kind="real",
            units="[]",
            intention="derived",
            description="Tipper magnitude",
            derived_from="tipper",
            see_also=("tipper_phase",),
        ),
        DataTypeDefinition(
            name="TIPPHS",
            tag="tipper_phase",
            data_kind="real",
            units="degrees",
            intention="derived",
            description="Tipper phase",
            derived_from="tipper",
            see_also=("tipper_magnitude",),
        ),
        DataTypeDefinition(
            name="INDMAG",
            tag="induction_arrow_magnitude",
            data_kind="complex",
            units="[]",
            intention="derived",
            description="Induction arrow magnitude",
            derived_from="tipper",
            see_also=("induction_arrow_angle",),
        ),
        DataTypeDefinition(
            name="INDANG",
            tag="induction_arrow_angle",
            data_kind="complex",
            units="degrees",
            intention="derived",
            description="Induction arrow angle (Parkinson convention)",
            derived_from="tipper",
            see_also=("induction_arrow_magnitude",),
        ),
        DataTypeDefinition(
            name="TSTRIKE",
            tag="tipper_strike",
            data_kind="real",
            units="degrees",
            intention="derived",
            description="Tipper strike",
            derived_from="tipper",
            see_also=("tipper_skew", "tipper_ellipticity"),
        ),
        DataTypeDefinition(
            name="TSKEW",
            tag="tipper_skew",
            data_kind="real",
            units="[]",
            intention="derived",
            description="Tipper skew",
            derived_from="tipper",
            see_also=("tipper_strike", "tipper_ellipticity"),
        ),
        DataTypeDefinition(
            name="TELLIP",
            tag="tipper_ellipticity",
            data_kind="real",
            units="[]",
            intention="derived",
            description="Tipper ellipticity",
            derived_from="tipper",
            see_also=("tipper_skew", "tipper_strike"),
        ),
    )

    for definition in (*primary, *derived):
        register_emtf_datatype(definition)


_register_defaults()
