.. _geology_concepts:

Package concepts
=================

Every class in :mod:`pycsamt.geology` -- rock entries, boreholes, fault
traces -- follows the same small set of conventions. They are worth
learning once, here, rather than re-explaining on every page that
follows.

Leaf values and containers
---------------------------

Classes in this package come in two shapes. A **leaf value object**
describes one fact -- one rock/resistivity range, one borehole depth
interval, one field measurement -- as a ``@dataclass`` that validates
itself the moment it is constructed:

.. code-block:: pycon

   >>> from pycsamt.geology import RockEntry
   >>> entry = RockEntry(
   ...     name="Granite (weathered)",
   ...     rho_min=50,
   ...     rho_max=2000,
   ...     source="Palacky (1988)",
   ... )
   >>> entry
   RockEntry(name='Granite (weathered)', rho_min=50, rho_max=2000, color='#AAAAAA', description='', code=0, source='Palacky (1988)')

A **container object** collects many leaf values that belong together
-- every rock entry in a table, every depth interval in one borehole,
every structural observation on one profile. Containers are plain
classes with an explicit ``__init__``, a ``from_csv`` classmethod, and
a query API, rather than dataclasses:

.. code-block:: pycon

   >>> from pycsamt.geology import RockDatabase
   >>> db = RockDatabase.default()
   >>> db
   RockDatabase(49 entries)
   >>> db.entries[0]
   RockEntry(name='Sulfide ore body', rho_min=0.001, rho_max=0.1, color='#2C3E50', description='Massive sulfides, pyrite, chalcopyrite', code=1, source='Telford, Geldart & Sheriff (1990); Slichter & Telkes (1942)')

The same split recurs throughout the package:
:class:`~pycsamt.geology.borehole.Interval` (leaf) inside
:class:`~pycsamt.geology.Borehole` (container);
:class:`~pycsamt.geology.StructuralMeasurement`,
:class:`~pycsamt.geology.LinearMeasurement`, and
:class:`~pycsamt.geology.FaultTrace` (leaves) inside
:class:`~pycsamt.geology.StructuralModel` (container). Once you know
which shape a class is, you already know most of its API.

Repr, dictionaries, and metadata
----------------------------------

Every class inherits :class:`~pycsamt.api.property.PyCSAMTObject`,
which is what produced the compact, readable ``repr()`` output above
without either class writing its own ``__repr__``. It also provides
:meth:`~pycsamt.api.property.PyCSAMTObject.to_dict`, a plain-``dict``
view that is easy to serialize or hand to :mod:`pandas`:

.. code-block:: pycon

   >>> entry.to_dict()
   {'name': 'Granite (weathered)', 'rho_min': 50, 'rho_max': 2000, 'color': '#AAAAAA', 'description': '', 'code': 0, 'source': 'Palacky (1988)'}

Container objects that carry provenance --
:class:`~pycsamt.geology.RockDatabase` and
:class:`~pycsamt.geology.StructuralModel` -- additionally mix in
:class:`~pycsamt.api.property.MetadataMixin` for a free-form
``metadata`` dictionary. ``RockDatabase.default()`` already stamps its
own origin into it, and :meth:`~pycsamt.api.property.MetadataMixin.update_metadata`
adds to it without disturbing what is already there:

.. code-block:: pycon

   >>> db.update_metadata(project="demo-survey", reviewer="D.K.")
   RockDatabase(49 entries)
   >>> db.metadata_dict()
   {'origin': 'default', 'project': 'demo-survey', 'reviewer': 'D.K.'}

Leaf value objects do not carry metadata -- there is one fact, and
provenance for the whole collection belongs on the container, not
repeated on every entry.

Editing an object: what ``clone()`` and ``update()`` actually check
----------------------------------------------------------------------

:meth:`~pycsamt.api.property.PyCSAMTObject.clone` returns a modified
copy and :meth:`~pycsamt.api.property.PyCSAMTObject.update` modifies
an object in place; both accept the same keyword overrides and both
call :meth:`~pycsamt.api.property.PyCSAMTObject.validate` afterwards.
The important detail is what ``validate()`` actually does: the base
implementation is a no-op, so whether an edit is checked depends on
whether the specific class *overrides* ``validate()`` -- and not every
class does, because not every class has the same amount to check.

:class:`~pycsamt.geology.StructuralMeasurement` has real cross-field
structure to protect (``dip_direction_deg`` must stay consistent with
``strike_deg``), so it overrides ``validate()`` to re-run the same
check its constructor uses. Editing it with an inconsistent value is
caught either way:

.. code-block:: pycon

   >>> from pycsamt.geology import StructuralMeasurement
   >>> bedding = StructuralMeasurement(
   ...     x=100.0, strike_deg=45.0, dip_deg=30.0,
   ...     dip_direction_deg=135.0, kind="bedding",
   ... )
   >>> bedding.clone(dip_direction_deg=0.0)
   Traceback (most recent call last):
       ...
   ValueError: dip_direction_deg (0.0) is not within 20.0 deg of strike_deg (45.0) +/- 90 -- check for a transposed strike/dip-direction reading.

:class:`~pycsamt.geology.borehole.Interval` validates a single
constraint (``bottom > top``) in ``__post_init__``, but does not
override ``validate()`` -- so construction is checked, and an edit via
``clone()``/``update()`` is not:

.. code-block:: pycon

   >>> from pycsamt.geology import Interval
   >>> Interval(top=20.0, bottom=10.0, lithology="Sand")
   Traceback (most recent call last):
       ...
   ValueError: Interval bottom (10.0) must be > top (20.0).
   >>> good = Interval(top=10.0, bottom=20.0, lithology="Sand")
   >>> edited = good.clone(bottom=5.0)
   >>> edited.thickness
   -5.0

``RockEntry`` goes further still and never validates at all, at
construction or otherwise -- ``rho_min`` greater than ``rho_max`` is
accepted silently, since nothing downstream assumes an order between
them beyond what :meth:`~pycsamt.geology.RockEntry.contains` computes
correctly regardless.

.. warning::

   Do not rely on ``clone()``/``update()`` to catch every mistake.
   Whether an edited object gets re-validated is a per-class property
   of whether that class overrides ``validate()``, not a guarantee of
   the base ``PyCSAMTObject`` API. When editing a value in place --
   particularly ``Interval`` bounds or any leaf field without an
   obvious cross-check -- inspect the result, or reconstruct the
   object directly so ``__post_init__`` runs again.

Position is always profile-relative
-------------------------------------

Every class that has a location in this package places it with a
single ``x`` in metres along a 2-D profile -- never a latitude/
longitude pair:

.. code-block:: pycon

   >>> from pycsamt.geology import Borehole, FaultTrace
   >>> bh = Borehole(
   ...     name="BH01", x=1050.0,
   ...     intervals=[
   ...         Interval(top=0, bottom=5, lithology="Topsoil"),
   ...         Interval(top=5, bottom=40, lithology="Sand"),
   ...     ],
   ... )
   >>> bh
   Borehole('BH01', x=1050.0 m, 2 intervals, depth=40.0 m)
   >>> FaultTrace(x=500.0, dip_deg=70.0, downthrown_side="right")
   FaultTrace(x=500.0 m, dip=70 deg, down=right, throw=?)

This matches how a :class:`~pycsamt.interp.ResistivityModel`'s
``x_centers`` are already expressed, so a borehole, fault trace, or
structural measurement lines up directly against a section plot
without a coordinate lookup. Real-world placement -- the latitude/
longitude or projected easting/northing behind that ``x`` -- is
:mod:`pycsamt.site`/:mod:`pycsamt.gis`'s responsibility, not this
package's; :mod:`pycsamt.geology` only ever needs to know where
something sits *along the line*.

Where to go next
-------------------

:doc:`rock_database` covers :class:`~pycsamt.geology.RockDatabase` in
full, including where :data:`~pycsamt.geology.rock_library.BUILTIN_ROCKS`
comes from and how to source a table from elsewhere.
:doc:`borehole` and :doc:`structural` cover the two other leaf/
container families introduced above. For applying any of this to a
real resistivity section, see :doc:`../interpretation/index`.
