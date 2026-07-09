.. _extending:

Extending pyCSAMT
=================

When you write your own class on top of pyCSAMT — a reader for a custom
instrument format, a processing result, a station container — inherit
from the package base objects instead of ``object``. One line of
inheritance buys readable ``repr``, light serialization, the vectorized
MT math, and interoperability with the EDI/transfer-function machinery
used everywhere else in the package.

The Object Hierarchy
--------------------

:class:`~pycsamt.api.PyCSAMTObject`
    The root. Dependency-light readable ``__repr__``/``__str__`` and
    shallow introspection. Not numerical — result and config
    dataclasses all over the package (agents, IoT, pipeline) subclass
    it. Customize the summary with ``__repr_fields__`` /
    ``__repr_exclude__``.

:class:`~pycsamt.core.CoreObject`
    Adds safe array-aware display (large arrays are summarized by
    shape, never printed), an overridable ``_repr_fields`` hook,
    ``as_dict()`` shallow serialization, and ``summary()``.

:class:`~pycsamt.core.MTBase`
    The MT computational layer: numerically safe, vectorized NumPy
    helpers with broadcast semantics —

    * ``omega(f)`` — angular frequency;
    * ``rho_phase_from_z(z, f, phase_unit="deg")`` — apparent
      resistivity and phase from impedance;
    * ``z_from_rho_phase(rho, phase, f)`` — the inverse;
    * ``determinant_z(z)`` / ``rho_phase_from_det(z, f)`` — rotation
      invariants for ``(..., 2, 2)`` tensors;
    * ``rotate_impedance(z, theta_deg)`` — tensor rotation;
    * ``skin_depth(f, rho)`` — penetration depth in metres.

    Operations are element-wise: pass one component
    (``z[:, 0, 1]`` with ``freq``) or broadcast the frequency yourself
    (``freq[:, None, None]``) for full-tensor input.

:class:`~pycsamt.core.TFBundle`
    The neutral transfer-function payload — ``freq``, ``z (n, 2, 2)``,
    ``z_err``, ``tipper``, ``rho``, station metadata. It is the lingua
    franca between formats, backends, and your code: anything that can
    produce or consume a bundle interoperates with the rest of pyCSAMT.

:class:`~pycsamt.core.mixins.BundleMixin` / :class:`~pycsamt.core.mixins.BundleContainerMixin`
    Contracts on top of the bundle: implement ``to_bundle()`` (and
    optionally ``from_bundle()``) and the mixin supplies
    ``ensure_station_name``, ``to_edi()``/``from_edi()`` dispatch, and
    — for containers — ``iter_bundles()`` and ``to_edi_collection()``.

Choosing A Base
---------------

* A config or result object that only needs to *print well* →
  ``PyCSAMTObject`` (or ``CoreObject`` when it holds arrays).
* A class that computes with impedances, frequencies, or resistivities
  → ``MTBase``.
* A class that represents station data other pyCSAMT tools should be
  able to consume → add ``BundleMixin`` and implement ``to_bundle()``.

A Worked Example
----------------

A single-station reader for a custom format, plugged into everything at
once:

.. code-block:: python

   import numpy as np
   from pycsamt.core import MTBase, TFBundle
   from pycsamt.core.mixins import BundleMixin


   class MySounding(MTBase, BundleMixin):
       """One station from a custom instrument format."""

       def __init__(self, station, freq, z, z_err=None):
           self.station = station
           self.freq = np.asarray(freq, float)
           self.z = np.asarray(z, complex).reshape(-1, 2, 2)
           self.z_err = z_err

       # BundleMixin contract ------------------------------------
       def to_bundle(self):
           return TFBundle(freq=self.freq, z=self.z, z_err=self.z_err,
                           station=self.station)

       @classmethod
       def from_bundle(cls, bundle):
           return cls(bundle.station, bundle.freq, bundle.z, bundle.z_err)

Everything below comes from the bases — none of it is defined on
``MySounding``:

.. code-block:: python

   s = MySounding("S01", freq, z)

   repr(s)               # arrays summarized by shape, never dumped
   s.as_dict()           # {'station': ..., 'freq': ..., 'z': ..., 'z_err': ...}

   # MT computations (MTBase)
   rho_xy, phase_xy = s.rho_phase_from_z(s.z[:, 0, 1], s.freq)
   depth = s.skin_depth(s.freq, rho_xy)
   z_rot = s.rotate_impedance(s.z, 30.0)

   # transfer-function interop (BundleMixin)
   bundle = s.to_bundle()
   again = MySounding.from_bundle(bundle)

Format adapters are an explicit extension point: ``to_edi()`` and
``from_edi()`` dispatch through the adapter registry, so register a
callable for your kind (or rely on
:func:`~pycsamt.core.base.pick_adapter_key`, which infers ``'edi'``,
``'avg'``, or ``'j'`` from class and module names):

.. code-block:: python

   from pycsamt.core import register_adapter

   register_adapter("mysounding", lambda obj, **kw: obj.to_bundle())
   out = s.to_edi(key="mysounding")

Unit Conventions
----------------

:class:`~pycsamt.core.MTBase` documents both apparent-resistivity
conventions and this matters when you extend it:

* **SI** (E in V/m, H in A/m): :math:`\rho_a = |Z|^2 / (\mu_0 \omega)`.
* **Field units** (E in mV/km, B in nT): :math:`\rho_a \approx (0.2 / f)\,|Z|^2`
  — the Zonge-style convention.

The EDI files bundled with pyCSAMT (and most field EDIs) carry
impedances in **field units**. If your subclass ingests raw instrument
data, decide the convention at the boundary and keep it consistent —
mixing the two inflates or deflates :math:`\rho_a` by orders of
magnitude.

Rules For New Public Classes
----------------------------

If your class ships inside pyCSAMT (rather than in your own project):

* follow the numpydoc conventions in :doc:`docstring_style`;
* respect the public-API rules in :doc:`api_policy` — public means
  exported, documented, and tested;
* put tests next to the subpackage (``pycsamt/<subpackage>/tests/``),
  including a ``to_bundle()``/``from_bundle()`` round-trip when the
  bundle contract is implemented.
