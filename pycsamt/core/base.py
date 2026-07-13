# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""Foundational objects and transfer-function containers for pyCSAMT."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field, is_dataclass
from dataclasses import fields as dc_fields
from typing import (
    Any,
    Protocol,
)

import numpy as np

from ..api.property import PyCSAMTObject
from . import config as _core_config
from .config import StationNamePolicy, get_config

__all__ = [
    "TFBundle",
    "SupportsToBundle",
    "SupportsFromBundle",
    "ensure_station",
    "pick_adapter_key",
    "to_edi",
]


class CoreObject(PyCSAMTObject):
    r"""
    Minimal base class for PyCSAMT objects and mixins.

    The goal is to provide stable, readable ``__repr__`` and
    ``__str__`` implementations, plus a couple of helpers for
    light introspection and serialization that do *not* pull
    heavy dependencies.

    Notes
    -----
    * ``__repr__`` shows a compact field summary. Large arrays
      or containers are summarized safely (shape/len only).
    * Subclasses may override :meth:`_repr_fields` to control
      which attributes appear in the summary.
    * :meth:`as_dict` is shallow by default and avoids importing
      external libs. It aims to be robust, not exhaustive.

    Examples
    --------
    >>> class P(CoreObject):
    ...     def __init__(self):
    ...         self.name = "S001"
    ...         self.data = [1, 2, 3]
    ...
    >>> P()  # doctest: +ELLIPSIS
    P(name='S001', data=[len=3])

    See Also
    --------
    _repr_fields : Hook to pick fields for ``__repr__``.
    as_dict : Lightweight serialization helper.
    summary : Short single-line description.
    """

    def __repr__(self) -> str:  # noqa: D401
        name = self.__class__.__name__
        items = []
        for k in self._repr_fields():
            try:
                v = getattr(self, k)
            except Exception:  # pragma: no cover
                continue
            items.append(f"{k}={self._short(v)}")
        inside = ", ".join(items)
        return f"{name}({inside})"

    def __str__(self) -> str:  # noqa: D401
        return self.summary()

    def summary(self, *, max_fields: int = 6) -> str:
        r"""
        Return a short one-line description.

        Parameters
        ----------
        max_fields : int, optional
            Maximum number of fields to include.

        Returns
        -------
        str
            A compact summary string.
        """
        name = self.__class__.__name__
        fs = list(self._repr_fields())[:max_fields]
        items = []
        for k in fs:
            try:
                v = getattr(self, k)
            except Exception:
                continue
            items.append(f"{k}={self._short(v)}")
        more = ""
        if len(list(self._repr_fields())) > max_fields:
            more = ", ..."
        return f"{name}({', '.join(items)}{more})"

    def as_dict(
        self,
        *,
        public_only: bool = True,
        max_depth: int = 1,
    ) -> dict:
        r"""
        Serialize to a lightweight ``dict`` safely.

        Parameters
        ----------
        public_only : bool, optional
            If ``True``, include only public attributes (no
            leading underscore). Dataclasses always use their
            declared fields.
        max_depth : int, optional
            Maximum recursion depth for nested mappings and
            sequences. Use small values to avoid deep walks.

        Returns
        -------
        dict
            A JSON-friendly mapping (best-effort).
        """
        return self._to_dict(self, public_only=public_only, depth=max_depth)

    def _repr_fields(self) -> Iterable[str]:
        r"""
        Select attribute names to show in ``__repr__``.

        Returns
        -------
        Iterable[str]
            Field names to display. By default:
            * dataclasses → declared field names
            * objects     → public attributes in ``__dict__``
        """
        exclude = set(getattr(self, "__repr_exclude__", set()))
        if is_dataclass(self):
            return (f.name for f in dc_fields(self) if f.name not in exclude)
        if hasattr(self, "__dict__"):
            return (
                k
                for k in self.__dict__.keys()
                if not k.startswith("_") and k not in exclude
            )
        return ()  # pragma: no cover

    @staticmethod
    def _short(v: object) -> str:
        """Summarize values safely for ``__repr__``."""
        # primitives
        if v is None or isinstance(v, (bool, int, float, complex)):
            return repr(v)
        if isinstance(v, str):
            if len(v) > 32:
                return repr(v[:29] + "…")
            return repr(v)
        # numpy arrays (optional)
        try:  # pragma: no cover for environments w/o numpy
            import numpy as _np  # local import, no hard dep

            if isinstance(v, _np.ndarray):
                return f"ndarray(shape={tuple(v.shape)}, dtype={v.dtype})"
        except Exception:
            pass
        # mappings
        if isinstance(v, Mapping):
            keys = list(v.keys())
            head = ", ".join(map(str, keys[:3]))
            if len(keys) > 3:
                head += ", ..."
            return f"dict(len={len(keys)}, keys=[{head}])"
        # sequences (but not strings/bytes)
        if isinstance(v, (bytes, bytearray)):
            return f"bytes(len={len(v)})"
        if isinstance(v, (list, tuple, set)):
            n = len(v)
            # sample a few items, summarized
            sample = ", ".join(CoreObject._short(x) for x in list(v)[:3])
            if n > 3:
                sample += ", ..."
            typ = type(v).__name__
            return f"{typ}([{sample}])" if n else f"{typ}([])"
        # objects exposing shape/dtype
        try:
            shp = getattr(v, "shape", None)
            dt = getattr(v, "dtype", None)
            if shp is not None:
                return f"array_like(shape={shp}, dtype={dt})"
        except Exception:
            pass
        # fallback
        return repr(v)

    @classmethod
    def _to_dict(
        cls,
        obj: object,
        *,
        public_only: bool,
        depth: int,
    ) -> dict:
        """Shallow, safe dict conversion."""
        if depth < 0:
            return {"_": "max_depth"}
        out: dict = {}
        # dataclass
        if is_dataclass(obj):
            for f in dc_fields(obj):  # type: ignore
                k = f.name
                v = getattr(obj, k, None)
                out[k] = cls._coerce(v, public_only, depth - 1)
            return out
        # generic object with __dict__
        if hasattr(obj, "__dict__"):
            for k, v in obj.__dict__.items():
                if public_only and k.startswith("_"):
                    continue
                out[k] = cls._coerce(v, public_only, depth - 1)
            return out
        # fallback: empty
        return out

    @classmethod
    def _coerce(
        cls,
        v: object,
        public_only: bool,
        depth: int,
    ):
        """Coerce nested values for :meth:`as_dict`."""
        # primitives
        if v is None or isinstance(v, (bool, int, float, complex, str)):
            return v
        # numpy arrays → summary dict (no heavy stats)
        try:  # pragma: no cover
            import numpy as _np

            if isinstance(v, _np.ndarray):
                return {
                    "type": "ndarray",
                    "shape": tuple(v.shape),
                    "dtype": str(v.dtype),
                }
        except Exception:
            pass
        # mapping
        if isinstance(v, Mapping):
            out = {}
            for k, x in list(v.items())[:32]:
                out[str(k)] = cls._coerce(x, public_only, depth - 1)
            return out
        # sequence (not bytes)
        if isinstance(v, (list, tuple, set)):
            head = list(v)[:32]
            return [cls._coerce(x, public_only, depth - 1) for x in head]
        # object → nested dict (one level)
        return cls._to_dict(v, public_only=public_only, depth=depth - 1)


@dataclass
class TFBundle(CoreObject):
    r"""
    Lightweight, neutral payload for transfer functions.

    The bundle holds frequency-indexed arrays and optional site
    metadata. It is intentionally permissive about array types
    (e.g., lists or NumPy arrays) so that backends can populate
    and consume it without hard dependencies.

    Parameters
    ----------
    freq : array_like or None
        Frequency vector, length ``n``.
    z : array_like or None
        Impedance tensor, shape ``(n, 2, 2)`` (complex).
    z_err : array_like or None
        Errors for ``z``, shape ``(n, 2, 2)`` (float).
    tipper : array_like or None
        Magnetic tipper. Common shapes are ``(n, 2)`` or
        ``(n, 1, 2)``; callers may normalize as needed.
    tipper_err : array_like or None
        Errors for the tipper, shape like ``tipper``.
    rho : array_like or None
        Apparent resistivity, length ``n`` (float).
    phase : array_like or None
        Phase, length ``n``. Units are backend-specific; the
        Zonge path often stores milliradians.
    station : str or None
        Preferred station/site name.
    station_id : str, int or None
        Identifier used to synthesize a name when missing.
    lat, lon, elev : float or None
        Optional geographic location and elevation.
    azimuth : float or None
        Sensor azimuth or site orientation (degrees).
    attrs : dict
        Free-form attributes that accompany the data.

    Notes
    -----
    Only one of ``z`` or ``(rho, phase)`` is required for most
    workflows. Transformers may fill the missing part using
    package settings (see :mod:`pycsamt.core.config`).

    Examples
    --------
    >>> from pycsamt.core.base import TFBundle
    >>> TFBundle(freq=[1.0], z=[[ [0+1j, 0], [0, 0+1j] ]],
    ...          station=\"S001\")
    TFBundle(...)
    """

    # neutral payload used across backends
    freq: Any | None = None
    z: Any | None = None
    z_err: Any | None = None
    tipper: Any | None = None
    tipper_err: Any | None = None
    rho: Any | None = None
    phase: Any | None = None
    station: str | None = None
    station_id: str | int | None = None
    lat: float | None = None
    lon: float | None = None
    elev: float | None = None
    azimuth: float | None = None
    attrs: dict[str, Any] = field(default_factory=dict)

    def is_empty(self) -> bool:
        r"""
        Return True if there is no usable TF content.

        A bundle is considered empty when neither ``z`` nor the pair
        ``(rho, phase)`` is present.

        Returns
        -------
        bool
            ``True`` if empty, ``False`` otherwise.

        Examples
        --------
        >>> TFBundle().is_empty()
        True
        >>> TFBundle(rho=[100], phase=[45]).is_empty()
        False
        """
        have_z = self.z is not None
        have_rp = (self.rho is not None) and (self.phase is not None)
        return not (have_z or have_rp)


class SupportsToBundle(Protocol):
    r"""
    Protocol for objects that can export a :class:`TFBundle`.

    Any class implementing:

    ``to_bundle(self) -> TFBundle``

    is considered compatible. This keeps interop light while
    avoiding a hard dependency on specific backends.
    """

    def to_bundle(self) -> TFBundle: ...  # noqa: E701


class SupportsFromBundle(Protocol):
    r"""
    Protocol for classes that can be built from a
    :class:`TFBundle`.

    Any class implementing the classmethod:

    ``from_bundle(cls, bundle: TFBundle) -> Any``

    is considered compatible. This is commonly used by test
    doubles or thin wrappers in frontends.
    """

    @classmethod
    def from_bundle(cls, bundle: TFBundle): ...  # noqa: E701


def ensure_station(
    name: str | None,
    station_id: str | int | None,
    *,
    policy: StationNamePolicy | None = None,
) -> str:
    r"""
    Return a valid station name using policy rules.

    Validation is performed by the active
    :class:`~pycsamt.core.config.StationNamePolicy`. If the given
    ``name`` is invalid or missing, a synthetic one is derived
    from ``station_id``.

    Parameters
    ----------
    name : str or None
        Preferred station/site name.
    station_id : str, int or None
        Identifier used to synthesize a fallback name.
    policy : StationNamePolicy, optional
        Custom policy. Defaults to the global policy from
        :func:`~pycsamt.core.config.get_config`.

    Returns
    -------
    str
        Validated or synthesized station name.

    See Also
    --------
    pycsamt.core.config.StationNamePolicy
        Normalization, synthesis and limits.

    Examples
    --------
    >>> ensure_station(None, 7)
    'S007'
    """
    pol = policy or get_config().station_policy
    return pol.ensure(name, station_id)


def pick_adapter_key(
    obj: Any,
    *,
    hint: str | None = None,
) -> str | None:
    r"""
    Infer an adapter key (``'avg'``, ``'j'``, ``'edi'``) from an
    object.

    Heuristics look at the object's module and class name. If a
    ``hint`` is provided, it is returned verbatim (lower-cased).

    Parameters
    ----------
    obj : Any
        Source object whose kind should be inferred.
    hint : str, optional
        Explicit override for the key.

    Returns
    -------
    str or None
        Inferred key or ``None`` if it cannot be decided.

    Notes
    -----
    This function does not validate that an adapter exists for
    the returned key. Use :func:`to_edi` to dispatch safely.

    Examples
    --------
    >>> class X: pass
    >>> X.__module__ = 'pycsamt.zonge.avg'
    >>> pick_adapter_key(X())  # zonge path
    'avg'
    """

    if hint:
        return hint.lower()

    try:
        mod = obj.__class__.__module__.lower()
        cls = obj.__class__.__name__.lower()
    except Exception:
        return None

    if ("zonge" in mod) or ("avg" in cls):
        return "avg"
    if ("jones" in mod) or (cls in {"jfile", "jcollection"}):
        return "j"
    if ("seg" in mod) and ("edi" in (mod + cls)):
        return "edi"
    return None


def to_edi(
    source: Any,
    *,
    key: str | None = None,
    **kwargs: Any,
) -> Any:
    r"""
    Dispatch ``source`` to a registered adapter and return EDI.

    This function consults the adapter registry managed by
    :mod:`pycsamt.core.config`. It attempts to infer the adapter
    key when not provided and calls the associated factory.

    Parameters
    ----------
    source : Any
        The object to convert (AVG, Jones, etc.).
    key : str, optional
        Adapter key (e.g., ``'avg'``, ``'j'``, ``'edi'``). If
        omitted, :func:`pick_adapter_key` is used.
    **kwargs : Any
        Extra keyword arguments forwarded to the adapter.

    Returns
    -------
    Any
        An EDI object or an EDI collection, depending on the
        adapter.

    Raises
    ------
    RuntimeError
        If the key cannot be inferred or no adapter is registered
        for it.

    See Also
    --------
    pycsamt.core.config.register_adapter
        Register new adapters at runtime.
    pick_adapter_key
        Lightweight key inference from object metadata.

    Examples
    --------
    Register a trivial adapter and convert an object:

    >>> from pycsamt.core.config import register_adapter
    >>> class Dummy: pass
    >>> Dummy.__module__ = 'pycsamt.zonge.avg'
    >>> register_adapter('avg', lambda src, **k: {'edi': True})
    >>> to_edi(Dummy())
    {'edi': True}
    """

    k = key or pick_adapter_key(source)
    if not k:
        raise RuntimeError("Cannot infer adapter key")

    factory = _core_config.get_adapter(k)
    if not factory:
        raise RuntimeError(f"No adapter registered for: {k}")

    return factory(source, **kwargs)


class MTBase(CoreObject):
    r"""
    Common electromagnetic and MT utilities.

    This mixin-like base provides numerically safe, vectorized
    helpers built on NumPy. All methods accept array-like inputs
    and return NumPy arrays with broadcast semantics.

    Notes
    -----
    Shapes are preserved and operations are element-wise unless
    stated otherwise. For impedance tensors shaped ``(..., 2, 2)``,
    "determinant" operations collapse the last two axes.

    **Unit conventions.** Two equivalent formulas are common for
    apparent resistivity from impedance magnitude ``|Z|``:

    * **SI (E in V/m, H in A/m):**
      ``ρₐ = |Z|² / (μ₀·ω) = |Z|² · RHO_FACTOR / f``, with
      ``ω = 2π f`` and ``RHO_FACTOR = 1/(μ₀·2π)``.
    * **Field units (E in mV/km, B in nT):**
      ``ρₐ ≈ (0.2 / f) · |E/H|²`` — a legacy Zonge-style form.
      This is reflected by :data:`ZONGE_RHO_FACTOR`, i.e. ``0.2``.

    Be consistent. Converting mV/km → V/m and nT → T, and B → H
    via ``H = B/μ₀``, recovers the SI formulation.

    Constants
    ---------
    MU0 : float
        Magnetic permeability of free space (H/m).
        ``4π × 10⁻⁷``.
    EPS0 : float
        Electric permittivity of free space (F/m).
        ``8.854187817 × 10⁻¹²``.
    C : float
        Speed of light in vacuum (m/s). From CODATA.
    C0 : float
        ``1/√(μ₀·ε₀)`` (m/s). Equal to :data:`C` by definition.
    ETA0 : float
        Wave impedance of free space (Ω). ``√(μ₀/ε₀)``.
    TWO_PI : float
        ``2π``. Handy for ``ω = 2π f`` conversions.
    DEG2RAD : float
        Degrees to radians. ``π/180``.
    RAD2DEG : float
        Radians to degrees. ``180/π``.
    RHO_FACTOR : float
        ``1/(μ₀·2π)``. For ``ρₐ = |Z|² · RHO_FACTOR / f`` (SI).
    ZONGE_RHO_FACTOR : float
        ``0.2``. Legacy factor in ``ρₐ ≈ (0.2/f)·|E/H|²`` when
        using E in mV/km and B/H in nT without explicit SI
        conversion in the formula.
    B_TO_H : float
        Conversion factor from magnetic flux density to field
        intensity (A/m per Tesla). ``1/μ₀``.
    # Unit scales (field practice)
    MICROVOLTS_TO_VOLTS : float
        ``1e-6``.
    PICOTESLA_TO_TESLA : float
        ``1e-12``.
    NANOTESLA_TO_TESLA : float
        ``1e-9``.
    MV_PER_KM_TO_V_PER_M : float
        ``1e-6``. Since ``1 mV/km = 1e-3 V / 1e3 m``.
    METERS_TO_KILOMETERS : float
        ``1e-3``.
    PERCENT_FACTOR : float
        ``100.0``. For percent amplitudes or coherence.
    Z_UNIT_MVK_NT_TO_SI : float
        ``(1e-6) / (1e-9) = 1e3``. Converts ``Z`` expressed as
        ``(mV/km)/nT`` into SI-consistent ``(V/m)/T``.

    See Also
    --------
    rho_phase_from_z : Apparent resistivity and phase from ``Z``.
    z_from_rho_phase : Build ``Z`` magnitude/phase from ``ρ, φ``.
    determinant_z : Rotation-invariant impedance.
    rho_phase_from_det : Invariants from determinant impedance.
    rotate_impedance : Rotate 2×2 impedance tensors.
    skin_depth : Diffusive skin depth.
    tipper_amp_phase : Amplitude and phase of tipper vector.

    Examples
    --------
    Compute SI apparent resistivity from ``|Z|`` and frequency::

        >>> from pycsamt.core.base import MTBase
        >>> import numpy as np
        >>> Zmag = np.asarray([5.0, 10.0])        # |Z| in ohms
        >>> f = np.asarray([1.0, 0.1])            # Hz
        >>> rho_a = (Zmag**2) * MTBase.RHO_FACTOR / f
        >>> rho_a.shape
        (2,)

    Convert field units to SI and compare formulas::

        >>> E_mVkm = 100.0      # mV/km
        >>> B_nT = 50.0         # nT (≈ μ0·H if not converted)
        >>> Z_field = (E_mVkm * MTBase.MV_PER_KM_TO_V_PER_M) / \
        ...           (B_nT * MTBase.NANOTESLA_TO_TESLA)
        >>> Z_si_like = Z_field  # (V/m)/T
        >>> f = 1.0  # Hz
        >>> rho_si = (Z_si_like**2) * MTBase.RHO_FACTOR / f
        >>> rho_zonge = (MTBase.ZONGE_RHO_FACTOR / f) * (E_mVkm / B_nT)**2
        >>> rho_si > 0 and rho_zonge > 0
        True

    Rotate an impedance tensor by a site azimuth::

        >>> Z = np.zeros((3, 2, 2), dtype=complex)  # demo only
        >>> phi = 30.0 * MTBase.DEG2RAD
        >>> # rotate_impedance would apply R(φ) Z Rᵀ(φ)
        >>> # Z_rot = MTBase.rotate_impedance(Z, phi)  # hypothetical

    References
    ----------
    .. [1] Chave & Jones (2012). *The Magnetotelluric Method*. CUP.
    .. [2] Simpson & Bahr (2005). *Practical MT*. Cambridge Univ.
    .. [3] Wait (1982). *Geo-Electromagnetism*. Academic Press.
    .. [4] Constable (2016). *Scaling relations in EM*. SEG Notes.
    """

    # --- Fundamental EM constants (CODATA) -----------------------------
    MU0: float = 4.0e-7 * float(np.pi)  # H/m  (μ0)
    EPS0: float = 8.854_187_817e-12  # F/m  (ε0)
    C: float = 299_792_458.0  # m/s  (speed of light)

    # Derived free-space constants
    C0: float = 1.0 / np.sqrt(MU0 * EPS0)  # m/s  (= C by definition)
    ETA0: float = np.sqrt(MU0 / EPS0)  # Ω    (wave impedance ~376.73 Ω)

    # Handy math constants
    TWO_PI: float = 2.0 * float(np.pi)  # 2π
    DEG2RAD: float = float(np.pi) / 180.0  # deg → rad
    RAD2DEG: float = 180.0 / float(np.pi)  # rad → deg

    # Apparent resistivity factors
    # ρa = |Z|² / (μ0·ω) = |Z|² * [1/(μ0·2π)] / f   (SI form)
    RHO_FACTOR: float = 1.0 / (MU0 * 2.0 * float(np.pi))  # SI factor

    # Zonge field-units convention (E in mV/km, B or H in nT):
    # ρa ≈ (0.2 / f) * |E/H|²  → 0.2 = 1/5
    ZONGE_RHO_FACTOR: float = 0.2  # legacy factor

    # --- Unit conversions ------------------------------------------------
    # Voltage
    MICROVOLTS_TO_VOLTS: float = 1e-6  # μV → V
    MILLIVOLTS_TO_VOLTS: float = 1e-3  # mV → V

    # Electric field
    MV_PER_KM_TO_V_PER_M: float = 1e-6  # (mV/km) → (V/m)

    # Magnetic field / flux density
    PICOTESLA_TO_TESLA: float = 1e-12  # pT → T
    NANOTESLA_TO_TESLA: float = 1e-9  # nT → T
    MICROTESLA_TO_TESLA: float = 1e-6  # μT → T
    H_TO_B: float = MU0  # H (A/m) → B (T): B = μ0·H
    B_TO_H: float = 1.0 / MU0  # B (T) → H (A/m): H = B/μ0

    # Length
    METERS_TO_KILOMETERS: float = 1e-3  # m → km
    KILOMETERS_TO_METERS: float = 1e3  # km → m

    # Percent scaling
    PERCENT_FACTOR: float = 100.0  # ×100 for %

    # Mixed-unit convenience:
    # Convert Z in (mV/km)/nT → (V/m)/T (SI-consistent)
    # factor = (mV/km→V/m) / (nT→T) = 1e-6 / 1e-9 = 1e3
    Z_UNIT_MVK_NT_TO_SI: float = 1e3  # (mV/km)/nT → (V/m)/T

    @staticmethod
    def _as_c(x: object) -> np.ndarray:
        return np.asarray(x, dtype=np.complex128)

    @staticmethod
    def _as_f(x: object) -> np.ndarray:
        return np.asarray(x, dtype=np.float64)

    @staticmethod
    def _ang(
        z: np.ndarray,
        *,
        unit: str = "deg",
    ) -> np.ndarray:
        phi = np.angle(z)
        if unit == "deg":
            return np.degrees(phi)
        if unit == "mrad":
            return phi * 1.0e3
        return phi  # "rad"

    @staticmethod
    def _to_rad(
        phi: object,
        *,
        unit: str,
    ) -> np.ndarray:
        p = np.asarray(phi, dtype=np.float64)
        if unit == "deg":
            return np.radians(p)
        if unit == "mrad":
            return p * 1.0e-3
        return p  # "rad"

    @staticmethod
    def omega(f: object) -> np.ndarray:
        r"""
        Angular frequency ``ω = 2π f``.
        """
        ff = np.asarray(f, dtype=np.float64)
        return 2.0 * float(np.pi) * ff

    #  Z ↔ ρa, φ

    def rho_phase_from_z(
        self,
        z: object,
        f: object,
        *,
        phase_unit: str = "deg",
    ) -> tuple[np.ndarray, np.ndarray]:
        r"""
        Apparent resistivity and phase from impedance.

        Parameters
        ----------
        z : array_like
            Impedance, scalar or tensor ``(..., 2, 2)``.
        f : array_like
            Frequency (Hz), broadcastable with ``z``.
        phase_unit : {'deg', 'rad', 'mrad'}, optional
            Output unit for phase. Default is degrees.

        Returns
        -------
        rho : ndarray
            Apparent resistivity (Ω·m), element-wise.
        phi : ndarray
            Phase in requested unit, element-wise.

        Notes
        -----
        Uses ``ρa = |Z|^2 / (μ0 ω)`` per component. For tensor
        invariants see :meth:`rho_phase_from_det`.
        """
        zz = self._as_c(z)
        w = self.omega(f)
        rho = (np.abs(zz) ** 2) / (self.MU0 * w)
        phi = self._ang(zz, unit=phase_unit)
        return rho, phi

    def z_from_rho_phase(
        self,
        rho: object,
        phi: object,
        f: object,
        *,
        phase_unit: str = "deg",
    ) -> np.ndarray:
        r"""
        Build impedance magnitude/phase from ``ρa`` and ``φ``.

        Parameters
        ----------
        rho : array_like
            Apparent resistivity (Ω·m).
        phi : array_like
            Phase in ``phase_unit``.
        f : array_like
            Frequency (Hz).
        phase_unit : {'deg', 'rad', 'mrad'}, optional
            Unit of ``phi``. Default is degrees.

        Returns
        -------
        Z : ndarray
            Complex impedance with ``|Z| = sqrt(μ0 ω ρa)`` and
            angle = ``φ``. Shape follows broadcasting rules.
        """
        w = self.omega(f)
        amp = np.sqrt(self.MU0 * w * self._as_f(rho))
        ang = self._to_rad(phi, unit=phase_unit)
        return amp * (np.cos(ang) + 1j * np.sin(ang))

    # determinant invariant

    @staticmethod
    def determinant_z(z: object) -> np.ndarray:
        r"""
        Determinant impedance ``Z_det = sqrt(-Z_xy Z_yx)``.

        Parameters
        ----------
        z : array_like, shape (..., 2, 2)
            Impedance tensor(s).

        Returns
        -------
        Z_det : ndarray, shape (...)
            Rotation-invariant complex impedance.

        Notes
        -----
        The minus sign assumes the common MT sign convention.
        """
        zz = np.asarray(z, dtype=complex)
        zxy = zz[..., 0, 1]
        zyx = zz[..., 1, 0]
        return np.sqrt(-zxy * zyx)

    def rho_phase_from_det(
        self,
        z: object,
        f: object,
        *,
        phase_unit: str = "deg",
    ) -> tuple[np.ndarray, np.ndarray]:
        r"""
        Invariant apparent resistivity and phase from ``Z``.

        Parameters
        ----------
        z : array_like, shape (..., 2, 2)
            Impedance tensor(s).
        f : array_like
            Frequency (Hz).
        phase_unit : {'deg', 'rad', 'mrad'}, optional
            Output unit for phase.

        Returns
        -------
        rho_det : ndarray
            Determinant apparent resistivity (Ω·m).
        phi_det : ndarray
            Determinant phase in requested unit.
        """
        zd = self.determinant_z(z)
        return self.rho_phase_from_z(zd, f, phase_unit=phase_unit)

    # rotation

    @staticmethod
    def rotate_impedance(
        z: object,
        theta_deg: float,
    ) -> np.ndarray:
        r"""
        Rotate 2×2 impedance tensor(s) by ``theta_deg``.

        Parameters
        ----------
        z : array_like, shape (..., 2, 2)
            Impedance tensor(s).
        theta_deg : float
            Rotation angle in degrees. Positive rotates
            the coordinate frame clockwise (convention).

        Returns
        -------
        Z_rot : ndarray, shape (..., 2, 2)
            Rotated impedance tensors.

        Notes
        -----
        Uses ``Z' = R Z Rᵀ`` with

        ``R = [[cos θ,  sin θ], [-sin θ, cos θ]]``.
        """
        zz = np.asarray(z, dtype=complex)
        th = np.deg2rad(theta_deg)
        c = np.cos(th)
        s = np.sin(th)
        r = np.array([[c, s], [-s, c]])
        rt = r.T
        return r @ zz @ rt

    # diffusion / scales

    def skin_depth(
        self,
        f: object,
        rho: object,
    ) -> np.ndarray:
        r"""
        Skin depth ``δ = sqrt(2ρ / (μ0 ω))`` in meters.

        Parameters
        ----------
        f : array_like
            Frequency (Hz).
        rho : array_like
            Resistivity (Ω·m).

        Returns
        -------
        δ : ndarray
            Skin depth in meters.

        Notes
        -----
        For quick estimates, ``δ ≈ 503 √(ρ/f)`` (ρ in Ω·m,
        f in Hz).
        """
        w = self.omega(f)
        return np.sqrt(2.0 * self._as_f(rho) / (self.MU0 * w))

    # time-frequency

    @staticmethod
    def freq_to_period(f: object) -> np.ndarray:
        r"""
        Convert frequency (Hz) to period (s), safe at zero.
        """
        ff = np.asarray(f, dtype=float)
        out = np.empty_like(ff)
        with np.errstate(divide="ignore", invalid="ignore"):
            out = np.where(ff == 0.0, np.inf, 1.0 / ff)
        return out

    @staticmethod
    def period_to_freq(t: object) -> np.ndarray:
        r"""
        Convert period (s) to frequency (Hz), safe at inf.
        """
        tt = np.asarray(t, dtype=float)
        out = np.empty_like(tt)
        with np.errstate(divide="ignore", invalid="ignore"):
            out = np.where(np.isinf(tt), 0.0, 1.0 / tt)
        return out

    # tipper

    @staticmethod
    def tipper_amp_phase(
        t: object,
        *,
        phase_unit: str = "deg",
    ) -> tuple[np.ndarray, np.ndarray]:
        r"""
        Amplitude and phase of tipper vectors.

        Parameters
        ----------
        t : array_like, shape (..., 2) or (..., 1, 2)
            Tipper components ``(Tx, Ty)``.
        phase_unit : {'deg', 'rad', 'mrad'}, optional
            Output unit for phase. Default is degrees.

        Returns
        -------
        amp : ndarray
            Vector magnitude ``sqrt(Tx² + Ty²)``.
        phi : ndarray
            Phase of the complex vector ``Tx + i Ty``.

        Notes
        -----
        If shape is ``(..., 1, 2)`` it is squeezed to
        ``(..., 2)`` prior to computation.
        """
        tt = np.asarray(t, dtype=complex)
        if tt.ndim >= 2 and tt.shape[-2:] == (1, 2):
            tt = tt[..., 0, :]
        tx = tt[..., 0]
        ty = tt[..., 1]
        vec = tx + 1j * ty
        amp = np.abs(vec)
        phi = np.angle(vec)
        if phase_unit == "deg":
            phi = np.degrees(phi)
        elif phase_unit == "mrad":
            phi = phi * 1.0e3
        return amp, phi

    def phase_tensor(self, z: object) -> np.ndarray:
        r"""
        Phase tensor ``Φ = X⁺ Y`` from impedance ``Z = X + iY``.

        This computes a real 2×2 phase tensor using a robust
        pseudoinverse ``X⁺`` of the real part of ``Z``. Shapes are
        broadcast over the leading axes.

        Parameters
        ----------
        z : array_like, shape (..., 2, 2)
            Impedance tensor(s) ``Z = X + iY``.

        Returns
        -------
        Phi : ndarray, shape (..., 2, 2)
            Real phase tensor.

        Notes
        -----
        If ``X`` is ill-conditioned, the pseudoinverse stabilizes the
        solution. ``Φ`` is invariant to site gain and captures the
        2-D/3-D character of the response.

        Examples
        --------
        >>> from pycsamt.core.base import MTBase
        >>> import numpy as np
        >>> Z = np.array([[[0+1j, 2+3j],[4+5j, 6+7j]]], dtype=complex)
        >>> Phi = MTBase().phase_tensor(Z)
        >>> Phi.shape
        (1, 2, 2)

        See Also
        --------
        phase_tensor_params
        phase_tensor_azimuth
        """

        zz = np.asarray(z, dtype=complex)
        X = zz.real
        Y = zz.imag
        # Use pseudoinverse for robustness
        X_pinv = np.linalg.pinv(X)
        return np.einsum("...ij,...jk->...ik", X_pinv, Y)

    def phase_tensor_params(
        self,
        z: object,
        *,
        angle_unit: str = "deg",
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        r"""
        Principal phase angles, azimuth, skew, and ellipticity.

        Derives quick-look parameters from the phase tensor ``Φ``:
        principal phase angles (``φ_max``, ``φ_min``), phase tensor
        azimuth ``α``, skew angle ``β``, and ellipticity.

        Parameters
        ----------
        z : array_like, shape (..., 2, 2)
            Impedance tensor(s).
        angle_unit : {'deg', 'rad', 'mrad'}, optional
            Output unit for angular quantities. Default is degrees.

        Returns
        -------
        phi_max : ndarray, shape (...)
            Larger principal phase angle derived from ``Φ``.
        phi_min : ndarray, shape (...)
            Smaller principal phase angle derived from ``Φ``.
        alpha : ndarray, shape (...)
            Phase tensor azimuth.
        beta : ndarray, shape (...)
            Phase tensor skew angle.
        ellipt : ndarray, shape (...)
            Ellipticity proxy ``(λ_max - λ_min)/(λ_max + λ_min)``,
            where ``λ`` are the eigenvalues used internally.

        Notes
        -----
        This implementation obtains parameters from an eigen
        decomposition of ``Φ`` and simple angle formulas (Caldwell-
        style azimuth/skew heuristics). It is suitable for quality
        control and strike reconnaissance.

        Examples
        --------
        >>> from pycsamt.core.base import MTBase
        >>> import numpy as np
        >>> Z = np.zeros((5, 2, 2), complex)
        >>> z = MTBase()
        >>> pmax, pmin, az, skew, eps = z.phase_tensor_params(Z)
        >>> pmax.shape == pmin.shape == az.shape == skew.shape == eps.shape
        True

        See Also
        --------
        phase_tensor
        phase_tensor_azimuth
        """

        Phi = self.phase_tensor(z).astype(float)
        a, b = Phi[..., 0, 0], Phi[..., 0, 1]
        c, d = Phi[..., 1, 0], Phi[..., 1, 1]
        # Caldwell et al. 2004 parameterization
        alpha = 0.5 * np.arctan2(b - c, a + d)  # azimuth
        beta = 0.5 * np.arctan2(b + c, a - d)  # skew
        # Eigenvalues (real for real Phi)
        # vals = np.linalg.eigvals(Phi)
        # vals = np.real_if_close(vals, tol=1e5)
        # lam_max = np.maximum(vals[..., 0], vals[..., 1])
        # lam_min = np.minimum(vals[..., 0], vals[..., 1])
        # phi_max = np.arctan(lam_max)
        # phi_min = np.arctan(lam_min)

        # existing azimuth/skew (alpha, beta) formulas may stay as-is
        # but for principal phases use singular values:
        s = np.linalg.svd(Phi, compute_uv=False)  # shape (..., 2)
        lam_max = np.maximum(s[..., 0], s[..., 1])
        lam_min = np.minimum(s[..., 0], s[..., 1])
        phi_max = np.arctan(lam_max)
        phi_min = np.arctan(lam_min)

        ellipt = (lam_max - lam_min) / (lam_max + lam_min + 1e-18)
        if angle_unit == "deg":
            phi_max = np.degrees(phi_max)
            phi_min = np.degrees(phi_min)
            alpha = np.degrees(alpha)
            beta = np.degrees(beta)
        elif angle_unit == "mrad":
            k = 1e3
            phi_max = phi_max * k
            phi_min = phi_min * k
            alpha = alpha * k
            beta = beta * k
        return phi_max, phi_min, alpha, beta, ellipt

    def phase_tensor_azimuth(
        self,
        z: object,
        *,
        unit: str = "deg",
    ) -> np.ndarray:
        r"""
        Azimuth of the phase tensor principal axis.

        Parameters
        ----------
        z : array_like, shape (..., 2, 2)
            Impedance tensor(s).
        unit : {'deg', 'rad', 'mrad'}, optional
            Output unit. Default is degrees.

        Returns
        -------
        alpha : ndarray, shape (...)
            Phase tensor azimuth.

        Examples
        --------
        >>> from pycsamt.core.base import MTBase
        >>> import numpy as np
        >>> Z = np.zeros((3, 2, 2), complex)
        >>> MTBase().phase_tensor_azimuth(Z).shape
        (3,)

        See Also
        --------
        phase_tensor
        phase_tensor_params
        """
        _, _, az, _, _ = self.phase_tensor_params(z, angle_unit=unit)
        return az

    @staticmethod
    def tipper_rotate(t: object, theta_deg: float) -> np.ndarray:
        r"""
        Rotate tipper vectors by ``theta_deg`` in the horizontal plane.

        Parameters
        ----------
        t : array_like, shape (..., 2) or (..., 1, 2)
            Tipper components ``(T_x, T_y)``. A trailing ``(1, 2)``
            is squeezed.
        theta_deg : float
            Rotation angle in degrees. Positive angles rotate the
            coordinates clockwise (same convention as
            :meth:`rotate_impedance`).

        Returns
        -------
        t_rot : ndarray, shape (..., 2)
            Rotated tipper vectors.

        Examples
        --------
        >>> from pycsamt.core.base import MTBase
        >>> import numpy as np
        >>> T = np.array([[0.2+0.1j, -0.1+0.3j]])
        >>> MTBase.tipper_rotate(T, 30.0).shape
        (1, 2)

        See Also
        --------
        rotate_impedance
        induction_arrows
        """

        tt = np.asarray(t, dtype=complex)
        if tt.ndim >= 2 and tt.shape[-2:] == (1, 2):
            tt = tt[..., 0, :]
        th = np.deg2rad(theta_deg)
        c, s = np.cos(th), np.sin(th)
        R = np.array([[c, s], [-s, c]])
        return np.einsum("ij,...j->...i", R, tt)

    @staticmethod
    def induction_arrows(
        t: object,
        *,
        convention: str = "wiese",
        use_imag: bool = False,
    ) -> tuple[np.ndarray, np.ndarray]:
        r"""
        Compute 2-D induction arrows from tipper components.

        Parameters
        ----------
        t : array_like, shape (..., 2) or (..., 1, 2)
            Tipper components ``(T_x, T_y)``.
        convention : {'wiese', 'parkinson'}, optional
            Arrow convention. ``'wiese'`` uses a 90° rotation from
            the real (or imaginary) tipper; ``'parkinson'`` uses the
            components directly. Default ``'wiese'``.
        use_imag : bool, optional
            If ``True``, use imaginary parts; else use real parts.
            Default ``False``.

        Returns
        -------
        ax : ndarray, shape (...)
            X-component of the arrow.
        ay : ndarray, shape (...)
            Y-component of the arrow.

        Notes
        -----
        Wiese arrows often highlight current channeling (2-D/3-D).
        Use imaginary arrows (``use_imag=True``) for deeper structure
        emphasis at higher periods.

        Examples
        --------
        >>> from pycsamt.core.base import MTBase
        >>> import numpy as np
        >>> T = np.array([[0.2+0.1j, -0.1+0.3j]])
        >>> ax, ay = MTBase.induction_arrows(T, convention='wiese')
        >>> ax.shape == ay.shape
        True

        See Also
        --------
        tipper_rotate
        """

        tt = np.asarray(t, dtype=complex)
        if tt.ndim >= 2 and tt.shape[-2:] == (1, 2):
            tt = tt[..., 0, :]
        tx = tt[..., 0]
        ty = tt[..., 1]
        rx = tx.imag if use_imag else tx.real
        ry = ty.imag if use_imag else ty.real
        if convention.lower() == "wiese":
            ax, ay = -ry, rx
        else:  # "parkinson" or fallback
            ax, ay = rx, ry
        return ax, ay

    def swift_skew(
        self,
        z: object,
        *,
        unit: str = "deg",
    ) -> tuple[np.ndarray, np.ndarray]:
        r"""
        Swift skew parameter from the impedance tensor.

        Computes the complex Swift skew ratio

        ``s = (Z_xx + Z_yy) / (Z_xy - Z_yx)``,

        returning its magnitude and angle.

        Parameters
        ----------
        z : array_like, shape (..., 2, 2)
            Impedance tensor(s).
        unit : {'deg', 'rad', 'mrad'}, optional
            Unit for the returned angle. Default is degrees.

        Returns
        -------
        amp : ndarray, shape (...)
            Magnitude ``|s|`` (dimensionless).
        ang : ndarray, shape (...)
            Argument of ``s`` in the requested unit.

        Notes
        -----
        Low ``|s|`` values are consistent with 1-D/2-D structure; an
        elevated skew suggests 3-D effects or distortion.

        Examples
        --------
        >>> from pycsamt.core.base import MTBase
        >>> import numpy as np
        >>> Z = np.ones((4, 2, 2), dtype=complex)
        >>> amp, ang = MTBase().swift_skew(Z)
        >>> amp.shape == ang.shape
        True
        """

        zz = np.asarray(z, dtype=complex)
        zxx, zxy = zz[..., 0, 0], zz[..., 0, 1]
        zyx, zyy = zz[..., 1, 0], zz[..., 1, 1]
        s = (zxx + zyy) / (zxy - zyx + 1e-30)
        amp = np.abs(s)
        ang = np.angle(s)
        if unit == "deg":
            ang = np.degrees(ang)
        elif unit == "mrad":
            ang = ang * 1e3
        return amp, ang

    def apparent_conductivity_from_z(
        self,
        z: object,
        f: object,
    ) -> np.ndarray:
        r"""
        Apparent conductivity ``σₐ = μ₀ ω / |Z|²`` (S/m).

        Parameters
        ----------
        z : array_like
            Impedance (per component or invariant).
        f : array_like
            Frequency (Hz), broadcastable to ``z``.

        Returns
        -------
        sigma_a : ndarray
            Apparent conductivity (S/m).

        Notes
        -----
        This is the reciprocal of the common SI resistivity formula
        ``ρₐ = |Z|² / (μ₀ ω)``.

        Examples
        --------
        >>> from pycsamt.core.base import MTBase
        >>> import numpy as np
        >>> Z = np.array([5.0, 10.0])       # Ω
        >>> f = np.array([1.0, 0.1])        # Hz
        >>> MTBase().apparent_conductivity_from_z(Z, f).shape
        (2,)

        See Also
        --------
        rho_phase_from_z
        rho_phase_from_det
        """

        zz = np.asarray(z, dtype=complex)
        w = self.omega(f)
        return (self.MU0 * w) / (np.abs(zz) ** 2 + 1e-30)

    def halfspace_impedance(
        self,
        f: object,
        rho: object,
    ) -> np.ndarray:
        r"""
        Plane-wave impedance over a uniform half-space.

        For a homogeneous half-space of resistivity ``ρ``, the
        surface impedance is

        ``Z = (1 + i) / √2 · √(μ₀ ω ρ)``,

        with a constant 45° phase and amplitude ``∝ √(ω ρ)``.

        Parameters
        ----------
        f : array_like
            Frequency (Hz).
        rho : array_like
            Resistivity (Ω·m).

        Returns
        -------
        Z : ndarray
            Complex half-space impedance.

        Examples
        --------
        >>> from pycsamt.core.base import MTBase
        >>> import numpy as np
        >>> f = np.logspace(-3, 3, 5)
        >>> Z = MTBase().halfspace_impedance(f, 100.0)
        >>> Z.shape
        (5,)
        """

        w = self.omega(f)
        amp = np.sqrt(self.MU0 * w * self._as_f(rho))
        return (1.0 + 1.0j) * (amp / np.sqrt(2.0))

    def z_mvk_nt_to_ohms(self, z: object) -> np.ndarray:
        r"""
        Convert ``Z`` from ``(mV/km)/nT`` to ohms (Ω).

        This helper assumes the denominator is magnetic flux density
        ``B`` in nT (not ``H``). It converts to SI ``(V/m)/T`` and
        then multiplies by ``μ₀`` to obtain ``V/A = Ω``.

        Parameters
        ----------
        z : array_like
            Impedance in mixed field units ``(mV/km)/nT``.

        Returns
        -------
        Z_ohms : ndarray
            Impedance in ohms (Ω).

        Notes
        -----
        If your denominator is already the magnetic field intensity
        ``H`` in A/m, you *should not* multiply by ``μ₀``. Use
        careful unit bookkeeping in mixed-unit workflows.

        Examples
        --------
        >>> from pycsamt.core.base import MTBase
        >>> import numpy as np
        >>> z_field = np.array([50.0])  # (mV/km)/nT
        >>> Z_ohm = MTBase().z_mvk_nt_to_ohms(z_field)
        >>> Z_ohm.shape
        (1,)
        """

        zz = np.asarray(z, dtype=complex)
        return zz * self.Z_UNIT_MVK_NT_TO_SI * self.MU0

    @staticmethod
    def rotate_fields(
        e: object,
        h: object,
        theta_deg: float,
    ) -> tuple[np.ndarray, np.ndarray]:
        r"""
        Rotate horizontal electric and magnetic field vectors.

        Parameters
        ----------
        e : array_like, shape (..., 2)
            Horizontal electric field components ``(E_x, E_y)``.
        h : array_like, shape (..., 2)
            Horizontal magnetic field components ``(H_x, H_y)``.
        theta_deg : float
            Rotation angle in degrees. Positive rotates the
            coordinate frame clockwise.

        Returns
        -------
        e_rot : ndarray, shape (..., 2)
            Rotated electric field components.
        h_rot : ndarray, shape (..., 2)
            Rotated magnetic field components.

        Examples
        --------
        >>> from pycsamt.core.base import MTBase
        >>> import numpy as np
        >>> e = np.array([[1.0, 0.0]])
        >>> h = np.array([[0.0, 1.0]])
        >>> e_r, h_r = MTBase.rotate_fields(e, h, 45.0)
        >>> e_r.shape == h_r.shape
        True

        See Also
        --------
        rotate_impedance
        tipper_rotate
        """

        ee = np.asarray(e)
        hh = np.asarray(h)
        th = np.deg2rad(theta_deg)
        c, s = np.cos(th), np.sin(th)
        R = np.array([[c, s], [-s, c]])
        e_rot = np.einsum("ij,...j->...i", R, ee)
        h_rot = np.einsum("ij,...j->...i", R, hh)
        return e_rot, h_rot
