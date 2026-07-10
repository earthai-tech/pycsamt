# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

from __future__ import annotations

from typing import Any

import numpy as np

from ..core.base import MTBase, TFBundle, ensure_station
from ..core.config import get_config
from ..seg.heads import Head, Info
from ..seg.meas import DefineMeas
from ..seg.mtemap import MTEMAP

__all__ = [
    "TransformerMixin",
]


class TransformerMixin(MTBase):
    r"""
    Template utilities for AVG/J -> EDI transformations.

    This module defines :class:`TransformerMixin`, a light
    template that standardizes the steps needed to convert
    foreign transfer-function containers (e.g. Zonge AVG or
    Jones J) into EDI files or collections.

    The mixin focuses on *finalization* of a neutral payload
    (:class:`pycsamt.core.base.TFBundle`): validating the
    station name, ordering and de-duplicating frequencies, and
    optionally filling missing pieces (``Z`` vs ``rho/phase``)
    according to user config.

    Notes
    -----
    Subclasses provide the data plumbing:

    * :meth:`TransformerMixin.extract` pulls a :class:`TFBundle`
      out of the source object or file.
    * :meth:`TransformerMixin.emit_edi` builds an EDI object
      from the finalized bundle.
    * :meth:`TransformerMixin.post_emit` can attach metadata or
      perform small last-mile fixes.

    The default hooks for *filling* missing parts are no-ops.
    Backends can override
    :meth:`TransformerMixin.compute_res_from_z` and
    :meth:`TransformerMixin.compute_z_from_res`.

    Subclasses implement the *data path* while this mixin
    handles standard *finalization* steps on a
    :class:`~pycsamt.core.base.TFBundle`:

    1. Validate or synthesize the station name using the
       global policy.
    2. Order frequencies per config (asc/desc).
    3. De-duplicate nearly equal frequencies using a relative
       tolerance.
    4. Optionally fill missing parts (``Z`` vs ``rho/phase``)
       via overridable hooks.

    Finalization obeys the global configuration from
    :mod:`pycsamt.core.config`. In particular:

    * ``freq_order`` controls sorting.
    * ``freq_tol`` controls de-duplication.
    * ``compute_res_from_z`` and ``compute_z_from_res`` toggle
      filling behavior.


    See Also
    --------
    pycsamt.core.transformers.AVGtoEDI
        Concrete AVG → EDI/EDICollection converter.
    pycsamt.core.transformers.JtoEDI
        Concrete J → EDI/EDICollection converter.
    pycsamt.core.base.TFBundle
        Neutral payload used across backends.
    pycsamt.core.config
        Global configuration and frequency policies.

    Examples
    --------
    A minimal custom transformer:

    >>> from pycsamt.core._transformers import TransformerMixin
    >>> from pycsamt.core.base import TFBundle
    >>>
    >>> class MyX(TransformerMixin):
    ...     def extract(self, source):
    ...         # pull TF data from source into a bundle
    ...         return TFBundle(freq=[1.0], rho=[100.], phase=[45.])
    ...     def emit_edi(self, bundle):
    ...         # return an EDI-like stub for illustration
    ...         return {'freq': bundle.freq, 'rho': bundle.rho}
    ...
    >>> out = MyX().transform(object())
    >>> isinstance(out, dict)
    True

    References
    ----------
    .. [1] Simpson, F. & Bahr, K. (2005). *Practical MT*.
    .. [2] Egbert, G. D. (1997). Robust MT processing.

    """


    def extract(self, source: Any) -> TFBundle:  # noqa: D401
        r"""
        Extract a :class:`TFBundle` from ``source``.

        This is the *ingest* step. Implementors parse the foreign
        object (or file path) and return a neutral bundle. Keep the
        method free of side effects.

        Parameters
        ----------
        source : Any
            Backend-specific object or path.

        Returns
        -------
        TFBundle
            A bundle containing frequencies and transfer functions.

        Raises
        ------
        NotImplementedError
            Must be provided by subclasses.

        Notes
        -----
        Prefer leaving station naming and frequency manipulation to
        the mixin. Provide raw content here.
        """

        raise NotImplementedError

    def emit_edi(self, bundle: TFBundle) -> Any:  # noqa: D401
        r"""
        Materialize an EDI object from a finalized bundle.

        This is the *emit* step. Implementors build and return an
        EDI object (or a compatible stub).

        Must be provided by subclasses.

        Notes
        -----
        Do not attempt to reorder frequencies here. The bundle has
        already been finalized by :meth:`_finalize`.
        """

        raise NotImplementedError


    def post_emit(
        self,
        edi_obj: Any,
        source: Any,
        bundle: TFBundle,
    ) -> Any:
        r"""
        Optional last-mile adjustments after EDI creation.

        Use this hook to enrich the EDI object with auxiliary
        metadata (e.g. site location, elevation) that is awkward
        to thread through earlier steps.

        Parameters
        ----------
        edi_obj : Any
            The EDI object returned by :meth:`emit_edi`.
        source : Any
            The original source object used in :meth:`extract`.
        bundle : TFBundle
            The finalized bundle.

        Returns
        -------
        Any
            The potentially modified EDI object.

        Notes
        -----
        The default implementation returns ``edi_obj`` unmodified.
        """

        return edi_obj

    def compute_res_from_z(self, b: TFBundle) -> TFBundle:
        r"""
        Compute apparent resistivity/phase from ``Z``.

        This hook is called by :meth:`_fill_missing` when
        ``compute_res_from_z`` is True in the global config, and
        the bundle holds ``Z`` but lacks ``(rho, phase)``.

        Parameters
        ----------
        b : TFBundle
            Input bundle.

        Returns
        -------
        TFBundle
            The updated bundle (may be the same instance).

        Notes
        -----
        The default implementation is a no-op. Subclasses may
        override to apply the standard MT relations [1]_ [2]_.
        """
        return b

    def compute_z_from_res(self, b: TFBundle) -> TFBundle:
        r"""
        Reconstruct ``Z`` from apparent resistivity/phase.

        This hook is called by :meth:`_fill_missing` when
        ``compute_z_from_res`` is True in the global config, and
        the bundle lacks ``Z`` but holds ``(rho, phase)``.

        Parameters
        ----------
        b : TFBundle
            Input bundle.

        Returns
        -------
        TFBundle
            The updated bundle (may be the same instance).

        Notes
        -----
        The default implementation is a no-op. Subclasses may
        override to apply physically consistent reconstruction.
        """

        return b


    def _order_freq(self, b: TFBundle) -> TFBundle:
        r"""
        Order frequencies and align arrays accordingly.

        The order is defined by ``freq_order`` in the global
        config (``'asc'`` or ``'desc'``). All frequency-indexed
        arrays in the bundle are re-indexed to match.

        Parameters
        ----------
        b : TFBundle
            Bundle to reorder.

        Returns
        -------
        TFBundle
            The same bundle, ordered in place.

        Notes
        -----
        If ``freq`` is ``None`` the bundle is returned unchanged.
        """

        if b.freq is None:
            return b
        order = get_config().freq_order
        if np is not None:
            arr = np.asarray(b.freq)
            idx = np.argsort(arr)
            if order == "desc":
                idx = idx[::-1]
            b.freq = arr[idx]
            if b.z is not None:
                b.z = b.z[idx]
            if b.z_err is not None:
                b.z_err = b.z_err[idx]
            if b.tipper is not None:
                b.tipper = b.tipper[idx]
            if b.tipper_err is not None:
                b.tipper_err = b.tipper_err[idx]
            if b.rho is not None:
                b.rho = b.rho[idx]
            if b.phase is not None:
                b.phase = b.phase[idx]
            return b
        seq = list(zip(b.freq, range(len(b.freq))))
        seq.sort(key=lambda x: x[0], reverse=(order == "desc"))
        take = [i for _, i in seq]
        b.freq = [b.freq[i] for i in take]
        if b.z is not None:
            b.z = [b.z[i] for i in take]
        if b.z_err is not None:
            b.z_err = [b.z_err[i] for i in take]
        if b.tipper is not None:
            b.tipper = [b.tipper[i] for i in take]
        if b.tipper_err is not None:
            b.tipper_err = [b.tipper_err[i] for i in take]
        if b.rho is not None:
            b.rho = [b.rho[i] for i in take]
        if b.phase is not None:
            b.phase = [b.phase[i] for i in take]
        return b

    def _dedup_freq(self, b: TFBundle) -> TFBundle:
        r"""
        De-duplicate nearly equal frequencies with a relative tol.

        Frequencies are scanned in their current order. A sample is
        kept when its relative difference to the last kept value is
        greater than ``freq_tol`` in the global config. All
        frequency-indexed arrays are subsetted consistently.

        Parameters
        ----------
        b : TFBundle
            Bundle to de-duplicate.

        Returns
        -------
        TFBundle
            The same bundle, with potential rows removed.

        Notes
        -----
        This step expects an already ordered ``freq``. Provide
        consistent types (NumPy or list) in the bundle.
        """

        if b.freq is None:
            return b
        tol = get_config().freq_tol
        if np is not None:
            arr = np.asarray(b.freq, dtype=float)
            keep = [0]
            for i in range(1, arr.size):
                prev = arr[keep[-1]]
                if abs(arr[i] - prev) / max(prev, 1.0) > tol:
                    keep.append(i)
            idx = np.asarray(keep)
            b.freq = arr[idx]
            if b.z is not None:
                b.z = b.z[idx]
            if b.z_err is not None:
                b.z_err = b.z_err[idx]
            if b.tipper is not None:
                b.tipper = b.tipper[idx]
            if b.tipper_err is not None:
                b.tipper_err = b.tipper_err[idx]
            if b.rho is not None:
                b.rho = b.rho[idx]
            if b.phase is not None:
                b.phase = b.phase[idx]
            return b
        # fallback list mode
        out_f = []
        take = []
        for i, f in enumerate(b.freq):
            if not out_f:
                out_f.append(f)
                take.append(i)
                continue
            prev = out_f[-1]
            denom = prev if prev else 1.0
            if abs(f - prev) / denom > tol:
                out_f.append(f)
                take.append(i)
        b.freq = out_f
        if b.z is not None:
            b.z = [b.z[i] for i in take]
        if b.z_err is not None:
            b.z_err = [b.z_err[i] for i in take]
        if b.tipper is not None:
            b.tipper = [b.tipper[i] for i in take]
        if b.tipper_err is not None:
            b.tipper_err = [b.tipper_err[i] for i in take]
        if b.rho is not None:
            b.rho = [b.rho[i] for i in take]
        if b.phase is not None:
            b.phase = [b.phase[i] for i in take]
        return b

    def _fill_missing(self, b: TFBundle) -> TFBundle:
        r"""
        Fill missing TF parts according to global configuration.

        If the bundle has ``Z`` but lacks ``(rho, phase)`` and
        ``compute_res_from_z`` is enabled, this calls
        :meth:`compute_res_from_z`. Conversely, if it has
        ``(rho, phase)`` but lacks ``Z`` and
        ``compute_z_from_res`` is enabled, this calls
        :meth:`compute_z_from_res`.

        Parameters
        ----------
        b : TFBundle
            Bundle to enrich.

        Returns
        -------
        TFBundle
            The updated bundle.
        """

        cfg = get_config()
        have_z = b.z is not None
        have_rp = (b.rho is not None) and (b.phase is not None)
        if cfg.compute_res_from_z and have_z and not have_rp:
            b = self.compute_res_from_z(b)
        if cfg.compute_z_from_res and not have_z and have_rp:
            b = self.compute_z_from_res(b)
        return b

    def _finalize(
        self,
        b: TFBundle,
        *,
        name: str | None = None,
        station_id: str | int | None = None,
    ) -> TFBundle:
        r"""
        Finalize a bundle: name, fill, order, and de-duplicate.

        Steps:

        1. Ensure a valid station name using the global policy.
        2. Fill missing parts using :meth:`_fill_missing`.
        3. Order frequencies per :meth:`_order_freq`.
        4. De-duplicate frequencies per :meth:`_dedup_freq`.

        Parameters
        ----------
        b : TFBundle
            Bundle to finalize.
        name : str, optional
            Preferred station name override.
        station_id : str or int, optional
            Identifier used to synthesize a name if needed.

        Returns
        -------
        TFBundle
            The finalized bundle.

        See Also
        --------
        pycsamt.core.config.StationNamePolicy
            Rules used to validate/synthesize names.
        """

        sid = b.station_id if station_id is None else station_id
        b.station = ensure_station(name or b.station, sid)
        b = self._fill_missing(b)
        b = self._order_freq(b)
        b = self._dedup_freq(b)
        return b

    def transform(
        self,
        source: Any,
        *,
        name: str | None = None,
        station_id: str | int | None = None,
    ) -> Any:
        r"""
        One-shot transform: extract → finalize → emit → post.

        This orchestrates the whole conversion pipeline:

        1. :meth:`extract` to build a :class:`TFBundle`.
        2. :meth:`_finalize` to normalize the bundle.
        3. :meth:`emit_edi` to create the EDI object.
        4. :meth:`post_emit` for optional enrichment.

        Parameters
        ----------
        source : Any
            Source object or file path.
        name : str, optional
            Preferred station name override.
        station_id : str or int, optional
            Identifier used if a name must be synthesized.

        Returns
        -------
        Any
            The resulting EDI object (or compatible subtype).

        Examples
        --------
        >>> from pycsamt.core._transformers import TransformerMixin
        >>> class Mini(TransformerMixin):
        ...     def extract(self, s):
        ...         from pycsamt.core.base import TFBundle
        ...         return TFBundle(freq=[1.], rho=[100.], phase=[45.])
        ...     def emit_edi(self, b):
        ...         return {'ok': True, 'n': len(b.freq)}
        ...
        >>> Mini().transform(object())['ok']
        True
        """

        b = self.extract(source)
        b = self._finalize(
            b,
            name=name,
            station_id=station_id,
        )
        edi = self.emit_edi(b)
        return self.post_emit(edi, source, b)

    def _ensure_head(self, ed, station: str, empty: float):
        h = ed.get_section("head")
        if h is None:
            h = Head()
            h.dataid = (station or "").strip()
            h.stdvers = "SEG 1.0"
            h.progvers = "PYCSAMT"
            h.empty = float(empty)
            h.lat = 0.0
            h.long = 0.0
            h.elev = 0.0
            ed.add_section("head", h)

        return h


    def _ensure_info(self, ed, survey_id: str):
        info = ed.get_section("info")
        if info is None:
            info = Info()
            info.update(
                info_text=[
                    f"  SURVEY ID:{survey_id or ''}",
                    "  ROTATION=FIX",
                ]
            )
            ed.add_section("info", info)
        return info


    def _ensure_definemeas(
        self,
        ed,
        *,
        units: str = "M",
        reftype: str = "CART",
    ):
        dm = ed.get_section("definemeas")
        if dm is None:
            dm = DefineMeas()
            dm.units = units
            dm.reftype = reftype
            h = ed.get_section("head")
            if h is not None:
                dm.reflat = getattr(h, "lat", 0.0)
                dm.reflong = getattr(h, "long", 0.0)
                dm.refelev = getattr(h, "elev", 0.0)
            ed.add_section("definemeas", dm)

        return dm

    def _ensure_mtsect(self, ed, sectid: str, nfreq: int):
        mt = ed.get_section("mtsect")
        if mt is None:
            mt = MTEMAP()
            ed.add_section("mtsect", mt)
        mt.sectid = sectid
        mt.nfreq = int(nfreq)
        return mt
