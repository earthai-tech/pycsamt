# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later

from __future__ import annotations

from collections.abc import Iterable
from typing import Any

__all__ = ["ComponentsMixin"]


class ComponentsMixin:
    r"""
    Lightweight helpers to access and manage structured
    components stored in ``self.sections``.  The host class
    is expected to provide a minimal mapping-like API:

    * ``add_section(key, obj)``
    * ``get_section(key)``
    * ``has_section(key)``

    Keys are normalized to lowercase to keep lookups stable
    across callers and file variants.

    Notes
    -----
    This mixin centralizes common section plumbing used by
    EDI readers and writers.  It exposes both generic CRUD
    helpers and a small set of typed accessors for frequent
    sections (e.g. ``head``, ``info``, spectra, and time
    series).

    The following utilities are provided:

    * ``cget(key, default=None)`` :
      Safe fetch from ``self.sections`` with key
      normalization and a default.

    * ``cset(key, obj)`` :
      Insert or replace a component under the normalized
      key.

    * ``chas(key)`` :
      Presence check under the normalized key.

    * ``cdrop(key)`` :
      Remove a component if present.  Missing keys are
      ignored.

    * ``snapshot(keys=None)`` :
      Return a dict containing a shallow snapshot of the
      requested components.  When ``keys`` is ``None`` all
      components are included.

    * ``compose_headers_from(keys_order=None)`` :
      Serialize a subset of sections by calling each
      component's ``write()`` method, if available.  A
      component may return either ``list[str]`` or
      ``str``.  Missing sections, components without a
      ``write()`` method, or components whose ``write()``
      raises are safely ignored.  The result is a single
      concatenated ``str``.  If ``keys_order`` is given,
      it controls the emission order; otherwise the mixin
      uses a sensible default order for EDI files.

    Typed accessors mirror common sections and include
    small quality-of-life fallbacks:

    * ``get_head()``, ``set_head(obj)``
    * ``get_info()``, ``set_info(obj)``
    * ``get_definemeasurement()``, ``set_definemeasurement(obj)``
      with a fallback that also checks the key
      ``'definemeas'`` when the canonical key is absent.
    * ``get_spectra_sect()``, ``set_spectra_sect(obj)``
      and ``get_spectra_io()``, ``set_spectra_io(obj)``,
      along with ``get_spectra()`` / ``set_spectra(obj)``.
      When a specific spectra header key is missing, the
      generic ``'spectra'`` key is consulted.
    * ``get_timeseries_sect()``, ``set_timeseries_sect(obj)``
      and ``get_timeseries_io()``, ``set_timeseries_io(obj)``,
      along with ``get_timeseries()`` / ``set_timeseries(obj)``.
      As with spectra, a generic ``'timeseries'`` key is used
      as a fallback for header lookups.

    These helpers keep parsing and writing code concise while
    remaining tolerant to small input variations.

    Examples
    --------
    Create a tiny host and manage sections::

        class Host(ComponentsMixin):
            def __init__(self):
                self.sections = {}
            def add_section(self, k, v):
                self.sections[k] = v
            def get_section(self, k):
                return self.sections.get(k)
            def has_section(self, k):
                return k in self.sections

        h = Host()
        h.cset("HEAD", {"dataid": "S01"})
        assert h.chas("head")
        assert h.cget("Head")["dataid"] == "S01"

        snap = h.snapshot(keys=["head", "info"])
        # {'head': {...}, 'info': None}

    Compose headers from components that implement a
    ``write()`` method::

        class _Head:
            def write(self):
                return [">HEAD\n", "  DATAID=S01\n"]
        class _Info:
            def write(self):
                return ">INFO\n  PROJECT=P\n"

        h.cset("head", _Head())
        h.cset("info", _Info())
        text = h.compose_headers_from(keys_order=["info", "head"])
        # INFO appears before HEAD in the output.

    See Also
    --------
    pycsamt.seg.cbase.ParseMixin :
        Source discovery and normalization utilities.
    pycsamt.seg.collection.EDICollection :
        Higher-level container that aggregates parsed EDI
        files and uses this mixin for section handling.

    References
    ----------
    .. [1] SEG EDI MT/EMAP standard (1987).  MTNet.
       https://www.mtnet.info/docs/seg_mt_emap_1987.pdf
    """

    # ---------- generic helpers ----------

    @staticmethod
    def _norm(k: str) -> str:
        return str(k).lower()

    def cget(self, key: str, default: Any = None) -> Any:
        k = self._norm(key)
        get = getattr(self, "get_section", None)
        if callable(get):
            val = get(k)
            return default if val is None else val
        # fallback if host uses `sections` directly
        sec = getattr(self, "sections", None)
        if isinstance(sec, dict):
            return sec.get(k, default)
        return default

    def cset(self, key: str, obj: Any) -> None:
        k = self._norm(key)
        add = getattr(self, "add_section", None)
        if callable(add):
            add(k, obj)
            return
        sec = getattr(self, "sections", None)
        if isinstance(sec, dict):
            sec[k] = obj

    def chas(self, key: str) -> bool:
        k = self._norm(key)
        has = getattr(self, "has_section", None)
        if callable(has):
            return bool(has(k))
        sec = getattr(self, "sections", None)
        if isinstance(sec, dict):
            return k in sec
        return False

    def cdrop(self, key: str) -> None:
        k = self._norm(key)
        sec = getattr(self, "sections", None)
        if isinstance(sec, dict) and k in sec:
            del sec[k]

    def snapshot(self, keys: Iterable[str] | None = None
                 ) -> dict[str, Any]:
        """
        Return a shallow snapshot mapping key -> component or
        ``None`` if absent.  If ``keys`` is omitted, use all
        current keys.
        """
        sec = getattr(self, "sections", None)
        if not isinstance(sec, dict):
            return {}
        if keys is None:
            keys = list(sec.keys())
        out: dict[str, Any] = {}
        for k in keys:
            kk = self._norm(k)
            out[kk] = sec.get(kk, None)
        return out

    # ---------- typed getters ----------

    # Core structural sections
    def get_head(self) -> Any:
        return self.cget("head")

    def set_head(self, obj: Any) -> None:
        self.cset("head", obj)

    def get_info(self) -> Any:
        return self.cget("info")

    def set_info(self, obj: Any) -> None:
        self.cset("info", obj)

    def get_definemeasurement(self) -> Any:
        # support both "definemeasurement" and "definemeas"
        dm = self.cget("definemeasurement")
        if dm is None:
            dm = self.cget("definemeas")
        return dm

    def set_definemeasurement(self, obj: Any) -> None:
        self.cset("definemeasurement", obj)

    # MT/EMAP section header
    def get_mtsect(self) -> Any:
        return self.cget("mtsect")

    def set_mtsect(self, obj: Any) -> None:
        self.cset("mtsect", obj)

    # OTHER section + IO
    def get_other(self) -> Any:
        return self.cget("other")

    def set_other(self, obj: Any) -> None:
        self.cset("other", obj)

    def get_other_io(self) -> Any:
        return self.cget("otherio")

    def set_other_io(self, obj: Any) -> None:
        self.cset("otherio", obj)

    # ---------- spectra components ----------

    def get_spectra_sect(self) -> Any:
        # when read via iter_sections: key is "spectra"
        # for clarity expose a separate accessor alias
        s = self.cget("spectra_sect")
        if s is None:
            s = self.cget("spectra")
        return s

    def set_spectra_sect(self, obj: Any) -> None:
        self.cset("spectra_sect", obj)

    def get_spectra_io(self) -> Any:
        return self.cget("spectra_io")

    def set_spectra_io(self, obj: Any) -> None:
        self.cset("spectra_io", obj)

    def get_spectra(self) -> Any:
        return self.cget("spectra")

    def set_spectra(self, obj: Any) -> None:
        self.cset("spectra", obj)

    # ---------- time-series components ----------

    def get_timeseries_sect(self) -> Any:
        s = self.cget("timeseries_sect")
        if s is None:
            s = self.cget("timeseries")
        return s

    def set_timeseries_sect(self, obj: Any) -> None:
        self.cset("timeseries_sect", obj)

    def get_timeseries_io(self) -> Any:
        return self.cget("timeseries_io")

    def set_timeseries_io(self, obj: Any) -> None:
        self.cset("timeseries_io", obj)

    def get_timeseries(self) -> Any:
        return self.cget("timeseries")

    def set_timeseries(self, obj: Any) -> None:
        self.cset("timeseries", obj)

    # ---------- composition helpers ----------

    def compose_headers_from(
        self, keys_order: list[str] | None = None
    ) -> str:
        """
        Compose textual headers by calling ``write()`` on
        registered section objects.  The order can be
        controlled with ``keys_order``; unknown keys are
        ignored; remaining keys are appended by insertion
        order.
        """
        sec = getattr(self, "sections", None)
        if not isinstance(sec, dict):
            return ""

        out: list[str] = []
        seen: set[str] = set()

        def _emit_one(k: str) -> None:
            kk = self._norm(k)
            if kk in seen:
                return
            obj = sec.get(kk, None)
            if obj is None:
                return
            write = getattr(obj, "write", None)
            if callable(write):
                try:
                    lines = write()
                    if isinstance(lines, list):
                        out.extend(lines)
                    elif isinstance(lines, str):
                        out.append(lines)
                except Exception:
                    # be tolerant in composition
                    pass
            seen.add(kk)

        if keys_order:
            for k in keys_order:
                _emit_one(k)

        # append any remaining known sections
        for k in sec.keys():
            _emit_one(k)

        return "".join(out)

