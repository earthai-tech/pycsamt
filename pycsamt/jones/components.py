# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later

from __future__ import annotations

from typing import Any, Dict, Iterable, List, Optional


__all__ = ["JComponentMixin"]


class JComponentMixin:
    r"""
    Lightweight facade to manage Jones (J-format) components on a
    host object.  The host is expected to expose a mapping-like
    container (``self.sections`` or ``self.components``).  Keys are
    normalized to lowercase for stable lookups.

    The mixin centralizes small CRUD helpers (``cget/cset/chas/
    cdrop/snapshot``) and offers typed accessors for frequent J
    parts such as ``banner``, ``info``, ``head``, ``heads``,
    ``blocks``, and derived data objects (``z``, ``tip``, ``res``).

    It also provides ``compose_headers`` which serializes the
    header portion of a J file by calling ``write()`` on available
    components.  If a combined ``heads`` component is present, it
    takes precedence over separate ``banner``/``head``/``info``
    parts to avoid duplication.
    """
    @staticmethod
    def _norm(key: str) -> str:
        return str(key).strip().lower()

    def _bag(self) -> Dict[str, Any]:
        # Prefer `sections`, fall back to `components`.
        bag = getattr(self, "sections", None)
        if not isinstance(bag, dict):
            bag = getattr(self, "components", None)
        if not isinstance(bag, dict):
            bag = {}
            setattr(self, "sections", bag)
        return bag

    def cget(self, key: str, default: Any = None) -> Any:
        return self._bag().get(self._norm(key), default)

    def cset(self, key: str, obj: Any) -> None:
        self._bag()[self._norm(key)] = obj

    def chas(self, key: str) -> bool:
        return self._norm(key) in self._bag()

    def cdrop(self, key: str) -> None:
        self._bag().pop(self._norm(key), None)

    def snapshot(
        self, keys: Optional[Iterable[str]] = None
    ) -> Dict[str, Any]:
        bag = self._bag()
        if keys is None:
            keys = list(bag.keys())
        out: Dict[str, Any] = {}
        for k in keys:
            kk = self._norm(k)
            out[kk] = bag.get(kk, None)
        return out

    # Banner
    def get_banner(self) -> Any:
        return self.cget("banner")

    def set_banner(self, obj: Any) -> None:
        self.cset("banner", obj)

    # Info (>KEY=VALUE)
    def get_info(self) -> Any:
        return self.cget("info")

    def set_info(self, obj: Any) -> None:
        self.cset("info", obj)

    # Single head (station / dtype / count)
    def get_head(self) -> Any:
        return self.cget("head")

    def set_head(self, obj: Any) -> None:
        self.cset("head", obj)

    # Combined Heads facade (banner + head + info)
    def get_heads(self) -> Any:
        return self.cget("heads")

    def set_heads(self, obj: Any) -> None:
        self.cset("heads", obj)

    # Data blocks collection
    def get_blocks(self) -> Any:
        return self.cget("blocks")

    def set_blocks(self, obj: Any) -> None:
        self.cset("blocks", obj)

    # Derived data objects (from pycsamt.z.*)
    def get_z(self) -> Any:
        # impedance tensor
        return self.cget("z")

    def set_z(self, obj: Any) -> None:
        self.cset("z", obj)

    def get_tip(self) -> Any:
        # tipper
        return self.cget("tip")

    def set_tip(self, obj: Any) -> None:
        self.cset("tip", obj)

    def get_res(self) -> Any:
        # resistivity/phase
        return self.cget("res")

    def set_res(self, obj: Any) -> None:
        self.cset("res", obj)

    def compose_headers(
        self, *, prefer_heads: bool = True, join: bool = False
    ) -> List[str] | str:
        """
        Serialize header lines by querying available components.
        If ``prefer_heads`` is True and a ``heads`` component is
        present, it is used exclusively.  Otherwise the mixin
        attempts ``banner`` + ``head`` + ``info`` in that order.

        When ``join`` is True, returns a single string joined
        with newlines.  Otherwise returns a list of lines.
        """
        out: List[str] = []

        def _extend(obj: Any) -> None:
            if obj is None:
                return
            write = getattr(obj, "write", None)
            if not callable(write):
                return
            try:
                lines = write()
            except Exception:
                return
            if isinstance(lines, list):
                out.extend([str(s) for s in lines])
            elif isinstance(lines, str):
                # accept already-joined text
                out.extend(str(lines).splitlines())

        bag = self._bag()

        if prefer_heads and "heads" in bag:
            _extend(bag.get("heads"))
        else:
            _extend(bag.get("banner"))
            _extend(bag.get("head"))
            _extend(bag.get("info"))

        if join:
            return "\n".join(out)
        return out
