# pycsamt/_session.py 
# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later

from __future__ import annotations 

from pathlib import Path 
from typing import Any, Optional  

from .core.transformers import AVGtoEDI, JtoEDI 
from .core.registry import RegistryAPI
from .core.base import to_edi 

from .jones.j import JFile
from .zonge.avg import AVG
from .seg.edi import EDIFile
from .seg.collection import EDICollection


__all__ = [
    "Session",
    "work_session",
    "Normalize", 
    "normalize_session"

]


class Session:
    def __init__(
        self,
        root: Path | str,
        *,
        manifest_name: str = "manifest.json",
        auto_capture: bool = True,
        capture_children: bool = False,
        max_children: int = 256,
    ) -> None:
        self.root = Path(root)
        self.auto_capture = bool(auto_capture)
        self.capture_children = bool(capture_children)
        self.max_children = int(max_children)
        self.reg = RegistryAPI(
            self.root,
            manifest_name=manifest_name,
        )
        self._orig_to_edi = None
        self._orig_t_x = None
        self._orig_t_j = None


    def _record(self, obj: Any, *, tag: str) -> None:
        try:
            self.reg.add_object(obj, tags=[tag])
        except Exception:
            pass
        if not self.capture_children:
            return
        try:
            if hasattr(obj, "__iter__") and not hasattr(obj, "Z"):
                n = 0
                for it in obj:  # type: ignore
                    self.reg.add_object(it, tags=[tag])
                    n += 1
                    if n >= self.max_children:
                        break
        except Exception:
            pass

    def _wrap_to_edi(self) -> None:
        
        from .core import base as b
        
        if not self.auto_capture:
            return

        if self._orig_to_edi is None:
            self._orig_to_edi = to_edi

        def _wrapped(source: Any, *a: Any, **k: Any) -> Any:
            out = self._orig_to_edi(source, *a, **k)
            self._record(out, tag="to_edi")
            return out

        b.to_edi = _wrapped  # type: ignore

    def _wrap_transformers(self) -> None:
        from .core import transformers as tr
        
        if not self.auto_capture:
            return

        def _wrap_method(func):
            def inner(self_, *a, **k):  # noqa: ANN001
                out = func(self_, *a, **k)
                try:
                    tag = f"{self_.__class__.__name__}.transform"
                except Exception:
                    tag = "transform"
                self._record(out, tag=tag)
                return out
            return inner

        if self._orig_t_x is None and hasattr(tr, "AVGtoEDI"):
            self._orig_t_x = tr.AVGtoEDI.transform
            tr.AVGtoEDI.transform = _wrap_method(self._orig_t_x)  # type: ignore
        if self._orig_t_j is None and hasattr(tr, "JtoEDI"):
            self._orig_t_j = tr.JtoEDI.transform
            tr.JtoEDI.transform = _wrap_method(self._orig_t_j)  # type: ignore

    def _restore(self) -> None:
        from .core import base as b
        from .core import transformers as tr
        try:
            if self._orig_to_edi is not None:
                b.to_edi = self._orig_to_edi  # type: ignore
        except Exception:
            pass
        try:
            
            if self._orig_t_x is not None:
                tr.AVGtoEDI.transform = self._orig_t_x  # type: ignore
            if self._orig_t_j is not None:
                tr.JtoEDI.transform = self._orig_t_j  # type: ignore
        except Exception:
            pass

    def __enter__(self) -> "Session":
        self._wrap_to_edi()
        self._wrap_transformers()
        return self

    def __exit__(
        self,
        exc_type,
        exc,
        tb,
    ) -> Optional[bool]:
        self._restore()
        try:
            self.reg.save()
        except Exception:
            pass
        return None


def work_session(
    root: Path | str,
    *,
    manifest_name: str = "manifest.json",
    auto_capture: bool = True,
    capture_children: bool = False,
    max_children: int = 256,
) -> Session:
    return Session(
        root,
        manifest_name=manifest_name,
        auto_capture=auto_capture,
        capture_children=capture_children,
        max_children=max_children,
    )


class Normalize:
    def __init__(
        self,
        root: Path | str,
        *,
        manifest_name: str = "manifest.json",
        topo_src: Any | None = None,
        auto_register: bool = True,
    ) -> None:
        self.root = Path(root)
        self.reg = RegistryAPI(
            self.root,
            manifest_name=manifest_name,
        )
        self.topo_src = topo_src
        self.auto_register = bool(auto_register)

    def _as_edi_coll(self, src: Any) -> Any:

        if isinstance(src, EDICollection):
            return src
        if isinstance(src, EDIFile):
            return EDICollection(items=[src], verbose=0)
        return None

    def _try_topo(self, avg: Any) -> None:
        if self.topo_src is None:
            return
        try:
            if hasattr(avg, "add_topography"):
                avg.add_topography(self.topo_src)
                return
        except Exception:
            pass
        try:
            if hasattr(self.topo_src, "frame"):
                avg.topo = self.topo_src
        except Exception:
            pass

    def _to_avg(self, src: Any) -> Any:
        
        
       
        if isinstance(src, AVG):
            return src
        if isinstance(src, (str, Path)):
            return AVG.from_file(src)
        return None

    def _to_j(self, src: Any) -> Any:
        
        

        if isinstance(src, JFile):
            return src
        if isinstance(src, (str, Path)):
            return JFile.from_file(src)
        return None

    def _normalize(self, source: Any) -> Any:
        out = self._as_edi_coll(source)
        if out is not None:
            return out
        a = self._to_avg(source)
        if a is not None:
            self._try_topo(a)
            
            return AVGtoEDI().transform(a)
        j = self._to_j(source)
        if j is not None:
            return JtoEDI().transform(j)
        try:
            
         
            if isinstance(source, (str, Path)):
                ed = EDIFile.from_file(source)
                return EDICollection(items=[ed], verbose=0)
        except Exception:
            pass
        
        obj = to_edi(source)
        return self._as_edi_coll(obj) or obj

    def __enter__(self) -> "Normalize":
        return self

    def __exit__(
        self,
        exc_type,
        exc,
        tb,
    ) -> bool | None:
        try:
            self.reg.save()
        except Exception:
            pass
        return None

    def load(self, source: Any) -> Any:
        out = self._normalize(source)
        if self.auto_register:
            try:
                self.reg.add_object(out, tags=["normalized"])
            except Exception:
                pass
        return out

def normalize_session(
    root: Path | str,
    *,
    manifest_name: str = "manifest.json",
    topo_src: Any | None = None,
    auto_register: bool = True,
):
    return Normalize(
        root,
        manifest_name=manifest_name,
        topo_src=topo_src,
        auto_register=auto_register,
    )
