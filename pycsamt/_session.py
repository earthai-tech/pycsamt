# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

from __future__ import annotations 

from pathlib import Path 
from typing import Any, Optional  

from .core.transformers import AVGtoEDI, JtoEDI 
from .core.registry import RegistryAPI
from .core.base import to_edi, CoreObject 
from .jones.j import JFile
from .jones.collection import JCollection
from .zonge.avg import AVG
from .seg.edi import EDIFile
from .seg.collection import EDICollection


__all__ = [
    "Session",
    "work_session",
    "Normalize", 
    "normalize_session"

]


class Session(CoreObject):
    r"""
    Context-managed capture of key transformation outputs.
    
    A :class:`Session` wraps selected conversion utilities and
    transformers so that their **results** are automatically
    registered in a small on-disk registry. This enables light
    provenance tracking of what you produced during a workflow
    (e.g., the output from :func:`to_edi` or from
    :class:`AVGtoEDI` / :class:`JtoEDI`).
    
    When the session starts, it temporarily monkey-patches the
    following call sites to intercept and register outputs:
    
    * :func:`pycsamt.core.base.to_edi`
    * :meth:`pycsamt.core.transformers.AVGtoEDI.transform`
    * :meth:`pycsamt.core.transformers.JtoEDI.transform`
    
    On exit, the original functions are restored and the manifest
    is saved.
    
    Parameters
    ----------
    root : path-like
        Working directory that will hold a registry manifest and
        any artifacts you add explicitly via the underlying
        :class:`~pycsamt.core.registry.RegistryAPI`.
    manifest_name : str, optional
        Manifest filename under ``root``. Defaults to
        ``"manifest.json"``.
    auto_capture : bool, optional
        If ``True`` (default), wrap the functions above so that
        their **outputs** are registered with tags describing the
        source call (e.g., ``"AVGtoEDI.transform"``).
    capture_children : bool, optional
        If ``True``, when a captured output is iterable (e.g., a
        collection), the session attempts to also register each
        child item, up to ``max_children``. Defaults to
        ``False``.
    max_children : int, optional
        Maximum number of iterable children to capture when
        ``capture_children=True``. Defaults to ``256``.
    
    Attributes
    ----------
    root : pathlib.Path
        Absolute path to the working directory.
    reg : pycsamt.core.registry.RegistryAPI
        High-level registry façade used to record artifacts.
    
    Notes
    -----
    The session captures **results**, not inputs. To store inputs
    explicitly, call :meth:`reg.add_file` or
    :meth:`reg.add_object` on the session's registry.
    
    The capture is process-local and uses monkey-patching. In
    multi-threaded or multi-process settings, consider isolate
    usage or add explicit locking to your code.
    
    Examples
    --------
    Basic usage with explicit object registration::
    
        >>> from pycsamt.session import Session
        >>> from pycsamt.seg.edi import EDIFile
        >>> with Session("work") as ses:
        ...     edi = EDIFile.from_file("data/site001.edi")
        ...     _ = ses.reg.add_object(edi, tags=["raw"])
        ...     out = ses.reg.list()
        >>> isinstance(out, list)
        True
    
    Automatic capture from transformers::
    
        >>> from pycsamt.session import Session
        >>> from pycsamt.zonge.avg import AVG
        >>> from pycsamt.core.transformers import AVGtoEDI
        >>> with Session("work", auto_capture=True) as ses:
        ...     avg = AVG.from_file("data/S01.AVG")  # doctest: +SKIP
        ...     edi = AVGtoEDI().transform(avg)      # doctest: +SKIP
        ...     # The result of transform is recorded automatically.
        ...     got = ses.reg.find(tag="AVGtoEDI.transform")  # doctest: +SKIP
    
    Automatic capture from :func:`to_edi`::
    
        >>> from pycsamt.session import Session
        >>> from pycsamt.core.base import to_edi
        >>> from pycsamt.jones.j import JFile
        >>> with Session("work", auto_capture=True) as ses:
        ...     j = JFile.from_file("data/site001.j")  # doctest: +SKIP
        ...     edi_like = to_edi(j)                   # doctest: +SKIP
        ...     hits = ses.reg.find(tag="to_edi")      # doctest: +SKIP
    
    Capture children from iterable outputs (e.g., collections)::
    
        >>> from pycsamt import session as pcs  # alias form
        >>> from pycsamt.seg.collection import EDICollection
        >>> with pcs.Session("work",
        ...                  auto_capture=True,
        ...                  capture_children=True,
        ...                  max_children=8) as ses:
        ...     # Suppose `coll` is an iterable of EDI-like items.
        ...     coll = EDICollection([])  # example only
        ...     ses.reg.add_object(coll, tags=["batch"])
        ...     # Children of `coll` would be recorded when captured.
    
    Inspect what was recorded and persist::
    
        >>> with Session("work") as ses:
        ...     _ = ses.reg.list()
        ...     # Saving also happens on context exit.
        >>> # The manifest `work/manifest.json` now reflects updates.
    
    See Also
    --------
    pycsamt.core.registry.RegistryAPI
    pycsamt.core.base.to_edi
    pycsamt.core.transformers.AVGtoEDI
    pycsamt.core.transformers.JtoEDI
    pycsamt.seg.edi.EDIFile
    pycsamt.seg.collection.EDICollection
    pycsamt.zonge.avg.AVG
    pycsamt.jones.j.JFile
    
    References
    ----------
    .. [1] Python Stdlib. *context managers* (`with` statement).
    .. [2] NumPy. *np.savez_compressed* used by registry packers.
    """

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
        r"""
        Register ``obj`` in the session registry with a tag.
        
        Optionally iterates over ``obj`` to capture a
        bounded number of children when ``capture_children=True``.
        """

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
        r"""
        Enable capture for :func:`to_edi` within this process.
 
        Monkey-patches the public ``to_edi`` so that
        its **outputs** are recorded with the tag ``"to_edi"``.
        """

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
        r"""
        Enable capture for AVG→EDI and J→EDI transformer outputs.

        Wraps ``AVGtoEDI.transform`` and ``JtoEDI.transform`` so 
        their return values are recorded with class-qualified tags.
        """

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
        r"""
        Restore original call sites and stop interception.
        
        Best-effort restoration; failures are
        silently ignored to keep teardown robust.
        """

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
    r"""
    Create and return a :class:`Session` with common defaults.
    
    This is a convenience wrapper around :class:`Session`
    mirroring its constructor. See :class:`Session` for details
    on automatic capture and child registration.
    
    Parameters
    ----------
    root : path-like
        Working directory for the session. Created if missing.
    manifest_name : str, optional
        Registry manifest filename. Defaults to ``"manifest.json"``.
    auto_capture : bool, optional
        Enable interception of selected transforms to register
        their outputs. Defaults to ``True``.
    capture_children : bool, optional
        When ``True``, iterable outputs are traversed and child
        items are registered up to ``max_children``. Defaults to
        ``False``.
    max_children : int, optional
        Max number of children to register from one iterable
        output. Defaults to ``256``.
    
    Returns
    -------
    Session
        A context-manageable session instance.
    
    Examples
    --------
    Use the helper to run a small workflow::
    
        >>> from pycsamt.session import work_session
        >>> with work_session("work") as ses:
        ...     # Use ses.reg.* to add artifacts
        ...     items = ses.reg.list()
    
    Alias style::
    
        >>> import pycsamt as pcs
        >>> with pcs.work_session("work") as ses:
        ...     _ = ses.reg.list()
    
    See Also
    --------
    Session
    pycsamt.core.registry.RegistryAPI
    """

    return Session(
        root,
        manifest_name=manifest_name,
        auto_capture=auto_capture,
        capture_children=capture_children,
        max_children=max_children,
    )


class Normalize(CoreObject):
    r"""
    Normalize heterogeneous EM sources into a common form.
    
    ``Normalize`` provides a small, file-system-friendly
    facade that tries to convert various inputs (files and
    objects) into a standard, iterable representation, favoring
    :class:`~pycsamt.seg.collection.EDICollection` when possible.
    
    It accepts AVG files or objects, Jones J files, single EDI
    files, and already-built collections. When provided, a
    ``topo_src`` is attached to AVG sources before conversion.
    
    Optionally, normalized outputs are registered in a small
    manifest via :class:`~pycsamt.core.registry.RegistryAPI`.
    
    Parameters
    ----------
    root : path-like
        Working directory where a manifest is stored.
    manifest_name : str, optional
        Manifest filename under ``root``. Defaults to
        ``"manifest.json"``.
    topo_src : Any or None, optional
        Optional topography source to enrich AVG data before
        conversion. It may be an object accepted by
        ``AVG.add_topography`` or another object exposing a
        ``frame`` attribute to assign as ``avg.topo``.
    auto_register : bool, optional
        If ``True`` (default), normalized outputs are recorded in
        the registry with the tag ``"normalized"``.
    
    Attributes
    ----------
    root : pathlib.Path
        Absolute path to the working directory.
    reg : pycsamt.core.registry.RegistryAPI
        Registry used to persist normalized outputs.
    topo_src : Any or None
        Topography source passed at construction.
    auto_register : bool
        Whether :meth:`load` registers outputs automatically.
    
    Notes
    -----
    The normalization rules are:
    
    1. If the source is an :class:`EDICollection`, return it.
    2. If the source is an :class:`EDIFile`, wrap it into a
       one-item :class:`EDICollection`.
    3. If the source is an :class:`AVG` (or AVG file path), try
       to attach ``topo_src`` then transform with
       :class:`AVGtoEDI`.
    4. If the source is a :class:`JFile` (or J file path),
       transform with :class:`JtoEDI`.
    5. If the source is a string or path, try loading an
       :class:`EDIFile` and wrap into a collection.
    6. Otherwise, delegate to :func:`to_edi` and return an EDI-
       like object or a collection when possible.
    
    Errors while attaching topography or performing transforms
    are handled defensively; the method falls back to the next
    strategy or returns the best-effort result.
    
    Examples
    --------
    Normalize from an AVG file into an EDI collection::
    
        >>> from pycsamt.session import Normalize
        >>> with Normalize("work") as nz:
        ...     coll = nz.load("data/S01.AVG")  # doctest: +SKIP
        ...     # 'coll' is typically an EDICollection
    
    Attach topography before AVG→EDI::
    
        >>> topo = object()  # placeholder for a real topo source
        >>> with Normalize("work", topo_src=topo) as nz:  # doctest: +SKIP
        ...     coll = nz.load("data/S02.AVG")           # doctest: +SKIP
    
    Normalize from a J file path::
    
        >>> with Normalize("work") as nz:               # doctest: +SKIP
        ...     coll = nz.load("data/site001.j")        # doctest: +SKIP
    
    Already-normalized inputs pass through unchanged::
    
        >>> from pycsamt.seg.collection import EDICollection
        >>> with Normalize("work") as nz:
        ...     empty = EDICollection([])               # example
        ...     same = nz.load(empty)
        >>> isinstance(same, EDICollection)
        True
    
    Alias form via the package facade::
    
        >>> import pycsamt as pcs
        >>> with pcs.Normalize("work") as nz:           # doctest: +SKIP
        ...     out = nz.load("data/site001.edi")       # doctest: +SKIP
    
    See Also
    --------
    pycsamt.session.normalize_session
    pycsamt.core.registry.RegistryAPI
    pycsamt.core.transformers.AVGtoEDI
    pycsamt.core.transformers.JtoEDI
    pycsamt.core.base.to_edi
    pycsamt.seg.edi.EDIFile
    pycsamt.seg.collection.EDICollection
    pycsamt.zonge.avg.AVG
    pycsamt.jones.j.JFile
    
    References
    ----------
    .. [1] D. Wight (1991). *SEG MT/EMAP EDI Standard*.
    .. [2] Zonge International. *AVG data format*, tech notes.
    """

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
        r"""
        Return an :class:`EDICollection` view of ``src`` when possible.
    
        Rules
        -----
        * Pass through existing **non-empty** :class:`EDICollection`.
        * Wrap a single :class:`EDIFile` into a one-item collection.
        * If ``src`` is a path-like:
          - If it ends with ``.edi``, read with :meth:`EDIFile.from_file`
            and wrap to a one-item collection.
          - Otherwise (folder / glob / mix), try to construct an
            :class:`EDICollection` directly from ``sources=src``.
        * Tolerate raw ``list``/``tuple`` of :class:`EDIFile`` by
          normalizing to :class:`EDICollection`.
        * Return ``None`` when nothing valid can be produced (including
          empty collections).
        """
        # Pass through a non-empty collection
        if isinstance(src, EDICollection):
            return src if len(src) > 0 else None
    
        # Single EDI file -> one-item collection
        if isinstance(src, EDIFile):
            return EDICollection(items=[src], verbose=0)
    
        # Path-like input
        if isinstance(src, (str, Path)):
            p = Path(src)
            suf = p.suffix.lower()
    
            # Single .edi file → read and wrap
            if suf == ".edi":
                try:
                    ed = EDIFile.from_file(p)
                    coll = EDICollection(items=[ed], verbose=0)
                    return coll if len(coll) > 0 else None
                except Exception:
                    return None
    
            # Folder / glob / non-.edi path → let EDICollection load sources
            # Try the modern API first: EDICollection(sources=src, ...)
            try:
                coll = EDICollection(sources=src, verbose=0)  # real implementation
                return coll if len(coll) > 0 else None
            except TypeError:
                # Fallback for simpler stubs that don't accept 'sources'
                try:
                    coll = EDICollection(src, verbose=0)       # positional fallback
                    return coll if len(coll) > 0 else None
                except Exception:
                    return None
            except Exception:
                return None
    
        # Raw list/tuple of EDIFile -> normalize to EDICollection
        if isinstance(src, (list, tuple)) and src:
            items = [it for it in src if isinstance(it, EDIFile)]
            if not items:
                return None
            return EDICollection(items=items, verbose=0)
    
        return None

    def _try_topo(self, avg: Any) -> None:
        r"""
        Attach :attr:`topo_src` to an AVG object if supported.
        
        Calls ``avg.add_topography`` when available,
        otherwise tries to assign ``avg.topo`` from a ``frame`` attr.
        Silently ignores failures.
        """

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
        r"""
        Return an :class:`AVG` object from ``src`` or ``None``.
        
        Accepts an :class:`AVG` instance or a file
        path, which is read via :meth:`AVG.from_file`.
        """
        # 1) keep the instance check
        if isinstance(src, AVG):
            return src
        # Duck-typing fallback for environments where multiple 
        # copies of the class may exist
        # or for subclasses:

        if hasattr(src, "add_topography"):  # AVG has this
            # and this attribute exists
            return src
 
        # 2) only treat paths with AVG-like suffixes as AVG files
        if isinstance(src, (str, Path)):
            suf = Path(src).suffix.lower()
            if suf in {".avg"}:  # add others if you really support them
                return AVG.from_file(src)
        return None

    def _to_j(self, src: Any) -> Any:
        r"""
        Return a Jones object from ``src`` if it looks like Jones data.
    
        Rules
        -----
        * Accept a :class:`JFile` instance directly.
        * Accept a **non-empty** :class:`JCollection` (return ``None`` if empty).
        * If ``src`` is a path-like:
          - If it ends with a Jones-like suffix (``.j``, ``.jones``,
            ``.dat``, ``.txt``), try :meth:`JFile.from_file`.
          - Otherwise (folder / glob / mix), try to construct a
            :class:`JCollection` directly from ``sources=src``.
        * Tolerate raw ``list``/``tuple`` of :class:`JFile`` by normalizing
          to :class:`JCollection` (if available).
        * Return ``None`` when nothing valid can be produced (including
          empty collections).
        """
 
        # JFile instance → pass through
        if isinstance(src, JFile):
            return src
    
        # JCollection instance → pass through if non-empty
        if isinstance(src, JCollection):
            return src if len(src) > 0 else None
    
        # Path-like
        if isinstance(src, (str, Path)):
            p = Path(src)
            suf = p.suffix.lower()
    
            # Jones-like single-file suffix → try JFile.from_file
            if suf in {".j", ".jones", ".dat", ".txt"}:
                try:
                    return JFile.from_file(p)
                except Exception:
                    return None
    
            # Folder / glob / non-Jones-suffix → try JCollection(sources=...)
            try:
                coll = JCollection(sources=src, verbose=0)  # real API
                return coll if len(coll) > 0 else None
            except TypeError:
                # Fallback for stubs that don't accept 'sources='
                try:
                    coll = JCollection(src, verbose=0)       # positional
                    return coll if len(coll) > 0 else None
                except Exception:
                    return None
            except Exception:
                return None
    
        # Raw list/tuple of JFile → normalize to JCollection when available
        if isinstance(src, (list, tuple)) and src:
            items = [it for it in src if isinstance(it, JFile)]
            if not items:
                return None
            try:
                coll = JCollection(items=items, verbose=0)
            except TypeError:
                coll = JCollection(items)  # very simple stub fallback
            return coll if len(coll) > 0 else None
    
        return None

    def _normalize(self, source: Any) -> Any:
        r"""
        Core normalization strategy used by :meth:`load`.
    
        Resolution order
        ----------------
        1) Try to view ``source`` as an :class:`EDICollection` via
           :meth:`_as_edi_coll` (handles ``EDICollection``, ``EDIFile``,
           ``.edi`` file paths, folders, and glob patterns).
        2) Try Jones route via :meth:`_to_j` and
           :class:`~pycsamt.core.transformers.JtoEDI`.
        3) Try Zonge/AVG route via :meth:`_to_avg` (inject topography if
           configured) and :class:`~pycsamt.core.transformers.AVGtoEDI`.
        4) As a last chance for path-like inputs that slipped through,
           attempt to build an :class:`EDICollection` directly from the
           path (folder/glob).
        5) Fallback to :func:`to_edi` and normalize its result to an
           :class:`EDICollection` when possible.
    
        Returns
        -------
        Any
            An :class:`EDICollection` when possible, otherwise the object
            returned by :func:`to_edi`.
        """
        
        # 1) Immediate EDI pass-through (or single-file wrap)
        out = self._as_edi_coll(source)
        # If this is already an EDI object (not a string/Path), return it.
        if out is not None and not isinstance(source, (str, Path)):
            return out
    
        # 2) Path-like routing by suffix
        if isinstance(source, (str, Path)):
            p = Path(source)
            suf = p.suffix.lower()
    
            if suf in {".j", ".jones", ".dat", ".txt"}:
                j = self._to_j(source)
                if j is not None:
                    return JtoEDI().transform(j)
    
            if suf == ".avg":
                a = self._to_avg(source)
                if a is not None:
                    self._try_topo(a)
                    return AVGtoEDI().transform(a)
    
            if suf == ".edi":
                try:
                    ed = EDIFile.from_file(p)
                    return EDICollection(items=[ed], verbose=0)
                except Exception:
                    pass
    
            # No suffix (folder/glob) → try EDICollection(sources=...)
            if suf == "":
                try:
                    try:
                        coll = EDICollection(sources=source, verbose=0)
                    except TypeError:
                        coll = EDICollection(source, verbose=0)
                    if len(coll) > 0:
                        return coll
                except Exception:
                    pass
    
            # If we got here, path-like wasn’t resolved → fall through
    
        # 3) Non-path objects: Jones/AVG instance routes
        j = self._to_j(source)
        if j is not None:
            return JtoEDI().transform(j)
    
        a = self._to_avg(source)
        if a is not None:
            self._try_topo(a)
            return AVGtoEDI().transform(a)
    
        # 4) Ultimate fallback: adapter registry
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
        r"""
        Normalize ``source`` into a common, iterable EM representation.
        
        This method attempts to coerce the input into an
        :class:`EDICollection` when feasible. If the class was
        constructed with ``auto_register=True``, the resulting object
        is registered in the manifest with the tag ``"normalized"``.
        
        Parameters
        ----------
        source : Any
            Input to normalize. May be an :class:`EDICollection`,
            :class:`EDIFile`, :class:`AVG`, :class:`JFile`, a file
            path to one of those, or another EM object supported by
            :func:`to_edi`.
        
        Returns
        -------
        Any
            Usually an :class:`EDICollection`. In fallback cases, an
            EDI-like object produced by :func:`to_edi`.
        
        Notes
        -----
        Topography is attached to AVG inputs when possible, using
        :attr:`topo_src`. Failures are ignored, and the method moves
        on to the next strategy.
        
        Examples
        --------
        Basic usage returning an EDI collection::
        
            >>> from pycsamt.session import Normalize
            >>> with Normalize("work") as nz:
            ...     coll = nz.load("data/site001.edi")  # doctest: +SKIP
        
        Normalization from AVG with auto-registration::
        
            >>> from pycsamt.session import Normalize
            >>> with Normalize("work", auto_register=True) as nz:  # doctest: +SKIP
            ...     coll = nz.load("data/S03.AVG")                 # doctest: +SKIP
            ...     hits = nz.reg.find(tag="normalized")           # doctest: +SKIP
        
        Alias style via the package facade::
        
            >>> import pycsamt as pcs
            >>> with pcs.Normalize("work") as nz:                  # doctest: +SKIP
            ...     out = nz.load("data/site001.j")                # doctest: +SKIP
        
        See Also
        --------
        pycsamt.core.transformers.AVGtoEDI
        pycsamt.core.transformers.JtoEDI
        pycsamt.core.base.to_edi
        pycsamt.seg.collection.EDICollection
        """

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
    r"""
    Create and return a :class:`Normalize` with common defaults.
    
    This is a convenience wrapper around :class:`Normalize`,
    mirroring its constructor.
    
    Parameters
    ----------
    root : path-like
        Working directory.
    manifest_name : str, optional
        Manifest filename under ``root``. Defaults to
        ``"manifest.json"``.
    topo_src : Any or None, optional
        Optional topography source for AVG inputs.
    auto_register : bool, optional
        Register outputs from :meth:`Normalize.load`. Defaults to
        ``True``.
    
    Returns
    -------
    Normalize
        A context-manageable normalizer instance.
    
    Examples
    --------
    Use the helper to normalize a file quickly::
    
        >>> from pycsamt.session import normalize_session
        >>> with normalize_session("work") as nz:  # doctest: +SKIP
        ...     coll = nz.load("data/site001.edi") # doctest: +SKIP
    
    Alias form::
    
        >>> import pycsamt as pcs
        >>> with pcs.normalize_session("work") as nz:          # doctest: +SKIP
        ...     coll = nz.load("data/S01.AVG")                 # doctest: +SKIP
    
    See Also
    --------
    Normalize
    pycsamt.core.registry.RegistryAPI
    """

    return Normalize(
        root,
        manifest_name=manifest_name,
        topo_src=topo_src,
        auto_register=auto_register,
    )
