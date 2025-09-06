# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: GPL-3.0

"""
pycsamt.seg.collection

Batch helper around many :class:`~pycsamt.seg.edi.Edi` objects.

* handles either file paths **or** pre-instantiated ``Edi``;
* keeps long/lat/elev as *numpy* arrays for fast maths;
* exposes impedance / ρ / φ blocks through dict-like views;
* ships a compact human-readable ``__str__`` & ``__repr__``;
* provides quick stats and bulk-rewrite utilities.

"""
from __future__ import annotations

from pathlib import Path
from typing import Iterable, List, Mapping, Sequence

import numpy as np
import time 

from ..exceptions import FileHandlingError
from ..gis.utils import scale_position
from ..utils.func_utils import make_ids, scale_position as _sp_to_dms
from ..log.logger import get_logger

from .base import CollectionBase
from .edi import Edi
from .utils import quick_edi_stats, show_edi_stats

__all__ = ["EdiCollection"]

logger = get_logger(__name__)

class EdiCollection(CollectionBase[Edi]):
    r"""
    Container and helper for a *survey* of many
    :class:`~pycsamt.seg.edi.Edi` files.

    ``EdiCollection`` loads each site, extracts its key
    metadata and transfer-function arrays, and exposes
    them through convenient mappings keyed by *station
    id*.  The object is **immutable** – all data are held
    in RAM and the class provides only *read-only* views.
    A lightweight wrapper around :class:`~numpy.ndarray`
    is used internally so that vectorised operations
    remain possible without copy-on-write penalties.

    Parameters
    ----------
    paths : iterable of str or Path, optional
        A sequence of ``.edi`` paths.  Mutually exclusive
        with *objs* – if both are given *objs* win.
    objs : iterable of :class:`Edi`, optional
        Pre-instantiated objects, already parsed with
        :pymeth:`Edi.read`.  Useful in multiprocessing
        pipelines.
    survey_name : str, optional
        Label used when auto-generating *station ids*.
        Purely cosmetic; leave *None* to rely on file
        names.
    verbose : int or bool, default ``0``
        * 0  → silent  
        * 1  → one-line recap (``quick_edi_stats``)  
        * ≥2 → full banner (``show_edi_stats``)

    Attributes
    ----------
    station_ids : tuple of str
        One id per site.  The default pattern is
        ``S00``, ``S01`` … unless overridden by
        *survey_name* or :pymeth:`rewrite`.
    lat, lon : ndarray of float, shape (n_sites,)
        Geographic coordinates in **decimal degrees**,
        west-to-east / south-to-north order corrected by
        :pyfunc:`pycsamt.gis.utils.scale_position`.
    elev : ndarray of float, shape (n_sites,)
        Site elevations (metres).
    z, z_err : dict[str, ndarray]
        Complex 2 × 2 impedance tensors and 1-σ errors.
    rho, rho_err : dict[str, ndarray]
        Apparent resistivities (Ω·m) and errors.
    phi, phi_err : dict[str, ndarray]
        Phases (°) and errors.

    Notes
    -----
    * **Impedance stacking**   All TF arrays are kept in
      their original frequency sampling.  You are
      responsible for resampling / interpolating before
      stacking across sites.
    * **Memory footprint**   Raw tensors for large MT
      surveys can easily reach hundreds of MB.  Consider
      using :pymeth:`rewrite` to produce trimmed
      light-weight EDI files whenever you just need
      meta-information.
    * **Verbosity**   Reading failures never abort the
      run – they are logged and counted in the summary
      banner so you can inspect them later.

    Methods
    -------
    summary() : str
        One-line recap used in :pycode{print(obj)}.
    rewrite(dst, *, prefix=None, correct_ll=False)
        Bulk-save every site to *dst*, optionally fixing
        coordinates and renaming stations.

    Examples
    --------
    >>> from pycsamt.seg.collection import EdiCollection
    >>> paths = sorted(Path("survey/edi").glob("*.edi"))
    >>> coll = EdiCollection(paths, verbose=1)
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    EDI read  :    12/12    —  success 100.00 %
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    >>> coll.lat.mean(), coll.lon.mean()
    (-22.507, 26.340)
    >>> abs(coll.z['S03']).max()
    14.2
    >>> # export every site under a new prefix
    >>> coll.rewrite(dst=Path("survey/out"), prefix="K1")

    References
    ----------
    * Wight D. E. & Bostick F. X., *SEG-EMAP Data
      Interchange Standard*, 1988.
    * Egbert G. D. (2017). “Impedance tensor
      interpolation”, *Computers & Geosciences*, 103,
      44-56.  https://doi.org/10.1016/j.cageo.2016.12.003
    """

    def __init__(                    
        self,
        paths: Iterable[str | Path] | None = None,
        *,
        objs: Iterable[Edi] | None = None,
        survey_name: str | None = None,
        verbose: int | bool = 0,      # ← NEW ARG
    ) -> None:   
        self.verbose = int(verbose)   # store once for later
        
        if objs is None and paths is None:
            raise FileHandlingError("need either *paths* or *objs*")

        #  Load objects while tracking successes / failures
        t0 = time.perf_counter()
        items: List[Edi] = []
        failures: list[str] = []

        if objs is not None:
            items = list(objs)
        else:
            for p in paths or []:
                try:
                    edi = Edi(path=Path(p)).read()
                    items.append(edi)
                except Exception as err:          # noqa: BLE001
                    logger.warning("could not read %s (%s)", p, err)
                    failures.append(str(p))
 

        super().__init__(items, name=survey_name or "edi")

        # --- compute quick look-ups -----------------------------------
        self._ids: List[str] = make_ids(self, prefix="S")
        self._lat, self._lon, self._elev = self._fetch_geo()

        # impedance, ρ, φ and associated σ
        self._z, self._z_err = self._stack("z", "z_err")
        self._rho, self._rho_err = self._stack(
            "resistivity", "resistivity_err"
        )
        self._phi, self._phi_err = self._stack("phase", "phase_err")


        #  Console feedback (verbose>0)
        if self.verbose:
            tot = len(paths or items) if objs is None else len(objs)
            ok = len(self)
            elapsed = time.perf_counter() - t0

            if self.verbose == 1:
                quick_edi_stats(total=tot, ok=ok, label=self.name.upper())
            else:  # verbose >= 2  → detailed banner
                show_edi_stats(
                    collected=tot,
                    succeeded=ok,
                    failed=len(failures),
                    elapsed=elapsed,
                    title=self.name.upper(),
            )
                
    @property
    def lat(self) -> np.ndarray:
        """Decimal degrees north."""
        return self._lat.copy()

    @property
    def lon(self) -> np.ndarray:
        """Decimal degrees east."""
        return self._lon.copy()

    @property
    def elev(self) -> np.ndarray:
        """Elevation (meters)."""
        return self._elev.copy()


    @property
    def station_ids(self) -> Sequence[str]:
        return tuple(self._ids)

    @property
    def z(self) -> Mapping[str, np.ndarray]:
        """Raw complex 2×2 Z tensor per site."""
        return self._z

    @property
    def z_err(self) -> Mapping[str, np.ndarray]:
        return self._z_err

    @property
    def rho(self) -> Mapping[str, np.ndarray]:
        return self._rho

    @property
    def rho_err(self) -> Mapping[str, np.ndarray]:
        return self._rho_err

    @property
    def phi(self) -> Mapping[str, np.ndarray]:
        return self._phi

    @property
    def phi_err(self) -> Mapping[str, np.ndarray]:
        return self._phi_err

    def summary(self) -> str:  # noqa: D401
        """
        One-line recap used by ``__str__``.
        """
        txt = (
            f"{len(self)} sites • "
            f"lat {self._lat.min():.2f}–{self._lat.max():.2f} °  "
            f"lon {self._lon.min():.2f}–{self._lon.max():.2f} °"
        )
        return txt

    def rewrite(
        self,
        *,
        dst: Path,
        prefix: str | None = None,
        correct_ll: bool = False,
    ) -> None:
        """
        Write every site to *dst* directory.

        Parameters
        ----------
        dst
            Destination folder (will be created).
        prefix
            Optional prefix replacing the default ``"S"`` in
            auto-generated ids.
        correct_ll
            If *True* converts scaled coordinates back to DMS
            strings and rewrites ``>=DEFINEMEAS`` values.
        """
        dst.mkdir(parents=True, exist_ok=True)
        for idx, edi in enumerate(self):
            new_id = (
                f"{prefix or 'S'}{idx:02d}"
                if prefix is not None
                else self._ids[idx]
            )
            edi = edi.copy()  # shallow, keeps Z/Tip arrays
            edi.head.dataid = new_id
            if correct_ll:
                # convert to DD:MM:SS
                lat_dms = _sp_to_dms(self._lat[idx], todms=True)
                lon_dms = _sp_to_dms(self._lon[idx], todms=True)
                edi.head.lat = self._lat[idx]
                edi.head.lon = self._lon[idx]
                dm = edi.definemeasurement
                dm.reflat = lat_dms
                dm.reflong = lon_dms
                dm.refelev = self._elev[idx]
            out = dst / f"{new_id}.edi"
            edi.write(out)
            logger.debug("wrote %s", out)


    def _fetch_geo(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        lat, lon, elev = [], [], []
        for o in self:
            lat.append(o.head.lat)
            lon.append(o.head.lon)
            elev.append(o.head.elev)
        # ensure same order (south-north, west-east)
        lon, *_ = scale_position(lon)
        lat, *_ = scale_position(lat)
        return np.asarray(lat), np.asarray(lon), np.asarray(elev)

    def _stack(
        self, attr: str, err_attr: str
    ) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
        mats, errs = {}, {}
        for sid, edi in zip(self._ids, self):
            mats[sid] = getattr(edi.Z, attr)
            errs[sid] = getattr(edi.Z, err_attr)
        return mats, errs


    def __str__(self) -> str:  # noqa: D401
        base = super().__str__().splitlines()[0]
        return f"{base}  •  {self.summary()}"
