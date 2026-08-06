# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""OccamModel — read and write the Occam2DModel file.

The model file links the mesh to the inversion parameterisation.  It
specifies which mesh cells share the same free parameter (the
``Num Layers`` / column-mapping block) and global model metadata.

Entry points
------------
``OccamModel.from_mesh(mesh, config)``
    Build a model definition consistent with a given mesh.
``OccamModel.read(path)``
    Parse an existing ``Occam2DModel`` file.
``OccamModel.write(path)``
    Serialise to disk in OCCAM2MTMOD_1.0 format.

File format (OCCAM2MTMOD_1.0)
------------------------------
::

    Format:           OCCAM2MTMOD_1.0
    Model Name:       <name>
    Description:      <description>
    Mesh File:        <mesh_filename>
    Mesh Type:        PW2D
    Statics File:     none
    Prejudice File:   none
    Binding Offset:   <float>
    Num Layers:       <N>
      <n_merge>  <n_cols>
      <n_cols integers — parameter codes for each model column>
      ... (repeated N times)
    NO. EXCEPTIONS:   <K>

Parameter codes
---------------
Each integer in a layer's data row encodes the cell type:

* **7** — boundary column (tied to the binding value, not a free parameter)
* **even ≥ 2** — active column; the exact value encodes how many mesh columns
  are merged into one model column (2 = 1:1, 4 = 2:1, …)

The total cell count across all layers equals the ``Param Count`` recorded in
the Startup / .iter file.
"""

from __future__ import annotations

from pathlib import Path
from typing import Union

import numpy as np

from .base import OccamBase
from .config import OccamConfig

PathLike = Union[str, Path]

__all__ = ["OccamModel"]

_FORMAT_TAG = "OCCAM2MTMOD_1.0"

# Header keyword → attribute name (upper-cased key for matching)
_HEADER_MAP: dict[str, str] = {
    "FORMAT": "format_str",
    "MODEL NAME": "name",
    "DESCRIPTION": "description",
    "MESH FILE": "mesh_file",
    "MESH TYPE": "mesh_type",
    "STATICS FILE": "statics_file",
    "PREJUDICE FILE": "prejudice_file",
    "BINDING OFFSET": "binding_offset",
    "NUM LAYERS": "n_layers",
    "NO. EXCEPTIONS": "n_exceptions",
}

_FLOAT_FIELDS = {"binding_offset"}
_INT_FIELDS = {"n_layers", "n_exceptions"}


# -----------------------------------------------------------------------
# Low-level parser
# -----------------------------------------------------------------------


def _parse_model(path: Path) -> dict:
    """Parse an OCCAM2MTMOD_1.0 file.

    Returns
    -------
    dict
        Keys are all ``_HEADER_MAP`` values plus ``"layers"``.
        Each layer has ``n_merge``, ``n_cols``, and
        ``params``.

    Raises
    ------
    FileNotFoundError
    ValueError
        If the Format tag is absent or wrong.
    """
    if not path.exists():
        raise FileNotFoundError(f"Occam2DModel file not found: {path}")

    result: dict = {v: None for v in _HEADER_MAP.values()}
    result["layers"] = []

    with path.open("r", errors="replace") as fh:
        lines = [line.rstrip("\n") for line in fh]

    i = 0
    N = len(lines)

    while i < N:
        raw = lines[i].strip()
        i += 1

        if not raw:
            continue

        # Header line?
        if ":" in raw:
            raw_key, _, raw_val = raw.partition(":")
            key = raw_key.strip().upper()
            val = raw_val.strip()
            attr = _HEADER_MAP.get(key)
            if attr is not None:
                if attr in _INT_FIELDS:
                    try:
                        result[attr] = int(float(val))
                    except ValueError:
                        result[attr] = 0
                elif attr in _FLOAT_FIELDS:
                    try:
                        result[attr] = float(val)
                    except ValueError:
                        result[attr] = 0.0
                elif attr == "format_str":
                    result[attr] = val
                    if val.upper() != _FORMAT_TAG:
                        raise ValueError(
                            f"Expected format '{_FORMAT_TAG}', got '{val}' in {path}"
                        )
                else:
                    result[attr] = val

                # After Num Layers we switch to reading layer blocks
                if attr == "n_layers" and result["n_layers"]:
                    n_layers = result["n_layers"]
                    for _ in range(n_layers):
                        # Control line: n_merge  n_cols
                        while i < N and not lines[i].strip():
                            i += 1
                        ctrl = lines[i].strip().split()
                        i += 1
                        n_merge = int(ctrl[0])
                        n_cols = int(ctrl[1])
                        # Data line: n_cols integers
                        while i < N and not lines[i].strip():
                            i += 1
                        params = np.array(
                            [int(x) for x in lines[i].strip().split()],
                            dtype=np.int32,
                        )
                        i += 1
                        result["layers"].append(
                            {
                                "n_merge": n_merge,
                                "n_cols": n_cols,
                                "params": params,
                            }
                        )

    if result["format_str"] is None:
        raise ValueError(
            f"File does not contain a valid OCCAM2MTMOD_1.0 header: {path}"
        )

    return result


# -----------------------------------------------------------------------
# OccamModel
# -----------------------------------------------------------------------


class OccamModel(OccamBase):
    r"""Represent the Occam2D model-parameter definition.

    ``OccamModel`` links a finite-element ``OccamMesh`` to the
    inversion parameter vector used by Occam2D. The model file
    does not store resistivity values. Instead, it defines how
    mesh cells are grouped into free or fixed parameters.
    The startup and iteration files then store one value for
    each parameter counted by :attr:`n_params`.

    Each model layer contains integer column codes. Boundary
    code ``7`` marks fixed edge columns tied to the binding
    value, while active even codes represent free inversion
    columns [OccamModel-1]_. If :math:`p_j` is the code for one model
    column, then the number of mesh cells represented by that
    column is encoded by the code value itself. The PyCSAMT
    builder uses ``2`` for interior columns and ``7`` for the
    two boundary columns:

    .. math::

        \mathbf{p}
        =
        [7,\;2,\;2,\;\ldots,\;2,\;7].

    Parameters
    ----------
    name : str, default "MODEL MADE BY PYCSAMT"
        Model name written to the ``Model Name`` header field.
        Use this for a short label that identifies the model
        family, processing run, or inversion setup.
    description : str, default "SMOOTH INVERSION"
        Description written to the ``Description`` header
        field. It is intended for human-readable provenance
        and is preserved when the model is written to disk.
    config : OccamConfig, optional
        Configuration object used for default file names and
        related Occam2D settings. If omitted, a default
        :class:`OccamConfig` is created.
    verbose : int or bool, default 0
        Verbosity level inherited from :class:`OccamBase`.
        Positive values enable progress messages through the
        instance logger.
    logger : logging.Logger, optional
        Logger used for progress and diagnostic messages. If
        omitted, a class-specific PyCSAMT logger is created.

    Attributes
    ----------
    format_str : str
        Occam model format tag. The writer uses
        ``"OCCAM2MTMOD_1.0"``.
    name : str
        Model name written to the file header.
    description : str
        Human-readable model description.
    config : OccamConfig
        Configuration used by this object.
    mesh_file : str
        Filename of the associated mesh used by the model
        file, usually ``"Occam2DMesh"``.
    mesh_type : str
        Mesh type string written to the header. Occam2D PW2D
        meshes use ``"PW2D"``.
    statics_file : str
        Optional static-shift file. The default ``"none"``
        means no static correction file is referenced.
    prejudice_file : str
        Optional prejudice model file. The default ``"none"``
        means no prejudice file is referenced.
    binding_offset : float, default 0.0
        Horizontal offset of the binding column. This is the
        reference value used by boundary columns.
    n_layers : int
        Number of active model layers after the header.
        It is normally the number of non-air mesh rows.
    layers : list of dict
        Per-layer parameter specification. Each entry has the
        following keys:

        ``n_merge`` : int
            Number of mesh z-rows merged into this layer.
        ``n_cols`` : int
            Number of model columns in this layer.
        ``params`` : numpy.ndarray of int, shape (n_cols,)
            Parameter codes for each model column. Code ``7``
            marks a boundary column; active even values mark
            free inversion parameters.

    n_exceptions : int
        Number of exception records at the end of the model
        file. PyCSAMT currently writes ``0``.

    Notes
    -----
    ``n_params`` is the sum of ``n_cols`` over all model
    layers. This value must match the ``Param Count`` in the
    startup and iteration files. ``n_free_params`` excludes
    boundary columns with code ``7``.

    See Also
    --------
    OccamMesh
        Defines finite-element cells grouped by this model.
    OccamStartup.from_model
        Creates an initial vector with matching size.
    InputBuilder
        Builds data, mesh, model, and startup files together.

    Examples
    --------
    Build a model definition from an existing mesh:

    >>> from pycsamt.models.occam2d import OccamMesh
    >>> from pycsamt.models.occam2d import OccamModel
    >>> mesh = OccamMesh.read("occam_run/Occam2DMesh")
    >>> model = OccamModel.from_mesh(mesh)
    >>> model.n_params

    Read and write an existing model file:

    >>> from pycsamt.models.occam2d import OccamModel
    >>> model = OccamModel.read("occam_run/Occam2DModel")
    >>> model.write("copy/Occam2DModel")

    Create a custom empty container for tests:

    >>> from pycsamt.models.occam2d import OccamModel
    >>> model = OccamModel(name="SYNTHETIC MODEL")
    >>> model.n_layers, model.n_params

    References
    ----------
    .. [OccamModel-1] deGroot-Hedlin, C., and Constable, S.,
       "Occam's inversion to generate smooth, two-dimensional
       models from magnetotelluric data", Geophysics, 55(12),
       1613-1624, 1990.
    .. [OccamModel-2] Constable, S. C., Parker, R. L., and Constable,
       C. G., "Occam's inversion: A practical algorithm for
       generating smooth models from electromagnetic sounding
       data", Geophysics, 52(3), 289-300, 1987.
    """

    def __init__(
        self,
        name: str = "MODEL MADE BY PYCSAMT",
        description: str = "SMOOTH INVERSION",
        config: OccamConfig | None = None,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.format_str: str = _FORMAT_TAG
        self.name: str = name
        self.description: str = description
        self.config: OccamConfig = config or OccamConfig()
        self.mesh_file: str = self.config.mesh_file
        self.mesh_type: str = "PW2D"
        self.statics_file: str = "none"
        self.prejudice_file: str = "none"
        self.binding_offset: float = 0.0
        self.n_layers: int = 0
        self.layers: list[dict] = []
        self.n_exceptions: int = 0

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------
    @classmethod
    def from_mesh(
        cls,
        mesh: OccamMesh,  # noqa: F821
        config: OccamConfig | None = None,
        **kwargs,
    ) -> OccamModel:
        r"""Build a model definition from a populated mesh.

        The method converts a finite-element mesh into the
        column mapping required by ``OCCAM2MTMOD_1.0``. Air
        rows are ignored. Each remaining earth row becomes one
        model layer with ``n_merge = 1``. Horizontally, the
        model uses fixed boundary columns and active interior
        columns:

        .. math::

            \mathbf{p}
            =
            [7,\;2,\;2,\;\ldots,\;2,\;7].

        The leading and trailing ``7`` codes represent seven
        mesh cells each at the profile boundaries. Interior
        ``2`` codes represent two mesh cells per free model
        column. Therefore the total horizontal cell count is

        .. math::

            n_x = 7 + 2n_i + 7,

        where :math:`n_i` is the number of interior model
        columns. Meshes built by :meth:`OccamMesh.from_data`
        are constructed to satisfy this layout.

        Parameters
        ----------
        mesh : OccamMesh
            Populated mesh object defining horizontal and
            vertical finite-element cells. It must provide
            ``n_xcells``, ``n_zcells``, and ``n_airlayers``.
            Air layers are excluded from the model; all other
            z-cells become inversion layers.
        config : OccamConfig, optional
            Configuration object used for file names and
            related Occam2D defaults. If omitted, a default
            :class:`OccamConfig` is created.
        **kwargs
            Additional keyword arguments forwarded to the
            ``OccamModel`` constructor. This is commonly used
            for ``name``, ``description``, ``verbose``, or
            ``logger``.

        Returns
        -------
        OccamModel
            Model-definition object ready to be written as an
            ``Occam2DModel`` file. The returned object has
            ``n_layers`` equal to the number of active earth
            rows and ``layers`` populated with parameter-code
            arrays.

        Raises
        ------
        ValueError
            Raised when the mesh has no active earth layers or
            fewer than fourteen horizontal cells. Fourteen
            cells are required for the two seven-cell boundary
            columns.

        See Also
        --------
        OccamMesh.from_data
            Builds meshes that match this parameterization.
        OccamModel.write
            Serializes the returned model definition.
        OccamStartup.from_model
            Creates a startup vector with matching parameter
            count.

        Examples
        --------
        Build a model from a mesh read from disk:

        >>> from pycsamt.models.occam2d import OccamMesh
        >>> from pycsamt.models.occam2d import OccamModel
        >>> mesh = OccamMesh.read("occam_run/Occam2DMesh")
        >>> model = OccamModel.from_mesh(mesh)
        >>> model.n_layers

        Pass metadata through to the model header:

        >>> model = OccamModel.from_mesh(
        ...     mesh,
        ...     name="PROFILE A MODEL",
        ...     description="smooth TE-TM inversion",
        ... )
        """
        cfg = config or OccamConfig()
        obj = cls(config=cfg, **kwargs)
        n_xcells = mesh.n_xcells
        n_air = mesh.n_airlayers
        n_active = mesh.n_zcells - n_air

        if n_active <= 0:
            raise ValueError(
                f"OccamModel.from_mesh: no active layers "
                f"(n_zcells={mesh.n_zcells}, n_airlayers={n_air})"
            )
        if n_xcells < 14:
            raise ValueError(
                f"OccamModel.from_mesh: n_xcells={n_xcells} < 14; "
                "mesh too narrow for boundary columns"
            )

        # Build per-layer code array: [7, 2, 2, ..., 2, 7].
        # Generated meshes normally have an even interior width, but accepting
        # an odd-width existing mesh makes preparation robust when users bring
        # in older or hand-edited Occam2DMesh files.
        interior_width = n_xcells - 14
        interior = [2] * (interior_width // 2)
        if interior_width % 2:
            if interior:
                interior[-1] += 1
            else:
                interior = [1]
        codes = np.array([7] + interior + [7], dtype=np.int32)

        assert int(codes.sum()) == n_xcells, (
            f"code sum {codes.sum()} != n_xcells {n_xcells}"
        )
        n_cols = len(codes)

        layers: list[dict] = []
        for _ in range(n_active):
            layers.append(
                {
                    "n_merge": 1,
                    "n_cols": n_cols,
                    "params": codes.copy(),
                }
            )

        obj.mesh_file = cfg.mesh_file
        obj.n_layers = n_active
        obj.layers = layers
        obj.n_exceptions = 0

        if obj.verbose:
            obj.logger.info(
                "OccamModel.from_mesh: %d layers, %d cols/layer, %d total params",
                n_active,
                n_cols,
                obj.n_params,
            )
        return obj

    # ------------------------------------------------------------------
    # I/O
    # ------------------------------------------------------------------
    @classmethod
    def read(cls, path: PathLike, **kwargs) -> OccamModel:
        """Read an existing ``OCCAM2MTMOD_1.0`` model file.

        The reader parses the model header and the per-layer
        parameter-code blocks. Numeric header values are cast
        to ``int`` or ``float`` where appropriate. Layer
        ``params`` arrays are stored as ``numpy.int32`` for
        comparison with generated model mappings.

        Parameters
        ----------
        path : path-like
            Path to an Occam2D model file. The value may be a
            string, :class:`pathlib.Path`, or any object
            accepted by :class:`pathlib.Path`.
        **kwargs
            Additional keyword arguments forwarded to the
            ``OccamModel`` constructor before parsed values
            are attached. Use this for ``config``,
            ``verbose``, or ``logger``.

        Returns
        -------
        OccamModel
            Parsed model-definition container with header
            fields and layer mappings populated from ``path``.

        Raises
        ------
        FileNotFoundError
            Raised when ``path`` does not exist.
        ValueError
            Raised when the format tag is missing or is not
            ``"OCCAM2MTMOD_1.0"``.

        See Also
        --------
        OccamModel.write
            Writes model definitions in the same format.
        OccamStartup.read
            Reads startup or iteration files that depend on
            the same parameter count.

        Examples
        --------
        >>> from pycsamt.models.occam2d import OccamModel
        >>> model = OccamModel.read("occam_run/Occam2DModel")
        >>> model.n_layers, model.n_params
        """
        p = Path(path)
        parsed = _parse_model(p)
        obj = cls(**kwargs)

        for attr in _HEADER_MAP.values():
            val = parsed.get(attr)
            if val is not None and hasattr(obj, attr):
                setattr(obj, attr, val)

        obj.layers = parsed["layers"]

        if obj.verbose:
            obj.logger.info(
                "OccamModel.read: %d layers, %d total cells from %s",
                obj.n_layers,
                obj.n_params,
                p,
            )
        return obj

    def write(self, path: PathLike) -> Path:
        """Write this model in ``OCCAM2MTMOD_1.0`` format.

        The writer serializes the current header fields and
        layer mappings to the Occam2D model layout. Parent
        directories are created before writing. The object is
        not modified, so the same instance can be written to
        multiple run directories.

        Parameters
        ----------
        path : path-like
            Destination path for the model file. The value may
            be a string, :class:`pathlib.Path`, or any object
            accepted by :class:`pathlib.Path`.

        Returns
        -------
        pathlib.Path
            Path to the file that was written.

        See Also
        --------
        OccamModel.read
            Parses model files written by this method.
        InputBuilder.build
            Calls this method during input generation.

        Examples
        --------
        >>> from pycsamt.models.occam2d import OccamModel
        >>> model = OccamModel.read("source/Occam2DModel")
        >>> written = model.write("copy/Occam2DModel")
        """
        p = Path(path)
        p.parent.mkdir(parents=True, exist_ok=True)

        _W = 18  # width for the keyword field including colon

        def _kv(kw: str, val) -> str:
            key = f"{kw}:"
            return f"{key:<{_W}}{val}\n"

        lines: list[str] = [
            _kv("Format", self.format_str),
            _kv("Model Name", self.name),
            _kv("Description", self.description),
            _kv("Mesh File", self.mesh_file),
            _kv("Mesh Type", self.mesh_type),
            _kv("Statics File", self.statics_file),
            _kv("Prejudice File", self.prejudice_file),
            _kv("Binding Offset", f"{self.binding_offset:.1f}"),
            _kv("Num Layers", self.n_layers),
        ]

        for layer in self.layers:
            n_merge = layer["n_merge"]
            n_cols = layer["n_cols"]
            params = layer["params"]
            lines.append(f"     {n_merge}   {n_cols}\n")
            lines.append("    " + "    ".join(str(v) for v in params) + "\n")

        lines.append(f"NO. EXCEPTIONS:   {self.n_exceptions}\n")

        with p.open("w") as fh:
            fh.writelines(lines)
        self.path = p
        return p

    # ------------------------------------------------------------------
    # Derived
    # ------------------------------------------------------------------
    @property
    def n_params(self) -> int:
        """Total model cells, equal to iter ``Param Count``."""
        return sum(layer["n_cols"] for layer in self.layers)

    @property
    def n_free_params(self) -> int:
        """Number of free (non-boundary) model cells."""
        return sum(int((layer["params"] != 7).sum()) for layer in self.layers)
