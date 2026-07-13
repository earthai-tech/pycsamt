# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Full reader and writer for MARE2DEM ``.resistivity`` files.

Port of:
  * ``m2d_readResistivity.m``
  * ``m2d_writeResistivity.m``

The ``.resistivity`` file is the primary input file passed to the
MARE2DEM binary.  It contains all inversion control parameters, file
references, regularization settings, and the region-by-region
resistivity table that defines the 2-D FEM model.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path

import numpy as np

__all__ = [
    "ResistivityFile",
    "read_resistivity",
    "write_resistivity",
]

# ---------------------------------------------------------------------------
# Anisotropy modes and their column count
# ---------------------------------------------------------------------------

_ANISO_NCOLS: dict[str, int] = {
    "isotropic": 1,
    "isotropic_ip": 4,  # Cole-Cole
    "isotropic_complex": 2,
    "triaxial": 3,
    "tix": 2,
    "tiy": 2,
    "tiz": 2,
    "tiz_ratio": 2,
}


def _nrho(anisotropy: str) -> int:
    return _ANISO_NCOLS.get(anisotropy.lower().strip(), 1)


# ---------------------------------------------------------------------------
# Data class
# ---------------------------------------------------------------------------


@dataclass
class ResistivityFile:
    """Contents of one MARE2DEM ``.resistivity`` file.

    Attributes
    ----------
    resistivity_file : str
        Output filename (used when writing).
    poly_file : str
        Triangle mesh (``.poly``) file stem referenced by this model.
    data_file : str
        ``.emdata`` file referenced by this model.
    settings_file : str
        ``.settings`` file referenced by this model.
    version : str
        Format version string (e.g. ``"mare2dem_1.1"``).
    anisotropy : str
        Anisotropy mode: ``"isotropic"`` (default), ``"triaxial"``,
        ``"tix"``, ``"tiy"``, ``"tiz"``, ``"tiz_ratio"``,
        ``"isotropic_ip"``, ``"isotropic_complex"``.
    target_misfit : float
        Target normalized RMS misfit.
    max_iterations : int
        Maximum inversion iterations.
    iteration : int
        Current iteration number (output from inversion).
    log10_lagrange : float
        Log10 Lagrange (regularization) trade-off value.
    roughness : float or None
        Model roughness at this iteration.
    misfit : float or None
        Normalized misfit at this iteration.
    date_time : str
        Date/time stamp written by MARE2DEM.
    bounds_transform : str
        Bounds transform type (``"bandpass"``).
    global_bounds : numpy.ndarray, shape (2,)
        Lower and upper log10-resistivity bounds.
    roughness_penalty_method : str
        Smoothing type (``"gradient"`` or ``"first_difference"``).
    yz_penalty_weights : numpy.ndarray, shape (2,)
        Smoothing weights (y, z).
    penalty_cut_weight : float
        Penalty cut weight.
    roughness_with_prejudice : bool
        Use prejudice in regularization.
    beta_mgs : float
        Minimum gradient support weight.
    anisotropy_penalty_weight : float or None
        Anisotropic penalty weight.
    anisotropy_ratio_roughness_weight : float or None
        Anisotropic ratio roughness weight.
    debug_level : int
        Verbosity / print level.
    inversion_method : str
        Inversion method string.
    rms_threshold : float
        Misfit decrease threshold.
    converge_slowly : str
        ``"yes"`` or ``"no"``.
    resistivity : numpy.ndarray, shape (n_regions, nrho)
        Log10-resistivity per region per component.
    free_parameter : numpy.ndarray, shape (n_regions, nrho)
        Free-parameter index (0 = fixed, >0 = free parameter number).
    bounds : numpy.ndarray, shape (n_regions, 2*nrho)
        Lower / upper bounds per region component.
    prejudice : numpy.ndarray, shape (n_regions, 2*nrho)
        Prejudice value and weight per region component.
    """

    resistivity_file: str = "mare2dem.resistivity"
    poly_file: str = "mare2dem.poly"
    data_file: str = "mare2dem.emdata"
    settings_file: str = "mare2dem.settings"
    version: str = "mare2dem_1.1"
    anisotropy: str = "isotropic"

    target_misfit: float = 1.0
    max_iterations: int = 100
    iteration: int = 0
    log10_lagrange: float = 5.0
    roughness: float | None = None
    misfit: float | None = None
    date_time: str = ""

    bounds_transform: str = "bandpass"
    global_bounds: np.ndarray = field(
        default_factory=lambda: np.array([-2.0, 5.0], dtype=float)
    )
    roughness_penalty_method: str = "gradient"
    yz_penalty_weights: np.ndarray = field(
        default_factory=lambda: np.array([1.0, 1.0], dtype=float)
    )
    penalty_cut_weight: float = 0.1
    roughness_with_prejudice: bool = False
    beta_mgs: float = 0.0
    anisotropy_penalty_weight: float | None = None
    anisotropy_ratio_roughness_weight: float | None = None
    debug_level: int = 1
    inversion_method: str = "occam"
    rms_threshold: float = 0.85
    converge_slowly: str = "no"

    # Region tables
    resistivity: np.ndarray = field(
        default_factory=lambda: np.empty((0, 1), dtype=float)
    )
    free_parameter: np.ndarray = field(
        default_factory=lambda: np.empty((0, 1), dtype=float)
    )
    bounds: np.ndarray = field(
        default_factory=lambda: np.empty((0, 2), dtype=float)
    )
    prejudice: np.ndarray = field(
        default_factory=lambda: np.empty((0, 2), dtype=float)
    )

    # ---- misc optional fields ----
    data_group_file: str = ""
    joint_inv_weight_type: str = ""
    penalty_file: str = ""
    fixed_mu_cut: float | None = None

    @property
    def num_regions(self) -> int:
        """Number of resistivity regions."""
        return len(self.resistivity)

    def __repr__(self) -> str:
        return (
            f"ResistivityFile(num_regions={self.num_regions}, "
            f"anisotropy={self.anisotropy!r}, "
            f"target_misfit={self.target_misfit})"
        )


# ---------------------------------------------------------------------------
# Reader
# ---------------------------------------------------------------------------


def _strtok(line: str) -> tuple[str, str]:
    idx = line.find(":")
    if idx < 0:
        return line.strip().lower(), ""
    key = line[:idx].strip().lower()
    val = line[idx + 1 :].strip()
    # strip trailing comment
    for ch in ("!", "%"):
        ci = val.find(ch)
        if ci >= 0:
            val = val[:ci].strip()
    return key, val


def read_resistivity(
    path: str | Path,
    *,
    no_data: bool = False,
) -> ResistivityFile:
    """Read a MARE2DEM ``.resistivity`` file.

    Port of ``m2d_readResistivity.m``.

    Parameters
    ----------
    path : path-like
        File to read.
    no_data : bool, default False
        When ``True``, stop reading before the region resistivity
        table. Useful for quickly inspecting header metadata in large
        iteration output files.

    Returns
    -------
    ResistivityFile
        Parsed file contents.

    Examples
    --------
    >>> from pycsamt.models.mare2dem.iotools.resistivity import read_resistivity
    >>> rf = read_resistivity("mare2dem.resistivity")
    >>> rf.num_regions
    1234
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")

    rf = ResistivityFile()
    rf.resistivity_file = path.name

    lines = path.read_text(errors="replace").splitlines()
    i = 0
    while i < len(lines):
        raw = lines[i]
        i += 1
        code, value = _strtok(raw)

        if not code or code.startswith("!") or code.startswith("%"):
            continue

        if code in ("format", "version"):
            rf.version = value
        elif code in ("model file", "poly file"):
            rf.poly_file = value.strip()
        elif code == "data file":
            rf.data_file = value.strip()
        elif code == "data group file":
            rf.data_group_file = value.strip().lower()
        elif code == "joint inversion weight":
            rf.joint_inv_weight_type = value.strip().lower()
        elif code == "settings file":
            rf.settings_file = value.strip()
        elif code == "penalty file":
            rf.penalty_file = value.strip()
        elif code in ("maximum interations", "maximum iterations"):
            try:
                rf.max_iterations = int(value)
            except ValueError:
                pass
        elif code == "bounds transform":
            rf.bounds_transform = value.strip()
        elif code == "global bounds":
            value = value.replace(",", " ")
            parts = value.split()
            if len(parts) >= 2:
                rf.global_bounds = np.array(
                    [float(parts[0]), float(parts[1])]
                )
        elif code == "roughness penalty method":
            rf.roughness_penalty_method = value.strip()
        elif code == "roughness weights (y,z)":
            value = value.replace(",", " ")
            parts = value.split()
            if len(parts) >= 2:
                rf.yz_penalty_weights = np.array(
                    [float(parts[0]), float(parts[1])]
                )
        elif code == "penalty cut weight":
            try:
                rf.penalty_cut_weight = float(value)
            except ValueError:
                pass
        elif code == "roughness with prejudice":
            rf.roughness_with_prejudice = value.strip().lower() == "yes"
        elif code == "min. gradient support weight":
            try:
                rf.beta_mgs = float(value)
            except ValueError:
                pass
        elif code == "aniso. penalty weight":
            try:
                rf.anisotropy_penalty_weight = float(value)
            except ValueError:
                pass
        elif code == "aniso. ratio roughness weight":
            try:
                rf.anisotropy_ratio_roughness_weight = float(value)
            except ValueError:
                pass
        elif code in ("debug level", "print level"):
            try:
                rf.debug_level = int(value)
            except ValueError:
                pass
        elif code == "target misfit":
            try:
                rf.target_misfit = float(value)
            except ValueError:
                pass
        elif code == "iteration":
            try:
                rf.iteration = int(value)
            except ValueError:
                pass
        elif code in ("lagrange value", "log10 lagrange value"):
            try:
                rf.log10_lagrange = float(value)
            except ValueError:
                pass
        elif code == "model roughness":
            try:
                rf.roughness = float(value)
            except ValueError:
                pass
        elif code == "model misfit":
            try:
                rf.misfit = float(value)
            except ValueError:
                pass
        elif code == "date/time":
            rf.date_time = value.strip()
        elif code == "inversion method":
            rf.inversion_method = value.strip().lower()
        elif code == "fixed mu cut":
            try:
                rf.fixed_mu_cut = float(value)
            except ValueError:
                pass
        elif code == "converge slowly":
            rf.converge_slowly = value.strip().lower()
        elif code == "misfit decrease threshold":
            try:
                rf.rms_threshold = float(value)
            except ValueError:
                pass
        elif code == "anisotropy":
            rf.anisotropy = value.strip()
        elif code == "number of regions":
            if no_data:
                break
            n_regions = int(value)
            nrho = _nrho(rf.anisotropy)
            # expected columns per row: index + rho*nrho + freeparam*nrho + bounds*(2*nrho) + prej*(2*nrho)
            # skip comment line
            if i < len(lines):
                i += 1  # comment header line

            # Read remaining numeric data fast
            remaining_text = "\n".join(lines[i:])
            vals = np.fromstring(
                remaining_text.replace("\n", " "), sep=" ", dtype=float
            )
            cols_per_row = (
                1 + nrho + nrho + 2 * nrho + 2 * nrho
            )  # idx + rho + param + bounds + prej
            expected = n_regions * cols_per_row
            if len(vals) >= expected:
                vals = vals[:expected].reshape(n_regions, cols_per_row)
                c = 1  # skip index column
                rf.resistivity = vals[:, c : c + nrho]
                c += nrho
                rf.free_parameter = vals[:, c : c + nrho]
                c += nrho
                rf.bounds = vals[:, c : c + 2 * nrho]
                c += 2 * nrho
                rf.prejudice = vals[:, c : c + 2 * nrho]
            else:
                # fallback slow read
                rf.resistivity = np.zeros((n_regions, nrho))
                rf.free_parameter = np.zeros((n_regions, nrho))
                rf.bounds = np.zeros((n_regions, 2 * nrho))
                rf.prejudice = np.zeros((n_regions, 2 * nrho))
                n_read = 0
                j = i
                while n_read < n_regions and j < len(lines):
                    line = lines[j].strip()
                    j += 1
                    if not line or line.startswith("!"):
                        continue
                    row = [float(v) for v in line.split()]
                    if len(row) < cols_per_row:
                        continue
                    c = 1
                    rf.resistivity[n_read, :] = row[c : c + nrho]
                    c += nrho
                    rf.free_parameter[n_read, :] = row[c : c + nrho]
                    c += nrho
                    rf.bounds[n_read, :] = row[c : c + 2 * nrho]
                    c += 2 * nrho
                    rf.prejudice[n_read, :] = row[c : c + 2 * nrho]
                    n_read += 1
            break

    return rf


# ---------------------------------------------------------------------------
# Writer
# ---------------------------------------------------------------------------


def write_resistivity(
    rf: ResistivityFile, path: str | Path | None = None
) -> Path:
    """Write a :class:`ResistivityFile` to *path*.

    Port of ``m2d_writeResistivity.m``.

    Parameters
    ----------
    rf : ResistivityFile
        Model to write.
    path : path-like or None
        Destination file. When ``None``, uses ``rf.resistivity_file``.

    Returns
    -------
    pathlib.Path
        Path of the written file.

    Examples
    --------
    >>> from pycsamt.models.mare2dem.iotools.resistivity import write_resistivity
    >>> write_resistivity(rf, "mare2dem_iter10.resistivity")
    PosixPath('mare2dem_iter10.resistivity')
    """
    dest = Path(path or rf.resistivity_file)
    if dest.suffix.lower() != ".resistivity":
        dest = dest.with_suffix(".resistivity")
    dest.parent.mkdir(parents=True, exist_ok=True)

    aniso = rf.anisotropy.lower().strip()
    nrho = _nrho(aniso)

    # column header parts
    _rho_labels = {
        "isotropic": ["Rho"],
        "isotropic_ip": ["Rho", "Eta", "Tau", "C"],
        "isotropic_complex": ["Rho Real", "Rho Imag"],
        "triaxial": ["Rho-x", "Rho-y", "Rho-z"],
        "tix": ["Rho-x", "Rho-yz"],
        "tiy": ["Rho-y", "Rho-xz"],
        "tiz": ["Rho-z", "Rho-h"],
        "tiz_ratio": ["Rho-z", "Rho z/h"],
    }
    rho_labels = _rho_labels.get(aniso, ["Rho"])

    date_str = rf.date_time or datetime.now().strftime("%a %b %d %H:%M:%S %Y")
    bounds_str = (
        f"{rf.global_bounds[0]:.10g}, {rf.global_bounds[1]:.10g}"
        if rf.global_bounds is not None and len(rf.global_bounds) >= 2
        else "-2, 5"
    )
    weights_str = (
        f"{rf.yz_penalty_weights[0]:.10g}, {rf.yz_penalty_weights[1]:.10g}"
        if rf.yz_penalty_weights is not None
        and len(rf.yz_penalty_weights) >= 2
        else "1.0, 1.0"
    )

    with dest.open("w") as fh:
        fh.write(
            f"Format:                         {'mare2dem_1.1':<32s} ! input \n"
        )
        fh.write(
            f"Model File:                     {rf.poly_file:<32s} ! input \n"
        )
        fh.write(
            f"Data File:                      {rf.data_file:<32s} ! input \n"
        )
        if rf.data_group_file:
            fh.write(
                f"Data Group File:                {rf.data_group_file:<32s} ! opt. input \n"
            )
        if rf.joint_inv_weight_type:
            fh.write(
                f"Joint inversion weight:         {rf.joint_inv_weight_type:<32s} ! opt. input\n"
            )
        fh.write(
            f"Settings File:                  {rf.settings_file:<32s} ! input \n"
        )
        fh.write(
            f"Maximum Iterations:             {rf.max_iterations:<32d} ! opt. input \n"
        )
        fh.write(
            f"Bounds Transform:               {'bandpass':<32s} ! opt. input \n"
        )
        fh.write(
            f"Global Bounds:                  {bounds_str:<32s} ! opt. input \n"
        )
        fh.write(
            f"Roughness Penalty Method:       {rf.roughness_penalty_method:<32s} ! opt. input \n"
        )
        fh.write(
            f"Roughness Weights (y,z):        {weights_str:<32s} ! opt. input \n"
        )
        fh.write(
            f"Penalty Cut Weight:             {rf.penalty_cut_weight:<32g} ! opt. input \n"
        )
        prej_str = "yes" if rf.roughness_with_prejudice else "no"
        fh.write(
            f"Roughness With Prejudice:       {prej_str:<32s} ! opt. input \n"
        )
        fh.write(
            f"Min. Gradient Support Weight:   {rf.beta_mgs:<32g} ! opt. input \n"
        )
        if aniso not in ("isotropic", "isotropic_ip", "isotropic_complex"):
            if rf.anisotropy_penalty_weight is not None:
                fh.write(
                    f"Aniso. Penalty Weight:          {rf.anisotropy_penalty_weight:<32g} ! opt. input \n"
                )
            if rf.anisotropy_ratio_roughness_weight is not None:
                fh.write(
                    f"Aniso. Ratio Roughness Weight:  {rf.anisotropy_ratio_roughness_weight:<32g} ! opt. input \n"
                )
        fh.write(
            f"Print Level:                    {rf.debug_level:<32d} ! opt. input  \n"
        )
        fh.write(
            f"Target Misfit:                  {rf.target_misfit:<32g} ! require for inversion\n"
        )
        fh.write(
            f"Misfit Decrease Threshold:      {rf.rms_threshold:<32g} ! opt. input\n"
        )
        fh.write(
            f"Converge Slowly:                {rf.converge_slowly:<32s} ! opt. input\n"
        )
        fh.write(
            f"Log10 Lagrange Value:           {rf.log10_lagrange:<32g} ! input/output\n"
        )
        roughness_str = (
            f"{rf.roughness:.10g}" if rf.roughness is not None else " "
        )
        misfit_str = f"{rf.misfit:.10g}" if rf.misfit is not None else " "
        fh.write(
            f"Model Roughness:                {roughness_str:<32s} ! output from inversion\n"
        )
        fh.write(
            f"Model Misfit:                   {misfit_str:<32s} ! output from inversion\n"
        )
        fh.write(
            f"Date/Time:                      {date_str:<32s} ! output from inversion\n"
        )
        fh.write(
            f"Anisotropy:                     {rf.anisotropy:<32s} ! input \n"
        )

        n_regions = rf.num_regions
        fh.write(
            f"Number of regions:              {n_regions:<32d} ! input \n"
        )

        # column header
        rho_hdr = " ".join(f"{'%-13s' % lb}" for lb in rho_labels)
        param_hdr = " ".join(
            f"{'%-10s' % ('Param ' + lb.split('-')[-1][:3])}"
            for lb in rho_labels
        )
        lower_upper = " ".join(
            f"{'%-13s' % ('Lower' + lb.split('-')[-1][:3])} {'%-13s' % ('Upper' + lb.split('-')[-1][:3])}"
            for lb in rho_labels
        )
        prej_hdr = " ".join(
            f"{'%-13s' % 'Prej'} {'%-13s' % 'Weight'}" for _ in rho_labels
        )
        fh.write(
            f"{'!#':<8s} {rho_hdr} {param_hdr} {lower_upper} {prej_hdr}\n"
        )

        # update free_parameter numbering
        if rf.free_parameter is not None and len(rf.free_parameter):
            fp = rf.free_parameter.copy()
            n_free = 0
            for ii in range(n_regions):
                row_free = fp[ii] > 0
                if np.any(row_free):
                    counts = np.cumsum(row_free)
                    fp[ii] = np.where(row_free, n_free + counts, 0)
                    n_free = int(fp[ii].max())
        else:
            fp = np.zeros((n_regions, nrho), dtype=float)

        for ii in range(n_regions):
            rho_str = " ".join(f"{v:<13.7g}" for v in rf.resistivity[ii])
            fp_str = " ".join(f"{int(v):<10d}" for v in fp[ii])
            bnds = (
                rf.bounds[ii]
                if rf.bounds is not None and len(rf.bounds) > ii
                else np.zeros(2 * nrho)
            )
            bnds_str = " ".join(f"{v:<13.7g}" for v in bnds)
            prej = (
                rf.prejudice[ii]
                if rf.prejudice is not None and len(rf.prejudice) > ii
                else np.zeros(2 * nrho)
            )
            prej_str = " ".join(f"{v:<13.7g}" for v in prej)
            fh.write(
                f"{ii + 1:<9d} {rho_str} {fp_str} {bnds_str} {prej_str}\n"
            )

    return dest
