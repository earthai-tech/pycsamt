# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Writer for the MARE2DEM ``.settings`` parallel-decomposition file.

Port of ``m2d_writeSettingsFile.m``.

The settings file controls the MPI parallel decomposition: how many
transmitters, receivers, and frequencies are grouped per worker
process. The defaults here match the MATLAB script defaults.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

__all__ = ["SettingsFile", "write_settings"]


@dataclass
class SettingsFile:
    """Parameters for one MARE2DEM ``.settings`` file.

    Attributes
    ----------
    tolerance : float
        Target solution accuracy in percent.
    tx_per_group : int
        Maximum transmitters per parallel group (≤ 10 recommended).
    csem_rx_per_group : int
        CSEM receivers per parallel group.
    csem_freq_per_group : int
        CSEM frequencies per group (usually 1).
    mt_rx_per_group : int
        MT receivers per parallel group.
    mt_freq_per_group : int
        MT frequencies per group (usually 1).
    use_mesh_coarsening : bool
        Enable moving-window mesh coarsening for long profiles.
    use_mt_scattered_field : bool
        Use scattered-field MT formulation (for deep-water scenarios).
    print_adaptive : bool
        Print adaptive refinement iteration stats.
    print_decomposition : bool
        Print parallel decomposition settings.
    """

    tolerance: float = 1.0
    tx_per_group: int = 10
    csem_rx_per_group: int = 40
    csem_freq_per_group: int = 1
    mt_rx_per_group: int = 40
    mt_freq_per_group: int = 1
    use_mesh_coarsening: bool = True
    use_mt_scattered_field: bool = False
    print_adaptive: bool = True
    print_decomposition: bool = True


def write_settings(
    sf: SettingsFile,
    path: str | Path,
    *,
    overwrite: bool = True,
) -> Path:
    """Write a MARE2DEM ``.settings`` file.

    Port of ``m2d_writeSettingsFile.m``.

    Parameters
    ----------
    sf : SettingsFile
        Settings to write.
    path : path-like
        Destination file.
    overwrite : bool, default True
        If ``False`` and the file already exists, return the existing path
        without writing.

    Returns
    -------
    pathlib.Path
        Path of the written (or existing) file.

    Examples
    --------
    >>> from pycsamt.models.mare2dem.iotools.settings import SettingsFile, write_settings
    >>> sf = SettingsFile(tx_per_group=5, csem_rx_per_group=20)
    >>> write_settings(sf, "mare2dem.settings")
    PosixPath('mare2dem.settings')
    """
    dest = Path(path)
    if dest.exists() and not overwrite:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)

    coarsen = "yes" if sf.use_mesh_coarsening else "no"
    scattered = "yes" if sf.use_mt_scattered_field else "no"
    p_adaptive = "yes" if sf.print_adaptive else "no"
    p_decomp = "yes" if sf.print_decomposition else "no"

    with dest.open("w") as fh:
        fh.write(
            f"Tolerance (%):                    {int(sf.tolerance)}"
            "  ! target solution accuracy\n"
        )
        fh.write("\n!\n! CSEM settings:\n!\n")
        fh.write(
            f"Transmitters per group:           {sf.tx_per_group}"
            "     ! set this <= 10\n"
        )
        fh.write(
            f"CSEM receivers per group:         {sf.csem_rx_per_group}"
            "    ! adjust to maximize cluster usage\n"
        )
        fh.write(
            f"CSEM frequencies per group:       {sf.csem_freq_per_group}"
            "     ! usually 1\n"
        )
        fh.write(
            f"Use mesh coarsening:              {coarsen}"
            "   ! moving-window mesh coarsening\n"
        )
        fh.write("\n!\n! MT settings:\n!\n")
        fh.write(
            f"MT receivers per group:           {sf.mt_rx_per_group}"
            "    ! adjust to maximize cluster usage\n"
        )
        fh.write(
            f"MT frequencies per group:         {sf.mt_freq_per_group}"
            "     ! usually 1\n"
        )
        fh.write(
            f"Use MT scattered field:           {scattered}"
            "    ! yes for deep-water resistive lithosphere\n"
        )
        fh.write("\n\n")
        fh.write(
            f"Print adaptive:                   {p_adaptive}"
            "   ! print adaptive refinement stats\n"
        )
        fh.write(
            f"Print decomposition:              {p_decomp}"
            "   ! print parallel decomposition settings\n"
        )
        fh.write("\n!\n! Advice:\n!\n")
        fh.write(
            "! Adjust receivers per group so that the total number of "
            "groups is >= number of MPI processors.\n"
        )

    return dest
