# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Reader and writer for MARE2DEM ``.emdata_group`` data-group files.

Port of:
  * ``m2d_readDataGroupFile.m``
  * ``m2d_writeDataGroupFile.m``

Data-group files are optional companions to ``.emdata`` files. They
assign each datum to a named group so MARE2DEM can apply per-group
data weights during joint inversion or when balancing heterogeneous
data volumes (e.g. towed CSEM vs seafloor CSEM vs MT).

Format (``EMDataGroup_1.0``)::

    Format:  EMDataGroup_1.0
    ! optional comment
    # groups:   3
    MT
    Seafloor CSEM
    Towed CSEM
    # data:   1200
    1
    1
    2
    ...
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

__all__ = ["DataGroupFile", "read_data_group", "write_data_group"]


@dataclass
class DataGroupFile:
    """Contents of one MARE2DEM ``.emdata_group`` file.

    Attributes
    ----------
    path : pathlib.Path or None
        Source file path.
    comment : str
        Optional free-text comment from the file header.
    group_names : list of str
        Ordered list of group name strings.
    group_indices : numpy.ndarray, shape (n_data,)
        1-based group index for each datum.  The value at position *i*
        selects the group name at ``group_names[group_indices[i] - 1]``.
    """

    path: Path | None = None
    comment: str = ""
    group_names: list[str] = field(default_factory=list)
    group_indices: np.ndarray = field(
        default_factory=lambda: np.empty(0, dtype=int)
    )

    @property
    def n_groups(self) -> int:
        return len(self.group_names)

    @property
    def n_data(self) -> int:
        return len(self.group_indices)

    def __repr__(self) -> str:
        return f"DataGroupFile(n_groups={self.n_groups}, n_data={self.n_data})"


def read_data_group(path: str | Path) -> DataGroupFile:
    """Read a MARE2DEM ``.emdata_group`` file.

    Port of ``m2d_readDataGroupFile.m``.

    Parameters
    ----------
    path : path-like
        File to read (``Format: EMDataGroup_1.0``).

    Returns
    -------
    DataGroupFile
        Parsed data-group file.

    Raises
    ------
    FileNotFoundError
        When *path* does not exist.
    ValueError
        When the format string is not ``EMDataGroup_1.0``.

    Examples
    --------
    >>> from pycsamt.models.mare2dem.iotools.data_group import read_data_group
    >>> dg = read_data_group("survey.emdata_group")
    >>> dg.group_names
    ['MT', 'Seafloor CSEM', 'Towed CSEM']
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")

    dg = DataGroupFile(path=path.resolve())
    lines = path.read_text(errors="replace").splitlines()

    i = 0
    n_groups = 0

    while i < len(lines):
        raw = lines[i].strip()
        i += 1
        if not raw or raw.startswith("!") or raw.startswith("%"):
            if raw.startswith("!") and not dg.comment:
                dg.comment = raw.lstrip("!")
            continue

        sep = raw.find(":")
        if sep < 0:
            continue
        code = raw[:sep].strip().lower()
        value = raw[sep + 1 :].strip()
        for ch in ("!", "%"):
            ci = value.find(ch)
            if ci >= 0:
                value = value[:ci].strip()

        if code in ("format", "version"):
            if value.lower().strip() not in ("emdatagroup_1.0",):
                raise ValueError(f"Unsupported format: {value!r}")
        elif code == "# groups":
            n_groups = int(value)
            names: list[str] = []
            while len(names) < n_groups and i < len(lines):
                s = lines[i].strip()
                i += 1
                if s and not s.startswith("!") and not s.startswith("%"):
                    names.append(s)
            dg.group_names = names
        elif code in ("# data", "#data"):
            n_data = int(value)
            remaining = "\n".join(lines[i:])
            vals = np.fromstring(
                remaining.replace("\n", " "), sep=" ", dtype=int
            )
            dg.group_indices = vals[:n_data]
            break

    return dg


def write_data_group(dg: DataGroupFile, path: str | Path) -> Path:
    """Write a :class:`DataGroupFile` to *path*.

    Port of ``m2d_writeDataGroupFile.m``.

    Parameters
    ----------
    dg : DataGroupFile
        Data to write.
    path : path-like
        Destination file.

    Returns
    -------
    pathlib.Path
        Path of the written file.

    Raises
    ------
    ValueError
        When :attr:`DataGroupFile.group_names` is empty or
        :attr:`DataGroupFile.group_indices` are out of range.

    Examples
    --------
    >>> from pycsamt.models.mare2dem.iotools.data_group import (
    ...     DataGroupFile,
    ...     write_data_group,
    ... )
    >>> import numpy as np
    >>> dg = DataGroupFile(
    ...     group_names=["MT", "CSEM"], group_indices=np.array([1, 1, 2, 2])
    ... )
    >>> write_data_group(dg, "survey.emdata_group")
    PosixPath('survey.emdata_group')
    """
    if not dg.group_names:
        raise ValueError("group_names must not be empty.")
    if len(dg.group_indices) > 0:
        if dg.group_indices.max() > len(dg.group_names):
            raise ValueError(
                f"group_indices max ({dg.group_indices.max()}) exceeds "
                f"number of group names ({len(dg.group_names)})."
            )
        if dg.group_indices.min() < 1:
            raise ValueError("group_indices must be >= 1 (1-based).")

    dest = Path(path)
    dest.parent.mkdir(parents=True, exist_ok=True)

    with dest.open("w") as fh:
        fh.write("Format:  EMDataGroup_1.0\n")
        if dg.comment:
            fh.write(f"!{dg.comment}\n")
        fh.write(f"# groups:    {len(dg.group_names)}\n")
        for name in dg.group_names:
            fh.write(f"{name}\n")
        fh.write(f"# data:    {len(dg.group_indices)}\n")
        for idx in dg.group_indices:
            fh.write(f"{int(idx)}\n")

    return dest
