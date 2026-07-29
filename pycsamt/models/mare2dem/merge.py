# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Merge two or more MARE2DEM ``.emdata`` files.

Port of ``m2d_mergeDataFiles.m``.

:func:`merge_data_files` reads multiple ``.emdata`` files and
combines them into one file suitable for joint inversion of MT and
CSEM data (or multiple CSEM profiles).

Merging rules
-------------
* MT and CSEM sections are merged independently.
* Common receivers and transmitters (same position / parameters) are
  deduplicated; indices in the DATA block are remapped accordingly.
* Duplicate receivers can optionally be **kept** for towed-array
  datasets where the same physical offset appears at different absolute
  positions as the transmitter moves.
* UTM origin must be identical across all input files.
* DC resistivity data are not yet supported (raises ``NotImplementedError``).
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import Union

import numpy as np

from .iotools.emdata import (
    CSEMConfig,
    EMDataFile,
    read_emdata,
    write_emdata,
)

__all__ = ["merge_data_files", "merge_emdata"]

PathLike = Union[str, Path]


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _merge_arrays(
    a: np.ndarray,
    b: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Merge rows of *b* into *a*, deduplicating identical rows.

    Returns
    -------
    combined : numpy.ndarray
        Union of *a* and *b* rows (unique rows, *a* order preserved).
    ib : numpy.ndarray
        Mapping: for each row index in *b*, its 1-based index in *combined*.
    """
    if len(a) == 0:
        ib = np.arange(1, len(b) + 1, dtype=int)
        return b.copy(), ib

    combined = list(a)
    ib = np.zeros(len(b), dtype=int)
    a_list = [tuple(row) for row in a]

    for i, row_b in enumerate(b):
        key = tuple(row_b)
        if key in a_list:
            ib[i] = a_list.index(key) + 1
        else:
            combined.append(row_b)
            a_list.append(key)
            ib[i] = len(combined)

    return np.array(combined, dtype=a.dtype), ib


def _merge_lists(
    a: list[str],
    b: list[str],
) -> tuple[list[str], np.ndarray]:
    """Merge string lists, deduplicating; return combined list and b-mapping."""
    combined = list(a)
    ib = np.zeros(len(b), dtype=int)
    for i, name in enumerate(b):
        if name in combined:
            ib[i] = combined.index(name) + 1
        else:
            combined.append(name)
            ib[i] = len(combined)
    return combined, ib


def _check_utm(em1: EMDataFile, em2: EMDataFile) -> None:
    u1, u2 = em1.utm, em2.utm
    if u1.north0 != u2.north0 or u1.east0 != u2.east0 or u1.theta != u2.theta:
        raise ValueError(
            "Cannot merge: UTM origins or 2-D strike angles differ between "
            "the input files.  All files must share the same UTM reference."
        )


def _merge_mt(
    out: EMDataFile,
    src: EMDataFile,
    mt_data: np.ndarray,
) -> EMDataFile:
    """Merge MT section from *src* into *out*."""
    if src.mt is None or not len(src.mt.frequencies):
        return out

    if out.mt is None:
        out.mt = src.mt
        # remap DATA rows: trivial since indices already 1-based from source
        out.data[:, 0].astype(int)
        out.data = np.vstack([out.data, mt_data]) if len(out.data) else mt_data
        return out

    # merge frequencies
    freq_combined, ifreq_map = _merge_arrays(
        out.mt.frequencies.reshape(-1, 1),
        src.mt.frequencies.reshape(-1, 1),
    )
    out.mt.frequencies = freq_combined.ravel()
    mt_data = mt_data.copy()
    mt_data[:, 1] = ifreq_map[mt_data[:, 1].astype(int) - 1]

    # merge receivers
    rx_combined, irx_map = _merge_arrays(out.mt.receivers, src.mt.receivers)
    out.mt.receivers = rx_combined
    mt_data[:, 3] = irx_map[mt_data[:, 3].astype(int) - 1]
    mt_data[:, 2] = irx_map[mt_data[:, 2].astype(int) - 1]  # Tx#==H-rx# for MT

    # merge receiver names
    names_combined, _ = _merge_lists(out.mt.receiver_name, src.mt.receiver_name)
    out.mt.receiver_name = names_combined

    out.data = np.vstack([out.data, mt_data]) if len(out.data) else mt_data
    return out


def _merge_csem(
    out: EMDataFile,
    src: EMDataFile,
    csem_data: np.ndarray,
    keep_duplicate_rx: bool,
) -> EMDataFile:
    """Merge CSEM section from *src* into *out*."""
    if src.csem is None or not len(src.csem.transmitters):
        return out

    if out.csem is None:
        out.csem = src.csem
        out.data = np.vstack([out.data, csem_data]) if len(out.data) else csem_data
        return out

    if out.csem.phase_convention.lower() != src.csem.phase_convention.lower():
        raise ValueError("Cannot merge CSEM sections with different phase conventions.")

    # frequencies
    freq_combined, ifreq_map = _merge_arrays(
        out.csem.frequencies.reshape(-1, 1),
        src.csem.frequencies.reshape(-1, 1),
    )
    out.csem.frequencies = freq_combined.ravel()
    csem_data = csem_data.copy()
    csem_data[:, 1] = ifreq_map[csem_data[:, 1].astype(int) - 1]

    # transmitters (include type as extra boolean column for uniqueness)
    def _tx_with_type(csem_cfg: CSEMConfig) -> np.ndarray:
        txs = csem_cfg.transmitters.copy()
        e_flags = np.array(
            [1 if t.lower() == "edipole" else 0 for t in csem_cfg.transmitter_type],
            dtype=float,
        ).reshape(-1, 1)
        return np.hstack([txs, e_flags])

    tx_a = _tx_with_type(out.csem)
    tx_b = _tx_with_type(src.csem)
    tx_combined, itx_map = _merge_arrays(tx_a, tx_b)

    # restore type list from flag column
    out.csem.transmitters = tx_combined[:, :-1]
    out.csem.transmitter_type = [
        "edipole" if f > 0.5 else "bdipole" for f in tx_combined[:, -1]
    ]
    names_combined_tx, _ = _merge_lists(
        out.csem.transmitter_name, src.csem.transmitter_name
    )
    out.csem.transmitter_name = names_combined_tx
    csem_data[:, 2] = itx_map[csem_data[:, 2].astype(int) - 1]

    # receivers
    if keep_duplicate_rx:
        n_existing = len(out.csem.receivers)
        out.csem.receivers = np.vstack([out.csem.receivers, src.csem.receivers])
        irx_map = np.arange(
            n_existing + 1,
            n_existing + 1 + len(src.csem.receivers),
            dtype=int,
        )
        out.csem.receiver_name = out.csem.receiver_name + src.csem.receiver_name
        csem_data[:, 3] = irx_map[csem_data[:, 3].astype(int) - 1]
    else:
        rx_combined, irx_map = _merge_arrays(out.csem.receivers, src.csem.receivers)
        out.csem.receivers = rx_combined
        names_combined_rx, _ = _merge_lists(
            out.csem.receiver_name, src.csem.receiver_name
        )
        out.csem.receiver_name = names_combined_rx
        csem_data[:, 3] = irx_map[csem_data[:, 3].astype(int) - 1]

    out.data = np.vstack([out.data, csem_data]) if len(out.data) else csem_data
    return out


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def merge_emdata(
    files: list[EMDataFile],
    *,
    keep_duplicate_rx: bool = False,
    comment: str = "",
) -> EMDataFile:
    """Merge a list of :class:`EMDataFile` objects into one.

    Parameters
    ----------
    files : list of EMDataFile
        At least two files to merge.
    keep_duplicate_rx : bool, default False
        When ``True``, identical receiver locations are kept separate
        (needed for towed CSEM arrays with identical offsets).
    comment : str
        Comment written into the merged output file header.

    Returns
    -------
    EMDataFile
        Merged data file (not yet written to disk).

    Raises
    ------
    ValueError
        When fewer than two files are given, or when UTM origins differ.
    """
    if len(files) < 2:
        raise ValueError("Need at least two files to merge.")

    import copy

    out = copy.deepcopy(files[0])
    if not comment:
        comment = (
            f"Merged {len(files)} data files with pycsamt "
            f"on {datetime.now().strftime('%Y-%m-%d %H:%M')}"
        )
    out.comment = comment

    for src in files[1:]:
        if src.data is None or len(src.data) == 0:
            continue

        _check_utm(out, src)

        codes = src.data[:, 0].astype(int)

        if src.dc is not None:
            raise NotImplementedError(
                "DC resistivity data merging is not yet implemented."
            )

        mt_mask = (codes > 100) & (codes < 200)
        csem_mask = codes < 100

        if mt_mask.any():
            out = _merge_mt(out, src, src.data[mt_mask])

        if csem_mask.any():
            out = _merge_csem(out, src, src.data[csem_mask], keep_duplicate_rx)

    return out


def merge_data_files(
    files_to_merge: list[PathLike],
    out_file: PathLike,
    *,
    keep_duplicate_rx: bool = False,
) -> EMDataFile:
    """Merge MARE2DEM ``.emdata`` files and write the result.

    Port of ``m2d_mergeDataFiles.m``.

    Parameters
    ----------
    files_to_merge : list of path-like
        At least two ``.emdata`` files to merge.  Include the full path
        if files are not in the current directory.
    out_file : path-like
        Output ``.emdata`` filename.
    keep_duplicate_rx : bool, default False
        Keep identical receiver locations (needed for towed CSEM arrays).

    Returns
    -------
    EMDataFile
        Merged file contents (also written to *out_file*).

    Examples
    --------
    >>> from pycsamt.models.mare2dem.merge import merge_data_files
    >>> em = merge_data_files(["mt.emdata", "csem.emdata"], "joint.emdata")
    >>> em.n_data
    3200
    """
    if len(files_to_merge) < 2:
        raise ValueError("Need at least two files to merge.")

    ems = [read_emdata(f) for f in files_to_merge]
    merged = merge_emdata(ems, keep_duplicate_rx=keep_duplicate_rx)
    write_emdata(merged, out_file)
    return merged
