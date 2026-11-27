# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Write round-trip regressions: MTSECT NFREQ and INFO free text."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np

from pycsamt.seg.edi import EDIFile
from pycsamt.z.z import Z


def _nfreq(text: str) -> int:
    sect = text.split(">=MTSECT", 1)[1]
    return int(re.search(r"NFREQ\s*=\s*(\d+)", sect).group(1))


def test_mtsect_nfreq_follows_trimmed_z(edi_imp_file: Path,
                                       tmp_path: Path) -> None:
    ed = EDIFile(edi_imp_file, verbose=0)
    f = np.asarray(ed.Z.freq)
    keep = slice(0, max(2, f.size // 2))
    ed.Z = Z(z_array=ed.Z.z[keep], z_err_array=ed.Z.z_err[keep],
             freq=f[keep])
    ed.write(new_edifn="trim.edi", savepath=str(tmp_path))
    text = (tmp_path / "trim.edi").read_text(encoding="utf-8")
    assert _nfreq(text) == f[keep].size
    back = EDIFile(tmp_path / "trim.edi", verbose=0)
    assert np.asarray(back.Z.freq).size == f[keep].size


def test_info_free_text_with_equals_survives(edi_imp_file: Path,
                                            tmp_path: Path) -> None:
    src = edi_imp_file.read_text(encoding="utf-8", errors="replace")
    note = ["    Station coordinates: CGCS2000 GK CM 111E",
            "    N = 2830301.633 m, E = 793059.577 m",
            "    UTM_zone = 49"]
    src = src.replace(">INFO", ">INFO\n" + "\n".join(note), 1)
    p = tmp_path / "info.edi"
    p.write_text(src, encoding="utf-8")
    ed = EDIFile(p, verbose=0)
    ed.write(new_edifn="out.edi", savepath=str(tmp_path))
    text = (tmp_path / "out.edi").read_text(encoding="utf-8")
    for line in note:
        assert line.strip() in text
    # known INFO keys are still parsed as fields
    info = ed.get_section("info")
    assert hasattr(info, "maxinfo")
