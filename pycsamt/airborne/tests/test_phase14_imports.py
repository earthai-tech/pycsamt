from __future__ import annotations

import subprocess
import sys


def test_common_airborne_import_does_not_load_native_backends():
    code = r'''
import sys
import pycsamt.airborne
print(int("pycsamt.seg.edi" in sys.modules))
print(int("pycsamt.emtf.xml.reader" in sys.modules))
print(int("pycsamt.airborne.mobilemt" in sys.modules))
print(int("pycsamt.airborne.ztem" in sys.modules))
print(int("pycsamt.airborne.afmag" in sys.modules))
'''
    result = subprocess.run(
        [sys.executable, "-c", code],
        check=True,
        capture_output=True,
        text=True,
    )
    assert result.stdout.strip().splitlines() == ["0", "0", "0", "0", "0"]
