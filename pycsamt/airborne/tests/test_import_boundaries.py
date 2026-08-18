from __future__ import annotations

import subprocess
import sys


def test_import_airborne_does_not_eagerly_import_seg_or_xml_reader():
    code = """
import sys
import pycsamt.airborne
print('pycsamt.seg' in sys.modules)
print('pycsamt.seg.edi' in sys.modules)
print('pycsamt.emtf.xml.reader' in sys.modules)
"""
    out = subprocess.check_output([sys.executable, "-c", code], text=True)
    assert out.strip().splitlines() == ["False", "False", "False"]
