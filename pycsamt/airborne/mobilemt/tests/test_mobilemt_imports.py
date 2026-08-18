from __future__ import annotations

import json
import subprocess
import sys


def test_mobilemt_import_does_not_load_historical_edi_backend():
    code = r'''
import json
import sys
import pycsamt.airborne.mobilemt
print(json.dumps({
    "seg": "pycsamt.seg" in sys.modules,
    "edi": "pycsamt.seg.edi" in sys.modules,
}))
'''
    result = subprocess.run(
        [sys.executable, "-c", code],
        check=True,
        capture_output=True,
        text=True,
    )
    state = json.loads(result.stdout.strip().splitlines()[-1])
    assert state == {"seg": False, "edi": False}
