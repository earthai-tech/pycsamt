# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Science API evidence is static, bounded and confined to source roots."""
from pycsamt.assistant.tools.api_evidence import inspect_science_api


def test_actual_loader_signature_and_source():
    card = inspect_science_api("pycsamt.emtools._core.ensure_sites")
    assert card["source_path"] == "pycsamt/emtools/_core.py"
    assert "recursive" in card["signature"]
    assert card["line"] > 0


def test_static_inspection_does_not_execute_module(tmp_path):
    module = tmp_path / "pycsamt" / "emtools" / "example.py"
    module.parent.mkdir(parents=True)
    module.write_text("raise RuntimeError('must never execute')\ndef measure(sites, *, units='m'):\n    '''Measurement documentation.'''\n    return sites\n")
    card = inspect_science_api("pycsamt.emtools.example.measure", root=tmp_path)
    assert "units='m'" in card["signature"]
    assert "Measurement" in card["doc"]
    assert inspect_science_api("pycsamt.app.secret.value", root=tmp_path) is None
    assert inspect_science_api("pycsamt.emtools...secret", root=tmp_path) is None
    assert inspect_science_api("pycsamt.emtools.example.missing", root=tmp_path) is None
