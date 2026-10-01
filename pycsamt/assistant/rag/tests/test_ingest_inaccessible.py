# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Unsupported or inaccessible files must not break index freshness."""
import logging
from pathlib import Path

from pycsamt.assistant.rag.ingest import iter_index_files, source_fingerprint


def test_specific_inversion_tags_win_over_generic_ai_path():
    from pycsamt.assistant.rag.config import infer_workflow

    assert infer_workflow("pycsamt/agents/inv2d_agent.py", "AI inversion U-Net") == "inv2d"
    assert infer_workflow("docs/ai_inversion/", "EMInverter2D profile") == "inv2d"
    assert infer_workflow("docs/ai_inversion/", "EMInverter3D model") == "inv3d"
    assert infer_workflow("pycsamt/ai/inversion/", "EMInverter1D") == "ai_inversion"


def test_external_solver_preparation_is_not_neural_inversion():
    from pycsamt.assistant.rag.config import infer_workflow

    assert infer_workflow("export inputs for external 2D inversion") == "pre_inversion"
    assert infer_workflow("prepare third-party inversion input") == "pre_inversion"
    assert infer_workflow("run 2d inversion") == "inv2d"


def test_filters_before_stat_and_reports_unreadable_sources(tmp_path, monkeypatch, caplog):
    package = tmp_path / "pycsamt"
    package.mkdir()
    good = package / "good.py"
    good.write_text("def load(): pass\n")
    html = package / "index.html"
    html.write_text("unsupported vendor link")
    bad = package / "unreadable.py"
    bad.write_text("x = 1")
    original = Path.is_file

    def is_file(path):
        if path == html:
            raise AssertionError("Unsupported HTML must be filtered before stat")
        if path == bad:
            raise OSError(1920, "unreadable reparse point")
        return original(path)

    monkeypatch.setattr(Path, "is_file", is_file)
    # The packaged logging config sets ``propagate: no`` on the ``pycsamt``
    # logger, so once another test has loaded it the warning never reaches
    # caplog's root handler -- attach the handler to the logger directly.
    ingest_logger = logging.getLogger("pycsamt.assistant.rag.ingest")
    ingest_logger.addHandler(caplog.handler)
    monkeypatch.setattr(ingest_logger, "disabled", False)
    caplog.set_level(logging.WARNING, logger=ingest_logger.name)
    try:
        assert list(iter_index_files(tmp_path)) == [good]
    finally:
        ingest_logger.removeHandler(caplog.handler)
    assert list(iter_index_files(tmp_path)) == [good]
    assert source_fingerprint(tmp_path) == source_fingerprint(tmp_path, files=[good])
    assert "unreadable.py" in caplog.text
