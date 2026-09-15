# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Additional coverage-focused tests for :mod:`pycsamt.agents.mare2dem_agent`.

Complements ``test_mare2dem_agent.py`` (which already covers the
``emdata``/``mt``/``csem`` prepare paths, run-mode-without-a-binary, and
report mode against real bundled sample data) with: the EDI-pathway
(``sites``/``path``) branch, the ``models.mare2dem`` import guard, the
``download_source`` auto-build branches (success / build-failed /
download-exception), the user-supplied ``resistivity`` override, the
``read_emdata`` exception guard, the ``InputBuilder`` exception guard, and
the report-mode ``InversionResult`` exception guard.

The MARE2DEM binary itself is never available in this environment (see
``models/_solver_build``), so ``SourceManager`` is mocked for the
auto-download branches rather than performing a real download/compile.
"""

from __future__ import annotations

import sys
import types

import pytest

from pycsamt.agents import Mare2DEMAgent

pytestmark = pytest.mark.filterwarnings("ignore::RuntimeWarning")


def test_models_mare2dem_import_error_fails(tmp_path, monkeypatch):
    fake_mod = types.ModuleType("pycsamt.models.mare2dem")
    monkeypatch.setitem(sys.modules, "pycsamt.models.mare2dem", fake_mod)
    agent = Mare2DEMAgent()
    result = agent.execute(
        {"emdata": "x.emdata", "output_dir": str(tmp_path)}
    )
    assert result.status == "failed"
    assert "pycsamt.models.mare2dem not available" in result.error


@pytest.fixture()
def bh_edi_dir(project_root):
    # The bundled 3edis toy dataset has too few/degenerate station
    # coordinates for MARE2DEM's real profile-geometry fit (SVD does not
    # converge); the real Broken Hill survey has proper along-profile
    # station spacing.
    d = project_root / "data" / "MT" / "broken-hill" / "edis"
    if not d.exists():
        pytest.skip(f"Broken Hill EDI dataset not found: {d}")
    return d


def test_edi_pathway_builds_emdata(bh_edi_dir, tmp_path):
    agent = Mare2DEMAgent()
    result = agent.execute(
        {"path": str(bh_edi_dir), "output_dir": str(tmp_path)}
    )
    assert result.status in ("success", "needs_review")
    assert result["data_path"] is not None
    assert result["data_path"].exists()
    assert result["n_mt_receivers"] > 0


def test_edi_pathway_confidence_weighting(bh_edi_dir, tmp_path):
    agent = Mare2DEMAgent()
    result = agent.execute(
        {
            "path": str(bh_edi_dir),
            "output_dir": str(tmp_path),
            "confidence_weighting": "true",
        }
    )
    assert result["data_path"].exists()


def test_auto_download_build_succeeds(tmp_path, monkeypatch):
    import pycsamt.models.mare2dem as m2d_pkg

    state = {"built": False}
    monkeypatch.setattr(
        m2d_pkg.SourceManager, "is_built", lambda self: state["built"]
    )
    monkeypatch.setattr(
        m2d_pkg.SourceManager, "download", lambda self: None
    )

    def _build(self):
        state["built"] = True

    monkeypatch.setattr(m2d_pkg.SourceManager, "build", _build)
    agent = Mare2DEMAgent()
    result = agent.execute(
        {
            "mt": {"frequencies": [1.0], "rx_y": [0.0], "lTE": True},
            "output_dir": str(tmp_path),
            "download_source": True,
        }
    )
    assert result["binary_found"] is True
    assert result["source_downloaded"] is True


def test_auto_download_build_reports_failure(tmp_path, monkeypatch):
    import pycsamt.models.mare2dem as m2d_pkg

    monkeypatch.setattr(m2d_pkg.SourceManager, "is_built", lambda self: False)
    monkeypatch.setattr(
        m2d_pkg.SourceManager, "download", lambda self: None
    )
    monkeypatch.setattr(m2d_pkg.SourceManager, "build", lambda self: None)
    agent = Mare2DEMAgent()
    result = agent.execute(
        {
            "mt": {"frequencies": [1.0], "rx_y": [0.0], "lTE": True},
            "output_dir": str(tmp_path),
            "download_source": True,
        }
    )
    assert any(
        "build may have failed" in w for w in result.warnings
    )
    assert result["binary_found"] is False


def test_auto_download_exception_is_recorded(tmp_path, monkeypatch):
    import pycsamt.models.mare2dem as m2d_pkg

    monkeypatch.setattr(m2d_pkg.SourceManager, "is_built", lambda self: False)
    monkeypatch.setattr(
        m2d_pkg.SourceManager,
        "download",
        lambda self: (_ for _ in ()).throw(RuntimeError("network down")),
    )
    agent = Mare2DEMAgent()
    result = agent.execute(
        {
            "mt": {"frequencies": [1.0], "rx_y": [0.0], "lTE": True},
            "output_dir": str(tmp_path),
            "download_source": True,
        }
    )
    assert any("Auto-build failed" in w for w in result.warnings)


def test_resistivity_override_is_copied(tmp_path):
    resist_src = tmp_path / "custom.resistivity"
    resist_src.write_text("fake resistivity content")
    out_dir = tmp_path / "run"
    agent = Mare2DEMAgent()
    result = agent.execute(
        {
            "mt": {"frequencies": [1.0], "rx_y": [0.0], "lTE": True},
            "resistivity": str(resist_src),
            "output_dir": str(out_dir),
        }
    )
    assert result["resistivity_path"].read_text() == "fake resistivity content"


def test_read_emdata_exception_is_swallowed(tmp_path, monkeypatch):
    import pycsamt.models.mare2dem as m2d_pkg

    monkeypatch.setattr(
        m2d_pkg,
        "read_emdata",
        lambda path: (_ for _ in ()).throw(RuntimeError("parse boom")),
    )
    agent = Mare2DEMAgent()
    result = agent.execute(
        {
            "mt": {"frequencies": [1.0], "rx_y": [0.0], "lTE": True},
            "output_dir": str(tmp_path),
        }
    )
    # statistics fall back to their zero defaults; no crash
    assert result["n_mt_receivers"] == 0
    assert result["n_data"] == 0


def test_input_builder_exception_fails(tmp_path, monkeypatch):
    import pycsamt.models.mare2dem as m2d_pkg

    monkeypatch.setattr(
        m2d_pkg.InputBuilder,
        "build",
        lambda self, *a, **k: (_ for _ in ()).throw(
            RuntimeError("builder boom")
        ),
    )
    agent = Mare2DEMAgent()
    result = agent.execute(
        {
            "mt": {"frequencies": [1.0], "rx_y": [0.0], "lTE": True},
            "output_dir": str(tmp_path),
        }
    )
    assert result.status == "failed"
    assert "InputBuilder failed" in result.error


def test_report_mode_inversion_result_exception_fails(tmp_path, monkeypatch):
    import pycsamt.models.mare2dem as m2d_pkg

    monkeypatch.setattr(
        m2d_pkg,
        "InversionResult",
        lambda output_dir, config: (_ for _ in ()).throw(
            RuntimeError("scan boom")
        ),
    )
    agent = Mare2DEMAgent()
    result = agent.execute(
        {"output_dir": str(tmp_path), "mode": "report"}
    )
    assert result.status == "failed"
    assert "Could not scan run directory" in result.error
