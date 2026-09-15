# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Coverage-focused tests for :mod:`pycsamt.agents.denoising`.

Targets: ``ensure_sites`` failure, the unknown-method fallback, each
denoising method dispatch branch (rpca / hampel / emap / pipeline / ai),
the SNR-gain exception guard, the comparison-figure exception guard, the
LLM interpretation path, and the private ``_apply_ai_denoiser`` helper's
no-data / import-error / generic-exception / happy-path branches.
"""

from __future__ import annotations

import sys
import types

import numpy as np
import pytest

from pycsamt.agents.denoising import DenoisingAgent, _apply_ai_denoiser

pytestmark = pytest.mark.filterwarnings("ignore::RuntimeWarning")


# ── fake site helpers (bypass ensure_sites for direct-helper tests) ─────────


class _FakeZ:
    def __init__(self, z, freq):
        self.z = z
        self.freq = freq


class _FakeSite:
    def __init__(self, z=None, freq=None, station="st"):
        self.station = station
        if z is not None:
            self.Z = _FakeZ(z, freq)


def _zblock(n=8, seed=0):
    rng = np.random.default_rng(seed)
    z = rng.normal(size=(n, 2, 2)) + 1j * rng.normal(size=(n, 2, 2))
    freq = np.linspace(1.0, 100.0, n)
    return z, freq


# ── DenoisingAgent.execute() tests ───────────────────────────────────────────


def test_no_sites_or_path(no_llm_kw):
    agent = DenoisingAgent(**no_llm_kw)
    result = agent.execute({})
    assert result.status == "failed"
    assert "No 'sites' or 'path'" in result.error


def test_ensure_sites_raises(no_llm_kw, monkeypatch):
    import pycsamt.emtools._core as core

    monkeypatch.setattr(
        core,
        "ensure_sites",
        lambda *a, **k: (_ for _ in ()).throw(ValueError("bad")),
    )
    agent = DenoisingAgent(**no_llm_kw)
    result = agent.execute({"path": "/nope"})
    assert result.status == "failed"
    assert "bad" in result.error


def test_unknown_method_falls_back_to_rpca(no_llm_kw, loaded_sites, monkeypatch):
    import pycsamt.emtools.remove_noise as rn

    monkeypatch.setattr(rn, "rpca_offdiag_denoise", lambda sites, **k: sites)
    agent = DenoisingAgent(**no_llm_kw, method="not_a_method")
    result = agent.execute({"sites": loaded_sites})
    assert result.status == "success"
    assert any("Unknown method" in w for w in result.warnings)
    assert result["method"] == "rpca"


def test_method_rpca(no_llm_kw, loaded_sites, monkeypatch):
    import pycsamt.emtools.remove_noise as rn

    monkeypatch.setattr(rn, "rpca_offdiag_denoise", lambda sites, **k: sites)
    agent = DenoisingAgent(**no_llm_kw, method="rpca")
    result = agent.execute({"sites": loaded_sites})
    assert result.status == "success"
    assert result["method"] == "rpca"


def test_method_hampel(no_llm_kw, loaded_sites, monkeypatch):
    import pycsamt.emtools.remove_noise as rn

    monkeypatch.setattr(rn, "hampel_filter_freq", lambda sites, **k: sites)
    agent = DenoisingAgent(**no_llm_kw, method="hampel")
    result = agent.execute({"sites": loaded_sites})
    assert result.status == "success"
    assert result["method"] == "hampel"


def test_method_emap_success(no_llm_kw, loaded_sites, monkeypatch):
    import pycsamt.emtools.remove_noise as rn

    monkeypatch.setattr(rn, "apply_emap_filter", lambda sites, **k: sites)
    agent = DenoisingAgent(**no_llm_kw, method="emap")
    result = agent.execute({"sites": loaded_sites})
    assert result.status == "success"


def test_method_emap_failure_is_recorded(no_llm_kw, loaded_sites, monkeypatch):
    import pycsamt.emtools.remove_noise as rn

    monkeypatch.setattr(
        rn,
        "apply_emap_filter",
        lambda *a, **k: (_ for _ in ()).throw(RuntimeError("emap boom")),
    )
    agent = DenoisingAgent(**no_llm_kw, method="emap")
    result = agent.execute({"sites": loaded_sites})
    assert result.status == "success"
    assert any("EMAP filter failed" in w for w in result.warnings)


def test_method_pipeline_success(no_llm_kw, loaded_sites, monkeypatch):
    import pycsamt.emtools.remove_noise as rn

    monkeypatch.setattr(rn, "remove_noise_pipeline", lambda sites, **k: sites)
    agent = DenoisingAgent(**no_llm_kw, method="pipeline")
    result = agent.execute({"sites": loaded_sites})
    assert result.status == "success"


def test_method_pipeline_failure_is_recorded(
    no_llm_kw, loaded_sites, monkeypatch
):
    import pycsamt.emtools.remove_noise as rn

    monkeypatch.setattr(
        rn,
        "remove_noise_pipeline",
        lambda *a, **k: (_ for _ in ()).throw(RuntimeError("pipe boom")),
    )
    agent = DenoisingAgent(**no_llm_kw, method="pipeline")
    result = agent.execute({"sites": loaded_sites})
    assert result.status == "success"
    assert any("Noise pipeline failed" in w for w in result.warnings)


def test_method_ai_dispatch(no_llm_kw, loaded_sites, monkeypatch):
    import pycsamt.agents.denoising as den

    monkeypatch.setattr(
        den, "_apply_ai_denoiser", lambda sites, warnings: sites
    )
    agent = DenoisingAgent(**no_llm_kw, method="ai")
    result = agent.execute({"sites": loaded_sites})
    assert result.status == "success"


def test_snr_gain_exception_is_swallowed(no_llm_kw, loaded_sites, monkeypatch):
    import pycsamt.agents.denoising as den
    import pycsamt.emtools.remove_noise as rn

    monkeypatch.setattr(rn, "rpca_offdiag_denoise", lambda sites, **k: sites)
    # a plain list can't be boolean-indexed by np.isfinite(...) -> TypeError,
    # caught by the surrounding `except Exception: pass`.
    monkeypatch.setattr(den, "_compute_snr_proxy", lambda sites: [1.0, 2.0])
    agent = DenoisingAgent(**no_llm_kw, method="rpca")
    result = agent.execute({"sites": loaded_sites})
    assert result.status == "success"
    assert result["snr_gain"] == 0.0


def test_comparison_figure_exception_is_captured(
    no_llm_kw, loaded_sites, monkeypatch
):
    import pycsamt.emtools.inspect as insp
    import pycsamt.emtools.remove_noise as rn

    monkeypatch.setattr(rn, "rpca_offdiag_denoise", lambda sites, **k: sites)
    monkeypatch.setattr(
        insp,
        "pseudosection",
        lambda *a, **k: (_ for _ in ()).throw(RuntimeError("plot boom")),
    )
    agent = DenoisingAgent(**no_llm_kw, method="rpca")
    result = agent.execute({"sites": loaded_sites})
    assert any("Comparison figure" in w for w in result.warnings)
    assert result["figures"] == {}


def test_ama_happy_path_saves_figure(no_llm_kw, loaded_sites, tmp_output, monkeypatch):
    import pycsamt.emtools.remove_noise as rn

    monkeypatch.setattr(rn, "rpca_offdiag_denoise", lambda sites, **k: sites)
    agent = DenoisingAgent(**no_llm_kw, method="rpca")
    result = agent.execute(
        {"sites": loaded_sites, "output_dir": str(tmp_output)}
    )
    assert result.status == "success"
    assert result["figures"]
    assert result["figure_paths"]


def test_llm_interpretation_called_when_api_key_set(loaded_sites, monkeypatch):
    import pycsamt.emtools.remove_noise as rn

    monkeypatch.setattr(rn, "rpca_offdiag_denoise", lambda sites, **k: sites)
    agent = DenoisingAgent(api_key="fake-key", method="rpca")
    agent.query_llm = lambda *a, **k: "mocked interpretation"
    result = agent.execute({"sites": loaded_sites})
    assert result.llm_interpretation == "mocked interpretation"


# ── _apply_ai_denoiser direct tests ──────────────────────────────────────────


def test_apply_ai_denoiser_no_valid_z_returns_sites_unchanged():
    warnings: list[str] = []
    sites: list = []
    result = _apply_ai_denoiser(sites, warnings)
    assert result is sites
    assert any("No valid Z data" in w for w in warnings)


def test_apply_ai_denoiser_import_error_falls_back(monkeypatch):
    fake_mod = types.ModuleType("pycsamt.ai.processing.denoise")
    monkeypatch.setitem(
        sys.modules, "pycsamt.ai.processing.denoise", fake_mod
    )
    z, fr = _zblock()
    sites = [_FakeSite(z=z, freq=fr, station="s1")]
    warnings: list[str] = []
    result = _apply_ai_denoiser(sites, warnings)
    assert result is sites
    assert any("requires PyTorch or TensorFlow" in w for w in warnings)


def test_apply_ai_denoiser_happy_path(monkeypatch):
    class _FakeDenoiser:
        def __init__(self, n_freqs):
            self.n_freqs = n_freqs

        def fit(self, X, epochs=20, verbose=False):
            return None

        def predict(self, X):
            return X

    fake_mod = types.ModuleType("pycsamt.ai.processing.denoise")
    fake_mod.EMDenoiser = _FakeDenoiser
    fake_mod.prepare_z_features = lambda arr: arr
    monkeypatch.setitem(
        sys.modules, "pycsamt.ai.processing.denoise", fake_mod
    )
    z, fr = _zblock()
    sites = [_FakeSite(z=z, freq=fr, station="s1")]
    warnings: list[str] = []
    result = _apply_ai_denoiser(sites, warnings)
    assert result is sites
    assert any("AI denoiser applied" in w for w in warnings)


def test_apply_ai_denoiser_generic_exception_falls_back(monkeypatch):
    class _BoomDenoiser:
        def __init__(self, n_freqs):
            raise RuntimeError("gpu unavailable")

    fake_mod = types.ModuleType("pycsamt.ai.processing.denoise")
    fake_mod.EMDenoiser = _BoomDenoiser
    fake_mod.prepare_z_features = lambda arr: arr
    monkeypatch.setitem(
        sys.modules, "pycsamt.ai.processing.denoise", fake_mod
    )
    z, fr = _zblock()
    sites = [_FakeSite(z=z, freq=fr, station="s1")]
    warnings: list[str] = []
    result = _apply_ai_denoiser(sites, warnings)
    assert result is sites
    assert any("AI denoiser failed" in w for w in warnings)
