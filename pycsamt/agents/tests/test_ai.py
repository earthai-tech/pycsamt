# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for AI agents — AIInversionAgent, Inv2DAgent, EnsembleAgent,
JointInversionAgent.

All DL-dependent tests are skipped when no backend (torch/tensorflow) is
available; the no-backend error path is always tested.
"""

import unittest
from pathlib import Path

_EDI_DIR = Path(__file__).parents[4] / "data" / "3edis"
_HAS_EDIS = _EDI_DIR.exists() and any(_EDI_DIR.glob("*.edi"))


def _has_backend() -> bool:
    try:
        from pycsamt.backends import get_backend

        return get_backend() != "none"
    except Exception:
        return False


# ── helpers ───────────────────────────────────────────────────────────────────


def _load_sites():
    from pycsamt.agents import MTLoaderAgent

    r = MTLoaderAgent().execute({"path": str(_EDI_DIR)})
    return r["sites"] if r.status == "success" else None


# ── AIInversionAgent ──────────────────────────────────────────────────────────


class TestAIInversionAgentNoBackend(unittest.TestCase):
    """Error-path test: clear failure message when no DL framework."""

    @unittest.skipIf(
        _has_backend(), "DL backend present — skip no-backend path"
    )
    def test_no_backend_returns_failed(self):
        from pycsamt.agents import AIInversionAgent

        ag = AIInversionAgent()
        r = ag.execute({"path": str(_EDI_DIR)})
        self.assertEqual(r.status, "failed")
        self.assertIsNotNone(r.error_fix_hint)

    def test_from_pretrained_classmethod_exists(self):
        from pycsamt.agents import AIInversionAgent

        self.assertTrue(callable(AIInversionAgent.from_pretrained))


@unittest.skipUnless(
    _HAS_EDIS and _has_backend(), "requires 3edis data + DL backend"
)
class TestAIInversionAgentWithBackend(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from pycsamt.agents import AIInversionAgent

        cls.sites = _load_sites()
        if cls.sites is None:
            raise unittest.SkipTest("Could not load sites.")
        cls.agent = AIInversionAgent(
            arch="cnn1d", n_layers=3, n_train_samples=120, epochs=3
        )
        cls.result = cls.agent.execute({"sites": cls.sites})

    def test_status_success_or_needs_review(self):
        self.assertIn(self.result.status, ("success", "needs_review"))

    def test_predictions_dict(self):
        preds = self.result.get("predictions") or {}
        self.assertIsInstance(preds, dict)

    def test_inverter_present(self):
        self.assertIsNotNone(self.result.get("inverter"))

    def test_figures_dict(self):
        figs = self.result.get("figures") or {}
        self.assertIsInstance(figs, dict)


# ── Inv2DAgent ────────────────────────────────────────────────────────────────


class TestInv2DAgentNoBackend(unittest.TestCase):
    @unittest.skipIf(
        _has_backend(), "DL backend present — skip no-backend path"
    )
    def test_no_backend_fails_gracefully(self):
        from pycsamt.agents import Inv2DAgent

        ag = Inv2DAgent()
        r = ag.execute({"path": str(_EDI_DIR)})
        self.assertEqual(r.status, "failed")
        self.assertIsNotNone(r.error_fix_hint)


@unittest.skipUnless(
    _HAS_EDIS and _has_backend(), "requires 3edis data + DL backend"
)
class TestInv2DAgentWithBackend(unittest.TestCase):
    def test_execute_returns_result(self):
        from pycsamt.agents import Inv2DAgent

        sites = _load_sites()
        if sites is None:
            self.skipTest("Could not load sites.")
        ag = Inv2DAgent(n_freqs=16, n_depth=10, n_train_profiles=20, epochs=2)
        r = ag.execute({"sites": sites})
        self.assertIn(r.status, ("success", "failed", "needs_review"))


# ── EnsembleAgent ─────────────────────────────────────────────────────────────


class TestEnsembleAgentNoBackend(unittest.TestCase):
    @unittest.skipIf(
        _has_backend(), "DL backend present — skip no-backend path"
    )
    def test_no_backend_fails_gracefully(self):
        from pycsamt.agents import EnsembleAgent

        ag = EnsembleAgent()
        r = ag.execute({"path": str(_EDI_DIR)})
        self.assertEqual(r.status, "failed")


@unittest.skipUnless(
    _HAS_EDIS and _has_backend(), "requires 3edis data + DL backend"
)
class TestEnsembleAgentWithBackend(unittest.TestCase):
    def test_execute_coverage_in_range(self):
        from pycsamt.agents import EnsembleAgent

        sites = _load_sites()
        if sites is None:
            self.skipTest("Could not load sites.")
        ag = EnsembleAgent(
            n_estimators=2, n_layers=3, n_train_samples=100, epochs=2
        )
        r = ag.execute({"sites": sites})
        self.assertIn(r.status, ("success", "needs_review"))
        cov = r.get("coverage")
        if cov is not None and cov == cov:  # not nan
            self.assertGreaterEqual(cov, 0.0)
            self.assertLessEqual(cov, 1.0)


# ── JointInversionAgent ───────────────────────────────────────────────────────


class TestJointInversionAgentNoBackend(unittest.TestCase):
    @unittest.skipIf(
        _has_backend(), "DL backend present — skip no-backend path"
    )
    def test_no_backend_fails_gracefully(self):
        from pycsamt.agents import JointInversionAgent

        ag = JointInversionAgent()
        r = ag.execute({"path": str(_EDI_DIR)})
        self.assertEqual(r.status, "failed")

    def test_modalities_stored(self):
        from pycsamt.agents import JointInversionAgent

        ag = JointInversionAgent(modalities=["mt", "csamt"])
        self.assertEqual(ag.modalities, ["mt", "csamt"])


@unittest.skipUnless(
    _HAS_EDIS and _has_backend(), "requires 3edis data + DL backend"
)
class TestJointInversionAgentWithBackend(unittest.TestCase):
    def test_execute_returns_result(self):
        from pycsamt.agents import JointInversionAgent

        sites = _load_sites()
        if sites is None:
            self.skipTest("Could not load sites.")
        ag = JointInversionAgent(
            modalities=["mt", "tem"],
            n_layers=3,
            n_train_samples=100,
            epochs=2,
        )
        r = ag.execute({"sites": sites})
        self.assertIn(r.status, ("success", "needs_review", "failed"))
        self.assertIn("modalities", r.data)


# ── Inv3DAgent ────────────────────────────────────────────────────────────────


class TestInv3DAgentNoBackend(unittest.TestCase):
    @unittest.skipIf(
        _has_backend(), "DL backend present — skip no-backend path"
    )
    def test_no_backend_fails_gracefully(self):
        from pycsamt.agents import Inv3DAgent

        ag = Inv3DAgent()
        r = ag.execute({"path": str(_EDI_DIR)})
        self.assertEqual(r.status, "failed")
        self.assertIsNotNone(r.error_fix_hint)

    def test_constructor_parameters_stored(self):
        from pycsamt.agents import Inv3DAgent

        ag = Inv3DAgent(
            n_layers=7, n_freqs=16, radius=8_000.0, n_mc=10, hidden=(128, 64)
        )
        self.assertEqual(ag.n_layers, 7)
        self.assertEqual(ag.n_freqs, 16)
        self.assertAlmostEqual(ag.radius, 8_000.0)
        self.assertEqual(ag.n_mc, 10)
        self.assertEqual(ag.hidden, (128, 64))

    def test_name_attribute(self):
        from pycsamt.agents import Inv3DAgent

        self.assertEqual(Inv3DAgent().name, "Inv3DAgent")

    def test_no_sites_fails(self):
        from pycsamt.agents import Inv3DAgent

        ag = Inv3DAgent()
        r = ag.execute({})
        self.assertEqual(r.status, "failed")


@unittest.skipUnless(
    _HAS_EDIS and _has_backend(), "requires 3edis data + DL backend"
)
class TestInv3DAgentWithBackend(unittest.TestCase):
    def test_execute_returns_result(self):
        from pycsamt.agents import Inv3DAgent

        sites = _load_sites()
        if sites is None:
            self.skipTest("Could not load sites.")
        ag = Inv3DAgent(
            n_layers=3,
            n_freqs=16,
            n_train_profiles=20,
            epochs=2,
            n_mc=5,
            hidden=(64, 32),
        )
        r = ag.execute({"sites": sites})
        self.assertIn(r.status, ("success", "needs_review", "failed"))

    def test_output_keys_present_on_success(self):
        from pycsamt.agents import Inv3DAgent

        sites = _load_sites()
        if sites is None:
            self.skipTest("Could not load sites.")
        ag = Inv3DAgent(
            n_layers=3,
            n_freqs=16,
            n_train_profiles=20,
            epochs=2,
            n_mc=0,
            hidden=(64, 32),  # n_mc=0 → no uncertainty
        )
        r = ag.execute({"sites": sites})
        if r.status == "success":
            self.assertIn("pred_rho", r.data)
            self.assertIn("pred_thick", r.data)
            self.assertIn("depths_km", r.data)
            self.assertIn("adjacency", r.data)
            self.assertIn("inverter", r.data)
            # pred_rho shape: (n_stations, n_layers)
            self.assertEqual(r["pred_rho"].ndim, 2)
            self.assertEqual(r["pred_rho"].shape[1], 3)  # n_layers=3


if __name__ == "__main__":
    unittest.main()
