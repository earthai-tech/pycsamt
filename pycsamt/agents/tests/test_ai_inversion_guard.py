# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
AI inversion must fail fast when no station has usable impedance.

Regression: the agent trained a network on synthetic data (minutes)
*before* touching the real EDIs, then predicted 0 stations and produced
an empty result. It now checks for at least one station with usable
impedance features up front and returns immediately.

Skipped if no DL backend (torch/tensorflow) is installed.
"""

from __future__ import annotations

import time
import unittest

import pycsamt.agents.ai_inversion as aimod
from pycsamt.agents.ai_inversion import AIInversionAgent
from pycsamt.agents.inv3d_agent import Inv3DAgent


def test_inv3d_agent_accepts_explicit_early_stopping_patience():
    agent = Inv3DAgent(epochs=50, patience=10, n_mc=0)
    assert agent.epochs == 50
    assert agent.patience == 10


def _backend_available() -> bool:
    try:
        from pycsamt.backends import get_backend_instance

        return get_backend_instance() is not None
    except Exception:
        return False


@unittest.skipUnless(_backend_available(), "no DL backend")
class TestAIInversionFastFail(unittest.TestCase):
    def setUp(self):
        self._orig = aimod._z_to_features

    def tearDown(self):
        aimod._z_to_features = self._orig

    def test_no_usable_impedance_fails_before_training(self):
        # Force every station to yield no features → guard must trip.
        aimod._z_to_features = lambda *a, **k: None

        agent = AIInversionAgent(
            arch="resnet",
            n_layers=5,
            epochs=30,
            n_train_samples=2000,
        )
        t0 = time.time()
        res = agent.execute({"path": "data/3edis"})
        elapsed = time.time() - t0

        self.assertEqual(res.status, "failed")
        self.assertIn("usable impedance", (res.error or "").lower())
        # Fast: must NOT have trained (training on 2000 samples would
        # take many seconds/minutes; the guard returns near-instantly).
        self.assertLess(elapsed, 20.0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
