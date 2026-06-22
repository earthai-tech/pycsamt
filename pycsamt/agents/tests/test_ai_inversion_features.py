# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
The observed feature vector must match the inverter's training layout.

Regression: _z_to_features returned |Z| over both off-diagonals as an
(n, 4) array, but EMInverter1D is trained on a flat
``[log10(rho_a_xy), phase_xy]`` vector of length 2 * n_freq
(see pycsamt.ai._base._z_list_to_array). The mismatch made every
prediction fail ("operands could not be broadcast … (1,40,4) (80,)"),
so the inversion produced 0 stations and no figure on valid data.
"""
from __future__ import annotations

import os
import unittest

import numpy as np

from pycsamt.agents.ai_inversion import _z_to_features
from pycsamt.ai._base import _z_list_to_array


class _FakeZ:
    """Minimal Z exposing resistivity_xy / phase_xy like the real class."""
    def __init__(self, freq, rho_xy, phase_xy):
        self.freq = np.asarray(freq, float)
        self.resistivity_xy = np.asarray(rho_xy, float)
        self.phase_xy = np.asarray(phase_xy, float)


class TestFeatureLayout(unittest.TestCase):

    def setUp(self):
        # observed grid deliberately differs from the target grid
        self.fr = np.logspace(0, 4, 53)
        self.rho = np.full(53, 100.0)
        self.pha = np.full(53, 45.0)
        self.zobj = _FakeZ(self.fr, self.rho, self.pha)
        self.z = np.zeros((53, 2, 2), complex)  # unused (resistivity_xy present)
        self.freqs_target = np.logspace(-4, 3, 40)

    def test_feature_length_matches_training_width(self):
        feat = _z_to_features(
            self.zobj, self.z, self.fr, self.freqs_target
        )
        self.assertIsNotNone(feat)
        # flat vector, length 2 * n_target (rho + phase)
        self.assertEqual(feat.shape, (2 * len(self.freqs_target),))
        self.assertTrue(np.all(np.isfinite(feat)))

    def test_matches_z_list_to_array_width(self):
        # A training row built by the inverter's own helper, on the
        # target grid, must have the same width as the observed feature.
        zt = _FakeZ(
            self.freqs_target,
            np.full(40, 100.0), np.full(40, 45.0),
        )
        train_row = _z_list_to_array(
            [zt], include_phase=True, log_rho=True
        )
        feat = _z_to_features(
            self.zobj, self.z, self.fr, self.freqs_target
        )
        self.assertEqual(feat.shape[0], train_row.shape[1])

    def test_too_few_points_returns_none(self):
        zobj = _FakeZ([10.0], [100.0], [45.0])
        feat = _z_to_features(
            zobj, np.zeros((1, 2, 2), complex),
            np.array([10.0]), self.freqs_target,
        )
        self.assertIsNone(feat)

    def test_fallback_from_impedance_when_no_rho(self):
        # Z object without resistivity_xy → compute from the tensor.
        class _Bare:
            pass
        z = np.zeros((20, 2, 2), complex)
        z[:, 0, 1] = (1 + 1j)
        fr = np.logspace(0, 3, 20)
        feat = _z_to_features(_Bare(), z, fr, self.freqs_target)
        self.assertIsNotNone(feat)
        self.assertEqual(feat.shape, (2 * len(self.freqs_target),))


if __name__ == "__main__":
    unittest.main(verbosity=2)
