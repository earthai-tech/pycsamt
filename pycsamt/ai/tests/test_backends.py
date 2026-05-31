# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.backends — detection, graceful degradation, API."""
import unittest


class TestBackendDetection(unittest.TestCase):

    def test_get_backend_returns_string(self):
        from pycsamt.backends import get_backend
        name = get_backend()
        self.assertIsInstance(name, str)
        self.assertIn(name, {"torch", "tensorflow", "none"})

    def test_list_backends_keys(self):
        from pycsamt.backends import list_backends
        avail = list_backends()
        self.assertIn("torch", avail)
        self.assertIn("tensorflow", avail)
        for v in avail.values():
            self.assertIsInstance(v, bool)

    def test_get_backend_instance_does_not_crash(self):
        from pycsamt.backends import get_backend_instance
        inst = get_backend_instance()
        # Either a NeuralBackend subclass or None (no framework)
        self.assertTrue(inst is None or hasattr(inst, "resolve_device"))

    def test_set_backend_invalid_raises(self):
        from pycsamt.backends import set_backend
        with self.assertRaises(ValueError):
            set_backend("foobar")

    def test_set_backend_auto_no_crash(self):
        """set_backend('auto') must not crash even when no framework is present."""
        from pycsamt.backends import set_backend, get_backend
        # Reset to auto — should select first available or 'none'
        set_backend("auto")
        name = get_backend()
        self.assertIn(name, {"torch", "tensorflow", "none"})

    def test_backend_config_singleton(self):
        from pycsamt.backends._config import BackendConfig
        a = BackendConfig()
        b = BackendConfig()
        self.assertIs(a, b)

    def test_active_backend_utility(self):
        from pycsamt.ai._backend_utils import active_backend
        name = active_backend()
        self.assertIn(name, {"torch", "tensorflow", "none"})

    def test_resolve_device_no_crash(self):
        from pycsamt.ai._backend_utils import resolve_device
        dev = resolve_device(None)
        self.assertIsInstance(dev, str)

    def test_probe_backend(self):
        from pycsamt.backends import probe_backend
        for name in ("torch", "tensorflow"):
            result = probe_backend(name)
            self.assertIsInstance(result, bool)


if __name__ == "__main__":
    unittest.main()
