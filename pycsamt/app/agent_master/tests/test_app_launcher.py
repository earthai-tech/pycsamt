"""Tests for the application-launch task (web / mapview / desktop).

Covers request detection routing and the launch card component.
Launcher functions that spawn real subprocesses are not exercised here.
"""
from __future__ import annotations

import unittest

import pycsamt.app.agent_master.callbacks.chat as C


class TestDetectAppRequest(unittest.TestCase):
    def test_desktop_requests(self):
        for text in (
            "open the desktop app",
            "launch the gui",
            "start the desktop application",
            "open pycsamt desktop",
        ):
            kind, reason = C._detect_app_request(text)
            self.assertEqual(kind, "desktop", text)
            self.assertEqual(reason, "")

    def test_mapview_requests(self):
        for text in (
            "open the map view",
            "launch mapview",
            "open the map workbench",
        ):
            kind, reason = C._detect_app_request(text)
            self.assertEqual(kind, "mapview", text)
            self.assertEqual(reason, "")

    def test_map_viz_redirects_to_mapview(self):
        for text in (
            "show me a 3d map of the survey",
            "I want an interactive pseudosection",
        ):
            kind, reason = C._detect_app_request(text)
            self.assertEqual(kind, "mapview", text)
            self.assertTrue(reason)

    def test_web_requests(self):
        kind, reason = C._detect_app_request("open the web app")
        self.assertEqual(kind, "web")
        self.assertEqual(reason, "")
        kind, reason = C._detect_app_request(
            "open the full pipeline editor"
        )
        self.assertEqual(kind, "web")
        self.assertTrue(reason)

    def test_plot_workflows_not_hijacked(self):
        # chat plot workflows must stay in the chat
        for text in (
            "plot the phase tensor map of L22",
            "plot the pseudosection",
            "run qc",
            "read the data",
        ):
            self.assertIsNone(C._detect_app_request(text), text)


class TestLaunchBubble(unittest.TestCase):
    def test_web_card_has_link(self):
        b = C._launch_bubble("web", url="http://127.0.0.1:8051")
        s = str(b)
        self.assertIn("http://127.0.0.1:8051", s)
        self.assertIn("am-webapp-link", s)
        self.assertIn("Launching pyCSAMT Web App", s)

    def test_mapview_card(self):
        b = C._launch_bubble(
            "mapview", url="http://127.0.0.1:8770"
        )
        s = str(b)
        self.assertIn("Launching MapView", s)
        self.assertIn("8770", s)

    def test_desktop_card_no_link(self):
        b = C._launch_bubble(
            "desktop", note="A native pyCSAMT window should appear."
        )
        s = str(b)
        self.assertIn("Launching pyCSAMT Desktop", s)
        self.assertNotIn("am-webapp-link", s)
        self.assertIn("native pyCSAMT window", s)

    def test_desktop_failure_card(self):
        b = C._launch_bubble(
            "desktop", note="PySide6 is probably not installed.",
            ok=False,
        )
        s = str(b)
        self.assertIn("Could not launch", s)
        self.assertIn("PySide6", s)


if __name__ == "__main__":
    unittest.main()
