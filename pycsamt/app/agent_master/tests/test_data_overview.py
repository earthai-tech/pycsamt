"""Tests for the inline data-overview task.

"read the EDI data / stations / sites" is answered deterministically from
the stored survey: detection gate, no-data guidance, payload shaping, and
the overview card component.
"""
from __future__ import annotations

import unittest

import pycsamt.app.agent_master.callbacks.chat as C


class TestDataReadGate(unittest.TestCase):
    def test_positive_phrasings(self):
        for text in (
            "read the edi data",
            "read edi data",
            "read stations",
            "read the sites",
            "show loaded data",
            "list the stations",
            "describe the survey",
            "check the dataset",
            "statistics of the data",
            "survey data overview",
            "what data is loaded?",
        ):
            self.assertTrue(
                C._looks_like_data_read(text), text
            )

    def test_negative_phrasings(self):
        for text in (
            "plot the stations on a map",
            "run qc on the data",
            "invert the data with occam2d",
            "generate code to read edi files",
            "what is an edi file?",
            "correct the static shift",
            "hello there",
            "export the sites to csv",
        ):
            self.assertFalse(
                C._looks_like_data_read(text), text
            )


class TestOverviewPayloadAndCard(unittest.TestCase):
    LINES = [
        {
            "label": "L22", "n_stations": 24,
            "stations": [f"S{i:02d}" for i in range(24)],
            "freq": (0.125, 8192.0), "max_nfreq": 40,
            "qc": 82.0, "flagged": 2, "length_km": 3.2,
            "tipper": True,
        },
        {
            "label": "K1", "n_stations": 18,
            "stations": [f"K{i:02d}" for i in range(18)],
            "freq": (1.0, 1024.0), "max_nfreq": 30,
            "qc": 71.0, "flagged": 0, "length_km": 2.1,
            "tipper": False,
        },
    ]

    def test_payload_totals(self):
        card = C._overview_payload(list(self.LINES), [])
        tiles = {t["label"]: t["value"] for t in card["tiles"]}
        self.assertEqual(tiles["stations"], "42")
        self.assertIn("0.125 Hz", tiles["frequency range"])
        self.assertIn("8.19 kHz", tiles["frequency range"])
        self.assertEqual(tiles["mean QC score"], "76/100")
        self.assertIn("5.3 km", tiles["total length"])
        self.assertIn("2 lines", card["scope"])
        self.assertIn("tipper present", card["scope"])
        # station chips only for single-line overviews
        self.assertEqual(card["stations"], [])

    def test_payload_single_line_has_chips(self):
        card = C._overview_payload([self.LINES[0]], [])
        self.assertEqual(len(card["stations"]), 24)
        self.assertEqual(
            card["tiles"][2]["label"], "profile length"
        )

    def test_card_component_builds(self):
        card = C._overview_payload(
            list(self.LINES), ["L22: 2 of 24 station(s) flagged"]
        )
        div = C._data_overview_card(card)
        self.assertEqual(div.className, "am-ov-card")
        # 2 lines → the per-line table is present
        rendered = str(div)
        self.assertIn("am-ov-table", rendered)
        self.assertIn("am-ov-warn", rendered)
        self.assertIn("am-ov-hint", rendered)

    def test_bubble_embeds_card(self):
        card = C._overview_payload([self.LINES[0]], [])
        bub = C._agent_bubble(
            "Read **L22**.", kind=C.KIND_ANSWER, card=card,
        )
        self.assertIn("am-ov-card", str(bub))


class TestLinesQuery(unittest.TestCase):
    GROUPS = {
        "L22PLT": [f"s{i}.edi" for i in range(24)],
        "K1": [f"k{i}.edi" for i in range(18)],
    }

    def test_gate(self):
        for text in (
            "what are the lines ?",
            "which lines are loaded",
            "list the lines",
            "how many lines do we have",
        ):
            self.assertTrue(
                C._looks_like_lines_query(text), text
            )
            self.assertTrue(
                C._looks_like_data_read(text), text
            )
        self.assertFalse(
            C._looks_like_lines_query("read the data")
        )

    def test_lines_reply_from_store(self):
        jid = C._new_job()
        C._dispatch_data_overview(
            jid, "what are the lines?",
            {"groups": self.GROUPS, "path": "x"}, {},
            step=lambda *a, **k: None,
        )
        job = C._get_job(jid)
        self.assertEqual(job["kind"], C.KIND_ANSWER)
        self.assertIn("2 survey lines are loaded", job["result"])
        self.assertIn("**K1** — 18 EDI files", job["result"])
        self.assertIn("read line K1", job["result"])

    def test_lines_reply_no_data(self):
        jid = C._new_job()
        C._dispatch_data_overview(
            jid, "what are the lines?", {}, {},
            step=lambda *a, **k: None,
        )
        job = C._get_job(jid)
        self.assertEqual(job["kind"], C.KIND_META)
        self.assertIn("Load EDI", job["result"])


class TestSmartUnknownReply(unittest.TestCase):
    STORE = {
        "groups": {
            "L22PLT": ["a.edi", "b.edi"],
            "K1": ["c.edi"],
        },
        "path": "x",
    }

    def test_unknown_line_proposes_loaded_lines(self):
        out = C._smart_unknown_reply(
            "plot only the line B9", self.STORE
        )
        self.assertIn("couldn't find a line", out)
        self.assertIn("**L22PLT**", out)
        self.assertIn("**K1**", out)

    def test_ordinal_line_gets_plot_menu(self):
        # "line 2" → the 2nd loaded line (sorted: K1, L22PLT)
        out = C._smart_unknown_reply(
            "plot only the line 2", self.STORE
        )
        self.assertIn("**L22PLT**", out)
        self.assertIn("phase pseudosection of L22PLT", out)
        self.assertIn("“line 2”", out)

    def test_valid_line_gets_plot_menu(self):
        out = C._smart_unknown_reply(
            "draw something for K1", self.STORE
        )
        self.assertIn("for **K1**", out)
        self.assertIn("plot the sounding curves of K1", out)

    def test_plotish_without_line(self):
        out = C._smart_unknown_reply(
            "plot the data nicely", {"path": "x"}
        )
        self.assertIn("which figure", out)
        self.assertIn("plot the phase pseudosection", out)

    def test_non_plot_request_returns_none(self):
        self.assertIsNone(
            C._smart_unknown_reply(
                "make me a coffee", self.STORE
            )
        )


class TestDispatchNoData(unittest.TestCase):
    def test_no_data_guidance(self):
        jid = C._new_job()
        steps: list = []

        def step(label, status="done"):
            steps.append((label, status))

        C._dispatch_data_overview(
            jid, "read the edi data", {}, {}, step=step,
        )
        job = C._get_job(jid)
        self.assertEqual(job["status"], "done")
        self.assertEqual(job["kind"], C.KIND_META)
        self.assertIn("Load EDI", job["result"])
        self.assertIsNone(job.get("card"))


if __name__ == "__main__":
    unittest.main()
