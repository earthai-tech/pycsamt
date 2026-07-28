# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Natural-language workflow routing benchmark — Experiment 1.

Measures the offline (regex-only) classifier accuracy across
all supported workflow types.  Each entry in ``_CASES`` is:

    (request_text, expected_workflow, expected_data_path)

``expected_data_path`` may be ``None`` when the request
intentionally contains no path.

Run::

    pytest pycsamt/agents/tests/test_nl_routing.py -v

Metrics reported at the end of the test session (via a
module-level fixture):
  * workflow classification accuracy
  * path extraction accuracy (when path is expected)
  * overall JSON / schema validity rate
"""

from __future__ import annotations

import unittest

# ── test cases ────────────────────────────────────────────────────────────────
# (request, expected_workflow, expected_path_or_None)
_CASES: list[tuple[str, str, str | None]] = [
    # ── qc ────────────────────────────────────────────────────────────────────
    (
        "Load EDI files from /data/willy and run quality control",
        "qc",
        "/data/willy",
    ),
    (
        "Check data quality for the WILLY survey at /mnt/edi/willy",
        "qc",
        "/mnt/edi/willy",
    ),
    (
        "Flag bad stations and compute SNR for /data/L22PLT",
        "qc",
        "/data/L22PLT",
    ),
    (
        "Clean noisy frequencies in /surveys/line1",
        "qc",
        "/surveys/line1",
    ),
    # ── phase_analysis ────────────────────────────────────────────────────────
    (
        "Run phase tensor analysis on /data/willy EDIs",
        "phase_analysis",
        "/data/willy",
    ),
    (
        "Compute strike and dimensionality for /data/L22PLT",
        "phase_analysis",
        "/data/L22PLT",
    ),
    (
        "Phase tensor, Mohr circles, and skew analysis /data/profile1",
        "phase_analysis",
        "/data/profile1",
    ),
    (
        "Argand diagram and PT analysis at /home/user/edi",
        "phase_analysis",
        "/home/user/edi",
    ),
    # ── pre_inversion ─────────────────────────────────────────────────────────
    (
        "Prepare Occam2D input files for /data/willy",
        "pre_inversion",
        "/data/willy",
    ),
    (
        "Write Occam2D mesh and startup file from /data/L22PLT",
        "pre_inversion",
        "/data/L22PLT",
    ),
    (
        "Pre-inversion preparation using Occam2D for /data/profile1",
        "pre_inversion",
        "/data/profile1",
    ),
    (
        "Prepare inversion data files for /surveys/csamt",
        "pre_inversion",
        "/surveys/csamt",
    ),
    # ── ai_inversion ──────────────────────────────────────────────────────────
    (
        "Run AI inversion on /data/willy EDI files",
        "ai_inversion",
        "/data/willy",
    ),
    (
        "Deep learning 1D inversion for /data/L22PLT",
        "ai_inversion",
        "/data/L22PLT",
    ),
    (
        "Use CNN inverter on /data/profile1",
        "ai_inversion",
        "/data/profile1",
    ),
    (
        "Run neural network inversion on /surveys/mt_data",
        "ai_inversion",
        "/surveys/mt_data",
    ),
    (
        "1D AI inversion with machine learning on /home/user/edi",
        "ai_inversion",
        "/home/user/edi",
    ),
    (
        "Apply inv1d deep learning model to /data/willy",
        "ai_inversion",
        "/data/willy",
    ),
    # ── inv2d ─────────────────────────────────────────────────────────────────
    (
        "Run U-Net 2D neural profile inversion on /data/willy",
        "inv2d",
        "/data/willy",
    ),
    (
        "2D AI inversion using deep learning for /data/L22PLT",
        "inv2d",
        "/data/L22PLT",
    ),
    (
        "Profile inversion with lateral continuity AI /data/profile1",
        "inv2d",
        "/data/profile1",
    ),
    (
        "u-net 2d inversion for /surveys/line2",
        "inv2d",
        "/surveys/line2",
    ),
    # ── inv3d ─────────────────────────────────────────────────────────────────
    (
        "Run GCN 3D spatial AI inversion on /data/willy",
        "inv3d",
        "/data/willy",
    ),
    (
        "Graph convolutional network inversion for /data/L22PLT",
        "inv3d",
        "/data/L22PLT",
    ),
    (
        "3D AI inversion using spatial graph model /data/3d_survey",
        "inv3d",
        "/data/3d_survey",
    ),
    (
        "Graph network inversion on /surveys/block1",
        "inv3d",
        "/surveys/block1",
    ),
    # ── ensemble_inversion ────────────────────────────────────────────────────
    (
        "Ensemble inversion with uncertainty bands for /data/willy",
        "ensemble_inversion",
        "/data/willy",
    ),
    (
        "Run Bayesian inversion with confidence intervals /data/L22PLT",
        "ensemble_inversion",
        "/data/L22PLT",
    ),
    (
        "Uncertainty quantification for AI inversion on /data/profile1",
        "ensemble_inversion",
        "/data/profile1",
    ),
    (
        "Calibrated inversion ensemble for /surveys/mt1",
        "ensemble_inversion",
        "/surveys/mt1",
    ),
    # ── joint_inversion ───────────────────────────────────────────────────────
    (
        "Joint multi-modal TEM+MT inversion for /data/willy",
        "joint_inversion",
        "/data/willy",
    ),
    (
        "Run joint inversion combining MT and TEM on /data/L22PLT",
        "joint_inversion",
        "/data/L22PLT",
    ),
    (
        "Multi-physics joint inversion /data/joint_survey",
        "joint_inversion",
        "/data/joint_survey",
    ),
    # ── modem ─────────────────────────────────────────────────────────────────
    (
        "Write ModEM 3D data file for /data/willy",
        "modem",
        "/data/willy",
    ),
    (
        "Prepare ModEM input for 3D inversion at /data/L22PLT",
        "modem",
        "/data/L22PLT",
    ),
    # ── tipper ────────────────────────────────────────────────────────────────
    (
        "Tipper analysis for /data/willy",
        "tipper",
        "/data/willy",
    ),
    (
        "Compute Parkinson vectors and Wiese arrows /data/L22PLT",
        "tipper",
        "/data/L22PLT",
    ),
    # ── tipper_plot ───────────────────────────────────────────────────────────
    # "induction arrows" is an exact tipper_plot keyword, so the mixed
    # phrasing routes to the plotting workflow by design
    (
        "Tipper analysis and induction arrows for /data/willy",
        "tipper_plot",
        "/data/willy",
    ),
    # ── sensitivity ───────────────────────────────────────────────────────────
    (
        "Compute Bostick depth of investigation for /data/willy",
        "sensitivity",
        "/data/willy",
    ),
    (
        "Depth of investigation and sensitivity kernels /data/L22PLT",
        "sensitivity",
        "/data/L22PLT",
    ),
    # ── rotation ──────────────────────────────────────────────────────────────
    (
        "Rotate tensors to principal strike axis for /data/willy",
        "rotation",
        "/data/willy",
    ),
    (
        "Tensor rotation and coordinate frame correction /data/L22PLT",
        "rotation",
        "/data/L22PLT",
    ),
    # ── freq_decimation ───────────────────────────────────────────────────────
    (
        "Decimate frequencies for /data/willy EDI files",
        "freq_decimation",
        "/data/willy",
    ),
    (
        "Period selection and dead band removal /data/L22PLT",
        "freq_decimation",
        "/data/L22PLT",
    ),
    # ── batch ─────────────────────────────────────────────────────────────────
    (
        "Batch process all profiles in /surveys/willy",
        "batch",
        "/surveys/willy",
    ),
    (
        "Run all profiles in parallel for /data/batch_survey",
        "batch",
        "/data/batch_survey",
    ),
    # ── full ──────────────────────────────────────────────────────────────────
    (
        "Run full pipeline on /data/willy EDIs",
        "full",
        "/data/willy",
    ),
    (
        "Complete processing workflow all steps for /data/L22PLT",
        "full",
        "/data/L22PLT",
    ),
    # ── full_ai_workflow ──────────────────────────────────────────────────────
    (
        "Full AI workflow on /data/willy including inversion",
        "full_ai_workflow",
        "/data/willy",
    ),
    (
        "Run full AI pipeline with all steps on /data/L22PLT",
        "full_ai_workflow",
        "/data/L22PLT",
    ),
    # ── path edge cases ───────────────────────────────────────────────────────
    (
        "Run AI inversion on /mnt/data/edi5 save to /out/results",
        "ai_inversion",
        "/mnt/data/edi5",
    ),
    (
        "Load /data/WILLY_DATA/L22PLT run phase tensor analysis",
        "phase_analysis",
        "/data/WILLY_DATA/L22PLT",
    ),
    (
        "QC analysis for /home/daniel/projects/data/willy",
        "qc",
        "/home/daniel/projects/data/willy",
    ),
    # ── no path cases ─────────────────────────────────────────────────────────
    (
        "Run quality control on the survey data",
        "qc",
        None,
    ),
    (
        "Perform phase tensor analysis",
        "phase_analysis",
        None,
    ),
]

# ── minimum acceptable accuracy thresholds ───────────────────────────────────
_MIN_WORKFLOW_ACCURACY = 0.85  # 85% workflow classification correct
_MIN_PATH_ACCURACY = 0.90  # 90% path extraction correct (when expected)


class TestNLRoutingBenchmark(unittest.TestCase):
    """
    Regression suite for the offline NL workflow router.

    Each test_case_* method checks one (request, workflow, path)
    triple.  The test_summary method checks aggregate accuracy.
    """

    @classmethod
    def setUpClass(cls):
        from pycsamt.agents import ContextInputAgent

        cls.agent = ContextInputAgent()  # no API key — regex only
        cls.results: list[dict] = []

        for req, exp_wf, exp_path in _CASES:
            r = cls.agent.execute({"request": req})
            cfg = r["config"]
            cls.results.append(
                {
                    "request": req,
                    "expected_wf": exp_wf,
                    "got_wf": cfg.get("workflow", ""),
                    "expected_path": exp_path,
                    "got_path": cfg.get("data_path", ""),
                    "wf_ok": cfg.get("workflow", "") == exp_wf,
                    "path_ok": (
                        exp_path is None or cfg.get("data_path", "") == exp_path
                    ),
                }
            )

    # ── per-workflow-type subtests ────────────────────────────────────────────

    def _check_group(self, expected_wf: str) -> None:
        cases = [r for r in self.results if r["expected_wf"] == expected_wf]
        if not cases:
            self.skipTest(f"No cases for workflow {expected_wf!r}")
        wrong = [c for c in cases if not c["wf_ok"]]
        self.assertEqual(
            wrong,
            [],
            msg=(
                f"Workflow {expected_wf!r}: "
                f"{len(wrong)}/{len(cases)} misclassified.\n"
                + "\n".join(f"  [{c['got_wf']!r}] {c['request'][:60]}" for c in wrong)
            ),
        )

    def test_qc_classification(self):
        self._check_group("qc")

    def test_phase_analysis_classification(self):
        self._check_group("phase_analysis")

    def test_pre_inversion_classification(self):
        self._check_group("pre_inversion")

    def test_ai_inversion_classification(self):
        self._check_group("ai_inversion")

    def test_inv2d_classification(self):
        self._check_group("inv2d")

    def test_inv3d_classification(self):
        self._check_group("inv3d")

    def test_ensemble_inversion_classification(self):
        self._check_group("ensemble_inversion")

    def test_joint_inversion_classification(self):
        self._check_group("joint_inversion")

    def test_modem_classification(self):
        self._check_group("modem")

    def test_tipper_classification(self):
        self._check_group("tipper")

    def test_tipper_plot_classification(self):
        self._check_group("tipper_plot")

    def test_sensitivity_classification(self):
        self._check_group("sensitivity")

    def test_rotation_classification(self):
        self._check_group("rotation")

    def test_freq_decimation_classification(self):
        self._check_group("freq_decimation")

    def test_batch_classification(self):
        self._check_group("batch")

    def test_full_classification(self):
        self._check_group("full")

    def test_full_ai_workflow_classification(self):
        self._check_group("full_ai_workflow")

    # ── path extraction ───────────────────────────────────────────────────────

    def test_path_extraction_accuracy(self):
        """At least 90% of expected paths are extracted correctly."""
        path_cases = [r for r in self.results if r["expected_path"]]
        n_ok = sum(1 for r in path_cases if r["path_ok"])
        accuracy = n_ok / len(path_cases) if path_cases else 0.0
        wrong = [r for r in path_cases if not r["path_ok"]]
        self.assertGreaterEqual(
            accuracy,
            _MIN_PATH_ACCURACY,
            msg=(
                f"Path accuracy {accuracy:.1%} < {_MIN_PATH_ACCURACY:.0%}.\n"
                "Failed extractions:\n"
                + "\n".join(
                    f"  got={r['got_path']!r} "
                    f"want={r['expected_path']!r}\n"
                    f"  req: {r['request'][:70]}"
                    for r in wrong
                )
            ),
        )

    # ── aggregate accuracy ────────────────────────────────────────────────────

    def test_overall_workflow_accuracy(self):
        """At least 85% of all NL requests are classified correctly."""
        n_ok = sum(1 for r in self.results if r["wf_ok"])
        accuracy = n_ok / len(self.results) if self.results else 0.0
        wrong = [r for r in self.results if not r["wf_ok"]]
        self.assertGreaterEqual(
            accuracy,
            _MIN_WORKFLOW_ACCURACY,
            msg=(
                f"Workflow accuracy {accuracy:.1%} < "
                f"{_MIN_WORKFLOW_ACCURACY:.0%}.\n"
                "Misclassified:\n"
                + "\n".join(
                    f"  [{r['got_wf']!r} != {r['expected_wf']!r}] "
                    f"{r['request'][:60]}"
                    for r in wrong
                )
            ),
        )

    # ── WorkflowPlan schema ───────────────────────────────────────────────────

    def test_workflow_plan_returned(self):
        """ContextInputAgent always returns a WorkflowPlan object."""
        from pycsamt.agents import (
            ContextInputAgent,
            WorkflowPlan,
        )

        ag = ContextInputAgent()
        r = ag.execute({"request": "AI inversion on /data/willy EDIs"})
        self.assertIn("workflow_plan", r.data)
        self.assertIsInstance(r["workflow_plan"], WorkflowPlan)

    def test_workflow_plan_workflow_type_matches_config(self):
        """WorkflowPlan.workflow_type must equal config['workflow']."""
        from pycsamt.agents import ContextInputAgent

        ag = ContextInputAgent()
        for req, _exp_wf, _ in _CASES[:10]:
            r = ag.execute({"request": req})
            plan = r["workflow_plan"]
            cfg = r["config"]
            self.assertEqual(
                plan.workflow_type,
                cfg["workflow"],
                msg=f"Mismatch for: {req!r}",
            )

    def test_workflow_plan_data_path_matches_config(self):
        """WorkflowPlan.data_path must equal config['data_path']."""
        from pycsamt.agents import ContextInputAgent

        ag = ContextInputAgent()
        r = ag.execute({"request": "QC on /data/willy save to /out"})
        self.assertEqual(
            r["workflow_plan"].data_path,
            r["config"]["data_path"],
        )

    def test_no_path_request_returns_needs_review(self):
        """Requests with no path get status='needs_review'."""
        from pycsamt.agents import ContextInputAgent

        ag = ContextInputAgent()
        r = ag.execute({"request": "Run quality control"})
        self.assertEqual(r.status, "needs_review")


# ── standalone accuracy reporter ─────────────────────────────────────────────


def _print_accuracy_report():
    """Print a summary table when this module is run directly."""
    from pycsamt.agents import ContextInputAgent

    ag = ContextInputAgent()

    wf_totals: dict[str, list[bool]] = {}
    path_hits: list[bool] = []

    for req, exp_wf, exp_path in _CASES:
        r = ag.execute({"request": req})
        cfg = r["config"]
        got_wf = cfg.get("workflow", "")
        ok = got_wf == exp_wf
        wf_totals.setdefault(exp_wf, []).append(ok)
        if exp_path:
            path_hits.append(cfg.get("data_path", "") == exp_path)

    print("\n=== NL Routing Benchmark ===")
    print(f"{'Workflow':<25} {'Correct':>8} {'Total':>7} {'Acc':>6}")
    print("-" * 50)
    overall = []
    for wf in sorted(wf_totals):
        hits = wf_totals[wf]
        acc = sum(hits) / len(hits)
        overall.extend(hits)
        mark = "" if acc >= _MIN_WORKFLOW_ACCURACY else " !"
        print(f"{wf:<25} {sum(hits):>8} {len(hits):>7} {acc:>5.1%}{mark}")
    n = len(overall)
    print("-" * 50)
    print(f"{'OVERALL':<25} {sum(overall):>8} {n:>7} {sum(overall) / n:>5.1%}")
    if path_hits:
        pa = sum(path_hits) / len(path_hits)
        print(
            f"\nPath extraction accuracy: "
            f"{sum(path_hits)}/{len(path_hits)} = {pa:.1%}"
        )


if __name__ == "__main__":
    _print_accuracy_report()
