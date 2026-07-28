# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for ToolAgent (strike / dimensionality / validator)."""

from __future__ import annotations

import os
import unittest

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from pycsamt.agents.tooling import (  # noqa: E402
    _EXPORT_BUNDLES,
    TOOL_KINDS,
    ToolAgent,
    _circular_strike_mean,
    _df_to_text,
    _get_latlon,
    _ll_to_utm,
    _safe_filename,
)

_DATA = os.path.join("data", "3edis")
_HAS_DATA = os.path.isdir(_DATA)

# A larger line, used for the Stratagem EDI-native pipeline (Wave E).
_WILLY = os.path.join("data", "AMT", "WILLY_DATA", "L18PLT")
_HAS_WILLY = os.path.isdir(_WILLY)


class TestToolHelpers(unittest.TestCase):
    def test_kinds(self):
        self.assertEqual(
            set(TOOL_KINDS),
            {
                "strike",
                "dimensionality",
                "validator",
                "coords",
                "elevation",
                "converter",
                "batch_export",
                "freq_editor",
                "layered_model",
                "correction",
            },
        )

    def test_correction_registry_coercion(self):
        # registry metadata + parameter coercion (no EDI data needed).
        from pycsamt.agents._corrections import (
            CORRECTION_METHODS,
            coerce_kwargs,
            fn_for,
            param_specs,
        )

        self.assertIn("corr_ss_ama", CORRECTION_METHODS)
        self.assertEqual(fn_for("corr_ss_ama"), "correct_ss_ama")
        names = {p.name for p in param_specs("corr_ss_ama")}
        self.assertTrue({"half_window", "max_skew"}.issubset(names))
        # str inputs from the modal are coerced to the right types; bad/blank
        # values fall back to the ParamSpec default.
        kw = coerce_kwargs(
            "corr_ss_ama",
            {"half_window": "5", "max_skew": "", "weights": "box"},
        )
        self.assertEqual(kw["half_window"], 5)
        self.assertIsInstance(kw["half_window"], int)
        self.assertEqual(kw["weights"], "box")
        self.assertEqual(kw["max_skew"], 6.0)  # default for blank

    def test_correction_missing_selection(self):
        # correction with no method selected → failed, not a crash.
        r = ToolAgent().execute({"path": _DATA, "kind": "correction"})
        self.assertEqual(r.status, "failed")

    def test_static_shift_method_routing(self):
        # Specific static-shift phrases route to their correction workflow,
        # while the bare phrase still routes to the legacy static_shift wf.
        from pycsamt.agents._workflows import (
            classify_workflow,
        )

        cases = {
            "apply ama static shift": "corr_ss_ama",
            "loess static shift correction": "corr_ss_loess",
            "bilateral static shift": "corr_ss_bilateral",
            "reference median static shift": "corr_ss_refmedian",
            "hanning emap static shift": "corr_ss_emap",
            "remove static shift": "static_shift",
            "apply a notch filter": "corr_notch",
            "log frequency smoothing": "corr_smooth_logfreq",
            "smooth rho phase trend": "corr_smooth_rho_phase",
            "rotate by 30 degrees": "corr_rotate_angle",
            "rotate to geoelectric strike": "corr_rotate_strike",
            "rotate to phase tensor strike": "corr_rotate_pt_strike",
            "rotate to profile azimuth": "corr_rotate_profile",
            "antisymmetrize the tensor": "corr_antisymmetrize",
            "project stations onto the line": "corr_coord_projection",
            "regularize station spacing": "corr_coord_spacing",
            "snap outlier stations": "corr_coord_snap",
            "smooth elevation along profile": "corr_coord_elevation",
            "shift coordinates": "corr_coord_shift",
            "interpolate missing coordinates": "corr_coord_interpolate",
            "near field correction": "corr_near_field",
            "stratagem qc report": "corr_strat_qc",
            "stratagem static shift": "corr_strat_static_shift",
            "stratagem noise removal": "corr_strat_noise",
            "stratagem frequency filter": "corr_strat_freq_filter",
            "stratagem full pipeline": "corr_strat_full",
        }
        for text, expected in cases.items():
            self.assertEqual(classify_workflow(text), expected, msg=text)

    def test_correction_schemas_generated(self):
        # Every registered correction method gets a param-modal schema whose
        # field keys match its catalogue ParamSpec names.
        from pycsamt.agents._corrections import (
            CORRECTION_METHODS,
            param_specs,
        )
        from pycsamt.app.agent_master.callbacks.params import (
            _SCHEMAS,
        )

        for wf in CORRECTION_METHODS:
            self.assertIn(wf, _SCHEMAS, msg=wf)
            keys = {f["key"] for f in _SCHEMAS[wf]["fields"]}
            for ps in param_specs(wf):
                self.assertIn(ps.name, keys, msg=f"{wf}:{ps.name}")

    def test_all_catalogue_methods_registered(self):
        # Every catalogue method (every category) is exposed as a workflow id.
        from pycsamt.agents._corrections import (
            CORRECTION_METHODS,
        )
        from pycsamt.app.desktop.controllers.correction_controller import (
            CATALOGUE,
        )

        catalogue_fns = {
            entry["fn"] for methods in CATALOGUE.values() for entry in methods.values()
        }
        registered_fns = {m["fn"] for m in CORRECTION_METHODS.values()}
        self.assertEqual(catalogue_fns, registered_fns)

    # ── Wave F: conversational Q&A vs commands ───────────────────────────────
    def test_correction_questions_route_to_question(self):
        # Conceptual questions answer (QUESTION) rather than triggering a run.
        from pycsamt.agents.router import (
            classify_intent_offline,
        )

        for q in [
            "what is static shift?",
            "when should I rotate to strike?",
            "what does the notch filter remove?",
            "explain near-field correction",
        ]:
            self.assertEqual(classify_intent_offline(q)[0], "question", msg=q)

    def test_correction_commands_route_to_workflow(self):
        from pycsamt.agents._workflows import (
            classify_workflow,
        )
        from pycsamt.agents.router import (
            classify_intent_offline,
        )

        cases = {
            "apply ama static shift": "corr_ss_ama",
            "rotate to geoelectric strike": "corr_rotate_strike",
            "run stratagem full pipeline": "corr_strat_full",
            "correct near field": "corr_near_field",
        }
        for text, wf in cases.items():
            self.assertEqual(classify_intent_offline(text)[0], "workflow", msg=text)
            self.assertEqual(classify_workflow(text), wf, msg=text)

    def test_capability_text_lists_corrections(self):
        from pycsamt.app.agent_master.callbacks.chat import (
            _capability_text,
        )

        txt = _capability_text()
        self.assertIn("Correct & condition", txt)
        for cat in (
            "Static Shift",
            "Noise Removal",
            "Tensor Rotation",
            "Coordinates",
            "Source Effects",
            "Stratagem",
        ):
            self.assertIn(cat, txt, msg=cat)

    def test_layered_model_dataless(self):
        # layered_model needs no EDI dataset.
        import matplotlib.pyplot as plt

        r = ToolAgent().execute(
            {
                "kind": "layered_model",
                "preset": "custom",
                "resistivities": "100, 10, 500",
                "thicknesses": "300, 800",
            }
        )
        self.assertEqual(r.status, "success", msg=r.summary)
        figs = r.data.get("figures") or {}
        self.assertTrue(figs)
        for f in figs.values():
            self.assertIsInstance(f, plt.Figure)
        plt.close("all")

    def test_layered_model_preset(self):
        import matplotlib.pyplot as plt

        r = ToolAgent().execute(
            {
                "kind": "layered_model",
                "preset": "random",
                "n_layers": 5,
            }
        )
        self.assertEqual(r.status, "success", msg=r.summary)
        self.assertIn("5 layers", r.data["table_text"] + r.summary)
        plt.close("all")

    def test_circular_strike_mean(self):
        import numpy as np

        self.assertAlmostEqual(_circular_strike_mean([10, 10, 10]), 10.0, 1)
        self.assertTrue(np.isnan(_circular_strike_mean([])))

    def test_df_to_text(self):
        import pandas as pd

        df = pd.DataFrame({"a": [1.23456, 2.0], "b": ["x", "y"]})
        txt = _df_to_text(df, columns=["a", "b"], ndigits=2)
        self.assertIn("1.23", txt)
        self.assertIn("x", txt)

    def test_ll_to_utm(self):
        # London ~ UTM zone 30N; easting/northing must be finite & in range.
        import numpy as np

        e, n, z = _ll_to_utm(51.5074, -0.1278, None, "N", "WGS84")
        self.assertEqual(int(z), 30)
        self.assertTrue(np.isfinite(e) and np.isfinite(n))
        self.assertGreater(n, 5e6)  # northing for ~51°N

    def test_get_latlon(self):
        class _Obj:
            lat, lon = 12.5, -3.25

        la, lo = _get_latlon(_Obj())
        self.assertAlmostEqual(la, 12.5)
        self.assertAlmostEqual(lo, -3.25)
        # missing coords → (None, None)
        self.assertEqual(_get_latlon(object()), (None, None))

    def test_safe_filename(self):
        self.assertEqual(_safe_filename("Rho/Phi (xy)"), "Rho_Phi__xy_")
        self.assertEqual(_safe_filename("   "), "figure")

    def test_export_bundles(self):
        # every bundle resolves to known PlotAgent kinds
        from pycsamt.agents.plotting import PLOT_KINDS

        for kinds in _EXPORT_BUNDLES.values():
            for k in kinds:
                self.assertIn(k, PLOT_KINDS)

    def test_unknown_kind(self):
        self.assertEqual(
            ToolAgent().execute({"path": _DATA, "kind": "nope"}).status,
            "failed",
        )

    def test_missing_data(self):
        self.assertEqual(ToolAgent().execute({"kind": "strike"}).status, "failed")


@unittest.skipUnless(_HAS_DATA, "sample EDI data not available")
class TestToolAgent(unittest.TestCase):
    def tearDown(self):
        plt.close("all")

    def _run(self, kind, **params):
        r = ToolAgent().execute({"path": _DATA, "kind": kind, **params})
        self.assertEqual(r.status, "success", msg=r.summary)
        self.assertTrue((r.data or {}).get("table_text"))
        for fig in (r.data.get("figures") or {}).values():
            self.assertIsInstance(fig, plt.Figure)
        return r

    def test_strike(self):
        self._run("strike", method="consensus")

    def test_dimensionality(self):
        self._run("dimensionality")

    def test_validator(self):
        r = self._run("validator")
        # validator returns a table, no figure
        self.assertEqual((r.data.get("figures") or {}), {})

    # ── Wave 0/A: static-shift correction methods ────────────────────────────
    def _run_correction(self, wf_id, path=None, expect_corrected=True, **params):
        from pycsamt.agents._corrections import fn_for

        r = ToolAgent().execute(
            {
                "path": path or _DATA,
                "kind": "correction",
                "corr_wf": wf_id,
                "fn_name": fn_for(wf_id),
                **params,
            }
        )
        self.assertEqual(r.status, "success", msg=r.summary)
        if expect_corrected:
            self.assertIsNotNone((r.data or {}).get("corrected_sites"))
        else:  # diagnostic step (e.g. Stratagem QC)
            self.assertIsNone((r.data or {}).get("corrected_sites"))
        self.assertTrue((r.data or {}).get("table_text"))
        for fig in (r.data.get("figures") or {}).values():
            self.assertIsInstance(fig, plt.Figure)
        return r

    def test_correction_ama(self):
        # The generic correction kind applies a catalogue method and hands the
        # corrected Sites back through corrected_sites for the postproc modal.
        self._run_correction("corr_ss_ama", half_window=3, max_skew=6.0)

    def test_correction_ss_bilateral(self):
        # Bilateral SS is NOT covered by the legacy StaticShiftAgent.
        self._run_correction("corr_ss_bilateral", half_window=4, max_skew=6.0)

    def test_correction_ss_emap(self):
        # Hanning EMAP (Torres-Verdín) — also not in the legacy agent.
        self._run_correction(
            "corr_ss_emap", window_m=1500.0, spacing_m=200.0, comp="det"
        )

    # ── Wave B: noise-removal corrections ────────────────────────────────────
    def test_correction_notch(self):
        self._run_correction(
            "corr_notch",
            mains_hz=50.0,
            n_harm=10,
            tol_hz=0.08,
            mode="interp",
            also="both",
        )

    def test_correction_smooth_rho_phase(self):
        # boolean (check) params arrive as real bools after coercion.
        self._run_correction(
            "corr_smooth_rho_phase",
            components="offdiag",
            degree=3,
            smooth_rho=True,
            smooth_phase=True,
            blend=1.0,
            robust=True,
        )

    # ── Wave C: tensor-rotation corrections ──────────────────────────────────
    def test_correction_rotate_angle(self):
        self._run_correction("corr_rotate_angle", angle=30.0)

    def test_correction_rotate_strike(self):
        self._run_correction("corr_rotate_strike", method="swift")

    def test_correction_antisymmetrize(self):
        self._run_correction("corr_antisymmetrize", how="rms")

    # ── Wave D: coordinate corrections (corrected Sites, position figure) ─────
    def test_correction_coord_projection(self):
        self._run_correction("corr_coord_projection", azimuth=-1.0, keep_elevation=True)

    def test_correction_coord_shift(self):
        # A pure offset must change the stored coordinates.
        from pycsamt.emtools._core import ensure_sites
        from pycsamt.gis.coord_correction import (
            _get_coords_df,
        )

        r = self._run_correction(
            "corr_coord_shift",
            delta_lat=0.001,
            delta_lon=0.002,
            delta_elev=5.0,
        )
        before = _get_coords_df(ensure_sites(_DATA, recursive=True, verbose=0))
        after = _get_coords_df(r.data["corrected_sites"])
        self.assertAlmostEqual(
            float(after["lat"].iloc[0]) - float(before["lat"].iloc[0]),
            0.001,
            places=4,
        )

    def test_correction_coord_interpolate(self):
        self._run_correction("corr_coord_interpolate", fill_nan_only=True)

    # ── Wave E: Source Effects + Stratagem ───────────────────────────────────
    def test_correction_near_field(self):
        self._run_correction("corr_near_field", source_offset=500.0)

    @unittest.skipUnless(_HAS_WILLY, "Stratagem EDI directory not available")
    def test_correction_strat_qc(self):
        # QC is diagnostic: a report table, no corrected_sites to apply.
        r = self._run_correction("corr_strat_qc", path=_WILLY, expect_corrected=False)
        self.assertIn("Stratagem", r.summary)

    @unittest.skipUnless(_HAS_WILLY, "Stratagem EDI directory not available")
    def test_correction_strat_static_shift(self):
        self._run_correction(
            "corr_strat_static_shift",
            path=_WILLY,
            half_window=3,
            max_skew=6.0,
        )

    @unittest.skipUnless(_HAS_WILLY, "Stratagem EDI directory not available")
    def test_correction_strat_full(self):
        self._run_correction("corr_strat_full", path=_WILLY, mains_hz=50.0)

    def test_correction_strat_requires_directory(self):
        # Stratagem without a directory path → graceful failure, not a crash.
        from pycsamt.agents._corrections import fn_for

        r = ToolAgent().execute(
            {
                "path": "definitely/not/a/dir",
                "kind": "correction",
                "corr_wf": "corr_strat_qc",
                "fn_name": fn_for("corr_strat_qc"),
            }
        )
        self.assertEqual(r.status, "failed")

    # ── Wave C: data / IO tools ──────────────────────────────────────────────
    def test_coords(self):
        # No figure; succeeds even when stations lack coordinates.
        r = self._run("coords", datum="WGS84")
        self.assertEqual((r.data.get("figures") or {}), {})

    def test_converter_csv(self):
        import tempfile

        out = tempfile.mkdtemp(prefix="wc_conv_")
        r = self._run("converter", format="csv", output_dir=out)
        self.assertTrue(os.path.isfile(os.path.join(out, "survey_stations.csv")))
        self.assertEqual((r.data.get("figures") or {}), {})

    def test_converter_json(self):
        import tempfile

        out = tempfile.mkdtemp(prefix="wc_conv_")
        self._run("converter", format="json", output_dir=out)
        self.assertTrue(os.path.isfile(os.path.join(out, "survey_stations.json")))

    def test_batch_export(self):
        import tempfile

        out = tempfile.mkdtemp(prefix="wc_fig_")
        r = self._run(
            "batch_export",
            plots="overview",
            format="png",
            dpi=80,
            output_dir=out,
        )
        saved = [f for f in os.listdir(out) if f.endswith(".png")]
        self.assertTrue(saved, msg="no figures written")
        self.assertTrue(r.data.get("figures") or {})

    def test_elevation_offline(self):
        # Stations without coordinates → graceful, no network call.
        r = ToolAgent().execute({"path": _DATA, "kind": "elevation"})
        self.assertEqual(r.status, "success")

    def test_elevation_fetch_mocked(self):
        # Inject synthetic coords + a stubbed API so no network is touched.
        import pycsamt.agents.tooling as T
        import pycsamt.gis.utils as G

        orig_coords = T._station_coords
        orig_api = G.get_elevation_from_api
        T._station_coords = lambda sites: [
            ("A", 12.0, -3.0),
            ("B", 12.1, -3.1),
            ("C", None, None),
        ]
        G.get_elevation_from_api = lambda lat, lon, api_name=None: ([100.0, 200.0])
        try:
            r = ToolAgent().execute({"path": _DATA, "kind": "elevation"})
        finally:
            T._station_coords = orig_coords
            G.get_elevation_from_api = orig_api
        self.assertEqual(r.status, "success", msg=r.summary)
        self.assertIn("100", r.data["table_text"])
        self.assertIn("200", r.data["table_text"])

    def test_freq_editor(self):
        # Mutating tool: succeeds and hands back an edited survey for the
        # post-processing modal.
        r = ToolAgent().execute(
            {
                "path": _DATA,
                "kind": "freq_editor",
                "mode": "recover",
                "method": "composite",
                "threshold": 0.5,
            }
        )
        self.assertEqual(r.status, "success", msg=r.summary)
        self.assertIn("corrected_sites", r.data)
        self.assertIsNotNone(r.data["corrected_sites"])
        for fig in (r.data.get("figures") or {}).values():
            self.assertIsInstance(fig, plt.Figure)


if __name__ == "__main__":
    unittest.main(verbosity=2)
