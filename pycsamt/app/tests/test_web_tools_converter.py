from __future__ import annotations

from pycsamt.app.web.callbacks.tools import (
    _batch_export_tool_figures,
    _converter_options_from_values,
    _converter_result_view,
    _render_stored_tool_result,
)


def test_converter_options_include_avg_topography_controls():
    options = _converter_options_from_values(
        output_dir="/tmp/out",
        stn_path="/tmp/K1.stn",
        utm_zone="49N",
        epsg="32649",
        convert_coords=True,
        e_labels="EX,EY",
        h_labels="HX,HY",
        suffix="",
        spectra_flags=[],
        freq_order="ascending",
        freq_tol="0.001",
        core_flags=["compute_z", "compute_rho_phi", "load_session"],
        station_name="K1",
    )

    assert options["output_dir"] == "/tmp/out"
    assert options["stn_path"] == "/tmp/K1.stn"
    assert options["utm_zone"] == "49N"
    assert options["epsg"] == "32649"
    assert options["convert_stn_coords"] is True
    assert options["freq_order"] == "ascending"
    assert options["freq_tol"] == 0.001
    assert options["compute_z"] is True
    assert options["compute_rho_phi"] is True
    assert options["station_name"] == "K1"


def test_converter_options_include_spectra_flags():
    options = _converter_options_from_values(
        output_dir="",
        stn_path="",
        utm_zone="",
        epsg="",
        convert_coords=False,
        e_labels="Ex,Ey",
        h_labels="Hx,Hy",
        suffix="_IMP",
        spectra_flags=["estimate_error", "use_remote", "skip_errors"],
        freq_order="descending",
        freq_tol=None,
        core_flags=[],
        station_name="",
    )

    assert options["e_labels"] == "Ex,Ey"
    assert options["h_labels"] == "Hx,Hy"
    assert options["station_suffix"] == "_IMP"
    assert options["estimate_error"] is True
    assert options["use_remote"] is True
    assert options["skip_errors"] is True
    assert options["compute_z"] is False


def test_converter_result_view_reports_coordinate_counts():
    view = _converter_result_view(
        stats={
            "n_total": 2,
            "n_failures": 1,
            "rows": [
                {
                    "station": "S01",
                    "n_freqs": 3,
                    "f_min": 1.0,
                    "f_max": 100.0,
                    "lat": 10.0,
                    "lon": 20.0,
                    "elev": 50.0,
                    "has_Z": True,
                    "has_tipper": False,
                },
                {
                    "station": "S02",
                    "n_freqs": 0,
                    "f_min": float("nan"),
                    "f_max": float("nan"),
                    "lat": float("nan"),
                    "lon": float("nan"),
                    "elev": float("nan"),
                    "has_Z": False,
                    "has_tipper": False,
                },
            ],
        },
        failures=[{"source": "bad.edi", "error": "not spectra"}],
        output_dir="/tmp/out",
        loaded=True,
    )

    assert view is not None


def test_stored_strike_result_renders_in_light_mode():
    payload = {
        "tool": "strike",
        "method_label": "Consensus strike",
        "scope_label": "L18",
        "n_stations": 2,
        "median_strike": 42.5,
        "median_iqr": 8.0,
        "records": [
            {"station": "S01", "line": "L18", "ang_axial": 40.0},
            {"station": "S02", "line": "L18", "ang_axial": 45.0},
        ],
        "table": [
            {"Station": "S01", "Line": "L18", "Strike (°)": 40.0},
            {"Station": "S02", "Line": "L18", "Strike (°)": 45.0},
        ],
        "columns": ["Station", "Line", "Strike (°)"],
        "page_size": 10,
    }

    view = _render_stored_tool_result("strike", {"strike": payload}, "light")

    assert view is not None


def test_batch_export_tool_saves_station_map_and_manifest(tmp_path):
    result = _batch_export_tool_figures(
        [],
        dest=tmp_path,
        fmt="png",
        dpi=72,
        items=["map"],
        flags={"manifest"},
        store_data={
            "station_records": [
                {
                    "ID": "S01",
                    "Line": "L18",
                    "Latitude": 10.0,
                    "Longitude": 20.0,
                },
                {
                    "ID": "S02",
                    "Line": "L18",
                    "Latitude": 10.1,
                    "Longitude": 20.2,
                },
            ],
        },
        active_lines_store={"active": ["L18"]},
    )

    assert result["saved"] == 1
    assert result["failed"] == 0
    assert (tmp_path / "station_map.png").exists()
    assert (tmp_path / "pycsamt_batch_export_manifest.csv").exists()
