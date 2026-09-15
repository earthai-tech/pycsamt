# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from pycsamt.constants import MU_0, PI
from pycsamt.exceptions import ProcessingError
from pycsamt.zonge.proc_utils import (
    ama,
    flma,
    get_reference_frequency,
    get_skew,
    get_strike,
    interpolate_to_log_space,
    prepare_strike_frame,
    smooth_rho_from_phase,
    tma,
)


@pytest.fixture(scope="module")
def base_df() -> pd.DataFrame:
    stations = np.repeat([100, 200, 300], 5)
    freqs = np.tile([1024, 512, 256, 128, 64], 3)
    data = {
        "station": stations,
        "freq": freqs,
        "comp": ["ExHy"] * 15,
        "rho": np.linspace(50, 200, 15),
        "phase": np.linspace(1000, 1400, 15),
    }
    df = pd.DataFrame(data)
    omega = 2 * PI * df["freq"]
    z_mag = np.sqrt(df["rho"] * omega * MU_0)
    phase_rad = df["phase"] * 1e-3
    df["z"] = z_mag * np.exp(1j * phase_rad)
    return df


# ─────────────────────────────────────────────────────────────────────────
# tma
# ─────────────────────────────────────────────────────────────────────────


def test_tma_raises_for_even_window_size(base_df):
    with pytest.raises(ValueError, match="odd integer"):
        tma(base_df["rho"], window_size=4)


def test_tma_raises_when_string_profile_without_data():
    with pytest.raises(ValueError, match="DataFrame must be provided"):
        tma("rho", window_size=3)


def test_tma_array_input_returns_ndarray(base_df):
    out = tma(base_df["rho"].to_numpy(), window_size=3)
    assert isinstance(out, np.ndarray)
    assert len(out) == len(base_df)


# ─────────────────────────────────────────────────────────────────────────
# flma
# ─────────────────────────────────────────────────────────────────────────


def test_flma_raises_when_string_profile_without_data():
    with pytest.raises(ValueError, match="DataFrame must be provided"):
        flma("z", "station", 100.0)


def test_flma_dataframe_input_adds_column(base_df):
    out = flma("z", "station", 100.0, data=base_df)
    assert isinstance(out, pd.DataFrame)
    assert "z_flma" in out.columns


def test_flma_array_input_returns_ndarray(base_df):
    out = flma(
        base_df["z"].to_numpy(),
        base_df["station"].to_numpy(),
        100.0,
    )
    assert isinstance(out, np.ndarray)


def test_flma_all_nan_window_produces_nan():
    z = pd.Series([np.nan, np.nan, np.nan])
    stn = pd.Series([0.0, 10.0, 20.0])
    out = flma(z, stn, dipole_length=1.0, filter_width_dipoles=0.1)
    assert out.isna().all()


# ─────────────────────────────────────────────────────────────────────────
# ama
# ─────────────────────────────────────────────────────────────────────────


def test_ama_raises_when_string_profile_without_data():
    with pytest.raises(ValueError, match="'data' when"):
        ama("z", "station", 100.0, 1024.0)


def test_ama_dataframe_input_adds_column(base_df):
    out = ama("z", "station", 100.0, 1024.0, data=base_df, iterations=1)
    assert isinstance(out, pd.DataFrame)
    assert "z_ama" in out.columns


def test_ama_series_input_returns_series(base_df):
    out = ama(base_df["z"], base_df["station"], 100.0, 1024.0, iterations=1)
    assert isinstance(out, pd.Series)


def test_ama_all_nan_window_produces_nan():
    z = pd.Series([np.nan, np.nan, np.nan])
    stn = pd.Series([0.0, 10.0, 20.0])
    out = ama(z, stn, dipole_length=1.0, frequency=1024.0, iterations=1)
    assert out.isna().all()


# ─────────────────────────────────────────────────────────────────────────
# interpolate_to_log_space
# ─────────────────────────────────────────────────────────────────────────


def test_interpolate_raises_when_columns_missing():
    df = pd.DataFrame({"station": [1], "freq": [1.0]})
    with pytest.raises(ProcessingError):
        interpolate_to_log_space(df)


def test_interpolate_skips_soundings_with_fewer_than_two_points(base_df):
    df = base_df[base_df["station"] == 100].iloc[:1]
    out = interpolate_to_log_space(df, num_points=5)
    assert out.empty


# ─────────────────────────────────────────────────────────────────────────
# smooth_rho_from_phase
# ─────────────────────────────────────────────────────────────────────────


def test_smooth_rho_raises_when_columns_missing():
    df = pd.DataFrame({"station": [1], "freq": [1.0]})
    with pytest.raises(ProcessingError):
        smooth_rho_from_phase(df)


def test_smooth_rho_skips_soundings_with_fewer_than_five_points(base_df):
    df = base_df[base_df["station"] == 100].iloc[:3]
    out = smooth_rho_from_phase(df)
    assert out["rho_smoothed"].isna().all()


# ─────────────────────────────────────────────────────────────────────────
# get_reference_frequency
# ─────────────────────────────────────────────────────────────────────────


def test_get_reference_frequency_raises_when_columns_missing():
    df = pd.DataFrame({"station": [1], "freq": [1.0]})
    with pytest.raises(ProcessingError):
        get_reference_frequency(df)


def test_get_reference_frequency_falls_back_when_no_clean_data():
    df = pd.DataFrame(
        {"station": [1, 2], "freq": [10.0, 20.0], "pc_rho": [50.0, 80.0]}
    )
    with pytest.warns(UserWarning, match="Falling back"):
        ref = get_reference_frequency(df, qc_threshold=20.0)
    assert ref == 20.0


# ─────────────────────────────────────────────────────────────────────────
# get_strike / get_skew
# ─────────────────────────────────────────────────────────────────────────


def test_get_strike_raises_when_columns_missing():
    df = pd.DataFrame({"station": [1], "freq": [1.0]})
    with pytest.raises(ProcessingError):
        get_strike(df)


def test_get_strike_fills_missing_components_with_zero():
    df = pd.DataFrame(
        {
            "station": [1, 1],
            "freq": [10.0, 10.0],
            "comp": ["ExHx", "ExHy"],
            "z": [1 + 1j, 2 + 2j],
        }
    )
    out = get_strike(df)
    assert len(out) == 1
    assert pd.notna(out["strike_angle"].iloc[0])


def test_get_skew_raises_when_columns_missing():
    df = pd.DataFrame({"station": [1], "freq": [1.0]})
    with pytest.raises(ProcessingError):
        get_skew(df)


def test_get_skew_fills_missing_components_with_zero():
    df = pd.DataFrame(
        {
            "station": [1, 1],
            "freq": [10.0, 10.0],
            "comp": ["ExHx", "ExHy"],
            "z": [1 + 1j, 2 + 2j],
        }
    )
    out = get_skew(df)
    assert len(out) == 1


# ─────────────────────────────────────────────────────────────────────────
# prepare_strike_frame
# ─────────────────────────────────────────────────────────────────────────


def test_prepare_strike_frame_raises_when_neither_source_given():
    with pytest.raises(ProcessingError, match="Provide 'z_frame'"):
        prepare_strike_frame()


def test_prepare_strike_frame_from_z_frame_prefer_z():
    z_frame = pd.DataFrame(
        {
            "station": [1, 2],
            "freq": [10.0, 20.0],
            "comp": ["ExHy", "ExHy"],
            "z": [1 + 1j, 2 + 2j],
        }
    )
    out = prepare_strike_frame(z_frame=z_frame)
    assert list(out.columns) == ["station", "freq", "comp", "z"]
    assert np.iscomplexobj(out["z"])
    assert len(out) == 2


def test_prepare_strike_frame_z_frame_missing_columns_falls_back_to_df():
    z_frame = pd.DataFrame({"station": [1], "freq": [10.0]})  # missing comp/z
    df = pd.DataFrame(
        {
            "station": [1],
            "freq": [10.0],
            "comp": ["ExHy"],
            "rho": [100.0],
            "phase": [45.0],
        }
    )
    out = prepare_strike_frame(z_frame=z_frame, df=df)
    assert "z" in out.columns
    assert np.iscomplexobj(out["z"])


def test_prepare_strike_frame_derives_z_from_df_deg_phase():
    df = pd.DataFrame(
        {
            "station": [1],
            "freq": [10.0],
            "comp": ["ExHy"],
            "rho": [100.0],
            "phase": [45.0],
        }
    )
    out = prepare_strike_frame(df=df, phase_unit="deg")
    assert np.isclose(np.angle(out["z"].iloc[0]), np.deg2rad(45.0))


def test_prepare_strike_frame_derives_z_from_df_mrad_phase():
    df = pd.DataFrame(
        {
            "station": [1],
            "freq": [10.0],
            "comp": ["ExHy"],
            "rho": [100.0],
            "phase": [450.0],  # mrad
        }
    )
    out = prepare_strike_frame(df=df, phase_unit="mrad")
    assert np.isclose(np.angle(out["z"].iloc[0]), 450.0 / 1000.0)


def test_prepare_strike_frame_derives_z_from_df_rad_phase():
    df = pd.DataFrame(
        {
            "station": [1],
            "freq": [10.0],
            "comp": ["ExHy"],
            "rho": [100.0],
            "phase": [0.5],  # already radians
        }
    )
    out = prepare_strike_frame(df=df, phase_unit="rad")
    assert np.isclose(np.angle(out["z"].iloc[0]), 0.5)


def test_prepare_strike_frame_auto_phase_detects_radians():
    df = pd.DataFrame(
        {
            "station": [1],
            "freq": [10.0],
            "comp": ["ExHy"],
            "rho": [100.0],
            "phase": [0.5],  # |0.5| <= 1.5*pi -> treated as radians
        }
    )
    out = prepare_strike_frame(df=df, phase_unit="auto")
    assert np.isclose(np.angle(out["z"].iloc[0]), 0.5)


def test_prepare_strike_frame_auto_phase_detects_degrees():
    df = pd.DataFrame(
        {
            "station": [1],
            "freq": [10.0],
            "comp": ["ExHy"],
            "rho": [100.0],
            "phase": [90.0],  # > 1.5*pi, <= 180 -> degrees
        }
    )
    out = prepare_strike_frame(df=df, phase_unit="auto")
    assert np.isclose(np.angle(out["z"].iloc[0]), np.deg2rad(90.0))


def test_prepare_strike_frame_auto_phase_detects_mrad():
    df = pd.DataFrame(
        {
            "station": [1],
            "freq": [10.0],
            "comp": ["ExHy"],
            "rho": [100.0],
            "phase": [900.0],  # > 180 -> mrad
        }
    )
    out = prepare_strike_frame(df=df, phase_unit="auto")
    assert np.isclose(np.angle(out["z"].iloc[0]), 900.0 / 1000.0)


def test_prepare_strike_frame_df_raises_when_missing_columns():
    df = pd.DataFrame({"station": [1], "freq": [10.0]})
    with pytest.raises(ProcessingError, match="Missing columns"):
        prepare_strike_frame(df=df)


def test_prepare_strike_frame_df_raises_when_empty_after_dropna():
    df = pd.DataFrame(
        {
            "station": [1],
            "freq": [10.0],
            "comp": ["ExHy"],
            "rho": [np.nan],
            "phase": [45.0],
        }
    )
    with pytest.raises(ProcessingError, match="No data after NA drop"):
        prepare_strike_frame(df=df)


def test_prepare_strike_frame_filters_components():
    df = pd.DataFrame(
        {
            "station": [1, 1],
            "freq": [10.0, 10.0],
            "comp": ["ExHy", "ExHx"],
            "rho": [100.0, 100.0],
            "phase": [45.0, 45.0],
        }
    )
    out = prepare_strike_frame(df=df, components=["ExHy"])
    assert set(out["comp"]) == {"ExHy"}


def test_prepare_strike_frame_z_frame_only_no_df_no_prefer_z():
    z_frame = pd.DataFrame(
        {
            "station": [2, 1],
            "freq": [20.0, 10.0],
            "comp": ["ExHy", "ExHy"],
            "z": [2 + 2j, 1 + 1j],
        }
    )
    out = prepare_strike_frame(z_frame=z_frame, prefer="df")
    # df is None, so falls through to the z_frame fallback branch
    assert list(out.columns) == ["station", "freq", "comp", "z"]
    assert out["station"].tolist() == [1, 2]  # ensure_sorted applied


def test_prepare_strike_frame_z_col_non_complex_dtype_is_coerced():
    z_frame = pd.DataFrame(
        {
            "station": [1],
            "freq": [10.0],
            "comp": ["ExHy"],
            "z": [1.5],  # plain float, not complex
        }
    )
    out = prepare_strike_frame(z_frame=z_frame)
    assert np.iscomplexobj(out["z"])
    assert out["z"].iloc[0] == 1.5 + 0j


def test_prepare_strike_frame_z_frame_fallback_non_complex_dtype_is_coerced():
    z_frame = pd.DataFrame(
        {
            "station": [1],
            "freq": [10.0],
            "comp": ["ExHy"],
            "z": [1.5],  # plain float, not complex
        }
    )
    out = prepare_strike_frame(z_frame=z_frame, prefer="df")
    assert np.iscomplexobj(out["z"])
    assert out["z"].iloc[0] == 1.5 + 0j


def test_prepare_strike_frame_z_col_already_complex_dtype():
    z_frame = pd.DataFrame(
        {
            "station": [1],
            "freq": [10.0],
            "comp": ["ExHy"],
            "z": np.array([1 + 1j], dtype=complex),
        }
    )
    out = prepare_strike_frame(z_frame=z_frame)
    assert np.iscomplexobj(out["z"])


def test_prepare_strike_frame_copy_false_does_not_copy():
    z_frame = pd.DataFrame(
        {
            "station": [1],
            "freq": [10.0],
            "comp": ["ExHy"],
            "z": [1 + 1j],
        }
    )
    out = prepare_strike_frame(z_frame=z_frame, copy=False)
    assert not out.empty
