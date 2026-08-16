# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.geology.structural — StructuralMeasurement,
LinearMeasurement, FaultTrace, StructuralModel."""

from __future__ import annotations

import pytest

from pycsamt.geology.structural import (
    FaultTrace,
    LinearMeasurement,
    StructuralMeasurement,
    StructuralModel,
)

# ─────────────────────────────────────────────────────────────────────────────
# StructuralMeasurement
# ─────────────────────────────────────────────────────────────────────────────


def test_structural_measurement_basic_construction():
    m = StructuralMeasurement(
        x=500.0, kind="bedding", strike_deg=45.0, dip_deg=30.0,
        dip_direction_deg=135.0,
    )
    assert m.x == 500.0
    assert m.kind == "bedding"
    assert m.dip_azimuth_ok is True
    assert m.station is None
    assert m.confidence == 1.0


def test_structural_measurement_accepts_minus_90_side():
    # strike 135, dip_direction 45 == strike - 90
    m = StructuralMeasurement(
        x=0, kind="joint", strike_deg=135.0, dip_deg=50.0,
        dip_direction_deg=45.0,
    )
    assert m.dip_azimuth_ok is True


def test_structural_measurement_normalizes_strike_and_dip_direction():
    m = StructuralMeasurement(
        x=0, kind="bedding", strike_deg=405.0, dip_deg=10.0,
        dip_direction_deg=495.0,
    )
    assert m.strike_deg == 45.0
    assert m.dip_direction_deg == 135.0


def test_structural_measurement_rejects_inconsistent_dip_direction():
    with pytest.raises(ValueError, match="transposed"):
        StructuralMeasurement(
            x=0, kind="bedding", strike_deg=45.0, dip_deg=30.0,
            dip_direction_deg=45.0,
        )


def test_structural_measurement_rejects_out_of_range_dip():
    with pytest.raises(ValueError, match="dip_deg"):
        StructuralMeasurement(
            x=0, kind="bedding", strike_deg=0.0, dip_deg=95.0,
            dip_direction_deg=90.0,
        )
    with pytest.raises(ValueError, match="dip_deg"):
        StructuralMeasurement(
            x=0, kind="bedding", strike_deg=0.0, dip_deg=-1.0,
            dip_direction_deg=90.0,
        )


def test_structural_measurement_rejects_out_of_range_confidence():
    with pytest.raises(ValueError, match="confidence"):
        StructuralMeasurement(
            x=0, kind="bedding", strike_deg=0.0, dip_deg=10.0,
            dip_direction_deg=90.0, confidence=1.5,
        )


def test_structural_measurement_tolerance_is_configurable():
    # 30 deg off from strike+90 (=135); default tolerance (20) rejects,
    # a wider explicit tolerance accepts.
    with pytest.raises(ValueError):
        StructuralMeasurement(
            x=0, kind="bedding", strike_deg=45.0, dip_deg=30.0,
            dip_direction_deg=165.0,
        )
    m = StructuralMeasurement(
        x=0, kind="bedding", strike_deg=45.0, dip_deg=30.0,
        dip_direction_deg=165.0, dip_direction_tolerance_deg=35.0,
    )
    assert m.dip_azimuth_ok is True


def test_structural_measurement_from_right_hand_rule():
    m = StructuralMeasurement.from_right_hand_rule(
        500.0, "foliation", dip_direction_deg=200.0, dip_deg=40.0,
    )
    assert m.strike_deg == 110.0
    assert m.dip_direction_deg == 200.0
    assert m.dip_deg == 40.0
    assert m.dip_azimuth_ok is True


def test_structural_measurement_repr():
    m = StructuralMeasurement(
        x=500.0, kind="bedding", strike_deg=45.0, dip_deg=30.0,
        dip_direction_deg=135.0,
    )
    assert repr(m) == "StructuralMeasurement(x=500.0 m, 'bedding', 45/30->135)"


def test_structural_measurement_clone_revalidates():
    m = StructuralMeasurement(
        x=500.0, kind="bedding", strike_deg=45.0, dip_deg=30.0,
        dip_direction_deg=135.0,
    )
    with pytest.raises(ValueError, match="transposed"):
        m.clone(dip_direction_deg=0.0)
    # original untouched
    assert m.dip_direction_deg == 135.0
    good = m.clone(dip_deg=40.0)
    assert good.dip_deg == 40.0
    assert good is not m


def test_structural_measurement_update_revalidates():
    m = StructuralMeasurement(
        x=500.0, kind="bedding", strike_deg=45.0, dip_deg=30.0,
        dip_direction_deg=135.0,
    )
    with pytest.raises(ValueError, match="dip_deg"):
        m.update(dip_deg=200.0)


# ─────────────────────────────────────────────────────────────────────────────
# LinearMeasurement
# ─────────────────────────────────────────────────────────────────────────────


def test_linear_measurement_basic_construction():
    lm = LinearMeasurement(x=500.0, kind="fold_axis", trend_deg=210.0, plunge_deg=15.0)
    assert lm.x == 500.0
    assert lm.trend_deg == 210.0
    assert lm.plunge_deg == 15.0


def test_linear_measurement_normalizes_trend():
    lm = LinearMeasurement(x=0, kind="lineation", trend_deg=370.0, plunge_deg=10.0)
    assert lm.trend_deg == 10.0


def test_linear_measurement_rejects_out_of_range_plunge():
    with pytest.raises(ValueError, match="plunge_deg"):
        LinearMeasurement(x=0, kind="lineation", trend_deg=0.0, plunge_deg=91.0)


def test_linear_measurement_rejects_out_of_range_confidence():
    with pytest.raises(ValueError, match="confidence"):
        LinearMeasurement(
            x=0, kind="lineation", trend_deg=0.0, plunge_deg=10.0,
            confidence=-0.1,
        )


def test_linear_measurement_repr():
    lm = LinearMeasurement(x=500.0, kind="fold_axis", trend_deg=210.0, plunge_deg=15.0)
    assert repr(lm) == "LinearMeasurement(x=500.0 m, 'fold_axis', 210/15)"


def test_linear_measurement_clone_revalidates():
    lm = LinearMeasurement(x=0, kind="fold_axis", trend_deg=200.0, plunge_deg=10.0)
    with pytest.raises(ValueError, match="plunge_deg"):
        lm.clone(plunge_deg=120.0)


# ─────────────────────────────────────────────────────────────────────────────
# FaultTrace
# ─────────────────────────────────────────────────────────────────────────────


def test_fault_trace_basic_construction():
    ft = FaultTrace(x=500.0, dip_deg=70.0, downthrown_side="right", throw_m=12.0)
    assert ft.x == 500.0
    assert ft.sense == "unknown"
    assert ft.throw_m == 12.0


def test_fault_trace_rejects_bad_downthrown_side():
    with pytest.raises(ValueError, match="downthrown_side"):
        FaultTrace(x=0, dip_deg=60.0, downthrown_side="up")


def test_fault_trace_rejects_bad_sense():
    with pytest.raises(ValueError, match="sense"):
        FaultTrace(x=0, dip_deg=60.0, downthrown_side="left", sense="oblique")


def test_fault_trace_rejects_negative_throw():
    with pytest.raises(ValueError, match="throw_m"):
        FaultTrace(x=0, dip_deg=60.0, downthrown_side="left", throw_m=-5.0)


def test_fault_trace_rejects_out_of_range_dip():
    with pytest.raises(ValueError, match="dip_deg"):
        FaultTrace(x=0, dip_deg=91.0, downthrown_side="left")


def test_fault_trace_normalizes_optional_strike():
    ft = FaultTrace(x=0, dip_deg=60.0, downthrown_side="left", strike_deg=400.0)
    assert ft.strike_deg == 40.0


def test_fault_trace_strike_none_by_default():
    ft = FaultTrace(x=0, dip_deg=60.0, downthrown_side="left")
    assert ft.strike_deg is None


def test_fault_trace_repr_unknown_throw():
    ft = FaultTrace(x=500.0, dip_deg=70.0, downthrown_side="right")
    assert repr(ft) == "FaultTrace(x=500.0 m, dip=70 deg, down=right, throw=?)"


def test_fault_trace_repr_known_throw():
    ft = FaultTrace(x=500.0, dip_deg=70.0, downthrown_side="right", throw_m=12.0)
    assert repr(ft) == "FaultTrace(x=500.0 m, dip=70 deg, down=right, throw=12.0 m)"


def test_fault_trace_clone_revalidates():
    ft = FaultTrace(x=0, dip_deg=60.0, downthrown_side="left")
    with pytest.raises(ValueError, match="downthrown_side"):
        ft.clone(downthrown_side="sideways")


def test_fault_trace_update_returns_self():
    ft = FaultTrace(x=0, dip_deg=60.0, downthrown_side="left")
    out = ft.update(throw_m=5.0)
    assert out is ft
    assert ft.throw_m == 5.0


# ─────────────────────────────────────────────────────────────────────────────
# StructuralModel — construction, mutation
# ─────────────────────────────────────────────────────────────────────────────


def test_structural_model_default_empty():
    model = StructuralModel()
    assert model.planar == []
    assert model.linear == []
    assert model.faults == []
    assert model.metadata == {}


def test_structural_model_repr():
    model = StructuralModel(
        faults=[FaultTrace(x=500.0, dip_deg=70.0, downthrown_side="right")],
    )
    assert repr(model) == "StructuralModel(0 planar, 0 linear, 1 faults)"


def test_structural_model_add_methods():
    model = StructuralModel()
    model.add_planar(
        StructuralMeasurement(
            x=0, kind="bedding", strike_deg=0.0, dip_deg=10.0,
            dip_direction_deg=90.0,
        )
    )
    model.add_linear(
        LinearMeasurement(x=0, kind="fold_axis", trend_deg=0.0, plunge_deg=5.0)
    )
    model.add_fault(FaultTrace(x=0, dip_deg=60.0, downthrown_side="left"))
    assert len(model.planar) == 1
    assert len(model.linear) == 1
    assert len(model.faults) == 1


def test_structural_model_metadata_mixin():
    model = StructuralModel()
    assert model.metadata_dict() == {}
    model.update_metadata(source="field_notebook_2024")
    assert model.metadata == {"source": "field_notebook_2024"}


# ─────────────────────────────────────────────────────────────────────────────
# StructuralModel — queries
# ─────────────────────────────────────────────────────────────────────────────


def _demo_model():
    return StructuralModel(
        faults=[
            FaultTrace(x=200.0, dip_deg=60.0, downthrown_side="left"),
            FaultTrace(x=800.0, dip_deg=45.0, downthrown_side="right"),
        ],
        planar=[
            StructuralMeasurement(
                x=100.0, kind="bedding", strike_deg=0.0, dip_deg=10.0,
                dip_direction_deg=90.0,
            ),
            StructuralMeasurement(
                x=900.0, kind="bedding", strike_deg=0.0, dip_deg=15.0,
                dip_direction_deg=90.0,
            ),
        ],
    )


def test_nearest_returns_closest_item():
    model = _demo_model()
    nearest = model.nearest(250.0, kind="faults")
    assert nearest.x == 200.0


def test_nearest_respects_max_distance():
    model = _demo_model()
    assert model.nearest(250.0, kind="faults", max_distance=10.0) is None
    assert model.nearest(250.0, kind="faults", max_distance=100.0) is not None


def test_nearest_on_empty_list_returns_none():
    assert StructuralModel().nearest(0.0, kind="faults") is None


def test_nearest_invalid_kind_raises():
    with pytest.raises(ValueError, match="kind must be one of"):
        _demo_model().nearest(0.0, kind="bogus")


def test_nearest_planar_kind():
    model = _demo_model()
    nearest = model.nearest(150.0, kind="planar")
    assert nearest.x == 100.0


def test_within_filters_by_range():
    model = _demo_model()
    sub = model.within(0.0, 500.0)
    assert len(sub.faults) == 1
    assert sub.faults[0].x == 200.0
    assert len(sub.planar) == 1
    assert sub.planar[0].x == 100.0


def test_within_returns_independent_model():
    model = _demo_model()
    sub = model.within(0.0, 500.0)
    sub.add_fault(FaultTrace(x=100.0, dip_deg=10.0, downthrown_side="left"))
    assert len(model.faults) == 2  # original unaffected


# ─────────────────────────────────────────────────────────────────────────────
# to_dict
# ─────────────────────────────────────────────────────────────────────────────


def test_structural_model_to_dict_shape():
    model = _demo_model()
    d = model.to_dict()
    assert set(d) == {"planar", "linear", "faults"}
    assert len(d["faults"]) == 2
    assert d["faults"][0]["x"] == 200.0
    assert len(d["planar"]) == 2


def test_leaf_objects_to_dict_uses_pycsamt_object():
    m = StructuralMeasurement(
        x=500.0, kind="bedding", strike_deg=45.0, dip_deg=30.0,
        dip_direction_deg=135.0,
    )
    d = m.to_dict()
    assert d["x"] == 500.0
    assert d["kind"] == "bedding"
    assert d["strike_deg"] == 45.0


# ─────────────────────────────────────────────────────────────────────────────
# CSV I/O
# ─────────────────────────────────────────────────────────────────────────────


def test_from_csv_all_three(tmp_path):
    planar_csv = tmp_path / "planar.csv"
    planar_csv.write_text(
        "x,kind,strike_deg,dip_deg,dip_direction_deg,confidence,notes\n"
        "500,bedding,45,30,135,0.9,near outcrop A\n"
        "750,foliation,10,60,100,,\n"
    )
    linear_csv = tmp_path / "linear.csv"
    linear_csv.write_text(
        "x,kind,trend_deg,plunge_deg\n500,fold_axis,210,15\n"
    )
    faults_csv = tmp_path / "faults.csv"
    faults_csv.write_text(
        "x,dip_deg,downthrown_side,sense,throw_m,evidence\n"
        "500,70,right,normal,12,resistivity offset\n"
    )

    model = StructuralModel.from_csv(
        planar_path=planar_csv, linear_path=linear_csv, faults_path=faults_csv,
    )
    assert len(model.planar) == 2
    assert model.planar[0].confidence == 0.9
    assert model.planar[0].notes == "near outcrop A"
    assert model.planar[1].confidence == 1.0  # default, blank column
    assert model.planar[1].notes == ""

    assert len(model.linear) == 1
    assert model.linear[0].trend_deg == 210.0

    assert len(model.faults) == 1
    assert model.faults[0].sense == "normal"
    assert model.faults[0].throw_m == 12.0
    assert model.faults[0].evidence == "resistivity offset"


def test_from_csv_optional_paths_yield_empty_lists(tmp_path):
    faults_csv = tmp_path / "faults.csv"
    faults_csv.write_text(
        "x,dip_deg,downthrown_side\n500,70,right\n"
    )
    model = StructuralModel.from_csv(faults_path=faults_csv)
    assert model.planar == []
    assert model.linear == []
    assert len(model.faults) == 1


def test_from_csv_faults_missing_optional_columns(tmp_path):
    faults_csv = tmp_path / "faults.csv"
    faults_csv.write_text("x,dip_deg,downthrown_side\n500,70,right\n")
    model = StructuralModel.from_csv(faults_path=faults_csv)
    ft = model.faults[0]
    assert ft.sense == "unknown"
    assert ft.throw_m is None
    assert ft.strike_deg is None
    assert ft.evidence == ""


def test_from_csv_raises_on_missing_header(tmp_path):
    empty_csv = tmp_path / "empty.csv"
    empty_csv.write_text("")
    with pytest.raises(ValueError, match="no header"):
        StructuralModel.from_csv(faults_path=empty_csv)
