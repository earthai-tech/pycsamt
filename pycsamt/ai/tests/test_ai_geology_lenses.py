"""Contracts for :mod:`pycsamt.ai.geology.lenses`."""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.ai.geology import (
    ElectricalLayer,
    EllipsoidalLens,
    GeologyGrid,
    LensGeology,
    generate_layered_geology,
    insert_lenses,
)


def _base_2d(nx=20, nz=12):
    grid = GeologyGrid.regular_2d(nx=nx, nz=nz, dx_m=100, dz_m=50)
    return generate_layered_geology(grid, [ElectricalLayer("earth", 100)], [], seed=0)


def test_2d_lens_geometry_rotation_weights_and_serialization():
    grid = _base_2d().grid
    lens = EllipsoidalLens(
        "conductor",
        1000,
        300,
        500,
        120,
        5,
        dip_deg=25,
        transition_fraction=0.2,
    )
    assert lens.dimension == 2
    lens.validate_grid(grid)
    radius = lens.normalized_radius(grid)
    weight = lens.blend_weight(grid)
    assert radius.shape == weight.shape == grid.shape
    assert np.any(radius <= 1)
    assert np.all((weight >= 0) & (weight <= 1))
    assert np.any((weight > 0) & (weight < 1))
    assert EllipsoidalLens.from_dict(lens.to_dict()) == lens


def test_insert_sharp_lens_changes_only_body_and_keeps_base_immutable():
    base = _base_2d()
    lens = EllipsoidalLens("conductor", 1000, 300, 400, 100, 5)
    model = insert_lenses(base, [lens])
    mask = model.lens_mask("conductor")
    assert np.any(mask)
    np.testing.assert_allclose(model.resistivity_ohm_m[mask], 5)
    np.testing.assert_allclose(model.resistivity_ohm_m[~mask], 100)
    np.testing.assert_array_equal(base.resistivity_ohm_m, 100)
    assert np.max(model.overlap_count) == 1
    assert not model.resistivity_ohm_m.flags.writeable


def test_smooth_transition_blends_in_log_resistivity_space():
    base = _base_2d(nx=40, nz=20)
    lens = EllipsoidalLens("smooth", 2000, 500, 800, 300, 1, transition_fraction=0.5)
    model = insert_lenses(base, [lens])
    values = model.resistivity_ohm_m[model.lens_mask(0)]
    assert np.min(values) == pytest.approx(1)
    assert np.any((values > 1) & (values < 100))
    assert np.max(values) <= 100


def test_overlap_error_and_precedence_policies():
    base = _base_2d()
    conductive = EllipsoidalLens("conductive", 900, 300, 400, 150, 5)
    resistive = EllipsoidalLens("resistive", 1100, 300, 400, 150, 1000)
    with pytest.raises(ValueError, match="overlap"):
        insert_lenses(base, [conductive, resistive], conflict_policy="error")

    first = insert_lenses(base, [conductive, resistive], conflict_policy="first")
    last = insert_lenses(base, [conductive, resistive], conflict_policy="last")
    overlap = first.overlap_count > 1
    assert np.any(overlap)
    np.testing.assert_allclose(first.resistivity_ohm_m[overlap], 5)
    np.testing.assert_allclose(last.resistivity_ohm_m[overlap], 1000)
    conductive_result = insert_lenses(
        base,
        [resistive, conductive],
        conflict_policy="most_conductive",
    )
    resistive_result = insert_lenses(
        base,
        [conductive, resistive],
        conflict_policy="most_resistive",
    )
    np.testing.assert_allclose(conductive_result.resistivity_ohm_m[overlap], 5)
    np.testing.assert_allclose(resistive_result.resistivity_ohm_m[overlap], 1000)


def test_3d_rotated_dipping_ellipsoid():
    grid = GeologyGrid.regular_3d(nx=12, ny=10, nz=8, dx_m=100, dy_m=100, dz_m=50)
    base = generate_layered_geology(grid, [ElectricalLayer("earth", 100)], [], seed=0)
    lens = EllipsoidalLens(
        "body",
        600,
        200,
        350,
        120,
        10,
        center_y_m=500,
        radius_y_m=250,
        azimuth_deg=30,
        dip_deg=15,
    )
    model = insert_lenses(base, [lens])
    assert model.resistivity_ohm_m.shape == (8, 10, 12)
    assert np.any(model.lens_mask(0))


def test_lens_npz_roundtrip_is_self_contained(tmp_path):
    base = _base_2d()
    model = insert_lenses(
        base,
        [EllipsoidalLens("body", 1000, 300, 400, 100, 10)],
    )
    restored = LensGeology.from_npz(model.to_npz(tmp_path / "lenses.npz"))
    assert restored.model_hash == model.model_hash
    assert restored.base.model_hash == base.model_hash
    np.testing.assert_array_equal(restored.lens_index, model.lens_index)


def test_lens_validation_rejects_dimension_nonintersection_and_bad_parameters():
    base = _base_2d()
    with pytest.raises(ValueError, match="supplied together"):
        EllipsoidalLens("bad", 0, 1, 2, 1, 10, center_y_m=0)
    with pytest.raises(ValueError, match="does not intersect"):
        insert_lenses(
            base,
            [EllipsoidalLens("outside", 1e6, 1e6, 10, 10, 10)],
        )
    with pytest.raises(ValueError, match="dimensions differ"):
        insert_lenses(
            base,
            [
                EllipsoidalLens(
                    "3d", 1000, 300, 100, 100, 10, center_y_m=0, radius_y_m=100
                )
            ],
        )
    with pytest.raises(ValueError, match="azimuth_deg"):
        insert_lenses(
            base,
            [EllipsoidalLens("bad azimuth", 1000, 300, 100, 100, 10, azimuth_deg=10)],
        )


def test_constructor_rejects_tampered_resistivity_and_overlap_maps():
    base = _base_2d()
    lens = EllipsoidalLens("body", 1000, 300, 400, 100, 10)
    model = insert_lenses(base, [lens])
    bad_resistivity = np.array(model.resistivity_ohm_m, copy=True)
    bad_resistivity[0, 0] = 99
    with pytest.raises(ValueError, match="inconsistent"):
        LensGeology(
            base,
            (lens,),
            bad_resistivity,
            model.lens_index,
            model.overlap_count,
        )
    bad_overlap = np.array(model.overlap_count, copy=True)
    bad_overlap[0, 0] += 1
    with pytest.raises(ValueError, match="overlap_count is inconsistent"):
        LensGeology(
            base,
            (lens,),
            model.resistivity_ohm_m,
            model.lens_index,
            bad_overlap,
        )
