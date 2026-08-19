"""Inspect and structurally QC an already decoded airborne dataset."""

from __future__ import annotations

import numpy as np

from pycsamt.airborne import (
    AirborneEMDataset,
    AirborneEMLine,
    NavigationTrack,
    assess_airborne_qc,
    inspect_airborne,
)
from pycsamt.airborne.ztem import build_ztem_emtf

navigation = NavigationTrack(
    sample_ids=("P001", "P002", "P003"),
    latitude=(5.0, 5.001, 5.002),
    longitude=(-3.0, -3.001, -3.002),
)
line = AirborneEMLine(
    line_id="L001",
    navigation=navigation,
    attrs={"technology": "ZTEM"},
)
for sample_id in ("P001", "P003"):
    document = build_ztem_emtf(
        np.ones((2, 2), dtype=complex),
        frequency=(30.0, 90.0),
    )
    line.add_emtf(sample_id, document)

survey = AirborneEMDataset(
    name="example",
    lines={"L001": line},
    attrs={"technology": "ZTEM"},
)

print(inspect_airborne(survey).to_dict(max_depth=2))
print(assess_airborne_qc(survey).to_dict(max_depth=3))
