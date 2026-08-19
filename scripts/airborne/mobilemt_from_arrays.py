"""Build a MobileMT line from decoded arrays without assuming a file schema."""

import numpy as np

from pycsamt.airborne import NavigationTrack
from pycsamt.airborne.mobilemt import (
    build_mobilemt_dataset,
    build_mobilemt_line,
)

frequency = np.array([20.0, 100.0, 1000.0])
sample_ids = ("P001", "P002", "P003")

navigation = NavigationTrack(
    sample_ids=sample_ids,
    latitude=(5.0000, 5.0005, 5.0010),
    longitude=(-3.0000, -3.0005, -3.0010),
    platform_elevation=(220.0, 221.0, 219.5),
    terrain_elevation=(120.0, 121.5, 120.5),
)

# samples x frequency x (Hx, Hy, Hz) x (Ex, Ey)
admittance = np.zeros((3, 3, 3, 2), dtype=complex)
admittance[:, :, 0, 0] = 1.0 + 0.1j
admittance[:, :, 1, 1] = 0.8 + 0.2j
admittance[:, :, 2, 0] = 0.2 + 0.05j

line = build_mobilemt_line(
    "L001",
    navigation,
    admittance,
    frequency=frequency,
    record_mask=(True, False, True),
)

dataset = build_mobilemt_dataset("demo_mobilemt", [line])

print(dataset.n_lines)
print(dataset.n_samples)
print(dataset.n_records)
print(dataset.transfer_function_names)
