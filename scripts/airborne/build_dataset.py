"""Build a format-neutral airborne EM dataset without a native parser."""

import numpy as np

from pycsamt.airborne import AirborneEMDataset, AirborneEMLine, NavigationTrack
from pycsamt.emtf import EMTF, TransferFunction

nav = NavigationTrack(
    sample_ids=("P001", "P002"),
    latitude=(5.0, 5.001),
    longitude=(-3.0, -3.001),
    terrain_elevation=(100.0, 102.0),
    platform_elevation=(180.0, 181.0),
)
line = AirborneEMLine(line_id="L001", navigation=nav)

periods = 1.0 / np.array([100.0, 1000.0])
y = np.zeros((2, 3, 2), dtype=complex)
tf = TransferFunction(
    name="admittance",
    data=y,
    input_channels=("Ex", "Ey"),
    output_channels=("Hx", "Hy", "Hz"),
    periods=periods,
)
doc = EMTF(periods=periods)
doc.add_transfer_function(tf)
line.add_emtf("P001", doc)

dataset = AirborneEMDataset(name="example", lines={"L001": line})
print(dataset)
print(dataset.transfer_function_names)
