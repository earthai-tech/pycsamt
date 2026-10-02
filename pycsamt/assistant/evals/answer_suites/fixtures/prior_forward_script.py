"""Prior artifact: 1-D MT forward response for a three-layer model."""
import numpy as np

from pycsamt.forward.em1d import MT1DForward
from pycsamt.forward.synthetic import LayeredModel

resistivity = [100.0, 10.0, 1000.0]  # ohm m; last value is the half-space
thickness = [500.0, 1000.0]  # m; finite layers only
freqs = np.logspace(-2, 2, 20)  # Hz

model = LayeredModel(resistivity=resistivity, thickness=thickness)
response = MT1DForward(freqs).run(model)
np.savetxt(
    "forward_response.csv",
    np.column_stack([response.freqs, response.rho_a, response.phase]),
    delimiter=",",
    header="freq_hz,rho_a_ohm_m,phase_deg",
)
