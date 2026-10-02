"""Prior artifact: estimate AMA factors, export them, then apply correction."""
import os

from pycsamt.emtools import ensure_sites
from pycsamt.emtools.ss import correct_ss_ama, estimate_ss_ama

output_dir = "results/static_shift"
os.makedirs(output_dir, exist_ok=True)

sites = ensure_sites("data/3edis", verbose=0)
factors = estimate_ss_ama(sites)
factors.to_csv(os.path.join(output_dir, "ss_factors.csv"))
if len(factors) == 0:
    raise SystemExit("AMA returned no usable factors; correction skipped.")
corrected = correct_ss_ama(sites, inplace=False)
print(f"Corrected {len(list(corrected))} stations")
