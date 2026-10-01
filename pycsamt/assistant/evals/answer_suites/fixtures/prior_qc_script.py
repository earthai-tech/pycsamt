"""Prior artifact: QC table for data/3edis saved as results/qc/qc_results.csv."""
import os

from pycsamt.emtools import ensure_sites
from pycsamt.emtools.qc import build_qc_table

output_dir = "results/qc"
os.makedirs(output_dir, exist_ok=True)

sites = ensure_sites("data/3edis", verbose=0)
qc_table = build_qc_table(sites)
qc_table.to_csv(os.path.join(output_dir, "qc_results.csv"), index=False)
print(qc_table.head())
