# Quality control (QC) recipe

## User intent phrases
- quality control, QC, check data quality
- flag bad stations, SNR, dead band, clean data

## Required inputs
- EDI directory / survey line
- output directory

## Preferred API
```python
from pycsamt.agents import DataQCAgent
result = DataQCAgent().execute({"path": edi_dir, "output_dir": output_dir})
```

Low-level helpers:
```python
from pycsamt.emtools._core import ensure_sites
from pycsamt.emtools.qc import build_qc_table, qc_flags, plot_frequency_confidence_psection

sites = ensure_sites(edi_dir)
qc_table = build_qc_table(sites)
flags = qc_flags(sites, min_frac_ok=0.6, min_snr_med=2.0)
```

## Key symbols
- `pycsamt.agents.qc.DataQCAgent`
- `pycsamt.emtools.qc.build_qc_table`
- `pycsamt.emtools.qc.qc_flags`
- `pycsamt.emtools.qc.plot_frequency_confidence_psection`

## Expected outputs
- per-station QC scores / flags
- a frequency-confidence pseudosection figure
- a QC report
