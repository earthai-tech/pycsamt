# Phase-tensor & strike / dimensionality recipe

## User intent phrases
- phase tensor, PT analysis
- strike, geoelectric strike, skew
- dimensionality, Mohr circle, Argand diagram

## Required inputs
- EDI directory / survey line
- output directory

## Preferred API
```python
from pycsamt.agents import PhaseAnalysisAgent
result = PhaseAnalysisAgent().execute({"path": edi_dir, "output_dir": output_dir})
```

## Key symbols
- `pycsamt.agents.phase_analysis.PhaseAnalysisAgent`

## Expected outputs
- phase-tensor pseudosection
- strike / dimensionality analysis (1-D / 2-D / 3-D fractions, strike angle)
- optional Mohr / Argand figures
- a phase-tensor report

## Notes
- Strongly 3-D data (high skew) limits 1-D/2-D interpretation; the analysis
  reports the dimensionality breakdown so you can judge this.
