# Static-shift correction recipe

## User intent phrases
- static shift, static-shift, galvanic shift/distortion
- correct / remove static shift
- shift correction, AMA correction
- "write static shift code", "process EDI line for static shift"

## Required inputs
- EDI directory, or a survey-line name resolved via the project registry
- output directory
- method: `ama` (default), `loess`, `refmedian`, `bilateral`
- sort_by: `lon` (default), `lat`, `name`
- half_window (spatial AMA window, default 3)
- optional period band `pband = [T_min, T_max]` in seconds

## Preferred high-level API (run the workflow)
Use `StaticShiftAgent` when the user wants it executed.

```python
from pycsamt.agents import StaticShiftAgent

agent = StaticShiftAgent(method="ama", half_window=3)
result = agent.execute({
    "path": edi_dir,
    "output_dir": output_dir,
})
corrected = result["corrected_sites"]
```

## Preferred low-level API (standalone code)
Use `pycsamt.emtools.ss` when the user asks for a plain script.

```python
from pycsamt.emtools._core import ensure_sites
from pycsamt.emtools.ss import (
    estimate_ss_ama,
    correct_ss_ama,
    plot_ss_summary,
    plot_ss_1d_curves,
)

sites = ensure_sites(edi_dir, recursive=True, verbose=1)
ss_table = estimate_ss_ama(sites, half_window=3, sort_by="lon")
sites_corr = correct_ss_ama(sites, half_window=3, sort_by="lon")
```

## Key symbols
- `pycsamt.agents.static_shift.StaticShiftAgent`
- `pycsamt.emtools.ss.estimate_ss_ama`
- `pycsamt.emtools.ss.correct_ss_ama`
- `pycsamt.emtools.ss.apply_ss_factors`
- `pycsamt.emtools.ss.plot_ss_summary`
- `pycsamt.emtools.ss.plot_ss_1d_curves`

## Expected outputs
- static-shift factor table (`station`, `delta_log10_rho`, `fac_rho`, `fac_z`, `n_used`)
- corrected `Sites` object
- summary figure (before / after / delta)
- per-station 1-D apparent-resistivity curves
- optional corrected EDI files

## Notes / gotchas
- AMA assumes ~1-D/2-D structure; on strongly 3-D data (high skew) it may
  decline (empty factor table, before == after). That is expected, not a bug.
- Static-shift factors must be finite & positive; non-finite factors are
  ignored so they never destroy the impedance.
