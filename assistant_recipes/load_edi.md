# Load EDI data recipe

## User intent phrases
- load EDI / load the data / read EDI files
- open my survey, ingest EDIs, build Sites

## Required inputs
- EDI directory or file, or a survey-line name (project registry)

## Preferred API
Use `ensure_sites` (the canonical loader) or `MTLoaderAgent`.

```python
from pycsamt.emtools._core import ensure_sites
sites = ensure_sites(edi_dir, recursive=True, on_dup="replace", verbose=1)
print(f"Loaded {len(list(sites))} stations")
```

```python
from pycsamt.agents import MTLoaderAgent
result = MTLoaderAgent().execute({"path": edi_dir})
```

## Key symbols
- `pycsamt.emtools._core.ensure_sites`
- `pycsamt.agents.loader.MTLoaderAgent`
- `pycsamt.site.Sites`

## Expected outputs
- a `Sites` collection of stations with impedance tensors
- a per-station QC scan (MTLoaderAgent)

## Notes
- `recursive=True` searches sub-directories.
- A station needs valid impedance (Z) data; files with no Z are skipped.
