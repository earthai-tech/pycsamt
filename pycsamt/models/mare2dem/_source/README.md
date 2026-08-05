# MARE2DEM source — populated on demand, not vendored

This directory is intentionally **empty in the pycsamt git repository** (this
`README.md` is the only tracked file). MARE2DEM is a third-party solver
hosted on Bitbucket, distributed under its own license — pycsamt does not
vendor its Fortran source the way it vendors ModEM/Occam2D under the sibling
`models/modem/_source/` and `models/occam2d/_source/` directories.

Instead, `pycsamt.models.mare2dem.source.SourceManager` downloads it here on
demand:

```python
from pycsamt.models.mare2dem import SourceManager

sm = SourceManager(verbose=1)
sm.download()   # clones/extracts the real ~5-6 MB source tree into this dir
sm.build()      # compiles it (Linux/macOS/WSL only)
```

or via the CLI: `pycsamt build mare2dem`.

Everything `download()` writes here (the real `.f90`/`.c`/`Makefile` sources,
and any nested `.git` directory if the `git` download method is used) is
covered by `.gitignore` and must never be committed. It is also excluded
from the PyPI sdist/wheel — see the `mare2dem/_source` note in `MANIFEST.in`
— so `pip install pycsamt` stays small; `SourceManager` downloads directly
into a user-data directory instead when the installed package tree is
read-only (see `resolve_source_dir()`'s docstring for the full fallback
order).
