# Gabbs Valley — USGS magnetotelluric survey, Nevada (2020)

Fifty-nine wideband magnetotelluric (MT) soundings collected by the U.S.
Geological Survey in July–August 2020 across the Gabbs Valley region,
Nevada, to characterize crustal control on the local geothermal system.
Full metadata (citation, processing lineage, coordinate reference, and the
shapefile attribute schema) is in `gv_survey_v2.xml` (FGDC Content Standard
for Digital Geospatial Metadata).

Only 3 of the 59 EDI files are tracked here — `gv100`, `gv130`, `gv163`,
spanning the low, middle, and high end of the station numbering — to keep
the repository small while giving `pycsamt.models.occam1d`'s worked example
(`docs/source/user_guide/models/occam1d.rst`) real field data instead of a
synthetic sounding. The full EDI set, the raw time series (`.h5`), and the
per-station quicklook images are available from the source DOI below if a
larger example is needed.

## Attribution — required if you use this data

Public-domain USGS data release. Cite:

> Peacock, J. R., Siler, D. L., Dean, B., Zielinski, L., 2021,
> Magnetotelluric Data from Gabbs Valley Region, Nevada, 2020: U.S.
> Geological Survey data release, <https://doi.org/10.5066/P9GZ9Z56>

Related publication:

> Peacock, J. R., and Siler, D. L., 2021, Bottom-up and top-down control on
> hydrothermal resources in the Great Basin: An example from Gabbs Valley,
> Nevada. *Geophysical Research Letters* (in review at time of the data
> release).

## Survey characteristics

- 59 sites, each with 3 orthogonal magnetic induction coils and 2 horizontal
  orthogonal electric dipoles (5 channels per station).
- Recorded an average of 18 hours per site (alternating 256 samples/s for
  7h50 and 4096 samples/s for 10 min), processed together for remote
  reference with BIRRP v5.2.1 (Chave and Thomson, 2004).
- EDI files rotated to geographic coordinates (`X` = geographic North),
  `STDVERS=SEG 1.0`, produced with `PROGNAME=MTpy`.
- Frequency coverage varies by station (~40–48 frequencies per sounding,
  roughly 0.0005–770 Hz); per the survey abstract, data quality is good
  over periods of 0.007–2048 s with some estimates less robust at the
  longest periods. Missing bands are marked with the EDI's own `EMPTY`
  sentinel (`1e32`) rather than omitted, which is why
  `pycsamt.models.occam1d.processing.extract_sounding` explicitly filters
  non-finite/non-positive observations before building an inversion sounding
  (see that function's docstring).
- Coordinates: WGS 84 (`gv_mt_station_locations.shp`); elevations from USGS
  3DEP, NAVD88.

## Usage

```python
from pycsamt.models.occam1d import Occam1DBatch, Occam1DConfig

config = Occam1DConfig(mode="determinant", n_layers=30, depth_max=5000.0)
batch = Occam1DBatch(
    "data/gv_data/gv_final_edi", "occam1d-inversion", config=config
).build_all()
```
