# CSAMT — Tongkeng example survey

A real controlled-source audio-frequency magnetotelluric (CSAMT) field
survey used throughout pyCSAMT as the built-in CSAMT example/demo dataset
(docstrings, `pycsamt.seg` true-file tests, and CSAMT theory/user-guide
examples).

## What's here

One profile line, 10 stations spaced 50 m apart (0–450 m), acquired
2020-04-23 by Zhejiang University (`ACQBY=ZJU UNIV`) near Zixing, Hunan
Province, China (~26.05°N, 110.49°E). Station `DATAID`s are `S00`–`S09`;
the filename encodes distance along the profile in metres.

| File         | Station | Distance (m) |
|--------------|:-------:|--------------:|
| `csa000.edi` | S00     | 0             |
| `csa050.edi` | S01     | 50            |
| `csa100.edi` | S02     | 100           |
| `csa150.edi` | S03     | 150           |
| `csa200.edi` | S04     | 200           |
| `csa250.edi` | S05     | 250           |
| `csa300.edi` | S06     | 300           |
| `csa350.edi` | S07     | 350           |
| `csa400.edi` | S08     | 400           |
| `csa450.edi` | S09     | 450           |

Each station carries 17 frequencies (8196.7 Hz down to 0.125 Hz).

## Citation

This survey is the field dataset behind the published CSAMT groundwater
study for the Tongkeng area. If you reuse this data, please cite:

> Kouadio, K.L., Xu, Y., Liu, C., Boukhalfa, Z. (2020). Two-dimensional
> inversion of CSAMT data and three-dimensional geological mapping for
> groundwater exploration in Tongkeng Area, Hunan Province, China.
> *Journal of Applied Geophysics*, 104204.
> https://doi.org/10.1016/j.jappgeo.2020.104204

See also [Kouadio2020] in `docs/source/references.rst`.

## Usage

Read with the public API:

```python
from pycsamt.api import read_edis

survey = read_edis("data/CSAMT", recursive=False)
```
