# WILLY_DATA — AMT example survey

A real audio-magnetotelluric (AMT) field survey used throughout pyCSAMT as
the built-in example/demo dataset (docstrings, agent examples, unit tests,
and the `docs/source/emtools/` gallery).

## What's here

Five station lines were recorded, but **only `L18PLT` and `L22PLT` are
tracked in this repository** to keep it small. `L26PLT`, `L30PLT`, and
`L34PLT` are excluded via `.gitignore` — if you already have them locally
(e.g. from the original delivery of this dataset), they stay usable on
disk and untouched by git; they are simply never committed. If you need
those additional lines and do not already have them, contact the
corresponding author of the paper cited below.

| Line     | Tracked in git | Stations |
|----------|:--------------:|---------:|
| L18PLT   | yes            | 28       |
| L22PLT   | yes            | 25       |
| L26PLT   | no (local only)| —        |
| L30PLT   | no (local only)| —        |
| L34PLT   | no (local only)| —        |

## Citation

This survey is the field dataset behind the published AMT study of the
Baohuashan Cu-Mo porphyry deposit (Jiangsu Province). If you reuse this
data, please cite:

> Kouabena, K.A.W., Zhou, J., Chen, R., Yin, L., Cai, H., Lu, Z., Gu, J.,
> Yu, W. (2025). Enhanced prediction of deep-seated Cu-Mo porphyry
> mineralization: A comprehensive interpretation based on 2D inversion of
> audio-magnetotelluric data. *Ore Geology Reviews*, 185, 106798.
> https://doi.org/10.1016/j.oregeorev.2025.106798

See also [Kouabena2025] in `docs/source/references.rst`.

## Survey characteristics (from the EDI headers)

- Format: SEG EDI 1.0, processed with `MTPROC V1.0.7`.
- Acquired late October–November 2023, near 32°N, 119°E.
- No vertical-field (`Hz`) channel — these lines have **no tipper**. For
  gallery examples that need real tipper data, see `data/MT/` (KAP03)
  instead; `docs/source/emtools/examples/_synthetic.py` can also attach a
  synthetic tipper response on top of this survey's real station geometry
  when useful.

## Usage

Read with the public API:

```python
from pycsamt.api import read_edis

survey = read_edis("data/AMT/WILLY_DATA/L18PLT", recursive=False)
```
