# Broken Hill — high-resolution 3-D magnetotelluric survey, NSW, Australia (2024)

Twenty-one "ultra-wide-band" magnetotelluric (MT) soundings collected by
the University of Sydney in June–July 2024 around the world-class
**Broken Hill Pb-Zn-Ag deposit**, on the south-eastern margin of the
Curnamona Province, New South Wales. Average station spacing ~2.5 km.

The survey and its 3-D ModEM resistivity model are published in:

> AlQahtani, Y., Özaydın, S., Chatzaras, V., Rey, P. F., & Passos, T.
> (2026). Why does the Broken Hill deposit sit in resistive crust?
> Magnetotelluric evidence for metamorphic decoupling of a world-class
> mineral system. *Journal of Geophysical Research: Solid Earth*, **131**,
> e2026JB035666. <https://doi.org/10.1029/2026JB035666>

The paper's central result: the crust around Broken Hill is
**predominantly resistive** (> 1000 Ω·m), with only a few discrete,
shallow (< 3 km) conductive anomalies (C1–C3) that sit directly on the
sulphide-rich Broken Hill Group (chiefly the Hores Gneiss). There is no
robust vertical electrical connection between those shallow conductors and
the deeper (> 25 km) Eastern Nackara Arc–Broken Hill Conductor
(ENAC-BHC) — interpreted as a source-to-sink pathway that was electrically
and structurally *decoupled* during granulite-facies metamorphism.

This dataset gives pyCSAMT a real, tipper-bearing 3-D broadband MT survey
with a matching published inversion, useful as example data for the
`pycsamt.emtools` (phase tensor, induction vectors) and
`pycsamt.models.modem` workflows.

## Attribution — required if you use this data

The MT data and 3-D model are an **open-access CC-BY** data release. If you
use them, cite **both** the paper (above) **and** the data release:

> AlQahtani, Y., Özaydın, S., Chatzaras, V., Rey, P. F., & Passos, T.
> (2026). Open data for the article titled "Why does the Broken Hill
> deposit sit in resistive crust? Magnetotelluric evidence for metamorphic
> decoupling of a world-class mineral system". *Zenodo*.
> <https://doi.org/10.5281/zenodo.21272106>

### Which Zenodo version

| Version | DOI | Published | Title |
|---|---|---|---|
| concept (always latest) | `10.5281/zenodo.21091923` | — | — |
| **v2** (this folder) | `10.5281/zenodo.21272106` | 2026-07-09 | "Why does the Broken Hill deposit sit in resistive crust? …" |
| v1 | `10.5281/zenodo.21091924` | 2026-07-01 | "Should electrically resistive regions be avoided in magnetotelluric mineral exploration? …" |

The article's own reference list cites the **v1** DOI (`…21091924`). v2 is
a **title/metadata-only** re-release — the `MT_Data_Models.zip` payload is
**byte-for-byte identical** between v1 and v2 (6,339,884 bytes, MD5
`9d5f7910da4e5dd8d3fc34376d9c3668`), and identical to the files in this
folder. Nothing in the data or the model changed; prefer the concept DOI,
or cite v2 to match the published article title.

The authors acknowledge the **Wilyakali/Wiljagali people**, Traditional
Custodians of the land on which the survey was conducted, and the
landowners of the Rupee, Clevedale, Nine Mile, K Tank, Donsandel and
Gum Paddock Stations and the Living Desert State Park. Acquisition was
supported by ARC grants DP22010070 and LP190100146 and the University of
Sydney.

The inversion used **ModEM** (Kelbert et al., 2014, *Computers &
Geosciences* 66, 40–53; <https://doi.org/10.5281/zenodo.21760502>); the
EDI files were processed with **mtpy-v2** (Kirkby et al., 2019, *JOSS*
4(37), 1358).

## What's here

| Path | Tracked | Contents |
|------|:-------:|----------|
| `edis/` | yes, all 21 | Phoenix EMpower EDI exports, `BH_1` … `BH_21`. Each carries a full impedance tensor **and** a vertical-field (tipper) transfer function. Two naming variants: `BH_N_imp.edi` and `BH_N_imp_rev.edi` (a revised re-export — see below). 96–104 frequencies per site, ~10 kHz down to ~7×10⁻⁴ Hz. |
| `final-models/BH_31.dat` | yes | ModEM data file for the published inversion — full impedance + full vertical components, 24 periods (333 s – 0.0033 s), 21 sites, with the real per-datum error floors. |
| `final-models/BH_31_NLCG_030.dat` | yes | Forward response of the recovered model at NLCG iteration 30 (error column set to ModEM's `1E13` "response" sentinel). |
| `final-models/BH_31_NLCG_030.res` | yes | Normalized residuals (observed − predicted) for the same iteration. |
| `final-models/BH_31_NLCG_030.rho` | **no** (11 MB) | Recovered resistivity model, ModEM WS format, 129 × 120 × 55 cells, `LOGE`. Get it from the Zenodo DOI above. |
| `final-models/BH_31_NLCG_030.prm` | **no** (11 MB) | Companion parameter/prior model in the same grid and format. Zenodo. |
| `JGR Solid Earth - 2026 - AlQahtani … .pdf` | **no** | The published article (open access; not redistributed here). Read it at the paper DOI. |
| `naser-comments.jpg` | **no** | Screenshot of a LinkedIn post by Dr. Naser Meqbel about an independent ModEM re-inversion of this dataset — kept locally as a pointer, not part of the release. |

## Survey characteristics

- **21 stations**, recorded 18 June – 3 July 2024, one sounding at a time,
  ~9 h to ~1 day per site (see each EDI's `>INFO DURATION`).
- Phoenix **MTU-5C** receivers, **MTC-155** induction coils, 100 m
  orthogonal electric dipoles, instrument azimuth 0°.
- Processing: Phoenix **EMpower v2.22.0.1**, 50 Hz comb filter, robust
  multiple-coherence stacking; EDI files written by mtpy-v2 (`FILEBY=MTpy`,
  `FILEDATE=2024/09/19`).
- Usable band, per the paper: **~0.008–1000 Hz**. The EDI headers carry a
  wider nominal grid (down to ~7×10⁻⁴ Hz); the longest periods are the
  least robust.
- Dimensionality: phase-tensor skew is high at high frequencies and the
  geoelectric strike varies spatially — the data are **3-D**, which is why
  the published model is a 3-D ModEM inversion.

## Coordinates and reference frame

- `LAT`/`LON` in the EDI `>HEAD` are geographic **WGS84**, spanning roughly
  **31.86–32.03° S, 141.44–141.63° E** (near the town of Broken Hill).
  `ELEV` is in metres, 239–322 m.
- The EDI impedances are in the **acquisition frame**: `ZROT = 0` at every
  frequency and EMpower reports `COORDINATE_SYSTEM=Geomagnetic North`
  (`>INFO` lists a station declination of 9–18°). Rotate to a common
  geographic frame before joint analysis.
- `final-models/BH_31.dat` is already the **inversion-ready** form:
  ModEM local `X` (north) / `Y` (east) metres about the model origin
  `-31.95556, 41.53481` (ModEM subtracts 100° from longitude), geographic
  orientation.

## A pyCSAMT bug found while adding this dataset

Every EDI here has `DATAID=None` in `>HEAD` — EMpower writes that literal
placeholder when the station id field is left unset (the real name is only
in `>INFO STATION NAME` and the filename). `pycsamt.api.read_edis` keyed
its collection on `EDIFile.station`, which returned the raw `"None"`
string, so **all 21 soundings collapsed onto one station** and the survey
reported a single site.

Fixed in `pycsamt/seg/edi.py` (`EDIFile.station`) and
`pycsamt/seg/cbase.py` (`CoreParser._fast_station`): placeholder DATAID
values (`None`, `NULL`, `NA`, `-`, …) are now treated as missing and the
filename stem is used instead. Regression test in
`pycsamt/seg/tests/test_seg_cbase.py`
(`test_coreparser_nullish_dataid_falls_back_to_stem`).

## Usage

```python
from pycsamt.api import read_edis

survey = read_edis("data/MT/broken-hill/edis", recursive=False)
# -> APISurvey with 21 sites (BH_1_imp, BH_2_imp_rev, ...)
```

```python
from pycsamt.models.modem import ModEmData

data = ModEmData.read("data/MT/broken-hill/final-models/BH_31.dat")
# data.n_sites == 21, data.n_periods == 24  (impedance + tipper)
```

```python
from pycsamt.models.modem import InversionResult
from pycsamt.models.modem.plot import PlotDepthMap

r = InversionResult("data/MT/broken-hill/final-models")
fig = PlotDepthMap(
    r,
    depths={"(a) 1 km": 1000, "(b) 2 km": 2000},
    origin_lat=-31.95556, origin_lon=141.53481,
    mask_outside_hull=True, contours=[1000],
).plot()
```

```python
from pycsamt.emtools import plot_phase_tensor_map_grid

fig = plot_phase_tensor_map_grid(
    "data/MT/broken-hill/edis",
    frequencies=[30, 3, 0.3, 0.03],
    c_by="skew",
    tipper_convention="parkinson",
)
```

## Tutorial

`docs/source/tutorials/model_broken_hill_mt_3d.rst` is a complete
walkthrough on this dataset — load, QC, dimensionality, ModEM input, and
the published inversion model read back and interpreted. Its figures are
regenerated (and its code verified) by
`docs/scripts/generate_tutorial_broken_hill.py`.
