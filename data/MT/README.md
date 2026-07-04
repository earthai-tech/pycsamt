# KAP03 — SAMTEX Kaapvaal Craton MT profile

Twenty-six long-period magnetotelluric (MT) EDI files from the **KAP03**
profile of **SAMTEX** (Southern African Magnetotelluric Experiment), a
NE–SW transect across the Kaapvaal Craton, South Africa, from the SW end
(`kap103`) to the NE end (`kap175`), roughly 60 km long.

Downloaded from <https://www.mtnet.info/data/kap03/kap03.html>. Unlike the
AMT lines in `data/AMT/WILLY_DATA/`, this survey has a real vertical-field
(tipper) channel, which makes it the reference dataset for the
`pycsamt.emtools.tf` (induction-vector) gallery example.

## Attribution — required if you use this data

Per the terms stated on the source page: *"If you download any of these
time series data, or the EDI MT impedance tensor estimates, and present or
publish them, you are obliged to acknowledge the SAMTEX Consortium as
their source."* Cite:

> Jones, A.G., R.L. Evans, M.R. Muller, M.P. Hamilton, M.P. Miensopust,
> X. Garcia, P. Cole, T. Ngwisanyi, D. Hutchins, C.J.S Fourie, H. Jelsma,
> T. Aravanis, W. Pettit, S. Webb, J. Wasborg, and The SAMTEX Team, 2009.
> Area selection for diamonds using magnetotellurics: Examples from
> southern Africa. *Lithos*, **112S**, 83–92.

## Survey characteristics

- 26 sites recorded October/November 2003 with GSC LIMS systems and
  Phoenix LRMT clones, ~1 month per site at 5-second sampling.
- EDI files BIRRP-processed, dated June 2004; `STDVERS=2.3`,
  `PROGVERS="JONES 2.3"`.
- Frequency band (from the EDI headers): ~5.9×10⁻⁵–4×10⁻² Hz (periods of
  ~25 s–4.7 h), 20 frequencies per site.

## Usage

```python
from pycsamt.api import read_edis

survey = read_edis("data/MT/kap03lmt_edis", recursive=False)
```
