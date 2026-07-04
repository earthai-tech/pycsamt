# SPECTRA — de-identified cross-power spectra EDI examples

Two example EDI files with a raw `>=SPECTRASECT` cross-power-spectra
block, used as the sample data for the `pycsamt.emtools.spectra` gallery
example. This is a different (richer) EDI structure than the standard
Z-tensor-only files used elsewhere in the gallery (`data/AMT/WILLY_DATA/`,
`data/MT/kap03lmt_edis/`), which only carry the final impedance/tipper
estimates rather than the underlying spectra.

- `spectra01.edi` — short-period/AMT band, 51 frequencies, ~10400–1.7 Hz.
- `spectra02.edi` — broadband/long-period, 73 frequencies, ~320–0.0004 Hz.

## Provenance and anonymization

Both files are derived from real field surveys. The originals are
private data collected at a protected/sensitive location and are **not**
included in this repository or tracked by git. These two files are
sanitized copies: every location-, operator-, and hardware-identifying
field (`LAT`/`LONG`/`ELEV`, `REFLAT`/`REFLONG`/`REFELEV`, `COUNTRY`/
`STATE`, `ACQBY`/`FILEBY`/`COMPANY`/`SURVEY`/`JOB`, station/section IDs,
acquisition dates, and MTU-box/reference serial numbers) has been
replaced with generic placeholders. The instrument model name
(`HARDWARE: MTU5A`), the `>HMEAS`/`>EMEAS` channel geometry, the
`NCHAN`/`NFREQ`/`MAXBLKS`/channel-ID structure, and every numeric
`>SPECTRA` cross-power block are untouched — these carry no location or
operator information and are exactly what the gallery example needs to
demonstrate coherence, PSD, and Z/tipper-from-spectra computations.

## Usage

```python
from pycsamt.seg.spectra import Spectra

sp = Spectra.from_file("data/MT/SPECTRA/spectra01.edi")
```
