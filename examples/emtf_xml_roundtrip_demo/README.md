# EDI &lt;-&gt; EMTF-XML round-trip demo

Converts three real field stations from `data/gv_data/gv_final_edi/`
(`gv100`, `gv101`, `gv102` by default — non-orthogonal channel geometry,
tipper present) to EMTF-XML and back to EDI, through **both** public
entry points:

- **Path A — `Site`/`Sites`**: `Sites([...])` -> `Sites.write_xml()` ->
  `to_sites()` (rediscovering the `*.xml` directory, the same coercion
  `ensure_sites()` uses) -> `Sites.write()`.
- **Path B — the raw `EMTF` document**: `EMTF.from_edi()` ->
  `EMTF.write_xml()` -> `EMTF.from_xml()` ->
  `pycsamt.emtf.converters.edi.write_edi()`.

Run from the repository root:

```bash
python examples/emtf_xml_roundtrip_demo/run_demo.py
```

Convert different (or more) stations:

```bash
python examples/emtf_xml_roundtrip_demo/run_demo.py --station gv103 --station gv104
```

## What gets written

```
output/
  xml/                        Path A: Sites.write_xml() output
  edi_roundtrip_via_site/     Path A: round-tripped EDI (Sites.write())
  xml_via_emtf/               Path B: EMTF.write_xml() output
  edi_roundtrip_via_emtf/     Path B: round-tripped EDI (write_edi())
  demo-summary.json           machine-readable verification report
```

Open any of the `.xml`/`.edi` files directly to inspect them.

## What gets verified

For each station, the script compares periods, impedance, tipper, and
variance (max absolute difference) across five pairs:

1. Path A's round-tripped `EMTF` vs. the original.
2. Path B's round-tripped `EMTF` vs. the original.
3. Path A's result vs. Path B's result (do both APIs agree?).
4. Path A's **on-disk** written EDI, re-parsed from disk, vs. the
   original.
5. Path B's **on-disk** written EDI, re-parsed from disk, vs. the
   original.

The last two re-read the actual files from disk rather than trusting
the in-memory objects that wrote them — that distinction matters: it's
exactly what caught a real bug in `Sites.write()` while building this
demo (it looked for a nonexistent `EDIFile.to_file()` method and
silently fell back to writing a one-line placeholder instead of real
EDI content; fixed to use `write_edi()`, the same writer Path B uses).

All five comparisons come back exact (`max|dZ|=0.0` etc.) on the
bundled `gv100`/`gv101`/`gv102` fixtures — full float64 precision is
preserved end-to-end, and there's no missing data in these three files
to trigger a `NaN`/data-loss warning. Any `DataLossWarning` raised
during a conversion is captured and listed in
`demo-summary.json`'s `data_loss_warnings_via_*` fields (empty for
these fixtures).
