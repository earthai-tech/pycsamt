# PCSM conversion examples

PCSM is the human-readable ASCII projection of the canonical HDF5 PCSF
container. Both encodings reconstruct the same `PCSFModel`; PCSM is not a
separate inversion-result schema.

From the repository root, run:

```console
python examples/pcsm_conversion/run_demo.py
```

The command converts the real Occam2D, MARE2DEM, and ModEM PCSF artifacts
created by `examples/pcsf_conversion_demo/`. It keeps the complete Occam2D
and MARE2DEM text files and short documentation excerpts for both ModEM and
MARE2DEM. Use `--full-modem` to also retain the roughly 19.7 MB plain ModEM
file and its gzip-compressed form.

It also converts `occam2d_with_bln_topo.pcsf` (`pycsamt.format.topo_source`'s
positional/UTM topo attachment — see `pcsf_conversion_demo/README.md`) to
`occam2d_with_topo.pcsm`, confirming `stations/lon`/`lat` survive the
PCSF→PCSM round trip, not just a raw `.pcsf` write.

The excerpts contain literal omission comments and are intentionally not
valid PCSM inputs. The complete Occam2D and MARE2DEM examples are parsed again
after writing as a round-trip verification.
