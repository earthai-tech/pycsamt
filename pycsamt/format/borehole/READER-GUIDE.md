# PCBH 0.1 independent reader guide

A PCBH file is UTF-8 JSON, conventionally named `*.pcbh.json`. Detect it from
the root `pcbh_version`, not its filename. Validate first against
`https://pycsamt.org/schemas/pcbh/0.1/schema.json`, then enforce the semantic
rules in `SPEC.md`; JSON Schema alone cannot detect every interval, trajectory,
or cross-reference error.

Minimum reading sequence:

1. Reject oversized input and excessive nesting before object construction.
2. Require version `0.1.x`; reject a version whose major/minor is unknown.
3. Read the document CRS, units, and conventions before any coordinates.
4. Resolve each collar `(x, y, z)` in that CRS. Measured depth (`md`) starts at
   the collar and increases down the borehole; elevation follows `z_positive`.
5. For a vertical trajectory, derive points from collar to total depth. For a
   surveyed trajectory, interpolate/desurvey ordered stations using its stated
   method, inclination reference, azimuth direction, and north reference.
6. Treat intervals as half-open `[from_md, to_md)` except their final endpoint.
   Do not invent material for gaps. Resolve codes through the relevant
   dictionary and retain `data_nature`.
7. Preserve unknown namespaced `extensions`; never interpret their names as
   core fields.

The public minimal example is
`https://pycsamt.org/examples/pcbh/0.1/minimal-vertical.pcbh.json`.
The Python reference reader is `pycsamt.format.borehole.read_pcbh`.
Implementations should report errors with a JSON-style path so users can locate
the bad value.

For 3-D display, compute trajectories in source CRS, transform the document as
a whole when alignment is needed, retain scientific radius separately from a
minimum visible display radius, and render observed logs distinctly from
interpreted or modeled data. Never use a missing vertical datum as evidence
that elevations already align.
