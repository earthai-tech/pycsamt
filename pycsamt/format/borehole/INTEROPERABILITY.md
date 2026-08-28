# PCBH interoperability matrix

PCBH is the lossless project exchange format. Exporters are views and must
report omissions through `PCBHExportReport` rather than imply equivalence.

| Target | Useful mapping | Expected loss or constraint |
|---|---|---|
| CSV directory | Collars, surveys, interval families, structures | Relational manifest is required to retain identity, units, and relationships |
| LAS 2.0 subset | One hole's depth-indexed logs and selected metadata | Multi-hole packaging, rich CRS, structures, dictionaries, and arbitrary extensions are reduced |
| GeoJSON | Collars and desurveyed trajectories as map features | Depth-log semantics and 3-D styling remain properties or are omitted |
| VTK XML PolyData | 3-D paths, interval segments, scalar attributes | Exchange metadata and domain dictionaries are not a complete PCBH round trip |
| glTF/GLB 2.0 | Presentation-ready tubes, colors, and scene geometry | It is a visualization asset, not a scientific borehole archive |
| PCSF | Embed or checksum-reference the unchanged PCBH document | Alignment must explicitly handle CRS and vertical datum |

Future mappings to geology-domain standards should be added only with a field
mapping, controlled-vocabulary policy, conformance fixtures, and explicit loss
report. A format name alone is not sufficient evidence of interoperability.

