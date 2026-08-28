# PCBH governance and extension registry

PCBH is an open pyCSAMT exchange contract. Changes are proposed through a
public issue and pull request, with schema, specification, migration notes,
fixtures, and tests in the same change. Maintainers review scientific meaning,
interoperability, security limits, and backward compatibility.

## Compatibility policy

- Readers accept compatible patch releases of the implemented `0.1` line and
  reject unknown major/minor contracts rather than guessing.
- Patch changes may clarify prose, fix a schema defect that contradicts the
  specification, or add non-normative examples. They cannot rename fields or
  change their meaning.
- A minor pre-1.0 version may add optional fields. A breaking change increments
  the minor version while PCBH remains pre-1.0; after 1.0 it increments major.
- Writers preserve unknown JSON-compatible, namespaced extension values.
- A version is frozen only after its schema, fixtures, parser, exporters, app
  rendering, and supported-Python CI matrix pass together.

## Extension registration

Extension keys use `owner:name`, for example `acme:core_recovery`. The prefix
identifies an organization or project that controls the meaning; `pcbh` is
reserved for extensions maintained by this project. Registration is optional
for private exchange but recommended for public interoperability.

A registry proposal must provide the key, owner/contact, purpose, JSON shape,
units and null rules, one valid example, validation rules, security/resource
considerations, and license. Registered keys are added below by pull request.
Conflicting semantics receive a different key; an existing key is never
silently redefined.

| Prefix/key | Owner | Status | Meaning |
|---|---|---|---|
| `pcbh:*` | pyCSAMT maintainers | reserved | First-party experimental fields |

Deprecation requires a replacement or rationale, reader support for at least
one published compatibility line, and a migration note. Security reports may
follow the repository security policy; malformed input must remain bounded by
the parser's byte, nesting, borehole, interval, and CSV-row limits.

