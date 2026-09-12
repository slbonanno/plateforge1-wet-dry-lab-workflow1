# emit

Files for external systems.

**Inputs:** `PCK`, `PLT`, `RUN` ids
**Outputs:** `XFR` transfer lists, `DOC` forms and reports

Every emitter is two layers: a neutral internal structure, then a formatter per
target system. Formatters live behind the `EMITTERS` registry.

Planned targets:

- INTEGRA VIALAB hit-picking worklist
- generic transfer CSV
- sequencing submission form
- run summary report

No emitter is written before a real example of its target format exists in
`fixtures/` and is documented in `docs/formats/`.
