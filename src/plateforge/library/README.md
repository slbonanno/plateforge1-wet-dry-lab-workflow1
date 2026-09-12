# library

Sequence acquisition and in-silico construct generation.

**Inputs:** external sequence resources (OAS), sampling criteria
**Outputs:** `SEQ` sequence records, `LIB` sequence sets, `FMT` reformatted constructs

Planned scope:

- interface with OAS: query, download, cache, parse
- classify and store a large pool of sequences in the `library` store
- diversity-aware sampling: draw N sequences spanning V-gene families and
  CDRH3 length/identity ranges that resemble a real discovery campaign
- reformat onto scFv, VHH, VH-only, CDRH3-only, and rescaffolded frameworks
- generate synthetic data tables by sampling the pool

Downstream modules receive `LIB` or `FMT` ids, never objects from this module.
