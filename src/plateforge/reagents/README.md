# reagents

Reagent catalog, lot tracking, and the rules layer that reasons about confounds.

**Inputs:** reagent and lot definitions, `RUN` ids and their reagent usage
**Outputs:** `RGT`, `LOT`, `CAV` artifacts

This is not a passive log. Catalog entries carry attributes — host species,
clonality, conjugate, tag specificity, cross-reactivity — and registered rules
evaluate combinations in use to raise caveats:

- anti-tag secondary against a tagged antigen
- secondary species matching the sample source species
- substrate incompatible with media components
- lot expired or changed mid-campaign

Caveats attach to the run and surface in the analysis summary. Adding a rule is
a new file plus one `@RULES.register(...)` decorator.
