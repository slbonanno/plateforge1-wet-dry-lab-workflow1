# Octet BLI exports — nothing here yet

Octet instruments read 96-well sample plates, so the plate-in/plate-out shape
in `decisions/0019-plate-lineage.md` holds. The export format still needs a
real file before a reader is written (rule 6, Q7).

## What to drop in

- one raw export from a **quantitation** run (ProteinA tips, IgG standard
  curve) — this is the one that feeds the normalisation step
- the plate map the run was set up with, so the mapping from sensor position
  back to well is checked rather than assumed

## What to note alongside it

- file type actually produced (`.frd`, `.csv`, a folder per run?)
- where the concentration result lives, and its units
- how the standard curve is reported, and whether R² travels with it
- whether well identifiers are `A1` or `A01`
