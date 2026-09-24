# Integra ASSIST / VIAFLO worklists — nothing here yet

Rule 6: no emitter until a real example is in this directory. A worklist
emitted against a guessed format is a wasted plate of reagent.

## What to drop in

One saved file per program shape we need, unmodified, straight out of the
software:

1. a **plate-to-plate transfer** (the 500 µL supernatant harvest)
2. a **reagent addition** — one source trough to all 96 wells (bead addition)
3. a **variable-volume dilution** — different volume per well (BLI normalisation)

## What to note alongside each

- exact column headers, including case and separators
- whether wells are `A1` or `A01` (this project is always `A01`, rule 3, and
  converts at the boundary)
- how the source labware is addressed (position number? named deck slot?)
- volume units, and the minimum and maximum the tip head accepts
- whether one file holds one step or a whole program

Once a file is here, the emitter goes in `emit/` as one registered handler
(rule 5) and its test reads this fixture.
