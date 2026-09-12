# External formats

Every external system we read from or write to gets a note in this folder and
at least one real, unmodified example file in `fixtures/`.

A note should record:

- where the file came from (instrument model, software version, who exported it)
- what the columns or blocks mean, including anything surprising
- what we verified by hand versus what we assumed
- the fixture filename it corresponds to

## Rule

No parser and no emitter is written before a real example file exists here.
Guessed formats are the main way this project can quietly produce wrong files.

## Status

| Format | Direction | Fixture | Verified |
|---|---|---|---|
| SpectraMax plate export | in | — | no |
| BioTek plate export | in | — | no |
| INTEGRA VIALAB hit-picking worklist | out | — | no |
| Sequencing submission form | out | — | no |
