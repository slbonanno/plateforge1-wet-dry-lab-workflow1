# 0021 — Instruments are registry entries, and each declares how well we know it

Date: 2026-09-23
Status: accepted

## Context

The goal is an agent a lab member can ask questions of. That agent will one
day ask "which liquid handler do you have?" and write the right file. So the
pipeline has to hold more than one lab's equipment.

The obvious worry is that this means a restructure. It does not. Rule 5 has
said since the beginning that extension happens through registries, and a
liquid handler is exactly that shape: one function, one decorator. What
*needed* deciding was not the mechanism but the honesty policy, because the
failure mode here is not a crash — it is a plausible-looking file that moves
the wrong volume to the wrong well.

## Decision

### Three levels of confidence, declared per emitter

| | meaning |
|---|---|
| `verified` | a real file from this instrument is in `fixtures/` and the emitter's test reads it |
| `documented` | the format is specified publicly and **at least two independent sources agree** on the field order; emitted on request, stamped, never silently |
| `sketch` | the shape is known, the details are not; refuses, and says what would unblock it |

`emit()` requires `allow_unverified=True` for anything below `verified`, and
every such file is written with a `.CAVEAT.txt` beside it naming the sources.
**Nothing is `verified` today** except our own format, and a test asserts that
— so the day a fixture lands, that test is what changes.

### Where each instrument stands

| Instrument | Confidence | Why |
|---|---|---|
| generic (ours) | verified | every column is one we defined |
| Tecan EVOware / Fluent | documented | `.gwl` is a simple semicolon record; robotools and Dioscuri agree on the ten-field layout |
| Opentrons | documented | the Python Protocol API is public and versioned; what is *not* known is the deck |
| INTEGRA VIALAB | sketch | CSV import is a D-ONE feature and the templates are behind a login (Q15) |
| Hamilton VENUS | sketch | methods are a proprietary project format, not a worklist; PyHamilton is an integration, not an emitter |
| Beckman Biomek | sketch | transfer-file format not publicly specified |
| Agilent Bravo | sketch | VWorks protocols are XML with an unpublished schema |

Opentrons is the interesting case: the "worklist" is a program, so the format
is not in doubt at all — but slot numbers and labware names are guesses, and
the caveat says exactly that rather than implying the whole file is suspect.

### The instrument is a property of the lab

`core.lab` reads `$PLATEFORGE_DATA/lab.json`. A lab has the instruments it
has; that does not change per run, so it is read once rather than threaded
through every call. A run may override it — someone borrows a robot for an
afternoon — and `Lab.source` records that an override happened, because a
silent override is how a worklist gets written for the wrong machine.

An unconfigured lab still runs, on defaults that say they are defaults. A
fresh clone should not need a config file to do anything.

Unknown keys in the config are kept, not dropped (rule 7): a lab will own
something this schema did not anticipate.

## Why not drive the instruments directly

Every emitter writes a file a human loads, reviews and runs. That is not a
gap to close later. A review step between a computed volume and a moving
pipette is worth keeping, and it is the only thing standing between an
arithmetic error and ninety-six wasted wells.
