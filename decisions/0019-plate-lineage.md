# 0019 — Plates are containers; the clone is the thing being tracked

Date: 2026-09-23
Status: accepted — built 2026-09-23, minus the instrument emitters

## The problem

96 clones are about to make a journey:

    ordered DNA → transfection → culture → harvest supernatant
      → ProteinA Dynabead IP → wash → low-pH elution → neutralise → magnet
      → IgG prep → BLI (concentration) → normalise → ELISA

Each arrow is a new physical plate with a new barcode holding the *same 96
samples*, and each arrow is something a human or a robot did, at a time, with
reagents from lots, sometimes wrongly. The system has to answer, months later:
"what is in well D07 of plate PF-000312, and everything that happened to it."

It also has to let steps be added, removed or reordered without a rewrite,
because the protocol will change.

## The decision

**The plate is not the sample.** A `CLN` clone id is the through-line; plates
are containers it passes through. Every well of every plate records which clone
it descends from and which well it came from. That single choice is what makes
steps addable and removable: a step is a function from plates to plates, and
nothing downstream depends on how many steps preceded it.

Three new artifact types, all already reserved in `core.ids`:

| Type | Is | Key fields |
|---|---|---|
| `PLT` | one physical plate | barcode, format, `content_kind`, parent plate(s) |
| `XFR` | one transfer list | source plate/well → dest plate/well, volume |
| `RUN` | one step that happened | step type, inputs, outputs, params, reagent lots, operator, timestamps, deviations |

and one table that carries the weight:

```sql
CREATE TABLE plate_wells (
    plate_id        TEXT NOT NULL,
    well            TEXT NOT NULL,          -- always A01 (rule 3)
    clone_id        TEXT,                   -- the through-line; NULL for controls
    content_kind    TEXT NOT NULL,          -- dna | cells | supernatant | beads
                                            -- | eluate | igg_prep | dilution
    volume_ul       REAL,
    concentration   REAL,
    conc_units      TEXT,
    source_plate_id TEXT,                   -- where this well came from
    source_well     TEXT,
    role            TEXT,                   -- sample | control | blank | empty
    flags           TEXT,                   -- per-well deviations
    meta_json       TEXT NOT NULL DEFAULT '{}',
    PRIMARY KEY (plate_id, well)
);
```

`source_plate_id`/`source_well` is deliberately per-well rather than per-plate.
Most steps are 1:1 and preserve the map — which is what gives "same
organization as the original plate" for free — but a rearray, a
cherry-pick or a plate split is then the same mechanism, not a special case.

**Steps are a registry** (rule 5). A step type is a new file plus a decorator,
exposing two halves:

```python
@STEPS.register("proteinA_capture", consumes="supernatant", produces="eluate")
class ProteinACapture:
    def plan(self, source: Plate, **params) -> Plan:
        """Output plate layout + the worklists to run it. No side effects."""
    def record(self, plan: Plan, **observed) -> Run:
        """What actually happened: lots, timings, deviations, operator."""
```

`plan` and `record` are separate because they happen hours apart and because
`plan` must be re-runnable without inventing history. A step that a human
performs (the IP washes) has a `plan` that emits instructions and a `record`
that takes the human's word for it — the same shape as a robot step, so the
protocol does not fork based on who did the work.

The planned sequence, as step types:

| Step | Robot | Emits |
|---|---|---|
| `transfect` | no | plate layout |
| `harvest_supernatant` | Integra | transfer worklist, 500 µL |
| `add_beads` | Integra | reagent-addition worklist (pre-washed ProteinA Dynabeads) |
| `immunoprecipitate` | no | bench instructions; records incubation, washes |
| `elute_neutralise` | Integra | elution + neutralisation worklist |
| `magnet_transfer` | Integra | transfer worklist → IgG prep plate |
| `quantify_bli` | Octet | sample plate map; ingests the export |
| `normalise` | Integra | dilution worklist computed from BLI |
| `elisa` | — | the existing `assay` path |

## What is blocked, and on what

Rule 6 is not negotiable here — a worklist emitted against a guessed format
is a wasted plate.

- **Integra ASSIST/VIAFLO worklists.** No real file has been seen. Nothing is
  emitted until one is in `fixtures/integra/`. (Q7, Q15)
- **Octet BLI export.** Octet instruments do read 96-well sample plates, so
  the plate-in/plate-out shape holds; the export format still needs a real
  file before a reader is written. (Q7)
- **Bead and buffer volumes, incubation times, wash counts.** Protocol
  parameters, not software decisions. They belong in a JSON protocol
  definition (rule 8) supplied per campaign, not hardcoded.

## Why not one big LIMS table

Because the protocol will change, and a schema that encodes today's nine steps
as columns has to be migrated every time one is inserted. Steps as registry
entries plus a generic well-lineage table means adding `size_exclusion` between
elution and BLI is a new file and a line in a protocol JSON.

## Consequences

- `core.ids` needs no new prefixes; `PLT`, `XFR` and `RUN` already exist.
- `assay` gains the step registry and `plate_wells`; `emit` gains the
  instrument worklist emitters, one per format, each gated on a fixture.
- The clone table's `well` and `plate_barcode` become the *first* row of
  `plate_wells` rather than the only record of position.


## Built, and what changed on contact with the code

`assay.plates`, `assay.steps` and `assay.protocol` implement this. Two things
the sketch got wrong, both found by running the chain end to end:

**`barcode` cannot be UNIQUE.** Several steps change what is in a container
without moving it to a new one — adding beads, washing, eluting. Each state is
its own `PLT` because artifacts are never modified (rule 10), and they share a
barcode because it is the same piece of plastic. With the column unique,
`INSERT OR REPLACE` silently deleted the earlier state and the trace stopped
halfway along the chain. `plates.states_of(barcode)` now returns every state of
one container, oldest first.

**`core.stores.migrate` was quietly broken**, and this schema is what exposed
it. It split a table body on commas *before* stripping `--` comments, so a
comment containing a comma — which a comment enumerating a column's allowed
values almost always does — split one column definition in two. The first half
ended inside the comment and was discarded as empty, taking the **next column
with it**; the second half was read as a column definition, which is how a
table came to be asked for a column named `not`. Silently skipping a declared
column is the dangerous half: the schema says it exists, the database never
gains it, and it surfaces much later as a missing column at write time.

## What exists

| | |
|---|---|
| `plates.create/fill/save/load` | a plate, every well present, positions normalised |
| `plates.trace(plate, well)` | every well this one descends from, back to the order |
| `plates.where_is(clone)` | every container a clone has been in |
| `plates.states_of(barcode)` | every recorded state of one physical plate |
| `plates.layout_matches(a, b)` | assert the map survived, rather than trusting it |
| `steps.STEPS` | the registry; `plan` / `commit` / `record` |
| `protocol` | transfect → harvest → beads → IP → elute → magnet → BLI → normalise |

Verified on 96 wells: the chain runs end to end, `layout_matches(dna, final)`
holds, and D07 traces back eight hops to the plate that was ordered.

`normalise` is the only step whose volumes differ per well, computed from the
BLI concentrations by c₁v₁ = c₂v₂. Wells that cannot reach the target — too
dilute, or below the tip's minimum — are flagged on the output well and
counted in the step's params rather than quietly transferred wrong.

## Still not built, deliberately

No instrument file is written. A step declares the transfers it needs; turning
those into a worklist is `emit`'s job and is gated on a real example of each
machine's format (rule 6, Q7, Q15). Protocol volumes and incubation times are
parameters with placeholder defaults, declared in each step's `assumed` dict,
and belong in a campaign's protocol JSON (rule 8, Q19).
