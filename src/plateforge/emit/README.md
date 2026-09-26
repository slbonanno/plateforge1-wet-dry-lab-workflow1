# emit

Files for external systems: instrument worklists, order forms, reports.

**Inputs:** planned transfers from `assay.steps`, clone tables from `library`
**Outputs:** files a human loads, reviews and runs

## Confidence, not just format

Every emitter declares how well its format is actually known (decision 0021),
because the failure here is not a crash -- it is a plausible-looking file that
moves the wrong volume to the wrong well.

| | meaning |
|---|---|
| `verified` | a real file from the instrument is in `fixtures/` and the test reads it |
| `documented` | publicly specified, two independent sources agree; stamped when emitted |
| `sketch` | refuses, and says what would unblock it |

```python
from plateforge.emit import worklists

worklists.available()                     # what exists, and how well we know it
worklists.emit(plan.transfer_frame, "generic")                    # always fine
worklists.emit(plan.transfer_frame, "tecan_evo", allow_unverified=True)
```

Anything below `verified` needs `allow_unverified=True` and is written with a
`.CAVEAT.txt` beside it naming the sources.

Among the **instrument worklists**, nothing is verified but our own format --
a test asserts that, so the day a real instrument file lands, that test is
what changes.

The **sequencing order forms** are the exception, and the first emitters in
this repo that are genuinely verified: all three vendors' real templates are
in `fixtures/sequencing/`, and the tests read them.

## Sequencing order forms

```python
from plateforge.emit import sequencing

sequencing.available()                              # and how each is known
form = sequencing.order(clones, "genewiz", primer="CMV_F")
form.write("out/")                                  # plus a .NOTES.txt
```

| form | vendor | wells | fill order |
|---|---|---|---|
| `genewiz` *(default)* | Azenta / GENEWIZ | `A01` | both, as two columns |
| `elim` | ELIM Biopharmaceuticals | `A1` | row-major |
| `ucberkeley` | UC Berkeley DNA Sequencing Facility | `A1` | column-major |
| `generic` | none | `A01` | ours |

Three vendors, three well spellings, two fill directions (decision 0025). A
column-major plate pasted into a row-major form is a 96-well transposition
that every later step preserves and nothing downstream can detect -- which is
rule 3's whole reason for existing. Each emitter converts once, on the way out.

Azenta's is the one that fails quietly: it gives `Well (H)` and `Well (V)` as
two orderings of the same 96 positions, paired row by row, and the first
version of the emitter had them the wrong way round. A test now reads their
template and asserts our pairing equals theirs for all 95 rows.

Each form's own stated limits are enforced *before* writing, not discovered at
upload: Azenta's 500-sample ceiling, ELIM's 50-character names and controlled
vocabularies, and Berkeley's requirement to leave at least one well empty --
which is how that facility confirms plate orientation, so a full 96-sample
plate is refused rather than submitted.

## Adding an instrument

One function plus a decorator (rule 5):

```python
@EMITTERS.register("my_robot", software="...", confidence=DOCUMENTED,
                   suffix=".csv")
def my_robot(transfers, **kwargs) -> Worklist:
    ...
```

## The lab

`core.lab` holds which instruments this lab has, in
`$PLATEFORGE_DATA/lab.json`. A run may override it, and the override is
recorded rather than assumed.

```python
from plateforge.core import lab
lab.write({"name": "...", "liquid_handler": "integra_vialab",
           "has_single_channel_module": True, "bli": "octet"})
lab.instrument_for("liquid_handler")
```

## Nothing here drives an instrument

Every output is a file a human reviews first. That is deliberate: a review
step between a computed volume and a moving pipette is the only thing between
an arithmetic error and ninety-six wasted wells.
