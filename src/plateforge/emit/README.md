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
`.CAVEAT.txt` beside it naming the sources. Nothing is verified yet but our
own format -- a test asserts that, so the day a fixture lands, that test is
what changes.

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
