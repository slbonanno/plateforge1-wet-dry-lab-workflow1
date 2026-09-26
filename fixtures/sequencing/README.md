# Sanger order forms — three vendors' real templates

Ground truth for `emit.sequencing` (rule 6). The tests read these files and
check that what we emit fits them; the GENEWIZ test compares our well pairing
against the template's own, row by row.

| file | vendor | wells | fill order | shape |
|---|---|---|---|---|
| `azenta_sanger_form_v2.xlsx` | Azenta / GENEWIZ **(default)** | `A01` | both, as two columns | one table, 500 rows |
| `elim_seq_orderform_96well.xls` | ELIM Biopharmaceuticals | `A1` | row-major | one table, 96 rows, dropdowns |
| `ucberkeley_full_plate_order_form.xlsx` | UC Berkeley DNA Sequencing Facility | `A1` | column-major | printed form, two side-by-side blocks |

Three vendors, three well spellings, two fill directions. A column-major plate
pasted into a row-major form is a 96-well transposition that every later step
preserves and nothing downstream can detect, which is why `core.wells`
normalises on the way in and each emitter converts once, on the way out.

## Per-form specifics the emitters enforce

**Azenta** — headers must not be renamed (the template says so in its own
notes); columns may be reordered and extra columns are ignored. Several
primers per row are separated by `;`, but pre-mixed submissions allow only
one. 500 rows. `Well (H)` and `Well (V)` are two orderings of the same 96
positions paired row by row — getting them the wrong way round transposes the
plate, and did, until a test compared our pairing to theirs.

**ELIM** — sample names are at most 50 characters and letters, numbers, dash
and underscore only. `Template Type`, `Premix?` and `GC Rich?` are dropdowns
with fixed option lists parked in columns P-R. One sheet per template plate.

**UC Berkeley** — "LEAVE AT LEAST ONE WELL EMPTY ON PLATE" is how the facility
confirms plate orientation, so it is a constraint and a full 96-sample plate
is refused. The form also needs a header block we do not fill: submitter, PI,
department, chartstring or PO, and the plate barcode.

## Still worth capturing

An **actual submitted form** from whichever vendor gets used, plus whatever
they send back on completion. The emitters are verified against the blank
templates; a round trip would also pin what the vendor does to the file.
