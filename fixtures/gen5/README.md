# BioTek / Agilent Gen5 exports — nothing here yet

The grid finder in `assay.readers` does **not** need a fixture to work: it
locates a plate by its shape (consecutive column numbers, `A`..`H` beneath,
a rectangle of numbers between them), which is true of every microplate
export. That part is tested against `readers.write_gen5_like`.

What needs a real file is the **metadata**, which is currently collected
verbatim and left unparsed.

## What to export

One workbook from a session where several plates were read at once — ideally
a real target/control pair — straight out of Gen5, unmodified.

## What to note alongside it

- which sheet corresponds to which physical plate, and how you knew
- where the wavelength appears, and its exact wording (`Read 1:450`?)
- where the plate name or barcode appears, if anywhere
- whether the timestamp is per plate or per session
- whether a saturated well reads as a number, `OVRFLW`, or something else
- whether one sheet ever holds more than one read

With that, the preamble becomes parsed fields instead of strings, and
`assay.assign` can match on barcode without being told the mapping.

## Two more plates worth capturing, once hit calling is on real data

Both are cheap if they happen during a run that was going to happen anyway.

**A plate with a known positive control well.** Any well where you know an
antibody binds the coated antigen. `hits.DEFAULTS["min_signal"] = 0.2` is the
weakest number in the repo — it is the figure people reach for, not one we
measured (Q21). One real positive, and its background, replaces it.

**The same plate read twice: at the normal development time, and again after
it has over-developed.** Decision 0024 found that saturation goes with
*better* calls, not worse, because a well only reaches the ceiling when
something really bound. That was measured on the simulator, whose optics are
our own model, so it is currently a claim about our assumptions rather than
about TMB. Two reads of one plate test it directly — and tell us whether a
short second read is the cheap way to rank saturated wells (Q22).

Note alongside: the two development times, and whether the reader was
re-blanked between reads.
