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
