# Vendor order tables — what is known and what is assumed

Rule 6: never write a file emitter without a real example of the template in
`fixtures/`. This file records the state of that evidence, so a future session
does not have to rediscover which layouts were guessed.

## Status

| Vendor | Product | Layout verified | Template in `fixtures/` | Limits source |
|---|---|---|---|---|
| (ours) | `generic` | n/a — every column is ours | n/a | none applied |
| IDT | eBlocks Gene Fragments | **no** | **no** | published product limits |
| GenScript | — | not written | no | — |
| Genewiz / Azenta | — | not written | no | — |
| Twist | — | not written | no | — |

## IDT eBlocks

**Assumed layout:** two columns, `Name` and `Sequence`, one row per fragment.

**Not verified.** No real IDT bulk-upload template has been seen. Until one is
in `fixtures/idt/`, `vendors.emit(..., "idt_eblocks")` refuses unless called
with `allow_unverified=True`, and any file it writes is accompanied by
`order_idt_eblocks.CAVEAT.txt`.

**Limits applied regardless**, because a fragment outside them fails whatever
the column headings say:

- length 300–1500 bp
- GC 25–65%
- no homopolymer run longer than 8

These are published eBlocks Gene Fragment limits. They have not been checked
against a current spec sheet or against a specific account's quote, and IDT's
own complexity screen is stricter than any of them — it evaluates repeats,
hairpins and local GC windows that are not modelled here. A fragment that
passes these checks can still be rejected at the quote stage.

## What to capture when a template arrives

1. Save the real file, unmodified, under `fixtures/<vendor>/`.
2. Note the exact column headers, including case and any trailing spaces.
3. Note whether the upload wants one row per well or one row per plate, and
   whether well positions are `A1` or `A01` — this project is always `A01`
   internally (rule 3) and converts at the boundary if a vendor differs.
4. Flip `template_verified=True` on that emitter, in the same commit.

## What was tried

- 2026-09-23: IDT, GenScript and Genewiz upload templates are behind account
  logins, so none could be retrieved from this environment. No unauthenticated
  URL was found that serves them.
