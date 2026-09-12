# 0005 — A clone is a sequence realised in a construct

Date: 2026-09-12
Status: accepted

## Context

A sequence can exist abstractly (a row in a catalog) or physically (cloned into
a phagemid, picked from a pan, rescaffolded onto human IgG in a new plasmid).
The same amino acid sequence may exist in several constructs at once, in
different plates, with different behaviour.

## Decision

`SEQ` is a catalog entry: an amino acid sequence and its annotations. It has no
physical existence.

`CLN` is created when a sequence is assigned to a construct. That is the point
at which the thing exists in the world — the sequence cloned out of an in-vivo
sample and living in a phagemid is the first clone.

When a clone is reformatted into a different construct — phagemid to IgG
expression plasmid — a **new** `CLN` is minted, linked to the previous one by a
`reformatted_from` edge, and both point at the same `SEQ`.

## Why a new CLN rather than a mutating one

The new construct is a physically different object: different plasmid,
different glycerol stock, different well, potentially different behaviour in
assay. Both often exist simultaneously and must be distinguishable. Minting a
new ID also keeps clones immutable, matching the rest of the ledger.

"Same molecule, different format" queries remain cheap: join on `seq_id`, or
walk the `reformatted_from` chain.

## Consequences

- Clone counts are per-construct, not per-molecule. Reports must say which.
- A rearray within the same construct does **not** mint a new clone; only a
  construct change does.
- This is a data convention, not a code structure, so it is inexpensive to
  revisit if the wet-lab reality turns out otherwise.
