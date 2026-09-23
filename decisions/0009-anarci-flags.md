# 0009 — Most ANARCI flags describe the read, not the molecule

Date: 2026-09-22
Status: accepted
Closes: Q10

## Context

`Filter(exclude_liabilities=True)` treated any non-empty `ANARCI_status` as a
liability. Run against Briney et al. 2019 it kept **0 of 125,743** sequences.

## Evidence

Flag counts over 19,864 sequences surviving gene and length filters:

| Flag | Count |
|---|---|
| `Shorter than IMGT defined: fw1` | 17,592 |
| `Deletions: 1..15, 73` | 11,471 |
| `Deletions: 1..14, 73` | 3,453 |
| `Shorter than IMGT defined: fw1, fw4` | 2,272 |
| `Missing Conserved Cysteine: 104` | 480 |

## Decision

`COVERAGE_PATTERNS` classifies `"shorter than imgt defined"` and `"deletions:"`
as descriptions of read coverage rather than the molecule. A sequence is a
liability only if it carries a flag outside that set.

Unrecognised flags are liabilities. A flag nobody has examined should fail
closed, and the funnel now makes an over-strict filter visible in one run.

## Why "Deletions" counts as coverage

Debatable, and recorded as such. A mid-domain deletion would be a real event.
But position 73 is deleted in roughly 90% of Briney sequences, and positions
1–15 in most — that is the amplicon starting inside FR1, not eighteen thousand
deletions. Rejecting them would reject a repertoire on the basis of primer
design. Revisit if a dataset ever shows sparse, isolated deletions.

`Missing Conserved Cysteine` stays a liability: that is a structural problem.

## Consequences

- The default filter is now usable on real repertoire data.
- `oas.flag_counts()` exists so this can be re-derived per dataset rather than
  inherited as folklore.
- A separate filter, `max_ambiguous`, was added in the same pass: real CDR3s
  come back as `XXXXXXXXXXYW`, and without it those reach the sampler and land
  in a picked library.
