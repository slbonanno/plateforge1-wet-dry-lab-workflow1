# 0008 — What "diverse enough" means

Date: 2026-09-22
Status: accepted
Closes: Q1

## Context

The sampler needs a definition of diversity that can be checked, not eyeballed.
"Enough similarity that it recapitulates a real campaign, enough difference that
it is not a few point mutants" has to become numbers.

## Decision

Three criteria, all evaluated after sampling by `diversity.Spec.check()`:

| Criterion | Default | Why |
|---|---|---|
| `max_pairwise_identity` | 0.80 | No two members are near-duplicates. Normalized edit distance on CDRH3, length-aware, so a short loop inside a long one is not counted as identical. |
| `min_cdr3_len_span` | 8 aa | The set covers a range of loop geometries rather than one length regime. |
| `min_genes` / `min_per_gene` | panel size / 1 | Every panel germline is represented. |

Selection is greedy farthest-point on CDRH3 within each germline's quota, run
per CDRH3 length bin. That maximises the *minimum* pairwise distance rather than
merely removing exact duplicates, which is the behaviour that matters when the
pool contains expanded clonal lineages.

The check returns a report rather than raising, so a set that misses a criterion
is still usable with the miss recorded.

## Duplicate loops across germlines

Added 2026-09-22 after a 96-pick on real Briney data returned
`max_pairwise_identity: 1.0`. `TKGVGSVRNDALHI` was selected twice, once as
IGHV3-23 and once as IGHV3-53.

Quotas are filled one germline at a time, and a real repertoire contains
identical CDRH3s assigned to different V genes — convergence, or an ambiguous
V call. Each gene's independent selection picked the same loop. At plate scale
that means two wells synthesising one molecule.

`sample(unique_cdr3=True)` is now the default: the candidate pool is
deduplicated on CDRH3 and each gene excludes loops already taken. When distinct
loops run out the sampler returns fewer than asked rather than repeating, and
the caller is told. Ordering a duplicate is worse than ordering 94 clones.

## Consequences

- The criteria are defaults on a dataclass, so a campaign can tighten or relax
  them without touching the sampler.
- Sampling is seeded and deterministic; the same seed gives the same set.
- Cost is O(k·N) distance computations, fine at the scale this project uses
  (≤ 1k input sequences, a few hundred picks).
