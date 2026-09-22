# 0007 — Scaffold panel and light chain

Date: 2026-09-22
Status: accepted
Closes: Q1 (partially), Q6

## Context

The library module needs a defined set of human heavy-chain germlines to sample
from, and one light chain to pair everything against during reformatting. The
goal is a synthetic set that plausibly recapitulates what an in-vivo discovery
campaign yields — not a portfolio of approved drugs.

## What the literature actually says

There is no clean ranked list of "most common therapeutic VH genes". The most
recent systematic pass mapped every Thera-SAbDab antibody in trials or approved
to its nearest human germline and found broad dispersion with no clear
preference for particular germline combinations; IGHV3-23 appears often, likely
because of its natural abundance rather than selection for it. IGHV1-46 is
numerically prominent (83 of 817) but many of those antibodies are not fully
human, so the count reflects an animal-derived germline mapped to its nearest
human match.

Three frequency distributions get conflated and should not be:

- therapeutic frequency (what got approved)
- naive repertoire frequency (what a healthy person carries)
- antigen-response frequency (what a campaign against a target yields)

We are simulating discovery output, so the latter two are the better model.

## Decision

Panel `default-v1`, human heavy chains only:

| Gene | Weight | Role | Why |
|---|---|---|---|
| IGHV3-23 | 0.60 | workhorse | Most common human VH; considered well-behaved; behind many approved therapeutics (bevacizumab, daratumumab, dupilumab, ranibizumab and others). The baseline. |
| IGHV1-69 | 0.25 | long hydrophobic CDRH3 | Dominant in antiviral responses (influenza, HIV, dengue, SARS-CoV-2); biased toward long hydrophobic loops. Supplies clones that behave unlike the baseline. |
| IGHV3-53 | 0.15 | short CDRH3 | Public response germline; ~10% of SARS-CoV-2 RBD antibodies against 0.5–2.6% in naive repertoires. Short loops, well-behaved; the third canonical shape. |

A second panel `single-v1` holds IGHV3-23 alone, for when downstream behavioral
variety is not needed yet.

Light chain: **IGKV1-39**, fixed, for all reformatting.

## Why these three rather than three popular ones

Germline is the only handle on *behavioral* difference available before real
assay data exists. A panel of three well-behaved germlines produces clones that
are interchangeable downstream, so every analysis passes for the wrong reason.
These three differ in CDRH3 length regime and hydrophobicity, which is what the
assay module needs to have something to discriminate.

## On IGKV1-39, honestly

It is the best-supported generic pairing partner, not a proven optimum. It sits
at roughly 13% of therapeutic repertoires and around 30% of validated common
light chains, and that enrichment has been linked to broad VH–VL compatibility;
separately, Fab light chains paired against a fixed VH showed five- to ten-fold
lower expression than the same VH with IGKV1-39, and unmutated IGKV1-39 was
argued to aid production stability across VH genes.

The authors of the common-light-chain work explicitly decline to infer intrinsic
superiority from those numbers, noting that structural databases over-represent
crystallizable antibodies and treating the enrichment as a working hypothesis.
So: chosen as the best-supported default, with known database bias in the
supporting evidence. Not chosen because it is optimal.

IGLV1-51 was considered. It is lambda, so it works as a specific pairing but not
as a single generic scaffold. Worth adding later as a second light chain
specifically to prove the reformatting code is not kappa-hardcoded.

## Consequences

- Panels are data (`Panel.as_dict()` / `from_dict()`), so changing the mix is a
  config edit, not a code change.
- Quotas are computed by weight with largest-remainder and a floor of one per
  gene, so small n still covers the panel.
- Paired-chain support remains open (Q6 partially). Unpaired heavy only for now.

## Sources

- Data-driven analyses of human antibody variable domain germlines, mAbs 2025 —
  https://www.tandfonline.com/doi/full/10.1080/19420862.2025.2507950
- Structural basis of a public antibody response to SARS-CoV-2 —
  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7302194/
- An allelic atlas of immunoglobulin heavy chain variable regions —
  https://pmc.ncbi.nlm.nih.gov/articles/PMC11746035/
- AI-guided design of common light chains, bioRxiv 2025 —
  https://www.biorxiv.org/content/10.1101/2025.10.11.681265v2.full
- Selecting for developability in eukaryotic cell display systems, US11499150 —
  https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/11499150
