# 0013 — Construct design and codon optimization

Date: 2026-09-22
Status: accepted

## Decisions taken

**Architecture: VH–(G4S)3–VL.** Chosen by the repository owner. Orientation and
linker are both parameters, because both change expression and neither has a
universally right answer. `LINKERS` holds G4S3, G4S4 and Whitlow; any literal
string also works.

**A clone is minted here.** `construct.build()` is where a `SEQ` becomes a
`CLN` (decisions/0005), because that is the point at which a sequence is
realised in a construct. Clone ids are minted, never derived, so designing the
same sequences twice gives two distinct sets rather than silently colliding.

**Assembly enzyme: BsaI.** Codon optimization removes `GGTCTC` and its reverse
complement `GAGACC` — a Type IIs enzyme cuts either strand, so both must go. An
internal site means the insert is cut during Golden Gate and the assembly
fails. Other enzymes are one entry in `ENZYME_SITES`.

## What the real output forced

**GC had to be balanced, not maximised.** Taking the highest-usage human codon
everywhere gives ~71% GC on an antibody V region, because the common human
codons are GC-rich. Synthesis vendors flag that. Codons are now swapped toward
the middle of the GC window, considering only synonymous alternatives above a
usage floor so balancing does not quietly introduce rare codons.

Aiming for *inside* the window was not enough: the balancer stopped at the
first acceptable value, parked on the boundary, and later repairs pushed it
back out — every clone in a 96-plate came back flagged for being 0.1% over.
It now targets mid-window, which is better DNA anyway.

**Homopolymers are as disqualifying as restriction sites.** Proline next to
other C-ending codons produced `CCCCCC`; vendors reject long runs. Sites and
runs are now the same problem — a forbidden substring — repaired the same way.

**Repairs are scored by severity, not by count.** A 12-base poly-A run split
into runs of 6 and 5 is still one violation, so counting violations rejected
the improving swap and gave up. This showed up on a poly-K tract, where lysine
has only `AAA` and `AAG`. Scoring by total residues in violating runs fixes it.

## The codon table

`HUMAN_USAGE` is the commonly cited Homo sapiens usage. It is a **default, not
a verified reference**: expression is sensitive to it and vendors have their
own. `usage=` takes a replacement. Everything else in the module is mechanical
and independent of the exact frequencies.

## Not done, and why

Vendor order forms (Q15) need a real example of each template in `fixtures/`
first. Gibson/HiFi/InFusion adapters (Q16) depend on the destination vector,
which is a wet-lab decision rather than a lookup. Both are deliberately absent
rather than guessed.
