"""Reverse translation with codon optimization and restriction site avoidance.

Two jobs. Turn a designed amino acid sequence into DNA a vendor can synthesise,
and keep that DNA free of the Type IIs sites the assembly depends on: a BsaI
site inside an insert is cut during Golden Gate and the assembly fails.

The codon table below is the commonly cited Homo sapiens usage. It is a
DEFAULT, not a verified reference -- expression is sensitive to it, vendors
have their own, and `usage=` takes a replacement. Everything else here is
mechanical and does not depend on the exact frequencies.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field

# aa -> [(codon, relative usage), ...], highest first.
HUMAN_USAGE: dict[str, list[tuple[str, float]]] = {
    "A": [("GCC", 0.40), ("GCT", 0.26), ("GCA", 0.23), ("GCG", 0.11)],
    "R": [("AGA", 0.21), ("AGG", 0.20), ("CGG", 0.20), ("CGC", 0.18), ("CGA", 0.11), ("CGT", 0.08)],
    "N": [("AAC", 0.53), ("AAT", 0.47)],
    "D": [("GAC", 0.54), ("GAT", 0.46)],
    "C": [("TGC", 0.55), ("TGT", 0.45)],
    "Q": [("CAG", 0.75), ("CAA", 0.25)],
    "E": [("GAG", 0.58), ("GAA", 0.42)],
    "G": [("GGC", 0.34), ("GGA", 0.25), ("GGG", 0.25), ("GGT", 0.16)],
    "H": [("CAC", 0.59), ("CAT", 0.41)],
    "I": [("ATC", 0.48), ("ATT", 0.36), ("ATA", 0.16)],
    "L": [("CTG", 0.41), ("CTC", 0.20), ("TTG", 0.13), ("CTT", 0.13), ("CTA", 0.07), ("TTA", 0.07)],
    "K": [("AAG", 0.57), ("AAA", 0.43)],
    "M": [("ATG", 1.00)],
    "F": [("TTC", 0.55), ("TTT", 0.45)],
    "P": [("CCC", 0.33), ("CCT", 0.28), ("CCA", 0.27), ("CCG", 0.11)],
    "S": [("AGC", 0.24), ("TCC", 0.22), ("TCT", 0.18), ("AGT", 0.15), ("TCA", 0.15), ("TCG", 0.06)],
    "T": [("ACC", 0.36), ("ACA", 0.28), ("ACT", 0.24), ("ACG", 0.12)],
    "W": [("TGG", 1.00)],
    "Y": [("TAC", 0.57), ("TAT", 0.43)],
    "V": [("GTG", 0.47), ("GTC", 0.24), ("GTT", 0.18), ("GTA", 0.11)],
    "*": [("TGA", 0.52), ("TAA", 0.28), ("TAG", 0.20)],
}

# Type IIs recognition sites that must not appear inside an insert.
ENZYME_SITES = {
    "BsaI": "GGTCTC",
    "BsmBI": "CGTCTC",
    "Esp3I": "CGTCTC",
    "BbsI": "GAAGAC",
    "SapI": "GCTCTTC",
}

_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def reverse_complement(dna: str) -> str:
    return dna.translate(_COMPLEMENT)[::-1]


def sites_to_avoid(enzymes: tuple[str, ...] = ("BsaI",)) -> list[str]:
    """Recognition sites and their reverse complements -- a Type IIs enzyme
    cuts either strand, so both orientations have to go."""
    out: list[str] = []
    for name in enzymes:
        site = ENZYME_SITES[name.strip()] if name.strip() in ENZYME_SITES else name.strip().upper()
        for variant in (site, reverse_complement(site)):
            if variant not in out:
                out.append(variant)
    return out


@dataclass
class Optimized:
    dna: str
    protein: str
    gc: float
    sites_removed: int = 0
    gc_swaps: int = 0
    unresolved: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    notes: str = ""

    @property
    def clean(self) -> bool:
        return not self.unresolved and not self.warnings


def homopolymers(dna: str, limit: int = 6) -> list[str]:
    """Runs a synthesis vendor will complain about."""
    return [m.group(0) for m in re.finditer(r"(A{%d,}|C{%d,}|G{%d,}|T{%d,})"
                                            % (limit, limit, limit, limit), dna)]


def _codon_gc(codon: str) -> int:
    return sum(1 for c in codon if c in "GC")


def gc_fraction(dna: str) -> float:
    return (dna.count("G") + dna.count("C")) / len(dna) if dna else 0.0


def translate(dna: str, usage: dict | None = None) -> str:
    """Back-translate DNA to protein, for verifying a round trip."""
    table = {}
    for aa, codons in (usage or HUMAN_USAGE).items():
        for codon, _w in codons:
            table[codon] = aa
    return "".join(table.get(dna[i:i + 3], "X") for i in range(0, len(dna) - 2, 3))


def optimize(protein: str, enzymes: tuple[str, ...] = ("BsaI",),
             usage: dict | None = None, max_passes: int = 12,
             gc_window: tuple[float, float] = (0.40, 0.60),
             min_usage: float = 0.10, homopolymer_limit: int = 6) -> Optimized:
    """Reverse translate, balance GC, then remove forbidden sites.

    Taking the highest-usage codon everywhere pushes GC to ~70% for an
    antibody V region, because the common human codons are GC-rich. Synthesis
    vendors flag that, so codons are swapped toward `gc_window` -- only for
    synonymous alternatives still above `min_usage`, so balancing GC does not
    quietly introduce rare codons.

    Only synonymous alternatives are ever considered, so the protein is
    unchanged by construction, which is checked rather than assumed.
    """
    usage = usage or HUMAN_USAGE
    protein = protein.strip().upper()
    unknown = sorted({aa for aa in protein if aa not in usage})
    codons = [usage[aa][0][0] if aa in usage else "NNN" for aa in protein]
    avoid = sites_to_avoid(enzymes)

    # --- GC balancing -----------------------------------------------------
    # Aim for the middle of the window, not merely inside it. Stopping at the
    # first acceptable value parks the result on the boundary, where later
    # site and homopolymer repairs push it back out and every clone comes back
    # flagged for being 0.1% over. Mid-window is also simply better DNA.
    lo, hi = gc_window
    target = (lo + hi) / 2
    tolerance = (hi - lo) / 4
    swaps = 0
    for _ in range(len(codons) * 2):
        gc = gc_fraction("".join(codons))
        if abs(gc - target) <= tolerance:
            break
        want_less = gc > target
        best = None
        for i, (aa, current) in enumerate(zip(protein, codons)):
            for alt, weight in usage.get(aa, []):
                if alt == current or weight < min_usage:
                    continue
                delta = _codon_gc(alt) - _codon_gc(current)
                if (want_less and delta < 0) or (not want_less and delta > 0):
                    score = (abs(delta), weight)
                    if best is None or score > best[0]:
                        best = (score, i, alt)
        if best is None:
            break
        _score, i, alt = best
        codons[i] = alt
        swaps += 1

    # --- repair: restriction sites and homopolymer runs --------------------
    # Both are the same problem -- a forbidden substring -- so both are fixed
    # the same way: re-code a codon the offender spans, preferring the
    # highest-usage synonymous alternative that removes it without creating
    # another. A vendor rejects a long run as readily as an assembly fails on
    # an internal BsaI site.
    def first_problem(dna: str):
        for site in avoid:
            at = dna.find(site)
            if at >= 0:
                return at, len(site)
        m = re.search(r"(A{%d,}|C{%d,}|G{%d,}|T{%d,})"
                      % ((homopolymer_limit,) * 4), dna)
        if m:
            return m.start(), len(m.group(0))
        return None

    def badness(dna: str, lo_i: int, hi_i: int) -> int:
        """Severity, not a count. A run of 12 A's split into 6 and 5 is still
        one violation but is strictly better, and scoring by count would
        reject the swap and give up -- which is what happened on a poly-K
        tract, where lysine's only alternative codon is AAG.
        """
        window = dna[max(0, lo_i - 12):hi_i + 12]
        score = 1000 * sum(window.count(site) for site in avoid)
        for m in re.finditer(r"(A{%d,}|C{%d,}|G{%d,}|T{%d,})"
                             % ((homopolymer_limit,) * 4), window):
            score += len(m.group(0))
        return score

    removed = 0
    for _ in range(max_passes * 4):
        dna = "".join(codons)
        found = first_problem(dna)
        if found is None:
            break
        at, span = found
        first, last = at // 3, (at + span - 1) // 3
        fixed = False
        for i in range(first, min(last + 1, len(codons))):
            aa = protein[i] if i < len(protein) else None
            best = None
            for alt, weight in usage.get(aa, []):
                if alt == codons[i]:
                    continue
                trial = codons[:i] + [alt] + codons[i + 1:]
                score = badness("".join(trial), at, at + span)
                if best is None or score < best[0]:
                    best = (score, alt, weight)
            if best and best[0] < badness(dna, at, at + span):
                codons[i] = best[1]
                removed += 1
                fixed = True
                break
        if not fixed:
            break

    dna = "".join(codons)
    remaining = [s for s in avoid if s in dna]
    back = translate(dna, usage)
    if back.replace("X", "") != protein.replace("".join(unknown), "") and not unknown:
        raise AssertionError("reverse translation changed the protein")

    notes = ""
    if unknown:
        notes = f"no codon for {', '.join(unknown)}; emitted as NNN"

    gc = gc_fraction(dna)
    warnings: list[str] = []
    if not (lo <= gc <= hi):
        warnings.append(f"GC {gc:.1%} outside {lo:.0%}-{hi:.0%}")
    runs = homopolymers(dna, homopolymer_limit)
    if runs:
        warnings.append(f"homopolymer runs: {', '.join(sorted(set(runs)))}")

    return Optimized(dna=dna, protein=protein, gc=gc,
                     sites_removed=removed, gc_swaps=swaps,
                     unresolved=remaining + ([f"unknown residues: {unknown}"] if unknown else []),
                     warnings=warnings, notes=notes)
