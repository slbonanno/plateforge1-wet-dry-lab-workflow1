"""How a coding sequence becomes an orderable fragment.

`library.vector` says what the insert lands between. This says how it is
joined. The two are separate so that changing cloning strategy does not mean
redefining the vector, and changing vector does not mean rewriting the
strategy.

A strategy is a function of (coding DNA, vector) -> Fragment, registered by
name. Adding one is a new function plus a decorator; nothing existing moves.

    golden_gate   Type IIs (BsaI by default). Fusion sites come from the
                  vector's constant flanks, so all 96 inserts share one
                  vector. Scarless.
    gibson        Homology arms copied from the vector context. No enzyme
                  sites, so nothing has to be removed from the insert, but
                  the arms are long and the fragment gets expensive.
    blunt         The coding sequence with flanking pad only -- for ordering
                  before the cloning route is decided. Explicitly not
                  cloneable as-is, and says so.

Every Fragment carries its own layout, so a figure or a check can point at
exactly which stretch of the ordered DNA is what.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from ..core import registry
from . import codon, vector as vectors

STRATEGIES = registry.Registry("cloning strategy")

# A neutral filler used to pad short fragments up to a vendor's minimum, and
# to give a Type IIs enzyme something to grip outside its recognition site.
# Fixed, not random: two runs must design the same fragment. Checked at import
# for BsaI/BsmBI sites and homopolymers.
FILLER = ("TGACTGCATCAGTACGATCTGACTAGCATGCTAGTCAGATCGTACTGATGCAGTCAGTAC"
          "GATCTAGCATGCTAGTCAGATCGTACAGTCATGCTAGATCGTACGATCAGTCATGCTAGC")

# Enzymes need flanking DNA beyond their recognition site to cut efficiently.
MIN_FLANK = 8


@dataclass
class Fragment:
    """One orderable piece of DNA, and what every part of it is."""
    dna: str
    strategy: str
    vector: str
    layout: list[tuple[str, int, int]] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    cloneable: bool = True

    @property
    def length(self) -> int:
        return len(self.dna)

    @property
    def gc(self) -> float:
        return codon.gc_fraction(self.dna)

    def part(self, name: str) -> str:
        for label, start, stop in self.layout:
            if label == name:
                return self.dna[start:stop]
        raise KeyError(f"no part {name!r}; have {[l for l, _, _ in self.layout]}")

    def describe(self) -> str:
        rows = [f"{self.strategy} into {self.vector} — {self.length} bp, "
                f"GC {self.gc:.0%}"]
        rows += [f"  {start:>4}-{stop:<4} {label}"
                 for label, start, stop in self.layout]
        rows += [f"  ! {w}" for w in self.warnings]
        return "\n".join(rows)


def _filler(n: int, offset: int = 0) -> str:
    if n <= 0:
        return ""
    doubled = FILLER * (2 + (n + offset) // len(FILLER))
    return doubled[offset:offset + n]


def _pad_to(core: str, minimum: int, flank: int) -> tuple[str, str]:
    """(5', 3') filler: at least `flank` each side, more to reach `minimum`."""
    five = three = flank
    short = minimum - (len(core) + five + three)
    if short > 0:
        five += short // 2 + short % 2
        three += short // 2
    return _filler(five), _filler(three, offset=31)


def _internal_sites(dna: str, enzymes: tuple[str, ...]) -> list[str]:
    return [s for s in codon.sites_to_avoid(enzymes) if s in dna]


def _assemble(parts: list[tuple[str, str]]) -> tuple[str, list]:
    """Concatenate (label, dna) pairs and record where each landed."""
    dna, layout, at = [], [], 0
    for label, piece in parts:
        if not piece:
            continue
        dna.append(piece)
        layout.append((label, at, at + len(piece)))
        at += len(piece)
    return "".join(dna), layout


@STRATEGIES.register("golden_gate", enzyme="BsaI", scarless=True)
def golden_gate(cds: str, vec: vectors.Vector, *, enzyme: str = "BsaI",
                min_length: int = 0, flank: int = MIN_FLANK,
                spacer: str = "A") -> Fragment:
    """Type IIs fragment: pad - BsaI - N - OH5 - CDS - OH3 - N - BsaI - pad.

    The two overhangs are copies of vector sequence, not additions to it: the
    assembled plasmid reads upstream + CDS + downstream with nothing inserted
    at either seam. That is what scarless means here.
    """
    site = codon.sites_to_avoid((enzyme,))[0]
    rc = codon.reverse_complement(site)
    five, three = vec.fusion_sites()

    core_parts = [
        (f"{enzyme} site", site), ("spacer", spacer),
        ("5' fusion site", five), ("CDS", cds),
        ("3' fusion site", three), ("spacer", spacer),
        (f"{enzyme} site (rev)", rc),
    ]
    core = "".join(p for _, p in core_parts)
    pad5, pad3 = _pad_to(core, min_length, flank)
    dna, layout = _assemble([("5' pad", pad5)] + core_parts + [("3' pad", pad3)])

    warnings = list(vec.check_fusion_sites())
    # The recognition sites we put there are expected; anything beyond two of
    # each orientation is an internal site that would cut the fragment apart.
    if dna.count(site) > 1 or dna.count(rc) > 1:
        warnings.append(f"internal {enzyme} site: the fragment would be cut "
                        "somewhere other than its ends")
    joined = five + cds + three
    for extra in _internal_sites(joined, ("BsmBI", "BbsI")):
        warnings.append(f"contains {extra}, a Type IIs site for another enzyme")
    if len(cds) % 3:
        warnings.append(f"CDS is {len(cds)} bp, not a multiple of 3")
    runs = codon.homopolymers(dna, limit=6)
    if runs:
        warnings.append(f"homopolymer run(s): {', '.join(sorted(set(runs)))}")

    return Fragment(dna=dna, strategy="golden_gate", vector=vec.name,
                    layout=layout, warnings=warnings)


@STRATEGIES.register("gibson", enzyme=None, scarless=True)
def gibson(cds: str, vec: vectors.Vector, *, arm: int = 25,
           min_length: int = 0, flank: int = 0, **_) -> Fragment:
    """Homology arms copied from the vector on either side of the insert.

    No enzyme sites, so nothing has to be removed from the coding sequence --
    but the arms are real vector sequence, so this needs the vector's actual
    upstream and downstream DNA, not just its fusion sites.
    """
    up, down = vec.upstream_dna(), vec.downstream_dna()
    warnings = []
    if len(up) < arm or len(down) < arm:
        warnings.append(f"only {min(len(up), len(down))} bp of flanking "
                        f"sequence is declared; asked for {arm} bp arms")
    parts = [("5' homology arm", up[-arm:]), ("CDS", cds),
             ("3' homology arm", down[:arm])]
    core = "".join(p for _, p in parts)
    pad5, pad3 = _pad_to(core, min_length, flank)
    dna, layout = _assemble([("5' pad", pad5)] + parts + [("3' pad", pad3)])
    if len(cds) % 3:
        warnings.append(f"CDS is {len(cds)} bp, not a multiple of 3")
    return Fragment(dna=dna, strategy="gibson", vector=vec.name,
                    layout=layout, warnings=warnings)


@STRATEGIES.register("blunt", enzyme=None, scarless=False)
def blunt(cds: str, vec: vectors.Vector, *, min_length: int = 0,
          flank: int = 0, **_) -> Fragment:
    """The coding sequence and nothing else. Not cloneable as ordered.

    For pricing a synthesis run, or ordering before the route is settled.
    Marked `cloneable=False` so no downstream step can treat it as finished.
    """
    pad5, pad3 = _pad_to(cds, min_length, flank)
    dna, layout = _assemble([("5' pad", pad5), ("CDS", cds), ("3' pad", pad3)])
    return Fragment(
        dna=dna, strategy="blunt", vector=vec.name, layout=layout,
        cloneable=False,
        warnings=["no cloning features: this fragment cannot be assembled "
                  "into the vector as ordered"])


def design(cds: str, vec: vectors.Vector | str = "pcdna-igg1-dual-vk-v1",
           strategy: str = "golden_gate", **kwargs) -> Fragment:
    """Build the orderable fragment for one coding sequence."""
    if isinstance(vec, str):
        vec = vectors.get(vec)
    return STRATEGIES.get(strategy)(cds, vec, **kwargs)


def verify(fragment: Fragment, cds: str, vec: vectors.Vector,
           expected_protein: str | None = None) -> list[str]:
    """Check the assembled plasmid, not just the fragment.

    The failure this catches is the expensive one: a fragment that looks fine
    on its own but shifts the frame or changes a residue once it is in the
    vector. Checked by translating the joined ORF, which is the only way to
    know.
    """
    problems = []
    if fragment.part("CDS") != cds:
        problems.append("the fragment's CDS is not the sequence it was built from")

    # Every cassette, not only the insert's. A second cassette cannot be
    # broken by the insert, but it can be broken by an edit to the vector,
    # and this is the check that would catch it.
    for name, orf in vec.orfs(cds).items():
        if not orf:
            continue
        if len(orf) % 3:
            problems.append(f"{name}: ORF is {len(orf)} bp, not a multiple of 3")
        protein = codon.translate(orf)
        if protein.count("*") > 1 or (protein.count("*") == 1
                                      and not protein.endswith("*")):
            problems.append(f"{name}: internal stop codon")

    if expected_protein:
        got = codon.translate(vec.orf(cds)).rstrip("*")
        if expected_protein not in got:
            problems.append("the designed protein is not in the assembled ORF")
    return problems
