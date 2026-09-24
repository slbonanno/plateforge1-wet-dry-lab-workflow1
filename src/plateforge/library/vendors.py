"""Order tables, one registered emitter per vendor.

Rule 6 of CLAUDE.md: never write a file emitter without a real example of the
vendor's template in `fixtures/`. Format guesses are how this project fails --
a plate of DNA ordered against a guessed column layout is a real cost, paid
weeks later.

So this module separates two things that are usually conflated:

  * the **generic** table, which is ours. Every column it has, we defined.
    It is always safe to emit and always the source of truth.
  * a **vendor** table, which is theirs. Column names, order and casing are
    whatever their upload form wants, and we do not know that until we have
    seen one.

Each vendor emitter declares `template_verified`. While that is False the
emitter still runs -- you need something to look at -- but it refuses unless
asked explicitly, stamps its output, and records the caveat so it travels with
the order rather than living in someone's memory. Flip it to True in the same
commit that adds the fixture, never before.

Constraints (length, GC) are published product limits and are checked either
way, because a fragment outside them will be rejected or fail synthesis no
matter what the column headings say.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import pandas as pd

from ..core import registry
from . import codon

VENDORS = registry.Registry("vendor order format")


@dataclass
class Limits:
    """Published product limits. `notes` says what was checked and when."""
    min_bp: int | None = None
    max_bp: int | None = None
    gc_range: tuple[float, float] | None = None
    max_homopolymer: int | None = None
    notes: str = ""

    def check(self, dna: str) -> list[str]:
        out = []
        n = len(dna)
        if self.min_bp and n < self.min_bp:
            out.append(f"{n} bp is below the {self.min_bp} bp minimum")
        if self.max_bp and n > self.max_bp:
            out.append(f"{n} bp is above the {self.max_bp} bp maximum")
        if self.gc_range:
            gc = codon.gc_fraction(dna)
            lo, hi = self.gc_range
            if not lo <= gc <= hi:
                out.append(f"GC {gc:.0%} is outside {lo:.0%}-{hi:.0%}")
        if self.max_homopolymer:
            runs = codon.homopolymers(dna, limit=self.max_homopolymer)
            if runs:
                out.append(f"homopolymer run(s) over {self.max_homopolymer}: "
                           f"{', '.join(sorted(set(runs)))}")
        return out


@dataclass
class Order:
    """A vendor-shaped table plus everything that qualifies it."""
    vendor: str
    product: str
    table: pd.DataFrame
    issues: pd.DataFrame
    template_verified: bool
    limits: Limits
    caveats: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return self.issues.empty

    def as_dict(self) -> dict:
        return {
            "vendor": self.vendor,
            "product": self.product,
            "rows": int(len(self.table)),
            "columns": list(self.table.columns),
            "template_verified": self.template_verified,
            "rows_with_issues": int(self.issues["name"].nunique()
                                    if len(self.issues) else 0),
            "limits": {"min_bp": self.limits.min_bp, "max_bp": self.limits.max_bp,
                       "gc_range": list(self.limits.gc_range or ()),
                       "notes": self.limits.notes},
            "caveats": self.caveats,
        }


def _issues(names, seqs, limits: Limits) -> pd.DataFrame:
    rows = []
    for name, dna in zip(names, seqs):
        for problem in limits.check(str(dna)):
            rows.append({"name": name, "issue": problem})
    return pd.DataFrame(rows, columns=["name", "issue"])


# --- ours -------------------------------------------------------------------

GENERIC_LIMITS = Limits(notes="no vendor limits applied")


@VENDORS.register("generic", template_verified=True, product="none")
def generic(clones: pd.DataFrame, *, name_column: str = "clone_id",
            sequence_column: str = "order_dna_seq", **_) -> Order:
    """Our own table. Every column here is one we defined, so it cannot drift."""
    table = pd.DataFrame({
        "name": clones[name_column].astype(str),
        "sequence": clones[sequence_column].astype(str),
        "length_bp": clones[sequence_column].astype(str).str.len(),
        "gc": [round(codon.gc_fraction(str(s)), 3) for s in clones[sequence_column]],
        "well": clones.get("well"),
        "plate_barcode": clones.get("plate_barcode"),
        "construct": clones.get("construct"),
        "v_gene": clones.get("v_gene"),
        "cdr3_aa": clones.get("cdr3_aa"),
        "cloning": clones.get("cloning_strategy"),
        "vector": clones.get("vector"),
        "notes": clones.get("synthesis_warnings"),
    })
    return Order(vendor="generic", product="none", table=table,
                 issues=pd.DataFrame(columns=["name", "issue"]),
                 template_verified=True, limits=GENERIC_LIMITS)


# --- theirs -----------------------------------------------------------------

IDT_EBLOCKS = Limits(
    min_bp=300, max_bp=1500, gc_range=(0.25, 0.65), max_homopolymer=8,
    notes="Published eBlocks Gene Fragment limits. Not checked against an "
          "account's own quote or a current spec sheet.")


@VENDORS.register("idt_eblocks", template_verified=False, product="eBlocks Gene Fragments")
def idt_eblocks(clones: pd.DataFrame, *, name_column: str = "clone_id",
                sequence_column: str = "order_dna_seq", **_) -> Order:
    """IDT eBlocks bulk upload.

    Believed to be two columns, Name and Sequence. UNVERIFIED -- no real IDT
    template is in fixtures/ yet (Q15). The length and GC limits are published
    product limits and are checked regardless.
    """
    names = clones[name_column].astype(str)
    seqs = clones[sequence_column].astype(str)
    table = pd.DataFrame({"Name": names, "Sequence": seqs})
    return Order(
        vendor="IDT", product="eBlocks Gene Fragments", table=table,
        issues=_issues(names, seqs, IDT_EBLOCKS),
        template_verified=False, limits=IDT_EBLOCKS,
        caveats=["Column layout is assumed, not taken from a real IDT upload "
                 "template. Confirm against the site before ordering, and add "
                 "the template to fixtures/ (Q15)."])


def emit(clones: pd.DataFrame, vendor: str = "generic",
         allow_unverified: bool = False, **kwargs) -> Order:
    """Build one vendor's order table.

    An unverified layout will not be produced by accident: ask for it.
    """
    meta = VENDORS.meta(vendor)
    if not meta.get("template_verified", False) and not allow_unverified:
        raise ValueError(
            f"{vendor!r} has no real upload template in fixtures/ yet, so its "
            "column layout is a guess. Pass allow_unverified=True to emit it "
            "anyway (it will be stamped), or use vendor='generic'.")
    return VENDORS.get(vendor)(clones, **kwargs)


def available() -> pd.DataFrame:
    """Every registered vendor and whether its layout has been verified."""
    return pd.DataFrame([
        {"vendor": name, "product": VENDORS.meta(name).get("product"),
         "template_verified": VENDORS.meta(name).get("template_verified")}
        for name in sorted(VENDORS)
    ])
