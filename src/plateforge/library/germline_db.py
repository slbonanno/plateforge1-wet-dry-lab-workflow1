"""Germline reference sequences, by gene name.

Needed for two jobs the construct step depends on: the fixed light chain an
scFv is built against, and the `Hu_germline` annotation on every ordered clone.

There is no offline germline database in any pip package, and IMGT/GENE-DB and
OGRDB are both plain HTTP downloads that may or may not be reachable from a
given machine -- OPIG taught us not to assume. So this is a registry of
sources, tried in order, with one backend that always works:

  pool    consensus of OAS's own germline_alignment_aa, per gene, mapped onto
          IMGT columns. No network, and it is IgBlast's IMGT-derived germline
          rather than a guess -- but it only covers genes present in the pool,
          and only the span the reads covered.
  fasta   an IMGT/GENE-DB FASTA already on disk. Authoritative and complete.
  url     the same, fetched. Works where the host is reachable.

Every result carries where it came from and how well supported it is, because
a germline assembled from 12 reads and one downloaded from IMGT should not be
used interchangeably without knowing which is which.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from ..core import registry
from . import germlines, imgt

SOURCES = registry.Registry("germline source")

IMGT_GENEDB_URL = (
    "https://www.imgt.org/download/GENE-DB/"
    "IMGTGENEDB-ReferenceSequences.fasta-AA-WithGaps-F+ORF+inframeP"
)


@dataclass
class Germline:
    """One germline reference, with its provenance."""
    gene: str
    sequence: str                       # ungapped amino acids
    source: str
    by_position: dict[str, str] = field(default_factory=dict)   # IMGT column -> residue
    support: int = 0                    # sequences the consensus drew on (pool only)
    allele: str | None = None
    functionality: str | None = None
    notes: str = ""

    @property
    def is_derived(self) -> bool:
        """True when this was inferred from observed data rather than looked up."""
        return self.source == "pool"

    def summary(self) -> str:
        where = f"{self.source}" + (f" (n={self.support})" if self.is_derived else "")
        return f"{self.gene} [{where}] {len(self.sequence)} aa"


# --- pool-derived ----------------------------------------------------------

@SOURCES.register("pool", needs_network=False, authoritative=False)
def from_pool(gene: str, rows: pd.DataFrame | None = None,
              min_support: int = 5) -> Germline | None:
    """Per-column consensus of the germline OAS reports for this gene.

    Positions no read covered are absent rather than guessed, so a germline
    built from 5' truncated reads is short and says so.
    """
    if rows is None:
        from . import pool as _pool
        idx = _pool.index(where="v_gene = ?", params=(gene,))
        if idx.empty:
            return None
        rows = _pool.fetch(list(idx["seq_id"]),
                           columns=["seq_id", "v_gene", "aa_gapped",
                                    "germline_aa", "anarci_numbering"])
    rows = rows[rows["v_gene"] == gene]
    if rows.empty:
        return None

    aln = imgt.build(rows)
    if aln is None:
        return None
    by_pos, used = imgt.germline_row(rows, list(aln.matrix.columns))
    if used < min_support or not by_pos:
        return None

    ordered = sorted(by_pos, key=imgt.position_key)
    return Germline(
        gene=gene,
        sequence="".join(by_pos[p] for p in ordered),
        source="pool",
        by_position=by_pos,
        support=used,
        notes=("consensus of observed germline calls; covers only the span the "
               "reads covered"),
    )


# --- IMGT FASTA ------------------------------------------------------------

_IMGT_HEADER = re.compile(r"^>([^|]*)\|([^|]+)\|([^|]*)\|([^|]*)\|")


def parse_imgt_fasta(text: str) -> dict[str, Germline]:
    """Parse IMGT/GENE-DB FASTA.

    Headers are pipe-delimited:
    >accession|IGHV3-23*01|Homo sapiens|F|V-REGION|...

    Keeps the first allele seen per gene, which for IMGT's ordering is *01,
    and drops IMGT's dot gaps from the stored sequence.
    """
    out: dict[str, Germline] = {}
    name = allele = func = None
    chunks: list[str] = []

    def flush():
        if not name or not chunks:
            return
        gene = name.split("*", 1)[0]
        if gene in out:
            return
        raw = "".join(chunks)
        out[gene] = Germline(
            gene=gene,
            sequence=raw.replace(".", "").replace("-", "").strip("*"),
            source="fasta",
            allele=name,
            functionality=func,
            notes="IMGT/GENE-DB reference",
        )

    for line in text.splitlines():
        if line.startswith(">"):
            flush()
            chunks = []
            m = _IMGT_HEADER.match(line)
            name, func = (m.group(2), m.group(4)) if m else (None, None)
            allele = name
        elif name:
            chunks.append(line.strip())
    flush()
    return out


@SOURCES.register("fasta", needs_network=False, authoritative=True)
def from_fasta(gene: str, path: str | Path) -> Germline | None:
    text = Path(path).read_text()
    return parse_imgt_fasta(text).get(gene)


@SOURCES.register("url", needs_network=True, authoritative=True)
def from_url(gene: str, url: str = IMGT_GENEDB_URL) -> Germline | None:
    import urllib.request
    req = urllib.request.Request(url, headers={"User-Agent": "plateforge"})
    with urllib.request.urlopen(req, timeout=60) as resp:
        text = resp.read().decode("utf-8", errors="replace")
    found = parse_imgt_fasta(text).get(gene)
    if found:
        found.source = "url"
    return found


# --- resolution ------------------------------------------------------------

def resolve(gene: str, order: tuple[str, ...] = ("fasta", "pool"),
            fasta_path: str | Path | None = None,
            rows: pd.DataFrame | None = None,
            url: str = IMGT_GENEDB_URL) -> Germline | None:
    """First source that yields a germline for this gene.

    Default order prefers an authoritative local FASTA and falls back to the
    pool, so a machine with IMGT data gets the real thing and one without still
    gets something usable and clearly labelled.
    """
    gene = germlines.gene(gene) or gene
    for name in order:
        try:
            if name == "fasta":
                if not fasta_path:
                    continue
                found = from_fasta(gene, fasta_path)
            elif name == "url":
                found = from_url(gene, url)
            elif name == "pool":
                found = from_pool(gene, rows=rows)
            else:
                found = SOURCES.get(name)(gene)
        except Exception:                                   # noqa: BLE001
            continue
        if found:
            return found
    return None


def coverage_report(genes: list[str], **kwargs) -> pd.DataFrame:
    """What each gene resolved to, and from where. Run before trusting a batch."""
    rows = []
    for gene in genes:
        found = resolve(gene, **kwargs)
        rows.append({
            "gene": gene,
            "resolved": found is not None,
            "source": found.source if found else None,
            "length": len(found.sequence) if found else 0,
            "support": found.support if found else 0,
            "derived": found.is_derived if found else None,
        })
    return pd.DataFrame(rows)
