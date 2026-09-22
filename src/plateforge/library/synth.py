"""Synthetic OAS-format data units.

Two jobs. It lets the parser be tested without network access or a 1.1 TB
download, and it is the seed of the dummy-data generator the later modules need.

The output is byte-for-byte in OAS layout: a CSV-quoted JSON metadata line,
then AIRR column names, then rows.
"""
from __future__ import annotations

import csv
import gzip
import json
import random
from pathlib import Path

# Framework context per germline, trimmed to what the pipeline actually reads.
FRAMEWORKS = {
    "IGHV3-23": ("EVQLLESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKG"
                 "RFTISRDNSKNTLYLQMNSLRAEDTAVYYCAK", "WGQGTLVTVSS", 12, 4),
    "IGHV1-69": ("QVQLVQSGAEVKKPGSSVKVSCKASGGTFSSYAISWVRQAPGQGLEWMGGIIPIFGTANYAQKFQG"
                 "RVTITADESTSTAYMELSSLRSEDTAVYYCAR", "WGQGTLVTVSS", 18, 5),
    "IGHV3-53": ("EVQLVESGGGLIQPGGSLRLSCAASGFTVSSNYMSWVRQAPGKGLEWVSVIYSGGSTYYADSVKG"
                 "RFTISRDNSKNTLYLQMNSLRAEDTAVYYCAR", "WGQGTLVTVSS", 9, 3),
}
_AA = "ACDEFGHIKLMNPQRSTVWY"
_HYDROPHOBIC = "AFILMVWY"


def _cdr3(rng: random.Random, mean_len: int, spread: int, hydrophobic: bool) -> str:
    n = max(5, int(rng.gauss(mean_len, spread)))
    alphabet = _HYDROPHOBIC + _AA if hydrophobic else _AA
    return "AR" + "".join(rng.choice(alphabet) for _ in range(n - 4)) + "DY"


def _mutate(rng: random.Random, cdr3: str, n_mut: int) -> str:
    """A clonal sibling: the same loop with a few point substitutions."""
    chars = list(cdr3)
    for _ in range(n_mut):
        i = rng.randrange(2, max(3, len(chars) - 2))
        chars[i] = rng.choice(_AA)
    return "".join(chars)


def _cdr3_population(rng: random.Random, n: int, mean_len: int, spread: int,
                     hydrophobic: bool, clonality: float) -> list[str]:
    """Draw n CDRH3s where a fraction are clonal siblings of a founder.

    Real repertoires and real panning outputs are dominated by expanded
    lineages: many near-identical loops plus a tail of singletons. Without
    this, every synthetic sequence is already maximally distinct and the
    diversity sampler has nothing to do.
    """
    out: list[str] = []
    while len(out) < n:
        founder = _cdr3(rng, mean_len, spread, hydrophobic)
        out.append(founder)
        if rng.random() < clonality:
            for _ in range(rng.randint(3, 25)):
                if len(out) >= n:
                    break
                out.append(_mutate(rng, founder, rng.randint(1, 3)))
    rng.shuffle(out)
    return out[:n]


def make_unit(path: str | Path, n: int = 500, genes: list[str] | None = None,
              species: str = "human", seed: int = 0, liability_rate: float = 0.05,
              clonality: float = 0.45, gzip_output: bool = True) -> Path:
    """Write a synthetic OAS data unit. Returns the path.

    `clonality` is the chance that a founder CDRH3 spawns an expanded lineage
    of near-identical siblings, as happens in a real campaign.
    """
    rng = random.Random(seed)
    genes = genes or list(FRAMEWORKS)
    path = Path(path)

    meta = {
        "Run": "SYNTHETIC", "Link": "n/a", "Author": "plateforge synthetic",
        "Species": species, "Strain": "none", "Age": "none", "BSource": "PBMC",
        "BType": "Unsorted-B-Cells", "Vaccine": "None", "Disease": "None",
        "Subject": "synthetic", "Longitudinal": "no", "Chain": "Heavy",
        "Unique sequences": n, "Isotype": "IGHG", "Total sequences": n,
        "Organism": species, "MiAIRR": "yes",
    }
    columns = ["sequence", "locus", "stop_codon", "vj_in_frame", "v_frameshift",
               "productive", "v_call", "d_call", "j_call", "sequence_alignment",
               "sequence_alignment_aa", "cdr3", "cdr3_aa", "junction_aa",
               "v_identity", "j_identity", "Redundancy", "ANARCI_status", "Isotype"]

    # Draw each germline's loops as a population so lineages stay within a gene.
    assignments = [rng.choice(genes) for _ in range(n)]
    loops: dict[str, list[str]] = {}
    for g in set(assignments):
        _, _, mean_len, spread = FRAMEWORKS[g]
        loops[g] = _cdr3_population(rng, assignments.count(g), mean_len, spread,
                                    hydrophobic=(g == "IGHV1-69"), clonality=clonality)

    opener = gzip.open if gzip_output else open
    with opener(path, "wt", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([json.dumps(meta)])
        w.writerow(columns)
        for g in assignments:
            fr, fr4, mean_len, spread = FRAMEWORKS[g]
            cdr3 = loops[g].pop()
            aa = fr + cdr3 + fr4
            liability = "|Unusual residue|" if rng.random() < liability_rate else "||"
            w.writerow([
                "N" * 30, "IGH", "F", "T", "F", "T",
                f"{g}*0{rng.randint(1, 3)}", "IGHD3-10*01",
                f"IGHJ{rng.randint(1, 6)}*0{rng.randint(1, 2)}",
                "N" * 30, aa, "NNN", cdr3, "C" + cdr3 + "W",
                round(rng.uniform(0.85, 1.0), 3), round(rng.uniform(0.85, 1.0), 3),
                rng.randint(1, 50), liability, "IGHG",
            ])
    return path
