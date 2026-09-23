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

# Framework context per germline: (FR1-FR3, FR4, mean CDRH3 length, spread).
#
# These framework strings are ILLUSTRATIVE, not authoritative IMGT references.
# They exist so synthetic sequences have realistic length and composition and
# so germlines differ from one another in CDRH3 regime. Nothing in the pipeline
# treats them as a germline standard: the alignment figure uses the observed
# consensus of the data itself as its reference row, which is also what works
# for real OAS data.
FRAMEWORKS = {
    "IGHV3-23": ("EVQLLESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKG"
                 "RFTISRDNSKNTLYLQMNSLRAEDTAVYYCAK", "WGQGTLVTVSS", 12, 4),
    "IGHV1-69": ("QVQLVQSGAEVKKPGSSVKVSCKASGGTFSSYAISWVRQAPGQGLEWMGGIIPIFGTANYAQKFQG"
                 "RVTITADESTSTAYMELSSLRSEDTAVYYCAR", "WGQGTLVTVSS", 18, 5),
    "IGHV3-53": ("EVQLVESGGGLIQPGGSLRLSCAASGFTVSSNYMSWVRQAPGKGLEWVSVIYSGGSTYYADSVKG"
                 "RFTISRDNSKNTLYLQMNSLRAEDTAVYYCAR", "WGQGTLVTVSS", 9, 3),
    "IGHV1-46": ("QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYYMHWVRQAPGQGLEWMGIINPSGGSTSYAQKFQG"
                 "RVTMTRDTSTSTVYMELSSLRSEDTAVYYCAR", "WGQGTLVTVSS", 14, 4),
    "IGHV4-34": ("QVQLQQWGAGLLKPSETLSLTCAVYGGSFSGYYWSWIRQPPGKGLEWIGEINHSGSTNYNPSLKS"
                 "RVTISVDTSKNQFSLKLSSVTAADTAVYYCAR", "WGQGTLVTVSS", 15, 5),
    "IGHV3-30": ("QVQLVESGGGVVQPGRSLRLSCAASGFTFSSYAMHWVRQAPGKGLEWVAVISYDGSNKYYADSVKG"
                 "RFTISRDNSKNTLYLQMNSLRAEDTAVYYCAR", "WGQGTLVTVSS", 13, 4),
}

# Padded width of the CDRH3 slot, so sequences from one germline come out
# column-aligned the way IMGT-gapped OAS output is.
CDR3_SLOT = 34

# Approximate IMGT region boundaries within the FR1-FR3 framework strings
# above. Illustrative, like the frameworks themselves -- real OAS units carry
# their own fwr1_aa / cdr1_aa / ... columns from IgBlast and those are used
# in preference to anything here.
REGION_SPLITS = ((0, 25), (25, 33), (33, 50), (50, 58), (58, None))

# The framework strings above end at the conserved cysteine plus the first two
# CDRH3 residues (…CAR / …CAK). _cdr3() already emits those two, so trim them
# here; otherwise the regions overlap and would not concatenate back into the
# sequence, which the alignment's boundary detection relies on.
FRAMEWORKS = {
    g: (fr[:-2] if fr.endswith(("CAR", "CAK")) else fr, fr4, mean, spread)
    for g, (fr, fr4, mean, spread) in FRAMEWORKS.items()
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


def _gap_pad(cdr3: str, width: int = CDR3_SLOT) -> str:
    """Centre-pad a loop with IMGT-style '.' gaps to a fixed column width."""
    if len(cdr3) >= width:
        return cdr3[:width]
    pad = width - len(cdr3)
    mid = (len(cdr3) + 1) // 2
    return cdr3[:mid] + "." * pad + cdr3[mid:]


def make_unit(path: str | Path, n: int = 500, genes: list[str] | None = None,
              species: str = "human", seed: int = 0, liability_rate: float = 0.05,
              clonality: float = 0.45, shm_rate: float = 0.03,
              ambiguous_rate: float = 0.04, gzip_output: bool = True) -> Path:
    """Write a synthetic OAS data unit. Returns the path.

    `clonality` is the chance that a founder CDRH3 spawns an expanded lineage
    of near-identical siblings, as happens in a real campaign.
    `shm_rate` is the per-residue rate of framework substitution, standing in
    for somatic hypermutation so sequences diverge from their germline.
    `ambiguous_rate` is the chance a sequence carries a run of X, as real
    repertoires do where the basecaller could not resolve positions.
    Sequences are emitted with the CDRH3 gap-padded to a fixed column width,
    matching the IMGT-gapped layout of real OAS `sequence_alignment_aa`.
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
               "sequence_alignment_aa", "germline_alignment_aa",
               "v_germline_alignment_aa", "fwr1_aa", "cdr1_aa", "fwr2_aa",
               "cdr2_aa", "fwr3_aa", "cdr3", "cdr3_aa", "junction_aa",
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
            aa = fr + _gap_pad(cdr3) + fr4
            # The germline reference, as IgBlast reports it: unmutated V
            # framework, gaps across the junction (which is not germline
            # encoded), then the J-derived FR4.
            germline_aa = fr + "." * CDR3_SLOT + fr4
            v_germline_aa = fr
            if shm_rate:
                chars = list(aa)
                for _ in range(max(0, int(rng.gauss(shm_rate * len(fr), 1.5)))):
                    i = rng.randrange(len(fr))
                    if chars[i] != ".":
                        chars[i] = rng.choice(_AA)
                aa = "".join(chars)
            # Regions are sliced from the OBSERVED sequence, after SHM, which
            # is what IgBlast reports for real data.
            regions = [aa[:len(fr)][a:b] for a, b in REGION_SPLITS] + [cdr3]
            if rng.random() < ambiguous_rate:
                chars = list(aa)
                start = rng.randrange(0, max(1, len(chars) - 10))
                for i in range(start, min(start + rng.randint(4, 10), len(chars))):
                    if chars[i] != ".":
                        chars[i] = "X"
                aa = "".join(chars)
            liability = "|Unusual residue|" if rng.random() < liability_rate else "||"
            w.writerow([
                "N" * 30, "IGH", "F", "T", "F", "T",
                f"{g}*0{rng.randint(1, 3)}", "IGHD3-10*01",
                f"IGHJ{rng.randint(1, 6)}*0{rng.randint(1, 2)}",
                "N" * 30, aa, germline_aa, v_germline_aa,
                *regions[:5], "NNN", regions[5], "C" + cdr3 + "W",
                round(rng.uniform(0.85, 1.0), 3), round(rng.uniform(0.85, 1.0), 3),
                rng.randint(1, 50), liability, "IGHG",
            ])
    return path
