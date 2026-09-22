"""Diversity-aware sampling from the pool.

Answers the question that blocked this module: what makes a sampled set
"diverse enough" to stand in for a real discovery campaign. The answer here is
three explicit criteria, checkable after the fact:

  1. germline spread   - every panel gene represented, at its quota
  2. CDRH3 length span - the set covers a range of loop lengths, not one
  3. pairwise identity - no two members are near-duplicates of each other

Selection is greedy farthest-point on CDRH3 within each germline's quota, which
maximises the minimum pairwise distance rather than merely avoiding exact
duplicates.
"""
from __future__ import annotations

import random
from dataclasses import dataclass, field

import pandas as pd

from . import germlines


def levenshtein(a: str, b: str) -> int:
    if a == b:
        return 0
    if len(a) < len(b):
        a, b = b, a
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb)))
        prev = cur
    return prev[-1]


def identity(a: str, b: str) -> float:
    """Normalized sequence identity in [0, 1]. Length-aware, so a short loop and
    a long one are far apart even if the short one is a substring."""
    if not a and not b:
        return 1.0
    longest = max(len(a), len(b))
    return 1.0 - (levenshtein(a, b) / longest) if longest else 1.0


def distance(a: str, b: str) -> float:
    return 1.0 - identity(a, b)


@dataclass(frozen=True)
class Spec:
    """Acceptance criteria for a sampled set."""
    max_pairwise_identity: float = 0.80
    min_cdr3_len_span: int = 8
    min_genes: int = 1
    min_per_gene: int = 1

    def check(self, df: pd.DataFrame) -> dict:
        """Evaluate a sampled set against the criteria. Never raises."""
        cdr3 = list(df["cdr3_aa"])
        worst = 0.0
        worst_pair: tuple[str, str] | None = None
        for i in range(len(cdr3)):
            for j in range(i + 1, len(cdr3)):
                ident = identity(cdr3[i], cdr3[j])
                if ident > worst:
                    worst, worst_pair = ident, (cdr3[i], cdr3[j])
        lens = df["cdr3_len"]
        span = int(lens.max() - lens.min()) if len(lens) else 0
        per_gene = df["v_gene"].value_counts().to_dict()
        results = {
            "max_pairwise_identity": {
                "value": round(worst, 3), "limit": self.max_pairwise_identity,
                "pass": worst <= self.max_pairwise_identity, "worst_pair": worst_pair,
            },
            "cdr3_len_span": {
                "value": span, "limit": self.min_cdr3_len_span,
                "pass": span >= self.min_cdr3_len_span,
            },
            "genes": {
                "value": len(per_gene), "limit": self.min_genes,
                "pass": len(per_gene) >= self.min_genes, "per_gene": per_gene,
            },
            "min_per_gene": {
                "value": min(per_gene.values()) if per_gene else 0,
                "limit": self.min_per_gene,
                "pass": bool(per_gene) and min(per_gene.values()) >= self.min_per_gene,
            },
        }
        results["pass"] = all(v["pass"] for k, v in results.items() if isinstance(v, dict))
        return results


def greedy_farthest(seqs: list[str], n: int, seed: int = 0) -> list[int]:
    """Indices of n sequences chosen to maximise the minimum pairwise distance."""
    if n >= len(seqs):
        return list(range(len(seqs)))
    rng = random.Random(seed)
    first = rng.randrange(len(seqs))
    chosen = [first]
    mind = [distance(seqs[first], s) for s in seqs]
    mind[first] = -1.0
    while len(chosen) < n:
        pick = max(range(len(seqs)), key=lambda i: mind[i])
        chosen.append(pick)
        mind[pick] = -1.0
        for i in range(len(seqs)):
            if mind[i] >= 0:
                mind[i] = min(mind[i], distance(seqs[pick], seqs[i]))
    return chosen


def sample(candidates: pd.DataFrame, n: int, panel: germlines.Panel | None = None,
           spec: Spec | None = None, seed: int = 0,
           length_bins: int = 4) -> pd.DataFrame:
    """Draw n sequences spanning the panel's germlines and a range of CDRH3 lengths.

    `candidates` is a pool index frame (needs v_gene, cdr3_aa, cdr3_len, seq_id).
    """
    panel = panel or germlines.DEFAULT
    spec = spec or Spec()
    if candidates.empty:
        raise ValueError("no candidate sequences")

    quotas = panel.quotas(n)
    picks: list[pd.DataFrame] = []
    shortfall = 0

    for gene, quota in quotas.items():
        pool = candidates[candidates["v_gene"] == gene]
        if pool.empty:
            shortfall += quota
            continue
        quota = min(quota, len(pool))
        # Spread across CDRH3 length bins first, then maximise distance within.
        pool = pool.copy()
        try:
            pool["_bin"] = pd.qcut(pool["cdr3_len"], q=min(length_bins, pool["cdr3_len"].nunique()),
                                   duplicates="drop", labels=False)
        except (ValueError, IndexError):
            pool["_bin"] = 0
        bins = sorted(pool["_bin"].dropna().unique())
        per_bin = [quota // len(bins)] * len(bins)
        for i in range(quota - sum(per_bin)):
            per_bin[i % len(bins)] += 1

        taken: list[pd.DataFrame] = []
        leftover = 0
        for b, k in zip(bins, per_bin):
            sub = pool[pool["_bin"] == b].reset_index(drop=True)
            k = min(k, len(sub))
            leftover += (per_bin[bins.index(b)] - k)
            if k:
                idx = greedy_farthest(list(sub["cdr3_aa"]), k, seed=seed + int(b))
                taken.append(sub.iloc[idx])
        got = pd.concat(taken, ignore_index=True) if taken else pool.head(0)
        if leftover > 0:
            rest = pool[~pool["seq_id"].isin(got["seq_id"])].reset_index(drop=True)
            if len(rest):
                extra = min(leftover, len(rest))
                idx = greedy_farthest(list(rest["cdr3_aa"]), extra, seed=seed)
                got = pd.concat([got, rest.iloc[idx]], ignore_index=True)
        shortfall += quota - len(got)
        picks.append(got)

    out = pd.concat(picks, ignore_index=True) if picks else candidates.head(0)

    if shortfall > 0:
        rest = candidates[~candidates["seq_id"].isin(out["seq_id"])].reset_index(drop=True)
        if len(rest):
            idx = greedy_farthest(list(rest["cdr3_aa"]), min(shortfall, len(rest)), seed=seed)
            out = pd.concat([out, rest.iloc[idx]], ignore_index=True)

    return out.drop(columns=[c for c in out.columns if c.startswith("_")]).reset_index(drop=True)
