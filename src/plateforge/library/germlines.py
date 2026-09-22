"""Germline gene call parsing and the working scaffold panel.

OAS reports V/J genes as IMGT allele calls, sometimes with several
comma-separated ties. Everything that needs a gene or family goes through the
parsers here so the string handling lives in one place.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field

_CALL = re.compile(r"^(IG[HKL][VDJ])(\d+)(?:[-/](\S+?))?(?:\*(\d+))?$")


def split_calls(call: str | None) -> list[str]:
    """OAS ties are comma-separated. Return the individual calls, in order."""
    if not call or not isinstance(call, str):
        return []
    return [c.strip() for c in call.split(",") if c.strip()]


def top_call(call: str | None) -> str | None:
    """The first call of a tie, which is what OAS ranks highest."""
    calls = split_calls(call)
    return calls[0] if calls else None


def gene(call: str | None) -> str | None:
    """'IGHV3-23*01' -> 'IGHV3-23'. Allele stripped, ties resolved to the first."""
    top = top_call(call)
    if not top:
        return None
    return top.split("*", 1)[0]


def family(call: str | None) -> str | None:
    """'IGHV3-23*01' -> 'IGHV3'."""
    top = top_call(call)
    if not top:
        return None
    m = _CALL.match(top)
    return f"{m.group(1)}{m.group(2)}" if m else None


def allele(call: str | None) -> str | None:
    """'IGHV3-23*01' -> '01'. None when the call carries no allele."""
    top = top_call(call)
    if not top or "*" not in top:
        return None
    return top.split("*", 1)[1]


def chain_of(call: str | None) -> str | None:
    """'IGHV3-23*01' -> 'H'; IGKV/IGLV -> 'L'."""
    top = top_call(call)
    if not top or len(top) < 3:
        return None
    return "H" if top[2] == "H" else "L"


@dataclass(frozen=True)
class Scaffold:
    """One germline in the working panel, with why it is there."""
    gene: str
    weight: float
    role: str
    rationale: str


@dataclass(frozen=True)
class Panel:
    name: str
    scaffolds: tuple[Scaffold, ...]
    light_chain: str | None = None
    notes: str = ""

    @property
    def genes(self) -> list[str]:
        return [s.gene for s in self.scaffolds]

    def weight_of(self, g: str) -> float:
        for s in self.scaffolds:
            if s.gene == g:
                return s.weight
        return 0.0

    def quotas(self, n: int) -> dict[str, int]:
        """Split n across the panel by weight, largest-remainder, every gene >= 1."""
        if n < len(self.scaffolds):
            raise ValueError(f"n={n} is smaller than the {len(self.scaffolds)}-gene panel")
        total = sum(s.weight for s in self.scaffolds)
        raw = {s.gene: (s.weight / total) * n for s in self.scaffolds}
        floors = {g: max(1, int(v)) for g, v in raw.items()}
        short = n - sum(floors.values())
        order = sorted(raw, key=lambda g: raw[g] - int(raw[g]), reverse=True)
        i = 0
        while short > 0:
            floors[order[i % len(order)]] += 1
            short -= 1
            i += 1
        while short < 0:
            g = max((g for g in floors if floors[g] > 1), key=lambda g: floors[g], default=None)
            if g is None:
                break
            floors[g] -= 1
            short += 1
        return floors

    def as_dict(self) -> dict:
        return {
            "name": self.name,
            "light_chain": self.light_chain,
            "notes": self.notes,
            "scaffolds": [
                {"gene": s.gene, "weight": s.weight, "role": s.role, "rationale": s.rationale}
                for s in self.scaffolds
            ],
        }

    @classmethod
    def from_dict(cls, d: dict) -> "Panel":
        return cls(
            name=d["name"],
            scaffolds=tuple(
                Scaffold(s["gene"], float(s["weight"]), s.get("role", ""), s.get("rationale", ""))
                for s in d["scaffolds"]
            ),
            light_chain=d.get("light_chain"),
            notes=d.get("notes", ""),
        )


PANELS: dict[str, Panel] = {}


def register_panel(panel: Panel) -> Panel:
    if panel.name in PANELS:
        raise KeyError(f"panel {panel.name!r} already registered")
    PANELS[panel.name] = panel
    return panel


def get_panel(name: str) -> Panel:
    try:
        return PANELS[name]
    except KeyError:
        raise KeyError(f"unknown panel {name!r}; known: {sorted(PANELS)}") from None


# --- the working panel. Reasoning is in decisions/0007-germline-panel.md ----
DEFAULT = register_panel(Panel(
    name="default-v1",
    light_chain="IGKV1-39",
    notes="Human heavy chains only. Fixed kappa light chain for reformatting.",
    scaffolds=(
        Scaffold("IGHV3-23", 0.60, "workhorse",
                 "Most common human VH and the germline behind many approved "
                 "therapeutics; the well-behaved baseline."),
        Scaffold("IGHV1-69", 0.25, "long-hydrophobic-cdrh3",
                 "Dominant in antiviral responses, biased toward long hydrophobic "
                 "CDRH3s; supplies clones that behave unlike the baseline."),
        Scaffold("IGHV3-53", 0.15, "short-cdrh3",
                 "Short-CDRH3 public response germline; the third canonical "
                 "shape, keeping the panel from being one loop geometry."),
    ),
))

SINGLE = register_panel(Panel(
    name="single-v1",
    light_chain="IGKV1-39",
    notes="Fastest path: one germline. Use when downstream behavioral variety is not needed yet.",
    scaffolds=(
        Scaffold("IGHV3-23", 1.0, "workhorse", "See default-v1."),
    ),
))
