"""Sanger sequencing order forms, one emitter per vendor.

These are the first emitters in this repo that are actually `verified`: the
real template from each vendor is in `fixtures/sequencing/`, and the tests
read the template and check that what we write fits it. Nothing here is a
guess at a layout.

Three vendors, and they disagree about nearly everything a plate map is:

| | well spelling | fill order | shape |
|---|---|---|---|
| GENEWIZ / Azenta | `A01` | gives you both, as two columns | one long table, 500 rows |
| ELIM Biopharm | `A1` | row-major (A1..A12, B1..) | one table, 96 rows, dropdowns |
| UC Berkeley | `A1` | column-major (A1, B1, .. H1, A2) | printed form, two side-by-side blocks |

That table is the entire argument for rule 3. A plate laid out column-major
and pasted into a row-major form is a 96-well transposition that every later
step faithfully preserves, and nothing downstream can detect. So every
emitter takes wells in our canonical `A01` form and converts on the way out,
once, in one place.

Azenta's form is the one that can go wrong quietly, because it offers `Well
(H)` and `Well (V)` as separate columns and the sample numbering runs down
the rows either way. Our plates are column-major (`library.diversity.
group_for_plate`), so `Well (V)` is the column that lines up -- and the
emitter writes both, so a human can see which it used.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from ..core import registry, wells as wellmod
from .worklists import SKETCH, VERIFIED

FORMS = registry.Registry("sequencing order form")

# What each vendor's own template says it requires. Checked before writing,
# because a form rejected at upload costs a day and the check costs nothing.
LIMITS = {
    "genewiz": {"max_samples": 500, "primers_per_row": 1,
                "primer_separator": ";"},
    "elim": {"max_samples": 96, "name_chars": 50,
             "allowed": "letters, numbers, dash and underscore only"},
    "ucberkeley": {"max_samples": 95, "must_leave_empty": 1},
}


@dataclass
class OrderForm:
    """A filled order form, and what qualifies it."""
    vendor: str
    frame: pd.DataFrame
    suffix: str
    confidence: str
    template: str = ""
    n_samples: int = 0
    notes: list[str] = field(default_factory=list)
    issues: list[str] = field(default_factory=list)

    @property
    def verified(self) -> bool:
        return self.confidence == VERIFIED

    def write(self, directory: str | Path, stem: str = "sequencing_order") -> Path:
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)
        path = directory / f"{stem}_{self.vendor}{self.suffix}"
        if self.suffix == ".csv":
            self.frame.to_csv(path, index=False)
        else:
            try:
                with pd.ExcelWriter(path, engine="openpyxl") as writer:
                    self.frame.to_excel(writer, index=False, sheet_name="order")
            except ImportError as exc:
                raise ImportError(
                    f"writing the {self.vendor} form needs openpyxl: "
                    "pip install -e '.[sequencing]'. The 'generic' form is a "
                    "CSV and needs nothing.") from exc
        if self.notes or self.issues:
            path.with_suffix(".NOTES.txt").write_text(
                f"{self.vendor} sequencing order — {self.n_samples} samples\n"
                f"format confidence: {self.confidence}\n"
                + (f"template: {self.template}\n" if self.template else "")
                + "\n" + "\n".join(f"- {n}" for n in self.notes)
                + ("\n\nissues:\n" + "\n".join(f"  ! {i}" for i in self.issues)
                   if self.issues else "") + "\n")
        return path


def _clean(frame: pd.DataFrame) -> pd.DataFrame:
    """Normalise wells once, on the way in. Nothing downstream re-parses."""
    out = frame.copy()
    out["well"] = out["well"].map(lambda w: wellmod.normalize(str(w)))
    if "sample_name" not in out.columns:
        out["sample_name"] = out.get("clone_id", pd.Series(dtype=str))
    return out.sort_values("well").reset_index(drop=True)


def _short(well: str) -> str:
    """`A01` -> `A1`, for the vendors whose forms use it."""
    return f"{well[0]}{int(well[1:])}"


# --- GENEWIZ / Azenta --------------------------------------------------------

@FORMS.register("genewiz", vendor="Azenta Life Sciences (GENEWIZ)",
                confidence=VERIFIED, suffix=".xlsx",
                template="fixtures/sequencing/azenta_sanger_form_v2.xlsx")
def genewiz(frame: pd.DataFrame, primer: str = "", **kwargs) -> OrderForm:
    """The Azenta "Excel Form Version 2" layout.

    Column headers are reproduced exactly, including the spacing, because the
    template says in its own notes: DO NOT change column headers. Columns may
    be reordered and extra columns are ignored on upload, so the order here is
    the template's for readability rather than out of necessity.

    `Well (V)` is filled with the actual well; `Well (H)` is written too, so
    that a human opening the file can see the plate both ways and catch a
    transposition before it is submitted.
    """
    data = _down_the_plate(_clean(frame))
    limits = LIMITS["genewiz"]
    issues, notes = [], []

    primers = data.get("primer", pd.Series([primer] * len(data))).fillna(primer)
    if any(primers.astype(str).str.contains(";")):
        notes.append("rows with several primers use ';' as the template "
                     "requires; pre-mixed submissions allow only one primer "
                     "per row.")
    if len(data) > limits["max_samples"]:
        issues.append(f"{len(data)} samples; the form holds "
                      f"{limits['max_samples']}")

    # Both columns are fixed orderings of all 96 wells, paired row by row:
    # row i is the ith well across (H) and the ith well down (V). They are
    # properties of the *form*, not of our samples. Our plate runs down the
    # columns, so the sample in form row i must be the one in `Well (V)` --
    # which is asserted rather than assumed, because getting this backwards
    # transposes the whole plate and nothing downstream could tell.
    across = [wellmod.from_rc(i // 12 + 1, i % 12 + 1) for i in range(len(data))]
    down = [wellmod.from_rc(i % 8 + 1, i // 8 + 1) for i in range(len(data))]
    misplaced = [(a, b) for a, b in zip(down, data["well"]) if a != b]
    if misplaced:
        issues.append(
            f"{len(misplaced)} sample(s) are not in the wells the form's "
            f"vertical ordering expects (first: form says {misplaced[0][0]}, "
            f"sample is in {misplaced[0][1]}). The plate is not contiguous "
            "column-major; fill the form by well rather than by row order.")

    out = pd.DataFrame({
        "Well (H)": across,
        "Well (V)": down,
        "Sample #": range(1, len(data) + 1),
        "DNA Name": data["sample_name"].astype(str),
        "Length (bp)": data.get("length_bp", pd.Series([""] * len(data))),
        "Concentration (ng/uL)": data.get("concentration_ng_ul",
                                          pd.Series([""] * len(data))),
        "Primer": primers.astype(str),
        "Difficult Template": data.get("difficult", pd.Series([""] * len(data))),
        "Notes": data.get("notes", pd.Series([""] * len(data))),
        "Is Empty": ["" for _ in range(len(data))],
    })
    notes.append("plate is column-major, so the sample order follows "
                 "'Well (V)'. Check that before uploading.")
    return OrderForm("genewiz", out, ".xlsx", VERIFIED,
                     template=FORMS.meta("genewiz")["template"],
                     n_samples=len(data), notes=notes, issues=issues)


def _down_the_plate(frame: pd.DataFrame) -> pd.DataFrame:
    """Order rows the way the plate is filled: down each column, then across."""
    out = frame.copy()
    out["_order"] = [wellmod.to_rc(w)[1] * 100 + wellmod.to_rc(w)[0]
                     for w in out["well"]]
    return out.sort_values("_order").drop(columns="_order").reset_index(drop=True)


# --- ELIM Biopharm -----------------------------------------------------------

@FORMS.register("elim", vendor="ELIM Biopharmaceuticals",
                confidence=VERIFIED, suffix=".xlsx",
                template="fixtures/sequencing/elim_seq_orderform_96well.xls")
def elim(frame: pd.DataFrame, primer: str = "", template_type: str = "Plasmid",
         premix: str = "Premix", **kwargs) -> OrderForm:
    """ELIM's 96-well form: row-major wells, `A1`, and controlled vocabularies.

    `Template Type`, `Premix?` and `GC Rich?` are dropdowns in the real file
    with fixed option lists (parked in columns P-R of the template). Writing
    anything outside them is rejected, so they are validated here rather than
    discovered at upload.
    """
    data = _clean(frame)
    limits = LIMITS["elim"]
    issues, notes = [], []

    if template_type not in TEMPLATE_TYPES:
        issues.append(f"template type {template_type!r} is not one of "
                      f"{', '.join(TEMPLATE_TYPES)}")
    if premix not in PREMIX:
        issues.append(f"premix {premix!r} is not one of {', '.join(PREMIX)}")
    if len(data) > limits["max_samples"]:
        issues.append(f"{len(data)} samples; one sheet holds "
                      f"{limits['max_samples']} — the template says to add a "
                      "sheet per extra plate")

    long_names = data[data["sample_name"].astype(str).str.len() > limits["name_chars"]]
    if len(long_names):
        issues.append(f"{len(long_names)} sample name(s) over "
                      f"{limits['name_chars']} characters")
    bad = data[~data["sample_name"].astype(str).str.fullmatch(r"[A-Za-z0-9_-]+")]
    if len(bad):
        issues.append(f"{len(bad)} sample name(s) use characters outside "
                      f"{limits['allowed']}: "
                      + ", ".join(bad["sample_name"].astype(str).head(3)))

    # Row-major: the form's rows run A1..A12, B1..B12.
    ordered = data.copy()
    ordered["_order"] = [wellmod.to_rc(w)[0] * 100 + wellmod.to_rc(w)[1]
                         for w in ordered["well"]]
    ordered = ordered.sort_values("_order").drop(columns="_order")

    out = pd.DataFrame({
        "Reaction #": range(1, len(ordered) + 1),
        "Well": [_short(w) for w in ordered["well"]],
        "Template Name*": ordered["sample_name"].astype(str),
        "Template Type": template_type,
        "Template Conc (ng/ul)": ordered.get("concentration_ng_ul",
                                             pd.Series([""] * len(ordered))),
        "Template Size (bp)*": ordered.get("length_bp",
                                           pd.Series([""] * len(ordered))),
        "Primer Name*": ordered.get("primer", pd.Series([primer] * len(ordered))
                                    ).fillna(primer).astype(str),
        "Primer Conc (uM)": ordered.get("primer_um", pd.Series([""] * len(ordered))),
        "Premix?*": premix,
        "Preferred Annealling Temp. (oC)": ordered.get(
            "anneal_c", pd.Series([""] * len(ordered))),
        "GC Rich?": ordered.get("gc_rich", pd.Series(["No"] * len(ordered))),
        "Notes": ordered.get("notes", pd.Series([""] * len(ordered))),
    })
    notes.append("wells are row-major (A1..A12, then B1) and spelled 'A1', "
                 "which is the template's convention, not ours.")
    return OrderForm("elim", out, ".xlsx", VERIFIED,
                     template=FORMS.meta("elim")["template"],
                     n_samples=len(ordered), notes=notes, issues=issues)


TEMPLATE_TYPES = ("Plasmid", "PCR Product", "BAC DNA", "Genomic DNA", "Cosmid",
                  "Other (e.g. Bacteria)")
PREMIX = ("Premix", "Non-Premix")


# --- UC Berkeley DNA Sequencing Facility -------------------------------------

@FORMS.register("ucberkeley", vendor="UC Berkeley DNA Sequencing Facility",
                confidence=VERIFIED, suffix=".xlsx",
                template="fixtures/sequencing/ucberkeley_full_plate_order_form.xlsx")
def ucberkeley(frame: pd.DataFrame, primer: str = "", **kwargs) -> OrderForm:
    """The Berkeley full-plate form: column-major, `A1`, two side-by-side blocks.

    The template says "LEAVE AT LEAST ONE WELL EMPTY ON PLATE" in the A1 cell,
    which is a constraint rather than a suggestion -- it is how the facility
    confirms plate orientation. A full 96-sample plate is therefore an error
    here, and this refuses rather than letting it be found on the bench.
    """
    data = _clean(frame)
    limits = LIMITS["ucberkeley"]
    issues, notes = [], []

    if len(data) > limits["max_samples"]:
        issues.append(
            f"{len(data)} samples on a 96-well plate, but the facility "
            "requires at least one well left empty so they can confirm the "
            "plate's orientation. Drop one sample or split the plate.")

    filled = {w: n for w, n in zip(data["well"], data["sample_name"].astype(str))}
    left, right = [], []
    for column in range(1, 13):
        for row in range(1, 9):
            well = wellmod.from_rc(row, column)
            entry = (_short(well), filled.get(well, ""), "")
            (left if column <= 6 else right).append(entry)

    out = pd.DataFrame({
        "well #": [a for a, _, _ in left],
        "SAMPLE NAME": [b for _, b, _ in left],
        "SPECIAL INSTRUCTION": [c for _, _, c in left],
        "well # (7-12)": [a for a, _, _ in right],
        "SAMPLE NAME (7-12)": [b for _, b, _ in right],
        "SPECIAL INSTRUCTION (7-12)": [c for _, _, c in right],
    })
    notes.append("two side-by-side blocks: columns 1-6 on the left, 7-12 on "
                 "the right, wells running down each column. Paste into the "
                 "facility's own form, which also needs the PI, chartstring "
                 "and plate barcode in its header block.")
    notes.append(f"{96 - len(data)} well(s) left empty.")
    if primer:
        notes.append(f"primer {primer} is supplied with the plate; this form "
                     "has no primer column.")
    return OrderForm("ucberkeley", out, ".xlsx", VERIFIED,
                     template=FORMS.meta("ucberkeley")["template"],
                     n_samples=len(data), notes=notes, issues=issues)


# --- a plain table, which cannot be wrong about anyone else's format ---------

@FORMS.register("generic", vendor="none", confidence=VERIFIED, suffix=".csv",
                template="")
def generic(frame: pd.DataFrame, primer: str = "", **kwargs) -> OrderForm:
    data = _clean(frame)
    out = pd.DataFrame({
        "well": data["well"],
        "sample_name": data["sample_name"].astype(str),
        "primer": data.get("primer", pd.Series([primer] * len(data))).fillna(primer),
        "length_bp": data.get("length_bp", pd.Series([""] * len(data))),
        "concentration_ng_ul": data.get("concentration_ng_ul",
                                        pd.Series([""] * len(data))),
        "clone_id": data.get("clone_id", pd.Series([""] * len(data))),
    })
    return OrderForm("generic", out, ".csv", VERIFIED, n_samples=len(data))


DEFAULT = "genewiz"


def order(frame: pd.DataFrame, vendor: str = DEFAULT, *, primer: str = "",
          allow_issues: bool = False, **kwargs) -> OrderForm:
    """Fill one vendor's order form from a table of wells and sample names.

    Refuses when the form's own stated limits are broken, unless
    `allow_issues=True` -- the point is to fail here rather than at upload.
    """
    form = FORMS.get(vendor)(frame, primer=primer, **kwargs)
    if form.issues and not allow_issues:
        raise ValueError(
            f"{vendor} order form: " + "; ".join(form.issues)
            + ". Fix these, or pass allow_issues=True to write it anyway.")
    return form


def available() -> pd.DataFrame:
    """Every form, and how well its layout is actually known."""
    rows = []
    for name in FORMS:
        meta = FORMS.meta(name)
        rows.append({"form": name, "vendor": meta.get("vendor", ""),
                     "confidence": meta.get("confidence", SKETCH),
                     "suffix": meta.get("suffix", ""),
                     "template": meta.get("template", ""),
                     "default": name == DEFAULT})
    return pd.DataFrame(rows).sort_values("form").reset_index(drop=True)
