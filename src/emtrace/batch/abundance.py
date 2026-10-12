"""Cross-sample OTU/feature abundance tables (rows = species, columns = samples).

Reads ``reports/abund_table.tsv`` from each sample subdirectory of a results
directory and pivots it into standard OTU-style tables: one row per species,
one column per sample. Tables are produced per division group:

    bacteria, bacteriophages, viruses,
    invertebrates_parasites, invertebrates_other

The NCBI ``Invertebrates`` division mixes eukaryotic parasites (protozoa and
helminths) with host/environmental animals (ticks, mosquitoes, marine taxa).
It is split here using a curated genus allowlist (``PARASITE_GENERA``) into a
``parasites`` group (veterinary/medical pathogens) and an ``other`` group.

Two value matrices are emitted per group:
  * counts — EM-adjusted read counts (``readn_em``, the EM-expected reads,
             rounded to whole reads); raw ``readn`` where ``readn_em`` is
             missing (tables written before the EM column existed, or taxa
             without an EM estimate)
  * prob   — EMtrace's EM assignment probability (``prob``); this is NOT a
             statistical p-value.

Reuses the per-sample loaders from :mod:`emtrace.batch.taxa`.
"""

from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional

import pandas as pd

from emtrace.batch.taxa import em_read_counts, find_sample_dirs, load_abund_table, _natural_key

logger = logging.getLogger("emtrace.batch.abundance")

# Batch (column) order: samples are named "<batch>_barcodeNN"; columns are
# ordered by this batch order, then natural-sorted by barcode within a batch.
BATCH_ORDER: List[str] = [
    "RUN1",
    "RUN2",
    "RUN3",
    "RUN4.Fecal",
    "RUN4.Wastewater",
    "RUN5.1.Mosquito",
    "RUN5.1.Ticks",
    "RUN5.2.Mosquito",
    "RUN5.2.Mosquito_BU",
    "RUN5.2.Ticks",
    "RUN6.Fecal",
    "RUN6.Wastewater",
    "RUN7.cDNA.Ticks",
    "RUN7.DNA.mosquitos",
    "RUN7.DNA.Ticks.first",
    "RUN7.DNA.Ticks.rerun",
]

# Division (output file / row-block) order.
DIVISION_ORDER: List[str] = [
    "bacteria",
    "bacteriophages",
    "viruses",
    "invertebrates_parasites",
    "invertebrates_other",
]

# Curated genus allowlist: veterinary/medical eukaryotic parasites. Invertebrate
# taxa whose genus is in this set go to "invertebrates_parasites"; everything
# else in the Invertebrates division goes to "invertebrates_other".
PARASITE_GENERA: frozenset = frozenset({
    # Apicomplexa / protozoa
    "Babesia", "Theileria", "Cytauxzoon", "Hepatozoon", "Plasmodium",
    "Hepatocystis", "Trypanosoma", "Leishmania", "Eimeria", "Cystoisospora",
    "Cryptosporidium", "Toxoplasma", "Neospora", "Sarcocystis", "Entamoeba",
    "Giardia",
    # Trematodes
    "Dicrocoelium", "Fasciola", "Schistosoma", "Paramphistomum", "Calicophoron",
    # Nematodes
    "Haemonchus", "Trichostrongylus", "Ostertagia", "Cooperia",
    "Oesophagostomum", "Strongyloides", "Ascaris", "Toxocara", "Trichuris",
    "Bunostomum", "Nematodirus", "Onchocerca", "Dirofilaria", "Setaria",
    "Nippostrongylus",
    # Cestodes
    "Taenia", "Echinococcus", "Moniezia", "Dipylidium", "Hymenolepis",
})


# Taxonomic ranks for the taxpath lineage string, high -> low. "division" is
# excluded on purpose: it is an NCBI GenBank grouping, not a taxonomic rank.
_RANK_COLS = ("class", "order", "family", "genus", "species")


def _first_nonempty(series) -> str:
    for v in series:
        if str(v).strip() not in ("", "nan", "None", "<NA>"):
            return v
    return ""


def _batch_of(sample: str) -> str:
    """Return the batch label of a sample name '<batch>_barcodeNN'."""
    return re.sub(r"_barcode\d+$", "", sample)


def _sample_sort_key(sample: str):
    """Order samples by BATCH_ORDER then natural barcode; unknown batches last."""
    batch = _batch_of(sample)
    try:
        bi = BATCH_ORDER.index(batch)
    except ValueError:
        bi = len(BATCH_ORDER)
    return (bi, _natural_key(batch), _natural_key(sample))


def _group_of(division: str, genus: str, parasite_genera: frozenset) -> Optional[str]:
    """Map a taxon's division (+ genus for Invertebrates) to an output group."""
    d = str(division).strip().lower()
    if d == "bacteria":
        return "bacteria"
    if "phage" in d:
        return "bacteriophages"
    if d == "viruses":
        return "viruses"
    if d == "invertebrates":
        if str(genus).strip() in parasite_genera:
            return "invertebrates_parasites"
        return "invertebrates_other"
    return None


def build_abundance(
    resdir: Path,
    max_samples: Optional[int] = None,
    verbose: bool = False,
    parasite_genera: Optional[Iterable[str]] = None,
) -> Dict[str, Dict[str, pd.DataFrame]]:
    """Build per-division OTU tables.

    Returns ``{group: {"counts": df, "prob": df}}`` where each df has a leading
    ``species`` column followed by one column per sample (ordered by batch).
    Rows are species, ordered by total read count (descending).
    """
    genera = frozenset(parasite_genera) if parasite_genera else PARASITE_GENERA

    sample_dirs = find_sample_dirs(resdir, max_samples)
    if not sample_dirs:
        logger.warning("No valid sample directories found in %s", resdir)
        return {}

    all_samples = [name for name, _ in sample_dirs]
    ordered_samples = sorted(all_samples, key=_sample_sort_key)

    try:
        from rich.progress import track as _track
        iter_samples = _track(sample_dirs, description="Reading samples…")
    except ImportError:
        iter_samples = sample_dirs

    records: List[pd.DataFrame] = []
    for sample_name, sample_dir in iter_samples:
        if verbose:
            logger.info("Loading %s", sample_name)
        df = load_abund_table(sample_dir)
        if df.empty:
            continue
        df = df.copy()
        df["group"] = [
            _group_of(d, g, genera)
            for d, g in zip(df.get("division", ""), df.get("genus", ""))
        ]
        df = df[df["group"].notna()]
        if df.empty:
            continue
        # taxon label: species, falling back to organism when species is blank
        species = df.get("species", "").astype(str).str.strip()
        organism = df.get("organism", "").astype(str).str.strip()
        df["taxon"] = species.where(species != "", organism)
        df = df[df["taxon"] != ""]
        if df.empty:
            continue
        df["Sample"] = sample_name
        df["readn"] = em_read_counts(df, sample_name)
        df["prob"] = pd.to_numeric(df["prob"], errors="coerce")
        for col in _RANK_COLS + ("taxid",):
            if col not in df.columns:
                df[col] = ""
        records.append(
            df[["group", "taxon", "Sample", "readn", "prob", "taxid", *_RANK_COLS]]
        )

    if not records:
        return {}

    long = pd.concat(records, ignore_index=True)

    # taxonomy per taxon: first non-empty taxid + ranks, and a lineage string.
    # "division" is intentionally excluded — it is not an NCBI taxonomic rank.
    taxmap = long.groupby("taxon").agg(
        {c: _first_nonempty for c in ("taxid",) + _RANK_COLS}
    )
    taxmap["taxpath"] = [
        ";".join(v for v in (row[r] for r in _RANK_COLS) if str(v).strip())
        for _, row in taxmap.iterrows()
    ]

    out: Dict[str, Dict[str, pd.DataFrame]] = {}
    for group in DIVISION_ORDER:
        sub = long[long["group"] == group]
        if sub.empty:
            continue

        counts = sub.pivot_table(
            index="taxon", columns="Sample", values="readn",
            aggfunc="sum", fill_value=0,
        )
        prob = sub.pivot_table(
            index="taxon", columns="Sample", values="prob",
            aggfunc="max",
        )

        # all samples as columns (0 / blank where the group is absent), batch order
        counts = counts.reindex(columns=ordered_samples, fill_value=0)
        prob = prob.reindex(columns=ordered_samples).fillna(0.0)

        # rows: species ordered by total read count (descending), then name
        order = (
            counts.sum(axis=1).sort_values(ascending=False, kind="stable").index
        )
        counts = counts.loc[order]
        prob = prob.loc[order]

        counts = counts.reset_index().rename(columns={"taxon": "species"})
        prob = prob.reset_index().rename(columns={"taxon": "species"})
        # taxid as 2nd column, taxpath as last column (both matrices)
        for mat in (counts, prob):
            mat.insert(1, "taxid", mat["species"].map(taxmap["taxid"]))
            mat["taxpath"] = mat["species"].map(taxmap["taxpath"])
        out[group] = {"counts": counts, "prob": prob}

    return out


def write_abundance_xlsx(
    tables: Dict[str, Dict[str, pd.DataFrame]],
    path: Path,
) -> int:
    """Write the read-count tables to a single styled .xlsx (one sheet/division).

    Styling per sheet: the header row is rotated up ("Rotate Text Up", 90°);
    the sample columns (3rd through second-to-last, i.e. excluding species,
    taxid and taxpath) are set to a compact width of 4.0. Returns the number of
    sheets written. Raises ImportError if openpyxl is unavailable.
    """
    from openpyxl.styles import Alignment
    from openpyxl.utils import get_column_letter

    n = 0
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        for group in DIVISION_ORDER:
            if group not in tables:
                continue
            df = tables[group]["counts"]
            sheet = f"{group}.counts"[:31]
            df.to_excel(writer, sheet_name=sheet, index=False)
            ws = writer.sheets[sheet]

            # header row: "Rotate Text Up" (90°)
            for cell in ws[1]:
                cell.alignment = Alignment(textRotation=90, vertical="bottom")
            ws.row_dimensions[1].height = 150

            ncol = df.shape[1]  # species, taxid, <samples...>, taxpath
            ws.column_dimensions[get_column_letter(1)].width = 30  # species
            ws.column_dimensions[get_column_letter(2)].width = 10  # taxid
            for c in range(3, ncol):  # 3rd .. second-last = sample columns
                ws.column_dimensions[get_column_letter(c)].width = 4.0
            ws.column_dimensions[get_column_letter(ncol)].width = 50  # taxpath

            # keep taxonomy columns + header visible while scrolling
            ws.freeze_panes = "C2"
            n += 1
    return n
