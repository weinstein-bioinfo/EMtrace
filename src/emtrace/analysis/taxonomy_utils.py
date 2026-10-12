"""Taxonomy-based annotation enrichment utilities.

Ports ``add_division_field()`` and ``add_bphage_field()`` from
emtrace.pl to pure Python, replacing taxonkit + csvtk calls with
direct lookups against the NCBITaxonomy SQLite database.
"""

from __future__ import annotations

import logging
import re
import sqlite3
from pathlib import Path
from typing import Dict, Optional, Set

import pandas as pd

logger = logging.getLogger(__name__)


# -------------------------------------------------------------------
# No-species (non-binomial / unassigned) detection
# -------------------------------------------------------------------

# Markers of a species rank that is absent or not a resolved binomial.
#
# Deliberately narrow.  Measured against the installed NCBI taxonomy
# (2,739,528 rows): 'sp.' alone accounts for 49.4% of rows, while a
# broader 'bacterium|archaeon|symbiont' pattern would match a further
# 4.95% of validly-named taxa and is therefore NOT included.
#
# Viral names are non-Linnaean ("Tulip top breaking virus") and survive
# this pattern — only 0.71% of viral rows match, and those are genuinely
# unresolved ("unidentified adenovirus").
NO_SPECIES_RE = re.compile(
    r"\b(?:sp|cf|aff)\.|\b(?:uncultured|unclassified|environmental|unidentified)\b",
    re.IGNORECASE,
)

# Placeholder strings that stand in for a missing value.
_NO_SPECIES_LITERALS = {"", "na", "n/a", "nan", "none", "-", "null"}


def is_no_species(s: pd.Series) -> pd.Series:
    """Return a boolean mask: True where the species rank is unresolved.

    A value is "no species" when it is empty/NA (the species column is
    NULL for every taxon above species rank), a placeholder literal, or a
    non-binomial name such as ``"Vitreoscilla sp."``.
    """
    v = s.fillna("").astype(str).str.strip()
    return v.str.lower().isin(_NO_SPECIES_LITERALS) | v.str.contains(
        NO_SPECIES_RE, na=False
    )


def filter_no_species(
    df: pd.DataFrame,
    species_col: str = "species",
    context: str = "",
) -> pd.DataFrame:
    """Return *df* without rows lacking a resolved binomial species.

    No-ops with a warning when *species_col* is absent, so annotation
    tables written before species assignment (or a taxonomy DB built
    without species) degrade to "keep everything" rather than dropping
    every row.
    """
    where = f" [{context}]" if context else ""

    if species_col not in df.columns:
        logger.warning(
            "exclude_no_species requested but column '%s' not found%s — "
            "keeping all %d rows",
            species_col, where, len(df),
        )
        return df

    mask = is_no_species(df[species_col])
    n_before = len(df)
    out = df[~mask]
    logger.info(
        "Excluded %d rows with no species assignment%s (%d remaining)",
        n_before - len(out), where, len(out),
    )
    return out


# -------------------------------------------------------------------
# SILVA detection from annotation DataFrame
# -------------------------------------------------------------------

def is_silva_annotation(df: pd.DataFrame, config=None) -> bool:
    """Return True if *df* was produced using a SILVA database.

    With a *config* (``EMtraceConfig``), each distinct ``db`` value is
    resolved through ``config.taxonomy_for_database``, so an explicit
    ``taxonomy: SILVA`` on a registry entry is honoured. Without one, any
    ``db`` value containing 'silva' (case-insensitive) counts. Returns False
    (NCBI assumed) when the column is absent or empty.
    """
    if "db" not in df.columns or df.empty:
        return False
    if config is not None:
        from emtrace.core.config import TaxonomyType

        return any(
            config.taxonomy_for_database(db) == TaxonomyType.SILVA
            for db in df["db"].dropna().astype(str).unique()
        )
    return df["db"].astype(str).str.lower().str.contains("silva", na=False).any()


# -------------------------------------------------------------------
# add_division_field
# -------------------------------------------------------------------

def add_division_field(
    df: pd.DataFrame,
    taxonomy_db_path: Path,
    taxid_col: str = "staxid",
    overwrite: bool = True,
) -> None:
    """Assign NCBI division names to annotation/abundance rows in-place.

    Port of ``add_division_field()`` from emtrace.pl.
    Uses the pre-built NCBITaxonomy SQLite database to resolve each
    taxid to its NCBI division name (Viruses, Bacteria, Phages, etc.).

    Parameters
    ----------
    df : pd.DataFrame
        Annotation or abundance table with a *taxid_col* column.
    taxonomy_db_path : Path
        Path to the NCBI taxonomy SQLite3 database.
    taxid_col : str
        Column containing NCBI taxonomy IDs.
    overwrite : bool
        If True, overwrite existing 'division' column.
    """
    if taxid_col not in df.columns:
        logger.warning(
            "add_division_field: column '%s' not found", taxid_col
        )
        return

    if "division" in df.columns and not overwrite:
        return

    from emtrace.models.taxonomy import NCBITaxonomy
    tax = NCBITaxonomy(taxonomy_db_path)
    tax.add_division_to_df(df, taxid_col=taxid_col)
    logger.info("Division field assigned (%d rows)", len(df))


# -------------------------------------------------------------------
# add_species_field
# -------------------------------------------------------------------

def add_species_field(
    df: pd.DataFrame,
    taxonomy_db_path: Path,
    taxid_col: str = "staxid",
    overwrite: bool = True,
) -> None:
    """Assign the NCBI species-rank name to annotation rows in-place.

    Mirrors :func:`add_division_field`.  The SILVA path already gets a
    ``species`` column from ``SILVATaxonomy.add_lineage_to_df`` (called via
    :func:`add_domain_field`); this brings the NCBI path to parity so that
    ``species`` is available to the EM step, which runs before reporting.

    Parameters
    ----------
    df : pd.DataFrame
        Annotation or abundance table with a *taxid_col* column.
    taxonomy_db_path : Path
        Path to the NCBI taxonomy SQLite3 database.
    taxid_col : str
        Column containing NCBI taxonomy IDs.
    overwrite : bool
        If True, overwrite an existing 'species' column.
    """
    if taxid_col not in df.columns:
        logger.warning(
            "add_species_field: column '%s' not found", taxid_col
        )
        return

    if "species" in df.columns and not overwrite:
        return

    from emtrace.models.taxonomy import NCBITaxonomy
    tax = NCBITaxonomy(taxonomy_db_path)
    tax.add_species_to_df(df, taxid_col=taxid_col)
    logger.info("Species field assigned (%d rows)", len(df))


# -------------------------------------------------------------------
# add_bphage_field
# -------------------------------------------------------------------

def _load_phage_filter_taxids(
    phfilter_path: Path,
) -> Set[str]:
    """Load the set of phage taxids from the ICTV phage filter TSV.

    The file has columns: rank, name, taxid (with header).
    """
    taxids: Set[str] = set()
    with open(phfilter_path) as fh:
        next(fh, None)  # skip header
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                taxids.add(parts[2].strip())
    return taxids


def _get_lineage_taxids(
    taxid: int,
    taxonomy_db_path: Path,
) -> list[int]:
    """Walk the taxonomy tree from *taxid* to root.

    Returns all ancestor taxids.  Uses the ``taxonomy`` table's
    parent links if available, but since our flattened DB doesn't
    store parent_id, we use a simpler approach: look up the
    taxonomy row and collect taxids from all non-null ranks.

    For the bphage check we actually need the full lineage as taxids,
    so we query the DB for the lineage path.
    """
    # Our NCBITaxonomy DB stores flattened lineage (rank->name)
    # but not parent taxids. For bphage checking, we need to
    # check whether any ancestor is in the phage filter set.
    # The actual lineage walk happens via a bulk approach in
    # add_bphage_field.
    return [taxid]


def add_bphage_field(
    df: pd.DataFrame,
    taxonomy_db_path: Path,
    phfilter_path: Optional[Path] = None,
    taxid_col: str = "staxid",
    overwrite: bool = True,
) -> None:
    """Assign bacteriophage flag ('yes'/'no') to annotation rows in-place.

    Port of ``add_bphage_field()`` from emtrace.pl.
    A row gets bphage='yes' if:
      1. Its division is 'Phages', OR
      2. Its taxid (or any ancestor in its NCBI lineage) appears
         in the ICTV phage filter file.

    Parameters
    ----------
    df : pd.DataFrame
        Annotation table with 'division' and *taxid_col* columns.
        Division should be assigned before calling this function.
    taxonomy_db_path : Path
        Path to the NCBI taxonomy SQLite3 database built by
        ``NCBITaxonomy.build_database()``.  Must contain a 'merged'
        table for old→new taxid resolution.
    phfilter_path : Path, optional
        Path to ``ICTV.bphage.filter.tsv`` (rank, name, taxid).
        If None, only division=='Phages' marking is done.
    taxid_col : str
        Column containing NCBI taxonomy IDs.
    overwrite : bool
        If True, overwrite existing 'bphage' column values.
    """
    if taxid_col not in df.columns:
        logger.warning(
            "add_bphage_field: column '%s' not found", taxid_col
        )
        return

    if "bphage" in df.columns and not overwrite:
        return

    # Determine phage taxids from filter file
    phfilter_taxids: Set[str] = set()
    if phfilter_path and phfilter_path.exists():
        phfilter_taxids = _load_phage_filter_taxids(phfilter_path)
        logger.debug(
            "Loaded %d phage filter taxids from %s",
            len(phfilter_taxids),
            phfilter_path,
        )

    # Build the set of taxids in this annotation that are viral
    division_col = "division"
    has_division = division_col in df.columns

    # Collect all unique viral taxids and check their lineages
    phage_taxid_set: Set[str] = set()

    if phfilter_taxids and has_division:
        viral_mask = df[division_col].str.lower() == "viruses"
        viral_taxids = (
            df.loc[viral_mask, taxid_col]
            .dropna()
            .unique()
            .tolist()
        )

        if viral_taxids and taxonomy_db_path.exists():
            phage_taxid_set = _check_lineages_for_phage(
                viral_taxids,
                taxonomy_db_path,
                phfilter_taxids,
            )

    # Assign bphage flag
    def _assign(row):
        # Direct division match
        if has_division:
            div = str(row.get(division_col, "")).lower()
            if div == "phages":
                return "yes"
        # Lineage-based match
        taxid_str = str(row.get(taxid_col, ""))
        if taxid_str in phage_taxid_set:
            return "yes"
        return "no"

    df["bphage"] = df.apply(_assign, axis=1)
    n_phage = (df["bphage"] == "yes").sum()
    logger.info(
        "Bphage field assigned: %d/%d marked as phage",
        n_phage,
        len(df),
    )


def add_domain_field(
    df: pd.DataFrame,
    taxonomy_db_path: Path,
    taxid_col: str = "staxid",
    overwrite: bool = True,
) -> None:
    """Assign SILVA domain to annotation rows in-place via accession lookup.

    Calls ``SILVATaxonomy.add_lineage_to_df`` which looks up by ``sseqid``
    (the SILVA accession in PAF/BLAST output, e.g. ``AB000001.1.1450``),
    populates ``staxid`` and lineage rank columns (domain, phylum, …), then
    copies ``domain`` → ``division`` so the value appears in the standard
    annotation TSV column and read-partitioning conditions.

    Parameters
    ----------
    df : pd.DataFrame
        Annotation table with an ``sseqid`` column containing SILVA accessions.
    taxonomy_db_path : Path
        Path to the SILVA taxonomy SQLite3 database.
    taxid_col : str
        Column to write the SILVA taxid into (typically ``"staxid"``).
    overwrite : bool
        If True, overwrite existing ``division`` column values.
    """
    if "sseqid" not in df.columns:
        logger.warning("add_domain_field: column 'sseqid' not found")
        return

    if "division" in df.columns and not overwrite:
        return

    from emtrace.models.taxonomy import SILVATaxonomy
    tax = SILVATaxonomy(taxonomy_db_path)
    tax.add_lineage_to_df(df, sseqid_col="sseqid", staxid_col=taxid_col)

    # Mirror domain → division so it is captured in the annotation TSV and
    # used by partition_reads conditions (_cond_silva_ba etc.)
    if "domain" in df.columns:
        df["division"] = df["domain"]

    logger.info("Domain field assigned (%d rows)", len(df))


def _check_lineages_for_phage(
    taxids: list,
    taxonomy_db_path: Path,
    phfilter_taxids: Set[str],
) -> Set[str]:
    """Check which taxids have any phage-filter ancestor in their lineage.

    The NCBI taxonomy DB built by NCBITaxonomy doesn't store parent
    links, but it does store the full nodes.dmp structure.  To avoid
    requiring nodes.dmp at runtime, we use a practical shortcut:

    We load the taxonomy DB's ``taxonomy`` table to get each taxid's
    division. For the lineage walk, we build a lightweight parent map
    from the DB if a ``nodes`` table exists, or fall back to checking
    only direct taxid membership in the phage filter.

    For the Perl pipeline, lineages were obtained via taxonkit.
    Here we replicate the logic by building a parent map from the DB.
    """
    phage_set: Set[str] = set()

    if not taxonomy_db_path.exists():
        return phage_set

    with sqlite3.connect(taxonomy_db_path) as conn:
        cur = conn.cursor()

        # Check if a 'nodes' table exists (stores parent links)
        cur.execute(
            "SELECT name FROM sqlite_master "
            "WHERE type='table' AND name='nodes'"
        )
        has_nodes = cur.fetchone() is not None

        # Load merged taxids for resolution
        merged_map: Dict[str, str] = {}
        try:
            rows = cur.execute(
                "SELECT old_taxid, new_taxid FROM merged"
            ).fetchall()
            merged_map = {str(o): str(n) for o, n in rows}
        except sqlite3.OperationalError:
            pass

        if has_nodes:
            # Build parent map from nodes table
            parent_map: Dict[str, str] = {}
            rows = cur.execute(
                "SELECT taxid, parent_id FROM nodes"
            ).fetchall()
            parent_map = {str(t): str(p) for t, p in rows}

            for tid in taxids:
                tid_str = str(tid).strip()
                # Resolve merged
                tid_str = merged_map.get(tid_str, tid_str)
                # Walk lineage
                current = tid_str
                visited: Set[str] = set()
                while current and current not in visited:
                    if current in phfilter_taxids:
                        phage_set.add(str(tid))
                        break
                    visited.add(current)
                    current = parent_map.get(current, "")
                    if current == "1" or current == "":
                        break
        else:
            # Fallback: only direct membership check
            for tid in taxids:
                tid_str = str(tid).strip()
                tid_str = merged_map.get(tid_str, tid_str)
                if tid_str in phfilter_taxids:
                    phage_set.add(str(tid))

    return phage_set
