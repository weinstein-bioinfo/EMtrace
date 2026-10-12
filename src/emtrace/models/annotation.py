"""Annotation table schema.

Defines the 18-field annotation table used throughout emtrace.
Both a Pydantic model (for row-level validation) and DataFrame helpers
(for bulk operations) are provided.

Fields (in order):
    search, db, qseqid, orf, qseqlen, sseqid, bitscore, evalue,
    alen, pident, qlen, qcov, slen, scov, staxid, sname, division
"""

from __future__ import annotations

import csv
from typing import List, Optional

import pandas as pd
from pydantic import BaseModel, Field, field_validator

# Ordered field list — kept in sync with ANNOT_TABLE_FIELDS in SeqAn.pm
ANNOT_FIELDS: List[str] = [
    "search",    # search program: minimap2, blastn, blastx, diamond.blastx, diamond.blastp
    "db",        # database basename
    "qseqid",    # query sequence ID
    "orf",       # ORF ID extracted from query header (or "NA")
    "qseqlen",   # full query sequence length
    "sseqid",    # subject/target sequence ID
    "bitscore",  # alignment bitscore
    "evalue",    # E-value
    "alen",      # alignment length
    "pident",    # percent identity
    "qlen",      # query alignment region length
    "qcov",      # query coverage (fraction, 5 decimal places)
    "slen",      # subject sequence length
    "scov",      # subject coverage (fraction, 5 decimal places)
    "staxid",    # subject taxonomy ID (NCBI or SILVA)
    "sname",     # subject name / title
    "division",  # NCBI division or SILVA domain
]

# Fields added later by reporting steps (not in raw annotation table)
ANNOT_LINEAGE_FIELDS: List[str] = [
    "organism", "species", "genus", "family", "order", "class", "kingdom",
]


class AnnotationRow(BaseModel):
    """Pydantic model for a single annotation table row.

    Used for row-level validation when building the table.
    For bulk operations use ``empty_annotation_df()`` and pandas.
    """

    search: str
    db: str
    qseqid: str
    orf: str = "NA"
    qseqlen: Optional[int] = None
    sseqid: str
    bitscore: float
    evalue: float = Field(default=0.0, ge=0.0)
    alen: int = Field(default=0, ge=0)
    pident: float = Field(default=0.0, ge=0.0, le=100.0)
    qlen: int = Field(default=0, ge=0)
    qcov: float = Field(default=0.0, ge=0.0, le=1.0)
    slen: int = Field(default=0, ge=0)
    scov: float = Field(default=0.0, ge=0.0, le=1.0)
    staxid: str = "NA"
    sname: str = ""
    division: str = ""


# ---------------------------------------------------------------------------
# DataFrame helpers
# ---------------------------------------------------------------------------

ANNOT_DTYPES = {
    "search": str,
    "db": str,
    "qseqid": str,
    "orf": str,
    "qseqlen": "Int64",   # nullable int
    "sseqid": str,
    "bitscore": float,
    "evalue": float,
    "alen": "Int64",
    "pident": float,
    "qlen": "Int64",
    "qcov": float,
    "slen": "Int64",
    "scov": float,
    "staxid": str,
    "sname": str,
    "division": str,
}


def empty_annotation_df() -> pd.DataFrame:
    """Return an empty DataFrame with the annotation table schema."""
    return pd.DataFrame(columns=ANNOT_FIELDS).astype(
        {k: v for k, v in ANNOT_DTYPES.items() if k in ANNOT_FIELDS}
    )


def read_annotation_tsv(path) -> pd.DataFrame:
    """Read a TSV annotation table from disk."""
    df = pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )
    # Ensure all expected columns are present
    for col in ANNOT_FIELDS:
        if col not in df.columns:
            df[col] = ""
    # Cast numeric columns
    for col in ("bitscore", "evalue", "pident", "qcov", "scov"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    for col in ("qseqlen", "alen", "qlen", "slen"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce").astype("Int64")
    return df


def write_annotation_tsv(df: pd.DataFrame, path) -> None:
    """Write annotation DataFrame to TSV."""
    df.to_csv(path, sep="\t", index=False)


# String columns that may contain commas or special characters —
# quoted on CSV write.  Includes raw annotation fields plus the
# lineage columns added by the reporting step.
_QUOTED_COLS: List[str] = [
    "qseqid", "orf", "sseqid", "sname", "bphage", "division",
] + ANNOT_LINEAGE_FIELDS


def read_annotation_csv(path) -> pd.DataFrame:
    """Read a CSV annotation table written by ``write_annotation_csv()``.

    Handles double-quoted string fields and applies the same type
    casts as ``read_annotation_tsv()``.
    """
    df = pd.read_csv(
        path,
        sep=",",
        dtype=str,
        keep_default_na=False,
        quotechar='"',
    )
    for col in ANNOT_FIELDS:
        if col not in df.columns:
            df[col] = ""
    for col in ("bitscore", "evalue", "pident", "qcov", "scov"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    for col in ("qseqlen", "alen", "qlen", "slen"):
        if col in df.columns:
            df[col] = pd.to_numeric(
                df[col], errors="coerce",
            ).astype("Int64")
    return df


def write_annotation_csv(df: pd.DataFrame, path) -> None:
    """Write annotation DataFrame to CSV.

    String columns (identifiers, names, lineage fields) are always
    double-quoted to handle commas and special characters in subject
    names and taxon names.  Numeric columns are written unquoted.
    """
    df = df.copy()
    for col in _QUOTED_COLS:
        if col in df.columns:
            df[col] = df[col].astype(str)

    df.to_csv(
        path,
        sep=",",
        index=False,
        quoting=csv.QUOTE_NONNUMERIC,
        quotechar='"',
    )


def validate_annotation_df(df: pd.DataFrame) -> pd.DataFrame:
    """Ensure required columns present and types are correct; return df."""
    missing = [c for c in ANNOT_FIELDS if c not in df.columns]
    if missing:
        raise ValueError(f"Annotation table missing columns: {missing}")
    return df


def blast_to_annotation(
    hits: pd.DataFrame,
    search: str,
    db: str,
    dbtype: str = "nucl",
) -> pd.DataFrame:
    """Convert a BLAST/Diamond hit DataFrame to annotation table format.

    Maps the 14/15-field BLAST/Diamond output to the 18-field annotation
    table.  Computes qcov and scov from alignment coordinates.

    Parameters
    ----------
    hits : pd.DataFrame
        Output from ``parse_blast_output()`` or ``parse_diamond_output()``.
    search : str
        Search program name (e.g. "blastn", "diamond.blastp").
    db : str
        Database basename.
    dbtype : str
        "nucl" or "prot".
    """
    if hits.empty:
        return empty_annotation_df()

    out = pd.DataFrame()

    # Query / subject IDs
    out["qseqid"] = hits["qseqid"]
    out["search"] = search
    out["db"] = db

    # Extract ORF tag from qseqid if present (ORF=start-end:strand)
    out["orf"] = hits["qseqid"].str.extract(r"ORF=([^\s]+)", expand=False).fillna("NA")

    out["qseqlen"] = pd.to_numeric(hits.get("qlen", pd.Series(dtype=object)), errors="coerce").astype("Int64")
    out["sseqid"] = hits.get("sseqid", hits.get("saccver", ""))
    out["bitscore"] = pd.to_numeric(hits["bitscore"], errors="coerce")
    out["evalue"] = pd.to_numeric(hits.get("evalue", pd.Series(0, index=hits.index)), errors="coerce").fillna(0)
    out["alen"] = pd.to_numeric(hits["length"], errors="coerce").astype("Int64")
    out["pident"] = pd.to_numeric(hits["pident"], errors="coerce")

    # Query alignment region length
    qstart = pd.to_numeric(hits.get("qstart", pd.Series(0, index=hits.index)), errors="coerce").fillna(0)
    qend = pd.to_numeric(hits.get("qend", pd.Series(0, index=hits.index)), errors="coerce").fillna(0)
    out["qlen"] = (qend - qstart).abs().astype("Int64")

    # Query coverage
    qlen_num = pd.to_numeric(hits.get("qlen", pd.Series(dtype=object)), errors="coerce")
    out["qcov"] = (out["qlen"] / qlen_num.replace(0, pd.NA)).round(5)

    # Subject
    slen_num = pd.to_numeric(hits.get("slen", pd.Series(dtype=object)), errors="coerce")
    out["slen"] = slen_num.astype("Int64")
    sstart = pd.to_numeric(hits.get("sstart", pd.Series(0, index=hits.index)), errors="coerce").fillna(0)
    send = pd.to_numeric(hits.get("send", pd.Series(0, index=hits.index)), errors="coerce").fillna(0)
    s_alen = (send - sstart).abs()
    out["scov"] = (s_alen / slen_num.replace(0, pd.NA)).round(5)

    staxid_col = hits.get("staxid", hits.get("staxids", pd.Series("NA", index=hits.index)))
    out["staxid"] = staxid_col.astype(str)
    out["sname"] = hits.get("stitle", "")
    out["division"] = ""

    return out[ANNOT_FIELDS].reset_index(drop=True)


def paf_to_annotation(
    hits: pd.DataFrame,
    search: str = "minimap2",
    db: str = "",
    c: bool = False,
) -> pd.DataFrame:
    """Convert a parsed PAF DataFrame to annotation table format.

    Parameters
    ----------
    hits : pd.DataFrame
        Output from ``parse_paf_text()`` / ``parse_paf()``.
    search : str
        Search program name.
    db : str
        Database basename.
    c : bool
        Whether minimap2 was run with ``-c`` (full DP alignment).

        ``c=True`` — DP mode: bitscore from ``AS:i:`` tag; accurate
        alignment coordinates available; qcov and scov are computed from
        ``qstart``/``qend`` and ``sstart``/``send``.

        ``c=False`` (default, fast) — chain-only mode: bitscore from
        ``s1:i:`` tag; alignment coordinates are approximate; qcov and
        scov are set to ``pd.NA`` because the reported span may not reflect
        the true aligned region.
    """
    if hits.empty:
        return empty_annotation_df()

    out = pd.DataFrame()
    out["qseqid"] = hits["qseqid"]
    out["search"] = search
    out["db"] = db
    out["orf"] = hits["qseqid"].str.extract(r"ORF=([^\s]+)", expand=False).fillna("NA")
    out["qseqlen"] = hits["qlen"].astype("Int64")
    out["sseqid"] = hits["sseqid"]
    out["bitscore"] = pd.to_numeric(hits["bitscore"], errors="coerce")
    out["evalue"] = 0.0   # minimap2 doesn't report E-values
    out["alen"] = hits["alen"].astype("Int64")

    # pident from matches/alen — valid in both modes
    out["pident"] = (hits["matches"] / hits["alen"].replace(0, pd.NA) * 100).round(2)

    if c:
        # Full DP alignment: coordinates reflect the true aligned region
        out["qlen"] = (hits["qend"] - hits["qstart"]).abs().astype("Int64")
        out["qcov"] = (out["qlen"] / hits["qlen"].replace(0, pd.NA)).round(5)
        out["slen"] = hits["slen"].astype("Int64")
        s_alen = (hits["send"] - hits["sstart"]).abs()
        out["scov"] = (s_alen / hits["slen"].replace(0, pd.NA)).round(5)
    else:
        # Chain-only: coordinates are approximate; mark coverage as NA
        out["qlen"] = pd.array([pd.NA] * len(hits), dtype="Int64")
        out["qcov"] = pd.NA
        out["slen"] = hits["slen"].astype("Int64")
        out["scov"] = pd.NA

    out["staxid"] = "NA"   # filled later from FASTA header or taxonomy lookup
    out["sname"] = hits["sseqid"]
    out["division"] = ""

    return out[ANNOT_FIELDS].reset_index(drop=True)
