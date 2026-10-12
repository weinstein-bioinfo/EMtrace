"""Abundance table schema.

Fields: taxid, readn, bitscore, prob, species, genus, family, order, class, division
"""

from __future__ import annotations

import csv
from typing import List

import pandas as pd
from pydantic import BaseModel, Field

ABUNDANCE_FIELDS: List[str] = [
    "taxid",     # NCBI or SILVA taxonomy ID
    "readn",     # read count mapped to this taxon
    "bitscore",  # sum or mean bitscore across mapped reads
    "prob",      # EM-estimated abundance probability
    "readn_em",  # EM expected read count: prob x reads in the EM input
    "organism",  # scientific name from NCBI names.dmp / SILVA organism_name
    "species",
    "genus",
    "family",
    "order",
    "class",
    "division",
]


class AbundanceRow(BaseModel):
    """Pydantic model for a single abundance table row."""

    taxid: str
    readn: int = Field(default=0, ge=0)
    bitscore: float = Field(default=0.0, ge=0.0)
    prob: float = Field(default=0.0, ge=0.0, le=1.0)
    organism: str = ""
    species: str = ""
    genus: str = ""
    family: str = ""
    order: str = ""
    class_: str = Field(default="", alias="class")

    model_config = {"populate_by_name": True}
    division: str = ""


def empty_abundance_df() -> pd.DataFrame:
    """Return an empty DataFrame with the abundance table schema."""
    return pd.DataFrame(columns=ABUNDANCE_FIELDS)


def read_abundance_tsv(path) -> pd.DataFrame:
    """Read a TSV abundance table from disk."""
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    for col in ABUNDANCE_FIELDS:
        if col not in df.columns:
            df[col] = ""
    for col in ("readn",):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce").astype("Int64")
    for col in ("bitscore", "prob"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def write_abundance_tsv(df: pd.DataFrame, path) -> None:
    """Write abundance DataFrame to TSV."""
    df.to_csv(path, sep="\t", index=False)


# Taxon name columns that may contain commas or other special characters
_TAXON_COLS: List[str] = [
    "species", "genus", "family", "order", "class", "division",
]


def read_abundance_csv(path) -> pd.DataFrame:
    """Read a CSV abundance table from disk.

    Handles quoted taxon name fields written by
    ``write_abundance_csv()``.
    """
    df = pd.read_csv(
        path,
        sep=",",
        dtype=str,
        keep_default_na=False,
        quotechar='"',
    )
    for col in ABUNDANCE_FIELDS:
        if col not in df.columns:
            df[col] = ""
    for col in ("readn",):
        if col in df.columns:
            df[col] = pd.to_numeric(
                df[col], errors="coerce",
            ).astype("Int64")
    for col in ("bitscore", "prob"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def write_abundance_csv(df: pd.DataFrame, path) -> None:
    """Write abundance DataFrame to CSV.

    Taxon name columns (species → division) are always
    double-quoted to handle commas and other special characters
    that appear in taxonomic names.
    """
    df = df.copy()

    # Pre-quote taxon name columns by wrapping values with the csv
    # module so pandas writes them with quoting=QUOTE_NONNUMERIC.
    # We force the columns to object dtype (strings) so pandas
    # treats them as non-numeric and quotes them.
    for col in _TAXON_COLS:
        if col in df.columns:
            df[col] = df[col].astype(str)

    df.to_csv(
        path,
        sep=",",
        index=False,
        quoting=csv.QUOTE_NONNUMERIC,
        quotechar='"',
    )
