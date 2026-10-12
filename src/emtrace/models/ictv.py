"""ICTV VMR genome composition, molecule type, and host model.

Parses the ICTV Virus Metadata Resource (VMR) Excel table to extract
per-species genome [composition] (dsDNA, ssRNA, etc.) and host information.
Stores the results in a small SQLite3 database keyed by ICTV species name,
providing fast lookups during annotation enrichment.

Adapted from virefset/models/ictv_builder.py parse_vmr() — but simplified:
we only need genome_composition/molecule_type and host, not full accession
expansion or GenBank fetching.
"""

from __future__ import annotations

import datetime
import logging
import sqlite3
from pathlib import Path
from typing import Dict, Optional

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# ICTV VMR parser
# ---------------------------------------------------------------------------

# VMR columns we care about (others are ignored)
_VMR_COLS = [
    "Species",
    "Genus",
    "Family",
    "Order",
    "Genome",               # primary molecule type field (ssRNA, dsDNA, …)
    "Genome composition",  # fallback field (legacy VMR releases)
    "Host source",
]


def parse_vmr(vmr_path: Path) -> "pandas.DataFrame":  # type: ignore[name-defined]
    """Parse an ICTV VMR Excel file into a per-species DataFrame.

    Handles multi-sheet workbooks by finding the sheet whose name
    contains 'VMR' (e.g. 'VMR MSL40').

    Returns
    -------
    pd.DataFrame
        Columns: species, genus, family, order, molecule_type, host
        One row per species (duplicates resolved by keeping last).
    """
    import pandas as pd

    vmr_path = Path(vmr_path)
    xl = pd.ExcelFile(vmr_path)

    # Find VMR sheet — prefer names starting with "VMR"
    sheet = None
    for name in xl.sheet_names:
        if "VMR" in str(name).upper():
            sheet = name
            break
    if sheet is None:
        sheet = xl.sheet_names[0]
        logger.warning("No VMR sheet found — using first sheet: %s", sheet)
    else:
        logger.info("Using VMR sheet: %s", sheet)

    df = pd.read_excel(vmr_path, sheet_name=sheet, dtype=str)

    # Keep only the columns we need (gracefully handle missing ones)
    available = [c for c in _VMR_COLS if c in df.columns]
    df = df[available].copy().fillna("")

    # Determine molecule type
    genome_col = "Genome" if "Genome" in df.columns else None
    composition_col = "Genome composition" if "Genome composition" in df.columns else None

    def _get_mol_type(row) -> str:
        # Priority: Genome field first
        if genome_col:
            mt = row.get(genome_col, "").strip()
            if mt:
                return mt
        if composition_col:
            return row.get(composition_col, "").strip()
        return ""

    result = pd.DataFrame({
        "species": df.get("Species", ""),
        "genus": df.get("Genus", ""),
        "family": df.get("Family", ""),
        "order": df.get("Order", ""),
        "molecule_type": df.apply(_get_mol_type, axis=1),
        "host": df.get("Host source", ""),
    })

    # Drop rows with no species
    result = result[result["species"].str.strip() != ""]
    # Deduplicate (last row wins, matching Perl behaviour)
    result = result.drop_duplicates(subset=["species"], keep="last").reset_index(drop=True)

    logger.info("Parsed %d ICTV species records from %s", len(result), vmr_path)
    return result


# ---------------------------------------------------------------------------
# SQLite3 database
# ---------------------------------------------------------------------------

class ICTVCompositionDB:
    """
    Small SQLite3 database storing ICTV species → molecule_type + host.

    Built from the ICTV VMR Excel table via ``build_from_vmr()``.
    Queried during annotation enrichment by species name.
    """

    def __init__(self, db_path: Path) -> None:
        self.db_path = Path(db_path)

    def build_from_vmr(
        self,
        vmr_path: Path,
        version: Optional[str] = None,
    ) -> None:
        """Parse VMR Excel and write to SQLite3.

        Parameters
        ----------
        vmr_path : Path
            Path to ICTV VMR .xlsx file.
        version : str, optional
            Version label written to db_info (e.g. 'MSL40').
        """
        import pandas as pd

        df = parse_vmr(vmr_path)
        self.db_path.parent.mkdir(parents=True, exist_ok=True)

        with sqlite3.connect(self.db_path) as conn:
            cur = conn.cursor()
            cur.execute("DROP TABLE IF EXISTS ictv")
            cur.execute(
                """CREATE TABLE ictv (
                    species TEXT PRIMARY KEY,
                    genus TEXT,
                    family TEXT,
                    "order" TEXT,
                    molecule_type TEXT,
                    host TEXT
                )"""
            )
            rows = [
                (row.species, row.genus, row.family, row.order,
                 row.molecule_type, row.host)
                for row in df.itertuples(index=False)
            ]
            cur.executemany("INSERT OR REPLACE INTO ictv VALUES (?,?,?,?,?,?)", rows)

            if version is not None:
                cur.execute("DROP TABLE IF EXISTS db_info")
                cur.execute(
                    "CREATE TABLE db_info (version TEXT, description TEXT, created_at TEXT)"
                )
                cur.execute(
                    "INSERT INTO db_info VALUES (?, ?, ?)",
                    (
                        version,
                        f"ICTV VMR genome composition/host database (version {version})",
                        datetime.datetime.now().isoformat(timespec="seconds"),
                    ),
                )

            conn.commit()

        logger.info(
            "ICTV composition DB written to %s (%d species)", self.db_path, len(rows)
        )

    def get_by_species(self, species: str) -> Optional[Dict[str, str]]:
        """Return the ICTV record for *species*, or None if not found."""
        if not self.db_path.exists():
            return None
        with sqlite3.connect(self.db_path) as conn:
            conn.row_factory = sqlite3.Row
            cur = conn.cursor()
            cur.execute("SELECT * FROM ictv WHERE species = ?", (species,))
            row = cur.fetchone()
        if row is None:
            return None
        return dict(row)

    def get_molecule_type(self, species: str) -> str:
        """Return molecule_type for *species*, or empty string."""
        rec = self.get_by_species(species)
        return rec["molecule_type"] if rec else ""

    def get_host(self, species: str) -> str:
        """Return host source for *species*, or empty string."""
        rec = self.get_by_species(species)
        return rec["host"] if rec else ""

    def add_molecule_type_to_df(
        self,
        df,
        species_col: str = "species",
        out_col: str = "molecule_type",
    ) -> None:
        """Add molecule_type column to annotation/abundance DataFrame in-place."""
        if not self.db_path.exists():
            df[out_col] = ""
            return
        with sqlite3.connect(self.db_path) as conn:
            mol_map = dict(
                conn.execute("SELECT species, molecule_type FROM ictv").fetchall()
            )
        df[out_col] = df[species_col].map(lambda s: mol_map.get(str(s), ""))

    def get_db_info(self) -> Optional[Dict[str, str]]:
        """Return db_info provenance dict."""
        if not self.db_path.exists():
            return None
        with sqlite3.connect(self.db_path) as conn:
            cur = conn.cursor()
            try:
                cur.execute("SELECT version, description, created_at FROM db_info")
                row = cur.fetchone()
            except sqlite3.OperationalError:
                return None
        if row is None:
            return None
        return {"version": row[0], "description": row[1], "created_at": row[2]}
