"""Taxonomy models for EMtrace.

Two taxonomy sources are supported:

1. **NCBITaxonomy** — adapted from virefset/models/taxonomy.py.
   Parses NCBI taxonomy dump files (nodes.dmp, names.dmp, division.dmp)
   into a flat SQLite3 database and provides division lookups.
   Key difference from virefset: virus_only=False by default, and an
   extra `get_division()` method for the add_division_field() logic.

2. **SILVATaxonomy**
   Parses the SILVA taxmap TSV (taxmap_slv_ssu_ref_138.2.txt) into an
   in-memory dict and optional SQLite3 database for fast accession →
   taxpath lookups used in SILVA-mode annotation reporting.

Examples
--------
NCBI taxonomy::

    from pathlib import Path
    import pandas as pd
    from emtrace.models.taxonomy import NCBITaxonomy

    tax = NCBITaxonomy(Path("/data/taxonomy/taxonomy.db"))

    # Build once from NCBI taxdump files
    tax.build_database(
        nodes_dmp=Path("/data/taxonomy/nodes.dmp"),
        names_dmp=Path("/data/taxonomy/names.dmp"),
        division_dmp=Path("/data/taxonomy/division.dmp"),
        merged_dmp=Path("/data/taxonomy/merged.dmp"),
        version="v20260101",
    )

    lineage = tax.get_taxonomy(11234)  # {"family": "Paramyxoviridae", ...}
    division = tax.get_division(11234) # "Viruses"

    # Add 'division' column to annotation DataFrame in-place
    tax.add_division_to_df(df, taxid_col="staxid")

SILVA taxonomy::

    from emtrace.models.taxonomy import SILVATaxonomy

    silva = SILVATaxonomy(Path("/data/taxonomy/SILVA/silva.db"))

    # Build once from SILVA taxmap file
    silva.build_database(
        taxmap_path=Path(
            "/data/taxonomy/SILVA/taxmap_slv_ssu_ref_138.2.txt"
        ),
        version="138.2",
    )

    rec = silva.get_by_accession("AB000001.1.1450")
    # {"domain": "Bacteria", "family": "Rhodobacteraceae", ...}

    # Add staxid + lineage columns (domain..genus) to df in-place
    # sseqid column must contain SILVA accessions (acc.start.stop)
    silva.add_lineage_to_df(df, sseqid_col="sseqid", staxid_col="staxid")
"""

from __future__ import annotations

import datetime
import logging
import re
import sqlite3
from pathlib import Path
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)


# =============================================================================
# NCBI Taxonomy
# =============================================================================

class NCBITaxonomy:
    """Parses NCBI taxonomy dumps into a flattened SQLite3 database.

    Adapted from virefset/models/taxonomy.py with two key differences:
    - ``virus_only`` defaults to **False** (emtrace classifies all domains)
    - Adds ``get_division()`` for direct taxid→NCBI division name lookup,
      which powers add_division_field() in the annotation step.

    The database schema stores one row per taxid with ranks as columns,
    plus a ``division`` column with the NCBI division name.
    """

    RANKS_OF_INTEREST: Tuple[str, ...] = (
        "realm", "kingdom", "phylum", "class", "order", "family", "genus", "species"
    )

    def __init__(self, db_path: Path) -> None:
        self.db_path = Path(db_path)

    # ------------------------------------------------------------------ #
    # Private parsing helpers
    # ------------------------------------------------------------------ #

    def _parse_nodes(self, path: Path) -> Dict[int, Tuple[int, str, int]]:
        """Parse nodes.dmp → {taxid: (parent_id, rank, division_id)}."""
        logger.info("Parsing %s", path)
        nodes: Dict[int, Tuple[int, str, int]] = {}
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                parts = [p.strip() for p in line.split("|")]
                taxid = int(parts[0])
                parent_id = int(parts[1])
                rank = parts[2]
                div_id = int(parts[4])
                nodes[taxid] = (parent_id, rank, div_id)
        return nodes

    def _parse_names(self, path: Path) -> Dict[int, str]:
        """Parse names.dmp → {taxid: scientific_name}."""
        logger.info("Parsing %s", path)
        names: Dict[int, str] = {}
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                parts = [p.strip() for p in line.split("|")]
                taxid = int(parts[0])
                name = parts[1]
                name_class = parts[3]
                if name_class == "scientific name":
                    names[taxid] = name
        return names

    def _parse_divisions(self, path: Path) -> Dict[int, str]:
        """Parse division.dmp → {division_id: division_name}."""
        divs: Dict[int, str] = {}
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                parts = [p.strip() for p in line.split("|")]
                if len(parts) >= 3:
                    divs[int(parts[0])] = parts[2]
        return divs

    def _parse_merged(self, path: Path) -> Dict[int, int]:
        """Parse merged.dmp → {old_taxid: new_taxid}."""
        merged: Dict[int, int] = {}
        if not path.exists():
            return merged
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                line = line.rstrip("\t|\n")
                parts = [p.strip() for p in line.split("\t|\t")]
                if len(parts) >= 2:
                    try:
                        merged[int(parts[0])] = int(parts[1])
                    except ValueError:
                        pass
        return merged

    # ------------------------------------------------------------------ #
    # Database build
    # ------------------------------------------------------------------ #

    def build_database(
        self,
        nodes_dmp: Path,
        names_dmp: Path,
        division_dmp: Path,
        merged_dmp: Optional[Path] = None,
        virus_only: bool = False,
        version: Optional[str] = None,
    ) -> None:
        """Build the taxonomy SQLite3 database from NCBI dump files.

        Parameters
        ----------
        nodes_dmp, names_dmp, division_dmp : Path
            NCBI taxdump files.
        merged_dmp : Path, optional
            merged.dmp for old→new taxid resolution.
        virus_only : bool
            Restrict to the NCBI Viruses (VRL) division only.
        version : str, optional
            Version string written to db_info.
        """
        nodes = self._parse_nodes(nodes_dmp)
        names = self._parse_names(names_dmp)
        divisions = self._parse_divisions(division_dmp)
        merged = self._parse_merged(merged_dmp) if merged_dmp else {}

        # Build {taxid: division_id} after merging
        taxid_div: Dict[int, int] = {
            taxid: node[2] for taxid, node in nodes.items()
        }

        virus_div_id: Optional[int] = None
        if virus_only:
            for div_id, name in divisions.items():
                if "Viruses" in name or name == "VRL":
                    virus_div_id = div_id
                    break
            logger.info("Virus division ID: %s", virus_div_id)

        self.db_path.parent.mkdir(parents=True, exist_ok=True)

        with sqlite3.connect(self.db_path) as conn:
            cur = conn.cursor()

            # Schema
            rank_cols = ", ".join(f'"{r}" TEXT' for r in self.RANKS_OF_INTEREST)
            cur.execute("DROP TABLE IF EXISTS taxonomy")
            cur.execute(
                f"CREATE TABLE taxonomy "
                f"(taxid INTEGER PRIMARY KEY, organism TEXT, {rank_cols}, division TEXT)"
            )
            cur.execute("DROP TABLE IF EXISTS merged")
            cur.execute("CREATE TABLE merged (old_taxid INTEGER, new_taxid INTEGER)")
            if merged:
                cur.executemany(
                    "INSERT INTO merged VALUES (?, ?)", list(merged.items())
                )

            # Build rows
            batch: List[List] = []
            taxid_count = 0
            for taxid, (_, _, div_id) in nodes.items():
                if virus_only and div_id != virus_div_id:
                    continue

                lineage: Dict[str, Optional[str]] = {r: None for r in self.RANKS_OF_INTEREST}
                curr = taxid
                while curr != 1 and curr in nodes:
                    p_id, rank, _ = nodes[curr]
                    if rank in self.RANKS_OF_INTEREST:
                        lineage[rank] = names.get(curr)
                    curr = p_id

                div_name = divisions.get(div_id, "")
                organism = names.get(taxid, "")
                row = [taxid, organism] + [lineage[r] for r in self.RANKS_OF_INTEREST] + [div_name]
                batch.append(row)
                taxid_count += 1

                if len(batch) >= 50_000:
                    self._insert_batch(cur, batch)
                    batch = []

            if batch:
                self._insert_batch(cur, batch)

            cur.execute("CREATE INDEX IF NOT EXISTS idx_taxid ON taxonomy(taxid)")

            # db_info
            if version is not None:
                self._write_db_info(cur, version, virus_only, taxid_count)

            conn.commit()

        logger.info("Taxonomy DB written to %s (%d taxids)", self.db_path, taxid_count)

    def _insert_batch(self, cur: sqlite3.Cursor, batch: List[List]) -> None:
        n = 2 + len(self.RANKS_OF_INTEREST) + 1  # taxid + organism + ranks + division
        placeholders = ", ".join(["?"] * n)
        cur.executemany(
            f"INSERT OR REPLACE INTO taxonomy VALUES ({placeholders})", batch
        )

    def _write_db_info(
        self,
        cur: sqlite3.Cursor,
        version: str,
        virus_only: bool,
        taxid_count: int,
    ) -> None:
        try:
            ver_int = int(version.lstrip("v"))
        except ValueError:
            ver_int = 0
        cur.execute(f"PRAGMA user_version = {ver_int}")
        scope = "Viruses only" if virus_only else "All taxa"
        desc = (
            f"NCBI taxonomy SQLite database\n"
            f"Input: NCBI taxonomy dump {version}\n"
            f"Stats: {taxid_count:,} taxids, {scope}"
        )
        cur.execute("DROP TABLE IF EXISTS db_info")
        cur.execute("CREATE TABLE db_info (version TEXT, description TEXT, created_at TEXT)")
        cur.execute(
            "INSERT INTO db_info VALUES (?, ?, ?)",
            (version, desc, datetime.datetime.now().isoformat(timespec="seconds")),
        )

    # ------------------------------------------------------------------ #
    # Query methods
    # ------------------------------------------------------------------ #

    def get_taxonomy(self, taxid: int) -> Optional[Dict[str, Optional[str]]]:
        """Return lineage dict for *taxid*, or None if not found."""
        if not self.db_path.exists():
            return None
        with sqlite3.connect(self.db_path) as conn:
            conn.row_factory = sqlite3.Row
            cur = conn.cursor()
            cur.execute("SELECT * FROM taxonomy WHERE taxid = ?", (taxid,))
            row = cur.fetchone()
        if row is None:
            return None
        return {rank: row[rank] for rank in self.RANKS_OF_INTEREST}

    def get_division(self, taxid: int) -> str:
        """Return the NCBI division name for *taxid* (empty string if unknown).

        Resolves merged taxids automatically.
        This is the core lookup for add_division_field().
        """
        if not self.db_path.exists():
            return ""
        with sqlite3.connect(self.db_path) as conn:
            cur = conn.cursor()
            # Resolve merged
            cur.execute("SELECT new_taxid FROM merged WHERE old_taxid = ?", (taxid,))
            row = cur.fetchone()
            resolved = row[0] if row else taxid
            cur.execute("SELECT division FROM taxonomy WHERE taxid = ?", (resolved,))
            row = cur.fetchone()
        return row[0] if row else ""

    def get_organism(self, taxid: int) -> str:
        """Return the scientific name (organism) for *taxid* (empty string if unknown).

        For NCBI/GenBank records this is the canonical "organism" field — the
        scientific name of the assigned taxid, which may be more granular than
        species (e.g. subspecies, strain, isolate).  Resolves merged taxids.
        """
        if not self.db_path.exists():
            return ""
        with sqlite3.connect(self.db_path) as conn:
            cur = conn.cursor()
            cur.execute("SELECT new_taxid FROM merged WHERE old_taxid = ?", (taxid,))
            row = cur.fetchone()
            resolved = row[0] if row else taxid
            try:
                cur.execute("SELECT organism FROM taxonomy WHERE taxid = ?", (resolved,))
                row = cur.fetchone()
            except sqlite3.OperationalError:
                return ""  # legacy DB without organism column
        return row[0] if row else ""

    def add_division_to_df(self, df, taxid_col: str = "staxid") -> None:
        """Add/overwrite the 'division' column in a DataFrame in-place.

        Resolves each staxid to its NCBI division name.
        """
        if not self.db_path.exists():
            df["division"] = ""
            return
        with sqlite3.connect(self.db_path) as conn:
            merged_map = dict(
                conn.execute("SELECT old_taxid, new_taxid FROM merged").fetchall()
            )
            div_map = dict(
                conn.execute("SELECT taxid, division FROM taxonomy").fetchall()
            )

        def _lookup(taxid_str: str) -> str:
            try:
                taxid = int(taxid_str)
            except (ValueError, TypeError):
                return ""
            resolved = merged_map.get(taxid, taxid)
            return div_map.get(resolved, "")

        df["division"] = df[taxid_col].map(_lookup)

    def add_species_to_df(self, df, taxid_col: str = "staxid") -> None:
        """Add/overwrite the 'species' column in a DataFrame in-place.

        Resolves each staxid to its species-rank name.  Taxa above species
        rank have no species name and resolve to "" (7.6% of the NCBI DB),
        as do taxids absent from the database.
        """
        if not self.db_path.exists():
            df["species"] = ""
            return
        with sqlite3.connect(self.db_path) as conn:
            merged_map = dict(
                conn.execute("SELECT old_taxid, new_taxid FROM merged").fetchall()
            )
            sp_map = dict(
                conn.execute("SELECT taxid, species FROM taxonomy").fetchall()
            )

        def _lookup(taxid_str: str) -> str:
            try:
                taxid = int(taxid_str)
            except (ValueError, TypeError):
                return ""
            resolved = merged_map.get(taxid, taxid)
            # species is NULL for above-species ranks — normalise to ""
            return sp_map.get(resolved) or ""

        df["species"] = df[taxid_col].map(_lookup)

    def get_db_info(self) -> Optional[Dict[str, str]]:
        """Return db_info provenance dict (version, description, created_at)."""
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


# =============================================================================
# SILVA Taxonomy
# =============================================================================

# SILVA taxmap TSV columns (taxmap_slv_ssu_ref_138.2.txt):
# primaryAccession  start  stop  path  organism_name  taxid
SILVA_TAXMAP_RANKS: Tuple[str, ...] = (
    "domain", "phylum", "class", "order", "family", "genus",
)

# Second words of an organism name that do not name a species.
_SILVA_NON_EPITHETS = frozenset({
    "sp", "sp.", "spp", "spp.", "cf", "cf.", "aff", "aff.", "gen", "gen.",
    "nov", "nov.", "str", "str.", "genomosp", "genomosp.", "bacterium",
    "archaeon", "endosymbiont", "symbiont", "group", "clone", "phage", "virus",
})


def clean_silva_species(organism_name: str, genus: str) -> str:
    """Reduce a SILVA ``organism_name`` to a binomial, or "" if it names no species.

    SILVA's taxid identifies a lineage down to genus; ``organism_name`` is the
    free-text source organism of the ENA record. It is only trusted as a species
    when it is a plain binomial whose genus word belongs to the SILVA genus —
    this rejects host names filed under a bacterial path ("Solea senegalensis"
    under *Escherichia-Shigella*), ``uncultured``/``metagenome`` labels and
    ``Genus sp.`` names. Strain designations are dropped. Composite SILVA genera
    ("Escherichia-Shigella") accept any of their component genus words.
    """
    if not organism_name or not genus:
        return ""
    words = organism_name.replace("[", "").replace("]", "").split()
    candidatus = bool(words) and words[0] == "Candidatus"
    if candidatus:
        words = words[1:]
    if len(words) < 2:
        return ""
    genus_word, epithet = words[0], words[1]
    genus_words = set(re.split(r"[-\s]+", genus.replace("[", "").replace("]", "")))
    if genus_word not in genus_words:
        return ""
    if epithet in _SILVA_NON_EPITHETS or not re.fullmatch(r"[a-z][a-z-]+", epithet):
        return ""
    return ("Candidatus " if candidatus else "") + f"{genus_word} {epithet}"


def assign_silva_species_taxids(
    pairs: "set[Tuple[str, str]]",
) -> Dict[Tuple[str, str], str]:
    """Map each (lineage taxid, species) pair to a species-level key.

    SILVA has no species taxids, so a synthetic key ``<taxid>.<k>`` is issued
    per distinct binomial under a lineage, numbered in sorted order so a rebuild
    of the same release reproduces the same keys. A pair with an empty species
    keeps the plain lineage taxid — reads that cannot be named to species are
    counted at genus, exactly as before.
    """
    out: Dict[Tuple[str, str], str] = {}
    by_taxid: Dict[str, List[str]] = {}
    for taxid, species in pairs:
        if species:
            by_taxid.setdefault(taxid, []).append(species)
        else:
            out[(taxid, "")] = taxid
    for taxid, names in by_taxid.items():
        for k, species in enumerate(sorted(names), start=1):
            out[(taxid, species)] = f"{taxid}.{k}"
    return out


class SILVATaxonomy:
    """Parses SILVA taxonomy files and provides accession→taxpath lookups.

    Port of the Perl ``parse_taxmap_silva()`` function.  Builds an SQLite3
    database for fast lookups — suitable for SILVA-mode annotation reporting.

    Input files (from config.yaml → taxonomy.SILVA):
        taxmap_slv  — taxmap_slv_ssu_ref_138.2.txt
                      Columns: primaryAccession, start, stop, path,
                               organism_name, taxid
        tax_slv     — tax_slv_ssu_138.2.txt
                      Columns: path, taxid, rank
    """

    def __init__(self, db_path: Path) -> None:
        self.db_path = Path(db_path)

    # ------------------------------------------------------------------ #
    # Parsing helpers
    # ------------------------------------------------------------------ #

    @staticmethod
    def _make_accession(fields: List[str]) -> str:
        """Build accession key from first three fields (acc.start.stop)."""
        return f"{fields[0]}.{fields[1]}.{fields[2]}"

    @staticmethod
    def _parse_taxpath(path: str) -> Dict[str, str]:
        """Parse a SILVA taxpath into rank→name by POSITION.

        Fallback used only when no ``tax_slv`` rank file is supplied. It assumes
        the path is exactly ``domain;phylum;class;order;family;genus``, which
        held through SILVA 138.2 but **not** from 144 onward — 144 inserted a
        ``kingdom`` rank, so every SILVA 144 bacterial lineage has seven levels
        and positional parsing shifts each rank down by one, landing the family
        name in the genus column. Prefer `_parse_taxpath_ranked`.
        """
        parts = [p.strip() for p in path.rstrip(";").split(";") if p.strip()]
        result: Dict[str, str] = {}
        for i, rank in enumerate(SILVA_TAXMAP_RANKS):
            if i < len(parts):
                result[rank] = parts[i]
        return result

    @staticmethod
    def _parse_taxpath_ranked(
        path: str, rank_of: Dict[str, str],
    ) -> Dict[str, str]:
        """Parse a SILVA taxpath using SILVA's own per-path rank declarations.

        `rank_of` maps a path prefix (``"Bacteria;Pseudomonadati;"``) to the
        rank SILVA assigns it. Each element is placed in the column SILVA says
        it belongs to, so the result is correct whatever depth a release uses
        and whatever intermediate ranks it introduces. Ranks SILVA declares that
        this schema does not carry (``kingdom``, ``subgenus``) are skipped
        rather than shifting everything after them.
        """
        parts = [p.strip() for p in path.rstrip(";").split(";") if p.strip()]
        result: Dict[str, str] = {}
        prefix = ""
        for part in parts:
            prefix = f"{prefix}{part};"
            rank = rank_of.get(prefix)
            if rank in SILVA_TAXMAP_RANKS:
                result[rank] = part
        return result

    def parse_taxmap(
        self,
        taxmap_path: Path,
        rank_of: Optional[Dict[str, str]] = None,
    ) -> Dict[str, Dict]:
        """Parse taxmap_slv_ssu_ref_*.txt into a dict keyed by accession.

        When `rank_of` is supplied (from ``tax_slv_ssu_*.txt``) each lineage
        element is assigned the rank SILVA declares for it. Without it the
        parser falls back to positional assignment, which is only correct for
        releases up to 138.2 — see `_parse_taxpath`.

        Distinct paths are cached, because a SILVA release has millions of
        accessions but only tens of thousands of distinct lineages.

        Returns
        -------
        dict
            {accession: {
                "path": str, "organism_name": str, "taxid": str,
                "domain": str, "phylum": str, ..., "genus": str
            }}
        """
        logger.info("Parsing SILVA taxmap: %s", taxmap_path)
        if rank_of is None:
            logger.warning(
                "No tax_slv rank file supplied — falling back to positional "
                "rank assignment, which is WRONG for SILVA 144 and later "
                "(their lineages carry an extra 'kingdom' level)."
            )
        result: Dict[str, Dict] = {}
        path_cache: Dict[str, Dict[str, str]] = {}
        with open(taxmap_path, encoding="utf-8") as fh:
            next(fh, None)  # skip header
            for ln, line in enumerate(fh, start=2):
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 6:
                    logger.debug("Skipping short line %d", ln)
                    continue
                acc = self._make_accession(fields)
                path = fields[3]
                ranks = path_cache.get(path)
                if ranks is None:
                    ranks = (
                        self._parse_taxpath_ranked(path, rank_of)
                        if rank_of is not None
                        else self._parse_taxpath(path)
                    )
                    path_cache[path] = ranks
                result[acc] = {
                    "path": path,
                    "organism_name": fields[4],
                    "taxid": fields[5],
                    **ranks,
                }
        logger.info(
            "Parsed %d SILVA accessions over %d distinct lineages",
            len(result), len(path_cache),
        )
        return result

    def parse_tax_slv(self, tax_path: Path) -> Dict[str, Tuple[str, str]]:
        """Parse tax_slv_ssu_*.txt → {path: (taxid, rank)}.

        Used to look up the taxonomic rank of a given SILVA path.
        """
        result: Dict[str, Tuple[str, str]] = {}
        with open(tax_path, encoding="utf-8") as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 3:
                    result[parts[0]] = (parts[1], parts[2])
        return result

    # ------------------------------------------------------------------ #
    # Database build
    # ------------------------------------------------------------------ #

    def build_database(
        self,
        taxmap_path: Path,
        tax_path: Optional[Path] = None,
        version: Optional[str] = None,
    ) -> None:
        """Build a SQLite3 database from SILVA taxonomy files.

        Parameters
        ----------
        taxmap_path : Path
            taxmap_slv_ssu_ref_*.txt
        tax_path : Path, optional
            tax_slv_ssu_*.txt (for rank information)
        version : str, optional
            Version string written to db_info.
        """
        rank_of = None
        if tax_path is not None and Path(tax_path).exists():
            rank_of = {
                p: rank for p, (_, rank) in self.parse_tax_slv(Path(tax_path)).items()
            }
            logger.info("Loaded %d SILVA rank declarations", len(rank_of))
        records = self.parse_taxmap(taxmap_path, rank_of=rank_of)
        self.db_path.parent.mkdir(parents=True, exist_ok=True)

        rank_cols = ", ".join(f'"{r}" TEXT' for r in SILVA_TAXMAP_RANKS)

        for rec in records.values():
            rec["species"] = clean_silva_species(rec["organism_name"], rec.get("genus", ""))
        species_taxid = assign_silva_species_taxids(
            {(rec["taxid"], rec["species"]) for rec in records.values()}
        )
        logger.info(
            "Assigned %d SILVA species keys",
            sum(1 for (_, sp) in species_taxid if sp),
        )

        with sqlite3.connect(self.db_path) as conn:
            cur = conn.cursor()
            cur.execute("DROP TABLE IF EXISTS silva")
            cur.execute(
                f"CREATE TABLE silva "
                f"(accession TEXT PRIMARY KEY, taxid TEXT, "
                f"path TEXT, organism_name TEXT, {rank_cols}, "
                f"species TEXT, species_taxid TEXT)"
            )

            # One row per reportable taxon: every species key plus the plain
            # lineage taxid, so abundance tables written by an older DB still link.
            cur.execute("DROP TABLE IF EXISTS silva_species")
            cur.execute(
                "CREATE TABLE silva_species "
                "(species_taxid TEXT PRIMARY KEY, taxid TEXT, species TEXT)"
            )
            lineage_rows = {(t, t, "") for (t, _) in species_taxid}
            cur.executemany(
                "INSERT OR REPLACE INTO silva_species VALUES (?, ?, ?)",
                sorted(lineage_rows | {(k, t, sp) for (t, sp), k in species_taxid.items()}),
            )

            batch = []
            for acc, rec in records.items():
                row = [acc, rec["taxid"], rec["path"], rec["organism_name"]] + [
                    rec.get(r, "") for r in SILVA_TAXMAP_RANKS
                ] + [rec["species"], species_taxid[(rec["taxid"], rec["species"])]]
                batch.append(row)
                if len(batch) >= 50_000:
                    self._insert_silva_batch(cur, batch)
                    batch = []
            if batch:
                self._insert_silva_batch(cur, batch)

            cur.execute("CREATE INDEX IF NOT EXISTS idx_acc ON silva(accession)")
            cur.execute("CREATE INDEX IF NOT EXISTS idx_taxid ON silva(taxid)")

            if version is not None:
                cur.execute("DROP TABLE IF EXISTS db_info")
                cur.execute(
                    "CREATE TABLE db_info (version TEXT, description TEXT, created_at TEXT)"
                )
                cur.execute(
                    "INSERT INTO db_info VALUES (?, ?, ?)",
                    (
                        version,
                        f"SILVA taxonomy SQLite database (version {version})",
                        datetime.datetime.now().isoformat(timespec="seconds"),
                    ),
                )

            conn.commit()

        logger.info("SILVA DB written to %s", self.db_path)

    def _insert_silva_batch(self, cur: sqlite3.Cursor, batch: List[List]) -> None:
        n = 6 + len(SILVA_TAXMAP_RANKS)
        placeholders = ", ".join(["?"] * n)
        cur.executemany(
            f"INSERT OR REPLACE INTO silva VALUES ({placeholders})", batch
        )

    # ------------------------------------------------------------------ #
    # Query methods
    # ------------------------------------------------------------------ #

    def get_by_accession(self, accession: str) -> Optional[Dict[str, str]]:
        """Look up SILVA taxonomy for a given accession string.

        The accession format used in SILVA FASTA headers is
        ``primaryAccession.start.stop`` (e.g. ``AB000001.1.1450``).
        """
        if not self.db_path.exists():
            return None
        with sqlite3.connect(self.db_path) as conn:
            conn.row_factory = sqlite3.Row
            cur = conn.cursor()
            cur.execute("SELECT * FROM silva WHERE accession = ?", (accession,))
            row = cur.fetchone()
        if row is None:
            return None
        return dict(row)

    def add_lineage_to_df(
        self,
        df,
        sseqid_col: str = "sseqid",
        staxid_col: str = "staxid",
    ) -> None:
        """Add SILVA lineage columns to annotation DataFrame in-place.

        Sets ``staxid`` from the accession lookup and adds rank columns plus
        ``species``. On a species-aware DB ``staxid`` is the species-level key
        (``<taxid>.<k>``), so EM and abundance resolve species; on an older DB
        it is the lineage taxid and ``species`` stays empty.
        """
        if not self.db_path.exists():
            for rank in SILVA_TAXMAP_RANKS:
                df[rank] = ""
            df["species"] = ""
            return

        has_species = self.has_species()
        taxid_expr = "species_taxid" if has_species else "taxid"
        species_expr = "species" if has_species else "''"
        with sqlite3.connect(self.db_path) as conn:
            rows = conn.execute(
                f"SELECT accession, {taxid_expr}, organism_name, {species_expr}, "
                + ", ".join(f'"{r}"' for r in SILVA_TAXMAP_RANKS)
                + " FROM silva"
            ).fetchall()

        cols = ["accession", "taxid", "organism_name", "species"] + list(SILVA_TAXMAP_RANKS)
        silva_df = __import__("pandas").DataFrame(rows, columns=cols)
        silva_map = silva_df.set_index("accession").to_dict("index")

        def _lookup(acc: str) -> Dict[str, str]:
            return silva_map.get(str(acc), {})

        for rank in ["taxid", "organism_name", "species"] + list(SILVA_TAXMAP_RANKS):
            df[rank if rank != "taxid" else staxid_col] = df[sseqid_col].map(
                lambda a, r=rank: _lookup(a).get(r, "")
            )

    def has_species(self) -> bool:
        """True if the DB was built with species columns (``silva_species`` table)."""
        with sqlite3.connect(self.db_path) as conn:
            return conn.execute(
                "SELECT 1 FROM sqlite_master WHERE type='table' AND name='silva_species'"
            ).fetchone() is not None

    def lineage_by_taxid(self) -> Dict[str, Dict[str, str]]:
        """Map each reportable taxid to ``species`` plus its lineage ranks.

        Keys are species keys and plain lineage taxids alike. Lineage ranks are
        constant within a lineage taxid, so collapsing ``silva`` by taxid is
        safe for them — unlike ``organism_name``, which is never returned here.
        """
        ranks = ", ".join(f'l."{r}"' for r in SILVA_TAXMAP_RANKS)
        lineage = (
            "(SELECT taxid, "
            + ", ".join(f'"{r}"' for r in SILVA_TAXMAP_RANKS)
            + " FROM silva GROUP BY taxid) l"
        )
        if self.has_species():
            sql = (
                f"SELECT s.species_taxid, s.species, {ranks} "
                f"FROM silva_species s JOIN {lineage} ON s.taxid = l.taxid"
            )
        else:
            sql = f"SELECT l.taxid, '', {ranks} FROM {lineage}"
        with sqlite3.connect(self.db_path) as conn:
            rows = conn.execute(sql).fetchall()
        cols = ["species"] + list(SILVA_TAXMAP_RANKS)
        return {str(r[0]): dict(zip(cols, r[1:])) for r in rows}
