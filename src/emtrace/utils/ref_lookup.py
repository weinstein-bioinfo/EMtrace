"""Local reference-sequence retrieval for the assembly step.

Resolves a set of accessions against the FASTA files that back the
configured annotation databases, so an assembly run can build its
combined reference without contacting NCBI. The annotation table tells
us which database produced each accession, so we only scan the relevant
DB(s) per accession instead of every configured database.

Two backends:
- FASTA-backed DBs (Minimap2): the original ``.fa``/``.fasta`` sits
  beside the ``.mmi`` index. Streamed via :func:`iter_fasta`.
- BLAST/Diamond DBs: queried via ``blastdbcmd -db <path> -entry <acc>``
  (native, or inside WSL on Windows). Skipped if ``blastdbcmd`` is not
  installed (warned once).
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, Optional, Set, Tuple

import pandas as pd

from emtrace.core.config import DatabaseEntry
from emtrace.io.fasta import fasta_id_from_header, iter_fasta
from emtrace.utils import wsl
from emtrace.utils.shell import run as shell_run


def resolve_db_fasta(db_entry: DatabaseEntry) -> Optional[Path]:
    """Return the local FASTA path for *db_entry*, or None if not available.

    For Minimap2 DBs: tries the configured path, then strips a trailing
    ``.mmi`` and looks for ``{base}``, ``{base}.fa``, ``{base}.fasta``.
    For BLAST/Diamond DBs: returns None (sequence extraction goes through
    ``blastdbcmd``).
    """
    if db_entry.search.lower() != "minimap2":
        return None

    raw = Path(db_entry.db)
    candidates = [raw]
    if raw.suffix == ".mmi":
        base = raw.with_suffix("")
        candidates.extend([base, base.with_suffix(".fa"), base.with_suffix(".fasta")])
    else:
        candidates.extend([raw.with_suffix(".fa"), raw.with_suffix(".fasta")])

    for c in candidates:
        if c.exists() and c.suffix.lower() in (".fa", ".fasta"):
            return c
    return None


def _extract_from_fasta(
    fasta_path: Path,
    wanted: Set[str],
    out_handle,
) -> Set[str]:
    """Stream *fasta_path*, write matching records to *out_handle*, return found IDs."""
    found: Set[str] = set()
    remaining = set(wanted)
    for header, seq in iter_fasta(fasta_path):
        if not remaining:
            break
        sid = fasta_id_from_header(header)
        if sid in remaining:
            out_handle.write(f"{header}\n{seq}\n")
            found.add(sid)
            remaining.discard(sid)
    return found


def _extract_from_blastdb(
    db_path: Path,
    accession: str,
    out_handle,
    logger: logging.Logger,
    tool: Tuple[str, bool] = ("blastdbcmd", False),
) -> bool:
    """Extract a single accession from a BLAST DB via ``blastdbcmd``.

    *tool* is ``(path, via_wsl)`` from :func:`emtrace.utils.wsl.locate_tool`.
    Returns True when the sequence was written to *out_handle*.
    """
    tool_path, via_wsl = tool
    try:
        result = shell_run(
            [tool_path, "-db", str(db_path), "-entry", accession, "-outfmt", "%f"],
            check=False,
            capture=True,
            tool_name="blastdbcmd",
            wsl=via_wsl,
        )
    except Exception as exc:  # noqa: BLE001
        logger.debug("blastdbcmd error for %s in %s: %s", accession, db_path, exc)
        return False

    if result.returncode != 0 or not (result.stdout or "").startswith(">"):
        return False

    text = result.stdout
    if not text.endswith("\n"):
        text += "\n"
    out_handle.write(text)
    return True


def lookup_local(
    accessions: Set[str],
    annotation_df: pd.DataFrame,
    databases: Dict[str, DatabaseEntry],
    output_fa: Path,
    logger: logging.Logger,
) -> Tuple[Set[str], Set[str]]:
    """Look up *accessions* in the DB(s) the annotation table says they came from.

    Sequences are appended to *output_fa* (which must already exist; the
    caller is expected to truncate it first). Returns ``(found, missing)``.
    """
    if not accessions:
        return set(), set()

    # Map: db_name -> set of accessions to search there
    by_db: Dict[str, Set[str]] = {}
    if not annotation_df.empty and {"db", "sseqid"}.issubset(annotation_df.columns):
        hits = annotation_df[annotation_df["sseqid"].isin(accessions)]
        for db_name, sub in hits.groupby("db"):
            by_db[str(db_name)] = set(sub["sseqid"].astype(str))

    # Accessions that the annotation table never saw — cannot resolve locally
    in_annot = {a for accs in by_db.values() for a in accs}
    no_origin = set(accessions) - in_annot
    if no_origin:
        logger.info(
            "Local lookup: %d accession(s) have no annotation hit, will fall back: %s",
            len(no_origin), sorted(no_origin),
        )

    found: Set[str] = set()
    blastdbcmd_warned = [False]  # nonlocal flag

    with open(output_fa, "a") as out:
        for db_name, wanted in by_db.items():
            db_entry = databases.get(db_name)
            if db_entry is None:
                logger.warning(
                    "Local lookup: db '%s' from annotation not in current config; skipping",
                    db_name,
                )
                continue

            fasta_path = resolve_db_fasta(db_entry)

            if fasta_path is not None:
                logger.info(
                    "Local lookup: scanning %s (%d accession(s)) -> %s",
                    db_name, len(wanted), fasta_path,
                )
                hit_ids = _extract_from_fasta(fasta_path, wanted, out)
                found.update(hit_ids)
                if hit_ids:
                    logger.info(
                        "Local lookup: extracted %d/%d from %s",
                        len(hit_ids), len(wanted), db_name,
                    )
                continue

            # BLAST/Diamond path
            blastdbcmd = wsl.locate_tool("blastdbcmd")
            if blastdbcmd is None:
                if not blastdbcmd_warned[0]:
                    logger.warning(
                        "blastdbcmd not found; cannot extract from BLAST DB '%s' (and any other BLAST DBs)",
                        db_name,
                    )
                    blastdbcmd_warned[0] = True
                continue

            db_path = Path(db_entry.db)
            extracted_here = 0
            for acc in sorted(wanted):
                if _extract_from_blastdb(db_path, acc, out, logger, blastdbcmd):
                    found.add(acc)
                    extracted_here += 1
            logger.info(
                "Local lookup: extracted %d/%d from %s via blastdbcmd",
                extracted_here, len(wanted), db_name,
            )

    missing = set(accessions) - found
    return found, missing
