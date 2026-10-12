"""BWA wrapper — short/long read aligner for host-read filtering.

Uses ``bwa mem`` with samtools-free pure-Python SAM parsing to extract
mapped read IDs.  Mirrors the Perl ``filter_host_reads()`` implementation
in ``Lazypipe::SeqAn`` (replacing sambamba with pure-Python SAM parsing
and seqkit with :mod:`emtrace.io.fastq`).

Index management
----------------
``bwa index -a bwtsw`` is run automatically when the expected index files
(``<db>.amb`` and ``<db>.sa``) are not found.

SAM flag reference (used for filtering)
----------------------------------------
  Bit 0x4 (4)  — read unmapped
  Bit 0x100 (256) — secondary alignment (skipped)
  Bit 0x800 (2048) — supplementary alignment (skipped)

The ``AS:i:`` optional tag carries the alignment score (bwa mem output).
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional, Set

from emtrace.wrappers.base import ToolWrapper
from emtrace.utils.shell import run as shell_run

logger = logging.getLogger("emtrace.wrappers.bwa")

# SAM flags
_FLAG_UNMAPPED      = 0x4
_FLAG_SECONDARY     = 0x100
_FLAG_SUPPLEMENTARY = 0x800
_SKIP_FLAGS         = _FLAG_UNMAPPED | _FLAG_SECONDARY | _FLAG_SUPPLEMENTARY


class BWA(ToolWrapper):
    """BWA aligner wrapper for host read filtering.

    Parameters
    ----------
    num_threads : int
        Threads passed to ``bwa mem -t``.
    """

    tool_name = "bwa"
    install_hint = "conda install -c bioconda bwa"

    def _build_command(self, **kwargs) -> list[str]:
        # Not used directly — index() and mem() build commands manually.
        raise NotImplementedError("Use index() or mem() directly")

    # ------------------------------------------------------------------
    # Index management
    # ------------------------------------------------------------------

    def is_indexed(self, db: str | Path) -> bool:
        """Return True if BWA index files exist for *db*."""
        db = Path(db)
        return (
            Path(str(db) + ".amb").exists()
            and Path(str(db) + ".sa").exists()
        )

    def index(
        self,
        db: str | Path,
        algorithm: str = "bwtsw",
        log_file: Optional[Path] = None,
    ) -> None:
        """Build BWA index for *db* if not already indexed.

        Parameters
        ----------
        db : Path
            Reference FASTA (plain or gzipped).
        algorithm : str
            Indexing algorithm: ``"bwtsw"`` (large genomes) or ``"is"``.
        log_file : Path, optional
            Append bwa index stderr to this file.
        """
        db = Path(db)
        if self.is_indexed(db):
            self.logger.debug("BWA index already exists for %s", db)
            return

        self.logger.info("Building BWA index for %s (algorithm=%s)", db, algorithm)
        cmd = [self._tool_path, "index", "-a", algorithm, str(db)]

        result = shell_run(
            cmd,
            check=True,
            capture=True,
            tool_name="bwa",
        )
        if log_file is not None:
            combined = (result.stdout or "") + (result.stderr or "")
            if combined:
                with open(log_file, "a") as fh:
                    fh.write(combined)

    # ------------------------------------------------------------------
    # Alignment
    # ------------------------------------------------------------------

    def mem(
        self,
        db: str | Path,
        query: str | Path,
        output_sam: Path,
        score_threshold: int = 0,
        log_file: Optional[Path] = None,
    ) -> Path:
        """Run ``bwa mem`` and write SAM to *output_sam*.

        Parameters
        ----------
        db : Path
            Reference FASTA / BWA index prefix.
        query : Path
            Query FASTQ (plain or gzipped; single-end only for Nanopore).
        output_sam : Path
            Destination for the SAM output.
        score_threshold : int
            Passed as ``-T`` to bwa mem (minimum alignment score to output).
        log_file : Path, optional
            Append bwa mem stderr to this file.

        Returns
        -------
        Path
            *output_sam*
        """
        cmd = [
            self._tool_path, "mem",
            "-t", str(self.num_threads),
            "-T", str(score_threshold),
            str(db),
            str(query),
        ]
        result = shell_run(
            cmd,
            check=True,
            capture=True,
            tool_name="bwa",
        )
        output_sam.write_text(result.stdout or "")
        if log_file is not None and result.stderr:
            with open(log_file, "a") as fh:
                fh.write(result.stderr)
        return output_sam

    # ------------------------------------------------------------------
    # Host filtering
    # ------------------------------------------------------------------

    def host_filter(
        self,
        host_db: str | Path,
        query: str | Path,
        score_threshold: int = 30,
        min_score: int = 100,
        min_mapq: int = 0,
        sam_tmp: Optional[Path] = None,
        log_file: Optional[Path] = None,
    ) -> Set[str]:
        """Align *query* reads against *host_db* and return mapped read IDs.

        Replaces the Perl pipeline::

            bwa mem -T score_threshold db query > hostgen.sam
            sambamba view -S -F "not(unmapped) and mapq>=MAPQ and AS>=score" sam

        with pure-Python SAM parsing (no sambamba dependency).

        Steps
        -----
        1. Ensure BWA index exists (runs ``bwa index`` if needed).
        2. Run ``bwa mem -T score_threshold`` → SAM written to *sam_tmp*.
        3. Parse SAM: keep primary mapped alignments with MAPQ >= min_mapq
           and AS >= min_score.
        4. Delete the SAM file; return the set of mapped read IDs.

        Parameters
        ----------
        host_db : Path
            Reference FASTA or BWA index prefix.
        query : Path
            Input reads (gzipped FASTQ).
        score_threshold : int
            Passed as ``-T`` to bwa mem (min score for bwa to output alignment).
        min_score : int
            Post-alignment filter: min AS tag value to call a read as host.
        min_mapq : int
            Post-alignment filter: min MAPQ to call a read as host.
        sam_tmp : Path, optional
            Temporary SAM path.  Deleted after parsing.
        log_file : Path, optional
            Append aligner output to this file.

        Returns
        -------
        set[str]
            Read IDs of reads that map to the host.
        """
        host_db = Path(host_db)
        query   = Path(query)

        # Ensure index
        self.index(host_db, log_file=log_file)

        # Temp SAM location
        if sam_tmp is None:
            sam_tmp = query.parent / "hostgen_bwa.sam"

        try:
            self.mem(
                db=host_db,
                query=query,
                output_sam=sam_tmp,
                score_threshold=score_threshold,
                log_file=log_file,
            )
            mapped = _parse_sam_mapped_ids(
                sam_tmp,
                min_score=min_score,
                min_mapq=min_mapq,
            )
        finally:
            if sam_tmp.exists():
                sam_tmp.unlink()

        self.logger.info(
            "bwa mem: %d reads mapped to %s", len(mapped), host_db.name
        )
        return mapped


# ---------------------------------------------------------------------------
# Pure-Python SAM parser
# ---------------------------------------------------------------------------

def _parse_sam_mapped_ids(
    sam_path: Path,
    min_score: int = 0,
    min_mapq: int = 0,
) -> Set[str]:
    """Parse a SAM file and return IDs of reads that pass mapping filters.

    Criteria (mirrors sambamba ``not(unmapped) and mapq>=MAPQ and AS>=score``):

    * Alignment flags: not unmapped (0x4), not secondary (0x100),
      not supplementary (0x800).
    * MAPQ >= *min_mapq*.
    * ``AS:i:`` tag >= *min_score* (when tag is present; reads without AS
      pass if other criteria are met).

    Parameters
    ----------
    sam_path : Path
        Path to SAM file (as written by bwa mem).
    min_score, min_mapq : int
        Thresholds (inclusive).

    Returns
    -------
    set[str]
        Qualifying read IDs (column 1 of SAM, without /1 /2 suffixes).
    """
    mapped: Set[str] = set()

    with open(sam_path, errors="replace") as fh:
        for line in fh:
            if line.startswith("@"):
                continue                          # SAM header
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                continue

            flag  = int(fields[1])
            mapq  = int(fields[4])
            qname = fields[0]

            # Skip unmapped / secondary / supplementary
            if flag & _SKIP_FLAGS:
                continue

            # MAPQ filter
            if mapq < min_mapq:
                continue

            # AS tag filter
            if min_score > 0:
                as_score = _extract_as_tag(fields[11:])
                if as_score is not None and as_score < min_score:
                    continue

            # Strip /1 /2 pair suffixes so IDs match FASTQ headers
            read_id = qname.split("/")[0].split()[0]
            mapped.add(read_id)

    return mapped


def _extract_as_tag(optional_fields: list[str]) -> Optional[int]:
    """Return the value of the ``AS:i:`` optional SAM field, or None."""
    for field in optional_fields:
        if field.startswith("AS:i:"):
            try:
                return int(field[5:])
            except ValueError:
                return None
    return None
