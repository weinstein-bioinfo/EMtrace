"""Pure-Python ORF detection using orfipy_core.

Port of the virefset OrfipyWrapper for emtrace.
Uses the orfipy_core Python API directly — no subprocess.

ORF FASTA IDs follow the format:
    {read_id}_{start}-{end}:{strand}
where start/end are 0-based BED coordinates (end exclusive).
This matches the ID format produced by ``_orf_fasta_id()`` and is
what Diamond blastp will report as ``qseqid``.

Typical usage::

    from emtrace.wrappers.orfipy import Orfipy
    orfipy = Orfipy(config)
    pep_path, orf_map = orfipy.detect_orfs(query_fa, output_dir)
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, List, Literal, Optional, Set

import pandas as pd
from Bio.Seq import Seq

from emtrace.core.exceptions import DependencyError
from emtrace.io.fasta import iter_fasta

if TYPE_CHECKING:
    from emtrace.core.config import OrfipyConfig

logger = logging.getLogger(__name__)

Strand = Literal["f", "r", "b"]
FilterMode = Literal[0, 1, 2]


# ---------------------------------------------------------------------------
# Codon tables
# ---------------------------------------------------------------------------

CODON_TABLES = {
    1: {   # Standard
        "starts": ["ATG", "CTG", "TTG"],
        "stops":  ["TAA", "TAG", "TGA"],
    },
    2: {   # Vertebrate Mitochondrial
        "starts": ["ATT", "ATC", "ATA", "ATG", "GTG"],
        "stops":  ["TAA", "TAG", "AGA", "AGG"],
    },
    3: {   # Yeast Mitochondrial
        "starts": ["ATA", "ATG", "GTG"],
        "stops":  ["TAA", "TAG"],
    },
    4: {   # Mold, Protozoan, Coelenterate Mitochondrial
        "starts": ["TTA", "TTG", "CTG", "ATT", "ATC", "ATA", "ATG", "GTG"],
        "stops":  ["TAA", "TAG"],
    },
    5: {   # Invertebrate Mitochondrial
        "starts": ["TTG", "ATT", "ATC", "ATA", "ATG", "GTG"],
        "stops":  ["TAA", "TAG"],
    },
    11: {  # Bacterial, Archaeal and Plant Plastid
        "starts": ["TTG", "CTG", "ATT", "ATC", "ATA", "ATG", "GTG"],
        "stops":  ["TAA", "TAG", "TGA"],
    },
}


# ---------------------------------------------------------------------------
# Data containers
# ---------------------------------------------------------------------------

@dataclass
class ORFRecord:
    """Single ORF prediction record (BED6 row)."""

    contig: str
    start: int     # 0-based, BED start
    end: int       # exclusive, BED end
    name: str      # internal ORF name (e.g. readID_ORF1)
    desc: str
    score: int = 0
    strand: str = "+"

    @property
    def length(self) -> int:
        return self.end - self.start

    def to_bed_line(self) -> str:
        return (
            f"{self.contig}\t{self.start}\t{self.end}\t"
            f"{self.name}\t{self.score}\t{self.strand}"
        )

    def orf_fasta_id(self) -> str:
        """FASTA record ID: ``{contig}_{start}-{end}:{strand}``.

        This format is what Diamond blastp reports as qseqid and what
        ``detect_orfs()`` uses as the ``orf_id`` key in the orf_map.
        """
        return f"{self.contig}_{self.start}-{self.end}:{self.strand}"

    def overlaps(self, other: "ORFRecord", same_strand: bool = False) -> bool:
        if self.contig != other.contig:
            return False
        if same_strand and self.strand != other.strand:
            return False
        return self.start < other.end and other.start < self.end

    def contains(self, other: "ORFRecord") -> bool:
        if self.contig != other.contig:
            return False
        return self.start <= other.start and self.end >= other.end


@dataclass
class OrfipyResult:
    """Output paths and statistics from an orfipy run."""

    bed_file: Path
    nucl_fasta: Optional[Path] = None
    prot_fasta: Optional[Path] = None
    num_orfs: int = 0
    num_filtered: int = 0


# ---------------------------------------------------------------------------
# Main class
# ---------------------------------------------------------------------------

class Orfipy:
    """Pure-Python ORF detection via orfipy_core API.

    Parameters
    ----------
    config : OrfipyConfig
        ORF detection settings from the pipeline config.

    Raises
    ------
    DependencyError
        If ``orfipy_core`` is not installed.
    """

    def __init__(self, config: "OrfipyConfig") -> None:
        self.min_len = config.min_len
        # orfipy_core expects an integer; treat inf as a very large value
        self.max_len = (
            10_000_000
            if isinstance(config.max_len, float) and math.isinf(config.max_len)
            else int(config.max_len)
        )
        self.table = config.table
        self.partial_5 = config.partial_5
        self.partial_3 = config.partial_3
        self.strand = config.strand
        self.include_stop = False
        self._check_dependency()

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def detect_orfs(
        self,
        query: Path,
        output_dir: Path,
    ) -> tuple[Optional[Path], pd.DataFrame]:
        """Detect ORFs in *query* FASTA and return peptide FASTA + map.

        Parameters
        ----------
        query : Path
            Nucleotide FASTA input.
        output_dir : Path
            Directory for output files (created if absent).

        Returns
        -------
        pep_path : Path or None
            Path to the protein FASTA, or None if no ORFs were found.
        orf_map : pd.DataFrame
            Columns: ``orf_id``, ``read_id``, ``orf_start``, ``orf_stop``,
            ``strand``.  ``orf_id`` == Diamond blastp ``qseqid``.
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        bed_path = output_dir / "orfs.bed"
        pep_path = output_dir / "orfs.pep.fa"

        logger.info(
            "Detecting ORFs in %s (min_len=%d, strand=%s)",
            query.name, self.min_len, self.strand,
        )

        result = self.find_orfs(
            contigs=query,
            orfs_bed=bed_path,
            orfs_prot=pep_path,
        )

        _EMPTY = pd.DataFrame(
            columns=["orf_id", "read_id", "orf_start", "orf_stop", "strand"]
        )
        if result.num_orfs == 0:
            logger.info("No ORFs detected in %s", query.name)
            return None, _EMPTY

        orf_map = self._parse_bed(bed_path)
        logger.info("Detected %d ORFs in %s", result.num_orfs, query.name)
        return pep_path, orf_map

    def find_orfs(
        self,
        contigs: Path,
        orfs_bed: Path,
        orfs_nucl: Optional[Path] = None,
        orfs_prot: Optional[Path] = None,
        filter_mode: FilterMode = 0,
    ) -> OrfipyResult:
        """Find ORFs in a FASTA file.

        Parameters
        ----------
        contigs : Path
            Nucleotide FASTA input.
        orfs_bed : Path
            Output BED6 file with ORF coordinates.
        orfs_nucl : Path, optional
            Output nucleotide FASTA.
        orfs_prot : Path, optional
            Output protein FASTA.
        filter_mode : int
            0 = no filtering, 1 = remove contained ORFs,
            2 = remove overlapping ORFs.

        Returns
        -------
        OrfipyResult
        """
        contigs = Path(contigs)
        if not contigs.exists():
            raise FileNotFoundError(f"Input FASTA not found: {contigs}")

        Path(orfs_bed).parent.mkdir(parents=True, exist_ok=True)

        all_orfs: List[ORFRecord] = []
        sequences: dict[str, str] = {}

        for header, seq in iter_fasta(contigs):
            seq_id = header.lstrip(">").split()[0]
            sequences[seq_id] = seq
            all_orfs.extend(self.find_orfs_in_sequence(seq_id, seq))

        logger.info("Found %d ORFs before filtering", len(all_orfs))

        num_filtered = 0
        if filter_mode > 0:
            original = len(all_orfs)
            all_orfs = self._filter_overlapping_orfs(all_orfs, mode=filter_mode)
            num_filtered = original - len(all_orfs)
            logger.info("Filtered %d overlapping ORFs", num_filtered)

        with open(orfs_bed, "w", encoding="utf-8", newline="\n") as fh:
            for orf in all_orfs:
                fh.write(orf.to_bed_line() + "\n")

        result = OrfipyResult(
            bed_file=Path(orfs_bed),
            num_orfs=len(all_orfs),
            num_filtered=num_filtered,
        )

        if orfs_nucl:
            Path(orfs_nucl).parent.mkdir(parents=True, exist_ok=True)
            self._write_nucl_fasta(all_orfs, sequences, Path(orfs_nucl))
            result.nucl_fasta = Path(orfs_nucl)

        if orfs_prot:
            Path(orfs_prot).parent.mkdir(parents=True, exist_ok=True)
            self._write_prot_fasta(all_orfs, sequences, Path(orfs_prot))
            result.prot_fasta = Path(orfs_prot)

        return result

    def find_orfs_in_sequence(
        self, seq_id: str, sequence: str
    ) -> List[ORFRecord]:
        """Find ORFs in a single sequence using orfipy_core."""
        import orfipy_core  # pylint: disable=import-outside-toplevel

        starts, stops = self._get_codons()
        try:
            orfs = orfipy_core.orfs(
                seq=sequence.upper(),
                name=seq_id,
                minlen=self.min_len,
                maxlen=self.max_len,
                strand=self.strand,
                starts=starts,
                stops=stops,
                include_stop=self.include_stop,
                partial3=self.partial_3,
                partial5=self.partial_5,
            )
        except (TypeError, ValueError, AttributeError) as exc:
            logger.warning("orfipy_core.orfs failed for %s: %s", seq_id, exc)
            orfs = []

        records = []
        for counter, orf in enumerate(orfs, start=1):
            start, end, strand, desc = orf
            records.append(ORFRecord(
                contig=seq_id,
                start=start,
                end=end,
                name=f"{seq_id}_ORF{counter}",
                desc=desc,
                score=0,
                strand=strand,
            ))
        return records

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _check_dependency(self) -> None:
        try:
            import orfipy_core  # noqa: F401
        except ImportError as exc:
            raise DependencyError(
                "orfipy_core",
                "pip install orfipy",
            ) from exc

    def _get_codons(self) -> tuple[List[str], List[str]]:
        codons = CODON_TABLES.get(self.table, CODON_TABLES[1])
        return codons["starts"], codons["stops"]

    def _write_nucl_fasta(
        self,
        orfs: List[ORFRecord],
        sequences: dict[str, str],
        path: Path,
    ) -> None:
        with open(path, "w", encoding="utf-8", newline="\n") as fh:
            for orf in orfs:
                seq = sequences.get(orf.contig, "")
                if not seq:
                    continue
                orf_seq = seq[orf.start:orf.end]
                if orf.strand == "-":
                    orf_seq = self._reverse_complement(orf_seq)
                fh.write(f">{orf.orf_fasta_id()}\n{orf_seq.upper()}\n")

    def _write_prot_fasta(
        self,
        orfs: List[ORFRecord],
        sequences: dict[str, str],
        path: Path,
    ) -> None:
        with open(path, "w", encoding="utf-8", newline="\n") as fh:
            for orf in orfs:
                seq = sequences.get(orf.contig, "")
                if not seq:
                    continue
                orf_seq = seq[orf.start:orf.end]
                if orf.strand == "-":
                    orf_seq = self._reverse_complement(orf_seq)
                prot = self._translate(orf_seq.upper())
                fh.write(f">{orf.orf_fasta_id()}\n{prot}\n")

    @staticmethod
    def _parse_bed(bed_path: Path) -> pd.DataFrame:
        """Parse orfipy BED output into an orf_map DataFrame.

        Derives ``orf_id`` as ``{seqname}_{start}-{end}:{strand}``,
        matching the ``orf_fasta_id()`` format used in the protein FASTA
        and therefore in Diamond blastp ``qseqid`` output.
        """
        empty = pd.DataFrame(
            columns=["orf_id", "read_id", "orf_start", "orf_stop", "strand"]
        )
        if not bed_path.exists() or bed_path.stat().st_size == 0:
            return empty
        try:
            df = pd.read_csv(
                bed_path,
                sep="\t",
                header=None,
                names=["read_id", "start", "end", "name", "score", "strand"],
                dtype=str,
                comment="#",
            )
        except Exception:
            return empty

        df["orf_start"] = pd.to_numeric(df["start"], errors="coerce")
        df["orf_stop"] = pd.to_numeric(df["end"], errors="coerce")
        # Build orf_id matching orf_fasta_id(): contig_start-end:strand
        df["orf_id"] = (
            df["read_id"] + "_"
            + df["start"] + "-"
            + df["end"] + ":"
            + df["strand"]
        )
        return df[["orf_id", "read_id", "orf_start", "orf_stop", "strand"]].copy()

    @staticmethod
    def _reverse_complement(seq: str) -> str:
        complement = str.maketrans("ACGTacgt", "TGCAtgca")
        return seq.translate(complement)[::-1]

    def _translate(self, seq: str) -> str:
        return str(Seq(seq).translate(table=str(self.table), to_stop=False))

    def _filter_overlapping_orfs(
        self,
        orfs: List[ORFRecord],
        mode: FilterMode,
    ) -> List[ORFRecord]:
        if mode == 0 or not orfs:
            return orfs
        sorted_orfs = sorted(orfs, key=lambda x: x.length, reverse=True)
        kept: List[ORFRecord] = []
        removed: Set[int] = set()
        for i, orf in enumerate(sorted_orfs):
            if i in removed:
                continue
            kept.append(orf)
            for j, other in enumerate(sorted_orfs[i + 1:], start=i + 1):
                if j in removed:
                    continue
                if mode == 1 and orf.contains(other):
                    removed.add(j)
                elif mode == 2 and orf.overlaps(other):
                    removed.add(j)
        return kept


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------

def check_orfipy_available() -> bool:
    """Return True if orfipy_core is installed."""
    try:
        import orfipy_core  # noqa: F401
        return True
    except ImportError:
        return False


def get_orfipy_version() -> Optional[str]:
    """Return orfipy_core version string, or None if not installed."""
    try:
        import orfipy_core
        return getattr(orfipy_core, "__version__", "unknown")
    except ImportError:
        return None
