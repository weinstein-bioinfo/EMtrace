"""Pure-Python FASTA/FASTQ I/O utilities.

Replaces seqkit fq2fa, seqkit grep, seqkit fx2tab calls from the Perl
pipeline with native Python equivalents.
"""

from __future__ import annotations

import gzip
import logging
from pathlib import Path
from typing import Dict, Generator, Iterable, Set, Tuple

from emtrace.io.fastq import iter_fastq, read_id_from_header

logger = logging.getLogger(__name__)


# -------------------------------------------------------------------
# FASTA iterators
# -------------------------------------------------------------------

def _open_fasta(path: Path | str):
    """Return an open file handle for a plain or gzipped FASTA."""
    path = Path(path)
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def iter_fasta(
    path: Path | str,
) -> Generator[Tuple[str, str], None, None]:
    """Yield (header, sequence) tuples from a FASTA file.

    Works with both plain and gzipped input.
    Multi-line sequences are concatenated.
    """
    with _open_fasta(path) as fh:
        header = ""
        seq_parts: list[str] = []
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header:
                    yield header, "".join(seq_parts)
                header = line
                seq_parts = []
            else:
                seq_parts.append(line)
        if header:
            yield header, "".join(seq_parts)


def fasta_id_from_header(header: str) -> str:
    """Extract the sequence ID from a FASTA header line."""
    return header.lstrip(">").split()[0]


# -------------------------------------------------------------------
# FASTQ -> FASTA conversion
# -------------------------------------------------------------------

def fastq_to_fasta(
    input_fq: Path | str,
    output_fa: Path | str,
) -> Path:
    """Convert a FASTQ file (plain or gz) to a plain FASTA file.

    Replaces ``seqkit fq2fa -w0``.
    """
    output_fa = Path(output_fa)
    with open(output_fa, "w") as out:
        for header, seq, _, _ in iter_fastq(input_fq):
            # Convert @header to >header
            fa_header = ">" + header.lstrip("@")
            out.write(f"{fa_header}\n{seq}\n")
    return output_fa


# -------------------------------------------------------------------
# Sequence info table (replaces seqkit fx2tab)
# -------------------------------------------------------------------

def reads_info_tsv(
    input_fq: Path | str,
    output_tsv: Path | str,
) -> Path:
    """Create a seqid<tab>length info table from a FASTQ file.

    Replaces ``seqkit fx2tab --length --name --only-id``.
    Output has columns: seqid, length (with header).
    """
    output_tsv = Path(output_tsv)
    with open(output_tsv, "w") as out:
        out.write("seqid\tlength\n")
        for header, seq, _, _ in iter_fastq(input_fq):
            rid = read_id_from_header(header)
            out.write(f"{rid}\t{len(seq)}\n")
    return output_tsv


# -------------------------------------------------------------------
# Sequence extraction by ID
# -------------------------------------------------------------------

def extract_fasta_by_ids(
    input_fa: Path | str,
    output_fa: Path | str,
    keep_ids: Iterable[str],
) -> int:
    """Write sequences whose IDs are in *keep_ids* to *output_fa*.

    Replaces ``seqkit grep -f``.

    Returns
    -------
    int
        Number of sequences written.
    """
    id_set: Set[str] = set(keep_ids)
    output_fa = Path(output_fa)
    n = 0
    with open(output_fa, "w") as out:
        for header, seq in iter_fasta(input_fa):
            if fasta_id_from_header(header) in id_set:
                out.write(f"{header}\n{seq}\n")
                n += 1
    return n


def exclude_fasta_by_ids(
    input_fa: Path | str,
    output_fa: Path | str,
    exclude_ids: Iterable[str],
) -> int:
    """Write sequences whose IDs are NOT in *exclude_ids*.

    Replaces ``seqkit grep -vf``.

    Returns
    -------
    int
        Number of sequences written.
    """
    id_set: Set[str] = set(exclude_ids)
    output_fa = Path(output_fa)
    n = 0
    with open(output_fa, "w") as out:
        for header, seq in iter_fasta(input_fa):
            if fasta_id_from_header(header) not in id_set:
                out.write(f"{header}\n{seq}\n")
                n += 1
    return n


# -------------------------------------------------------------------
# Sequence length map
# -------------------------------------------------------------------

def fasta_lengths(path: Path | str) -> Dict[str, int]:
    """Return {seqid: length} dict from a FASTA file."""
    return {
        fasta_id_from_header(h): len(s)
        for h, s in iter_fasta(path)
    }
