"""Pure-Python FASTQ reading and partitioning utilities.

All functions handle plain and gzipped (.gz) files transparently.
"""

from __future__ import annotations

import gzip
import logging
import shutil
from pathlib import Path
from typing import Generator, Iterable, Set, Tuple

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Low-level iterators
# ---------------------------------------------------------------------------

def _open_fastq(path: Path | str):
    """Return an open file handle for a plain or gzipped FASTQ."""
    path = Path(path)
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def iter_fastq(
    path: Path | str,
) -> Generator[Tuple[str, str, str, str], None, None]:
    """Yield (header, seq, plus, qual) tuples from a FASTQ file.

    Works with both plain and gzipped input.
    """
    with _open_fastq(path) as fh:
        while True:
            h = fh.readline()
            if not h:
                break
            s = fh.readline()
            p = fh.readline()
            q = fh.readline()
            if not (h and s and p and q):
                break
            yield (
                h.rstrip("\n"),
                s.rstrip("\n"),
                p.rstrip("\n"),
                q.rstrip("\n"),
            )


def read_id_from_header(header: str) -> str:
    """Extract the read ID from a FASTQ header line.

    Strips the leading '@' and returns the part before the first space,
    matching the Perl regex ``^([^/\\s]+)\\s?``.
    """
    return header.lstrip("@").split()[0]


# ---------------------------------------------------------------------------
# ID-set extraction from alignment files
# ---------------------------------------------------------------------------

def mapped_ids_from_sam(
    sam_path: Path | str,
    min_mapq: int = 0,
) -> Set[str]:
    """Return read IDs with at least one primary mapped alignment in a SAM.

    Parameters
    ----------
    sam_path : Path
        Path to SAM file (uncompressed).
    min_mapq : int
        Minimum mapping quality; alignments below this are skipped.

    Returns
    -------
    set[str]
        Read IDs (without '@') of mapped reads.
    """
    sam_path = Path(sam_path)
    mapped: Set[str] = set()
    with open(sam_path) as fh:
        for line in fh:
            if line.startswith("@"):
                continue
            fields = line.split("\t", 12)
            if len(fields) < 5:
                continue
            flag = int(fields[1])
            # Skip unmapped (0x4), secondary (0x100), supplementary (0x800)
            if flag & 0x4:
                continue
            if flag & 0x100 or flag & 0x800:
                continue
            if int(fields[4]) < min_mapq:
                continue
            mapped.add(fields[0])
    return mapped


def mapped_ids_from_paf(
    paf_path: Path | str,
    min_score: int = 0,
    min_mapq: int = 0,
) -> Set[str]:
    """Return query read IDs from a PAF file, optionally filtered by AS score and MAPQ.

    Parameters
    ----------
    paf_path : Path
        Path to PAF file.
    min_score : int
        Minimum AS:i: bitscore; alignments below are skipped.
    min_mapq : int
        Minimum mapping quality (PAF column 12); alignments below are skipped.

    Returns
    -------
    set[str]
        Query read IDs of mapped reads.
    """
    paf_path = Path(paf_path)
    mapped: Set[str] = set()
    with open(paf_path) as fh:
        for line in fh:
            if not line.strip():
                continue
            fields = line.split("\t")
            if len(fields) < 12:
                continue
            if min_mapq > 0 and int(fields[11]) < min_mapq:
                continue
            if min_score > 0:
                score = 0
                for tag in fields[12:]:
                    if tag.startswith("AS:i:"):
                        score = int(tag[5:])
                        break
                if score < min_score:
                    continue
            mapped.add(fields[0])
    return mapped


# ---------------------------------------------------------------------------
# Partitioning
# ---------------------------------------------------------------------------

def partition_fastq(
    input_fq: Path | str,
    pass_fq: Path | str,
    flt_fq: Path | str,
    filter_ids: Iterable[str],
    *,
    gzip_output: bool = True,
    compresslevel: int = 1,
) -> Tuple[int, int]:
    """Partition a FASTQ into two files based on a set of read IDs.

    Reads whose IDs are in *filter_ids* go to *flt_fq*; the rest to
    *pass_fq*.  Replaces ``seqkit grep -vf / -f`` calls from the Perl
    pipeline.

    Parameters
    ----------
    input_fq : Path
        Input FASTQ (plain or gzipped).
    pass_fq : Path
        Output for reads NOT in filter_ids (passing / host-free reads).
    flt_fq : Path
        Output for reads IN filter_ids (host-matched / removed reads).
    filter_ids : iterable of str
        Read IDs to route to *flt_fq*.
    gzip_output : bool
        If True, gzip both output files.
    compresslevel : int
        gzip compression level (1 = fastest).

    Returns
    -------
    (n_pass, n_flt) : tuple[int, int]
        Read counts for passing and filtered outputs.
    """
    id_set: Set[str] = set(filter_ids)

    def _open_out(path: Path):
        if gzip_output:
            return gzip.open(path, "wt", compresslevel=compresslevel, newline="\n")
        return open(path, "w", newline="\n")

    pass_fq = Path(pass_fq)
    flt_fq = Path(flt_fq)

    n_pass = n_flt = 0
    with _open_out(pass_fq) as fh_pass, _open_out(flt_fq) as fh_flt:
        for header, seq, plus, qual in iter_fastq(input_fq):
            rid = read_id_from_header(header)
            record = f"{header}\n{seq}\n{plus}\n{qual}\n"
            if rid in id_set:
                fh_flt.write(record)
                n_flt += 1
            else:
                fh_pass.write(record)
                n_pass += 1

    return n_pass, n_flt


def extract_fastq(
    input_fq: Path | str,
    output_fq: Path | str,
    keep_ids: Iterable[str],
    *,
    gzip_output: bool = True,
    compresslevel: int = 1,
) -> int:
    """Write only reads whose IDs are in *keep_ids* to *output_fq*.

    Returns
    -------
    int
        Number of reads written.
    """
    id_set: Set[str] = set(keep_ids)

    def _open_out(path: Path):
        if gzip_output:
            return gzip.open(path, "wt", compresslevel=compresslevel, newline="\n")
        return open(path, "w", newline="\n")

    output_fq = Path(output_fq)
    n = 0
    with _open_out(output_fq) as fh:
        for header, seq, plus, qual in iter_fastq(input_fq):
            if read_id_from_header(header) in id_set:
                fh.write(f"{header}\n{seq}\n{plus}\n{qual}\n")
                n += 1
    return n


def copy_fastq_gz(
    input_fq: Path | str,
    output_fq: Path | str,
    compresslevel: int = 1,
) -> Path:
    """Copy a FASTQ (plain or gz) to an always-gzipped output file.

    Used when no reads matched any host genome (pass-through copy).
    """
    input_fq = Path(input_fq)
    output_fq = Path(output_fq)

    if str(input_fq).endswith(".gz"):
        # Raw byte copy — preserve existing gzip stream
        with open(input_fq, "rb") as src, open(output_fq, "wb") as dst:
            shutil.copyfileobj(src, dst)
    else:
        with open(input_fq, "rb") as src:
            with gzip.open(
                output_fq, "wb", compresslevel=compresslevel
            ) as dst:
                shutil.copyfileobj(src, dst)

    return output_fq


def count_reads(path: Path | str) -> int:
    """Count reads in a FASTQ file (plain or gzipped). 4 lines per record."""
    n = 0
    with _open_fastq(path) as fh:
        for _ in fh:
            n += 1
    return n // 4


def subsample_fastq(
    input_path: Path | str,
    output_path: Path | str,
    size: float,
    seed: int = 42,
) -> int:
    """Write a random subsample of *input_path* to *output_path*.

    Parameters
    ----------
    size : float
        ``0 < size <= 1``: fraction of total reads to keep.
        ``size > 1``: absolute number of reads to keep (rounded to int).
    seed : int
        Random seed for reproducibility.

    Returns the number of reads written.
    """
    import random
    rng = random.Random(seed)

    total = count_reads(input_path)
    if size <= 1.0:
        n_keep = max(1, round(total * size))
    else:
        n_keep = min(int(size), total)

    # Reservoir sampling — O(total) time, O(n_keep) memory
    kept: list = []
    for i, record in enumerate(iter_fastq(input_path)):
        if len(kept) < n_keep:
            kept.append(record)
        else:
            j = rng.randint(0, i)
            if j < n_keep:
                kept[j] = record

    with open(output_path, "w", newline="\n") as fh:
        for header, seq, plus, qual in kept:
            fh.write(f"{header}\n{seq}\n{plus}\n{qual}\n")

    return len(kept)
