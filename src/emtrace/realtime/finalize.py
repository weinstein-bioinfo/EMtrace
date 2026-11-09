"""Connect watch-mode batches to the job's standard layout.

Watch mode processes each batch in ``<job>/watch/batches/batch_NNNN``. Report,
EM, assembly, plot and pack read fixed paths under ``<job>`` (``reads/``,
``annotation/``, ``results/``), so at Stop the completed batches are merged
there and those steps work as they do after ``emtrace run``.

A batch is complete once it has a ``batch.json``:

* Reads move: each batch file is appended to its counterpart in ``<job>/reads``
  and then deleted, so reads are not stored twice. A read file still in a batch
  has not been merged yet, which makes repeating a merge safe.
* Annotation tables are rebuilt from every batch, so batch tables stay in place
  and a resumed watch can reload them.
"""

from __future__ import annotations

import json
import logging
import shutil
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import pandas as pd

from emtrace.core.pipelinejob import WATCH_DIRNAME, PipelineJob
from emtrace.models.annotation import read_annotation_tsv, write_annotation_tsv
from emtrace.utils.paths import ensure_dir

#: FastqWatcher file identity: (device, inode, size, mtime).
FileId = Tuple[int, ...]

logger = logging.getLogger("emtrace.realtime.finalize")

BATCH_MANIFEST = "batch.json"
BATCH_INPUT = "batch_input.fq.gz"
#: Merged raw input in ``<job>/reads``; becomes the job's ``input_reads``.
MERGED_INPUT = "reads.input.fq.gz"
ASSEMBLY_CANDIDATES = "assembly_candidates.tsv"

_READ_PATTERNS = ("reads.*.fq", "reads.*.fq.gz", "hflt.*.reads.fq.gz")
_ANNOT_PATTERN = "annot*.tsv"
#: job.files keys for merged read stages, least to most processed.
_READ_STAGE_KEYS = (
    ("reads.porechop.fq.gz", "porechop_reads"),
    ("reads.chopped.fq.gz", "chopped_reads"),
    ("reads.trim.hflt.fq.gz", "host_filtered_reads"),
)


# ------------------------------------------------------------------
# Batch bookkeeping
# ------------------------------------------------------------------

def batch_dirs(watch_dir: Path) -> List[Path]:
    """All batch directories, complete or not, in batch order."""
    return sorted(
        p for p in (Path(watch_dir) / "batches").glob("batch_*") if p.is_dir()
    )


def completed_batches(watch_dir: Path) -> List[Path]:
    """Batch directories of batches that finished successfully."""
    return [b for b in batch_dirs(watch_dir) if (b / BATCH_MANIFEST).is_file()]


def next_batch_number(watch_dir: Path) -> int:
    """Number for the next batch: never reuses a directory, even a failed one."""
    numbers = []
    for batch in batch_dirs(watch_dir):
        suffix = batch.name[len("batch_"):]
        if suffix.isdigit():
            numbers.append(int(suffix))
    return max(numbers, default=0) + 1


def write_batch_manifest(
    batch_dir: Path,
    batch_num: int,
    files: Iterable[Tuple[Path, Optional[FileId]]],
    n_reads: int,
) -> None:
    """Mark a batch complete, recording its source files and read count."""
    payload = {
        "batch": batch_num,
        "n_reads": n_reads,
        "files": [
            {"path": str(path), "id": list(fid) if fid else None}
            for path, fid in files
        ],
    }
    tmp = batch_dir / f"{BATCH_MANIFEST}.tmp"
    tmp.write_text(json.dumps(payload, indent=1), encoding="utf-8")
    tmp.replace(batch_dir / BATCH_MANIFEST)


def read_batch_manifest(batch_dir: Path) -> dict:
    """The ``batch.json`` of a completed batch."""
    return json.loads((batch_dir / BATCH_MANIFEST).read_text(encoding="utf-8"))


def has_unmerged_batches(watch_dir: Path) -> bool:
    """True if a completed batch still holds reads not merged into the job."""
    return any((b / BATCH_INPUT).exists() for b in completed_batches(watch_dir))


def resume_tables(batch_dir: Path) -> List[Path]:
    """The batch tables a later merge or resume needs (they may not all exist)."""
    annotation = batch_dir / "annotation"
    return [*annotation.glob(_ANNOT_PATTERN), annotation / ASSEMBLY_CANDIDATES]


# ------------------------------------------------------------------
# Which jobs can host a watch
# ------------------------------------------------------------------

def is_watch_job(job_dir: Path) -> bool:
    """True for a job watched before, or created for watching (no input reads)."""
    job_dir = Path(job_dir)
    if (job_dir / WATCH_DIRNAME).is_dir():
        return True
    manifest = PipelineJob.read_job_manifest(job_dir)
    return manifest is not None and not manifest.get("input_reads")


def merged_input(job_dir: Path) -> Optional[Path]:
    """The merged input reads of a watch job, once a watch has been stopped."""
    path = Path(job_dir) / "reads" / MERGED_INPUT
    return path if path.is_file() else None


def watch_conflict(job_dir: Path) -> Optional[str]:
    """Why *job_dir* cannot host a watch, or None if it can.

    A watch merges its results into the job's standard folders, so a job with
    its own input reads or results is refused rather than overwritten. Jobs
    that were watched before can always continue.
    """
    job_dir = Path(job_dir)
    if (job_dir / WATCH_DIRNAME).is_dir():
        return None
    manifest = PipelineJob.read_job_manifest(job_dir) or {}
    if manifest.get("input_reads"):
        return f"it has its own input reads ({manifest['input_reads']})"
    if any((job_dir / "annotation").glob("*.tsv")):
        return "it already has annotation results"
    reads_dir = job_dir / "reads"
    if any(p for pattern in _READ_PATTERNS for p in reads_dir.glob(pattern)):
        return "it already has processed reads"
    return None


def dirs_overlap(watch_dir: Path, job_dir: Path) -> bool:
    """True if either directory contains the other.

    The watcher is recursive, so a job folder inside the watched folder would
    feed its own batch files back in as new input.
    """
    a, b = Path(watch_dir).resolve(), Path(job_dir).resolve()
    return a.is_relative_to(b) or b.is_relative_to(a)


# ------------------------------------------------------------------
# Merging
# ------------------------------------------------------------------

def _append_and_remove(src: Path, dest: Path) -> None:
    """Append *src* to *dest* (gzip members concatenate validly), delete *src*."""
    with open(dest, "ab") as out_fh, open(src, "rb") as in_fh:
        shutil.copyfileobj(in_fh, out_fh, length=1024 * 1024)
    src.unlink()


def _merge_reads(job: PipelineJob, batches: List[Path]) -> int:
    """Move batch reads into ``<job>/reads``; returns the number of files moved."""
    ensure_dir(job.reads_dir)
    moved = 0
    for batch in batches:
        batch_input = batch / BATCH_INPUT
        if batch_input.exists():
            _append_and_remove(batch_input, job.reads_dir / MERGED_INPUT)
            moved += 1
        for pattern in _READ_PATTERNS:
            for src in sorted((batch / "reads").glob(pattern)):
                if ".tmp." in src.name or not src.is_file():
                    continue
                _append_and_remove(src, job.reads_dir / src.name)
                moved += 1
    return moved


def _merge_table(batches: List[Path], name: str, dest: Path) -> bool:
    """Rebuild *dest* from every batch's table *name*.

    An empty result removes *dest* rather than leaving a header-only table,
    which report and EM would prefer over a usable ``annot1.tsv``.
    """
    frames: List[pd.DataFrame] = []
    for batch in batches:
        path = batch / "annotation" / name
        if path.is_file() and path.stat().st_size > 0:
            df = read_annotation_tsv(path)
            if not df.empty:
                frames.append(df)
    if not frames:
        dest.unlink(missing_ok=True)
        return False
    write_annotation_tsv(pd.concat(frames, ignore_index=True), dest)
    return True


def publish_annotations(job: PipelineJob, watch_dir: Path) -> None:
    """Rebuild ``<job>/annotation`` from all completed batches.

    Each table name is merged on its own, so ``annot1`` and ``annot2`` stay
    separate as in a single ``run`` on the same reads.
    ``assembly_candidates.tsv`` comes last: ``candidates_are_current`` requires
    it to be at least as new as the reads.
    """
    batches = completed_batches(watch_dir)
    ensure_dir(job.annotation_dir)
    names = sorted({
        p.name for b in batches for p in (b / "annotation").glob(_ANNOT_PATTERN)
    })
    keys: Dict[str, str] = {"annot1.tsv": "annot1", "annot2.tsv": "annot2"}
    for name in names:
        dest = job.annotation_dir / name
        merged = _merge_table(batches, name, dest)
        if name in keys:
            if merged:
                job.set_file(keys[name], dest)
            else:
                job.files.pop(keys[name], None)

    candidates = job.annotation_dir / ASSEMBLY_CANDIDATES
    if _merge_table(batches, ASSEMBLY_CANDIDATES, candidates):
        job.set_file("assembly_candidates", candidates)
    else:
        job.files.pop("assembly_candidates", None)


def merge_batches_into_job(job: PipelineJob, watch_dir: Path) -> int:
    """Merge all completed batches into the job's standard layout.

    Sets ``job.input_reads`` to the merged input and rewrites ``job.yaml``.
    Returns the number of completed batches.
    """
    batches = completed_batches(watch_dir)
    moved = _merge_reads(job, batches)
    publish_annotations(job, watch_dir)

    merged_input = job.reads_dir / MERGED_INPUT
    if merged_input.exists():
        job.input_reads = merged_input
        job.set_file("input_reads", merged_input)
        job.set_file("current_reads", merged_input)
    for fname, key in _READ_STAGE_KEYS:
        path = job.reads_dir / fname
        if path.exists():
            job.set_file(key, path)
            job.set_file("current_reads", path)
    job._write_job_manifest()

    logger.info(
        "Merged %d batch(es) into %s (%d read file(s) moved)",
        len(batches), job.output_dir, moved,
    )
    return len(batches)
