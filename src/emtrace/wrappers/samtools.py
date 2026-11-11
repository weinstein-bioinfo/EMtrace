"""Samtools wrapper — BAM sorting, fixmate, markdup, indexing, depth coverage.

Provides both individual subcommand methods and the composite
``sort_fixmate_markdup()`` pipeline used in ref-based assembly.
"""

from __future__ import annotations

from pathlib import Path

from emtrace.wrappers.base import ToolWrapper
from emtrace.utils.shell import run_shell, shell_path


class Samtools(ToolWrapper):
    """Samtools BAM processing wrapper."""

    tool_name = "samtools"
    install_hint = "conda install -c bioconda samtools"

    def _build_command(self, *, subcmd: str, args: list[str]) -> list[str]:
        return [self._tool_path, subcmd] + args

    # ------------------------------------------------------------------ #
    # Individual subcommands
    # ------------------------------------------------------------------ #

    def sort(
        self,
        input_file: str | Path,
        output_file: str | Path,
        by_name: bool = False,
        output_fmt: str = "bam",
    ) -> Path:
        """Sort SAM/BAM file."""
        output_file = Path(output_file)
        args = ["-O", output_fmt, "-@", str(self.num_threads)]
        if by_name:
            args.append("-n")
        args.extend(["-o", str(output_file), str(input_file)])
        self.run(subcmd="sort", args=args)
        self.logger.info("Sorted -> %s", output_file)
        return output_file

    def fixmate(
        self,
        input_file: str | Path,
        output_file: str | Path,
    ) -> Path:
        """Fill in mate coordinates and insert size fields."""
        output_file = Path(output_file)
        self.run(
            subcmd="fixmate",
            args=["-@", str(self.num_threads), "-m",
                  str(input_file), str(output_file)],
        )
        return output_file

    def markdup(
        self,
        input_file: str | Path,
        output_file: str | Path,
        remove: bool = True,
        stats: bool = True,
    ) -> Path:
        """Mark (and optionally remove) duplicate reads."""
        output_file = Path(output_file)
        args = ["-@", str(self.num_threads)]
        if remove:
            args.append("-r")
        if stats:
            args.append("-s")
        args.extend([
            "-T", f"{output_file}.tmp",
            str(input_file),
            str(output_file),
        ])
        self.run(subcmd="markdup", args=args)
        self.logger.info("Marked duplicates -> %s", output_file)
        return output_file

    def index(self, bam_file: str | Path) -> Path:
        """Index a BAM file (creates .bai)."""
        bam_file = Path(bam_file)
        self.run(
            subcmd="index",
            args=["-@", str(self.num_threads), str(bam_file)],
        )
        self.logger.info("Indexed %s", bam_file)
        return bam_file

    def view(
        self,
        input_file: str | Path,
        output_file: str | Path,
        flags_include: int | None = None,
        flags_exclude: int | None = None,
        output_fmt: str = "bam",
    ) -> Path:
        """Convert/filter SAM/BAM with samtools view."""
        output_file = Path(output_file)
        args = [
            "-@", str(self.num_threads),
            "-O", output_fmt,
            "-o", str(output_file),
        ]
        if flags_include is not None:
            args.extend(["-f", str(flags_include)])
        if flags_exclude is not None:
            args.extend(["-F", str(flags_exclude)])
        args.append(str(input_file))
        self.run(subcmd="view", args=args)
        return output_file

    def flagstat(self, bam_file: str | Path) -> str:
        """Run samtools flagstat and return stdout."""
        return self.run(
            subcmd="flagstat",
            args=["-@", str(self.num_threads), str(bam_file)],
        )

    def count_mapped(self, bam_file: str | Path) -> int:
        """Count primary mapped reads (``samtools view -c -F 0x904``).

        Excludes unmapped (0x4), secondary (0x100) and supplementary (0x800)
        records, so the result is the number of distinct reads aligned to the
        reference. Returns 0 for a BAM with no alignments.
        """
        out = self.run(
            subcmd="view",
            args=["-c", "-F", "2308", str(bam_file)],
        )
        return int(out.strip() or 0)

    def coverage(self, bam_file: str | Path) -> str:
        """Run samtools coverage and return TSV stdout.

        Columns: #rname, startpos, endpos, numreads, covbases,
                 coverage, meandepth, meanbaseq, meanmapq
        """
        return self.run(
            subcmd="coverage",
            args=[str(bam_file)],
        )

    # ------------------------------------------------------------------ #
    # Composite pipeline
    # ------------------------------------------------------------------ #

    def depth(
        self,
        bam: str | Path,
        output: str | Path,
        min_base_quality: int | None = None,
        min_map_quality: int | None = None,
        all_positions: bool = True,
    ) -> Path:
        """Run ``samtools depth`` and write per-base depth to *output*.

        Output is a 3-column TSV (no header): chrom, pos (1-based), depth.

        Parameters
        ----------
        bam : Path
            Input BAM (must be sorted and indexed).
        output : Path
            Output depth file.
        min_base_quality : int, optional
            Skip bases below this base quality (``-q``).
        min_map_quality : int, optional
            Skip alignments below this mapping quality (``-Q``).
        all_positions : bool
            Output all positions including zero-depth (``-a``), default True.
        """
        output = Path(output)
        args: list[str] = []
        if all_positions:
            args.append("-a")
        if min_base_quality is not None:
            args.extend(["-q", str(min_base_quality)])
        if min_map_quality is not None:
            args.extend(["-Q", str(min_map_quality)])
        args.append(str(bam))
        stdout = self.run(subcmd="depth", args=args)
        output.write_text(stdout, newline="\n")
        self.logger.info("Depth -> %s", output)
        return output

    def depth_low_coverage_mask(
        self,
        bam: str | Path,
        output_bed: str | Path,
        min_coverage: int = 5,
        min_base_quality: int | None = None,
        min_map_quality: int | None = None,
        count_deletions: bool = True,
    ) -> Path:
        """Identify low-coverage regions and write them as a BED file.

        Replaces the bedtools pipeline::

            bedtools genomecov -ibam BAM -bga | awk '$4 < cov' | bedtools merge

        Uses ``samtools depth -aa`` then merges consecutive low-depth
        positions in pure Python — no bedtools required.  ``-aa`` (not ``-a``)
        is essential: a reference sequence with *zero* mapped reads is treated
        by ``-a`` as an "unused" reference and emits no rows at all, which would
        leave the BED empty, mask nothing, and let ``bcftools consensus`` return
        the full unmasked reference (mis-scored as 100% coverage).  ``-aa`` emits
        all-zero-depth rows for such references so they are fully masked.

        Parameters
        ----------
        bam : Path
            Input BAM (must be sorted and indexed).
        output_bed : Path
            Output BED file with low-coverage intervals (0-based, half-open).
        min_coverage : int
            Positions with depth strictly below this threshold are masked.
        min_base_quality, min_map_quality : int, optional
            Forwarded to ``samtools depth`` (``-q``/``-Q``).
        count_deletions : bool
            Count reads with a deletion at a position as covering it
            (``samtools depth -J``).  Without it a position that the reads
            delete looks uncovered, is masked with N, and the called deletion
            can never reach the consensus.
        """
        bam = Path(bam)
        output_bed = Path(output_bed)

        # Run samtools depth -aa (include zero-depth positions AND references
        # with no mapped reads — see docstring)
        args: list[str] = ["-aa"]
        if count_deletions:
            args.append("-J")
        if min_base_quality is not None:
            args.extend(["-q", str(min_base_quality)])
        if min_map_quality is not None:
            args.extend(["-Q", str(min_map_quality)])
        args.append(str(bam))
        depth_text = self.run(subcmd="depth", args=args)

        # Merge consecutive low-depth positions into BED intervals
        intervals: list[tuple[str, int, int]] = []
        cur_chrom: str | None = None
        cur_start: int | None = None
        cur_end: int | None = None

        for line in depth_text.splitlines():
            if not line:
                continue
            parts = line.split("\t")
            chrom, pos, dep = parts[0], int(parts[1]), int(parts[2])
            bed_start = pos - 1  # convert 1-based pos to 0-based BED start

            if dep < min_coverage:
                if cur_chrom == chrom and cur_end == bed_start:
                    cur_end = bed_start + 1          # extend current interval
                else:
                    if cur_chrom is not None:        # emit previous interval
                        intervals.append((cur_chrom, cur_start, cur_end))  # type: ignore[arg-type]
                    cur_chrom, cur_start, cur_end = chrom, bed_start, bed_start + 1
            else:
                if cur_chrom is not None:            # close open interval
                    intervals.append((cur_chrom, cur_start, cur_end))  # type: ignore[arg-type]
                    cur_chrom = cur_start = cur_end = None  # type: ignore[assignment]

        if cur_chrom is not None:                    # flush final interval
            intervals.append((cur_chrom, cur_start, cur_end))  # type: ignore[arg-type]

        with open(output_bed, "w", newline="\n") as fh:
            for chrom, start, end in intervals:
                fh.write(f"{chrom}\t{start}\t{end}\n")

        self.logger.info(
            "Low-coverage mask (%d intervals) -> %s",
            len(intervals), output_bed,
        )
        return output_bed

    @staticmethod
    def maskfasta(
        fasta: str | Path,
        bed: str | Path,
        output: str | Path,
        mask_char: str = "N",
    ) -> Path:
        """Mask FASTA regions in *bed* with *mask_char* (default 'N').

        Pure-Python replacement for ``bedtools maskfasta``.
        BED intervals are 0-based, half-open (standard BED format).

        Parameters
        ----------
        fasta : Path
            Input FASTA (plain or gzipped).
        bed : Path
            BED file with regions to mask.
        output : Path
            Output masked FASTA (plain text).
        mask_char : str
            Character to substitute; default ``'N'``.
        """
        import gzip  # pylint: disable=import-outside-toplevel

        fasta = Path(fasta)
        bed = Path(bed)
        output = Path(output)

        # Load BED intervals per chromosome
        mask: dict[str, list[tuple[int, int]]] = {}
        if bed.exists() and bed.stat().st_size > 0:
            with open(bed) as fh:
                for line in fh:
                    line = line.rstrip("\n")
                    if not line or line.startswith("#"):
                        continue
                    parts = line.split("\t")
                    chrom, start, end = parts[0], int(parts[1]), int(parts[2])
                    mask.setdefault(chrom, []).append((start, end))

        def _open(p: Path):
            return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p)

        with _open(fasta) as fh_in, open(output, "w", newline="\n") as fh_out:
            header = ""
            seq_parts: list[str] = []

            def _flush(hdr: str, parts: list[str]) -> None:
                if not hdr:
                    return
                seq_id = hdr.lstrip(">").split()[0]
                seq = bytearray("".join(parts).encode())
                for start, end in mask.get(seq_id, []):
                    seq[start:end] = (mask_char * (end - start)).encode()
                fh_out.write(f"{hdr}\n{seq.decode()}\n")

            for line in fh_in:
                line = line.rstrip("\n")
                if line.startswith(">"):
                    _flush(header, seq_parts)
                    header = line
                    seq_parts = []
                else:
                    seq_parts.append(line)
            _flush(header, seq_parts)

        return output

    def sort_fixmate_markdup(
        self,
        input_sam: str | Path,
        output_bam: str | Path,
    ) -> Path:
        """Run the full SAM->BAM pipeline: sort -> fixmate -> markdup.

        This mirrors the Perl pipeline::

            cat SAM | samtools sort -n | samtools fixmate -m \\
                | samtools sort > BAM
            samtools markdup -r -s BAM MD.BAM
            samtools index MD.BAM

        Parameters
        ----------
        input_sam : Path
            Input SAM file from aligner.
        output_bam : Path
            Output deduplicated, indexed BAM.

        Returns
        -------
        Path
            Path to the output BAM file.
        """
        output_bam = Path(output_bam)
        t = str(self.num_threads)

        # Piped sort -> fixmate -> sort -> output
        tool = shell_path(self._tool_path)
        pipe_cmd = (
            f"cat {shell_path(input_sam)} | "
            f"{tool} sort -O bam -@ {t} -n | "
            f"{tool} fixmate -@ {t} -m - - | "
            f"{tool} sort -O bam -@ {t} -o {shell_path(output_bam)}"
        )
        run_shell(pipe_cmd, tool_name="samtools")

        # markdup
        md_bam = output_bam.with_suffix(".md.bam")
        self.markdup(output_bam, md_bam)

        # index
        self.index(md_bam)

        self.logger.info("BAM pipeline complete -> %s", md_bam)
        return md_bam
