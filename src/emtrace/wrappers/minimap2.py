"""Minimap2 wrapper — the most-used aligner in emtrace.

Supports long-read alignment (map-ont), all-vs-all overlap (ava-ont),
and assembly-level mapping (asm20).  Parses PAF output into DataFrames.
"""

from __future__ import annotations

import os
import secrets
from pathlib import Path
from typing import Optional

import pandas as pd

from emtrace.utils.filelock import file_lock
from emtrace.utils.shell import run as shell_run
from emtrace.wrappers.base import ToolWrapper

# Standard PAF columns (first 12 mandatory fields)
PAF_COLUMNS = [
    "qseqid",   # 0  query name
    "qlen",     # 1  query length
    "qstart",   # 2  query start (0-based)
    "qend",     # 3  query end
    "strand",   # 4  + or -
    "sseqid",   # 5  target name
    "slen",     # 6  target length
    "sstart",   # 7  target start (0-based)
    "send",     # 8  target end
    "matches",  # 9  number of matching bases
    "alen",     # 10 alignment length
    "mapq",     # 11 mapping quality
]


def _extract_paf_score(fields: list[str], c: bool) -> int:
    """Extract bitscore from PAF optional tags.

    When *c* is True (DP alignment mode, ``-c`` flag), the DP alignment
    score ``AS:i:`` is used.  When False (chain-only mode, default), the
    chaining score ``s1:i:`` is used.  Falls back to 0 if the tag is absent.
    """
    tag_prefix = "AS:i:" if c else "s1:i:"
    for tag in fields[12:]:
        if tag.startswith(tag_prefix):
            return int(tag[len(tag_prefix):])
    return 0


def parse_paf(
    paf_path: Path | str,
    min_score: int = 0,
    c: bool = False,
) -> pd.DataFrame:
    """Parse a PAF file into a DataFrame with bitscore from the appropriate tag.

    Parameters
    ----------
    paf_path : Path
        Path to PAF file.
    min_score : int
        Minimum score to keep. Applied to ``s1:i:`` (chain score, default)
        or ``AS:i:`` (DP alignment score, when *c* is True).
    c : bool
        When True, read DP alignment score from ``AS:i:`` tag (minimap2 was
        run with ``-c``).  When False (default), read chaining score from
        ``s1:i:`` tag.

    Returns
    -------
    pd.DataFrame
        Columns: PAF_COLUMNS + ``bitscore``.
    """
    rows = []
    with open(paf_path) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue

            bitscore = _extract_paf_score(fields, c)
            if bitscore < min_score:
                continue

            rows.append({
                "qseqid": fields[0],
                "qlen": int(fields[1]),
                "qstart": int(fields[2]),
                "qend": int(fields[3]),
                "strand": fields[4],
                "sseqid": fields[5],
                "slen": int(fields[6]),
                "sstart": int(fields[7]),
                "send": int(fields[8]),
                "matches": int(fields[9]),
                "alen": int(fields[10]),
                "mapq": int(fields[11]),
                "bitscore": bitscore,
            })

    if not rows:
        return pd.DataFrame(columns=PAF_COLUMNS + ["bitscore"])
    return pd.DataFrame(rows)


def parse_paf_text(text: str, c: bool = False) -> pd.DataFrame:
    """Parse PAF text (stdout) into a DataFrame.

    Parameters
    ----------
    c : bool
        When True, read DP alignment score from ``AS:i:`` (``-c`` mode).
        When False (default), read chaining score from ``s1:i:``.
    """
    rows = []
    for line in text.splitlines():
        fields = line.split("\t")
        if len(fields) < 12:
            continue
        bitscore = _extract_paf_score(fields, c)
        rows.append({
            "qseqid": fields[0],
            "qlen": int(fields[1]),
            "qstart": int(fields[2]),
            "qend": int(fields[3]),
            "strand": fields[4],
            "sseqid": fields[5],
            "slen": int(fields[6]),
            "sstart": int(fields[7]),
            "send": int(fields[8]),
            "matches": int(fields[9]),
            "alen": int(fields[10]),
            "mapq": int(fields[11]),
            "bitscore": bitscore,
        })
    if not rows:
        return pd.DataFrame(columns=PAF_COLUMNS + ["bitscore"])
    return pd.DataFrame(rows)


class Minimap2(ToolWrapper):
    """Minimap2 long-read aligner wrapper."""

    tool_name = "minimap2"
    install_hint = "conda install -c bioconda minimap2"

    def _build_command(
        self,
        *,
        db: str | Path,
        query: str | Path,
        preset: str = "map-ont",
        output: str | Path | None = None,
        s: int = 0,
        m: int = 0,
        c: bool = False,
        f: float = 0.0002,
        secondary: str = "no",
        n: int = 5,
        k: int = 0,
        w: int = 0,
        index_size: str = "6G",
        K: str = "500M",
        sam_output: bool = False,
        query2: str | Path | None = None,
    ) -> list[str]:
        cmd = [
            self._tool_path,
            "-I", index_size,
            "-t", str(self.num_threads),
            "-x", preset,
            "--secondary", secondary,
            "-f", str(f),
        ]
        if K:
            cmd.extend(["-K", str(K)])
        if k > 0:
            cmd.extend(["-k", str(k)])
        if w > 0:
            cmd.extend(["-w", str(w)])
        if m > 0:
            cmd.extend(["-m", str(m)])
        if c:
            cmd.append("-c")
            if s > 0:
                cmd.extend(["-s", str(s)])
        # -N 0 is a real setting (no secondary alignments), not "unset":
        # only a negative n leaves -N to minimap2's default.
        if n >= 0:
            cmd.extend(["-N", str(n)])
        if sam_output:
            cmd.append("-a")
        cmd.extend([str(db), str(query)])
        if query2 is not None:
            cmd.append(str(query2))
        return cmd

    def align(
        self,
        db: str | Path,
        query: str | Path,
        output: str | Path | None = None,
        preset: str = "map-ont",
        s: int = 0,
        m: int = 0,
        c: bool = False,
        f: float = 0.0002,
        secondary: str = "no",
        n: int = 5,
        k: int = 0,
        w: int = 0,
        index_size: str = "6G",
        K: str = "500M",
        sam_output: bool = False,
        query2: str | Path | None = None,
        log_file: Optional[Path] = None,
    ) -> Path | pd.DataFrame:
        """Run minimap2 alignment.

        Parameters
        ----------
        s : int
            DP Aligment Score
        m : int
            Chain score threshold (``-m``).  Applied regardless of *c*.
        c : bool
            Compute full DP alignment and CIGAR (``-c``).  Enables ``AS:i:``
            score tag and accurate qcov/scov; 2–5× slower than chain-only mode.
        f : float
            Ignore top *f* fraction of most frequent minimizers (``-f``).
            Range 0.0002–1.0.  Minimap2 default is 0.0002.
        n : int
            Secondary alignments to retain (``-N``).  Emitted for any
            ``n >= 0``, so ``0`` keeps primary alignments only; a negative
            value omits the flag and leaves minimap2's default.
        K : str
            Minibatch size for mapping (``-K``), e.g. ``"500M"`` (minimap2
            default).  All reads in a minibatch are aligned and their results
            buffered in RAM before ordered output, so peak memory scales with
            *K*.  Lower it (e.g. ``"100M"``) to reduce peak RAM on large-index
            searches at the cost of slightly more I/O.

        If *output* is given, write PAF/SAM to that file and return its Path.
        Otherwise, capture stdout and return a parsed DataFrame.
        """
        kwargs = dict(
            db=db, query=query, preset=preset, output=output,
            s=s, m=m, c=c, f=f,
            secondary=secondary, n=n, k=k, w=w,
            index_size=index_size, K=K, sam_output=sam_output,
            query2=query2,
        )

        if output is not None:
            out_path = Path(output)
            stdout = self.run(**kwargs, log_file=log_file)
            out_path.write_text(stdout, newline="\n")
            self.logger.info("Wrote %s", out_path)
            return out_path

        stdout = self.run(**kwargs, log_file=log_file)
        if sam_output:
            return stdout  # type: ignore[return-value]
        return parse_paf_text(stdout, c=c)

    def align_to_file(
        self,
        db: str | Path,
        query: str | Path,
        output: str | Path,
        **kwargs,
    ) -> Path:
        """Run minimap2 and write output to file."""
        return self.align(db=db, query=query, output=output, **kwargs)  # type: ignore[return-value]

    def build_index(
        self,
        db: str | Path,
        mmi_path: str | Path,
        k: int = 0,
        w: int = 0,
    ) -> None:
        """Build a minimap2 binary index (.mmi) from *db* and save to *mmi_path*.

        The index is built with ``map-ont`` preset to match host-filter usage.
        ``k`` / ``w`` (k-mer size / minimizer window) are baked into the index
        when ``> 0``; minimap2 ignores ``-k`` / ``-w`` given at map time against a
        prebuilt index, so a non-default k/w index must be built here.

        Built atomically: minimap2 writes to a per-process temporary file which
        is renamed onto *mmi_path* only after a successful build.  If the build
        fails or is interrupted, the partial file is removed so a truncated /
        corrupt index never lands at *mmi_path* (where a later run would silently
        reuse it).  The temp name is unique per process (pid + random token) so
        that concurrent builders (e.g. HPC array tasks targeting the same shared
        index) never clobber one another's partial file — the final atomic
        rename then leaves exactly one complete index.
        """
        mmi_path = Path(mmi_path)
        tmp_path = mmi_path.with_name(
            f"{mmi_path.name}.tmp.{os.getpid()}.{secrets.token_hex(4)}"
        )
        self.logger.info("Building minimap2 index: %s -> %s", db, mmi_path)
        cmd = [
            self._tool_path,
            "-t", str(self.num_threads),
            "-x", "map-ont",
        ]
        if k > 0:
            cmd.extend(["-k", str(k)])
        if w > 0:
            cmd.extend(["-w", str(w)])
        cmd.extend([
            "-d", str(tmp_path),
            str(db),
        ])
        try:
            shell_run(cmd, check=True, tool_name="minimap2", wsl=self.via_wsl)
        except BaseException:
            # Remove the partial index on any failure (including interruption)
            tmp_path.unlink(missing_ok=True)
            raise
        os.replace(tmp_path, mmi_path)  # atomic: only a complete index appears
        self.logger.info("Index saved: %s", mmi_path)

    def host_filter(
        self,
        host_db: str | Path,
        query: str | Path,
        output: str | Path,
        m: int = 40,
        s: int = 80,
        c: bool = False,
        k: int = 0,
        w: int = 0,
        mmi_path: str | Path | None = None,
        save_mmi: bool = False,
    ) -> Path:
        """Align reads against host genome for filtering.

        Two operating modes controlled by *c*:

        ``c=False`` (default, fast)
            minimap2 is called **without** ``-c``; DP alignment is skipped.
            Only the chaining score threshold ``-m {m}`` is applied.
            PAF output contains ``s1:i:`` (chaining score) but **not** ``AS:i:``.

        ``c=True`` (slower, ~2-5×)
            minimap2 is called with ``-c``; full DP alignment is computed.
            Both ``-m {m}`` and ``-s {s}`` thresholds are applied.
            PAF output contains both ``s1:i:`` and ``AS:i:``.

        Index handling: *mmi_path* is used when it exists on disk (faster).
        When *save_mmi* is True and the index is absent it is built first
        (carrying *k*/*w* into the index, since minimap2 ignores ``-k``/``-w``
        supplied at map time against a prebuilt index).

        Parameters
        ----------
        m : int
            Chaining score threshold (``-m``); always applied.
        s : int
            DP alignment score threshold (``-s``); applied only when ``c=True``.
        c : bool
            Compute full DP alignment and CIGAR (``-c``).
        k : int
            k-mer size (``-k``); applied when ``> 0``. Baked into the prebuilt
            index when one is built; also passed at map time for the FASTA path.
        w : int
            Minimizer window (``-w``); applied when ``> 0`` (see *k*).
        mmi_path : Path, optional
            Expected path for the prebuilt ``.mmi`` index.
        save_mmi : bool
            When True, build the index at *mmi_path* if it does not exist.

        Returns path to the PAF output file.
        """
        effective_db: str | Path = host_db
        using_prebuilt_index = False
        if mmi_path is not None:
            mmi = Path(mmi_path)
            if mmi.exists():
                effective_db = mmi          # always use existing index
                using_prebuilt_index = True
            elif save_mmi:
                # Serialize concurrent builders (e.g. HPC array tasks sharing one
                # index dir): hold an advisory lock while building, and re-check
                # existence inside it so only the first task builds and the rest
                # reuse the result. Locking is best-effort — if unsupported, fall
                # back to building (the per-process temp name keeps that safe).
                with file_lock(mmi.with_name(mmi.name + ".lock")):
                    if not mmi.exists():
                        self.build_index(host_db, mmi, k=k, w=w)
                effective_db = mmi          # use freshly built (or just-built) index
                using_prebuilt_index = True

        cmd = [
            self._tool_path,
            "-I", "6G",
            "-t", str(self.num_threads),
            "-x", "map-ont",
            "--secondary", "no",
            "-m", str(m),
        ]
        # -k/-w are baked into a prebuilt index; passing them again is ignored by
        # minimap2 (and warns). Only add them when aligning against a raw FASTA.
        if not using_prebuilt_index:
            if k > 0:
                cmd.extend(["-k", str(k)])
            if w > 0:
                cmd.extend(["-w", str(w)])
        if c:
            cmd.append("-c")
            cmd.extend(["-s", str(s)])
        cmd.extend([str(effective_db), str(query)])

        result = shell_run(
            cmd, check=True, capture=True, tool_name="minimap2", wsl=self.via_wsl,
        )
        out_path = Path(output)
        out_path.write_text(result.stdout or "", newline="\n")
        return out_path

    def self_align(
        self,
        query: str | Path,
        output: str | Path,
    ) -> Path:
        """All-vs-all overlap alignment for chimera detection."""
        return self.align(
            db=query,
            query=query,
            output=output,
            preset="ava-ont",
            secondary="no",
        )  # type: ignore[return-value]
