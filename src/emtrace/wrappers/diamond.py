"""Diamond wrapper — accelerated blastx/blastp searches.

Output format matches the Perl pipeline's 13-field tabular format:
qseqid, sseqid, bitscore, length, pident, qlen, qstart, qend,
slen, sstart, send, staxids, stitle
"""

from __future__ import annotations

import io
from pathlib import Path

import pandas as pd

from emtrace.wrappers.base import ToolWrapper

# Diamond output fields (note: staxids not staxid, and no qcovs)
DIAMOND_FIELDS = [
    "qseqid", "sseqid", "bitscore", "evalue", "length", "pident",
    "qlen", "qstart", "qend", "slen", "sstart", "send",
    "staxids", "stitle",
]


def parse_diamond_output(
    text: str, fields: list[str] | None = None
) -> pd.DataFrame:
    """Parse Diamond tabular output text into a DataFrame."""
    fields = fields or DIAMOND_FIELDS
    if not text.strip():
        return pd.DataFrame(columns=fields)
    return pd.read_csv(
        io.StringIO(text),
        sep="\t",
        header=None,
        names=fields,
        dtype={"staxids": str},
        na_values=["N/A"],
    )


def parse_diamond_file(
    path: Path | str, fields: list[str] | None = None
) -> pd.DataFrame:
    """Parse a Diamond tabular output file into a DataFrame."""
    fields = fields or DIAMOND_FIELDS
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame(columns=fields)
    return pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=fields,
        dtype={"staxids": str},
        na_values=["N/A"],
    )


class _DiamondBase(ToolWrapper):
    """Shared logic for Diamond blastx/blastp."""

    tool_name = "diamond"
    install_hint = "conda install -c bioconda diamond"

    diamond_mode: str = ""  # "blastx" or "blastp"
    default_min_orf: int = 150
    # blastx-only flags; subclasses can disable
    _use_min_orf: bool = True
    _use_query_gencode: bool = True

    def _build_command(
        self,
        *,
        query: str | Path,
        db: str | Path,
        output: str | Path | None = None,
        evalue: float = 10,
        max_target_seqs: int = 5,
        min_score: int | None = None,
        min_orf: int | None = None,
        query_gencode: int = 1,
        sensitivity: str = "very-sensitive",
        taxonlist: str | None = None,
        tmpdir: str | Path | None = None,
        fields: list[str] | None = None,
    ) -> list[str]:
        fields = fields or DIAMOND_FIELDS

        cmd = [
            self._tool_path,
            self.diamond_mode,
            "--threads", str(self.num_threads),
            "--evalue", str(evalue),
            "--max-target-seqs", str(max_target_seqs),
            f"--{sensitivity}",
            "--outfmt", "6", *fields,
            "--db", str(db),
            "--query", str(query),
        ]
        if self._use_min_orf:
            cmd.extend(["--min-orf", str(min_orf or self.default_min_orf)])
        if self._use_query_gencode:
            cmd.extend(["--query-gencode", str(query_gencode)])
        if min_score is not None:
            cmd.extend(["--min-score", str(min_score)])
        if taxonlist:
            cmd.extend(["--taxonlist", taxonlist])
        if tmpdir:
            cmd.extend(["--tmpdir", str(tmpdir)])
        if output:
            cmd.extend(["--out", str(output)])
        return cmd

    def search(
        self,
        query: str | Path,
        db: str | Path,
        output: str | Path | None = None,
        evalue: float = 10,
        max_target_seqs: int = 5,
        min_score: int | None = None,
        min_orf: int | None = None,
        query_gencode: int = 1,
        sensitivity: str = "very-sensitive",
        taxonlist: str | None = None,
        tmpdir: str | Path | None = None,
        fields: list[str] | None = None,
    ) -> pd.DataFrame | Path:
        """Run Diamond search.

        If *output* is given, Diamond writes directly to that file (--out)
        and we return the Path.  Otherwise, capture stdout and return
        a parsed DataFrame.
        """
        fields = fields or DIAMOND_FIELDS

        stdout = self.run(
            query=query, db=db, output=output, evalue=evalue,
            max_target_seqs=max_target_seqs, min_score=min_score,
            min_orf=min_orf, query_gencode=query_gencode,
            sensitivity=sensitivity, taxonlist=taxonlist,
            tmpdir=tmpdir, fields=fields,
        )

        if output is not None:
            out_path = Path(output)
            self.logger.info("Wrote %s", out_path)
            return out_path

        return parse_diamond_output(stdout, fields=fields)


class DiamondBlastx(_DiamondBase):
    """Diamond blastx (translated nucleotide vs protein DB)."""

    diamond_mode = "blastx"
    default_min_orf = 150


class DiamondBlastp(_DiamondBase):
    """Diamond blastp (protein vs protein DB).

    Query sequences must already be protein (e.g. produced by orfipy).
    ``--min-orf`` and ``--query-gencode`` are blastx-only and are omitted.
    """

    diamond_mode = "blastp"
    _use_min_orf = False
    _use_query_gencode = False
