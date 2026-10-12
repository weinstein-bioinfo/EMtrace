"""BLAST wrapper — blastn, blastp, blastx with emtrace output format.

Output format matches the Perl pipeline's 14-field tabular format:
qseqid, sseqid, bitscore, length, pident, qlen, qstart, qend,
slen, sstart, send, qcovs, staxid, stitle
"""

from __future__ import annotations

import io
from pathlib import Path

import pandas as pd

from emtrace.wrappers.base import ToolWrapper

# Default BLAST output fields (outfmt 6)
BLAST_FIELDS = [
    "qseqid", "sseqid", "bitscore", "evalue", "length", "pident",
    "qlen", "qstart", "qend", "slen", "sstart", "send",
    "qcovs", "staxid", "stitle",
]


def parse_blast_output(text: str, fields: list[str] | None = None) -> pd.DataFrame:
    """Parse BLAST tabular output text into a DataFrame."""
    fields = fields or BLAST_FIELDS
    if not text.strip():
        return pd.DataFrame(columns=fields)
    return pd.read_csv(
        io.StringIO(text),
        sep="\t",
        header=None,
        names=fields,
        dtype={"staxid": str},
        na_values=["N/A"],
    )


def parse_blast_file(path: Path | str, fields: list[str] | None = None) -> pd.DataFrame:
    """Parse a BLAST tabular output file into a DataFrame."""
    fields = fields or BLAST_FIELDS
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame(columns=fields)
    return pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=fields,
        dtype={"staxid": str},
        na_values=["N/A"],
    )


class _BlastBase(ToolWrapper):
    """Shared logic for BLAST family tools (blastn, blastp, blastx)."""

    blast_program: str = ""
    default_fields: list[str] = BLAST_FIELDS

    def _build_command(
        self,
        *,
        query: str | Path,
        db: str | Path,
        evalue: float = 10,
        max_target_seqs: int = 5,
        taxids: str | None = None,
        fields: list[str] | None = None,
        query_gencode: int | None = None,
    ) -> list[str]:
        fields = fields or self.default_fields
        outfmt = "6 " + " ".join(fields)

        cmd = [
            self._tool_path,
            "-num_threads", str(self.num_threads),
            "-evalue", str(evalue),
            "-max_target_seqs", str(max_target_seqs),
            "-outfmt", outfmt,
            "-db", str(db),
            "-query", str(query),
        ]
        if taxids:
            cmd.extend(["-taxids", taxids])
        if query_gencode is not None:
            cmd.extend(["-query_gencode", str(query_gencode)])
        return cmd

    def search(
        self,
        query: str | Path,
        db: str | Path,
        output: str | Path | None = None,
        evalue: float = 10,
        max_target_seqs: int = 5,
        taxids: str | None = None,
        fields: list[str] | None = None,
        query_gencode: int | None = None,
    ) -> pd.DataFrame | Path:
        """Run BLAST search.

        If *output* is given, write to file and return Path.
        Otherwise return parsed DataFrame.
        """
        fields = fields or self.default_fields
        kwargs = dict(
            query=query, db=db, evalue=evalue,
            max_target_seqs=max_target_seqs, taxids=taxids,
            fields=fields, query_gencode=query_gencode,
        )
        stdout = self.run(**kwargs)

        if output is not None:
            out_path = Path(output)
            out_path.write_text(stdout)
            self.logger.info("Wrote %s", out_path)
            return out_path

        return parse_blast_output(stdout, fields=fields)


class Blastn(_BlastBase):
    """BLASTn nucleotide search wrapper."""

    tool_name = "blastn"
    blast_program = "blastn"
    install_hint = "conda install -c bioconda blast"


class Blastp(_BlastBase):
    """BLASTp protein search wrapper."""

    tool_name = "blastp"
    blast_program = "blastp"
    install_hint = "conda install -c bioconda blast"

    # blastp uses saccver instead of sseqid in the Perl pipeline
    default_fields: list[str] = [
        "qseqid", "saccver", "bitscore", "evalue", "length", "pident",
        "qlen", "qstart", "qend", "slen", "sstart", "send",
        "qcovs", "staxid", "stitle",
    ]


class Blastx(_BlastBase):
    """BLASTx translated nucleotide search wrapper."""

    tool_name = "blastx"
    blast_program = "blastx"
    install_hint = "conda install -c bioconda blast"

    default_fields: list[str] = [
        "qseqid", "saccver", "bitscore", "evalue", "length", "pident",
        "qlen", "qstart", "qend", "slen", "sstart", "send",
        "qcovs", "staxid", "stitle",
    ]

    def search(
        self,
        query: str | Path,
        db: str | Path,
        output: str | Path | None = None,
        evalue: float = 10,
        max_target_seqs: int = 5,
        taxids: str | None = None,
        fields: list[str] | None = None,
        query_gencode: int = 1,
    ) -> pd.DataFrame | Path:
        """Run blastx with default genetic code 1."""
        return super().search(
            query=query, db=db, output=output, evalue=evalue,
            max_target_seqs=max_target_seqs, taxids=taxids,
            fields=fields, query_gencode=query_gencode,
        )
