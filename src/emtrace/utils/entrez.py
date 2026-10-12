"""NCBI Entrez utilities — fetch nucleotide sequences by accession."""

from __future__ import annotations

import logging
import time
from io import StringIO
from pathlib import Path
from typing import List

from Bio import Entrez, SeqIO

logger = logging.getLogger(__name__)

# NCBI rate limits: 3 req/s without key, 10 req/s with key
_DELAY_NO_KEY = 0.34
_DELAY_WITH_KEY = 0.11


class EntrezFetch:
    """Fetch nucleotide sequences from NCBI Entrez as FASTA.

    Example::

        fetcher = EntrezFetch(email="user@example.com")
        fetcher.fetch_fasta(["NC_001560", "NC_001542"], Path("refs.fa"))
    """

    def __init__(self, email: str, api_key: str = "", batch_size: int = 200):
        if not email:
            raise ValueError("NCBI Entrez requires an email address (set runtime.entrez_email)")
        Entrez.email = email
        if api_key:
            Entrez.api_key = api_key
        self.batch_size = batch_size
        self._delay = _DELAY_WITH_KEY if api_key else _DELAY_NO_KEY

    def fetch_fasta(self, accessions: List[str], output: Path) -> Path:
        """Fetch *accessions* from NCBI nucleotide db and write FASTA to *output*.

        Returns *output*. Skips accessions that fail (logs a warning per batch).
        """
        output.parent.mkdir(parents=True, exist_ok=True)
        written = 0
        with open(output, "w") as fh:
            for i in range(0, len(accessions), self.batch_size):
                batch = accessions[i:i + self.batch_size]
                try:
                    handle = Entrez.efetch(db="nucleotide", id=batch, rettype="fasta", retmode="text")
                    text = handle.read()
                    handle.close()
                    for record in SeqIO.parse(StringIO(text), "fasta"):
                        SeqIO.write(record, fh, "fasta")
                        written += 1
                    logger.debug("Fetched batch %d-%d (%d records)", i, i + len(batch), len(batch))
                except Exception as exc:
                    logger.warning("Entrez fetch failed for batch %d: %s", i // self.batch_size, exc)
                if i + self.batch_size < len(accessions):
                    time.sleep(self._delay)
        logger.info("Fetched %d sequence(s) from NCBI -> %s", written, output)
        return output
