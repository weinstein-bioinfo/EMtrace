"""Statistics helpers — merged into ReportStep.

All helper functions and constants are consumed by
``ReportStep`` in ``steps/report.py``.  The former ``StatsStep``
class has been removed; this module now serves as a pure utility
library.

Note: matplotlib is an optional dependency (``pip install emtrace[plots]``).
If not available, histogram generation is skipped with a warning.
"""

from __future__ import annotations

import csv
import io
import statistics
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd

from emtrace.core.config import TaxonomyType
from emtrace.io.fastq import iter_fastq

try:
    import yaml as _yaml
    _YAML_AVAILABLE = True
except ImportError:
    _YAML_AVAILABLE = False


# ------------------------------------------------------------------
# FASTQ stats helper
# ------------------------------------------------------------------

def _fastq_stats(path: Path) -> Dict[str, object]:
    """Compute read-level statistics for a FASTQ file.

    Returns dict with keys:
        file, read_count, min_len, max_len, median_len, mean_len, total_bp, n50
    """
    lengths: List[int] = []
    for _, seq, _, _ in iter_fastq(path):
        lengths.append(len(seq))

    if not lengths:
        return {
            "file": path.name,
            "read_count": 0,
            "min_len": 0,
            "max_len": 0,
            "median_len": 0,
            "mean_len": 0.0,
            "total_bp": 0,
            "n50": 0,
        }

    lengths.sort(reverse=True)
    total = sum(lengths)
    half = total / 2
    cumsum = 0
    n50 = 0
    for ln in lengths:
        cumsum += ln
        if cumsum >= half:
            n50 = ln
            break

    return {
        "file": path.name,
        "read_count": len(lengths),
        "min_len": lengths[-1],
        "max_len": lengths[0],
        "median_len": int(statistics.median(lengths)),
        "mean_len": round(statistics.mean(lengths), 1),
        "total_bp": total,
        "n50": n50,
    }


def _read_lengths(path: Path) -> List[int]:
    """Extract read lengths from a FASTQ file."""
    return [len(seq) for _, seq, _, _ in iter_fastq(path)]


# ------------------------------------------------------------------
# BAM stats helper (samtools flagstat + coverage)
# ------------------------------------------------------------------

def _bam_stats(
    bam_path: Path,
    num_threads: int = 1,
) -> Dict[str, object]:
    """Compute assembly-level stats from a BAM file.

    Returns dict with keys:
        bam, mapped_reads, mean_depth, coverage_pct
    Uses samtools flagstat and samtools coverage.
    """
    from emtrace.wrappers.samtools import Samtools

    st = Samtools(num_threads=num_threads)

    # flagstat → mapped reads count
    flagstat_out = st.flagstat(bam_path)
    mapped = 0
    for line in flagstat_out.splitlines():
        if "mapped (" in line and "primary" not in line:
            mapped = int(line.split()[0])
            break

    # coverage → per-contig stats
    cov_out = st.coverage(bam_path)
    mean_depth = 0.0
    coverage_pct = 0.0
    if cov_out.strip():
        reader = csv.DictReader(
            io.StringIO(cov_out), delimiter="\t",
        )
        total_bases = 0
        covered_bases = 0
        weighted_depth = 0.0
        for row in reader:
            rlen = ( int(row["endpos"]) - int(row["startpos"]) + 1 )
            total_bases += rlen
            covered_bases += int(row["covbases"])
            weighted_depth += float(row["meandepth"]) * rlen
        if total_bases > 0:
            mean_depth = round( weighted_depth / total_bases, 2 )
            coverage_pct = round( covered_bases / total_bases * 100, 2 )

    return {
        "bam": bam_path.name,
        "mapped_reads": mapped,
        "mean_depth": mean_depth,
        "coverage_pct": coverage_pct,
    }


# ------------------------------------------------------------------
# Division read counts from abundance table
# ------------------------------------------------------------------

def _division_counts(abund_path: Path) -> List[Dict[str, object]]:
    """Read division read counts from abund_table.tsv (NCBI taxonomy)."""
    if not abund_path.exists():
        return []
    df = pd.read_csv(abund_path, sep="\t", dtype=str, keep_default_na=False)
    if "division" not in df.columns or "readn" not in df.columns:
        return []
    df["readn"] = pd.to_numeric(df["readn"], errors="coerce").fillna(0)
    grouped = (
        df.groupby("division")["readn"]
        .sum()
        .sort_values(ascending=False)
    )
    return [
        {"division": div, "read_count": int(n)}
        for div, n in grouped.items()
        if str(div).strip()
    ]


# ------------------------------------------------------------------
# Histogram plotting
# ------------------------------------------------------------------

def _plot_histograms(panels: List[Tuple[str, List[int]]], output_path: Path) -> bool:
    """Create a multi-panel read-length histogram.

    Parameters
    ----------
    panels : list of (label, lengths)
        Each entry produces one subplot.
    output_path : Path
        Output PNG path.

    Returns True on success, False if matplotlib unavailable.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return False

    n = len(panels)
    if n == 0:
        return False

    fig, axes = plt.subplots(
        1, n, figsize=(5 * n, 4), squeeze=False,
    )
    for i, (label, lengths) in enumerate(panels):
        ax = axes[0][i]
        if lengths:
            ax.hist(
                lengths, bins=50,
                color="#4a86c8", edgecolor="white", linewidth=0.3,
            )
            ax.set_xlabel("Read length (bp)")
            ax.set_ylabel("Count")
        ax.set_title(label)
    plt.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=150)
    plt.close(fig)
    return True


# ------------------------------------------------------------------
# HTML report
# ------------------------------------------------------------------

def _generate_html_report(
    fastq_rows: List[Dict[str, object]],
    division_rows: List[Dict[str, object]],
    assembly_rows: List[Dict[str, object]],
    fig_path: Optional[Path],
    output_path: Path,
    stacked_bar_path: Optional[Path] = None,
) -> None:
    """Write a self-contained HTML stats report."""
    import base64

    def _embed_png(path: Optional[Path], heading: str) -> str:
        if not path or not path.exists():
            return ""
        b64 = base64.b64encode(path.read_bytes()).decode("ascii")
        return (
            f'<h2>{heading}</h2>\n'
            '<img src="data:image/png;base64,'
            + b64
            + '" style="max-width:100%;">\n'
        )

    fig_html = _embed_png(fig_path, "Read Length Distribution")
    stacked_bar_html = _embed_png(stacked_bar_path, "Taxonomic Composition")

    # Build HTML tables
    def _table(headers, rows):
        lines = ["<table>", "<tr>"]
        for h in headers:
            lines.append(f"<th>{h}</th>")
        lines.append("</tr>")
        for row in rows:
            lines.append("<tr>")
            for h in headers:
                val = row.get(h, row.get(h.lower(), ""))
                lines.append(f"<td>{val}</td>")
            lines.append("</tr>")
        lines.append("</table>")
        return "\n".join(lines)

    fastq_headers = [
        "file", "read_count", "min_len", "max_len",
        "median_len", "mean_len", "total_bp", "n50",
    ]
    div_headers = ["division", "read_count"]
    asm_headers = [
        "bam", "mapped_reads", "mean_depth", "coverage_pct",
    ]

    parts = [
        "<!DOCTYPE html>",
        "<html><head>",
        "<meta charset='utf-8'>",
        "<title>EMtrace Stats Report</title>",
        "<style>",
        "body { font-family: sans-serif; margin: 2em; }",
        "table { border-collapse: collapse; margin: 1em 0; }",
        "th, td { border: 1px solid #ccc; padding: 4px 10px;"
        " text-align: right; }",
        "th { background: #f0f0f0; }",
        "h1 { color: #333; }",
        "h2 { color: #555; margin-top: 1.5em; }",
        "</style>",
        "</head><body>",
        "<h1>EMtrace Stats Report</h1>",
    ]

    if fastq_rows:
        parts.append("<h2>FASTQ Statistics</h2>")
        parts.append(_table(fastq_headers, fastq_rows))

    if division_rows:
        parts.append(
            "<h2>Read Counts by Taxonomic Division</h2>"
        )
        parts.append(_table(div_headers, division_rows))

    if stacked_bar_html:
        parts.append(stacked_bar_html)

    if assembly_rows:
        parts.append("<h2>Assembly Statistics</h2>")
        parts.append(_table(asm_headers, assembly_rows))

    if fig_html:
        parts.append(fig_html)

    parts.extend(["</body>", "</html>"])

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(parts), encoding="utf-8")


_FASTQ_HEADERS = ["stage", "file", "read_count", "min_len", "max_len", "median_len", "mean_len", "total_bp", "n50"]
_DIV_HEADERS   = ["division", "read_count"]
_ASM_HEADERS   = ["bam", "mapped_reads", "mean_depth", "coverage_pct"]
