"""Abundance visualization step.

Generates standard metagenomic abundance visualizations from abund_table.tsv:

  - Krona text data file (tab-separated hierarchy for ktImportText)
  - Sankey diagram  (Plotly — multi-level taxonomy flow, interactive HTML)
  - Stacked bar chart — side-by-side bars per top-level group (division),
    stacked by species/genus within each group (matplotlib PNG + Plotly HTML)
  - Sunburst (Plotly — hierarchical radial drill-down, interactive HTML)
  - Treemap (Plotly — optional, off by default)
  - Heatmaps at species and genus level (PDF only, PyComplexHeatmap / seaborn)

Design decisions:
  - Plotly is the primary engine for interactive HTML outputs (sankey,
    sunburst, stacked bar).  These are displayed inside the app via
    QWebEngineView (PyQt6-WebEngine); falling back to the system browser
    when QWebEngineView is unavailable.
  - pySankey / pySankey2 are NOT used: unmaintained and limited to 2 levels.
  - Static PNG/PDF export uses kaleido for Plotly figures.
  - matplotlib handles the publication-quality stacked bar PNG.
  - BioSankey (Xu et al., 2019) pattern: rank-pair aggregation for Sankey.

Stacked bar design:
  One bar per top-level group (e.g. division), stacked segments show species
  (or genus) composition within that group — equivalent to PhyloSeq's
  plot_bar() composition plot.

Taxonomy rank chain (NCBI and SILVA):
  division → class → order → family → genus → species

Output (under ``{sample.reports_dir}/plots/``):
    ``krona.data.txt``          — Krona text input (ktImportText)
    ``sankey.html``             — interactive Sankey diagram
    ``sankey.pdf``              — static Sankey (requires kaleido)
    ``stacked_bar.html``        — interactive stacked bar (Plotly)
    ``stacked_bar.png``         — static stacked bar (matplotlib)
    ``sunburst.html``           — interactive sunburst (Plotly)
    ``sunburst.pdf``            — static sunburst (requires kaleido)
    ``heatmap_species.pdf``     — species-level clustered heatmap
    ``heatmap_genus.pdf``       — genus-level clustered heatmap
    ``coverage.png/.pdf``       — per-reference read-depth (pysam + matplotlib)
    ``coverage.html``           — interactive coverage (Plotly)
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING, List, Optional

import pandas as pd

from emtrace.steps import StepName, register_step
from emtrace.steps.base import PipelineStep
from emtrace.utils.paths import ensure_dir

if TYPE_CHECKING:
    from emtrace.core.pipelinejob import PipelineJob

logger = logging.getLogger("emtrace.steps.plots")

# ---------------------------------------------------------------------------
# Taxonomy rank chain — order matters (broad → narrow)
# ---------------------------------------------------------------------------
_RANK_CHAIN: List[str] = [
    "division", "class", "order", "family", "genus", "species",
]

# ---------------------------------------------------------------------------
# Color palette (ColorBrewer Set2 + Tableau10, colorblind-friendly)
# ---------------------------------------------------------------------------
_RANK_COLORS = [
    "#4e79a7", "#f28e2b", "#e15759", "#76b7b2",
    "#59a14f", "#edc948", "#b07aa1", "#ff9da7",
    "#9c755f", "#bab0ac",
]
_LINK_COLOR = "rgba(180,180,180,0.35)"


def _hex_to_rgba(hex_color: str, alpha: float) -> str:
    h = hex_color.lstrip("#")
    r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
    return f"rgba({r},{g},{b},{alpha:.3f})"


def _hex_to_rgba_mpl(hex_color: str, alpha: float) -> tuple:
    h = hex_color.lstrip("#")
    return (int(h[0:2], 16) / 255, int(h[2:4], 16) / 255, int(h[4:6], 16) / 255, alpha)


# ===========================================================================
# Krona data (moved from report.py)
# ===========================================================================

def generate_krona_data(
    abund: pd.DataFrame,
    output: Path,
    toptaxrank: str = "division",
) -> None:
    """Write Krona text input from an abundance DataFrame.

    Format (port of ``abundtable2krona.pl``):
        readn<TAB>rank1<TAB>rank2<TAB>...

    Compatible with ``ktImportText`` from the KronaTools suite.
    """
    rank_cols = [toptaxrank, "class", "order", "family", "genus", "species"]
    rank_cols = [c for c in rank_cols if c in abund.columns]

    output.parent.mkdir(parents=True, exist_ok=True)
    with open(output, "w", encoding="utf-8") as fh:
        for _, row in abund.iterrows():
            readn = row.get("readn", 0)
            parts: List[str] = []
            for col in rank_cols:
                val = str(row.get(col, "")).strip()
                if val and val not in ("NA", "nan", ""):
                    parts.append(val)
            if parts:
                fh.write(f"{readn}\t" + "\t".join(parts) + "\n")
    logger.info("Krona data -> %s", output)


# ===========================================================================
# Sankey diagram
# ===========================================================================

def _build_sankey_links(
    abund: pd.DataFrame,
    rank_chain: List[str],
    value_col: str = "readn",
):
    """Build node list and link arrays for a Plotly Sankey figure.

    Returns (node_labels, display_labels, node_colors, src_idxs, tgt_idxs, values).
    Nodes are namespaced as ``rank:value`` internally to avoid collisions when
    the same name appears at two different taxonomy levels.
    """
    # Keep only ranks present in the DataFrame with non-trivial values
    active_ranks: List[str] = []
    for rank in rank_chain:
        if rank in abund.columns and abund[rank].replace("", pd.NA).dropna().shape[0] > 0:
            active_ranks.append(rank)
    if len(active_ranks) < 2:
        return None

    def _clean(val) -> str:
        s = str(val).strip()
        return s if s and s not in ("NA", "nan", "") else "unknown"

    source_ns: List[str] = []
    target_ns: List[str] = []
    link_values: List[float] = []

    for i in range(len(active_ranks) - 1):
        src_col, tgt_col = active_ranks[i], active_ranks[i + 1]
        grp = (
            abund
            .assign(
                _src=abund[src_col].apply(_clean),
                _tgt=abund[tgt_col].apply(_clean),
            )
            .groupby(["_src", "_tgt"], as_index=False)[value_col]
            .sum()
        )
        for _, row in grp.iterrows():
            source_ns.append(f"{src_col}:{row['_src']}")
            target_ns.append(f"{tgt_col}:{row['_tgt']}")
            link_values.append(float(row[value_col]))

    # Ordered unique node list
    seen: dict = {}
    for label in source_ns + target_ns:
        if label not in seen:
            seen[label] = len(seen)
    node_ns = list(seen.keys())
    node_idx = seen

    # Display labels: strip rank prefix
    display_labels = [n.split(":", 1)[1] for n in node_ns]

    # Color each node by its rank level
    rank_color = {rank: _RANK_COLORS[i % len(_RANK_COLORS)] for i, rank in enumerate(active_ranks)}
    node_colors = [rank_color.get(n.split(":")[0], "#aaaaaa") for n in node_ns]

    return (
        node_ns,
        display_labels,
        node_colors,
        [node_idx[s] for s in source_ns],
        [node_idx[t] for t in target_ns],
        link_values,
    )


def plot_sankey(
    abund: pd.DataFrame,
    output_dir: Path,
    toptaxrank: str = "division",
    value_col: str = "readn",
    title: str = "Taxonomic Abundance — Sankey",
    top_n: int = 50,
) -> None:
    """Generate a multi-level Sankey diagram from an abundance DataFrame.

    Outputs:
        sankey.html  — interactive Plotly HTML
        sankey.pdf   — static PDF (requires kaleido)

    Parameters
    ----------
    abund : pd.DataFrame
        Abundance table with taxonomy columns (division … species) and
        a read-count column (readn).
    output_dir : Path
        Directory to write output files.
    toptaxrank : str
        Broadest taxonomy rank to start the flow (default: "division").
    value_col : str
        Column to use as flow magnitude (default: "readn").
    top_n : int
        Keep only the top N taxa by readn at the broadest rank to avoid
        overcrowded diagrams (default: 50).
    """
    try:
        import plotly.graph_objects as go
    except ImportError:
        logger.warning("plotly not installed — skipping Sankey diagram")
        return

    # Trim to top_n to keep diagram readable
    if toptaxrank in abund.columns and len(abund) > top_n:
        top_taxa = abund.nlargest(top_n, value_col)[toptaxrank].unique()
        abund = abund[abund[toptaxrank].isin(top_taxa)]

    chain = [r for r in _RANK_CHAIN if r != toptaxrank]
    chain = [toptaxrank] + chain  # ensure toptaxrank is first

    result = _build_sankey_links(abund, chain, value_col)
    if result is None:
        logger.warning("Sankey: need at least 2 non-empty taxonomy ranks — skipping")
        return

    _, display_labels, node_colors, src_idxs, tgt_idxs, values = result

    fig = go.Figure(go.Sankey(
        arrangement="snap",
        node=dict(
            pad=12,
            thickness=18,
            line=dict(color="white", width=0.4),
            label=display_labels,
            color=node_colors,
        ),
        link=dict(
            source=src_idxs,
            target=tgt_idxs,
            value=values,
            color=_LINK_COLOR,
        ),
    ))
    fig.update_layout(
        title_text=title,
        font=dict(size=11, family="Arial"),
        height=max(600, len(display_labels) * 14),
        margin=dict(l=20, r=20, t=50, b=20),
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    html_path = output_dir / "sankey.html"
    fig.write_html(str(html_path))
    logger.info("Sankey HTML -> %s", html_path)

    _write_static(fig, pdf_path=output_dir / "sankey.pdf")


# ===========================================================================
# Stacked bar chart
# ===========================================================================

def plot_stacked_bar(
    abund: pd.DataFrame,
    output_dir: Path,
    toptaxrank: str = "division",
    value_col: str = "readn",
    label: str = "Sample",
    top_n_groups: int = 10,
    top_n_inner: int = 15,
    normalize: bool = False,
    inner_rank: Optional[str] = None,
    prob_col: Optional[str] = "prob",
) -> None:
    """Generate a multi-level stacked bar chart of taxonomic composition.

    Shows one bar per top-level group (e.g. division), stacked by the inner
    rank (species → genus → family, first present column).  This is the
    standard metagenomics composition plot equivalent to PhyloSeq's plot_bar().

    Outputs:
        stacked_bar.png  — matplotlib static PNG
        stacked_bar.html — interactive Plotly HTML

    Parameters
    ----------
    top_n_groups : int
        Number of top-level groups (X-axis bars) to show.
    top_n_inner : int
        Number of inner-rank taxa shown as stacked segments; rest → "Other".
    inner_rank : str, optional
        Column to stack within each group.  Auto-detected: species → genus →
        family → order (first present, different from toptaxrank).
    normalize : bool
        Show relative abundance (%) per group instead of read counts.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    if toptaxrank not in abund.columns:
        logger.warning("Stacked bar: column '%s' not found — skipping", toptaxrank)
        return

    # Auto-detect inner rank
    if inner_rank is None:
        for candidate in ("species", "genus", "family", "order", "class"):
            if candidate in abund.columns and candidate != toptaxrank:
                inner_rank = candidate
                break
    if inner_rank is None or inner_rank not in abund.columns:
        logger.warning("Stacked bar: no inner rank column found — skipping")
        return

    has_prob = bool(
        prob_col and prob_col in abund.columns and abund[prob_col].gt(0).any()
    )

    def _clean(s) -> str:
        v = str(s).strip()
        return v if v and v not in ("NA", "nan", "") else "unknown"

    def _agg_with_prob(df, keys):
        if has_prob:
            df = df.copy()
            df["_pw"] = df[value_col] * df[prob_col]
            out = df.groupby(keys, as_index=False).agg({value_col: "sum", "_pw": "sum"})
            out[prob_col] = out["_pw"] / out[value_col].replace(0, 1)
            return out.drop(columns=["_pw"])
        return df.groupby(keys, as_index=False)[value_col].sum()

    work = abund.copy()
    work["_grp"]   = work[toptaxrank].apply(_clean)
    work["_inner"] = work[inner_rank].apply(_clean)

    # Sum by (group, inner)
    agg = _agg_with_prob(work, ["_grp", "_inner"])

    # Keep top N groups by total readn
    grp_totals = agg.groupby("_grp")[value_col].sum().nlargest(top_n_groups)
    top_grps   = grp_totals.index.tolist()
    agg = agg[agg["_grp"].isin(top_grps)].copy()

    # Keep top N inner taxa globally; collapse rest to "Other"
    inner_totals = agg.groupby("_inner")[value_col].sum().nlargest(top_n_inner)
    top_inner    = set(inner_totals.index)
    agg["_inner"] = agg["_inner"].where(agg["_inner"].isin(top_inner), "Other")
    agg = _agg_with_prob(agg, ["_grp", "_inner"])

    # Pivot: rows = inner taxa, columns = groups
    pivot = agg.pivot_table(
        index="_inner", columns="_grp", values=value_col, fill_value=0,
    )
    pivot = pivot[[g for g in top_grps if g in pivot.columns]]

    # "Other" at bottom of legend/stack
    main_idx  = [i for i in pivot.index if i != "Other"]
    other_idx = [i for i in pivot.index if i == "Other"]
    pivot = pivot.loc[main_idx + other_idx]

    # Build prob pivot before normalising readn counts
    if has_prob:
        pivot_prob = agg.pivot_table(
            index="_inner", columns="_grp", values=prob_col, fill_value=0,
        )
        pivot_prob = pivot_prob[[g for g in top_grps if g in pivot_prob.columns]]
        pivot_prob = pivot_prob.reindex(pivot.index, fill_value=0)
    else:
        pivot_prob = None

    if normalize:
        col_sums = pivot.sum(axis=0).replace(0, 1)
        pivot = pivot.div(col_sums, axis=1) * 100
        ylabel = "Relative abundance (%)"
    else:
        ylabel = "Read count"

    # Threshold for showing p= text labels (7% of max column total)
    segment_threshold = float(pivot.sum(axis=0).max()) * 0.07 if pivot_prob is not None else 0.0

    taxa_list = pivot.index.tolist()
    grp_list  = pivot.columns.tolist()
    colors    = [_RANK_COLORS[i % len(_RANK_COLORS)] for i in range(len(taxa_list))]

    # --- matplotlib stacked bar ---
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np

        x       = np.arange(len(grp_list))
        bottoms = np.zeros(len(grp_list))
        fig, ax = plt.subplots(figsize=(max(5, len(grp_list) * 1.2 + 3), 6))
        for i, taxon in enumerate(taxa_list):
            vals = pivot.loc[taxon].values.astype(float)
            ax.bar(x, vals, bottom=bottoms, color=colors[i],
                   label=taxon, edgecolor="white", linewidth=0.4)
            if pivot_prob is not None:
                probs = pivot_prob.loc[taxon].values.astype(float)
                for j, (v, p) in enumerate(zip(vals, probs)):
                    if v >= segment_threshold and p > 0:
                        ax.text(
                            x[j], bottoms[j] + v / 2, f"p={p:.2f}",
                            ha="center", va="center",
                            fontsize=6, color="white", fontweight="bold",
                        )
            bottoms += vals
        ax.set_xticks(x)
        ax.set_xticklabels(grp_list, rotation=35, ha="right", fontsize=9)
        ax.set_ylabel(ylabel, fontsize=11)
        ax.set_xlabel(toptaxrank.capitalize(), fontsize=11)
        ax.set_title(f"Taxonomic Composition — {label}", fontsize=12)
        ax.legend(
            title=inner_rank.capitalize(),
            bbox_to_anchor=(1.01, 1), loc="upper left",
            fontsize=8, title_fontsize=9,
        )
        plt.tight_layout()
        png_path = output_dir / "stacked_bar.png"
        plt.savefig(png_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        logger.info("Stacked bar PNG -> %s", png_path)
    except ImportError:
        logger.warning("matplotlib not available — skipping static stacked bar")

    # --- Plotly interactive stacked bar ---
    try:
        import plotly.graph_objects as go

        fig = go.Figure()
        y_unit = "%" if normalize else "reads"
        for i, taxon in enumerate(taxa_list):
            y_vals = pivot.loc[taxon].tolist()
            if pivot_prob is not None:
                probs = pivot_prob.loc[taxon].tolist()
                text_labels = [
                    f"p={p:.2f}" if y_vals[j] >= segment_threshold and p > 0 else ""
                    for j, p in enumerate(probs)
                ]
                fig.add_trace(go.Bar(
                    name=taxon,
                    x=grp_list,
                    y=y_vals,
                    marker_color=colors[i],
                    customdata=[[p] for p in probs],
                    text=text_labels,
                    textposition="inside",
                    textfont=dict(size=9, color="white"),
                    hovertemplate=(
                        f"<b>%{{x}}</b> — {taxon}<br>"
                        f"%{{y:,.1f}} {y_unit}<br>"
                        "P(EM): %{customdata[0]:.3f}"
                        "<extra></extra>"
                    ),
                ))
            else:
                fig.add_trace(go.Bar(
                    name=taxon,
                    x=grp_list,
                    y=y_vals,
                    marker_color=colors[i],
                ))
        fig.update_layout(
            barmode="stack",
            title=f"Taxonomic Composition — {label}",
            xaxis_title=toptaxrank.capitalize(),
            yaxis_title=ylabel,
            legend_title=inner_rank.capitalize(),
            font=dict(size=11, family="Arial"),
            margin=dict(l=60, r=20, t=50, b=120),
        )
        html_path = output_dir / "stacked_bar.html"
        fig.write_html(str(html_path))
        logger.info("Stacked bar HTML -> %s", html_path)
    except ImportError:
        logger.warning("plotly not available — skipping interactive stacked bar")


# ===========================================================================
# Sunburst and treemap
# ===========================================================================

def plot_sunburst(
    abund: pd.DataFrame,
    output_dir: Path,
    toptaxrank: str = "division",
    value_col: str = "readn",
    title: str = "Taxonomic Abundance — Sunburst",
    top_n: int = 100,
    prob_col: Optional[str] = "prob",
) -> None:
    """Generate a Plotly sunburst chart (radial hierarchical drill-down).

    This is the closest Python-native equivalent to Krona's interactive
    radial chart.

    Outputs:
        sunburst.html — interactive HTML
        sunburst.pdf  — static PDF (requires kaleido)
    """
    try:
        import plotly.express as px
    except ImportError:
        logger.warning("plotly not installed — skipping sunburst diagram")
        return

    if top_n and len(abund) > top_n:
        abund = abund.nlargest(top_n, value_col)

    chain = [r for r in _RANK_CHAIN if r == toptaxrank or _RANK_CHAIN.index(r) > _RANK_CHAIN.index(toptaxrank)]
    chain = [r for r in chain if r in abund.columns]
    if not chain:
        logger.warning("Sunburst: no taxonomy columns found — skipping")
        return

    has_prob = bool(
        prob_col and prob_col in abund.columns and abund[prob_col].gt(0).any()
    )

    # Clean up NA / empty values
    extra_cols = [prob_col] if has_prob else []
    work = abund[chain + [value_col] + extra_cols].copy()
    for col in chain:
        work[col] = work[col].fillna("unknown").replace("", "unknown")
    if has_prob:
        work[prob_col] = pd.to_numeric(work[prob_col], errors="coerce").fillna(0)

    output_dir.mkdir(parents=True, exist_ok=True)

    fig = px.sunburst(
        work,
        path=chain,
        values=value_col,
        title=title,
        color_discrete_sequence=_RANK_COLORS,
        hover_data={prob_col: ":.3f"} if has_prob else None,
    )
    fig.update_layout(font=dict(size=11, family="Arial"), margin=dict(t=50, l=0, r=0, b=0))

    html_path = output_dir / "sunburst.html"
    fig.write_html(str(html_path))
    logger.info("Sunburst HTML -> %s", html_path)

    _write_static(fig, pdf_path=output_dir / "sunburst.pdf")


def plot_treemap(
    abund: pd.DataFrame,
    output_dir: Path,
    toptaxrank: str = "division",
    value_col: str = "readn",
    title: str = "Taxonomic Abundance — Treemap",
    top_n: int = 100,
) -> None:
    """Generate a Plotly treemap (space-filling hierarchical chart).

    Complementary to the sunburst — shows the same hierarchy with
    area proportional to abundance, which makes relative sizes easier
    to compare at a glance.

    Outputs:
        treemap.html — interactive HTML
        treemap.png  — static PNG (requires kaleido)
    """
    try:
        import plotly.express as px
    except ImportError:
        logger.warning("plotly not installed — skipping treemap")
        return

    if top_n and len(abund) > top_n:
        abund = abund.nlargest(top_n, value_col)

    chain = [r for r in _RANK_CHAIN if r == toptaxrank or _RANK_CHAIN.index(r) > _RANK_CHAIN.index(toptaxrank)]
    chain = [r for r in chain if r in abund.columns]
    if not chain:
        logger.warning("Treemap: no taxonomy columns found — skipping")
        return

    work = abund[chain + [value_col]].copy()
    for col in chain:
        work[col] = work[col].fillna("unknown").replace("", "unknown")

    output_dir.mkdir(parents=True, exist_ok=True)

    fig = px.treemap(
        work,
        path=chain,
        values=value_col,
        title=title,
        color_discrete_sequence=_RANK_COLORS,
    )
    fig.update_layout(font=dict(size=11, family="Arial"), margin=dict(t=50, l=0, r=0, b=0))

    html_path = output_dir / "treemap.html"
    fig.write_html(str(html_path))
    logger.info("Treemap HTML -> %s", html_path)

    _write_static(fig, output_dir / "treemap.png", output_dir / "treemap.pdf")


# ===========================================================================
# Annotated heatmap (PyComplexHeatmap, with seaborn/matplotlib fallback)
# ===========================================================================

def plot_heatmap(
    abund: pd.DataFrame,
    output_dir: Path,
    value_col: str = "readn",
    rank_col: Optional[str] = None,
    annot_ranks: Optional[List[str]] = None,
    top_n: int = 40,
    log_scale: bool = True,
    sample_label: str = "Sample",
    title: str = "Taxonomic Abundance — Heatmap",
    output_filename: str = "heatmap",
) -> None:
    """Clustered heatmap of top-N taxa with lineage annotation tracks.

    Primary backend: ``PyComplexHeatmap.ClusterMapPlotter`` (Python port of
    R ComplexHeatmap).  Falls back to ``seaborn.clustermap`` when
    PyComplexHeatmap is unavailable, and to ``matplotlib.imshow`` as a
    last resort.

    Output: ``{output_filename}.pdf`` (PDF only — use output_filename to
    differentiate species-level and genus-level heatmaps).
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Pick row label column: species > genus > the first present rank
    if rank_col is None:
        for candidate in ("species", "genus", "family", "division"):
            if candidate in abund.columns:
                rank_col = candidate
                break
    if not rank_col or rank_col not in abund.columns:
        logger.warning("Heatmap: no usable row rank column — skipping")
        return

    if annot_ranks is None:
        annot_ranks = [r for r in ("division", "family") if r in abund.columns and r != rank_col]

    # Numeric values
    vals = pd.to_numeric(abund[value_col], errors="coerce").fillna(0)
    work = abund.assign(_val=vals)
    work = work[work["_val"] > 0]
    if work.empty:
        logger.warning("Heatmap: no non-zero values in '%s' — skipping", value_col)
        return

    # Deduplicate by row label (sum values, take first of annotation ranks)
    agg_map = {"_val": "sum"}
    for r in annot_ranks:
        agg_map[r] = "first"
    grouped = work.groupby(rank_col, as_index=False).agg(agg_map)
    grouped = grouped.sort_values("_val", ascending=False).head(top_n)
    if grouped.empty:
        logger.warning("Heatmap: nothing to plot after aggregation — skipping")
        return

    labels = grouped[rank_col].astype(str).tolist()
    heat_vals = grouped["_val"].to_numpy(dtype=float)
    if log_scale:
        import numpy as np
        heat_vals = np.log10(heat_vals + 1.0)
        cbar_label = f"log10({value_col}+1)"
    else:
        cbar_label = value_col

    heat_df = pd.DataFrame({sample_label: heat_vals}, index=labels)

    annot_df = None
    if annot_ranks:
        annot_df = grouped[annot_ranks].copy()
        annot_df.index = labels
        annot_df = annot_df.fillna("unknown").replace("", "unknown")

    pdf_path = output_dir / f"{output_filename}.pdf"

    # --- Primary: PyComplexHeatmap ---
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from PyComplexHeatmap import (
            ClusterMapPlotter,
            HeatmapAnnotation,
            anno_simple,
        )

        row_annotation = None
        if annot_df is not None and not annot_df.empty:
            annos = {}
            cmaps = ["Set2", "tab20", "Set3", "Paired"]
            for i, rank in enumerate(annot_ranks):
                annos[rank.capitalize()] = anno_simple(
                    annot_df[rank], cmap=cmaps[i % len(cmaps)], add_text=False,
                )
            row_annotation = HeatmapAnnotation(
                **annos, axis=0, verbose=0,
            )

        ClusterMapPlotter(
            data=heat_df,
            left_annotation=row_annotation,
            row_cluster=len(labels) > 1,
            col_cluster=False,
            cmap="YlOrRd",
            show_rownames=True,
            show_colnames=True,
            label=cbar_label,
            row_dendrogram=False,
            figsize=(5.0, max(6.0, len(labels) * 0.22)),
        )
        plt.suptitle(title, fontsize=12)
        plt.savefig(pdf_path, bbox_inches="tight")
        plt.close("all")
        logger.info("Heatmap (PyComplexHeatmap) -> %s", pdf_path)
        return
    except ImportError:
        logger.warning(
            "PyComplexHeatmap not installed — falling back to seaborn.clustermap",
        )
    except Exception as exc:
        logger.warning("PyComplexHeatmap failed (%s) — falling back", exc)
        plt.close("all")

    # --- Fallback: seaborn.clustermap ---
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import seaborn as sns

        g = sns.clustermap(
            heat_df,
            cmap="YlOrRd",
            row_cluster=len(labels) > 1,
            col_cluster=False,
            figsize=(5.0, max(6.0, len(labels) * 0.22)),
            cbar_kws={"label": cbar_label},
        )
        g.fig.suptitle(title, y=1.02, fontsize=12)
        g.savefig(pdf_path, bbox_inches="tight")
        plt.close(g.fig)
        logger.info("Heatmap (seaborn) -> %s", pdf_path)
        return
    except ImportError:
        logger.warning(
            "seaborn not installed — falling back to matplotlib.imshow",
        )
    except Exception as exc:
        logger.warning("seaborn.clustermap failed (%s) — falling back", exc)
        plt.close("all")

    # --- Last resort: matplotlib.imshow ---
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(4.0, max(6.0, len(labels) * 0.22)))
        im = ax.imshow(heat_df.values, aspect="auto", cmap="YlOrRd")
        ax.set_yticks(range(len(labels)))
        ax.set_yticklabels(labels, fontsize=8)
        ax.set_xticks([0])
        ax.set_xticklabels([sample_label])
        ax.set_title(title, fontsize=11)
        fig.colorbar(im, ax=ax, label=cbar_label)
        plt.tight_layout()
        plt.savefig(pdf_path)
        plt.close(fig)
        logger.info("Heatmap (matplotlib) -> %s", pdf_path)
    except ImportError:
        logger.warning("matplotlib not installed — skipping heatmap")


# ===========================================================================
# Per-reference read-depth coverage plot (BAM → pysam)
# ===========================================================================

def plot_coverage(
    bam_path: Path,
    output_dir: Path,
    top_n: int = 12,
    min_reads: int = 5,
    title: str = "Per-reference read depth",
) -> None:
    """Plot per-base read depth across the top-N references in a BAM.

    The single most diagnostic plot for viral confirmation: a flat
    coverage profile across a reference is strong evidence of a real
    hit, while a single spike is often artefactual.

    Outputs:
        coverage.png  (matplotlib small-multiples)
        coverage.pdf
        coverage.html (Plotly, if available)
    """
    try:
        import pysam
    except ImportError:
        logger.warning("pysam not installed — skipping coverage plot")
        return

    if not bam_path.exists():
        logger.warning("Coverage: BAM not found at %s — skipping", bam_path)
        return

    output_dir.mkdir(parents=True, exist_ok=True)

    # Pick top-N references by mapped-read count
    try:
        bam = pysam.AlignmentFile(str(bam_path), "rb")
    except Exception as exc:
        logger.warning("Coverage: cannot open %s (%s) — skipping", bam_path, exc)
        return

    try:
        stats = [
            (s.contig, bam.get_reference_length(s.contig), s.mapped)
            for s in bam.get_index_statistics()
        ]
    except Exception as exc:
        logger.warning(
            "Coverage: no BAM index (%s) — run `samtools index` first; skipping",
            exc,
        )
        bam.close()
        return

    stats = [s for s in stats if s[2] >= min_reads]
    stats.sort(key=lambda s: s[2], reverse=True)
    stats = stats[:top_n]
    if not stats:
        logger.info("Coverage: no references with >= %d reads — skipping", min_reads)
        bam.close()
        return

    # Per-reference depth arrays
    import numpy as np
    profiles: List[tuple] = []
    for ref, length, mapped in stats:
        try:
            cov = bam.count_coverage(contig=ref, quality_threshold=0)
            depth = np.asarray(cov[0], dtype=np.int64) \
                + np.asarray(cov[1], dtype=np.int64) \
                + np.asarray(cov[2], dtype=np.int64) \
                + np.asarray(cov[3], dtype=np.int64)
        except Exception as exc:
            logger.warning("Coverage: %s failed (%s) — skipping ref", ref, exc)
            continue
        covered = float((depth > 0).sum()) / max(len(depth), 1) * 100
        mean_depth = float(depth.mean()) if len(depth) else 0.0
        profiles.append((ref, length, mapped, depth, covered, mean_depth))
    bam.close()

    if not profiles:
        logger.warning("Coverage: no usable references — skipping")
        return

    # --- matplotlib small-multiples ---
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        n = len(profiles)
        fig, axes = plt.subplots(
            n, 1, figsize=(8.0, max(2.0, n * 1.4)), sharex=False,
        )
        if n == 1:
            axes = [axes]
        for ax, (ref, length, mapped, depth, covered, mean_depth) in zip(axes, profiles):
            ax.fill_between(range(len(depth)), depth, step="pre", color="#4e79a7", linewidth=0)
            ax.set_ylabel("Depth", fontsize=9)
            ax.set_title(
                f"{ref}  ({mapped} reads, {covered:.1f}% covered, mean {mean_depth:.1f}×)",
                fontsize=10, loc="left",
            )
            ax.tick_params(axis="both", labelsize=8)
        axes[-1].set_xlabel("Position (bp)")
        fig.suptitle(title, fontsize=12)
        plt.tight_layout()
        png_path = output_dir / "coverage.png"
        pdf_path = output_dir / "coverage.pdf"
        plt.savefig(png_path, dpi=150)
        plt.savefig(pdf_path)
        plt.close(fig)
        logger.info("Coverage plot -> %s", png_path)
    except ImportError:
        logger.warning("matplotlib not installed — skipping static coverage plot")

    # --- Plotly interactive ---
    try:
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots

        n = len(profiles)
        fig = make_subplots(
            rows=n, cols=1,
            subplot_titles=[
                f"{ref} ({mapped} reads, {cov:.1f}% covered, mean {md:.1f}×)"
                for ref, _, mapped, _, cov, md in profiles
            ],
            vertical_spacing=0.04,
        )
        for i, (ref, _, _, depth, _, _) in enumerate(profiles, start=1):
            fig.add_trace(
                go.Scatter(
                    x=list(range(len(depth))),
                    y=depth.tolist(),
                    mode="lines",
                    fill="tozeroy",
                    name=ref,
                    line=dict(color="#4e79a7", width=1),
                ),
                row=i, col=1,
            )
        fig.update_layout(
            title_text=title,
            height=max(300, 200 * n),
            showlegend=False,
            font=dict(size=11, family="Arial"),
            margin=dict(l=60, r=20, t=60, b=40),
        )
        fig.update_xaxes(title_text="Position (bp)", row=n, col=1)
        html_path = output_dir / "coverage.html"
        fig.write_html(str(html_path))
        logger.info("Coverage HTML -> %s", html_path)
    except ImportError:
        logger.warning("plotly not installed — skipping interactive coverage plot")


# ===========================================================================
# Static export helper
# ===========================================================================

def _write_static(fig, png_path: Optional[Path] = None, pdf_path: Optional[Path] = None) -> None:
    """Write Plotly figure to PNG and/or PDF using kaleido."""
    try:
        if png_path:
            fig.write_image(str(png_path), scale=2)
            logger.info("Static PNG -> %s", png_path)
        if pdf_path:
            fig.write_image(str(pdf_path))
            logger.info("Static PDF -> %s", pdf_path)
    except Exception as exc:
        # kaleido not installed, or other rendering error
        logger.warning(
            "Static export failed (%s) — install kaleido for PNG/PDF output", exc,
        )


# ===========================================================================
# Pipeline step
# ===========================================================================

@register_step(StepName.PLOT)
class PlotStep(PipelineStep):
    """Generate abundance visualizations (Krona, Sankey, bar, sunburst, treemap).

    Reads ``abund_table.tsv`` from the reports directory and writes all
    plot outputs to ``{reports_dir}/plots/``.

    Requires plotly for interactive/static outputs and matplotlib for
    the static stacked bar chart.  All outputs are skipped gracefully
    if their dependencies are missing.
    """

    step_name = "plot"
    step_description = "Abundance visualization (Sankey, sunburst, treemap, bar)"

    def validate_inputs(self, job: PipelineJob, **kwargs) -> bool:
        # Plots are optional output: a missing/empty abund_table.tsv (e.g. a run
        # with no annotation hits) is not a pipeline failure. Always validate;
        # run() skips gracefully so later steps (pack/clean) still execute.
        return True

    def run(self, job: PipelineJob, **kwargs) -> PipelineJob:
        abund_path = self._abund_path(job)
        if not abund_path.exists():
            self.logger.warning(
                "abund_table.tsv not found at %s — skipping plots", abund_path,
            )
            return job
        abund = pd.read_csv(abund_path, sep="\t", dtype=str, keep_default_na=False)
        if abund.empty:
            self.logger.warning("abund_table.tsv is empty — skipping plots")
            return job

        # Numeric columns
        for col in ("readn", "bitscore", "prob"):
            if col in abund.columns:
                abund[col] = pd.to_numeric(abund[col], errors="coerce").fillna(0)

        toptaxrank = self._get_toptaxrank()
        plots_dir = job.reports_dir / "plots"
        ensure_dir(plots_dir)

        label = job.name or "Sample"
        report = self.config.report

        # Krona data (disabled — superseded by sunburst.html; kept for future use)
        # krona_path = plots_dir / "krona.data.txt"
        # generate_krona_data(abund, krona_path, toptaxrank)
        # job.set_file("krona_data", krona_path)

        # 1. Sankey
        if report.sankey:
            self.logger.info("Generating Sankey diagram...")
            plot_sankey(abund, plots_dir, toptaxrank=toptaxrank)

        # 2. Stacked bar
        if report.stacked_bar:
            self.logger.info("Generating stacked bar chart...")
            plot_stacked_bar(abund, plots_dir, toptaxrank=toptaxrank, label=label)

        # 3. Sunburst
        if report.sunburst:
            self.logger.info("Generating sunburst chart...")
            plot_sunburst(abund, plots_dir, toptaxrank=toptaxrank)

        # 5. Heatmaps — species level and genus level (PDF only)
        if report.heatmap_species:
            self.logger.info("Generating species-level heatmap...")
            plot_heatmap(
                abund, plots_dir, sample_label=label,
                rank_col="species",
                annot_ranks=[r for r in ("division", "genus") if r in abund.columns],
                title="Taxonomic Abundance — Species",
                output_filename="heatmap_species",
            )
        if report.heatmap_genus:
            self.logger.info("Generating genus-level heatmap...")
            plot_heatmap(
                abund, plots_dir, sample_label=label,
                rank_col="genus",
                annot_ranks=[r for r in ("division", "family") if r in abund.columns],
                title="Taxonomic Abundance — Genus",
                output_filename="heatmap_genus",
            )

        # 7. Per-reference coverage profile (requires BAM)
        bam_path = job.get_file("assembly_bam")
        if bam_path and Path(bam_path).exists():
            self.logger.info("Generating coverage plot...")
            plot_coverage(Path(bam_path), plots_dir)
        else:
            self.logger.info("No BAM in job — skipping coverage plot")

        self.logger.info("Plots -> %s", plots_dir)
        sunburst = plots_dir / "sunburst.html"
        if sunburst.exists():
            job.set_file("plots_sunburst", sunburst)
        return job

    def _abund_path(self, job: PipelineJob) -> Path:
        p = job.get_file("abund_table")
        if p:
            return p
        return job.reports_dir / "reports" / "abund_table.tsv"

    def _get_toptaxrank(self) -> str:
        # Plots read the abundance table, where both taxonomies expose their
        # top rank as "division" (SILVA's domain is copied into it). The
        # taxonomy sources' toptaxrank ("domain" for SILVA) names an
        # annotation-table column and does not apply here.
        return "division"
