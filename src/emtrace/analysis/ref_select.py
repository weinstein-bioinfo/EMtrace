"""Automatic reference selection for reference-based assembly (``assembly.ref_auto``).

Two hit sources are combined:

  primary    the assembly reads searched against ``assembly.database`` (RefSeq
             viral genomes by default) — ``annotation/assembly_candidates.tsv``,
             written by :func:`search_reference_db`, the same file the GUI's
             "Pick from viral hits" picker reads.
  secondary  the annotation hits (``annot2.tsv``, else ``annot1.tsv``), viral rows.

Only rows from nucleotide searches (``NUCLEOTIDE_SEARCHES``: minimap2, blastn)
are candidates: a protein search (Diamond, blastp) hits protein accessions,
which cannot be fetched as nucleotide references.  A chained round
(e.g. ``minimap.refseq.abv.1g,diamond.uniref.vi``) mixes both in one table;
its protein rows are dropped and recorded in the selection table
(:func:`split_nucleotide_hits`).

Picks are made in each source independently (:func:`pick_references`):

  1. hits are grouped by species (``species``, else ``staxid``, else accession),
  2. species are kept when they reach ``min_reads`` and ``min_read_frac``,
  3. accessions are chosen by greedy read cover over all kept species: each
     round takes the accession with the highest alignment quality summed over
     the reads not yet covered, as long as it recruits at least ``min_reads``
     of them.  Segments recruit different reads, so each is picked; strains of
     one segment share reads, so only the best is.  The cover is global, so a
     separately-named taxon recruiting the same reads (e.g. a reassortant) is
     not picked either.  Caps: ``max_refs_per_species``, ``max_species``.

The sources are then combined by one rule (:func:`combine_picks`): a primary
pick is replaced by the secondary pick sharing the most reads with it only if
the secondary aligns better, on those shared reads, by more than
``secondary_margin``.  Secondary picks that share fewer than ``min_reads`` reads
with every primary pick are added as they are.

Alignment quality (``quality_metric``), per read, best hit per accession:

  qcov_identity  matched bases / read length × 100 = pident × alen / qseqlen.
                 Rewards identity *and* the whole read aligned; comparable
                 across aligners.  Margin in percentage points.
  pident_qcov    pident × qcov.  Margin in percentage points.
  bitscore       summed bitscore.  Margin in percent (relative).
  bits_per_base  bitscore / read length.  Margin in percent (relative).
  reads_won      percent of shared reads on which a side has the higher
                 qcov_identity.  Margin in percentage points.

pident alone is deliberately not offered: a short local alignment can have high
identity without the target being the closer homolog.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import pandas as pd

logger = logging.getLogger(__name__)

QUALITY_METRICS = ("qcov_identity", "pident_qcov", "bitscore", "bits_per_base", "reads_won")
_RELATIVE_METRICS = ("bitscore", "bits_per_base")

# Search programs whose subjects are nucleotide sequences (annotation ``search``
# column, lower-cased).  An allowlist: an unknown search type is not a candidate.
NUCLEOTIDE_SEARCHES = ("minimap2", "blastn")

SELECTION_COLUMNS = [
    "species", "staxid", "accession", "source", "db", "reads", "new_reads",
    "shared_reads", "quality_primary", "quality_secondary", "selected", "reason",
]


# -------------------------------------------------------------------
# Search of assembly.database
# -------------------------------------------------------------------

def search_reference_db(
    reads: Path,
    db_entry,
    num_threads: int,
    k: int = 0,
    w: int = 0,
    log_file: Optional[Path] = None,
) -> pd.DataFrame:
    """Map *reads* to the reference DB *db_entry* and return annotation-format hits.

    Runs minimap2 ``map-ont`` with base-level alignment (``-c``) and secondary
    alignments kept, so ``pident``, ``qcov`` and an ``AS``-based bitscore are
    real and every hit per read is available to the greedy read cover.  Taxids
    come from the ``.acc2taxid`` sidecar beside the DB FASTA.  All rows are
    marked ``division = Viruses`` (the assembly DB is viral by convention).

    A pre-built ``.mmi`` index is used when present, unless *k*/*w* are set:
    an index fixes k and w, so those need the FASTA (indexed on the fly, slower).
    """
    from emtrace.models.annotation import empty_annotation_df, paf_to_annotation
    from emtrace.steps.annotate import _apply_acc2taxid
    from emtrace.wrappers.minimap2 import Minimap2, parse_paf

    db_path = Path(db_entry.db)
    fasta = Path(str(db_path)[:-4]) if str(db_path).endswith(".mmi") else db_path
    mmi = fasta.with_name(fasta.name + ".mmi")
    if k or w:
        effective_db = fasta
    elif str(db_path).endswith(".mmi") and db_path.exists():
        effective_db = db_path
    elif mmi.exists():
        effective_db = mmi
    else:
        effective_db = fasta
    if not effective_db.exists():
        raise FileNotFoundError(f"Reference DB not found: {effective_db}")

    hits = Minimap2(num_threads=num_threads).align(
        db=effective_db, query=reads, preset="map-ont",
        c=True, secondary="yes", n=5, k=k, w=w, log_file=log_file,
    )
    if isinstance(hits, Path):
        hits = parse_paf(hits, c=True)
    if hits.empty:
        return empty_annotation_df()

    annot = paf_to_annotation(hits, search="Minimap2", db=db_entry.name, c=True)
    acc2taxid = fasta.with_name(fasta.name + ".acc2taxid")
    if acc2taxid.exists():
        _apply_acc2taxid(annot, acc2taxid)
    annot["division"] = "Viruses"
    return annot


def candidates_are_current(path: Path, reads: Path) -> bool:
    """True if *path* can be reused as the primary hits for *reads*.

    It must be at least as new as the reads and carry base-level alignment
    fields: files written before the search used ``-c`` have an empty ``qcov``
    and chain-score bitscores, and must be regenerated.
    """
    if not path.exists() or path.stat().st_mtime < reads.stat().st_mtime:
        return False
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, nrows=50)
    if df.empty:
        return True  # a search that found nothing
    return "qcov" in df.columns and (df["qcov"].str.strip() != "").any()


# -------------------------------------------------------------------
# Per-read alignment quality
# -------------------------------------------------------------------

def _read_len(df: pd.DataFrame) -> pd.Series:
    """Read length: qseqlen, else aligned query length / qcov."""
    qseqlen = pd.to_numeric(df.get("qseqlen"), errors="coerce")
    alt = pd.to_numeric(df.get("qlen"), errors="coerce") / pd.to_numeric(df.get("qcov"), errors="coerce")
    return qseqlen.where(qseqlen > 0, alt)


def row_quality(df: pd.DataFrame, metric: str) -> pd.Series:
    """Per-hit alignment quality for *metric* (``reads_won`` scores as qcov_identity)."""
    pident = pd.to_numeric(df["pident"], errors="coerce")
    bits = pd.to_numeric(df["bitscore"], errors="coerce")
    if metric in ("qcov_identity", "reads_won"):
        alen = pd.to_numeric(df["alen"], errors="coerce")
        q = (pident * alen / _read_len(df)).clip(upper=100)
    elif metric == "pident_qcov":
        q = pident * pd.to_numeric(df["qcov"], errors="coerce")
    elif metric == "bitscore":
        q = bits
    elif metric == "bits_per_base":
        q = bits / _read_len(df)
    else:
        raise ValueError(f"unknown quality_metric: {metric}")
    return q.fillna(0.0)


def split_nucleotide_hits(df: Optional[pd.DataFrame],
                          source: str) -> Tuple[Optional[pd.DataFrame], List[dict]]:
    """Keep the rows of nucleotide searches; report the others.

    Rows with an empty ``search`` (or a table without the column) are kept.
    Returns ``(kept, rejected_rows)``: one rejected row (SELECTION_COLUMNS) per
    dropped ``(db, search)`` pair, with its number of distinct reads.
    """
    if df is None or df.empty or "search" not in df.columns:
        return df, []
    search = df["search"].fillna("").astype(str).str.strip()
    keep = (search == "") | search.str.lower().isin(NUCLEOTIDE_SEARCHES)
    if keep.all():
        return df, []
    dropped = df[~keep].assign(search=search[~keep])
    db = dropped["db"] if "db" in dropped.columns else pd.Series("", index=dropped.index)
    dropped = dropped.assign(db=db.fillna("").astype(str).str.strip())
    rejected: List[dict] = []
    for (db_name, prog), grp in dropped.groupby(["db", "search"]):
        n = grp["qseqid"].nunique()
        logger.info("ref_auto: ignored protein hits of %d read(s) from %s (%s)",
                    n, db_name or "?", prog)
        rejected.append(_row("", "", "", source, db_name, n,
                             reason=f"protein search ({prog}): accessions are not "
                                    "nucleotide references"))
    return df[keep], rejected


def prepare_hits(df: Optional[pd.DataFrame], metric: str) -> pd.DataFrame:
    """Reduce raw hits to one row per (read, accession) with a group key and quality."""
    cols = ["qseqid", "sseqid", "db", "species", "staxid", "group", "bitscore", "quality"]
    if df is None or df.empty:
        return pd.DataFrame(columns=cols)
    df = df.copy()
    if "division" in df.columns and (df["division"].astype(str).str.strip() != "").any():
        df = df[df["division"].astype(str).str.strip().str.lower() == "viruses"]
    df = df[df["sseqid"].astype(str).str.strip() != ""]
    if df.empty:
        return pd.DataFrame(columns=cols)
    for c in ("species", "staxid", "db"):
        if c not in df.columns:
            df[c] = ""
        df[c] = df[c].fillna("").astype(str).str.strip()
    taxid = df["staxid"].where(~df["staxid"].isin(["", "NA", "0"]), "")
    df["group"] = df["species"].where(df["species"] != "",
                                      ("taxid:" + taxid).where(taxid != "", "acc:" + df["sseqid"]))
    df["bitscore"] = pd.to_numeric(df["bitscore"], errors="coerce").fillna(0.0)
    df["quality"] = row_quality(df, metric)
    # best hit per (read, accession): secondaries to the same target collapse
    df = df.sort_values("quality", ascending=False).drop_duplicates(["qseqid", "sseqid"])
    return df[cols].reset_index(drop=True)


# -------------------------------------------------------------------
# Picking within one source
# -------------------------------------------------------------------

@dataclass
class Pick:
    """One chosen accession."""
    accession: str
    group: str
    species: str
    staxid: str
    db: str
    source: str
    reads: Set[str]
    new_reads: int


def pick_references(hits: pd.DataFrame, cfg, source: str) -> Tuple[List[Pick], List[dict]]:
    """Species ranking + greedy read cover over prepared *hits*.

    Returns ``(picks, rejected_rows)``; rejected rows follow SELECTION_COLUMNS.
    """
    picks: List[Pick] = []
    rejected: List[dict] = []
    if hits.empty:
        return picks, rejected
    n_total = hits["qseqid"].nunique()

    sp = (hits.groupby("group")
              .agg(reads=("qseqid", "nunique"), bits=("bitscore", "sum"))
              .sort_values(["reads", "bits"], ascending=False))
    kept = sp[(sp["reads"] >= cfg.min_reads) & (sp["reads"] / n_total >= cfg.min_read_frac)]
    for g, row in sp.drop(kept.index).iterrows():
        first = hits[hits["group"] == g].iloc[0]
        rejected.append(_row(first["species"], first["staxid"], "", source, "", int(row["reads"]),
                             reason=f"species below min_reads={cfg.min_reads} / "
                                    f"min_read_frac={cfg.min_read_frac}"))
    # Greedy read cover, global across species: an accession of another
    # species that recruits the same reads (e.g. a reassortant with its own
    # taxid) adds no new reads and is not picked.
    kh = hits[hits["group"].isin(kept.index)]
    reads_by_acc = kh.groupby("sseqid")["qseqid"].apply(set).to_dict()
    qual = kh.set_index(["sseqid", "qseqid"])["quality"]
    meta = kh.drop_duplicates("sseqid").set_index("sseqid")
    covered: Set[str] = set()
    chosen: List[str] = []
    per_species: Dict[str, int] = {}
    while True:
        best, best_score, best_new = None, 0.0, 0
        for acc, rs in reads_by_acc.items():
            if acc in chosen:
                continue
            g = meta.loc[acc, "group"]
            if per_species.get(g, 0) >= cfg.max_refs_per_species:
                continue
            if g not in per_species and len(per_species) >= cfg.max_species:
                continue
            new = rs - covered
            if len(new) < cfg.min_reads:
                continue
            score = float(qual.loc[acc].loc[list(new)].sum())
            if score > best_score:
                best, best_score, best_new = acc, score, len(new)
        if best is None:
            break
        m = meta.loc[best]
        chosen.append(best)
        covered |= reads_by_acc[best]
        per_species[m["group"]] = per_species.get(m["group"], 0) + 1
        picks.append(Pick(best, m["group"], m["species"], m["staxid"], m["db"], source,
                          reads_by_acc[best], best_new))
    for acc, rs in reads_by_acc.items():
        if acc in chosen:
            continue
        m = meta.loc[acc]
        new = len(rs - covered)
        if new >= cfg.min_reads:
            reason = (f"beyond max_refs_per_species={cfg.max_refs_per_species} "
                      f"or max_species={cfg.max_species}")
        elif chosen:
            reason = "reads already covered by a chosen accession"
        else:
            reason = "no accession reached min_reads"
        rejected.append(_row(m["species"], m["staxid"], acc, source, m["db"], len(rs),
                             new_reads=new, reason=reason))
    return picks, rejected


# -------------------------------------------------------------------
# Combining primary and secondary picks
# -------------------------------------------------------------------

def compare_on_shared(hits_p: pd.DataFrame, acc_p: str, hits_s: pd.DataFrame, acc_s: str,
                      shared: Set[str], metric: str) -> Tuple[float, float]:
    """Quality of each side on the *shared* reads, in the metric's comparison units."""
    qp = hits_p[hits_p["sseqid"] == acc_p].set_index("qseqid")["quality"].reindex(list(shared)).fillna(0.0)
    qs = hits_s[hits_s["sseqid"] == acc_s].set_index("qseqid")["quality"].reindex(list(shared)).fillna(0.0)
    if metric == "reads_won":
        n = max(len(shared), 1)
        return 100.0 * float((qp > qs).sum()) / n, 100.0 * float((qs > qp).sum()) / n
    if metric == "bitscore":
        return float(qp.sum()), float(qs.sum())
    return float(qp.mean()), float(qs.mean())


def secondary_wins(q_p: float, q_s: float, metric: str, margin: float) -> bool:
    """True if the secondary side is better by more than *margin*."""
    if metric in _RELATIVE_METRICS:
        return q_p > 0 and 100.0 * (q_s - q_p) / q_p > margin or (q_p <= 0 < q_s)
    return q_s - q_p > margin


def combine_picks(primary: List[Pick], secondary: List[Pick],
                  hits_p: pd.DataFrame, hits_s: pd.DataFrame, cfg) -> Tuple[List[Pick], List[dict]]:
    """Apply the margin rule; returns ``(final_picks, rows)`` with a row per pick considered."""
    metric, margin = cfg.quality_metric, cfg.secondary_margin
    final: List[Pick] = []
    rows: List[dict] = []
    matched: Set[int] = set()
    for p in primary:
        best_i, best_shared = None, set()
        for i, s in enumerate(secondary):
            shared = p.reads & s.reads
            if len(shared) > len(best_shared):
                best_i, best_shared = i, shared
        if best_i is None or len(best_shared) < cfg.min_reads:
            final.append(p)
            rows.append(_row_pick(p, selected=True, reason="primary pick; no matching secondary pick"))
            continue
        s = secondary[best_i]
        matched.add(best_i)
        if s.accession == p.accession:
            final.append(p)
            rows.append(_row_pick(p, selected=True, shared=len(best_shared),
                                  reason="primary pick; secondary picked the same accession"))
            continue
        q_p, q_s = compare_on_shared(hits_p, p.accession, hits_s, s.accession, best_shared, metric)
        win = secondary_wins(q_p, q_s, metric, margin)
        verdict = (f"{metric}: secondary {q_s:.2f} vs primary {q_p:.2f}, margin {margin}")
        rows.append(_row_pick(p, selected=not win, shared=len(best_shared), q_p=q_p, q_s=q_s,
                              reason=("replaced; " if win else "kept; ") + verdict))
        rows.append(_row_pick(s, selected=win, shared=len(best_shared), q_p=q_p, q_s=q_s,
                              reason=("replaces primary " + p.accession + "; " if win
                                      else "not better by the margin; ") + verdict))
        final.append(s if win else p)

    n_species = len({f.group for f in final})
    for i, s in enumerate(secondary):
        if i in matched:
            continue
        overlaps = [p for p in primary if len(p.reads & s.reads) >= cfg.min_reads]
        if overlaps:
            rows.append(_row_pick(s, selected=False,
                                  reason=f"overlaps primary pick {overlaps[0].accession}"))
            continue
        new_species = s.group not in {f.group for f in final}
        if new_species and n_species >= cfg.max_species:
            rows.append(_row_pick(s, selected=False, reason=f"beyond max_species={cfg.max_species}"))
            continue
        final.append(s)
        n_species += new_species
        rows.append(_row_pick(s, selected=True,
                              reason="secondary pick; not found by the primary search"))
    return final, rows


# -------------------------------------------------------------------
# Entry point
# -------------------------------------------------------------------

def select_references(primary: Optional[pd.DataFrame], secondary: Optional[pd.DataFrame],
                      cfg) -> Tuple[List[Pick], pd.DataFrame]:
    """Select reference accessions from the two hit tables (either may be None).

    *cfg* is an ``assembly.ref_auto`` config.  Returns the final picks and the
    selection table (SELECTION_COLUMNS) to write as ``ref_selection.tsv``.
    """
    metric = cfg.quality_metric
    primary, prot_p = split_nucleotide_hits(primary, "primary")
    secondary, prot_s = split_nucleotide_hits(secondary, "secondary")
    hp = prepare_hits(primary, metric)
    hs = prepare_hits(secondary, metric)
    picks_p, rej_p = pick_references(hp, cfg, "primary")
    picks_s, rej_s = pick_references(hs, cfg, "secondary")
    if picks_p:
        final, rows = combine_picks(picks_p, picks_s, hp, hs, cfg)
    else:
        final = picks_s
        rows = [_row_pick(s, selected=True,
                          reason="secondary pick; no primary picks") for s in picks_s]
    table = pd.DataFrame(rows + rej_p + rej_s + prot_p + prot_s, columns=SELECTION_COLUMNS)
    logger.info("Reference selection: %s", ", ".join(p.accession for p in final) or "none")
    return final, table


def _row(species, staxid, accession, source, db, reads, new_reads="", shared="",
         q_p="", q_s="", selected=False, reason="") -> Dict[str, object]:
    return dict(zip(SELECTION_COLUMNS, [species, staxid, accession, source, db, reads, new_reads,
                                        shared, q_p, q_s, selected, reason]))


def _row_pick(p: Pick, selected: bool, reason: str, shared="", q_p="", q_s="") -> Dict[str, object]:
    fmt = lambda v: round(v, 3) if isinstance(v, float) else v  # noqa: E731
    return _row(p.species, p.staxid, p.accession, p.source, p.db, len(p.reads), p.new_reads,
                shared, fmt(q_p), fmt(q_s), selected, reason)
