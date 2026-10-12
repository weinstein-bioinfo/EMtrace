"""Expectation-Maximization abundance estimation.

Numpy-vectorized port of the Perl EM algorithm from
``Lazypipe::SeqAn`` (``EM_loop``, ``EM_step``, ``EM_logL``,
``EM_q_t_logprob``).

The EM estimates taxon abundance F(t) from alignment bitscores.
Bitscores are converted to alignment probabilities using
temperature-scaled softmax (Frith 2020,
https://doi.org/10.1093/bioinformatics/btz576).

Main entry point::

    from emtrace.analysis.em import em_loop
    result = em_loop(annot_df, temp=1.0)
    # result.F        — dict taxid → probability
    # result.log_L    — final log-likelihood
    # result.F_iter   — dict iteration → {taxid: prob}
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Dict, Optional

import numpy as np
import pandas as pd

from emtrace.analysis.taxonomy_utils import filter_no_species

logger = logging.getLogger(__name__)


# ------------------------------------------------------------------
# Result container
# ------------------------------------------------------------------

@dataclass
class EMResult:
    """Container for EM algorithm output."""

    F: Dict[str, float]
    log_L: float
    n_iterations: int
    F_iter: Dict[int, Dict[str, float]] = field(
        default_factory=dict,
    )
    # Hit structure of the EM input, per read (query):
    n_queries: int = 0   # reads entering EM
    n_single: int = 0    # one candidate taxon: contributes a plain count
    n_tied: int = 0      # several taxa, all equal score: EM can only split ties
    n_ranked: int = 0    # several taxa, unequal scores: EM weighs alternatives

    def readn_em(self) -> Dict[str, float]:
        """Expected reads per taxon, F(t) x reads in the EM input.

        Exact for the uniform read prior (``q_prob=None``).
        """
        return {t: f * self.n_queries for t, f in self.F.items()}


def hit_structure(
    q_t_logprob: Dict[str, Dict[str, float]],
) -> Dict[str, int]:
    """Count reads with one, tied, or ranked candidate taxa.

    Taxa within a read are tied when their log P(q|t) are equal, i.e. their
    best bitscores are equal (temperature scaling preserves equality).
    """
    n_single = n_tied = n_ranked = 0
    for t_probs in q_t_logprob.values():
        if len(t_probs) <= 1:
            n_single += 1
            continue
        vals = list(t_probs.values())
        if max(vals) - min(vals) <= 1e-12:
            n_tied += 1
        else:
            n_ranked += 1
    return {
        "n_queries": len(q_t_logprob),
        "n_single": n_single,
        "n_tied": n_tied,
        "n_ranked": n_ranked,
    }


def log_hit_structure(
    result: EMResult,
    log: logging.Logger,
    filter_tophits: bool,
) -> bool:
    """Log the EM input's hit structure; warn when EM cannot weigh alternatives.

    EM resolves ambiguity only through reads with *ranked* candidate taxa.
    When no read has one, its estimates are plain read counts (single-taxon
    reads) or equal splits of tied reads, and ``em.temp`` has no effect.
    Returns True if a warning was emitted.
    """
    n = result.n_queries
    if n == 0:
        return False
    log.info(
        "EM input: %d reads - single-taxon %.1f%%, tied-only %.1f%%, "
        "ranked alternatives %.1f%%",
        n, 100.0 * result.n_single / n, 100.0 * result.n_tied / n,
        100.0 * result.n_ranked / n,
    )
    if result.n_ranked > 0:
        return False
    if filter_tophits:
        log.warning(
            "annotation.filter_tophits is on (--tophits): EM input holds only "
            "top-scoring hits, so EM only splits exact ties and em.temp has no "
            "effect. Drop --tophits / set annotation.filter_tophits: false to let "
            "EM weigh alternatives."
        )
    else:
        log.warning(
            "No read has ranked alternative hits (single or tied only): EM "
            "estimates reduce to read counting and tie splitting. Check "
            "max_target_seqs / minimap2 -N."
        )
    return True


# ------------------------------------------------------------------
# Score → log-probability conversion
# ------------------------------------------------------------------

def build_q_t_logprob(
    df: pd.DataFrame,
    *,
    qseqid_col: str = "qseqid",
    staxid_col: str = "staxid",
    bitscore_col: str = "bitscore",
    temp: float = 1.0,
) -> dict:
    """Build log P(q|t) matrix from annotation bitscores.

    For each query *q* and taxon *t*, keeps the maximum bitscore
    across all alignments, then converts to probabilities via
    temperature-scaled softmax in log-space.

    Parameters
    ----------
    df : pd.DataFrame
        Annotation table with qseqid, staxid, bitscore columns.
    temp : float
        Temperature for softmax (default 1.0).

    Returns
    -------
    dict
        ``{qseqid: {staxid: log_prob, ...}, ...}``
        Also sets attributes on the returned dict:
        - ``taxa``: sorted list of unique taxids
        - ``queries``: sorted list of unique query ids
    """
    required = [qseqid_col, staxid_col, bitscore_col]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(
            f"Missing columns: {missing}"
        )

    # Keep max bitscore per (query, taxid) pair
    sub = df[[qseqid_col, staxid_col, bitscore_col]].copy()
    sub[bitscore_col] = pd.to_numeric(
        sub[bitscore_col], errors="coerce",
    )
    sub = sub.dropna(subset=[bitscore_col])
    sub[staxid_col] = sub[staxid_col].astype(str)

    max_bits = (
        sub.groupby([qseqid_col, staxid_col])[bitscore_col]
        .max()
        .reset_index()
    )

    # Temperature-scale the bitscores
    max_bits["scaled"] = max_bits[bitscore_col] / temp

    # Build nested dict and compute log-softmax per query
    q_t_logprob: Dict[str, Dict[str, float]] = {}

    for q, grp in max_bits.groupby(qseqid_col):
        scores = grp["scaled"].values
        taxids = grp[staxid_col].values

        # Log-sum-exp trick for numerical stability
        s_max = scores.max()
        log_denom = s_max + np.log(
            np.sum(np.exp(scores - s_max))
        )

        q_dict: Dict[str, float] = {}
        for t, s in zip(taxids, scores):
            q_dict[str(t)] = float(s - log_denom)
        q_t_logprob[str(q)] = q_dict

    logger.info(
        "Built log P(q|t): %d queries, %d taxa",
        len(q_t_logprob),
        len(
            {t for qd in q_t_logprob.values() for t in qd}
        ),
    )
    return q_t_logprob


# ------------------------------------------------------------------
# Log-likelihood
# ------------------------------------------------------------------

def em_log_likelihood(
    q_t_logprob: Dict[str, Dict[str, float]],
    F: Dict[str, float],
) -> float:
    """Compute log-likelihood: sum_q log sum_t P(q|t)*F(t).

    Uses log-sum-exp trick for numerical stability.
    """
    log_L = 0.0
    for q, t_probs in q_t_logprob.items():
        taxids = sorted(t_probs.keys())
        log_p = np.array([t_probs[t] for t in taxids])
        f_vals = np.array([F.get(t, 0.0) for t in taxids])

        # log(F) with -inf for zero
        log_f = np.where(
            f_vals > 0, np.log(f_vals), -np.inf,
        )
        log_pf = log_p + log_f

        # log-sum-exp
        max_lpf = log_pf.max()
        if np.isfinite(max_lpf):
            log_L += max_lpf + np.log(
                np.sum(np.exp(log_pf - max_lpf))
            )

    return float(log_L)


# ------------------------------------------------------------------
# Single EM step (E + M)
# ------------------------------------------------------------------

def em_step(
    q_t_logprob: Dict[str, Dict[str, float]],
    q_logprob: Dict[str, float],
    F: Dict[str, float],
) -> Dict[str, float]:
    """One EM iteration: compute posterior P(t|q), update F(t).

    E-step: P(t|q) = P(q|t)*F(t) / sum_t' P(q|t')*F(t')
    M-step: F(t)   = sum_q P(t|q)*P(q)
    """
    # E-step: compute log P(t|q) for all (t, q) pairs
    t_q_logprob: Dict[str, Dict[str, float]] = {}

    for q, t_probs in q_t_logprob.items():
        taxids = sorted(t_probs.keys())
        log_p = np.array([t_probs[t] for t in taxids])
        f_vals = np.array([F.get(t, 0.0) for t in taxids])
        log_f = np.where(
            f_vals > 0, np.log(f_vals), -np.inf,
        )
        log_pf = log_p + log_f

        # Bayes denominator via log-sum-exp
        max_lpf = log_pf.max()
        log_denom = max_lpf + np.log(
            np.sum(np.exp(log_pf - max_lpf))
        )

        for i, t in enumerate(taxids):
            log_post = log_pf[i] - log_denom
            if t not in t_q_logprob:
                t_q_logprob[t] = {}
            t_q_logprob[t][q] = float(log_post)

    # M-step: F(t) = sum_q exp(log P(t|q) + log P(q))
    F_new: Dict[str, float] = {}
    for t, q_posts in t_q_logprob.items():
        f_update = 0.0
        for q, log_post in q_posts.items():
            log_pq = q_logprob.get(q, 0.0)
            f_update += np.exp(log_post + log_pq)
        F_new[t] = float(f_update)

    return F_new


# ------------------------------------------------------------------
# Main EM loop
# ------------------------------------------------------------------

def em_loop(
    df: pd.DataFrame,
    *,
    qseqid_col: str = "qseqid",
    staxid_col: str = "staxid",
    bitscore_col: str = "bitscore",
    temp: float = 1.0,
    max_iter: int = 1000,
    log_L_threshold: float = 0.01,
    exclude_no_species: bool = False,
    q_prob: Optional[Dict[str, float]] = None,
    F_init: Optional[Dict[str, float]] = None,
    prune_min_reads: float = 0.0,
) -> EMResult:
    """Run the full EM algorithm on an annotation table.

    Parameters
    ----------
    df : pd.DataFrame
        Annotation table (must have qseqid, staxid, bitscore).
    temp : float
        Temperature for bitscore → probability (default 1.0).
    max_iter : int
        Maximum EM iterations (default 1000).
    log_L_threshold : float
        Convergence: stop when logL improvement < threshold.
    exclude_no_species : bool
        If True, exclude hits whose ``species`` rank is unresolved —
        empty/NA, or a non-binomial name such as ``"Vitreoscilla sp."``.
        Requires a ``species`` column; without one the filter no-ops
        with a warning.
    q_prob : dict, optional
        Prior P(q) for each query. If None, uniform 1/|Q|.
    F_init : dict, optional
        Warm-start F(t) from a previous EM run. New taxa
        are assigned a small initial probability. If None,
        F(t) is initialised uniformly to 1/|T|.
    prune_min_reads : float
        Two-pass EM, as in Emu (Curry et al. 2022). When > 0, after the
        first EM converges every taxon with fewer expected reads
        (F(t) x reads) is removed, each read's candidates are restricted to
        the surviving taxa (reads left with none leave the EM input), and
        the EM is run again to convergence from the pruned F(t), so the
        removed taxa's reads are re-assigned. 0 (default) = one pass.

    Returns
    -------
    EMResult
        Final F(t) abundances, log-likelihood, and per-iteration
        history.
    """
    work = df.copy()

    # Optionally drop hits with no resolved binomial species.
    #
    # This tests the `species` rank of the assigned taxpath, NOT `sname`.
    # For minimap2 hits `sname` is the subject accession, not an organism
    # name (see paf_to_annotation), so the old sname-based regex silently
    # matched nothing on every minimap2-based run.
    if exclude_no_species:
        work = filter_no_species(work, context="EM input")

    if work.empty:
        logger.warning("No annotation data for EM")
        return EMResult(F={}, log_L=0.0, n_iterations=0)

    # Build log P(q|t) matrix
    q_t_lp = build_q_t_logprob(
        work,
        qseqid_col=qseqid_col,
        staxid_col=staxid_col,
        bitscore_col=bitscore_col,
        temp=temp,
    )

    if not q_t_lp:
        return EMResult(F={}, log_L=0.0, n_iterations=0)

    # Collect all taxa
    all_taxa = sorted(
        {t for qd in q_t_lp.values() for t in qd}
    )
    n_taxa = len(all_taxa)

    # Init F(t)
    if F_init is not None:
        # Warm-start: reuse previous F(t), assign small
        # probability to new taxa not seen before
        F: Dict[str, float] = {}
        new_taxa = [
            t for t in all_taxa if t not in F_init
        ]
        if new_taxa:
            # Redistribute a small fraction to new taxa
            frac = min(
                0.01 * len(new_taxa) / n_taxa, 0.1,
            )
            scale = 1.0 - frac
            new_p = frac / len(new_taxa)
            for t in all_taxa:
                if t in F_init:
                    F[t] = F_init[t] * scale
                else:
                    F[t] = new_p
        else:
            F = {t: F_init.get(t, 0.0) for t in all_taxa}
        # Normalise
        total_f = sum(F.values())
        if total_f > 0:
            F = {t: v / total_f for t, v in F.items()}
        logger.info(
            "EM warm-start from %d prior taxa "
            "(%d new)",
            len(F_init), len(new_taxa),
        )
    else:
        F = {t: 1.0 / n_taxa for t in all_taxa}
    F_iter: Dict[int, Dict[str, float]] = {0: dict(F)}

    # Init P(q) — uniform if not provided
    if q_prob is None:
        n_q = len(q_t_lp)
        q_logprob = {
            q: np.log(1.0 / n_q) for q in q_t_lp
        }
    else:
        q_logprob = {
            q: np.log(p) for q, p in q_prob.items()
            if p > 0
        }

    # Initial log-likelihood
    log_L = em_log_likelihood(q_t_lp, F)
    logger.info(
        "EM init: %d queries, %d taxa, logL=%.4f",
        len(q_t_lp), n_taxa, log_L,
    )

    F, log_L, iteration = _iterate(
        q_t_lp, q_logprob, F, log_L, F_iter, 0, max_iter, log_L_threshold,
    )

    if prune_min_reads > 0:
        n_q = len(q_t_lp)
        keep = {t for t, f in F.items() if f * n_q >= prune_min_reads}
        q_t_lp = {
            q: {t: lp for t, lp in tp.items() if t in keep}
            for q, tp in q_t_lp.items()
        }
        q_t_lp = {q: tp for q, tp in q_t_lp.items() if tp}
        F = {t: F[t] for t in keep}
        total_f = sum(F.values())
        F = {t: v / total_f for t, v in F.items()} if total_f > 0 else F
        if q_prob is None:
            q_logprob = {q: np.log(1.0 / len(q_t_lp)) for q in q_t_lp}
        else:
            z = sum(q_prob.get(q, 0.0) for q in q_t_lp)
            q_logprob = {
                q: np.log(q_prob[q] / z) for q in q_t_lp
                if q_prob.get(q, 0.0) > 0
            }
        logger.info(
            "EM second pass: kept %d of %d taxa (>= %g expected reads); "
            "%d of %d reads keep a candidate",
            len(keep), n_taxa, prune_min_reads, len(q_t_lp), n_q,
        )
        log_L = em_log_likelihood(q_t_lp, F)
        F, log_L, iteration = _iterate(
            q_t_lp, q_logprob, F, log_L, F_iter, iteration, max_iter,
            log_L_threshold,
        )

    return EMResult(
        F=F,
        log_L=log_L,
        n_iterations=iteration,
        F_iter=F_iter,
        **hit_structure(q_t_lp),
    )


def _iterate(
    q_t_lp: Dict[str, Dict[str, float]],
    q_logprob: Dict[str, float],
    F: Dict[str, float],
    log_L: float,
    F_iter: Dict[int, Dict[str, float]],
    start: int,
    max_iter: int,
    log_L_threshold: float,
) -> tuple[Dict[str, float], float, int]:
    """EM iterations until the log-likelihood gain < threshold (or max_iter).

    Appends each iteration's F(t) to ``F_iter``, numbered from ``start + 1``.
    Returns ``(F, log_L, last iteration number)``.
    """
    iteration = start
    for i in range(max_iter):
        iteration += 1

        F_new = em_step(q_t_lp, q_logprob, F)
        log_L_new = em_log_likelihood(q_t_lp, F_new)
        log_L_diff = log_L_new - log_L

        logger.debug(
            "EM iter %d: logL=%.6f, diff=%.6f",
            iteration, log_L_new, log_L_diff,
        )

        log_L = log_L_new
        F = F_new
        F_iter[iteration] = dict(F)

        if log_L_diff < log_L_threshold and i > 1:
            logger.info(
                "EM converged at iteration %d "
                "(logL=%.4f, diff=%.6f)",
                iteration, log_L, log_L_diff,
            )
            break
    else:
        logger.info(
            "EM reached max_iter=%d (logL=%.4f)",
            max_iter, log_L,
        )
    return F, log_L, iteration
