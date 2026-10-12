"""EM abundance estimation step.

Runs the Expectation-Maximization algorithm on annotation tables
to estimate per-taxon abundance probabilities F(t).

Input:
    ``annot2.tsv`` (preferred) or ``annot1.tsv`` from the annotation
    step — whichever has data.

Output (under ``{sample.reports_dir}/``):
    ``EM.taxid_prob.tsv``          — final F(t) by taxid
    ``EM.taxid_prob.byiter.tsv``   — F(t) history per iteration
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import pandas as pd

from emtrace.analysis.em import em_loop, log_hit_structure
from emtrace.models.annotation import read_annotation_tsv
from emtrace.steps import StepName, register_step
from emtrace.steps.base import PipelineStep
from emtrace.utils.paths import ensure_dir

if TYPE_CHECKING:
    from emtrace.core.pipelinejob import PipelineJob


@register_step(StepName.EM)
class EMAbundanceStep(PipelineStep):
    """Estimate taxon abundances via EM on alignment bitscores."""

    step_name = "em"
    step_description = "EM abundance estimation"

    def validate_inputs(
        self, sample: PipelineJob, **kwargs,
    ) -> bool:
        """EM has no hard prerequisites.

        A missing or empty annotation table is a valid edge case (e.g.
        no reads were identified) and is handled as a graceful skip in
        ``run()`` rather than a validation failure.
        """
        return True

    def run(self, sample: PipelineJob, **kwargs) -> PipelineJob:
        reports_dir = sample.reports_dir
        ensure_dir(reports_dir)

        # Pick best annotation table — absence is a valid no-op skip
        annot_path = self._pick_annotation(sample)
        if annot_path is None:
            self.logger.warning(
                "No annotation table found — skipping EM abundance "
                "(no taxa identified)"
            )
            return sample
        self.logger.info(
            "EM input: %s", annot_path,
        )

        annot_df = read_annotation_tsv(annot_path)
        if annot_df.empty:
            self.logger.warning(
                "Annotation table is empty — skipping EM"
            )
            return sample

        # EM parameters from config or kwargs
        temp = float(
            kwargs.get("score2prob_temp", 1.0)
        )
        max_iter = int(kwargs.get("maxiter", 1000))
        log_L_diff = float(
            kwargs.get("logLdiff", 0.01)
        )
        # `--nospnorep` is applied to cfg.em.exclude_no_species by the CLI,
        # so it arrives through config below rather than through kwargs.
        no_sp = False

        # If config has EM parameters, prefer those
        em_cfg = getattr(self.config, "em", None)
        if em_cfg is not None:
            temp = getattr(em_cfg, "temp", temp)
            max_iter = getattr(
                em_cfg, "max_iter", max_iter,
            )
            log_L_diff = getattr(
                em_cfg, "log_L_threshold", log_L_diff,
            )
            no_sp = getattr(
                em_cfg, "exclude_no_species", no_sp,
            )
        prune = getattr(em_cfg, "prune_min_reads", 0.0) if em_cfg is not None else 0.0

        self.logger.info(
            "EM params: temp=%.2f, max_iter=%d, "
            "logL_threshold=%.4f, exclude_sp=%s, prune_min_reads=%g",
            temp, max_iter, log_L_diff, no_sp, prune,
        )

        # Run EM
        result = em_loop(
            annot_df,
            temp=temp,
            max_iter=max_iter,
            log_L_threshold=log_L_diff,
            exclude_no_species=no_sp,
            prune_min_reads=prune,
        )

        self.logger.info(
            "EM finished: %d iterations, logL=%.4f, "
            "%d taxa",
            result.n_iterations,
            result.log_L,
            len(result.F),
        )

        log_hit_structure(
            result, self.logger, self.config.annotation.filter_tophits,
        )

        # Write EM.taxid_prob.tsv (sorted by prob descending)
        prob_path = reports_dir / "EM.taxid_prob.tsv"
        self._write_taxid_prob(result.F, prob_path, n_queries=result.n_queries)
        sample.set_file("em_taxid_prob", prob_path)

        # Write EM.taxid_prob.byiter.tsv
        byiter_path = (
            reports_dir / "EM.taxid_prob.byiter.tsv"
        )
        self._write_byiter(
            result.F, result.F_iter, byiter_path,
        )
        sample.set_file("em_taxid_prob_byiter", byiter_path)

        return sample

    # -----------------------------------------------------------
    # Helpers
    # -----------------------------------------------------------

    def _pick_annotation(self, sample: PipelineJob):
        """Return the best annotation table path, or None.

        Checks {output_dir}/annotation/annot2.tsv first (non-empty),
        then annot1.tsv.  Emits a warning for each missing candidate.
        """
        ann_dir = sample.annotation_dir

        for name in ("annot2.tsv", "annot1.tsv"):
            p = ann_dir / name
            if not p.exists():
                continue
            try:
                df = pd.read_csv(p, sep="\t", nrows=1)
                if len(df) > 0:
                    return p
            except Exception:
                continue

        self.logger.warning(
            "No annotation table found under %s — skipping EM abundance", ann_dir
        )
        return None

    @staticmethod
    def _write_taxid_prob(
        F: dict, path, n_queries: int = 0,
    ) -> None:
        """Write taxid, prob, readn_em TSV sorted by prob desc.

        ``readn_em`` = prob x reads in the EM input (expected reads); "NA"
        when the read count is unknown.
        """
        sorted_taxa = sorted(
            F.items(), key=lambda x: x[1], reverse=True,
        )
        with open(path, "w", encoding="utf-8") as fh:
            fh.write("taxid\tprob\treadn_em\n")
            for taxid, prob in sorted_taxa:
                readn_em = prob * n_queries if n_queries > 0 else "NA"
                fh.write(f"{taxid}\t{prob}\t{readn_em}\n")

    @staticmethod
    def _write_byiter(
        F_final: dict,
        F_iter: dict,
        path,
    ) -> None:
        """Write per-iteration F(t) table.

        Rows = taxa (sorted by final prob descending).
        Columns = iterations (0, 1, 2, …).
        """
        iters = sorted(F_iter.keys())
        sorted_taxa = sorted(
            F_final.items(),
            key=lambda x: x[1],
            reverse=True,
        )
        with open(path, "w", encoding="utf-8") as fh:
            header = "taxid\t" + "\t".join(
                str(i) for i in iters
            )
            fh.write(header + "\n")
            for taxid, _ in sorted_taxa:
                vals = [
                    str(F_iter[i].get(taxid, 0.0))
                    for i in iters
                ]
                fh.write(
                    f"{taxid}\t" + "\t".join(vals) + "\n"
                )
