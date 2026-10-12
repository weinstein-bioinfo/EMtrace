"""Tests for analysis/em.py — EM abundance estimation algorithm.

This is the most algorithmically critical component per the ROADMAP.
Tests validate:
  - Score-to-probability conversion (softmax with temperature)
  - Single EM step (E + M)
  - Log-likelihood computation
  - Full EM loop convergence
  - Warm-start from prior F(t)
  - Edge cases (single taxon, empty input, etc.)
"""

import math

import numpy as np
import pandas as pd
import pytest

from emtrace.analysis.em import (
    EMResult,
    build_q_t_logprob,
    em_log_likelihood,
    em_loop,
    em_step,
)


# ------------------------------------------------------------------
# Fixtures
# ------------------------------------------------------------------

@pytest.fixture
def simple_annot_df():
    """Annotation table with clear taxon signal.

    read1 → taxon A (high score) and B (low score)
    read2 → taxon A (medium) and C (medium)
    read3 → taxon B (high)
    """
    return pd.DataFrame({
        "qseqid": [
            "r1", "r1", "r2", "r2", "r3",
        ],
        "staxid": [
            "A", "B", "A", "C", "B",
        ],
        "bitscore": [
            500.0, 100.0, 300.0, 280.0, 450.0,
        ],
    })


@pytest.fixture
def uniform_annot_df():
    """Annotation table where all taxa have equal scores."""
    return pd.DataFrame({
        "qseqid": ["r1", "r1", "r2", "r2"],
        "staxid": ["X", "Y", "X", "Y"],
        "bitscore": [100.0, 100.0, 100.0, 100.0],
    })


# ------------------------------------------------------------------
# Score → probability
# ------------------------------------------------------------------

class TestBuildQtLogprob:
    """Test bitscore → log P(q|t) conversion."""

    def test_probabilities_sum_to_one(self, simple_annot_df):
        lp = build_q_t_logprob(simple_annot_df)
        for q, t_dict in lp.items():
            probs = [math.exp(v) for v in t_dict.values()]
            assert abs(sum(probs) - 1.0) < 1e-10

    def test_higher_score_higher_prob(
        self, simple_annot_df,
    ):
        lp = build_q_t_logprob(simple_annot_df)
        # r1: A has bitscore 500, B has 100
        assert lp["r1"]["A"] > lp["r1"]["B"]

    def test_temperature_effect(self, simple_annot_df):
        lp_t1 = build_q_t_logprob(
            simple_annot_df, temp=1.0,
        )
        lp_t10 = build_q_t_logprob(
            simple_annot_df, temp=10.0,
        )
        # Higher temperature → more uniform probs
        # Difference between A and B should be smaller
        diff_t1 = lp_t1["r1"]["A"] - lp_t1["r1"]["B"]
        diff_t10 = lp_t10["r1"]["A"] - lp_t10["r1"]["B"]
        assert diff_t10 < diff_t1

    def test_uniform_scores(self, uniform_annot_df):
        lp = build_q_t_logprob(uniform_annot_df)
        # Equal scores → equal probabilities
        for q, t_dict in lp.items():
            probs = list(t_dict.values())
            assert abs(probs[0] - probs[1]) < 1e-10

    def test_max_bitscore_per_pair(self):
        """Only max bitscore per (q, t) should be kept."""
        df = pd.DataFrame({
            "qseqid": ["r1", "r1"],
            "staxid": ["A", "A"],
            "bitscore": [100.0, 500.0],
        })
        lp = build_q_t_logprob(df)
        # Only one taxon → probability = 1.0
        assert abs(math.exp(lp["r1"]["A"]) - 1.0) < 1e-10

    def test_missing_columns_raises(self):
        df = pd.DataFrame({"x": [1]})
        with pytest.raises(ValueError, match="Missing"):
            build_q_t_logprob(df)


# ------------------------------------------------------------------
# Log-likelihood
# ------------------------------------------------------------------

class TestEmLogLikelihood:
    """Test log-likelihood computation."""

    def test_logL_finite(self, simple_annot_df):
        lp = build_q_t_logprob(simple_annot_df)
        F = {"A": 0.5, "B": 0.3, "C": 0.2}
        ll = em_log_likelihood(lp, F)
        assert math.isfinite(ll)
        assert ll < 0  # log-likelihood is negative

    def test_logL_increases_after_step(
        self, simple_annot_df,
    ):
        lp = build_q_t_logprob(simple_annot_df)
        F = {"A": 1 / 3, "B": 1 / 3, "C": 1 / 3}
        q_logprob = {
            q: math.log(1 / 3) for q in lp
        }
        ll_before = em_log_likelihood(lp, F)
        F_new = em_step(lp, q_logprob, F)
        ll_after = em_log_likelihood(lp, F_new)
        assert ll_after >= ll_before


# ------------------------------------------------------------------
# Single EM step
# ------------------------------------------------------------------

class TestEmStep:
    """Test single E-M step."""

    def test_F_sums_approximately_to_one(
        self, simple_annot_df,
    ):
        lp = build_q_t_logprob(simple_annot_df)
        F = {"A": 1 / 3, "B": 1 / 3, "C": 1 / 3}
        q_logprob = {
            q: math.log(1 / 3) for q in lp
        }
        F_new = em_step(lp, q_logprob, F)
        total = sum(F_new.values())
        assert abs(total - 1.0) < 0.01

    def test_dominant_taxon_increases(
        self, simple_annot_df,
    ):
        lp = build_q_t_logprob(simple_annot_df)
        F = {"A": 1 / 3, "B": 1 / 3, "C": 1 / 3}
        q_logprob = {
            q: math.log(1 / 3) for q in lp
        }
        F_new = em_step(lp, q_logprob, F)
        # A and B have strong signal, C is weaker
        assert F_new["A"] > F["A"] or F_new["B"] > F["B"]


# ------------------------------------------------------------------
# Full EM loop
# ------------------------------------------------------------------

class TestEmLoop:
    """Test full EM algorithm."""

    def test_converges(self, simple_annot_df):
        result = em_loop(simple_annot_df)
        assert isinstance(result, EMResult)
        assert result.n_iterations > 0
        assert math.isfinite(result.log_L)
        assert len(result.F) == 3  # A, B, C

    def test_F_sums_to_approximately_one(
        self, simple_annot_df,
    ):
        result = em_loop(simple_annot_df)
        total = sum(result.F.values())
        assert abs(total - 1.0) < 0.05

    def test_dominant_taxon_gets_highest_prob(
        self, simple_annot_df,
    ):
        result = em_loop(simple_annot_df)
        # A and B both have strong signal
        top = max(result.F, key=result.F.get)
        assert top in ("A", "B")

    def test_iteration_history(self, simple_annot_df):
        result = em_loop(simple_annot_df, max_iter=5)
        assert 0 in result.F_iter
        assert len(result.F_iter) >= 2

    def test_max_iter_respected(self, simple_annot_df):
        result = em_loop(
            simple_annot_df, max_iter=3,
        )
        assert result.n_iterations <= 3

    def test_temperature_parameter(self, simple_annot_df):
        r1 = em_loop(simple_annot_df, temp=1.0)
        r10 = em_loop(simple_annot_df, temp=10.0)
        # Higher temp → more uniform F
        std_1 = np.std(list(r1.F.values()))
        std_10 = np.std(list(r10.F.values()))
        assert std_10 < std_1

    def test_uniform_scores_equal_F(
        self, uniform_annot_df,
    ):
        result = em_loop(uniform_annot_df)
        vals = list(result.F.values())
        assert abs(vals[0] - vals[1]) < 0.01

    def test_empty_input(self):
        df = pd.DataFrame(
            columns=["qseqid", "staxid", "bitscore"],
        )
        result = em_loop(df)
        assert result.F == {}
        assert result.n_iterations == 0

    def test_single_taxon(self):
        df = pd.DataFrame({
            "qseqid": ["r1", "r2"],
            "staxid": ["A", "A"],
            "bitscore": [100.0, 200.0],
        })
        result = em_loop(df)
        assert len(result.F) == 1
        assert abs(result.F["A"] - 1.0) < 1e-10

    def test_exclude_no_species(self):
        """Filter reads the `species` rank, not `sname`.

        `sname` is the subject accession for minimap2 hits, so it is
        deliberately unhelpful here — only `species` decides.
        """
        df = pd.DataFrame({
            "qseqid": ["r1", "r1", "r2"],
            "staxid": ["A", "B", "A"],
            "bitscore": [100.0, 200.0, 300.0],
            "sname": ["ACC1.1.500", "ACC2.1.500", "ACC1.1.500"],
            "species": [
                "Virus A", "Salmonella sp.", "Virus A",
            ],
        })
        result = em_loop(
            df, exclude_no_species=True,
        )
        # B ("sp.") should be excluded
        assert "B" not in result.F
        assert "A" in result.F


# ------------------------------------------------------------------
# Warm-start EM
# ------------------------------------------------------------------

class TestEmWarmStart:
    """Test warm-start from previous F(t)."""

    def test_warm_start_converges(self, simple_annot_df):
        # First run
        r1 = em_loop(simple_annot_df)

        # Warm-start from r1.F
        r2 = em_loop(
            simple_annot_df, F_init=r1.F,
        )
        assert r2.n_iterations > 0
        assert math.isfinite(r2.log_L)

    def test_warm_start_fewer_iterations(
        self, simple_annot_df,
    ):
        # Cold start
        r_cold = em_loop(
            simple_annot_df,
            log_L_threshold=0.001,
        )

        # Warm start should converge faster
        r_warm = em_loop(
            simple_annot_df,
            F_init=r_cold.F,
            log_L_threshold=0.001,
        )
        assert r_warm.n_iterations <= r_cold.n_iterations

    def test_warm_start_new_taxa(self):
        """Warm-start with new taxa not in F_init."""
        df1 = pd.DataFrame({
            "qseqid": ["r1", "r1"],
            "staxid": ["A", "B"],
            "bitscore": [300.0, 200.0],
        })
        r1 = em_loop(df1)

        # Add new taxon C
        df2 = pd.DataFrame({
            "qseqid": ["r1", "r1", "r2"],
            "staxid": ["A", "B", "C"],
            "bitscore": [300.0, 200.0, 400.0],
        })
        r2 = em_loop(df2, F_init=r1.F)
        assert "C" in r2.F
        assert r2.F["C"] > 0


# ------------------------------------------------------------------
# Two-pass EM (prune_min_reads)
# ------------------------------------------------------------------

def _two_taxa_with_noise() -> pd.DataFrame:
    """40 reads that hit taxon A best and also a noise taxon N a little worse,
    plus 3 reads that hit only N: N gets a few expected reads in pass one."""
    rows = []
    for i in range(40):
        rows += [{"qseqid": f"r{i}", "staxid": "A", "bitscore": 100.0},
                 {"qseqid": f"r{i}", "staxid": "N", "bitscore": 98.0}]
    for i in range(3):
        rows.append({"qseqid": f"n{i}", "staxid": "N", "bitscore": 100.0})
    return pd.DataFrame(rows)


def test_two_pass_off_by_default_keeps_noise() -> None:
    res = em_loop(_two_taxa_with_noise(), temp=1.0)
    assert "N" in res.F and res.F["N"] > 0


def test_two_pass_prunes_and_reassigns() -> None:
    res = em_loop(_two_taxa_with_noise(), temp=1.0, prune_min_reads=10)
    assert set(res.F) == {"A"}
    assert res.F["A"] == pytest.approx(1.0)
    # the 3 reads that could only be N leave the EM input
    assert res.n_queries == 40
