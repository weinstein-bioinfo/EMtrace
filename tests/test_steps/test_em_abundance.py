"""Tests for steps/em.py — EM abundance pipeline step."""

from pathlib import Path

import pandas as pd
import pytest

from emtrace.core.config import EMtraceConfig
from emtrace.core.pipelinejob import PipelineJob as Sample
from emtrace.models.annotation import write_annotation_tsv
from emtrace.steps.em import EMAbundanceStep


@pytest.fixture
def em_sample(
    tmp_dir: Path, sample_annotation_df: pd.DataFrame,
):
    """Sample with annotation table ready for EM."""
    s = Sample(
        name="em_test",
        output_dir=tmp_dir / "output",
        num_threads=1,
    )
    s.init_dirs()

    # Write annotation table where the step looks for it
    annot_path = s.annotation_dir / "annot1.tsv"
    write_annotation_tsv(sample_annotation_df, annot_path)
    s.set_file("annot1", annot_path)
    return s


class TestEMAbundanceStep:
    """Test the EM abundance pipeline step."""

    def test_validate_inputs(self, em_sample: Sample):
        cfg = EMtraceConfig()
        step = EMAbundanceStep(cfg)
        assert step.validate_inputs(em_sample) is True

    def test_no_annotation_is_skipped(self, tmp_dir: Path):
        # No annotation table is a valid edge case: validation passes and
        # run() skips EM without writing output.
        cfg = EMtraceConfig()
        step = EMAbundanceStep(cfg)
        s = Sample(
            name="empty",
            output_dir=tmp_dir / "empty_out",
        )
        s.init_dirs()
        assert step.validate_inputs(s) is True
        result = step.run(s)
        assert result.get_file("em_taxid_prob") is None

    def test_run_produces_output(
        self, em_sample: Sample,
    ):
        cfg = EMtraceConfig()
        step = EMAbundanceStep(cfg)
        result = step.run(em_sample)

        # Check output files
        prob_path = result.get_file("em_taxid_prob")
        assert prob_path is not None
        assert prob_path.exists()

        byiter_path = result.get_file(
            "em_taxid_prob_byiter",
        )
        assert byiter_path is not None
        assert byiter_path.exists()

    def test_output_format(self, em_sample: Sample):
        cfg = EMtraceConfig()
        step = EMAbundanceStep(cfg)
        result = step.run(em_sample)

        prob_path = result.get_file("em_taxid_prob")
        df = pd.read_csv(prob_path, sep="\t")
        assert "taxid" in df.columns
        assert "prob" in df.columns
        assert len(df) > 0
        # Probs should be positive
        assert (df["prob"] > 0).all()
