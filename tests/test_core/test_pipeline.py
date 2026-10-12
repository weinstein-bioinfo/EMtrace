"""Tests for core/pipeline.py — Pipeline orchestrator."""

import pytest

from emtrace.core.config import EMtraceConfig
from emtrace.core.exceptions import StepNotFoundError
from emtrace.core.pipeline import (
    Pipeline,
    WORKFLOWS,
)
from emtrace.steps import StepName


class TestParseSteps:
    """Test Pipeline.parse_steps()."""

    def test_workflow_name(self):
        steps = Pipeline.parse_steps("main")
        assert StepName.PORE in steps
        assert StepName.CHOP in steps
        assert StepName.ANN in steps
        assert StepName.EM in steps

    def test_individual_steps(self):
        steps = Pipeline.parse_steps("pore,flt,ann")
        assert steps == [
            StepName.PORE, StepName.FLT, StepName.ANN,
        ]

    def test_aliases(self):
        steps = Pipeline.parse_steps("porechop,filter,annot")
        assert steps == [
            StepName.PORE, StepName.FLT, StepName.ANN,
        ]

    def test_unknown_raises(self):
        with pytest.raises(StepNotFoundError):
            Pipeline.parse_steps("nonexistent")

    def test_deduplication(self):
        steps = Pipeline.parse_steps("pore,pore,ann")
        assert steps.count(StepName.PORE) == 1

    def test_realtime_workflow(self):
        steps = Pipeline.parse_steps("realtime")
        assert StepName.PORE in steps
        assert StepName.ANN in steps
        assert StepName.EM in steps
        assert StepName.PACK not in steps


class TestWorkflows:
    """Test WORKFLOWS dict."""

    def test_main_workflow_exists(self):
        assert "main" in WORKFLOWS

    def test_realtime_workflow_exists(self):
        assert "realtime" in WORKFLOWS

    def test_all_workflows_have_steps(self):
        for name, steps in WORKFLOWS.items():
            assert len(steps) > 0, f"Workflow '{name}' is empty"


class TestPipelineInit:
    """Test Pipeline initialization."""

    def test_creates_step_instances(self):
        cfg = EMtraceConfig()
        pipe = Pipeline(cfg)
        # Should have at least the registered steps
        assert StepName.PORE in pipe._steps
        assert StepName.CHOP in pipe._steps
        assert StepName.ANN in pipe._steps
        assert StepName.EM in pipe._steps

    def test_run_workflow_unknown_raises(self):
        cfg = EMtraceConfig()
        pipe = Pipeline(cfg)
        from emtrace.core.pipelinejob import PipelineJob as Sample
        s = Sample(name="test")
        with pytest.raises(StepNotFoundError):
            pipe.run_workflow(s, "nonexistent_workflow")
