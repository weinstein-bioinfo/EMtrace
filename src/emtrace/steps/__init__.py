"""Pipeline steps registry.

Defines StepName enum, STEP_REGISTRY dict, and register_step
decorator.  Step modules are imported at the bottom to trigger
registration — keep new step imports there.
"""

from enum import Enum
from typing import Dict, Type

from emtrace.steps.base import PipelineStep


class StepName(str, Enum):
    """Available pipeline steps with CLI-friendly names."""

    PORE = "pore"        # Adapter removal (Porechop_ABI / Porechop)
    CHOP = "chop"        # Quality/length filtering (Chopper)
    FLT = "flt"          # Host filtering
    ANN = "ann"          # Annotation (rounds 1 & 2)
    EM = "em"            # EM abundance estimation
    ASS = "ass"          # Reference-based assembly
    REP = "rep"          # Reporting (includes statistics)
    PLOT = "plot"        # Abundance visualizations
    PACK = "pack"        # Pack output files
    CLEAN = "clean"      # Clean intermediate files




STEP_REGISTRY: Dict[StepName, Type[PipelineStep]] = {}


def register_step(step_name: StepName):
    """Decorator to register a PipelineStep subclass."""

    def decorator(cls: Type[PipelineStep]):
        STEP_REGISTRY[step_name] = cls
        return cls

    return decorator


# -----------------------------------------------------------------
# Import step modules to trigger registration.
# These MUST come after STEP_REGISTRY and register_step are defined.
# -----------------------------------------------------------------
# flake8: noqa: E402
from emtrace.steps.porechop import PorechopStep  # noqa: E402, F401
from emtrace.steps.chopper import ChopperStep  # noqa: E402, F401
from emtrace.steps.bgfilter import BgFilterStep  # noqa: E402, F401
from emtrace.steps.annotate import AnnotateStep  # noqa: E402, F401
from emtrace.steps.assembly import AssemblyStep  # noqa: E402, F401
from emtrace.steps.report import ReportStep  # noqa: E402, F401
from emtrace.steps.em import EMAbundanceStep  # noqa: E402, F401
from emtrace.steps.pack import PackStep, CleanStep  # noqa: E402, F401
from emtrace.steps.plots import PlotStep  # noqa: E402, F401


__all__ = [
    "StepName",
    "STEP_REGISTRY",
    "register_step",
    "PipelineStep",
    "PorechopStep",
    "ChopperStep",
    "BgFilterStep",
    "AnnotateStep",
    "EMAbundanceStep",
    "AssemblyStep",
    "ReportStep",
    "PackStep",
    "CleanStep",
    "PlotStep",
]
