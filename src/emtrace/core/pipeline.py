"""Pipeline orchestrator — runs sequences of pipeline steps.

Follows the ViRefSet pattern: StepName enum, STEP_REGISTRY, predefined
WORKFLOWS, and a Pipeline class that initialises and executes steps.

Also ports the Perl ``options_format()`` logic for ``--anns``,
``--ann1``, ``--ann2``, and ``--hostgen`` into
``resolve_annotation_strategy()`` and ``resolve_host_databases()``.
"""

from __future__ import annotations

import logging
import os
import re
from pathlib import Path
from typing import Dict, List, Sequence

from emtrace.core.config import EMtraceConfig
from emtrace.core.exceptions import (
    ConfigurationError,
    DatabaseNotConfiguredError,
    StepNotFoundError,
)
from emtrace.core.pipelinejob import AnnotationRound, PipelineJob
from emtrace.steps import STEP_REGISTRY, StepName
from emtrace.steps.base import PipelineStep

logger = logging.getLogger("emtrace.core.pipeline")

# ------------------------------------------------------------------
# Predefined step sequences
# ------------------------------------------------------------------

WORKFLOWS: Dict[str, List[StepName]] = {
    "main": [
        StepName.PORE,
        StepName.CHOP,
        StepName.FLT,
        StepName.ANN,
        StepName.EM,
        StepName.REP,
        StepName.PLOT,
        StepName.PACK,
        StepName.CLEAN,
    ],
    "annotate": [
        StepName.ANN,
    ],
    "report": [
        StepName.REP,
        StepName.PLOT,
    ],
    "realtime": [
        StepName.PORE,
        StepName.CHOP,
        StepName.FLT,
        StepName.ANN,
        StepName.EM,
        StepName.REP,
    ],
    "preprocess": [
        StepName.PORE,
        StepName.CHOP,
        StepName.FLT,
    ],
}

# Maps CLI step abbreviations to StepName (Perl options_format)
_STEP_ALIASES: Dict[str, StepName] = {
    "pore": StepName.PORE,
    "porechop": StepName.PORE,
    "chop": StepName.CHOP,
    "chopper": StepName.CHOP,
    "flt": StepName.FLT,
    "filter": StepName.FLT,
    "ann": StepName.ANN,
    "annot": StepName.ANN,
    "em": StepName.EM,
    "emabund": StepName.EM,
    "ass": StepName.ASS,
    "assemble": StepName.ASS,
    "assembly": StepName.ASS,
    "rep": StepName.REP,
    "report": StepName.REP,
    "pack": StepName.PACK,
    "clean": StepName.CLEAN,
    "plot": StepName.PLOT,
    "plots": StepName.PLOT,
}


class Pipeline:
    """Instantiates registered steps and runs them in sequence.

    Example::

        config = EMtraceConfig.from_yaml("config.yaml")
        pipeline = Pipeline(config)
        sample = PipelineJob(name="s1", input_reads=Path("reads.fq.gz"))
        sample = pipeline.setup_job(sample, ann1="minimap.nt.abv")
        sample = pipeline.run_workflow(sample, "main")
    """

    def __init__(self, config: EMtraceConfig):
        self.config = config
        self._steps: Dict[StepName, PipelineStep] = {}
        self._initialize_steps()

    def _initialize_steps(self) -> None:
        """Create one instance of every registered step."""
        for name, cls in STEP_REGISTRY.items():
            self._steps[name] = cls(self.config)

    # ==============================================================
    # Step execution
    # ==============================================================

    def run(
        self,
        sample: PipelineJob,
        steps: Sequence[StepName],
        **kwargs,
    ) -> PipelineJob:
        """Execute a sequence of steps on a sample.

        Parameters
        ----------
        sample : PipelineJob
            The sample state to process.
        steps : sequence of StepName
            Steps to execute, in order.

        Returns
        -------
        PipelineJob
            The updated sample after all steps complete.
        """
        # Filter to only registered steps
        available = []
        for s in steps:
            if s in self._steps:
                available.append(s)
            else:
                logger.warning(
                    "Step '%s' not registered — skipping",
                    s.value,
                )

        logger.info(
            "Running %d steps: %s",
            len(available),
            " -> ".join(s.value for s in available),
        )
        for step_name in available:
            step = self._steps[step_name]
            sample = step.execute(sample, **kwargs)
        logger.info("Pipeline finished")
        return sample

    def run_workflow(
        self,
        sample: PipelineJob,
        workflow: str,
        **kwargs,
    ) -> PipelineJob:
        """Run a predefined workflow by name."""
        if workflow not in WORKFLOWS:
            raise StepNotFoundError(
                workflow,
                available=list(WORKFLOWS),
            )
        return self.run(
            sample, WORKFLOWS[workflow], **kwargs,
        )

    # ==============================================================
    # Parsing
    # ==============================================================

    @staticmethod
    def parse_steps(step_string: str) -> List[StepName]:
        """Parse a comma-separated step string into StepName list.

        Supports workflow names (``"main"``) and individual step
        names/aliases (``"pre,flt,ann"``).  Mirrors the Perl
        ``options_format()`` ``--pipe`` logic.

        Parameters
        ----------
        step_string : str
            e.g. ``"pre,flt,ann"`` or ``"main"``.

        Returns
        -------
        list of StepName
        """
        result: List[StepName] = []

        def _add(sn: StepName) -> None:
            if sn not in result:
                result.append(sn)

        for token in step_string.split(","):
            token = token.strip().lower()
            if not token:
                continue
            # Check workflow name
            if token in WORKFLOWS:
                for sn in WORKFLOWS[token]:
                    _add(sn)
                continue
            # Check alias
            if token in _STEP_ALIASES:
                _add(_STEP_ALIASES[token])
                continue
            # Check direct StepName value
            valid = {s.value: s for s in StepName}
            if token in valid:
                _add(valid[token])
                continue
            raise StepNotFoundError(
                token,
                available=sorted(
                    set(list(_STEP_ALIASES) + [s.value for s in StepName])
                ),
            )

        # Reorder into canonical pipeline order (StepName declaration order) so
        # combining a workflow with extra steps (e.g. "main,assembly") runs each
        # step at its correct position — assembly before reporting, not after.
        canonical = {s: i for i, s in enumerate(StepName)}
        result.sort(key=lambda s: canonical[s])
        return result

    # ==============================================================
    # PipelineJob setup — annotation strategy + bgfilter databases
    # ==============================================================

    def setup_job(
        self,
        sample: PipelineJob,
        *,
        anns: str = "",
        ann1: str = "",
        ann2: str = "",
        bgfilter: str = "",
    ) -> PipelineJob:
        """Configure annotation rounds and bgfilter databases on *sample*.

        Ports the ``options_format()`` logic for ``--anns``,
        ``--ann1``, ``--ann2``, and ``--bgfilter``.

        Parameters
        ----------
        sample : PipelineJob
            The sample to configure (mutated in-place).
        anns : str
            Named annotation strategy (e.g. ``"vi.nt"``).
            Overrides *ann1* / *ann2* with the strategy
            definition from config.
        ann1 : str
            Comma-separated annotation database keys for
            round 1 (e.g. ``"minimap.nt.abv"``).
        ann2 : str
            Comma-separated ``target:db`` pairs for round 2
            (e.g. ``"vi:blastn.nt,ab:diamondx.nr"``).
        bgfilter : str
            Comma-separated bgfilter keys or FASTA paths.
        """
        # Resolve annotation strategy → ann1/ann2
        if anns:
            ann1, ann2 = self._resolve_strategy(anns)
        elif not ann1:
            ann1 = self.config.annotation.ann1
            ann2 = ann2 or self.config.annotation.ann2

        # Parse ann1
        if ann1:
            round1 = self._parse_ann1(ann1)
            if round1.databases:
                sample.ann1 = round1

        # Parse ann2
        if ann2:
            round2 = self._parse_ann2(ann2)
            if round2.databases:
                sample.ann2 = round2

        # Resolve host databases
        bgfilter = bgfilter or self.config.bgfilter.bgfilter
        if bgfilter and bgfilter != "0":
            sample.bgfilter_databases = (
                self._resolve_host_databases(bgfilter)
            )

        return sample

    def _resolve_strategy(
        self, anns: str,
    ) -> tuple[str, str]:
        """Resolve a named strategy to (ann1, ann2) strings.

        Port of Perl lines 1294-1304.
        """
        strat_str = self.config.resolve_strategy(anns)
        ann1_match = re.search(
            r"--ann1\s+([\w.,]+)", strat_str, re.I,
        )
        ann2_match = re.search(
            r"--ann2\s+([\w.,:]+)", strat_str, re.I,
        )
        if not ann1_match:
            raise ConfigurationError(
                f"Invalid annotation strategy '{anns}'"
                f": no --ann1 in '{strat_str}'"
            )
        ann1 = ann1_match.group(1).lower()
        ann2 = (
            ann2_match.group(1).lower()
            if ann2_match else ""
        )
        return ann1, ann2

    def _parse_ann1(
        self, ann1_str: str,
    ) -> AnnotationRound:
        """Parse ``--ann1`` string into AnnotationRound.

        Port of Perl lines 1380-1406.
        Round 1: first db targets 'all' reads, subsequent
        dbs target 'un' (unannotated) reads.
        """
        databases = []
        target = "all"
        for dbname in ann1_str.split(","):
            dbname = dbname.strip().lower()
            if not dbname:
                continue
            entry = self._lookup_ann_db(dbname)
            if entry is None:
                continue
            databases.append({
                "db": entry.db,
                "name": entry.name,
                "search": entry.search,
                "target": target,
            })
            target = "un"
        if not databases and ann1_str.strip():
            raise ConfigurationError(
                f"--ann1 '{ann1_str}' resolved to no configured database. "
                f"Available: {self.config.list_databases()}"
            )
        return AnnotationRound(databases=databases)

    def _parse_ann2(
        self, ann2_str: str,
    ) -> AnnotationRound:
        """Parse ``--ann2`` string into AnnotationRound.

        Port of Perl lines 1409-1441.
        Format: ``target:dbname,...`` where target is
        vi/ph/ab/un.
        """
        valid_targets = {"vi", "ph", "ab", "un", "euk"}
        databases = []
        for token in ann2_str.split(","):
            token = token.strip().lower()
            if not token:
                continue
            if ":" in token:
                target, dbname = token.split(":", 1)
            else:
                target, dbname = "vi", token
            if target not in valid_targets:
                raise ConfigurationError(
                    f"Invalid ann2 target '{target}'"
                    f" — must be one of {valid_targets}"
                )
            entry = self._lookup_ann_db(dbname)
            if entry is None:
                continue
            databases.append({
                "db": entry.db,
                "name": entry.name,
                "search": entry.search,
                "target": target,
            })
        return AnnotationRound(databases=databases)

    def _lookup_ann_db(self, dbname: str):
        """Look up an annotation database entry by key.

        Returns ``None`` and logs a warning when the database is not
        configured, so callers can skip missing entries rather than abort.
        """
        try:
            return self.config.get_database(dbname)
        except DatabaseNotConfiguredError:
            logger.warning(
                "Annotation database '%s' not found in config — skipping. "
                "Available: %s",
                dbname, self.config.list_databases(),
            )
            return None

    def _resolve_host_databases(
        self, hostgen: str,
    ) -> List[Dict[str, str]]:
        """Resolve ``--hostgen`` to list of host database dicts.

        Port of Perl lines 1320-1365.  Supports:
          - Config keys: ``"Homo_sapiens,Mus_musculus"``
          - Direct FASTA paths: ``"/path/to/genome.fa"``
        """
        result = []
        for key in hostgen.split(","):
            key = key.strip()
            if not key or key == "0":
                continue

            # Check if it's a direct FASTA path
            if re.search(
                r"\.(fasta|faa|fna|fa)(\.gz)?$",
                key, re.I,
            ):
                name = Path(key).stem
                name = re.sub(r"\.(fasta|faa|fna|fa)$", "", name, flags=re.I)
                result.append({
                    "db": key,
                    "name": name,
                })
            else:
                # Look up in config
                try:
                    entry = self.config.get_host_database(key)
                    db_path = os.path.expandvars(
                        os.path.expanduser(entry.db)
                    )
                    result.append({
                        "db": db_path,
                        "name": entry.name or key,
                        "search": entry.search,
                    })
                except DatabaseNotConfiguredError:
                    raise ConfigurationError(
                        f"Host database '{key}' not found"
                        f" in config. Available: "
                        f"{self.config.list_host_databases()}"
                    )
        return result
