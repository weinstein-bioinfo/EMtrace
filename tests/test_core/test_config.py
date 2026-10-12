"""Tests for core/config.py — Pydantic configuration models."""

from pathlib import Path

import pandas as pd
import pytest

from emtrace.analysis.taxonomy_utils import is_silva_annotation
from emtrace.core.config import (
    AnnotationConfig,
    BlastnConfig,
    EMConfig,
    EMtraceConfig,
    RuntimeConfig,
    TaxonomyConfig,
    TaxonomyType,
)


class TestAnnotationConfig:
    """Test AnnotationConfig defaults and overrides."""

    def test_defaults(self):
        cfg = AnnotationConfig()
        assert cfg.blastn.bitscore == 200
        assert cfg.minimap.s == 200
        assert cfg.diamond_blastp.bitscore == 200
        assert cfg.orfipy.min_len == 300

    def test_override_blastn(self):
        cfg = AnnotationConfig(blastn=BlastnConfig(bitscore=100))
        assert cfg.blastn.bitscore == 100


class TestEMConfig:
    """Test EMConfig defaults."""

    def test_defaults(self):
        cfg = EMConfig()
        assert cfg.temp == 1.0
        assert cfg.max_iter == 1000
        assert cfg.log_L_threshold == 0.01
        assert cfg.exclude_no_species is False


class TestRuntimeConfig:
    """Test RuntimeConfig."""

    def test_defaults(self):
        cfg = RuntimeConfig()
        assert cfg.num_threads == 8


class TestEMtraceConfig:
    """Test top-level config construction."""

    def test_default_construction(self):
        cfg = EMtraceConfig()
        assert cfg.annotation.ann1 == "minimap.nt.abv"
        assert cfg.pipe == "main"
        assert isinstance(cfg.em, EMConfig)
        assert isinstance(cfg.taxonomy, TaxonomyConfig)

    def test_from_yaml(self, sample_config_yaml: Path, caplog):
        cfg = EMtraceConfig.from_yaml(sample_config_yaml)
        assert cfg.runtime.num_threads == 4
        assert cfg.annotation.blastn.bitscore == 200
        assert cfg.em.max_iter == 500
        # The fixture carries a legacy `taxonomy.selected`: loaded, but flagged.
        assert cfg.taxonomy.NCBI.db == "/tmp/taxonomy"
        assert "Ignoring taxonomy.selected" in caplog.text

    def test_to_yaml_roundtrip(
        self, tmp_dir: Path, sample_config_yaml: Path,
    ):
        cfg = EMtraceConfig.from_yaml(sample_config_yaml)
        out = tmp_dir / "out.yaml"
        cfg.to_yaml(out)
        assert out.exists()
        # to_yaml dumps model_dump, not the custom _from_raw
        # format, so re-read and check a stable field
        cfg2 = EMtraceConfig.from_yaml(out)
        assert cfg2.em.max_iter == cfg.em.max_iter

    def test_database_lookup_missing_raises(self):
        cfg = EMtraceConfig()
        with pytest.raises(Exception):
            cfg.get_host_database("nonexistent")

    def test_list_databases_empty(self):
        cfg = EMtraceConfig()
        assert cfg.list_host_databases() == []
        assert cfg.databases != {}  # defaults are always populated


class TestTaxonomyFollowsDatabase:
    """Taxonomy is not selectable: each database reports through its own."""

    def test_no_global_switch(self):
        tax = EMtraceConfig().taxonomy
        assert not hasattr(tax, "selected")
        assert not hasattr(tax, "active")
        assert tax.source(TaxonomyType.SILVA) is tax.SILVA

    @pytest.mark.parametrize("name, expected", [
        ("minimap.silva", TaxonomyType.SILVA),
        ("blastn.silva", TaxonomyType.SILVA),
        ("minimap.rrnDB", TaxonomyType.NCBI),
        ("minimap.nt.abv", TaxonomyType.NCBI),
        # Annotation tables may carry a path rather than a registry key.
        ("/db/minimap.SILVA_SSURef.fasta", TaxonomyType.SILVA),
        ("not.in.registry", TaxonomyType.NCBI),
    ])
    def test_taxonomy_for_database(self, name, expected):
        cfg = EMtraceConfig.create_default()
        assert cfg.taxonomy_for_database(name) == expected

    def test_explicit_taxonomy_beats_name(self, tmp_dir: Path):
        # A SILVA-derived database whose name does not say "silva".
        p = tmp_dir / "cfg.yaml"
        p.write_text(
            "databases:\n"
            "  minimap.ssu138:\n"
            "    search: Minimap2\n"
            "    db: /db/ssu138.fasta\n"
            "    taxonomy: silva\n"
        )
        cfg = EMtraceConfig.from_yaml(p)
        assert cfg.databases["minimap.ssu138"].taxonomy == TaxonomyType.SILVA
        assert cfg.taxonomy_for_database("minimap.ssu138") == TaxonomyType.SILVA

    def test_is_silva_annotation_uses_config(self, tmp_dir: Path):
        p = tmp_dir / "cfg.yaml"
        p.write_text(
            "databases:\n"
            "  minimap.ssu138:\n"
            "    search: Minimap2\n"
            "    db: /db/ssu138.fasta\n"
            "    taxonomy: SILVA\n"
        )
        cfg = EMtraceConfig.from_yaml(p)
        annot = pd.DataFrame({"db": ["minimap.ssu138"], "qseqid": ["r1"]})
        assert is_silva_annotation(annot, cfg)
        # Without a config only the name rule applies.
        assert not is_silva_annotation(annot)

    def test_legacy_selected_keys_warn(self, tmp_dir: Path, caplog):
        p = tmp_dir / "cfg.yaml"
        p.write_text(
            "taxonomy:\n"
            "  ncbi:\n"
            "    selected: false\n"
            "  silva:\n"
            "    selected: true\n"
            "    db: /db/silva.db\n"
        )
        cfg = EMtraceConfig.from_yaml(p)
        assert cfg.taxonomy.SILVA.db == "/db/silva.db"
        assert "taxonomy.ncbi.selected" in caplog.text
        assert "taxonomy.silva.selected" in caplog.text
