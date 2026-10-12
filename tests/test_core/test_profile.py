"""Bundled settings profiles (`--profile`) and the deep merge behind them.

`--config` replaces the bundled defaults outright; `--profile` overlays them.
The distinction matters because frozen config files (the benchmark's, for
instance) rely on replace semantics — anything they omit must keep falling
back to the model defaults, not silently pick up default.yaml.
"""

import click
import pytest
from click.testing import CliRunner

from emtrace.__main__ import cli
from emtrace.core.config import EMtraceConfig, _deep_merge
from emtrace.core.exceptions import ConfigurationError


# ------------------------------------------------------------------ #
# _deep_merge
# ------------------------------------------------------------------ #

def test_mappings_merge_key_by_key():
    assert _deep_merge({"a": {"x": 1, "y": 2}}, {"a": {"y": 9}}) == {"a": {"x": 1, "y": 9}}


def test_scalars_replace():
    assert _deep_merge({"a": 1}, {"a": 2}) == {"a": 2}


def test_lists_replace_wholesale_rather_than_append():
    # The rule that keeps parasite_genera and anns.options predictable.
    assert _deep_merge({"l": [1, 2, 3]}, {"l": [9]}) == {"l": [9]}


def test_new_keys_are_added():
    assert _deep_merge({"a": 1}, {"b": 2}) == {"a": 1, "b": 2}


def test_base_is_not_mutated():
    base = {"a": {"x": 1}}
    _deep_merge(base, {"a": {"x": 2}})
    assert base == {"a": {"x": 1}}


def test_mapping_replaces_a_scalar():
    assert _deep_merge({"a": 1}, {"a": {"x": 2}}) == {"a": {"x": 2}}


# ------------------------------------------------------------------ #
# The 16S profile
# ------------------------------------------------------------------ #

@pytest.fixture(scope="module")
def s16():
    return EMtraceConfig.from_profile("16S")


@pytest.fixture(scope="module")
def default():
    return EMtraceConfig.create_default()


def test_16S_profile_turns_the_tail_trim_on(s16, default):
    assert s16.report.tail == 0.01
    assert s16.report.tail_viruses is True
    assert default.report.tail == 0.0          # and the defaults stay off


def test_16S_profile_sets_the_benchmarked_cell(s16):
    assert s16.annotation.anns == "16S.rrnDB.dedup"
    assert s16.annotation.minimap.max_target_seqs == 25
    assert s16.em.temp == 10.0
    assert (s16.preprocess.chopper.quality,
            s16.preprocess.chopper.minlength,
            s16.preprocess.chopper.maxlength) == (10, 1200, 1800)


def test_16S_strategy_resolves_against_the_shipped_databases(s16):
    # A profile naming a strategy default.yaml does not define would fail only
    # at run time, deep inside the pipeline.
    assert s16.resolve_strategy("16S.rrnDB.dedup") == "--ann1 minimap.rrnDB.dedup"
    assert s16.get_database("minimap.rrnDB.dedup").group == "16S"


def test_profile_inherits_everything_it_does_not_name(s16, default):
    # Proves the overlay merges onto default.yaml rather than replacing it.
    assert len(s16.databases) == len(default.databases)
    assert s16.report.min_qcov == default.report.min_qcov
    assert s16.report.sort_reads == default.report.sort_reads
    assert s16.taxonomy.NCBI.db == default.taxonomy.NCBI.db
    assert s16.preprocess.bgfilter.bgfilter == default.preprocess.bgfilter.bgfilter


def test_profile_merge_is_onto_default_yaml_not_model_defaults(s16):
    # sort_reads is False in default.yaml and True in the model. A profile
    # that does not mention it must inherit the YAML value.
    from emtrace.core.config import ReportConfig
    assert ReportConfig().sort_reads is True
    assert s16.report.sort_reads is False


def test_unknown_profile_names_the_available_ones():
    with pytest.raises(ConfigurationError, match="Unknown profile 'nope'"):
        EMtraceConfig.from_profile("nope")
    try:
        EMtraceConfig.from_profile("nope")
    except ConfigurationError as exc:
        assert "16S" in str(exc)


# ------------------------------------------------------------------ #
# --config keeps replace semantics
# ------------------------------------------------------------------ #

def test_config_file_does_not_inherit_default_yaml(tmp_path):
    # The benchmark's frozen configs depend on this: a key they omit falls
    # back to the MODEL default, not to default.yaml.
    cfg_file = tmp_path / "minimal.yaml"
    cfg_file.write_text("report:\n  min_qcov: 0.4\n")
    cfg = EMtraceConfig.from_yaml(cfg_file)
    assert cfg.report.min_qcov == 0.4
    assert cfg.report.sort_reads is True       # model default, not YAML's False
    assert cfg.report.tail == 0.0


def test_config_file_without_tail_does_not_trim(tmp_path):
    # What every frozen benchmark config gets.
    cfg_file = tmp_path / "frozen.yaml"
    cfg_file.write_text("report:\n  min_qcov: 0.4\n")
    assert EMtraceConfig.from_yaml(cfg_file).report.tail == 0.0


# ------------------------------------------------------------------ #
# CLI wiring
# ------------------------------------------------------------------ #

def test_load_config_rejects_config_plus_profile(tmp_path):
    from emtrace.__main__ import _load_config

    cfg_file = tmp_path / "c.yaml"
    cfg_file.write_text("report: {}\n")
    with pytest.raises(click.UsageError, match="mutually exclusive"):
        _load_config(str(cfg_file), "16S")


def test_load_config_picks_the_right_source(tmp_path):
    from emtrace.__main__ import _load_config

    cfg_file = tmp_path / "c.yaml"
    cfg_file.write_text("report:\n  min_qcov: 0.25\n")
    assert _load_config(str(cfg_file), None).report.min_qcov == 0.25
    assert _load_config(None, "16S").report.tail == 0.01
    assert _load_config(None, None).report.tail == 0.0


def test_config_and_profile_are_mutually_exclusive(tmp_path):
    cfg_file = tmp_path / "c.yaml"
    cfg_file.write_text("report: {}\n")
    reads = tmp_path / "x.fq"
    reads.write_text("@r\nACGT\n+\nIIII\n")
    result = CliRunner().invoke(
        cli, ["run", "--config", str(cfg_file), "--profile", "16S",
              "-1", str(reads), "-s", "smoke"],
    )
    assert result.exit_code != 0
    assert "mutually exclusive" in result.output


@pytest.mark.parametrize("command", ["run", "watch"])
def test_profile_is_offered_on_run_and_watch(command):
    out = CliRunner().invoke(cli, [command, "--help"]).output
    assert "--profile" in out


def test_tail_flags_are_offered_on_run():
    out = CliRunner().invoke(cli, ["run", "--help"]).output
    assert "--tail " in out or "--tail FLOAT" in out
    assert "--tail-viruses" in out
    # The help must steer people away from the unrelated read-end trimmer.
    assert "--tailcrop" in out
