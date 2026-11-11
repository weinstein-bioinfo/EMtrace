"""Tests for utils/wsl.py — WSL routing of Linux-only tools on Windows.

WSL itself is not available on the CI/macOS test machines. These tests cover
the parts that do not need it: path conversion, the wsl.exe command line, and
the bash runner/killer scripts (run here with plain bash, a new session
standing in for ``setsid``).
"""

import subprocess
import sys
import time
from pathlib import Path

import pytest

from emtrace.utils import shell, wsl


@pytest.mark.parametrize(
    "win, expected",
    [
        (r"C:\Users\me\run 1\reads.fq", "/mnt/c/Users/me/run 1/reads.fq"),
        ("D:/data/db.mmi", "/mnt/d/data/db.mmi"),
        ("C:\\", "/mnt/c/"),
        (r"\\wsl.localhost\Ubuntu\home\me\db.mmi", "/home/me/db.mmi"),
        (r"\\wsl$\Ubuntu-24.04\home\me", "/home/me"),
        (r"\\wsl.localhost\Ubuntu", "/"),
        ("/home/me/miniforge3/envs/emtrace/bin/minimap2", "/home/me/miniforge3/envs/emtrace/bin/minimap2"),
        (r"out\sample\x.bam", "out/sample/x.bam"),
    ],
)
def test_to_wsl_path(win, expected):
    assert wsl.to_wsl_path(win) == expected


@pytest.mark.parametrize(
    "arg, expected",
    [
        ("-x", "-x"),
        ("map-ont", "map-ont"),
        ("6G", "6G"),
        ("INFO/AD", "INFO/AD"),
        ("QUAL>=20 && INFO/DP>5", "QUAL>=20 && INFO/DP>5"),
        (r"C:\db\ref.fa", "/mnt/c/db/ref.fa"),
        (r"results\s1\reads.fq.gz", "results/s1/reads.fq.gz"),
    ],
)
def test_to_wsl_arg_only_rewrites_paths(arg, expected):
    assert wsl.to_wsl_arg(arg) == expected


def test_wrap_command(monkeypatch):
    monkeypatch.setenv("EMTRACE_WSL_DISTRO", "Ubuntu")
    monkeypatch.setattr(wsl, "search_path", lambda: "/home/me/miniforge3/envs/emtrace/bin")
    argv, pidfile = wsl.wrap_command(
        ["/home/me/bin/minimap2", "-t", "4", r"C:\db\ref.mmi", r"C:\run\reads.fq"],
        cwd=r"C:\run" if sys.platform == "win32" else "/tmp",
    )
    assert argv[:3] == ["wsl.exe", "-d", "Ubuntu"]
    assert argv[3] == "--cd"
    i = argv.index("-e")
    assert argv[i + 1:i + 4] == ["setsid", "-w", "bash"]
    assert pidfile.startswith("/tmp/emtrace-") and pidfile.endswith(".pid")
    assert argv[i + 7:i + 9] == [pidfile, "/home/me/miniforge3/envs/emtrace/bin"]
    assert argv[-3:] == ["4", "/mnt/c/db/ref.mmi", "/mnt/c/run/reads.fq"]


def test_wrap_command_shell_keeps_pipeline(monkeypatch):
    monkeypatch.setenv("EMTRACE_WSL_DISTRO", "")
    monkeypatch.setattr(wsl, "search_path", lambda: "")
    pipeline = "cat '/mnt/c/a b.sam' | samtools sort -o /mnt/c/out.bam"
    argv, _ = wsl.wrap_command([pipeline], shell=True)
    assert argv[0] == "wsl.exe" and argv[1] == "-e"
    assert argv[-3:] == ["bash", "-c", pipeline]


def test_wsl_helpers_are_inert_off_windows(monkeypatch):
    if wsl.IS_WINDOWS:
        pytest.skip("Windows")
    assert wsl.is_available() is False
    assert wsl.find_tool("minimap2") is None
    assert wsl.search_path() == ""
    assert wsl.locate_tool("nonexistent_tool_xyz_999") is None
    assert wsl.tool_argv("/usr/bin/samtools", False, ["--version"]) == ["/usr/bin/samtools", "--version"]


def test_shell_path_quotes_spaces():
    if wsl.IS_WINDOWS:
        assert shell.shell_path(r"C:\my data\x.fq") == "'/mnt/c/my data/x.fq'"
    else:
        assert shell.shell_path("/Users/me/my data/x.fq") == "'/Users/me/my data/x.fq'"
        assert shell.shell_path(Path("/tmp/x.fq")) == "/tmp/x.fq"


# ---------------------------------------------------------------------------
# The bash scripts that run inside WSL
# ---------------------------------------------------------------------------

posix_only = pytest.mark.skipif(sys.platform == "win32", reason="needs bash")


def _runner(tmp_path, *cmd, prefix=""):
    pidfile = tmp_path / "run.pid"
    proc = subprocess.Popen(
        ["bash", "-c", wsl._RUNNER, "emtrace", str(pidfile), prefix, *cmd],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        start_new_session=True,  # what `setsid -w` does inside WSL
    )
    return proc, pidfile


@posix_only
def test_runner_passes_exit_status_and_cleans_up(tmp_path):
    proc, pidfile = _runner(tmp_path, "bash", "-c", "echo out; echo err >&2; exit 3")
    out, err = proc.communicate(timeout=10)
    assert (proc.returncode, out, err) == (3, "out\n", "err\n")
    assert not pidfile.exists()


@posix_only
def test_runner_puts_tool_dirs_first_on_path(tmp_path):
    bindir = tmp_path / "envbin"
    bindir.mkdir()
    tool = bindir / "mytool"
    tool.write_text("#!/bin/sh\necho from-env\n")
    tool.chmod(0o755)
    proc, _ = _runner(tmp_path, "mytool", prefix=str(bindir))
    out, _ = proc.communicate(timeout=10)
    assert out == "from-env\n"


@posix_only
def test_killer_stops_the_whole_pipeline(tmp_path):
    proc, pidfile = _runner(tmp_path, "bash", "-c", "sleep 30 | sleep 30")
    deadline = time.monotonic() + 5
    while not pidfile.exists() and time.monotonic() < deadline:
        time.sleep(0.05)
    pgid = int(pidfile.read_text())
    subprocess.run(["bash", "-c", wsl._KILLER, "emtrace", str(pidfile), "10"], timeout=15, check=True)
    proc.communicate(timeout=5)
    assert proc.returncode != 0
    assert not pidfile.exists()
    # No process left in the group
    left = subprocess.run(["pgrep", "-g", str(pgid)], capture_output=True, text=True)
    assert left.stdout.strip() == ""


@posix_only
def test_killer_without_pidfile_is_a_noop(tmp_path):
    res = subprocess.run(
        ["bash", "-c", wsl._KILLER, "emtrace", str(tmp_path / "missing.pid"), "1"], timeout=10,
    )
    assert res.returncode == 0
