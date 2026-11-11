"""Run Linux-only tools through WSL (Windows Subsystem for Linux) on Windows.

minimap2, samtools, bcftools, chopper and porechop_abi have no maintained
native Windows builds. On Windows, EMtrace looks for each tool natively first
and, failing that, inside WSL, where the user installs them from Bioconda (see
docs/windows_install.md).

Tool lookup
    ``wsl -e`` starts a non-interactive shell that does not read ``.bashrc``,
    so a conda environment is never on its PATH. :func:`search_path` therefore
    probes the usual conda locations for an ``emtrace`` environment, and the
    tool paths are resolved against those directories first.

Settings (environment variable, else key in ``user_config.yaml``)
    ``EMTRACE_WSL_DISTRO`` / ``wsl_distro``
        WSL distribution to use; empty = the WSL default distribution.
    ``EMTRACE_WSL_BIN`` / ``wsl_bin_dir``
        Linux directory (or ``:``-separated directories) holding the tools,
        searched before the built-in conda locations.

Paths
    Windows paths are rewritten to their WSL mounts (``C:\\x`` -> ``/mnt/c/x``,
    ``\\\\wsl.localhost\\Ubuntu\\home\\u`` -> ``/home/u``). List-form commands
    are rewritten argument by argument; shell pipelines must quote their paths
    with :func:`emtrace.utils.shell.shell_path`.

Force Stop
    Terminating ``wsl.exe`` does not reliably stop the Linux process behind it.
    Every WSL command therefore runs in its own Linux process group whose id is
    written to a pid file; :func:`kill_group` signals that group from a second
    ``wsl.exe`` call.

All functions are no-ops (return ``None`` / ``False`` / ``""``) off Windows.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
import threading
import uuid
from pathlib import Path
from typing import Optional, Sequence

IS_WINDOWS = sys.platform == "win32"

#: ``creationflags`` that keep a console window from flashing up for each
#: child process of the windowed (``console=False``) GUI bundle.
NO_WINDOW: int = getattr(subprocess, "CREATE_NO_WINDOW", 0)

#: Name of the conda environment the install guide creates inside WSL.
ENV_NAME = "emtrace"

# Conda installations probed for ``envs/<ENV_NAME>/bin`` and ``bin``, relative
# to the WSL user's $HOME.
_CONDA_ROOTS = (
    "miniforge3",
    "mambaforge",
    "miniconda3",
    "anaconda3",
    "micromamba",
    ".local/share/mamba",
)

# Seconds to wait for a probe. WSL can take several seconds to boot its VM on
# the first call after login.
_PROBE_TIMEOUT = 60

_lock = threading.Lock()
_available: Optional[bool] = None
_search_path: Optional[str] = None
_tool_cache: dict[str, Optional[str]] = {}


# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------

def _setting(env_var: str, config_key: str) -> str:
    value = os.environ.get(env_var)
    if value is not None:
        return value.strip()
    from emtrace.utils.user_config import load_user_config

    return str(load_user_config().get(config_key) or "").strip()


def distro() -> str:
    """Configured WSL distribution, or ``""`` for the WSL default."""
    return _setting("EMTRACE_WSL_DISTRO", "wsl_distro")


def configured_bin_dirs() -> list[str]:
    """Tool directories configured by the user (Linux paths)."""
    raw = _setting("EMTRACE_WSL_BIN", "wsl_bin_dir")
    return [d for d in re.split(r"[:;]", raw) if d.strip()]


def wsl_exe() -> list[str]:
    """``wsl.exe`` plus the distribution option, if one is configured."""
    d = distro()
    return ["wsl.exe", "-d", d] if d else ["wsl.exe"]


def reset_cache() -> None:
    """Forget probe results (after the user changes the WSL settings)."""
    global _available, _search_path
    with _lock:
        _available = None
        _search_path = None
        _tool_cache.clear()


# ---------------------------------------------------------------------------
# Path conversion
# ---------------------------------------------------------------------------

_DRIVE = re.compile(r"^([A-Za-z]):(?:[\\/]|$)")
_WSL_UNC = re.compile(r"^[\\/]{2}(?:wsl\.localhost|wsl\$)[\\/][^\\/]+(.*)$", re.IGNORECASE)


def to_wsl_path(path: str | os.PathLike) -> str:
    """Convert a Windows path to the path WSL sees.

    >>> to_wsl_path(r"C:\\Users\\me\\run 1\\reads.fq")
    '/mnt/c/Users/me/run 1/reads.fq'
    >>> to_wsl_path(r"\\\\wsl.localhost\\Ubuntu\\home\\me\\db.mmi")
    '/home/me/db.mmi'
    >>> to_wsl_path("/home/me/db.mmi")
    '/home/me/db.mmi'

    Relative paths keep their meaning (WSL starts in the same working
    directory) and only have their separators flipped. Assumes the default
    ``/mnt`` automount root.
    """
    p = os.fspath(path)
    m = _WSL_UNC.match(p)
    if m:
        return m.group(1).replace("\\", "/") or "/"
    m = _DRIVE.match(p)
    if m:
        rest = p[2:].replace("\\", "/").lstrip("/")
        return f"/mnt/{m.group(1).lower()}/{rest}"
    return p.replace("\\", "/")


def to_wsl_arg(arg: str) -> str:
    """Convert one command-line argument if it is a Windows path.

    Used to rewrite list-form commands generically: drive-letter and
    ``\\\\wsl.localhost`` paths are converted; any other argument containing a
    backslash (a relative Windows path) has its separators flipped. No
    EMtrace tool option uses a literal backslash.
    """
    if _WSL_UNC.match(arg) or _DRIVE.match(arg) or "\\" in arg:
        return to_wsl_path(arg)
    return arg


# ---------------------------------------------------------------------------
# Probes
# ---------------------------------------------------------------------------

def _bash(script: str, *args: str, timeout: float = _PROBE_TIMEOUT) -> Optional[subprocess.CompletedProcess]:
    """Run a bash *script* inside WSL; ``None`` if wsl.exe failed to run."""
    try:
        return subprocess.run(
            [*wsl_exe(), "-e", "bash", "-c", script, "emtrace", *args],
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=timeout,
            creationflags=NO_WINDOW,
        )
    except (OSError, subprocess.SubprocessError):
        return None


def is_available() -> bool:
    """True if WSL is installed with a working distribution (cached)."""
    global _available
    if not IS_WINDOWS:
        return False
    with _lock:
        if _available is None:
            res = _bash("exit 0")
            _available = res is not None and res.returncode == 0
        return _available


def search_path() -> str:
    """``:``-joined WSL directories put before PATH for every WSL command.

    The configured directories first, then each conda installation's
    ``envs/emtrace/bin``, then each installation's base ``bin``; only
    directories that exist are kept. Cached.
    """
    global _search_path
    if not IS_WINDOWS:
        return ""
    with _lock:
        if _search_path is None:
            candidates = [
                *configured_bin_dirs(),
                *(f"{root}/envs/{ENV_NAME}/bin" for root in _CONDA_ROOTS),
                *(f"{root}/bin" for root in _CONDA_ROOTS),
            ]
            script = (
                'for d in "$@"; do '
                'case $d in /*) p=$d;; *) p=$HOME/$d;; esac; '
                '[ -d "$p" ] && printf "%s\\n" "$p"; '
                "done; exit 0"
            )
            res = _bash(script, *candidates)
            dirs = res.stdout.splitlines() if res is not None and res.returncode == 0 else []
            _search_path = ":".join(d for d in dirs if d.startswith("/"))
        return _search_path


def find_tool(tool: str) -> Optional[str]:
    """Absolute path of *tool* inside WSL, or ``None`` (cached)."""
    if not IS_WINDOWS:
        return None
    with _lock:
        if tool in _tool_cache:
            return _tool_cache[tool]
    prefix = search_path()
    res = _bash('PATH="$1${1:+:}$PATH"; command -v -- "$2"', prefix, tool)
    path = None
    if res is not None and res.returncode == 0:
        out = res.stdout.strip().splitlines()
        if out and out[0].startswith("/"):
            path = out[0]
    with _lock:
        _tool_cache[tool] = path
    return path


def is_executable(path: str) -> bool:
    """True if the Linux *path* is an executable file inside WSL."""
    if not IS_WINDOWS:
        return False
    res = _bash('[ -f "$1" ] && [ -x "$1" ]', path)
    return res is not None and res.returncode == 0


def locate_tool(tool: str) -> Optional[tuple[str, bool]]:
    """Find *tool*: ``(path, via_wsl)``, or ``None`` if it is not installed.

    The native PATH wins; on Windows WSL is tried next.
    """
    native = shutil.which(tool)
    if native:
        return native, False
    if IS_WINDOWS:
        wsl_path = find_tool(tool)
        if wsl_path:
            return wsl_path, True
    return None


# ---------------------------------------------------------------------------
# Running commands
# ---------------------------------------------------------------------------

# Runs inside ``setsid -w``, so this bash leads a new process group whose id
# is $$. Records it for kill_group(), puts the tool directories on PATH, runs
# the command, and passes its exit status through.
_RUNNER = (
    'f=$1; p=$2; shift 2; '
    'printf "%s" "$$" > "$f"; '
    'if [ -n "$p" ]; then PATH="$p:$PATH"; export PATH; fi; '
    '"$@"; rc=$?; rm -f "$f"; exit $rc'
)

_KILLER = (
    'p=$(cat "$1" 2>/dev/null) || exit 0; '
    'kill -TERM -- "-$p" 2>/dev/null; '
    'i=0; while [ "$i" -lt "$2" ] && kill -0 -- "-$p" 2>/dev/null; '
    'do sleep 0.1; i=$((i+1)); done; '
    'kill -KILL -- "-$p" 2>/dev/null; rm -f "$1"; exit 0'
)


def wrap_command(
    cmd: Sequence[str],
    *,
    cwd: Optional[str | os.PathLike] = None,
    shell: bool = False,
) -> tuple[list[str], str]:
    """Build the ``wsl.exe`` command line that runs *cmd* inside WSL.

    *cmd* is an argument list (Windows paths in it are converted), or with
    ``shell=True`` a one-element list holding a bash command string whose
    paths are already WSL paths.

    Returns ``(argv, pidfile)``; *pidfile* is the Linux path that will hold
    the command's process-group id while it runs (see :func:`kill_group`).
    """
    pidfile = f"/tmp/emtrace-{uuid.uuid4().hex}.pid"
    if shell:
        inner = ["bash", "-c", cmd[0]]
    else:
        inner = [to_wsl_arg(str(a)) for a in cmd]
    argv = wsl_exe()
    if cwd is not None:
        argv += ["--cd", to_wsl_path(Path(cwd).absolute())]
    argv += ["-e", "setsid", "-w", "bash", "-c", _RUNNER, "emtrace",
             pidfile, search_path(), *inner]
    return argv, pidfile


def kill_group(pidfile: str, grace_s: float = 2.0) -> None:
    """SIGTERM the WSL process group recorded in *pidfile*, then SIGKILL.

    Blocks for at most *grace_s* (plus WSL start-up) seconds.
    """
    ticks = max(1, int(grace_s * 10))
    _bash(_KILLER, pidfile, str(ticks), timeout=grace_s + 30)


def tool_argv(path: str, via_wsl: bool, args: Sequence[str]) -> list[str]:
    """Command line for a one-off tool call outside the pipeline runner.

    For short probes (version checks, GUI helpers) that must not take over
    the Force Stop slot in :mod:`emtrace.utils.shell`.
    """
    if via_wsl:
        argv, _ = wrap_command([path, *args])
        return argv
    return [path, *args]
