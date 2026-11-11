"""Subprocess runner with logging, error handling, and force-kill support.

External tools are spawned via :func:`run` (list form) or :func:`run_shell`
(``shell=True`` for piped commands).  Both use :class:`subprocess.Popen`
under the hood and register the live process in a small module-level slot
so that :func:`kill_active` can terminate it from another thread (used by
the GUI's Force Stop action).

On POSIX, every child is started in its own session/process group
(``start_new_session=True``) so that signals reach the whole tree —
critical for ``run_shell`` pipelines that go through ``/bin/sh``.

On Windows, tools found only inside WSL run through ``wsl.exe``
(``run(..., wsl=True)``), and every ``run_shell`` pipeline runs in WSL's bash,
because the pipelines rely on ``cat``, ``pigz`` and POSIX redirection.
Paths embedded in a pipeline string must be quoted with :func:`shell_path`.
See :mod:`emtrace.utils.wsl`.
"""

import logging
import os
import shlex
import signal
import subprocess
import sys
import threading
import time
from pathlib import Path
from typing import Optional, Sequence

from emtrace.core.exceptions import ExternalToolError
from emtrace.utils import wsl as _wsl

logger = logging.getLogger("emtrace.utils.shell")

_IS_WINDOWS = sys.platform == "win32"

# Active-process registry.  Only one pipeline worker runs at a time, so a
# single slot suffices.  Guarded by a lock so the GUI thread can read it
# while the worker thread is mid-Popen.
_active_lock = threading.Lock()
_active_proc: Optional[subprocess.Popen] = None
# Pid file of the active process's Linux process group when it runs in WSL.
_active_wsl_pidfile: Optional[str] = None


def _register_active(proc: subprocess.Popen, wsl_pidfile: Optional[str] = None) -> None:
    global _active_proc, _active_wsl_pidfile
    with _active_lock:
        _active_proc = proc
        _active_wsl_pidfile = wsl_pidfile


def _clear_active(proc: subprocess.Popen) -> None:
    global _active_proc, _active_wsl_pidfile
    with _active_lock:
        if _active_proc is proc:
            _active_proc = None
            _active_wsl_pidfile = None


def shell_path(path: str | Path) -> str:
    """Quote *path* for a :func:`run_shell` command string.

    On Windows the pipeline runs in WSL, so the path is converted to its WSL
    form first. Also protects paths containing spaces on macOS/Linux.
    """
    p = _wsl.to_wsl_path(path) if _IS_WINDOWS else str(path)
    return shlex.quote(p)


def kill_active(grace_s: float = 2.0) -> bool:
    """Terminate the currently-running subprocess (and its children).

    Sends ``SIGTERM`` to the child's process group, waits up to
    ``grace_s`` seconds, then escalates to ``SIGKILL`` for anything still
    alive.  On Windows, sends ``CTRL_BREAK_EVENT`` then falls back to
    ``Popen.kill``; a command running in WSL has its Linux process group
    signalled first (killing ``wsl.exe`` alone may leave it running).

    Safe to call from any thread.  Returns ``True`` if a process was
    signalled, ``False`` if no subprocess is currently active.
    """
    with _active_lock:
        proc = _active_proc
        wsl_pidfile = _active_wsl_pidfile
    if proc is None or proc.poll() is not None:
        return False

    pid = proc.pid
    logger.warning("Force-killing subprocess pid=%s", pid)

    if wsl_pidfile:
        _wsl.kill_group(wsl_pidfile, grace_s)

    if _IS_WINDOWS:
        try:
            proc.send_signal(signal.CTRL_BREAK_EVENT)
        except (OSError, ValueError):
            pass
    else:
        try:
            os.killpg(os.getpgid(pid), signal.SIGTERM)
        except (ProcessLookupError, PermissionError):
            pass

    deadline = time.monotonic() + grace_s
    while time.monotonic() < deadline:
        if proc.poll() is not None:
            return True
        time.sleep(0.05)

    if _IS_WINDOWS:
        try:
            proc.kill()
        except OSError:
            pass
    else:
        try:
            os.killpg(os.getpgid(pid), signal.SIGKILL)
        except (ProcessLookupError, PermissionError):
            pass
    return True


def _popen_kwargs(extra: Optional[dict] = None) -> dict:
    kw: dict = dict(extra or {})
    if _IS_WINDOWS:
        kw["creationflags"] = (
            kw.get("creationflags", 0)
            | subprocess.CREATE_NEW_PROCESS_GROUP
            | _wsl.NO_WINDOW
        )
    else:
        kw["start_new_session"] = True
    return kw


def _run_popen(
    proc_factory,
    *,
    name: str,
    stdin_data: Optional[str],
    timeout: Optional[int],
    check: bool,
    wsl_pidfile: Optional[str] = None,
) -> subprocess.CompletedProcess:
    """Spawn a Popen via ``proc_factory``, register it, drive communicate."""
    try:
        proc = proc_factory()
    except FileNotFoundError as exc:
        raise ExternalToolError(
            f"Command not found: {exc.filename or name}",
            tool=name,
            returncode=-1,
        )

    _register_active(proc, wsl_pidfile)
    killed_by_user = False
    try:
        try:
            stdout, stderr = proc.communicate(input=stdin_data, timeout=timeout)
        except subprocess.TimeoutExpired:
            kill_active(grace_s=1.0)
            try:
                stdout, stderr = proc.communicate(timeout=2)
            except subprocess.TimeoutExpired:
                stdout = stderr = ""
            raise ExternalToolError(
                f"{name} timed out after {timeout}s",
                tool=name,
                returncode=-1,
            )
    finally:
        _clear_active(proc)
        # Returncode < 0 on POSIX = killed by signal (e.g. -SIGTERM = -15).
        # Treat that as a user force-quit so the message is helpful.
        if proc.returncode is not None and proc.returncode < 0:
            killed_by_user = True

    if check and proc.returncode != 0:
        if killed_by_user:
            raise ExternalToolError(
                f"{name} was force-killed (signal {-proc.returncode})",
                tool=name,
                returncode=proc.returncode,
                stderr=(stderr or ""),
            )
        stderr_tail = (stderr or "")[-500:]
        raise ExternalToolError(
            f"{name} failed (exit {proc.returncode}): {stderr_tail}",
            tool=name,
            returncode=proc.returncode,
            stderr=stderr or "",
        )

    return subprocess.CompletedProcess(
        args=proc.args,
        returncode=proc.returncode,
        stdout=stdout,
        stderr=stderr,
    )


def run(
    cmd: Sequence[str | Path],
    *,
    cwd: Optional[Path] = None,
    check: bool = True,
    capture: bool = True,
    stdin_data: Optional[str] = None,
    timeout: Optional[int] = None,
    tool_name: Optional[str] = None,
    wsl: bool = False,
) -> subprocess.CompletedProcess:
    """Run a command and return the CompletedProcess.

    Parameters
    ----------
    cmd : sequence of str/Path
        Command and arguments.
    cwd : Path, optional
        Working directory.
    check : bool
        Raise ExternalToolError on non-zero exit.
    capture : bool
        Capture stdout/stderr (True) or inherit terminal (False).
    stdin_data : str, optional
        Data to pipe into stdin.
    timeout : int, optional
        Timeout in seconds.
    tool_name : str, optional
        Tool name for error messages (defaults to cmd[0]).
    wsl : bool
        Run the command inside WSL (Windows only): ``cmd[0]`` is a Linux
        path and Windows paths among the arguments are converted.

    Returns
    -------
    subprocess.CompletedProcess
    """
    cmd_list = [str(c) for c in cmd]
    name = tool_name or Path(cmd_list[0]).name

    pidfile: Optional[str] = None
    popen_cwd = str(cwd) if cwd else None
    encoding: Optional[str] = None
    if wsl:
        logger.debug("Running (WSL): %s", " ".join(cmd_list))
        cmd_list, pidfile = _wsl.wrap_command(cmd_list, cwd=cwd or None)
        popen_cwd = None    # passed to wsl.exe as --cd
        encoding = "utf-8"  # Linux tools write UTF-8, not the Windows code page
    else:
        logger.debug("Running: %s", " ".join(cmd_list))

    stdout_pipe = subprocess.PIPE if capture else None
    stderr_pipe = subprocess.PIPE if capture else None
    stdin_pipe = subprocess.PIPE if stdin_data is not None else None

    def factory() -> subprocess.Popen:
        return subprocess.Popen(
            cmd_list,
            cwd=popen_cwd,
            stdin=stdin_pipe,
            stdout=stdout_pipe,
            stderr=stderr_pipe,
            text=True,
            encoding=encoding,
            errors="replace",  # a stray non-UTF-8 byte must not crash the runner
            **_popen_kwargs(),
        )

    return _run_popen(
        factory,
        name=name,
        stdin_data=stdin_data,
        timeout=timeout,
        check=check,
        wsl_pidfile=pidfile,
    )


def run_shell(
    cmd: str,
    *,
    cwd: Optional[Path] = None,
    check: bool = True,
    capture: bool = True,
    timeout: Optional[int] = None,
    tool_name: str = "shell",
) -> subprocess.CompletedProcess:
    """Run a shell command string (supports pipes, redirects).

    Use this for multi-tool pipelines like
    ``samtools sort | samtools fixmate | samtools sort``.

    Quote every path in *cmd* with :func:`shell_path`. On Windows the
    command runs in WSL's bash, which must be installed.
    """
    stdout_pipe = subprocess.PIPE if capture else None
    stderr_pipe = subprocess.PIPE if capture else None

    if _IS_WINDOWS:
        if not _wsl.is_available():
            raise ExternalToolError(
                f"{tool_name} needs WSL (Windows Subsystem for Linux), which "
                "was not found. See docs/windows_install.md.",
                tool=tool_name,
                returncode=-1,
            )
        logger.debug("Running (WSL shell): %s", cmd)
        argv, pidfile = _wsl.wrap_command([cmd], cwd=cwd or None, shell=True)

        def factory() -> subprocess.Popen:
            return subprocess.Popen(
                argv,
                stdout=stdout_pipe,
                stderr=stderr_pipe,
                text=True,
                encoding="utf-8",
                errors="replace",
                **_popen_kwargs(),
            )
    else:
        logger.debug("Running (shell): %s", cmd)
        pidfile = None

        def factory() -> subprocess.Popen:
            return subprocess.Popen(
                cmd,
                shell=True,
                cwd=str(cwd) if cwd else None,
                stdout=stdout_pipe,
                stderr=stderr_pipe,
                text=True,
                errors="replace",  # a stray non-UTF-8 byte must not crash the runner
                **_popen_kwargs(),
            )

    return _run_popen(
        factory,
        name=tool_name,
        stdin_data=None,
        timeout=timeout,
        check=check,
        wsl_pidfile=pidfile,
    )
