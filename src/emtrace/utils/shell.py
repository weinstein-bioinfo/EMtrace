"""Subprocess runner with logging, error handling, and force-kill support.

External tools are spawned via :func:`run` (list form) or :func:`run_shell`
(``shell=True`` for piped commands).  Both use :class:`subprocess.Popen`
under the hood and register the live process in a small module-level slot
so that :func:`kill_active` can terminate it from another thread (used by
the GUI's Force Stop action).

On POSIX, every child is started in its own session/process group
(``start_new_session=True``) so that signals reach the whole tree —
critical for ``run_shell`` pipelines that go through ``/bin/sh``.
"""

import logging
import os
import signal
import subprocess
import sys
import threading
import time
from pathlib import Path
from typing import Optional, Sequence

from emtrace.core.exceptions import ExternalToolError

logger = logging.getLogger("emtrace.utils.shell")

_IS_WINDOWS = sys.platform == "win32"

# Active-process registry.  Only one pipeline worker runs at a time, so a
# single slot suffices.  Guarded by a lock so the GUI thread can read it
# while the worker thread is mid-Popen.
_active_lock = threading.Lock()
_active_proc: Optional[subprocess.Popen] = None


def _register_active(proc: subprocess.Popen) -> None:
    global _active_proc
    with _active_lock:
        _active_proc = proc


def _clear_active(proc: subprocess.Popen) -> None:
    global _active_proc
    with _active_lock:
        if _active_proc is proc:
            _active_proc = None


def kill_active(grace_s: float = 2.0) -> bool:
    """Terminate the currently-running subprocess (and its children).

    Sends ``SIGTERM`` to the child's process group, waits up to
    ``grace_s`` seconds, then escalates to ``SIGKILL`` for anything still
    alive.  On Windows, sends ``CTRL_BREAK_EVENT`` then falls back to
    ``Popen.kill``.

    Safe to call from any thread.  Returns ``True`` if a process was
    signalled, ``False`` if no subprocess is currently active.
    """
    with _active_lock:
        proc = _active_proc
    if proc is None or proc.poll() is not None:
        return False

    pid = proc.pid
    logger.warning("Force-killing subprocess pid=%s", pid)

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
            kw.get("creationflags", 0) | subprocess.CREATE_NEW_PROCESS_GROUP
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

    _register_active(proc)
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

    Returns
    -------
    subprocess.CompletedProcess
    """
    cmd_list = [str(c) for c in cmd]
    name = tool_name or Path(cmd_list[0]).name

    logger.debug("Running: %s", " ".join(cmd_list))

    stdout_pipe = subprocess.PIPE if capture else None
    stderr_pipe = subprocess.PIPE if capture else None
    stdin_pipe = subprocess.PIPE if stdin_data is not None else None

    def factory() -> subprocess.Popen:
        return subprocess.Popen(
            cmd_list,
            cwd=str(cwd) if cwd else None,
            stdin=stdin_pipe,
            stdout=stdout_pipe,
            stderr=stderr_pipe,
            text=True,
            errors="replace",  # a stray non-UTF-8 byte must not crash the runner
            **_popen_kwargs(),
        )

    return _run_popen(
        factory,
        name=name,
        stdin_data=stdin_data,
        timeout=timeout,
        check=check,
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
    """
    logger.debug("Running (shell): %s", cmd)

    stdout_pipe = subprocess.PIPE if capture else None
    stderr_pipe = subprocess.PIPE if capture else None

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
    )
