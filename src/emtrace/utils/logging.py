"""Logging setup with colored console output and optional file logging."""

import logging
from contextlib import contextmanager
from pathlib import Path
from typing import Optional

try:
    import colorlog

    HAS_COLORLOG = True
except ImportError:
    HAS_COLORLOG = False

_LOG_COLORS = {
    "DEBUG": "cyan",
    "INFO": "green",
    "WARNING": "yellow",
    "ERROR": "red",
    "CRITICAL": "bold_red",
}

DEFAULT_FORMAT = "%(asctime)s | %(name)-24s | %(levelname)-8s | %(message)s"
DEFAULT_DATE_FORMAT = "%Y-%m-%d %H:%M:%S"


def setup_logging(
    level: int = logging.INFO,
    log_file: Optional[Path] = None,
    log_file_level: Optional[int] = None,
    colored: bool = True,
    format_string: Optional[str] = None,
    date_format: Optional[str] = None,
) -> logging.Logger:
    """Configure the root 'emtrace' logger.

    Parameters
    ----------
    level : int
        Console log level.
    log_file : Path, optional
        If given, also log to this file.
    log_file_level : int, optional
        File log level (defaults to *level*).
    colored : bool
        Use colorlog for console output when available.
    format_string : str, optional
        Custom format string.
    date_format : str, optional
        Custom date format.

    Returns
    -------
    logging.Logger
        The configured 'emtrace' logger.
    """
    logger = logging.getLogger("emtrace")
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()

    fmt = format_string or DEFAULT_FORMAT
    datefmt = date_format or DEFAULT_DATE_FORMAT

    # Console handler
    console = logging.StreamHandler()
    console.setLevel(level)
    if colored and HAS_COLORLOG:
        color_fmt = fmt.replace("%(levelname)", "%(log_color)s%(levelname)")
        formatter = colorlog.ColoredFormatter(
            color_fmt, datefmt=datefmt, log_colors=_LOG_COLORS
        )
    else:
        formatter = logging.Formatter(fmt, datefmt=datefmt)
    console.setFormatter(formatter)
    logger.addHandler(console)

    # File handler
    if log_file:
        log_file = Path(log_file)
        log_file.parent.mkdir(parents=True, exist_ok=True)
        fh = logging.FileHandler(log_file)
        fh.setLevel(log_file_level or level)
        fh.setFormatter(logging.Formatter(fmt, datefmt=datefmt))
        logger.addHandler(fh)

    logger.propagate = False
    return logger


def get_logger(name: str) -> logging.Logger:
    """Get a child logger under the 'emtrace' namespace."""
    if name.startswith("emtrace"):
        return logging.getLogger(name)
    return logging.getLogger(f"emtrace.{name}")


def set_level(level: int) -> None:
    """Change the console handler log level."""
    logger = logging.getLogger("emtrace")
    for handler in logger.handlers:
        if isinstance(handler, logging.StreamHandler) and not isinstance(
            handler, logging.FileHandler
        ):
            handler.setLevel(level)


@contextmanager
def log_level(level: int):
    """Temporarily change console log level."""
    logger = logging.getLogger("emtrace")
    previous = []
    for handler in logger.handlers:
        if isinstance(handler, logging.StreamHandler) and not isinstance(
            handler, logging.FileHandler
        ):
            previous.append((handler, handler.level))
            handler.setLevel(level)
    try:
        yield
    finally:
        for handler, prev_level in previous:
            handler.setLevel(prev_level)


def quiet_external_loggers(level: int = logging.WARNING) -> None:
    """Suppress noisy third-party loggers."""
    for name in ("urllib3", "Bio", "matplotlib", "PIL"):
        logging.getLogger(name).setLevel(level)
