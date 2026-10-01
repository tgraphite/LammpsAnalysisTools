"""Error hierarchy and leveled logging to stderr.

Decision: no exit codes for success; ERROR interrupts. All diagnostics go to
stderr so stdout/CSV are never polluted.
"""

from __future__ import annotations

import sys

_LEVELS = {"info": 10, "warning": 20, "error": 30}


class LatError(Exception):
    """Fatal problem: interrupts the run."""


class LatWarning(UserWarning):
    """Non-fatal problem: surfaced in header and stderr, never alters results."""


class Logger:
    def __init__(self, level: str = "info", stream=None):
        if level not in _LEVELS:
            raise LatError(f"unknown log level {level!r}; expected one of {list(_LEVELS)}")
        self.level = level
        self.stream = stream if stream is not None else sys.stderr
        self.warnings: list[str] = []

    def _emit(self, tag: str, msg: str) -> None:
        print(f"[lat {tag}] {msg}", file=self.stream)

    def info(self, msg: str) -> None:
        if _LEVELS[self.level] <= _LEVELS["info"]:
            self._emit("INFO", msg)

    def notice(self, msg: str) -> None:
        """Always shown on stderr (used for produced file paths), never on stdout."""
        self._emit("INFO", msg)

    def warning(self, msg: str) -> None:
        self.warnings.append(msg)
        if _LEVELS[self.level] <= _LEVELS["warning"]:
            self._emit("WARNING", msg)

    def error(self, msg: str) -> None:
        self._emit("ERROR", msg)
        raise LatError(msg)


_logger = Logger()


def get_logger() -> Logger:
    return _logger


def set_logger(logger: Logger) -> None:
    global _logger
    _logger = logger
