"""ANSI colors for the CLI. No extra dependencies.

Honors NO_COLOR, --no-color, and non-TTY streams.
"""

from __future__ import annotations

import os
import sys

RESET = "\033[0m"
BOLD = "\033[1m"
DIM = "\033[2m"
RED = "\033[31m"
GREEN = "\033[32m"
YELLOW = "\033[33m"
CYAN = "\033[36m"
MAGENTA = "\033[35m"

_disabled = False


def reset_color() -> None:
    global _disabled
    _disabled = False


def disable_color() -> None:
    global _disabled
    _disabled = True


def use_color(stream=None) -> bool:
    if _disabled:
        return False
    if os.environ.get("NO_COLOR"):
        return False
    if os.environ.get("TERM") == "dumb":
        return False
    target = stream if stream is not None else sys.stderr
    return bool(getattr(target, "isatty", lambda: False)())


def paint(text: str, *codes: str, stream=None) -> str:
    if not codes or not use_color(stream):
        return text
    return f"{''.join(codes)}{text}{RESET}"


def enable_windows_vt() -> None:
    """Turn on virtual terminal processing so ANSI works in Windows consoles."""
    if sys.platform != "win32":
        return
    try:
        import ctypes

        handle = ctypes.windll.kernel32.GetStdHandle(-11)
        mode = ctypes.c_ulong()
        if ctypes.windll.kernel32.GetConsoleMode(handle, ctypes.byref(mode)):
            ctypes.windll.kernel32.SetConsoleMode(handle, mode.value | 0x0004)
    except Exception:
        return
