"""Entry point: parse arguments, dispatch, report errors."""

from __future__ import annotations

import logging
import sys

from biocol import BlastError, DiamondError, FastaError, HmmError, KofamError, MetadataError

from biocol.cli.commands import run_from_blast, run_from_fasta
from biocol.cli.parser import build_parser
from biocol.cli.style import (
    BOLD,
    CYAN,
    GREEN,
    RED,
    YELLOW,
    disable_color,
    enable_windows_vt,
    paint,
    reset_color,
)

_COMMANDS = {
    "run": run_from_fasta,
    "from-blast": run_from_blast,
}


class _CliFormatter(logging.Formatter):
    """One readable line per stage, colored when the terminal allows it."""

    def format(self, record: logging.LogRecord) -> str:
        message = record.getMessage()
        label = paint("biocol", BOLD, CYAN, stream=sys.stderr)
        if record.levelno >= logging.ERROR:
            body = paint(message, BOLD, RED, stream=sys.stderr)
        elif record.levelno >= logging.WARNING:
            body = paint(message, YELLOW, stream=sys.stderr)
        else:
            body = message
        return f"{label}  {body}"


def _configure_cli_logging() -> None:
    logger = logging.getLogger("biocol")
    logger.setLevel(logging.INFO)
    logger.propagate = False
    logger.handlers.clear()
    handler = logging.StreamHandler(sys.stderr)
    handler.setLevel(logging.INFO)
    handler.setFormatter(_CliFormatter())
    logger.addHandler(handler)


def main(argv: list[str] | None = None) -> int:
    argv_list = list(sys.argv[1:] if argv is None else argv)
    reset_color()
    if "--no-color" in argv_list:
        disable_color()
    enable_windows_vt()
    _configure_cli_logging()

    parser = build_parser()
    stripped = [item for item in argv_list if item != "--no-color"]
    if not stripped:
        parser.print_help()
        print("error: a command is required (run or from-blast)", file=sys.stderr)
        return 2
    args = parser.parse_args(stripped)
    try:
        output = _COMMANDS[args.command](args)
    except (FileNotFoundError, FastaError, BlastError, DiamondError, HmmError, KofamError, MetadataError, OSError, ValueError) as exc:
        print(paint(f"error: {exc}", BOLD, RED, stream=sys.stderr), file=sys.stderr)
        return 1
    print(paint("Done.", BOLD, GREEN, stream=sys.stderr), file=sys.stderr)
    print(paint(str(output), GREEN, stream=sys.stdout))
    return 0


def console() -> None:
    raise SystemExit(main())
