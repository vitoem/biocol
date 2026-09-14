"""Run the Diamond executable."""

from __future__ import annotations

import logging
import shutil
import subprocess

from biocol.exceptions import DiamondExecutionError

logger = logging.getLogger(__name__)


def run_diamond_command(command: list[str]) -> None:
    executable = command[0]
    if shutil.which(executable) is None:
        raise DiamondExecutionError(
            f"'{executable}' was not found in PATH. Install Diamond in the conda environment."
        )
    logger.debug("Diamond command: %s", " ".join(command))
    completed = subprocess.run(
        command,
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        raise DiamondExecutionError(
            f"{executable} failed (exit code {completed.returncode}): {completed.stderr.strip()}"
        )
